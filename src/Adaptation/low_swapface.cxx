//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_swap.hxx"
#include "msh_swap.hxx" // for swapOptions

#include "../Mesh/Mesh.hxx"

#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../low_topo.hxx"
#include "../low_geo/measure.hxx"
#include "../low_geo/misc.hxx"
#include "../utils/aux_misc.hxx"
#include "../low_geo/lenedg.hxx"
#include "../utils/mprintf.hxx"
#include "../low_geo/normal.hxx"
#include "../io_libmeshb.hxx"
#include "../quality/low_metqua.hxx"
#include "../linalg/det.hxx"

#include "../BezierOffsets/low_gaps.hxx"

#include "../msh_checktopo.hxx"

namespace Metris{



// Swap edge between two triangles (including surface w/ tets)
// Return 0 if nothing done, 1 if error, -1 if swap done
// Compute using norm specified in opt: if 0, take max.
// If norm is -1, use edge length instead.
template<class MFT, int gdim, int ideg>
int swapface(Mesh<MFT>& msh, int iface, swapOptions opt,
             MshCavity &cav, CavWrkArrs &work,
             double *qnrm0_, double *qnrm1_,
             StepDistanceObjectiveState *globalObjective,
             int ithread){

  #ifdef STEPDISTANCE
  constexpr QuaFun iquaf = QuaFun::StepDistance;
  #else
  constexpr QuaFun iquaf = QuaFun::SizeShape;
  #endif

  #if 0
  bool istop = false;
  static int nwarn = 0;
  if(nwarn++ < 10)printf("## DEBUG REMOVE THIS \n");
  int iverb0 = msh.param->iverb;
  int ivdepth0 = msh.param->ivdepth;
  if(msh.fac2poi(iface,0) == 2148 && msh.fac2poi(iface,1) == 2152 && msh.fac2poi(iface,2) == 2139
  || msh.fac2poi(iface,0) == 2148 && msh.fac2poi(iface,2) == 2152 && msh.fac2poi(iface,1) == 2139
  || msh.fac2poi(iface,1) == 2148 && msh.fac2poi(iface,2) == 2152 && msh.fac2poi(iface,0) == 2139
  || msh.fac2poi(iface,1) == 2148 && msh.fac2poi(iface,0) == 2152 && msh.fac2poi(iface,2) == 2139
  || msh.fac2poi(iface,2) == 2148 && msh.fac2poi(iface,0) == 2152 && msh.fac2poi(iface,1) == 2139
  || msh.fac2poi(iface,2) == 2148 && msh.fac2poi(iface,1) == 2152 && msh.fac2poi(iface,0) == 2139
  //
  || msh.fac2poi(iface,0) == 2137 && msh.fac2poi(iface,1) == 2152 && msh.fac2poi(iface,2) == 2139
  || msh.fac2poi(iface,0) == 2137 && msh.fac2poi(iface,2) == 2152 && msh.fac2poi(iface,1) == 2139
  || msh.fac2poi(iface,1) == 2137 && msh.fac2poi(iface,2) == 2152 && msh.fac2poi(iface,0) == 2139
  || msh.fac2poi(iface,1) == 2137 && msh.fac2poi(iface,0) == 2152 && msh.fac2poi(iface,2) == 2139
  || msh.fac2poi(iface,2) == 2137 && msh.fac2poi(iface,0) == 2152 && msh.fac2poi(iface,1) == 2139
  || msh.fac2poi(iface,2) == 2137 && msh.fac2poi(iface,1) == 2152 && msh.fac2poi(iface,0) == 2139
     ){
    printf("\n\n\n## DEBUG SET MAX PRINTS iface = {} vertices {} {} {}\n",iface,
      msh.fac2poi(iface,0),msh.fac2poi(iface,1),msh.fac2poi(iface,2));
    msh.param->iverb = 5;
    msh.param->ivdepth = 5;
    istop = true;
    writeMesh("debug_swapface0", msh);
  }
  #endif

  INCVDEPTH(msh.param);

  double &qnrm0 = *qnrm0_;
  double &qnrm1 = *qnrm1_;
  #ifdef TESTQUALITYALGO
  const int spnorm = 0;
  #else
  const int spnorm = opt.swap_norm; // swap norm -> < 0 is length, >= 0 is Lp over the 2
  #endif

  if(spnorm > 0) METRIS_ASSERT(spnorm == msh.param->opt_pnorm);


  constexpr int tdim = 2;
  constexpr AsDeg asdmet = AsDeg::P1;

  if(isdeadent(iface,msh.fac2poi)) return 0;

  CavOprOpt opts;
  CavOprInfo info;
  opts.allow_topological_correction = false;
  opts.skip_topo_checks = false;
  opts.allow_remove_points = false;
  opts.allow_remove_points_superdim = false;
  opts.dryrun = false;
  cav.lcedg.set_n(0);
  cav.lcfac.set_n(0);
  cav.lctet.set_n(0);
  cav.lcedg.allocate(1);
  cav.lcfac.allocate(2);
  cav.inewp = 0;
  if(msh.get_tdim() >= 3) cav.lctet.allocate(10);

  int iele0 = -1;
  if(msh.get_tdim() == 3){
    iele0 = msh.fac2tet(iface,0);
    if(iele0 < 0) iele0 = msh.fac2tet(iface,1);
    METRIS_ASSERT(iele0 >= 0);
  }


  double quae1; // quality of current face iface

  if(spnorm >= 0){
    #ifdef TESTQUALITYALGO
    if (msh.get_tdim() == 2){
      quae1 = metqua<MFT,gdim,tdim,iquaf>(msh,AsDeg::P1,asdmet,iface,1.0);
    }
    else{
      quae1 = 1e30;
    }
    #else
    quae1 = metqua<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,iface,1.0);
    METRIS_ASSERT_MSG(quae1 > -1.0e-16, "Negative quae1 {:e} iface {}",quae1,iface);
    #endif
  }

  CPRINTF1("-- START swapface iface = {} verts {} {} {}",iface
    ,msh.fac2poi(iface,0),msh.fac2poi(iface,1),msh.fac2poi(iface,2));
  if(spnorm >= 0 && DOPRINTS1()){
    PRINTF(" initial quality = {}\n",quae1);
  }else if(DOPRINTS1()){
    PRINTF("\n");
  }

  int iref = msh.fac2ref[iface];

  // Old qualities associated to each possible swap.
  // If pnorm >= 0, this is the p-norm of quality accross
  // Otherwise it is the length of the edge.
  double quaol[3], norfac[3], eval[18], norCAD[3][3];

  // If using regular quality (spnorm >= 0), the nordev is already taken into account
  // Otherwise, we add it here. Reuse normal of this face.
  if(gdim >= 3 && spnorm < 0){
    getnorfacP1(msh.fac2poi[iface], msh.coord, norfac);
    if(normalize_vec<3>(norfac)){
      CPRINTF1(" # face {} normal {} {} {} vanishes\n",iface,norfac[0],norfac[1],norfac[2]);
      return 1;
    }
  }
  for(int ied = 0; ied < 3; ied++){
    quaol[ied] = -1; // In bounds 0, 1, -1 is disregarded

    int ifac2 = msh.fac2fac(iface,ied);
    if(ifac2 < 0) continue; // Can't swap across nm edge or bdry

    // Different refs, skip
    if(iref != msh.fac2ref[ifac2]) continue;

    // Note: manifold but edge in-between is ineig >= 0
    int iedge = msh.fac2edg(iface, ied);
    if(iedge >= 0) continue;


    if(spnorm >= 0){
      #ifdef TESTQUALITYALGO
      if (msh.get_tdim() == 2){
        quaol[ied] = metqua<MFT,gdim,tdim,iquaf>(msh,AsDeg::P1,asdmet,ifac2,1.0);
      }
      else{
        quaol[ied] = 1e30;
      }
      #else
      quaol[ied] = metqua<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,ifac2,1.0);
      #endif
    }else{
      double sz[2], len;

      if constexpr(gdim >= 3){
        // In surface case, compute length only in tangent plane.
        // This is necessary e.g. in curved boundary layer cases.
        // One case (Sandia bump 1e-6 BL spacing) had:
        // - 2000 edge length using getlenedg_geosz
        // - 0.86 using getlenedg_geosz_plane
        len = getlenedg_geosz_plane<MFT,gdim,ideg>(msh, iface, 2, ied, sz);
      }else{
        len = getlenedg_geosz<MFT,gdim,ideg>(msh, iface, 2, ied, sz);
      }

      // Attribute quality between 0 and 1, multiplicatively symmetric:
      // i.e. q(sqrt2) = q(1/sqrt2).
      quaol[ied] = len < 1.0 ? 1.0 - len
                             : 1.0 - 1.0 / len;


      CPRINTF1(" - ied {} faces {} {} len score {}\n",ied,iface,ifac2,
               quaol[ied]);

      // To this, we want to add a normal deviation term in 3D
      if(gdim >= 3 && msh.CAD()){
        // Get normals of the two faces
        double norfa2[3];
        getnorfacP1(msh.fac2poi[ifac2], msh.coord, norfa2);
        if(normalize_vec<3>(norfa2)){
          CPRINTF1(" # face {} normal {} {} {} vanishes\n",ifac2,norfa2[0],norfa2[1],norfa2[2]);
          return 1;
        }

        // Get the cad normal at the mid edge
        double ibpo1 = msh.poi2ebp(msh.fac2poi(iface,lnoed2[ied][0]), 2, iface, -1);
        double ibpo2 = msh.poi2ebp(msh.fac2poi(iface,lnoed2[ied][1]), 2, iface, -1);
        double uv[2];
        for(int ii = 0; ii < 2; ii++)
          uv[ii] = (msh.bpo2rbi(ibpo1,ii) + msh.bpo2rbi(ibpo2,ii)) / 2.0;
        ego obj = msh.CAD.cad2fac[iref];
        int ierro = EG_evaluate(obj, uv, eval);
        if(ierro > 0) {
          CPRINTF1(" # EG_evaluate failed for face {} at uv {} {} ierro = {}\n",iref,uv[0],uv[1],ierro);
          quaol[ied] = -1;
          continue;
        }

        vecprod(&eval[3], &eval[6], norCAD[ied]);
        if(normalize_vec<3>(norCAD[ied])){
          CPRINTF1(" # CAD normal {} {} {} vanishes\n",norCAD[ied][0],norCAD[ied][1],norCAD[ied][2]);
          return 1;
        }

        double dtpr1 = abs(getprdl2<3>(norfac, norCAD[ied]));
        double dtpr2 = abs(getprdl2<3>(norfa2, norCAD[ied]));
        double qndev = 1 - MIN(dtpr1, dtpr2);

        CPRINTF1(" - ied {} faces {} {} nordevs {} {} nordev qual = {}\n",ied,iface,ifac2,
                 1-dtpr1,1-dtpr2,qndev);

        quaol[ied] = msh.param->qua_surf_wt_quality*quaol[ied]
                   + msh.param->qua_surf_wt_normal*qndev;
      }

      CPRINTF1(" - edge {} length {} quality {:15.7e}\n",ied,len, quaol[ied]);

    }
    CPRINTF1(" - candidate neighbour {} has qual {}\n",ifac2,quaol[ied]);
  }


  int idx[3] = {0,1,2};
  sortupto8_dec<double>(quaol,idx,3);

  // Simulate swaps as P1
  // Improve when curvature added to cavity
  //intAr2 fac2pol(2,3);
  int edg2pol[2]; // only for length based
  //fac2pol.set_n(2);
  cav.lcfac.set_n(2);
  cav.lcfac[0] = iface;
  for(int iix = 0; iix < 3; iix++){
    double acceptedOldNumerator = 0.;
    double acceptedNewNumerator = 0.;
    double acceptedOldTargetWeight = 0.;
    double acceptedNewTargetWeight = 0.;
    int acceptedOldCount = 0;
    int acceptedNewCount = 0;
    int ied = idx[iix];
    double quae2 = quaol[ied];
    if(quae2 < 0) continue;

    if(spnorm >= 0){
      CPRINTF1(" - consider swap qface = {} qneigh = {} \n", quae1,quae2);
    }else{
      CPRINTF1(" - consider swap lenqua init = {}\n",quae2);
    }

    // Quality of previous configuration
    #ifdef TESTQUALITYALGO
    const double qsum0 = quae1 + quae2; // current config
    const double qmax0 = MAX(quae1,quae2);
    #endif
    if(spnorm == 0){
      qnrm0 = MAX(quae1,quae2);
    }else if (spnorm > 0){
      qnrm0 = pow(pow(quae1,spnorm) + pow(quae2,spnorm), 1.0/spnorm);
    }else{
      qnrm0 = quae2;
    }
    int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
    int ip2 = msh.fac2poi(iface,lnoed2[ied][1]);

    int ifac2 = msh.fac2fac(iface,ied);
    METRIS_ASSERT(ifac2 >= 0);
    int ie2 = getedgfac(msh,ifac2,ip1,ip2);

    int nfac0 = msh.nface;
    msh.set_nface(msh.nface+2);

    msh.fac2ref[nfac0+0] = iref;
    msh.fac2ref[nfac0+1] = iref;

    msh.fac2poi(nfac0+0,0) = msh.fac2poi(iface,ied);
    msh.fac2poi(nfac0+0,1) = ip1;
    msh.fac2poi(nfac0+0,2) = msh.fac2poi(ifac2,ie2);

    msh.fac2poi(nfac0+1,0) = msh.fac2poi(iface,ied);
    msh.fac2poi(nfac0+1,1) = msh.fac2poi(ifac2,ie2);
    msh.fac2poi(nfac0+1,2) = ip2;

    bool iflat = false;
    for(int ii = 0; ii <= 1; ii++){
      iflat = !isvalideltP1<gdim,2>(msh,nfac0+ii);
      if(iflat){
        CPRINTF1(" - new face {}: {} {} {} would be flat",ii+1,msh.fac2poi(nfac0+ii,0),
          msh.fac2poi(nfac0+ii,1),msh.fac2poi(nfac0+ii,2));
        break;
      }
    }

    if(iflat){
      msh.set_nface(nfac0);
      continue;
    }

    // Need to create new ibpoi entries as well for quality function nordev.
    if constexpr (gdim >= 3){
      int ibpon, ibpoo;
      // ip1
      ibpoo = msh.poi2ebp(ip1, 2, iface, iref);
      ibpon = msh.newbpotopo(Vertex{ip1}, 2, nfac0+0);
      for(int ii = 0; ii < nrbi; ii++)
        msh.bpo2rbi(ibpon, ii) = msh.bpo2rbi(ibpoo, ii);
      // ip2
      ibpoo = msh.poi2ebp(ip2, 2, iface, iref);
      ibpon = msh.newbpotopo(Vertex{ip2}, 2, nfac0+1);
      for(int ii = 0; ii < nrbi; ii++)
        msh.bpo2rbi(ibpon, ii) = msh.bpo2rbi(ibpoo, ii);
      // ip3 (msh.fac2poi(iface,ied))
      ibpoo = msh.poi2ebp(msh.fac2poi(iface,ied), 2, iface, iref);
      ibpon = msh.newbpotopo(Vertex{msh.fac2poi(iface,ied)}, 2, nfac0+1);
      for(int ii = 0; ii < nrbi; ii++)
        msh.bpo2rbi(ibpon, ii) = msh.bpo2rbi(ibpoo, ii);
      // ip4 (msh.fac2poi(ifac2,ie2))
      ibpoo = msh.poi2ebp(msh.fac2poi(ifac2,ie2), 2, ifac2, iref);
      METRIS_ASSERT(ibpoo >= 0);
      ibpon = msh.newbpotopo(Vertex{msh.fac2poi(ifac2,ie2)}, 2, nfac0+1);
      METRIS_ASSERT(ibpon >= 0);
      for(int ii = 0; ii < nrbi; ii++)
        msh.bpo2rbi(ibpon, ii) = msh.bpo2rbi(ibpoo, ii);

      // For rembpotags.
      msh.tag[ithread]++;
      msh.fac2tag(ithread,nfac0+0) = msh.tag[ithread];
      msh.fac2tag(ithread,nfac0+1) = msh.tag[ithread];
    }

    bool skipswap = false;
    double qunw[2];
    #ifdef TESTQUALITYALGO
    double qsum1 = 0;
    double qmax1 = -1;
    #endif
    if(spnorm >= 0){

      qnrm1 = spnorm == 0 ? -1 : 0;
      for(int ifanw = nfac0+0; ifanw < nfac0+2; ifanw++){
        //qunw1 = metqua0<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,fac2pol[0],1.0);
        #ifdef TESTQUALITYALGO
        if (msh.get_tdim() == 2){
          qunw[ifanw-nfac0] = metqua<MFT,gdim,tdim,iquaf>(msh,AsDeg::P1,asdmet,ifanw,1.0);
        }
        else{
          qunw[ifanw-nfac0] = 1e30;
        }
        qsum1 += qunw[ifanw-nfac0];
        if (qunw[ifanw-nfac0] > qmax1) qmax1 = qunw[ifanw-nfac0];
        #else
        qunw[ifanw-nfac0] = metqua<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,ifanw,1.0);
        #endif
        CPRINTF1(" - new face {} quality = {} \n",ifanw-nfac0,qunw[ifanw-nfac0]);
        // Can skip already if using max
        #ifndef TESTQUALITYALGO
        if(spnorm == 0 && qunw[ifanw-nfac0] + opt.swap_thres > qnrm0){
          skipswap = true;
          break;
        }
        #endif

        if(spnorm == 0){
          qnrm1 = MAX(qunw[ifanw-nfac0], qnrm1);
        }else{
          qnrm1 += pow(qunw[ifanw-nfac0], spnorm);
        }
      }

      #ifdef TESTQUALITYALGO
      if(msh.get_tdim() == 2){
        acceptedOldNumerator = qsum0;
        acceptedNewNumerator = qsum1;
        acceptedOldCount = 2;
        acceptedNewCount = 2;
        if constexpr(iquaf == QuaFun::StepDistance){
          if(msh.param->step_distance_cavity_target_average){
            const double weightOld0 =
                step_distance_element_target_weight<MFT,gdim,tdim>(
                    msh,asdmet,iface);
            const double weightOld1 =
                step_distance_element_target_weight<MFT,gdim,tdim>(
                    msh,asdmet,ifac2);
            const double weightNew0 =
                step_distance_element_target_weight<MFT,gdim,tdim>(
                    msh,asdmet,nfac0);
            const double weightNew1 =
                step_distance_element_target_weight<MFT,gdim,tdim>(
                    msh,asdmet,nfac0+1);
            acceptedOldNumerator = weightOld0*quae1 + weightOld1*quae2;
            acceptedNewNumerator = weightNew0*qunw[0] + weightNew1*qunw[1];
            acceptedOldTargetWeight = weightOld0 + weightOld1;
            acceptedNewTargetWeight = weightNew0 + weightNew1;
          }
        }
        const double objectiveOld = step_distance_region_objective(
            acceptedOldNumerator,acceptedOldTargetWeight,
            iquaf == QuaFun::StepDistance
                && msh.param->step_distance_cavity_target_average);
        const double objectiveNew = step_distance_region_objective(
            acceptedNewNumerator,acceptedNewTargetWeight,
            iquaf == QuaFun::StepDistance
                && msh.param->step_distance_cavity_target_average);
        if(globalObjective == nullptr
           && !objective_strictly_improves(objectiveNew,objectiveOld)){
          skipswap = true;
        }
        if(!skipswap && globalObjective != nullptr
           && !globalObjective->accepts_replacement(
                acceptedOldNumerator,acceptedOldCount,
                acceptedOldTargetWeight,
                acceptedNewNumerator,acceptedNewCount,
                acceptedNewTargetWeight)){
          skipswap = true;
        }
      }
      #endif

      if(skipswap) goto cleanup;

      if(spnorm > 0){
        qnrm1 = pow(qnrm1, 1.0/spnorm);
        if constexpr(iquaf == QuaFun::StepDistance){
          if(msh.param->step_distance_cavity_target_average){
            const double weightOld0 =
                step_distance_element_target_weight<MFT,gdim,tdim>(
                    msh,asdmet,iface);
            const double weightOld1 =
                step_distance_element_target_weight<MFT,gdim,tdim>(
                    msh,asdmet,ifac2);
            const double weightNew0 =
                step_distance_element_target_weight<MFT,gdim,tdim>(
                    msh,asdmet,nfac0);
            const double weightNew1 =
                step_distance_element_target_weight<MFT,gdim,tdim>(
                    msh,asdmet,nfac0+1);
            qnrm0 = (weightOld0*quae1 + weightOld1*quae2)
                    /(weightOld0 + weightOld1);
            qnrm1 = (weightNew0*qunw[0] + weightNew1*qunw[1])
                    /(weightNew0 + weightNew1);
          }
        }
      }

    }else{

      edg2pol[0] = msh.fac2poi(iface,ied);
      edg2pol[1] = msh.fac2poi(ifac2,ie2);
      double sz[2], len;
      len = getlenedg_geosz<MFT,gdim,ideg>(msh, edg2pol, sz);
      double qulen = len < 1 ? 1.0 - len
                             : 1.0 - 1.0 / len;
      CPRINTF1(" - new config ied {} len = {} quality = {} \n",ied,len,qulen);

      if(msh.idim == 2){
        qnrm1 = qulen;
      }else{
        // Compute nordev related quality
        double norfa2[3];
        double qudev = -1;
        for(int ifanw = nfac0 + 0; ifanw <= nfac0 + 1; ifanw++){
          getnorfacP1(msh.fac2poi[ifanw], msh.coord, norfa2);
          if(normalize_vec<3>(norfa2)){
            CPRINTF1(" # face {} normal {} {} {} vanishes\n",ifanw,norfa2[0],norfa2[1],norfa2[2]);
            skipswap = true;
            break;
          }
          double dtprd = getprdl2<3>(norCAD[ied], norfa2);
          CPRINTF1(" - new face {} nordev quality = {} \n",ifanw-nfac0,1-dtprd);
          qudev = MAX(qudev,1-dtprd);
        }
        if(skipswap) goto cleanup;
        qnrm1 = qulen * msh.param->qua_surf_wt_quality
              + qudev * msh.param->qua_surf_wt_normal;
      }
    }

    cleanup:
    // Clean up before cavity operator
    msh.set_nface(nfac0);
    if constexpr (gdim >= 3){
      msh.rembpotag(ip1, ithread);
      msh.rembpotag(ip2, ithread);
      msh.rembpotag(msh.fac2poi(iface,ied), ithread);
      msh.rembpotag(msh.fac2poi(ifac2,ie2), ithread);
    }

    if(skipswap) continue;

    #ifndef TESTQUALITYALGO
    if(qnrm1 + opt.swap_thres > qnrm0) continue;
    #endif

    cav.lcfac[1] = ifac2;
    cav.ipins = msh.fac2poi(iface,ied);

    cav.lctet.set_n(0);

    // #ifdef TESTQUALITYALGO
    // std::cout << "skipswap before cheking tets qual = " << skipswap << std::endl;
    // #endif
    if(msh.get_tdim() == 3){
      // Fill tet cavity with edge shell
      intAr1 dum;
      int iopen;
      shell(msh, ip1, ip2, 3, iele0, dum, dum, cav.lctet, &iopen);

      #ifdef TESTQUALITYALGO

      // here it is important to compare the quality of the tets before and after reconnection
      // so far we have that swapping the faces improves the quality of those faces, but if the tets quality go down I think we should reject the op

      msh.tag[ithread]++;
      double qtetsum0 = 0;
      double qtetweight0 = 0;
      double qtetmax0 = -1;
      int ntetqual0 = 0;
      // compute quality of original cavity (also tag tets in there for later)
      for (int itet : cav.lctet){
        msh.tet2tag(ithread,itet) = msh.tag[ithread];
        double qua = metqua<MFT,3,3,iquaf>(msh,AsDeg::P1,asdmet,itet,1.0);
        if constexpr(iquaf == QuaFun::StepDistance){
          if(msh.param->step_distance_cavity_target_average){
            const double weight =
                step_distance_element_target_weight<MFT,3,3>(
                    msh,asdmet,itet);
            qtetsum0 += step_distance_region_contribution(qua,weight,true);
            qtetweight0 += weight;
          }else{
            qtetsum0 += qua;
          }
        }else{
          qtetsum0 += qua;
        }
        ntetqual0++;
        if (qua > qtetmax0) qtetmax0 = qua;
      }

      // now put together reconnected config and probe quality

      // add one temporary tet in mesh
      const int ntet0 = msh.nentt(3);
      msh.set_nentt(3,ntet0+1);
      const int tmpEntt = ntet0;
      const int ipins = cav.ipins;

      double qtetsum1 = 0;
      double qtetweight1 = 0;
      double qtetmax1 = -1;
      int ninvalid = 0;
      int ntetqual1 = 0;
      for (int itet : cav.lctet){

        // loop over tet faces and fetch neighbors
        for (int jj = 0; jj < 4; jj++){

          const int itetnei = msh.tet2tet(itet,jj);
          if (itetnei >= 0 && msh.tet2tag(ithread,itetnei) == msh.tag[ithread]) continue; // neighbor tet in cavity

          // here we have either a neighbor not in cavity or no neighbor, so face jj is boundary face of cavity

          // put together new tet
          const int ip1 = msh.tet2poi(itet,lnofa3[jj][0]);
          const int ip2 = msh.tet2poi(itet,lnofa3[jj][1]);
          const int ip3 = msh.tet2poi(itet,lnofa3[jj][2]);

          if (ip1 == ipins || ip2 == ipins || ip3 == ipins) continue;

          msh.tet2poi(tmpEntt,jj) = ipins;
          msh.tet2poi(tmpEntt,lnofa3[jj][0]) = ip1;
          msh.tet2poi(tmpEntt,lnofa3[jj][1]) = ip2;
          msh.tet2poi(tmpEntt,lnofa3[jj][2]) = ip3;

          msh.tet2ref[tmpEntt] = msh.tet2ref[itet];

          double meas;
          if(!isvalideltP1<3,3>(msh, tmpEntt, nullptr, &meas)){
            ninvalid++;
            continue;
          }
          double qua = metqua<MFT,3,3,iquaf>(msh, AsDeg::P1, AsDeg::P1, tmpEntt, 1.0);
          if constexpr(iquaf == QuaFun::StepDistance){
            if(msh.param->step_distance_cavity_target_average){
              const double weight =
                  step_distance_element_target_weight<MFT,3,3>(
                      msh,AsDeg::P1,tmpEntt);
              qtetsum1 += step_distance_region_contribution(qua,weight,true);
              qtetweight1 += weight;
            }else{
              qtetsum1 += qua;
            }
          }else{
            qtetsum1 += qua;
          }
          ntetqual1++;
          if (qua > qtetmax1) qtetmax1 = qua;
        }
      }
      // revert mesh back to original state
      msh.set_nentt(3,ntet0);

      if(ntetqual1 == 0) skipswap = true;
      acceptedOldNumerator = qtetsum0;
      acceptedNewNumerator = qtetsum1;
      acceptedOldTargetWeight = qtetweight0;
      acceptedNewTargetWeight = qtetweight1;
      acceptedOldCount = ntetqual0;
      acceptedNewCount = ntetqual1;
      if constexpr(iquaf == QuaFun::StepDistance){
        if(ntetqual1 > 0){
          qtetsum0 = step_distance_region_objective(
              qtetsum0,qtetweight0,
              msh.param->step_distance_cavity_target_average);
          qtetsum1 = step_distance_region_objective(
              qtetsum1,qtetweight1,
              msh.param->step_distance_cavity_target_average);
        }
      }
      if(globalObjective == nullptr
         && !objective_strictly_improves(qtetsum1,qtetsum0)) skipswap = true;
      if(!skipswap && globalObjective != nullptr
         && !globalObjective->accepts_replacement(
              acceptedOldNumerator,acceptedOldCount,
              acceptedOldTargetWeight,
              acceptedNewNumerator,acceptedNewCount,
              acceptedNewTargetWeight)){
        skipswap = true;
      }

      // std::cout << "tet related qual info: " << std::endl;
      // std::cout << "qtetmax1 = " << qtetmax1 << ", qtetmax0 = " << qtetmax0 << ", qtetsum1 = " << qtetsum1 << ", qtetsum0 = " << qtetsum0 << std::endl;

      if (ninvalid > 1) skipswap = true;
      if (qtetmax1 < 0) skipswap = true; // just to be safe...
      #endif
    }

    #ifdef TESTQUALITYALGO
    if (skipswap) continue;
    #endif

    if(spnorm >= 0){
      CPRINTF1(" - enact swap ||({},{})|| = {} -> ||({},{})|| = {} \n ",
                                             quae1,quae2,qnrm0,qunw[0],qunw[1],qnrm1);
    }else{
      CPRINTF1(" - enact swap {} -> {} improvement {}\n ",qnrm0,qnrm1,qnrm1 - qnrm0);
    }
    if(tdim == 3) CPRINTF1(" - cavity ntetr {} \n",cav.lctet.get_n());

    // #ifdef TESTQUALITYALGO
    // std::cout << "call cavity_operator " << std::endl;
    // #endif

    int ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithread);

    // #ifdef TESTQUALITYALGO
    // std::cout << "ierro = " << ierro << std::endl;
    // #endif


    #if 0
    if(istop){
      istop = false;
      msh.param->iverb = iverb0;
      msh.param->ivdepth = ivdepth0;
      if(info.done && ierro == 0){
        writeMesh("debug_swapface1", msh);
        printf("## DEBUG STOP HERE \n");
        wait();
      }
    }
    #endif
    if(info.done && ierro == 0){
      if(globalObjective != nullptr){
        globalObjective->replace(
            acceptedOldNumerator,acceptedOldCount,
            acceptedOldTargetWeight,
            acceptedNewNumerator,acceptedNewCount,
            acceptedNewTargetWeight);
      }
      CPRINTF1("-- END swapface did {} - {} -> {} - {} \n",iface,
                                                 ifac2,msh.nface-2,msh.nface-1);
      return -1; // Return did op
    }
  }

  #if 0
  if(istop){
    msh.param->iverb = iverb0;
    msh.param->ivdepth = ivdepth0;
    //writeMesh("debug_swapface1", msh);
    //printf("## DEBUG STOP HERE \n");
    //wait();
  }
  #endif
  return 0;
}


#define BOOST_PP_LOCAL_MACRO(n)\
template int swapface<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical>& msh, \
                                    int iface, swapOptions opt, \
                                    MshCavity &cav, CavWrkArrs &work, \
                                    double *qumx0, double *qnrm1, \
                                    StepDistanceObjectiveState*, int ithread);\
template int swapface<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        >& msh, \
                                    int iface, swapOptions opt, \
                                    MshCavity &cav, CavWrkArrs &work, \
                                    double *qumx0, double *qnrm1, \
                                    StepDistanceObjectiveState*, int ithread);\
template int swapface<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical>& msh, \
                                    int iface, swapOptions opt, \
                                    MshCavity &cav, CavWrkArrs &work, \
                                    double *qumx0, double *qnrm1,\
                                    StepDistanceObjectiveState*, int ithread);\
template int swapface<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        >& msh, \
                                    int iface, swapOptions opt, \
                                    MshCavity &cav, CavWrkArrs &work, \
                                    double *qumx0, double *qnrm1, \
                                    StepDistanceObjectiveState*, int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

} // end namespace
