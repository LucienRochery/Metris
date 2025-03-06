//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_lineadapt.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../utils/aux_misc.hxx"
#include "../low_lenedg.hxx"
#include "../io_libmeshb.hxx"
#include "../adapt/low_increasecav.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../Localization/msh_localization.hxx"
#include "../linalg/det.hxx"
#include "../low_geo.hxx"
#include "../utils/mprintf.hxx"
#include "../msh_checktopo.hxx"


namespace Metris{


// crvlen is the length of the CAD curve computed in the metric field. 
// icor0 is seed corner on edge 
// lnewt: t coord of new points
// iseed0: an edge on this curve close to first t 
template<class MFT>
void insPointsCurve(Mesh<MFT>& msh, int iref, const double* range, const int* lcorn,
                    const dblAr1 &lnewt, const intAr1 &ledge, int ithrd1, int ithrd2){
  GETVDEPTH(msh);

  const int miter = 10; 
  const int mdseed = 2; // How many t's appart do we create the cavities from

  CavOprOpt opts;
  CavOprInfo info;
  CavWrkArrs work;
  opts.allow_topological_correction = true;
  opts.skip_topo_checks = true;
  opts.allow_remove_points = true;
  opts.dryrun = false;
  opts.geodev1 = 1.0; // lax
  int mcfac = 100, mcedg = 10;
  MshCavity cav(0,mcfac,mcedg);



  int ninsp = lnewt.get_n();
  intAr1 t2mrk(ninsp);  // Mark t coords inserted
  intAr2 t2lnk(ninsp,2);// Keep track of cavity left and right t bounds
  intAr1 t2sed(ninsp+2);// Cavity seed. Offset by two for range[0] and range[1]
  intAr1 t2poi(ninsp);  // As vertices are created, store them to help seeding
  t2mrk.fill(0);
  int tmark = 1;
  t2sed.fill(-1);
  t2poi.fill(-1);

  // Assumptions used throughout: ninsp is ordered in increasing order, as is range
  if(DOPRINTS1()){
    CPRINTF1("-- START insPointsCurve iref %d print first 100 t:",iref);
    dblAr1(MIN(100,lnewt.get_n()), &lnewt[0]).print();
  }
  if(DOPRINTS2()){
    double result[18];
    ego obj = msh.CAD.cad2edg[iref]; 
    int npoi0 = msh.npoin;
    for(int inewt = 0; inewt < ninsp; inewt++){
      double tcur = lnewt[inewt];
      METRIS_ENFORCE(EG_evaluate(obj, &tcur, result) == EGADS_SUCCESS);
      int ipnew = msh.newpoitopo(0, -1);
      msh.newbpotopo(ipnew,0,ipnew);
      for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipnew,ii) = result[ii];
    }
    writeMesh("genPoints_ref"+std::to_string(iref),msh);
    for(int ipoin = npoi0; ipoin < msh.npoin;ipoin++){
      msh.killpoint(ipoin);
    }
  }
  METRIS_ENFORCE(lnewt[0] > range[0]);
  METRIS_ENFORCE(ninsp == 1 || lnewt[1] > lnewt[0]);


  // Fill edges seeds for extremities
  double dt = abs(range[1] - range[0]);
  for(int icor = 0; icor <= 1; icor++){
    // Do it manually because we may have two inciding edges (loop) and we 
    // need to go by t similarity
    bool ifnd = false;
    for(int ib = msh.poi2bpo[lcorn[icor]]; ib >= 0; ib = msh.bpo2ibi(ib,3)){
      if(msh.bpo2ibi(ib,1) != 1) continue;
      int iedge = msh.bpo2ibi(ib,2);
      if(msh.edg2ref[iedge] != iref) continue;

      double t = msh.bpo2rbi(ib,0);
      // This ib entry corresponds to range[0] 
      if(abs(t - range[0]) < 1.0e-6 * dt){
        t2sed[0] = iedge;
        ifnd = true;
        break;
      }else if(abs(t - range[1]) < 1.0e-6 * dt){
        t2sed[ninsp+1] = iedge;
        ifnd = true;
        break;
      }
    }
    METRIS_ASSERT_MSG(ifnd,"Missing ibpoi entries for corners in line");
  }
  for(int ii = 1; ii <= ninsp; ii++){
    t2sed[ii] = ledge[ii-1];
  }

  #ifndef NDEBUG
    try{
      for(int inewt = 0; inewt < ninsp+2; inewt++){
        int iseed = t2sed[inewt];
        METRIS_ASSERT_MSG(iseed >= 0 && iseed < msh.nedge && !isdeadent(iseed,msh.edg2poi),
          "Invalid edge seed "<<iseed);
      }
      for(int inewt = 0; inewt < ninsp; inewt++){
        int iseed = t2sed[inewt+1];
        double tedg[2];
        for(int ii = 0; ii < 2; ii++){
          int ipoin = msh.edg2poi(iseed,ii);
          int ibpoi = msh.poi2ebp(ipoin,1,iseed,iref);
          METRIS_ASSERT(ibpoi >= 0);
          tedg[ii] = msh.bpo2rbi(ibpoi,0);
        }
        double tcur = lnewt[inewt];
        METRIS_ASSERT_MSG(tcur >= tedg[0] && tcur <= tedg[1]
                       || tcur <= tedg[0] && tcur >= tedg[1],
                       "Current t = "<<tcur<<" not in edge bounds = "
                       <<tedg[0]<<"  "<<tedg[1]);
      }
    }catch(const MetrisExcept &e){
      printf("## WRONG SEEDS, DUMP ALL\n");
      for(int inewt = 0; inewt < ninsp+2; inewt++){
        int iseed = t2sed[inewt];
        MPRINTF("iseed %d ref %d vertices ",iseed,msh.edg2ref[iseed]);
        intAr1(getnnod1(msh.curdeg),msh.edg2poi[iseed]).print();
      }
      int nobd = 0;
      for(int inewt = 0; inewt < ninsp; inewt++){
        int iseed = t2sed[inewt+1];
        double tedg[2];
        for(int ii = 0; ii < 2; ii++){
          int ipoin = msh.edg2poi(iseed,ii);
          int ibpoi = msh.poi2ebp(ipoin,1,iseed,iref);
          METRIS_ASSERT(ibpoi >= 0);
          tedg[ii] = msh.bpo2rbi(ibpoi,0);
        }
        double tcur = lnewt[inewt];
        bool inbd = tcur >= tedg[0] && tcur <= tedg[1]
                 || tcur <= tedg[0] && tcur >= tedg[1];
        nobd += !inbd;
        MPRINTF("inewt %d/%d seed %d tcur %f tedg %f %f inbd %d \n",inewt,ninsp+2,iseed,tcur,
                tedg[0], tedg[1], inbd);
      }
      MPRINTF(" out of bounds %d / %d \n",nobd,ninsp);
      throw(e);
    }
  #endif


  // Each t will create a cavity from the left to the right t coordinate
  // As points are inserted, these links are updated to avoid collapses.
  for(int ii = 0; ii < ninsp; ii++){
    t2lnk(ii,0) = MAX(ii-mdseed,-1);    // -1    means range[0]
    t2lnk(ii,1) = MIN(ii+mdseed,ninsp); // ninsp means range[1]
  }


  double result[18];
  ego obj = msh.CAD.cad2edg[iref]; 
  int nedg0;

  // Attempts to correct errors
  for(int niter = 0; niter < miter; niter++){
    int nerro = 0;
    int nsucc = 0;
    for(int inewt = 0; inewt < ninsp; inewt++){
      if(t2mrk[inewt] >= tmark) continue;
      INCVDEPTH(msh);

      CPRINTF1(" - insert newt %d / %d \n", inewt, ninsp);

      double tcur = lnewt[inewt];

      int ipnew = msh.newpoitopo(1, -1);
      int ibnew = msh.newbpotopo(ipnew,1,t2sed[inewt+1]);
      cav.ipins = ipnew;
      int ierro = EG_evaluate(obj, &tcur, result);
      METRIS_ASSERT(ierro == 0);
      for(int ii = 0; ii < msh.idim; ii++) msh.coord(cav.ipins,ii) = result[ii];
      msh.bpo2rbi(ibnew,0) = tcur;

      // Build cavity by going from prv to nxt
      int itprv = t2lnk(inewt,0);
      METRIS_ASSERT(itprv >= -1);
      int ieprv = t2sed[itprv+1];
      METRIS_ASSERT(ieprv >= 0 && !isdeadent(ieprv,msh.edg2poi));

      int itnxt = t2lnk(inewt,1);
      METRIS_ASSERT(itnxt >= -1);
      int ienxt = t2sed[itnxt+1];
      METRIS_ASSERT_MSG(ienxt >= 0 && !isdeadent(ienxt,msh.edg2poi),
        "Invalid ienxt = "<<ienxt<<" w/ itnxt = "<<itnxt<<" ninsp = "<<ninsp);

      #ifndef NDEBUG
      // Check the seed contains the t value 
      {
        int iedge = t2sed[inewt+1];
        double tedg[2];
        for(int ii = 0; ii < 2; ii++){
          int ipoin = msh.edg2poi(iedge,ii);
          int ibpoi = msh.poi2ebp(ipoin,1,iedge,iref);
          METRIS_ASSERT(ibpoi >= 0);
          tedg[ii] = msh.bpo2rbi(ibpoi,0);
        }
        METRIS_ASSERT_MSG(tcur >= tedg[0] && tcur <= tedg[1]
                       || tcur <= tedg[0] && tcur >= tedg[1],
                       "Current t = "<<tcur<<" not in edge bounds = "
                       <<tedg[0]<<"  "<<tedg[1]);
      }
      #endif

      ierro = msh.interpMetBack(cav.ipins,1,t2sed[inewt+1],iref,&result[3]);
      if(ierro != 0){
        CPRINTF1(" - failed interpMetBack ierro = %d \n",ierro);
        nerro++;
        goto cleanup;
      }

      msh.tag[ithrd1]++;
      cav.lcedg.set_n(0);
      cav.lcfac.set_n(0);
      cav.lctet.set_n(0);
      cav.lcedg.stack(t2sed[inewt+1]);
      msh.edg2tag(ithrd1,t2sed[inewt+1]) = msh.tag[ithrd1];

      // Seed cavity from both sides
      for(int it = 0; it <= 1; it++){
        int itext = it == 0 ? itprv : itnxt;
        int ieext = it == 0 ? ieprv : ienxt;

        if(itext == -1 || itext == ninsp || t2mrk[itext] < tmark){
          if(msh.edg2tag(ithrd1,ieext) < msh.tag[ithrd1]){
            cav.lcedg.stack(ieext);
            msh.edg2tag(ithrd1,ieext) = msh.tag[ithrd1];
          }
        }else if(t2mrk[itext] >= tmark){
          // Point has been inserted. Then the seed must contain the mesh vertex
          // with that t. We need to know if that's to the "left" or "right" 
          // (t-wise) of it, in order not to collapse the vertex. Hence
          // do not add this edge if both t coordinates are lesser than tprv
          // or higher than tnxt.  
          double text = lnewt[inewt]; 
          bool iadded = false;
          for(int ii = 0; ii < 2; ii++){
            int ipoin = msh.edg2poi(ieext,ii);
            int ibpoi = msh.poi2ebp(ipoin, 1, ieext, iref);
            METRIS_ASSERT(ibpoi >= 0);
            CPRINTF1(" - iext %d ipoin %d t = %f against text %f \n",
                     it, ipoin, msh.bpo2rbi(ibpoi,0),text);
            if(it == 0 && msh.bpo2rbi(ibpoi,0) >= text
            || it == 1 && msh.bpo2rbi(ibpoi,0) <= text){
              if(msh.edg2tag(ithrd1,ieext) < msh.tag[ithrd1]){
                cav.lcedg.stack(ieext);
                msh.edg2tag(ithrd1,ieext) = msh.tag[ithrd1];
              }
              iadded = true;
              break;
            }
          }

          if(!iadded){
            // If we didn't add that edge, we should add the neighbour
            // We have t2poi to help us. 
            int ipext = t2poi[itext];
            int iver  = msh.template getveredg<1>(ieext, ipext);
            METRIS_ASSERT(iver >= 0);
            int iedge = msh.edg2edg(ieext,1-iver);
            METRIS_ASSERT_MSG(iedge >= 0, 
              "Left/right link neighbour has no right/left neighbour");
            if(msh.edg2tag(ithrd1,iedge) < msh.tag[ithrd1]){
              cav.lcedg.stack(iedge);
              msh.edg2tag(ithrd1,iedge) = msh.tag[ithrd1];
            }
          }
        }//if itext
      }// for it 

      if(DOPRINTS2()) writeMeshCavity("linecav0",msh,cav,ithrd2);


      // fill the cavity adding edges against the t coordinate going
      // towards this t. The edge that contains t is already in the cavity
      // from t2sed. Skip this one. 
      for(int icav = 0; icav < cav.lcedg.get_n(); icav++){
        int iedge = cav.lcedg[icav];
        CPRINTF1(" - debug iedge %d = %d %d , seed skip %d \n",
          iedge,msh.edg2poi(iedge,0),msh.edg2poi(iedge,1),t2sed[inewt+1]);
        if(iedge == t2sed[inewt+1]) continue;
        int imin = -1;
        double dmin = 1.0e30;
        for(int ii = 0; ii < 2; ii++){
          int ipoin = msh.edg2poi(iedge,ii);
          int ibpoi = msh.poi2ebp(ipoin,1,iedge,iref);
          METRIS_ASSERT(ibpoi >= 0);
          double t = msh.bpo2rbi(ibpoi,0);
          CPRINTF1(" - test edge %d ver %d = %d t = %f cur = %f\n",
                   iedge,ii,ipoin,t,tcur);
          double dist = abs(t-tcur);
          if(dist < dmin){
            imin = ii;
            dmin = dist;
          }
        }// for int ii 
        METRIS_ASSERT(imin >= 0);
        // Neighbour against that t coordinate
        int ienei = msh.edg2edg(iedge,1-imin);
        CPRINTF1(" - add neighb %d against %d = %d\n",ienei,1-imin,
                 msh.edg2poi(iedge,imin));
        if(msh.edg2tag(ithrd1,ienei) < msh.tag[ithrd1]){
          cav.lcedg.stack(ienei);
          msh.edg2tag(ithrd1,ienei) = msh.tag[ithrd1];
        }
      }// for int icav

      if(DOPRINTS2()) writeMeshCavity("linecav1",msh,cav,ithrd2);

      // Add faces to the cavity. 
      cav.lcfac.set_n(0);
      for(int iedge : cav.lcedg){
        int ifac1 = msh.edg2fac[iedge];
        METRIS_ASSERT(ifac1 >= 0);
        if(msh.fac2tag(ithrd1,ifac1) >= msh.tag[ithrd1]) continue;
        msh.fac2tag(ithrd1,ifac1) = msh.tag[ithrd1];
        cav.lcfac.stack(ifac1);
        int iedl = getedgfac(msh,ifac1,msh.edg2poi(iedge,0),msh.edg2poi(iedge,1));
        METRIS_ASSERT(iedl >= 0);
        int ifac2 = msh.fac2fac(ifac1,iedl);
        if(ifac2 == -1){
          continue;
        }else if(ifac2 < -1){ // Non manifold, add all
          int ied;
          ifac2 = - ifac2 - 2;
          // starting with ifac2
          if(msh.fac2tag(ithrd1,ifac2) < msh.tag[ithrd1]){
            msh.fac2tag(ithrd1,ifac2) = msh.tag[ithrd1];
            cav.lcfac.stack(ifac2);
          }
          while(getnextfacnm(msh,ifac1,
                             msh.edg2poi(iedge,0),msh.edg2poi(iedge,1),
                             &ifac2, &ied)){
            if(msh.fac2tag(ithrd1,ifac2) >= msh.tag[ithrd1]) continue;
            msh.fac2tag(ithrd1,ifac2) = msh.tag[ithrd1];
            cav.lcfac.stack(ifac2);
          }
        }else{
          if(msh.fac2tag(ithrd1,ifac2) >= msh.tag[ithrd1]) continue;
          msh.fac2tag(ithrd1,ifac2) = msh.tag[ithrd1];
          cav.lcfac.stack(ifac2);
        }// if ifac2 == -1
      }// for iedge : cav.lcedg

      if(DOPRINTS2()) writeMeshCavity("linecav2",msh,cav,ithrd2);

      double nrmal[3];
      if(msh.idim == 3){
        // Compute normal at ipins. 
        // Use the cavity because we're lazy. It's also not too bad yet as 
        // it shouldn't be skewed, just a ribbon of triangles along a few edges.
        getnorballref<1>(msh,cav.lcfac,-1,nrmal);
        cav.nrmal = nrmal;
      }
      ierro = increase_cavity_Delaunay(msh, cav, ithrd2, nrmal);
      if(DOPRINTS2()) writeMeshCavity("linecav4",msh,cav);
      if(ierro != 0){
        CPRINTF1(" - failed increase_cavity_Delaunay ierro = %d \n",ierro);
        if(DOPRINTS2()) writeMeshCavity("linecav4",msh,cav,ithrd2);
        if(msh.param->interactive && DOPRINTS1()) wait();
        nerro++;
        goto cleanup;
      }

      ierro = increase_cavity2D(msh,cav,ithrd2);
      if(ierro != 0){
        CPRINTF1(" - failed increase_cavity2D ierro = %d \n",ierro);
        if(DOPRINTS2()) writeMeshCavity("linecav3",msh,cav,ithrd2);
        if(msh.param->interactive && DOPRINTS1()) wait();
        nerro++;
        goto cleanup;
      }
      if(DOPRINTS2()) writeMeshCavity("linecav3",msh,cav,ithrd2);



      nedg0 = msh.nedge;
      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
        ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd2);
      }}CT_FOR1(ideg);

      if(ierro != 0){
        CPRINTF1(" - failed cavity_operator ierro = %d \n",ierro);
        nerro++;
        if(msh.param->interactive && DOPRINTS2()) wait();
        goto cleanup;
      }

      nsucc++;
      t2mrk[inewt] = tmark;
      t2poi[inewt] = cav.ipins;
      t2sed[inewt+1] = nedg0;

      if(DOPRINTS2()) writeMesh("linestep",msh);

      CPRINTF1(" - Debug us = %d left link = %d right = %d \n",
               inewt,t2lnk(inewt,0),t2lnk(inewt,1));
      // Rework the links. People who left link to further left than us should
      // link to us instead. Idem right. 
      for(int ii = 1; ii <= mdseed; ii++){
        // The right link of people to our left:
        if(inewt-ii >= 0 && t2lnk(inewt-ii,1) > inewt){
          CPRINTF2(" - update right link of %d: %d -> %d \n",
                   inewt-ii,t2lnk(inewt-ii,1),inewt);
          t2lnk(inewt-ii,1) = inewt;
        }
        // The left link of people to our right:
        if(inewt+ii < ninsp  && t2lnk(inewt+ii,0) < inewt){
          CPRINTF2(" - update left link of %d: %d -> %d \n",
                   inewt+ii,t2lnk(inewt+ii,0),inewt);
          t2lnk(inewt+ii,0) = inewt;
        }
      }// for ii < mdseed

      // Update seeds. Unfortunately, we don't know how many people were 
      // affected, as potentially even the whole line mesh may have been
      // destroyed (not likely, but a possibility)
      // Walk left and right, as long as we encounter dead seeds. 
      for(int ipm = -1; ipm <= 1; ipm += 2){
        for(int kk = 1; inewt + ipm*kk >= -1 && inewt + ipm*kk <= ninsp; kk++){
          int iother = inewt + ipm*kk;
          int iseed  = t2sed[iother + 1];
          METRIS_ASSERT(iseed >= 0 && iseed < msh.nedge);
          if(!isdeadent(iseed,msh.edg2poi)) break;

          double tother = iother >= 0 && iother < ninsp ? lnewt[iother]
                        : iother == -1 ? range[0] : range[1];
          CPRINTF1(" - update seed for %d range = %f %f \n",iother,range[0],range[1]);
          bool ifnd = false;
          // This is in fact a very small loop (most often 2 elts)
          for(int iedge = nedg0; iedge < msh.nedge; iedge++){
            METRIS_ASSERT(!isdeadent(iedge,msh.edg2poi));
            if(msh.edg2ref[iedge] != iref) continue;
            double tedg[2];
            for(int jj = 0; jj < 2; jj++){
              int ipoin = msh.edg2poi(iedge,jj);
              int ibpoi = msh.poi2ebp(ipoin,1,iedge,iref);
              METRIS_ASSERT(ibpoi >= 0);
              tedg[jj] = msh.bpo2rbi(ibpoi,0);
            }// for jj
            CPRINTF1(" - seek t = %f tedg = %f %f \n",tother,tedg[0],tedg[1]);
            if(tother >= tedg[0] && tother <= tedg[1]
            || tother >= tedg[1] && tother <= tedg[0]){
              ifnd = true;
              CPRINTF1(" - updated seed for it = %d: %d -> %d \n",iother,
                       t2sed[iother+1],iedge);
              t2sed[iother+1] = iedge;
            }
          }// for iedge
          METRIS_ASSERT_MSG(ifnd,
                       "Failed to update seed because t not found in new edges")

        }// for int kk
      }// for int ipm

      continue;
      cleanup:
      msh.killpoint(cav.ipins);
    }// for inewt

    CPRINTF1(" - iter %d inserted %d nerro %d\n",niter,nsucc,nerro);
    #ifndef NDEBUG
    if(DOPRINTS2()){
      writeMesh("line"+std::to_string(iref)+"iter"+std::to_string(niter),msh);
      if(msh.param->interactive){
        CPRINTF1("## WAIT HERE 1 \n");
        wait();
      }
    }
    #endif
    if(nerro == 0) break;
  }// for niter

  if(DOPRINTS2()){
    writeMesh("line"+std::to_string(iref),msh);
    if(msh.param->interactive){
      CPRINTF1("## WAIT HERE 2 \n");
      wait();
    }
  }
}


template void insPointsCurve<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, 
   int iref, const double* range, const int*lcorn, const dblAr1 &lnewt, 
   const intAr1 &ledge, int ithrd1, int ithrd2);
template void insPointsCurve<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, 
   int iref, const double* range, const int*lcorn, const dblAr1 &lnewt, 
   const intAr1 &ledge, int ithrd1, int ithrd2);





#if 0


// crvlen is the length of the CAD curve computed in the metric field. 
// icor0 is seed corner on edge 
// lnewt: t coord of new points
// iseed0: an edge on this curve close to first t 
template<class MFT>
void insPointsCurve(Mesh<MFT>& msh, int iref, const double* range, const int* lcorn,
                    const dblAr1 &lnewt, int ithrd1, int ithrd2){
  GETVDEPTH(msh);

  int npoi0 = msh.npoin;

  CavOprOpt opts;
  CavOprInfo info;
  CavWrkArrs work;
  opts.allow_topological_correction = true;
  opts.skip_topo_checks = true;
  opts.allow_remove_points = true;
  opts.dryrun = false;
  opts.geodev1 = 1.0; // lax
  int mcfac = 100, mcedg = 10;
  MshCavity cav(0,mcfac,mcedg);


  int ninsp = lnewt.get_n();
  intAr1 t2mrk(ninsp+1);
  t2mrk.set_n(ninsp+1);
  t2mrk.fill(0);
  int tmark = 1;
  double result[18], bary[2], coopr[3], normal[3];
  ego obj = msh.CAD.cad2edg[iref]; 

  // Iterate in case some errors require restart 
  for(int niter = 0; niter < 10; niter++){

    int nerro = 0;
    int nsucc = 0;
    double tprev = range[0];
    int iseed1 = -1;
    {
      double dt = abs(range[1] - range[0]);
      // Do it manually because we may have two inciding edges (loop) and we 
      // need to go by t similarity
      for(int ib = msh.poi2bpo[lcorn[0]]; ib >= 0; ib = msh.bpo2ibi(ib,3)){
        if(msh.bpo2ibi(ib,1) != 1) continue;
        int itmp = msh.bpo2ibi(ib,2);
        if(msh.edg2ref[itmp] != iref) continue;

        double t = msh.bpo2rbi(ib,0);
        //printf("Test ib %d t %f range %f %f \n",ib,t,range[0],range[1]);
        if(abs(t - range[0]) < 1.0e-6 * dt){
          iseed1 = itmp;
          break;
        }
      }
      METRIS_ASSERT_MSG(iseed1 >= 0, "COULD NOT FIND t CORNER IN BPOIS")
      METRIS_ASSERT(msh.edg2ref[iseed1] == iref);
    }
    // Loop one over, last is range[1]. We still need to keep track of that. 
    for(int inewt = 0; inewt < ninsp + 1; inewt++){
      INCVDEPTH(msh);

      double tcoor;
      int ierro;
      if(inewt < ninsp){
        if(t2mrk[inewt] >= tmark) continue;
        tcoor = lnewt[inewt]; 
        cav.ipins = msh.newpoitopo(1,-1);

        ierro = EG_evaluate(obj, &tcoor, result);
        METRIS_ASSERT(ierro == 0);
        for(int ii = 0; ii < msh.idim; ii++) msh.coord(cav.ipins,ii) = result[ii];
      }else{
        tcoor = range[1];
        cav.ipins = lcorn[1];
      }
      cav.lcedg.set_n(0);
      cav.lcfac.set_n(0);
      cav.nrmal = NULL;

      CPRINTF1(" - inewt %d / %d tcoor = %f \n",inewt,ninsp,tcoor);



      int iseed2 = iseed1;


      { // namespace preservation (goto)
        // get the edge seed by localizing in t space on the front mesh 
        int itag0 = msh.tag[ithrd1];
        if(msh.idim == 2){
          ierro = locMesh<2,1,1>(msh, &iseed2, msh.coord[cav.ipins], 1, &tcoor, iref, &result[3],
                                 coopr, bary, 1.0e-6, ithrd1);
        }else{
          ierro = locMesh<3,1,1>(msh, &iseed2, msh.coord[cav.ipins], 1, &tcoor, iref, &result[3],
                                 coopr, bary, 1.0e-6, ithrd1);
        }
        if(ierro != 0){
          CPRINTF1(" - failed locMesh ierro = %d \n",ierro);
          nerro++;
          goto cleanup;
        }
        int itag1 = msh.tag[ithrd1];
        if(itag1 != itag0+1) METRIS_THROW_MSG(TODOExcept(), 
                                        "Revise getting edges from iedg1 to iedg2")

        // Now we're going to be a little too clever for our own good, probably. 
        // Since we're calling locMesh with tdim = 1, no lower dims will be called
        // Hence the tag is only prev tag + 1
        // And the tagged elements are all those between iseed1 and iseed2.
        // These are thus the guys we want in our cavity. 
        int iecur = iseed1;
        int ieprv = -1;
        while(true){
          cav.lcedg.stack(iecur);
          METRIS_ASSERT(msh.edg2tag(ithrd1,iecur) == msh.tag[ithrd1]);

          bool ifnd = false;
          for(int ii = 0; ii < 2; ii++){
            int ienei = msh.edg2edg(iecur,ii);
            if(ienei < 0) continue;
            if(ienei == ieprv) continue;
            if(msh.edg2tag(ithrd1,ienei) < msh.tag[ithrd1]) continue;
            ifnd = true;
            ieprv = iecur;
            iecur = ienei;
          }
          if(!ifnd) break;
          if(iecur == iseed1){
            METRIS_THROW_MSG(TODOExcept(), "Investigate this looping after loc");
            break; // could be a loop
          }
        }
        CPRINTF1(" - step 1: edge cavity n = %d\n",cav.lcedg.get_n());
        msh.tag[ithrd2]++;
        for(int iedgl = 0; iedgl < cav.lcedg.get_n(); iedgl++){
          INCVDEPTH(msh);
          int iedge = cav.lcedg[iedgl];
          CPRINTF1(" - iedgl %d / %d \n",iedgl,cav.lcedg.get_n());

          // Skip edges that contain only new points. Remove the edge
          if(msh.edg2poi(iedge,0) >= npoi0 
          && msh.edg2poi(iedge,1) >= npoi0 ){
            cav.lcedg[iedgl] = cav.lcedg[cav.lcedg.get_n() - 1];
            cav.lcedg.pop();
            CPRINTF1(" - removed edge i %d = %d \n",iedgl,iedge);
            iedgl--;
            continue;
          }

          int ifac1 = msh.edg2fac[iedge];
          CPRINTF1(" - iedgl %d / %d consider ifac1 = %d tad = %d <? %d \n", 
                   iedgl,cav.lcedg.get_n(),ifac1,
                   msh.fac2tag(ithrd2,ifac1),msh.tag[ithrd2]);
          if(msh.fac2tag(ithrd2,ifac1) >= msh.tag[ithrd2]) continue;
          msh.fac2tag(ithrd2,ifac1) = msh.tag[ithrd2];
          cav.lcfac.stack(ifac1);

          int iedl = getedgfac(msh,ifac1,msh.edg2poi(iedge,0),msh.edg2poi(iedge,1));
          METRIS_ASSERT(iedl >= 0);
          int ifac2 = msh.fac2fac(ifac1,iedl);
          CPRINTF1(" - consider ifac2 = %d tad = %d <? %d \n", ifac2,
                   msh.fac2tag(ithrd2,MAX(ifac2,0)),msh.tag[ithrd2]);
          if(ifac2 >= 0){
            if(msh.fac2tag(ithrd2,ifac2) < msh.tag[ithrd2]){
              msh.fac2tag(ithrd2,ifac2) = msh.tag[ithrd2];
              cav.lcfac.stack(ifac2); 
            }
          }

        }
      } // local namespace
      CPRINTF1(" - step 2: edge cav n = %d face cav n = %d\n",
                            cav.lcedg.get_n(),cav.lcfac.get_n());

      if(DOPRINTS2()){
        writeMeshCavity("inscurvecav0",msh,cav);
        writeMesh("inscurvepreop",msh);
      }

      // Now forget about previous seed.
      iseed1 = iseed2;
      if(inewt < ninsp){
        int ibins = msh.newbpotopo(cav.ipins,1,iseed2);
        msh.bpo2rbi(ibins,0) = tcoor;
      }

      // Compute a normal if needed.  
      if(msh.idim == 3){
        cav.nrmal = normal;
        getnorballref<1>(msh, cav.lcfac, -1, normal);
      }



      // Then interpolate metric at this point 
      if(inewt < ninsp){
        ierro = msh.interpMetBack(cav.ipins,1,iseed1,iref,&result[3]);
        if(ierro != 0){
          CPRINTF1(" - failed interpMetBack ierro = %d \n",ierro);
          nerro++;
          goto cleanup;
        }
      }


      // Now increase cavity
      ierro = increase_cavity2D(msh,cav,ithrd1);
      if(ierro != 0){
        CPRINTF1(" - failed increase_cavity2D ierro = %d \n",ierro);
        if(msh.param->dbgfull && DOPRINTS1()) wait();
        nerro++;
        goto cleanup;
      }

      if(DOPRINTS2()) writeMeshCavity("inscurvecav1",msh,cav);

      if(msh.idim == 2){
        increase_cavity_Delaunay(msh, cav, ithrd1);
        if(DOPRINTS2()) writeMeshCavity("inscurvecav2",msh,cav);
      }

      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
        ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd1);
      }}CT_FOR1(ideg);

      if(ierro != 0){
        CPRINTF1(" - failed cavity_operator ierro = %d \n",ierro);
        nerro++;
        if(msh.param->dbgfull && DOPRINTS1()) wait();
        goto cleanup;
      }

      if(msh.param->dbgfull)  check_topo(msh);

      nsucc++;
      // Restart seed from one of the new edges. The one with t going towards next
      t2mrk[inewt] = tmark;

      if(inewt < ninsp){ // Last will throw, also useless
        bool ifnd = false;
        for(int ii = 1; ii <= 2; ii++){
          int iedge = msh.nedge - ii;
          int iver = msh.template getverent<1>(iedge, 1, cav.ipins);
          METRIS_ASSERT(iver >= 0);
          int ipother = msh.edg2poi(iedge,1-iver);
          METRIS_ASSERT(ipother >= 0);
          int ibother = msh.poi2ebp(ipother,1,iedge,iref);
          METRIS_ASSERT(ibother >= 0);
          double tother = msh.bpo2rbi(ibother,0);
          // If this t coord is closer to previous t, then we're going backwards
          if(abs(tother - tprev) <= abs(tcoor - tprev)) continue;
          ifnd = true;
          iseed1 = iedge;
        }
        METRIS_ASSERT(ifnd);
        tprev = tcoor;
      }

      continue;
      cleanup:
      msh.killpoint(cav.ipins);
    }// for int inewt


    CPRINTF1(" - iter %d inserted %d nerro %d\n",niter,nsucc,nerro);
    //#ifndef NDEBUG
    //if(nsucc == 1 && nerro == 1 && msh.param->dbgfull){
    //  MPRINTF("WAIT AT NSUCC 1 NERRO 1 \n");
    //  wait();
    //}
    //#endif
    if(DOPRINTS2()) writeMesh("insPoints_ref"+std::to_string(iref)+"_iter"+std::to_string(niter),msh);
    if(nerro == 0) break;

  }// for int niter


}


template void insPointsCurve<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, 
   int iref, const double* range, const int*lcorn, const dblAr1 &lnewt, int ithrd1, int ithrd2);
template void insPointsCurve<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, 
   int iref, const double* range, const int*lcorn, const dblAr1 &lnewt, int ithrd1, int ithrd2);

#endif







} //namespace Metris
