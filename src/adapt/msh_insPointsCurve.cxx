//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_lineadapt.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../utils/aux_misc.hxx"
#include "../low_lenedg.hxx"
#include "../low_topo.hxx"
#include "../io_libmeshb.hxx"
#include "../adapt/low_increasecav.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../Localization/msh_localization.hxx"
#include "../linalg/det.hxx"
#include "../low_normal.hxx"
#include "../utils/mprintf.hxx"
#include "../msh_checktopo.hxx"


namespace Metris{


// crvlen is the length of the CAD curve computed in the metric field. 
// icor0 is seed corner on edge 
// lnewt: t coord of new points
// iseed0: an edge on this curve close to first t 
template<class MFT>
void insPointsCurve(Mesh<MFT>& msh, int iref, const double* range, const int* lcorn,
                    const dblAr1 &lnewt, const intAr1 &ledge, 
                    int ithrd1, int ithrd2, int ithrd3){
  GETVDEPTH(msh.param);

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
  int mctet = 100, mcfac = 100, mcedg = 10;
  MshCavity cav(mctet,mcfac,mcedg);



  int ninsp = lnewt.get_n();
  intAr1 t2mrk(ninsp);  // Mark t coords inserted
  intAr2 t2lnk(ninsp,2);// Keep track of cavity left and right t bounds
  intAr1 t2sed(ninsp+2);// Cavity seed. Offset by two for range[0] and range[1]
  intAr1 t2poi(ninsp);  // As vertices are created, store them to help seeding
  t2mrk.fill(0);
  int tmark = 1;
  t2sed.fill(-1);
  t2poi.fill(-1);

  intAr1 lshell(10);

  // Assumptions used throughout: ninsp is ordered in increasing order, as is range
  if(DOPRINTS1()){
    CPRINTF1("-- START insPointsCurve iref %d\n",iref);
    //dblAr1(MIN(100,lnewt.get_n()), &lnewt[0]).print();
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
      INCVDEPTH(msh.param);


      double tcur = lnewt[inewt];

      CPRINTF1(" - insert newt %d / %d at t = %f \n", inewt, ninsp, tcur);

      if(msh.param->dbgfull) check_topo(msh,ithrd2);

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

      if(DOPRINTS2()) writeMeshCavity("linecav0",msh,cav);


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

      if(DOPRINTS2()) writeMeshCavity("linecav1",msh,cav);

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

      if(DOPRINTS2()) writeMeshCavity("linecav2",msh,cav);


      // Add tetrahedra. We need the edge shells but also the face supported tets
      // These can be distinct as we can have more faces than edge supported due 
      // to cav increase
      if(msh.nelem > 0){
        for(int iedge : cav.lcedg){
          int iopen;
          int ipoi1 = msh.edg2poi(iedge,0);
          int ipoi2 = msh.edg2poi(iedge,1);
          int iface = msh.edg2fac[iedge];
          int iele0 = msh.fac2tet(iface,0);
          if(iele0 < 0) iele0 = msh.fac2tet(iface,1);
          METRIS_ASSERT(iele0 >= 0);

          intAr1 dum;
          shell3(msh, ipoi1, ipoi2, iele0, lshell, dum, &iopen);
          METRIS_ASSERT(iopen != 0); // should all be open shells.

          for(int ielem : lshell){
            if(msh.tet2tag(ithrd1,ielem) >= msh.tag[ithrd1]) continue;
            msh.tet2tag(ithrd1,ielem) = msh.tag[ithrd1];
            cav.lctet.stack(ielem);
          }
        }// for iedge
        for(int iface : cav.lcfac){
          for(int ii = 0; ii < 2; ii++){
            int ielem = msh.fac2tet(iface,ii);
            if(ielem < 0) continue;
            if(msh.tet2tag(ithrd1, ielem) >= msh.tag[ithrd1]) continue;
            msh.tet2tag(ithrd1, ielem) = msh.tag[ithrd1];
            cav.lctet.stack(ielem);
          }
        }
      }//if msh.nelem > 0

      if(DOPRINTS2()) writeMeshCavity("linecav3",msh,cav);

      double nrmal[3];
      if(msh.idim == 3){
        // Compute normal at ipins. 
        // Use the cavity because we're lazy. It's also not too bad yet as 
        // it shouldn't be skewed, just a ribbon of triangles along a few edges.
        getnorballref<1>(msh,cav.lcfac,-1,nrmal);
      }

      // 
      ierro = increase_cavity(msh, cav, true, ithrd2, ithrd3);
      if(DOPRINTS2()) writeMeshCavity("linecav4",msh,cav);
      if(ierro != 0){
        CPRINTF1(" - failed increase_cavity ierro = %d \n",ierro);
        if(DOPRINTS2()) writeMeshCavity("linecav4",msh,cav);
        if(msh.param->interactive && DOPRINTS1()) wait();
        nerro++;
        goto cleanup;
      }

      /*
      ierro = increase_cavity_Delaunay(msh, cav, ithrd2, nrmal);
      if(DOPRINTS2()) writeMeshCavity("linecav4",msh,cav);
      if(ierro != 0){
        CPRINTF1(" - failed increase_cavity_Delaunay ierro = %d \n",ierro);
        if(DOPRINTS2()) writeMeshCavity("linecav4",msh,cav);
        if(msh.param->interactive && DOPRINTS1()) wait();
        nerro++;
        goto cleanup;
      }

      ierro = increase_cavity_validity(msh,cav,ithrd2);
      if(ierro != 0){
        CPRINTF1(" - failed increase_cavity_validity ierro = %d \n",ierro);
        if(DOPRINTS2()) writeMeshCavity("linecav3",msh,cav);
        if(msh.param->interactive && DOPRINTS1()) wait();
        nerro++;
        goto cleanup;
      }
      if(DOPRINTS2()) writeMeshCavity("linecav5",msh,cav);
      */



      nedg0 = msh.nedge;
      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
        ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd2);
      }}CT_FOR1(ideg);

      //if(ierro != 0){
      //  printf("## WAIT: in insPointsCurve cavity ierro %d \n",ierro);
      //  wait();
      //}

      if(ierro != 0){
        CPRINTF1(" - failed cavity_operator ierro = %d \n",ierro);
        nerro++;
        if(msh.param->interactive && DOPRINTS2()) wait();
        goto cleanup;
      }

      // After ierro check to avoid dangling ipins
      if(msh.param->dbgfull) check_topo(msh,ithrd2);

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

    CPRINTF1(" - insPoint ref %d iter %d inserted %d nerro %d\n",iref,niter,nsucc,nerro);
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
   const intAr1 &ledge, int ithrd1, int ithrd2, int ithrd3);
template void insPointsCurve<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, 
   int iref, const double* range, const int*lcorn, const dblAr1 &lnewt, 
   const intAr1 &ledge, int ithrd1, int ithrd2, int ithrd3);


} //namespace Metris
