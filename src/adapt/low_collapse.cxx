//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../adapt/low_collapse.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../low_topo.hxx"
#include "../low_normal.hxx"
#include "../utils/mprintf.hxx"
#include "../msh_structs.hxx"
#include "../io_libmeshb.hxx"
#include "../msh_checktopo.hxx"
#include "../adapt/low_increasecav.hxx"
#include "../low_lenedg.hxx"
#include "../low_geo.hxx"

#include <unordered_set>

namespace Metris{

/*
Collapse a vertex in short edge 
Cavity is passed in to reuse allocations
*/
template<class MFT>
int colledgsurf(Mesh<MFT>& msh, int tdim, int ientt, int iedl, double qmax_suf, 
                MshCavity &cav, CavWrkArrs &work, 
                intAr1 &lerro, int ithrd1, int ithrd2, int ithrd3){

  GETVDEPTH(msh.param);


  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  //METRIS_ASSERT(ithrd3 >= 0 && ithrd3 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);
  //METRIS_ASSERT(ithrd1 != ithrd3);
  //METRIS_ASSERT(ithrd2 != ithrd3);

  intAr2& ent2poi = msh.ent2poi(tdim);

 
  CavOprOpt  opts;
  CavOprInfo info;
  opts.allow_topological_correction = true; // To fetch missing edges
  opts.skip_topo_checks = false;
  opts.allow_remove_points = true;
  //opts.dryrun   = true;
  opts.dryrun   = false;
  opts.qmax_suf = qmax_suf;
  cav.lcedg.set_n(0);
  cav.lcfac.set_n(0);
  cav.lctet.set_n(0);

  int ierro = 0;

  auto lnoed = tdim == 2 ? lnoed2 : lnoed3;

  int ip1 = ent2poi(ientt,lnoed[iedl][0]);
  int ip2 = ent2poi(ientt,lnoed[iedl][1]);


  CPRINTF1("-- START colledgsurf ientt = %d tdim %d iedl = %d = (%d,%d) \n",
           ientt,tdim,iedl,ip1,ip2);
  if(DOPRINTS2()) writeMesh("debug_collapse0.meshb",msh);

  // Track best insertion dry run to run it at the end
  double qmabest = 1.0e30;
  int ivebest = -1; // p collapse vertex of 
  int ipibest = -1; // p insert 

  int tdimp[2], tdimc;
  for(int ii = 0; ii < 2 ; ii ++){
    int ipoin = ent2poi(ientt,lnoed[iedl][ii]);
    tdimp[ii] = msh.getpoitdim(ipoin);
    if(DOPRINTS2() && tdimp[ii] == 0){
      CPRINTF2("Topo dim 0 point %d \n",ipoin);
      int ib = msh.poi2bpo[ipoin];
      CPRINTF2("poi2bpo = %d bpo2ibi = ",ib);
      intAr1(nibi,msh.bpo2ibi[ib]).print();
    }
  }
  CPRINTF1(" - topo dims %d %d \n",tdimp[0],tdimp[1]);

  tdimc = MAX(tdimp[0], tdimp[1]);

  bool imani;
  int  iopen;
  for(int iver = 0; iver < 2; iver++){
    // Collapse this one 
    int ipcol = ent2poi(ientt,lnoed[iedl][iver]);

    // Collapse only among the highest-dimensional points
    if(tdimp[iver] < tdimc) continue;

    int ibcol = msh.poi2bpo[ipcol];
    if(ibcol >= 0){
      if(msh.bpo2ibi(ibcol,1) == 0) continue; // Skip corners
    }

    ierro = ball(msh, ipcol, cav.lcedg, cav.lcfac, cav.lctet,
                 &iopen, ithrd1);

    //if(tdim == 2){
    //  ierro = ball2(msh,ipcol,ientt,cav.lcfac,cav.lcedg,&iopen,&imani,ithrd1);
    //  if(!imani) METRIS_THROW_MSG(TODOExcept(), "Non manifold not handled");
    //}else{
    //  ierro = ball3(msh,ipcol,ientt,cav.lctet,&iopen,ithrd1);
    //  METRIS_ASSERT(ierro == 0);
    //  if(iopen){
    //    int tdimp = msh.poi2ent(ipcol,1);
    //    METRIS_ASSERT(tdimp == 1 || tdimp == 2);
    //    int iface = -1;
    //    if(tdimp == 1){
    //      int iedge = msh.poi2ent(ipcol,0);
    //      iface = msh.edg2fac[iedge];
    //    }else{
    //      iface = msh.poi2ent(ipcol,0);
    //    }
    //    ierro = ball2(msh,ipcol,iface,cav.lcfac,cav.lcedg,&iopen,&imani,ithrd1);
    //  }
    //}
    CPRINTF1(" - try collapse poi = %d seed ntetr = %d nface = %d nedge = %d \n",
                    ipcol,cav.lctet.get_n(),cav.lcfac.get_n(),cav.lcedg.get_n());
    //METRIS_ASSERT(iopen == 0); // open won't collapse
    METRIS_ASSERT(ierro == 0);


    // Try the cavity call with different ipins in neighbours of ipcol
    msh.tag[ithrd1]++;
    int tag0 = msh.tag[ithrd1]; // We'll reuse the tag for elements in subroutines
    // We'll stack on top of the ball, so we need to be able to prune to nbalf
    // with each attempt, as well as restrict search of ipins to ball 
    int nbalt = cav.lctet.get_n();
    int nbalf = cav.lcfac.get_n();
    int nbale = cav.lcedg.get_n();
    intAr1& lcent = cav.lcent(msh.get_tdim());
    for(int icent : lcent){
      METRIS_ASSERT(!isdeadent(icent,ent2poi));

      // Doesn't change but easy to get it here 
      for(int ive2 = 0; ive2 < msh.get_tdim() + 1; ive2 ++){
        int ipins = ent2poi(icent,ive2);
        if(ipins == ipcol) continue;
        if(msh.poi2tag(ithrd1,ipins) >= tag0) continue;
        msh.poi2tag(ithrd1,ipins) = tag0;

        // Check ipins has same (or lower) topological dimension as ipcol
        // e.g. a triangle point can be collapsed and reconnection done to an edge
        // point, but not to a volume point. 
        if(msh.getpoitdim(ipins) > tdimp[iver]){
          CPRINTF1(" - point %d dim %d > ipins dim %d -> reject reconnection\n",
            ipins,msh.getpoitdim(ipins),tdimp[iver]);
          continue;
        }

        // Check ipins has same ref in topo dim of ipcol.
        // Counter-example: a boundary point, ball's boundary hits other surface refs. 
        // This does not happen in the volume. Indeed, if a point's ball
        // has two domain refs, then the point itself has two domain refs.
        // Hence it was actually a boundary point. 
        if(tdimp[iver] < msh.get_tdim()){
          // As we never collapse a corner, the point's tdim is >= 1.
          // Hence it has a unique ref to the lowest dim entity group. 
          METRIS_ASSERT(cav.lcent(tdimp[iver]).get_n() > 0);
          int iref = msh.ent2ref(tdimp[iver])[cav.lcent(tdimp[iver])[0]];
          #ifndef NDEBUG
          for(int ient1 : cav.lcent(tdimp[iver])){
            METRIS_ASSERT(iref == msh.ent2ref(tdimp[iver])[ient1]);
          }
          #endif

          if(msh.getpoitdim(ipins) == tdimp[iver]){
            // Easy, just check seed reference
            int ient1 = msh.poi2ent(ipins,0);
            METRIS_ASSERT(msh.poi2ent(ipins,1) == tdimp[iver]);
            int iref1 = msh.ent2ref(tdimp[iver])[ient1];
            if(iref1 != iref){
              CPRINTF1(" - point %d dim %d = dim ipins but ref %d != %d\n",
                       ipins,tdimp[iver],iref1,iref);
              continue;
            }
          }else{
            // Use boundary info
            bool ifnd = false;
            for(int ibpoi = msh.poi2bpo[ipins]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
              int tdim1 = msh.bpo2ibi(ibpoi,1);
              if(tdim1 != tdimp[iver]) continue;
              int ient1 = msh.bpo2ibi(ibpoi,2);
              int iref1 = msh.ent2ref(tdimp[iver])[ient1];
              CPRINTF1(" - check entity %d dim %d ref %d ipcol ref = %d\n",
                       ient1,tdim1,iref1,iref);
              if(iref1 != iref) continue;
              CPRINTF1(" -> found ref\n");
              ifnd = true;
              break;
            }
            if(!ifnd){
              CPRINTF1(" - did not find ref %d dim %d of ipcol in ipins %d refs\n",
                       iref,tdimp[iver],ipins);
              continue;
            }
          }
        }// if tdimp[iver] < msh.get_tdim()

        cav.ipins = ipins;
        cav.lctet.set_n(nbalt);
        cav.lcfac.set_n(nbalf); // Revert to simple ball 
        cav.lcedg.set_n(nbale); // Revert 
        CPRINTF1(" - try reinsert point %d tag = %d vs %d \n",
                             ipins,msh.poi2tag(ithrd1,ipins),tag0);


        if(DOPRINTS2()) writeMeshCavity("collapse_cavity0.meshb", msh, cav);

        //static int nwarnprt = 0;
        //if(nwarnprt++ < 10) printf("## PUT BACK DELAUNAY IN LOW COLLAPSE\n");
        ierro = increase_cavity(msh, cav, true, ithrd1, ithrd2);
        if(ierro != 0){
          CPRINTF1("# increase_cavity error %d \n",ierro);
          continue;
        }
        //// Increase cavity with Delaunay criterion
        //ierro = increase_cavity_Delaunay(msh, cav, ithrd2, nrmal);
        //if(ierro != 0) continue;
        ////int nprem = increase_cavity_lenedg(msh,cav,ipins,ithrd2,ithrd3);
        //ierro = increase_cavity_validity(msh,cav,ithrd2);
        //if(ierro > 0) continue;

        if(DOPRINTS2()) writeMeshCavity("collapse_cavity1.meshb", msh, cav);

        //ierro = collrejcav_dens(msh, cav, ithrd2, ithrd3);
        //if(ierro > 0){
        //  CPRINTF1(" # reject cavity\n");
        //  continue;
        //}


        CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
          ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd2);
        }}CT_FOR1(ideg);

        if(ierro > 0){
          lerro[ierro-1]++;
        }
        //printf("Debug wait;\n");
        // If operation was done, out
        if(info.done){
          CPRINTF1("-- END colledgsurf successful using ipcol = %d ipins = %d \n",ipcol,ipins);
          if(DOPRINTS2()) writeMesh("debug_collapse1.meshb",msh);
          msh.poi2ent(ipcol,0) = -1;
          msh.poi2ent(ipcol,1) = -1;
        }
        if(msh.param->dbgfull){
          check_topo(msh, msh.nbpoi, msh.npoin, msh.nedge, msh.nface, msh.nelem, ithrd2); 
        }
        if(info.done) return 0;

        CPRINTF1(" - return qmax = %f \n",info.qmax_end);
        if(info.qmax_end < qmabest && ierro == 0){
          CPRINTF1(" - new best quality !\n");
          qmabest = info.qmax_end;
          ivebest = iver;
          ipibest = ipins;
        }
      }
    }
  }

  if(ivebest == -1){
    // Failed to find a reinsertion candidate. 
    // Let's not rely on exceptions as failure here is not invalid behaviour 
    // and could happen quite frequently 
    return 1;
  }

  cav.ipins = ipibest;
  opts.dryrun = false;

  if(ivebest == 0){ // Otherwise ball hasn't changed !
    int ipcol = ent2poi(ientt,lnoed[iedl][ivebest]);

    ierro = ball(msh, ipcol, cav.lcedg, cav.lcfac, cav.lctet,
                 &iopen, ithrd1);

    ierro = increase_cavity(msh, cav, true, ithrd1, ithrd2);
  }

  if(DOPRINTS2()) writeMeshCavity("collapse_cavity0.meshb", msh, cav);


  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd2);
  }}CT_FOR1(ideg);

  METRIS_ASSERT(ierro == 0);
  METRIS_ASSERT(abs(info.qmax_end - qmabest) < 1.0e-6);
  METRIS_ASSERT(info.done);

  return 0;
}

template int colledgsurf<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
                                          int tdim, int ientt, int iedl, double qmax_suf, 
                                          MshCavity &cav, CavWrkArrs &work, 
                                          intAr1 &lerro, int ithrd1, int ithrd2, int ithrd3);
template int colledgsurf<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
                                          int tdim, int ientt, int iedl, double qmax_suf, 
                                          MshCavity &cav, CavWrkArrs &work, 
                                          intAr1 &lerro, int ithrd1, int ithrd2, int ithrd3);


template<class MFT>
int collversurf(Mesh<MFT>& msh, int iface, int iver, double qmax_suf, 
                MshCavity &cav, CavWrkArrs &work, 
                intAr1 &lerro, int ithrd1, int ithrd2){
  if(msh.nelem > 0) METRIS_THROW_MSG(TODOExcept(), "Implement + tet nelem = "<<msh.nelem)
  GETVDEPTH(msh.param);

  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  //METRIS_ASSERT(ithrd3 >= 0 && ithrd3 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);
  //METRIS_ASSERT(ithrd1 != ithrd3);
  //METRIS_ASSERT(ithrd2 != ithrd3);

  CavOprOpt  opts;
  CavOprInfo info;
  opts.allow_topological_correction = true; // To fetch missing edges
  opts.skip_topo_checks = false;
  opts.allow_remove_points = true;
  //opts.dryrun   = true;
  opts.dryrun   = false;
  opts.qmax_suf = qmax_suf;

  int ierro = 0;

  // Track best insertion dry run to run it at the end
  double qmabest = 1.0e30;
  int ivebest = -1; // p collapse vertex of 
  //int ipibest = -1; // p insert 

  // Collapse this one 
  int ipcol = msh.fac2poi(iface,iver);

  CPRINTF1("-- START colledgsurf iface = %d iver = %d ipoin = %d \n",iface,iver,ipcol);
  if(DOPRINTS2()) writeMesh("debug_collapse0.meshb",msh);


  int ibcol = msh.poi2bpo[ipcol];
  if(ibcol >= 0){
    if(msh.bpo2ibi(ibcol,1) == 0) return 1;
  }

  int iopen;
  ierro = ball(msh,ipcol,
               cav.lcedg,cav.lcfac,cav.lctet,
               &iopen,ithrd1);
  METRIS_ASSERT(ierro == 0);


  CPRINTF1(" - try collapse poi = %d ball nface = %d nedge = %d \n",
                              ipcol,cav.lcfac.get_n(),cav.lcedg.get_n());

  // Try the cavity call with different ipins in neighbours of ipcol
  msh.tag[ithrd1]++;
  int tag0 = msh.tag[ithrd1];
  // We'll stack on top of the ball, so we need to be able to prune to nbalf
  // with each attempt, as well as restrict search of ipins to ball 
  int nbalf = cav.lcfac.get_n();
  int nbale = cav.lcedg.get_n();
  for(int iicfc = 0; iicfc < nbalf; iicfc++){
    int icfac = cav.lcfac[iicfc]; 
    METRIS_ASSERT(icfac >= 0 && icfac < msh.nface);
    METRIS_ASSERT(!isdeadent(icfac,msh.fac2poi));

    // Doesn't change but easy to get it here 

    for(int ive2 = 0; ive2 < 3; ive2 ++){
      int ipins = msh.fac2poi(icfac,ive2);
      if(ipins == ipcol) continue;
      if(msh.poi2tag(ithrd1,ipins) >= tag0) continue;
      msh.poi2tag(ithrd1,ipins) = tag0;

      //if(msh.isboundary_edges() && cav.lcedg.get_n() > 0 
      //|| msh.isboundary_faces() && cav.lcfac.get_n() > 0){
      //  int ibpoi = msh.poi2bpo[ipins];
      //  if(ibpoi < 0) continue;
      //}

      cav.ipins = ipins;
      cav.lctet.set_n(0);
      cav.lcfac.set_n(nbalf); // Revert to simple ball 
      cav.lcedg.set_n(nbale); // Revert 
      CPRINTF1(" - try reinsert point %d tag = %d vs %d \n",
                           ipins,msh.poi2tag(ithrd1,ipins),tag0);


      if(DOPRINTS2()) writeMeshCavity("collapse_cavity0.meshb", msh, cav);

      ierro = increase_cavity(msh, cav, true, ithrd1, ithrd2);
      if(ierro > 0){
        CPRINTF1(" - failed to increase cavity ierro = %d \n",ierro);
        continue;
      }


      if(DOPRINTS2()) writeMeshCavity("collapse_cavity1.meshb", msh, cav);


      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
        ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd2);
      }}CT_FOR1(ideg);
      if(ierro > 0){
        lerro[ierro-1]++;
      }
      //printf("Debug wait;\n");
      // If operation was done, out
      if(info.done){
        CPRINTF1("-- END colledgsurf successful using ipcol = %d ipins = %d \n",ipcol,ipins);
        if(DOPRINTS2()) writeMesh("debug_collapse1.meshb",msh);
        msh.poi2ent(ipcol,0) = -1;
        msh.poi2ent(ipcol,1) = -1;
      }
      if(msh.param->dbgfull){
        check_topo(msh, msh.nbpoi, msh.npoin, msh.nedge, msh.nface, msh.nelem, ithrd2); 
      }
      if(info.done) return 0;

      CPRINTF1(" - return qmax = %f \n",info.qmax_end);
      if(info.qmax_end < qmabest && ierro == 0){
        CPRINTF1(" - new best quality !\n");
        qmabest = info.qmax_end;
        ivebest = iver;
      }
    }
  }

  if(ivebest == -1){
    // Failed to find a reinsertion candidate. 
    // Let's not rely on exceptions as failure here is not invalid behaviour 
    // and could happen quite frequently 
    return 1;
  }

  if(ivebest == 0){ // Otherwise ball hasn't changed !
    ierro = ball(msh, ipcol, cav.lcedg, cav.lcfac, cav.lctet,
                 &iopen, ithrd1);
    
    ierro = increase_cavity(msh, cav, true, ithrd1, ithrd2);
    METRIS_ASSERT(ierro == 0);
  }

  if(DOPRINTS2()) writeMeshCavity("collapse_cavity0.meshb", msh, cav);

  opts.dryrun = false;

  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd2);
  }}CT_FOR1(ideg);

  METRIS_ASSERT(ierro == 0);
  METRIS_ASSERT(abs(info.qmax_end - qmabest) < 1.0e-6);
  METRIS_ASSERT(info.done);

  return 0;
}



template int collversurf<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
                   int iface, int iver, double qmax_suf, 
                   MshCavity &cav, CavWrkArrs &work, 
                   intAr1 &lerro, int ithrd1, int ithrd2);
template int collversurf<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
                   int iface, int iver, double qmax_suf, 
                   MshCavity &cav, CavWrkArrs &work, 
                   intAr1 &lerro, int ithrd1, int ithrd2);



// Reject proposed cavity based on density changes
// Each point counts in the cavity as the ratio of the ball volume that is contained
// in the cavity

// Accelerate by using nentt sized work array
// Implement getmeasent for surface
template<class MFT>
int collrejcav_dens(Mesh<MFT>& msh, MshCavity &cav, int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);

  //printf("## DEBUG forced iverb = 5 ivdepth = 5\n");
  //iverb__ = 5;
  //ivdepth__ = 5;

  const int tdim = cav.lctet.get_n() > 0 ? 3 
                 : cav.lcfac.get_n() > 0 ? 2 
                                         : 1;
  METRIS_ASSERT(tdim != 1);

  METRIS_ENFORCE(tdim == msh.idim);

  const intAr1& lcent = cav.lcent(tdim);
  const int ncent = lcent.get_n();
  const intAr2& ent2poi = msh.ent2poi(tdim);
  const intAr2& ent2ent = msh.ent2ent(tdim);
  intAr2& ent2tag = msh.ent2tag(tdim);


  // To precompute element volumes
  cav.rwrk1.set_n(ncent);
  // To accumulate volumes per point
  dblWrkAr1 volpoc = msh.get_rwork(msh.npoin);

  double met[6];

  msh.tag[ithrd1]++;
  int etag = msh.tag[ithrd1];

  // Tag cavity elements, compute volumes and zero out their point volumes
  double voltot = 0;
  for(int icent = 0; icent < ncent; icent++){
    int ientt = lcent[icent];
    ent2tag(ithrd1,ientt) = etag;

    double volM;
    MSH_DIM_DEG0(msh)
    volM = getmeasent<MFT,gdim,ideg>(msh, ientt);
    MSH_DIM_DEG1()

    METRIS_ASSERT(volM > 0);
    cav.rwrk1[icent] = volM;

    voltot += volM;

    // zero out point volumes, we'll accumulate in a later loop
    for(int iver = 0; iver < tdim + 1; iver++){
      int ipoin = ent2poi(ientt, iver);
      volpoc[ipoin] = 0;
    }
  }

  // Tag points on cavity boundary, these are not internal 
  // Also compute their cavity-internal impinging volume
  for(int icent = 0; icent < ncent; icent++){
    int ientt = lcent[icent];
    // Only consider cav boundary facets
    for(int ifa = 0; ifa < tdim + 1; ifa++){
      int ienei = ent2ent(ientt, ifa);
      if(ienei >= 0 && ent2tag(ithrd1,ienei) >= etag) continue;
      // Accumulate volume at facet vertices
      for(int iver = 0; iver < tdim + 1; iver++){
        if(iver == ifa) continue;
        int ipoin = ent2poi(ientt, iver);
        if(msh.poi2tag(ithrd1, ipoin) == msh.tag[ithrd1]) continue;
        msh.poi2tag(ithrd1, ipoin) = msh.tag[ithrd1];
        volpoc[ipoin] += cav.rwrk1[icent];
      }// for iver
    }// for ifa
  }// for icent

  // Loop again, this time compute their full ball volumes and count contributions
  double nwpoi = 0;
  msh.tag[ithrd1]++;
  intAr1 dum;
  intAr1 &lbfac = tdim == 2 ? cav.iwrk1 : dum;
  intAr1 &lbtet = tdim == 3 ? cav.iwrk1 : dum;
  intAr1 &lbent = tdim == 2 ? lbfac : lbtet;
  cav.iwrk1.allocate(10);
  for(int ientt : lcent){
    INCVDEPTH(msh.param);
    // Only consider cav boundary facets
    for(int ifa = 0; ifa < tdim + 1; ifa++){
      int ienei = ent2ent(ientt, ifa);
      if(ienei >= 0 && ent2tag(ithrd1,ienei) >= etag) continue;
      // Accumulate volume at facet vertices
      for(int iver = 0; iver < tdim + 1; iver++){
        if(iver == ifa) continue;
        int ipoin = ent2poi(ientt, iver);
        if(msh.poi2tag(ithrd1, ipoin) == msh.tag[ithrd1]) continue;
        msh.poi2tag(ithrd1, ipoin) = msh.tag[ithrd1];

        // Compute the ball and its volume
        int iopen;
        int ierro = ball(msh, ipoin, dum, lbfac, lbtet, &iopen, ithrd2);
        METRIS_ASSERT(ierro == 0);
        if(ierro != 0) return 1+ierro;

        double volB = 0;
        for(int ientt : lbent){
          MSH_DIM_DEG0(msh)
          volB += getmeasent<MFT,gdim,ideg>(msh, ientt);
          MSH_DIM_DEG1()
        }
        CPRINTF1(" - cav bdry pt %d has tot vol %e, internal %e, counts as %e\n",
                 ipoin, volB, volpoc[ipoin], volpoc[ipoin] / volB);
        nwpoi += volpoc[ipoin] / volB;
      }// for iver
    }// for ifa
  }// for icent

  // Now count the internal points. Boundary points are tagged.
  int nrempt = 0;
  for(int icent = 0; icent < ncent; icent++){
    int ientt = lcent[icent];
    for(int iver = 0; iver < tdim + 1; iver++){
      int ipoin = ent2poi(ientt,iver);
      if(msh.poi2tag(ithrd1, ipoin) == msh.tag[ithrd1]) continue;
      msh.poi2tag(ithrd1, ipoin) = msh.tag[ithrd1];
      // An internal point, and one not seen yet.
      nrempt++;
    }
  }

  double dens0 = (nrempt + nwpoi) / voltot;
  double dens1 = nwpoi / voltot; 

  CPRINTF1(" - counted nrempt = %d boundary points %f\n",nrempt,nwpoi);
  CPRINTF1(" - initial cavity density = %e final = %e, nentt = %d\n",
           dens0, dens1, ncent);

  //if(nrempt > 5){
  //  writeMeshCavity("collapse_cavity0.meshb", msh, cav);
  //  printf("## DEBUG WAIT\n");
  //  wait();
  //}

  // If initial is closer to optimal density than final, return.
  double pi = 3.141592653589793238462643383279502884;
  double opt_dens = msh.get_tdim() == 2 ? pi / 4 : 0.54;
  if(abs(dens0 - opt_dens) < abs(dens1 - opt_dens)){
    CPRINTF1(" # %e > %e -> reject\n",abs(dens1 - opt_dens), abs(dens0 - opt_dens));
    return 1;
  }

  return 0;
}



template int collrejcav_dens<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, MshCavity &cav, int ithrd1, int ithrd2);
template int collrejcav_dens<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, MshCavity &cav, int ithrd1, int ithrd2);

#if 0
// This idea is probably doomed to fail: 2 long edges do not mean 2 new points...
// Reject proposed cavity based on edge length:
// if more long edges are created than short edges are destroyed, reject
// Return 1 if reject, 0 otherwise
template<class MFT>
int collrejcav_len(Mesh<MFT>& msh, MshCavity &cav, int ithrd1){


  GETVDEPTH(msh.param);

  printf("## DEBUG forced iverb = 5 ivdepth = 5\n");
  iverb__ = 5;
  ivdepth__ = 5;

  // Tag points that won't be deleted: there is at least one elt outside
  // the cavity that has the point. 
  int tdim = cav.lctet.get_n() > 0 ? 3 
           : cav.lcfac.get_n() > 0 ? 2 
                                   : 1;
  const intAr1& lcent = cav.lcent(tdim);
  const intAr2& ent2poi = msh.ent2poi(tdim);
  const intAr2& ent2ent = msh.ent2ent(tdim);
  intAr2& ent2tag = msh.ent2tag(tdim);

  // Store here the edges whose length is not to be computed
  std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp;

  // Tag cavity elements
  for(int ientt : lcent){
    ent2tag(ithrd1,ientt) = msh.tag[ithrd1];
  }


  // Start by adding all edges on the cavity boundary to the set. 
  // In 3D, this doesn't mean they're on a boundary face... so we need to separate
  //auto ledfa = tdim == 2 ? ledfa2 : ledfa3;
  //int nedfa = tdim == 2 ? 1 : 3;
  if(tdim == 3){
    for(int itetr : lcent){
      for(int ifa = 0; ifa < tdim + 1; ifa++){
        int itnei = msh.tet2tet(itetr,ifa);
        if(itnei >= 0 && msh.tet2tag(ithrd1,itnei) == msh.tag[ithrd1]) continue;
        // Edges on the cavity boundary:
        for(int iedf = 0; iedf < 3; iedf++){
          int ied = ledfa3[ifa][iedf];
          int ipoi1 = msh.tet2poi(itetr, lnoed3[ied][0]);
          int ipoi2 = msh.tet2poi(itetr, lnoed3[ied][1]);
          auto key = stup2(ipoi1, ipoi2);
          nocomp.insert(key);
        }// for ied
      }// for ifa
    }// for itetr
  }

  const int nedl = (tdim*(tdim+1))/2;
  int nshort0 = 0, nlong0 = 0;
  const auto lnoed = tdim == 2 ? lnoed2 : lnoed3;
  double len, sz[2];
  for(int ientt : lcent){
    for(int ied = 0; ied < nedl; ied++){
      int ipoi1 = ent2poi(ientt, lnoed[ied][0]);
      int ipoi2 = ent2poi(ientt, lnoed[ied][1]);

      // In this case, we haven't added to nocomp
      // Also seize opportunity to tag the points
      if(tdim == 2){
        int ifnei = msh.fac2fac(ientt,ied);
        if(ifnei < 0 || msh.fac2tag(ithrd1,ifnei) < msh.tag[ithrd1]) continue;
      }

      auto key = stup2(ipoi1, ipoi2);

      if(nocomp.find(key) != nocomp.end()) continue;

      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
        len = msh.idim == 2 ? 
          getlenedg_geosz<MFT,2,ideg>(msh,ientt,tdim,ied,sz) :
          getlenedg_geosz<MFT,3,ideg>(msh,ientt,tdim,ied,sz);
      }}CT_FOR1(ideg);

      if(len > sqrt(2)) nlong0++;
      if(len < 1/sqrt(2)) nshort0++;

      nocomp.insert(key);
    }// for ied
  }// for ientt

  // Compute nshort and nlong in final cavity
  int nshort1 = 0, nlong1 = 0;
  int edg2pol[2] = {cav.ipins, -1};
  for(int ientt : lcent){
    for(int ifa = 0; ifa < tdim + 1; ifa++){
      int ienei = ent2ent(ientt,ifa);
      if(ienei >= 0 && ent2tag(ithrd1,ienei) == msh.tag[ithrd1]) continue;
      // Get points on face (3D) / edge (2D)
      for(int ipfa = 0; ipfa < tdim; ipfa++){
        int ipoin = tdim == 2 ? lnoed2[ifa][ipfa] : lnofa3[ifa][ipfa];
        if(ipoin == cav.ipins) continue;
        if(msh.poi2tag(ithrd1, ipoin) == msh.tag[ithrd1]) continue;
        msh.poi2tag(ithrd1, ipoin) = msh.tag[ithrd1];

        edg2pol[1] = ipoin;
        CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
          len = msh.idim == 2 ? 
            getlenedg_geosz<MFT,2,ideg>(msh,edg2pol,sz) :
            getlenedg_geosz<MFT,3,ideg>(msh,edg2pol,sz);
        }}CT_FOR1(ideg);
        if(len > sqrt(2)) nlong1++;
        if(len < 1/sqrt(2)) nshort1++;

      }// for ipfa
    }// for ifa
  }// for ientt

  CPRINTF1(" - collrejcav_len got short %d -> %d, long %d -> %d\n",
           nshort0,nshort1,nlong0,nlong1);

  writeMeshCavity("collapse_cavity0.meshb", msh, cav);
  printf("Debug wait\n");
  wait();

  METRIS_ASSERT(nshort1 <= nshort0);
  METRIS_ASSERT(nlong1 >= nlong0);

  // More long are created than short are destroyed
  if(nlong1 - nlong0 >= nshort0 - nshort1) return 1;
  return 0;
}

template int collrejcav_len<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, MshCavity &cav, int ithrd1);
template int collrejcav_len<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, MshCavity &cav, int ithrd1);
#endif


} // end namespace

