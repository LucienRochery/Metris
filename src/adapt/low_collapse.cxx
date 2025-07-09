//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "low_collapse.hxx"
#include "low_insert.hxx"
#include "low_cavqual.hxx"

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
int collapseEdge(Mesh<MFT>& msh, int tdim, int ientt, int iedl, double qmax_suf, 
                 MshCavity &cav, CavWrkArrs &work, 
                 intAr1 &lerro, int ithrd1, int ithrd2, int ithrd3){
  return insertEdge(msh, tdim, ientt, iedl, -1, true, cav, work, lerro, ithrd1, ithrd2);
}

template int collapseEdge<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
                                          int tdim, int ientt, int iedl, double qmax_suf, 
                                          MshCavity &cav, CavWrkArrs &work, 
                                          intAr1 &lerro, int ithrd1, int ithrd2, int ithrd3);
template int collapseEdge<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
                                          int tdim, int ientt, int iedl, double qmax_suf, 
                                          MshCavity &cav, CavWrkArrs &work, 
                                          intAr1 &lerro, int ithrd1, int ithrd2, int ithrd3);




/*
Collapse a vertex in short edge 
Cavity is passed in to reuse allocations
*/
template<class MFT>
int collapseEdge2(Mesh<MFT>& msh, int tdim, int ientt, int iedl, double qmax_suf, 
                  MshCavity &cav, CavWrkArrs &work, 
                  intAr1 &lerro, int ithrd1, int ithrd2, int ithrd3){
  METRIS_ASSERT(ientt >= 0);
  GETVDEPTH(msh.param);


  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  //METRIS_ASSERT(ithrd3 >= 0 && ithrd3 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);
  //METRIS_ASSERT(ithrd1 != ithrd3);
  //METRIS_ASSERT(ithrd2 != ithrd3);

  intAr2& ent2poi = msh.ent2poi(tdim);

 
  auto lnoed = tdim == 2 ? lnoed2 : lnoed3;

  int ip1 = ent2poi(ientt,lnoed[iedl][0]);
  int ip2 = ent2poi(ientt,lnoed[iedl][1]);


  CPRINTF1("-- START collapseEdge ientt = %d tdim %d iedl = %d = (%d,%d) \n",
           ientt,tdim,iedl,ip1,ip2);
  if(DOPRINTS2()) writeMesh("debug_collapse0.meshb",msh);


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

    int ierro = collapseVertex<MFT>(msh, ipcol, qmax_suf, cav, work, lerro, ithrd1, ithrd2);
    if(ierro <= 0) return 0;

    #if 0

    ierro = ball(msh, ipcol, cav.lcedg, cav.lcfac, cav.lctet,
                 &iopen, ithrd1);
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

        ierro = collrejcav_lenqua(msh, cav, false, false, false, -1, nocomp, ithrd2);
        if(ierro > 0){
          CPRINTF1(" # reject cavity\n");
          continue;
        }

        CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
          ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd2);
        }}CT_FOR1(ideg);

        if(ierro > 0){
          lerro[ierro-1]++;
        }
        //printf("Debug wait;\n");
        // If operation was done, out
        if(info.done){
          CPRINTF1("-- END collapseEdge successful using ipcol = %d ipins = %d \n",ipcol,ipins);
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
      }// for ive2
    }// for icent
    #endif
  }// for iver


  return 1;
}

template int collapseEdge2<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
                                          int tdim, int ientt, int iedl, double qmax_suf, 
                                          MshCavity &cav, CavWrkArrs &work, 
                                          intAr1 &lerro, int ithrd1, int ithrd2, int ithrd3);
template int collapseEdge2<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
                                          int tdim, int ientt, int iedl, double qmax_suf, 
                                          MshCavity &cav, CavWrkArrs &work, 
                                          intAr1 &lerro, int ithrd1, int ithrd2, int ithrd3);


template<class MFT>
int collapseVertex(Mesh<MFT>& msh, int ipcol, double qmax_suf, 
                   MshCavity &cav, CavWrkArrs &work, 
                   intAr1 &lerro, int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);

  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);

  cav.reset();

  CavOprOpt  opts;
  CavOprInfo info;
  opts.allow_topological_correction = true; // To fetch missing edges
  opts.skip_topo_checks = false;
  opts.allow_remove_points = true;
  opts.dryrun   = false;
  opts.qmax_suf = qmax_suf;

  int ierro = 0;

  // work for collrejcav_lenqua
  std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp;

  int tdimp = msh.getpoitdim(ipcol);
  if(tdimp == 0) return 1;

  CPRINTF1("-- START collapseVertex ipcol = %d tdimp = %d\n",ipcol,tdimp);


  int iopen;
  ierro = ball(msh, ipcol, cav.lcedg, cav.lcfac, cav.lctet,
               &iopen, false, ithrd1);
  CPRINTF1(" - try collapse poi = %d ball nface = %d nedge = %d \n",
                              ipcol,cav.lcfac.get_n(),cav.lcedg.get_n());
  METRIS_ASSERT(ierro == 0);


  // Try the cavity call with different ipins in neighbours of ipcol
  msh.tag[ithrd1]++;
  int tag0 = msh.tag[ithrd1]; // We'll reuse the tag for elements in subroutines
  // We'll stack on top of the ball, so we need to be able to prune to nbalf
  // with each attempt, as well as restrict search of ipins to ball 
  int nbalt = cav.lctet.get_n();
  int nbalf = cav.lcfac.get_n();
  int nbale = cav.lcedg.get_n();
  int tdim = msh.get_tdim();
  intAr2& ent2poi = msh.ent2poi(tdim);
  intAr1& lcent = cav.lcent(tdim);
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
      if(msh.getpoitdim(ipins) > tdimp){
        CPRINTF1(" - point %d dim %d > ipins dim %d -> reject reconnection\n",
          ipins,msh.getpoitdim(ipins),tdimp);
        continue;
      }

      // Check ipins has same ref in topo dim of ipcol.
      // Counter-example: a boundary point, ball's boundary hits other surface refs. 
      // This does not happen in the volume. Indeed, if a point's ball
      // has two domain refs, then the point itself has two domain refs.
      // Hence it was actually a boundary point. 
      if(tdimp < msh.get_tdim()){
        // As we never collapse a corner, the point's tdim is >= 1.
        // Hence it has a unique ref to the lowest dim entity group. 
        METRIS_ASSERT(cav.lcent(tdimp).get_n() > 0);
        int iref = msh.ent2ref(tdimp)[cav.lcent(tdimp)[0]];
        #ifndef NDEBUG
        for(int ient1 : cav.lcent(tdimp)){
          METRIS_ASSERT(iref == msh.ent2ref(tdimp)[ient1]);
        }
        #endif

        if(msh.getpoitdim(ipins) == tdimp){
          // Easy, just check seed reference
          int ient1 = msh.poi2ent(ipins,0);
          METRIS_ASSERT(msh.poi2ent(ipins,1) == tdimp);
          int iref1 = msh.ent2ref(tdimp)[ient1];
          if(iref1 != iref){
            CPRINTF1(" - point %d dim %d = dim ipins but ref %d != %d\n",
                     ipins,tdimp,iref1,iref);
            continue;
          }
        }else{
          // Use boundary info
          bool ifnd = false;
          for(int ibpoi = msh.poi2bpo[ipins]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
            int tdim1 = msh.bpo2ibi(ibpoi,1);
            if(tdim1 != tdimp) continue;
            int ient1 = msh.bpo2ibi(ibpoi,2);
            int iref1 = msh.ent2ref(tdimp)[ient1];
            CPRINTF1(" - check entity %d dim %d ref %d ipcol ref = %d\n",
                     ient1,tdim1,iref1,iref);
            if(iref1 != iref) continue;
            CPRINTF1(" -> found ref\n");
            ifnd = true;
            break;
          }
          if(!ifnd){
            CPRINTF1(" - did not find ref %d dim %d of ipcol in ipins %d refs\n",
                     iref,tdimp,ipins);
            continue;
          }
        }
      }// if tdimp < msh.get_tdim()

      cav.ipins = ipins;
      cav.lctet.set_n(nbalt);
      cav.lcfac.set_n(nbalf); // Revert to simple ball 
      cav.lcedg.set_n(nbale); // Revert 
      CPRINTF1(" - try reinsert point %d tag = %d vs %d \n",
                           ipins,msh.poi2tag(ithrd1,ipins),tag0);

      if(DOPRINTS2()) writeMeshCavity("collapse_cavity0.meshb", msh, cav);

      ierro = collrejcav_lenqua(msh, cav, false, false, false, -1, nocomp, ithrd2);
      if(ierro > 0){
        CPRINTF1(" # reject cavity\n");
        continue;
      }

      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
        ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithrd2);
      }}CT_FOR1(ideg);

      if(ierro > 0){
        lerro[ierro-1]++;
      }
      //printf("Debug wait;\n");
      // If operation was done, out
      if(info.done){
        CPRINTF1("-- END collapseEdge successful using ipcol = %d ipins = %d \n",ipcol,ipins);
        if(DOPRINTS2()) writeMesh("debug_collapse1.meshb",msh);
        msh.killpoint(ipcol);
      }
      if(info.done) return 0;

      CPRINTF1(" - cavity return qmax = %f \n",info.qmax_end);
    }// for ive2
  }// for icent
  return 1;
}



template int collapseVertex<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh, 
                   int ipcol, double qmax_suf, 
                   MshCavity &cav, CavWrkArrs &work, 
                   intAr1 &lerro, int ithrd1, int ithrd2);
template int collapseVertex<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh, 
                   int ipcol, double qmax_suf, 
                   MshCavity &cav, CavWrkArrs &work, 
                   intAr1 &lerro, int ithrd1, int ithrd2);



} // end namespace

