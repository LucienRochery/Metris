//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "low_cavqual.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../utils/mprintf.hxx"
#include "../cavity/msh_cavity.hxx"

#include "../low_topo.hxx"
#include "../low_lenedg.hxx"
#include "../low_geo.hxx"

#ifdef TRACY_ENABLE
#include "Tracy.hpp"
#endif

namespace Metris{


// Reject proposed cavity based on edge length score same as in swaps. 
// If filter_long/short is set, then the output edges that are long/short
// are ignored. In the case of an insertion, new long edges can be ignored. 
// Set grow_check to true to inspect elements outside the cavity. Only useful
// for insertions. 
// nocomp is work
template<class MFT>
int collrejcav_lenqua(Mesh<MFT>& msh, MshCavity &cav, 
                      bool filter_long, bool filter_short, bool grow_check, 
                      double lenqua_short_max, 
                      std::unordered_set<std::tuple<int,int>,tup2_hash::hash> &nocomp,
                      int ithrd1){

  #ifdef TRACY_ENABLE
  // Tracy profiler:
  ZoneScoped; 
  #endif

  METRIS_ASSERT(!filter_long || !filter_short);
  GETVDEPTH(msh.param);

  //printf("## DEBUG forced iverb = 5 ivdepth = 5\n");
  //iverb__ = 5;
  //ivdepth__ = 5;

  CPRINTF1("-- START collrejcav_lenqua filter long %d short %d grow %d\n",
           filter_long,filter_short,grow_check);


  //// Tag points that won't be deleted: there is at least one elt outside
  //// the cavity that has the point. 
  //int tdim = cav.lctet.get_n() > 0 ? 3 
  //         : cav.lcfac.get_n() > 0 ? 2 
  //                                 : 1;
  int tdim = MAX(1,msh.getpoitdim(cav.ipins));
  intAr1& lcent = cav.lcent(tdim);
  const intAr2& ent2poi = msh.ent2poi(tdim);
  const intAr2& ent2ent = msh.ent2ent(tdim);
  intAr2& ent2tag = msh.ent2tag(tdim);

  // Store here the edges whose length is not to be computed
  nocomp.clear();
  nocomp.reserve(lcent.get_n());

  msh.tag[ithrd1]++;
  // Tag cavity elements
  for(int ientt : lcent){
    ent2tag(ithrd1,ientt) = msh.tag[ithrd1];
  }

  const int nedl = (tdim*(tdim+1))/2;
  double qua0 = -1;
  const auto lnoed = tdim == 1 ? lnoed1 : 
                     tdim == 2 ? lnoed2 : lnoed3;
  double len, sz[2];

  {// tracy scope
  #ifdef TRACY_ENABLE
  ZoneScopedN("Initial cavity");
  #endif
  // Start by adding all edges on the cavity boundary to nocomp. 
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
          int ipoi1 = msh.tet2poi(itetr, lnoed[ied][0]);
          int ipoi2 = msh.tet2poi(itetr, lnoed[ied][1]);
          auto key = stup2(ipoi1, ipoi2);
          {
          #ifdef TRACY_ENABLE
          ZoneScopedN("nocomp insert 1");
          #endif
          nocomp.insert(key);
          }
        }// for ied
      }// for ifa
    }// for itetr
  }

  // Compute lengths of internal edges in initial cavity 
  for(int ientt : lcent){
    for(int ied = 0; ied < nedl; ied++){
      #ifdef TRACY_ENABLE
      ZoneScopedN("initial cav inner loop 2");
      #endif

      int ipoi1 = ent2poi(ientt, lnoed[ied][0]);
      int ipoi2 = ent2poi(ientt, lnoed[ied][1]);

      // In this case, we haven't added to nocomp
      // Also seize opportunity to tag the points
      if(tdim == 2){
        int ifnei = msh.fac2fac(ientt,ied);
        if(ifnei < 0 || msh.fac2tag(ithrd1,ifnei) < msh.tag[ithrd1]) continue;
      }

      auto key = stup2(ipoi1, ipoi2);
      {
      #ifdef TRACY_ENABLE
      ZoneScopedN("nocomp find 1");
      #endif
      if(nocomp.find(key) != nocomp.end()) continue;
      }   

      {
      #ifdef TRACY_ENABLE
      ZoneScopedN("initial cav lenedg");
      #endif
      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
        len = msh.idim == 2 ? 
          getlenedg_geosz<MFT,2,ideg>(msh,ientt,tdim,ied,sz) :
          getlenedg_geosz<MFT,3,ideg>(msh,ientt,tdim,ied,sz);
      }}CT_FOR1(ideg);
      }

      double quaed = len < 1.0 ? 1.0 - len 
                               : 1.0 - 1.0 / len;

      CPRINTF1(" - 0 orig len = %e score %e \n", len, quaed);
      //CPRINTF1(" met 1: %f %f %f met 2 : %f %f %f\n",msh.met(ipoi1,0),msh.met(ipoi1,1),msh.met(ipoi1,2)
      //  ,msh.met(ipoi2,0),msh.met(ipoi2,1),msh.met(ipoi2,2));
      qua0 = MAX(qua0, quaed);

      {
      #ifdef TRACY_ENABLE
      ZoneScopedN("nocomp insert 2");
      #endif
      nocomp.insert(key);
      }
    }// for ied
  }// for ientt
  }// tracy scope


  // Compute lengths of internal edges in final cavity 
  double qua1 = -1;
  int edg2pol[2] = {cav.ipins, -1};
  int ncent0 = lcent.get_n();
  int maxtag = msh.tag[ithrd1];
  // Potentially grow the cavity (and then revert) to catch short edges. 
  const int mgrow = 5;
  int icen0 = 0, icen1 = lcent.get_n();
  {
  #ifdef TRACY_ENABLE
  ZoneScopedN("Final cavity");
  #endif
  for(int ngrow = 0; ngrow < mgrow; ngrow++){


    int nadded = 0;
    for(int icent = icen0; icent < icen1; icent++){
      INCVDEPTH(msh.param);
      int ientt = lcent[icent];
      for(int ifa = 0; ifa < tdim + 1; ifa++){
        int ienei = ent2ent(ientt,ifa);
        // If no growth is allowed, ent2tag(ithrd1,ientt) is simply tag[ithrd1]
        // Otherwise, it is tag + generation. 
        if(ienei >= 0 && ent2tag(ithrd1,ienei) >= ent2tag(ithrd1,ientt)) continue;
        // Get points on face (3D) / edge (2D)
        for(int ipfa = 0; ipfa < tdim; ipfa++){
          int iver = tdim == 1 ? ipfa : 
                     tdim == 2 ? lnoed2[ifa][ipfa] : lnofa3[ifa][ipfa];
          int ipoin = ent2poi(ientt, iver);
          if(ipoin == cav.ipins) continue;
          if(msh.poi2tag(ithrd1, ipoin) == msh.tag[ithrd1]) continue;
          msh.poi2tag(ithrd1, ipoin) = msh.tag[ithrd1];

          edg2pol[1] = ipoin;
          //CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
            len = msh.idim == 2 ? 
              getlenedg_geosz<MFT,2,1>(msh,edg2pol,sz) :
              getlenedg_geosz<MFT,3,1>(msh,edg2pol,sz);
          //}}CT_FOR1(ideg);

          double quaed = len < 1.0 ? 1.0 - len 
                                   : 1.0 - 1.0 / len;

          CPRINTF1(" - %d new edg2pol = %d %d len = %e score %e \n",ngrow+1,edg2pol[0],edg2pol[1],
            len, quaed);
          if(len > sqrt(2) && filter_long){
            CPRINTF1(" -> skip long edge (filter_long)\n");
            continue;
          }else if(len < 1/sqrt(2) && filter_short){
            CPRINTF1(" -> skip short edge (filter_short)\n");
            continue;
          }
        //CPRINTF1(" met 1: %f %f %f met 2 : %f %f %f\n",msh.met(edg2pol[0],0),msh.met(edg2pol[0],1),msh.met(edg2pol[0],2)
        //  ,msh.met(edg2pol[1],0),msh.met(edg2pol[1],1),msh.met(edg2pol[1],2));
          qua1 = MAX(qua1, quaed);

        }// for ipfa


        // If cavity is allowed to grow, add the neighbour to the cavity.
        if(grow_check && ienei >= 0){
          // We tag these outside elements at tag + 1. 
          // Otherwise, we risk adding all the elements surrounding a cavity
          // boundary point as if it were to be collapsed, and then not check 
          // its length. 
          if(ent2tag(ithrd1,ienei) < msh.tag[ithrd1]){
            lcent.stack(ienei);
            ent2tag(ithrd1,ienei) = ent2tag(ithrd1,ientt) + 1;
            maxtag = MAX(maxtag, ent2tag(ithrd1, ienei));
            nadded++;
          }
        }
      }// for ifa
    }// for ientt
    icen0 = icen1;
    icen1 = lcent.get_n();
    // Also captures case flag is off
    CPRINTF1(" - cavity grown by %d new qua1 = %e\n",nadded,qua1);
    //if(abs(qua1 - qua1_prev) > 1.0e-6 && qua1_prev > 0){
    //  printf("## DEBUG AN UPDATE IN QUA1 %e -> %e iter %d\n",qua1_prev, qua1,ngrow);
    //  wait();
    //}
    if(nadded == 0) break;
  }// for ngrow

  }

  msh.tag[ithrd1] = maxtag;
  lcent.set_n(ncent0);

  CPRINTF1(" - END collrejcav_lenqua got lenqua %e -> %e\n", qua0, qua1);

  if(qua1 >= 0.99*qua0) return 1;
  if(filter_long && qua1 >= lenqua_short_max){
    CPRINTF1(" # reject due to qua1 >= lenqua_short_max: %e >= %e\n",qua1,lenqua_short_max);
    return 2;
  }
  return 0;
}

template int collrejcav_lenqua<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical>& msh,
 MshCavity &cav, bool filter1, bool filter2, bool igrow, double lenqua_short_max,
 std::unordered_set<std::tuple<int,int>,tup2_hash::hash> &nocomp,
 int ithrd1);
template int collrejcav_lenqua<MetricFieldFE        >(Mesh<MetricFieldFE        >& msh,
 MshCavity &cav, bool filter1, bool filter2, bool igrow, double lenqua_short_max,
 std::unordered_set<std::tuple<int,int>,tup2_hash::hash> &nocomp,
 int ithrd1);





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


  msh.tag[ithrd1]++;
  int etag = msh.tag[ithrd1];

  // Tag cavity elements, compute volumes and zero out their point volumes
  double voltot = 0;
  for(int icent = 0; icent < ncent; icent++){
    int ientt = lcent[icent];
    ent2tag(ithrd1,ientt) = etag;

    double volM;
    MSH_DIM_DEG0(msh){
    volM = getmeasent<MFT,gdim,ideg>(msh, ientt);
    }MSH_DIM_DEG1();

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
        int ierro = ball(msh, ipoin, dum, lbfac, lbtet, &iopen, false, ithrd2);
        METRIS_ASSERT(ierro == 0);
        if(ierro != 0) return 1+ierro;

        double volB = 0;
        for(int ientt : lbent){
          MSH_DIM_DEG0(msh){
          volB += getmeasent<MFT,gdim,ideg>(msh, ientt);
          }MSH_DIM_DEG1();
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




}// namespace
