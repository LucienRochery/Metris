//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_insert.hxx"
#include "low_insert.hxx"

#include "../low_geo/lenedg.hxx"
#include "../aux_topo.hxx"
#include "../io_libmeshb.hxx"
#include "../utils/aux_timer.hxx"
#include "../utils/mprintf.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../adapt/msh_swap.hxx"
#include "../BezierOffsets/low_gaps.hxx"
#include "../low_geo/misc.hxx"
#include "../Mesh/Mesh.hxx"
#include "../msh_checktopo.hxx"
#include "../aux_histogram.hxx"
#include "../msh_lenedg.hxx"

#include <cmath>


namespace Metris{

// This version tags elements and also edges via a hash table. 
// Initially, all elements are untagged (active) and the hash table is empty. 
// When an insertion fails on an edge, the edge is added to the hashtable. 
// When insertions fail on all edges of an element, the element becomes tagged (inactive)
// When an insertion is carried out, neighbouring edges are not untagged because
// the majority of rejections are due to error INS2D_ERR_SHORTEDG which is only
// made worse by an insertion. 
// Edges eliminated by other errors can be removed from the hashtable in a second
// time. 
// Note their errors might depend on greater than 1 neighbourhood; hence there 
// is no perfect solution other than to allow a full run at some point (in the end). 

// lpins work array sized dynamically by this routine ; it's an argument solely because this will be called several times, save on alloc
// also: as iterations go, fewer and fewer edges are long, no use allocating more than once to maximum needed size (first iter)

// insertLongEdges is called.
template<class MFT, int gdim, int ideg>
double insertLongEdges(Mesh<MFT> &msh, int *ninser, int ithrd1, int ithrd2){
  //printf("## DEBUG SET MAX PRINTS\n");
  //wait();
  //msh.param->iverb = 5;
  //msh.param->ivdepth = 15;
  GETVDEPTH(msh.param);
  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);


  // Swap norm -1: length-based. 
  //swapOptions swapOpt(100, -1, 0.0);
  swapOptions swapOpt(100, 0, 0.005);

  //msh.met.setSpace(MetSpace::Log);
  msh.met.setSpace(MetSpace::Exp);

  CPRINTF2("-- insertLongEdges start \n");
  #ifndef NDEBUG
  CPRINTF2(" - Note: improve by generating several points per edge. Generated but not used cf loop nn/2 \n");
  CPRINTF2(" - Note: improve by filtering point propositions \n");
  #endif

  lenStat lenstat0;
  intAr2 ilned;
  dblAr1 rlned;
  dblAr1 lenbds = {1.0/sqrt(2), sqrt(2)};
  getLengthEdges<MFT>(msh,msh.get_tdim(),-1,ilned,rlned,lenstat0);
  CPRINTF1(" - Length qua short = %e\n",lenstat0.qua_short);
  CPRINTF1(" -            long  = %e\n",lenstat0.qua_long);
  double lenqua_short_max = (lenstat0.qua_short + lenstat0.qua_long)/2;
  CPRINTF1(" - pct unit %e using qua threshold %e\n",lenstat0.prop_unit*100,lenqua_short_max);


  double stat = 0;

  const int tdim = msh.get_tdim();

  const int nedgl = (tdim*(tdim+1))/2;

  const intAr2 lnoed(nedgl,2,tdim == 1 ? lnoed1[0] :
                             tdim == 2 ? lnoed2[0] : lnoed3[0]);
  const intAr2 &ent2poi = msh.ent2poi(tdim);

  // Error counting:
  const int mcaverr = CAV_ERR_NERROR;
  intAr1 lcaverr(mcaverr);
  const int minserr = 100;
  intAr1 linserr(minserr);

  // Outer loop iterations:
  const int miter = 10;

  int ierro;
  *ninser = 0;

  MshCavity cav(100,100,1);
  CavWrkArrs work;

  // This can be improvied, but it'll do for now:
  const int medge = tdim == 2 ? msh.nface : msh.nelem;
  HshTab_I2I ledge;
  ledge.reserve(medge);

  double t0s = get_wall_time();
  int nedge_tot = 0;
  for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
    INCVDEPTH(msh.param);
    if(isdeadent(ientt,ent2poi)) continue;
    for(int ied = 0; ied < nedgl; ied++){
      INCVDEPTH(msh.param);

      // Check edge already seen
      int ip1 = ent2poi(ientt, lnoed(ied,0));
      int ip2 = ent2poi(ientt, lnoed(ied,1));
      auto key = stup2(ip1,ip2);
      if(ledge.find(key) != ledge.end()) continue;

      nedge_tot++;

      double sz[2];
      double len = getlenedg_geosz<MFT,gdim,ideg>(msh,ientt,tdim,ied,sz);
      if(len < sqrt(2)) continue;
      ledge[key] = ientt;

      // New long edge, add to stack
      CPRINTF1(" + long edge %d %d seed %d len = %e\n",ip1,ip2,ientt,len);

    }// for ied
  }// for ientt
  double t1s = get_wall_time();
  printf(" - init time %f nlong = %d\n",t1s-t0s,(int)ledge.size());

  if(ledge.size() == 0){
    CPRINTF1(" - END insertLongEdges: no long edges\n");
    return 0;
  }


  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh.param);
    int nlong = ledge.size();
    CPRINTF1(" - START ins loop %d/%d nlong = %d\n",niter+1,miter,nlong);

    int nskip = 0, ntry = 0, ninser1 = 0, nerro = 0, nadded = 0;
    lcaverr.fill(0);
    linserr.fill(0);

    // Don't let the loop grow beyond initial size.
    const auto edge_it1 = ledge.end();
    double t0 = get_wall_time();
    for(auto edge_it = ledge.begin(); edge_it != edge_it1;){
      ntry++;
      INCVDEPTH(msh.param);
      std::tuple<int,int> key = edge_it->first;
      int ientt = edge_it->second;
      METRIS_ASSERT(ientt >= 0);
      if(isdeadent(ientt,ent2poi)){
        CPRINTF1(" - ientt %d dead, skip\n",ientt);
        nskip++;
        // Remove the edge from the hash table.
        edge_it = ledge.erase(edge_it);
        continue;
      }
      ntry++;
      int ip1 = std::get<0>(key);
      int ip2 = std::get<1>(key);
      int ied = getedgent(msh, tdim, ientt, ip1, ip2);
      METRIS_ASSERT(ied >= 0);

      CPRINTF1(" - enact ins ientt = %d ied = %d edg %d %d\n",ientt,ied,ip1,ip2);
      int nent0 = msh.nentt(tdim); 
      ierro = insertEdge(msh,tdim,ientt,ied,lenqua_short_max,false,
                         cav,work,lcaverr,ithrd1,ithrd2);

      if(ierro <= 0){
        ninser1++;
        // constrain point
        msh.poicstr[cav.ipins] = true;
        // Remove the edge from the edge hash table.
        edge_it = ledge.erase(edge_it);

        //// Try to swap new elements
        int nent1 = msh.nentt(tdim);
        //MSH_DIM_DEG0(msh){
        //  swapMesh<MFT,gdim,ideg>(msh,opt_swap,&nswap,ithrd1,ithrd2,ithrd3,nent0);
        //}MSH_DIM_DEG1();





        // Look for long edges in new elements.
        nent1 = msh.nentt(tdim);
        for(int ientn = nent0; ientn < nent1; ientn++){
          INCVDEPTH(msh.param);

          for(int iedn = 0; iedn < nedgl; iedn++){

            // If the edge exists and points to a dead element, we must update the seed
            // If it points to a live element, there is no harm in updating the seed
            // Hence we always update the seed if the edge exists.
            int jp1 = ent2poi(ientn, lnoed(iedn,0));
            int jp2 = ent2poi(ientn, lnoed(iedn,1));
            auto key = stup2(jp1,jp2);
            if(ledge.find(key) != ledge.end()){
              ledge[key] = ientn;
              continue;
            }

            double sz[2];
            double len = getlenedg_geosz<MFT,gdim,ideg>(msh,ientn,tdim,iedn,sz);
            if(len < sqrt(2)) continue;

            nadded++;

            // New long edge, add to stack
            CPRINTF1(" - + long edge %d %d seed %d len = %e\n",jp1,jp2,ientn,len);
            ledge[key] = ientn;
          }// for iedn
        }// for ientn

      }else{
        // Remove the edge from the edge hash table.
        edge_it = ledge.erase(edge_it);
        //edge_it++;
        CPRINTF2(" # insertion failed ierro = %d \n",ierro);
        linserr[ierro - 1] ++;
        nerro++;
      }
    }// for edge_it

    double stat0 = ninser1 / (double) nedge_tot;

    double t1 = get_wall_time();
    int ncallps = 1000*(int)((ninser1 / (t1-t0)) / 1000);
    CPRINTF2(" - END t = %f nlong %d ntry %d nskip %d nadded %d ninser %d = %d /s; nerro %d stat %f\n",
              t1-t0,nlong,ntry,nskip,nadded,ninser1,ncallps,nerro,stat0);
    if(DOPRINTS2() && nerro > 0){
      CPRINTF2(" - cavity ierro list:\n");
      for(int ii = 0; ii < mcaverr; ii++){
        if(lcaverr[ii] == 0) continue;
        CPRINTF2("   ierro = %d : %d \n",ii+1,lcaverr[ii]);
      }
      CPRINTF2(" - inspoi ierro list:\n");
      for(int ii = 0; ii < minserr; ii++){
        if(linserr[ii] == 0) continue;
        CPRINTF2("   ierro = %d : %d \n",ii+1,linserr[ii]);
      }
    }

    *ninser += ninser1;



    if(stat0 < msh.param->adp_stagn_stop) break;

    if(DOPRINTS3()){
      writeMesh("insert"+std::to_string(niter),msh);
      msh.met.writeMetricFile("insert"+std::to_string(niter));
    }

    //printf("## WAIT HERE \n");
    //wait();

    #ifndef NDEBUG
    check_topo(msh,ithrd1);
    #endif

  }// for niter
  msh.cleanup();
  stat = (double) *ninser /(double) nedge_tot;

  return stat;
}


#define BOOST_PP_LOCAL_MACRO(n)\
template double insertLongEdges<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh,\
                                int* ninser, int ithrd1, int ithrd2);\
template double insertLongEdges<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh,\
                                int* ninser, int ithrd1, int ithrd2);\
template double insertLongEdges<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh,\
                                int* ninser, int ithrd1, int ithrd2);\
template double insertLongEdges<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh,\
                                int* ninser, int ithrd1, int ithrd2);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()





} // end namespace