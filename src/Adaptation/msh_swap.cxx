//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_swap.hxx"
#include "low_swap.hxx"

#include "../Mesh/Mesh.hxx"

#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../low_topo.hxx"
#include "../io_libmeshb.hxx"
#include "../utils/aux_timer.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/aux_misc.hxx"
#include "../msh_checktopo.hxx"




namespace Metris{


// Greedy swaps: if a swap improves, do it
template<class MFT, int gdim, int ideg>
double swapMesh(Mesh<MFT> &msh, swapOptions swapOpt, int *nswap, int ithrd1, int ithrd2, int ithrd3, int ient0){
  GETVDEPTH(msh.param);

  if(msh.param->opt_swap_niter <= 0){
    CPRINTF1("-- Provided swap_niter == 0: skip swapping\n");
    return 0;
  }

  //printf("## DEBUG SET MAX PRINTS HERE \n");
  //wait();
  //msh.param->ivdepth = 5;
  //msh.param->iverb = 5;

  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd3 >= 0 && ithrd3 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);
  METRIS_ASSERT(ithrd1 != ithrd3);
  METRIS_ASSERT(ithrd2 != ithrd3);

  //MetSpace ispac0 = msh.met.getSpace();
  //msh.met.setSpace(MetSpace::Log);

  //if(msh.get_tdim() != 2){
  //  printf("## SKIP SWAP2D FOR TETRA \n");
  //  return 0;
  //}
  double stat = 0;


  MshCavity cav(64,2,0);
  CavWrkArrs work;
  intAr1 lshell(100);


  *nswap = 0;
  int miter = msh.param->opt_swap_niter;
  msh.tag[ithrd1]++;

  #ifdef TESTQUALITYALGO
  int nswapface = 0;
  int nswaptet = 0;
  #endif
  for(int tdim = 2; tdim <= msh.get_tdim(); tdim++){
    INCVDEPTH(msh.param);

    int nerro_tdim = 0;
    int nswap_tdim = 0;
    int nswap3fa = 0;
    int nswap3ed = 0;
    double t02 = get_cpu_time();

    for(int niter = 0; niter < miter; niter++){
      int   ntry  = 0;
      //if(tdim == 3){
      //  static int nwarnprt1 = 0;
      //  if(nwarnprt1++ < 10) printf("\n\n## WARNING1: no 3D swaps\n\n\n");
      //  continue;
      //}
      int nent0 = msh.nentt(tdim);
      intAr2& ent2tag = msh.ent2tag(tdim);

      int nerro_niter = 0;
      int nswap_niter = 0;

      double t01 = get_cpu_time();

      for(int ientt = ient0; ientt < nent0; ientt++){
        INCVDEPTH(msh.param);
        if(isdeadent(ientt,msh.ent2poi(tdim))) continue;
        //ntry++;

        if(ent2tag(ithrd1,ientt) >= msh.tag[ithrd1]) continue;

        int info;
        int nent1 = msh.nentt(tdim);
        double qumx0,qumx1;
        #ifndef NDEBUG
          try{
        #endif
          if(tdim == 2){
            info = swapface<MFT,gdim,ideg>(msh, ientt, swapOpt, cav, work, &qumx0, &qumx1, ithrd2);
            #ifdef TESTQUALITYALGO
            if (info < 0) nswapface++;
            #endif
          }else{
            if constexpr(gdim == 3){
              info = swaptetra<MFT,ideg>(msh, ientt, swapOpt, cav, work, &qumx0, &qumx1, ithrd2, ithrd3);
              //printf("## WAIT AFTER SWAPTETRA\n");
              //wait();
              #ifdef TESTQUALITYALGO
              if (info < 0) nswaptet++;
              #endif
            }else{
              METRIS_THROW_MSG("in dim < 3, ntetra > 0");
            }
          }
        #ifndef NDEBUG
          if(msh.param->dbgfull)  check_topo(msh,ithrd2);
          }catch(MetrisExcept &e){
            fmt::print(stderr,"Caught exception in swapface, writing mesh:\n");
            writeMesh("swap_except.meshb",msh);
            msh.met.writeMetricFile("swap_except.solb");
            writeMesh("swap_except_back", *(msh.bak));
            msh.bak->met.writeMetricFile("swap_except_back.solb");
            std::string CADname = msh.param->outmPrefix + "swap_except_CAD.egads";
            EG_saveModel(msh.CAD.EGADS_model, CADname.c_str());
            fmt::print(stderr,"Wrote CAD file {}\n",CADname);
            throw(e);
          }
        #endif

        //CPRINTF2(" - swap try entity {} info = {} \n",ientt, info);


        if(info == 0){ // Nothing done
          // Tag entity as inert
          ent2tag(ithrd1,ientt) = msh.tag[ithrd1];
        }else if(info > 0){ // Error
          nerro_niter++;
        }else if(info < 0){ // Successful swap
          CPRINTF2(" - swap successful\n");
          METRIS_ASSERT(nent1 == msh.nface - 2 || tdim == 3);
          if(tdim == 2){
            for(int ient1 = nent1; ient1 < msh.nentt(tdim); ient1++){
              for(int ii = 0; ii < 3; ii++){
                int ifnei = msh.fac2fac(ient1,ii);
                if(ifnei < 0) continue; // nm not eligible to swap w/ this either
                // Unmark as inert if tagged
                msh.fac2tag(ithrd1,ifnei) = msh.tag[ithrd1] - 1;
              }
            }
          }else{
            if(info == -1){
              nswap3ed++;
            }else{
              nswap3fa++;
            }
            // In case 3D, we need to use the shells.
            intAr1 dum1;
            for(int ient1 = nent1; ient1 < msh.nentt(tdim); ient1++){
              for(int iedgl = 0; iedgl < 6; iedgl++){
                int ip1 = msh.tet2poi(ient1,lnoed3[iedgl][0]);
                int ip2 = msh.tet2poi(ient1,lnoed3[iedgl][1]);
                int iopen;
                shell3(msh, ip1, ip2, ient1, lshell, dum1, &iopen);
                for(int ient2 : lshell)
                  msh.tet2tag(ithrd1,ient2) = msh.tag[ithrd1] - 1;
              }// for iedgl
            }// for ient1
          }// if tdim == 2
          nswap_niter++;
        }
        ntry++;
      }
      double t11 = get_cpu_time();

      double stat0 =(double)nswap_niter / (double)msh.nentt(tdim);
      if(ntry == 0) stat  = MAX(stat, 0);
      else          stat  = MAX(stat, stat0);

      int ncallps_niter = 1000*(int)((nswap_niter / (t11-t01)) / 1000);
      CPRINTF2(" - swaps full iter {} ntry = {} nswap {} = {} /s; nerro {} stat {:.2e}",niter,
              ntry, nswap_niter, ncallps_niter,nerro_niter, stat0);
      if(stat0 < msh.param->adp_stagn_stop && DOPRINTS2())
        fmt::print(LOGFILE__," < adp_stagn_stop = {} -> break.\n",msh.param->adp_stagn_stop);
      else if(DOPRINTS2()) fmt::print(LOGFILE__,"\n");
      nswap_tdim += nswap_niter;
      nerro_tdim += nerro_niter;

      if(stat0 < msh.param->adp_stagn_stop) break;

    }// for niter


    double t12 = get_cpu_time();
    int ncallps_tdim = 1000*(int)((nswap_tdim / (t12-t02)) / 1000);
    if(tdim == 2){
      CPRINTF2(" - swaps dim {} time {:.2e}s swapped {} = {} /s nerro {}\n",tdim,t12-t02,nswap_tdim,ncallps_tdim,nerro_tdim)
    }else{
      CPRINTF2(" - swaps dim {} time {:.2e}s swapped {} = {} /s nerro {} edge = {} face = {}\n",
               tdim,t12-t02,nswap_tdim,ncallps_tdim,nerro_tdim,nswap3ed,nswap3fa);
    }
    *nswap += nswap_tdim;

  }// for tdim


  //msh.met.setSpace(ispac0);

  #ifdef TESTQUALITYALGO
  std::cout << "End swapMesh" << std::endl;
  std::cout << "nswapface = " << nswapface << std::endl;
  std::cout << "nswaptet = " << nswaptet << std::endl;
  #endif

  return stat;
}


// ---------- Forward declarations
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double swapMesh<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh,\
                    swapOptions swapOpt,int *nswap, int ithrd1, int ithrd2, int ithrd3, int ient0);\
template double swapMesh<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh,\
                    swapOptions swapOpt,int *nswap, int ithrd1, int ithrd2, int ithrd3, int ient0);\
template double swapMesh<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh,\
                    swapOptions swapOpt,int *nswap, int ithrd1, int ithrd2, int ithrd3, int ient0);\
template double swapMesh<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh,\
                    swapOptions swapOpt,int *nswap, int ithrd1, int ithrd2, int ithrd3, int ient0);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()





} // end namespace