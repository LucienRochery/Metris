//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_swap2D.hxx"
#include "low_swap2D.hxx"

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
double swap2D(Mesh<MFT> &msh, swapOptions swapOpt, int *nswap, int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);

  if(msh.param->opt_swap_niter <= 0){
    CPRINTF1("-- Provided swap_niter == 0: skip swapping\n");
    return 0;
  }

  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);

  MetSpace ispac0 = msh.met.getSpace();
  msh.met.setSpace(MetSpace::Log);

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
  int ntry0 = -1;
  msh.tag[ithrd1]++;
  
  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh.param);
    
    bool onebad = false;
    int ntry  = 0;
    int nerro2 = 0;
    int nswap2 = 0;
    double t02 = get_wall_time();

    for(int tdim = 2; tdim <= msh.get_tdim(); tdim++){
      if(tdim == 3){
        static int nwarnprt1 = 0;
        if(nwarnprt1++ < 10) printf("\n\n## WARNING1: no 3D swaps\n\n\n");
        continue;
      }
      int nent0 = msh.nentt(tdim);
      intAr2& ent2tag = msh.ent2tag(tdim);

      int nerro1 = 0;
      int nswap1 = 0;

      double t01 = get_wall_time();

      for(int ientt = 0; ientt < nent0; ientt++){
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
          }else{
            info = 0;
            static int nwarnprt = 0;
            if(nwarnprt++ < 10) printf("\n\n## WARNING: no 3D swaps\n\n\n");
          }
        #ifndef NDEBUG
          if(msh.param->dbgfull)  check_topo(msh,ithrd2);
          }catch(MetrisExcept &e){
            printf("Caught exception in swapface, writing mesh:\n");
            writeMesh("swap_except.meshb",msh);
            msh.met.writeMetricFile("swap_except.solb");
            writeMesh("swap_except_back", *(msh.bak));
            msh.bak->met.writeMetricFile("swap_except_back.solb");
            std::string CADname = msh.param->outmPrefix + "swap_except_CAD.egads";
            EG_saveModel(msh.CAD.EGADS_model, CADname.c_str());
            std::cout<<"Wrote CAD file "<<CADname<<"\n";
            throw(e);
          }
        #endif

        //CPRINTF2(" - swap try entity %d info = %d \n",ientt, info);


        if(info == 0){ // Nothing done 
          // Tag entity as inert 
          ent2tag(ithrd1,ientt) = msh.tag[ithrd1];
        }else if(info > 0){ // Error 
          nerro1++;
        }else if(info < 0){ // Successful swap
          CPRINTF1(" - swap successful\n");
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
          nswap1++;
          onebad = true;
        }
        ntry++; 
      }
      double t11 = get_wall_time();
      int ncallps1 = 1000*(int)((nswap1 / (t11-t01)) / 1000);
      CPRINTF1(" - swaps dim %d swapped %d = %d /s nerro %d\n",
               tdim,nswap1,ncallps1,nerro1)
      nswap2 += nswap1;
      nerro2 += nerro1;
    }// for tdim

    if(ntry0 < 0)  ntry0 = ntry;
    if(ntry0 == 0) stat = MAX(stat,0);
    else           stat = MAX(stat, (double)nswap2 / (double)ntry0);

    double t12 = get_wall_time();
    int ncallps2 = 1000*(int)((nswap2 / (t12-t02)) / 1000);
    CPRINTF1(" - swaps full iter ntry = %d nswap %d = %d /s; nerro %d stat %f \n",
            ntry, nswap2, ncallps2,nerro2, (double)nswap2 / (double)ntry0);
    *nswap += nswap2;

    if(!onebad) break;
  }// for niter


  msh.met.setSpace(ispac0);
  
  return stat;
}


// ---------- Forward declarations
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double swap2D<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh,\
                              swapOptions swapOpt,int *nswap, int ithrd1, int ithrd2);\
template double swap2D<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh,\
                              swapOptions swapOpt,int *nswap, int ithrd1, int ithrd2);\
template double swap2D<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh,\
                              swapOptions swapOpt,int *nswap, int ithrd1, int ithrd2);\
template double swap2D<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh,\
                              swapOptions swapOpt,int *nswap, int ithrd1, int ithrd2);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()





} // end namespace