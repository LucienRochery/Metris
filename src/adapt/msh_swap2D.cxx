//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_swap2D.hxx"
#include "low_swap2D.hxx"

#include "../Mesh/Mesh.hxx"

#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../io_libmeshb.hxx"
#include "../aux_timer.hxx"
#include "../mprintf.hxx"
#include "../aux_utils.hxx"
#include "../msh_checktopo.hxx"




namespace Metris{


// Greedy swaps: if a swap improves, do it
template<class MFT, int gdim, int ideg>
double swap2D(Mesh<MFT> &msh, swapOptions swapOpt, int *nswap, int ithrd1, int ithrd2){
  GETVDEPTH(msh);

  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);

  MetSpace ispac0 = msh.met.getSpace();
  msh.met.setSpace(MetSpace::Log);

  if(msh.get_tdim() != 2) METRIS_THROW_MSG(TODOExcept(), 
    "Implement collapseShortEdges on tdim != 2, got tdim = "<<msh.get_tdim());
  double stat = 0;

  int mcfac = 2; // more than 2 is a corner collapse: no!
  MshCavity cav(0,mcfac,0);
  CavWrkArrs work;

  *nswap = 0;
  int miter = 100;
  // Untaged elements are to be considered 
  #ifndef GLOFRO
  msh.tag[ithrd1]++;
  #endif
  int ntry0 = -1;
  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh);
    
    bool onebad = false;
    int nerro = 0;
    int nswap1 = 0;
    int ntry  = 0;
    int nfac0 = msh.nface;

    double t0 = get_wall_time();

    for(int iface = 0; iface < nfac0; iface++){
      INCVDEPTH(msh);
      if(isdeadent(iface,msh.fac2poi)) continue;
      if(msh.fac2tag(ithrd1,iface) >= msh.tag[ithrd1]) continue;

      int info; 
      int nfac1 = msh.nface;
      double qumx0,qumx1;
      #ifndef NDEBUG
        try{
      #endif
        info = swapface<MFT,gdim,ideg>(msh, iface, swapOpt, cav, work, &qumx0, &qumx1, ithrd2);
      #ifndef NDEBUG
        if(msh.param->dbgfull)  check_topo(msh);
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


      if(info == 0){ // Nothing done 
        // Tag face as inert 
        msh.fac2tag(ithrd1,iface) = msh.tag[ithrd1];
      }else if(info > 0){ // Error 
        nerro ++;
      }else if(info < 0){ // Successful swap
        CPRINTF1(" - swap successful\n");
        METRIS_ASSERT(nfac1 == msh.nface - 2);
        for(int ifanw = nfac1; ifanw < msh.nface; ifanw++){
          for(int ii = 0; ii < 3; ii++){
            int ifnei = msh.fac2fac(ifanw,ii);
            if(ifnei < 0) continue; // nm not eligible to swap w/ this either 
            // Unmark as inert if tagged 
            msh.fac2tag(ithrd1,ifnei) = msh.tag[ithrd1] - 1;
          }
        }
        nswap1++;
        onebad = true;
      }
      ntry++; 
    }
    if(ntry0 < 0) ntry0 = ntry;
    if(ntry0 == 0) stat = 0;
    else stat = MAX(stat, (double)nswap1 / (double)ntry0);

    double t1 = get_wall_time();
    int ncallps = 1000*(int)((nswap1 / (t1-t0)) / 1000);
    CPRINTF1("Loop end ntry = %d  nswap %d = %d /s; nerro %d coll \n",
                                                    ntry, nswap1, ncallps,nerro);
    *nswap += nswap1;

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