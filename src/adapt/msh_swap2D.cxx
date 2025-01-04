//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../adapt/msh_swap2D.hxx"
#include "../adapt/low_swap2D.hxx"

#include "../Mesh/Mesh.hxx"

#include "../aux_topo.hxx"
#include "../io_libmeshb.hxx"
#include "../aux_timer.hxx"




namespace Metris{


// Greedy swaps: if a swap improves, do it
template<class MFT, int gdim, int ideg>
double swap2D(Mesh<MFT> &msh, swapOptions swapOpt, int ithrd1, int ithrd2){
  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);

  MetSpace ispac0 = msh.met.getSpace();
  msh.met.setSpace(MetSpace::Log);

  int iverb = msh.param->iverb;

  if(msh.get_tdim() != 2) METRIS_THROW_MSG(TODOExcept(), 
    "Implement collapseShortEdges on tdim != 2, got tdim = "<<msh.get_tdim());
  double stat = 0;

  int nswap = 0;
  int miter = 100, niter = 0;
  msh.tag[ithrd2]++;
  int ntry0 = -1;
  do{
    niter++;
    if(niter >= miter) break;
    
    bool onebad = false;
    int nerro = 0;
    nswap = 0;
    int ntry  = 0;
    int nfac0 = msh.nface;

    double t0 = get_wall_time();

    for(int iface = 0; iface < nfac0; iface++){
      if(isdeadent(iface,msh.fac2poi)) continue;
      if(msh.fac2tag(ithrd2,iface) >= msh.tag[ithrd2]) continue;

      int info; 
      int nfac1 = msh.nface;
      double qumx0,qumx1;
      #ifndef NDEBUG
        try{
      #endif
        info = swapface<MFT,gdim,ideg>(msh, iface, swapOpt, &qumx0, &qumx1, ithrd1);
      #ifndef NDEBUG
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
        msh.fac2tag(ithrd2,iface) = msh.tag[ithrd2];
      }else if(info > 0){ // Error 
        nerro ++;
      }else if(info < 0){ // Successful swap
        if(iverb >= 4) printf(" - swap successful\n");
        METRIS_ASSERT(nfac1 == msh.nface - 2);
        for(int ifanw = nfac1; ifanw < msh.nface; ifanw++){
          for(int ii = 0; ii < 3; ii++){
            int ifnei = msh.fac2fac(ifanw,ii);
            if(ifnei < 0) continue; // nm not eligible to swap w/ this either 
            // Unmark as inert if tagged 
            msh.fac2tag(ithrd2,ifnei) = msh.tag[ithrd2] - 1;
          }
        }
        nswap++;
        onebad = true;
      }
      ntry++; 
    }
    if(ntry0 < 0) ntry0 = ntry;
    stat = MAX(stat, (double)nswap / (double)ntry0);

    double t1 = get_wall_time();
    int ncallps = 1000*(int)((nswap / (t1-t0)) / 1000);
    if(iverb >= 2) printf("Loop end ntry = %d  nswap %d = %d /s; nerro %d coll \n",
                                                    ntry, nswap, ncallps,nerro);

    if(!onebad) break;

  }while(true);


  msh.met.setSpace(ispac0);
  
  return stat;
}


// ---------- Forward declarations
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double swap2D<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh,\
                              swapOptions swapOpt, int ithrd1, int ithrd2);\
template double swap2D<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh,\
                              swapOptions swapOpt, int ithrd1, int ithrd2);\
template double swap2D<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh,\
                              swapOptions swapOpt, int ithrd1, int ithrd2);\
template double swap2D<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh,\
                              swapOptions swapOpt, int ithrd1, int ithrd2);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()





} // end namespace