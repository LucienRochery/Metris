//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "lenedg.hxx"
#include "normal.hxx"
#include "misc.hxx"

#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"

#include "../Mesh/MeshMetric.hxx"
#include "../linalg/symidx.hxx"
#include "../linalg/matprods.hxx"


#ifdef TRACY_ENABLE
#include "Tracy.hpp"
#endif


namespace Metris{

// -----------------------------------------------------------------------------
template<int gdim>
double getlenedg(const double x1[], const double x2[], const double metl[]){
  if constexpr(gdim == 2){
    return sqrt((x1[0] - x2[0]) * (x1[0] - x2[0]) * metl[0]
              + (x1[0] - x2[0]) * (x1[1] - x2[1]) * metl[1] * 2
              + (x1[1] - x2[1]) * (x1[1] - x2[1]) * metl[2]);
  }else if(gdim == 3){
    return sqrt((x1[0] - x2[0]) * (x1[0] - x2[0]) * metl[0]
              + (x1[0] - x2[0]) * (x1[1] - x2[1]) * metl[1] * 2
              + (x1[1] - x2[1]) * (x1[1] - x2[1]) * metl[2]
              + (x1[0] - x2[0]) * (x1[2] - x2[2]) * metl[3] * 2
              + (x1[1] - x2[1]) * (x1[2] - x2[2]) * metl[4] * 2
              + (x1[2] - x2[2]) * (x1[2] - x2[2]) * metl[5]);
  }
}
template double getlenedg<2>(const double x1[], const double x2[], const double metl[]);
template double getlenedg<3>(const double x1[], const double x2[], const double metl[]);

// -----------------------------------------------------------------------------
template<int gdim>
double getlenedg(const double dx[], const double metl[]){
  if constexpr(gdim == 2){
    return sqrt(dx[0] * dx[0] * metl[0]
              + dx[0] * dx[1] * metl[1] * 2
              + dx[1] * dx[1] * metl[2]);
  }else if(gdim == 3){
    return sqrt(dx[0] * dx[0] * metl[0]
              + dx[0] * dx[1] * metl[1] * 2
              + dx[1] * dx[1] * metl[2]
              + dx[0] * dx[2] * metl[3] * 2
              + dx[1] * dx[2] * metl[4] * 2
              + dx[2] * dx[2] * metl[5]);
  }
}
template double getlenedg<2>(const double dx[], const double metl[]);
template double getlenedg<3>(const double dx[], const double metl[]);

// -----------------------------------------------------------------------------
template<int gdim>
double getlenedgsq(const double x1[], const double x2[], const double metl[]){
  if constexpr(gdim == 2){
    return (x1[0] - x2[0]) * (x1[0] - x2[0]) * metl[0]
         + (x1[0] - x2[0]) * (x1[1] - x2[1]) * metl[1] * 2
         + (x1[1] - x2[1]) * (x1[1] - x2[1]) * metl[2];
  }else if(gdim == 3){
    return (x1[0] - x2[0]) * (x1[0] - x2[0]) * metl[0]
         + (x1[0] - x2[0]) * (x1[1] - x2[1]) * metl[1] * 2
         + (x1[1] - x2[1]) * (x1[1] - x2[1]) * metl[2]
         + (x1[0] - x2[0]) * (x1[2] - x2[2]) * metl[3] * 2
         + (x1[1] - x2[1]) * (x1[2] - x2[2]) * metl[4] * 2
         + (x1[2] - x2[2]) * (x1[2] - x2[2]) * metl[5];
  }
}
template double getlenedgsq<2>(const double x1[], const double x2[], const double metl[]);
template double getlenedgsq<3>(const double x1[], const double x2[], const double metl[]);
// -----------------------------------------------------------------------------
template<int gdim>
double getlenedgsq(const double dx[], const double metl[]){
  if constexpr(gdim == 2){
    return dx[0] * dx[0] * metl[0]
         + dx[0] * dx[1] * metl[1] * 2
         + dx[1] * dx[1] * metl[2];
  }else if(gdim == 3){
    return dx[0] * dx[0] * metl[0]
         + dx[0] * dx[1] * metl[1] * 2
         + dx[1] * dx[1] * metl[2]
         + dx[0] * dx[2] * metl[3] * 2
         + dx[1] * dx[2] * metl[4] * 2
         + dx[2] * dx[2] * metl[5];
  }
}
template double getlenedgsq<2>(const double dx[], const double metl[]);
template double getlenedgsq<3>(const double dx[], const double metl[]);

// -----------------------------------------------------------------------------
// The metric is in log format.
template<int gdim>
double getlenedg_log(const double dx[], const double metl[], int miter, double tol){
  constexpr int nnmet = (gdim*(gdim+1))/2;
  double metacc[2][nnmet];

  // Start with identity
  double len = getnrml2<gdim>(dx);
  // + met
  len += getlenedgsq<gdim>(dx,metl);

  int iwhich = 0;
  for(int ii = 0; ii < gdim; ii++)
    for(int jj = ii; jj < gdim; jj++)
      metacc[iwhich][sym2idx(ii,jj)] = metl[sym2idx(ii,jj)];

  int niter = 2;
  double err;
  do{
    symXsymsub_fac<gdim>(metl,metacc[iwhich],1.0 / niter,metacc[1-iwhich]);
    // Accumulates into metacc[1 - iwhich]
    iwhich = 1 - iwhich;

    err = getlenedgsq<gdim>(dx,metacc[iwhich]);
    len += err;

    niter++;
    // length is squared
  }while(err*err > tol*tol*len && niter < miter);

  return sqrt(len);
}

template double getlenedg_log<2>(const double dx[], const double metl[],int miter, double tol);
template double getlenedg_log<3>(const double dx[], const double metl[],int miter, double tol);


// Geometric size interpolation 
// metric given in metSpac format
template<class MetricFieldType, int gdim, int ideg>
double getlenedg_geosz(const MeshMetric<MetricFieldType> &msh,
                       int ientt, int tdimn, int iedg){
  #ifdef TRACY_ENABLE
  ZoneScopedN("getlenedg_geosz2");
  #endif
  double sz[2];
  double len = getlenedg_geosz<MetricFieldType,gdim,ideg>(msh,ientt,tdimn,iedg,sz);
  return len;
}


// Same but also return the sizes (e.g. for insertion)
template<class MetricFieldType, int gdim, int ideg>
double getlenedg_geosz(const MeshMetric<MetricFieldType> &msh,
                       int ientt, int tdimn, int iedg, 
                       double *sz){
  #ifdef TRACY_ENABLE
  ZoneScopedN("getlenedg_geosz1");
  #endif
  const int nedgl = (tdimn*(tdimn+1))/2;

  int lnoed1[1][2] = {{0, 1}};
  const intAr2 lnoed(nedgl,2,tdimn == 1 ? lnoed1[0] :
                             tdimn == 2 ? lnoed2[0] : lnoed3[0]);

  const intAr2 &ent2poi = tdimn == 1 ? msh.edg2poi : 
                          tdimn == 2 ? msh.fac2poi : msh.tet2poi;

  int edg2pol[getnnod1(ideg)];

  edg2pol[0] = ent2poi(ientt,lnoed[iedg][0]);
  edg2pol[1] = ent2poi(ientt,lnoed[iedg][1]);
  int idx0 = tdimn + 1 + iedg*(ideg-1);
  for(int ii = 0; ii < ideg-1; ii++){
    edg2pol[2+ii] = ent2poi[ientt][idx0+ii];
  }


  return getlenedg_geosz<MetricFieldType,gdim,ideg>(msh,&edg2pol[0],sz);
}


// Same but also return the sizes (e.g. for insertion)
template<class MetricFieldType, int gdim, int ideg>
double getlenedg_geosz(const MeshMetric<MetricFieldType> &msh,const int *edg2pol, double *sz){
  #ifdef TRACY_ENABLE
  ZoneScopedN("getlenedg_geosz0");
  #endif
  constexpr int nnmet = (gdim*(gdim+1))/2;
  double dum[nnmet],tang[gdim];
  double bar1[2];//,bary[tdimn+1];
 
  if constexpr (ideg > 1){
    for(int ii = 0; ii < 2; ii++){
      bar1[0] = 1.0 - (double) ii;
      bar1[1] = (double) ii;
      eval1<gdim,ideg>(msh.coord, edg2pol, msh.getBasis(), 
                       DifVar::Bary, DifVar::None, 
                       bar1, dum, tang, NULL);

   
      if(msh.met.getSpace() == MetSpace::Log){
        sz[ii] = getlenedg_log<gdim>(tang,msh.met[edg2pol[ii]],100,1.0e-6);
      }else{
        sz[ii] = getlenedg<gdim>(tang,msh.met[edg2pol[ii]]);
      }

    }
  }else{// ideg == 1
    //GETVDEPTH(msh.param);
    for(int ii = 0; ii < gdim; ii++) tang[ii] = msh.coord(edg2pol[1],ii)
                                              - msh.coord(edg2pol[0],ii);
    if(msh.met.getSpace() == MetSpace::Log){
      sz[0] = getlenedg_log<gdim>(tang,msh.met[edg2pol[0]],100,1.0e-6);
      sz[1] = getlenedg_log<gdim>(tang,msh.met[edg2pol[1]],100,1.0e-6);
    }else{
      sz[0] = getlenedg<gdim>(tang,msh.met[edg2pol[0]]);
      sz[1] = getlenedg<gdim>(tang,msh.met[edg2pol[1]]);
    }
    //CPRINTF2(" - getlenedg_geosz tang = {}, sz = {} {}\n",
    //         dblAr1(gdim,tang), sz[0], sz[1]);
  }

  if(abs(sz[1]) < Defaults::ltol) return 0;

  double a = sz[0]/sz[1];
  double len = -1;
  if(abs(a-1.0) < 1.0e-12){
    len = sz[1];
  }else{
    len = sz[1] * (a-1.0)/log(a);
  }

  return len;
}



// This version uses surface tangent directions to compute the metric size at 
// the edge extremities. 
template<class MetricFieldType, int gdim, int ideg>
double getlenedg_geosz_plane(const MeshMetric<MetricFieldType> &msh,
                             int ientt, int tdimn, int iedg, 
                             double *sz){
  static_assert(gdim == 3);
  METRIS_ASSERT(tdimn < gdim);
  METRIS_ASSERT(tdimn == 2); // we need a surface ref

  const int nedgl = (tdimn*(tdimn+1))/2;

  int lnoed1[1][2] = {{0, 1}};
  const intAr2 lnoed(nedgl,2,tdimn == 1 ? lnoed1[0] :
                             tdimn == 2 ? lnoed2[0] : lnoed3[0]);

  const intAr2 &ent2poi = tdimn == 1 ? msh.edg2poi : 
                          tdimn == 2 ? msh.fac2poi : msh.tet2poi;

  int edg2pol[getnnod1(ideg)];

  edg2pol[0] = ent2poi(ientt,lnoed[iedg][0]);
  edg2pol[1] = ent2poi(ientt,lnoed[iedg][1]);
  int idx0 = tdimn + 1 + iedg*(ideg-1);
  for(int ii = 0; ii < ideg-1; ii++){
    edg2pol[2+ii] = ent2poi[ientt][idx0+ii];
  }

  double nrmals[2][3];
  int iref = msh.fac2ref[ientt];
  getnorpoiref(msh, edg2pol[0], iref, nrmals[0]);
  getnorpoiref(msh, edg2pol[1], iref, nrmals[1]);

  return getlenedg_geosz_plane<MetricFieldType,gdim,ideg>(msh,&edg2pol[0],nrmals[0],sz);
}

// Provide unit normals to project on plane. 
template<class MetricFieldType, int gdim, int ideg>
double getlenedg_geosz_plane(const MeshMetric<MetricFieldType> &msh,
                             const int *edg2pol, double *nrmals, double *sz){
  constexpr int nnmet = (gdim*(gdim+1))/2;
  double dum[nnmet],tang[gdim];
  double bar1[2];//,bary[tdimn+1];
  static_assert(gdim == 3);
 
  if constexpr (ideg > 1){
    for(int ii = 0; ii < 2; ii++){
      bar1[0] = 1.0 - (double) ii;
      bar1[1] = (double) ii;
      eval1<gdim,ideg>(msh.coord, edg2pol, msh.getBasis(), 
                       DifVar::Bary, DifVar::None, 
                       bar1, dum, tang, NULL);

      double dtprd = getprdl2<gdim>(tang, &nrmals[gdim*ii]);
      for(int jj = 0; jj < gdim; jj++) tang[jj] -= dtprd*nrmals[gdim*ii + jj];
   
      if(msh.met.getSpace() == MetSpace::Log){
        sz[ii] = getlenedg_log<gdim>(tang,msh.met[edg2pol[ii]],100,1.0e-6);
      }else{
        sz[ii] = getlenedg<gdim>(tang,msh.met[edg2pol[ii]]);
      }

    }
  }else{// ideg == 1
    for(int ii = 0; ii < gdim; ii++) tang[ii] = msh.coord(edg2pol[1],ii)
                                              - msh.coord(edg2pol[0],ii);

    double dtpr1 = getprdl2<gdim>(tang, &nrmals[gdim*0]);
    double dtpr2 = getprdl2<gdim>(tang, &nrmals[gdim*1]);
    double tan1[3], tan2[3];
    for(int jj = 0; jj < gdim; jj++) tan1[jj] = tang[jj] - dtpr1*nrmals[gdim*0 + jj];
    for(int jj = 0; jj < gdim; jj++) tan2[jj] = tang[jj] - dtpr2*nrmals[gdim*1 + jj];

    //printf("Debug dtpr1 = {} dtpr2 = {}\n",dtpr1, dtpr2);
    //printf("Debug tang = {} {} {}\n",tang[0],tang[1],tang[2]);
    //printf("Debug nrmal1 = {} {} {}\n",nrmals[gdim*0+0],nrmals[gdim*0+1],nrmals[gdim*0+2]);
    //printf("Debug nrmal2 = {} {} {}\n",nrmals[gdim*1+0],nrmals[gdim*1+1],nrmals[gdim*1+2]);

    if(msh.met.getSpace() == MetSpace::Log){
      sz[0] = getlenedg_log<gdim>(tan1,msh.met[edg2pol[0]],100,1.0e-6);
      sz[1] = getlenedg_log<gdim>(tan2,msh.met[edg2pol[1]],100,1.0e-6);
    }else{
      sz[0] = getlenedg<gdim>(tan1,msh.met[edg2pol[0]]);
      sz[1] = getlenedg<gdim>(tan2,msh.met[edg2pol[1]]);
    }
  }

  double a = sz[0]/sz[1];
  double len = -1;
  if(abs(a-1.0) < 1.0e-12){
    len = sz[1];
  }else{
    len = sz[1] * (a-1.0)/log(a);
  }

  return len;
}
 
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double getlenedg_geosz<MetricFieldAnalytical, 2, n >(\
        const MeshMetric<MetricFieldAnalytical> &msh,int ientt, int tdimn, int iedg);\
template double getlenedg_geosz<MetricFieldFE        , 2, n >(\
        const MeshMetric<MetricFieldFE        > &msh,int ientt, int tdimn, int iedg);\
template double getlenedg_geosz<MetricFieldAnalytical, 3, n >(\
        const MeshMetric<MetricFieldAnalytical> &msh,int ientt, int tdimn, int iedg);\
template double getlenedg_geosz<MetricFieldFE        , 3, n >(\
        const MeshMetric<MetricFieldFE        > &msh,int ientt, int tdimn, int iedg);\
template double getlenedg_geosz<MetricFieldAnalytical, 2, n >(\
        const MeshMetric<MetricFieldAnalytical> &msh,int ientt, int tdimn, int iedg, double* sz);\
template double getlenedg_geosz<MetricFieldFE        , 2, n >(\
        const MeshMetric<MetricFieldFE        > &msh,int ientt, int tdimn, int iedg, double* sz);\
template double getlenedg_geosz<MetricFieldAnalytical, 3, n >(\
        const MeshMetric<MetricFieldAnalytical> &msh,int ientt, int tdimn, int iedg, double* sz);\
template double getlenedg_geosz<MetricFieldFE        , 3, n >(\
        const MeshMetric<MetricFieldFE        > &msh,int ientt, int tdimn, int iedg, double* sz);\
template double getlenedg_geosz<MetricFieldAnalytical, 2, n >(\
        const MeshMetric<MetricFieldAnalytical> &msh,const int *edg2pol, double *sz);\
template double getlenedg_geosz<MetricFieldFE        , 2, n >(\
        const MeshMetric<MetricFieldFE        > &msh,const int *edg2pol, double *sz);\
template double getlenedg_geosz<MetricFieldAnalytical, 3, n >(\
        const MeshMetric<MetricFieldAnalytical> &msh,const int *edg2pol, double *sz);\
template double getlenedg_geosz<MetricFieldFE        , 3, n >(\
        const MeshMetric<MetricFieldFE        > &msh,const int *edg2pol, double *sz);\
template double getlenedg_geosz_plane<MetricFieldAnalytical, 3, n >(\
        const MeshMetric<MetricFieldAnalytical> &msh,int ientt, int tdimn, int iedg, double* sz);\
template double getlenedg_geosz_plane<MetricFieldFE        , 3, n >(\
        const MeshMetric<MetricFieldFE        > &msh,int ientt, int tdimn, int iedg, double* sz);\
template double getlenedg_geosz_plane<MetricFieldAnalytical, 3, n >(\
        const MeshMetric<MetricFieldAnalytical> &msh,const int *edg2pol, double *nrmals, double *sz);\
template double getlenedg_geosz_plane<MetricFieldFE        , 3, n >(\
        const MeshMetric<MetricFieldFE        > &msh,const int *edg2pol, double *nrmals, double *sz);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



///*
//Compute length by uniform quadrature. 
//Metric is log-Euclidean interpolated along edge
//Improve this with better quadrature schemes in the future.
//*/
//template<int gdim, int ideg>
//double getlenedg_quad(const int* __restrict__ edg2pol, 
//                      const dblAr2& __restrict__ coord, 
//                      const dblAr2& __restrict__ met, 
//                      int nquad){
//  constexpr int nnmet = (gdim*(gdim+1))/2;
//
//  double len = 0.0;
//  double dum[nnmet],tang[gdim],metl[nnmet];
//  double bary[2];
//  double dx = 1.0/(nquad - 1.0);
//  for(int iquad = 0; iquad < nquad; iquad++){
//    bary[0] = iquad/(nquad-1.0);
//    bary[1] = 1.0 - bary[0];
//    eval1_bezier<gdim,ideg>(coord, edg2pol,DifVar::Bary,bary,  dum, tang);
//    eval1_bezier<gdim,ideg>(met  , edg2pol,DifVar::None,bary, metl, dum );
//    len += getlenedg<gdim>(tang,metl)*dx;
//  }
//
//  return len;
//}

// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
//#define BOOST_PP_LOCAL_MACRO(n)\
//template double getlenedg_quad<2, n >(const int* __restrict__ edg2pol, const dblAr2& __restrict__ coord, const dblAr2& __restrict__ met, int nquad);\
//template double getlenedg_quad<3, n >(const int* __restrict__ edg2pol, const dblAr2& __restrict__ coord, const dblAr2& __restrict__ met, int nquad);
//#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
//#include BOOST_PP_LOCAL_ITERATE()


template<class MetricFieldType, int gdim, int ideg>
double getlenedg_quad(MeshMetric<MetricFieldType> &msh,
                      int ientt, int tdimn, int iedg, int nquad){

  constexpr int nnmet = (gdim*(gdim+1))/2;
  const int nedgl = (tdimn*(tdimn+1))/2;
  METRIS_ASSERT(0 <= iedg && iedg < nedgl); // as many metric comps as edges

  double len = 0.0;
  double dum[nnmet],tang[gdim],metl[nnmet];
  //double bary[tdimn+1];
  double dx = 1.0/nquad;

  //for(int ii = 0; ii < tdimn + 1 ; ii ++) bary[ii] = 0;

  int lnoed1[2] = {0, 1};
  const int* lnoed = tdimn == 1 ? lnoed1 :
                     tdimn == 2 ? lnoed2[iedg] : lnoed3[iedg];
  intAr2 &ent2poi = msh.ent2poi(tdimn);

  int edg2pol[getnnod1(ideg)];

  // In this case, metric does not need reinterpolating. 
  bool atnodes = nquad == msh.nnode(1);
  METRIS_ASSERT(!(atnodes && (  msh.met.getBasis() != FEBasis::Lagrange
                             || msh.met.getSpace() != MetSpace::Exp)));

  if constexpr (ideg == 1){
    int ipoi1 = ent2poi(ientt,lnoed[0]);
    int ipoi2 = ent2poi(ientt,lnoed[1]);
    for(int ii = 0; ii < gdim; ii++){
      tang[ii] = msh.coord(ipoi1, ii) - msh.coord(ipoi2, ii);
    }
  }

  if(ideg > 1 || !atnodes){
    edg2pol[0] = ent2poi(ientt,lnoed[0]);
    edg2pol[1] = ent2poi(ientt,lnoed[1]);
    int idx0 = tdimn + 1 + iedg*(ideg-1);
    for(int ii = 0; ii < ideg-1; ii++){
      edg2pol[2+ii] = ent2poi[ientt][idx0+ii];
    }
  }

  double bar1[2];
  for(int iquad = 0; iquad < nquad; iquad++){

    if constexpr(ideg > 1){
      if(nquad > 1){
        bar1[0] = iquad/(nquad-1.0);
        bar1[1] = 1.0 - bar1[0];
      }else{
        bar1[0] = 0.5;
        bar1[1] = 0.5;
      }
      eval1<gdim,ideg>(msh.coord, edg2pol, msh.getBasis(), DifVar::Bary, DifVar::None, 
                       bar1,  dum, tang, NULL);
    }

    if(!atnodes){
      if constexpr(ideg == 1){ // case where not previously init
        if(nquad > 1){
          bar1[0] = iquad/(nquad-1.0);
          bar1[1] = 1.0 - bar1[0];
        }else{
          bar1[0] = 0.5;
          bar1[1] = 0.5;
        }
      }
      msh.met.getMetBary(AsDeg::Pk, DifVar::None, 
                         MetSpace::Exp, edg2pol, 1, bar1, metl, NULL);
      len += getlenedg<gdim>(tang,metl)*dx;
    }else{
      int ipoin = ent2poi(ientt, lnoed[iquad]);
      len += getlenedg<gdim>(tang,msh.met[ipoin])*dx;
    }

    METRIS_ASSERT_MSG(!std::isnan(len),"nan len increment in getlenedg_quad\n"
      "ientt {} tdim {} iedg {} nquad {} gdim {} ideg {} \n"
      "dx = {}, tang = {}\nmetl = {}",
      ientt, tdimn, iedg, nquad, gdim, ideg, dx, dblAr1(gdim,tang), dblAr1(nnmet,metl));
  }

  return len;
}

// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double getlenedg_quad<MetricFieldAnalytical, 2, n >(MeshMetric<MetricFieldAnalytical> &msh,int ientt, int tdimn, int iedg, int nquad);\
template double getlenedg_quad<MetricFieldFE        , 2, n >(MeshMetric<MetricFieldFE        > &msh,int ientt, int tdimn, int iedg, int nquad);\
template double getlenedg_quad<MetricFieldAnalytical, 3, n >(MeshMetric<MetricFieldAnalytical> &msh,int ientt, int tdimn, int iedg, int nquad);\
template double getlenedg_quad<MetricFieldFE         , 3, n >(MeshMetric<MetricFieldFE         > &msh,int ientt, int tdimn, int iedg, int nquad);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()




}// end namesapce