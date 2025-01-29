//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "low_lenedg.hxx"

#include "Mesh/MeshMetric.hxx"
#include "linalg/symidx.hxx"
#include "linalg/matprods.hxx"
#include "low_geo.hxx"


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

  //printf("## DEBUG log: init len %f metl = ",len);
  //dblAr1(nnmet,metl).print();


  int niter = 2;
  double err;
  do{
    symXsymsub_fac<gdim>(metl,metacc[iwhich],1.0 / niter,metacc[1-iwhich]);
    // Accumulates into metacc[1 - iwhich]
    iwhich = 1 - iwhich;

    err = getlenedgsq<gdim>(dx,metacc[iwhich]);
    len += err;


    //printf("## DEBUG log: niter %d/%d len %f err %15.7e met = ",
    //      niter,miter,len,err);
    //dblAr1(nnmet,metacc[iwhich]).print();


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
double getlenedg_geosz(MeshMetric<MetricFieldType> &msh,
                       int ientt, int tdimn, int iedg,ExpTyp iexptyp){
  double sz[2];
  double len = getlenedg_geosz<MetricFieldType,gdim,ideg>(msh,ientt,tdimn,iedg,sz,iexptyp);
  return len;
}


// Same but also return the sizes (e.g. for insertion)
template<class MetricFieldType, int gdim, int ideg>
double getlenedg_geosz(MeshMetric<MetricFieldType> &msh,
                       int ientt, int tdimn, int iedg, 
                       double *sz,ExpTyp iexptyp){

  const int nedgl = (tdimn*(tdimn+1))/2;

  int lnoed1[1][2] = {{0, 1}};
  const intAr2 lnoed(nedgl,2,tdimn == 1 ? lnoed1[0] :
                             tdimn == 2 ? lnoed2[0] : lnoed3[0]);

  intAr2 &ent2poi = tdimn == 1 ? msh.edg2poi : 
                    tdimn == 2 ? msh.fac2poi : msh.tet2poi;

  int edg2pol[edgnpps[ideg]];

  edg2pol[0] = ent2poi(ientt,lnoed[iedg][0]);
  edg2pol[1] = ent2poi(ientt,lnoed[iedg][1]);
  int idx0 = tdimn + 1 + iedg*(ideg-1);
  for(int ii = 0; ii < ideg-1; ii++){
    edg2pol[2+ii] = ent2poi[ientt][idx0+ii];
  }


  return getlenedg_geosz<MetricFieldType,gdim,ideg>(msh,&edg2pol[0],sz,iexptyp);
}



// Same but also return the sizes (e.g. for insertion)
template<class MetricFieldType, int gdim, int ideg>
double getlenedg_geosz(MeshMetric<MetricFieldType> &msh,int *edg2pol, double *sz,
                       ExpTyp iexptyp){
  constexpr int nnmet = (gdim*(gdim+1))/2;
  double dum[nnmet],tang[gdim];
  double bar1[2];//,bary[tdimn+1];
 
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
        MeshMetric<MetricFieldAnalytical> &msh,int ientt, int tdimn, int iedg, ExpTyp iexptyp);\
template double getlenedg_geosz<MetricFieldFE        , 2, n >(\
        MeshMetric<MetricFieldFE        > &msh,int ientt, int tdimn, int iedg, ExpTyp iexptyp);\
template double getlenedg_geosz<MetricFieldAnalytical, 3, n >(\
        MeshMetric<MetricFieldAnalytical> &msh,int ientt, int tdimn, int iedg, ExpTyp iexptyp);\
template double getlenedg_geosz<MetricFieldFE        , 3, n >(\
        MeshMetric<MetricFieldFE        > &msh,int ientt, int tdimn, int iedg, ExpTyp iexptyp);\
template double getlenedg_geosz<MetricFieldAnalytical, 2, n >(\
        MeshMetric<MetricFieldAnalytical> &msh,int ientt, int tdimn, int iedg, double* sz,ExpTyp iexptyp);\
template double getlenedg_geosz<MetricFieldFE        , 2, n >(\
        MeshMetric<MetricFieldFE        > &msh,int ientt, int tdimn, int iedg, double* sz,ExpTyp iexptyp);\
template double getlenedg_geosz<MetricFieldAnalytical, 3, n >(\
        MeshMetric<MetricFieldAnalytical> &msh,int ientt, int tdimn, int iedg, double* sz,ExpTyp iexptyp);\
template double getlenedg_geosz<MetricFieldFE        , 3, n >(\
        MeshMetric<MetricFieldFE        > &msh,int ientt, int tdimn, int iedg, double* sz,ExpTyp iexptyp);\
template double getlenedg_geosz<MetricFieldAnalytical, 2, n >(\
        MeshMetric<MetricFieldAnalytical> &msh,int *edg2pol, double *sz,ExpTyp iexptyp);\
template double getlenedg_geosz<MetricFieldFE        , 2, n >(\
        MeshMetric<MetricFieldFE        > &msh,int *edg2pol, double *sz,ExpTyp iexptyp);\
template double getlenedg_geosz<MetricFieldAnalytical, 3, n >(\
        MeshMetric<MetricFieldAnalytical> &msh,int *edg2pol, double *sz,ExpTyp iexptyp);\
template double getlenedg_geosz<MetricFieldFE        , 3, n >(\
        MeshMetric<MetricFieldFE        > &msh,int *edg2pol, double *sz,ExpTyp iexptyp);
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

  int lnoed1[1][2] = {{0, 1}};
  const intAr2 lnoed(nedgl,2,tdimn == 1 ? lnoed1[0] :
                             tdimn == 2 ? lnoed2[0] : lnoed3[0]);
//  const int **lnoed = tdimn == 1 ? &lnoed1[0] :
//                      tdimn == 2 ? &lnoed2[0] : &lnoed3[0];
  intAr2 &ent2poi = tdimn == 1 ? msh.edg2poi : 
                    tdimn == 2 ? msh.fac2poi : msh.tet2poi;

  int edg2pol[edgnpps[ideg]];
  edg2pol[0] = ent2poi(ientt,lnoed[iedg][0]);
  edg2pol[1] = ent2poi(ientt,lnoed[iedg][1]);
  int idx0 = tdimn + 1 + iedg*(ideg-1);
  for(int ii = 0; ii < ideg-1; ii++){
    edg2pol[2+ii] = ent2poi[ientt][idx0+ii];
  }


  double bar1[2];
  for(int iquad = 0; iquad < nquad; iquad++){

    bar1[0] = iquad/(nquad-1.0);
    bar1[1] = 1.0 - bar1[0];
    eval1<gdim,ideg>(msh.coord, edg2pol, msh.getBasis(), DifVar::Bary, DifVar::None, 
                     bar1,  dum, tang, NULL);

    //bary[lnoed[iedg][0]] = bar1[0];
    //bary[lnoed[iedg][1]] = bar1[1];
    msh.met.getMetBary(AsDeg::Pk, DifVar::None, 
                       MetSpace::Exp, edg2pol, 1, bar1, metl, NULL);
//    eval1_bezier<gdim,ideg>(met  , edg2pol,DifVar::None,bary, metl, dum );

    len += getlenedg<gdim>(tang,metl)*dx;
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