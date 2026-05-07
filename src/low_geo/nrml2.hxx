//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_LOW_GEO_NRML2__
#define __METRIS_LOW_GEO_NRML2__

#include "../libs/SANS/Surreal/SurrealS.h"
#include "../types_scalar.hxx"
#include "../Arrays/aux_msharrays.hxx"

namespace Metris{

template<int n, typename ftype = double>
inline ftype geterrl2(const ftype x[], const ftype y[]){
  static_assert(n > 0);
  return (x[0]-y[0])*(x[0]-y[0]) + geterrl2<n-1,ftype>(&x[1],&y[1]);
}

template<> inline double geterrl2<1,double>(const double x[], const double y[]){
  return (x[0]-y[0])*(x[0]-y[0]);
}

#ifdef USE_MULTIPRECISION
  template<> inline float8 geterrl2<1,float8>(const float8 x[], const float8 y[]){
    return (x[0]-y[0])*(x[0]-y[0]);
  }
#endif

template<int n, int m, typename ftype = double>
inline ftype geterrl2(const SANS::SurrealS<m,ftype> x[], const SANS::SurrealS<m,ftype> y[]){
  static_assert(n > 0);
  ftype ret = 0;
  for(int i = 0; i < n ;i++){
    ret += (x[i].value()-y[i].value())*(x[i].value()-y[i].value());
  }
  return ret;
}


inline double geterrl2(int n, const double* x, const double *y){
  double ret = 0;
  for(int ii = 0; ii < n; ii++) ret += (x[ii] - y[ii])*(x[ii] - y[ii]);
  return ret;
}


template<int n, typename ftype = double>
inline ftype geterrl2(const MeshArray1D<ftype,int32_t> &x,
                      const MeshArray1D<ftype,int32_t> &y){
  static_assert(n > 0);
  ftype ret = 0;
  for(int ii = 0; ii < n; ii++){
    ret += (x[ii]-y[ii])*(x[ii]-y[ii]);
  }
  return ret;
}

template<int n, typename ftype = double>
inline ftype geterrl2(const MeshArray1D<ftype,int32_t> &x, const ftype* y){
  static_assert(n > 0);
  ftype ret = 0;
  for(int ii = 0; ii < n; ii++){
    ret += (x[ii]-y[ii])*(x[ii]-y[ii]);
  }
  return ret;
}

template<int n, typename ftype = double>
inline ftype geterrl2(const ftype* x, const MeshArray1D<ftype,int32_t> &y){
  static_assert(n > 0);
  ftype ret = 0;
  for(int ii = 0; ii < n; ii++){
    ret += (x[ii]-y[ii])*(x[ii]-y[ii]);
  }
  return ret;
}

template<int n, typename ftype = double>
inline ftype getprdl2(const ftype* __restrict__ x,const ftype* __restrict__ y){
  static_assert(n > 0);
  return x[0]*y[0] + getprdl2<n-1,ftype>(&x[1],&y[1]);
}
template<> inline double getprdl2<1,double>(const double* __restrict__ x,
                                             const double* __restrict__ y){
  return x[0]*y[0];
}

#ifdef USE_MULTIPRECISION
  template<> inline float8 getprdl2<1,float8>(const float8* __restrict__ x,
                                               const float8* __restrict__ y){
    return x[0]*y[0];
  }
#endif



// (x1-x2) . (y1-y2)
template<int n, typename ftype = double>
inline ftype getprdl2(const ftype* __restrict__ x1,const ftype* __restrict__ x2,
                      const ftype* __restrict__ y1,const ftype* __restrict__ y2){
  static_assert(n > 0);
  return (x1[0]-x2[0])*(y1[0]-y2[0])
         + getprdl2<n-1,ftype>(&x1[1],&x2[1],&y1[1],&y2[1]);
}
template<> inline double getprdl2<1,double>(
                   const double* __restrict__ x1,const double* __restrict__ x2,
                   const double* __restrict__ y1,const double* __restrict__ y2){
  return (x1[0]-x2[0])*(y1[0]-y2[0]);
}

#ifdef USE_MULTIPRECISION
  template<> inline float8 getprdl2<1,float8>(
                   const float8* __restrict__ x1,const float8* __restrict__ x2,
                   const float8* __restrict__ y1,const float8* __restrict__ y2){
    return (x1[0]-x2[0])*(y1[0]-y2[0]);
  }
#endif


template<int n, typename ftype = double>
inline ftype getnrml2(const ftype x[]){
  static_assert(n > 0);
  return getprdl2<n,ftype>(x,x);
}

}// namespace Metris

#endif
