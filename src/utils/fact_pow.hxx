//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_FACT_POW__
#define __METRIS_FACT_POW__

namespace Metris{

template<int n>
inline double idpow(double x){
  if constexpr(n%2 == 0){
    return idpow<n/2>(x) * idpow<n/2>(x);
  }else{
    return idpow<n/2>(x) * idpow<n/2>(x) * x;
  }
}
template<> inline double idpow<0>(double){return 1.0;}

template<int n, class ftype>
inline ftype irpow(ftype x){
  if constexpr(n%2 == 0){
    return irpow<n/2,ftype>(x) * irpow<n/2,ftype>(x);
  }else{
    return irpow<n/2,ftype>(x) * irpow<n/2,ftype>(x) * x;
  }
}
template<> inline double irpow<0,double>(double){return 1.0;}


#ifdef USE_MULTIPRECISION
  template<> inline float4 irpow<0,float4>(float4){return 1.0;}
  template<> inline float8 irpow<0,float8>(float8){return 1.0;}

  template<int n>
  inline float4 id4pow(float4 x){
    if(n%2 == 0){
      return id4pow<n/2>(x) * id4pow<n/2>(x);
    }else{
      return id4pow<n/2>(x) * id4pow<n/2>(x) * x;
    }
  }
  template<> inline float4 id4pow<0>(float4){return 1.0;}


  template<int n>
  inline float8 id8pow(float8 x){
    if(n%2 == 0){
      return id8pow<n/2>(x) * id8pow<n/2>(x);
    }else{
      return id8pow<n/2>(x) * id8pow<n/2>(x) * x;
    }
  }
  template<> inline float8 id8pow<0>(float8){return 1.0;}
#endif


template<int n>
inline int iipow(int x){
  if(n%2 == 0){
    return iipow<n/2>(x) * iipow<n/2>(x);
  }else{
    return iipow<n/2>(x) * iipow<n/2>(x) * x;
  }
}
template<> inline int iipow<0>(int ){return 1;}


template<int n> 
inline constexpr int ifact(){
  return n*ifact<n-1>();
}
template<> inline constexpr int ifact<0>(){return 1;}

inline int ifact(int n){
  if(n <= 1) return 1;
  return n*ifact(n-1);
}

} // End namespace

#endif