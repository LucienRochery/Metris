//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_EXPLOGMET__
#define __METRIS_EXPLOGMET__


#include "../linalg/eigen.hxx"
#include "../metris_constants.hxx"
#include <cmath>


namespace Metris{


// -----------------------------------------------------------------------------
template<int ndim, typename T>
void getlogmet_inp(T *met);

template <int ndim, typename T>
void getexpmet_inp(T* met);

template <int gdim, typename T>
void getspacmet_inp(T* met, MetSpace tarspac){
  if(tarspac == MetSpace::Log) getlogmet_inp<gdim,T>(met);
  else                         getexpmet_inp<gdim,T>(met);
}


// -----------------------------------------------------------------------------
// using eigendecomposition
template <int ndim, typename T>
void getexpmet_dsyevq(T* met);
// inp can == out. using Eigen's .exp()
template <int ndim, typename T>
void getexpmet_Eigen(T* inp, T* out);


template <int n>
void getexpmet_cpy(const double* met ,double*  __restrict__ expm);



} // End namespace

#endif
