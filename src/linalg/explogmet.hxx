//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_EXPLOGMET__
#define __METRIS_EXPLOGMET__


#include "../linalg/eigen.hxx"
#include "../metris_constants.hxx"
#include <cmath>


namespace Metris{


// -----------------------------------------------------------------------------
// Write log(met) in lmet
template<int ndim, typename T>
void getlogmet_cpy(const T* __restrict__ met, T* __restrict__ lmet);

// -----------------------------------------------------------------------------
// Write log(met) as sum to lmet already initialized
template<int ndim, typename T>
void getlogmet_sum(const T* __restrict__ met, T* __restrict__ lmet);

// -----------------------------------------------------------------------------
// Replace met with lmet
template<int ndim, typename T>
void getlogmet_inp(T *met);


// -----------------------------------------------------------------------------
// Three derivatives
// met[6], dmet[3][6]
// expm[6],dexp[3][6]
// This templated version is for benchmarking purposes. 
//template<int idif>
//int getexpmet_cpy_d0(const double* met ,const double* dmet, double tol, 
//                    double*  __restrict__ expm,double*  __restrict__ dexp, int iscal = 1);

// Using Surreals. About 10% slower but consistent with direct approach. 
void getexpmet_cpy_dS(const double* met ,const double* dmet,  
                      double*  __restrict__ expm,double*  __restrict__ dexp, 
                      double tol = 1.0e-12, int iscal = 1);

// -----------------------------------------------------------------------------
template <int n>
void getexpmet_cpy(const double* met ,double*  __restrict__ expm, double tol = 1.0e-12, int iscal = 1);
void getexpmet_cpy_d(const double* met ,const double* dmet,  
                     double*  __restrict__ expm,double*  __restrict__ dexp, 
                     double tol = 1.0e-12, int iscal = 1);

// -----------------------------------------------------------------------------
// lwork as required by LAPACK_dsyev
inline void getexpmet_cpy_LAPACK(const double lmet[], double met[]){
	double eigval[3], eigvec[9],rwork[10];
	geteigsym<3,double>(lmet,10,rwork,eigval,eigvec);
	eigval[0] = exp(eigval[0]);
	eigval[1] = exp(eigval[1]);
	eigval[2] = exp(eigval[2]);
	eig2met<3,double>(eigval,eigvec,met);
}

template <int ndim, typename T>
void getexpmet_inp(T* met);



template <int gdim, typename T>
void getspacmet_inp(T* met, MetSpace tarspac){
	if(tarspac == MetSpace::Log) getlogmet_inp<gdim,T>(met);
	else                         getexpmet_inp<gdim,T>(met);
}


} // End namespace

#endif
