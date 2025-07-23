//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __LOW_GEO__
#define __LOW_GEO__

#include "nrml2.hxx"

#include "../types.hxx"
#include "../Mesh/MeshFwd.hxx"
#include "../metris_constants.hxx"

#include "../../SANS/Surreal/SurrealS.h"


namespace Metris{

struct MetrisParameters;

enum class FEBasis;

#ifndef ALWAYS_INLINE
	// ALWAYS_INLINE is a macro to further encourage the compiler to inline a function
	#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__clang__)
		#define ALWAYS_INLINE inline __attribute__((always_inline))
	#elif defined(_MSC_VER)
		#define ALWAYS_INLINE __forceinline
	#else
		#warning Not forcing inline with this compiler... (Please add this compiler to tools/always_inline.h)
		#define ALWAYS_INLINE inline
	#endif
#endif



// 
// Symmetric matrix indices: 0 1 3
//                             2 4
//                               5

/*
	LAPACK ROUTINES
*/
//extern "C" {
//		extern void dpptrf_(char*,int*,double*,int*);
//		extern void dpptri_(char*,int*,double*,int*);
//		extern void dsyev(_char*,char*,int*,double*,int*,double*,double*,int*,int*);
//}


void vdiff_perp(const double* a, const double* b, double *res);
void vdiff_perp(const double* a, const double* b, int up, int lo, double *res);
void vdiff_perp_sum(const double* a, const double* b, int up, int lo, double *res);



bool isintetP1(const double *p1, const double *p2,
               const double *p3, const double *p4,
               const double *pp, double tol = 1.0e-16);
bool isinfacP1(const double *p1, const double *p2,
               const double *p3, const double *pp, double tol = 1.0e-16);

template<int gdim>
void inventP1(const int*__restrict__ ent2pol, const dblAr2 &coord, const double*__restrict__ coor0, 
              double*__restrict__ bary);

// -----------------------------------------------------------------------------





//// Discrete element quality defined as l^p sum of control polygon qualities
//template<int ideg>
//double geteltqua_disc(dblAr2 &coord, dblAr2 &met, int* __restrict__ tet2pol, int pnorm, int *ierro);



template<int n>
int normalize_vec(double *vec){
  double nrm = getnrml2<n>(vec);
  if(nrm < Constants::vecNrmTol) return 1;
  for(int ii = 0; ii < n; ii++) vec[ii] /= sqrt(nrm);
  return 0;
}



// Characteristic element length for tolerance scaling
// Minimum edge length for now
template<int gdim>
double getepsent(MeshBase &msh, int tdimn, int ientt);

} // End namespace





#endif