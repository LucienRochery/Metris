//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __LOW_GEO__
#define __LOW_GEO__


#include "types.hxx"
#include "Mesh/MeshFwd.hxx"

#include "../SANS/Surreal/SurrealS.h"


namespace Metris{

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

// -----------------------------------------------------------------------------
template <typename T=double>
T det3_vdif(const T* x1,const T* x2,
            const T* y1,const T* y2,
            const T* z1,const T* z2);

template<typename T=double>
T det2_vdif(const T* x1,const T* x2, 
            const T* y1,const T* y2);

template <typename T=double>
T* vdiff(const T* a, const T* b);

template <typename T=double>
T* vproduct(const T* a, const T* b);

template <typename T=double>
T* vdiff_perp(const T* a, const T* b);


double det3_vdif(const double* x1,const double* x2
								,const double* y1,const double* y2
								,const double* z1,const double* z2);

double det2_vdif(const double* x1,const double* x2
								,const double* y1,const double* y2);

double* vdiff(const double* a, const double* b);
double* vproduct(const double* a, const double* b);
double* vdiff_perp(const double* a, const double* b);


bool isintetP1(const double *p1, const double *p2,
               const double *p3, const double *p4,
               const double *pp, double tol = 1.0e-16);
bool isinfacP1(const double *p1, const double *p2,
               const double *p3, const double *pp, double tol = 1.0e-16);

template<int gdim>
void inventP1(const int*__restrict__ ent2pol, const dblAr2 &coord, const double*__restrict__ coor0, 
              double*__restrict__ bary);

// -----------------------------------------------------------------------------
// Length, area or volume.
template <int gdim>
double getmeasentP1(const int *ent2pol,const dblAr2& coord);
// Variant that checks tolerance, returns iflat = true if negative or flat
// nrmal can be NULL if not surface
template <int gdim, int tdim>
double getmeasentP1(const int *ent2pol,const dblAr2& coord, double vtol, 
                    double *nrmal, bool *iflat, int iverb = 0);
template <int idim>
void getheightentP1_aniso(const int *ent2pol,const dblAr2 &coord, double *metl, double *height);



template <int gdim>
void getmeasentP1grad(const int *ent2pol, const dblAr2& coord, int idof, double *grad);

void getnorfacP1(const int *fac2pol, const dblAr2 &coord, double *nrmal);

// Return outgoing normal 
int getnorpoiCAD(const MeshBase &msh, int ipoin, std::map<ego,int> &edgorient, 
                 double *norpoi);

template <int ideg>
void getnorpoi(const MeshBase &msh, int ipoin, const intAr1 &lball, double* norpoi);

template <int ideg>
void getnorpoiref(const MeshBase &msh, int ipoin, int iref, const intAr1 &lball, double* norpoi);


// -----------------------------------------------------------------------------
template<int gdim, int tdim, int ideg>
void getintmetxi(const dblAr2 &coord, const int* __restrict__ tet2pol, FEBasis ibasis,
	               const double* bary,double* __restrict__ met);
//template<int gdim, int tdim, int ideg>
//void getintmetxi(const dblAr2 &coord, const int* __restrict__ tet2pol, FEBasis ibasis,
//	               const double* bary,SANS::SurrealS<3,double>* __restrict__ metS);

//// Discrete element quality defined as l^p sum of control polygon qualities
//template<int ideg>
//double geteltqua_disc(dblAr2 &coord, dblAr2 &met, int* __restrict__ tet2pol, int pnorm, int *ierro);

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


template<int n, typename ftype = double>
inline ftype getnrml2(const ftype x[]){
  static_assert(n > 0);
  return getprdl2<n,ftype>(x,x);
}


template<typename ftype = double>
inline void getvecprod3(const ftype* x,const ftype* y,ftype* z){
	z[0] = x[1]*y[2] - x[2]*y[1];
	z[1] = x[2]*y[0] - x[0]*y[2];
	z[2] = x[0]*y[1] - x[1]*y[0];
}



// Characteristic element length for tolerance scaling
// Minimum edge length for now
template<int gdim>
double getepsent(MeshBase &msh, int tdimn, int ientt);

} // End namespace





#endif