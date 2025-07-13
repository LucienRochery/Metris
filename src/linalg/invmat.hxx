//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LINALG_INVMAT__
#define __METRIS_LINALG_INVMAT__


namespace Metris{
// -----------------------------------------------------------------------------
// met must be positive definite. Otherwise use inv3sym. 
template<int ndim, typename T>
int invspd(T* met);

template<int n>
int invmat(double mat[]);



// ----------------------------------------------------------------------- AUX
// inp and out can be the same
template<int ndim, typename T>
int invspd_Eigen(T *inp, T *out);
#ifdef METRIS_USE_LAPACK
  int invspd_LAPACK(int n, double met[]);
#endif


// Ad-hoc one (dims 2 and 3)
template<int n>
int invmat_naive(double mat[]);
#ifdef METRIS_USE_LAPACK
  int invmat_LAPACK(int n, double *mat);
#endif
template<int ndim, typename T>
int invmat_EigenLUPP(T* mat, T* inv);
template<int ndim, typename T>
int invmat_EigenLUFP(T* mat, T* inv);

} // End namespace

#endif
