//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LINALG_INVMAT__
#define __METRIS_LINALG_INVMAT__


namespace Metris{
/*
	LAPACK WRAPPERS & AUX
*/
// -----------------------------------------------------------------------------
// met must be positive definite. Otherwise use inv3sym. 
void invspd(int n, double met[]);
// -----------------------------------------------------------------------------
// met must be positive definite. Otherwise use inv3sym. 
void inv2spd(double met[]);

// 
template <int n>
void invsym(double* met);

template<typename T>
void inv3sym(T *met, T *inv);

// Matrix stored line first in C fashion
void invmat(int n, double mat[]);

// Ad-hoc one (dims 2 and 3)
template<int n>
void invmat(double mat[]);


} // End namespace

#endif
