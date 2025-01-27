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
int invspd(int n, double met[]);

// 
template <int n>
int invsym(double* met);

template<typename T>
int inv3sym(T *met, T *inv);

// Matrix stored line first in C fashion
int invmat(int n, double mat[]);

// Ad-hoc one (dims 2 and 3)
template<int n>
int invmat(double mat[]);


} // End namespace

#endif
