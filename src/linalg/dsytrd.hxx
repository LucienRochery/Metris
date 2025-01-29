//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



// ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
// ----------------------------------------------------------------------------
#ifndef __DSYTRD3_H__
#define __DSYTRD3_H__

namespace Metris{

// d size 3, e size 2
// A in sym2idx[] format cf sym2idx
// 1 2 4 
//   3 5 
//     6


template<typename T>
void dsytrd3(const T* __restrict__ A, T* __restrict__  Q, T* __restrict__  d, T* __restrict__  e);
template<typename T>
void dsytrd2(const T* __restrict__ A, T* __restrict__  Q, T* __restrict__  d, T* __restrict__  e);

template<typename T, int ndim>
inline void dsytrd(const T* __restrict__ A, T* __restrict__  Q, T* __restrict__  d, T* __restrict__  e){
	static_assert(ndim == 2 || ndim == 3);
	if constexpr(ndim == 2) dsytrd2(A, Q, d, e);
	else                    dsytrd3(A, Q, d, e);
}

} // End namespace
#endif
