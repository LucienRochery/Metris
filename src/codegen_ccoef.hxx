//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __CODEGEN_CCOEF__
#define __CODEGEN_CCOEF__

#include "types.hxx"
namespace Metris{

template<int ideg>
void ccoef_genbez3(const intAr2 & __restrict__ tet2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef);
template<int ideg>
void ccoef_genbez2(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef);
}//End namespace
#endif
