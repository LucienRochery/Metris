//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __CODEGEN_CCOEF_D__
#define __CODEGEN_CCOEF_D__

#include "types.hxx"
namespace Metris{

template<int ideg> void d_ccoef_genbez2(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,
                                        int ielem, int icoor,
                                        dblAr2& __restrict__ d_ccoef);
template<int ideg> void d_pt_ccoef_genbez2(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,
                                        int ielem, int inode,
                                        dblAr2& __restrict__ d_ccoef);
template<int ideg>
void d_ccoef_genbez3(const intAr2 & __restrict__ tet2poi,
                     const dblAr2& __restrict__ coord,
                     int ielem, int icoor,
                     dblAr2& __restrict__ d_ccoef);

}//End namespace
#endif
