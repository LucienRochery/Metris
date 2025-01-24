#ifndef __LAG2BEZ__
#define __LAG2BEZ__

//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "types.hxx"

namespace Metris{

template<int ideg, int szfld>
void lag2bez1(const int* __restrict__ lfld,
 								const dblAr2& __restrict__ rfld0,
 								dblAr2& __restrict__ rfld1);
template<int ideg, int szfld>
void lag2bez2(const int* __restrict__ lfld,
 								const dblAr2& __restrict__ rfld0,
 								dblAr2& __restrict__ rfld1);
template<int ideg, int szfld>
void lag2bez3(const int* __restrict__ lfld,
 								const dblAr2& __restrict__ rfld0,
 								dblAr2& __restrict__ rfld1);

} // End namespace


#endif
