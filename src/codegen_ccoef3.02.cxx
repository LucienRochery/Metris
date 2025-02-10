//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "codegen_ccoef.hxx"
#include "types.hxx"

namespace Metris{

template<int ideg>
void ccoef_genbez3([[maybe_unused]] const intAr2&__restrict__ tet2poi,[[maybe_unused]] const dblAr2& __restrict__ coord,[[maybe_unused]] int ielem,[[maybe_unused]] double*__restrict__ ccoef){}

double det3_vdif(const double* x1,const double* x2
                ,const double* y1,const double* y2
                ,const double* z1,const double* z2);

template<> void ccoef_genbez3<2>(const intAr2 & __restrict__ tet2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef){

  ccoef[  0] =   8*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  1;

  ccoef[  1] =   8*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  1;

  ccoef[  2] =   8*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  1;

  ccoef[  3] =   8*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  1;

  ccoef[  4] =   8*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3;

  ccoef[  5] =   8*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3;

  ccoef[  6] =   8*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3;

  ccoef[  7] =   8*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3;

  ccoef[  8] =   8*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3;

  ccoef[  9] =   8*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3;

  ccoef[ 10] =   8*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3;

  ccoef[ 11] =   8*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3;

  ccoef[ 12] =   8*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3;

  ccoef[ 13] =   8*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3;

  ccoef[ 14] =   8*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3;

  ccoef[ 15] =   8*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3
             +   8*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3;

  ccoef[ 16] =   4*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3;

  ccoef[ 17] =   4*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3;

  ccoef[ 18] =   4*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   7]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   3]],coord[tet2poi[ielem][   7]])/  3;

  ccoef[ 19] =   4*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   7]],coord[tet2poi[ielem][   0]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   1]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   6]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   2]],coord[tet2poi[ielem][   6]]
                            ,coord[tet2poi[ielem][   8]],coord[tet2poi[ielem][   4]])/  3
             +   4*det3_vdif(coord[tet2poi[ielem][   4]],coord[tet2poi[ielem][   0]]
                            ,coord[tet2poi[ielem][   5]],coord[tet2poi[ielem][   4]]
                            ,coord[tet2poi[ielem][   9]],coord[tet2poi[ielem][   6]])/  3;
}


} // End namespace
