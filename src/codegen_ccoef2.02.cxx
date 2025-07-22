//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "codegen_ccoef.hxx"
#include "linalg/det.hxx"
#include "types.hxx"

namespace Metris{

template<int ideg>
void ccoef_genbez2([[maybe_unused]]const intAr2 & __restrict__ fac2poi,[[maybe_unused]]const dblAr2& __restrict__ coord,[[maybe_unused]]int ielem,[[maybe_unused]]double* __restrict__ ccoef){}


template<> void ccoef_genbez2<2>(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef){

  ccoef[  0] =   4*detvdif2(coord[fac2poi[ielem][   5]],coord[fac2poi[ielem][   0]]
                            ,coord[fac2poi[ielem][   4]],coord[fac2poi[ielem][   0]])/  1;

  ccoef[  1] =   4*detvdif2(coord[fac2poi[ielem][   1]],coord[fac2poi[ielem][   5]]
                            ,coord[fac2poi[ielem][   3]],coord[fac2poi[ielem][   5]])/  1;

  ccoef[  2] =   4*detvdif2(coord[fac2poi[ielem][   3]],coord[fac2poi[ielem][   4]]
                            ,coord[fac2poi[ielem][   2]],coord[fac2poi[ielem][   4]])/  1;

  ccoef[  3] =   2*detvdif2(coord[fac2poi[ielem][   3]],coord[fac2poi[ielem][   4]]
                            ,coord[fac2poi[ielem][   3]],coord[fac2poi[ielem][   5]])/  1
             +   2*detvdif2(coord[fac2poi[ielem][   1]],coord[fac2poi[ielem][   5]]
                            ,coord[fac2poi[ielem][   2]],coord[fac2poi[ielem][   4]])/  1;

  ccoef[  4] =   2*detvdif2(coord[fac2poi[ielem][   3]],coord[fac2poi[ielem][   4]]
                            ,coord[fac2poi[ielem][   4]],coord[fac2poi[ielem][   0]])/  1
             +   2*detvdif2(coord[fac2poi[ielem][   5]],coord[fac2poi[ielem][   0]]
                            ,coord[fac2poi[ielem][   2]],coord[fac2poi[ielem][   4]])/  1;

  ccoef[  5] =   2*detvdif2(coord[fac2poi[ielem][   1]],coord[fac2poi[ielem][   5]]
                            ,coord[fac2poi[ielem][   4]],coord[fac2poi[ielem][   0]])/  1
             +   2*detvdif2(coord[fac2poi[ielem][   5]],coord[fac2poi[ielem][   0]]
                            ,coord[fac2poi[ielem][   3]],coord[fac2poi[ielem][   5]])/  1;
}


} // End namespace
