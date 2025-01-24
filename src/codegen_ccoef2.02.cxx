//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "codegen_ccoef.hxx"

#include "types.hxx"

namespace Metris{

template<int ideg>
void ccoef_genbez2(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef){}

double det2_vdif(const double* x1,const double* x2
                ,const double* y1,const double* y2);

template<> void ccoef_genbez2<2>(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord, int ielem, double* __restrict__ ccoef){

  ccoef[  0] =   8*det2_vdif(coord[fac2poi[ielem][   5]],coord[fac2poi[ielem][   0]]
                            ,coord[fac2poi[ielem][   4]],coord[fac2poi[ielem][   0]])/  2;

  ccoef[  1] =   8*det2_vdif(coord[fac2poi[ielem][   1]],coord[fac2poi[ielem][   5]]
                            ,coord[fac2poi[ielem][   3]],coord[fac2poi[ielem][   5]])/  2;

  ccoef[  2] =   8*det2_vdif(coord[fac2poi[ielem][   3]],coord[fac2poi[ielem][   4]]
                            ,coord[fac2poi[ielem][   2]],coord[fac2poi[ielem][   4]])/  2;

  ccoef[  3] =   4*det2_vdif(coord[fac2poi[ielem][   3]],coord[fac2poi[ielem][   4]]
                            ,coord[fac2poi[ielem][   3]],coord[fac2poi[ielem][   5]])/  2
             +   4*det2_vdif(coord[fac2poi[ielem][   1]],coord[fac2poi[ielem][   5]]
                            ,coord[fac2poi[ielem][   2]],coord[fac2poi[ielem][   4]])/  2;

  ccoef[  4] =   4*det2_vdif(coord[fac2poi[ielem][   3]],coord[fac2poi[ielem][   4]]
                            ,coord[fac2poi[ielem][   4]],coord[fac2poi[ielem][   0]])/  2
             +   4*det2_vdif(coord[fac2poi[ielem][   5]],coord[fac2poi[ielem][   0]]
                            ,coord[fac2poi[ielem][   2]],coord[fac2poi[ielem][   4]])/  2;

  ccoef[  5] =   4*det2_vdif(coord[fac2poi[ielem][   1]],coord[fac2poi[ielem][   5]]
                            ,coord[fac2poi[ielem][   4]],coord[fac2poi[ielem][   0]])/  2
             +   4*det2_vdif(coord[fac2poi[ielem][   5]],coord[fac2poi[ielem][   0]]
                            ,coord[fac2poi[ielem][   3]],coord[fac2poi[ielem][   5]])/  2;
}


} // End namespace
