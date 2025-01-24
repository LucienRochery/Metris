//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "codegen_ccoef_d.hxx"

#include "types.hxx"

namespace Metris{

template<int ideg> void d_ccoef_genbez2(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,
                                        int ielem, double*__restrict__ ccoef, 
                                        int icoor, dblAr2& __restrict__ d_ccoef){}
double det2_vdif(const double* x1,const double* x2
                ,const double* y1,const double* y2);

double vdiff_perp_x(const double* a,const double* b);

double vdiff_perp_y(const double* a,const double* b);

template<> void d_ccoef_genbez2<2>(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,
                                   int ielem, double*__restrict__ ccoef, 
                                   int icoor, dblAr2& __restrict__ d_ccoef){
  METRIS_ASSERT(icoor == 0 || icoor == 1);

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

  if(icoor == 0){

    d_ccoef[0][0] =   4*(vdiff_perp_x(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   4]]) + vdiff_perp_x(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   0]]));

    d_ccoef[0][4] =   4*(vdiff_perp_x(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   5]]));

    d_ccoef[0][5] =   4*(vdiff_perp_x(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   0]]));

    d_ccoef[1][1] =   4*(vdiff_perp_x(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   5]]));

    d_ccoef[1][3] =   4*(vdiff_perp_x(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   1]]));

    d_ccoef[1][5] =   4*(vdiff_perp_x(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   3]]) + vdiff_perp_x(coord[fac2poi[ielem][   1]], coord[fac2poi[ielem][   5]]));

    d_ccoef[2][2] =   4*(vdiff_perp_x(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   3]]));

    d_ccoef[2][3] =   4*(vdiff_perp_x(coord[fac2poi[ielem][   2]], coord[fac2poi[ielem][   4]]));

    d_ccoef[2][4] =   4*(vdiff_perp_x(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   2]]) + vdiff_perp_x(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   4]]));

    d_ccoef[3][1] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   2]], coord[fac2poi[ielem][   4]]));

    d_ccoef[3][2] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   1]]));

    d_ccoef[3][3] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   5]]) + vdiff_perp_x(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   3]]));

    d_ccoef[3][4] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   3]]) + vdiff_perp_x(coord[fac2poi[ielem][   1]], coord[fac2poi[ielem][   5]]));

    d_ccoef[3][5] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   4]]) + vdiff_perp_x(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   2]]));

    d_ccoef[4][0] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   4]]) + vdiff_perp_x(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   2]]));

    d_ccoef[4][2] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   5]]));

    d_ccoef[4][3] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   0]]));

    d_ccoef[4][4] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   4]]) + vdiff_perp_x(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   3]]) + vdiff_perp_x(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   0]]));

    d_ccoef[4][5] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   2]], coord[fac2poi[ielem][   4]]));

    d_ccoef[5][0] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   1]], coord[fac2poi[ielem][   5]]) + vdiff_perp_x(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   3]]));

    d_ccoef[5][1] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   0]]));

    d_ccoef[5][3] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   5]]));

    d_ccoef[5][4] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   1]]));

    d_ccoef[5][5] =   2*(vdiff_perp_x(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   4]]) + vdiff_perp_x(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   5]]) + vdiff_perp_x(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   0]]));

  }else{

    d_ccoef[0][0] =   4*(vdiff_perp_y(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   4]]) + vdiff_perp_y(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   0]]));

    d_ccoef[0][4] =   4*(vdiff_perp_y(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   5]]));

    d_ccoef[0][5] =   4*(vdiff_perp_y(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   0]]));

    d_ccoef[1][1] =   4*(vdiff_perp_y(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   5]]));

    d_ccoef[1][3] =   4*(vdiff_perp_y(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   1]]));

    d_ccoef[1][5] =   4*(vdiff_perp_y(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   3]]) + vdiff_perp_y(coord[fac2poi[ielem][   1]], coord[fac2poi[ielem][   5]]));

    d_ccoef[2][2] =   4*(vdiff_perp_y(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   3]]));

    d_ccoef[2][3] =   4*(vdiff_perp_y(coord[fac2poi[ielem][   2]], coord[fac2poi[ielem][   4]]));

    d_ccoef[2][4] =   4*(vdiff_perp_y(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   2]]) + vdiff_perp_y(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   4]]));

    d_ccoef[3][1] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   2]], coord[fac2poi[ielem][   4]]));

    d_ccoef[3][2] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   1]]));

    d_ccoef[3][3] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   5]]) + vdiff_perp_y(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   3]]));

    d_ccoef[3][4] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   3]]) + vdiff_perp_y(coord[fac2poi[ielem][   1]], coord[fac2poi[ielem][   5]]));

    d_ccoef[3][5] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   4]]) + vdiff_perp_y(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   2]]));

    d_ccoef[4][0] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   4]]) + vdiff_perp_y(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   2]]));

    d_ccoef[4][2] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   5]]));

    d_ccoef[4][3] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   0]]));

    d_ccoef[4][4] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   4]]) + vdiff_perp_y(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   3]]) + vdiff_perp_y(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   0]]));

    d_ccoef[4][5] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   2]], coord[fac2poi[ielem][   4]]));

    d_ccoef[5][0] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   1]], coord[fac2poi[ielem][   5]]) + vdiff_perp_y(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   3]]));

    d_ccoef[5][1] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   0]]));

    d_ccoef[5][3] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   5]]));

    d_ccoef[5][4] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   1]]));

    d_ccoef[5][5] =   2*(vdiff_perp_y(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   4]]) + vdiff_perp_y(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   5]]) + vdiff_perp_y(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   0]]));
  }
}


} // End namespace
