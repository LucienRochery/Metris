//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "codegen_ccoef_d.hxx"

#include "types.hxx"

namespace Metris{

template<int ideg> void d_pt_ccoef_genbez2(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,
                                        int ielem,
                                        int icoor,
                                        int inode,
                                        double*__restrict__ ccoef,
                                        double*__restrict__ d_ccoef){}

double det2_vdif(const double* x1,const double* x2
                ,const double* y1,const double* y2);

double* vdiff_perp(const double* a,const double* b);

template<> void d_pt_ccoef_genbez2<2>(const intAr2 & __restrict__ fac2poi, const dblAr2& __restrict__ coord,
                                        int ielem,
                                        int icoor,
                                        int inode,
                                        double*__restrict__ ccoef,
                                        double*__restrict__ d_ccoef){

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

    for(int i = 0; i < 6; i++) {
      d_ccoef[i] = 0;
    }
    if(inode  ==   0){

      d_ccoef[0] =   8*vdiff_perp(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   4]])[icoor] /   2 +   8*vdiff_perp(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   0]])[icoor] /   2;
      d_ccoef[4] =   8*vdiff_perp(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   5]])[icoor] /   2;
      d_ccoef[5] =   8*vdiff_perp(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   0]])[icoor] /   2;
    }
    else if(inode  ==   1){

      d_ccoef[1] =   8*vdiff_perp(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   5]])[icoor] /   2;
      d_ccoef[3] =   8*vdiff_perp(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   1]])[icoor] /   2;
      d_ccoef[5] =   8*vdiff_perp(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   3]])[icoor] /   2 +   8*vdiff_perp(coord[fac2poi[ielem][   1]], coord[fac2poi[ielem][   5]])[icoor] /   2;
    }
    else if(inode  ==   2){

      d_ccoef[2] =   8*vdiff_perp(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   3]])[icoor] /   2;
      d_ccoef[3] =   8*vdiff_perp(coord[fac2poi[ielem][   2]], coord[fac2poi[ielem][   4]])[icoor] /   2;
      d_ccoef[4] =   8*vdiff_perp(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   2]])[icoor] /   2 +   8*vdiff_perp(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   4]])[icoor] /   2;
    }
    else if(inode  ==   3){

      d_ccoef[1] =   4*vdiff_perp(coord[fac2poi[ielem][   2]], coord[fac2poi[ielem][   4]])[icoor] /   2;
      d_ccoef[2] =   4*vdiff_perp(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   1]])[icoor] /   2;
      d_ccoef[3] =   4*vdiff_perp(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   5]])[icoor] /   2 +   4*vdiff_perp(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   3]])[icoor] /   2;
      d_ccoef[4] =   4*vdiff_perp(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   3]])[icoor] /   2 +   4*vdiff_perp(coord[fac2poi[ielem][   1]], coord[fac2poi[ielem][   5]])[icoor] /   2;
      d_ccoef[5] =   4*vdiff_perp(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   4]])[icoor] /   2 +   4*vdiff_perp(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   2]])[icoor] /   2;
    }
    else if(inode  ==   4){

      d_ccoef[0] =   4*vdiff_perp(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   4]])[icoor] /   2 +   4*vdiff_perp(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   2]])[icoor] /   2;
      d_ccoef[2] =   4*vdiff_perp(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   5]])[icoor] /   2;
      d_ccoef[3] =   4*vdiff_perp(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   0]])[icoor] /   2;
      d_ccoef[4] =   4*vdiff_perp(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   4]])[icoor] /   2 +   4*vdiff_perp(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   3]])[icoor] /   2 +   4*vdiff_perp(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   0]])[icoor] /   2;
      d_ccoef[5] =   4*vdiff_perp(coord[fac2poi[ielem][   2]], coord[fac2poi[ielem][   4]])[icoor] /   2;
    }
    else if(inode  ==   5){

      d_ccoef[0] =   4*vdiff_perp(coord[fac2poi[ielem][   1]], coord[fac2poi[ielem][   5]])[icoor] /   2 +   4*vdiff_perp(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   3]])[icoor] /   2;
      d_ccoef[1] =   4*vdiff_perp(coord[fac2poi[ielem][   4]], coord[fac2poi[ielem][   0]])[icoor] /   2;
      d_ccoef[3] =   4*vdiff_perp(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   5]])[icoor] /   2;
      d_ccoef[4] =   4*vdiff_perp(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   1]])[icoor] /   2;
      d_ccoef[5] =   4*vdiff_perp(coord[fac2poi[ielem][   0]], coord[fac2poi[ielem][   4]])[icoor] /   2 +   4*vdiff_perp(coord[fac2poi[ielem][   3]], coord[fac2poi[ielem][   5]])[icoor] /   2 +   4*vdiff_perp(coord[fac2poi[ielem][   5]], coord[fac2poi[ielem][   0]])[icoor] /   2;
    }

  } 
 
} // End namespace
