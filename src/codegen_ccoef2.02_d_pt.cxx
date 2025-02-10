//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "codegen_ccoef_d.hxx"

#include "types.hxx"

namespace Metris{

template<int ideg>
void d_pt_ccoef_genbez2([[maybe_unused]] const intAr2&__restrict__ fac2poi, [[maybe_unused]] const dblAr2&__restrict__ coord,
                        [[maybe_unused]] int ielem, [[maybe_unused]] int inode,
                        [[maybe_unused]] dblAr2&__restrict__ d_ccoef){}

void vdiff_perp(const double* a, const double* b, int up, int lo, double *res);

void vdiff_perp_sum(const double* a, const double* b, int up, int lo, double *res);

template<> void d_pt_ccoef_genbez2<2>(const intAr2&__restrict__ fac2poi, const dblAr2&__restrict__ coord,
                                        int ielem, int inode,
                                        dblAr2&__restrict__ d_ccoef){

  d_ccoef.fill(6,2,0);

  if(inode ==   0){
    vdiff_perp(coord[fac2poi(ielem,5)], coord[fac2poi(ielem,4)],4,1,d_ccoef[0]);
    vdiff_perp(coord[fac2poi(ielem,3)], coord[fac2poi(ielem,2)],2,1,d_ccoef[4]);
    vdiff_perp(coord[fac2poi(ielem,1)], coord[fac2poi(ielem,3)],2,1,d_ccoef[5]);
  }
  else if(inode ==   1){
    vdiff_perp(coord[fac2poi(ielem,3)], coord[fac2poi(ielem,5)],4,1,d_ccoef[1]);
    vdiff_perp(coord[fac2poi(ielem,2)], coord[fac2poi(ielem,4)],2,1,d_ccoef[3]);
    vdiff_perp(coord[fac2poi(ielem,4)], coord[fac2poi(ielem,0)],2,1,d_ccoef[5]);
  }
  else if(inode ==   2){
    vdiff_perp(coord[fac2poi(ielem,4)], coord[fac2poi(ielem,3)],4,1,d_ccoef[2]);
    vdiff_perp(coord[fac2poi(ielem,5)], coord[fac2poi(ielem,1)],2,1,d_ccoef[3]);
    vdiff_perp(coord[fac2poi(ielem,0)], coord[fac2poi(ielem,5)],2,1,d_ccoef[4]);
  }
  else if(inode ==   3){
    vdiff_perp(coord[fac2poi(ielem,5)], coord[fac2poi(ielem,1)],4,1,d_ccoef[1]);
    vdiff_perp(coord[fac2poi(ielem,2)], coord[fac2poi(ielem,4)],4,1,d_ccoef[2]);
    vdiff_perp(coord[fac2poi(ielem,4)], coord[fac2poi(ielem,5)],2,1,d_ccoef[3]);
    vdiff_perp(coord[fac2poi(ielem,4)], coord[fac2poi(ielem,0)],2,1,d_ccoef[4]);
    vdiff_perp(coord[fac2poi(ielem,0)], coord[fac2poi(ielem,5)],2,1,d_ccoef[5]);
  }
  else if(inode ==   4){
    vdiff_perp(coord[fac2poi(ielem,0)], coord[fac2poi(ielem,5)],4,1,d_ccoef[0]);
    vdiff_perp(coord[fac2poi(ielem,3)], coord[fac2poi(ielem,2)],4,1,d_ccoef[2]);
    vdiff_perp(coord[fac2poi(ielem,1)], coord[fac2poi(ielem,3)],2,1,d_ccoef[3]);
    vdiff_perp(coord[fac2poi(ielem,5)], coord[fac2poi(ielem,3)],2,1,d_ccoef[4]);
    vdiff_perp(coord[fac2poi(ielem,5)], coord[fac2poi(ielem,1)],2,1,d_ccoef[5]);
  }
  else if(inode ==   5){
    vdiff_perp(coord[fac2poi(ielem,4)], coord[fac2poi(ielem,0)],4,1,d_ccoef[0]);
    vdiff_perp(coord[fac2poi(ielem,1)], coord[fac2poi(ielem,3)],4,1,d_ccoef[1]);
    vdiff_perp(coord[fac2poi(ielem,3)], coord[fac2poi(ielem,2)],2,1,d_ccoef[3]);
    vdiff_perp(coord[fac2poi(ielem,2)], coord[fac2poi(ielem,4)],2,1,d_ccoef[4]);
    vdiff_perp(coord[fac2poi(ielem,3)], coord[fac2poi(ielem,4)],2,1,d_ccoef[5]);
  }

} 
 
} // End namespace
