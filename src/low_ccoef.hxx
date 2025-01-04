//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __SRC_LOWCCOEF__
#define __SRC_LOWCCOEF__

#include "ho_constants.hxx"
#include "types.hxx"
#include <array>


namespace Metris{


class MeshBase;

/* 
Bézier coefficients of the Jacobian determinant
*/


template<int gdim, int tdim, int ideg>
void getccoef(const MeshBase &msh, int ientt, double *nrmal, double *ccoef);
template<int gdim, int tdim, int ideg>
void getsclccoef(const MeshBase &msh, int ientt, double *nrmal, 
                 double *ccoef, bool *iinva);

template<int gdim, int tdim, int ideg>
void ccoef_eval(FEBasis ibasis, const intAr2& ent2poi, const dblAr2& coord, int ientt, double *nrmal, double* ccoef);


// Sized for highest (tet) so applicable to all
namespace ccoef_eval_lfld{
  constexpr std::array<int,tetnpps[METRIS_MAX_DEG_JACOBIAN]> lfld{[]() constexpr{
    std::array<int,tetnpps[METRIS_MAX_DEG_JACOBIAN]> ret{};
    for(int i=0;i<tetnpps[METRIS_MAX_DEG_JACOBIAN];i++){
      ret[i] = i;
    }
    return ret;
  }()};
}


//template<int irnk1_1,int irnk1_2,int irnk2_1,int irnk2_2,int irnk3_1,int irnk3_2>
//double det3_vdif_poi(int ielem, intAr2& __restrict__ tet2poi, dblAr2& __restrict__ coord);

} // End namespace

#endif