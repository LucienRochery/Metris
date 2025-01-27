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


// Wrappers that may call the codegen routines 
template<int gdim, int tdim, int ideg>
void getccoef(const MeshBase &msh, int ientt, double *nrmal, double *ccoef);
template<int gdim, int tdim, int ideg>
void getsclccoef(const MeshBase &msh, int ientt, double *nrmal, 
                 double *ccoef, bool *iinva);

template<int idim, int ideg>
void getccoef_dcoord(const MeshBase &msh, int ientt, int icoor, double *ccoef, dblAr2& d_ccoef);

// Ccoefs if mesh is Lagrange
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


} // End namespace

#endif