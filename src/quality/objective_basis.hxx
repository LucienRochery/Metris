// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#ifndef __METRIS_OBJECTIVE_BASIS__
#define __METRIS_OBJECTIVE_BASIS__

#include "../aux_exceptions.hxx"
#include "../ho_constants.hxx"
#include "../low_eval.hxx"

#include <array>
#include <boost/hana/at_key.hpp>

namespace Metris
{

// Gradient of one geometry basis function in the regular reference frame.
// It is constant for P1 and varies with the quadrature point for higher-order
// Lagrange nodes and Bezier control points.
template<int tdim, int mshdeg>
std::array<double,tdim> objective_regular_basis_gradient(
    FEBasis dofbas,
    int local_control_point,
    const double *barycentric_coordinates)
{
  static_assert(tdim == 2 || tdim == 3);
  static_assert(mshdeg >= 1 && mshdeg <= METRIS_MAX_DEG);
  METRIS_ENFORCE_MSG(
      dofbas == FEBasis::Lagrange || dofbas == FEBasis::Bezier,
      "Objective derivative requires a Lagrange or Bezier control point");
  METRIS_ENFORCE_MSG(
      local_control_point >= 0
      && local_control_point < getnnode(tdim,mshdeg),
      "Objective derivative control point {} outside degree-{} element",
      local_control_point,mshdeg);

  constexpr auto ordering = ORDELT(tdim);
  const int *multi_index = ordering[mshdeg][local_control_point];
  double canonical_gradient[tdim];
  if(dofbas == FEBasis::Bezier){
    eval_bezierfunc<mshdeg,tdim>(
        multi_index,barycentric_coordinates,1,canonical_gradient);
  }else{
    eval_lagrangefunc<mshdeg,tdim>(
        multi_index,barycentric_coordinates,1,canonical_gradient);
  }

  std::array<double,tdim> regular_gradient{};
  for(int iregular = 0; iregular < tdim; iregular++){
    for(int icanonical = 0; icanonical < tdim; icanonical++){
      regular_gradient[iregular]
          += canonical_gradient[icanonical]
           * Constants::invtJ_0[hana::type_c<double>][tdim]
                              [iregular*tdim + icanonical];
    }
  }
  return regular_gradient;
}

} // namespace Metris

#endif // __METRIS_OBJECTIVE_BASIS__
