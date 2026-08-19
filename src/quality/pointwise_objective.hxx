// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_POINTWISE_OBJECTIVE__
#define __METRIS_POINTWISE_OBJECTIVE__

#include "quafun.hxx"

#include "../Mesh/MeshFwd.hxx"
#include "../aux_exceptions.hxx"

#include <array>

namespace Metris
{

enum class PointwiseDerivativeOrder
{
  Value,
  Gradient,
  Hessian
};

template<int gdim, typename ftype = double>
struct PointwiseObjectiveResult
{
  static_assert(gdim == 2 || gdim == 3);
  static constexpr int nhessian = gdim*(gdim + 1)/2;

  ftype psi{};
  std::array<ftype,gdim> gradient{};
  std::array<ftype,nhessian> hessian{};
  PointwiseDerivativeOrder derivative_order = PointwiseDerivativeOrder::Value;

  bool has_gradient() const
  {
    return derivative_order == PointwiseDerivativeOrder::Gradient
        || derivative_order == PointwiseDerivativeOrder::Hessian;
  }

  bool has_hessian() const
  {
    return derivative_order == PointwiseDerivativeOrder::Hessian;
  }
};

template<class MFT, int gdim, int tdim, QuaFun iquaf,
         typename ftype = double>
PointwiseObjectiveResult<gdim,ftype> evaluate_pointwise_objective_value(
    Mesh<MFT> &msh,
    AsDeg asdmsh,
    AsDeg asdmet,
    const int *ent2poi,
    const double *bary,
    const double *metric)
{
  static_assert(iquaf == QuaFun::SizeShape
                || iquaf == QuaFun::StepDistance);

  PointwiseObjectiveResult<gdim,ftype> result;
  result.psi = get_quafun_xi<MFT,gdim,tdim,iquaf,ftype>()(
      msh,asdmsh,asdmet,ent2poi,bary,metric);
  return result;
}

template<class MFT, int gdim, int tdim, QuaFun iquaf,
         typename ftype = double>
PointwiseObjectiveResult<gdim,ftype>
evaluate_pointwise_objective_derivatives(
    Mesh<MFT> &msh,
    AsDeg asdmsh,
    AsDeg asdmet,
    const int *ent2poi,
    const double *bary,
    const double *metric,
    int ivar,
    FEBasis dofbas,
    DifVar idifmet,
    PointwiseDerivativeOrder derivative_order)
{
  static_assert(iquaf == QuaFun::SizeShape
                || iquaf == QuaFun::StepDistance);
  METRIS_ENFORCE_MSG(
      derivative_order == PointwiseDerivativeOrder::Gradient
      || derivative_order == PointwiseDerivativeOrder::Hessian,
      "Differentiated pointwise objective requires gradient or Hessian order");

  PointwiseObjectiveResult<gdim,ftype> result;
  result.derivative_order = derivative_order;
  ftype *hessian_output
      = derivative_order == PointwiseDerivativeOrder::Hessian
      ? result.hessian.data()
      : NULL;
  result.psi = get_d_quafun_xi<MFT,gdim,tdim,iquaf,ftype>()(
      msh,asdmsh,asdmet,ent2poi,bary,metric,
      ivar,dofbas,idifmet,result.gradient.data(),hessian_output);
  return result;
}

} // namespace Metris

#endif
