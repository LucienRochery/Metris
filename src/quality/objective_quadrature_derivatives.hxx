// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_OBJECTIVE_QUADRATURE_DERIVATIVES__
#define __METRIS_OBJECTIVE_QUADRATURE_DERIVATIVES__

#include "objective_quadrature_value.hxx"

#include "../ho_constants.hxx"
#include "../linalg/symidx.hxx"
#include "../low_eval.hxx"

#include <array>

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

template<int gdim, typename ftype>
struct ObjectiveQuadratureDerivativeEvaluation
{
  static_assert(gdim == 2 || gdim == 3);
  static constexpr int nhessian = gdim*(gdim + 1)/2;

  PointwiseObjectiveResult<gdim,ftype> pointwise{};
  ftype additive_value{};
  std::array<ftype,gdim> additive_gradient{};
  std::array<ftype,nhessian> additive_hessian{};
  bool has_additive_value = false;
  bool differentiate_theta = false;
  bool reject_element = false;
  ftype rejection_value{};
};

template<class MFT, int gdim, int tdim, int mshdeg,
         QuaFun iquaf, typename ftype>
struct ObjectiveQuadratureDerivativePolicy;

template<class MFT, int gdim, int tdim, int mshdeg, typename ftype>
struct ObjectiveQuadratureDerivativePolicy<
    MFT,gdim,tdim,mshdeg,QuaFun::SizeShape,ftype>
{
  using ValuePolicy = ObjectiveQuadratureValuePolicy<
      MFT,gdim,tdim,mshdeg,QuaFun::SizeShape,ftype>;
  using Sample = ObjectiveQuadratureSample<gdim,tdim,mshdeg>;
  using Evaluation = ObjectiveQuadratureDerivativeEvaluation<gdim,ftype>;

  static Evaluation evaluate(
      Mesh<MFT> &msh,
      AsDeg asdmsh,
      AsDeg asdmet,
      const int *ent2poi,
      const Sample &sample,
      int local_control_point,
      FEBasis dofbas,
      PointwiseDerivativeOrder derivative_order,
      const double *regular_basis_gradient)
  {
    (void)regular_basis_gradient;
    Evaluation evaluation;
    if(derivative_order == PointwiseDerivativeOrder::Value){
      evaluation.pointwise = evaluate_pointwise_objective_value<
          MFT,gdim,tdim,QuaFun::SizeShape,ftype>(
              msh,asdmsh,asdmet,ent2poi,
              sample.barycentric_coordinates.data(),sample.metric.data());
    }else{
      evaluation.pointwise = evaluate_pointwise_objective_derivatives<
          MFT,gdim,tdim,QuaFun::SizeShape,ftype>(
              msh,asdmsh,asdmet,ent2poi,
              sample.barycentric_coordinates.data(),sample.metric.data(),
              local_control_point,dofbas,DifVar::None,derivative_order);
    }
    return evaluation;
  }

  static void enforce_theta_is_valid(bool theta_is_valid)
  {
    ValuePolicy::enforce_theta_is_valid(theta_is_valid);
  }
};

template<class MFT, int gdim, int tdim, int mshdeg, typename ftype>
struct ObjectiveQuadratureDerivativePolicy<
    MFT,gdim,tdim,mshdeg,QuaFun::StepDistance,ftype>
{
  using ValuePolicy = ObjectiveQuadratureValuePolicy<
      MFT,gdim,tdim,mshdeg,QuaFun::StepDistance,ftype>;
  using Sample = ObjectiveQuadratureSample<gdim,tdim,mshdeg>;
  using Evaluation = ObjectiveQuadratureDerivativeEvaluation<gdim,ftype>;

  static Evaluation evaluate(
      Mesh<MFT> &msh,
      AsDeg asdmsh,
      AsDeg asdmet,
      const int *ent2poi,
      const Sample &sample,
      int local_control_point,
      FEBasis dofbas,
      PointwiseDerivativeOrder derivative_order,
      const double *regular_basis_gradient)
  {
    Evaluation evaluation;
    if(derivative_order == PointwiseDerivativeOrder::Value){
      evaluation.pointwise = evaluate_pointwise_objective_value<
          MFT,gdim,tdim,QuaFun::StepDistance,ftype>(
              msh,asdmsh,asdmet,ent2poi,
              sample.barycentric_coordinates.data(),sample.metric.data());
    }else{
      evaluation.pointwise = evaluate_pointwise_objective_derivatives<
          MFT,gdim,tdim,QuaFun::StepDistance,ftype>(
              msh,asdmsh,asdmet,ent2poi,
              sample.barycentric_coordinates.data(),sample.metric.data(),
              local_control_point,dofbas,DifVar::None,derivative_order);
    }

    if(msh.param->step_distance_shape_volume
       && evaluation.pointwise.psi
          >= ftype(0.5*step_distance_shape_volume_rejection_quality)){
      evaluation.reject_element = true;
      evaluation.rejection_value
          = ftype(step_distance_shape_volume_rejection_quality);
      return evaluation;
    }

    if(msh.param->step_distance_cavity_target_average){
      return evaluation;
    }

    double metric_volume;
    double barrier;
    double barrier_gradient[gdim];
    constexpr int nhessian = gdim*(gdim + 1)/2;
    double barrier_hessian[nhessian];
    const double barrier_beta = msh.param->step_distance_shape_volume
                              ? 0.0
                              : msh.param->step_distance_barrier_beta;
    double *barrier_gradient_output
        = derivative_order == PointwiseDerivativeOrder::Value
        ? NULL
        : barrier_gradient;
    VolumeMeasureHelpers::eval_metric_volume_barrier_fixed_metric_grad<
        gdim,tdim,double>(
            sample.regular_jacobian_transpose.data(),sample.metric.data(),
            regular_basis_gradient,
            msh.param->step_distance_barrier_rho0,barrier_beta,
            &metric_volume,&barrier,barrier_gradient_output);

    evaluation.additive_value = static_cast<ftype>(barrier);
    evaluation.has_additive_value = true;
    if(derivative_order != PointwiseDerivativeOrder::Value){
      for(int icomponent = 0; icomponent < gdim; icomponent++){
        evaluation.additive_gradient[icomponent]
            = static_cast<ftype>(barrier_gradient[icomponent]);
      }
    }
    if(derivative_order == PointwiseDerivativeOrder::Hessian){
      VolumeMeasureHelpers::
          eval_metric_volume_barrier_fixed_metric_hess_by_surreal<
              gdim,tdim>(
                  sample.regular_jacobian_transpose.data(),
                  sample.metric.data(),regular_basis_gradient,
                  msh.param->step_distance_barrier_rho0,barrier_beta,
                  barrier_hessian);
      for(int ihessian = 0; ihessian < nhessian; ihessian++){
        evaluation.additive_hessian[ihessian]
            = static_cast<ftype>(barrier_hessian[ihessian]);
      }
    }

    #ifdef STEPDISTANCE_INCLUDE_GEOM_THETA_DERIV
    evaluation.differentiate_theta
        = derivative_order != PointwiseDerivativeOrder::Value
       && !msh.param->step_distance_shape_volume;
    #endif
    return evaluation;
  }

  static void enforce_theta_is_valid(bool theta_is_valid)
  {
    ValuePolicy::enforce_theta_is_valid(theta_is_valid);
  }
};

// Common differentiated objective traversal. The only geometry derivative
// associated with the active control point is its basis gradient in the
// regular reference frame. Pointwise objective policies remain responsible
// for their complete psi derivatives and any additive/rejection behavior.
template<class MFT, int gdim, int tdim, int mshdeg,
         QuaFun iquaf, typename ftype>
ftype integrate_objective_quadrature_derivatives(
    Mesh<MFT> &msh,
    AsDeg asdmsh,
    AsDeg asdmet,
    const int *ent2poi,
    const SimplexQuadratureView<tdim> quadrature,
    int local_control_point,
    FEBasis dofbas,
    ftype *derivative,
    ftype *hessian)
{
  static_assert(iquaf == QuaFun::SizeShape
                || iquaf == QuaFun::StepDistance);
  constexpr int nhessian = gdim*(gdim + 1)/2;

  using Policy = ObjectiveQuadratureDerivativePolicy<
      MFT,gdim,tdim,mshdeg,iquaf,ftype>;
  using Sample = ObjectiveQuadratureSample<gdim,tdim,mshdeg>;
  using Evaluation = ObjectiveQuadratureDerivativeEvaluation<gdim,ftype>;

  const bool differentiate = local_control_point >= 0;
  const PointwiseDerivativeOrder derivative_order
      = !differentiate
      ? PointwiseDerivativeOrder::Value
      : hessian == NULL
        ? PointwiseDerivativeOrder::Gradient
        : PointwiseDerivativeOrder::Hessian;
  if(differentiate){
    METRIS_ENFORCE(derivative != NULL);
    for(int icomponent = 0; icomponent < gdim; icomponent++){
      derivative[icomponent] = ftype(0);
    }
    if(hessian != NULL){
      for(int ihessian = 0; ihessian < nhessian; ihessian++){
        hessian[ihessian] = ftype(0);
      }
    }
  }

  const ObjectiveQuadratureTheta theta_mode
      = objective_quadrature_theta_mode<iquaf>(msh);
  ftype integral = ftype(0);

  for(int iquad = 0; iquad < quadrature.size(); iquad++){
    const SimplexQuadraturePointView<tdim> quadrature_point
        = quadrature[iquad];
    const Sample sample
        = prepare_objective_quadrature_sample<MFT,gdim,tdim,mshdeg>(
              msh,asdmet,ent2poi,quadrature_point,theta_mode);

    std::array<double,tdim> regular_basis_gradient{};
    const double *regular_basis_gradient_data = NULL;
    if(differentiate){
      regular_basis_gradient
          = objective_regular_basis_gradient<tdim,mshdeg>(
                dofbas,local_control_point,
                sample.barycentric_coordinates.data());
      regular_basis_gradient_data = regular_basis_gradient.data();
    }

    const Evaluation evaluation = Policy::evaluate(
        msh,asdmsh,asdmet,ent2poi,sample,
        local_control_point,dofbas,derivative_order,
        regular_basis_gradient_data);

    // ShapeVolume must retain its finite rejection sentinel even when the
    // rejected map also makes theta invalid.
    if(evaluation.reject_element){
      return evaluation.rejection_value;
    }
    Policy::enforce_theta_is_valid(sample.theta_is_valid);

    double theta_value = sample.theta;
    double theta_gradient[gdim]{};
    double theta_hessian[nhessian]{};
    if(evaluation.differentiate_theta){
      VolumeMeasureHelpers::eval_theta_fixed_metric_grad<gdim,tdim,double>(
          sample.regular_jacobian_transpose.data(),sample.metric.data(),
          regular_basis_gradient_data,&theta_value,theta_gradient);
      if(derivative_order == PointwiseDerivativeOrder::Hessian){
        VolumeMeasureHelpers::
            eval_theta_fixed_metric_hess_by_surreal<gdim,tdim>(
                sample.regular_jacobian_transpose.data(),sample.metric.data(),
                regular_basis_gradient_data,theta_hessian);
      }
    }

    const ftype quadrature_weight
        = static_cast<ftype>(sample.quadrature_weight);
    const ftype theta = static_cast<ftype>(theta_value);
    const ftype objective_weight
        = static_cast<ftype>(sample.quadrature_weight*theta_value);

    if(evaluation.has_additive_value){
      integral += quadrature_weight
                *(theta*evaluation.pointwise.psi
                  + evaluation.additive_value);
    }else{
      integral += objective_weight*evaluation.pointwise.psi;
    }
    if(!differentiate) continue;

    if(evaluation.has_additive_value
       || evaluation.differentiate_theta){
      for(int icomponent = 0; icomponent < gdim; icomponent++){
        derivative[icomponent] += quadrature_weight*(
            theta*evaluation.pointwise.gradient[icomponent]
          + evaluation.pointwise.psi
           * static_cast<ftype>(theta_gradient[icomponent])
          + evaluation.additive_gradient[icomponent]);
      }
    }else{
      for(int icomponent = 0; icomponent < gdim; icomponent++){
        derivative[icomponent]
            += objective_weight*evaluation.pointwise.gradient[icomponent];
      }
    }
    if(derivative_order != PointwiseDerivativeOrder::Hessian) continue;

    if(evaluation.has_additive_value
       || evaluation.differentiate_theta){
      for(int i = 0; i < gdim; i++){
        for(int j = i; j < gdim; j++){
          const int ihessian = sym2idx(i,j);
          hessian[ihessian] += quadrature_weight*(
              theta*evaluation.pointwise.hessian[ihessian]
            + evaluation.pointwise.psi
             * static_cast<ftype>(theta_hessian[ihessian])
            + static_cast<ftype>(theta_gradient[i])
             * evaluation.pointwise.gradient[j]
            + evaluation.pointwise.gradient[i]
             * static_cast<ftype>(theta_gradient[j])
            + evaluation.additive_hessian[ihessian]);
        }
      }
    }else{
      for(int ihessian = 0; ihessian < nhessian; ihessian++){
        hessian[ihessian]
            += objective_weight*evaluation.pointwise.hessian[ihessian];
      }
    }
  }

  return integral;
}

} // namespace Metris

#endif // __METRIS_OBJECTIVE_QUADRATURE_DERIVATIVES__
