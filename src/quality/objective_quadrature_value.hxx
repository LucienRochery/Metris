// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_OBJECTIVE_QUADRATURE_VALUE__
#define __METRIS_OBJECTIVE_QUADRATURE_VALUE__

#include "aux_volumeMeasure.hxx"
#include "objective_quadrature_sample.hxx"
#include "pointwise_objective.hxx"
#include "quafun_stepDistance.hxx"
#include "simplex_quadrature.hxx"

namespace Metris
{

template<int gdim, typename ftype>
struct ObjectiveQuadratureValueEvaluation
{
  static_assert(gdim == 2 || gdim == 3);

  ftype psi{};
  ftype additive_value{};
  bool has_additive_value = false;
  bool reject_element = false;
  ftype rejection_value{};
};

// All paper-aligned objective variants share the same integration measure.
// CavityTargetAverage is deliberately a normalized reference-space average,
// not an instance of that integral, and therefore remains the sole exception.
template<QuaFun iquaf, class MFT>
ObjectiveQuadratureTheta objective_quadrature_theta_mode(Mesh<MFT> &msh)
{
  static_assert(iquaf == QuaFun::SizeShape
                || iquaf == QuaFun::StepDistance);
  if constexpr(iquaf == QuaFun::StepDistance){
    const bool use_target_average
        = msh.param->step_distance_cavity_target_average;
    METRIS_ENFORCE_MSG(
        !(msh.param->step_distance_shape_volume && use_target_average),
        "Step Distance Shape Volume is a distinct integration variant");
    if(use_target_average){
      return ObjectiveQuadratureTheta::ReferenceAverage;
    }
  }
  return standard_objective_quadrature_theta_mode();
}

template<class MFT, int gdim, int tdim, int mshdeg,
         QuaFun iquaf, typename ftype>
struct ObjectiveQuadratureValuePolicy;

template<class MFT, int gdim, int tdim, int mshdeg, typename ftype>
struct ObjectiveQuadratureValuePolicy<
    MFT,gdim,tdim,mshdeg,QuaFun::SizeShape,ftype>
{
  using Sample = ObjectiveQuadratureSample<gdim,tdim,mshdeg>;
  using Evaluation = ObjectiveQuadratureValueEvaluation<gdim,ftype>;

  static Evaluation evaluate(Mesh<MFT> &msh,
                             AsDeg asdmsh,
                             AsDeg asdmet,
                             const int *ent2poi,
                             const Sample &sample)
  {
    const PointwiseObjectiveResult<gdim,ftype> pointwise_result
        = evaluate_pointwise_objective_value<
              MFT,gdim,tdim,QuaFun::SizeShape,ftype>(
                  msh,asdmsh,asdmet,ent2poi,
                  sample.barycentric_coordinates.data(),
                  sample.metric.data());

    Evaluation evaluation;
    evaluation.psi = pointwise_result.psi;
    return evaluation;
  }

  static void enforce_theta_is_valid(bool theta_is_valid)
  {
    METRIS_ENFORCE_MSG(theta_is_valid,
                       "Invalid SizeShape quadrature theta");
  }
};

template<class MFT, int gdim, int tdim, int mshdeg, typename ftype>
struct ObjectiveQuadratureValuePolicy<
    MFT,gdim,tdim,mshdeg,QuaFun::StepDistance,ftype>
{
  using Sample = ObjectiveQuadratureSample<gdim,tdim,mshdeg>;
  using Evaluation = ObjectiveQuadratureValueEvaluation<gdim,ftype>;

  static Evaluation evaluate(Mesh<MFT> &msh,
                             AsDeg asdmsh,
                             AsDeg asdmet,
                             const int *ent2poi,
                             const Sample &sample)
  {
    const PointwiseObjectiveResult<gdim,ftype> pointwise_result
        = evaluate_pointwise_objective_value<
              MFT,gdim,tdim,QuaFun::StepDistance,ftype>(
                  msh,asdmsh,asdmet,ent2poi,
                  sample.barycentric_coordinates.data(),
                  sample.metric.data());

    Evaluation evaluation;
    evaluation.psi = pointwise_result.psi;

    if(msh.param->step_distance_shape_volume
       && evaluation.psi
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
    const double barrier_beta = msh.param->step_distance_shape_volume
                              ? 0.0
                              : msh.param->step_distance_barrier_beta;
    VolumeMeasureHelpers::eval_metric_volume_barrier_fixed_metric_grad<
        gdim,tdim,double>(
            sample.regular_jacobian_transpose.data(),
            sample.metric.data(),NULL,
            msh.param->step_distance_barrier_rho0,
            barrier_beta,
            &metric_volume,&barrier,NULL);
    evaluation.additive_value = static_cast<ftype>(barrier);
    evaluation.has_additive_value = true;
    return evaluation;
  }

  static void enforce_theta_is_valid(bool theta_is_valid)
  {
    METRIS_ENFORCE_MSG(theta_is_valid,
                       "Invalid StepDistance quadrature theta");
  }
};

// Common value-only objective traversal. Objective policies provide the
// complete pointwise psi and any integration-level rejection or additive
// contribution. This function alone owns quadrature iteration, geometry and
// metric sample preparation, theta validation, and weighted accumulation.
template<class MFT, int gdim, int tdim, int mshdeg,
         QuaFun iquaf, typename ftype>
ftype integrate_objective_quadrature_value(
    Mesh<MFT> &msh,
    AsDeg asdmsh,
    AsDeg asdmet,
    const int *ent2poi,
    const SimplexQuadratureView<tdim> quadrature)
{
  static_assert(iquaf == QuaFun::SizeShape
                || iquaf == QuaFun::StepDistance);

  using Policy = ObjectiveQuadratureValuePolicy<
      MFT,gdim,tdim,mshdeg,iquaf,ftype>;
  using Sample = ObjectiveQuadratureSample<gdim,tdim,mshdeg>;
  using Evaluation = ObjectiveQuadratureValueEvaluation<gdim,ftype>;

  const ObjectiveQuadratureTheta theta_mode
      = objective_quadrature_theta_mode<iquaf>(msh);
  ftype integral = ftype(0);

  for(int iquad = 0; iquad < quadrature.size(); iquad++){
    const SimplexQuadraturePointView<tdim> quadrature_point
        = quadrature[iquad];
    const Sample sample
        = prepare_objective_quadrature_sample<MFT,gdim,tdim,mshdeg>(
              msh,asdmet,ent2poi,quadrature_point,theta_mode);
    const Evaluation evaluation
        = Policy::evaluate(msh,asdmsh,asdmet,ent2poi,sample);

    // ShapeVolume must retain its finite rejection sentinel even when the
    // rejected map also makes theta invalid.
    if(evaluation.reject_element){
      return evaluation.rejection_value;
    }
    Policy::enforce_theta_is_valid(sample.theta_is_valid);

    const ftype quadrature_weight
        = static_cast<ftype>(sample.quadrature_weight);
    const ftype theta = static_cast<ftype>(sample.theta);
    if(evaluation.has_additive_value){
      integral += quadrature_weight
                *(theta*evaluation.psi + evaluation.additive_value);
    }else{
      const ftype objective_weight
          = static_cast<ftype>(sample.quadrature_weight*sample.theta);
      integral += objective_weight*evaluation.psi;
    }
  }

  return integral;
}

} // namespace Metris

#endif // __METRIS_OBJECTIVE_QUADRATURE_VALUE__
