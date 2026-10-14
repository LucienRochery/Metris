// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_OBJECTIVE_QUADRATURE_SAMPLE__
#define __METRIS_OBJECTIVE_QUADRATURE_SAMPLE__

#include "../Mesh/Mesh.hxx"
#include "../aux_exceptions.hxx"
#include "../linalg/det.hxx"
#include "../low_eval.hxx"
#include "../metris_constants.hxx"

#include "aux_volumeMeasure.hxx"
#include "simplex_quadrature.hxx"

#include <array>
#include <cmath>
#include <type_traits>

namespace Metris
{

// The sample-preparation layer supports the integration conventions needed by
// focused tests and compatibility code while exposing one explicit theta to
// its consumer. Production objective traversal selects its mode through the
// common helper below, independently of the pointwise objective.
enum class ObjectiveQuadratureTheta
{
  ReferenceAverage,
  PhysicalMeasure,
  PhysicalMetricMeasure,
  RegularMetricMeasure
};

constexpr ObjectiveQuadratureTheta
standard_objective_quadrature_theta_mode() noexcept
{
  #ifdef INTQUALINRIEMSPACE
  return ObjectiveQuadratureTheta::PhysicalMetricMeasure;
  #else
  return ObjectiveQuadratureTheta::PhysicalMeasure;
  #endif
}

template<int gdim, int tdim, int mshdeg>
struct ObjectiveQuadratureSample
{
  static_assert(gdim == 2 || gdim == 3);
  static_assert(tdim == 2 || tdim == 3);
  static_assert(tdim <= gdim);
  static_assert(mshdeg >= 1 && mshdeg <= METRIS_MAX_DEG);

  static constexpr int geometry_degree = mshdeg;
  static constexpr int nmetric = gdim*(gdim + 1)/2;

  std::array<double,tdim + 1> barycentric_coordinates{};
  double quadrature_weight = 0.0;
  std::array<double,gdim> physical_coordinates{};
  // Both Jacobians use the transposed storage convention: tdim x gdim.
  std::array<double,tdim*gdim> canonical_jacobian_transpose{};
  std::array<double,tdim*gdim> regular_jacobian_transpose{};
  std::array<double,nmetric> metric{};
  double theta = 1.0;
  bool theta_is_valid = true;
};

namespace detail
{

template<int tdim>
int objective_quadrature_vertex(const double *barycentric_coordinates)
{
  for(int candidate = 0; candidate < tdim + 1; candidate++){
    bool is_vertex = barycentric_coordinates[candidate] == 1.0;
    for(int ibary = 0; ibary < tdim + 1 && is_vertex; ibary++){
      if(ibary == candidate) continue;
      is_vertex = barycentric_coordinates[ibary] == 0.0;
    }
    if(is_vertex) return candidate;
  }
  return -1;
}

template<int gdim, int tdim>
bool objective_jacobian_measure(
    const double *jacobian_transpose,
    double reference_measure,
    double *measure)
{
  // For a full-dimensional map, compute |det(J)| directly. Forming J*J^T
  // first is mathematically equivalent, but squares the condition number and
  // can round a positive determinant to zero for valid anisotropic elements.
  if constexpr(gdim == tdim){
    const double jacobian_determinant
        = VolumeMeasureHelpers::det_full<tdim,double>(jacobian_transpose);
    const double absolute_jacobian_determinant
        = std::abs(jacobian_determinant);
    if(!(absolute_jacobian_determinant > 0.0)
       || !std::isfinite(absolute_jacobian_determinant)){
      *measure = 0.0;
      return false;
    }
    *measure = reference_measure*absolute_jacobian_determinant;
    if(!(*measure > 0.0) || !std::isfinite(*measure)){
      *measure = 0.0;
      return false;
    }
    return true;
  }

  double gram[tdim*tdim];
  for(int i = 0; i < tdim; i++){
    for(int j = 0; j < tdim; j++){
      gram[i*tdim + j] = 0.0;
      for(int component = 0; component < gdim; component++){
        gram[i*tdim + j]
            += jacobian_transpose[i*gdim + component]
             * jacobian_transpose[j*gdim + component];
      }
    }
  }

  const double gram_determinant
      = VolumeMeasureHelpers::det_full<tdim,double>(gram);
  if(!(gram_determinant > 0.0) || !std::isfinite(gram_determinant)){
    *measure = 0.0;
    return false;
  }
  *measure = reference_measure*std::sqrt(gram_determinant);
  return true;
}

} // namespace detail

// Prepare one frozen objective sample for a degree-mshdeg entity map. Metric
// derivatives are deliberately disabled. Exact reference vertices retain the
// historical direct nodal metric sampling; every other FE sample is evaluated
// barycentrically and every other analytical sample at its physical point.
// An explicit owner selects entity-aware sampling in SPD space, including at
// vertices. Owner-free calls retain the historical continuous-field behavior.
template<class MFT, int gdim, int tdim, int mshdeg>
ObjectiveQuadratureSample<gdim,tdim,mshdeg>
prepare_objective_quadrature_sample(
    Mesh<MFT> &msh,
    AsDeg asdmet,
    const int *ent2poi,
    const SimplexQuadraturePointView<tdim> quadrature_point,
    ObjectiveQuadratureTheta theta_mode,
    int ientt = -1)
{
  METRIS_ENFORCE(ientt >= -1);
  if(ientt >= 0){
    METRIS_ENFORCE(gdim == tdim && msh.idim == gdim && msh.get_tdim() == tdim);
    METRIS_ENFORCE(msh.curdeg == mshdeg && asdmet == AsDeg::Pk);
    METRIS_ENFORCE(ientt < msh.nentt(tdim) && msh.ent2poi(tdim)(ientt,0) >= 0);
    for(int node = 0; node < getnnode(tdim,mshdeg); node++)
      METRIS_ENFORCE(ent2poi[node] == msh.ent2poi(tdim)(ientt,node));
  }else if constexpr(std::is_same_v<MFT,MetricFieldFE>){
    METRIS_ENFORCE_MSG(msh.met.getLayout() == MetricFieldLayout::Continuous,
                      "DG objective quadrature requires an owning entity");
  }
  ObjectiveQuadratureSample<gdim,tdim,mshdeg> sample;
  for(int ibary = 0; ibary < tdim + 1; ibary++){
    sample.barycentric_coordinates[ibary] = quadrature_point.bary[ibary];
  }
  sample.quadrature_weight = quadrature_point.weight;

  if constexpr(tdim == 2){
    eval2<gdim,mshdeg>(
        msh.coord,ent2poi,msh.getBasis(),
        DifVar::Bary,DifVar::None,
        sample.barycentric_coordinates.data(),
        sample.physical_coordinates.data(),
        sample.canonical_jacobian_transpose.data(),NULL);
  }else{
    eval3<gdim,mshdeg>(
        msh.coord,ent2poi,msh.getBasis(),
        DifVar::Bary,DifVar::None,
        sample.barycentric_coordinates.data(),
        sample.physical_coordinates.data(),
        sample.canonical_jacobian_transpose.data(),NULL);
  }

  for(int i = 0; i < tdim; i++){
    for(int component = 0; component < gdim; component++){
      sample.regular_jacobian_transpose[i*gdim + component] = 0.0;
      for(int k = 0; k < tdim; k++){
        sample.regular_jacobian_transpose[i*gdim + component]
            += Constants::invtJ_0[hana::type_c<double>][tdim][i*tdim + k]
             * sample.canonical_jacobian_transpose[k*gdim + component];
      }
    }
  }

  const int reference_vertex = detail::objective_quadrature_vertex<tdim>(
      sample.barycentric_coordinates.data());
  if(ientt >= 0){
    if constexpr(std::is_same_v<MFT,MetricFieldAnalytical>){
      msh.met.getMetPhys(DifVar::None,MetSpace::Exp,
          sample.physical_coordinates.data(),sample.metric.data(),nullptr);
    }else{
      msh.met.getMetOnEntity(ientt,DifVar::None,
          sample.barycentric_coordinates.data(),sample.metric.data(),nullptr);
    }
  }else if(reference_vertex >= 0){
    const int ipoin = ent2poi[reference_vertex];
    for(int imetric = 0; imetric < sample.nmetric; imetric++){
      sample.metric[imetric] = msh.met(ipoin,imetric);
    }
  }else if constexpr(
      std::is_same<MFT,MetricFieldAnalytical>::value){
    msh.met.getMetPhys(
        DifVar::None,msh.met.getSpace(),
        sample.physical_coordinates.data(),sample.metric.data(),NULL);
  }else{
    msh.met.getMetBary(
        asdmet,DifVar::None,msh.met.getSpace(),
        ent2poi,tdim,sample.barycentric_coordinates.data(),
        sample.metric.data(),NULL);
  }

  if(theta_mode == ObjectiveQuadratureTheta::ReferenceAverage){
    sample.theta = 1.0;
    return sample;
  }

  if(theta_mode == ObjectiveQuadratureTheta::RegularMetricMeasure){
    sample.theta_is_valid = detail::objective_jacobian_measure<gdim,tdim>(
        sample.regular_jacobian_transpose.data(),1.0,&sample.theta);
    const double metric_determinant = detsym<gdim>(sample.metric.data());
    if(!(metric_determinant > 0.0)
       || !std::isfinite(metric_determinant)){
      sample.theta = 0.0;
      sample.theta_is_valid = false;
      return sample;
    }
    sample.theta *= std::sqrt(metric_determinant);
    return sample;
  }

  constexpr double canonical_reference_measure
      = tdim == 2 ? 0.5 : 1.0/6.0;
  sample.theta_is_valid = detail::objective_jacobian_measure<gdim,tdim>(
      sample.canonical_jacobian_transpose.data(),
      canonical_reference_measure,&sample.theta);
  if(theta_mode == ObjectiveQuadratureTheta::PhysicalMeasure){
    return sample;
  }

  METRIS_ASSERT(theta_mode
                == ObjectiveQuadratureTheta::PhysicalMetricMeasure);
  const double metric_determinant = detsym<gdim>(sample.metric.data());
  if(!(metric_determinant > 0.0)
     || !std::isfinite(metric_determinant)){
    sample.theta = 0.0;
    sample.theta_is_valid = false;
    return sample;
  }
  sample.theta *= std::sqrt(metric_determinant);
  return sample;
}

// Entity-aware objective samples always contain an SPD metric. The shared
// point buffer is not used to select a trace or determine the metric measure.
template<class MFT, int gdim, int tdim, int mshdeg>
ObjectiveQuadratureSample<gdim,tdim,mshdeg>
prepare_objective_quadrature_sample(
    Mesh<MFT> &msh, AsDeg asdmet, int ientt,
    const SimplexQuadraturePointView<tdim> quadrature_point,
    ObjectiveQuadratureTheta theta_mode)
{
  METRIS_ENFORCE(ientt >= 0 && ientt < msh.nentt(tdim));
  return prepare_objective_quadrature_sample<MFT,gdim,tdim,mshdeg>(
      msh,asdmet,msh.ent2poi(tdim)[ientt],quadrature_point,theta_mode,ientt);
}

} // namespace Metris

#endif // __METRIS_OBJECTIVE_QUADRATURE_SAMPLE__
