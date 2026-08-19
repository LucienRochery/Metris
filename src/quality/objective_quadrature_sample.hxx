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

// The sample-preparation layer supports the integration conventions that are
// already present while exposing one explicit theta to its consumer. The
// objective paths can converge on one convention later without duplicating
// geometry or metric evaluation again.
enum class ObjectiveQuadratureTheta
{
  ReferenceAverage,
  PhysicalMeasure,
  PhysicalMetricMeasure,
  RegularMetricMeasure
};

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

// Prepare one frozen objective sample for a degree-mshdeg element map. Metric
// derivatives are deliberately disabled. Exact reference vertices retain the
// historical direct nodal metric sampling; every other FE sample is evaluated
// barycentrically and every other analytical sample at its physical point.
template<class MFT, int gdim, int tdim, int mshdeg>
ObjectiveQuadratureSample<gdim,tdim,mshdeg>
prepare_objective_quadrature_sample(
    Mesh<MFT> &msh,
    AsDeg asdmet,
    const int *ent2poi,
    const SimplexQuadraturePointView<tdim> quadrature_point,
    ObjectiveQuadratureTheta theta_mode)
{
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
  if(reference_vertex >= 0){
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

} // namespace Metris

#endif // __METRIS_OBJECTIVE_QUADRATURE_SAMPLE__
