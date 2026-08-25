// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_high_order_stepdistance_value

#include <boost/test/included/unit_test.hpp>

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "ho_constants.hxx"
#include "quality/low_metqua.hxx"
#include "quality/low_metqua_d.hxx"
#include "quality/objective_quadrature_derivatives.hxx"
#include "quality/pointwise_objective.hxx"
#include "quality/objective_quadrature_value.hxx"
#include "quality/simplex_quadrature.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>

namespace Metris
{

namespace
{

template<int gdim>
void identity_analytical_metric(
    const AnaMetCtx *, const double *, double scale,
    int derivative_order, double *metric, double *metric_derivative)
{
  constexpr int metric_count = gdim*(gdim + 1)/2;
  std::fill(metric,metric + metric_count,0.0);
  for(int component = 0; component < gdim; component++){
    metric[component*(component + 1)/2 + component] = scale;
  }
  if(derivative_order == 0) return;
  std::fill(metric_derivative,
            metric_derivative + gdim*metric_count,0.0);
}

template<class MFT, int gdim>
void initialize_curved_p2_stepdistance_element(
    Mesh<MFT> &mesh,
    MetrisParameters &parameters)
{
  static_assert(gdim == 2 || gdim == 3);
  constexpr int node_count = getnnode(gdim,2);
  constexpr int metric_count = gdim*(gdim + 1)/2;
  constexpr auto ordering = ORDELT(gdim);

  mesh.idim = gdim;
  mesh.curdeg = 2;
  mesh.strdeg = 2;
  mesh.forceBasisFlag(FEBasis::Lagrange);
  mesh.param = &parameters;
  mesh.set_npoin(node_count);
  mesh.set_nentt(gdim,1);
  mesh.met.forceBasisFlag(FEBasis::Lagrange);
  mesh.met.forceSpaceFlag(MetSpace::Exp);

  constexpr double vertices_2d[3][2] = {
      {0.05,0.08},{1.12,0.14},{0.18,0.96}};
  constexpr double offsets_2d[6][2] = {
      {0.0,0.0},{0.0,0.0},{0.0,0.0},
      {0.035,-0.025},{-0.020,0.040},{0.025,0.050}};
  constexpr double vertices_3d[4][3] = {
      {0.04,0.06,0.03},{1.08,0.12,0.02},
      {0.16,0.98,0.10},{0.10,0.20,0.93}};
  constexpr double offsets_3d[10][3] = {
      {0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0},
      {0.018,-0.012,0.010},{-0.012,0.022,0.015},
      {0.014,0.010,-0.016},{-0.018,0.012,0.020},
      {0.010,-0.016,0.012},{-0.014,0.018,-0.010}};

  intAr2 &element_to_point = mesh.ent2poi(gdim);
  for(int inode = 0; inode < node_count; inode++){
    element_to_point(0,inode) = inode;
    for(int component = 0; component < gdim; component++){
      mesh.coord(inode,component) = 0.0;
      for(int vertex = 0; vertex < gdim + 1; vertex++){
        const double barycentric_coordinate
            = ordering[2][inode][vertex]/2.0;
        if constexpr(gdim == 2){
          mesh.coord(inode,component)
              += barycentric_coordinate*vertices_2d[vertex][component];
        }else{
          mesh.coord(inode,component)
              += barycentric_coordinate*vertices_3d[vertex][component];
        }
      }
      if constexpr(gdim == 2){
        mesh.coord(inode,component) += offsets_2d[inode][component];
      }else{
        mesh.coord(inode,component) += offsets_3d[inode][component];
      }
    }
  }

  if constexpr(std::is_same<MFT,MetricFieldAnalytical>::value){
    if constexpr(gdim == 2){
      parameters.setAnalyticalMetric(
          AnaMetFun(identity_analytical_metric<2>));
    }else{
      parameters.setAnalyticalMetric(
          AnaMetFun(identity_analytical_metric<3>));
    }
    mesh.met.setAnalyticalMetric(parameters);
  }

  for(int ipoin = 0; ipoin < node_count; ipoin++){
    for(int imetric = 0; imetric < metric_count; imetric++){
      mesh.met(ipoin,imetric) = 0.0;
    }
    for(int component = 0; component < gdim; component++){
      mesh.met(ipoin,component*(component + 1)/2 + component) = 1.0;
    }
  }
}

template<class MFT, int gdim>
void check_stepdistance_value_geometry_dispatch()
{
  constexpr int metric_count = gdim*(gdim + 1)/2;
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_quadrature_order = 5;
  parameters.objective_p = 1.25;
  parameters.step_distance_regularization = 1.e-6;
  parameters.step_distance_shape_volume = false;
  parameters.step_distance_cavity_target_average = false;
  parameters.step_distance_barrier_beta = 0.0;

  Mesh<MFT> mesh;
  initialize_curved_p2_stepdistance_element<MFT,gdim>(mesh,parameters);
  const int *nodes = mesh.ent2poi(gdim)[0];
  const SimplexQuadratureView<gdim> quadrature
      = get_objective_quadrature<gdim>(
            parameters.objective_quadrature_order);

  const double expected_p2
      = integrate_objective_quadrature_value<
            MFT,gdim,gdim,2,QuaFun::StepDistance,double>(
                mesh,AsDeg::Pk,AsDeg::P1,nodes,quadrature);
  const double dispatched_p2
      = metqua<MFT,gdim,gdim,QuaFun::StepDistance,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,1.0);
  BOOST_CHECK_EQUAL(dispatched_p2,expected_p2);

  const double expected_p1
      = integrate_objective_quadrature_value<
            MFT,gdim,gdim,1,QuaFun::StepDistance,double>(
                mesh,AsDeg::P1,AsDeg::P1,nodes,quadrature);
  const double dispatched_p1
      = metqua<MFT,gdim,gdim,QuaFun::StepDistance,double>(
            mesh,AsDeg::P1,AsDeg::P1,0,1.0);
  BOOST_CHECK_EQUAL(dispatched_p1,expected_p1);
  BOOST_CHECK_GT(
      std::abs(dispatched_p2 - dispatched_p1),
      1.e-5*std::max(1.0,std::abs(dispatched_p2)));

  std::array<double,gdim + 1> barycentric_coordinates{};
  if constexpr(gdim == 2){
    barycentric_coordinates = {0.20,0.30,0.50};
  }else{
    barycentric_coordinates = {0.10,0.20,0.30,0.40};
  }
  std::array<double,metric_count> frozen_metric{};
  for(int component = 0; component < gdim; component++){
    frozen_metric[component*(component + 1)/2 + component] = 1.0;
  }

  const double pointwise_p1_before
      = evaluate_pointwise_objective_value<
            MFT,gdim,gdim,QuaFun::StepDistance,double>(
                mesh,AsDeg::P1,AsDeg::P1,nodes,
                barycentric_coordinates.data(),frozen_metric.data()).psi;
  const double pointwise_p2_before
      = evaluate_pointwise_objective_value<
            MFT,gdim,gdim,QuaFun::StepDistance,double>(
                mesh,AsDeg::Pk,AsDeg::P1,nodes,
                barycentric_coordinates.data(),frozen_metric.data()).psi;

  const int first_edge_node = gdim + 1;
  mesh.coord(first_edge_node,0) += 0.04;
  mesh.coord(first_edge_node,gdim - 1) -= 0.025;

  const double pointwise_p1_after
      = evaluate_pointwise_objective_value<
            MFT,gdim,gdim,QuaFun::StepDistance,double>(
                mesh,AsDeg::P1,AsDeg::P1,nodes,
                barycentric_coordinates.data(),frozen_metric.data()).psi;
  const double pointwise_p2_after
      = evaluate_pointwise_objective_value<
            MFT,gdim,gdim,QuaFun::StepDistance,double>(
                mesh,AsDeg::Pk,AsDeg::P1,nodes,
                barycentric_coordinates.data(),frozen_metric.data()).psi;

  BOOST_CHECK_EQUAL(pointwise_p1_after,pointwise_p1_before);
  BOOST_CHECK_GT(
      std::abs(pointwise_p2_after - pointwise_p2_before),
      1.e-6*std::max(1.0,std::abs(pointwise_p2_before)));
}

template<class MFT, int gdim>
void check_stepdistance_edge_derivative_dispatch()
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_quadrature_order = 5;
  parameters.objective_p = 1.25;
  parameters.step_distance_regularization = 1.e-6;
  parameters.step_distance_shape_volume = false;
  parameters.step_distance_cavity_target_average = false;
  parameters.step_distance_barrier_beta = 0.0;

  Mesh<MFT> mesh;
  initialize_curved_p2_stepdistance_element<MFT,gdim>(mesh,parameters);
  const int *nodes = mesh.ent2poi(gdim)[0];
  const SimplexQuadratureView<gdim> quadrature
      = get_objective_quadrature<gdim>(
            parameters.objective_quadrature_order);
  constexpr int active_edge_control_point = gdim + 1;

  double expected_gradient[gdim];
  const double expected_value
      = integrate_objective_quadrature_derivatives<
            MFT,gdim,gdim,2,QuaFun::StepDistance,double>(
                mesh,AsDeg::Pk,AsDeg::P1,nodes,quadrature,
                active_edge_control_point,FEBasis::Lagrange,
                expected_gradient,NULL);

  double dispatched_gradient[gdim];
  const double dispatched_value
      = d_metqua<MFT,gdim,gdim,QuaFun::StepDistance,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,
            active_edge_control_point,FEBasis::Lagrange,DifVar::None,
            dispatched_gradient,NULL,1.0);
  BOOST_CHECK_EQUAL(dispatched_value,expected_value);
  for(int component = 0; component < gdim; component++){
    BOOST_CHECK_EQUAL(dispatched_gradient[component],
                      expected_gradient[component]);
    BOOST_CHECK(std::isfinite(dispatched_gradient[component]));
  }

  const double value_only
      = metqua<MFT,gdim,gdim,QuaFun::StepDistance,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,1.0);
  BOOST_CHECK_CLOSE_FRACTION(dispatched_value,value_only,3.e-14);

  double gradient_norm_squared = 0.0;
  for(int component = 0; component < gdim; component++){
    gradient_norm_squared
        += dispatched_gradient[component]*dispatched_gradient[component];
  }
  BOOST_CHECK_GT(gradient_norm_squared,1.e-8);
}

} // namespace

BOOST_AUTO_TEST_CASE(p2_stepdistance_value_uses_degree_aware_geometry)
{
  check_stepdistance_value_geometry_dispatch<MetricFieldFE,2>();
  check_stepdistance_value_geometry_dispatch<MetricFieldAnalytical,2>();
  check_stepdistance_value_geometry_dispatch<MetricFieldFE,3>();
  check_stepdistance_value_geometry_dispatch<MetricFieldAnalytical,3>();
}

BOOST_AUTO_TEST_CASE(p2_edge_stepdistance_derivative_uses_degree_aware_dispatch)
{
  check_stepdistance_edge_derivative_dispatch<MetricFieldFE,2>();
  check_stepdistance_edge_derivative_dispatch<MetricFieldAnalytical,2>();
  check_stepdistance_edge_derivative_dispatch<MetricFieldFE,3>();
  check_stepdistance_edge_derivative_dispatch<MetricFieldAnalytical,3>();
}

} // namespace Metris
