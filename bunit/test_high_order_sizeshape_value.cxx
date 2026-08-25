// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_high_order_sizeshape_value

#include <boost/test/included/unit_test.hpp>

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "ho_constants.hxx"
#include "linalg/det.hxx"
#include "linalg/explogmet.hxx"
#include "low_eval.hxx"
#include "quality/low_metqua.hxx"
#include "quality/objective_quadrature_sample.hxx"
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
void evaluate_test_metric(const double *coordinate,
                          double scale,
                          double *metric)
{
  if constexpr(gdim == 2){
    metric[0] = scale*(1.40 + 0.20*coordinate[0]);
    metric[1] = scale*(0.05 + 0.02*coordinate[1]);
    metric[2] = scale*(0.90 + 0.15*coordinate[1]);
  }else{
    metric[0] = scale*(1.45 + 0.12*coordinate[0]);
    metric[1] = scale*(0.04 + 0.01*coordinate[1]);
    metric[2] = scale*(1.10 + 0.10*coordinate[1]);
    metric[3] = scale*(-0.02 + 0.01*coordinate[2]);
    metric[4] = scale*(0.03 + 0.01*coordinate[0]);
    metric[5] = scale*(0.88 + 0.09*coordinate[2]);
  }
}

template<int gdim>
void varying_analytical_metric(
    const AnaMetCtx *, const double *coordinate, double scale,
    int derivative_order, double *metric, double *metric_derivative)
{
  constexpr int metric_count = gdim*(gdim + 1)/2;
  evaluate_test_metric<gdim>(coordinate,scale,metric);
  if(derivative_order == 0) return;
  std::fill(metric_derivative,
            metric_derivative + gdim*metric_count,0.0);
}

template<class MFT, int gdim>
void initialize_curved_p2_element(Mesh<MFT> &mesh,
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
          AnaMetFun(varying_analytical_metric<2>));
    }else{
      parameters.setAnalyticalMetric(
          AnaMetFun(varying_analytical_metric<3>));
    }
    mesh.met.setAnalyticalMetric(parameters);
  }

  for(int ipoin = 0; ipoin < node_count; ipoin++){
    double metric[metric_count];
    evaluate_test_metric<gdim>(mesh.coord[ipoin],1.0,metric);
    for(int imetric = 0; imetric < metric_count; imetric++){
      mesh.met(ipoin,imetric) = metric[imetric];
    }
  }
}

template<int gdim>
double determinant_from_rows(const double *matrix)
{
  if constexpr(gdim == 2){
    return matrix[0]*matrix[3] - matrix[1]*matrix[2];
  }else{
    return matrix[0]*(matrix[4]*matrix[8] - matrix[5]*matrix[7])
         - matrix[1]*(matrix[3]*matrix[8] - matrix[5]*matrix[6])
         + matrix[2]*(matrix[3]*matrix[7] - matrix[4]*matrix[6]);
  }
}

template<int gdim>
std::array<double,gdim*(gdim + 1)/2>
independent_p1_fe_metric(const Mesh<MetricFieldFE> &mesh,
                         const int *nodes,
                         const double *barycentric_coordinates)
{
  constexpr int metric_count = gdim*(gdim + 1)/2;
  std::array<double,metric_count> expected_metric{};
  for(int vertex = 0; vertex < gdim + 1; vertex++){
    std::array<double,metric_count> logarithmic_metric{};
    for(int imetric = 0; imetric < metric_count; imetric++){
      logarithmic_metric[imetric] = mesh.met(nodes[vertex],imetric);
    }
    getlogmet_inp<gdim,double>(logarithmic_metric.data());
    for(int imetric = 0; imetric < metric_count; imetric++){
      expected_metric[imetric]
          += barycentric_coordinates[vertex]*logarithmic_metric[imetric];
    }
  }
  getexpmet_inp<gdim,double>(expected_metric.data());
  return expected_metric;
}

template<class MFT, int gdim>
void check_p2_quadrature_samples()
{
  constexpr int metric_count = gdim*(gdim + 1)/2;
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_quadrature_order = 5;

  Mesh<MFT> mesh;
  initialize_curved_p2_element<MFT,gdim>(mesh,parameters);
  const int *nodes = mesh.ent2poi(gdim)[0];
  const SimplexQuadratureView<gdim> quadrature
      = get_objective_quadrature<gdim>(
            parameters.objective_quadrature_order);
  const ObjectiveQuadratureTheta theta_mode
      = objective_quadrature_theta_mode<QuaFun::SizeShape>(mesh);
  constexpr auto ordering = ORDELT(gdim);
  constexpr int node_count = getnnode(gdim,2);

  for(int iquad = 0; iquad < quadrature.size(); iquad++){
    const SimplexQuadraturePointView<gdim> quadrature_point
        = quadrature[iquad];
    const ObjectiveQuadratureSample<gdim,gdim,2> sample
        = prepare_objective_quadrature_sample<MFT,gdim,gdim,2>(
              mesh,AsDeg::P1,nodes,quadrature_point,theta_mode);

    double expected_coordinates[gdim];
    double expected_jacobian_transpose[gdim*gdim];
    std::fill(expected_coordinates,expected_coordinates + gdim,0.0);
    std::fill(expected_jacobian_transpose,
              expected_jacobian_transpose + gdim*gdim,0.0);
    for(int inode = 0; inode < node_count; inode++){
      double basis_gradient[gdim];
      const double basis_value = eval_lagrangefunc<2,gdim>(
          ordering[2][inode],quadrature_point.bary,1,basis_gradient);
      for(int component = 0; component < gdim; component++){
        expected_coordinates[component]
            += basis_value*mesh.coord(nodes[inode],component);
        for(int derivative = 0; derivative < gdim; derivative++){
          expected_jacobian_transpose[derivative*gdim + component]
              += basis_gradient[derivative]
               * mesh.coord(nodes[inode],component);
        }
      }
    }

    for(int component = 0; component < gdim; component++){
      BOOST_CHECK_CLOSE_FRACTION(
          sample.physical_coordinates[component],
          expected_coordinates[component],2.e-14);
    }
    for(int ientry = 0; ientry < gdim*gdim; ientry++){
      BOOST_CHECK_CLOSE_FRACTION(
          sample.canonical_jacobian_transpose[ientry],
          expected_jacobian_transpose[ientry],2.e-14);
    }

    std::array<double,metric_count> expected_metric{};
    if constexpr(std::is_same<MFT,MetricFieldAnalytical>::value){
      evaluate_test_metric<gdim>(expected_coordinates,1.0,
                                 expected_metric.data());
    }else{
      expected_metric = independent_p1_fe_metric<gdim>(
          mesh,nodes,quadrature_point.bary);
    }
    for(int imetric = 0; imetric < metric_count; imetric++){
      BOOST_CHECK_CLOSE_FRACTION(
          sample.metric[imetric],expected_metric[imetric],3.e-13);
    }

    constexpr double reference_measure = gdim == 2 ? 0.5 : 1.0/6.0;
    double expected_theta
        = reference_measure*std::abs(
              determinant_from_rows<gdim>(expected_jacobian_transpose));
    #ifdef INTQUALINRIEMSPACE
    expected_theta *= std::sqrt(detsym<gdim>(expected_metric.data()));
    #endif
    BOOST_REQUIRE(sample.theta_is_valid);
    BOOST_CHECK_CLOSE_FRACTION(sample.theta,expected_theta,4.e-13);
  }

  const SimplexQuadraturePointView<gdim> probe = quadrature[0];
  const ObjectiveQuadratureSample<gdim,gdim,2> original_sample
      = prepare_objective_quadrature_sample<MFT,gdim,gdim,2>(
            mesh,AsDeg::P1,nodes,probe,theta_mode);
  const int first_ignored_node
      = std::is_same<MFT,MetricFieldAnalytical>::value ? 0 : gdim + 1;
  for(int inode = first_ignored_node; inode < node_count; inode++){
    for(int imetric = 0; imetric < metric_count; imetric++){
      mesh.met(nodes[inode],imetric) = 0.0;
    }
    for(int component = 0; component < gdim; component++){
      mesh.met(nodes[inode],component*(component + 1)/2 + component)
          = 50.0 + inode;
    }
  }
  const ObjectiveQuadratureSample<gdim,gdim,2> changed_nodal_sample
      = prepare_objective_quadrature_sample<MFT,gdim,gdim,2>(
            mesh,AsDeg::P1,nodes,probe,theta_mode);
  for(int imetric = 0; imetric < metric_count; imetric++){
    BOOST_CHECK_CLOSE_FRACTION(
        changed_nodal_sample.metric[imetric],
        original_sample.metric[imetric],2.e-14);
  }

  if constexpr(std::is_same<MFT,MetricFieldFE>::value){
    mesh.met(nodes[0],0) *= 1.8;
    const ObjectiveQuadratureSample<gdim,gdim,2> changed_corner_sample
        = prepare_objective_quadrature_sample<MFT,gdim,gdim,2>(
              mesh,AsDeg::P1,nodes,probe,theta_mode);
    double metric_change = 0.0;
    for(int imetric = 0; imetric < metric_count; imetric++){
      metric_change = std::max(
          metric_change,
          std::abs(changed_corner_sample.metric[imetric]
                   - original_sample.metric[imetric]));
    }
    BOOST_CHECK_GT(metric_change,1.e-4);
  }
}

template<class MFT, int gdim>
void check_sizeshape_geometry_dispatch()
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_quadrature_order = 5;
  parameters.objective_p = 1.0;

  Mesh<MFT> mesh;
  initialize_curved_p2_element<MFT,gdim>(mesh,parameters);
  const int *nodes = mesh.ent2poi(gdim)[0];
  const SimplexQuadratureView<gdim> quadrature
      = get_objective_quadrature<gdim>(
            parameters.objective_quadrature_order);

  const double expected_p2
      = integrate_objective_quadrature_value<
            MFT,gdim,gdim,2,QuaFun::SizeShape,double>(
                mesh,AsDeg::Pk,AsDeg::P1,nodes,quadrature);
  const double dispatched_p2
      = metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,1.0);
  BOOST_CHECK_EQUAL(dispatched_p2,expected_p2);

  const double expected_p1
      = integrate_objective_quadrature_value<
            MFT,gdim,gdim,1,QuaFun::SizeShape,double>(
                mesh,AsDeg::P1,AsDeg::P1,nodes,quadrature);
  const double dispatched_p1
      = metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
            mesh,AsDeg::P1,AsDeg::P1,0,1.0);
  BOOST_CHECK_EQUAL(dispatched_p1,expected_p1);

  BOOST_CHECK_GT(
      std::abs(dispatched_p2 - dispatched_p1),
      1.e-5*std::max(1.0,std::abs(dispatched_p2)));
}

} // namespace

BOOST_AUTO_TEST_CASE(p2_sizeshape_value_uses_degree_aware_dispatch)
{
  check_sizeshape_geometry_dispatch<MetricFieldFE,2>();
  check_sizeshape_geometry_dispatch<MetricFieldAnalytical,2>();
  check_sizeshape_geometry_dispatch<MetricFieldFE,3>();
  check_sizeshape_geometry_dispatch<MetricFieldAnalytical,3>();
}

BOOST_AUTO_TEST_CASE(p2_shared_samples_honor_geometry_and_metric_contract)
{
  check_p2_quadrature_samples<MetricFieldFE,2>();
  check_p2_quadrature_samples<MetricFieldAnalytical,2>();
  check_p2_quadrature_samples<MetricFieldFE,3>();
  check_p2_quadrature_samples<MetricFieldAnalytical,3>();
}

} // namespace Metris
