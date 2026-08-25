// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_high_order_sizeshape_value

#include <boost/test/included/unit_test.hpp>

#include "frozen_objective_value_ad.hxx"
#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "ho_constants.hxx"
#include "linalg/det.hxx"
#include "linalg/explogmet.hxx"
#include "linalg/symidx.hxx"
#include "low_eval.hxx"
#include "quality/low_metqua.hxx"
#include "quality/low_metqua_d.hxx"
#include "quality/objective_quadrature_derivatives.hxx"
#include "quality/objective_quadrature_sample.hxx"
#include "quality/objective_quadrature_value.hxx"
#include "quality/quafun_sizeshape.hxx"
#include "quality/simplex_quadrature.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

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
void independent_p2_geometry(const Mesh<MFT> &mesh,
                             const int *nodes,
                             const double *barycentric_coordinates,
                             double *physical_coordinates,
                             double *jacobian_transpose)
{
  constexpr int node_count = getnnode(gdim,2);
  constexpr auto ordering = ORDELT(gdim);
  std::fill(physical_coordinates,physical_coordinates + gdim,0.0);
  std::fill(jacobian_transpose,
            jacobian_transpose + gdim*gdim,0.0);

  for(int inode = 0; inode < node_count; inode++){
    int first_barycentric_index = -1;
    int second_barycentric_index = -1;
    double basis_value = 0.0;
    double basis_gradient[gdim] = {};

    for(int ibary = 0; ibary < gdim + 1; ibary++){
      if(ordering[2][inode][ibary] == 2){
        first_barycentric_index = ibary;
        second_barycentric_index = ibary;
      }else if(ordering[2][inode][ibary] == 1){
        if(first_barycentric_index < 0){
          first_barycentric_index = ibary;
        }else{
          second_barycentric_index = ibary;
        }
      }
    }
    BOOST_REQUIRE(first_barycentric_index >= 0);
    BOOST_REQUIRE(second_barycentric_index >= 0);

    const auto barycentric_gradient = [](int ibary, int derivative){
      if(ibary == 0) return -1.0;
      return ibary == derivative + 1 ? 1.0 : 0.0;
    };
    if(first_barycentric_index == second_barycentric_index){
      const double lambda
          = barycentric_coordinates[first_barycentric_index];
      basis_value = lambda*(2.0*lambda - 1.0);
      for(int derivative = 0; derivative < gdim; derivative++){
        basis_gradient[derivative]
            = (4.0*lambda - 1.0)
            * barycentric_gradient(first_barycentric_index,derivative);
      }
    }else{
      const double first_lambda
          = barycentric_coordinates[first_barycentric_index];
      const double second_lambda
          = barycentric_coordinates[second_barycentric_index];
      basis_value = 4.0*first_lambda*second_lambda;
      for(int derivative = 0; derivative < gdim; derivative++){
        basis_gradient[derivative]
            = 4.0*(second_lambda
                    *barycentric_gradient(first_barycentric_index,derivative)
                  + first_lambda
                    *barycentric_gradient(second_barycentric_index,derivative));
      }
    }

    for(int component = 0; component < gdim; component++){
      physical_coordinates[component]
          += basis_value*mesh.coord(nodes[inode],component);
      for(int derivative = 0; derivative < gdim; derivative++){
        jacobian_transpose[derivative*gdim + component]
            += basis_gradient[derivative]
             * mesh.coord(nodes[inode],component);
      }
    }
  }
}

template<int gdim>
double independent_metric_determinant(const double *metric)
{
  if constexpr(gdim == 2){
    return metric[0]*metric[2] - metric[1]*metric[1];
  }else{
    return metric[0]*(metric[2]*metric[5] - metric[4]*metric[4])
         - metric[1]*(metric[1]*metric[5] - metric[3]*metric[4])
         + metric[3]*(metric[1]*metric[4] - metric[2]*metric[3]);
  }
}

template<class MFT, int gdim>
double independent_p2_sizeshape_dense_reference(Mesh<MFT> &mesh)
{
  constexpr int metric_count = gdim*(gdim + 1)/2;
  constexpr int ngauss = 8;
  constexpr double gauss_points[ngauss] = {
      0.019855071751231884,0.10166676129318663,
      0.23723379504183550,0.40828267875217510,
      0.59171732124782490,0.76276620495816450,
      0.89833323870681337,0.98014492824876812};
  constexpr double gauss_weights[ngauss] = {
      0.05061426814518813,0.11119051722668724,
      0.15685332293894364,0.18134189168918099,
      0.18134189168918099,0.15685332293894364,
      0.11119051722668724,0.05061426814518813};
  constexpr double reference_measure = gdim == 2 ? 0.5 : 1.0/6.0;
  constexpr double quality_denominator = gdim == 2 ? 8.0 : 54.0;
  const int *nodes = mesh.ent2poi(gdim)[0];
  double integral = 0.0;

  const auto accumulate_sample = [&](const double *barycentric_coordinates,
                                     double normalized_weight){
    // The reference oracle owns its P2 basis, matrix algebra, pointwise
    // SizeShape formula, and quadrature loop. It intentionally does not call
    // eval2/eval3, sample preparation, or the production objective policy.
    double physical_coordinates[gdim];
    double jacobian_transpose[gdim*gdim];
    independent_p2_geometry<MFT,gdim>(
        mesh,nodes,barycentric_coordinates,
        physical_coordinates,jacobian_transpose);

    std::array<double,metric_count> metric{};
    if constexpr(std::is_same<MFT,MetricFieldAnalytical>::value){
      evaluate_test_metric<gdim>(physical_coordinates,1.0,metric.data());
    }else{
      metric = independent_p1_fe_metric<gdim>(
          mesh,nodes,barycentric_coordinates);
    }

    double regular_jacobian_transpose[gdim*gdim] = {};
    for(int row = 0; row < gdim; row++){
      for(int component = 0; component < gdim; component++){
        for(int derivative = 0; derivative < gdim; derivative++){
          regular_jacobian_transpose[row*gdim + component]
              += Constants::invtJ_0[hana::type_c<double>][gdim]
                    [row*gdim + derivative]
               * jacobian_transpose[derivative*gdim + component];
        }
      }
    }

    const auto metric_entry = [&](int row, int column){
      const int upper = std::max(row,column);
      const int lower = std::min(row,column);
      return metric[upper*(upper + 1)/2 + lower];
    };
    double trace = 0.0;
    for(int row = 0; row < gdim; row++){
      for(int first_component = 0;
          first_component < gdim; first_component++){
        for(int second_component = 0;
            second_component < gdim; second_component++){
          trace += regular_jacobian_transpose[row*gdim + first_component]
                 * metric_entry(first_component,second_component)
                 * regular_jacobian_transpose[row*gdim + second_component];
        }
      }
    }
    const double metric_determinant
        = independent_metric_determinant<gdim>(metric.data());
    const double regular_jacobian_determinant
        = determinant_from_rows<gdim>(regular_jacobian_transpose);
    const double transformed_determinant
        = regular_jacobian_determinant*regular_jacobian_determinant
        * metric_determinant;
    const double size_shape_quality
        = std::pow(trace,gdim)
        * (1.0 + 1.0/(transformed_determinant*transformed_determinant))
        / quality_denominator;
    double size_shape_error = size_shape_quality - 1.0;
    if(std::abs(size_shape_error)
       <= 32.0*std::numeric_limits<double>::epsilon()){
      size_shape_error = 0.0;
    }
    BOOST_REQUIRE_GE(size_shape_error,0.0);
    const double psi = mesh.param->objective_p == 1.0
                     ? size_shape_error
                     : std::pow(size_shape_error,mesh.param->objective_p);

    double theta = reference_measure
                 * std::abs(determinant_from_rows<gdim>(
                       jacobian_transpose));
    #ifdef INTQUALINRIEMSPACE
    theta *= std::sqrt(metric_determinant);
    #endif
    integral += normalized_weight*theta*psi;
  };

  for(int iu = 0; iu < ngauss; iu++){
    const double u = gauss_points[iu];
    for(int iv = 0; iv < ngauss; iv++){
      const double v = gauss_points[iv];
      if constexpr(gdim == 2){
        const double barycentric_coordinates[3]
            = {(1.0 - u)*(1.0 - v),u,(1.0 - u)*v};
        const double normalized_weight
            = 2.0*gauss_weights[iu]*gauss_weights[iv]*(1.0 - u);
        accumulate_sample(barycentric_coordinates,normalized_weight);
      }else{
        for(int iw = 0; iw < ngauss; iw++){
          const double w = gauss_points[iw];
          const double barycentric_coordinates[4]
              = {(1.0 - u)*(1.0 - v)*(1.0 - w),u,
                 (1.0 - u)*v,(1.0 - u)*(1.0 - v)*w};
          const double normalized_weight
              = 6.0*gauss_weights[iu]*gauss_weights[iv]*gauss_weights[iw]
              * (1.0 - u)*(1.0 - u)*(1.0 - v);
          accumulate_sample(barycentric_coordinates,normalized_weight);
        }
      }
    }
  }
  return integral;
}

template<class MFT, int gdim>
void check_p2_sizeshape_against_dense_reference()
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_p = 1.5;

  Mesh<MFT> mesh;
  initialize_curved_p2_element<MFT,gdim>(mesh,parameters);
  const double reference
      = independent_p2_sizeshape_dense_reference<MFT,gdim>(mesh);
  double values[4];
  double errors[4];
  for(int order = 2; order <= 5; order++){
    parameters.objective_quadrature_order = order;
    values[order - 2]
        = metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
              mesh,AsDeg::Pk,AsDeg::P1,0,1.0);
    errors[order - 2] = std::abs(values[order - 2] - reference);
    BOOST_CHECK(std::isfinite(values[order - 2]));
  }

  BOOST_TEST_MESSAGE(
      "P2 dense reference=" << reference
      << ", q2=" << values[0] << " (error " << errors[0] << ")"
      << ", q3=" << values[1] << " (error " << errors[1] << ")"
      << ", q4=" << values[2] << " (error " << errors[2] << ")"
      << ", q5=" << values[3] << " (error " << errors[3] << ")");
  BOOST_CHECK_LE(errors[3],errors[0]);
  BOOST_CHECK_LE(
      errors[3],2.e-5*std::max(1.0,std::abs(reference)));
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

template<class MFT, int gdim>
void check_sizeshape_edge_derivative_dispatch()
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_quadrature_order = 5;
  parameters.objective_p = 1.5;

  Mesh<MFT> mesh;
  initialize_curved_p2_element<MFT,gdim>(mesh,parameters);
  const int *nodes = mesh.ent2poi(gdim)[0];
  const SimplexQuadratureView<gdim> quadrature
      = get_objective_quadrature<gdim>(
            parameters.objective_quadrature_order);
  constexpr int active_edge_control_point = gdim + 1;

  double expected_gradient[gdim];
  const double expected_value
      = integrate_objective_quadrature_derivatives<
            MFT,gdim,gdim,2,QuaFun::SizeShape,double>(
                mesh,AsDeg::Pk,AsDeg::P1,nodes,quadrature,
                active_edge_control_point,FEBasis::Lagrange,
                expected_gradient,NULL);

  double dispatched_gradient[gdim];
  const double dispatched_value
      = d_metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
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
      = metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,1.0);
  BOOST_CHECK_CLOSE_FRACTION(dispatched_value,value_only,3.e-14);

  double gradient_norm_squared = 0.0;
  for(int component = 0; component < gdim; component++){
    gradient_norm_squared
        += dispatched_gradient[component]*dispatched_gradient[component];
  }
  BOOST_CHECK_GT(gradient_norm_squared,1.e-8);
}

template<int gdim>
struct FrozenP2SizeShapeSamples
{
  static constexpr int metric_count = gdim*(gdim + 1)/2;
  std::vector<std::array<double,gdim + 1>> barycentric_points;
  std::vector<std::array<double,gdim*gdim>>
      regular_jacobian_transposes;
  std::vector<std::array<double,metric_count>> metrics;
  std::vector<double> objective_weights;
};

template<class MFT, int gdim>
FrozenP2SizeShapeSamples<gdim>
capture_frozen_p2_sizeshape_samples(Mesh<MFT> &mesh)
{
  FrozenP2SizeShapeSamples<gdim> samples;
  const int *nodes = mesh.ent2poi(gdim)[0];
  const SimplexQuadratureView<gdim> quadrature
      = get_objective_quadrature<gdim>(
            mesh.param->objective_quadrature_order);
  const ObjectiveQuadratureTheta theta_mode
      = objective_quadrature_theta_mode<QuaFun::SizeShape>(mesh);
  samples.barycentric_points.reserve(quadrature.size());
  samples.regular_jacobian_transposes.reserve(quadrature.size());
  samples.metrics.reserve(quadrature.size());
  samples.objective_weights.reserve(quadrature.size());

  for(int iquad = 0; iquad < quadrature.size(); iquad++){
    const ObjectiveQuadratureSample<gdim,gdim,2> sample
        = prepare_objective_quadrature_sample<MFT,gdim,gdim,2>(
              mesh,AsDeg::P1,nodes,quadrature[iquad],theta_mode);
    BOOST_REQUIRE(sample.theta_is_valid);
    samples.barycentric_points.push_back(
        sample.barycentric_coordinates);
    samples.regular_jacobian_transposes.push_back(
        sample.regular_jacobian_transpose);
    samples.metrics.push_back(sample.metric);
    samples.objective_weights.push_back(
        sample.quadrature_weight*sample.theta);
  }
  return samples;
}

template<class MFT, int gdim>
double evaluate_frozen_p2_sizeshape_value(
    Mesh<MFT> &mesh,
    const FrozenP2SizeShapeSamples<gdim> &samples)
{
  const int *nodes = mesh.ent2poi(gdim)[0];
  double value = 0.0;
  for(std::size_t isample = 0;
      isample < samples.objective_weights.size(); isample++){
    value += samples.objective_weights[isample]
           * quafun_sizeshape<MFT,gdim,gdim,double>(
                 mesh,AsDeg::Pk,AsDeg::P1,nodes,
                 samples.barycentric_points[isample].data(),
                 samples.metrics[isample].data());
  }
  return value;
}

template<class MFT, int gdim>
std::array<double,gdim> evaluate_frozen_p2_sizeshape_gradient_by_ad(
    Mesh<MFT> &mesh,
    const FrozenP2SizeShapeSamples<gdim> &samples,
    const int active_control_point,
    double *primal_value)
{
  using S = SANS::SurrealS<gdim,double>;
  S value = S(0.0);
  for(std::size_t isample = 0;
      isample < samples.objective_weights.size(); isample++){
    const std::array<double,gdim> regular_basis_gradient
        = objective_regular_basis_gradient<gdim,2>(
              mesh.getBasis(),active_control_point,
              samples.barycentric_points[isample].data());
    value += S(samples.objective_weights[isample])
           * FrozenObjectiveValueAD::sizeshape_pointwise_value<gdim>(
                 samples.regular_jacobian_transposes[isample].data(),
                 samples.metrics[isample].data(),
                 regular_basis_gradient.data(),
                 mesh.param->objective_p);
  }

  std::array<double,gdim> gradient{};
  *primal_value = value.value();
  for(int component = 0; component < gdim; component++){
    gradient[component] = value.deriv(component);
  }
  return gradient;
}

template<class MFT, int gdim>
double evaluate_frozen_p2_sizeshape_derivatives(
    Mesh<MFT> &mesh,
    const FrozenP2SizeShapeSamples<gdim> &samples,
    const int active_control_point,
    double *gradient,
    double *hessian)
{
  constexpr int hessian_count = gdim*(gdim + 1)/2;
  for(int component = 0; component < gdim; component++){
    gradient[component] = 0.0;
  }
  if(hessian != NULL){
    for(int entry = 0; entry < hessian_count; entry++){
      hessian[entry] = 0.0;
    }
  }

  const int *nodes = mesh.ent2poi(gdim)[0];
  double value = 0.0;
  for(std::size_t isample = 0;
      isample < samples.objective_weights.size(); isample++){
    double pointwise_gradient[gdim];
    double pointwise_hessian[hessian_count];
    value += samples.objective_weights[isample]
           * d_quafun_sizeshape<MFT,gdim,gdim,double>(
                 mesh,AsDeg::Pk,AsDeg::P1,nodes,
                 samples.barycentric_points[isample].data(),
                 samples.metrics[isample].data(),active_control_point,
                 mesh.getBasis(),DifVar::None,pointwise_gradient,
                 hessian == NULL ? NULL : pointwise_hessian);
    for(int component = 0; component < gdim; component++){
      gradient[component] += samples.objective_weights[isample]
                           * pointwise_gradient[component];
    }
    if(hessian != NULL){
      for(int entry = 0; entry < hessian_count; entry++){
        hessian[entry] += samples.objective_weights[isample]
                        * pointwise_hessian[entry];
      }
    }
  }
  return value;
}

template<int gdim>
void evaluate_frozen_p2_sizeshape_hessian_by_gradient_ad(
    const FrozenP2SizeShapeSamples<gdim> &samples,
    const int active_control_point,
    FEBasis basis,
    double objective_power,
    double *gradient,
    double *hessian)
{
  using S = SANS::SurrealS<gdim,double>;
  S differentiated_gradient[gdim];
  for(int component = 0; component < gdim; component++){
    differentiated_gradient[component] = S(0.0);
  }

  for(std::size_t isample = 0;
      isample < samples.objective_weights.size(); isample++){
    const std::array<double,gdim> regular_basis_gradient
        = objective_regular_basis_gradient<gdim,2>(
              basis,active_control_point,
              samples.barycentric_points[isample].data());
    S pointwise_gradient[gdim];
    FrozenObjectiveValueAD::sizeshape_pointwise_gradient<gdim>(
        samples.regular_jacobian_transposes[isample].data(),
        samples.metrics[isample].data(),regular_basis_gradient.data(),
        objective_power,pointwise_gradient);
    for(int component = 0; component < gdim; component++){
      differentiated_gradient[component]
          += S(samples.objective_weights[isample])
           * pointwise_gradient[component];
    }
  }

  for(int component = 0; component < gdim; component++){
    gradient[component] = differentiated_gradient[component].value();
    for(int derivative = component;
        derivative < gdim; derivative++){
      hessian[sym2idx(component,derivative)]
          = differentiated_gradient[component].deriv(derivative);
    }
  }
}

template<class MFT, int gdim>
void check_sizeshape_edge_gradient(FEBasis geometry_basis)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_quadrature_order = 5;
  parameters.objective_p = 1.5;

  Mesh<MFT> mesh;
  initialize_curved_p2_element<MFT,gdim>(mesh,parameters);
  mesh.setBasis(geometry_basis);
  constexpr int active_edge_control_point = gdim + 1;
  const int active_point
      = mesh.ent2poi(gdim)(0,active_edge_control_point);
  const FrozenP2SizeShapeSamples<gdim> samples
      = capture_frozen_p2_sizeshape_samples<MFT,gdim>(mesh);

  double gradient[gdim];
  const double differentiated_value
      = d_metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,
            active_edge_control_point,geometry_basis,DifVar::None,
            gradient,NULL,1.0);
  const double frozen_value
      = evaluate_frozen_p2_sizeshape_value<MFT,gdim>(mesh,samples);
  BOOST_CHECK_CLOSE_FRACTION(
      differentiated_value,frozen_value,3.e-14);

  double automatic_value;
  const std::array<double,gdim> automatic_gradient
      = evaluate_frozen_p2_sizeshape_gradient_by_ad<MFT,gdim>(
            mesh,samples,active_edge_control_point,&automatic_value);
  BOOST_CHECK_CLOSE_FRACTION(
      automatic_value,frozen_value,3.e-14);
  for(int component = 0; component < gdim; component++){
    BOOST_CHECK_SMALL(
        gradient[component] - automatic_gradient[component],
        5.e-14*(1.0 + std::abs(automatic_gradient[component])));
  }

  constexpr double finite_difference_step = 2.e-6;
  for(int component = 0; component < gdim; component++){
    const double coordinate = mesh.coord(active_point,component);
    mesh.coord(active_point,component)
        = coordinate + finite_difference_step;
    const double plus_value
        = evaluate_frozen_p2_sizeshape_value<MFT,gdim>(mesh,samples);
    mesh.coord(active_point,component)
        = coordinate - finite_difference_step;
    const double minus_value
        = evaluate_frozen_p2_sizeshape_value<MFT,gdim>(mesh,samples);
    mesh.coord(active_point,component) = coordinate;

    const double finite_difference_gradient
        = (plus_value - minus_value)/(2.0*finite_difference_step);
    const double gradient_error
        = std::abs(gradient[component] - finite_difference_gradient);
    BOOST_TEST_MESSAGE(
        "P2 frozen edge gradient component " << component
        << ": analytical=" << gradient[component]
        << ", finite difference=" << finite_difference_gradient
        << ", error=" << gradient_error);
    BOOST_CHECK_SMALL(
        gradient_error,
        1.e-7*(1.0 + std::abs(finite_difference_gradient)));
  }
}

template<class MFT, int gdim>
void check_sizeshape_edge_hessian(FEBasis geometry_basis)
{
  constexpr int hessian_count = gdim*(gdim + 1)/2;
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_quadrature_order = 5;
  parameters.objective_p = 1.5;

  Mesh<MFT> mesh;
  initialize_curved_p2_element<MFT,gdim>(mesh,parameters);
  mesh.setBasis(geometry_basis);
  constexpr int active_edge_control_point = gdim + 1;
  const int active_point
      = mesh.ent2poi(gdim)(0,active_edge_control_point);
  const FrozenP2SizeShapeSamples<gdim> samples
      = capture_frozen_p2_sizeshape_samples<MFT,gdim>(mesh);

  double production_gradient[gdim];
  double production_hessian[hessian_count];
  const double production_value
      = d_metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,
            active_edge_control_point,geometry_basis,DifVar::None,
            production_gradient,production_hessian,1.0);

  double reconstructed_gradient[gdim];
  double reconstructed_hessian[hessian_count];
  const double reconstructed_value
      = evaluate_frozen_p2_sizeshape_derivatives<MFT,gdim>(
            mesh,samples,active_edge_control_point,
            reconstructed_gradient,reconstructed_hessian);
  BOOST_CHECK_CLOSE_FRACTION(
      production_value,reconstructed_value,3.e-14);
  for(int component = 0; component < gdim; component++){
    BOOST_CHECK_SMALL(
        production_gradient[component] - reconstructed_gradient[component],
        5.e-14*(1.0 + std::abs(reconstructed_gradient[component])));
  }
  for(int entry = 0; entry < hessian_count; entry++){
    BOOST_CHECK_SMALL(
        production_hessian[entry] - reconstructed_hessian[entry],
        5.e-14*(1.0 + std::abs(reconstructed_hessian[entry])));
  }

  double automatic_gradient[gdim];
  double automatic_hessian[hessian_count];
  evaluate_frozen_p2_sizeshape_hessian_by_gradient_ad<gdim>(
      samples,active_edge_control_point,mesh.getBasis(),
      parameters.objective_p,automatic_gradient,automatic_hessian);

  double maximum_automatic_hessian_error = 0.0;
  for(int component = 0; component < gdim; component++){
    BOOST_CHECK_SMALL(
        production_gradient[component] - automatic_gradient[component],
        5.e-14*(1.0 + std::abs(automatic_gradient[component])));
  }
  for(int entry = 0; entry < hessian_count; entry++){
    const double error
        = std::abs(production_hessian[entry] - automatic_hessian[entry]);
    maximum_automatic_hessian_error
        = std::max(maximum_automatic_hessian_error,error);
    BOOST_CHECK_SMALL(
        error,5.e-14*(1.0 + std::abs(automatic_hessian[entry])));
  }

  constexpr double finite_difference_step = 2.e-6;
  double maximum_finite_difference_hessian_error = 0.0;
  for(int derivative = 0; derivative < gdim; derivative++){
    const double coordinate = mesh.coord(active_point,derivative);
    double plus_gradient[gdim];
    double minus_gradient[gdim];

    mesh.coord(active_point,derivative)
        = coordinate + finite_difference_step;
    evaluate_frozen_p2_sizeshape_derivatives<MFT,gdim>(
        mesh,samples,active_edge_control_point,plus_gradient,NULL);
    mesh.coord(active_point,derivative)
        = coordinate - finite_difference_step;
    evaluate_frozen_p2_sizeshape_derivatives<MFT,gdim>(
        mesh,samples,active_edge_control_point,minus_gradient,NULL);
    mesh.coord(active_point,derivative) = coordinate;

    for(int component = 0; component <= derivative; component++){
      const int entry = sym2idx(component,derivative);
      const double finite_difference_hessian
          = (plus_gradient[component] - minus_gradient[component])
           /(2.0*finite_difference_step);
      const double error
          = std::abs(production_hessian[entry]
                     - finite_difference_hessian);
      maximum_finite_difference_hessian_error
          = std::max(maximum_finite_difference_hessian_error,error);
      BOOST_CHECK_SMALL(
          error,1.e-7*(1.0 + std::abs(finite_difference_hessian)));
    }
  }

  BOOST_TEST_MESSAGE(
      "P2 frozen edge Hessian: AD error="
      << maximum_automatic_hessian_error
      << ", finite-difference error="
      << maximum_finite_difference_hessian_error);
}

template<class MFT, int gdim>
void check_sizeshape_lagrange_bezier_equivalence()
{
  constexpr int hessian_count = gdim*(gdim + 1)/2;
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_quadrature_order = 5;
  parameters.objective_p = 1.5;

  Mesh<MFT> mesh;
  initialize_curved_p2_element<MFT,gdim>(mesh,parameters);
  constexpr int active_edge_control_point = gdim + 1;

  const FrozenP2SizeShapeSamples<gdim> lagrange_samples
      = capture_frozen_p2_sizeshape_samples<MFT,gdim>(mesh);
  const double lagrange_value
      = metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,1.0);
  double lagrange_gradient[gdim];
  double lagrange_hessian[hessian_count];
  const double lagrange_differentiated_value
      = d_metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,active_edge_control_point,
            FEBasis::Lagrange,DifVar::None,
            lagrange_gradient,lagrange_hessian,1.0);

  mesh.setBasis(FEBasis::Bezier);
  const FrozenP2SizeShapeSamples<gdim> bezier_samples
      = capture_frozen_p2_sizeshape_samples<MFT,gdim>(mesh);
  const double bezier_value
      = metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,1.0);
  double bezier_gradient[gdim];
  double bezier_hessian[hessian_count];
  const double bezier_differentiated_value
      = d_metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,active_edge_control_point,
            FEBasis::Bezier,DifVar::None,
            bezier_gradient,bezier_hessian,1.0);

  double maximum_value_error = std::max(
      std::abs(lagrange_value - bezier_value),
      std::abs(lagrange_differentiated_value
               - bezier_differentiated_value));
  double maximum_sample_error = 0.0;
  double maximum_gradient_transform_error = 0.0;
  double maximum_hessian_transform_error = 0.0;
  BOOST_CHECK_CLOSE_FRACTION(lagrange_value,bezier_value,5.e-14);
  BOOST_CHECK_CLOSE_FRACTION(
      lagrange_differentiated_value,bezier_differentiated_value,5.e-14);
  BOOST_REQUIRE_EQUAL(
      lagrange_samples.objective_weights.size(),
      bezier_samples.objective_weights.size());
  for(std::size_t isample = 0;
      isample < lagrange_samples.objective_weights.size(); isample++){
    for(int entry = 0; entry < gdim*gdim; entry++){
      maximum_sample_error = std::max(
          maximum_sample_error,
          std::abs(
              lagrange_samples.regular_jacobian_transposes[isample][entry]
              - bezier_samples.regular_jacobian_transposes[isample][entry]));
      BOOST_CHECK_SMALL(
          lagrange_samples.regular_jacobian_transposes[isample][entry]
              - bezier_samples.regular_jacobian_transposes[isample][entry],
          5.e-14*(1.0 + std::abs(
              lagrange_samples.regular_jacobian_transposes[isample][entry])));
    }
    for(int entry = 0;
        entry < FrozenP2SizeShapeSamples<gdim>::metric_count; entry++){
      maximum_sample_error = std::max(
          maximum_sample_error,
          std::abs(lagrange_samples.metrics[isample][entry]
                   - bezier_samples.metrics[isample][entry]));
      BOOST_CHECK_SMALL(
          lagrange_samples.metrics[isample][entry]
              - bezier_samples.metrics[isample][entry],
          5.e-14*(1.0 + std::abs(
              lagrange_samples.metrics[isample][entry])));
    }
    maximum_sample_error = std::max(
        maximum_sample_error,
        std::abs(lagrange_samples.objective_weights[isample]
                 - bezier_samples.objective_weights[isample]));
    BOOST_CHECK_SMALL(
        lagrange_samples.objective_weights[isample]
            - bezier_samples.objective_weights[isample],
        5.e-14*(1.0 + std::abs(
            lagrange_samples.objective_weights[isample])));
  }

  for(int component = 0; component < gdim; component++){
    maximum_gradient_transform_error = std::max(
        maximum_gradient_transform_error,
        std::abs(lagrange_gradient[component]
                 - 2.0*bezier_gradient[component]));
    BOOST_CHECK_SMALL(
        lagrange_gradient[component] - 2.0*bezier_gradient[component],
        5.e-14*(1.0 + std::abs(lagrange_gradient[component])));
  }
  for(int entry = 0; entry < hessian_count; entry++){
    maximum_hessian_transform_error = std::max(
        maximum_hessian_transform_error,
        std::abs(lagrange_hessian[entry]
                 - 4.0*bezier_hessian[entry]));
    BOOST_CHECK_SMALL(
        lagrange_hessian[entry] - 4.0*bezier_hessian[entry],
        5.e-14*(1.0 + std::abs(lagrange_hessian[entry])));
  }
  BOOST_TEST_MESSAGE(
      "P2 Lagrange-Bezier equivalence: value error="
      << maximum_value_error
      << ", sample error=" << maximum_sample_error
      << ", gradient-transform error="
      << maximum_gradient_transform_error
      << ", Hessian-transform error="
      << maximum_hessian_transform_error);
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

BOOST_AUTO_TEST_CASE(p2_sizeshape_matches_independent_dense_reference)
{
  check_p2_sizeshape_against_dense_reference<MetricFieldFE,2>();
  check_p2_sizeshape_against_dense_reference<MetricFieldAnalytical,2>();
  check_p2_sizeshape_against_dense_reference<MetricFieldFE,3>();
  check_p2_sizeshape_against_dense_reference<MetricFieldAnalytical,3>();
}

BOOST_AUTO_TEST_CASE(p2_edge_sizeshape_derivative_uses_degree_aware_dispatch)
{
  check_sizeshape_edge_derivative_dispatch<MetricFieldFE,2>();
  check_sizeshape_edge_derivative_dispatch<MetricFieldAnalytical,2>();
  check_sizeshape_edge_derivative_dispatch<MetricFieldFE,3>();
  check_sizeshape_edge_derivative_dispatch<MetricFieldAnalytical,3>();
}

BOOST_AUTO_TEST_CASE(p2_edge_sizeshape_gradient_matches_frozen_value_oracles)
{
  for(FEBasis basis : {FEBasis::Lagrange,FEBasis::Bezier}){
    check_sizeshape_edge_gradient<MetricFieldFE,2>(basis);
    check_sizeshape_edge_gradient<MetricFieldAnalytical,2>(basis);
    check_sizeshape_edge_gradient<MetricFieldFE,3>(basis);
    check_sizeshape_edge_gradient<MetricFieldAnalytical,3>(basis);
  }
}

BOOST_AUTO_TEST_CASE(p2_edge_sizeshape_hessian_matches_frozen_gradient_oracles)
{
  for(FEBasis basis : {FEBasis::Lagrange,FEBasis::Bezier}){
    check_sizeshape_edge_hessian<MetricFieldFE,2>(basis);
    check_sizeshape_edge_hessian<MetricFieldAnalytical,2>(basis);
    check_sizeshape_edge_hessian<MetricFieldFE,3>(basis);
    check_sizeshape_edge_hessian<MetricFieldAnalytical,3>(basis);
  }
}

BOOST_AUTO_TEST_CASE(p2_sizeshape_lagrange_and_bezier_are_equivalent)
{
  check_sizeshape_lagrange_bezier_equivalence<MetricFieldFE,2>();
  check_sizeshape_lagrange_bezier_equivalence<MetricFieldAnalytical,2>();
  check_sizeshape_lagrange_bezier_equivalence<MetricFieldFE,3>();
  check_sizeshape_lagrange_bezier_equivalence<MetricFieldAnalytical,3>();
}

} // namespace Metris
