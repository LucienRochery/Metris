// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_high_order_stepdistance_value

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
#include "quality/pointwise_objective.hxx"
#include "quality/objective_quadrature_value.hxx"
#include "quality/quafun_stepDistance.hxx"
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
  std::array<double,metric_count> metric{};
  for(int vertex = 0; vertex < gdim + 1; vertex++){
    std::array<double,metric_count> logarithmic_metric{};
    for(int imetric = 0; imetric < metric_count; imetric++){
      logarithmic_metric[imetric] = mesh.met(nodes[vertex],imetric);
    }
    getlogmet_inp<gdim,double>(logarithmic_metric.data());
    for(int imetric = 0; imetric < metric_count; imetric++){
      metric[imetric]
          += barycentric_coordinates[vertex]*logarithmic_metric[imetric];
    }
  }
  getexpmet_inp<gdim,double>(metric.data());
  return metric;
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
    double basis_value;
    double basis_gradient[gdim] = {};
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
double independent_p2_stepdistance_dense_reference(Mesh<MFT> &mesh)
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
  const int *nodes = mesh.ent2poi(gdim)[0];
  double integral = 0.0;

  const auto accumulate_sample = [&](const double *barycentric_coordinates,
                                     double normalized_weight){
    // This oracle owns the P2 basis, geometry, metric interpolation, matrix
    // assembly, StepDistance formula, and quadrature loop. It deliberately
    // bypasses the production objective integration and sample preparation.
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
    double gram_matrix[gdim*gdim] = {};
    for(int row = 0; row < gdim; row++){
      for(int column = 0; column < gdim; column++){
        for(int first_component = 0;
            first_component < gdim; first_component++){
          for(int second_component = 0;
              second_component < gdim; second_component++){
            gram_matrix[row*gdim + column]
                += regular_jacobian_transpose[row*gdim + first_component]
                 * metric_entry(first_component,second_component)
                 * regular_jacobian_transpose[
                       column*gdim + second_component];
          }
        }
      }
    }
    double packed_gram_matrix[metric_count];
    packed_gram_matrix[0] = gram_matrix[0];
    packed_gram_matrix[1] = gram_matrix[1];
    packed_gram_matrix[2] = gram_matrix[gdim + 1];
    if constexpr(gdim == 3){
      packed_gram_matrix[3] = gram_matrix[2];
      packed_gram_matrix[4] = gram_matrix[5];
      packed_gram_matrix[5] = gram_matrix[8];
    }
    double eigenvalues[gdim];
    double eigenvectors[gdim*gdim];
    geteigsym<gdim,double>(
        packed_gram_matrix,eigenvalues,eigenvectors);
    double determinant = 1.0;
    double mean_logarithm = 0.0;
    double logarithms[gdim];
    for(int eigenvalue = 0; eigenvalue < gdim; eigenvalue++){
      BOOST_REQUIRE_GT(eigenvalues[eigenvalue],0.0);
      determinant *= eigenvalues[eigenvalue];
      logarithms[eigenvalue] = std::log(eigenvalues[eigenvalue]);
      mean_logarithm += logarithms[eigenvalue]/gdim;
    }
    double squared_distance = 0.0;
    if(mesh.param->step_distance_shape_volume){
      for(int eigenvalue = 0; eigenvalue < gdim; eigenvalue++){
        const double centered_logarithm
            = logarithms[eigenvalue] - mean_logarithm;
        squared_distance += centered_logarithm*centered_logarithm;
      }
      const double volume_coordinate
          = determinant - 1.0/determinant;
      squared_distance += volume_coordinate*volume_coordinate/(4.0*gdim);
    }else{
      for(int eigenvalue = 0; eigenvalue < gdim; eigenvalue++){
        squared_distance += logarithms[eigenvalue]*logarithms[eigenvalue];
      }
    }
    const double regularization
        = mesh.param->step_distance_regularization;
    const double psi
        = std::pow(squared_distance + regularization*regularization,
                   mesh.param->objective_p/2.0)
        - std::pow(regularization,mesh.param->objective_p);

    const double metric_determinant
        = independent_metric_determinant<gdim>(metric.data());
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
void check_p2_stepdistance_against_dense_reference(
    bool shape_volume = false)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_p = 1.25;
  parameters.step_distance_regularization = 1.e-6;
  parameters.step_distance_shape_volume = shape_volume;
  parameters.step_distance_cavity_target_average = false;
  parameters.step_distance_barrier_rho0 = 10.0;
  parameters.step_distance_barrier_beta = shape_volume ? 1.e6 : 0.0;

  Mesh<MFT> mesh;
  initialize_curved_p2_stepdistance_element<MFT,gdim>(mesh,parameters);
  const double reference
      = independent_p2_stepdistance_dense_reference<MFT,gdim>(mesh);
  double values[4];
  double errors[4];
  for(int order = 2; order <= 5; order++){
    parameters.objective_quadrature_order = order;
    values[order - 2]
        = metqua<MFT,gdim,gdim,QuaFun::StepDistance,double>(
              mesh,AsDeg::Pk,AsDeg::P1,0,1.0);
    errors[order - 2] = std::abs(values[order - 2] - reference);
    BOOST_CHECK(std::isfinite(values[order - 2]));
  }

  BOOST_TEST_MESSAGE(
      "P2 StepDistance "
      << (shape_volume ? "ShapeVolume " : "basic ")
      << "dense reference=" << reference
      << ", q2=" << values[0] << " (error " << errors[0] << ")"
      << ", q3=" << values[1] << " (error " << errors[1] << ")"
      << ", q4=" << values[2] << " (error " << errors[2] << ")"
      << ", q5=" << values[3] << " (error " << errors[3] << ")");
  BOOST_CHECK_LE(errors[3],errors[0]);
  BOOST_CHECK_LE(
      errors[3],2.e-5*std::max(1.0,std::abs(reference)));
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

template<int gdim>
struct FrozenP2StepDistanceSamples
{
  static constexpr int metric_count = gdim*(gdim + 1)/2;
  std::vector<std::array<double,gdim + 1>> barycentric_points;
  std::vector<std::array<double,gdim*gdim>>
      regular_jacobian_transposes;
  std::vector<std::array<double,metric_count>> metrics;
  std::vector<double> quadrature_weights;
  std::vector<double> objective_weights;
};

template<class MFT, int gdim>
FrozenP2StepDistanceSamples<gdim>
capture_frozen_p2_stepdistance_samples(Mesh<MFT> &mesh)
{
  BOOST_REQUIRE(!mesh.param->step_distance_cavity_target_average);

  FrozenP2StepDistanceSamples<gdim> samples;
  const int *nodes = mesh.ent2poi(gdim)[0];
  const SimplexQuadratureView<gdim> quadrature
      = get_objective_quadrature<gdim>(
            mesh.param->objective_quadrature_order);
  const ObjectiveQuadratureTheta theta_mode
      = objective_quadrature_theta_mode<QuaFun::StepDistance>(mesh);
  samples.barycentric_points.reserve(quadrature.size());
  samples.regular_jacobian_transposes.reserve(quadrature.size());
  samples.metrics.reserve(quadrature.size());
  samples.quadrature_weights.reserve(quadrature.size());
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
    samples.quadrature_weights.push_back(sample.quadrature_weight);
    samples.objective_weights.push_back(
        sample.quadrature_weight*sample.theta);
  }
  return samples;
}

template<int gdim>
std::array<double,gdim> independent_p2_regular_basis_gradient(
    FEBasis basis,
    int local_control_point,
    const double *barycentric_coordinates)
{
  constexpr auto ordering = ORDELT(gdim);
  const int *multi_index = ordering[2][local_control_point];
  int first_barycentric_index = -1;
  int second_barycentric_index = -1;
  for(int ibary = 0; ibary < gdim + 1; ibary++){
    if(multi_index[ibary] == 2){
      first_barycentric_index = ibary;
      second_barycentric_index = ibary;
    }else if(multi_index[ibary] == 1){
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
  double canonical_gradient[gdim] = {};
  if(first_barycentric_index == second_barycentric_index){
    const double lambda
        = barycentric_coordinates[first_barycentric_index];
    const double scale = basis == FEBasis::Lagrange
                       ? 4.0*lambda - 1.0
                       : 2.0*lambda;
    for(int derivative = 0; derivative < gdim; derivative++){
      canonical_gradient[derivative]
          = scale*barycentric_gradient(
                first_barycentric_index,derivative);
    }
  }else{
    const double first_lambda
        = barycentric_coordinates[first_barycentric_index];
    const double second_lambda
        = barycentric_coordinates[second_barycentric_index];
    const double scale = basis == FEBasis::Lagrange ? 4.0 : 2.0;
    for(int derivative = 0; derivative < gdim; derivative++){
      canonical_gradient[derivative]
          = scale*(second_lambda
                    *barycentric_gradient(
                        first_barycentric_index,derivative)
                  + first_lambda
                    *barycentric_gradient(
                        second_barycentric_index,derivative));
    }
  }

  std::array<double,gdim> regular_gradient{};
  for(int regular_component = 0;
      regular_component < gdim; regular_component++){
    for(int canonical_component = 0;
        canonical_component < gdim; canonical_component++){
      regular_gradient[regular_component]
          += canonical_gradient[canonical_component]
           * Constants::invtJ_0[hana::type_c<double>][gdim]
                              [regular_component*gdim
                               + canonical_component];
    }
  }
  return regular_gradient;
}

template<int gdim>
std::array<double,gdim> evaluate_frozen_p2_collapse_barrier_by_ad(
    const FrozenP2StepDistanceSamples<gdim> &samples,
    int active_control_point,
    FEBasis basis,
    double rho0,
    double beta,
    const double *displacement,
    double *value)
{
  using S = SANS::SurrealS<gdim,double>;
  S integrated_barrier = S(0.0);
  for(std::size_t isample = 0;
      isample < samples.quadrature_weights.size(); isample++){
    const std::array<double,gdim> regular_basis_gradient
        = independent_p2_regular_basis_gradient<gdim>(
              basis,active_control_point,
              samples.barycentric_points[isample].data());
    std::array<double,gdim*gdim> current_jacobian
        = samples.regular_jacobian_transposes[isample];
    if(displacement != NULL){
      for(int regular_component = 0;
          regular_component < gdim; regular_component++){
        for(int physical_component = 0;
            physical_component < gdim; physical_component++){
          current_jacobian[regular_component*gdim + physical_component]
              += regular_basis_gradient[regular_component]
               * displacement[physical_component];
        }
      }
    }
    integrated_barrier
        += S(samples.quadrature_weights[isample])
         * FrozenObjectiveValueAD::metric_volume_barrier_value<gdim>(
               current_jacobian.data(),samples.metrics[isample].data(),
               regular_basis_gradient.data(),rho0,beta);
  }

  std::array<double,gdim> gradient{};
  *value = integrated_barrier.value();
  for(int component = 0; component < gdim; component++){
    gradient[component] = integrated_barrier.deriv(component);
  }
  return gradient;
}

template<int gdim>
void evaluate_frozen_p2_collapse_barrier_hessian_by_gradient_ad(
    const FrozenP2StepDistanceSamples<gdim> &samples,
    int active_control_point,
    FEBasis basis,
    double rho0,
    double beta,
    double *gradient,
    double *hessian)
{
  using S = SANS::SurrealS<gdim,double>;
  S differentiated_gradient[gdim];
  for(int component = 0; component < gdim; component++){
    differentiated_gradient[component] = S(0.0);
  }

  for(std::size_t isample = 0;
      isample < samples.quadrature_weights.size(); isample++){
    const std::array<double,gdim> regular_basis_gradient
        = independent_p2_regular_basis_gradient<gdim>(
              basis,active_control_point,
              samples.barycentric_points[isample].data());
    S pointwise_gradient[gdim];
    FrozenObjectiveValueAD::metric_volume_barrier_gradient<gdim>(
        samples.regular_jacobian_transposes[isample].data(),
        samples.metrics[isample].data(),regular_basis_gradient.data(),
        rho0,beta,pointwise_gradient);
    for(int component = 0; component < gdim; component++){
      differentiated_gradient[component]
          += S(samples.quadrature_weights[isample])
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
double evaluate_frozen_p2_stepdistance_value(
    Mesh<MFT> &mesh,
    const FrozenP2StepDistanceSamples<gdim> &samples)
{
  const int *nodes = mesh.ent2poi(gdim)[0];
  double value = 0.0;
  for(std::size_t isample = 0;
      isample < samples.objective_weights.size(); isample++){
    value += samples.objective_weights[isample]
           * quafun_stepDistance<MFT,gdim,gdim,double>(
                 mesh,AsDeg::Pk,AsDeg::P1,nodes,
                 samples.barycentric_points[isample].data(),
                 samples.metrics[isample].data());
  }
  return value;
}

template<class MFT, int gdim>
std::array<double,gdim> evaluate_frozen_p2_stepdistance_gradient_by_ad(
    Mesh<MFT> &mesh,
    const FrozenP2StepDistanceSamples<gdim> &samples,
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
           * FrozenObjectiveValueAD::step_distance_pointwise_value<gdim>(
                 samples.regular_jacobian_transposes[isample].data(),
                 samples.metrics[isample].data(),
                 regular_basis_gradient.data(),
                 mesh.param->objective_p,
                 mesh.param->step_distance_regularization,
                 mesh.param->step_distance_shape_volume);
  }

  std::array<double,gdim> gradient{};
  *primal_value = value.value();
  for(int component = 0; component < gdim; component++){
    gradient[component] = value.deriv(component);
  }
  return gradient;
}

template<class MFT, int gdim>
double evaluate_frozen_p2_stepdistance_derivatives(
    Mesh<MFT> &mesh,
    const FrozenP2StepDistanceSamples<gdim> &samples,
    const int active_control_point,
    double *gradient,
    double *hessian)
{
  constexpr int hessian_count = gdim*(gdim + 1)/2;
  std::fill(gradient,gradient + gdim,0.0);
  if(hessian != NULL){
    std::fill(hessian,hessian + hessian_count,0.0);
  }

  const int *nodes = mesh.ent2poi(gdim)[0];
  double value = 0.0;
  for(std::size_t isample = 0;
      isample < samples.objective_weights.size(); isample++){
    double pointwise_gradient[gdim];
    double pointwise_hessian[hessian_count];
    value += samples.objective_weights[isample]
           * d_quafun_stepDistance<MFT,gdim,gdim,double>(
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
void evaluate_frozen_p2_stepdistance_hessian_by_gradient_ad(
    const FrozenP2StepDistanceSamples<gdim> &samples,
    const int active_control_point,
    FEBasis basis,
    double objective_power,
    double regularization,
    bool shape_volume,
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
    if(shape_volume){
      FrozenObjectiveValueAD::shape_volume_pointwise_gradient<gdim>(
          samples.regular_jacobian_transposes[isample].data(),
          samples.metrics[isample].data(),regular_basis_gradient.data(),
          objective_power,regularization,pointwise_gradient);
    }else{
      FrozenObjectiveValueAD::step_distance_pointwise_gradient<gdim>(
          samples.regular_jacobian_transposes[isample].data(),
          samples.metrics[isample].data(),regular_basis_gradient.data(),
          objective_power,regularization,pointwise_gradient);
    }
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
void check_stepdistance_edge_gradient(
    FEBasis geometry_basis,
    bool shape_volume = false)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_quadrature_order = 5;
  parameters.objective_p = 1.25;
  parameters.step_distance_regularization = 1.e-6;
  parameters.step_distance_shape_volume = shape_volume;
  parameters.step_distance_cavity_target_average = false;
  parameters.step_distance_barrier_rho0 = 10.0;
  parameters.step_distance_barrier_beta = shape_volume ? 1.e6 : 0.0;

  Mesh<MFT> mesh;
  initialize_curved_p2_stepdistance_element<MFT,gdim>(mesh,parameters);
  mesh.setBasis(geometry_basis);
  constexpr int active_edge_control_point = gdim + 1;
  const int active_point
      = mesh.ent2poi(gdim)(0,active_edge_control_point);
  const FrozenP2StepDistanceSamples<gdim> samples
      = capture_frozen_p2_stepdistance_samples<MFT,gdim>(mesh);

  double gradient[gdim];
  const double differentiated_value
      = d_metqua<MFT,gdim,gdim,QuaFun::StepDistance,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,
            active_edge_control_point,geometry_basis,DifVar::None,
            gradient,NULL,1.0);
  const double frozen_value
      = evaluate_frozen_p2_stepdistance_value<MFT,gdim>(mesh,samples);
  BOOST_CHECK_CLOSE_FRACTION(
      differentiated_value,frozen_value,3.e-14);

  double automatic_value;
  const std::array<double,gdim> automatic_gradient
      = evaluate_frozen_p2_stepdistance_gradient_by_ad<MFT,gdim>(
            mesh,samples,active_edge_control_point,&automatic_value);
  BOOST_CHECK_CLOSE_FRACTION(
      automatic_value,frozen_value,3.e-14);
  double maximum_automatic_gradient_error = 0.0;
  for(int component = 0; component < gdim; component++){
    const double error
        = std::abs(gradient[component] - automatic_gradient[component]);
    maximum_automatic_gradient_error
        = std::max(maximum_automatic_gradient_error,error);
    BOOST_CHECK_SMALL(
        error,5.e-14*(1.0 + std::abs(automatic_gradient[component])));
  }

  constexpr double finite_difference_step = 2.e-6;
  double maximum_finite_difference_gradient_error = 0.0;
  for(int component = 0; component < gdim; component++){
    const double coordinate = mesh.coord(active_point,component);
    mesh.coord(active_point,component)
        = coordinate + finite_difference_step;
    const double plus_value
        = evaluate_frozen_p2_stepdistance_value<MFT,gdim>(mesh,samples);
    mesh.coord(active_point,component)
        = coordinate - finite_difference_step;
    const double minus_value
        = evaluate_frozen_p2_stepdistance_value<MFT,gdim>(mesh,samples);
    mesh.coord(active_point,component) = coordinate;

    const double finite_difference_gradient
        = (plus_value - minus_value)/(2.0*finite_difference_step);
    const double error
        = std::abs(gradient[component] - finite_difference_gradient);
    maximum_finite_difference_gradient_error
        = std::max(maximum_finite_difference_gradient_error,error);
    BOOST_CHECK_SMALL(
        error,1.e-7*(1.0 + std::abs(finite_difference_gradient)));
  }
  BOOST_TEST_MESSAGE(
      "P2 frozen StepDistance "
      << (shape_volume ? "ShapeVolume " : "basic ")
      << "edge gradient: AD error="
      << maximum_automatic_gradient_error
      << ", finite-difference error="
      << maximum_finite_difference_gradient_error);
}

template<class MFT, int gdim>
void check_stepdistance_edge_hessian(
    FEBasis geometry_basis,
    bool shape_volume = false)
{
  constexpr int hessian_count = gdim*(gdim + 1)/2;
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_quadrature_order = 5;
  parameters.objective_p = 1.25;
  parameters.step_distance_regularization = 1.e-6;
  parameters.step_distance_shape_volume = shape_volume;
  parameters.step_distance_cavity_target_average = false;
  parameters.step_distance_barrier_rho0 = 10.0;
  parameters.step_distance_barrier_beta = shape_volume ? 1.e6 : 0.0;

  Mesh<MFT> mesh;
  initialize_curved_p2_stepdistance_element<MFT,gdim>(mesh,parameters);
  mesh.setBasis(geometry_basis);
  constexpr int active_edge_control_point = gdim + 1;
  const int active_point
      = mesh.ent2poi(gdim)(0,active_edge_control_point);
  const FrozenP2StepDistanceSamples<gdim> samples
      = capture_frozen_p2_stepdistance_samples<MFT,gdim>(mesh);

  double production_gradient[gdim];
  double production_hessian[hessian_count];
  const double production_value
      = d_metqua<MFT,gdim,gdim,QuaFun::StepDistance,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,
            active_edge_control_point,geometry_basis,DifVar::None,
            production_gradient,production_hessian,1.0);

  double reconstructed_gradient[gdim];
  double reconstructed_hessian[hessian_count];
  const double reconstructed_value
      = evaluate_frozen_p2_stepdistance_derivatives<MFT,gdim>(
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
  evaluate_frozen_p2_stepdistance_hessian_by_gradient_ad<gdim>(
      samples,active_edge_control_point,mesh.getBasis(),
      parameters.objective_p,parameters.step_distance_regularization,
      shape_volume,
      automatic_gradient,automatic_hessian);

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
    evaluate_frozen_p2_stepdistance_derivatives<MFT,gdim>(
        mesh,samples,active_edge_control_point,plus_gradient,NULL);
    mesh.coord(active_point,derivative)
        = coordinate - finite_difference_step;
    evaluate_frozen_p2_stepdistance_derivatives<MFT,gdim>(
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
      "P2 frozen StepDistance "
      << (shape_volume ? "ShapeVolume " : "basic ")
      << "edge Hessian: AD-of-gradient error="
      << maximum_automatic_hessian_error
      << ", finite-difference error="
      << maximum_finite_difference_hessian_error);
}

template<class MFT, int gdim>
void check_stepdistance_collapse_barrier(FEBasis geometry_basis)
{
  constexpr int hessian_count = gdim*(gdim + 1)/2;
  constexpr int active_edge_control_point = gdim + 1;
  constexpr double barrier_rho0 = 2.0;
  constexpr double barrier_beta = 0.4;
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_quadrature_order = 5;
  parameters.objective_p = 1.25;
  parameters.step_distance_regularization = 1.e-6;
  parameters.step_distance_shape_volume = false;
  parameters.step_distance_cavity_target_average = false;
  parameters.step_distance_barrier_rho0 = barrier_rho0;
  parameters.step_distance_barrier_beta = 0.0;

  Mesh<MFT> mesh;
  initialize_curved_p2_stepdistance_element<MFT,gdim>(mesh,parameters);
  mesh.setBasis(geometry_basis);
  const FrozenP2StepDistanceSamples<gdim> samples
      = capture_frozen_p2_stepdistance_samples<MFT,gdim>(mesh);

  // Qualify the degree-aware basis gradient independently of the production
  // helper used by the shared differentiated traversal.
  for(std::size_t isample = 0;
      isample < samples.barycentric_points.size(); isample++){
    const std::array<double,gdim> independent_gradient
        = independent_p2_regular_basis_gradient<gdim>(
              geometry_basis,active_edge_control_point,
              samples.barycentric_points[isample].data());
    const std::array<double,gdim> production_gradient
        = objective_regular_basis_gradient<gdim,2>(
              geometry_basis,active_edge_control_point,
              samples.barycentric_points[isample].data());
    for(int component = 0; component < gdim; component++){
      BOOST_CHECK_SMALL(
          production_gradient[component] - independent_gradient[component],
          5.e-14*(1.0 + std::abs(independent_gradient[component])));
    }
  }

  const auto evaluate_production = [&](double beta,
                                       double *gradient,
                                       double *hessian){
    mesh.param->step_distance_barrier_beta = beta;
    return d_metqua<MFT,gdim,gdim,QuaFun::StepDistance,double>(
        mesh,AsDeg::Pk,AsDeg::P1,0,active_edge_control_point,
        geometry_basis,DifVar::None,gradient,hessian,1.0);
  };

  double active_gradient[gdim];
  double inactive_gradient[gdim];
  double active_hessian[hessian_count];
  double inactive_hessian[hessian_count];
  const double active_value = evaluate_production(
      barrier_beta,active_gradient,active_hessian);
  const double inactive_value = evaluate_production(
      0.0,inactive_gradient,inactive_hessian);
  const double production_barrier_value = active_value - inactive_value;
  BOOST_REQUIRE_GT(production_barrier_value,1.e-6);

  mesh.param->step_distance_barrier_beta = barrier_beta;
  const double active_value_only
      = metqua<MFT,gdim,gdim,QuaFun::StepDistance,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,1.0);
  mesh.param->step_distance_barrier_beta = 0.0;
  const double inactive_value_only
      = metqua<MFT,gdim,gdim,QuaFun::StepDistance,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,1.0);
  BOOST_CHECK_SMALL(
      (active_value_only - inactive_value_only)
          - production_barrier_value,
      5.e-14*(1.0 + std::abs(production_barrier_value)));

  double automatic_value;
  const std::array<double,gdim> automatic_gradient
      = evaluate_frozen_p2_collapse_barrier_by_ad<gdim>(
            samples,active_edge_control_point,geometry_basis,
            barrier_rho0,barrier_beta,NULL,&automatic_value);
  BOOST_CHECK_SMALL(
      production_barrier_value - automatic_value,
      5.e-14*(1.0 + std::abs(automatic_value)));

  double automatic_gradient_from_formula[gdim];
  double automatic_hessian[hessian_count];
  evaluate_frozen_p2_collapse_barrier_hessian_by_gradient_ad<gdim>(
      samples,active_edge_control_point,geometry_basis,
      barrier_rho0,barrier_beta,
      automatic_gradient_from_formula,automatic_hessian);

  double maximum_automatic_gradient_error = 0.0;
  double maximum_automatic_hessian_error = 0.0;
  double production_barrier_gradient[gdim];
  double production_barrier_hessian[hessian_count];
  for(int component = 0; component < gdim; component++){
    production_barrier_gradient[component]
        = active_gradient[component] - inactive_gradient[component];
    const double value_ad_error
        = std::abs(production_barrier_gradient[component]
                   - automatic_gradient[component]);
    const double formula_ad_error
        = std::abs(production_barrier_gradient[component]
                   - automatic_gradient_from_formula[component]);
    maximum_automatic_gradient_error = std::max(
        maximum_automatic_gradient_error,
        std::max(value_ad_error,formula_ad_error));
    BOOST_CHECK_SMALL(
        value_ad_error,
        5.e-14*(1.0 + std::abs(automatic_gradient[component])));
    BOOST_CHECK_SMALL(
        formula_ad_error,
        5.e-14*(1.0
                + std::abs(automatic_gradient_from_formula[component])));
  }
  for(int entry = 0; entry < hessian_count; entry++){
    production_barrier_hessian[entry]
        = active_hessian[entry] - inactive_hessian[entry];
    const double error
        = std::abs(production_barrier_hessian[entry]
                   - automatic_hessian[entry]);
    maximum_automatic_hessian_error
        = std::max(maximum_automatic_hessian_error,error);
    BOOST_CHECK_SMALL(
        error,5.e-14*(1.0 + std::abs(automatic_hessian[entry])));
  }

  constexpr double finite_difference_step = 2.e-6;
  double maximum_finite_difference_gradient_error = 0.0;
  double maximum_finite_difference_hessian_error = 0.0;
  for(int derivative = 0; derivative < gdim; derivative++){
    std::array<double,gdim> plus_displacement{};
    std::array<double,gdim> minus_displacement{};
    plus_displacement[derivative] = finite_difference_step;
    minus_displacement[derivative] = -finite_difference_step;
    double plus_value;
    double minus_value;
    const std::array<double,gdim> plus_gradient
        = evaluate_frozen_p2_collapse_barrier_by_ad<gdim>(
              samples,active_edge_control_point,geometry_basis,
              barrier_rho0,barrier_beta,plus_displacement.data(),
              &plus_value);
    const std::array<double,gdim> minus_gradient
        = evaluate_frozen_p2_collapse_barrier_by_ad<gdim>(
              samples,active_edge_control_point,geometry_basis,
              barrier_rho0,barrier_beta,minus_displacement.data(),
              &minus_value);

    const double finite_difference_gradient
        = (plus_value - minus_value)/(2.0*finite_difference_step);
    const double gradient_error
        = std::abs(production_barrier_gradient[derivative]
                   - finite_difference_gradient);
    maximum_finite_difference_gradient_error = std::max(
        maximum_finite_difference_gradient_error,gradient_error);
    BOOST_CHECK_SMALL(
        gradient_error,
        1.e-7*(1.0 + std::abs(finite_difference_gradient)));

    for(int component = 0; component <= derivative; component++){
      const int entry = sym2idx(component,derivative);
      const double finite_difference_hessian
          = (plus_gradient[component] - minus_gradient[component])
           /(2.0*finite_difference_step);
      const double hessian_error
          = std::abs(production_barrier_hessian[entry]
                     - finite_difference_hessian);
      maximum_finite_difference_hessian_error = std::max(
          maximum_finite_difference_hessian_error,hessian_error);
      BOOST_CHECK_SMALL(
          hessian_error,
          1.e-7*(1.0 + std::abs(finite_difference_hessian)));
    }
  }

  BOOST_TEST_MESSAGE(
      "P2 frozen StepDistance collapse barrier: gradient AD error="
      << maximum_automatic_gradient_error
      << ", Hessian AD-of-gradient error="
      << maximum_automatic_hessian_error
      << ", gradient finite-difference error="
      << maximum_finite_difference_gradient_error
      << ", Hessian finite-difference error="
      << maximum_finite_difference_hessian_error);
}

template<class MFT, int gdim>
void check_shape_volume_rejects_singular_p2_trial(FEBasis geometry_basis)
{
  constexpr int hessian_count = gdim*(gdim + 1)/2;
  constexpr int node_count = getnnode(gdim,2);
  constexpr int active_edge_control_point = gdim + 1;
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_quadrature_order = 5;
  parameters.objective_p = 1.25;
  parameters.step_distance_regularization = 1.e-6;
  parameters.step_distance_shape_volume = true;
  parameters.step_distance_cavity_target_average = false;
  parameters.step_distance_barrier_rho0 = 10.0;
  parameters.step_distance_barrier_beta = 1.e6;

  Mesh<MFT> mesh;
  initialize_curved_p2_stepdistance_element<MFT,gdim>(mesh,parameters);
  mesh.setBasis(geometry_basis);
  for(int point = 0; point < node_count; point++){
    for(int component = 0; component < gdim; component++){
      mesh.coord(point,component) *= 1.e-4;
    }
  }

  const double value
      = metqua<MFT,gdim,gdim,QuaFun::StepDistance,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,1.0);
  BOOST_CHECK_EQUAL(
      value,step_distance_shape_volume_rejection_quality);

  double gradient[gdim];
  double hessian[hessian_count];
  const double differentiated_value
      = d_metqua<MFT,gdim,gdim,QuaFun::StepDistance,double>(
            mesh,AsDeg::Pk,AsDeg::P1,0,
            active_edge_control_point,geometry_basis,DifVar::None,
            gradient,hessian,1.0);
  BOOST_CHECK_EQUAL(differentiated_value,value);
  for(int component = 0; component < gdim; component++){
    BOOST_CHECK_EQUAL(gradient[component],0.0);
  }
  for(int entry = 0; entry < hessian_count; entry++){
    BOOST_CHECK_EQUAL(hessian[entry],0.0);
  }
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

BOOST_AUTO_TEST_CASE(p2_stepdistance_matches_independent_dense_reference)
{
  check_p2_stepdistance_against_dense_reference<MetricFieldFE,2>();
  check_p2_stepdistance_against_dense_reference<MetricFieldAnalytical,2>();
  check_p2_stepdistance_against_dense_reference<MetricFieldFE,3>();
  check_p2_stepdistance_against_dense_reference<MetricFieldAnalytical,3>();
}

BOOST_AUTO_TEST_CASE(p2_edge_stepdistance_gradient_matches_frozen_value_oracles)
{
  for(FEBasis basis : {FEBasis::Lagrange,FEBasis::Bezier}){
    check_stepdistance_edge_gradient<MetricFieldFE,2>(basis);
    check_stepdistance_edge_gradient<MetricFieldAnalytical,2>(basis);
    check_stepdistance_edge_gradient<MetricFieldFE,3>(basis);
    check_stepdistance_edge_gradient<MetricFieldAnalytical,3>(basis);
  }
}

BOOST_AUTO_TEST_CASE(p2_edge_stepdistance_hessian_matches_frozen_gradient_oracles)
{
  for(FEBasis basis : {FEBasis::Lagrange,FEBasis::Bezier}){
    check_stepdistance_edge_hessian<MetricFieldFE,2>(basis);
    check_stepdistance_edge_hessian<MetricFieldAnalytical,2>(basis);
    check_stepdistance_edge_hessian<MetricFieldFE,3>(basis);
    check_stepdistance_edge_hessian<MetricFieldAnalytical,3>(basis);
  }
}

BOOST_AUTO_TEST_CASE(p2_stepdistance_collapse_barrier_is_degree_aware)
{
  for(FEBasis basis : {FEBasis::Lagrange,FEBasis::Bezier}){
    check_stepdistance_collapse_barrier<MetricFieldFE,2>(basis);
    check_stepdistance_collapse_barrier<MetricFieldAnalytical,2>(basis);
    check_stepdistance_collapse_barrier<MetricFieldFE,3>(basis);
    check_stepdistance_collapse_barrier<MetricFieldAnalytical,3>(basis);
  }
}

BOOST_AUTO_TEST_CASE(p2_stepdistance_shape_volume_matches_dense_reference)
{
  check_p2_stepdistance_against_dense_reference<MetricFieldFE,2>(true);
  check_p2_stepdistance_against_dense_reference<MetricFieldAnalytical,2>(true);
  check_p2_stepdistance_against_dense_reference<MetricFieldFE,3>(true);
  check_p2_stepdistance_against_dense_reference<MetricFieldAnalytical,3>(true);
}

BOOST_AUTO_TEST_CASE(p2_edge_stepdistance_shape_volume_gradient_is_validated)
{
  for(FEBasis basis : {FEBasis::Lagrange,FEBasis::Bezier}){
    check_stepdistance_edge_gradient<MetricFieldFE,2>(basis,true);
    check_stepdistance_edge_gradient<MetricFieldAnalytical,2>(basis,true);
    check_stepdistance_edge_gradient<MetricFieldFE,3>(basis,true);
    check_stepdistance_edge_gradient<MetricFieldAnalytical,3>(basis,true);
  }
}

BOOST_AUTO_TEST_CASE(p2_edge_stepdistance_shape_volume_hessian_is_validated)
{
  for(FEBasis basis : {FEBasis::Lagrange,FEBasis::Bezier}){
    check_stepdistance_edge_hessian<MetricFieldFE,2>(basis,true);
    check_stepdistance_edge_hessian<MetricFieldAnalytical,2>(basis,true);
    check_stepdistance_edge_hessian<MetricFieldFE,3>(basis,true);
    check_stepdistance_edge_hessian<MetricFieldAnalytical,3>(basis,true);
  }
}

BOOST_AUTO_TEST_CASE(p2_stepdistance_shape_volume_rejects_singular_trials)
{
  for(FEBasis basis : {FEBasis::Lagrange,FEBasis::Bezier}){
    check_shape_volume_rejects_singular_p2_trial<MetricFieldFE,2>(basis);
    check_shape_volume_rejects_singular_p2_trial<MetricFieldAnalytical,2>(basis);
    check_shape_volume_rejects_singular_p2_trial<MetricFieldFE,3>(basis);
    check_shape_volume_rejects_singular_p2_trial<MetricFieldAnalytical,3>(basis);
  }
}

} // namespace Metris
