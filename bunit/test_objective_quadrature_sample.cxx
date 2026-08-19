// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_objective_quadrature_sample

#include <boost/test/included/unit_test.hpp>

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "ho_constants.hxx"
#include "linalg/det.hxx"
#include "low_eval.hxx"
#include "low_geo/measure.hxx"
#include "quality/aux_volumeMeasure.hxx"
#include "quality/objective_quadrature_derivatives.hxx"
#include "quality/objective_quadrature_sample.hxx"
#include "quality/simplex_quadrature.hxx"

#include <algorithm>
#include <cmath>

namespace Metris
{

namespace
{

template<int gdim>
void initialize_p1_fe_element(Mesh<MetricFieldFE> &msh,
                              MetrisParameters &parameters)
{
  msh.idim = gdim;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &parameters;
  msh.set_npoin(gdim + 1);
  msh.set_nentt(gdim,1);
  msh.met.forceBasisFlag(FEBasis::Lagrange);
  msh.met.forceSpaceFlag(MetSpace::Exp);

  if constexpr(gdim == 2){
    msh.coord(0,0) = 0.05; msh.coord(0,1) = 0.10;
    msh.coord(1,0) = 1.15; msh.coord(1,1) = 0.18;
    msh.coord(2,0) = 0.25; msh.coord(2,1) = 0.92;
  }else{
    msh.coord(0,0) = 0.05; msh.coord(0,1) = 0.10;
    msh.coord(0,2) = 0.08;
    msh.coord(1,0) = 1.05; msh.coord(1,1) = 0.18;
    msh.coord(1,2) = 0.04;
    msh.coord(2,0) = 0.22; msh.coord(2,1) = 0.98;
    msh.coord(2,2) = 0.16;
    msh.coord(3,0) = 0.14; msh.coord(3,1) = 0.28;
    msh.coord(3,2) = 0.91;
  }

  intAr2 &element_to_point = msh.ent2poi(gdim);
  for(int inode = 0; inode < gdim + 1; inode++){
    element_to_point(0,inode) = inode;
  }

  constexpr int nmetric = gdim*(gdim + 1)/2;
  for(int ipoin = 0; ipoin < gdim + 1; ipoin++){
    for(int imetric = 0; imetric < nmetric; imetric++){
      msh.met(ipoin,imetric) = 0.0;
    }
    if constexpr(gdim == 2){
      msh.met(ipoin,0) = 1.2 + 0.2*ipoin;
      msh.met(ipoin,1) = 0.04 - 0.01*ipoin;
      msh.met(ipoin,2) = 0.9 + 0.15*ipoin;
    }else{
      msh.met(ipoin,0) = 1.3 + 0.15*ipoin;
      msh.met(ipoin,1) = 0.03;
      msh.met(ipoin,2) = 1.0 + 0.12*ipoin;
      msh.met(ipoin,3) = -0.02;
      msh.met(ipoin,4) = 0.025;
      msh.met(ipoin,5) = 0.85 + 0.10*ipoin;
    }
  }
}

template<int gdim, int mshdeg>
void check_regular_basis_gradient_against_geometry_perturbation(
    Mesh<MetricFieldFE> &msh,
    FEBasis dofbas,
    int local_control_point,
    const double *barycentric_coordinates)
{
  const std::array<double,gdim> regular_basis_gradient
      = objective_regular_basis_gradient<gdim,mshdeg>(
            dofbas,local_control_point,barycentric_coordinates);
  const SimplexQuadraturePointView<gdim> quadrature_point
      = {barycentric_coordinates,1.0};
  const int *nodes = msh.ent2poi(gdim)[0];
  constexpr double perturbation = 1.e-7;

  for(int perturbed_component = 0;
      perturbed_component < gdim;
      perturbed_component++){
    const int ipoin = nodes[local_control_point];
    const double coordinate = msh.coord(ipoin,perturbed_component);

    msh.coord(ipoin,perturbed_component) = coordinate + perturbation;
    const ObjectiveQuadratureSample<gdim,gdim,mshdeg> plus_sample
        = prepare_objective_quadrature_sample<
              MetricFieldFE,gdim,gdim,mshdeg>(
                  msh,AsDeg::P1,nodes,quadrature_point,
                  ObjectiveQuadratureTheta::ReferenceAverage);

    msh.coord(ipoin,perturbed_component) = coordinate - perturbation;
    const ObjectiveQuadratureSample<gdim,gdim,mshdeg> minus_sample
        = prepare_objective_quadrature_sample<
              MetricFieldFE,gdim,gdim,mshdeg>(
                  msh,AsDeg::P1,nodes,quadrature_point,
                  ObjectiveQuadratureTheta::ReferenceAverage);
    msh.coord(ipoin,perturbed_component) = coordinate;

    for(int iregular = 0; iregular < gdim; iregular++){
      for(int component = 0; component < gdim; component++){
        const int ientry = iregular*gdim + component;
        const double finite_difference
            = (plus_sample.regular_jacobian_transpose[ientry]
               - minus_sample.regular_jacobian_transpose[ientry])
             /(2.0*perturbation);
        const double expected_derivative
            = component == perturbed_component
            ? regular_basis_gradient[iregular]
            : 0.0;
        BOOST_CHECK_SMALL(
            finite_difference - expected_derivative,2.e-9);
      }
    }
  }
}

template<int gdim>
void check_p1_fe_sample()
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  Mesh<MetricFieldFE> msh;
  initialize_p1_fe_element<gdim>(msh,parameters);

  const SimplexQuadratureView<gdim> quadrature
      = get_vertex_barycenter_quadrature<gdim>();
  const SimplexQuadraturePointView<gdim> quadrature_point
      = quadrature[quadrature.size() - 1];
  const int *nodes = msh.ent2poi(gdim)[0];

  const ObjectiveQuadratureSample<gdim,gdim,1> reference_sample
      = prepare_objective_quadrature_sample<
            MetricFieldFE,gdim,gdim,1>(
                msh,AsDeg::P1,nodes,quadrature_point,
                ObjectiveQuadratureTheta::ReferenceAverage);
  const ObjectiveQuadratureSample<gdim,gdim,1> physical_sample
      = prepare_objective_quadrature_sample<
            MetricFieldFE,gdim,gdim,1>(
                msh,AsDeg::P1,nodes,quadrature_point,
                ObjectiveQuadratureTheta::PhysicalMeasure);
  const ObjectiveQuadratureSample<gdim,gdim,1> physical_metric_sample
      = prepare_objective_quadrature_sample<
            MetricFieldFE,gdim,gdim,1>(
                msh,AsDeg::P1,nodes,quadrature_point,
                ObjectiveQuadratureTheta::PhysicalMetricMeasure);
  const ObjectiveQuadratureSample<gdim,gdim,1> regular_metric_sample
      = prepare_objective_quadrature_sample<
            MetricFieldFE,gdim,gdim,1>(
                msh,AsDeg::P1,nodes,quadrature_point,
                ObjectiveQuadratureTheta::RegularMetricMeasure);

  BOOST_CHECK_EQUAL(reference_sample.geometry_degree,1);
  BOOST_CHECK(reference_sample.theta_is_valid);
  BOOST_CHECK(physical_sample.theta_is_valid);
  BOOST_CHECK(physical_metric_sample.theta_is_valid);
  BOOST_CHECK(regular_metric_sample.theta_is_valid);
  BOOST_CHECK_CLOSE_FRACTION(
      reference_sample.quadrature_weight,quadrature_point.weight,1.e-15);
  BOOST_CHECK_CLOSE_FRACTION(reference_sample.theta,1.0,1.e-15);

  for(int ibary = 0; ibary < gdim + 1; ibary++){
    BOOST_CHECK_CLOSE_FRACTION(
        reference_sample.barycentric_coordinates[ibary],
        quadrature_point.bary[ibary],1.e-15);
  }

  for(int component = 0; component < gdim; component++){
    double expected_coordinate = 0.0;
    for(int inode = 0; inode < gdim + 1; inode++){
      expected_coordinate += quadrature_point.bary[inode]
                           * msh.coord(nodes[inode],component);
    }
    BOOST_CHECK_CLOSE_FRACTION(
        reference_sample.physical_coordinates[component],
        expected_coordinate,1.e-14);
  }

  double expected_jacobian[gdim*gdim];
  double ignored_coordinates[gdim];
  if constexpr(gdim == 2){
    eval2<gdim,1>(msh.coord,nodes,msh.getBasis(),
                  DifVar::Bary,DifVar::None,quadrature_point.bary,
                  ignored_coordinates,expected_jacobian,NULL);
  }else{
    eval3<gdim,1>(msh.coord,nodes,msh.getBasis(),
                  DifVar::Bary,DifVar::None,quadrature_point.bary,
                  ignored_coordinates,expected_jacobian,NULL);
  }
  for(int ientry = 0; ientry < gdim*gdim; ientry++){
    BOOST_CHECK_SMALL(
        reference_sample.canonical_jacobian_transpose[ientry]
        - expected_jacobian[ientry],1.e-14);
  }

  constexpr int nmetric = gdim*(gdim + 1)/2;
  double expected_metric[nmetric];
  msh.met.getMetBary(
      AsDeg::P1,DifVar::None,MetSpace::Exp,nodes,gdim,
      quadrature_point.bary,expected_metric,NULL);
  for(int imetric = 0; imetric < nmetric; imetric++){
    BOOST_CHECK_CLOSE_FRACTION(
        reference_sample.metric[imetric],expected_metric[imetric],1.e-14);
  }

  const double expected_physical_measure
      = std::abs(getmeasentP1<gdim>(nodes,msh.coord));
  BOOST_CHECK_CLOSE_FRACTION(
      physical_sample.theta,expected_physical_measure,1.e-14);
  BOOST_CHECK_CLOSE_FRACTION(
      physical_metric_sample.theta,
      expected_physical_measure*std::sqrt(detsym<gdim>(expected_metric)),
      1.e-14);

  double expected_regular_theta;
  VolumeMeasureHelpers::eval_theta_fixed_metric_grad<gdim,gdim,double>(
      regular_metric_sample.regular_jacobian_transpose.data(),
      expected_metric,NULL,&expected_regular_theta,NULL);
  BOOST_CHECK_CLOSE_FRACTION(
      regular_metric_sample.theta,expected_regular_theta,1.e-14);

  const SimplexQuadraturePointView<gdim> vertex_point = quadrature[1];
  const ObjectiveQuadratureSample<gdim,gdim,1> vertex_sample
      = prepare_objective_quadrature_sample<
            MetricFieldFE,gdim,gdim,1>(
                msh,AsDeg::P1,nodes,vertex_point,
                ObjectiveQuadratureTheta::ReferenceAverage);
  for(int imetric = 0; imetric < nmetric; imetric++){
    BOOST_CHECK_EQUAL(vertex_sample.metric[imetric],msh.met(nodes[1],imetric));
  }

  check_regular_basis_gradient_against_geometry_perturbation<gdim,1>(
      msh,FEBasis::Lagrange,1,quadrature_point.bary);
}

void analytical_metric_2d(
    const AnaMetCtx *,
    const double *coordinate,
    double scale,
    int derivative_order,
    double *metric,
    double *metric_derivative)
{
  metric[0] = scale*(1.4 + 0.2*coordinate[0]);
  metric[1] = scale*(0.05 + 0.03*coordinate[1]);
  metric[2] = scale*(0.9 + 0.1*coordinate[0] + 0.15*coordinate[1]);
  if(derivative_order == 0) return;
  for(int ientry = 0; ientry < 6; ientry++){
    metric_derivative[ientry] = 0.0;
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(test_p1_triangle_fe_sample)
{
  check_p1_fe_sample<2>();
}

BOOST_AUTO_TEST_CASE(test_p1_tetrahedron_fe_sample)
{
  check_p1_fe_sample<3>();
}

BOOST_AUTO_TEST_CASE(test_analytical_metric_sampling)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.setAnalyticalMetric(AnaMetFun(analytical_metric_2d));

  Mesh<MetricFieldAnalytical> msh;
  msh.idim = 2;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &parameters;
  msh.set_npoin(3);
  msh.set_nface(1);
  msh.met.forceBasisFlag(FEBasis::Lagrange);
  msh.met.forceSpaceFlag(MetSpace::Exp);
  msh.met.setAnalyticalMetric(parameters);

  msh.coord(0,0) = 0.0; msh.coord(0,1) = 0.0;
  msh.coord(1,0) = 1.0; msh.coord(1,1) = 0.1;
  msh.coord(2,0) = 0.2; msh.coord(2,1) = 0.9;
  for(int inode = 0; inode < 3; inode++){
    msh.fac2poi(0,inode) = inode;
    msh.met.getMetPhys(DifVar::None,MetSpace::Exp,
                       msh.coord[inode],msh.met[inode],NULL);
  }

  // Exact rule vertices retain the historical nodal sample, even if it has
  // intentionally been made different from a fresh analytical evaluation.
  msh.met(0,0) += 0.25;
  const SimplexQuadratureView<2> quadrature
      = get_vertex_barycenter_quadrature<2>();
  const ObjectiveQuadratureSample<2,2,1> vertex_sample
      = prepare_objective_quadrature_sample<
            MetricFieldAnalytical,2,2,1>(
                msh,AsDeg::P1,msh.fac2poi[0],quadrature[0],
                ObjectiveQuadratureTheta::ReferenceAverage);
  BOOST_CHECK_EQUAL(vertex_sample.metric[0],msh.met(0,0));

  const ObjectiveQuadratureSample<2,2,1> barycenter_sample
      = prepare_objective_quadrature_sample<
            MetricFieldAnalytical,2,2,1>(
                msh,AsDeg::P1,msh.fac2poi[0],quadrature[3],
                ObjectiveQuadratureTheta::ReferenceAverage);
  double expected_metric[3];
  msh.met.getMetPhys(
      DifVar::None,MetSpace::Exp,
      barycenter_sample.physical_coordinates.data(),expected_metric,NULL);
  for(int imetric = 0; imetric < 3; imetric++){
    BOOST_CHECK_CLOSE_FRACTION(
        barycenter_sample.metric[imetric],expected_metric[imetric],1.e-14);
  }
}

BOOST_AUTO_TEST_CASE(test_quadratic_geometry_degree_is_not_baked_to_p1)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  Mesh<MetricFieldFE> msh;
  msh.idim = 2;
  msh.curdeg = 2;
  msh.strdeg = 2;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &parameters;
  msh.set_npoin(6);
  msh.set_nface(1);
  msh.met.forceBasisFlag(FEBasis::Lagrange);
  msh.met.forceSpaceFlag(MetSpace::Exp);

  constexpr auto ordering = ORDELT(2);
  for(int inode = 0; inode < 6; inode++){
    msh.fac2poi(0,inode) = inode;
    const double bary1 = ordering[2][inode][1]/2.0;
    const double bary2 = ordering[2][inode][2]/2.0;
    msh.coord(inode,0) = bary1;
    msh.coord(inode,1) = bary2;
    msh.met(inode,0) = 1.0;
    msh.met(inode,1) = 0.0;
    msh.met(inode,2) = 1.0;
  }
  msh.coord(3,1) += 0.2;

  const double barycentric_coordinates[3] = {0.2,0.3,0.5};
  const SimplexQuadraturePointView<2> quadrature_point
      = {barycentric_coordinates,0.37};
  const ObjectiveQuadratureSample<2,2,2> quadratic_sample
      = prepare_objective_quadrature_sample<MetricFieldFE,2,2,2>(
            msh,AsDeg::P1,msh.fac2poi[0],quadrature_point,
            ObjectiveQuadratureTheta::PhysicalMeasure);
  const ObjectiveQuadratureSample<2,2,1> linear_sample
      = prepare_objective_quadrature_sample<MetricFieldFE,2,2,1>(
            msh,AsDeg::P1,msh.fac2poi[0],quadrature_point,
            ObjectiveQuadratureTheta::PhysicalMeasure);

  BOOST_CHECK_EQUAL(quadratic_sample.geometry_degree,2);
  double expected_coordinates[2];
  double expected_jacobian[4];
  eval2<2,2>(msh.coord,msh.fac2poi[0],msh.getBasis(),
             DifVar::Bary,DifVar::None,barycentric_coordinates,
             expected_coordinates,expected_jacobian,NULL);
  for(int component = 0; component < 2; component++){
    BOOST_CHECK_CLOSE_FRACTION(
        quadratic_sample.physical_coordinates[component],
        expected_coordinates[component],1.e-14);
  }
  for(int ientry = 0; ientry < 4; ientry++){
    BOOST_CHECK_SMALL(
        quadratic_sample.canonical_jacobian_transpose[ientry]
        - expected_jacobian[ientry],1.e-14);
  }

  double geometry_difference = 0.0;
  for(int ientry = 0; ientry < 4; ientry++){
    geometry_difference = std::max(
        geometry_difference,
        std::abs(quadratic_sample.canonical_jacobian_transpose[ientry]
                 - linear_sample.canonical_jacobian_transpose[ientry]));
  }
  BOOST_TEST(geometry_difference > 1.e-3);

  check_regular_basis_gradient_against_geometry_perturbation<2,2>(
      msh,FEBasis::Lagrange,3,barycentric_coordinates);
  msh.forceBasisFlag(FEBasis::Bezier);
  check_regular_basis_gradient_against_geometry_perturbation<2,2>(
      msh,FEBasis::Bezier,3,barycentric_coordinates);
}

BOOST_AUTO_TEST_CASE(test_invalid_theta_is_reported_without_early_throw)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  Mesh<MetricFieldFE> msh;
  initialize_p1_fe_element<2>(msh,parameters);
  msh.coord(2,0) = 2.0;
  msh.coord(2,1) = 0.0;
  msh.coord(0,0) = 0.0;
  msh.coord(0,1) = 0.0;
  msh.coord(1,0) = 1.0;
  msh.coord(1,1) = 0.0;

  const SimplexQuadratureView<2> quadrature
      = get_vertex_barycenter_quadrature<2>();
  const ObjectiveQuadratureSample<2,2,1> sample
      = prepare_objective_quadrature_sample<MetricFieldFE,2,2,1>(
            msh,AsDeg::P1,msh.fac2poi[0],quadrature[3],
            ObjectiveQuadratureTheta::RegularMetricMeasure);
  BOOST_CHECK(!sample.theta_is_valid);
  BOOST_CHECK_EQUAL(sample.theta,0.0);
}

} // namespace Metris
