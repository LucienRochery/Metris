// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_quality_sizeshape_pointwise_value

#include <boost/test/included/unit_test.hpp>

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "aux_exceptions.hxx"
#include "quality/quafun_distortion.hxx"
#include "quality/quafun_sizeshape.hxx"

#include <cmath>

namespace Metris
{

namespace
{

void initialize_ideal_triangle(Mesh<MetricFieldFE> &msh,
                               MetrisParameters &parameters)
{
  msh.idim = 2;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &parameters;
  msh.set_npoin(3);
  msh.set_nface(1);

  msh.coord(0,0) = 0.0;
  msh.coord(0,1) = 0.0;
  msh.coord(1,0) = 1.0;
  msh.coord(1,1) = 0.0;
  msh.coord(2,0) = 0.5;
  msh.coord(2,1) = std::sqrt(3.0)/2.0;
  for(int inode = 0; inode < 3; inode++){
    msh.fac2poi(0,inode) = inode;
  }
}

void initialize_ideal_tetrahedron(Mesh<MetricFieldFE> &msh,
                                  MetrisParameters &parameters)
{
  msh.idim = 3;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &parameters;
  msh.set_npoin(4);
  msh.set_nelem(1);

  msh.coord(0,0) = 0.0;
  msh.coord(0,1) = 0.0;
  msh.coord(0,2) = 0.0;
  msh.coord(1,0) = -std::sqrt(3.0)/2.0;
  msh.coord(1,1) = 0.5;
  msh.coord(1,2) = 0.0;
  msh.coord(2,0) = -std::sqrt(3.0)/2.0;
  msh.coord(2,1) = -0.5;
  msh.coord(2,2) = 0.0;
  msh.coord(3,0) = -1.0/std::sqrt(3.0);
  msh.coord(3,1) = 0.0;
  msh.coord(3,2) = std::sqrt(2.0/3.0);
  for(int inode = 0; inode < 4; inode++){
    msh.tet2poi(0,inode) = inode;
  }
}

double evaluate_triangle_sizeshape(Mesh<MetricFieldFE> &msh,
                                   const double *metric)
{
  const double barycenter[3] = {1.0/3.0,1.0/3.0,1.0/3.0};
  return quafun_sizeshape<MetricFieldFE,2,2,double>(
      msh,AsDeg::P1,AsDeg::P1,msh.fac2poi[0],barycenter,metric);
}

double evaluate_tetrahedron_sizeshape(Mesh<MetricFieldFE> &msh,
                                      const double *metric)
{
  const double barycenter[4] = {0.25,0.25,0.25,0.25};
  return quafun_sizeshape<MetricFieldFE,3,3,double>(
      msh,AsDeg::P1,AsDeg::P1,msh.tet2poi[0],barycenter,metric);
}

} // namespace

BOOST_AUTO_TEST_CASE(test_ideal_state_is_zero_for_supported_powers)
{
  MetrisParameters parameters;
  parameters.iverb = 0;

  Mesh<MetricFieldFE> triangle;
  initialize_ideal_triangle(triangle,parameters);
  const double triangle_identity_metric[3] = {1.0,0.0,1.0};

  Mesh<MetricFieldFE> tetrahedron;
  initialize_ideal_tetrahedron(tetrahedron,parameters);
  const double tetrahedron_identity_metric[6]
      = {1.0,0.0,1.0,0.0,0.0,1.0};

  const double powers[3] = {1.0,2.0,1.5};
  for(int ipower = 0; ipower < 3; ipower++){
    parameters.objective_p = powers[ipower];
    BOOST_CHECK_SMALL(
        evaluate_triangle_sizeshape(triangle,triangle_identity_metric),
        1.0e-30);
    BOOST_CHECK_SMALL(
        evaluate_tetrahedron_sizeshape(tetrahedron,
                                       tetrahedron_identity_metric),
        1.0e-30);
  }
}

BOOST_AUTO_TEST_CASE(test_anisotropic_values_and_power_are_explicit)
{
  MetrisParameters parameters;
  parameters.iverb = 0;

  Mesh<MetricFieldFE> triangle;
  initialize_ideal_triangle(triangle,parameters);
  const double triangle_metric[3] = {4.0,0.0,1.0};
  const double triangle_error = 25.0*(1.0 + 1.0/16.0)/8.0 - 1.0;

  Mesh<MetricFieldFE> tetrahedron;
  initialize_ideal_tetrahedron(tetrahedron,parameters);
  const double tetrahedron_metric[6] = {4.0,0.0,1.0,0.0,0.0,1.0};
  const double tetrahedron_error
      = 216.0*(1.0 + 1.0/16.0)/54.0 - 1.0;

  const double powers[3] = {1.0,2.0,1.5};
  for(int ipower = 0; ipower < 3; ipower++){
    parameters.objective_p = powers[ipower];
    const double triangle_value
        = evaluate_triangle_sizeshape(triangle,triangle_metric);
    const double tetrahedron_value
        = evaluate_tetrahedron_sizeshape(tetrahedron,tetrahedron_metric);

    BOOST_CHECK(triangle_value >= 0.0);
    BOOST_CHECK(tetrahedron_value >= 0.0);
    BOOST_CHECK_CLOSE_FRACTION(
        triangle_value,std::pow(triangle_error,powers[ipower]),5.0e-14);
    BOOST_CHECK_CLOSE_FRACTION(
        tetrahedron_value,std::pow(tetrahedron_error,powers[ipower]),
        5.0e-14);
  }
}

BOOST_AUTO_TEST_CASE(test_sizeshape_is_independent_of_classical_opt_power)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_p = 1.5;

  Mesh<MetricFieldFE> msh;
  initialize_ideal_triangle(msh,parameters);
  const double metric[3] = {4.0,0.0,1.0};

  parameters.opt_power = 1;
  const double direct_setting_value = evaluate_triangle_sizeshape(msh,metric);
  parameters.opt_power = -1;
  const double inverse_setting_value = evaluate_triangle_sizeshape(msh,metric);

  BOOST_CHECK_EQUAL(inverse_setting_value,direct_setting_value);
}

BOOST_AUTO_TEST_CASE(test_invalid_states_are_not_hidden_by_ideal_clamping)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_p = 1.5;

  Mesh<MetricFieldFE> msh;
  initialize_ideal_triangle(msh,parameters);
  const double barycenter[3] = {1.0/3.0,1.0/3.0,1.0/3.0};

  const double indefinite_metric[3] = {2.0,0.0,-1.0};
  BOOST_CHECK_THROW(
      (quafun_sizeshape<MetricFieldFE,2,2,double>(
          msh,AsDeg::P1,AsDeg::P1,msh.fac2poi[0],barycenter,
          indefinite_metric)),
      MetrisExcept);

  msh.coord(2,0) = 0.0;
  msh.coord(2,1) = 0.0;
  const double identity_metric[3] = {1.0,0.0,1.0};
  parameters.opt_power = -1;
  BOOST_CHECK_THROW(
      (quafun_sizeshape<MetricFieldFE,2,2,double>(
          msh,AsDeg::P1,AsDeg::P1,msh.fac2poi[0],barycenter,
          identity_metric)),
      MetrisExcept);

  const double classical_inverse
      = quafun_distortion<MetricFieldFE,2,2,double>(
          msh,AsDeg::P1,AsDeg::P1,msh.fac2poi[0],barycenter,
          identity_metric);
  BOOST_CHECK_EQUAL(classical_inverse,0.0);
}

} // namespace Metris
