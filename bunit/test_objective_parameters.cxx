// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_objective_parameters

#include <boost/test/included/unit_test.hpp>

#include "MetrisRunner/MetrisParameters.hxx"
#include "aux_exceptions.hxx"
#include "metris_options.hxx"

namespace Metris
{

namespace
{

void construct_parameters(MetrisOptions &options)
{
  MetrisParameters parameters(options);
  (void)parameters;
}

} // namespace

BOOST_AUTO_TEST_CASE(test_objective_power_default_and_command_line_names)
{
  MetrisParameters default_parameters;
  BOOST_CHECK_EQUAL(default_parameters.objective_p,1.0);
  BOOST_CHECK_EQUAL(default_parameters.objective_quadrature_order,-1);
  BOOST_CHECK_NO_THROW(default_parameters.checkParameters());

  cargHandler automatic_arguments(
      "--objective-quadrature-order -1 "
      "--out objective_parameter_test.meshb --verb 0");
  MetrisOptions automatic_options(automatic_arguments.c,
                                  automatic_arguments.v);
  MetrisParameters automatic_parameters(automatic_options);
  BOOST_CHECK_EQUAL(automatic_parameters.objective_quadrature_order,-1);
  BOOST_CHECK_NO_THROW(automatic_parameters.checkParameters());

  cargHandler historical_arguments(
      "--objective-quadrature-order 0 "
      "--out objective_parameter_test.meshb --verb 0");
  MetrisOptions historical_options(historical_arguments.c,
                                   historical_arguments.v);
  MetrisParameters historical_parameters(historical_options);
  BOOST_CHECK_EQUAL(historical_parameters.objective_quadrature_order,0);
  BOOST_CHECK_NO_THROW(historical_parameters.checkParameters());

  cargHandler canonical_arguments(
      "--objective-p 1.75 --objective-quadrature-order 5 "
      "--out objective_parameter_test.meshb --verb 0");
  MetrisOptions canonical_options(canonical_arguments.c,
                                  canonical_arguments.v);
  MetrisParameters canonical_parameters(canonical_options);
  BOOST_CHECK_EQUAL(canonical_parameters.objective_p,1.75);
  BOOST_CHECK_EQUAL(canonical_parameters.objective_quadrature_order,5);
  BOOST_CHECK_NO_THROW(canonical_parameters.checkParameters());

  MetrisParameters degree_two_parameters;
  degree_two_parameters.objective_quadrature_order = 2;
  BOOST_CHECK_NO_THROW(degree_two_parameters.checkParameters());

  MetrisParameters degree_three_parameters;
  degree_three_parameters.objective_quadrature_order = 3;
  BOOST_CHECK_NO_THROW(degree_three_parameters.checkParameters());

  MetrisParameters degree_four_parameters;
  degree_four_parameters.objective_quadrature_order = 4;
  BOOST_CHECK_NO_THROW(degree_four_parameters.checkParameters());

  cargHandler compatibility_arguments(
      "--step-distance-p 1.75 --out objective_parameter_test.meshb --verb 0");
  MetrisOptions compatibility_options(compatibility_arguments.c,
                                      compatibility_arguments.v);
  MetrisParameters compatibility_parameters(compatibility_options);
  BOOST_CHECK_EQUAL(compatibility_parameters.objective_p,1.75);
  BOOST_CHECK_NO_THROW(compatibility_parameters.checkParameters());

  cargHandler identical_arguments(
      "--objective-p 2.5 --step-distance-p 2.5 "
      "--out objective_parameter_test.meshb --verb 0");
  MetrisOptions identical_options(identical_arguments.c,
                                  identical_arguments.v);
  MetrisParameters identical_parameters(identical_options);
  BOOST_CHECK_EQUAL(identical_parameters.objective_p,2.5);
}

BOOST_AUTO_TEST_CASE(test_objective_power_rejects_invalid_or_conflicting_values)
{
  MetrisParameters invalid_parameters;
  invalid_parameters.objective_p = 0.999;
  BOOST_CHECK_THROW(invalid_parameters.checkParameters(),MetrisExcept);

  MetrisParameters negative_quadrature_order;
  negative_quadrature_order.objective_quadrature_order = -2;
  BOOST_CHECK_THROW(
      negative_quadrature_order.checkParameters(),MetrisExcept);

  MetrisParameters unsupported_quadrature_order;
  unsupported_quadrature_order.objective_quadrature_order = 1;
  BOOST_CHECK_THROW(
      unsupported_quadrature_order.checkParameters(),MetrisExcept);

  unsupported_quadrature_order.objective_quadrature_order = 6;
  BOOST_CHECK_THROW(
      unsupported_quadrature_order.checkParameters(),MetrisExcept);

  cargHandler conflicting_arguments(
      "--objective-p 1.5 --step-distance-p 2.0 --verb 0");
  MetrisOptions conflicting_options(conflicting_arguments.c,
                                    conflicting_arguments.v);
  BOOST_CHECK_THROW(construct_parameters(conflicting_options),MetrisExcept);
}

} // namespace Metris
