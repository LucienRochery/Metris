// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_smoothing_progress

#include <boost/test/included/unit_test.hpp>

#include "MetrisRunner/MetrisParameters.hxx"
#include "aux_exceptions.hxx"
#include "metris_options.hxx"
#include "smoothing/smoothing_progress.hxx"

#include <array>

namespace Metris
{

BOOST_AUTO_TEST_CASE(relative_progress_is_invariant_under_objective_scaling)
{
  constexpr double value_before = 2.0;
  constexpr double substantive_value_after = 1.98;
  constexpr double negligible_value_after = 1.995;
  constexpr double relative_tolerance = 0.005;
  constexpr std::array<double,5> scales = {
      1.0e-12,1.0e-6,1.0,1.0e6,1.0e12};

  for (const double scale : scales)
  {
    BOOST_CHECK(smoothing_reduction_is_substantive(
        scale * value_before,
        scale * substantive_value_after,
        relative_tolerance));
    BOOST_CHECK(!smoothing_reduction_is_substantive(
        scale * value_before,
        scale * negligible_value_after,
        relative_tolerance));
  }
}

BOOST_AUTO_TEST_CASE(either_local_statistic_can_reactivate_the_neighborhood)
{
  constexpr double relative_tolerance = 0.005;

  BOOST_CHECK(smoothing_neighborhood_should_be_reactivated(
      10.0,9.9,4.0,3.999,relative_tolerance,relative_tolerance));
  BOOST_CHECK(smoothing_neighborhood_should_be_reactivated(
      10.0,9.999,4.0,3.9,relative_tolerance,relative_tolerance));
  BOOST_CHECK(!smoothing_neighborhood_should_be_reactivated(
      10.0,9.999,4.0,3.999,relative_tolerance,relative_tolerance));
}

BOOST_AUTO_TEST_CASE(invalid_or_non_improving_values_do_not_reactivate)
{
  BOOST_CHECK(!smoothing_reduction_is_substantive(0.0,0.0,0.005));
  BOOST_CHECK(!smoothing_reduction_is_substantive(1.0,1.0,0.005));
  BOOST_CHECK(!smoothing_reduction_is_substantive(1.0,1.1,0.005));
  BOOST_CHECK(!smoothing_reduction_is_substantive(1.0,-0.1,0.005));
}

BOOST_AUTO_TEST_CASE(adaptation_fallback_is_opt_in)
{
  MetrisParameters default_parameters;
  BOOST_CHECK(!default_parameters.adp_quality_smoothing);

  cargHandler arguments(
      "--adp-quality-smoothing "
      "--out smoothing_progress_test.meshb --verb 0");
  MetrisOptions options(arguments.c,arguments.v);
  MetrisParameters enabled_parameters(options);
  BOOST_CHECK(enabled_parameters.adp_quality_smoothing);
  BOOST_CHECK_NO_THROW(enabled_parameters.checkParameters());
}

BOOST_AUTO_TEST_CASE(relative_smoothing_tolerance_is_validated)
{
  MetrisParameters parameters;
  parameters.opt_smoo_tol = -1.0e-3;
  BOOST_CHECK_THROW(parameters.checkParameters(),MetrisExcept);

  parameters.opt_smoo_tol = 1.0;
  BOOST_CHECK_THROW(parameters.checkParameters(),MetrisExcept);

  parameters.opt_smoo_tol = 0.005;
  BOOST_CHECK_NO_THROW(parameters.checkParameters());
}

} // namespace Metris
