// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_high_order_validity_contract

#include "common_setup.hxx"

#include <cmath>
#include <limits>

using namespace Metris;

BOOST_AUTO_TEST_CASE(explicit_status_and_conservative_policy)
{
  ElementValidityResult result;

  BOOST_CHECK(result.is_uncertified());
  BOOST_CHECK(!result.is_certified());
  BOOST_CHECK(!result.is_invalid());
  BOOST_CHECK(!result.accepted_conservatively());
  BOOST_CHECK(std::isnan(result.normalized_lower_bound));
  BOOST_CHECK(std::isnan(result.normalized_witness));
  BOOST_CHECK_EQUAL(result.lower_bound_coefficient_index, -1);
  BOOST_CHECK_EQUAL(result.witness_sample_index, -1);

  result.status = ElementValidityStatus::Certified;
  result.normalized_lower_bound = 2.0e-4;
  result.lower_bound_coefficient_index = 4;
  BOOST_CHECK(result.is_certified());
  BOOST_CHECK(result.accepted_conservatively());

  result.status = ElementValidityStatus::Invalid;
  result.normalized_witness = -3.0e-3;
  result.witness_sample_index = 7;
  BOOST_CHECK(result.is_invalid());
  BOOST_CHECK(!result.accepted_conservatively());

  result.status = ElementValidityStatus::Uncertified;
  result.normalized_witness = std::numeric_limits<double>::quiet_NaN();
  result.witness_sample_index = -1;
  BOOST_CHECK(result.is_uncertified());
  BOOST_CHECK(!result.accepted_conservatively());
}
