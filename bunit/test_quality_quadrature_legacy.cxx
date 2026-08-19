// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_quality_quadrature_legacy

#include <boost/test/included/unit_test.hpp>

#include <array>
#include <cmath>

namespace
{

// Test-local specification of the legacy integration scheme currently embedded
// in low_metqua.cxx and low_metqua_d.cxx. Step 3 of the quadrature plan will
// replace this generator with the shared production SimplexQuadrature view.
template <int tdim>
struct LegacyReferenceQuadrature
{
  static_assert(tdim == 2 || tdim == 3);

  static constexpr int nquad = tdim + 2;
  using BarycentricPoint = std::array<double, tdim + 1>;

  static constexpr double weight()
  {
    return 1.0 / static_cast<double>(nquad);
  }

  static BarycentricPoint point(int iquad)
  {
    BarycentricPoint bary{};

    if (iquad < tdim + 1)
    {
      bary[iquad] = 1.0;
    }
    else
    {
      for (double &lambda : bary)
      {
        lambda = 1.0 / static_cast<double>(tdim + 1);
      }
    }

    return bary;
  }
};

template <int tdim>
double integrate_barycentric_coordinate(int icoord)
{
  using Rule = LegacyReferenceQuadrature<tdim>;

  double integral = 0.0;
  for (int iquad = 0; iquad < Rule::nquad; ++iquad)
  {
    integral += Rule::weight() * Rule::point(iquad)[icoord];
  }
  return integral;
}

template <int tdim>
void check_legacy_rule()
{
  using Rule = LegacyReferenceQuadrature<tdim>;
  constexpr double tolerance = 1.0e-15;

  BOOST_TEST(Rule::nquad == tdim + 2);

  double weight_sum = 0.0;
  for (int iquad = 0; iquad < Rule::nquad; ++iquad)
  {
    const auto bary = Rule::point(iquad);
    const double weight = Rule::weight();

    BOOST_TEST(weight > 0.0);
    weight_sum += weight;

    double bary_sum = 0.0;
    for (double lambda : bary)
    {
      BOOST_TEST(lambda >= 0.0);
      BOOST_TEST(lambda <= 1.0);
      bary_sum += lambda;
    }
    BOOST_TEST(bary_sum == 1.0, boost::test_tools::tolerance(tolerance));

    if (iquad < tdim + 1)
    {
      for (int icoord = 0; icoord < tdim + 1; ++icoord)
      {
        const double expected = icoord == iquad ? 1.0 : 0.0;
        BOOST_TEST(bary[icoord] == expected);
      }
    }
    else
    {
      for (double lambda : bary)
      {
        BOOST_TEST(lambda == 1.0 / static_cast<double>(tdim + 1),
                   boost::test_tools::tolerance(tolerance));
      }
    }
  }

  // Normalized reference-simplex integral of the constant one.
  BOOST_TEST(weight_sum == 1.0, boost::test_tools::tolerance(tolerance));

  // Exactness for the barycentric affine basis.
  for (int icoord = 0; icoord < tdim + 1; ++icoord)
  {
    BOOST_TEST(integrate_barycentric_coordinate<tdim>(icoord)
                   == 1.0 / static_cast<double>(tdim + 1),
               boost::test_tools::tolerance(tolerance));
  }

  // One nontrivial affine function, independent of the basis checks above.
  double computed = 0.0;
  constexpr double constant = -0.75;
  for (int iquad = 0; iquad < Rule::nquad; ++iquad)
  {
    const auto bary = Rule::point(iquad);
    double value = constant;
    for (int icoord = 0; icoord < tdim + 1; ++icoord)
    {
      value += static_cast<double>(2 * icoord + 1) * bary[icoord];
    }
    computed += Rule::weight() * value;
  }

  double expected = constant;
  for (int icoord = 0; icoord < tdim + 1; ++icoord)
  {
    expected += static_cast<double>(2 * icoord + 1)
              / static_cast<double>(tdim + 1);
  }
  BOOST_TEST(computed == expected,
             boost::test_tools::tolerance(tolerance));
}

} // namespace

BOOST_AUTO_TEST_CASE(legacy_triangle_rule)
{
  check_legacy_rule<2>();
}

BOOST_AUTO_TEST_CASE(legacy_tetrahedron_rule)
{
  check_legacy_rule<3>();
}
