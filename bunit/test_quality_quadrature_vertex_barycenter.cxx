// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_quality_quadrature_vertex_barycenter

#include <boost/test/included/unit_test.hpp>

#include "quality/simplex_quadrature.hxx"

using Metris::SimplexQuadraturePointView;
using Metris::SimplexQuadratureView;
using Metris::get_objective_quadrature;
using Metris::get_positive_simplex_quadrature;
using Metris::get_vertex_barycenter_quadrature;

namespace
{

int factorial(const int value)
{
  int result = 1;
  for (int factor = 2; factor <= value; factor++)
  {
    result *= factor;
  }
  return result;
}

double integer_power(const double base, const int exponent)
{
  double result = 1.0;
  for (int ipower = 0; ipower < exponent; ipower++)
  {
    result *= base;
  }
  return result;
}

template <int tdim>
double integrate_barycentric_coordinate(const SimplexQuadratureView<tdim> rule,
                                        int icoord)
{
  double integral = 0.0;
  for (int iquad = 0; iquad < rule.size(); iquad++)
  {
    const SimplexQuadraturePointView<tdim> point = rule[iquad];
    integral += point.weight * point.bary[icoord];
  }
  return integral;
}

template <int tdim>
void check_vertex_barycenter_rule()
{
  constexpr double tolerance = 1.0e-15;
  const SimplexQuadratureView<tdim> rule
      = get_vertex_barycenter_quadrature<tdim>();

  BOOST_TEST(rule.size() == tdim + 2);

  double weight_sum = 0.0;
  for (int iquad = 0; iquad < rule.size(); iquad++)
  {
    const SimplexQuadraturePointView<tdim> point = rule[iquad];

    BOOST_TEST(point.weight > 0.0);
    weight_sum += point.weight;

    double bary_sum = 0.0;
    for (int icoord = 0; icoord < tdim + 1; icoord++)
    {
      const double lambda = point.bary[icoord];
      BOOST_TEST(lambda >= 0.0);
      BOOST_TEST(lambda <= 1.0);
      bary_sum += lambda;
    }
    BOOST_TEST(bary_sum == 1.0, boost::test_tools::tolerance(tolerance));

    if (iquad < tdim + 1)
    {
      for (int icoord = 0; icoord < tdim + 1; icoord++)
      {
        const double expected = icoord == iquad ? 1.0 : 0.0;
        BOOST_TEST(point.bary[icoord] == expected);
      }
    }
    else
    {
      for (int icoord = 0; icoord < tdim + 1; icoord++)
      {
        BOOST_TEST(point.bary[icoord]
                       == 1.0 / static_cast<double>(tdim + 1),
                   boost::test_tools::tolerance(tolerance));
      }
    }
  }

  // Normalized reference-simplex integral of the constant one.
  BOOST_TEST(weight_sum == 1.0, boost::test_tools::tolerance(tolerance));

  // Exactness for the barycentric affine basis.
  for (int icoord = 0; icoord < tdim + 1; icoord++)
  {
    BOOST_TEST(integrate_barycentric_coordinate(rule, icoord)
                   == 1.0 / static_cast<double>(tdim + 1),
               boost::test_tools::tolerance(tolerance));
  }

  // One nontrivial affine function, independent of the basis checks above.
  double computed = 0.0;
  constexpr double constant = -0.75;
  for (int iquad = 0; iquad < rule.size(); iquad++)
  {
    const SimplexQuadraturePointView<tdim> point = rule[iquad];
    double value = constant;
    for (int icoord = 0; icoord < tdim + 1; icoord++)
    {
      value += static_cast<double>(2 * icoord + 1) * point.bary[icoord];
    }
    computed += point.weight * value;
  }

  double expected = constant;
  for (int icoord = 0; icoord < tdim + 1; icoord++)
  {
    expected += static_cast<double>(2 * icoord + 1)
              / static_cast<double>(tdim + 1);
  }
  BOOST_TEST(computed == expected,
             boost::test_tools::tolerance(tolerance));
}

template <int tdim>
void check_same_rule(const SimplexQuadratureView<tdim> selected_rule,
                     const SimplexQuadratureView<tdim> expected_rule)
{
  BOOST_TEST(selected_rule.size() == expected_rule.size());
  for(int iquad = 0; iquad < expected_rule.size(); iquad++)
  {
    const SimplexQuadraturePointView<tdim> selected_point
        = selected_rule[iquad];
    const SimplexQuadraturePointView<tdim> expected_point
        = expected_rule[iquad];
    BOOST_TEST(selected_point.weight == expected_point.weight);
    for(int ibary = 0; ibary < tdim + 1; ibary++)
    {
      BOOST_TEST(selected_point.bary[ibary]
                 == expected_point.bary[ibary]);
    }
  }
}

template <int tdim>
void check_objective_quadrature_selection()
{
  if constexpr(tdim == 2)
  {
    check_same_rule(
        get_objective_quadrature<tdim>(-1),
        get_positive_simplex_quadrature<tdim, 4>());
  }
  else
  {
    static_assert(tdim == 3);
    check_same_rule(
        get_objective_quadrature<tdim>(-1),
        get_positive_simplex_quadrature<tdim, 3>());
  }
  check_same_rule(
      get_objective_quadrature<tdim>(0),
      get_vertex_barycenter_quadrature<tdim>());
  check_same_rule(
      get_objective_quadrature<tdim>(2),
      get_positive_simplex_quadrature<tdim, 2>());
  check_same_rule(
      get_objective_quadrature<tdim>(3),
      get_positive_simplex_quadrature<tdim, 3>());
  check_same_rule(
      get_objective_quadrature<tdim>(4),
      get_positive_simplex_quadrature<tdim, 4>());
  check_same_rule(
      get_objective_quadrature<tdim>(5),
      get_positive_simplex_quadrature<tdim, 5>());
}

template <int tdim>
void check_positive_rule_geometry(const SimplexQuadratureView<tdim> rule)
{
  constexpr double tolerance = 1.0e-14;
  double weight_sum = 0.0;

  for (int iquad = 0; iquad < rule.size(); iquad++)
  {
    const SimplexQuadraturePointView<tdim> point = rule[iquad];
    BOOST_TEST(point.weight > 0.0);
    weight_sum += point.weight;

    double bary_sum = 0.0;
    for (int ibary = 0; ibary < tdim + 1; ibary++)
    {
      BOOST_TEST(point.bary[ibary] >= 0.0);
      BOOST_TEST(point.bary[ibary] <= 1.0);
      bary_sum += point.bary[ibary];
    }
    BOOST_TEST(bary_sum == 1.0,
               boost::test_tools::tolerance(tolerance));
  }

  BOOST_TEST(weight_sum == 1.0,
             boost::test_tools::tolerance(tolerance));
}

double integrate_triangle_monomial(const SimplexQuadratureView<2> rule,
                                   const int x_power,
                                   const int y_power)
{
  double integral = 0.0;
  for (int iquad = 0; iquad < rule.size(); iquad++)
  {
    const SimplexQuadraturePointView<2> point = rule[iquad];
    integral += point.weight
              * integer_power(point.bary[1], x_power)
              * integer_power(point.bary[2], y_power);
  }
  return integral;
}

double exact_normalized_triangle_monomial(const int x_power,
                                          const int y_power)
{
  const int total_degree = x_power + y_power;
  return 2.0 * static_cast<double>(factorial(x_power))
             * static_cast<double>(factorial(y_power))
       / static_cast<double>(factorial(total_degree + 2));
}

template <int degree>
void check_positive_triangle_rule()
{
  constexpr double tolerance = 1.0e-13;
  const SimplexQuadratureView<2> rule
      = get_positive_simplex_quadrature<2, degree>();
  check_positive_rule_geometry(rule);

  for (int x_power = 0; x_power <= degree; x_power++)
  {
    for (int y_power = 0; y_power <= degree - x_power; y_power++)
    {
      BOOST_TEST_CONTEXT("x power = " << x_power
                         << ", y power = " << y_power)
      {
        const double computed
            = integrate_triangle_monomial(rule, x_power, y_power);
        const double expected
            = exact_normalized_triangle_monomial(x_power, y_power);
        BOOST_TEST(computed == expected,
                   boost::test_tools::tolerance(tolerance));
      }
    }
  }
}

double integrate_tetrahedron_monomial(const SimplexQuadratureView<3> rule,
                                      const int x_power,
                                      const int y_power,
                                      const int z_power)
{
  double integral = 0.0;
  for (int iquad = 0; iquad < rule.size(); iquad++)
  {
    const SimplexQuadraturePointView<3> point = rule[iquad];
    integral += point.weight
              * integer_power(point.bary[1], x_power)
              * integer_power(point.bary[2], y_power)
              * integer_power(point.bary[3], z_power);
  }
  return integral;
}

double exact_normalized_tetrahedron_monomial(const int x_power,
                                             const int y_power,
                                             const int z_power)
{
  const int total_degree = x_power + y_power + z_power;
  return 6.0 * static_cast<double>(factorial(x_power))
             * static_cast<double>(factorial(y_power))
             * static_cast<double>(factorial(z_power))
       / static_cast<double>(factorial(total_degree + 3));
}

template <int degree>
void check_positive_tetrahedron_rule()
{
  constexpr double tolerance = 1.0e-13;
  const SimplexQuadratureView<3> rule
      = get_positive_simplex_quadrature<3, degree>();
  check_positive_rule_geometry(rule);

  for (int x_power = 0; x_power <= degree; x_power++)
  {
    for (int y_power = 0; y_power <= degree - x_power; y_power++)
    {
      for (int z_power = 0;
           z_power <= degree - x_power - y_power;
           z_power++)
      {
        BOOST_TEST_CONTEXT("x power = " << x_power
                           << ", y power = " << y_power
                           << ", z power = " << z_power)
        {
          const double computed = integrate_tetrahedron_monomial(
              rule, x_power, y_power, z_power);
          const double expected = exact_normalized_tetrahedron_monomial(
              x_power, y_power, z_power);
          BOOST_TEST(computed == expected,
                     boost::test_tools::tolerance(tolerance));
        }
      }
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(vertex_barycenter_triangle_rule)
{
  check_vertex_barycenter_rule<2>();
}

BOOST_AUTO_TEST_CASE(vertex_barycenter_tetrahedron_rule)
{
  check_vertex_barycenter_rule<3>();
}

BOOST_AUTO_TEST_CASE(objective_quadrature_order_selection)
{
  check_objective_quadrature_selection<2>();
  check_objective_quadrature_selection<3>();
  BOOST_CHECK_THROW(
      get_objective_quadrature<2>(-2),Metris::MetrisExcept);
  BOOST_CHECK_THROW(
      get_objective_quadrature<3>(1),Metris::MetrisExcept);
  BOOST_CHECK_THROW(
      get_objective_quadrature<2>(6),Metris::MetrisExcept);
}

BOOST_AUTO_TEST_CASE(imported_positive_quadrature_tables_are_available)
{
  const SimplexQuadratureView<2> triangle_degree_two
      = get_positive_simplex_quadrature<2, 2>();
  const SimplexQuadratureView<2> triangle_degree_three
      = get_positive_simplex_quadrature<2, 3>();
  const SimplexQuadratureView<3> tetrahedron_degree_two
      = get_positive_simplex_quadrature<3, 2>();
  const SimplexQuadratureView<3> tetrahedron_degree_three
      = get_positive_simplex_quadrature<3, 3>();
  const SimplexQuadratureView<2> triangle_degree_four
      = get_positive_simplex_quadrature<2, 4>();
  const SimplexQuadratureView<3> tetrahedron_degree_four
      = get_positive_simplex_quadrature<3, 4>();
  const SimplexQuadratureView<2> triangle_degree_five
      = get_positive_simplex_quadrature<2, 5>();
  const SimplexQuadratureView<3> tetrahedron_degree_five
      = get_positive_simplex_quadrature<3, 5>();

  BOOST_TEST(triangle_degree_two.size() == 3);
  BOOST_TEST(triangle_degree_three.size() == 6);
  BOOST_TEST(tetrahedron_degree_two.size() == 4);
  BOOST_TEST(tetrahedron_degree_three.size() == 8);
  BOOST_TEST(triangle_degree_four.size() == 6);
  BOOST_TEST(tetrahedron_degree_four.size() == 14);
  BOOST_TEST(triangle_degree_five.size() == 7);
  BOOST_TEST(tetrahedron_degree_five.size() == 14);
  check_same_rule(tetrahedron_degree_four,tetrahedron_degree_five);
}

BOOST_AUTO_TEST_CASE(positive_triangle_degree_two_rule_contract)
{
  check_positive_triangle_rule<2>();
}

BOOST_AUTO_TEST_CASE(positive_triangle_degree_three_rule_contract)
{
  check_positive_triangle_rule<3>();
}

BOOST_AUTO_TEST_CASE(positive_tetrahedron_degree_two_rule_contract)
{
  check_positive_tetrahedron_rule<2>();
}

BOOST_AUTO_TEST_CASE(positive_tetrahedron_degree_three_rule_contract)
{
  check_positive_tetrahedron_rule<3>();
}

BOOST_AUTO_TEST_CASE(positive_triangle_degree_four_rule_contract)
{
  check_positive_triangle_rule<4>();
}

BOOST_AUTO_TEST_CASE(positive_tetrahedron_degree_four_rule_contract)
{
  check_positive_tetrahedron_rule<4>();
}

BOOST_AUTO_TEST_CASE(positive_triangle_degree_five_rule_contract)
{
  check_positive_triangle_rule<5>();
}

BOOST_AUTO_TEST_CASE(positive_tetrahedron_degree_five_rule_contract)
{
  check_positive_tetrahedron_rule<5>();
}
