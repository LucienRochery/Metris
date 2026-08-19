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
void check_objective_order_zero_selects_vertex_barycenter()
{
  const SimplexQuadratureView<tdim> selected_rule
      = get_objective_quadrature<tdim>(0);
  const SimplexQuadratureView<tdim> historical_rule
      = get_vertex_barycenter_quadrature<tdim>();

  BOOST_TEST(selected_rule.size() == historical_rule.size());
  for(int iquad = 0; iquad < historical_rule.size(); iquad++)
  {
    const SimplexQuadraturePointView<tdim> selected_point
        = selected_rule[iquad];
    const SimplexQuadraturePointView<tdim> historical_point
        = historical_rule[iquad];
    BOOST_TEST(selected_point.weight == historical_point.weight);
    for(int ibary = 0; ibary < tdim + 1; ibary++)
    {
      BOOST_TEST(selected_point.bary[ibary]
                 == historical_point.bary[ibary]);
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
  check_objective_order_zero_selects_vertex_barycenter<2>();
  check_objective_order_zero_selects_vertex_barycenter<3>();
  BOOST_CHECK_THROW(
      get_objective_quadrature<2>(-1),Metris::MetrisExcept);
  BOOST_CHECK_THROW(
      get_objective_quadrature<3>(1),Metris::MetrisExcept);
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

  BOOST_TEST(triangle_degree_two.size() == 3);
  BOOST_TEST(triangle_degree_three.size() == 6);
  BOOST_TEST(tetrahedron_degree_two.size() == 4);
  BOOST_TEST(tetrahedron_degree_three.size() == 8);

  // Importing the tables does not activate them in the runtime selector.
  BOOST_CHECK_THROW(
      get_objective_quadrature<2>(2), Metris::MetrisExcept);
  BOOST_CHECK_THROW(
      get_objective_quadrature<3>(3), Metris::MetrisExcept);
}
