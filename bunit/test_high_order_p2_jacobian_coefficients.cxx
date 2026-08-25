// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_high_order_p2_jacobian_coefficients

#include <boost/test/included/unit_test.hpp>

#include "codegen_ccoef.hxx"
#include "codegen_lag2bez.hxx"
#include "ho_constants.hxx"
#include "linalg/det.hxx"
#include "low_eval.hxx"
#include "low_geo/ccoef.hxx"
#include "types.hxx"

#include <algorithm>
#include <array>
#include <cmath>

namespace Metris
{

namespace
{

constexpr int geometry_degree = 2;

void check_close(double expected, double computed, double tolerance = 2.e-12)
{
  const double scale = std::max({1.0,std::abs(expected),std::abs(computed)});
  BOOST_CHECK_SMALL(computed - expected,tolerance*scale);
}

template<int gdim>
double vertex_coordinate(int vertex, int component)
{
  if constexpr(gdim == 2){
    constexpr double vertices[3][2] = {
        { 0.10,-0.05},
        { 1.20, 0.10},
        {-0.05, 1.05}};
    return vertices[vertex][component];
  }else{
    constexpr double vertices[4][3] = {
        {0.05,0.02,0.03},
        {1.10,0.08,0.00},
        {0.12,1.00,0.13},
        {0.08,0.18,1.05}};
    return vertices[vertex][component];
  }
}

template<int gdim>
double curved_control_point_offset(int local_node, int component)
{
  if constexpr(gdim == 2){
    constexpr double offsets[6][2] = {
        {0.00, 0.00}, { 0.00,0.00}, {0.00,0.00},
        {0.04,-0.03}, {-0.02,0.05}, {0.03,0.06}};
    return offsets[local_node][component];
  }else{
    constexpr double offsets[10][3] = {
        { 0.000, 0.000, 0.000}, { 0.000, 0.000, 0.000},
        { 0.000, 0.000, 0.000}, { 0.000, 0.000, 0.000},
        { 0.020,-0.015, 0.010}, {-0.010, 0.025, 0.020},
        { 0.015, 0.010,-0.020}, {-0.020, 0.015, 0.025},
        { 0.010,-0.020, 0.015}, { 0.025, 0.005,-0.015}};
    return offsets[local_node][component];
  }
}

template<int gdim>
void initialize_p2_lagrange_element(bool curved,
                                    intAr2 &element_to_point,
                                    dblAr2 &coordinates)
{
  constexpr int nnode = getnnode(gdim,geometry_degree);
  constexpr auto ordering = ORDELT(gdim);

  for(int local_node = 0; local_node < nnode; local_node++){
    element_to_point(0,local_node) = local_node;
    for(int component = 0; component < gdim; component++){
      double coordinate = 0.0;
      for(int vertex = 0; vertex < gdim + 1; vertex++){
        const double barycentric_coordinate
            = ordering[geometry_degree][local_node][vertex]
             /(double)geometry_degree;
        coordinate += barycentric_coordinate
                    * vertex_coordinate<gdim>(vertex,component);
      }
      if(curved){
        coordinate += curved_control_point_offset<gdim>(local_node,component);
      }
      coordinates(local_node,component) = coordinate;
    }
  }
}

template<int gdim>
void check_off_lattice_point(const intAr2 &element_to_point,
                             const dblAr2 &lagrange_coordinates,
                             const dblAr2 &bezier_coordinates,
                             const double *lagrange_coefficients,
                             const double *bezier_coefficients,
                             const double *generated_coefficients,
                             const double *barycentric_coordinates)
{
  constexpr int jacobian_degree = gdim*(geometry_degree - 1);
  constexpr int ncoefficient = getnnode(gdim,jacobian_degree);
  constexpr auto evaluate_geometry
      = gdim == 2 ? eval2<gdim,geometry_degree>
                  : eval3<gdim,geometry_degree>;
  constexpr auto evaluate_determinant
      = gdim == 2 ? eval2<1,jacobian_degree>
                  : eval3<1,jacobian_degree>;

  double lagrange_point[gdim];
  double bezier_point[gdim];
  double lagrange_jacobian[gdim*gdim];
  double bezier_jacobian[gdim*gdim];
  evaluate_geometry(
      lagrange_coordinates,element_to_point[0],FEBasis::Lagrange,
      DifVar::Bary,DifVar::None,barycentric_coordinates,
      lagrange_point,lagrange_jacobian,nullptr);
  evaluate_geometry(
      bezier_coordinates,element_to_point[0],FEBasis::Bezier,
      DifVar::Bary,DifVar::None,barycentric_coordinates,
      bezier_point,bezier_jacobian,nullptr);

  for(int component = 0; component < gdim; component++){
    check_close(lagrange_point[component],bezier_point[component]);
  }
  for(int entry = 0; entry < gdim*gdim; entry++){
    check_close(lagrange_jacobian[entry],bezier_jacobian[entry]);
  }

  const double direct_determinant = detmat<gdim>(lagrange_jacobian);
  const double bezier_direct_determinant = detmat<gdim>(bezier_jacobian);
  check_close(direct_determinant,bezier_direct_determinant);

  int coefficient_nodes[ncoefficient];
  for(int coefficient = 0; coefficient < ncoefficient; coefficient++){
    coefficient_nodes[coefficient] = coefficient;
  }

  auto evaluate_coefficient_field = [&](const double *coefficients){
    double determinant;
    double coefficient_storage[ncoefficient];
    std::copy_n(coefficients,ncoefficient,coefficient_storage);
    dblAr2 coefficient_field(ncoefficient,1,coefficient_storage);
    evaluate_determinant(
        coefficient_field,coefficient_nodes,FEBasis::Bezier,
        DifVar::None,DifVar::None,barycentric_coordinates,
        &determinant,nullptr,nullptr);
    return determinant;
  };

  check_close(direct_determinant,
              evaluate_coefficient_field(lagrange_coefficients));
  check_close(direct_determinant,
              evaluate_coefficient_field(bezier_coefficients));
  check_close(direct_determinant,
              evaluate_coefficient_field(generated_coefficients));
}

template<int gdim>
void check_p2_element(bool curved)
{
  constexpr int nnode = getnnode(gdim,geometry_degree);
  constexpr int jacobian_degree = gdim*(geometry_degree - 1);
  constexpr int ncoefficient = getnnode(gdim,jacobian_degree);
  constexpr auto convert_lagrange_to_bezier
      = gdim == 2 ? lag2bez2<geometry_degree,gdim>
                  : lag2bez3<geometry_degree,gdim>;
  constexpr auto generate_coefficients
      = gdim == 2 ? ccoef_genbez2<geometry_degree>
                  : ccoef_genbez3<geometry_degree>;

  intAr2 element_to_point(1,nnode);
  dblAr2 lagrange_coordinates(nnode,gdim);
  dblAr2 bezier_coordinates(nnode,gdim);
  initialize_p2_lagrange_element<gdim>(
      curved,element_to_point,lagrange_coordinates);
  convert_lagrange_to_bezier(
      element_to_point[0],lagrange_coordinates,bezier_coordinates);

  double lagrange_coefficients[ncoefficient];
  double bezier_coefficients[ncoefficient];
  double generated_coefficients[ncoefficient];
  ccoef_eval<gdim,gdim,geometry_degree>(
      FEBasis::Lagrange,element_to_point,lagrange_coordinates,
      0,nullptr,lagrange_coefficients);
  ccoef_eval<gdim,gdim,geometry_degree>(
      FEBasis::Bezier,element_to_point,bezier_coordinates,
      0,nullptr,bezier_coefficients);
  generate_coefficients(
      element_to_point,bezier_coordinates,0,generated_coefficients);

  double minimum_coefficient = generated_coefficients[0];
  double maximum_coefficient = generated_coefficients[0];
  for(int coefficient = 0; coefficient < ncoefficient; coefficient++){
    check_close(generated_coefficients[coefficient],
                bezier_coefficients[coefficient]);
    check_close(generated_coefficients[coefficient],
                lagrange_coefficients[coefficient]);
    minimum_coefficient
        = std::min(minimum_coefficient,generated_coefficients[coefficient]);
    maximum_coefficient
        = std::max(maximum_coefficient,generated_coefficients[coefficient]);
  }
  BOOST_CHECK_GT(minimum_coefficient,0.0);
  if(curved){
    BOOST_CHECK_GT(maximum_coefficient - minimum_coefficient,1.e-3);
  }else{
    check_close(minimum_coefficient,maximum_coefficient);
  }

  if constexpr(gdim == 2){
    constexpr double off_lattice_points[][3] = {
        {0.17,0.31,0.52},
        {0.61,0.13,0.26},
        {0.23,0.69,0.08},
        {0.37,0.22,0.41}};
    for(const auto &barycentric_coordinates : off_lattice_points){
      check_off_lattice_point<gdim>(
          element_to_point,lagrange_coordinates,bezier_coordinates,
          lagrange_coefficients,bezier_coefficients,generated_coefficients,
          barycentric_coordinates);
    }
  }else{
    constexpr double off_lattice_points[][4] = {
        {0.11,0.22,0.29,0.38},
        {0.41,0.07,0.19,0.33},
        {0.09,0.53,0.16,0.22},
        {0.26,0.18,0.47,0.09},
        {0.31,0.24,0.12,0.33}};
    for(const auto &barycentric_coordinates : off_lattice_points){
      check_off_lattice_point<gdim>(
          element_to_point,lagrange_coordinates,bezier_coordinates,
          lagrange_coefficients,bezier_coefficients,generated_coefficients,
          barycentric_coordinates);
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(p2_triangle_coefficients_and_off_lattice_values)
{
  check_p2_element<2>(false);
  check_p2_element<2>(true);
}

BOOST_AUTO_TEST_CASE(p2_tetrahedron_coefficients_and_off_lattice_values)
{
  check_p2_element<3>(false);
  check_p2_element<3>(true);
}

} // namespace Metris
