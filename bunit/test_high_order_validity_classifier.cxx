// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_high_order_validity_classifier

#include <boost/test/included/unit_test.hpp>

#include "Mesh/MeshBase.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "ho_constants.hxx"
#include "low_geo/validity.hxx"

#include <algorithm>
#include <cmath>

namespace Metris {

namespace {

enum class GeometryCase {
  Certified,
  HighOrderInvalid,
  PositiveButUncertified,
  DegenerateP1
};

template<int gdim>
void initialize_p2_element(MeshBase &mesh,
                           MetrisParameters &parameters,
                           GeometryCase geometry_case)
{
  constexpr int degree = 2;
  constexpr int node_count = getnnode(gdim,degree);
  constexpr auto ordering = ORDELT(gdim);

  mesh.idim = gdim;
  mesh.curdeg = degree;
  mesh.strdeg = degree;
  mesh.forceBasisFlag(FEBasis::Lagrange);
  mesh.param = &parameters;
  mesh.set_npoin(node_count);
  if constexpr(gdim == 2){
    mesh.set_nface(1);
  }else{
    mesh.set_nelem(1);
  }

  intAr2 &element_to_point = mesh.ent2poi(gdim);
  const double coupling = std::sqrt(0.75);
  for(int node = 0; node < node_count; node++){
    element_to_point(0,node) = node;
    const double u = ordering[degree][node][1]/(double)degree;
    const double v = ordering[degree][node][2]/(double)degree;
    const double w = gdim == 3
                   ? ordering[degree][node][3]/(double)degree : 0.0;

    double x = u;
    double y = v;
    double z = w;
    switch(geometry_case){
    case GeometryCase::Certified:
      break;
    case GeometryCase::HighOrderInvalid:
      if constexpr(gdim == 2){
        y = v - 2.0*u*v;
      }else{
        z = w - 2.0*u*w;
      }
      break;
    case GeometryCase::PositiveButUncertified:
      x = u + coupling*v*v;
      y = v + coupling*u*u;
      break;
    case GeometryCase::DegenerateP1:
      if constexpr(gdim == 2){
        y = 0.0;
      }else{
        z = 0.0;
      }
      break;
    }

    mesh.coord(node,0) = x;
    mesh.coord(node,1) = y;
    if constexpr(gdim == 3) mesh.coord(node,2) = z;
  }
}

void check_result_diagnostics_equal(const ElementValidityResult &lagrange,
                                    const ElementValidityResult &bezier)
{
  BOOST_CHECK(lagrange.status == bezier.status);
  BOOST_CHECK_EQUAL(lagrange.lower_bound_coefficient_index < 0,
                    bezier.lower_bound_coefficient_index < 0);
  BOOST_CHECK_EQUAL(lagrange.witness_sample_index,
                    bezier.witness_sample_index);
  if(std::isnan(lagrange.normalized_lower_bound)){
    BOOST_CHECK(std::isnan(bezier.normalized_lower_bound));
  }else{
    const double scale = std::max(
        {1.0,std::abs(lagrange.normalized_lower_bound),
             std::abs(bezier.normalized_lower_bound)});
    BOOST_CHECK_SMALL(lagrange.normalized_lower_bound
                      - bezier.normalized_lower_bound,2.e-12*scale);
  }
  if(std::isnan(lagrange.normalized_witness)){
    BOOST_CHECK(std::isnan(bezier.normalized_witness));
  }else{
    const double scale = std::max(
        {1.0,std::abs(lagrange.normalized_witness),
             std::abs(bezier.normalized_witness)});
    BOOST_CHECK_SMALL(lagrange.normalized_witness
                      - bezier.normalized_witness,2.e-12*scale);
  }
}

template<int gdim>
ElementValidityResult classify_in_both_bases(GeometryCase geometry_case,
                                             double orientation_tolerance,
                                             ElementValidityResult *bezier)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.jtol = 1.e-6;
  parameters.vtol = orientation_tolerance;

  MeshBase mesh;
  initialize_p2_element<gdim>(mesh,parameters,geometry_case);
  const ElementValidityResult lagrange
      = classify_element_validity<gdim,2>(mesh,0);
  mesh.setBasis(FEBasis::Bezier);
  *bezier = classify_element_validity<gdim,2>(mesh,0);
  check_result_diagnostics_equal(lagrange,*bezier);
  return lagrange;
}

template<int gdim>
void check_certified_case()
{
  ElementValidityResult bezier;
  const ElementValidityResult result = classify_in_both_bases<gdim>(
      GeometryCase::Certified,1.e-12,&bezier);
  BOOST_CHECK(result.is_certified());
  BOOST_CHECK(result.accepted_conservatively());
  BOOST_CHECK_CLOSE_FRACTION(result.normalized_lower_bound,1.0,2.e-12);
  BOOST_CHECK_GE(result.lower_bound_coefficient_index,0);
  BOOST_CHECK(std::isnan(result.normalized_witness));
  BOOST_CHECK_EQUAL(result.witness_sample_index,-1);
}

template<int gdim>
void check_invalid_case()
{
  ElementValidityResult bezier;
  const ElementValidityResult result = classify_in_both_bases<gdim>(
      GeometryCase::HighOrderInvalid,1.e-12,&bezier);
  BOOST_CHECK(result.is_invalid());
  BOOST_CHECK(!result.accepted_conservatively());
  BOOST_CHECK_LT(result.normalized_lower_bound,1.e-6);
  BOOST_CHECK_GE(result.lower_bound_coefficient_index,0);
  BOOST_CHECK_LT(result.normalized_witness,1.e-6);
  BOOST_CHECK_GE(result.witness_sample_index,0);
}

template<int gdim>
void check_positive_uncertified_case()
{
  ElementValidityResult bezier;
  const ElementValidityResult result = classify_in_both_bases<gdim>(
      GeometryCase::PositiveButUncertified,1.e-12,&bezier);
  BOOST_CHECK(result.is_uncertified());
  BOOST_CHECK(!result.accepted_conservatively());
  BOOST_CHECK_LT(result.normalized_lower_bound,1.e-6);
  BOOST_CHECK_GE(result.lower_bound_coefficient_index,0);
  BOOST_CHECK(std::isnan(result.normalized_witness));
  BOOST_CHECK_EQUAL(result.witness_sample_index,-1);
}

template<int gdim>
void check_p1_orientation_prerequisite()
{
  ElementValidityResult bezier;
  const ElementValidityResult result = classify_in_both_bases<gdim>(
      GeometryCase::Certified,1.0,&bezier);
  BOOST_CHECK(result.is_uncertified());
  BOOST_CHECK(!result.accepted_conservatively());
  BOOST_CHECK_CLOSE_FRACTION(result.normalized_lower_bound,1.0,2.e-12);
  BOOST_CHECK(std::isnan(result.normalized_witness));
  BOOST_CHECK_EQUAL(result.witness_sample_index,-1);
}

template<int gdim>
void check_degenerate_normalization()
{
  ElementValidityResult bezier;
  const ElementValidityResult result = classify_in_both_bases<gdim>(
      GeometryCase::DegenerateP1,1.e-12,&bezier);
  BOOST_CHECK(result.is_uncertified());
  BOOST_CHECK(!result.accepted_conservatively());
  BOOST_CHECK(std::isnan(result.normalized_lower_bound));
  BOOST_CHECK_EQUAL(result.lower_bound_coefficient_index,-1);
  BOOST_CHECK(std::isnan(result.normalized_witness));
  BOOST_CHECK_EQUAL(result.witness_sample_index,-1);
}

} // namespace

BOOST_AUTO_TEST_CASE(p2_triangle_validity_outcomes)
{
  check_certified_case<2>();
  check_invalid_case<2>();
  check_positive_uncertified_case<2>();
  check_p1_orientation_prerequisite<2>();
  check_degenerate_normalization<2>();
}

BOOST_AUTO_TEST_CASE(p2_tetrahedron_validity_outcomes)
{
  check_certified_case<3>();
  check_invalid_case<3>();
  check_positive_uncertified_case<3>();
  check_p1_orientation_prerequisite<3>();
  check_degenerate_normalization<3>();
}

} // namespace Metris
