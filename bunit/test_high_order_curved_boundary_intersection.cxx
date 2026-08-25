// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_high_order_curved_boundary_intersection

#include "common_setup.hxx"

#include "API/MetrisAPI.hxx"
#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "MetrisRunner/MetrisRunner.hxx"
#include "low_geo/validity.hxx"
#include "low_topo.hxx"
#include "smoothing/msh_smooball.hxx"

#include <filesystem>
#include <memory>

using namespace Metris;

namespace
{

std::unique_ptr<MetrisRunner> make_two_component_planar_p2_mesh()
{
  constexpr double coordinates[6][2] = {
      {0.0, 0.0},{2.0, 0.0},{1.0, 1.0},
      {0.3,-0.2},{1.0,-0.6},{1.7,-0.2}};
  constexpr int edges[6][2] = {
      {0,1},{1,2},{2,0},
      {3,4},{4,5},{5,3}};
  constexpr int faces[2][3] = {{0,1,2},{3,4,5}};

  MetrisAPI input(
      2,1,6,0,0,6,false,6,2,0,
      FEBasis::Lagrange,FEBasis::Undefined,MetSpace::Exp);
  for(int point = 0; point < 6; point++){
    input.setCoord(point,coordinates[point]);
    input.setCorner(point,point);
  }
  for(int edge = 0; edge < 6; edge++){
    input.setElement(1,edge,edges[edge],edge < 3 ? 0 : 1);
  }
  input.setElement(2,0,faces[0],0);
  input.setElement(2,1,faces[1],1);

  MetrisParameters parameters;
  parameters.setAnalyticalMetric(1);
  parameters.usrTarDeg = 2;
  parameters.adp_niter = 0;
  parameters.opt_niter = 0;
  parameters.iverb = 0;
  parameters.outmPrefix
      = (std::filesystem::temp_directory_path()
         / "metris_high_order_curved_boundary_intersection").string()
        + "/";
  auto runner = std::make_unique<MetrisRunner>(
      &input,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner->degElevate(),1);
  return runner;
}

} // namespace

BOOST_AUTO_TEST_CASE(
    p2_curved_boundary_crossing_is_detected_despite_local_validity)
{
  auto runner = make_two_component_planar_p2_mesh();
  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner->msh_g);
  mesh.setBasis(FEBasis::Lagrange);

  const int curvedEdge = getedgglo(mesh,0,1);
  BOOST_REQUIRE_GE(curvedEdge,0);
  const int controlPoint = mesh.edg2poi(curvedEdge,2);
  BOOST_REQUIRE_GE(controlPoint,0);

  // Straight, disjoint components are initially accepted. Contacts at the
  // two endpoints shared with adjacent edges of the first triangle are legal.
  BOOST_CHECK(
      planarP2BoundaryIsIntersectionFreeAroundPoint(mesh,controlPoint));

  // Bow the first triangle's bottom edge into the separate lower component.
  // Both polynomial triangles retain positive Bernstein Jacobian bounds, so
  // element-local validity alone cannot identify this global boundary contact.
  mesh.coord(controlPoint,1) = -0.3;
  BOOST_REQUIRE((
      classify_element_validity<2,2>(mesh,0).accepted_conservatively()));
  BOOST_REQUIRE((
      classify_element_validity<2,2>(mesh,1).accepted_conservatively()));
  BOOST_CHECK(!planarP2BoundaryIsIntersectionFreeAroundPoint(
      mesh,controlPoint));

  mesh.coord(controlPoint,1) = 0.0;
  BOOST_CHECK(
      planarP2BoundaryIsIntersectionFreeAroundPoint(mesh,controlPoint));
}
