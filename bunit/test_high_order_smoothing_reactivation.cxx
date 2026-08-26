// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_high_order_smoothing_reactivation

#include "common_setup.hxx"

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "MetrisRunner/MetrisRunner.hxx"
#include "aux_badEntHandler.hxx"
#include "low_geo/validity.hxx"
#include "quality/low_metqua.hxx"
#include "smoothing/msh_smooball.hxx"

#include <cmath>
#include <filesystem>
#include <memory>
#include <vector>

using namespace Metris;

namespace
{

std::unique_ptr<MetrisRunner> make_two_triangle_runner(
    int degree,
    const char* case_name)
{
  MetrisParameters parameters;
  parameters.setAnalyticalMetric(1);
  parameters.setMeshIn(
      METRIS_ROOT_DIR "/examples/2D/misc/2tri2D.mesh");
  parameters.usrTarDeg = degree;
  parameters.adp_niter = 0;
  parameters.opt_niter = 0;
  parameters.opt_smoo_niter = 1;
  parameters.iflag1 = 1;
  parameters.iverb = 0;
  parameters.outmPrefix
      = (std::filesystem::temp_directory_path()
         / "metris_high_order_smoothing_reactivation"
         / case_name).string() + "/";

  auto runner = std::make_unique<MetrisRunner>(
      nullptr,nullptr,parameters);
  if(degree == 2) BOOST_REQUIRE_EQUAL(runner->degElevate(),1);
  return runner;
}

std::unique_ptr<MetrisRunner> make_cube_runner(const char* case_name)
{
  MetrisParameters parameters;
  parameters.setAnalyticalMetric(1);
  parameters.setMeshIn(
      METRIS_ROOT_DIR "/examples/3D/cube/cube.meshb");
  parameters.usrTarDeg = 2;
  parameters.adp_niter = 0;
  parameters.opt_niter = 0;
  parameters.iverb = 0;
  parameters.outmPrefix
      = (std::filesystem::temp_directory_path()
         / "metris_high_order_smoothing_reactivation"
         / case_name).string() + "/";

  auto runner = std::make_unique<MetrisRunner>(
      nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner->degElevate(),1);
  return runner;
}

bool point_belongs_to_face(const MeshBase& mesh,
                           int face,
                           int degree,
                           int point)
{
  const int node_count = getnnode(2,degree);
  for(int local = 0; local < node_count; local++){
    if(mesh.fac2poi(face,local) == point) return true;
  }
  return false;
}

int first_point_outside_face(const MeshBase& mesh,
                             int face,
                             int degree)
{
  for(int point = 0; point < mesh.npoin; point++){
    if(mesh.isdeadpoint(point)) continue;
    if(!point_belongs_to_face(mesh,face,degree,point)) return point;
  }
  return -1;
}

void set_all_live_point_tags(MeshBase& mesh, int ithread, int value)
{
  for(int point = 0; point < mesh.npoin; point++){
    if(!mesh.isdeadpoint(point)) mesh.poi2tag(ithread,point) = value;
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(p1_reactivation_retains_corner_only_contract)
{
  auto runner = make_two_triangle_runner(1,"p1_tags");
  auto& mesh
      = static_cast<Mesh<MetricFieldAnalytical>&>(*runner->msh_g);
  BOOST_REQUIRE_EQUAL(mesh.curdeg,1);

  constexpr int tag_thread = 1;
  mesh.tag[tag_thread]++;
  const int current_tag = mesh.tag[tag_thread];
  set_all_live_point_tags(mesh,tag_thread,current_tag);

  intAr1 region(1);
  region.stack(0);
  reactivateSmoothingRegionGeometry(
      mesh,2,1,region,tag_thread);

  for(int local = 0; local < 3; local++){
    BOOST_CHECK_EQUAL(
        mesh.poi2tag(tag_thread,mesh.fac2poi(0,local)),
        current_tag - 1);
  }
  const int unaffected = first_point_outside_face(mesh,0,1);
  BOOST_REQUIRE_GE(unaffected,0);
  BOOST_CHECK_EQUAL(mesh.poi2tag(tag_thread,unaffected),current_tag);
}

BOOST_AUTO_TEST_CASE(p2_reactivation_includes_every_edge_control_point)
{
  auto runner = make_two_triangle_runner(2,"p2_tags");
  auto& mesh
      = static_cast<Mesh<MetricFieldAnalytical>&>(*runner->msh_g);
  BOOST_REQUIRE_EQUAL(mesh.curdeg,2);

  constexpr int tag_thread = 1;
  mesh.tag[tag_thread]++;
  const int current_tag = mesh.tag[tag_thread];
  set_all_live_point_tags(mesh,tag_thread,current_tag);

  intAr1 region(1);
  region.stack(0);
  reactivateSmoothingRegionGeometry(
      mesh,2,2,region,tag_thread);

  for(int local = 0; local < 6; local++){
    const int point = mesh.fac2poi(0,local);
    BOOST_CHECK_EQUAL(
        mesh.poi2tag(tag_thread,point),current_tag - 1);
  }
  const int unaffected = first_point_outside_face(mesh,0,2);
  BOOST_REQUIRE_GE(unaffected,0);
  BOOST_CHECK_EQUAL(mesh.poi2tag(tag_thread,unaffected),current_tag);
}

BOOST_AUTO_TEST_CASE(p2_tetrahedral_reactivation_includes_all_geometry_nodes)
{
  auto runner = make_cube_runner("p2_tetrahedral_tags");
  auto& mesh
      = static_cast<Mesh<MetricFieldAnalytical>&>(*runner->msh_g);
  BOOST_REQUIRE_EQUAL(mesh.curdeg,2);

  constexpr int tag_thread = 1;
  mesh.tag[tag_thread]++;
  const int current_tag = mesh.tag[tag_thread];
  set_all_live_point_tags(mesh,tag_thread,current_tag);

  intAr1 region(1);
  region.stack(0);
  reactivateSmoothingRegionGeometry(
      mesh,3,2,region,tag_thread);

  for(int local = 0; local < 10; local++){
    BOOST_CHECK_EQUAL(
        mesh.poi2tag(tag_thread,mesh.tet2poi(0,local)),
        current_tag - 1);
  }
}

BOOST_AUTO_TEST_CASE(targeted_p2_control_move_reactivates_geometry_and_element)
{
  auto runner = make_two_triangle_runner(2,"targeted_p2");
  auto& mesh
      = static_cast<Mesh<MetricFieldAnalytical>&>(*runner->msh_g);
  mesh.met.setSpace(MetSpace::Exp);
  mesh.setBasis(FEBasis::Lagrange);

  constexpr int seed_face = 0;
  int interior_edge = -1;
  for(int edge = 0; edge < 3; edge++){
    if(mesh.fac2fac(seed_face,edge) >= 0){
      interior_edge = edge;
      break;
    }
  }
  BOOST_REQUIRE_GE(interior_edge,0);
  const int control_point
      = mesh.fac2poi(seed_face,3 + interior_edge);

  int boundary_edge = -1;
  for(int edge = 0; edge < 3; edge++){
    if(mesh.fac2fac(seed_face,edge) < 0){
      boundary_edge = edge;
      break;
    }
  }
  BOOST_REQUIRE_GE(boundary_edge,0);
  const int fixed_boundary_control
      = mesh.fac2poi(seed_face,3 + boundary_edge);
  const double fixed_boundary_coordinates[2] = {
      mesh.coord(fixed_boundary_control,0),
      mesh.coord(fixed_boundary_control,1)};

  mesh.coord(control_point,0) += 0.12;
  mesh.coord(control_point,1) += 0.12;
  const double perturbed[2] = {
      mesh.coord(control_point,0),mesh.coord(control_point,1)};

  constexpr int tag_thread = 1;
  const int next_tag = mesh.tag[tag_thread] + 1;
  set_all_live_point_tags(mesh,tag_thread,next_tag);
  mesh.poi2tag(tag_thread,control_point) = mesh.tag[tag_thread];

  BadEntHandler handler(2,100.,1.e-8);
  const double operations = smoothElement_Ball(
      mesh,seed_face,handler,QuaFun::SizeShape,tag_thread,2);
  BOOST_REQUIRE_GT(operations,0.);
  BOOST_CHECK_GT(
      std::hypot(mesh.coord(control_point,0) - perturbed[0],
                 mesh.coord(control_point,1) - perturbed[1]),
      1.e-10);

  const int active_tag = mesh.tag[tag_thread] - 1;
  for(int face = 0; face < mesh.nface; face++){
    for(int local = 0; local < 6; local++){
      BOOST_CHECK_EQUAL(
          mesh.poi2tag(tag_thread,mesh.fac2poi(face,local)),
          active_tag);
    }
  }
  BOOST_CHECK_EQUAL(
      mesh.coord(fixed_boundary_control,0),fixed_boundary_coordinates[0]);
  BOOST_CHECK_EQUAL(
      mesh.coord(fixed_boundary_control,1),fixed_boundary_coordinates[1]);

  BOOST_REQUIRE_EQUAL(handler.affectedEnttsAlive.size(),2);
  const auto objective
      = get_quafun<MetricFieldAnalytical,2,2>(QuaFun::SizeShape);
  for(int face = 0; face < mesh.nface; face++){
    const auto affected = handler.affectedEnttsAlive.find(face);
    BOOST_REQUIRE(affected != handler.affectedEnttsAlive.end());
    const double independently_recomputed = objective(
        mesh,AsDeg::Pk,AsDeg::Pk,face,1.0);
    BOOST_CHECK_CLOSE_FRACTION(
        affected->second,independently_recomputed,5.e-14);
    BOOST_CHECK((classify_element_validity<2,2>(
        mesh,face).is_certified()));
  }
}
