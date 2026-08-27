// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_high_order_classic_smoothing

#include "common_setup.hxx"

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "MetrisRunner/MetrisRunner.hxx"
#include "low_geo/validity.hxx"
#include "quality/low_metqua.hxx"
#include "smoothing/low_smooballdiff.hxx"
#include "smoothing/msh_smooball.hxx"
#include "smoothing/msh_smoolen.hxx"

#include <cmath>
#include <filesystem>
#include <memory>

using namespace Metris;

namespace
{

std::unique_ptr<MetrisRunner> make_linear_two_triangle_runner()
{
  MetrisParameters parameters;
  parameters.setAnalyticalMetric(1);
  parameters.setMeshIn(METRIS_ROOT_DIR "/examples/2D/misc/2tri2D.mesh");
  parameters.adp_niter = 0;
  parameters.opt_niter = 0;
  parameters.opt_smoo_niter = 1;
  parameters.opt_power = -1;
  parameters.iverb = 0;
  parameters.outmPrefix
      = (std::filesystem::temp_directory_path()
         / "metris_linear_classic_length_smoothing").string() + "/";

  return std::make_unique<MetrisRunner>(nullptr,nullptr,parameters);
}

std::unique_ptr<MetrisRunner> make_elevated_two_triangle_runner()
{
  MetrisParameters parameters;
  parameters.setAnalyticalMetric(1);
  parameters.setMeshIn(METRIS_ROOT_DIR "/examples/2D/misc/2tri2D.mesh");
  parameters.usrTarDeg = 2;
  parameters.adp_niter = 0;
  parameters.opt_niter = 0;
  parameters.iverb = 0;
  parameters.outmPrefix
      = (std::filesystem::temp_directory_path()
         / "metris_high_order_classic_smoothing").string() + "/";

  auto runner = std::make_unique<MetrisRunner>(
      nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner->degElevate(),1);
  return runner;
}

std::unique_ptr<MetrisRunner> make_elevated_square_cad_runner()
{
  MetrisParameters parameters;
  parameters.setAnalyticalMetric(1);
  parameters.setMeshIn(
      METRIS_ROOT_DIR "/examples/2D/square/square.meshb");
  parameters.setCAD(
      METRIS_ROOT_DIR "/examples/2D/square/square.egads");
  parameters.usrTarDeg = 2;
  parameters.adp_niter = 0;
  parameters.opt_niter = 0;
  parameters.iverb = 0;
  parameters.outmPrefix
      = (std::filesystem::temp_directory_path()
         / "metris_high_order_classic_smoothing_cad").string() + "/";

  auto runner = std::make_unique<MetrisRunner>(
      nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner->degElevate(),1);
  return runner;
}

std::unique_ptr<MetrisRunner> make_elevated_cube_runner()
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
         / "metris_high_order_classic_smoothing_3d").string() + "/";

  auto runner = std::make_unique<MetrisRunner>(
      nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner->degElevate(),1);
  return runner;
}

bool region_contains(const intAr1 &region, int entity)
{
  for(const int candidate : region){
    if(candidate == entity) return true;
  }
  return false;
}

template<int gdim>
void check_p2_region_is_certified(const MeshBase &mesh,
                                  const intAr1 &region)
{
  for(const int element : region){
    const ElementValidityResult validity
        = classify_element_validity<gdim,2>(mesh,element);
    BOOST_CHECK_MESSAGE(
        validity.is_certified(),
        "element " << element << " was not certified after smoothing");
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(classic_length_smoothing_remains_p1_only)
{
  auto linear_runner = make_linear_two_triangle_runner();
  auto &linear_mesh
      = static_cast<Mesh<MetricFieldAnalytical>&>(*linear_runner->msh_g);
  BOOST_REQUIRE_EQUAL(linear_mesh.curdeg,1);
  const double linear_stat = smoothMeshLength(linear_mesh,2,1,2);
  BOOST_CHECK(std::isfinite(linear_stat));
  BOOST_CHECK_GE(linear_stat,0.0);

  auto quadratic_runner = make_elevated_two_triangle_runner();
  auto &quadratic_mesh
      = static_cast<Mesh<MetricFieldAnalytical>&>(*quadratic_runner->msh_g);
  BOOST_REQUIRE_EQUAL(quadratic_mesh.curdeg,2);

  dblAr2 coordinates_before(quadratic_mesh.npoin,quadratic_mesh.idim);
  for(int ipoin = 0; ipoin < quadratic_mesh.npoin; ipoin++){
    for(int idim = 0; idim < quadratic_mesh.idim; idim++){
      coordinates_before(ipoin,idim) = quadratic_mesh.coord(ipoin,idim);
    }
  }

  BOOST_CHECK_EQUAL(smoothMeshLength(quadratic_mesh,2,1,2),0.0);
  for(int ipoin = 0; ipoin < quadratic_mesh.npoin; ipoin++){
    for(int idim = 0; idim < quadratic_mesh.idim; idim++){
      BOOST_CHECK_EQUAL(
          quadratic_mesh.coord(ipoin,idim),coordinates_before(ipoin,idim));
    }
  }
}

BOOST_AUTO_TEST_CASE(scaled_analytical_metric_remains_exact_after_degree_elevation)
{
  MetrisParameters parameters;
  parameters.setAnalyticalMetric(2);
  parameters.setMetricScale(0.164973429592);
  parameters.setMeshIn(
      METRIS_ROOT_DIR "/examples/2D/square/square.meshb");
  parameters.setCAD(
      METRIS_ROOT_DIR "/examples/2D/square/square.egads");
  parameters.usrTarDeg = 2;
  parameters.adp_niter = 0;
  parameters.opt_niter = 0;
  parameters.iverb = 0;
  parameters.outmPrefix
      = (std::filesystem::temp_directory_path()
         / "metris_scaled_analytical_degree_elevation").string() + "/";

  MetrisRunner runner(nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner.degElevate(),1);
  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);

  constexpr int metric_components = 3;
  for(int point = 0; point < mesh.npoin; point++){
    if(mesh.isdeadpoint(point)) continue;
    double expected[metric_components];
    mesh.met.getMetPhys(DifVar::None,mesh.met.getSpace(),
                        mesh.coord[point],expected,nullptr);
    for(int component = 0; component < metric_components; component++){
      BOOST_CHECK_SMALL(
          mesh.met(point,component)-expected[component],
          2.e-13*std::max(1.0,std::abs(expected[component])));
    }
  }
}

BOOST_AUTO_TEST_CASE(p2_edge_control_point_regions_in_two_triangles)
{
  auto runner = make_elevated_two_triangle_runner();
  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner->msh_g);
  BOOST_REQUIRE_EQUAL(mesh.curdeg,2);
  BOOST_REQUIRE_EQUAL(mesh.nface,2);

  int boundary_edge = -1;
  int interior_edge = -1;
  for(int local_edge = 0; local_edge < 3; local_edge++){
    if(mesh.fac2fac(0,local_edge) < 0){
      if(boundary_edge < 0) boundary_edge = local_edge;
    }else{
      interior_edge = local_edge;
    }
  }
  BOOST_REQUIRE_GE(boundary_edge,0);
  BOOST_REQUIRE_GE(interior_edge,0);

  intAr1 region(4);
  buildEdgeControlPointSmoothingRegion(
      mesh,2,0,boundary_edge,region);
  BOOST_REQUIRE_EQUAL(region.get_n(),1);
  BOOST_CHECK_EQUAL(region[0],0);

  const int neighbor = mesh.fac2fac(0,interior_edge);
  BOOST_REQUIRE_GE(neighbor,0);
  buildEdgeControlPointSmoothingRegion(
      mesh,2,0,interior_edge,region);
  BOOST_REQUIRE_EQUAL(region.get_n(),2);
  BOOST_CHECK(region_contains(region,0));
  BOOST_CHECK(region_contains(region,neighbor));

  const int boundary_control_point
      = mesh.fac2poi(0,3 + boundary_edge);
  const int interior_control_point
      = mesh.fac2poi(0,3 + interior_edge);
  BOOST_CHECK_GE(mesh.getverfac<2>(0,boundary_control_point),3);
  BOOST_CHECK_GE(mesh.getverfac<2>(0,interior_control_point),3);
  BOOST_CHECK_GE(
      mesh.getverfac<2>(neighbor,interior_control_point),3);
}

BOOST_AUTO_TEST_CASE(production_smoothing_objective_is_stepdistance_for_p2)
{
  BOOST_CHECK(
      productionSmoothingObjective(2,1) == DefaultQualityFunction);
  BOOST_CHECK(
      productionSmoothingObjective(2,2) == QuaFun::StepDistance);
  BOOST_CHECK(
      productionSmoothingObjective(3,2) == DefaultQualityFunction);
}

BOOST_AUTO_TEST_CASE(legacy_worst_element_guard_excludes_objective_paths)
{
  BOOST_CHECK(!isObjectiveDrivenSmoothing(QuaFun::Distortion));
  BOOST_CHECK(!isObjectiveDrivenSmoothing(QuaFun::Unit));
  BOOST_CHECK(isObjectiveDrivenSmoothing(QuaFun::SizeShape));
  BOOST_CHECK(isObjectiveDrivenSmoothing(QuaFun::StepDistance));
}

BOOST_AUTO_TEST_CASE(sizeshape_p2_interior_edge_control_point_is_smoothed)
{
  auto runner = make_elevated_two_triangle_runner();
  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner->msh_g);
  mesh.met.setSpace(MetSpace::Exp);
  mesh.setBasis(FEBasis::Lagrange);

  int interior_edge = -1;
  for(int local_edge = 0; local_edge < 3; local_edge++){
    if(mesh.fac2fac(0,local_edge) >= 0){
      interior_edge = local_edge;
      break;
    }
  }
  BOOST_REQUIRE_GE(interior_edge,0);
  const int control_point = mesh.fac2poi(0,3 + interior_edge);
  BOOST_REQUIRE_GE(control_point,0);

  mesh.coord(control_point,0) += 0.12;
  mesh.coord(control_point,1) += 0.12;
  const double perturbed_coordinates[2]
      = {mesh.coord(control_point,0),mesh.coord(control_point,1)};

  const auto sizeshape
      = get_quafun<MetricFieldAnalytical,2,2>(QuaFun::SizeShape);
  double objective_before = 0.0;
  for(int iface = 0; iface < mesh.nface; iface++){
    objective_before += sizeshape(
        mesh,AsDeg::Pk,AsDeg::Pk,iface,1.0);
  }

  const int point_tag = mesh.tag[1] + 1;
  for(int ipoin = 0; ipoin < mesh.npoin; ipoin++){
    mesh.poi2tag(1,ipoin) = point_tag;
  }
  mesh.poi2tag(1,control_point) = mesh.tag[1];
  runner->param->opt_smoo_niter = 1;

  BOOST_REQUIRE_NO_THROW(
      smoothInterior_Ball(
          mesh,QuaFun::SizeShape,1,2));

  double objective_after = 0.0;
  for(int iface = 0; iface < mesh.nface; iface++){
    objective_after += sizeshape(
        mesh,AsDeg::Pk,AsDeg::Pk,iface,1.0);
  }
  BOOST_CHECK_LE(objective_after,objective_before);
  const double displacement
      = std::hypot(
            mesh.coord(control_point,0) - perturbed_coordinates[0],
            mesh.coord(control_point,1) - perturbed_coordinates[1]);
  BOOST_CHECK_GT(displacement,1.e-10);

  intAr1 region(2);
  buildEdgeControlPointSmoothingRegion(
      mesh,2,0,interior_edge,region);
  check_p2_region_is_certified<2>(mesh,region);
}

BOOST_AUTO_TEST_CASE(sizeshape_p2_rejects_noncertified_final_candidate)
{
  auto runner = make_elevated_two_triangle_runner();
  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner->msh_g);
  mesh.met.setSpace(MetSpace::Exp);
  mesh.setBasis(FEBasis::Lagrange);

  int interior_edge = -1;
  for(int local_edge = 0; local_edge < 3; local_edge++){
    if(mesh.fac2fac(0,local_edge) >= 0){
      interior_edge = local_edge;
      break;
    }
  }
  BOOST_REQUIRE_GE(interior_edge,0);
  const int control_point = mesh.fac2poi(0,3 + interior_edge);

  intAr1 region(2);
  buildEdgeControlPointSmoothingRegion(
      mesh,2,0,interior_edge,region);
  BOOST_REQUIRE_EQUAL(region.get_n(),2);

  const double initial_coordinates[2]
      = {mesh.coord(control_point,0),mesh.coord(control_point,1)};
  bool found_noncertified_position = false;
  for(const double displacement : {0.5,1.0,2.0,4.0}){
    mesh.coord(control_point,0) = initial_coordinates[0] + displacement;
    mesh.coord(control_point,1) = initial_coordinates[1] + displacement;
    for(const int element : region){
      if(!classify_element_validity<2,2>(mesh,element)
              .accepted_conservatively()){
        found_noncertified_position = true;
        break;
      }
    }
    if(found_noncertified_position) break;
  }
  BOOST_REQUIRE(found_noncertified_position);

  const double rejected_coordinates[2]
      = {mesh.coord(control_point,0),mesh.coord(control_point,1)};
  double qavg0 = 0.0;
  double qmax0 = 0.0;
  double qavg1 = 0.0;
  double qmax1 = 0.0;
  const int status = smooballdiff<MetricFieldAnalytical,2,2>(
      mesh,control_point,region,&qavg0,&qmax0,&qavg1,&qmax1,
      QuaFun::SizeShape);

  BOOST_CHECK_NE(status,0);
  BOOST_CHECK_EQUAL(mesh.coord(control_point,0),rejected_coordinates[0]);
  BOOST_CHECK_EQUAL(mesh.coord(control_point,1),rejected_coordinates[1]);
}

BOOST_AUTO_TEST_CASE(sizeshape_p2_boundary_control_point_uses_cad_edge)
{
  auto runner = make_elevated_square_cad_runner();
  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner->msh_g);
  mesh.met.setSpace(MetSpace::Exp);
  mesh.setBasis(FEBasis::Lagrange);

  int seed_face = -1;
  int boundary_edge = -1;
  int control_point = -1;
  for(int iface = 0; iface < mesh.nface && control_point < 0; iface++){
    if(isdeadent(iface,mesh.fac2poi)) continue;
    for(int local_edge = 0; local_edge < 3; local_edge++){
      if(mesh.fac2fac(iface,local_edge) >= 0) continue;
      const int candidate = mesh.fac2poi(iface,3 + local_edge);
      const int boundary_record = mesh.poi2bpo[candidate];
      if(boundary_record < 0) continue;
      if(mesh.bpo2ibi(boundary_record,1) != 1) continue;
      seed_face = iface;
      boundary_edge = local_edge;
      control_point = candidate;
      break;
    }
  }
  BOOST_REQUIRE_GE(seed_face,0);
  BOOST_REQUIRE_GE(boundary_edge,0);
  BOOST_REQUIRE_GE(control_point,0);

  intAr1 region(2);
  buildEdgeControlPointSmoothingRegion(
      mesh,2,seed_face,boundary_edge,region);
  BOOST_REQUIRE_EQUAL(region.get_n(),1);
  BOOST_CHECK_EQUAL(region[0],seed_face);

  const int point_tag = mesh.tag[1] + 1;
  for(int ipoin = 0; ipoin < mesh.npoin; ipoin++){
    mesh.poi2tag(1,ipoin) = point_tag;
  }
  mesh.poi2tag(1,control_point) = mesh.tag[1];
  runner->param->opt_smoo_niter = 1;

  BOOST_REQUIRE_NO_THROW(
      smoothInterior_Ball(
          mesh,QuaFun::SizeShape,1,2));
  BOOST_CHECK_EQUAL(mesh.poi2tag(1,control_point),mesh.tag[1]);

  const int boundary_record = mesh.poi2bpo[control_point];
  BOOST_REQUIRE_GE(boundary_record,0);
  const int global_edge = mesh.bpo2ibi(boundary_record,2);
  BOOST_REQUIRE_GE(global_edge,0);
  const int reference = mesh.edg2ref[global_edge];
  const ego cad_edge = mesh.CAD.cad2edg[reference];
  const double parameter[2]
      = {mesh.bpo2rbi(boundary_record,0),0.0};
  double evaluation[18];
  BOOST_REQUIRE_EQUAL(
      EG_evaluate(cad_edge,parameter,evaluation),EGADS_SUCCESS);
  BOOST_CHECK_SMALL(mesh.coord(control_point,0) - evaluation[0],1.e-12);
  BOOST_CHECK_SMALL(mesh.coord(control_point,1) - evaluation[1],1.e-12);
  check_p2_region_is_certified<2>(mesh,region);
}

BOOST_AUTO_TEST_CASE(sizeshape_p2_complete_small_2d_mesh_smoothing)
{
  auto runner = make_elevated_square_cad_runner();
  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner->msh_g);
  mesh.met.setSpace(MetSpace::Exp);
  mesh.setBasis(FEBasis::Lagrange);
  runner->param->opt_smoo_niter = 1;

  BOOST_REQUIRE_NO_THROW(
      smoothInterior_Ball(
          mesh,QuaFun::SizeShape,1,2));
  BOOST_CHECK_EQUAL(mesh.curdeg,2);

  intAr1 all_faces(mesh.nface);
  for(int iface = 0; iface < mesh.nface; iface++){
    if(!isdeadent(iface,mesh.fac2poi)) all_faces.stack(iface);
  }
  check_p2_region_is_certified<2>(mesh,all_faces);
}

BOOST_AUTO_TEST_CASE(production_planar_p2_optimizer_runs_stepdistance_smoothing)
{
  auto runner = make_elevated_square_cad_runner();
  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner->msh_g);
  mesh.met.setSpace(MetSpace::Exp);
  mesh.setBasis(FEBasis::Lagrange);

  int control_point = -1;
  for(int iface = 0; iface < mesh.nface && control_point < 0; iface++){
    if(isdeadent(iface,mesh.fac2poi)) continue;
    for(int local_edge = 0; local_edge < 3; local_edge++){
      if(mesh.fac2fac(iface,local_edge) < 0) continue;
      const int candidate = mesh.fac2poi(iface,3 + local_edge);
      if(mesh.poi2bpo[candidate] >= 0) continue;
      control_point = candidate;
      break;
    }
  }
  BOOST_REQUIRE_GE(control_point,0);
  mesh.coord(control_point,0) += 0.01;
  mesh.coord(control_point,1) += 0.01;

  runner->param->opt_niter = 1;
  runner->param->opt_smoo_niter = 1;
  BOOST_REQUIRE_NO_THROW(runner->optimMesh());
  BOOST_CHECK_EQUAL(mesh.curdeg,2);

  intAr1 all_faces(mesh.nface);
  for(int iface = 0; iface < mesh.nface; iface++){
    if(!isdeadent(iface,mesh.fac2poi)) all_faces.stack(iface);
  }
  check_p2_region_is_certified<2>(mesh,all_faces);
}

BOOST_AUTO_TEST_CASE(sizeshape_p2_interior_3d_edge_control_point_and_shell)
{
  auto runner = make_elevated_cube_runner();
  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner->msh_g);
  mesh.met.setSpace(MetSpace::Exp);
  mesh.setBasis(FEBasis::Lagrange);

  int seed_tetrahedron = -1;
  int local_edge = -1;
  int control_point = -1;
  intAr1 region(10);
  for(int itetra = 0;
      itetra < mesh.nelem && control_point < 0;
      itetra++){
    if(isdeadent(itetra,mesh.tet2poi)) continue;
    for(int candidate_edge = 0; candidate_edge < 6; candidate_edge++){
      const int candidate = mesh.tet2poi(itetra,4 + candidate_edge);
      if(mesh.poi2bpo[candidate] >= 0) continue;
      buildEdgeControlPointSmoothingRegion(
          mesh,3,itetra,candidate_edge,region);
      if(region.get_n() < 2) continue;
      seed_tetrahedron = itetra;
      local_edge = candidate_edge;
      control_point = candidate;
      break;
    }
  }
  BOOST_REQUIRE_GE(seed_tetrahedron,0);
  BOOST_REQUIRE_GE(local_edge,0);
  BOOST_REQUIRE_GE(control_point,0);
  BOOST_REQUIRE_GE(region.get_n(),2);
  BOOST_CHECK(region_contains(region,seed_tetrahedron));

  mesh.coord(control_point,0) += 0.005;
  mesh.coord(control_point,1) += 0.004;
  mesh.coord(control_point,2) -= 0.003;
  const double perturbed_coordinates[3]
      = {mesh.coord(control_point,0),mesh.coord(control_point,1),
         mesh.coord(control_point,2)};
  check_p2_region_is_certified<3>(mesh,region);
  const auto sizeshape
      = get_quafun<MetricFieldAnalytical,3,3>(QuaFun::SizeShape);
  double objective_before = 0.0;
  for(const int itetra : region){
    objective_before += sizeshape(
        mesh,AsDeg::Pk,AsDeg::Pk,itetra,1.0);
  }

  const int point_tag = mesh.tag[1] + 1;
  for(int ipoin = 0; ipoin < mesh.npoin; ipoin++){
    mesh.poi2tag(1,ipoin) = point_tag;
  }
  mesh.poi2tag(1,control_point) = mesh.tag[1];
  runner->param->opt_smoo_niter = 1;

  BOOST_REQUIRE_NO_THROW(
      smoothInterior_Ball(
          mesh,QuaFun::SizeShape,1,2));

  double objective_after = 0.0;
  for(const int itetra : region){
    objective_after += sizeshape(
        mesh,AsDeg::Pk,AsDeg::Pk,itetra,1.0);
  }
  BOOST_CHECK_LE(objective_after,objective_before);
  const double displacement = std::sqrt(
      std::pow(mesh.coord(control_point,0) - perturbed_coordinates[0],2)
    + std::pow(mesh.coord(control_point,1) - perturbed_coordinates[1],2)
    + std::pow(mesh.coord(control_point,2) - perturbed_coordinates[2],2));
  BOOST_CHECK_GT(displacement,1.e-10);
  check_p2_region_is_certified<3>(mesh,region);
}
