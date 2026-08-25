// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_high_order_stepdistance_smoothing

#include "common_setup.hxx"

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "MetrisRunner/MetrisRunner.hxx"
#include "low_geo/validity.hxx"
#include "quality/msh_metqua.hxx"
#include "smoothing/low_smooballdiff.hxx"
#include "smoothing/msh_smooball.hxx"

#include <array>
#include <cmath>
#include <filesystem>
#include <memory>

using namespace Metris;

namespace
{

enum class StepDistanceVariant {
  Basic,
  CollapseBarrier,
  ShapeVolume,
  CavityTargetAverage
};

void configure_stepdistance_variant(MetrisParameters &parameters,
                                    StepDistanceVariant variant)
{
  parameters.iverb = 0;
  parameters.objective_quadrature_order = 5;
  parameters.objective_p = 1.25;
  parameters.step_distance_regularization = 1.e-6;
  parameters.step_distance_shape_volume
      = variant == StepDistanceVariant::ShapeVolume;
  parameters.step_distance_cavity_target_average
      = variant == StepDistanceVariant::CavityTargetAverage;
  parameters.step_distance_barrier_rho0 = 10.0;
  parameters.step_distance_barrier_beta
      = variant == StepDistanceVariant::CollapseBarrier ? 0.1
      : variant == StepDistanceVariant::ShapeVolume
        || variant == StepDistanceVariant::CavityTargetAverage ? 1.e6
                                                               : 0.0;
  parameters.opt_smoo_niter = 1;
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
         / "metris_high_order_stepdistance_smoothing_2d").string() + "/";

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
         / "metris_high_order_stepdistance_smoothing_3d").string() + "/";

  auto runner = std::make_unique<MetrisRunner>(
      nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner->degElevate(),1);
  return runner;
}

template<int gdim>
void require_certified_region(const MeshBase &mesh,
                              const intAr1 &region)
{
  for(const int element : region){
    const ElementValidityResult validity
        = classify_element_validity<gdim,2>(mesh,element);
    BOOST_REQUIRE_MESSAGE(
        validity.is_certified(),
        "element " << element << " is not Bernstein-certified");
  }
}

template<int gdim>
double stepdistance_mesh_objective(
    Mesh<MetricFieldAnalytical> &mesh)
{
  bool invalid_mesh = false;
  double minimum;
  double maximum;
  double average;
  const double objective
      = getmetquamesh<MetricFieldAnalytical,QuaFun::StepDistance>(
            mesh,gdim,AsDeg::Pk,AsDeg::Pk,
            &invalid_mesh,&minimum,&maximum,&average,NULL);
  BOOST_REQUIRE(!invalid_mesh);
  return objective;
}

template<int gdim>
double coordinate_displacement(
    const Mesh<MetricFieldAnalytical> &mesh,
    int point,
    const std::array<double,gdim> &coordinates)
{
  double squared_displacement = 0.0;
  for(int component = 0; component < gdim; component++){
    const double difference
        = mesh.coord(point,component) - coordinates[component];
    squared_displacement += difference*difference;
  }
  return std::sqrt(squared_displacement);
}

template<int gdim>
void require_noncertified_candidate_is_rejected(
    Mesh<MetricFieldAnalytical> &mesh,
    int control_point,
    const intAr1 &region)
{
  std::array<double,gdim> certified_coordinates{};
  for(int component = 0; component < gdim; component++){
    certified_coordinates[component] = mesh.coord(control_point,component);
  }

  bool found_noncertified_position = false;
  for(const double displacement : {0.5,1.0,2.0,4.0,8.0}){
    for(int component = 0; component < gdim; component++){
      const double direction = component == gdim - 1 ? -1.0 : 1.0;
      mesh.coord(control_point,component)
          = certified_coordinates[component] + direction*displacement;
    }
    for(const int element : region){
      if(!classify_element_validity<gdim,2>(mesh,element)
              .accepted_conservatively()){
        found_noncertified_position = true;
        break;
      }
    }
    if(found_noncertified_position) break;
  }
  BOOST_REQUIRE(found_noncertified_position);

  std::array<double,gdim> rejected_coordinates{};
  for(int component = 0; component < gdim; component++){
    rejected_coordinates[component] = mesh.coord(control_point,component);
  }
  double average_before = 0.0;
  double maximum_before = 0.0;
  double average_after = 0.0;
  double maximum_after = 0.0;
  const int status
      = smooballdiff<MetricFieldAnalytical,gdim,2>(
            mesh,control_point,region,
            &average_before,&maximum_before,
            &average_after,&maximum_after,
            QuaFun::StepDistance);
  BOOST_CHECK_NE(status,0);
  for(int component = 0; component < gdim; component++){
    BOOST_CHECK_EQUAL(
        mesh.coord(control_point,component),
        rejected_coordinates[component]);
  }
}

void check_planar_variant_smoothing(StepDistanceVariant variant)
{
  auto runner = make_elevated_two_triangle_runner();
  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner->msh_g);
  configure_stepdistance_variant(*runner->param,variant);
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

  intAr1 region(2);
  buildEdgeControlPointSmoothingRegion(
      mesh,2,0,interior_edge,region);
  BOOST_REQUIRE_EQUAL(region.get_n(),2);

  mesh.coord(control_point,0) += 0.12;
  mesh.coord(control_point,1) += 0.12;
  require_certified_region<2>(mesh,region);
  const std::array<double,2> perturbed_coordinates
      = {mesh.coord(control_point,0),mesh.coord(control_point,1)};
  const double objective_before = stepdistance_mesh_objective<2>(mesh);

  const int point_tag = mesh.tag[1] + 1;
  for(int point = 0; point < mesh.npoin; point++){
    mesh.poi2tag(1,point) = point_tag;
  }
  mesh.poi2tag(1,control_point) = mesh.tag[1];
  BOOST_REQUIRE_NO_THROW(
      smoothInterior_Ball(mesh,QuaFun::StepDistance,1,2));

  const double objective_after = stepdistance_mesh_objective<2>(mesh);
  BOOST_CHECK_LE(
      objective_after,
      objective_before + 1.e-12*(1.0 + std::abs(objective_before)));
  BOOST_CHECK_GT(
      coordinate_displacement<2>(
          mesh,control_point,perturbed_coordinates),
      1.e-10);
  require_certified_region<2>(mesh,region);
  require_noncertified_candidate_is_rejected<2>(
      mesh,control_point,region);
}

void check_tetrahedral_variant_smoothing(StepDistanceVariant variant)
{
  auto runner = make_elevated_cube_runner();
  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner->msh_g);
  configure_stepdistance_variant(*runner->param,variant);
  mesh.met.setSpace(MetSpace::Exp);
  mesh.setBasis(FEBasis::Lagrange);

  int seed_tetrahedron = -1;
  int local_edge = -1;
  int control_point = -1;
  intAr1 region(10);
  for(int tetrahedron = 0;
      tetrahedron < mesh.nelem && control_point < 0;
      tetrahedron++){
    if(isdeadent(tetrahedron,mesh.tet2poi)) continue;
    for(int candidate_edge = 0; candidate_edge < 6; candidate_edge++){
      const int candidate
          = mesh.tet2poi(tetrahedron,4 + candidate_edge);
      if(mesh.poi2bpo[candidate] >= 0) continue;
      buildEdgeControlPointSmoothingRegion(
          mesh,3,tetrahedron,candidate_edge,region);
      if(region.get_n() < 2) continue;
      seed_tetrahedron = tetrahedron;
      local_edge = candidate_edge;
      control_point = candidate;
      break;
    }
  }
  BOOST_REQUIRE_GE(seed_tetrahedron,0);
  BOOST_REQUIRE_GE(local_edge,0);
  BOOST_REQUIRE_GE(control_point,0);
  BOOST_REQUIRE_GE(region.get_n(),2);

  mesh.coord(control_point,0) += 0.02;
  mesh.coord(control_point,1) += 0.015;
  mesh.coord(control_point,2) -= 0.01;
  require_certified_region<3>(mesh,region);
  const std::array<double,3> perturbed_coordinates
      = {mesh.coord(control_point,0),mesh.coord(control_point,1),
         mesh.coord(control_point,2)};
  const double objective_before = stepdistance_mesh_objective<3>(mesh);

  const int point_tag = mesh.tag[1] + 1;
  for(int point = 0; point < mesh.npoin; point++){
    mesh.poi2tag(1,point) = point_tag;
  }
  mesh.poi2tag(1,control_point) = mesh.tag[1];

  BOOST_REQUIRE_NO_THROW(
      smoothInterior_Ball(mesh,QuaFun::StepDistance,1,2));

  const double objective_after = stepdistance_mesh_objective<3>(mesh);
  BOOST_CHECK_LE(
      objective_after,
      objective_before + 1.e-12*(1.0 + std::abs(objective_before)));
  BOOST_CHECK_GT(
      coordinate_displacement<3>(
          mesh,control_point,perturbed_coordinates),
      1.e-10);
  require_certified_region<3>(mesh,region);
  require_noncertified_candidate_is_rejected<3>(
      mesh,control_point,region);
}

} // namespace

BOOST_AUTO_TEST_SUITE(p2_stepdistance_2d)

BOOST_AUTO_TEST_CASE(basic_smoothing_and_validity)
{
  check_planar_variant_smoothing(StepDistanceVariant::Basic);
}

BOOST_AUTO_TEST_CASE(collapse_barrier_smoothing_and_validity)
{
  check_planar_variant_smoothing(StepDistanceVariant::CollapseBarrier);
}

BOOST_AUTO_TEST_CASE(shape_volume_smoothing_and_validity)
{
  check_planar_variant_smoothing(StepDistanceVariant::ShapeVolume);
}

BOOST_AUTO_TEST_CASE(cavity_target_average_smoothing_and_validity)
{
  check_planar_variant_smoothing(StepDistanceVariant::CavityTargetAverage);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(p2_stepdistance_3d)

BOOST_AUTO_TEST_CASE(basic_smoothing_and_validity)
{
  check_tetrahedral_variant_smoothing(StepDistanceVariant::Basic);
}

BOOST_AUTO_TEST_CASE(collapse_barrier_smoothing_and_validity)
{
  check_tetrahedral_variant_smoothing(StepDistanceVariant::CollapseBarrier);
}

BOOST_AUTO_TEST_CASE(shape_volume_smoothing_and_validity)
{
  check_tetrahedral_variant_smoothing(StepDistanceVariant::ShapeVolume);
}

BOOST_AUTO_TEST_CASE(cavity_target_average_smoothing_and_validity)
{
  check_tetrahedral_variant_smoothing(StepDistanceVariant::CavityTargetAverage);
}

BOOST_AUTO_TEST_SUITE_END()
