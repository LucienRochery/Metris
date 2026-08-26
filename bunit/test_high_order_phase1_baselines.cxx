// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_high_order_phase1_baselines

#include "common_setup.hxx"

#include "API/MetrisAPI.hxx"
#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "MetrisRunner/MetrisRunner.hxx"
#include "Adaptation/Insertion/EdgeSeed.hxx"
#include "Adaptation/Insertion/aux_insert.hxx"
#include "Adaptation/low_collapse.hxx"
#include "Adaptation/low_increasecav.hxx"
#include "aux_badEntHandler.hxx"
#include "cavity/msh_cavity.hxx"
#include "cavity/reconnect_geometry.hxx"
#include "aux_topo.hxx"
#include "ho_constants.hxx"
#include "low_eval.hxx"
#include "low_geo/validity.hxx"
#include "msh_checktopo.hxx"
#include "msh_metricCost.hxx"
#include "quality/low_metqua.hxx"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <numeric>
#include <string>
#include <vector>

using namespace Metris;

namespace
{

struct HighOrderBaselineCase
{
  const char *name;
  const char *mesh_path;
  int dimension;
  int expected_points;
  int expected_edges;
  int expected_faces;
  int expected_tetrahedra;
};

void metric_cost_quadrature_metric_2d(
    const AnaMetCtx *,
    const double *coordinate,
    double scale,
    int derivative_order,
    double *metric,
    double *metric_derivative)
{
  metric[0] = scale*(1.0 + coordinate[0]*coordinate[0]);
  metric[1] = 0.0;
  metric[2] = scale*(1.0 + coordinate[1]*coordinate[1]);
  if(derivative_order == 0) return;
  for(int ientry = 0; ientry < 6; ientry++){
    metric_derivative[ientry] = 0.0;
  }
}

std::string baseline_output_directory(const char *name, const char *stage)
{
  const std::filesystem::path directory
      = std::filesystem::temp_directory_path()
      / "metris_high_order_phase1"
      / name
      / stage;
  std::filesystem::create_directories(directory);
  return directory.string() + "/";
}

MetrisParameters baseline_parameters(const HighOrderBaselineCase &baseline_case,
                                     const char *stage)
{
  MetrisParameters parameters;
  parameters.setAnalyticalMetric(1);
  parameters.usrTarDeg = 2;
  parameters.adp_niter = 0;
  parameters.opt_niter = 0;
  parameters.iverb = 0;
  parameters.outmPrefix = baseline_output_directory(baseline_case.name,stage);
  return parameters;
}

void check_p2_mesh(const Mesh<MetricFieldAnalytical> &mesh,
                   const HighOrderBaselineCase &baseline_case)
{
  BOOST_REQUIRE_EQUAL(mesh.curdeg,2);
  BOOST_CHECK_EQUAL(mesh.idim,baseline_case.dimension);
  BOOST_CHECK_EQUAL(mesh.get_tdim(),baseline_case.dimension);
  BOOST_CHECK_EQUAL(mesh.npoin,baseline_case.expected_points);
  BOOST_CHECK_EQUAL(mesh.nedge,baseline_case.expected_edges);
  BOOST_CHECK_EQUAL(mesh.nface,baseline_case.expected_faces);
  BOOST_CHECK_EQUAL(mesh.nelem,baseline_case.expected_tetrahedra);
}

void run_baseline_case(const HighOrderBaselineCase &baseline_case)
{
  MetrisParameters elevation_parameters
      = baseline_parameters(baseline_case,"elevated");
  elevation_parameters.setMeshIn(baseline_case.mesh_path);

  MetrisRunner elevation_runner(nullptr,nullptr,elevation_parameters);
  auto &p1_mesh = static_cast<Mesh<MetricFieldAnalytical>&>(
      *elevation_runner.msh_g);
  const double p1_metric_cost
      = baseline_case.dimension == 2
      ? getMetricCost<MetricFieldAnalytical,2,2>(p1_mesh)
      : getMetricCost<MetricFieldAnalytical,3,3>(p1_mesh);

  BOOST_REQUIRE_NO_THROW(elevation_runner.runMetris());
  auto &elevated_mesh = static_cast<Mesh<MetricFieldAnalytical>&>(
      *elevation_runner.msh_g);
  check_p2_mesh(elevated_mesh,baseline_case);

  const double p2_metric_cost
      = baseline_case.dimension == 2
      ? getMetricCost<MetricFieldAnalytical,2,2>(elevated_mesh)
      : getMetricCost<MetricFieldAnalytical,3,3>(elevated_mesh);
  BOOST_CHECK(std::isfinite(p2_metric_cost));
  BOOST_CHECK_CLOSE(p2_metric_cost,p1_metric_cost,1.0e-10);

  MetrisAPI native_p2_data(elevation_runner);
  MetrisParameters native_parameters
      = baseline_parameters(baseline_case,"native");
  MetrisRunner native_runner(&native_p2_data,nullptr,native_parameters);
  BOOST_REQUIRE_NO_THROW(native_runner.runMetris());
  const auto &native_mesh = static_cast<const Mesh<MetricFieldAnalytical>&>(
      *native_runner.msh_g);
  check_p2_mesh(native_mesh,baseline_case);
}

template<int gdim>
void run_p2_cavity_insertion_case(
    const HighOrderBaselineCase &baseline_case)
{
  MetrisParameters parameters
      = baseline_parameters(baseline_case,"cavity");
  parameters.setMeshIn(baseline_case.mesh_path);
  MetrisRunner runner(nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner.degElevate(),1);

  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);
  intAr2 &element_to_point = mesh.ent2poi(gdim);
  const int seed_element = 0;
  BOOST_REQUIRE(!isdeadent(seed_element,element_to_point));

  const int inserted_point = mesh.npoin;
  mesh.set_npoin(inserted_point + 1);
  for(int coordinate = 0; coordinate < gdim; coordinate++){
    mesh.coord(inserted_point,coordinate) = 0.0;
    for(int vertex = 0; vertex < gdim + 1; vertex++){
      mesh.coord(inserted_point,coordinate)
          += mesh.coord(element_to_point(seed_element,vertex),coordinate)
           /(gdim + 1.0);
    }
  }
  constexpr int metric_entries = gdim*(gdim + 1)/2;
  for(int entry = 0; entry < metric_entries; entry++){
    mesh.met(inserted_point,entry) = 0.0;
    for(int vertex = 0; vertex < gdim + 1; vertex++){
      mesh.met(inserted_point,entry)
          += mesh.met(element_to_point(seed_element,vertex),entry)
           /(gdim + 1.0);
    }
  }
  mesh.poi2bpo[inserted_point] = -1;
  mesh.poi2bak[inserted_point] = -1;
  mesh.set_poi2ent(Vertex{inserted_point},0,-1);

  const int element_count_before = mesh.nentt(gdim);
  MshCavity cavity(
      gdim == 3 ? 10 : 0,
      gdim == 2 ? 10 : 0,
      0);
  cavity.ipins = inserted_point;
  cavity.inewp = 1;
  cavity.lcent(gdim).stack(seed_element);

  CavOprOpt options;
  options.allow_topological_correction = false;
  options.skip_topo_checks = false;
  options.qmax_nec = -1;
  options.qmax_suf = -1;
  options.qmax_iff = -1;
  CavWrkArrs work;
  CavOprInfo info;
  int acceptance_calls = 0;
  options.accept_completed_elements =
      [&](int tdim, int first_new, int end_new){
        BOOST_CHECK_EQUAL(tdim,gdim);
        acceptance_calls++;
        int candidate_count = 0;
        for(int element = first_new; element < end_new; element++){
          if(isdeadent(element,element_to_point)) continue;
          candidate_count++;
          BOOST_CHECK((
              classify_element_validity<gdim,2>(mesh,element).is_certified()));
        }
        BOOST_CHECK_EQUAL(candidate_count,gdim + 1);
        return true;
      };

  const int error = cavity_operator<MetricFieldAnalytical,2>(
      mesh,cavity,options,work,info,0);
  BOOST_REQUIRE_EQUAL(error,CAV_NOERR);
  BOOST_REQUIRE(info.done);
  BOOST_CHECK_EQUAL(acceptance_calls,1);
  BOOST_CHECK_EQUAL(work.lbad.get_n(),0);

  int completed_elements = 0;
  for(int element = element_count_before;
      element < mesh.nentt(gdim); element++){
    if(isdeadent(element,element_to_point)) continue;
    completed_elements++;
    const ElementValidityResult validity
        = classify_element_validity<gdim,2>(mesh,element);
    BOOST_CHECK(validity.is_certified());
  }
  BOOST_CHECK_EQUAL(completed_elements,gdim + 1);
}

void initialize_inserted_point_2d(Mesh<MetricFieldAnalytical>& mesh,
                                  int inserted_point,
                                  const double* coordinate,
                                  int seed_face){
  mesh.set_npoin(inserted_point + 1);
  for(int component = 0; component < 2; component++){
    mesh.coord(inserted_point,component) = coordinate[component];
  }
  for(int entry = 0; entry < 3; entry++){
    mesh.met(inserted_point,entry)
        = 0.5*(mesh.met(mesh.fac2poi(seed_face,0),entry)
             + mesh.met(mesh.fac2poi(seed_face,1),entry));
  }
  mesh.poi2bpo[inserted_point] = -1;
  mesh.poi2bak[inserted_point] = -1;
  mesh.set_poi2ent(Vertex{inserted_point},0,-1);
}

int find_live_face_with_edge(const Mesh<MetricFieldAnalytical>& mesh,
                             int first_face,
                             int endpoint0,
                             int endpoint1){
  for(int face = first_face; face < mesh.nface; face++){
    if(isdeadent(face,mesh.fac2poi)) continue;
    bool has_endpoint0 = false;
    bool has_endpoint1 = false;
    for(int vertex = 0; vertex < 3; vertex++){
      has_endpoint0 = has_endpoint0
                   || mesh.fac2poi(face,vertex) == endpoint0;
      has_endpoint1 = has_endpoint1
                   || mesh.fac2poi(face,vertex) == endpoint1;
    }
    if(has_endpoint0 && has_endpoint1) return face;
  }
  return -1;
}

void check_close_coordinate(const Mesh<MetricFieldAnalytical>& mesh,
                            int point,
                            const double* expected,
                            double tolerance = 2.e-13){
  for(int component = 0; component < 2; component++){
    BOOST_CHECK_SMALL(mesh.coord(point,component) - expected[component],
                      tolerance);
  }
}

double p1_triangle_objective(Mesh<MetricFieldAnalytical>& mesh,
                             const int* vertices){
  const int face_count = mesh.nface;
  mesh.set_nface(face_count + 1);
  for(int vertex = 0; vertex < 3; vertex++){
    mesh.fac2poi(face_count,vertex) = vertices[vertex];
  }
  const double objective = metqua<MetricFieldAnalytical,2,2,
      QuaFun::SizeShape>(mesh,AsDeg::P1,AsDeg::P1,face_count,1.0);
  mesh.set_nface(face_count);
  return objective;
}

template<QuaFun objective>
double live_planar_objective_sum(Mesh<MetricFieldAnalytical>& mesh,
                                 AsDeg geometry_degree){
  double sum = 0.;
  for(int face = 0; face < mesh.nface; face++){
    if(isdeadent(face,mesh.fac2poi)) continue;
    sum += metqua<MetricFieldAnalytical,2,2,objective>(
        mesh,geometry_degree,AsDeg::P1,face,1.0);
  }
  return sum;
}

} // namespace

BOOST_AUTO_TEST_CASE(high_order_phase1_four_baselines)
{
  const HighOrderBaselineCase baseline_cases[] = {
      {"triangle_p1_to_p2_and_native_p2",
       METRIS_ROOT_DIR "/examples/2D/misc/2tri2D.mesh",
       2,9,4,2,0},
      {"tetrahedron_p1_to_p2_and_native_p2",
       METRIS_ROOT_DIR
       "/bunit/meshes/high_order_phase1/one_tetrahedron_p1.mesh",
       3,10,6,4,1}};

  for(const HighOrderBaselineCase &baseline_case : baseline_cases){
    BOOST_TEST_CONTEXT(baseline_case.name){
      run_baseline_case(baseline_case);
    }
  }
}

BOOST_AUTO_TEST_CASE(completed_p2_cavity_insertions_are_certified)
{
  const HighOrderBaselineCase triangle_case = {
      "triangle_p2_cavity",
      METRIS_ROOT_DIR "/examples/2D/misc/2tri2D.mesh",
      2,9,4,2,0};
  const HighOrderBaselineCase tetrahedron_case = {
      "tetrahedron_p2_cavity",
      METRIS_ROOT_DIR
      "/bunit/meshes/high_order_phase1/one_tetrahedron_p1.mesh",
      3,10,6,4,1};

  run_p2_cavity_insertion_case<2>(triangle_case);
  run_p2_cavity_insertion_case<3>(tetrahedron_case);
}

BOOST_AUTO_TEST_CASE(completed_p2_cavity_rejection_rolls_back)
{
  const HighOrderBaselineCase baseline_case = {
      "triangle_p2_cavity_rejection",
      METRIS_ROOT_DIR "/examples/2D/misc/2tri2D.mesh",
      2,9,4,2,0};
  MetrisParameters parameters = baseline_parameters(baseline_case,"cavity");
  parameters.setMeshIn(baseline_case.mesh_path);
  MetrisRunner runner(nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner.degElevate(),1);
  auto& mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);

  const int seed_face = 0;
  double centroid[2] = {0.,0.};
  for(int vertex = 0; vertex < 3; vertex++){
    for(int component = 0; component < 2; component++){
      centroid[component] += mesh.coord(mesh.fac2poi(seed_face,vertex),component)/3.;
    }
  }
  const int inserted_point = mesh.npoin;
  initialize_inserted_point_2d(mesh,inserted_point,centroid,seed_face);
  const int points_before = mesh.npoin;
  const int faces_before = mesh.nface;

  MshCavity cavity(0,10,0);
  cavity.ipins = inserted_point;
  cavity.inewp = 1;
  cavity.lcfac.stack(seed_face);
  CavOprOpt options;
  options.allow_topological_correction = false;
  int acceptance_calls = 0;
  options.accept_completed_elements =
      [&](int, int first_new, int end_new){
        acceptance_calls++;
        BOOST_CHECK_EQUAL(end_new - first_new,3);
        return false;
      };
  CavWrkArrs work;
  CavOprInfo info;

  const int error = cavity_operator<MetricFieldAnalytical,2>(
      mesh,cavity,options,work,info,0);
  BOOST_CHECK_EQUAL(error,CAV_ERR_OBJECTIVE);
  BOOST_CHECK_EQUAL(acceptance_calls,1);
  BOOST_CHECK(!info.done);
  BOOST_CHECK_EQUAL(mesh.npoin,points_before);
  BOOST_CHECK_EQUAL(mesh.nface,faces_before);
  BOOST_CHECK(!isdeadent(seed_face,mesh.fac2poi));
}

BOOST_AUTO_TEST_CASE(curved_p2_edge_split_preserves_exact_child_geometry)
{
  const HighOrderBaselineCase baseline_case = {
      "curved_triangle_p2_edge_split",
      METRIS_ROOT_DIR "/examples/2D/misc/2tri2D.mesh",
      2,9,4,2,0};
  MetrisParameters parameters = baseline_parameters(baseline_case,"cavity");
  parameters.setMeshIn(baseline_case.mesh_path);
  MetrisRunner runner(nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner.degElevate(),1);
  auto& mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);
  BOOST_REQUIRE(mesh.getBasis() == FEBasis::Lagrange);

  int parent_face = -1;
  int parent_local_edge = -1;
  int neighbor_face = -1;
  for(int face = 0; face < mesh.nface && parent_face < 0; face++){
    for(int local_edge = 0; local_edge < 3; local_edge++){
      if(mesh.fac2fac(face,local_edge) >= 0){
        parent_face = face;
        parent_local_edge = local_edge;
        neighbor_face = mesh.fac2fac(face,local_edge);
        break;
      }
    }
  }
  BOOST_REQUIRE(parent_face >= 0);
  BOOST_REQUIRE(neighbor_face >= 0);

  const int endpoint0
      = mesh.fac2poi(parent_face,lnoed2[parent_local_edge][0]);
  const int endpoint1
      = mesh.fac2poi(parent_face,lnoed2[parent_local_edge][1]);
  const int parent_control
      = mesh.fac2poi(parent_face,3 + parent_local_edge);
  int parent_points[3] = {endpoint0,endpoint1,parent_control};

  const double dx = mesh.coord(endpoint1,0) - mesh.coord(endpoint0,0);
  const double dy = mesh.coord(endpoint1,1) - mesh.coord(endpoint0,1);
  const double length = std::sqrt(dx*dx + dy*dy);
  mesh.coord(parent_control,0)
      = 0.5*(mesh.coord(endpoint0,0) + mesh.coord(endpoint1,0))
      - 0.06*dy/length;
  mesh.coord(parent_control,1)
      = 0.5*(mesh.coord(endpoint0,1) + mesh.coord(endpoint1,1))
      + 0.06*dx/length;
  BOOST_REQUIRE((classify_element_validity<2,2>(
      mesh,parent_face).is_certified()));
  BOOST_REQUIRE((classify_element_validity<2,2>(
      mesh,neighbor_face).is_certified()));

  const double split_bary[2] = {0.5,0.5};
  const double left_bary[2] = {0.75,0.25};
  const double right_bary[2] = {0.25,0.75};
  double split_coordinate[2], left_coordinate[2], right_coordinate[2];
  eval1<2,2>(mesh.coord,parent_points,FEBasis::Lagrange,
             DifVar::None,DifVar::None,split_bary,
             split_coordinate,nullptr,nullptr);
  eval1<2,2>(mesh.coord,parent_points,FEBasis::Lagrange,
             DifVar::None,DifVar::None,left_bary,
             left_coordinate,nullptr,nullptr);
  eval1<2,2>(mesh.coord,parent_points,FEBasis::Lagrange,
             DifVar::None,DifVar::None,right_bary,
             right_coordinate,nullptr,nullptr);

  int inherited_controls[4];
  int inherited_count = 0;
  const int old_faces[2] = {parent_face,neighbor_face};
  for(int old_face : old_faces){
    for(int local_edge = 0; local_edge < 3; local_edge++){
      const int control = mesh.fac2poi(old_face,3 + local_edge);
      if(control != parent_control){
        inherited_controls[inherited_count++] = control;
      }
    }
  }
  BOOST_REQUIRE_EQUAL(inherited_count,4);

  int opposite_vertices[2] = {-1,-1};
  for(int old_rank = 0; old_rank < 2; old_rank++){
    for(int vertex = 0; vertex < 3; vertex++){
      const int point = mesh.fac2poi(old_faces[old_rank],vertex);
      if(point != endpoint0 && point != endpoint1){
        opposite_vertices[old_rank] = point;
        break;
      }
    }
    BOOST_REQUIRE(opposite_vertices[old_rank] >= 0);
  }

  const int inserted_point = mesh.npoin;
  initialize_inserted_point_2d(
      mesh,inserted_point,split_coordinate,parent_face);
  const int first_new_face = mesh.nface;
  MshCavity cavity(0,10,0);
  cavity.ipins = inserted_point;
  cavity.inewp = 1;
  cavity.lcfac.stack(parent_face);
  cavity.lcfac.stack(neighbor_face);
  cavity.split_edge_points.set_n(3);
  for(int inode = 0; inode < 3; inode++){
    cavity.split_edge_points[inode] = parent_points[inode];
  }
  cavity.split_edge_barycentric[0] = split_bary[0];
  cavity.split_edge_barycentric[1] = split_bary[1];
  cavity.preserve_split_edge_geometry = true;

  CavOprOpt options;
  options.allow_topological_correction = false;
  CavWrkArrs work;
  CavOprInfo info;
  const int error = cavity_operator<MetricFieldAnalytical,2>(
      mesh,cavity,options,work,info,0);
  BOOST_REQUIRE_EQUAL(error,CAV_NOERR);
  BOOST_REQUIRE(info.done);

  const int left_face = find_live_face_with_edge(
      mesh,first_new_face,endpoint0,inserted_point);
  const int right_face = find_live_face_with_edge(
      mesh,first_new_face,inserted_point,endpoint1);
  BOOST_REQUIRE(left_face >= 0);
  BOOST_REQUIRE(right_face >= 0);
  const int left_edge = getedgfac(mesh,left_face,endpoint0,inserted_point);
  const int right_edge = getedgfac(mesh,right_face,inserted_point,endpoint1);
  BOOST_TEST_CONTEXT("left restricted child"){
    check_close_coordinate(mesh,mesh.fac2poi(left_face,3 + left_edge),
                           left_coordinate);
  }
  BOOST_TEST_CONTEXT("right restricted child"){
    check_close_coordinate(mesh,mesh.fac2poi(right_face,3 + right_edge),
                           right_coordinate);
  }

  for(int opposite_vertex : opposite_vertices){
    const int spoke_face = find_live_face_with_edge(
        mesh,first_new_face,inserted_point,opposite_vertex);
    BOOST_REQUIRE(spoke_face >= 0);
    const int spoke_edge = getedgfac(
        mesh,spoke_face,inserted_point,opposite_vertex);
    const int spoke_control = mesh.fac2poi(spoke_face,3 + spoke_edge);
    const double affine_midpoint[2] = {
        0.5*(mesh.coord(inserted_point,0)
           + mesh.coord(opposite_vertex,0)),
        0.5*(mesh.coord(inserted_point,1)
           + mesh.coord(opposite_vertex,1))};
    BOOST_TEST_CONTEXT("affine new spoke to " << opposite_vertex){
      check_close_coordinate(mesh,spoke_control,affine_midpoint);
    }
  }

  for(int inherited_control : inherited_controls){
    bool reused = false;
    for(int face = first_new_face; face < mesh.nface && !reused; face++){
      if(isdeadent(face,mesh.fac2poi)) continue;
      for(int inode = 3; inode < getnnode(2,2); inode++){
        reused = reused || mesh.fac2poi(face,inode) == inherited_control;
      }
    }
    BOOST_CHECK(reused);
  }
  for(int face = first_new_face; face < mesh.nface; face++){
    if(isdeadent(face,mesh.fac2poi)) continue;
    BOOST_CHECK((classify_element_validity<2,2>(mesh,face).is_certified()));
  }
}

BOOST_AUTO_TEST_CASE(bezier_p2_split_coefficients_use_de_casteljau)
{
  const HighOrderBaselineCase baseline_case = {
      "bezier_triangle_p2_edge_split",
      METRIS_ROOT_DIR "/examples/2D/misc/2tri2D.mesh",
      2,9,4,2,0};
  MetrisParameters parameters = baseline_parameters(baseline_case,"cavity");
  parameters.setMeshIn(baseline_case.mesh_path);
  MetrisRunner runner(nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner.degElevate(),1);
  auto& mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);

  int parent_face = -1;
  int parent_local_edge = -1;
  for(int face = 0; face < mesh.nface && parent_face < 0; face++){
    for(int local_edge = 0; local_edge < 3; local_edge++){
      if(mesh.fac2fac(face,local_edge) >= 0){
        parent_face = face;
        parent_local_edge = local_edge;
        break;
      }
    }
  }
  BOOST_REQUIRE(parent_face >= 0);
  const int parent0
      = mesh.fac2poi(parent_face,lnoed2[parent_local_edge][0]);
  const int parent2
      = mesh.fac2poi(parent_face,lnoed2[parent_local_edge][1]);
  const int parent1 = mesh.fac2poi(parent_face,3 + parent_local_edge);

  mesh.coord(parent1,0)
      = 0.5*(mesh.coord(parent0,0) + mesh.coord(parent2,0)) + 0.04;
  mesh.coord(parent1,1)
      = 0.5*(mesh.coord(parent0,1) + mesh.coord(parent2,1)) - 0.03;
  mesh.setBasis(FEBasis::Bezier);
  BOOST_REQUIRE(mesh.getBasis() == FEBasis::Bezier);

  const int inserted_point = mesh.npoin;
  mesh.set_npoin(inserted_point + 1);
  MshCavity cavity(0,1,0);
  cavity.ipins = inserted_point;
  cavity.split_edge_points.set_n(3);
  cavity.split_edge_points[0] = parent0;
  cavity.split_edge_points[1] = parent2;
  cavity.split_edge_points[2] = parent1;
  cavity.split_edge_barycentric[0] = 0.63;
  cavity.split_edge_barycentric[1] = 0.37;
  cavity.preserve_split_edge_geometry = true;

  double left_coefficient[2];
  double right_coefficient[2];
  BOOST_REQUIRE((initialize_quadratic_split_child_coefficient<
      MetricFieldAnalytical,2>(
          mesh,cavity,inserted_point,parent0,left_coefficient)));
  BOOST_REQUIRE((initialize_quadratic_split_child_coefficient<
      MetricFieldAnalytical,2>(
          mesh,cavity,inserted_point,parent2,right_coefficient)));
  for(int component = 0; component < 2; component++){
    const double expected_left
        = 0.63*mesh.coord(parent0,component)
        + 0.37*mesh.coord(parent1,component);
    const double expected_right
        = 0.63*mesh.coord(parent1,component)
        + 0.37*mesh.coord(parent2,component);
    BOOST_CHECK_SMALL(left_coefficient[component] - expected_left,2.e-14);
    BOOST_CHECK_SMALL(right_coefficient[component] - expected_right,2.e-14);
  }
}

BOOST_AUTO_TEST_CASE(p2_growth_probe_is_local_and_uses_completed_geometry)
{
  const HighOrderBaselineCase baseline_case = {
      "p2_local_growth_probe",
      METRIS_ROOT_DIR
      "/bunit/meshes/high_order_phase1/four_triangle_strip_p1.mesh",
      2,15,9,4,0};
  MetrisParameters parameters = baseline_parameters(baseline_case,"growth");
  parameters.setMeshIn(baseline_case.mesh_path);
  MetrisRunner runner(nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner.degElevate(),1);
  auto& mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);
  mesh.met.setSpace(MetSpace::Exp);

  const int seed_face = 1;
  int seed_edge = -1;
  for(int edge = 0; edge < 3; edge++){
    const int point0 = mesh.fac2poi(seed_face,lnoed2[edge][0]);
    const int point1 = mesh.fac2poi(seed_face,lnoed2[edge][1]);
    if((point0 == 1 && point1 == 4) || (point0 == 4 && point1 == 1)){
      seed_edge = edge;
      break;
    }
  }
  BOOST_REQUIRE(seed_edge >= 0);
  const int parent_control = mesh.fac2poi(seed_face,3 + seed_edge);
  mesh.coord(parent_control,0) += 0.08;
  BOOST_REQUIRE((classify_element_validity<2,2>(mesh,1).is_certified()));
  BOOST_REQUIRE((classify_element_validity<2,2>(mesh,2).is_certified()));

  MshCavity cavity;
  EdgeSeed insertion_seed(mesh,cavity,2,2,seed_face,seed_edge);
  BOOST_REQUIRE_EQUAL(cavity.lcfac.get_n(),2);
  cavity.ipins = mesh.newpoint(
      PointType::Vertex,insertion_seed.tdimp,insertion_seed.iseed);
  cavity.inewp = 1;
  BOOST_REQUIRE_EQUAL(aux_bisecPointLen(
      mesh,insertion_seed,mesh.poi2bpo[cavity.ipins],false,cavity),0);
  BOOST_REQUIRE_EQUAL(increase_cavity(mesh,cavity,false,0,1),0);
  BOOST_REQUIRE_EQUAL(cavity.lcfac.get_n(),2);

  std::vector<bool> in_cavity(mesh.nface,false);
  for(int face : cavity.lcfac) in_cavity[face] = true;
  int outside = -1;
  for(int face : cavity.lcfac){
    for(int edge = 0; edge < 3; edge++){
      const int neighbor = mesh.fac2fac(face,edge);
      if(neighbor >= 0 && !in_cavity[neighbor]){
        outside = neighbor;
        break;
      }
    }
    if(outside >= 0) break;
  }
  BOOST_REQUIRE(outside >= 0);

  double current_p2 = metqua<MetricFieldAnalytical,2,2,
      QuaFun::SizeShape>(mesh,AsDeg::Pk,AsDeg::P1,outside,1.0);
  double enlarged_p2 = 0.;
  double current_p1 = metqua<MetricFieldAnalytical,2,2,
      QuaFun::SizeShape>(mesh,AsDeg::P1,AsDeg::P1,outside,1.0);
  double enlarged_p1 = 0.;
  int incident_count = 0;
  int enlarged_count = 0;
  for(int edge = 0; edge < 3; edge++){
    const int neighbor = mesh.fac2fac(outside,edge);
    if(neighbor >= 0 && in_cavity[neighbor]){
      incident_count++;
      int inside_edge = -1;
      for(int candidate_edge = 0; candidate_edge < 3; candidate_edge++){
        if(mesh.fac2fac(neighbor,candidate_edge) == outside){
          inside_edge = candidate_edge;
          break;
        }
      }
      BOOST_REQUIRE(inside_edge >= 0);
      const int vertices[3] = {
          cavity.ipins,
          mesh.fac2poi(neighbor,lnoed2[inside_edge][0]),
          mesh.fac2poi(neighbor,lnoed2[inside_edge][1])};
      bool valid = false;
      current_p2 += evaluate_completed_p2_cavity_cone<
          MetricFieldAnalytical,2,QuaFun::SizeShape>(
              mesh,cavity,neighbor,vertices,&valid);
      BOOST_REQUIRE(valid);
      current_p1 += p1_triangle_objective(mesh,vertices);
    }else{
      const int vertices[3] = {
          cavity.ipins,
          mesh.fac2poi(outside,lnoed2[edge][0]),
          mesh.fac2poi(outside,lnoed2[edge][1])};
      bool valid = false;
      enlarged_p2 += evaluate_completed_p2_cavity_cone<
          MetricFieldAnalytical,2,QuaFun::SizeShape>(
              mesh,cavity,outside,vertices,&valid);
      BOOST_REQUIRE(valid);
      enlarged_p1 += p1_triangle_objective(mesh,vertices);
      enlarged_count++;
    }
  }
  BOOST_REQUIRE_EQUAL(incident_count,1);
  BOOST_REQUIRE_EQUAL(enlarged_count,2);
  BOOST_CHECK(std::abs(current_p2 - current_p1) > 1.e-8
           || std::abs(enlarged_p2 - enlarged_p1) > 1.e-8);

  std::vector<double> qualities(mesh.nface);
  std::vector<int> sorted_faces(mesh.nface);
  std::iota(sorted_faces.begin(),sorted_faces.end(),0);
  for(int face = 0; face < mesh.nface; face++){
    qualities[face] = metqua<MetricFieldAnalytical,2,2,
        QuaFun::SizeShape>(mesh,AsDeg::Pk,AsDeg::P1,face,1.0);
  }
  std::sort(sorted_faces.begin(),sorted_faces.end(),
      [&](int left, int right){ return qualities[left] > qualities[right]; });
  BadEntHandler handler(2,100.,0.00001);
  handler.setCallbacks(
      [&](int face){ return qualities[face]; },
      [&](int face){ return isdeadent(face,mesh.fac2poi); });
  handler.seedFromSortedIDs(sorted_faces);

  int matching_probes = 0;
  cavity.inspect_growth_probe = [&](const CavityGrowthProbeInfo& probe){
    if(probe.outside_element != outside) return;
    matching_probes++;
    BOOST_CHECK_EQUAL(probe.topological_dimension,2);
    BOOST_CHECK_EQUAL(probe.geometry_degree,2);
    BOOST_CHECK_EQUAL(probe.current_cavity_element_count,2);
    BOOST_CHECK_EQUAL(probe.incident_cavity_element_count,incident_count);
    BOOST_CHECK_EQUAL(probe.current_local_element_count,incident_count + 1);
    BOOST_CHECK_EQUAL(probe.enlarged_local_element_count,enlarged_count);
    BOOST_CHECK(probe.current_configuration_valid);
    BOOST_CHECK(probe.enlarged_configuration_valid);
    BOOST_CHECK_CLOSE(probe.current_objective,current_p2,2.e-10);
    BOOST_CHECK_CLOSE(probe.enlarged_objective,enlarged_p2,2.e-10);
  };
  const int growth_result = increase_cavity_quality<
      MetricFieldAnalytical,QuaFun::SizeShape>(
          mesh,cavity,2,1,handler,0);
  BOOST_CHECK(growth_result == 0 || growth_result == -1);
  BOOST_CHECK_EQUAL(matching_probes,1);
}

BOOST_AUTO_TEST_CASE(initial_validity_growth_uses_completed_p2_geometry)
{
  const HighOrderBaselineCase baseline_case = {
      "p2_initial_validity_growth",
      METRIS_ROOT_DIR
      "/bunit/meshes/high_order_phase1/four_triangle_strip_p1.mesh",
      2,15,9,4,0};
  MetrisParameters parameters = baseline_parameters(baseline_case,"validity");
  parameters.setMeshIn(baseline_case.mesh_path);
  MetrisRunner runner(nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner.degElevate(),1);
  auto& mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);
  mesh.met.setSpace(MetSpace::Exp);

  const int seed_face = 1;
  int seed_edge = -1;
  for(int edge = 0; edge < 3; edge++){
    const int point0 = mesh.fac2poi(seed_face,lnoed2[edge][0]);
    const int point1 = mesh.fac2poi(seed_face,lnoed2[edge][1]);
    if((point0 == 1 && point1 == 4) || (point0 == 4 && point1 == 1)){
      seed_edge = edge;
      break;
    }
  }
  BOOST_REQUIRE(seed_edge >= 0);

  MshCavity cavity;
  EdgeSeed insertion_seed(mesh,cavity,2,2,seed_face,seed_edge);
  cavity.ipins = mesh.newpoint(
      PointType::Vertex,insertion_seed.tdimp,insertion_seed.iseed);
  cavity.inewp = 1;
  BOOST_REQUIRE_EQUAL(aux_bisecPointLen(
      mesh,insertion_seed,mesh.poi2bpo[cavity.ipins],false,cavity),0);
  BOOST_REQUIRE_EQUAL(cavity.lcfac.get_n(),2);

  std::vector<bool> initially_in_cavity(mesh.nface,false);
  for(int face : cavity.lcfac) initially_in_cavity[face] = true;
  int outside = -1;
  int inside = -1;
  int inside_edge = -1;
  for(int face : cavity.lcfac){
    for(int edge = 0; edge < 3; edge++){
      const int neighbor = mesh.fac2fac(face,edge);
      if(neighbor >= 0 && !initially_in_cavity[neighbor]){
        outside = neighbor;
        inside = face;
        inside_edge = edge;
        break;
      }
    }
    if(outside >= 0) break;
  }
  BOOST_REQUIRE(outside >= 0);
  BOOST_REQUIRE(inside >= 0);
  BOOST_REQUIRE(inside_edge >= 0);

  const int vertices[3] = {
      cavity.ipins,
      mesh.fac2poi(inside,lnoed2[inside_edge][0]),
      mesh.fac2poi(inside,lnoed2[inside_edge][1])};
  double p1_measure = 0.;
  BOOST_REQUIRE((isvalideltP1<2,2>(
      mesh,vertices,NULL,NULL,&p1_measure,-1.)));

  const int inherited_control = mesh.fac2poi(inside,3 + inside_edge);
  const double base[2] = {
      mesh.coord(inherited_control,0),
      mesh.coord(inherited_control,1)};
  const double edge_dx = mesh.coord(vertices[2],0) - mesh.coord(vertices[1],0);
  const double edge_dy = mesh.coord(vertices[2],1) - mesh.coord(vertices[1],1);
  const double edge_length = std::sqrt(edge_dx*edge_dx + edge_dy*edge_dy);
  bool found_p2_only_invalidity = false;
  for(int step = 1; step <= 80 && !found_p2_only_invalidity; step++){
    for(int sign : {-1,1}){
      const double displacement = sign*0.01*step;
      mesh.coord(inherited_control,0) = base[0] - displacement*edge_dy/edge_length;
      mesh.coord(inherited_control,1) = base[1] + displacement*edge_dx/edge_length;
      if(!classify_element_validity<2,2>(mesh,inside).is_certified()
      || !classify_element_validity<2,2>(mesh,outside).is_certified()){
        continue;
      }
      bool completed_valid = true;
      (void)evaluate_completed_p2_cavity_cone<
          MetricFieldAnalytical,2,QuaFun::SizeShape>(
              mesh,cavity,inside,vertices,&completed_valid);
      if(!completed_valid){
        found_p2_only_invalidity = true;
        break;
      }
    }
  }
  BOOST_REQUIRE_MESSAGE(found_p2_only_invalidity,
      "Could not construct a P1-valid but P2-invalid seed reconnection");
  BOOST_REQUIRE((isvalideltP1<2,2>(
      mesh,vertices,NULL,NULL,&p1_measure,-1.)));

  BOOST_REQUIRE_EQUAL(increase_cavity(mesh,cavity,false,0,1),0);
  bool outside_was_absorbed = false;
  for(int face : cavity.lcfac){
    outside_was_absorbed = outside_was_absorbed || face == outside;
  }
  BOOST_CHECK(outside_was_absorbed);
  BOOST_CHECK(cavity.lcfac.get_n() > 2);

  std::vector<bool> finally_in_cavity(mesh.nface,false);
  for(int face : cavity.lcfac) finally_in_cavity[face] = true;
  for(int face : cavity.lcfac){
    for(int edge = 0; edge < 3; edge++){
      const int neighbor = mesh.fac2fac(face,edge);
      if(neighbor >= 0 && finally_in_cavity[neighbor]) continue;
      const int cone_vertices[3] = {
          cavity.ipins,
          mesh.fac2poi(face,lnoed2[edge][0]),
          mesh.fac2poi(face,lnoed2[edge][1])};
      bool completed_valid = false;
      (void)evaluate_completed_p2_cavity_cone<
          MetricFieldAnalytical,2,QuaFun::SizeShape>(
              mesh,cavity,face,cone_vertices,&completed_valid);
      BOOST_CHECK(completed_valid);
    }
  }
}

BOOST_AUTO_TEST_CASE(p2_collapse_uses_completed_geometry_without_split_children)
{
  const HighOrderBaselineCase baseline_case = {
      "p2_completed_collapse",
      METRIS_ROOT_DIR
      "/bunit/meshes/high_order_phase1/interior_short_edge_p1.mesh",
      2,35,23,12,0};
  MetrisParameters parameters = baseline_parameters(baseline_case,"collapse");
  parameters.setMeshIn(baseline_case.mesh_path);
  MetrisRunner runner(nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner.degElevate(),1);
  auto& mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);
  mesh.met.setSpace(MetSpace::Exp);
  for(int point = 0; point < mesh.npoin; point++){
    for(int component = 0; component < 2; component++){
      mesh.coord(point,component) *= 0.05;
    }
  }
  for(int point = 0; point < mesh.bak->npoin; point++){
    for(int component = 0; component < 2; component++){
      mesh.bak->coord(point,component) *= 0.05;
    }
  }

  const int endpoint0 = 5;
  const int endpoint1 = 6;
  int seed_face = -1;
  int seed_edge = -1;
  for(int face = 0; face < mesh.nface && seed_face < 0; face++){
    for(int edge = 0; edge < 3; edge++){
      const int point0 = mesh.fac2poi(face,lnoed2[edge][0]);
      const int point1 = mesh.fac2poi(face,lnoed2[edge][1]);
      if((point0 == endpoint0 && point1 == endpoint1)
      || (point0 == endpoint1 && point1 == endpoint0)){
        seed_face = face;
        seed_edge = edge;
        break;
      }
    }
  }
  BOOST_REQUIRE(seed_face >= 0);
  BOOST_REQUIRE(seed_edge >= 0);
  BOOST_REQUIRE_EQUAL(mesh.getpoitdim(endpoint0),2);
  BOOST_REQUIRE_EQUAL(mesh.getpoitdim(endpoint1),2);

  // Make the disappearing edge truly curved. Its midpoint still determines
  // the replacement vertex, but no child of this trace may survive collapse.
  const int disappearing_control = mesh.fac2poi(seed_face,3 + seed_edge);
  mesh.coord(disappearing_control,1) += 0.00025;
  for(int face = 0; face < mesh.nface; face++){
    BOOST_REQUIRE((classify_element_validity<2,2>(mesh,face).is_certified()));
  }

  const double objective_before =
      live_planar_objective_sum<QuaFun::SizeShape>(mesh,AsDeg::Pk);
  std::vector<double> qualities(mesh.nface);
  std::vector<int> sorted_faces(mesh.nface);
  std::iota(sorted_faces.begin(),sorted_faces.end(),0);
  for(int face = 0; face < mesh.nface; face++){
    qualities[face] = metqua<MetricFieldAnalytical,2,2,
        QuaFun::SizeShape>(mesh,AsDeg::Pk,AsDeg::P1,face,1.0);
  }
  std::sort(sorted_faces.begin(),sorted_faces.end(),
      [&](int left, int right){ return qualities[left] > qualities[right]; });
  BadEntHandler handler(2,100.,0.00001);
  handler.setCallbacks(
      [&](int face){ return qualities[face]; },
      [&](int face){ return isdeadent(face,mesh.fac2poi); });
  handler.seedFromSortedIDs(sorted_faces);

  MshCavity cavity(0,100,20);
  CavWrkArrs work;
  intAr1 cavity_errors(CAV_ERR_NERROR);
  cavity_errors.set_n(CAV_ERR_NERROR);
  for(int error = 0; error < CAV_ERR_NERROR; error++){
    cavity_errors[error] = 0;
  }
  const int first_new_face = mesh.nface;
  const int result = collapseEdge<MetricFieldAnalytical>(
      mesh,2,seed_face,seed_edge,0.,cavity,work,cavity_errors,
      handler,0,1,2);
  BOOST_REQUIRE_EQUAL(result,0);
  BOOST_CHECK_EQUAL(cavity.split_edge_points.get_n(),0);
  BOOST_CHECK(!cavity.preserve_split_edge_geometry);
  BOOST_CHECK(mesh.isdeadpoint(endpoint0));
  BOOST_CHECK(mesh.isdeadpoint(endpoint1));
  BOOST_REQUIRE(cavity.ipins >= 0);
  BOOST_CHECK(!mesh.isdeadpoint(cavity.ipins));

  int new_live_faces = 0;
  int affine_spokes = 0;
  for(int face = first_new_face; face < mesh.nface; face++){
    if(isdeadent(face,mesh.fac2poi)) continue;
    new_live_faces++;
    BOOST_CHECK((classify_element_validity<2,2>(mesh,face).is_certified()));
    for(int edge = 0; edge < 3; edge++){
      const int point0 = mesh.fac2poi(face,lnoed2[edge][0]);
      const int point1 = mesh.fac2poi(face,lnoed2[edge][1]);
      if(point0 != cavity.ipins && point1 != cavity.ipins) continue;
      const int control = mesh.fac2poi(face,3 + edge);
      for(int component = 0; component < 2; component++){
        BOOST_CHECK_SMALL(
            mesh.coord(control,component)
            - 0.5*(mesh.coord(point0,component)
                 + mesh.coord(point1,component)),
            5.e-14);
      }
      affine_spokes++;
    }
  }
  BOOST_CHECK(new_live_faces > 0);
  BOOST_CHECK(affine_spokes > 0);

  const double objective_after =
      live_planar_objective_sum<QuaFun::SizeShape>(mesh,AsDeg::Pk);
  BOOST_CHECK(objective_after < objective_before);
}

BOOST_AUTO_TEST_CASE(p1_collapse_compatibility_is_retained)
{
  const HighOrderBaselineCase baseline_case = {
      "p1_collapse_compatibility",
      METRIS_ROOT_DIR
      "/bunit/meshes/high_order_phase1/interior_short_edge_p1.mesh",
      2,12,10,12,0};
  MetrisParameters parameters = baseline_parameters(baseline_case,"collapse");
  parameters.usrTarDeg = 1;
  parameters.setMeshIn(baseline_case.mesh_path);
  MetrisRunner runner(nullptr,nullptr,parameters);
  auto& mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);
  BOOST_REQUIRE_EQUAL(mesh.curdeg,1);
  mesh.met.setSpace(MetSpace::Exp);
  for(int point = 0; point < mesh.npoin; point++){
    for(int component = 0; component < 2; component++){
      mesh.coord(point,component) *= 0.05;
    }
  }
  for(int point = 0; point < mesh.bak->npoin; point++){
    for(int component = 0; component < 2; component++){
      mesh.bak->coord(point,component) *= 0.05;
    }
  }

  const int endpoint0 = 5;
  const int endpoint1 = 6;
  int seed_face = -1;
  int seed_edge = -1;
  for(int face = 0; face < mesh.nface && seed_face < 0; face++){
    for(int edge = 0; edge < 3; edge++){
      const int point0 = mesh.fac2poi(face,lnoed2[edge][0]);
      const int point1 = mesh.fac2poi(face,lnoed2[edge][1]);
      if((point0 == endpoint0 && point1 == endpoint1)
      || (point0 == endpoint1 && point1 == endpoint0)){
        seed_face = face;
        seed_edge = edge;
        break;
      }
    }
  }
  BOOST_REQUIRE(seed_face >= 0);
  BOOST_REQUIRE(seed_edge >= 0);

  const double objective_before =
      live_planar_objective_sum<QuaFun::SizeShape>(mesh,AsDeg::P1);
  std::vector<double> qualities(mesh.nface);
  std::vector<int> sorted_faces(mesh.nface);
  std::iota(sorted_faces.begin(),sorted_faces.end(),0);
  for(int face = 0; face < mesh.nface; face++){
    qualities[face] = metqua<MetricFieldAnalytical,2,2,
        QuaFun::SizeShape>(mesh,AsDeg::P1,AsDeg::P1,face,1.0);
  }
  std::sort(sorted_faces.begin(),sorted_faces.end(),
      [&](int left, int right){ return qualities[left] > qualities[right]; });
  BadEntHandler handler(2,100.,0.00001);
  handler.setCallbacks(
      [&](int face){ return qualities[face]; },
      [&](int face){ return isdeadent(face,mesh.fac2poi); });
  handler.seedFromSortedIDs(sorted_faces);

  MshCavity cavity(0,100,20);
  CavWrkArrs work;
  intAr1 cavity_errors(CAV_ERR_NERROR);
  cavity_errors.set_n(CAV_ERR_NERROR);
  for(int error = 0; error < CAV_ERR_NERROR; error++){
    cavity_errors[error] = 0;
  }
  const int first_new_face = mesh.nface;
  const int result = collapseEdge<MetricFieldAnalytical>(
      mesh,2,seed_face,seed_edge,0.,cavity,work,cavity_errors,
      handler,0,1,2);
  BOOST_REQUIRE_EQUAL(result,0);
  BOOST_CHECK_EQUAL(cavity.split_edge_points.get_n(),0);
  BOOST_CHECK(!cavity.preserve_split_edge_geometry);
  BOOST_CHECK(mesh.isdeadpoint(endpoint0));
  BOOST_CHECK(mesh.isdeadpoint(endpoint1));

  int new_live_faces = 0;
  for(int face = first_new_face; face < mesh.nface; face++){
    if(isdeadent(face,mesh.fac2poi)) continue;
    new_live_faces++;
    BOOST_CHECK((isvalideltP1<2,2>(mesh,face)));
  }
  BOOST_CHECK(new_live_faces > 0);
  const double objective_after =
      live_planar_objective_sum<QuaFun::SizeShape>(mesh,AsDeg::P1);
  BOOST_CHECK(objective_after < objective_before);
}

BOOST_AUTO_TEST_CASE(p1_insertion_compatibility_keeps_optional_final_hook_empty)
{
  const HighOrderBaselineCase baseline_case = {
      "triangle_p1_cavity_compatibility",
      METRIS_ROOT_DIR "/examples/2D/misc/2tri2D.mesh",
      2,4,4,2,0};
  MetrisParameters parameters = baseline_parameters(baseline_case,"cavity");
  parameters.usrTarDeg = 1;
  parameters.setMeshIn(baseline_case.mesh_path);
  MetrisRunner runner(nullptr,nullptr,parameters);
  auto& mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);
  BOOST_REQUIRE_EQUAL(mesh.curdeg,1);

  double centroid[2] = {0.,0.};
  for(int vertex = 0; vertex < 3; vertex++){
    for(int component = 0; component < 2; component++){
      centroid[component] += mesh.coord(mesh.fac2poi(0,vertex),component)/3.;
    }
  }
  const int inserted_point = mesh.npoin;
  initialize_inserted_point_2d(mesh,inserted_point,centroid,0);
  MshCavity cavity(0,10,0);
  cavity.ipins = inserted_point;
  cavity.inewp = 1;
  cavity.lcfac.stack(0);
  CavOprOpt options;
  BOOST_CHECK(!options.accept_completed_elements);
  CavWrkArrs work;
  CavOprInfo info;
  const int error = cavity_operator<MetricFieldAnalytical,1>(
      mesh,cavity,options,work,info,0);
  BOOST_CHECK_EQUAL(error,CAV_NOERR);
  BOOST_CHECK(info.done);
}

BOOST_AUTO_TEST_CASE(high_order_phase1_final_validity_is_unconditional)
{
  const HighOrderBaselineCase baseline_case = {
      "triangle_invalid_p2",
      METRIS_ROOT_DIR "/examples/2D/misc/2tri2D.mesh",
      2,9,4,2,0};

  MetrisParameters elevation_parameters
      = baseline_parameters(baseline_case,"elevated");
  elevation_parameters.setMeshIn(baseline_case.mesh_path);
  MetrisRunner elevation_runner(nullptr,nullptr,elevation_parameters);
  BOOST_REQUIRE_EQUAL(elevation_runner.degElevate(),1);

  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(
      *elevation_runner.msh_g);
  BOOST_REQUIRE_EQUAL(mesh.curdeg,2);
  const int control_point = mesh.fac2poi(0,3);
  BOOST_REQUIRE(control_point >= 0);
  mesh.coord(control_point,0) = 10.0;
  mesh.coord(control_point,1) = 10.0;

  MetrisAPI invalid_p2_data(elevation_runner);
  MetrisParameters native_parameters
      = baseline_parameters(baseline_case,"native");
  MetrisRunner native_runner(&invalid_p2_data,nullptr,native_parameters);
  BOOST_CHECK_THROW(native_runner.runMetris(),MetrisExcept);
}

BOOST_AUTO_TEST_CASE(high_order_final_topology_reports_uncertified)
{
  const HighOrderBaselineCase baseline_case = {
      "triangle_uncertified_p2",
      METRIS_ROOT_DIR "/examples/2D/misc/2tri2D.mesh",
      2,9,4,2,0};

  MetrisParameters parameters
      = baseline_parameters(baseline_case,"elevated");
  parameters.setMeshIn(baseline_case.mesh_path);
  MetrisRunner runner(nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner.degElevate(),1);

  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);
  constexpr auto ordering = ORDELT(2);
  const double coupling = std::sqrt(0.75);
  for(int local_node = 0; local_node < getnnod2(2); local_node++){
    const int point = mesh.fac2poi(0,local_node);
    const double u = ordering[2][local_node][1]/2.0;
    const double v = ordering[2][local_node][2]/2.0;
    mesh.coord(point,0) = u + coupling*v*v;
    mesh.coord(point,1) = v + coupling*u*u;
  }

  const ElementValidityResult validity
      = classify_element_validity<2,2>(mesh,0);
  BOOST_REQUIRE(validity.is_uncertified());

  bool threw_uncertified = false;
  try{
    check_topo(mesh,0);
  }catch(const MetrisExcept &exception){
    threw_uncertified
        = std::string(exception.what()).find("Uncertified")
          != std::string::npos;
  }
  BOOST_CHECK(threw_uncertified);
}

BOOST_AUTO_TEST_CASE(metric_cost_uses_objective_quadrature_order)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.setAnalyticalMetric(
      AnaMetFun(metric_cost_quadrature_metric_2d));

  Mesh<MetricFieldAnalytical> mesh;
  mesh.idim = 2;
  mesh.curdeg = 1;
  mesh.strdeg = 1;
  mesh.forceBasisFlag(FEBasis::Lagrange);
  mesh.param = &parameters;
  mesh.set_npoin(3);
  mesh.set_nface(1);
  mesh.met.forceBasisFlag(FEBasis::Lagrange);
  mesh.met.forceSpaceFlag(MetSpace::Exp);
  mesh.met.setAnalyticalMetric(parameters);

  mesh.coord(0,0) = 0.0; mesh.coord(0,1) = 0.0;
  mesh.coord(1,0) = 1.0; mesh.coord(1,1) = 0.0;
  mesh.coord(2,0) = 0.0; mesh.coord(2,1) = 1.0;
  for(int inode = 0; inode < 3; inode++){
    mesh.fac2poi(0,inode) = inode;
    mesh.met.getMetPhys(DifVar::None,MetSpace::Exp,
                        mesh.coord[inode],mesh.met[inode],NULL);
  }

  parameters.objective_quadrature_order = 0;
  const double historical_cost
      = getMetricCost<MetricFieldAnalytical,2,2>(mesh);
  parameters.objective_quadrature_order = 5;
  const double fifth_order_cost
      = getMetricCost<MetricFieldAnalytical,2,2>(mesh);

  BOOST_CHECK(std::isfinite(historical_cost));
  BOOST_CHECK(std::isfinite(fifth_order_cost));
  BOOST_CHECK_GT(std::abs(historical_cost - fifth_order_cost),1.e-6);
}
