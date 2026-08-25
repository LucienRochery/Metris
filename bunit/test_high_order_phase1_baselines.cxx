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
#include "cavity/msh_cavity.hxx"
#include "ho_constants.hxx"
#include "low_geo/validity.hxx"
#include "msh_checktopo.hxx"
#include "msh_metricCost.hxx"

#include <cmath>
#include <filesystem>
#include <string>

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

  const int error = cavity_operator<MetricFieldAnalytical,2>(
      mesh,cavity,options,work,info,0);
  BOOST_REQUIRE_EQUAL(error,CAV_NOERR);
  BOOST_REQUIRE(info.done);
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
