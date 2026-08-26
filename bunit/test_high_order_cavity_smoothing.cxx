// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_high_order_cavity_smoothing

#include "common_setup.hxx"

#include "Adaptation/Insertion/EdgeSeed.hxx"
#include "Adaptation/Insertion/aux_insert.hxx"
#include "Adaptation/low_increasecav.hxx"
#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "MetrisRunner/MetrisRunner.hxx"
#include "aux_badEntHandler.hxx"
#include "cavity/msh_cavity.hxx"
#include "cavity/reconnect_geometry.hxx"
#include "ho_constants.hxx"
#include "low_eval.hxx"
#include "low_geo/validity.hxx"
#include "quality/low_metqua.hxx"
#include "quality/objective_quadrature_sample.hxx"
#include "quality/objective_quadrature_value.hxx"
#include "quality/quafun_sizeshape.hxx"
#include "quality/simplex_quadrature.hxx"
#include "smoothing/low_smooballdiff.hxx"
#include "smoothing/msh_smooball.hxx"

#include <array>
#include <cmath>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

using namespace Metris;

namespace
{

std::string output_directory(const char* name)
{
  const std::filesystem::path directory
      = std::filesystem::temp_directory_path()
      / "metris_high_order_cavity_smoothing" / name;
  std::filesystem::create_directories(directory);
  return directory.string() + "/";
}

MetrisParameters parameters_for(const char* name, const char* mesh_path)
{
  MetrisParameters parameters;
  parameters.setAnalyticalMetric(1);
  parameters.setMeshIn(mesh_path);
  parameters.usrTarDeg = 2;
  parameters.adp_niter = 0;
  parameters.opt_niter = 0;
  parameters.iflag1 = 1;
  parameters.iverb = 0;
  parameters.outmPrefix = output_directory(name);
  return parameters;
}

template<int dimension>
class TemporaryQuadraticConeWorkspace
{
public:
  explicit TemporaryQuadraticConeWorkspace(
      Mesh<MetricFieldAnalytical>& mesh)
    : msh(mesh), element_count(mesh.nentt(dimension)),
      point_count(mesh.npoin), temporary_element(element_count)
  {
    msh.set_nentt(dimension,element_count + 1);
    for(int edge = 0;
        edge < quadratic_simplex_edge_count<dimension>(); edge++){
      const int point = msh.newpoitopo(
          PointType::CtrlPt,dimension,temporary_element);
      msh.ent2poi(dimension)(
          temporary_element,quadratic_simplex_edge_node<dimension>(edge))
          = point;
    }
  }

  ~TemporaryQuadraticConeWorkspace()
  {
    msh.set_nentt(dimension,element_count);
    msh.set_npoin(point_count);
  }

  int element() const noexcept { return temporary_element; }

private:
  Mesh<MetricFieldAnalytical>& msh;
  int element_count;
  int point_count;
  int temporary_element;
};

struct PlanarSplitCase
{
  int source_face = -1;
  int neighbor_face = -1;
  int local_edge = -1;
  int endpoint0 = -1;
  int endpoint1 = -1;
  int parent_control = -1;
};

PlanarSplitCase find_and_curve_interior_edge(
    Mesh<MetricFieldAnalytical>& mesh)
{
  PlanarSplitCase result;
  for(int face = 0; face < mesh.nface && result.source_face < 0; face++){
    for(int edge = 0; edge < 3; edge++){
      const int neighbor = mesh.fac2fac(face,edge);
      if(neighbor < 0) continue;
      result.source_face = face;
      result.neighbor_face = neighbor;
      result.local_edge = edge;
      break;
    }
  }
  BOOST_REQUIRE(result.source_face >= 0);
  result.endpoint0
      = mesh.fac2poi(result.source_face,lnoed2[result.local_edge][0]);
  result.endpoint1
      = mesh.fac2poi(result.source_face,lnoed2[result.local_edge][1]);
  result.parent_control
      = mesh.fac2poi(result.source_face,3 + result.local_edge);

  const double dx = mesh.coord(result.endpoint1,0)
                  - mesh.coord(result.endpoint0,0);
  const double dy = mesh.coord(result.endpoint1,1)
                  - mesh.coord(result.endpoint0,1);
  const double length = std::sqrt(dx*dx + dy*dy);
  mesh.coord(result.parent_control,0) -= 0.045*dy/length;
  mesh.coord(result.parent_control,1) += 0.045*dx/length;
  BOOST_REQUIRE((classify_element_validity<2,2>(
      mesh,result.source_face).is_certified()));
  BOOST_REQUIRE((classify_element_validity<2,2>(
      mesh,result.neighbor_face).is_certified()));
  return result;
}

void place_split_point(
    Mesh<MetricFieldAnalytical>& mesh,
    MshCavity& cavity,
    const PlanarSplitCase& split,
    double endpoint0_weight)
{
  const int parent_points[3] = {
      split.endpoint0,split.endpoint1,split.parent_control};
  const double barycentric[2] = {
      endpoint0_weight,1. - endpoint0_weight};
  double coordinate[2];
  eval1<2,2>(mesh.coord,parent_points,FEBasis::Lagrange,
             DifVar::None,DifVar::None,barycentric,
             coordinate,nullptr,nullptr);

  cavity.ipins = mesh.newpoint(PointType::Vertex,2,split.source_face);
  cavity.inewp = 1;
  for(int component = 0; component < 2; component++){
    mesh.coord(cavity.ipins,component) = coordinate[component];
  }
  mesh.poi2bpo[cavity.ipins] = -1;
  mesh.set_poi2ent(
      Vertex{cavity.ipins},2,split.source_face);
  BOOST_REQUIRE_EQUAL(mesh.interpMetBack(cavity.ipins),0);

  cavity.lcfac.stack(split.source_face);
  cavity.lcfac.stack(split.neighbor_face);
  cavity.split_edge_points.set_n(3);
  cavity.split_edge_points[0] = split.endpoint0;
  cavity.split_edge_points[1] = split.endpoint1;
  cavity.split_edge_points[2] = split.parent_control;
  cavity.split_edge_barycentric[0] = endpoint0_weight;
  cavity.split_edge_barycentric[1] = 1. - endpoint0_weight;
  cavity.preserve_split_edge_geometry = true;
}

template<int dimension>
void install_cone_vertices(
    Mesh<MetricFieldAnalytical>& mesh,
    const MshCavity& cavity,
    int source_element,
    int opposite,
    int candidate)
{
  mesh.ent2poi(dimension)(candidate,0) = cavity.ipins;
  if constexpr(dimension == 2){
    mesh.fac2poi(candidate,lnoed2[0][0])
        = mesh.fac2poi(source_element,lnoed2[opposite][0]);
    mesh.fac2poi(candidate,lnoed2[0][1])
        = mesh.fac2poi(source_element,lnoed2[opposite][1]);
  }else{
    mesh.tet2poi(candidate,lnofa3[0][0])
        = mesh.tet2poi(source_element,lnofa3[opposite][0]);
    mesh.tet2poi(candidate,lnofa3[0][1])
        = mesh.tet2poi(source_element,lnofa3[opposite][1]);
    mesh.tet2poi(candidate,lnofa3[0][2])
        = mesh.tet2poi(source_element,lnofa3[opposite][2]);
  }
}

template<int dimension>
void check_released_map_variation(
    Mesh<MetricFieldAnalytical>& mesh,
    MshCavity& cavity,
    int source_element,
    int opposite,
    int temporary_element)
{
  install_cone_vertices<dimension>(
      mesh,cavity,source_element,opposite,temporary_element);
  complete_quadratic_cone_element<MetricFieldAnalytical,
      dimension,dimension>(
          mesh,cavity,source_element,temporary_element,
          QuadraticConeSpokePolicy::ReleasedAffine);

  const int apex_local = mesh.getverent<2>(
      temporary_element,dimension,cavity.ipins);
  BOOST_REQUIRE_EQUAL(apex_local,0);
  for(int edge = 0;
      edge < quadratic_simplex_edge_count<dimension>(); edge++){
    const int local0 = quadratic_simplex_edge_vertex<dimension>(edge,0);
    const int local1 = quadratic_simplex_edge_vertex<dimension>(edge,1);
    if(local0 != apex_local && local1 != apex_local) continue;
    const int midpoint = mesh.ent2poi(dimension)(
        temporary_element,quadratic_simplex_edge_node<dimension>(edge));
    for(int component = 0; component < dimension; component++){
      BOOST_CHECK_SMALL(
          mesh.coord(midpoint,component)
          - 0.5*(mesh.coord(mesh.ent2poi(dimension)(temporary_element,local0),component)
                 + mesh.coord(mesh.ent2poi(dimension)(temporary_element,local1),component)),
          5.e-14);
    }
  }

  std::array<double,dimension + 1> barycentric{};
  barycentric[0] = 0.23;
  for(int vertex = 1; vertex < dimension + 1; vertex++){
    barycentric[vertex] = (1. - barycentric[0])/dimension;
  }
  double original_image[dimension];
  if constexpr(dimension == 2){
    eval2<dimension,2>(mesh.coord,mesh.ent2poi(dimension)[temporary_element],
        mesh.getBasis(),DifVar::None,DifVar::None,barycentric.data(),
        original_image,nullptr,nullptr);
  }else{
    eval3<dimension,2>(mesh.coord,mesh.ent2poi(dimension)[temporary_element],
        mesh.getBasis(),DifVar::None,DifVar::None,barycentric.data(),
        original_image,nullptr,nullptr);
  }

  const double displacement[3] = {1.7e-4,-2.3e-4,1.1e-4};
  for(int component = 0; component < dimension; component++){
    mesh.coord(cavity.ipins,component) += displacement[component];
  }
  complete_quadratic_cone_element<MetricFieldAnalytical,
      dimension,dimension>(
          mesh,cavity,source_element,temporary_element,
          QuadraticConeSpokePolicy::ReleasedAffine);
  double displaced_image[dimension];
  if constexpr(dimension == 2){
    eval2<dimension,2>(mesh.coord,mesh.ent2poi(dimension)[temporary_element],
        mesh.getBasis(),DifVar::None,DifVar::None,barycentric.data(),
        displaced_image,nullptr,nullptr);
  }else{
    eval3<dimension,2>(mesh.coord,mesh.ent2poi(dimension)[temporary_element],
        mesh.getBasis(),DifVar::None,DifVar::None,barycentric.data(),
        displaced_image,nullptr,nullptr);
  }
  for(int component = 0; component < dimension; component++){
    BOOST_CHECK_SMALL(
        displaced_image[component] - original_image[component]
        - barycentric[0]*displacement[component],
        2.e-13);
    mesh.coord(cavity.ipins,component) -= displacement[component];
  }
}

template<int dimension>
void check_objective_gradient_by_finite_difference(
    Mesh<MetricFieldAnalytical>& mesh,
    MshCavity& cavity)
{
  constexpr int metric_count = dimension*(dimension + 1)/2;
  struct FrozenSample
  {
    std::array<double,dimension + 1> barycentric{};
    std::array<double,metric_count> metric{};
    double objective_weight = 0.;
  };
  struct FrozenCone
  {
    int source = -1;
    int opposite = -1;
    std::vector<FrozenSample> samples;
  };

  std::vector<FrozenCone> frozen_cones;
  const int temporary_element = mesh.nentt(dimension) - 1;
  std::vector<bool> in_cavity(mesh.nentt(dimension),false);
  for(int element : cavity.lcent(dimension)) in_cavity[element] = true;
  const SimplexQuadratureView<dimension> quadrature
      = get_objective_quadrature<dimension>(
            mesh.param->objective_quadrature_order);
  const ObjectiveQuadratureTheta theta_mode
      = objective_quadrature_theta_mode<QuaFun::SizeShape>(mesh);
  for(int source : cavity.lcent(dimension)){
    for(int opposite = 0; opposite < dimension + 1; opposite++){
      const int neighbor = mesh.ent2ent(dimension)(source,opposite);
      if(neighbor >= 0 && in_cavity[neighbor]) continue;
      install_cone_vertices<dimension>(
          mesh,cavity,source,opposite,temporary_element);
      complete_quadratic_cone_element<MetricFieldAnalytical,
          dimension,dimension>(
              mesh,cavity,source,temporary_element,
              QuadraticConeSpokePolicy::ReleasedAffine);
      FrozenCone cone;
      cone.source = source;
      cone.opposite = opposite;
      for(int iquad = 0; iquad < quadrature.size(); iquad++){
        const ObjectiveQuadratureSample<dimension,dimension,2> sample
            = prepare_objective_quadrature_sample<
                MetricFieldAnalytical,dimension,dimension,2>(
                    mesh,AsDeg::P1,
                    mesh.ent2poi(dimension)[temporary_element],
                    quadrature[iquad],theta_mode);
        BOOST_REQUIRE(sample.theta_is_valid);
        FrozenSample frozen;
        frozen.barycentric = sample.barycentric_coordinates;
        frozen.metric = sample.metric;
        frozen.objective_weight
            = sample.quadrature_weight*sample.theta;
        cone.samples.push_back(frozen);
      }
      frozen_cones.push_back(std::move(cone));
    }
  }

  auto frozen_value = [&](){
    double sum = 0.;
    for(const FrozenCone& cone : frozen_cones){
      install_cone_vertices<dimension>(
          mesh,cavity,cone.source,cone.opposite,temporary_element);
      complete_quadratic_cone_element<MetricFieldAnalytical,
          dimension,dimension>(
              mesh,cavity,cone.source,temporary_element,
              QuadraticConeSpokePolicy::ReleasedAffine);
      for(const FrozenSample& sample : cone.samples){
        sum += sample.objective_weight
             * quafun_sizeshape<MetricFieldAnalytical,
                   dimension,dimension,double>(
                       mesh,AsDeg::Pk,AsDeg::P1,
                       mesh.ent2poi(dimension)[temporary_element],
                       sample.barycentric.data(),sample.metric.data());
      }
    }
    return sum;
  };

  double value;
  double gradient[dimension];
  double hessian[dimension*(dimension + 1)/2];
  BOOST_REQUIRE((evaluate_released_p2_cavity_objective<
      MetricFieldAnalytical,dimension>(
          mesh,cavity,QuaFun::SizeShape,0,
          value,gradient,hessian)));
  BOOST_CHECK(std::isfinite(value));
  BOOST_CHECK_CLOSE_FRACTION(value,frozen_value(),5.e-14);
  for(double entry : hessian) BOOST_CHECK(std::isfinite(entry));

  for(int variable = 0; variable < dimension; variable++){
    const double original = mesh.coord(cavity.ipins,variable);
    const double step = 2.e-6*MAX(1.0,std::abs(original));
    mesh.coord(cavity.ipins,variable) = original - step;
    const double minus_value = frozen_value();
    mesh.coord(cavity.ipins,variable) = original + step;
    const double plus_value = frozen_value();
    mesh.coord(cavity.ipins,variable) = original;
    const double finite_difference = (plus_value - minus_value)/(2.*step);
    BOOST_TEST_MESSAGE(
        "released P2 dim=" << dimension
        << " component=" << variable
        << " analytic=" << gradient[variable]
        << " finite_difference=" << finite_difference);
    BOOST_CHECK_SMALL(
        gradient[variable] - finite_difference,
        3.e-7*(1. + std::abs(finite_difference)));
  }
}

template<int dimension>
void check_stepdistance_released_evaluation_is_finite(
    Mesh<MetricFieldAnalytical>& mesh,
    MshCavity& cavity)
{
  double value;
  double gradient[dimension];
  double hessian[dimension*(dimension + 1)/2];
  BOOST_REQUIRE((evaluate_released_p2_cavity_objective<
      MetricFieldAnalytical,dimension>(
          mesh,cavity,QuaFun::StepDistance,0,
          value,gradient,hessian)));
  BOOST_CHECK(std::isfinite(value));
  for(double component : gradient) BOOST_CHECK(std::isfinite(component));
  for(double entry : hessian) BOOST_CHECK(std::isfinite(entry));
}

template<int dimension>
double preserved_cavity_objective(
    Mesh<MetricFieldAnalytical>& mesh,
    MshCavity& cavity,
    double& maximum,
    int& count)
{
  double sum = 0.;
  maximum = -1.0e30;
  count = 0;
  const intAr1& elements = cavity.lcent(dimension);
  std::vector<bool> in_cavity(mesh.nentt(dimension),false);
  for(int element : elements) in_cavity[element] = true;
  for(int source : elements){
    for(int opposite = 0; opposite < dimension + 1; opposite++){
      const int neighbor = mesh.ent2ent(dimension)(source,opposite);
      if(neighbor >= 0 && in_cavity[neighbor]) continue;
      int vertices[dimension + 1];
      vertices[0] = cavity.ipins;
      if constexpr(dimension == 2){
        vertices[1] = mesh.fac2poi(source,lnoed2[opposite][0]);
        vertices[2] = mesh.fac2poi(source,lnoed2[opposite][1]);
      }else{
        vertices[1] = mesh.tet2poi(source,lnofa3[opposite][0]);
        vertices[2] = mesh.tet2poi(source,lnofa3[opposite][1]);
        vertices[3] = mesh.tet2poi(source,lnofa3[opposite][2]);
      }
      bool valid = false;
      const double objective = evaluate_completed_p2_cavity_cone<
          MetricFieldAnalytical,dimension,QuaFun::SizeShape>(
              mesh,cavity,source,vertices,&valid);
      BOOST_REQUIRE(valid);
      sum += objective;
      maximum = MAX(maximum,objective);
      count++;
    }
  }
  return sum;
}

} // namespace

BOOST_AUTO_TEST_CASE(released_planar_p2_spokes_have_complete_chain_rule)
{
  MetrisParameters parameters = parameters_for(
      "released_planar_chain_rule",
      METRIS_ROOT_DIR "/examples/2D/misc/2tri2D.mesh");
  MetrisRunner runner(nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner.degElevate(),1);
  auto& mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);
  mesh.met.setSpace(MetSpace::Exp);
  const PlanarSplitCase split = find_and_curve_interior_edge(mesh);
  MshCavity cavity(0,10,0);
  place_split_point(mesh,cavity,split,0.64);

  TemporaryQuadraticConeWorkspace<2> workspace(mesh);
  int boundary_opposite = -1;
  for(int opposite = 0; opposite < 3; opposite++){
    if(mesh.fac2fac(split.source_face,opposite) != split.neighbor_face){
      boundary_opposite = opposite;
      break;
    }
  }
  BOOST_REQUIRE(boundary_opposite >= 0);
  check_released_map_variation<2>(
      mesh,cavity,split.source_face,boundary_opposite,workspace.element());
  check_objective_gradient_by_finite_difference<2>(mesh,cavity);
  check_stepdistance_released_evaluation_is_finite<2>(mesh,cavity);
}

BOOST_AUTO_TEST_CASE(released_tetrahedral_p2_spokes_have_complete_chain_rule)
{
  MetrisParameters parameters = parameters_for(
      "released_tetrahedral_chain_rule",
      METRIS_ROOT_DIR
      "/bunit/meshes/high_order_phase1/one_tetrahedron_p1.mesh");
  MetrisRunner runner(nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner.degElevate(),1);
  auto& mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);
  mesh.met.setSpace(MetSpace::Exp);

  MshCavity cavity(10,0,0);
  cavity.lctet.stack(0);
  cavity.ipins = mesh.newpoint(PointType::Vertex,3,0);
  cavity.inewp = 1;
  for(int component = 0; component < 3; component++){
    mesh.coord(cavity.ipins,component) = 0.;
    for(int vertex = 0; vertex < 4; vertex++){
      mesh.coord(cavity.ipins,component)
          += mesh.coord(mesh.tet2poi(0,vertex),component)
           * (vertex == 0 ? 0.34 : 0.22);
    }
  }
  mesh.poi2bpo[cavity.ipins] = -1;
  mesh.set_poi2ent(Vertex{cavity.ipins},3,0);
  BOOST_REQUIRE_EQUAL(mesh.interpMetBack(cavity.ipins),0);

  TemporaryQuadraticConeWorkspace<3> workspace(mesh);
  check_released_map_variation<3>(mesh,cavity,0,0,workspace.element());
  check_objective_gradient_by_finite_difference<3>(mesh,cavity);
  check_stepdistance_released_evaluation_is_finite<3>(mesh,cavity);
}

BOOST_AUTO_TEST_CASE(accepted_planar_p2_cavity_smoothing_releases_children)
{
  MetrisParameters parameters = parameters_for(
      "accepted_planar_release",
      METRIS_ROOT_DIR "/examples/2D/misc/2tri2D.mesh");
  MetrisRunner runner(nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner.degElevate(),1);
  auto& mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);
  mesh.met.setSpace(MetSpace::Exp);
  const PlanarSplitCase split = find_and_curve_interior_edge(mesh);
  MshCavity cavity(0,10,0);
  place_split_point(mesh,cavity,split,0.78);

  double initial_maximum;
  int candidate_count;
  const double initial_objective = preserved_cavity_objective<2>(
      mesh,cavity,initial_maximum,candidate_count);
  BOOST_REQUIRE_EQUAL(candidate_count,4);
  const double initial_coordinate[2] = {
      mesh.coord(cavity.ipins,0),mesh.coord(cavity.ipins,1)};

  TemporaryQuadraticConeWorkspace<2> workspace(mesh);
  BadEntHandler handler(2,100.,0.0);
  double final_objective;
  double final_maximum;
  double final_target_weight;
  const double smoothing_count = smoothCavity(
      mesh,cavity,handler,QuaFun::SizeShape,
      initial_objective,initial_maximum,candidate_count,
      final_objective,final_maximum,final_target_weight,0,1);
  BOOST_REQUIRE(smoothing_count > 0.);
  BOOST_CHECK(final_objective <= initial_objective);
  BOOST_CHECK_EQUAL(cavity.split_edge_points.get_n(),0);
  BOOST_CHECK(!cavity.preserve_split_edge_geometry);
  const double displacement = std::hypot(
      mesh.coord(cavity.ipins,0) - initial_coordinate[0],
      mesh.coord(cavity.ipins,1) - initial_coordinate[1]);
  BOOST_CHECK(displacement > 1.e-12);

  for(int source : cavity.lcfac){
    for(int opposite = 0; opposite < 3; opposite++){
      const int neighbor = mesh.fac2fac(source,opposite);
      bool neighbor_in_cavity = false;
      for(int candidate : cavity.lcfac){
        neighbor_in_cavity = neighbor_in_cavity || neighbor == candidate;
      }
      if(neighbor_in_cavity) continue;
      install_cone_vertices<2>(
          mesh,cavity,source,opposite,workspace.element());
      complete_quadratic_cone_element<MetricFieldAnalytical,2,2>(
          mesh,cavity,source,workspace.element(),
          QuadraticConeSpokePolicy::ReleasedAffine);
      BOOST_CHECK((classify_element_validity<2,2>(
          mesh,workspace.element()).is_certified()));
      for(int edge = 0; edge < 3; edge++){
        const int local0 = quadratic_simplex_edge_vertex<2>(edge,0);
        const int local1 = quadratic_simplex_edge_vertex<2>(edge,1);
        if(local0 != 0 && local1 != 0) continue;
        const int midpoint = mesh.fac2poi(
            workspace.element(),quadratic_simplex_edge_node<2>(edge));
        for(int component = 0; component < 2; component++){
          BOOST_CHECK_SMALL(
              mesh.coord(midpoint,component)
              - 0.5*(mesh.coord(mesh.fac2poi(workspace.element(),local0),component)
                   + mesh.coord(mesh.fac2poi(workspace.element(),local1),component)),
              5.e-14);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(rejected_planar_p2_cavity_smoothing_restores_provenance)
{
  MetrisParameters parameters = parameters_for(
      "rejected_planar_release",
      METRIS_ROOT_DIR "/examples/2D/misc/2tri2D.mesh");
  MetrisRunner runner(nullptr,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(runner.degElevate(),1);
  auto& mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*runner.msh_g);
  mesh.met.setSpace(MetSpace::Exp);
  const PlanarSplitCase split = find_and_curve_interior_edge(mesh);
  MshCavity cavity(0,10,0);
  place_split_point(mesh,cavity,split,0.78);

  double initialMaximum;
  int candidateCount;
  const double initialObjective = preserved_cavity_objective<2>(
      mesh,cavity,initialMaximum,candidateCount);
  const double coordinate[2] = {
      mesh.coord(cavity.ipins,0),mesh.coord(cavity.ipins,1)};
  const int splitPoints[3] = {
      cavity.split_edge_points[0],
      cavity.split_edge_points[1],
      cavity.split_edge_points[2]};

  TemporaryQuadraticConeWorkspace<2> workspace(mesh);
  // Requiring a 100% reduction makes every positive objective candidate fail
  // final acceptance even when the low-level optimizer finds a useful move.
  BadEntHandler handler(2,100.,100.0);
  double finalObjective;
  double finalMaximum;
  double finalTargetWeight;
  const double smoothingCount = smoothCavity(
      mesh,cavity,handler,QuaFun::SizeShape,
      initialObjective,initialMaximum,candidateCount,
      finalObjective,finalMaximum,finalTargetWeight,0,1);
  BOOST_CHECK_EQUAL(smoothingCount,0.);
  BOOST_CHECK_EQUAL(finalObjective,initialObjective);
  BOOST_CHECK_EQUAL(cavity.split_edge_points.get_n(),3);
  BOOST_CHECK(cavity.preserve_split_edge_geometry);
  for(int index = 0; index < 3; index++){
    BOOST_CHECK_EQUAL(cavity.split_edge_points[index],splitPoints[index]);
  }
  BOOST_CHECK_EQUAL(mesh.coord(cavity.ipins,0),coordinate[0]);
  BOOST_CHECK_EQUAL(mesh.coord(cavity.ipins,1),coordinate[1]);
}
