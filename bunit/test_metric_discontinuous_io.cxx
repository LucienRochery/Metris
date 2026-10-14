// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2026, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_metric_discontinuous_io
#include <boost/test/included/unit_test.hpp>
#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisRunner.hxx"
#include "quality/quafun_connection.hxx"
#include "quality/objective_quadrature_sample.hxx"
#include "smoothing/low_smooballdiff.hxx"
#include "smoothing/msh_smooball.hxx"
#include "low_geo/validity.hxx"
#include "io_libmeshb.hxx"
#include "nlohmann/json.hpp"

#include <array>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>

namespace Metris{
namespace{
const std::string fixture = METRIS_ROOT_DIR "/bunit/fixtures/dg_metric_v1/two-triangles";
constexpr int connectivity[2][6] = {{0,1,2,4,5,6},{1,3,2,7,4,8}};
constexpr double points[9][2] = {{0,0},{1,0},{0,1},{1,1},{.5,.5},{0,.5},{.5,0},{.5,1},{1,.5}};
constexpr double values[2][3][3] = {{{.2,.1,-.3},{.6,-.1,.4},{-.3,.4,-.1}},
                                  {{-.1,.8,.5},{.7,.4,1.4},{.4,-.1,1.5}}};
struct Support{
  MetrisParameters param;
  Mesh<MetricFieldFE> mesh;
  Support(){
    param.iverb = 0; mesh.param = &param;
    mesh.idim = 2; mesh.curdeg = mesh.strdeg = 2;
    mesh.forceBasisFlag(FEBasis::Lagrange);
    mesh.set_npoin(9); mesh.set_nentt(2,2);
    mesh.met.forceBasisFlag(FEBasis::Lagrange);
    for(int ipoin = 0; ipoin < 9; ipoin++){
      for(int dim = 0; dim < 2; dim++) mesh.coord(ipoin,dim) = points[ipoin][dim];
      mesh.met(ipoin,0) = mesh.met(ipoin,2) = 1.; mesh.met(ipoin,1) = 0.;
    }
    for(int ientt = 0; ientt < 2; ientt++){
      mesh.fac2ref[ientt] = 0;
      for(int node = 0; node < 6; node++) mesh.fac2poi(ientt,node) = connectivity[ientt][node];
    }
  }
};
struct Temporary{
  std::filesystem::path path = std::filesystem::temp_directory_path()
      / ("metris-dg-io-"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
  Temporary(){ std::filesystem::create_directory(path); }
  ~Temporary(){ std::error_code ignored; std::filesystem::remove_all(path,ignored); }
  std::string write(const std::string &text){
    const auto file = (path/"input.dgmet").string();
    std::ofstream output(file); output << text; output.close();
    BOOST_REQUIRE(output.good()); return file;
  }
};
std::string contents(const std::string &path){
  std::ifstream input(path); BOOST_REQUIRE(input.good());
  return {std::istreambuf_iterator<char>(input),std::istreambuf_iterator<char>()};
}
std::string changed(std::string text,const std::string &oldValue,const std::string &newValue){
  const auto position = text.find(oldValue); BOOST_REQUIRE(position != std::string::npos);
  text.replace(position,oldValue.size(),newValue); return text;
}

double checkCheckpointReload(MetrisRunner &run,const std::filesystem::path &base){
  auto &mesh = static_cast<Mesh<MetricFieldFE>&>(*run.msh_g);
  const auto space = mesh.met.getSpace();
  const auto basis = mesh.met.getBasis();
  const auto geometryBasis = mesh.getBasis(), outputBasis = run.param->outbasis;
  std::vector<double> coordinates(2*mesh.npoin), coefficients(9*mesh.nface), contributions(mesh.nface);
  for(int ipoin = 0; ipoin < mesh.npoin; ipoin++) for(int dim = 0; dim < 2; dim++)
    coordinates[2*ipoin+dim] = mesh.coord(ipoin,dim);
  double total = 0.;
  for(int ientt = 0; ientt < mesh.nface; ientt++){
    contributions[ientt] = metqua_connection<MetricFieldFE,2>(mesh,AsDeg::Pk,AsDeg::Pk,ientt);
    total += contributions[ientt];
    for(int node = 0; node < 3; node++) for(int k = 0; k < 3; k++)
      coefficients[9*ientt+3*node+k] = mesh.met.elementCoefficient(ientt,node)[k];
  }
  run.param->setMeshOut(base.string()+".meshb");
  run.writeOutputs();
  BOOST_CHECK(mesh.met.getSpace() == space && mesh.met.getBasis() == basis);
  BOOST_CHECK(mesh.getBasis() == geometryBasis && run.param->outbasis == outputBasis);
  for(int ipoin = 0; ipoin < mesh.npoin; ipoin++) for(int dim = 0; dim < 2; dim++)
    BOOST_CHECK_EQUAL(mesh.coord(ipoin,dim),coordinates[2*ipoin+dim]);
  for(int ientt = 0; ientt < mesh.nface; ientt++) for(int node = 0; node < 3; node++)
    for(int k = 0; k < 3; k++)
      BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(ientt,node)[k],coefficients[9*ientt+3*node+k]);
  BOOST_CHECK(!std::filesystem::exists(base.string()+".solb"));
  std::ifstream input(base.string()+".connection.json"); BOOST_REQUIRE(input.good());
  nlohmann::json record; input >> record;
  BOOST_REQUIRE_EQUAL(record.at("version").get<int>(),1);
  BOOST_REQUIRE_EQUAL(record.at("objective_name").get<std::string>(),"ConnectionNormSquared");
  const std::string measure = standard_objective_quadrature_theta_mode() == ObjectiveQuadratureTheta::PhysicalMetricMeasure
      ? "physical_metric" : "physical";
  BOOST_CHECK_EQUAL(record.at("integration_measure").get<std::string>(),measure);
  BOOST_CHECK_EQUAL(record.at("value").get<double>(),total);
  MetrisParameters param; param.iverb = 0; param.usrTarDeg = 2;
  param.setMeshIn((base.parent_path()/record.at("mesh_file").get<std::string>()).string());
  param.dgMetFileName = (base.parent_path()/record.at("met_file").get<std::string>()).string();
  param.objective_quadrature_order = record.at("quadrature_order").get<int>();
  BOOST_CHECK_EQUAL(param.objective_quadrature_order,run.param->objective_quadrature_order);
  MetrisRunner reloaded(nullptr,nullptr,param);
  auto &restored = static_cast<Mesh<MetricFieldFE>&>(*reloaded.msh_g);
  BOOST_REQUIRE_EQUAL(restored.npoin,mesh.npoin);
  BOOST_REQUIRE_EQUAL(restored.nface,mesh.nface);
  double restoredTotal = 0.;
  const double bary[3] = {.2,.3,.5};
  for(int ientt = 0; ientt < mesh.nface; ientt++){
    for(int node = 0; node < 6; node++) BOOST_CHECK_EQUAL(restored.fac2poi(ientt,node),mesh.fac2poi(ientt,node));
    for(int node = 0; node < 3; node++){
      double logarithm[3]; std::copy(mesh.met.elementCoefficient(ientt,node),mesh.met.elementCoefficient(ientt,node)+3,logarithm);
      if(space == MetSpace::Exp) getlogmet_inp<2,double>(logarithm);
      for(int k = 0; k < 3; k++) BOOST_CHECK_EQUAL(restored.met.elementCoefficient(ientt,node)[k],logarithm[k]);
    }
    double met[3], dmet[6], restoredMet[3], restoredD[6];
    mesh.met.getMetOnEntity(ientt,DifVar::Phys,bary,met,dmet);
    restored.met.getMetOnEntity(ientt,DifVar::Phys,bary,restoredMet,restoredD);
    for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(restoredMet[k]-met[k],2e-12*(1.+std::abs(met[k])));
    for(int k = 0; k < 6; k++) BOOST_CHECK_SMALL(restoredD[k]-dmet[k],2e-12*(1.+std::abs(dmet[k])));
    const double value = metqua_connection<MetricFieldFE,2>(restored,AsDeg::Pk,AsDeg::Pk,ientt);
    BOOST_CHECK_EQUAL(record.at("entities").at(ientt).at("entity").get<int>(),ientt+1);
    BOOST_CHECK_EQUAL(record.at("entities").at(ientt).at("value").get<double>(),contributions[ientt]);
    BOOST_CHECK_SMALL(value-contributions[ientt],2e-12*(1.+std::abs(contributions[ientt])));
    restoredTotal += value;
  }
  BOOST_CHECK_SMALL(restoredTotal-total,2e-12*(1.+std::abs(total)));
  return total;
}
}

BOOST_FIXTURE_TEST_CASE(import_preserves_log_coefficients_and_interface_traces,Support){
  mesh.met.readDiscontinuousMetricFile(fixture+".dgmet");
  BOOST_CHECK_EQUAL(mesh.met.getCoefficientCount(),6);
  for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++)
    for(int k = 0; k < 3; k++) BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(ientt,node)[k],values[ientt][node][k]);
  BOOST_CHECK(mesh.met.elementCoefficient(0,1) != mesh.met.elementCoefficient(1,0));
  const double bary[2][3] = {{0,.4,.6},{.4,0,.6}};
  double traces[2][3], derivatives[6];
  for(int ientt = 0; ientt < 2; ientt++)
    mesh.met.getMetOnEntity(ientt,DifVar::Phys,bary[ientt],traces[ientt],derivatives);
  BOOST_CHECK_GT(std::abs(traces[0][0]-traces[1][0]),1e-3);
}

BOOST_AUTO_TEST_CASE(import_remaps_entity_rows_corners_and_explicit_point_ids){
  std::array<int,3> permutation = {{0,1,2}};
  do{
    for(bool renumber : {false,true}){
      Support support; auto &mesh = support.mesh;
      std::vector<int> map(9);
      for(int ipoin = 0; ipoin < 9; ipoin++){
        map[ipoin] = renumber ? 8-ipoin : ipoin;
        for(int dim = 0; dim < 2; dim++) mesh.coord(map[ipoin],dim) = points[ipoin][dim];
      }
      for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++){
        mesh.fac2poi(ientt,node) = map[connectivity[1-ientt][permutation[node]]];
        mesh.fac2poi(ientt,3+node) = map[connectivity[1-ientt][3+permutation[node]]];
      }
      mesh.met.readDiscontinuousMetricFile(fixture+".dgmet",map);
      for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++)
        for(int k = 0; k < 3; k++)
          BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(ientt,node)[k],values[1-ientt][permutation[node]][k]);
    }
  }while(std::next_permutation(permutation.begin(),permutation.end()));
}

BOOST_FIXTURE_TEST_CASE(import_spd_and_bernstein_p1_reproduces_log_interpolant,Support){
  mesh.met.readDiscontinuousMetricFile(fixture+".dgmet");
  std::string text = contents(fixture+".dgmet");
  text = changed(text,"MetricSpace Log","MetricSpace SPD");
  text = changed(text,"MetricBasis Lagrange","MetricBasis Bernstein");
  text.erase(text.find("Coefficients 6")+std::string("Coefficients 6\n").size());
  std::ostringstream records; records << std::setprecision(17);
  for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++){
    double met[3]; getexpmet_frechet<2>(values[ientt][node],nullptr,0,met,nullptr);
    records << ientt+1 << ' ' << node+1 << ' ' << connectivity[ientt][node]+1
            << ' ' << met[0] << ' ' << met[1] << ' ' << met[2] << '\n';
  }
  Temporary files;
  mesh.met.readDiscontinuousMetricFile(files.write(text+records.str()+"End\n"));
  BOOST_CHECK(mesh.met.getBasis() == FEBasis::Bezier);
  BOOST_CHECK(mesh.met.getSpace() == MetSpace::Log);
  for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++)
    for(int k = 0; k < 3; k++)
      BOOST_CHECK_SMALL(mesh.met.elementCoefficient(ientt,node)[k]-values[ientt][node][k],2e-14);
}

BOOST_FIXTURE_TEST_CASE(invalid_import_is_transactional,Support){
  Temporary files;
  const auto valid = contents(fixture+".dgmet");
  const std::vector<std::string> invalid = {
    changed(valid,"MetrisDGMetric 1","MetrisDGMetric 2"),
    changed(valid,"MetricDegree 1","MetricDegree 2"),
    changed(valid,"TensorComponents xx xy yy","TensorComponents xx yy xy"),
    changed(valid,"5 0.5 0.5","5 0.6 0.5"),
    changed(valid,"1 1 2 3 5 6 7 1","1 1 2 3 8 6 7 1"),
    changed(valid,"1 1 2 3 5 6 7 1","1 1 2 3 5 6 7 2"),
    changed(valid,"2 2 4 3 8 5 9 1","1 2 4 3 8 5 9 1"),
    changed(valid,"1 2 2 0.6 -0.1 0.4","1 1 1 0.6 -0.1 0.4"),
    changed(valid,"1 2 2 0.6 -0.1 0.4","1 2 3 0.6 -0.1 0.4"),
    changed(valid,"1 1 1 0.2 0.1 -0.3","1 1 1 1000 0.1 -0.3"),
    changed(valid,"End","Extra\nEnd"), valid+"Extra\n",
    changed(changed(valid,"MetricSpace Log","MetricSpace SPD"),"1 1 1 0.2 0.1 -0.3","1 1 1 -1 0 -1")
  };
  for(const auto &text : invalid){
    BOOST_CHECK_THROW(mesh.met.readDiscontinuousMetricFile(files.write(text)),MetrisExcept);
    BOOST_CHECK(mesh.met.getLayout() == MetricFieldLayout::Continuous);
    BOOST_CHECK_EQUAL(mesh.met(0,0),1.);
  }
  mesh.met.readDiscontinuousMetricFile(fixture+".dgmet");
  for(const auto &text : invalid){
    BOOST_CHECK_THROW(mesh.met.readDiscontinuousMetricFile(files.write(text)),MetrisExcept);
    for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++)
      for(int k = 0; k < 3; k++) BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(ientt,node)[k],values[ientt][node][k]);
  }
}

BOOST_AUTO_TEST_CASE(runner_copies_dg_background_without_interpolation){
  MetrisParameters param; param.iverb = 0; param.usrTarDeg = 2;
  param.setMeshIn(fixture+".mesh"); param.dgMetFileName = fixture+".dgmet";
  MetrisRunner run(nullptr,nullptr,param);
  auto &mesh = static_cast<Mesh<MetricFieldFE>&>(*run.msh_g);
  BOOST_CHECK(mesh.met.getLayout() == MetricFieldLayout::Discontinuous);
  BOOST_CHECK_EQUAL(mesh.npoin,run.bak.npoin);
  for(int ipoin = 0; ipoin < mesh.npoin; ipoin++) for(int dim = 0; dim < 2; dim++)
    BOOST_CHECK_EQUAL(mesh.coord(ipoin,dim),run.bak.coord(ipoin,dim));
  for(int ientt = 0; ientt < mesh.nface; ientt++){
    for(int node = 0; node < 6; node++) BOOST_CHECK_EQUAL(mesh.fac2poi(ientt,node),run.bak.fac2poi(ientt,node));
    for(int node = 0; node < 3; node++) for(int k = 0; k < 3; k++){
      BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(ientt,node)[k],values[ientt][node][k]);
      BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(ientt,node)[k],run.bak.met.elementCoefficient(ientt,node)[k]);
    }
  }
  mesh.met.elementCoefficient(0,1)[0] += 1.;
  BOOST_CHECK_EQUAL(run.bak.met.elementCoefficient(0,1)[0],values[0][1][0]);
}

BOOST_AUTO_TEST_CASE(background_samples_preserve_entity_specific_interface_traces){
  for(MetSpace space : {MetSpace::Log,MetSpace::Exp}){
    MetrisParameters param; param.iverb = 0; param.usrTarDeg = 2;
    param.setMeshIn(fixture+".mesh"); param.dgMetFileName = fixture+".dgmet";
    MetrisRunner run(nullptr,nullptr,param);
    auto &mesh = static_cast<Mesh<MetricFieldFE>&>(*run.msh_g);
    mesh.met.setSpace(space);
    std::array<double,18> originalCoefficients;
    for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++)
      for(int k = 0; k < 3; k++)
        originalCoefficients[9*ientt+3*node+k] = mesh.met.elementCoefficient(ientt,node)[k];
    for(int ipoin = 0; ipoin < mesh.npoin; ipoin++) for(int k = 0; k < 3; k++){
      mesh.met.rfld(ipoin,k) = std::numeric_limits<double>::quiet_NaN();
      run.bak.met.rfld(ipoin,k) = std::numeric_limits<double>::quiet_NaN();
    }
    const int oldSeed = mesh.poi2bak[4];
    const double bary[2][3] = {{0,.5,.5},{.5,0,.5}};
    BackgroundMetricSample traces[2];
    for(int ientt = 0; ientt < 2; ientt++){
      double expected[3];
      if(space == MetSpace::Log) run.bak.met.getLogMetOnEntity(ientt,DifVar::None,bary[ientt],expected,nullptr);
      else run.bak.met.getMetOnEntity(ientt,DifVar::None,bary[ientt],expected,nullptr);
      BOOST_REQUIRE_EQUAL(mesh.sampleMetBackOnEntity(4,ientt,traces[ientt]),0);
      BOOST_CHECK_EQUAL(traces[ientt].location.ientt,ientt);
      BOOST_CHECK_EQUAL(traces[ientt].ienttMet,ientt);
      BOOST_CHECK(traces[ientt].space == space);
      BackgroundMetricSample automatic;
      BOOST_REQUIRE_EQUAL(mesh.sampleMetBack(4,automatic,ientt),0);
      for(int k = 0; k < 3; k++){
        BOOST_CHECK_SMALL(traces[ientt].met[k]-expected[k],2e-13);
        BOOST_CHECK_EQUAL(automatic.met[k],traces[ientt].met[k]);
      }
    }
    for(int ipoin : {1,2}) for(int ientt = 0; ientt < 2; ientt++){
      const int node = ipoin == 1 ? (ientt == 0 ? 1 : 0) : 2;
      double vertexBary[3] = {}; vertexBary[node] = 1.;
      double expected[3];
      if(space == MetSpace::Log) run.bak.met.getLogMetOnEntity(ientt,DifVar::None,vertexBary,expected,nullptr);
      else run.bak.met.getMetOnEntity(ientt,DifVar::None,vertexBary,expected,nullptr);
      BackgroundMetricSample vertex;
      const int seed = mesh.poi2bak[ipoin];
      BOOST_REQUIRE_EQUAL(mesh.sampleMetBackOnEntity(ipoin,ientt,vertex),0);
      BOOST_CHECK_EQUAL(vertex.ienttMet,ientt);
      for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(vertex.met[k]-expected[k],2e-13);
      BOOST_CHECK_EQUAL(mesh.poi2bak[ipoin],seed);
    }
    BOOST_CHECK_GT(std::abs(traces[0].met[2]-traces[1].met[2]),.1);
    BOOST_CHECK_EQUAL(mesh.poi2bak[4],oldSeed);
    BOOST_CHECK_THROW(mesh.sampleMetBack(4,traces[0]),MetrisExcept);
    BOOST_CHECK_THROW(mesh.interpMetBack(4),MetrisExcept);
    const double bottomBary[3] = {.5,.5,0};
    double expected[3];
    if(space == MetSpace::Log) run.bak.met.getLogMetOnEntity(0,DifVar::None,bottomBary,expected,nullptr);
    else run.bak.met.getMetOnEntity(0,DifVar::None,bottomBary,expected,nullptr);
    BackgroundMetricSample boundary;
    const int oldBoundarySeed = mesh.poi2bak[6];
    BOOST_REQUIRE_EQUAL(mesh.sampleMetBack(6,boundary,0),0);
    BOOST_CHECK_EQUAL(boundary.location.tdim,1);
    BOOST_CHECK_EQUAL(boundary.tdimMet,2);
    BOOST_CHECK_EQUAL(boundary.ienttMet,0);
    for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(boundary.met[k]-expected[k],2e-13);
    BOOST_REQUIRE_EQUAL(mesh.sampleMetBack(6,boundary,1),0);
    BOOST_CHECK_EQUAL(boundary.ienttMet,0);
    for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(boundary.met[k]-expected[k],2e-13);
    BOOST_CHECK_EQUAL(mesh.poi2bak[6],oldBoundarySeed);
    for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++)
      for(int k = 0; k < 3; k++){
        BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(ientt,node)[k],originalCoefficients[9*ientt+3*node+k]);
        BOOST_CHECK_EQUAL(run.bak.met.elementCoefficient(ientt,node)[k],values[ientt][node][k]);
      }
    for(int ipoin = 0; ipoin < mesh.npoin; ipoin++) for(int dim = 0; dim < 2; dim++){
      BOOST_CHECK_EQUAL(mesh.coord(ipoin,dim),points[ipoin][dim]);
      BOOST_CHECK_EQUAL(run.bak.coord(ipoin,dim),points[ipoin][dim]);
    }
  }
}

BOOST_AUTO_TEST_CASE(background_sample_localizes_moved_points_and_curved_maps_transactionally){
  MetrisParameters param; param.iverb = 0; param.usrTarDeg = 2;
  param.setMeshIn(fixture+".mesh"); param.dgMetFileName = fixture+".dgmet";
  MetrisRunner run(nullptr,nullptr,param);
  auto &mesh = static_cast<Mesh<MetricFieldFE>&>(*run.msh_g);
  const int oldSeed = mesh.poi2bak[4];
  mesh.coord(4,0) = .8; mesh.coord(4,1) = .7;
  BackgroundMetricSample sample;
  BOOST_REQUIRE_EQUAL(mesh.sampleMetBackOnEntity(4,0,sample),0);
  BOOST_CHECK_EQUAL(sample.ienttMet,1);
  const double crossedBary[3] = {.3,.5,.2};
  double expected[3]; run.bak.met.getLogMetOnEntity(1,DifVar::None,crossedBary,expected,nullptr);
  for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(sample.met[k]-expected[k],2e-13);
  run.bak.setBasis(FEBasis::Lagrange);
  run.bak.coord(4,0) += .1; run.bak.coord(4,1) -= .07;
  run.bak.setBasis(FEBasis::Bezier);
  const double curvedBary[3] = {.2,.3,.5};
  double point[2];
  eval2<2,2>(run.bak.coord,run.bak.fac2poi[0],run.bak.getBasis(),DifVar::None,DifVar::None,
             curvedBary,point,nullptr,nullptr);
  std::copy(point,point+2,mesh.coord[4]);
  BOOST_REQUIRE_EQUAL(mesh.sampleMetBackOnEntity(4,0,sample),0);
  BOOST_CHECK_EQUAL(sample.ienttMet,0);
  run.bak.met.getLogMetOnEntity(0,DifVar::None,curvedBary,expected,nullptr);
  for(int k = 0; k < 3; k++){
    BOOST_CHECK_SMALL(sample.met[k]-expected[k],2e-6);
    BOOST_CHECK_SMALL(sample.baryMet[k]-curvedBary[k],2e-6);
  }
  const auto previous = sample;
  mesh.coord(4,0) = 20.; mesh.coord(4,1) = 20.;
  BOOST_CHECK_NE(mesh.sampleMetBackOnEntity(4,0,sample),0);
  BOOST_CHECK(sample.met == previous.met);
  for(int k = 0; k < 3; k++) BOOST_CHECK_EQUAL(sample.baryMet[k],previous.baryMet[k]);
  BOOST_CHECK_EQUAL(sample.ienttMet,previous.ienttMet);
  BOOST_CHECK_EQUAL(mesh.poi2bak[4],oldSeed);
  for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++)
    for(int k = 0; k < 3; k++){
      BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(ientt,node)[k],values[ientt][node][k]);
      BOOST_CHECK_EQUAL(run.bak.met.elementCoefficient(ientt,node)[k],values[ientt][node][k]);
    }
}

BOOST_AUTO_TEST_CASE(dg_background_owner_fallback_is_stable_and_reference_constrained){
  MetrisParameters param; param.iverb = 0; param.usrTarDeg = 2;
  param.setMeshIn(fixture+".mesh"); param.dgMetFileName = fixture+".dgmet";
  MetrisRunner run(nullptr,nullptr,param);
  auto &mesh = static_cast<Mesh<MetricFieldFE>&>(*run.msh_g);
  for(FEBasis basis : {FEBasis::Lagrange,FEBasis::Bezier}){
    run.bak.setBasis(FEBasis::Lagrange);
    run.bak.coord(4,0) = .6; run.bak.coord(4,1) = .43;
    run.bak.setBasis(basis);
    for(double t : {0.,.23,.5,.79,1.}){
      const double bary[2][3] = {{0,1.-t,t},{1.-t,0,t}};
      double point[2];
      eval2<2,2>(run.bak.coord,run.bak.fac2poi[0],basis,DifVar::None,DifVar::None,
                 bary[0],point,nullptr,nullptr);
      std::copy(point,point+2,mesh.coord[4]);
      for(int repetition = 0; repetition < 3; repetition++){
        for(int ientt : {1,0,run.bak.nface+3,-1}){
          BackgroundMetricSample sample;
          BOOST_REQUIRE_EQUAL(mesh.sampleMetBackOnEntity(4,ientt,sample),0);
          const int expectedOwner = ientt == 1 ? 1 : 0;
          BOOST_CHECK_EQUAL(sample.ienttMet,expectedOwner);
          double expected[3];
          run.bak.met.getLogMetOnEntity(expectedOwner,DifVar::None,bary[expectedOwner],expected,nullptr);
          for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(sample.met[k]-expected[k],2e-10);
        }
      }
    }
  }
  run.bak.setBasis(FEBasis::Lagrange);
  run.bak.coord(4,0) = run.bak.coord(4,1) = .5;
  run.bak.setBasis(FEBasis::Bezier);
  for(double side : {-1.,1.}){
    mesh.coord(4,0) = .5+side*1e-7; mesh.coord(4,1) = .5;
    for(int ientt : {0,1,99}){
      BackgroundMetricSample sample;
      BOOST_REQUIRE_EQUAL(mesh.sampleMetBackOnEntity(4,ientt,sample),0);
      BOOST_CHECK_EQUAL(sample.ienttMet,side < 0 ? 0 : 1);
    }
  }
  mesh.coord(4,0) = mesh.coord(4,1) = .5;
  const int oldRef = run.bak.fac2ref[0];
  run.bak.fac2ref[1] = oldRef+1;
  for(int ientt : {-1,0,1,99}){
    BackgroundLocation location;
    BOOST_REQUIRE_EQUAL(mesh.locateMetBack(4,2,ientt,oldRef+1,nullptr,location),0);
    BOOST_CHECK_EQUAL(location.ientt,1);
  }
  BackgroundLocation unchanged; unchanged.ientt = 91;
  BOOST_CHECK_NE(mesh.locateMetBack(4,2,0,oldRef+2,nullptr,unchanged),0);
  BOOST_CHECK_EQUAL(unchanged.ientt,91);
}

BOOST_AUTO_TEST_CASE(dg_interface_tie_break_is_independent_of_a_stale_live_owner){
  Temporary files;
  auto geometry = contents(fixture+".mesh");
  geometry = changed(geometry,"Vertices\n9\n","Vertices\n12\n");
  geometry = changed(geometry,"1 0.5 0\nEdgesP2Ordering","1 0.5 0\n2 0 0\n1.5 0 0\n1.5 0.5 0\nEdgesP2Ordering");
  geometry = changed(geometry,"EdgesP2\n4\n","EdgesP2\n5\n");
  geometry = changed(geometry,"2 4 8 2\n","2 10 11 2\n10 4 12 2\n");
  geometry = changed(geometry,"Corners\n4\n","Corners\n5\n10\n");
  geometry = changed(geometry,"TrianglesP2\n2\n","TrianglesP2\n3\n");
  geometry = changed(geometry,"End","2 10 4 12 9 11 1\nEnd");
  const auto meshFile = files.path/"three-triangles.mesh";
  std::ofstream output(meshFile); output << geometry; output.close(); BOOST_REQUIRE(output.good());
  auto met = contents(fixture+".dgmet");
  met = changed(met,"Points 9","Points 12");
  met = changed(met,"Entities 2","10 2 0\n11 1.5 0\n12 1.5 0.5\nEntities 3");
  met = changed(met,"Coefficients 6","3 2 10 4 12 9 11 1\nCoefficients 9");
  met = changed(met,"End","3 1 2 0 0 0\n3 2 10 0 0 0\n3 3 4 0 0 0\nEnd");
  MetrisParameters param; param.iverb = 0; param.usrTarDeg = 2;
  param.setMeshIn(meshFile.string()); param.dgMetFileName = files.write(met);
  MetrisRunner run(nullptr,nullptr,param);
  auto &mesh = static_cast<Mesh<MetricFieldFE>&>(*run.msh_g);
  BOOST_REQUIRE_EQUAL(run.bak.nface,3);
  for(int repetition = 0; repetition < 3; repetition++) for(int ientt : {2,-1,99,1,0}){
    BackgroundMetricSample sample;
    BOOST_REQUIRE_EQUAL(mesh.sampleMetBackOnEntity(4,ientt,sample),0);
    BOOST_CHECK_EQUAL(sample.ienttMet,ientt == 1 ? 1 : 0);
  }
}

BOOST_AUTO_TEST_CASE(dg_boundary_owner_follows_the_located_reference_and_preserves_vertex_traces){
  MetrisParameters param; param.iverb = 0; param.usrTarDeg = 2;
  param.setMeshIn(fixture+".mesh"); param.dgMetFileName = fixture+".dgmet";
  MetrisRunner run(nullptr,nullptr,param);
  auto &mesh = static_cast<Mesh<MetricFieldFE>&>(*run.msh_g);
  const int oldSeed = mesh.poi2bak[6];
  for(double x : {.2,.6,.9}){
    mesh.coord(6,0) = x;
    for(int ientt : {1,0,99}){
      BackgroundMetricSample sample;
      BOOST_REQUIRE_EQUAL(mesh.sampleMetBack(6,sample,ientt),0);
      BOOST_CHECK_EQUAL(sample.location.tdim,1);
      BOOST_CHECK_EQUAL(sample.ienttMet,0);
      const double bary[3] = {1.-x,x,0.};
      double expected[3]; run.bak.met.getLogMetOnEntity(0,DifVar::None,bary,expected,nullptr);
      for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(sample.met[k]-expected[k],2e-13);
    }
  }
  mesh.coord(6,0) = 1.;
  for(int ientt : {1,0,99}){
    BackgroundMetricSample sample;
    BOOST_REQUIRE_EQUAL(mesh.sampleMetBack(6,sample,ientt),0);
    const int expectedOwner = ientt == 1 ? 1 : 0;
    const int node = expectedOwner == 0 ? 1 : 0;
    BOOST_CHECK_EQUAL(sample.ienttMet,expectedOwner);
    for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(sample.met[k]-values[expectedOwner][node][k],2e-13);
  }
  BOOST_CHECK_EQUAL(mesh.poi2bak[6],oldSeed);
  mesh.coord(6,0) = .5;
  BackgroundMetricSample sample;
  const int ienttBottom = oldSeed;
  const auto edgeNode = run.bak.edg2poi(ienttBottom,2);
  run.bak.edg2poi(ienttBottom,2) = run.bak.edg2poi((ienttBottom+1)%run.bak.nedge,2);
  sample.ienttMet = 91; sample.met.fill(7.);
  BOOST_CHECK_NE(mesh.sampleMetBack(6,sample,0),0);
  BOOST_CHECK_EQUAL(sample.ienttMet,91);
  for(double value : sample.met) BOOST_CHECK_EQUAL(value,7.);
  run.bak.edg2poi(ienttBottom,2) = edgeNode;
  for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++)
    for(int k = 0; k < 3; k++){
      BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(ientt,node)[k],values[ientt][node][k]);
      BOOST_CHECK_EQUAL(run.bak.met.elementCoefficient(ientt,node)[k],values[ientt][node][k]);
    }
}

BOOST_AUTO_TEST_CASE(dg_vertex_refresh_preserves_traces_and_updates_every_incident_copy){
  for(MetSpace space : {MetSpace::Log,MetSpace::Exp}){
    MetrisParameters param; param.iverb = 0; param.usrTarDeg = 2;
    param.setMeshIn(fixture+".mesh"); param.dgMetFileName = fixture+".dgmet";
    MetrisRunner run(nullptr,nullptr,param);
    auto &mesh = static_cast<Mesh<MetricFieldFE>&>(*run.msh_g);
    mesh.met.setSpace(space);
    std::array<double,18> original;
    for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++)
      for(int k = 0; k < 3; k++) original[9*ientt+3*node+k] = mesh.met.elementCoefficient(ientt,node)[k];
    const auto initial = mesh.captureDGMetAtVertex(1);
    BOOST_REQUIRE_EQUAL(initial.size(),2);
    for(const auto &coefficient : initial) BOOST_CHECK_EQUAL(coefficient.ienttBack,coefficient.ientt);
    for(int ipoin = 0; ipoin < mesh.npoin; ipoin++) for(int k = 0; k < 3; k++){
      mesh.met.rfld(ipoin,k) = std::numeric_limits<double>::quiet_NaN();
      run.bak.met.rfld(ipoin,k) = std::numeric_limits<double>::quiet_NaN();
    }
    mesh.coord(1,0) = .6; mesh.coord(1,1) = .4;
    const int sharedSeed = mesh.poi2bak[1];
    BOOST_REQUIRE_EQUAL(mesh.refreshDGMetAtVertex(1),0);
    const double bary[2][3] = {{0,.6,.4},{.6,0,.4}};
    for(const auto &coefficient : mesh.captureDGMetAtVertex(1)){
      double expected[3] = {};
      for(int node = 0; node < 3; node++) for(int k = 0; k < 3; k++)
        expected[k] += bary[coefficient.ientt][node]*values[coefficient.ientt][node][k];
      if(space == MetSpace::Exp) getexpmet_inp<2,double>(expected);
      BOOST_CHECK_EQUAL(coefficient.ienttBack,coefficient.ientt);
      for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(coefficient.met[k]-expected[k],2e-13);
    }
    BOOST_CHECK_GT(std::abs(mesh.met.elementCoefficient(0,1)[2]-mesh.met.elementCoefficient(1,0)[2]),.1);
    mesh.coord(1,0) = .9; mesh.coord(1,1) = .4;
    BOOST_REQUIRE_EQUAL(mesh.refreshDGMetAtVertex(1),0);
    double expected[3] = {};
    const double inside[3] = {.6,.3,.1};
    for(int node = 0; node < 3; node++) for(int k = 0; k < 3; k++) expected[k] += inside[node]*values[1][node][k];
    if(space == MetSpace::Exp) getexpmet_inp<2,double>(expected);
    for(const auto &coefficient : mesh.captureDGMetAtVertex(1)){
      BOOST_CHECK_EQUAL(coefficient.ienttBack,1);
      for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(coefficient.met[k]-expected[k],2e-13);
    }
    for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++){
      if(connectivity[ientt][node] == 1) continue;
      for(int k = 0; k < 3; k++) BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(ientt,node)[k],original[9*ientt+3*node+k]);
      BOOST_CHECK_EQUAL(mesh.met2bak[mesh.met.getElementCoefficientIndex(ientt,node)],ientt);
    }
    BOOST_CHECK_EQUAL(mesh.poi2bak[1],sharedSeed);
    mesh.restoreDGMetAtVertex(initial);
    for(const auto &coefficient : initial){
      const auto restored = mesh.met.elementCoefficient(coefficient.ientt,coefficient.node);
      for(int k = 0; k < 3; k++) BOOST_CHECK_EQUAL(restored[k],coefficient.met[k]);
      BOOST_CHECK_EQUAL(mesh.met2bak[mesh.met.getElementCoefficientIndex(coefficient.ientt,coefficient.node)],coefficient.ienttBack);
    }
  }
}

BOOST_AUTO_TEST_CASE(dg_vertex_refresh_commits_only_after_all_samples_succeed){
  MetrisParameters param; param.iverb = 0; param.usrTarDeg = 2;
  param.setMeshIn(fixture+".mesh"); param.dgMetFileName = fixture+".dgmet";
  MetrisRunner run(nullptr,nullptr,param);
  auto &mesh = static_cast<Mesh<MetricFieldFE>&>(*run.msh_g);
  const auto initial = mesh.captureDGMetAtVertex(1);
  mesh.coord(1,0) = .6; mesh.coord(1,1) = .4;
  const double old = run.bak.met.elementCoefficient(1,0)[0];
  run.bak.met.elementCoefficient(1,0)[0] = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK_THROW(mesh.refreshDGMetAtVertex(1),MetrisExcept);
  run.bak.met.elementCoefficient(1,0)[0] = old;
  mesh.coord(1,0) = mesh.coord(1,1) = 20.;
  BOOST_CHECK_NE(mesh.refreshDGMetAtVertex(1),0);
  const auto current = mesh.captureDGMetAtVertex(1);
  for(size_t i = 0; i < initial.size(); i++){
    BOOST_CHECK(current[i].met == initial[i].met);
    BOOST_CHECK_EQUAL(current[i].ienttBack,initial[i].ienttBack);
  }
}

BOOST_AUTO_TEST_CASE(dg_vertex_smoothing_checks_the_refreshed_candidate_objective){
  Temporary files;
  std::vector<std::array<double,2>> coordinates = {{{0,0}},{{1,0}},{{1,1}},{{0,1}},{{.5,.5}}};
  const int corners[4][3] = {{0,1,4},{1,2,4},{2,3,4},{3,0,4}};
  int nodes[4][6];
  std::map<std::pair<int,int>,int> middle;
  for(int ientt = 0; ientt < 4; ientt++){
    for(int node = 0; node < 3; node++) nodes[ientt][node] = corners[ientt][node];
    for(int edge = 0; edge < 3; edge++){
      const int first = corners[ientt][(edge+1)%3], second = corners[ientt][(edge+2)%3];
      const auto key = std::minmax(first,second);
      auto position = middle.find(key);
      if(position == middle.end()){
        const int ipoin = int(coordinates.size());
        coordinates.push_back({{.5*(coordinates[first][0]+coordinates[second][0]),
                                .5*(coordinates[first][1]+coordinates[second][1])}});
        position = middle.emplace(key,ipoin).first;
      }
      nodes[ientt][3+edge] = position->second;
    }
  }
  std::ostringstream geometry, met; geometry << std::setprecision(17); met << std::setprecision(17);
  geometry << "MeshVersionFormatted 2\nDimension 2\nVertices\n" << coordinates.size() << '\n';
  met << "MetrisDGMetric 1\nDimension 2\nTopology Triangle\nGeometryDegree 2\nGeometryBasis Lagrange\n"
         "GeometryNodeOrder v0 v1 v2 e12 e20 e01\nMetricDegree 1\nMetricBasis Lagrange\nMetricSpace Log\n"
         "TensorComponents xx xy yy\nIndexBase 1\nPoints " << coordinates.size() << '\n';
  for(size_t ipoin = 0; ipoin < coordinates.size(); ipoin++){
    geometry << coordinates[ipoin][0] << ' ' << coordinates[ipoin][1] << " 0\n";
    met << ipoin+1 << ' ' << coordinates[ipoin][0] << ' ' << coordinates[ipoin][1] << '\n';
  }
  geometry << "EdgesP2Ordering\n3\n2 0\n0 2\n1 1\nEdgesP2\n4\n";
  for(int ientt = 0; ientt < 4; ientt++) geometry << corners[ientt][0]+1 << ' ' << corners[ientt][1]+1 << ' ' << nodes[ientt][5]+1 << ' ' << ientt+1 << '\n';
  geometry << "Corners\n4\n1\n2\n3\n4\nTrianglesP2Ordering\n6\n2 0 0\n0 2 0\n0 0 2\n0 1 1\n1 0 1\n1 1 0\nTrianglesP2\n4\n";
  met << "Entities 4\n";
  for(int ientt = 0; ientt < 4; ientt++){
    met << ientt+1 << ' ';
    for(int node = 0; node < 6; node++){
      geometry << nodes[ientt][node]+1 << ' '; met << nodes[ientt][node]+1 << ' ';
    }
    geometry << "1\n"; met << "1\n";
  }
  geometry << "End\n"; met << "Coefficients 12\n";
  for(int ientt = 0; ientt < 4; ientt++) for(int node = 0; node < 3; node++){
    const int ipoin = nodes[ientt][node];
    const double x = coordinates[ipoin][0], y = coordinates[ipoin][1];
    met << ientt+1 << ' ' << node+1 << ' ' << ipoin+1 << ' ' << .1+.03*x << ' ' << .005*y << ' ' << -.1+.02*y << '\n';
  }
  met << "End\n";
  const auto meshFile = files.path/"fan.mesh";
  std::ofstream output(meshFile); output << geometry.str(); output.close(); BOOST_REQUIRE(output.good());
  MetrisParameters param; param.iverb = 0; param.usrTarDeg = 2; param.iflag1 = 1;
  param.setMeshIn(meshFile.string()); param.dgMetFileName = files.write(met.str());
  param.objective_quadrature_order = 5;
  MetrisRunner run(nullptr,nullptr,param);
  auto &mesh = static_cast<Mesh<MetricFieldFE>&>(*run.msh_g);
  mesh.setBasis(FEBasis::Lagrange);
  mesh.coord(4,0) = .54; mesh.coord(4,1) = .47;
  for(int ipoin = 0; ipoin < mesh.npoin; ipoin++) for(int k = 0; k < 3; k++) mesh.met.rfld(ipoin,k) = std::numeric_limits<double>::quiet_NaN();
  intAr1 region(4); region.set_n(0);
  double before = 0.;
  for(int ientt = 0; ientt < 4; ientt++){
    region.stack(ientt); before += metqua_connection<MetricFieldFE,2>(mesh,AsDeg::Pk,AsDeg::Pk,ientt);
  }
  double initial, maximum0, final, maximum1;
  const int status = smooballdiff<MetricFieldFE,2,2>(mesh,4,region,&initial,&maximum0,&final,&maximum1,QuaFun::Connection);
  BOOST_REQUIRE_EQUAL(status,0);
  BOOST_CHECK_SMALL(initial-before,2e-12);
  BOOST_CHECK_LT(final,before);
  double actual = 0.;
  for(int ientt = 0; ientt < 4; ientt++){
    actual += metqua_connection<MetricFieldFE,2>(mesh,AsDeg::Pk,AsDeg::Pk,ientt);
    BOOST_CHECK((classify_element_validity<2,2>(mesh,ientt).accepted_conservatively()));
    for(int node = 0; node < 3; node++){
      const int ipoin = nodes[ientt][node];
      const double x = mesh.coord(ipoin,0), y = mesh.coord(ipoin,1);
      const double expected[3] = {.1+.03*x,.005*y,-.1+.02*y};
      for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(mesh.met.elementCoefficient(ientt,node)[k]-expected[k],2e-12);
    }
  }
  BOOST_CHECK_SMALL(actual-final,2e-12);
  BOOST_CHECK_EQUAL(run.bak.coord(4,0),.5); BOOST_CHECK_EQUAL(run.bak.coord(4,1),.5);
  mesh.coord(4,0) = mesh.coord(4,1) = .5;
  for(int ientt = 0; ientt < 4; ientt++) for(int node = 0; node < 3; node++){
    double *coefficient = mesh.met.elementCoefficient(ientt,node);
    coefficient[0] = .02*ientt; coefficient[1] = 0.; coefficient[2] = -.02*ientt;
  }
  const auto unchanged = mesh.captureDGMetAtVertex(4);
  BOOST_REQUIRE_EQUAL((smooballdiff<MetricFieldFE,2,2>(mesh,4,region,&initial,&maximum0,&final,&maximum1,QuaFun::Connection)),0);
  const auto stillUnchanged = mesh.captureDGMetAtVertex(4);
  for(size_t i = 0; i < unchanged.size(); i++){
    BOOST_CHECK(unchanged[i].met == stillUnchanged[i].met);
    BOOST_CHECK_EQUAL(unchanged[i].ienttBack,stillUnchanged[i].ienttBack);
  }
  mesh.coord(4,0) = .54; mesh.coord(4,1) = .47;
  for(int ientt = 0; ientt < 4; ientt++) for(int node = 0; node < 3; node++){
    const int ipoin = nodes[ientt][node];
    double *coefficient = run.bak.met.elementCoefficient(ientt,node);
    coefficient[0] = .1+3.*coordinates[ipoin][0];
    coefficient[1] = 0.; coefficient[2] = -.1+2.*coordinates[ipoin][1];
  }
  const auto rejected = mesh.captureDGMetAtVertex(4);
  BOOST_REQUIRE_EQUAL((smooballdiff<MetricFieldFE,2,2>(mesh,4,region,&initial,&maximum0,&final,&maximum1,QuaFun::Connection)),2);
  BOOST_CHECK_EQUAL(mesh.coord(4,0),.54); BOOST_CHECK_EQUAL(mesh.coord(4,1),.47);
  const auto restored = mesh.captureDGMetAtVertex(4);
  for(size_t i = 0; i < rejected.size(); i++){
    BOOST_CHECK(rejected[i].met == restored[i].met);
    BOOST_CHECK_EQUAL(rejected[i].ienttBack,restored[i].ienttBack);
  }
  BOOST_CHECK_EQUAL(initial,final);
  for(double invalid : {std::numeric_limits<double>::quiet_NaN(),2000.}){
    for(int ientt = 0; ientt < 4; ientt++) for(int node = 0; node < 3; node++)
      run.bak.met.elementCoefficient(ientt,node)[0] = invalid;
    const auto saved = mesh.captureDGMetAtVertex(4);
    mesh.poi2bak[4] = 17;
    BOOST_CHECK_THROW((smooballdiff<MetricFieldFE,2,2>(mesh,4,region,&initial,&maximum0,&final,&maximum1,QuaFun::Connection)),MetrisExcept);
    BOOST_CHECK_EQUAL(mesh.coord(4,0),.54); BOOST_CHECK_EQUAL(mesh.coord(4,1),.47);
    BOOST_CHECK_EQUAL(mesh.poi2bak[4],17);
    const auto rolledBack = mesh.captureDGMetAtVertex(4);
    for(size_t i = 0; i < saved.size(); i++){
      BOOST_CHECK(saved[i].met == rolledBack[i].met);
      BOOST_CHECK_EQUAL(saved[i].ienttBack,rolledBack[i].ienttBack);
    }
    double restoredObjective = 0.;
    for(int ientt = 0; ientt < 4; ientt++)
      restoredObjective += metqua_connection<MetricFieldFE,2>(mesh,AsDeg::Pk,AsDeg::Pk,ientt);
    BOOST_CHECK_SMALL(restoredObjective-initial,2e-12);
  }
  MetrisParameters cgparam; cgparam.iverb = 0; cgparam.usrTarDeg = 2;
  cgparam.setMeshIn(meshFile.string()); cgparam.objective_quadrature_order = 5;
  MetrisRunner cgRun(nullptr,nullptr,cgparam);
  auto &cg = static_cast<Mesh<MetricFieldFE>&>(*cgRun.msh_g);
  cg.setBasis(FEBasis::Lagrange); cg.met.setSpace(MetSpace::Log);
  cgRun.bak.met.setSpace(MetSpace::Log);
  for(int ipoin = 0; ipoin < cg.npoin; ipoin++) for(int k = 0; k < 3; k++){
    cg.met(ipoin,k) = 0.; cgRun.bak.met(ipoin,k) = k == 1 ? 0. : 2000.;
  }
  cg.coord(4,0) = .54; cg.coord(4,1) = .47; cg.poi2bak[4] = -1;
  BOOST_CHECK_THROW((smooballdiff<MetricFieldFE,2,2>(cg,4,region,&initial,&maximum0,&final,&maximum1,QuaFun::Connection)),MetrisExcept);
  BOOST_CHECK_EQUAL(cg.coord(4,0),.54); BOOST_CHECK_EQUAL(cg.coord(4,1),.47);
  BOOST_CHECK_EQUAL(cg.poi2bak[4],-1);
  for(int ipoin = 0; ipoin < cg.npoin; ipoin++) for(int k = 0; k < 3; k++) BOOST_CHECK_EQUAL(cg.met(ipoin,k),0.);
}

BOOST_AUTO_TEST_CASE(dg_edge_node_smoothing_keeps_coefficients_and_avoids_background_sampling){
  for(MetSpace space : {MetSpace::Log,MetSpace::Exp}){
    MetrisParameters param; param.iverb = 0; param.usrTarDeg = 2;
    param.setMeshIn(fixture+".mesh"); param.dgMetFileName = fixture+".dgmet";
    param.objective_quadrature_order = 5;
    MetrisRunner run(nullptr,nullptr,param);
    auto &mesh = static_cast<Mesh<MetricFieldFE>&>(*run.msh_g);
    mesh.setBasis(FEBasis::Lagrange); mesh.met.setSpace(space);
    mesh.coord(4,0) += .08; mesh.coord(4,1) += .06;
    std::array<double,18> coefficients;
    std::array<int,6> owners;
    for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 3; node++){
      owners[3*ientt+node] = mesh.met2bak[mesh.met.getElementCoefficientIndex(ientt,node)];
      for(int k = 0; k < 3; k++) coefficients[9*ientt+3*node+k] = mesh.met.elementCoefficient(ientt,node)[k];
    }
    for(int ipoin = 0; ipoin < mesh.npoin; ipoin++) for(int k = 0; k < 3; k++)
      mesh.met.rfld(ipoin,k) = std::numeric_limits<double>::quiet_NaN();
    mesh.poi2bak[4] = -1; mesh.bak = nullptr;
    const double bary[3] = {.2,.3,.5};
    double met0[3], deriv0[6], met1[3], deriv1[6];
    mesh.met.getMetOnEntity(0,DifVar::Phys,bary,met0,deriv0);
    const std::array<double,2> position0 = {{mesh.coord(4,0),mesh.coord(4,1)}};
    double before, after;
    BOOST_REQUIRE((smoothSingleQuadraticEdgeControlPoint<MetricFieldFE>(
        mesh,4,0,0,QuaFun::Connection,&before,&after)));
    BOOST_CHECK_LT(after,before);
    BOOST_CHECK_GT(std::hypot(mesh.coord(4,0)-position0[0],mesh.coord(4,1)-position0[1]),1e-8);
    mesh.met.getMetOnEntity(0,DifVar::Phys,bary,met1,deriv1);
    double derivativeChange = 0.;
    for(int k = 0; k < 3; k++) BOOST_CHECK_EQUAL(met1[k],met0[k]);
    for(int k = 0; k < 6; k++) derivativeChange += std::abs(deriv1[k]-deriv0[k]);
    BOOST_CHECK_GT(derivativeChange,1e-5);
    double actual = 0.;
    for(int ientt = 0; ientt < 2; ientt++){
      BOOST_CHECK((classify_element_validity<2,2>(mesh,ientt).accepted_conservatively()));
      actual += metqua_connection<MetricFieldFE,2>(mesh,AsDeg::Pk,AsDeg::Pk,ientt);
      for(int node = 0; node < 3; node++){
        BOOST_CHECK_EQUAL(mesh.met2bak[mesh.met.getElementCoefficientIndex(ientt,node)],owners[3*ientt+node]);
        for(int k = 0; k < 3; k++) BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(ientt,node)[k],coefficients[9*ientt+3*node+k]);
      }
    }
    BOOST_CHECK_SMALL(actual-after,2e-12);
    BOOST_CHECK_EQUAL(mesh.poi2bak[4],-1);
    for(int ipoin = 0; ipoin < mesh.npoin; ipoin++){
      if(ipoin == 4) continue;
      for(int dim = 0; dim < 2; dim++) BOOST_CHECK_EQUAL(mesh.coord(ipoin,dim),points[ipoin][dim]);
    }
  }
}

BOOST_AUTO_TEST_CASE(continuous_sampling_matches_legacy_interpolation_without_committing){
  for(MetSpace space : {MetSpace::Log,MetSpace::Exp}){
    MetrisParameters param; param.iverb = 0; param.usrTarDeg = 2;
    param.setMeshIn(fixture+".mesh");
    MetrisRunner run(nullptr,nullptr,param);
    auto &mesh = static_cast<Mesh<MetricFieldFE>&>(*run.msh_g);
    run.bak.setBasis(FEBasis::Lagrange);
    run.bak.met.setBasis(FEBasis::Lagrange); run.bak.met.setSpace(MetSpace::Log);
    for(int ipoin = 0; ipoin < run.bak.npoin; ipoin++){
      run.bak.met(ipoin,0) = .2+run.bak.coord(ipoin,0);
      run.bak.met(ipoin,1) = .1*run.bak.coord(ipoin,1);
      run.bak.met(ipoin,2) = -.3+run.bak.coord(ipoin,1);
    }
    mesh.met.setSpace(space);
    mesh.coord(4,0) = .2; mesh.coord(4,1) = .3;
    const int oldSeed = mesh.poi2bak[4];
    const std::array<double,3> oldMet = {{mesh.met(4,0),mesh.met(4,1),mesh.met(4,2)}};
    BackgroundMetricSample sample;
    BOOST_REQUIRE_EQUAL(mesh.sampleMetBack(4,sample),0);
    BOOST_CHECK_EQUAL(mesh.poi2bak[4],oldSeed);
    for(int k = 0; k < 3; k++) BOOST_CHECK_EQUAL(mesh.met(4,k),oldMet[k]);
    double expected[3] = {.4,.03,0.};
    if(space == MetSpace::Exp) getexpmet_inp<2,double>(expected);
    for(int k = 0; k < 3; k++) BOOST_CHECK_SMALL(sample.met[k]-expected[k],2e-13);
    BOOST_REQUIRE_EQUAL(mesh.interpMetBack(4),0);
    for(int k = 0; k < 3; k++) BOOST_CHECK_EQUAL(mesh.met(4,k),sample.met[k]);
    BOOST_CHECK_EQUAL(mesh.poi2bak[4],sample.location.ientt);
    const auto accepted = sample;
    const int acceptedSeed = mesh.poi2bak[4];
    mesh.coord(4,0) = mesh.coord(4,1) = 20.;
    BOOST_CHECK_NE(mesh.sampleMetBack(4,sample),0);
    BOOST_CHECK(sample.met == accepted.met);
    BOOST_CHECK_NE(mesh.interpMetBack(4),0);
    BOOST_CHECK_EQUAL(mesh.poi2bak[4],acceptedSeed);
    for(int k = 0; k < 3; k++) BOOST_CHECK_EQUAL(mesh.met(4,k),accepted.met[k]);
  }
}

BOOST_AUTO_TEST_CASE(dgmet_command_line_loads_the_paired_input){
  std::vector<std::string> arguments = {"dg-io-test","--in",fixture+".mesh","--dgmet",fixture+".dgmet",
                                        "--tardeg","2","--verb","0"};
  std::vector<char*> argv;
  for(auto &argument : arguments) argv.push_back(argument.data());
  MetrisRunner run(int(argv.size()),argv.data());
  BOOST_CHECK(run.bak.met.getLayout() == MetricFieldLayout::Discontinuous);
  BOOST_CHECK_EQUAL(run.param->dgMetFileName,fixture+".dgmet");
}

BOOST_AUTO_TEST_CASE(checkpoint_reloads_changed_curved_geometry_and_local_coefficients){
  Temporary files;
  for(FEBasis geometryBasis : {FEBasis::Lagrange,FEBasis::Bezier})
    for(MetSpace space : {MetSpace::Log,MetSpace::Exp}){
      MetrisParameters param; param.iverb = 0; param.usrTarDeg = 2;
      param.setMeshIn(fixture+".mesh"); param.dgMetFileName = fixture+".dgmet";
      param.objective_quadrature_order = 5; param.outbasis = FEBasis::Bezier;
      MetrisRunner run(nullptr,nullptr,param);
      auto &mesh = static_cast<Mesh<MetricFieldFE>&>(*run.msh_g);
      mesh.setBasis(FEBasis::Lagrange);
      // A changed state with curvature and an independent local coefficient.
      mesh.coord(4,0) += .03; mesh.coord(4,1) -= .02;
      mesh.met.elementCoefficient(0,1)[0] += .17;
      mesh.setBasis(geometryBasis);
      mesh.met.setBasis(FEBasis::Bezier); mesh.met.setSpace(space);
      const std::string name = (geometryBasis == FEBasis::Lagrange ? "lag-" : "bez-")
                            + std::string(space == MetSpace::Log ? "log" : "spd");
      BOOST_CHECK_GT(checkCheckpointReload(run,files.path/name),0.);
      BOOST_CHECK_EQUAL(run.bak.met.elementCoefficient(0,1)[0],values[0][1][0]);
      for(int ipoin = 0; ipoin < run.bak.npoin; ipoin++) for(int dim = 0; dim < 2; dim++)
        BOOST_CHECK_EQUAL(run.bak.coord(ipoin,dim),points[ipoin][dim]);
    }
}

BOOST_AUTO_TEST_CASE(import_sans_target_values_and_derivatives){
  const char *directory = std::getenv("SANS_DG_TARGET_DIR");
  if(!directory){ BOOST_TEST_MESSAGE("SANS_DG_TARGET_DIR selects the exported SANS bundle"); return; }
  for(const auto &stem : {"target","target-spd"}){
    MetrisParameters param; param.iverb = 0; param.usrTarDeg = 2;
    const auto root = std::filesystem::path(directory);
    param.setMeshIn((root/(std::string(stem)+".meshb")).string());
    param.dgMetFileName = (root/(std::string(stem)+".dgmet")).string();
    MetrisRunner run(nullptr,nullptr,param);
    auto &mesh = static_cast<Mesh<MetricFieldFE>&>(*run.msh_g);
    BOOST_REQUIRE_EQUAL(mesh.met.getCoefficientCount(),54);
    std::ifstream coefficients(root/"target-coefficients.csv"); BOOST_REQUIRE(coefficients.good());
    std::string line; BOOST_REQUIRE(std::getline(coefficients,line));
    int coefficientCount = 0;
    while(std::getline(coefficients,line)){
      std::replace(line.begin(),line.end(),',',' ');
      std::istringstream record(line);
      int ientt,node,id; double x,y,met[3],logmet[3];
      BOOST_REQUIRE(record >> ientt >> node >> id >> x >> y);
      for(double &value : met) BOOST_REQUIRE(record >> value);
      for(double &value : logmet) BOOST_REQUIRE(record >> value);
      for(int k = 0; k < 3; k++){
        const double value = mesh.met.elementCoefficient(ientt,node)[k];
        if(std::string(stem) == "target") BOOST_CHECK_EQUAL(value,logmet[k]);
        else BOOST_CHECK_SMALL(value-logmet[k],2e-13);
        BOOST_CHECK_EQUAL(value,run.bak.met.elementCoefficient(ientt,node)[k]);
      }
      coefficientCount++;
    }
    BOOST_CHECK_EQUAL(coefficientCount,54);
    std::ifstream samples(root/"target-samples.csv"); BOOST_REQUIRE(samples.good());
    BOOST_REQUIRE(std::getline(samples,line));
    int count = 0;
    while(std::getline(samples,line)){
      std::replace(line.begin(),line.end(),',',' ');
      std::istringstream input(line);
      int id; double bary[3], point[2], expected[3], expectedD[6];
      BOOST_REQUIRE(input >> id >> bary[0] >> bary[1] >> bary[2] >> point[0] >> point[1]);
      for(double &value : expected) BOOST_REQUIRE(input >> value);
      for(double &value : expectedD) BOOST_REQUIRE(input >> value);
      double met[3], deriv[6], back[3], backD[6];
      double currentPoint[2];
      eval2<2,2>(mesh.coord,mesh.fac2poi[id-1],mesh.getBasis(),DifVar::None,DifVar::None,
                 bary,currentPoint,nullptr,nullptr);
      for(int dim = 0; dim < 2; dim++) BOOST_CHECK_SMALL(currentPoint[dim]-point[dim],2e-14);
      mesh.met.getMetOnEntity(id-1,DifVar::Phys,bary,met,deriv);
      run.bak.met.getMetOnEntity(id-1,DifVar::Phys,bary,back,backD);
      for(int k = 0; k < 3; k++){
        BOOST_CHECK_SMALL(met[k]-expected[k],3e-11*(1.+std::abs(expected[k])));
        BOOST_CHECK_EQUAL(met[k],back[k]);
      }
      for(int k = 0; k < 6; k++){
        BOOST_CHECK_SMALL(deriv[k]-expectedD[k],3e-11*(1.+std::abs(expectedD[k])));
        BOOST_CHECK_EQUAL(deriv[k],backD[k]);
      }
      count++;
    }
    BOOST_CHECK_EQUAL(count,108);
    if(const char *outputDirectory = std::getenv("METRIS_DG_CHECKPOINT_DIR")){
      std::filesystem::create_directories(outputDirectory);
      run.param->objective_quadrature_order = 5;
      const double value = checkCheckpointReload(run,std::filesystem::path(outputDirectory)/stem);
      std::cout << std::setprecision(17) << "DG_CHECKPOINT " << stem << " connection objective: " << value << '\n';
    }
  }
}
} // namespace Metris
