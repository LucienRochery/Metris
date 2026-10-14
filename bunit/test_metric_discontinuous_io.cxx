// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2026, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_metric_discontinuous_io
#include <boost/test/included/unit_test.hpp>
#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisRunner.hxx"
#include "quality/quafun_connection.hxx"
#include "quality/objective_quadrature_sample.hxx"
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
    const auto oldBoundary = boundary;
    BOOST_CHECK_NE(mesh.sampleMetBack(6,boundary,1),0);
    BOOST_CHECK(boundary.met == oldBoundary.met);
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
