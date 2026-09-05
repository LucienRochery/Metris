// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2026, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#include <src/MetrisRunner/MetrisRunner.hxx>
#include <src/Mesh/Mesh.hxx>
#include <src/io_libmeshb.hxx>
#include <src/low_topo.hxx>
#include <src/quality/msh_metqua.hxx>
#include <src/smoothing/msh_smooball.hxx>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Metris;

namespace
{

struct CustomOptions
{
  std::filesystem::path elevatedOutput;
  std::filesystem::path historyOutput;
  int maximumSweeps = 20;
  QuaFun objective = QuaFun::SizeShape;
  std::string objectiveName = "SizeShape";
  std::vector<char*> metrisArgv;
};

struct ControlPointSite
{
  int point;
  int seedTriangle;
  int localEdge;
};

CustomOptions extractCustomOptions(int argc, char** argv)
{
  CustomOptions options;
  options.metrisArgv.push_back(argv[0]);
  for(int argument = 1; argument < argc; argument++)
  {
    const std::string name(argv[argument]);
    if(name == "--elevated-output" || name == "--history-output"
       || name == "--maximum-sweeps" || name == "--smoothing-objective")
    {
      if(argument + 1 >= argc)
        throw std::runtime_error(name + " requires one value");
      const std::string value(argv[++argument]);
      if(name == "--elevated-output") options.elevatedOutput = value;
      else if(name == "--history-output") options.historyOutput = value;
      else if(name == "--maximum-sweeps") options.maximumSweeps = std::stoi(value);
      else if(value == "sizeshape")
      {
        options.objective = QuaFun::SizeShape;
        options.objectiveName = "SizeShape";
      }
      else if(value == "step-distance-sv")
      {
        options.objective = QuaFun::StepDistance;
        options.objectiveName = "StepDistanceSV";
      }
      else
      {
        throw std::runtime_error(
            "--smoothing-objective must be sizeshape or step-distance-sv");
      }
      continue;
    }
    options.metrisArgv.push_back(argv[argument]);
  }
  if(options.elevatedOutput.empty())
    throw std::runtime_error("--elevated-output is required");
  if(options.historyOutput.empty())
    throw std::runtime_error("--history-output is required");
  if(options.maximumSweeps < 1)
    throw std::runtime_error("--maximum-sweeps must be positive");
  return options;
}

double objectiveEnergy(Mesh<MetricFieldAnalytical>& mesh, QuaFun objective)
{
  bool invalid = false;
  double minimum = 0.;
  double maximum = 0.;
  double average = 0.;
  double energy = 0.;
  if(objective == QuaFun::SizeShape)
  {
    energy = getmetquamesh<MetricFieldAnalytical,QuaFun::SizeShape>(
        mesh,2,AsDeg::Pk,AsDeg::Pk,
        &invalid,&minimum,&maximum,&average,nullptr);
  }
  else if(objective == QuaFun::StepDistance)
  {
    energy = getmetquamesh<MetricFieldAnalytical,QuaFun::StepDistance>(
        mesh,2,AsDeg::Pk,AsDeg::Pk,
        &invalid,&minimum,&maximum,&average,nullptr);
  }
  else
  {
    throw std::runtime_error("Unsupported Q2-node smoothing objective");
  }
  if(invalid || !std::isfinite(energy))
    throw std::runtime_error("Q2-node smoothing objective evaluation is invalid");
  return energy;
}

std::vector<ControlPointSite> collectQuadraticControlPoints(
    const Mesh<MetricFieldAnalytical>& mesh,
    int originalPointCount)
{
  std::map<int,ControlPointSite> uniqueSites;
  for(int triangle = 0; triangle < mesh.nface; triangle++)
  {
    if(isdeadent(triangle,mesh.fac2poi)) continue;
    for(int localEdge = 0; localEdge < 3; localEdge++)
    {
      const int point = mesh.fac2poi(triangle,3 + localEdge);
      if(point < originalPointCount)
        throw std::runtime_error("Elevated edge node aliases a P1 point");
      uniqueSites.emplace(
          point,ControlPointSite{point,triangle,localEdge});
    }
  }

  std::vector<ControlPointSite> sites;
  sites.reserve(uniqueSites.size());
  for(const auto& [point,site] : uniqueSites) sites.push_back(site);
  return sites;
}

} // namespace

int main(int argc, char** argv)
{
  try
  {
    if(argc == 1)
    {
      std::cerr
          << "Usage: smoothQ2NodesOnly --in input.meshb --cad model.egads "
             "--anamet ID --sclmet SCALE --tardeg 2 --out final.meshb "
             "--elevated-output elevated.meshb --history-output history.csv "
             "[--maximum-sweeps 20] "
             "[--smoothing-objective sizeshape|step-distance-sv] "
             "[Metris options]\n";
      return 2;
    }

    CustomOptions options = extractCustomOptions(argc,argv);
    for(const auto& output : {options.elevatedOutput,options.historyOutput})
    {
      if(std::filesystem::exists(output))
        throw std::runtime_error("Refusing to overwrite: " + output.string());
      if(!output.parent_path().empty())
        std::filesystem::create_directories(output.parent_path());
    }

    MetrisRunner run(
        static_cast<int>(options.metrisArgv.size()),
        options.metrisArgv.data());
    if(run.metricFE)
      throw std::runtime_error("This focused test requires an analytic metric");

    auto& mesh = static_cast<Mesh<MetricFieldAnalytical>&>(*run.msh_g);
    if(mesh.idim != 2 || mesh.get_tdim() != 2 || mesh.curdeg != 1)
      throw std::runtime_error("Input must be a full-dimensional planar P1 mesh");
    if(run.param->usrTarDeg != 2)
      throw std::runtime_error("This focused test requires --tardeg 2");
    if(options.objective == QuaFun::StepDistance
       && !run.param->step_distance_shape_volume)
      throw std::runtime_error(
          "step-distance-sv smoothing requires --step-distance-shape-volume");

    const int originalPointCount = mesh.npoin;
    const int originalTriangleCount = mesh.nface;
    run.degElevate();
    mesh.setBasis(FEBasis::Lagrange);
    if(mesh.curdeg != 2 || mesh.nface != originalTriangleCount)
      throw std::runtime_error("Degree elevation changed the mesh topology");

    const std::vector<ControlPointSite> sites
        = collectQuadraticControlPoints(mesh,originalPointCount);
    if(static_cast<int>(sites.size()) != mesh.npoin - originalPointCount)
      throw std::runtime_error("Not every elevated point is an edge control point");

    writeMesh(options.elevatedOutput.string(),mesh);

    std::ofstream history(options.historyOutput);
    if(!history.good())
      throw std::runtime_error(
          "Cannot open history: " + options.historyOutput.string());
    history << std::setprecision(17)
            << "sweep,accepted_nodes,objective_energy\n";

    double energy = objectiveEnergy(mesh,options.objective);
    history << 0 << ',' << 0 << ',' << energy << '\n';
    std::cout << std::setprecision(17)
              << "Q2_NODE_SMOOTHING {\"stage\":\"elevated\","
              << "\"objective\":\"" << options.objectiveName << "\","
              << "\"p1_points\":" << originalPointCount << ','
              << "\"q2_nodes\":" << sites.size() << ','
              << "\"triangles\":" << originalTriangleCount << ','
              << "\"objective_energy\":" << energy << "}\n";

    for(int sweep = 1; sweep <= options.maximumSweeps; sweep++)
    {
      const double energyBefore = energy;
      int acceptedNodes = 0;
      for(const ControlPointSite& site : sites)
      {
        double regionalBefore = 0.;
        double regionalAfter = 0.;
        if(smoothSingleQuadraticEdgeControlPoint(
               mesh,site.point,site.seedTriangle,site.localEdge,
               options.objective,&regionalBefore,&regionalAfter))
        {
          acceptedNodes++;
        }
      }

      energy = objectiveEnergy(mesh,options.objective);
      if(energy > energyBefore
         + 1.e-12*std::max(1.,std::abs(energyBefore)))
        throw std::runtime_error("A local sweep increased the global objective");

      history << sweep << ',' << acceptedNodes << ',' << energy << '\n';
      std::cout << "Q2_NODE_SMOOTHING {\"stage\":\"sweep\","
                << "\"sweep\":" << sweep << ','
                << "\"accepted_nodes\":" << acceptedNodes << ','
                << "\"objective\":\"" << options.objectiveName << "\","
                << "\"objective_energy\":" << energy << "}\n";

      const double reduction = energyBefore - energy;
      if(acceptedNodes == 0
         || reduction <= 1.e-12*std::max(1.,std::abs(energyBefore)))
        break;
    }

    run.writeOutputs();
  }
  catch(const std::exception& exception)
  {
    std::cerr << "smoothQ2NodesOnly: " << exception.what() << '\n';
    return 1;
  }
  return 0;
}
