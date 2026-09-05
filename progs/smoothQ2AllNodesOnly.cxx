// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2026, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#include <src/MetrisRunner/MetrisRunner.hxx>
#include <src/Mesh/Mesh.hxx>
#include <src/io_libmeshb.hxx>
#include <src/quality/msh_metqua.hxx>
#include <src/smoothing/msh_smooball.hxx>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
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
  int maximumSweeps = 50;
  QuaFun objective = QuaFun::SizeShape;
  std::string objectiveName = "SizeShape";
  std::vector<char*> metrisArgv;
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
      else if(name == "--maximum-sweeps")
        options.maximumSweeps = std::stoi(value);
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
        throw std::runtime_error(
            "--smoothing-objective must be sizeshape or step-distance-sv");
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
    energy = getmetquamesh<MetricFieldAnalytical,QuaFun::SizeShape>(
          mesh,2,AsDeg::Pk,AsDeg::Pk,
          &invalid,&minimum,&maximum,&average,nullptr);
  else if(objective == QuaFun::StepDistance)
    energy = getmetquamesh<MetricFieldAnalytical,QuaFun::StepDistance>(
          mesh,2,AsDeg::Pk,AsDeg::Pk,
          &invalid,&minimum,&maximum,&average,nullptr);
  else
    throw std::runtime_error("Unsupported all-node smoothing objective");
  if(invalid || !std::isfinite(energy))
    throw std::runtime_error("All-node smoothing objective is invalid");
  return energy;
}

} // namespace

int main(int argc, char** argv)
{
  try
  {
    if(argc == 1)
    {
      std::cerr
          << "Usage: smoothQ2AllNodesOnly --in input.meshb --cad model.egads "
             "--anamet ID --sclmet SCALE --tardeg 2 --out final.meshb "
             "--elevated-output elevated.meshb --history-output history.csv "
             "[--maximum-sweeps 50] "
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
    if(run.param->opt_smoo_niter != 1)
      throw std::runtime_error(
          "This focused test requires --opt-smoo-niter 1 so each recorded "
          "outer iteration is one all-node sweep");

    const int originalPointCount = mesh.npoin;
    const int originalTriangleCount = mesh.nface;
    run.degElevate();
    mesh.setBasis(FEBasis::Lagrange);
    if(mesh.curdeg != 2 || mesh.nface != originalTriangleCount)
      throw std::runtime_error("Degree elevation changed the mesh topology");

    writeMesh(options.elevatedOutput.string(),mesh);

    std::ofstream history(options.historyOutput);
    if(!history.good())
      throw std::runtime_error(
          "Cannot open history: " + options.historyOutput.string());
    history << std::setprecision(17)
            << "sweep,operation_ratio,objective_energy\n";

    double energy = objectiveEnergy(mesh,options.objective);
    history << 0 << ',' << 0 << ',' << energy << '\n';
    std::cout << std::setprecision(17)
              << "Q2_ALL_NODE_SMOOTHING {\"stage\":\"elevated\","
              << "\"objective\":\"" << options.objectiveName << "\","
              << "\"p1_points\":" << originalPointCount << ','
              << "\"q2_nodes\":" << mesh.npoin - originalPointCount << ','
              << "\"triangles\":" << originalTriangleCount << ','
              << "\"objective_energy\":" << energy << "}\n";

    for(int sweep = 1; sweep <= options.maximumSweeps; sweep++)
    {
      const double energyBefore = energy;
      const double operationRatio =
          smoothInterior_Ball(mesh,options.objective,0,1);
      energy = objectiveEnergy(mesh,options.objective);
      if(energy > energyBefore
         + 1.e-12*std::max(1.,std::abs(energyBefore)))
        throw std::runtime_error(
            "An all-node sweep increased the global objective");

      history << sweep << ',' << operationRatio << ',' << energy << '\n';
      std::cout << "Q2_ALL_NODE_SMOOTHING {\"stage\":\"sweep\","
                << "\"sweep\":" << sweep << ','
                << "\"operation_ratio\":" << operationRatio << ','
                << "\"objective\":\"" << options.objectiveName << "\","
                << "\"objective_energy\":" << energy << "}\n";

      const double reduction = energyBefore - energy;
      if(operationRatio <= 0.
         || reduction <= 1.e-12*std::max(1.,std::abs(energyBefore)))
        break;
    }

    run.writeOutputs();
  }
  catch(const std::exception& exception)
  {
    std::cerr << "smoothQ2AllNodesOnly: " << exception.what() << '\n';
    return 1;
  }
  return 0;
}
