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
#include <limits>
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
  bool untilStalled = false;
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
    if(name == "--until-stalled"){
      options.untilStalled = true;
      continue;
    }
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
      else if(value == "connection")
      {
        options.objective = QuaFun::Connection;
        options.objectiveName = "ConnectionNormSquared";
      }
      else if(value == "distortion")
      {
        options.objective = QuaFun::Distortion;
        options.objectiveName = "LegacyDistortion";
      }
      else
        throw std::runtime_error(
            "--smoothing-objective must be sizeshape, step-distance-sv, connection, or distortion");
      continue;
    }
    options.metrisArgv.push_back(argv[argument]);
  }
  if(options.elevatedOutput.empty())
    throw std::runtime_error("--elevated-output is required");
  if(options.historyOutput.empty())
    throw std::runtime_error("--history-output is required");
  if(options.maximumSweeps < 0)
    throw std::runtime_error("--maximum-sweeps must be nonnegative (0 = elevation only)");
  return options;
}

template<class MFT>
double objectiveEnergy(Mesh<MFT>& mesh, QuaFun objective)
{
  bool invalid = false;
  double minimum = 0.;
  double maximum = 0.;
  double average = 0.;
  double energy = 0.;
  if(objective == QuaFun::SizeShape)
    energy = getmetquamesh<MFT,QuaFun::SizeShape>(
          mesh,2,AsDeg::Pk,AsDeg::Pk,
          &invalid,&minimum,&maximum,&average,nullptr);
  else if(objective == QuaFun::StepDistance)
    energy = getmetquamesh<MFT,QuaFun::StepDistance>(
          mesh,2,AsDeg::Pk,AsDeg::Pk,
          &invalid,&minimum,&maximum,&average,nullptr);
  else if(objective == QuaFun::Distortion)
    energy = getmetquamesh<MFT,QuaFun::Distortion>(
          mesh,2,AsDeg::Pk,AsDeg::Pk,
          &invalid,&minimum,&maximum,&average,nullptr);
  else if(objective == QuaFun::Connection)
    for(int face=0;face<mesh.nface;face++){
      if(mesh.fac2poi(face,0)<0) continue;
      energy += metqua_connection<MFT,2>(
          mesh,AsDeg::Pk,AsDeg::Pk,face);
    }
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
          << "Usage: smoothQ2AllNodes --in input.meshb --cad model.egads "
             "[--anamet ID --sclmet SCALE | --back background.meshb --met metric.solb] "
             "--tardeg 2 --out final.meshb "
             "--elevated-output elevated.meshb --history-output history.csv "
             "[--maximum-sweeps 50] "
             "[--until-stalled] "
             "[--smoothing-objective sizeshape|step-distance-sv|connection|distortion] "
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
    auto smooth = [&](auto& mesh){
    if(mesh.idim != 2 || mesh.get_tdim() != 2 || (mesh.curdeg != 1 && mesh.curdeg != 2))
      throw std::runtime_error("Input must be a full-dimensional planar P1 or P2 mesh");
    if(run.param->usrTarDeg != 2)
      throw std::runtime_error("This focused test requires --tardeg 2");
    if(options.objective == QuaFun::StepDistance
       && !run.param->step_distance_shape_volume)
      throw std::runtime_error(
          "step-distance-sv smoothing requires --step-distance-shape-volume");
    // A single production call owns the active set across all its sweeps.
    // Repeated one-sweep calls would reset tags and reconsider every node.
    run.param->opt_smoo_niter = options.maximumSweeps;
    std::vector<char> vertices(mesh.npoin,0);
    for(int face=0;face<mesh.nface;face++){
      if(mesh.fac2poi(face,0)<0) continue;
      for(int i=0;i<3;i++) vertices[mesh.fac2poi(face,i)]=1;
    }
    const int originalPointCount = std::count(vertices.begin(),vertices.end(),1);
    const int originalTriangleCount = mesh.nface;
    if(mesh.curdeg == 1) run.degElevate();
    mesh.setBasis(FEBasis::Lagrange);
    if(mesh.curdeg != 2 || mesh.nface != originalTriangleCount)
      throw std::runtime_error("Degree elevation changed the mesh topology");

    writeMesh(options.elevatedOutput.string(),mesh);
    const int geometryPointCount=mesh.npoin;
    std::vector<int> connectivity;
    for(int face=0;face<mesh.nface;face++)
      for(int i=0;i<6;i++) connectivity.push_back(mesh.fac2poi(face,i));

    std::ofstream history(options.historyOutput);
    if(!history.good())
      throw std::runtime_error(
          "Cannot open history: " + options.historyOutput.string());
    history << std::setprecision(17)
            << "batch,operation_ratio,objective_energy,max_coordinate_change\n";

    double energy = objectiveEnergy(mesh,options.objective);
    history << 0 << ',' << 0 << ',' << energy << ',' << 0 << '\n';
    std::cout << std::setprecision(17)
              << "Q2_ALL_NODE_SMOOTHING {\"stage\":\"elevated\","
              << "\"objective\":\"" << options.objectiveName << "\","
              << "\"p1_points\":" << originalPointCount << ','
              << "\"q2_nodes\":" << mesh.npoin - originalPointCount << ','
              << "\"triangles\":" << originalTriangleCount << ','
              << "\"objective_energy\":" << energy << "}\n";

    for(int batch=1;options.maximumSweeps>0;batch++)
    {
      const double energyBefore = energy;
      std::vector<double> coordinatesBefore(2*mesh.npoin);
      double coordinateScale=1.;
      for(int point=0;point<mesh.npoin;point++) for(int k=0;k<2;k++){
        coordinatesBefore[2*point+k]=mesh.coord(point,k);
        coordinateScale=std::max(coordinateScale,std::abs(mesh.coord(point,k)));
      }
      const double operationRatio =
          smoothInterior_Ball(mesh,options.objective,0,1);
      energy = objectiveEnergy(mesh,options.objective);
      if(energy > energyBefore
         + 1.e-12*std::max(1.,std::abs(energyBefore)))
        throw std::runtime_error(
            "An all-node sweep increased the global objective");

      double coordinateChange=0.;
      for(int point=0;point<mesh.npoin;point++) for(int k=0;k<2;k++)
        coordinateChange=std::max(coordinateChange,
            std::abs(mesh.coord(point,k)-coordinatesBefore[2*point+k]));
      const double roundingTolerance=64.*std::numeric_limits<double>::epsilon();
      const bool stalled=coordinateChange<=roundingTolerance*coordinateScale
          && energyBefore-energy<=roundingTolerance*std::max(1.,std::abs(energyBefore));
      history << batch << ',' << operationRatio << ',' << energy << ',' << coordinateChange << std::endl;
      std::cout << "Q2_ALL_NODE_SMOOTHING {\"stage\":\"production_smoothing\","
                << "\"batch\":" << batch << ','
                << "\"maximum_sweeps\":" << options.maximumSweeps << ','
                << "\"operation_ratio\":" << operationRatio << ','
                << "\"objective\":\"" << options.objectiveName << "\","
                << "\"objective_energy\":" << energy << ','
                << "\"max_coordinate_change\":" << coordinateChange << ','
                << "\"stalled\":" << (stalled ? "true" : "false") << "}" << std::endl;
      // Save recoverable checkpoints, then restart the active set if needed.
      // Every new batch scans all eligible nodes. A budget alone is not a
      // stopping condition in until-stalled mode.
      run.writeOutputs();
      if(!options.untilStalled || stalled) break;
    }

    if(mesh.nface != originalTriangleCount || mesh.npoin != geometryPointCount)
      throw std::runtime_error("Smoothing changed mesh counts");
    for(int face=0;face<mesh.nface;face++)
      for(int i=0;i<6;i++)
        if(mesh.fac2poi(face,i)!=connectivity[6*face+i])
          throw std::runtime_error("Smoothing changed connectivity");

    run.writeOutputs();
    };
    if(run.metricFE)
      smooth(static_cast<Mesh<MetricFieldFE>&>(*run.msh_g));
    else
      smooth(static_cast<Mesh<MetricFieldAnalytical>&>(*run.msh_g));
  }
  catch(const std::exception& exception)
  {
    std::cerr << "smoothQ2AllNodes: " << exception.what() << '\n';
    return 1;
  }
  return 0;
}
