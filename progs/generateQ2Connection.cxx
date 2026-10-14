// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2026, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#include <src/MetrisRunner/MetrisRunner.hxx>
#include <src/Mesh/Mesh.hxx>
#include <src/quality/msh_metqua.hxx>
#include <src/smoothing/msh_smooball.hxx>
#include <src/low_topo.hxx>
#include <src/msh_checktopo.hxx>

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
std::vector<char*> argumentPointers(std::vector<std::string>& arguments)
{
  std::vector<char*> pointers;
  for(auto& argument : arguments) pointers.push_back(argument.data());
  return pointers;
}

double connectionEnergy(Mesh<MetricFieldFE>& mesh)
{
  double energy = 0;
  for(int face = 0; face < mesh.nface; face++){
    if(isdeadent(face,mesh.fac2poi)) continue;
    energy += metqua_connection<MetricFieldFE,2>(mesh,AsDeg::Pk,AsDeg::Pk,face);
  }
  return energy;
}
}

int main(int argc, char** argv)
{
  try{
    std::string input, metric, cad, outputDirectory;
    for(int argument = 1; argument < argc; argument++){
      const std::string option = argv[argument];
      if(argument + 1 == argc) throw std::runtime_error("Missing value for " + option);
      const std::string value = argv[argument + 1];
      argument++;
      if(option == "--in") input = value;
      else if(option == "--met") metric = value;
      else if(option == "--cad") cad = value;
      else if(option == "--out-dir") outputDirectory = value;
      else throw std::runtime_error("Unknown option " + option);
    }
    if(input.empty() || metric.empty() || cad.empty() || outputDirectory.empty())
      throw std::runtime_error("Usage: generateQ2Connection --in support.meshb --met metric.solb --cad square.egads --out-dir directory");

    input = std::filesystem::absolute(input).string();
    metric = std::filesystem::absolute(metric).string();
    cad = std::filesystem::absolute(cad).string();
    const auto directory = std::filesystem::absolute(outputDirectory);
    std::filesystem::create_directories(directory / "tmp");
    const std::string generated = (directory / "generated.meshb").string();
    const std::string remeshed = (directory / "remeshed.meshb").string();
    if(std::filesystem::exists(generated) || std::filesystem::exists(remeshed))
      throw std::runtime_error("Output directory already contains a mesh; use a fresh directory");

    const std::vector<std::string> common = {
      argv[0], "--back", input, "--met", metric, "--cad", cad, "--tardeg", "2",
      "--objective-quadrature-order", "5", "--objective-p", "1",
      "--step-distance-shape-volume", "--step-distance-regularization", "1e-8",
      "--step-distance-barrier-beta", "0", "--opt-pnorm", "1", "--opt-power", "1",
      "--opt-smoo-tol", ".005", "--jtol", "1e-6", "--vtol", "1e-9",
      "--prefix", (directory / "tmp").string() + "/", "--verb", "1", "--vdepth", "1"
    };
    auto generationArguments = common;
    const std::vector<std::string> generationOptions = {
      "--in", input, "--out", generated, "--adapt", "20", "--adp-unit-stop", "99.9",
      "--adp-stagn-stop", ".01", "--adp-opt-niter", "1", "--opt-niter", "5",
      "--opt-smoo-niter", "10", "--opt-smoo-vertices-only",
      "--opt-swap-pnorm", "1", "--opt-swap-niter", "100",
      "--qua-surf-wt-quality", "1", "--qua-surf-wt-normal", "1"
    };
    generationArguments.insert(generationArguments.end(),generationOptions.begin(),generationOptions.end());
    auto generationPointers = argumentPointers(generationArguments);
    {
      MetrisRunner generation(generationPointers.size(),generationPointers.data());
      if(!generation.metricFE || generation.msh_g->idim != 2 || generation.msh_g->curdeg != 2)
        throw std::runtime_error("Expected a Q2 planar mesh with an FE metric");
      generation.runMetris();
    }

    // Reload the generated front, retaining the original FE background for smoothing.
    auto smoothingArguments = common;
    const std::vector<std::string> smoothingOptions = {
      "--in", generated, "--out", remeshed, "--iflag1", "1"
    };
    smoothingArguments.insert(smoothingArguments.end(),smoothingOptions.begin(),smoothingOptions.end());
    auto smoothingPointers = argumentPointers(smoothingArguments);
    MetrisRunner smoothing(smoothingPointers.size(),smoothingPointers.data());
    auto& mesh = static_cast<Mesh<MetricFieldFE>&>(*smoothing.msh_g);
    mesh.setBasis(FEBasis::Lagrange);
    constexpr int maximumSweeps = 100000;
    smoothing.param->opt_smoo_niter = maximumSweeps;
    smoothing.param->opt_smoo_vertices_only = false;
    const double initialEnergy = connectionEnergy(mesh);
    const bool freezeVertices = false;
    smoothInterior_Ball(mesh,QuaFun::Connection,0,1,freezeVertices);
    check_topo(mesh,0);
    smoothing.writeOutputs();
    std::cout << std::setprecision(17)
              << "CONNECTION_ENERGY " << initialEnergy << " -> " << connectionEnergy(mesh) << '\n';
  }
  catch(const std::exception& exception){
    std::cerr << "generateQ2Connection: " << exception.what() << '\n';
    return 1;
  }
  return 0;
}
