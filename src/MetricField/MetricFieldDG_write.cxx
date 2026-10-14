// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2026, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#include "MetricField.hxx"
#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../io_libmeshb.hxx"
#include "../quality/quafun_connection.hxx"
#include "../quality/objective_quadrature_sample.hxx"
#include "nlohmann/json.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>

namespace Metris{
namespace{
std::filesystem::path outputPath(const std::string &name){
  return std::filesystem::absolute(name).lexically_normal();
}
}

void MetricFieldFE::writeDiscontinuousMetricFiles(const std::string &meshFilename,
                                                 const std::string &metFilename) const{
  METRIS_ENFORCE_MSG(layout == MetricFieldLayout::Discontinuous,"DG output requires a discontinuous metric");
  validateDiscontinuous();
  METRIS_ENFORCE_MSG(msh.param != nullptr,"DG output requires mesh parameters");
  METRIS_ENFORCE_MSG(outputPath(meshFilename) != outputPath(metFilename),"DG output filenames must differ");
  METRIS_ENFORCE_MSG(std::filesystem::path(meshFilename).extension() == ".mesh"
                    || std::filesystem::path(meshFilename).extension() == ".meshb",
                    "DG geometry output requires a .mesh or .meshb filename");
  for(int ientt = 0; ientt < msh.nface; ientt++)
    METRIS_ENFORCE_MSG(msh.fac2poi(ientt,0) >= 0,"DG output requires live entity rows");
  for(int ientt = 0; ientt < msh.nedge; ientt++)
    METRIS_ENFORCE_MSG(msh.edg2poi(ientt,0) >= 0,"DG output requires live boundary entity rows");

  // Convert a geometry copy so file coordinates and the manifest agree exactly.
  MetrisParameters parameters(*msh.param);
  parameters.outbasis = FEBasis::Lagrange;
  MeshBase geometry;
  geometry.param = &parameters;
  geometry = msh;
  geometry.setBasis(FEBasis::Lagrange);
  for(int ipoin = 0; ipoin < geometry.npoin; ipoin++) for(int dim = 0; dim < 2; dim++)
    METRIS_ENFORCE_MSG(std::isfinite(geometry.coord(ipoin,dim)),"Nonfinite DG output geometry");
  for(int ientt = 0; ientt < msh.nface; ientt++) for(int node = 0; node < 3; node++){
    std::array<double,3> logarithm;
    const double *met = elementCoefficient(ientt,node);
    for(int k = 0; k < 3; k++){
      METRIS_ENFORCE_MSG(std::isfinite(met[k]),"Nonfinite DG output coefficient");
      logarithm[k] = met[k];
    }
    if(ispace == MetSpace::Exp) getlogmet_inp<2,double>(logarithm.data());
    double expmet[3];
    getexpmet_frechet<2>(logarithm.data(),nullptr,0,expmet,nullptr);
  }

  std::ofstream output(metFilename);
  METRIS_ENFORCE_MSG(output.good(),"Cannot open DG output {}",metFilename);
  output.imbue(std::locale::classic());
  output << std::setprecision(std::numeric_limits<double>::max_digits10)
         << "MetrisDGMetric 1\nDimension 2\nTopology Triangle\nGeometryDegree 2\nGeometryBasis Lagrange\n"
         << "GeometryNodeOrder v0 v1 v2 e12 e20 e01\nMetricDegree 1\nMetricBasis "
         << (ibasis == FEBasis::Lagrange ? "Lagrange" : "Bernstein")
         << "\nMetricSpace " << (ispace == MetSpace::Log ? "Log" : "SPD")
         << "\nTensorComponents xx xy yy\nIndexBase 1\nPoints " << geometry.npoin << '\n';
  for(int ipoin = 0; ipoin < geometry.npoin; ipoin++)
    output << ipoin+1 << ' ' << geometry.coord(ipoin,0) << ' ' << geometry.coord(ipoin,1) << '\n';
  output << "Entities " << geometry.nface << '\n';
  for(int ientt = 0; ientt < geometry.nface; ientt++){
    output << ientt+1;
    for(int node = 0; node < 6; node++) output << ' ' << geometry.fac2poi(ientt,node)+1;
    output << ' ' << geometry.fac2ref[ientt]+1 << '\n';
  }
  output << "Coefficients " << 3*geometry.nface << '\n';
  for(int ientt = 0; ientt < geometry.nface; ientt++) for(int node = 0; node < 3; node++){
    const double *met = elementCoefficient(ientt,node);
    output << ientt+1 << ' ' << node+1 << ' ' << geometry.fac2poi(ientt,node)+1
           << ' ' << met[0] << ' ' << met[1] << ' ' << met[2] << '\n';
  }
  output << "End\n";
  output.close();
  METRIS_ENFORCE_MSG(output.good(),"Error writing DG output {}",metFilename);
  writeMesh(meshFilename,geometry,false);
}

double writeDGConnectionCheckpoint(Mesh<MetricFieldFE> &msh,
    const std::string &meshFilename, const std::string &metFilename,
    const std::string &objectiveFilename){
  METRIS_ENFORCE_MSG(msh.met.getLayout() == MetricFieldLayout::Discontinuous,
                    "DG checkpoint requires a discontinuous metric");
  METRIS_ENFORCE_MSG(msh.param != nullptr,"DG checkpoint requires mesh parameters");
  METRIS_ENFORCE_MSG(outputPath(objectiveFilename) != outputPath(meshFilename)
                    && outputPath(objectiveFilename) != outputPath(metFilename),
                    "DG checkpoint filenames must differ");
  nlohmann::json contributions = nlohmann::json::array();
  double total = 0.;
  for(int ientt = 0; ientt < msh.nface; ientt++){
    METRIS_ENFORCE_MSG(msh.fac2poi(ientt,0) >= 0,"DG checkpoint requires live entity rows");
    const double value = metqua_connection<MetricFieldFE,2>(msh,AsDeg::Pk,AsDeg::Pk,ientt);
    total += value;
    contributions.push_back({{"entity",ientt+1},{"value",value}});
  }
  METRIS_ENFORCE(std::isfinite(total));
  const auto root = outputPath(objectiveFilename).parent_path();
  const char *measure = standard_objective_quadrature_theta_mode() == ObjectiveQuadratureTheta::PhysicalMetricMeasure
      ? "physical_metric" : "physical";
  const nlohmann::json record = {
    {"version",1},{"objective_name","ConnectionNormSquared"},{"connection_norm","target_metric_squared"},
    {"dimension",2},{"geometry_degree",2},{"metric_degree",1},
    {"mesh_file",outputPath(meshFilename).lexically_relative(root).generic_string()},
    {"met_file",outputPath(metFilename).lexically_relative(root).generic_string()},
    {"quadrature_order",msh.param->objective_quadrature_order},{"integration_measure",measure},
    {"value",total},{"entities",contributions}
  };
  msh.met.writeDiscontinuousMetricFiles(meshFilename,metFilename);
  std::ofstream output(objectiveFilename);
  METRIS_ENFORCE_MSG(output.good(),"Cannot open DG objective output {}",objectiveFilename);
  output << record.dump(2) << '\n';
  output.close();
  METRIS_ENFORCE_MSG(output.good(),"Error writing DG objective output {}",objectiveFilename);
  return total;
}
} // namespace Metris
