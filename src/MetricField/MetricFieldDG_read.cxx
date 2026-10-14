// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2026, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#include "MetricField.hxx"
#include "../Mesh/MeshBase.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <locale>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <utility>

namespace Metris{
namespace{
class DGRecords{
  std::ifstream input;
  const std::string filename;
  int lineNumber = 0;
public:
  explicit DGRecords(const std::string &path) : input(path), filename(path){
    METRIS_ENFORCE_MSG(input.good(),"Cannot open DG metric file {}",path);
  }
  bool next(std::istringstream &record){
    std::string line;
    while(std::getline(input,line)){
      lineNumber++;
      const auto start = line.find_first_not_of(" \t\r");
      if(start == std::string::npos || line[start] == '#') continue;
      record.clear(); record.str(line);
      record.imbue(std::locale::classic());
      return true;
    }
    METRIS_ENFORCE_MSG(input.eof(),"Error reading DG metric file {}",filename);
    return false;
  }
  template<class... T> void read(T&... values){
    std::istringstream record;
    METRIS_ENFORCE_MSG(next(record),"Missing DG record in {}",filename);
    METRIS_ENFORCE_MSG(bool((record >> ... >> values)),
                      "Invalid DG record at {}:{}",filename,lineNumber);
    record >> std::ws;
    METRIS_ENFORCE_MSG(record.eof(),"Extra DG values at {}:{}",filename,lineNumber);
  }
  void expect(const std::string &expected){
    std::istringstream record;
    METRIS_ENFORCE_MSG(next(record),"Missing DG header {}",expected);
    std::string actual, token;
    while(record >> token){ if(!actual.empty()) actual += ' '; actual += token; }
    METRIS_ENFORCE_MSG(actual == expected,"Expected DG header '{}', got '{}'",expected,actual);
  }
  int count(const std::string &label){
    std::string actual; int value;
    read(actual,value);
    METRIS_ENFORCE_MSG(actual == label && value > 0,"Invalid DG count for {}",label);
    return value;
  }
  std::string choice(const std::string &label,const std::string &a,const std::string &b){
    std::string actual,value; read(actual,value);
    METRIS_ENFORCE_MSG(actual == label && (value == a || value == b),
                      "Unsupported DG {}: {}",label,value);
    return value;
  }
  void end(){
    expect("End");
    std::istringstream record;
    METRIS_ENFORCE_MSG(!next(record),"Trailing records in DG metric file {}",filename);
  }
};
using Vertices = std::array<int,3>;
Vertices key(Vertices vertices){ std::sort(vertices.begin(),vertices.end()); return vertices; }
}

void MetricFieldFE::readDiscontinuousMetricFile(const std::string &filename,
                                                const std::vector<int> &sourceToPoint){
  METRIS_ENFORCE_MSG(metricClass() == MetricClass::MetricFieldFE,
                    "DG input requires an FE metric");
  METRIS_ENFORCE_MSG(msh.idim == 2 && msh.get_tdim() == 2 && msh.curdeg == 2
                    && msh.getBasis() == FEBasis::Lagrange,
                    "DG input requires planar Q2 Lagrange geometry");
  DGRecords records(filename);
  records.expect("MetrisDGMetric 1"); records.expect("Dimension 2");
  records.expect("Topology Triangle"); records.expect("GeometryDegree 2");
  records.expect("GeometryBasis Lagrange");
  records.expect("GeometryNodeOrder v0 v1 v2 e12 e20 e01");
  records.expect("MetricDegree 1");
  const auto basis = records.choice("MetricBasis","Lagrange","Bernstein");
  const auto space = records.choice("MetricSpace","Log","SPD");
  records.expect("TensorComponents xx xy yy"); records.expect("IndexBase 1");
  const int npoin = records.count("Points");
  METRIS_ENFORCE_MSG(npoin == msh.npoin,"DG point count does not match geometry");
  METRIS_ENFORCE_MSG(sourceToPoint.empty() || int(sourceToPoint.size()) == npoin,
                    "DG point map has the wrong size");
  std::vector<int> pointMap(npoin);
  std::vector<bool> seen(npoin,false), mapped(npoin,false);
  for(int i = 0; i < npoin; i++){
    const int ipoin = sourceToPoint.empty() ? i : sourceToPoint[i];
    METRIS_ENFORCE_MSG(ipoin >= 0 && ipoin < npoin && !mapped[ipoin],
                      "DG point map must be a bijection");
    mapped[ipoin] = true; pointMap[i] = ipoin;
  }
  for(int i = 0; i < npoin; i++){
    int id; double x,y; records.read(id,x,y);
    METRIS_ENFORCE_MSG(id >= 1 && id <= npoin && !seen[id-1],"Invalid or duplicate DG point ID");
    seen[id-1] = true;
    const int ipoin = pointMap[id-1];
    METRIS_ENFORCE_MSG(std::isfinite(x) && std::isfinite(y)
                      && x == msh.coord(ipoin,0) && y == msh.coord(ipoin,1),
                      "DG geometry coordinates differ at source point {}",id);
  }
  std::map<Vertices,int> owners;
  for(int ientt = 0; ientt < msh.nface; ientt++){
    if(msh.fac2poi(ientt,0) < 0) continue;
    const Vertices vertices = {{msh.fac2poi(ientt,0),msh.fac2poi(ientt,1),msh.fac2poi(ientt,2)}};
    METRIS_ENFORCE_MSG(key(vertices)[0] != key(vertices)[1] && key(vertices)[1] != key(vertices)[2],
                      "Degenerate DG owner vertices");
    METRIS_ENFORCE_MSG(owners.emplace(key(vertices),ientt).second,"Ambiguous DG owner connectivity");
  }
  const int nentt = records.count("Entities");
  METRIS_ENFORCE_MSG(nentt == int(owners.size()),"DG entity count does not match geometry");
  struct Owner { int ientt = -1; std::array<int,3> local; Vertices vertices; };
  std::vector<Owner> association(nentt);
  std::set<int> usedOwners;
  for(int i = 0; i < nentt; i++){
    int id, ref; std::array<int,6> points;
    records.read(id,points[0],points[1],points[2],points[3],points[4],points[5],ref);
    METRIS_ENFORCE_MSG(id >= 1 && id <= nentt && association[id-1].ientt < 0,
                      "Invalid or duplicate DG entity ID");
    auto &owner = association[id-1];
    for(int node = 0; node < 6; node++){
      METRIS_ENFORCE_MSG(points[node] >= 1 && points[node] <= npoin,"Invalid DG connectivity point");
      points[node] = pointMap[points[node]-1];
    }
    owner.vertices = {{points[0],points[1],points[2]}};
    const auto found = owners.find(key(owner.vertices));
    METRIS_ENFORCE_MSG(found != owners.end(),"DG entity {} is missing from geometry",id);
    owner.ientt = found->second;
    METRIS_ENFORCE_MSG(usedOwners.insert(owner.ientt).second,"Duplicate DG owner entity");
    METRIS_ENFORCE_MSG(ref == msh.fac2ref[owner.ientt]+1,"DG entity reference differs");
    for(int node = 0; node < 3; node++){
      const int *current = msh.fac2poi[owner.ientt];
      owner.local[node] = int(std::find(current,current+3,points[node])-current);
      METRIS_ENFORCE(owner.local[node] < 3);
    }
    for(int edge = 0; edge < 3; edge++){
      // The edge opposite corner edge joins the other two local vertices.
      const int currentEdge = owner.local[edge];
      METRIS_ENFORCE_MSG(points[3+edge] == msh.fac2poi(owner.ientt,3+currentEdge),
                        "DG Q2 edge-node association differs for entity {}",id);
    }
  }
  METRIS_ENFORCE(nentt <= std::numeric_limits<int>::max()/3);
  METRIS_ENFORCE_MSG(records.count("Coefficients") == 3*nentt,"Wrong DG coefficient count");
  MetricFieldFE staged(msh);
  staged.initializeDiscontinuous();
  std::vector<bool> coefficients(3*nentt,false);
  for(int i = 0; i < 3*nentt; i++){
    int id,node,point; double met[3];
    records.read(id,node,point,met[0],met[1],met[2]);
    METRIS_ENFORCE_MSG(id >= 1 && id <= nentt && node >= 1 && node <= 3
                      && !coefficients[3*(id-1)+node-1],"Invalid or duplicate DG coefficient ID");
    coefficients[3*(id-1)+node-1] = true;
    const auto &owner = association[id-1];
    METRIS_ENFORCE_MSG(point >= 1 && point <= npoin && pointMap[point-1] == owner.vertices[node-1],
                      "DG coefficient vertex does not match its owning entity");
    for(double value : met) METRIS_ENFORCE_MSG(std::isfinite(value),"Nonfinite DG coefficient");
    if(space == "SPD") getlogmet_inp<2,double>(met);
    double expmet[3];
    getexpmet_frechet<2>(met,nullptr,0,expmet,nullptr);
    std::copy(met,met+3,staged.elementCoefficient(owner.ientt,owner.local[node-1]));
  }
  records.end();
  staged.setBasis(basis == "Lagrange" ? FEBasis::Lagrange : FEBasis::Bezier);
  elementToCoefficient = std::move(staged.elementToCoefficient);
  localCoefficients = std::move(staged.localCoefficients);
  layout = staged.layout; metricDegree = staged.metricDegree;
  ispace = staged.ispace; ibasis = staged.ibasis;
}
} // namespace Metris
