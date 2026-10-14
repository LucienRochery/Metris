// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2026, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#ifndef METRIS_NODE_SMOOTHING_STATE_HXX
#define METRIS_NODE_SMOOTHING_STATE_HXX

#include "../Mesh/Mesh.hxx"
#include <algorithm>
#include <array>
#include <type_traits>
#include <vector>

namespace Metris{

template<class MFT>
class NodeSmoothingState{
  struct LocalMet{
    double *met;
    int *ienttBack;
    std::array<double,3> value;
    int ientt0;
  };
  struct BoundaryValue{
    double *parameters;
    std::array<double,2> value;
  };
  double *coordinate;
  std::array<double,3> coordinate0{};
  int idim, nnmet;
  double *pointMet = nullptr;
  std::array<double,6> met0{};
  int *ienttBack = nullptr, ienttBack0 = -1;
  std::vector<LocalMet> localMet;
  std::vector<BoundaryValue> boundary;
  bool active = true;
public:
  NodeSmoothingState(Mesh<MFT> &msh,int ipoin,const intAr1 &region)
      : coordinate(nullptr), idim(msh.idim), nnmet(msh.met.getnnmet()){
    METRIS_ENFORCE(ipoin >= 0 && ipoin < msh.npoin && region.get_n() > 0);
    METRIS_ENFORCE(idim >= 2 && idim <= 3);
    coordinate = msh.coord[ipoin];
    std::copy(coordinate,coordinate+idim,coordinate0.begin());
    if(msh.met.getLayout() == MetricFieldLayout::Continuous){
      pointMet = msh.met[ipoin];
      std::copy(pointMet,pointMet+nnmet,met0.begin());
    }else{
      const int node = msh.getverent(region[0],msh.get_tdim(),ipoin);
      METRIS_ENFORCE(node >= 0);
      if(node < msh.get_tdim()+1){
        for(const auto &coefficient : msh.captureDGMetAtVertex(ipoin)){
          const int index = msh.met.getElementCoefficientIndex(coefficient.ientt,coefficient.node);
          localMet.push_back({msh.met.elementCoefficient(coefficient.ientt,coefficient.node),
                              &msh.met2bak[index],coefficient.met,coefficient.ienttBack});
        }
      }
    }
    if constexpr(std::is_same_v<MFT,MetricFieldFE>){
      METRIS_ENFORCE(msh.poi2bak.get_n() > ipoin);
      ienttBack = &msh.poi2bak[ipoin]; ienttBack0 = *ienttBack;
    }
    int record = msh.poi2bpo[ipoin], iteration = 0;
    while(record >= 0 && msh.bpo2ibi(record,0) == ipoin){
      METRIS_ENFORCE(iteration++ < METRIS_MAX_WHILE);
      boundary.push_back({msh.bpo2rbi[record],{{msh.bpo2rbi(record,0),msh.bpo2rbi(record,1)}}});
      record = msh.bpo2ibi(record,3);
    }
  }
  NodeSmoothingState(const NodeSmoothingState&) = delete;
  NodeSmoothingState &operator=(const NodeSmoothingState&) = delete;
  ~NodeSmoothingState(){ if(active) rollback(); }
  void commit() noexcept { active = false; }
  void rollback() noexcept {
    std::copy(coordinate0.begin(),coordinate0.begin()+idim,coordinate);
    if(pointMet) std::copy(met0.begin(),met0.begin()+nnmet,pointMet);
    if(ienttBack) *ienttBack = ienttBack0;
    for(const auto &local : localMet){
      std::copy(local.value.begin(),local.value.end(),local.met);
      *local.ienttBack = local.ientt0;
    }
    for(const auto &value : boundary)
      std::copy(value.value.begin(),value.value.end(),value.parameters);
    active = false;
  }
};

} // namespace Metris
#endif
