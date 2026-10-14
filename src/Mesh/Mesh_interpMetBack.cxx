// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2026, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#include "Mesh.hxx"
#include <algorithm>
#include <type_traits>

namespace Metris{
namespace{
template<class MFT>
void commitSample(Mesh<MFT> &msh,int ipoin,const BackgroundMetricSample &sample){
  std::copy(sample.met.begin(),sample.met.begin()+msh.met.getnnmet(),msh.met[ipoin]);
  if constexpr(std::is_same_v<MFT,MetricFieldFE>) msh.poi2bak[ipoin] = sample.location.ientt;
}
}

template<class MFT>
int Mesh<MFT>::interpMetBack(int ipoin){
  this->met.requireContinuous("Pointwise background interpolation");
  BackgroundMetricSample sample;
  const int status = sampleMetBack(ipoin,sample);
  if(status == 0) commitSample(*this,ipoin,sample);
  return status;
}

template<class MFT>
int Mesh<MFT>::interpMetBack(int ipoin,int tdim,int iseed,int iref,const double *algnd){
  this->met.requireContinuous("Pointwise background interpolation");
  BackgroundMetricSample sample;
  const int status = sampleMetBack(ipoin,tdim,iseed,iref,algnd,sample);
  if(status == 0) commitSample(*this,ipoin,sample);
  return status;
}

template<class MFT>
int Mesh<MFT>::interpMetBack0(int ipoin,int tdim,int iseed,int iref,const double *algnd){
  this->met.requireContinuous("Pointwise background interpolation");
  BackgroundMetricSample sample;
  const int status = sampleMetBack0(ipoin,tdim,iseed,iref,algnd,sample,-1);
  if(status == 0) commitSample(*this,ipoin,sample);
  return status;
}

template<class MFT>
int Mesh<MFT>::interpMetBack00(int ipoin,int tdim,int iref,int ipseed,const double *algnd){
  this->met.requireContinuous("Pointwise background interpolation");
  BackgroundMetricSample sample;
  const int status = sampleMetBack00(ipoin,tdim,iref,ipseed,algnd,sample);
  if(status == 0) commitSample(*this,ipoin,sample);
  return status;
}

#define INSTANTIATE(MFT) \
template int Mesh<MFT>::interpMetBack(int); \
template int Mesh<MFT>::interpMetBack(int,int,int,int,const double*); \
template int Mesh<MFT>::interpMetBack0(int,int,int,int,const double*); \
template int Mesh<MFT>::interpMetBack00(int,int,int,int,const double*);
INSTANTIATE(MetricFieldFE)
INSTANTIATE(MetricFieldAnalytical)
#undef INSTANTIATE
} // namespace Metris
