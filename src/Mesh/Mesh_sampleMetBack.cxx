// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2026, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#include "Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../Localization/msh_localization.hxx"
#include "../Localization/low_localization.hxx"
#include "../low_geo/lenedg.hxx"
#include "../utils/CT_loop.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace Metris{
namespace{
template<class MFT>
int analyticalSample(Mesh<MFT> &msh,int ipoin,BackgroundMetricSample &sample){
  if constexpr(std::is_same_v<MFT,MetricFieldAnalytical>){
    BackgroundMetricSample result;
    result.space = msh.met.getSpace(); result.location.tdim = msh.idim;
    METRIS_ENFORCE(result.space == MetSpace::Log || result.space == MetSpace::Exp);
    result.location.status = 0;
    std::copy(msh.coord[ipoin],msh.coord[ipoin]+msh.idim,result.location.projected.begin());
    msh.met.getMetPhys(DifVar::None,result.space,msh.coord[ipoin],result.met.data(),nullptr);
    for(int k = 0; k < msh.met.getnnmet(); k++) METRIS_ENFORCE(std::isfinite(result.met[k]));
    sample = result;
    return 0;
  }
  return -1;
}
}

template<class MFT>
int Mesh<MFT>::locateMetBack(int ipoin,int tdim,int ientt,int iref,const double *algnd,
                             BackgroundLocation &location){
  METRIS_ENFORCE(ipoin >= 0 && ipoin < this->npoin && bak != nullptr);
  METRIS_ENFORCE(this->idim == bak->idim && (this->idim == 2 || this->idim == 3));
  if(bak->met.getLayout() == MetricFieldLayout::Discontinuous) bak->met.validateDiscontinuous();
  METRIS_ENFORCE(tdim >= 1 && tdim <= bak->get_tdim());
  METRIS_ENFORCE(ientt >= 0 && ientt < bak->nentt(tdim) && bak->ent2poi(tdim)(ientt,0) >= 0);
  if(iref >= 0 && bak->ent2ref(tdim)[ientt] != iref) return 4;
  BackgroundLocation result;
  result.tdim = tdim; result.ientt = ientt;
  result.bary.fill(std::numeric_limits<double>::quiet_NaN());
  result.projected.fill(std::numeric_limits<double>::quiet_NaN());
  const double *uvsrf = nullptr;
  if(tdim < this->get_tdim()){
    METRIS_ENFORCE_MSG(algnd != nullptr && iref >= 0,"Boundary localization needs reference and alignment");
    int ibpoi = this->poi2bpo[ipoin];
    METRIS_ENFORCE(ibpoi >= 0);
    if(this->bpo2ibi(ibpoi,1) == 0) ibpoi = this->bpo2ibi(ibpoi,3);
    METRIS_ENFORCE(ibpoi >= 0);
    uvsrf = this->bpo2rbi[ibpoi];
  }
  // Preserve a supplied DG owner whenever its full-dimensional map contains the point.
  if(bak->met.getLayout() == MetricFieldLayout::Discontinuous && tdim == bak->get_tdim()){
    if(locMeshQuick<2>(*bak,this->coord[ipoin])){
      result.status = LOC_ERR_OUTBB; location = result; return result.status;
    }
    if(inveval<2,2>(*bak,ientt,this->coord[ipoin],result.projected.data(),result.bary.data(),1e-6) == 0){
      result.status = 0; location = result; return 0;
    }
    result.bary.fill(std::numeric_limits<double>::quiet_NaN());
    result.projected.fill(std::numeric_limits<double>::quiet_NaN());
  }
  result.status = LOC_ERR_FAILP1;
  CT_FOR0_INC(1,METRIS_MAX_DEG,bdeg){if(bak->curdeg == bdeg){
    if(this->idim == 2){
      if(tdim == 1) result.status = locMesh<2,1,bdeg>(*bak,&result.ientt,this->coord[ipoin],tdim,uvsrf,iref,algnd,result.projected.data(),result.bary.data());
      else result.status = locMesh<2,2,bdeg>(*bak,&result.ientt,this->coord[ipoin],tdim,uvsrf,iref,nullptr,result.projected.data(),result.bary.data());
    }else{
      if(tdim == 1) result.status = locMesh<3,1,bdeg>(*bak,&result.ientt,this->coord[ipoin],tdim,uvsrf,iref,algnd,result.projected.data(),result.bary.data());
      else if(tdim == 2) result.status = locMesh<3,2,bdeg>(*bak,&result.ientt,this->coord[ipoin],tdim,uvsrf,iref,algnd,result.projected.data(),result.bary.data());
      else result.status = locMesh<3,3,bdeg>(*bak,&result.ientt,this->coord[ipoin],tdim,uvsrf,iref,nullptr,result.projected.data(),result.bary.data());
    }
  }}CT_FOR1(bdeg);
  location = result;
  return result.status;
}

template<class MFT>
int Mesh<MFT>::evaluateMetBackLocation(int ipoin,const BackgroundLocation &location,
                                       BackgroundMetricSample &sample,int ienttBack){
  if(location.tdim < 1 || location.ientt < 0) return 10+std::max(1,location.status);
  for(int i = 0; i <= location.tdim; i++) if(!std::isfinite(location.bary[i])) return 10+std::max(1,location.status);
  for(int i = 0; i < this->idim; i++) if(!std::isfinite(location.projected[i])) return 10+std::max(1,location.status);
  if(location.status == LOC_ERR_OUTBB) return 10+location.status;
  BackgroundMetricSample result;
  result.location = location; result.space = this->met.getSpace();
  METRIS_ENFORCE(result.space == MetSpace::Log || result.space == MetSpace::Exp);
  result.ienttMet = location.ientt; result.tdimMet = location.tdim;
  result.baryMet = location.bary;
  if(bak->met.getLayout() == MetricFieldLayout::Discontinuous){
    bak->met.validateDiscontinuous();
    double bary[3] = {};
    if(location.tdim == 2){
      std::copy(location.bary.begin(),location.bary.begin()+3,bary);
    }else{
      METRIS_ENFORCE_MSG(ienttBack >= 0 && ienttBack < bak->nface && bak->fac2poi(ienttBack,0) >= 0,
                        "DG boundary sampling requires an explicit owning triangle");
      result.ienttMet = ienttBack;
      const int *triangle = bak->fac2poi[ienttBack];
      const int *edge = bak->edg2poi[location.ientt];
      for(int node = 0; node < 2; node++){
        const int local = int(std::find(triangle,triangle+3,edge[node])-triangle);
        if(local == 3) return 6;
        bary[local] = location.bary[node];
      }
      // Match the Q2 edge node as well as its endpoints.
      int oppositeNode = 0;
      while(oppositeNode < 3 && (triangle[oppositeNode] == edge[0] || triangle[oppositeNode] == edge[1])) oppositeNode++;
      if(oppositeNode == 3 || bak->fac2poi(ienttBack,3+oppositeNode) != edge[2]) return 6;
    }
    if(result.space == MetSpace::Log)
      bak->met.getLogMetOnEntity(result.ienttMet,DifVar::None,bary,result.met.data(),nullptr);
    else bak->met.getMetOnEntity(result.ienttMet,DifVar::None,bary,result.met.data(),nullptr);
    result.tdimMet = 2;
    std::copy(bary,bary+3,result.baryMet.begin());
  }else{
    bak->met.getMetBary(AsDeg::Pk,DifVar::None,result.space,bak->ent2poi(location.tdim)[location.ientt],
                        location.tdim,location.bary.data(),result.met.data(),nullptr);
  }
  for(int k = 0; k < this->met.getnnmet(); k++) METRIS_ENFORCE(std::isfinite(result.met[k]));
  // Retain the existing acceptance threshold for a projected localization.
  bool checkDistance = location.status != 0;
  #ifndef NDEBUG
  checkDistance = true;
  #endif
  if(checkDistance){
    double tang[3];
    for(int dim = 0; dim < this->idim; dim++) tang[dim] = this->coord(ipoin,dim)-location.projected[dim];
    double length;
    if(this->idim == 2) length = result.space == MetSpace::Log
        ? getlenedg_log<2>(tang,result.met.data(),100,1e-6) : getlenedg<2>(tang,result.met.data());
    else length = result.space == MetSpace::Log
        ? getlenedg_log<3>(tang,result.met.data(),100,1e-6) : getlenedg<3>(tang,result.met.data());
    if(!(length < .5)) return 11;
  }
  sample = result;
  return 0;
}

template<class MFT>
int Mesh<MFT>::sampleMetBackOnEntity(int ipoin,int ienttBack,BackgroundMetricSample &sample){
  METRIS_ENFORCE(ipoin >= 0 && ipoin < this->npoin && bak != nullptr);
  if constexpr(std::is_same_v<MFT,MetricFieldAnalytical>) return analyticalSample(*this,ipoin,sample);
  BackgroundLocation location;
  locateMetBack(ipoin,bak->get_tdim(),ienttBack,-1,nullptr,location);
  return evaluateMetBackLocation(ipoin,location,sample,ienttBack);
}

template<class MFT>
int Mesh<MFT>::sampleMetBack(int ipoin,BackgroundMetricSample &sample,int ienttBack){
  METRIS_ENFORCE(ipoin >= 0 && ipoin < this->npoin);
  if constexpr(std::is_same_v<MFT,MetricFieldAnalytical>) return analyticalSample(*this,ipoin,sample);
  METRIS_ENFORCE(bak != nullptr);
  METRIS_ENFORCE_MSG(bak->met.getLayout() == MetricFieldLayout::Continuous || ienttBack >= 0,
                    "DG sampling requires an entity-specific background owner");
  const int tdim = this->getpoitdim(ipoin);
  METRIS_ENFORCE_MSG(tdim > 0,"Sample corner points with an explicit background entity");
  int iseed = this->poi2ent(ipoin,0);
  if(iseed < -1) iseed = -iseed-2;
  METRIS_ENFORCE(iseed >= 0 && iseed < this->nentt(tdim));
  const int iref = this->ent2ref(tdim)[iseed];
  double algnd[3]; this->get_algnd(ipoin,algnd);
  if(this->poi2bak[ipoin] >= 0 || (ienttBack >= 0 && tdim == bak->get_tdim()))
    if(sampleMetBack00(ipoin,tdim,iref,ipoin,algnd,sample,ienttBack) == 0) return 0;
  return sampleMetBack(ipoin,tdim,iseed,iref,algnd,sample,ienttBack);
}

template<class MFT>
int Mesh<MFT>::sampleMetBack(int ipoin,int tdim,int iseed,int iref,const double *algnd,
                             BackgroundMetricSample &sample,int ienttBack){
  METRIS_ENFORCE(ipoin >= 0 && ipoin < this->npoin);
  if constexpr(std::is_same_v<MFT,MetricFieldAnalytical>) return analyticalSample(*this,ipoin,sample);
  METRIS_ENFORCE(bak != nullptr && tdim >= 1 && tdim <= this->get_tdim());
  METRIS_ENFORCE(iseed >= 0 && iseed < this->nentt(tdim) && this->ent2poi(tdim)(iseed,0) >= 0);
  METRIS_ENFORCE(this->getpoitdim(ipoin) == tdim);
  METRIS_ENFORCE_MSG(bak->met.getLayout() == MetricFieldLayout::Continuous || ienttBack >= 0,
                    "DG sampling requires an entity-specific background owner");
  return sampleMetBack0(ipoin,tdim,iseed,iref,algnd,sample,ienttBack);
}

template<class MFT>
int Mesh<MFT>::sampleMetBack0(int ipoin,int tdim,int iseed,int iref,const double *algnd,
                              BackgroundMetricSample &sample,int ienttBack){
  const int nnode = getnnode(tdim,this->curdeg);
  const auto &ent2poi = this->ent2poi(tdim);
  for(int attempt = 0; attempt < 2; attempt++){
    bool skipped = false;
    for(int includeLower = 0; includeLower < 2; includeLower++){
      if(includeLower && !skipped) break;
      for(int node = 0; node < nnode; node++){
        const int ipseed = ent2poi(iseed,node), pdim = this->getpoitdim(ipseed);
        if(pdim > tdim || (pdim < tdim && !includeLower)){
          if(pdim < tdim) skipped = true;
          continue;
        }
        if(sampleMetBack00(ipoin,tdim,iref,ipseed,algnd,sample,ienttBack) == 0) return 0;
      }
    }
  }
  return 1;
}


template<class MFT>
int Mesh<MFT>::sampleMetBack00(int ipoin,int tdim,int iref,int ipseed,const double *algnd,
                               BackgroundMetricSample &sample,int ienttBack){
  METRIS_ENFORCE(ipoin >= 0 && ipoin < this->npoin);
  if constexpr(std::is_same_v<MFT,MetricFieldAnalytical>) return analyticalSample(*this,ipoin,sample);
  METRIS_ENFORCE(bak != nullptr && ipseed >= 0 && ipseed < this->npoin);
  METRIS_ENFORCE(tdim >= 1 && tdim <= bak->get_tdim());
  GETVDEPTH(this->param);
  if(bak->met.getLayout() == MetricFieldLayout::Discontinuous){
    METRIS_ENFORCE_MSG(ienttBack >= 0 && ienttBack < bak->nface && bak->fac2poi(ienttBack,0) >= 0,
                      "DG sampling requires an entity-specific background owner");
    if(tdim == bak->get_tdim()){
      BackgroundLocation location;
      const int status = locateMetBack(ipoin,tdim,ienttBack,iref,algnd,location);
      if(location.ientt < 0) return status;
      return evaluateMetBackLocation(ipoin,location,sample,ienttBack);
    }
  }
  int ientt = poi2bak[ipseed];
  if(ientt < 0){
    CPRINTF2(" - skip seed {} < 0\n",ientt);
    return 2; // Happens when called from the cavity operator.
  }

  CPRINTF2(" - init ientt = {}\n",ientt);

  int pdim_seed;

  // A point can be self-seeded; in that case, the ientt is already correct dim
  if(ipseed == ipoin) goto ientt_initialized;

  pdim_seed = this->getpoitdim(ipseed);

  if(pdim_seed > tdim){
    CPRINTF2(" - seed has dim {} > {} = point dim, skip\n",pdim_seed,tdim);
    return 3;
  }

  CPRINTF1(" - seed point dim {} initial bak seed {}\n",pdim_seed,ientt);

  if(pdim_seed == tdim){
    // Easiest case, ientt is already the one, and we just need to check ref.
    if(iref >= 0 && this->bak->ent2ref(tdim)[ientt] != iref){
      CPRINTF1(" # provided seed has ref {} != {} \n",this->bak->ent2ref(tdim)[ientt],iref);
      return 4;
    }
    goto ientt_initialized;
  }

  // Otherwise, the seed is dim < tdim: it can have several entity references of dimension tdim.
  if(pdim_seed == 0){
    // This is the easiest case, because we get a link to a vertex (corner)
    // in the back mesh.
    int ibpob = this->bak->poi2ebp(ientt, tdim, -1, iref);
    if(ibpob < 0){
      CPRINTF1(" # could not find seed entity dim {} ref {} from point {}\n",
               tdim, iref, poi2bak[ipseed]);
      return 5;
    }
    ientt = this->bak->bpo2ibi(ibpob,2);
    METRIS_ASSERT(this->bak->ent2ref(tdim)[ientt] == iref);
    CPRINTF1(" - corner seed tdim {} entity {}\n",tdim, ientt);
    goto ientt_initialized;
  }
  // If the seed is ref unconstrained, we will simply use edg2fac and fac2tet.
  // This is either if no ref constraint, or if single ref of target dim.
  if(iref < 0 || (tdim == 1 && this->CAD.ncaded == 1)
              || (tdim == 2 && this->CAD.ncadfa == 1)
              || (tdim == 3 && this->ndomn      == 1)){
    if(pdim_seed == 1){
      ientt = this->bak->edg2fac[ientt];
      CPRINTF1(" - seed edge -> face = {}\n",ientt);
    }
    if(pdim_seed <= 2 && tdim == 3){
      ientt = this->bak->fac2tet(ientt,0);
      if(ientt < 0) ientt = this->bak->fac2tet(ientt,1);
      CPRINTF1(" - seed face -> tetra = {}\n",ientt);
    }
    METRIS_ASSERT(ientt >= 0);
    goto ientt_initialized;
  }

  // Last case is pdim_seed >= 1, iref >= 0. If tdim is boundary, we have
  // poi2ebp. We can use this on a vertex of the back seed.
  if(this->isboundary_tdim(tdim)){
    int nnode = getnnode(pdim_seed, this->bak->curdeg);
    for(int inode = 0; inode < nnode; inode++){
      int ipoib = this->bak->ent2poi(pdim_seed)(ientt, inode);
      int ibpob = this->bak->poi2ebp(ipoib, tdim, -1, iref);
      if(ibpob < 0) continue;
      CPRINTF1(" - back seed entity {} vertex {} = {} points to back entity {} of ref {}\n",
               ientt, inode, ipoib, this->bak->bpo2ibi(ibpob, 2), iref);
      ientt = this->bak->bpo2ibi(ibpob, 2);
      METRIS_ASSERT(ientt >= 0 && ientt < this->bak->nentt(this->get_tdim()));
      goto ientt_initialized;
    }
    return 6;
  }else{
    static int nwarnprt = 0;
    if(nwarnprt++ < 10) printf("## IMPLEMENT BACK SEED INIT IN MULTI DOMAIN MESH");
    return 7;
  }


ientt_initialized:
  METRIS_ENFORCE(ientt >= 0 && ientt < bak->nentt(tdim));
  BackgroundLocation location;
  const int status = locateMetBack(ipoin,tdim,ientt,iref,algnd,location);
  if(location.ientt < 0) return status;
  return evaluateMetBackLocation(ipoin,location,sample,ienttBack);
}

#define INSTANTIATE(MFT) \
template int Mesh<MFT>::locateMetBack(int,int,int,int,const double*,BackgroundLocation&); \
template int Mesh<MFT>::evaluateMetBackLocation(int,const BackgroundLocation&,BackgroundMetricSample&,int); \
template int Mesh<MFT>::sampleMetBackOnEntity(int,int,BackgroundMetricSample&); \
template int Mesh<MFT>::sampleMetBack(int,BackgroundMetricSample&,int); \
template int Mesh<MFT>::sampleMetBack(int,int,int,int,const double*,BackgroundMetricSample&,int); \
template int Mesh<MFT>::sampleMetBack0(int,int,int,int,const double*,BackgroundMetricSample&,int); \
template int Mesh<MFT>::sampleMetBack00(int,int,int,int,const double*,BackgroundMetricSample&,int);
INSTANTIATE(MetricFieldFE)
INSTANTIATE(MetricFieldAnalytical)
#undef INSTANTIATE
} // namespace Metris
