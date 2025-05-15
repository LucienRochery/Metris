//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "Mesh.hxx"

#include "../MetrisRunner/MetrisParameters.hxx"
#include "../metris_constants.hxx"
#include "../ho_constants.hxx"
#include "../io_libmeshb.hxx"
#include "../aux_topo.hxx"
#include "../low_lenedg.hxx"
#include "../low_geo.hxx"
#include "../low_normal.hxx"
#include "../low_ccoef.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/aux_misc.hxx"
#include "../linalg/det.hxx"

#include "../Localization/msh_localization.hxx"
#include "../Localization/msh_interpFrontBack.hxx"


namespace Metris{

template<class MFT>
int Mesh<MFT>::newpoitopo(int tdimn, int ientt){
  int ipoin = MeshBase::newpoitopo(tdimn, ientt);

  // No need to add guesses here, that'll be done when calling interpMetBack
  poi2bak[ipoin] = -1;

  return ipoin;
}


template<>
void Mesh<MetricFieldAnalytical>::initialize(MetrisAPI *data, MeshBack &bak, 
  MetrisParameters &param){

  this->initializeCommon(data,bak,param);

  GETVDEPTH(this->param);

  if(param.ianamet >= 0){
    met.setAnalyticalMetric(param.ianamet);
  }else{
    METRIS_ASSERT(param.anamet_ptr != NULL);
    met.setAnalyticalMetric(param.anamet_ptr);
  }
  
  if(param.scaleMet){
    CPRINTF1("-- Front scaling metric by %15.7e\n", param.metScale);
    met.normalize(param.metScale);
  }

  met.forceBasisFlag(FEBasis::Lagrange);
  met.forceSpaceFlag(MetSpace::Exp);

  FEBasis mshbas0 = this->ibasis;
  this->setBasis(FEBasis::Lagrange);

  #if 0
  int tdim = this->get_tdim();
  const intAr2 &ent2poi = this->ent2poi(tdim);
  int nentt  = this->nentt(tdim);
  const int nnode = this->nnode(tdim);
  this->tag[0]++;
  double bary[4];
  const auto ordent = tdim == 1 ? ordedg.s[this->curdeg]
                    : tdim == 2 ? ordfac.s[this->curdeg] : ordtet.s[this->curdeg];

  for(int ientt = 0; ientt < nentt; ientt++){
    if(isdeadent(ientt,ent2poi)) continue;
    for(int iver = 0; iver < nnode; iver++){
      int ipoin = ent2poi(ientt,iver);
      if(this->poi2tag(0,ipoin) >= this->tag[0]) continue;
      this->poi2tag(0,ipoin) = this->tag[0];
      for(int ii = 0; ii < tdim+1; ii++)
        bary[ii] = ordent[iver][ii] / (double) this->curdeg;
      met.getMetBary(AsDeg::Pk, DifVar::None, MetSpace::Exp, ent2poi[ientt],
                     tdim, bary, met[ipoin], NULL);
    }
  }
  #endif
  for(int ipoin = 0; ipoin < npoin; ipoin++){
    met.getMetPhys(DifVar::None,MetSpace::Exp,
                   coord[ipoin],met[ipoin],NULL);
  }
  this->setBasis(mshbas0);

}



template<>
void Mesh<MetricFieldFE>::initialize(MetrisAPI *data, MeshBack &bak, 
  MetrisParameters &param){
  this->initializeCommon(data,bak,param);

  poi2bak.fill(-1);

  if(data == NULL && !param.inpBack){ // Copy case

    met = bak.met;

    for(int tdim = 1; tdim <= get_tdim(); tdim++){
      int nentt = this->nentt(tdim); 
      const intAr2 &ent2poi = this->ent2poi(tdim);
      int nnode = this->nnode(tdim);

      for(int ii = 0 ; ii < nentt; ii++){
        if(isdeadent(ii,ent2poi)) continue;
        for(int jj = 0; jj < nnode; jj++){
          int ip = ent2poi(ii,jj);
          int pdim = this->getpoitdim(ip);
          if(pdim == 0){
            // corner case, points to itself.
            poi2bak[ip] = ip;
            continue;
          }
          if(pdim != tdim) continue;
          poi2bak[ip] = ii;
        }
      }
    }


  }else{

    if(param.iverb >= 1) std::cout<<"-- Back metric interpolation\n";
    CT_FOR0_INC(1,METRIS_MAX_DEG,bdeg){if(bak.curdeg == bdeg){
      interpFrontBack<MetricFieldFE,bdeg>(*this,bak);
    }}CT_FOR1(bdeg);

  }


}


template<class MFT>
void Mesh<MFT>::initializeCommon(MetrisAPI *data, MeshBack &bak, 
                                 MetrisParameters &param){
  this->param = &param;
  this->bak   = &bak;

  // The back mesh has curdeg = strdeg regardless of target degree 
  this->strdeg = MAX(param.usrTarDeg,bak.curdeg);

  if(data == NULL && !param.inpBack){
    // In this case, simply copy from back mesh 
    MeshBase &dum = *this;
    dum = (MeshBase &) bak; // MeshBase::operator= 
  }else{
    MeshBase::initialize(data,param);
  }
  
  //this->cfa2tag.allocate(METRIS_MAXTAGS, this->ncadfa, true);
  //this->ced2tag.allocate(METRIS_MAXTAGS, this->ncaded, true);
  //this->cno2tag.allocate(METRIS_MAXTAGS, this->ncadno, true);

  //this->cfa2tag.fill(METRIS_MAXTAGS, this->ncadfa,0);
  //this->ced2tag.fill(METRIS_MAXTAGS, this->ncaded,0);
  //this->cno2tag.fill(METRIS_MAXTAGS, this->ncadno,0);

}



template<class MFT>
void Mesh<MFT>::set_npoin(int npoin, bool skipallocf){
  MeshMetric<MFT>::set_npoin(npoin, skipallocf);
  // Only allocate needed dimensions. Say dim = 2, allocate for edges and tri. 
  poi2bak.allocate(this->mpoin);
  poi2bak.set_n(this->npoin); 
}


template<class MFT>
void Mesh<MFT>::set_nentt(int tdim, int nentt, bool skipallocf){
  switch(tdim){
  case(-1):
    MeshBase::set_nbpoi(nentt);
    break;
  case(0):
    set_npoin(nentt,skipallocf);
    break;
  case(1):
    MeshBase::set_nedge(nentt,skipallocf);
    break;
  case(2):
    MeshBase::set_nface(nentt,skipallocf);
    break;
  case(3):
    MeshBase::set_nelem(nentt,skipallocf);
    break;
  }
}

template class Mesh<MetricFieldFE>;
template class Mesh<MetricFieldAnalytical>;


} // End namespace
