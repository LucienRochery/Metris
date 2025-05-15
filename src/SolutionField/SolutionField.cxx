//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "SolutionField.hxx"
#include "../low_eval.hxx"
#include "../utils/CT_loop.hxx"
#include "../Mesh/MeshBase.hxx"

namespace Metris{

void SolutionFieldAnalytical::setAnalyticalSolution(int ianasol_){
  int idim = msh->idim;
  METRIS_ASSERT(idim == 2 || idim == 3);
  
  if(ianasol_ <= 0 || ianasol_ > MAX_ANASOL_DEFINED(idim) )
    METRIS_THROW_MSG(WArgExcept(),"Invalid index: 1 - "<<MAX_ANASOL_DEFINED(idim )<<" accepted");
  
  this->ianasol = ianasol_;
  this->anasol  = (idim == 2 ? __ANASOL2D[this->ianasol-1] : __ANASOL3D[this->ianasol-1]);
}

void SolutionFieldAnalytical::setAnalyticalSolution(anasol_proto anasol_ptr){
  METRIS_ASSERT(anasol_ptr != NULL);
  
  this->ianasol = -1;
  this->anasol  = anasol_ptr;
}

SolutionFieldAnalytical &SolutionFieldAnalytical::operator=(const SolutionFieldAnalytical& inp){
  SolutionFieldBase::operator=(inp);
  this->ianasol = inp.ianasol;
  this->anasol  = inp.anasol;
  return *this;
}


double SolutionFieldAnalytical::getSolBary(int tdim, int ielem, 
                                           const double* __restrict__ bary, 
                                           std::initializer_list<double*> dfun) const{
  double coop[3];
  CT_FOR0_INC(2,3,gdim){if(gdim == msh->idim){
  CT_FOR0_INC(2,gdim,tdim_){if(tdim_ == tdim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh->curdeg){
      if constexpr(tdim_ == 1){
        eval1<gdim,ideg>(msh->coord, msh->edg2poi[ielem], msh->getBasis(), 
                         DifVar::None, DifVar::None, bary, coop, NULL, NULL);
      }else if(tdim == 2){
        eval2<gdim,ideg>(msh->coord, msh->fac2poi[ielem], msh->getBasis(), 
                         DifVar::None, DifVar::None, bary, coop, NULL, NULL);
      }else if(tdim_ == 3){
        eval3<gdim,ideg>(msh->coord, msh->tet2poi[ielem], msh->getBasis(), 
                         DifVar::None, DifVar::None, bary, coop, NULL, NULL);
      }
    }}CT_FOR1(ideg);
  }}CT_FOR1(tdim_);
  }}CT_FOR1(gdim);

  return anasol(NULL, coop, dfun);
}


} //namespace
