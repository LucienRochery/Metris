//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "MeshMetric.hxx"


#include "../linalg/det.hxx"
#include "../ho_constants.hxx"
#include "../aux_topo.hxx"


namespace Metris{


template<int gdim, int ideg,class MetricFieldType>
double getDomainVolume0(MeshMetric<MetricFieldType> &msh){
  static_assert(gdim == 2 || gdim == 3);

  int nentt       = gdim == 3 ? msh.nelem   : msh.nface;
  intAr2 &ent2poi = gdim == 3 ? msh.tet2poi : msh.fac2poi;
  constexpr auto ordent = ORDELT(gdim);

  MetSpace ispac0 = msh.met.getSpace();
  msh.met.setSpace(MetSpace::Log);

  // Nothing to do with the Jacobian, this is the highest supported degree.
  constexpr int idegq = METRIS_MAX_DEG_JACOBIAN;
  constexpr int nnodq = gdim == 3 ? tetnpps[idegq] : facnpps[idegq];

  constexpr auto eval = gdim == 2 ? eval2<gdim,ideg> : eval3<gdim,ideg>;

  MetSpace tarSpace = MetSpace::Log;
  constexpr int nnmet = (gdim*(gdim + 1))/2;
  double metl[nnmet];
  double dum[gdim];
  double jmat[gdim*gdim];
  double bary[gdim+1];

  // Integrate M^{1/2} over the domain, put into volM
  double volM = 0;
  const double dq = 1.0 / nnodq;
  for(int ientt = 0; ientt < nentt ; ientt++){
    if(isdeadent(ientt,ent2poi)) continue;

    for(int iquad = 0; iquad < nnodq; iquad++){
      for(int ii = 0; ii < gdim + 1; ii++){
        bary[ii] = ordent[idegq][iquad][ii] / (double) idegq;
      }
        
      // At quadrature point, compute metric 
      msh.met.getMetBary(AsDeg::Pk, DifVar::None, 
                     tarSpace, ent2poi[ientt], gdim, bary, metl, NULL);

      // and Jacobian 
      eval(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::Bary,DifVar::None,
           bary, dum, jmat, NULL);

      // Determinant of M^{1/2}
      double detJ = detmat<gdim>(jmat) / gdim / (gdim - 1); // Works in 2D and 3D
      double det12;
      if(tarSpace == MetSpace::Exp){
        det12 = detsym<gdim>(metl);
      }else{ // In log format, simply the trace
        det12 = metl[0] + metl[2];
        if constexpr(gdim == 3) det12 += metl[5];
        det12 = exp(det12);
      } 
      METRIS_ASSERT(det12 > 0.0); // Not a tolerance per se but suspecting sqrt should be valid for the same cases 
      // where this test is true (hoping)
      det12 = sqrt(det12);
      volM += det12 * detJ * dq ;
    }
  }
  msh.met.setSpace(ispac0);

  return volM;
}

// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double getDomainVolume0<2, n,  MetricFieldAnalytical >(MeshMetric<MetricFieldAnalytical> &msh);\
template double getDomainVolume0<3, n,  MetricFieldAnalytical >(MeshMetric<MetricFieldAnalytical> &msh);\
template double getDomainVolume0<2, n,  MetricFieldFE         >(MeshMetric<MetricFieldFE        > &msh);\
template double getDomainVolume0<3, n,  MetricFieldFE         >(MeshMetric<MetricFieldFE        > &msh);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()




template<class MetricFieldType>
void MeshMetric<MetricFieldType>::set_npoin(int npoin, bool skip_allocf){

  MeshBase::set_npoin(npoin, skip_allocf);
  int nnmet = (idim*(idim+1))/2;
  if(!skip_allocf) met.rfld.allocate(mpoin,nnmet);
  if(!skip_allocf) met.rfld.set_n(npoin);   
  
}


template<class MetricFieldType>
void MeshMetric<MetricFieldType>::set_nentt(int tdimn, int nentt, bool skip_allocf){
  switch(tdimn){
  case(-1):
    set_nbpoi(nentt);
    break;
  case(0):
    set_npoin(nentt,skip_allocf);
    break;
  case(1):
    set_nedge(nentt,skip_allocf);
    break;
  case(2):
    set_nface(nentt,skip_allocf);
    break;
  case(3):
    set_nelem(nentt,skip_allocf);
    break;
  }
}


template class MeshMetric<MetricFieldFE>;
template class MeshMetric<MetricFieldAnalytical>;


} // End Namespace
