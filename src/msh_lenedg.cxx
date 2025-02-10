//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "msh_lenedg.hxx"
#include "low_lenedg.hxx"

#include "aux_topo.hxx"
#include "aux_utils.hxx"
#include "CT_loop.hxx"
#include "Mesh/MeshMetric.hxx"

namespace Metris{

template<class MFT>
double getLengthEdges(MeshMetric<MFT> &msh, intAr2 &ilned, dblAr1 &rlned, LenTyp itype){
  double pct_unit = 0; 
  CT_FOR0_INC(2,3,gdim){if(msh.idim == gdim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
      pct_unit = getLengthEdges0<MFT,gdim,ideg>(msh,ilned,rlned,itype);
    }}CT_FOR1(ideg);
  }}CT_FOR1(gdim);
  return pct_unit;
}
template double getLengthEdges<MetricFieldAnalytical>(MeshMetric<MetricFieldAnalytical> &msh, 
                                      intAr2 &ilned, dblAr1 &rlned,LenTyp itype);
template double getLengthEdges<MetricFieldFE        >(MeshMetric<MetricFieldFE        > &msh, 
                                      intAr2 &ilned, dblAr1 &rlned,LenTyp itype);

template<class MFT>
double getLengthEdges_Bdry(MeshMetric<MFT> &msh, intAr2 &ilned, dblAr1 &rlned, LenTyp itype){
  double pct_unit = 0; 
  CT_FOR0_INC(2,3,gdim){if(msh.idim == gdim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
      pct_unit = getLengthEdges_Bdry0<MFT,gdim,ideg>(msh,ilned,rlned,itype);
    }}CT_FOR1(ideg);
  }}CT_FOR1(gdim);
  return pct_unit;
}
template double getLengthEdges_Bdry<MetricFieldAnalytical>(MeshMetric<MetricFieldAnalytical> &msh, 
                                      intAr2 &ilned, dblAr1 &rlned,LenTyp itype);
template double getLengthEdges_Bdry<MetricFieldFE        >(MeshMetric<MetricFieldFE        > &msh, 
                                      intAr2 &ilned, dblAr1 &rlned,LenTyp itype);


template<class MFT, int gdim, int ideg>
double getLengthEdges0(MeshMetric<MFT> &msh, intAr2 &ilned, dblAr1 &rlned, LenTyp itype){
  int tdimn = msh.get_tdim();
  METRIS_ASSERT(tdimn > 1); // implement lnoed1

  //msh.met.setSpace(MetSpace::Log);
  MetSpace ispac0 = msh.met.getSpace();

  if(itype == LenTyp::GeoSiz){    
    msh.met.setSpace(MetSpace::Exp);
  }else{
    msh.met.setSpace(MetSpace::Log);
  }
  int ned_unit = 0;
  int ned_totl = 0;

  // Only using quadrature to stat mesh: ok to go overboard
  const int nquad = 100;

  int nentt = msh.nentt(tdimn);
  intAr2 &ent2poi = msh.ent2poi(tdimn);
  int nedgl = (tdimn * (tdimn + 1)) / 2;
  const intAr2 lnoed(nedgl,2,tdimn == 2 ? lnoed2[0] : lnoed3[0]);
  HshTabInt2 hshTab; 
  hshTab.reserve(2*nentt);

  ilned.allocate(nentt,2); 
  rlned.allocate(nentt);
  ilned.set_n(0);
  rlned.set_n(0);

  for(int ientt = 0; ientt < nentt; ientt++){
    if(isdeadent(ientt,ent2poi)) continue;

    for(int iedgl = 0; iedgl < nedgl; iedgl++){
      int ip1 = ent2poi(ientt,lnoed(iedgl,0));
      int ip2 = ent2poi(ientt,lnoed(iedgl,1));
      auto key = stup2(ip1,ip2);
      auto t = hshTab.find(key);
      if(t != hshTab.end()) continue;
      hshTab.insert({key,ientt}); // dummy value

      double len;
      if(itype == LenTyp::GeoSiz){
        len = getlenedg_geosz<MFT,gdim,ideg>(msh,ientt,tdimn,iedgl,ExpTyp::Decomp);
      }else if(itype == LenTyp::Quad){
        len = getlenedg_quad<MFT,gdim,ideg>(msh,ientt,tdimn,iedgl,nquad);
      }else{
        METRIS_THROW_MSG(TODOExcept(),"Size interp scheme not implemented");
      }
      int iedgg = ilned.get_n();
      ilned.inc_n();
      ilned(iedgg,0) = std::get<0>(key);
      ilned(iedgg,1) = std::get<1>(key);

      if(len >= 1.0/sqrt(2) && len <= sqrt(2)) ned_unit++;
      ned_totl++;
      rlned.stack(len);
    }
  }

  msh.met.setSpace(ispac0);
  double pct_unit = ned_unit / (double) ned_totl;
  return pct_unit;
}


// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double getLengthEdges0<MetricFieldAnalytical, 2, n >(\
        MeshMetric<MetricFieldAnalytical> &msh,intAr2 &ilned, dblAr1 &rlned,LenTyp itype);\
template double getLengthEdges0<MetricFieldFE        , 2, n >(\
        MeshMetric<MetricFieldFE        > &msh,intAr2 &ilned, dblAr1 &rlned,LenTyp itype);\
template double getLengthEdges0<MetricFieldAnalytical, 3, n >(\
        MeshMetric<MetricFieldAnalytical> &msh,intAr2 &ilned, dblAr1 &rlned,LenTyp itype);\
template double getLengthEdges0<MetricFieldFE        , 3, n >(\
        MeshMetric<MetricFieldFE        > &msh,intAr2 &ilned, dblAr1 &rlned,LenTyp itype);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()




template<class MFT, int gdim, int ideg>
double getLengthEdges_Bdry0(MeshMetric<MFT> &msh, intAr2 &ilned, dblAr1 &rlned, 
                            LenTyp itype){
  constexpr int tdimn = 1;
  MetSpace ispac0 = msh.met.getSpace();

  if(itype == LenTyp::GeoSiz){    
    msh.met.setSpace(MetSpace::Exp);
  }else{
    msh.met.setSpace(MetSpace::Log);
  }

  int ned_unit = 0;
  int ned_totl = 0;

  // Only using quadrature to stat mesh: ok to go overboard
  const int nquad = 100;

  intAr2 &ent2poi = msh.ent2poi(1);

  ilned.allocate(msh.nedge,2); 
  rlned.allocate(msh.nedge);
  ilned.set_n(0);
  rlned.set_n(0);

  for(int iedge = 0; iedge < msh.nedge; iedge++){
    if(isdeadent(iedge,ent2poi)) continue;
    int ip1 = msh.edg2poi(iedge,0);
    int ip2 = msh.edg2poi(iedge,1);

    double len;
    if(itype == LenTyp::GeoSiz){
      len = getlenedg_geosz<MFT,gdim,ideg>(msh,iedge,tdimn,0,ExpTyp::Decomp);
    }else if(itype == LenTyp::Quad){
      len = getlenedg_quad<MFT,gdim,ideg>(msh,iedge,tdimn,0,nquad);
    }else{
      METRIS_THROW_MSG(TODOExcept(),"Size interp scheme not implemented");
    }
    int iedgg = ilned.get_n();
    ilned.inc_n();
    ilned(iedgg,0) = ip1;
    ilned(iedgg,1) = ip2;


    if(len >= 1.0/sqrt(2) && len <= sqrt(2)) ned_unit++;
    ned_totl++;
    rlned.stack(len);
  }

  double pct_unit = ned_unit / (double) ned_totl;
  msh.met.setSpace(ispac0);
  return pct_unit;
}


// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double getLengthEdges_Bdry0<MetricFieldAnalytical, 2, n >(\
        MeshMetric<MetricFieldAnalytical> &msh,intAr2 &ilned, dblAr1 &rlned,LenTyp itype);\
template double getLengthEdges_Bdry0<MetricFieldFE        , 2, n >(\
        MeshMetric<MetricFieldFE        > &msh,intAr2 &ilned, dblAr1 &rlned,LenTyp itype);\
template double getLengthEdges_Bdry0<MetricFieldAnalytical, 3, n >(\
        MeshMetric<MetricFieldAnalytical> &msh,intAr2 &ilned, dblAr1 &rlned,LenTyp itype);\
template double getLengthEdges_Bdry0<MetricFieldFE        , 3, n >(\
        MeshMetric<MetricFieldFE        > &msh,intAr2 &ilned, dblAr1 &rlned,LenTyp itype);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



}// end namespace
