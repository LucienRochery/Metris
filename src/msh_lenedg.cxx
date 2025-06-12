//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "msh_lenedg.hxx"
#include "low_lenedg.hxx"

#include "aux_topo.hxx"
#include "utils/aux_misc.hxx"
#include "utils/CT_loop.hxx"
#include "Mesh/MeshMetric.hxx"

namespace Metris{

template<class MFT>
double getLengthEdges(MeshMetric<MFT> &msh, int tdim, int iref, intAr2 &ilned, dblAr1 &rlned, LenTyp itype){
  
  METRIS_ASSERT(tdim >= 1); // implement lnoed1

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

  int nentt = msh.nentt(tdim);
  const intAr2 &ent2poi = msh.ent2poi(tdim);
  const intAr1 &ent2ref = msh.ent2ref(tdim);

  int nedgl = (tdim * (tdim + 1)) / 2;
  const intAr2 lnoed(nedgl,2,tdim == 1 ? lnoed1[0] :
                             tdim == 2 ? lnoed2[0] : lnoed3[0]);
  HshTab_I2I hshTab; 
  hshTab.reserve(2*nentt);

  ilned.allocate(nentt,2); 
  rlned.allocate(nentt);
  ilned.set_n(0);
  rlned.set_n(0);

  for(int ientt = 0; ientt < nentt; ientt++){
    if(isdeadent(ientt,ent2poi)) continue;
    if(iref >= 0 && iref != ent2ref[ientt]) continue;


    for(int iedgl = 0; iedgl < nedgl; iedgl++){

      int ip1 = ent2poi(ientt,lnoed(iedgl,0));
      int ip2 = ent2poi(ientt,lnoed(iedgl,1));
      auto key = stup2(ip1,ip2);
      if(tdim != 1){
        auto t = hshTab.find(key);
        if(t != hshTab.end()) continue;
        hshTab.insert({key,ientt}); // dummy value
      }

      double len;
      CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
        CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
          if(itype == LenTyp::GeoSiz){
            len = getlenedg_geosz<MFT,gdim,ideg>(msh,ientt,tdim,iedgl);
          }else if(itype == LenTyp::Quad){
            len = getlenedg_quad<MFT,gdim,ideg>(msh,ientt,tdim,iedgl,nquad);
          }else if(itype == LenTyp::LogIntrp){
            len = getlenedg_quad<MFT,gdim,ideg>(msh,ientt,tdim,iedgl,msh.nnode(1)-1);
          }else if(itype == LenTyp::BdryCor){
            METRIS_ENFORCE_MSG(tdim == 2 && gdim == 3, "BdryCor only available for triangles")
            static double sz[2];
            if constexpr(gdim == 3){
              len = getlenedg_geosz_plane<MFT,gdim,ideg>(msh,ientt,tdim,iedgl,sz);
            }
          }else{
            METRIS_THROW_MSG(TODOExcept(),"Size interp scheme not implemented");
          }
        }}CT_FOR1(ideg);
      }}CT_FOR1(gdim);
      if(std::isnan(len)){
        printf("## DEBUG NAN LEN EDGE !\n");
        printf("ientt = %d tdim = %d edge %d itype == GeoSiz? %d\n",
               ientt,tdim,iedgl,itype == LenTyp::GeoSiz);
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

template double getLengthEdges<MetricFieldAnalytical>(MeshMetric<MetricFieldAnalytical> &msh, 
                                int tdim,int iref,intAr2 &ilned, dblAr1 &rlned,LenTyp itype);
template double getLengthEdges<MetricFieldFE        >(MeshMetric<MetricFieldFE        > &msh, 
                                int tdim,int iref,intAr2 &ilned, dblAr1 &rlned,LenTyp itype);


}// end namespace
