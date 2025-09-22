//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "msh_lenedg.hxx"
#include "low_geo/lenedg.hxx"

#include "aux_topo.hxx"
#include "utils/aux_misc.hxx"
#include "utils/CT_loop.hxx"
#include "utils/mprintf.hxx"
#include "Mesh/MeshMetric.hxx"
#include "BezierOffsets/low_gaps.hxx"

namespace Metris{

template<class MFT>
void getLengthEdges(MeshMetric<MFT> &msh, int tdim, int iref, 
                    intAr2 &ilned, dblAr1 &rlned, lenStat& stat, LenTyp itype){
  GETVDEPTH(msh.param);
  METRIS_ASSERT(tdim >= 1); // implement lnoed1

  stat.qua_short = 0;
  stat.qua_long = 0;

  //msh.met.setSpace(MetSpace::Log);
  //MetSpace ispac0 = msh.met.getSpace();

  //if(itype == LenTyp::GeoSiz){    
  //  msh.met.setSpace(MetSpace::Exp);
  //}else{
  //  msh.met.setSpace(MetSpace::Log);
  //}
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

  int ipdum = -1;
  if(itype == LenTyp::MetCrv){
    ipdum = msh.newpoitopo(PointType::Vertex,msh.get_tdim(), -1);
    if(msh.getBasis() == FEBasis::Undefined || msh.curdeg == 1){
      METRIS_ASSERT(msh.curdeg == 1);
      msh.forceBasisFlag(FEBasis::Bezier);
    }else{
      METRIS_THROW_MSG("TODO: LenType::MetCrv not implemented for Pk");
    }
  }

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
          }else if(itype == LenTyp::MetCrv){

            #if METRIS_MAX_DEG > 1
              METRIS_ASSERT(msh.get_tdim() == gdim);
              double offset[gdim];
              int ients = -1, iedgs = -1;
              if(tdim == gdim){
                ients = ientt;
                iedgs = iedgl;
              }else if(tdim == 1){
                ients = msh.edg2fac[ientt];
                if(gdim == 2){
                  iedgs = getedgfac(msh, ients, ip1, ip2);
                }else{
                  ients = msh.fac2tet(ients,0);
                  iedgs = getedgtet(msh, ients, ip1, ip2);
                }
              }else{
                ients = msh.fac2tet(ientt,0);
                iedgs = getedgtet(msh, ients, ip1, ip2);
              }
              METRIS_ASSERT(iedgs >= 0);
              METRIS_ASSERT(ients >= 0);

              getBezOffsetsEdge<MFT, gdim, ideg>(msh, gdim, msh.ent2poi(gdim)[ients], iedgs, offset);

              for(int ii = 0; ii < gdim; ii++){
                msh.coord(ipdum, ii) = (msh.coord(ip1,ii) + msh.coord(ip2,ii))/2
                                    + offset[ii];
              }

              double sz[2];
              int edg2pol[3] = {ip1, ip2, ipdum};
              len = getlenedg_geosz<MFT, gdim, 2>(msh, edg2pol, sz);
            #else
              METRIS_THROW_MSG("LenTyp::MetCrv not available for P1");
            #endif

          }else{
            METRIS_THROW_MSG("TODO: Size interp scheme not implemented");
          }
        }}CT_FOR1(ideg);
      }}CT_FOR1(gdim);
      if(std::isnan(len)){
        MPRINTF("## DEBUG NAN LEN EDGE !\n");
        MPRINTF("ientt = {} tdim = {} edge {} itype == GeoSiz? {}\n",
               ientt,tdim,iedgl,itype == LenTyp::GeoSiz);
      }
      if(len > 1000){
        fmt::print("## DEBUG WAIT HERE LEN = {} ip1 = {} ip2 = {} ientt = {}\n",len,ip1,ip2,ientt);
      }
      int iedgg = ilned.get_n();
      ilned.inc_n();
      ilned(iedgg,0) = std::get<0>(key);
      ilned(iedgg,1) = std::get<1>(key);

      if(len >= 1.0/sqrt(2) && len <= sqrt(2)) ned_unit++;
      ned_totl++;
      rlned.stack(len);

      if(len < 1.0){
        stat.qua_short = MAX(stat.qua_short, 1 - len);
      }else{
        stat.qua_long  = MAX(stat.qua_long, 1 - 1/len);
      }
    }
  }

  stat.prop_unit = ned_unit / (double) ned_totl;
  stat.qua_glo = MAX(stat.qua_short, stat.qua_long);
}

template void getLengthEdges<MetricFieldAnalytical>(MeshMetric<MetricFieldAnalytical> &msh, 
                                int tdim,int iref,intAr2 &ilned, dblAr1 &rlned,
                                lenStat& stat, LenTyp itype);
template void getLengthEdges<MetricFieldFE        >(MeshMetric<MetricFieldFE        > &msh, 
                                int tdim,int iref,intAr2 &ilned, dblAr1 &rlned,
                                lenStat& stat, LenTyp itype);


}// end namespace
