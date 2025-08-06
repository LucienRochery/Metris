//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "misc.hxx"
#include "height.hxx"
#include "lenedg.hxx"

#include "../metris_constants.hxx"
#include "../Boundary/low_projsurf.hxx"
#include "../Mesh/Mesh.hxx"

namespace Metris{


template<class MFT, int tdim>
void getheightentP1_aniso(const Mesh<MFT> &msh, int ientt,
                          double *height){

  METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);
  constexpr int gdim = tdim;

  if constexpr(tdim == 2){

    for(int ied = 0; ied < tdim + 1; ied++){

      int ipoi1 = msh.fac2poi(ientt,lnoed2[ied][0]);
      int ipoi2 = msh.fac2poi(ientt,lnoed2[ied][1]);
      double tan[gdim];
      for(int ii = 0; ii < gdim; ii++) tan[ii] = msh.coord(ipoi2,ii)
                                               - msh.coord(ipoi1,ii);

      int ipoin = msh.fac2poi(ientt,ied);
      double x0 =  getprdl2<gdim>(msh.coord[ipoi1], tan);
      double x1 =  getprdl2<gdim>(msh.coord[ipoi2], tan);
      double tp = (getprdl2<gdim>(msh.coord[ipoin], tan) - x0) / (x1 - x0);

      //printf("Debug ied {}  x0 {} x1 {} xp {} tp {} ipoi1 {} ipoi2 {} ipoin {}\n",
      //  ied,x0,x1,getprdl2<gdim>(msh.coord[ipoin], tan),tp, ipoi1, ipoi2, ipoin);

      tp = MAX(0.0,MIN(1.0,tp));
      double dp[2];
      for(int ii = 0; ii < gdim; ii++) dp[ii] = (1.0 - tp) * msh.coord(ipoi1,ii)
                                              +        tp  * msh.coord(ipoi2,ii)
                                              -              msh.coord(ipoin,ii);
      double len = getlenedgsq<gdim>(dp, msh.met[ipoin]);
      height[ied] = sqrt(len);
    }
  }else{
    for(int ifa = 0; ifa < tdim + 1; ifa++){

      int ipoin = msh.tet2poi(ientt,ifa);
      int ipoi1 = msh.tet2poi(ientt,lnofa3[ifa][0]);
      int ipoi2 = msh.tet2poi(ientt,lnofa3[ifa][1]);
      int ipoi3 = msh.tet2poi(ientt,lnofa3[ifa][2]);

      double coopr[3];
      double bary[4];
      // ierro is whether point is inside face, but this doesn't matter. 
      //int ierro = 
      projptfacP1(msh.coord[ipoin], msh.coord[ipoi1], msh.coord[ipoi2], msh.coord[ipoi3],
                  bary, coopr);
      //METRIS_ENFORCE_MSG(ierro == 0,"projptfacP1 failed");


      double dp[3];
      for(int ii = 0; ii < gdim; ii++) dp[ii] = msh.coord(ipoin,ii) - coopr[ii];

      double len = getlenedgsq<gdim>(dp, msh.met[ipoin]);
      height[ifa] = sqrt(len);
    }
  }
}
template void getheightentP1_aniso<MetricFieldAnalytical, 2>(
         const Mesh<MetricFieldAnalytical> &msh, int ientt, double *height);
template void getheightentP1_aniso<MetricFieldFE        , 2>(
         const Mesh<MetricFieldFE        > &msh, int ientt, double *height);
template void getheightentP1_aniso<MetricFieldAnalytical, 3>(
         const Mesh<MetricFieldAnalytical> &msh, int ientt, double *height);
template void getheightentP1_aniso<MetricFieldFE        , 3>(
         const Mesh<MetricFieldFE        > &msh, int ientt, double *height);

}// namespace
