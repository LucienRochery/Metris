//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "msh_metricCost.hxx"

namespace Metris{

template<int gdim,int tdim>
double getMetricCost(MeshMetric<MetricFieldAnalytical>& msh){

  // TODO: HO, only linear case for now
  const int ideg = msh.curdeg;
  METRIS_ENFORCE_MSG(ideg == 1, "Implement metric cost for HO");
  static_assert(tdim == gdim, "Implement metric cost for surface mesh in 3D");

  const intAr2& ent2poi = msh.ent2poi(tdim);

  // helper to compute metric determinant
  auto evalMetDet = [&](const double* met) -> double {

    double det = -1.;
    if constexpr (tdim == 2){
      det = met[0]*met[2] - met[1]*met[1];
    }else {

      const double m11 = met[0];
      const double m12 = met[1];
      const double m22 = met[2];
      const double m13 = met[3];
      const double m23 = met[4];
      const double m33 = met[5];

      det =    m11*(m22*m33 - m23*m23)
             - m12*(m12*m33 - m13*m23)
             + m13*(m12*m23 - m13*m22);
    }

    METRIS_ENFORCE_MSG(det > 0.,"Non-positive metric determinant!");
    return det;
  };

  double cost = 0.;

  for (int ientt = 0; ientt < msh.nentt(tdim); ientt++){

    if(isdeadent(ientt,ent2poi)) continue;

    double measEntt = 0.;
    isvalideltP1<gdim,tdim>(msh, ientt, NULL, &measEntt);

    double costEntt = 0.;

    int nsample = 0;
    // compute integrand at vertices of the element
    for (int iver = 0; iver < tdim + 1; iver++){

      const int ipoin = ent2poi(ientt,iver);
      costEntt += sqrt(evalMetDet(msh.met[ipoin]));
      nsample++;
    }

    // compute integrand at barycenter of element
    double coordBary[gdim];
    for(int icoord = 0; icoord < gdim; icoord++){
      coordBary[icoord] = 0.;
      for(int ibary = 0; ibary < tdim+1; ibary++){
        coordBary[icoord] += msh.coord(ent2poi(ientt,ibary),icoord)/((double)tdim+1.);
      }
    }

    constexpr int nnmet = (gdim*(gdim+1))/2;
    double metBary[nnmet];

    msh.met.getMetPhys(DifVar::None,msh.met.getSpace(),coordBary,metBary,NULL);

    costEntt += sqrt(evalMetDet(metBary));
    nsample++;

    cost += measEntt * costEntt/(double)nsample;
  }

  return cost;
}

template double getMetricCost<2,2>(MeshMetric<MetricFieldAnalytical>& msh);
template double getMetricCost<3,3>(MeshMetric<MetricFieldAnalytical>& msh);

template<int gdim,int tdim>
double getMetricCost(MeshMetric<MetricFieldFE>& msh){

  METRIS_THROW_MSG("Implement metric cost for MFT = MetricFieldFE");
}

template double getMetricCost<2,2>(MeshMetric<MetricFieldFE>& msh);
template double getMetricCost<3,3>(MeshMetric<MetricFieldFE>& msh);

} // end namespace
