//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "msh_metricCost.hxx"

#include "quality/objective_quadrature_sample.hxx"
#include "quality/simplex_quadrature.hxx"

#include "utils/CT_loop.hxx"

namespace Metris{

namespace{

template<class MFT,int gdim,int tdim,int mshdeg>
double getMetricCost0(Mesh<MFT>& msh){
  static_assert(tdim == gdim, "Implement metric cost for surface mesh in 3D");

  const intAr2& ent2poi = msh.ent2poi(tdim);
  const SimplexQuadratureView<tdim> quadrature
      = get_objective_quadrature<tdim>(
            msh.param->objective_quadrature_order);

  double cost = 0.;

  for (int ientt = 0; ientt < msh.nentt(tdim); ientt++){

    if(isdeadent(ientt,ent2poi)) continue;

    for(int iquad = 0; iquad < quadrature.size(); iquad++){
      const auto sample
          = prepare_objective_quadrature_sample<MFT,gdim,tdim,mshdeg>(
                msh,AsDeg::P1,ent2poi[ientt],quadrature[iquad],
                ObjectiveQuadratureTheta::PhysicalMetricMeasure);
      METRIS_ENFORCE_MSG(sample.theta_is_valid,
                         "Invalid metric-cost quadrature sample");
      cost += sample.quadrature_weight*sample.theta;
    }
  }

  return cost;
}

} // namespace

template<class MFT,int gdim,int tdim>
double getMetricCost(MeshMetric<MFT>& mesh_metric){
  auto &msh = static_cast<Mesh<MFT>&>(mesh_metric);
  const MetSpace original_metric_space = msh.met.getSpace();
  msh.met.setSpace(MetSpace::Exp);

  double cost = 0.;
  bool degree_found = false;
  try{
    CT_FOR0_INC(1,METRIS_MAX_DEG,mshdeg){if(msh.curdeg == mshdeg){
      cost = getMetricCost0<MFT,gdim,tdim,mshdeg>(msh);
      degree_found = true;
    }}CT_FOR1(mshdeg);
    METRIS_ENFORCE_MSG(degree_found,
                       "Unsupported mesh degree {} in metric-cost integration",
                       msh.curdeg);
  }catch(...){
    msh.met.setSpace(original_metric_space);
    throw;
  }

  msh.met.setSpace(original_metric_space);
  return cost;
}

template double getMetricCost<MetricFieldAnalytical,2,2>(MeshMetric<MetricFieldAnalytical>& msh);
template double getMetricCost<MetricFieldAnalytical,3,3>(MeshMetric<MetricFieldAnalytical>& msh);
template double getMetricCost<MetricFieldFE,2,2>(MeshMetric<MetricFieldFE>& msh);
template double getMetricCost<MetricFieldFE,3,3>(MeshMetric<MetricFieldFE>& msh);

} // end namespace
