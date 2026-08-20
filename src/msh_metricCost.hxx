//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php
#ifndef __METRIS_MSH_METRICCOST__
#define __METRIS_MSH_METRICCOST__

#include "types.hxx"
#include "MetricField/MetricField.hxx"
#include "Mesh/MeshMetric.hxx"

namespace Metris{

template<class MFT,int gdim, int tdim>
double getMetricCost(MeshMetric<MFT> &msh);

}// end namespace
#endif
