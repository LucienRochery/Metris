//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_CHECKMET__
#define __METRIS_MSH_CHECKMET__


#include "../Mesh/MeshFwd.hxx"


namespace Metris{
class MetricFieldFE;


void checkMet(const MeshMetric<MetricFieldFE>& msh);


} // end namespace

#endif