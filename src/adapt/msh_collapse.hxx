//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_COLLAPSE__
#define __METRIS_MSH_COLLAPSE__


#include "../Mesh/MeshFwd.hxx"


namespace Metris{


template<class MetricFieldType, int gdim, int ideg>
double collapseShortEdges(Mesh<MetricFieldType> &msh, double qmax_suf,
                          int ithrd1 = 0, int ithrd2 = 1, int ithrd3 = 2);


}// end Namespace

#endif