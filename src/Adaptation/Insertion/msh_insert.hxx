//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_INSERT2D__
#define __METRIS_MSH_INSERT2D__


#include "../../Mesh/MeshFwd.hxx"


namespace Metris{


template<class MetricFieldType, int gdim, int ideg>
double insertLongEdges(Mesh<MetricFieldType> &msh, int tdim, int *ninser,
                       int ithrd1, int ithrd2);


}// end Namespace

#endif