//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_INTERPFRONTBACK__
#define __METRIS_MSH_INTERPFRONTBACK__


#include "../Mesh/MeshFwd.hxx"


namespace Metris{
template<class MetricFieldType, int bdeg>
void interpFrontBack(Mesh<MetricFieldType> &msh, MeshBack &bak, int ipoi0 = 0);



} // End namespace

#endif
