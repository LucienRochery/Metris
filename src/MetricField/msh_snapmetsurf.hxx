//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_SNAPMETSURF__
#define __METRIS_MSH_SNAPMETSURF__



namespace Metris{

template<class MFT> class MeshMetric;

template<class MetricFieldType>
void snapMetSurf(MeshMetric<MetricFieldType> &msh,
                 int ithread);




}//namespace
#endif