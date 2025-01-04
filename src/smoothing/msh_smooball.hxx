//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

/*
Ball-based as-P1 smoothing. 
Low-level drivers:
 - low_smoodirect -> facets propose points, average
 - low_smooballnewton -> Newton 
*/

#ifndef __METRIS_MSH_SMOOBAL__
#define __METRIS_MSH_SMOOBAL__



#include "../Mesh/MeshFwd.hxx" 
#include "../quality/quafun.hxx"

namespace Metris{

// Returns number of operations as double; this is because it may exceed element
// count, something that cannot be anticipated easily by selecting integer
// types. 
template<class MetricFieldType>
double smoothInterior_Ball(Mesh<MetricFieldType> &msh, 
             QuaFun iquaf = QuaFun::Distortion, int ithrd1 = 0, int ithrd2 = 1);


// idim: gdim = tdim
template<class MetricFieldType, int idim, int ideg>
double smoothInterior_Ball0(Mesh<MetricFieldType> &msh, 
             QuaFun iquaf = QuaFun::Distortion, int ithrd1 = 0, int ithrd2 = 1);




} // end namespace
#endif