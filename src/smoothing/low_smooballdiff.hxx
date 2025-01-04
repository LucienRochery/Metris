//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

/*
Low level routine for "direct" P1 ball smoothing. 
From each (facet, metric) pair, generate remaining vertex to be unit. Then average over ball. 
Simplest possible approach. 
*/

#ifndef __METRIS_LOW_SMOOBALLDIFF__
#define __METRIS_LOW_SMOOBALLDIFF__

#include "../Mesh/MeshFwd.hxx"
#include "../quality/quafun.hxx"
#include "../types.hxx"

namespace Metris{



// qball is work size nball
// same as smooballdirect but gradient descent 
template<class MetricFieldType, int idim, int ideg>
int smooballdiff(Mesh<MetricFieldType>& msh, int ipoin,
                 const intAr1 &lball, dblAr1 &qball,
                 double*__restrict__ qavg0, double*__restrict__ qmax0, 
                 double*__restrict__ qavg1, double*__restrict__ qmax1,
                 QuaFun iquaf = QuaFun::Distortion);



} // end namespace
#endif