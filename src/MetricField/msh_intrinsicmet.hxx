//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_INTRINSICMET__
#define __METRIS_MSH_INTRINSICMET__


#include "../types_arrays.hxx"

#include "../Mesh/MeshFwd.hxx"

namespace Metris{

struct MetrisParameters;
enum class FEBasis;

// -----------------------------------------------------------------------------
template<class MetricFieldType,int ideg>
void getMetMesh(const MetrisParameters &param, MeshMetric<MetricFieldType> &msh);


template<class MetricFieldType, int gdim, int tdim, int ideg>
void getMetMesh0_lplib(int ient0, int ient1,int ithread, 
                       MeshMetric<MetricFieldType> *msh_, dblWrkAr1 *rwork, int poitag);

// -----------------------------------------------------------------------------
template<int gdim, int tdim, int ideg>
int getintmetxi(const dblAr2 &coord, const int* __restrict__ tet2pol, FEBasis ibasis,
	               const double* bary,double* __restrict__ met);



} // End namespace


#endif

