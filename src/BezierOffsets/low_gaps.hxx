//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LOW_GAPS__
#define __METRIS_LOW_GAPS__


#include "../Mesh/MeshFwd.hxx"


namespace Metris{
//void scalrot3(const Mesh &msh, int ielem , 
//              double* __restrict__ metbar,
//              double* __restrict__ dmetbar,
//              double* __restrict__ scale ,
//              double* __restrict__ rotmat);

template<class MetricFieldType, int gdim, int ideg>
void getBezOffsetsEdge(Mesh<MetricFieldType> &msh, int tdimn, 
                       const int* ent2poi, int iedgl, double* offsets);

//template<int gdim, int ideg>
//void getBezOffsetsEdge0(const Mesh<MetricFieldAnalytical> &msh, MeshBack &bak, 
//  int ientt, int tdimn, double *__restrict__ dedg, double*__restrict__ bary, double*__restrict__ offset);
//template<int gdim, int ideg>
//void getBezOffsetsEdge0(const Mesh<MetricFieldFE> &msh, MeshBack &bak, 
//  int ientt, int tdimn, double *__restrict__ dedg, double*__restrict__ bary, double*__restrict__ offset);


//void d2unittensor(const Mesh<MetricFieldAnalytical> &msh, int ielem, double *tens3sym);
//void d2unittensor2(const MeshBase &msh, int ielem, double *tens3sym_);


// Common aux 
template <int gdim>
void scalrotJ0(const MeshBase &msh, int ielem  , 
               const double* __restrict__ srmet ,
               double* __restrict__ scale  ,
               double* __restrict__ rotJ0  );


} // End namespace

#endif