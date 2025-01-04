//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_QUAFUN_TRADET__
#define __METRIS_QUAFUN_TRADET__
/* 
  Auxiliary functions to compute:
   - tr(J_K^T J_0^{-T} M J_0^{-1} J_K)
   - det(J_K^T M J_K)^(1/n)
*/


#include "../Mesh/MeshFwd.hxx"


namespace Metris{

enum class AsDeg;
enum class FEBasis;
enum class DifVar;

/* ---- Scale invariant distortion measure form tr / det ---- */
// Pointwise 
template <class MetricFieldType, int gdim, int tdim, 
          typename ftype = double>
void quafun_tradet(Mesh<MetricFieldType> &msh,
                   AsDeg asdmsh, AsDeg asdmet,
                   const int*__restrict__ ent2poi,  
                   const double*__restrict__ bary, int power, 
                   double*__restrict__ met,
                   ftype*__restrict__ tra,
                   ftype*__restrict__ det);

// Differentiated w.r.t. ielem's ivar-th control point/node. 
template <class MetricFieldType, int gdim, 
           typename ftype = double>
void d_quafun_tradet(Mesh<MetricFieldType> &msh,
                     AsDeg asdmsh, AsDeg asdmet,
                     const int* ent2poi, 
                     const double*__restrict__ bary, 
                     int power, 
                     int ivar,
                     FEBasis dofbas, 
                     DifVar idifmet, 
                     ftype*__restrict__ tra, 
                     ftype*__restrict__ dtra, 
                     ftype*__restrict__ htra, 
                     ftype*__restrict__ det, 
                     ftype*__restrict__ ddet, 
                     ftype*__restrict__ hdet);


} // End namespace

#endif
