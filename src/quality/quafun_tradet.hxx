//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
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
#include "../../SANS/Surreal/SurrealS_fwd.h"

namespace SANS::DLA{
  template<int n, int m, typename T> class MatrixS;
}

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
                   const double*__restrict__ bary,
                   const double*__restrict__ met_,
                   ftype*__restrict__ tra,
                   ftype*__restrict__ det);

// Differentiated w.r.t. ielem's ivar-th control point/node. 
template <class MetricFieldType, int gdim, int tdim,
           typename ftype = double>
void d_quafun_tradet(Mesh<MetricFieldType> &msh,
                     AsDeg asdmsh, AsDeg asdmet,
                     const int* ent2poi, 
                     const double*__restrict__ bary, 
                     int ivar,
                     FEBasis dofbas, 
                     DifVar idifmet, 
                     const double*__restrict__ met_,
                     ftype*__restrict__ tra, 
                     ftype*__restrict__ dtra, 
                     ftype*__restrict__ htra, 
                     ftype*__restrict__ det, 
                     ftype*__restrict__ ddet, 
                     ftype*__restrict__ hdet);


template <class MFT, int gdim, int tdim, int nvar>
void d_quafun_tradet_SurrealS(Mesh<MFT> &msh, AsDeg asdmsh, AsDeg asdmet,
                              const int* ent2pol,
                              const double*__restrict__ bary, 
                              int ivar, 
                              FEBasis dofbas, 
                              DifVar idifmet, 
                              const double*__restrict__ met_,
                              SANS::SurrealS<nvar, double>&__restrict__ tra, 
                              SANS::SurrealS<nvar, double>&__restrict__ det,
                              const SANS::DLA::MatrixS<gdim,nvar,double> *dpoint);


} // End namespace

#endif
