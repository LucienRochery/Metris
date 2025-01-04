//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_METQUA__
#define __METRIS_MSH_METQUA__
/* 
  Functions relating to Q_M(K) = 1/n tr(J_K^T J_0^{-T} M J_0^{-1} J_K) / det(J_K^T M J_K)^(1/n)
*/

#include "../Mesh/MeshFwd.hxx"
#include "../metris_constants.hxx"
#include "../types.hxx"

namespace Metris{

  
// -- Whole mesh qualities
// Prefer calling this one. Computes L^2 conformity error 
template <class MFT, int tdim, AsDeg = AsDeg::Pk>
double getmetquamesh(Mesh<MFT> &msh,
                     bool *iinva, double *qmin, double *qmax, 
                     double *qavg, dblAr1 *lquae);
// More options:
template <class MFT, int tdim, AsDeg = AsDeg::Pk>
double getmetquamesh(Mesh<MFT> &msh, int power, 
                     int pnorm, double difto, 
                     bool *iinva, double *qmin, double *qmax, 
                     double *qavg, dblAr1 *lquae);

} // End namespace

#endif