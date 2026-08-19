//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
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
#include "../quality/low_metqua.hxx"

namespace Metris{


// -- Whole mesh qualities
// SizeShape and StepDistance return their complete additive objective (or the
// explicitly selected StepDistance mean). Classical Distortion and Unit retain
// the historical opt_pnorm aggregation.
template <class MFT, QuaFun iquaf = QuaFun::Distortion>
double getmetquamesh(Mesh<MFT> &msh, int tdim, AsDeg asdmsh, AsDeg asdmet,
                     bool *iinva, double *qmin, double *qmax,
                     double *qavg, dblAr1 *lquae);
} // End namespace

#endif
