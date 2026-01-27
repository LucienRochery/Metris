
//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_QUAFUN_HACKCOSTFUNCTIONAL__
#define __METRIS_QUAFUN_HACKCOSTFUNCTIONAL__
/*
  Functions relating to Q_M(K)
   = 1/n tr(J_0^{-T} J_K^T  M  J_K J_0^{-1}) / det(J_0^{-T} J_K^T  M  J_K J_0^{-1})^(1/n) * (1/2 * ( det(J_K^T M J_K) + 1/det(J_K^T M J_K) ) )^(2/n)
*/

#include "../Mesh/MeshFwd.hxx"

namespace Metris{

enum class AsDeg;
enum class FEBasis;
enum class DifVar;

/* ---- Just returns 2 so that the integrand in metqua is 1 and so it integrates just sqrt(det) ---- */
template <class MetricFieldType, int gdim, int tdim,
          typename ftype = double>
ftype quafun_hackCostFunctional(Mesh<MetricFieldType> &msh,
                                AsDeg asdmsh, AsDeg asdmet,
                                const int*__restrict__ ent2poi,
                                const double*__restrict__ bary,
                                const double*__restrict__ met);

} // End namespace

#endif // __METRIS_QUAFUN_HACKCOSTFUNCTIONAL__
