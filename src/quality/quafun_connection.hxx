// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#ifndef __METRIS_QUAFUN_CONNECTION__
#define __METRIS_QUAFUN_CONNECTION__

#include "../Mesh/MeshFwd.hxx"

namespace Metris{

enum class FEBasis;
enum class AsDeg;
enum class DifVar;

// Pointwise squared target-metric norm of the connection difference:
//   T^k_ij = Gamma(M_K)^k_ij - Gamma(M)^k_ij,
//   ||T||_M^2 = M_kl M^{ip} M^{jq} T^k_ij T^l_pq.
// M_K is induced by a constant reference metric and the element map F.
// Full-dimensional elements only (gdim = 2 or 3).
//
// Array conventions, with nmet = gdim*(gdim+1)/2:
//   jmat[a*gdim+k]              = d_a F^k (canonical reference coordinates)
//   hmat[sym2idx(a,b)*gdim+k]   = d_a d_b F^k
//   met[sym2idx(i,j)]           = M_ij (SPD, NOT its logarithm)
//   dmet[r*nmet+sym2idx(i,j)]   = d M_ij / d x^r (physical coordinates).
// Rejects singular Jacobians, non-SPD metrics, and nonfinite input/results
// by throwing. This is a value evaluator, not a trial-rejection policy.
template<int gdim>
double connection_difference_norm2(const double *jmat, const double *hmat,
                                   const double *met, const double *dmet);

// This uses the above, after using some helpers so compute its arguments at the evaluation point
template<class MetricFieldType, int gdim>
double quafun_connection(Mesh<MetricFieldType> &msh, const int *ent2poi,
                         const double *bary);

template<int gdim>
double d_connection_difference_norm2(
    const double *jmat, const double *hmat,
    const double *met, const double *dmet,
    const double *d_jmat, const double *d_hmat, double *gradient);

template<class MetricFieldType, int gdim>
double d_quafun_connection(Mesh<MetricFieldType> &msh, const int *ent2poi,
                           const double *bary, int ivar, FEBasis dofbas,
                           double *gradient);

template<int gdim>
double d2_connection_difference_norm2_frozen(
    const double *jmat, const double *hmat, const double *met, const double *dmet,
    const double *dphi, const double *hphi, double *gradient, double *hessian);

// Integrated wrappers for production node-by-node smoothing.
// These are needed because connection objective does not follow
// same policies as the distance objectives as described in IMR paper.
//
// Both gradient and Hessian freeze the full
// quadrature weight, M, and dM/dx; trial values re-evaluate all three.
// Full-dimensional Pk geometry/metric, double precision, raw squared-norm
// integral only (no distance aggregation or exponent transformation).
template<class MFT, int gdim>
double metqua_connection(Mesh<MFT>& msh, AsDeg asdmsh, AsDeg asdmet,
                         int element, double difto = 1.);
template<class MFT, int gdim>
double d_metqua_connection(Mesh<MFT>& msh, AsDeg asdmsh, AsDeg asdmet,
                           int element, int ivar, FEBasis dofbas, DifVar idifmet,
                           double *gradient, double *hessian, double difto = 1.);

} // namespace Metris

#endif
