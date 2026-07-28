
//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "quafun_stepDistance.hxx"

#include "../Mesh/Mesh.hxx"
#include "../metris_constants.hxx"
#include "../utils/aux_misc.hxx"

#include "../linalg/symidx.hxx"

#include "../utils/aux_pp_inc.hxx"
#include "../linalg/eigen.hxx"
#include "../low_eval.hxx"

#include "../libs/SANS/Surreal/SurrealS.h"

namespace Metris{

template<typename T>
T step_distance_power_from_squared(const T& distance2,
                                   double power,
                                   double regularization){
  const T eps2 = T(regularization*regularization);
  return pow(distance2 + eps2, power/2.0)
       - T(std::pow(regularization,power));
}

// For some special barys (nodes), met is already known -> pass it in
template <class MFT, int gdim, int tdim, typename ftype>
ftype quafun_stepDistance(Mesh<MFT> &msh,
                          AsDeg asdmsh, AsDeg asdmet,
                          const int*__restrict__ ent2pol,
                          const double*__restrict__ bary,
                          const double*__restrict__ met_){

  static_assert(gdim == 2 || gdim == 3);
  static_assert(tdim <= gdim);

  METRIS_ASSERT(gdim == msh.idim);
  METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);

  constexpr int nnmet_tdim = (tdim*(tdim+1))/2;
  constexpr int nnmet_gdim = (gdim*(gdim+1))/2;

  // ------------------------------------------------------------
  // Geometry: canonical-reference Jacobian.
  // jmat is tdim x gdim, i.e. Jcanonical^T.
  // ------------------------------------------------------------
  double coopr[gdim];
  double jmat[tdim*gdim];

  if constexpr(tdim == 2){
    eval2<gdim,1>(msh.coord, ent2pol, msh.getBasis(),
                  DifVar::Bary, DifVar::None,
                  bary, coopr, jmat, NULL);
  }else{
    eval3<gdim,1>(msh.coord, ent2pol, msh.getBasis(),
                  DifVar::Bary, DifVar::None,
                  bary, coopr, jmat, NULL);
  }

  // ------------------------------------------------------------
  // Jreg_T = Jreg^T = J0^{-T} Jcanonical^T.
  // Shape: tdim x gdim.
  // ------------------------------------------------------------
  ftype Jreg_T[tdim*gdim];

  for(int i = 0; i < tdim; i++){
    for(int a = 0; a < gdim; a++){
      Jreg_T[i*gdim+a] = ftype(0);
      for(int k = 0; k < tdim; k++){
        Jreg_T[i*gdim+a] +=
          Constants::invtJ_0[hana::type_c<ftype>][tdim][i*tdim+k]
          * (ftype)jmat[k*gdim+a];
      }
    }
  }

  // ------------------------------------------------------------
  // Metric.
  // In the frozen-metric metqua branch, met_ is passed in.
  // Keep fallback for standalone use.
  // ------------------------------------------------------------
  double met_local[nnmet_gdim];
  const double* met = met_;

  if(met == NULL){
    msh.met.getMetFullinfo(asdmet, DifVar::None, MetSpace::Exp,
                           ent2pol, tdim, bary, coopr,
                           met_local, NULL);
    met = met_local;
  }

  // ------------------------------------------------------------
  // A = J^T M J.
  //
  // Since Jreg_T = J^T, A = Jreg_T * M * Jreg_T^T.
  // Store A full tdim x tdim.
  // ------------------------------------------------------------
  ftype A[tdim*tdim];

  for(int i = 0; i < tdim; i++){
    for(int j = 0; j < tdim; j++){
      A[i*tdim+j] = ftype(0);

      for(int a = 0; a < gdim; a++){
        for(int b = 0; b < gdim; b++){

          int imet;
          if constexpr(gdim == 2){
            imet = (a == 0 && b == 0) ? 0 :
                   (a != b)           ? 1 : 2;
          }else{
            imet = (a == 0 && b == 0) ? 0 :
                   (a == 0 && b == 1) ? 1 :
                   (a == 1 && b == 0) ? 1 :
                   (a == 1 && b == 1) ? 2 :
                   (a == 0 && b == 2) ? 3 :
                   (a == 2 && b == 0) ? 3 :
                   (a == 1 && b == 2) ? 4 :
                   (a == 2 && b == 1) ? 4 : 5;
          }

          A[i*tdim+j] += Jreg_T[i*gdim+a] * (ftype)met[imet] * Jreg_T[j*gdim+b];
        }
      }
    }
  }

  // ------------------------------------------------------------
  // Pack A for geteigsym.
  // Packing:
  // 2D: [a00, a01, a11]
  // 3D: [a00, a01, a11, a02, a12, a22]
  // ------------------------------------------------------------
  ftype Ap[nnmet_tdim];

  if constexpr(tdim == 2){
    Ap[0] = A[0];
    Ap[1] = A[1];
    Ap[2] = A[3];
  }else{
    Ap[0] = A[0];
    Ap[1] = A[1];
    Ap[2] = A[4];
    Ap[3] = A[2];
    Ap[4] = A[5];
    Ap[5] = A[8];
  }

  ftype eigval[tdim];
  ftype eigvec[tdim*tdim];

  geteigsym<tdim,ftype>(Ap, eigval, eigvec);

  // ------------------------------------------------------------
  // phi = ||log(A)||_F^p. A small smooth norm regularization is used so
  // p < 2 remains differentiable at the ideal state A = I.
  // ------------------------------------------------------------
  ftype distance2 = ftype(0);

  for(int i = 0; i < tdim; i++){
    METRIS_ENFORCE(eigval[i] > ftype(0));
    ftype li = log(eigval[i]);
    distance2 += li*li;
  }

  return step_distance_power_from_squared(
      distance2,
      msh.param->step_distance_p,
      msh.param->step_distance_regularization);
}


#define EXPAND_TEMPLATE(r,SEQ) \
                  INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ),\
                              BOOST_PP_SEQ_ELEM(1, SEQ))
#define MFT_SEQ (MetricFieldFE)(MetricFieldAnalytical)




#define INSTANTIATE(MFT_VAL,FTYPE)\
template FTYPE quafun_stepDistance< MFT_VAL , 2, 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary,\
                   const double*__restrict__ met); \
template FTYPE quafun_stepDistance< MFT_VAL , 3, 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary,\
                   const double*__restrict__ met); \
template FTYPE quafun_stepDistance< MFT_VAL , 3, 3,  FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary,\
                   const double*__restrict__ met);
BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TEMPLATE,\
                              (MFT_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE
#undef EXPAND_TEMPLATE

// =======================================//
// Some helpers for derivatives
// =======================================//

template<int gdim, int tdim>
void get_gradN_P1_regular(int ivar, double* gradN){
  for(int i = 0; i < tdim; i++) gradN[i] = 0.0;

  if(ivar == 0){
    for(int i = 0; i < tdim; i++){
      for(int k = 0; k < tdim; k++){
        gradN[i] -= Constants::invtJ_0[hana::type_c<double>][tdim][i*tdim+k];
      }
    }
  }else{
    const int column = ivar - 1;
    for(int i = 0; i < tdim; i++){
      gradN[i] =
          Constants::invtJ_0[hana::type_c<double>][tdim][i*tdim+column];
    }
  }
}

template<int gdim, typename T>
int symidx_met(int a, int b){
  if(a > b) std::swap(a,b);

  if constexpr(gdim == 2){
    if(a == 0 && b == 0) return 0;
    if(a == 0 && b == 1) return 1;
    return 2;
  }else{
    if(a == 0 && b == 0) return 0;
    if(a == 0 && b == 1) return 1;
    if(a == 1 && b == 1) return 2;
    if(a == 0 && b == 2) return 3;
    if(a == 1 && b == 2) return 4;
    return 5;
  }
}

template<int tdim, typename T>
void sym_full_to_packed_tdim(const T* A, T* Ap){
  if constexpr(tdim == 2){
    Ap[0] = A[0];
    Ap[1] = A[1];
    Ap[2] = A[3];
  }else{
    Ap[0] = A[0];
    Ap[1] = A[1];
    Ap[2] = A[4];
    Ap[3] = A[2];
    Ap[4] = A[5];
    Ap[5] = A[8];
  }
}

template<int tdim, typename T>
void spectral_matrix_from_eigs(const T* eigval,
                               const T* eigvec,
                               const T* fval,
                               T* F){
  for(int i = 0; i < tdim*tdim; i++) F[i] = T(0);

  for(int a = 0; a < tdim; a++){
    for(int i = 0; i < tdim; i++){
      for(int j = 0; j < tdim; j++){
        F[i*tdim+j] += fval[a]*eigvec[a*tdim+i]*eigvec[a*tdim+j];
      }
    }
  }
}

template<int gdim, int tdim, typename T>
void eval_phi_grad_impl(const T* Jreg_T,
                        const T* met,
                        const T* gradN,
                        double power,
                        double regularization,
                        T* phi,
                        T* dphi){

  constexpr int nnmet_tdim = (tdim*(tdim+1))/2;

  // A = J^T M J = Jreg_T M Jreg_T^T
  T A[tdim*tdim];

  for(int i = 0; i < tdim; i++){
    for(int j = 0; j < tdim; j++){
      A[i*tdim+j] = T(0);

      for(int a = 0; a < gdim; a++){
        for(int b = 0; b < gdim; b++){
          const int imet = symidx_met<gdim,T>(a,b);
          A[i*tdim+j] += Jreg_T[i*gdim+a]*met[imet]*Jreg_T[j*gdim+b];
        }
      }
    }
  }

  T Ap[nnmet_tdim];
  sym_full_to_packed_tdim<tdim,T>(A, Ap);

  T eigval[tdim];
  T eigvec[tdim*tdim];

  geteigsym<tdim,T>(Ap, eigval, eigvec);

  T loglam[tdim];
  T invloglam[tdim];

  T distance2 = T(0);

  for(int i = 0; i < tdim; i++){
    loglam[i] = log(eigval[i]);
    invloglam[i] = loglam[i]/eigval[i];
    distance2 += loglam[i]*loglam[i];
  }

  *phi = step_distance_power_from_squared(
      distance2,power,regularization);

  if(dphi == NULL) return;

  const T eps2 = T(regularization*regularization);
  const T power_scale = T(power/2.0)
      *pow(distance2 + eps2,power/2.0 - 1.0);

  // A^{-1} log(A)
  T AinvlogA[tdim*tdim];
  spectral_matrix_from_eigs<tdim,T>(eigval, eigvec, invloglam, AinvlogA);

  // v = A^{-1}log(A) gradN
  T v[tdim];
  for(int i = 0; i < tdim; i++){
    v[i] = T(0);
    for(int j = 0; j < tdim; j++){
      v[i] += AinvlogA[i*tdim+j]*gradN[j];
    }
  }

  // u = J v, using Jreg_T = J^T
  T u[gdim];
  for(int a = 0; a < gdim; a++){
    u[a] = T(0);
    for(int i = 0; i < tdim; i++){
      u[a] += Jreg_T[i*gdim+a]*v[i];
    }
  }

  // dphi = 4 M u
  for(int a = 0; a < gdim; a++){
    dphi[a] = T(0);
    for(int b = 0; b < gdim; b++){
      const int imet = symidx_met<gdim,T>(a,b);
      dphi[a] += met[imet]*u[b];
    }
    // The final factor is (p/2)(distance2 + eps^2)^(p/2 - 1)
    // times the p=2 derivative of distance2.
    dphi[a] *= T(4)*power_scale;
  }
}

/*
NOTE: This does not yet use the metric field derivatives.
That is something we might not want, it'll be costlier as we can no longer
provide a metric, but instead must interpolate (which dominates CPU time
because of the matrix exponential).
Compute quality function and derivative with respect to node ivar
- gdim is geometric dimension: also topological dimension !
- mshdeg is mesh degree
- asdmet is AsDeg::P1 or AsDeg::Pk: MetricField handles its degree
- ftype is arithmetic precision (debug): double, float8...
- ivar is the DoF, specified Bézier or Lagrange by:
- dofbas is either FEBasis::Lagrange or FEBasis::Bezier -> whether the DoF is a
  Lagrange node or Bézier control point
- idifmet is either DifVar::None ("static" metric approximation) or DifVar::Phys
- quael is output quality and derivatives
*/
template <class MFT, int gdim, int tdim, typename ftype>
ftype d_quafun_stepDistance(Mesh<MFT> &msh,
                  AsDeg asdmsh, AsDeg asdmet,
                  const int* ent2pol,
                  const double*__restrict__ bary,
                  const double*__restrict__ met_,
                  int ivar,
                  FEBasis dofbas,
                  DifVar idifmet,
                  ftype*__restrict__ dquael,
                  ftype*__restrict__ hquael){

  static_assert(gdim == 2 || gdim == 3);
  static_assert(tdim <= gdim);

  METRIS_ASSERT(gdim == msh.idim);
  METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);

  constexpr int nnmet = (gdim*(gdim+1))/2;
  constexpr int nhess = (gdim*(gdim+1))/2;

  // For now: P1 only.
  METRIS_ASSERT(msh.curdeg == 1 || asdmsh == AsDeg::P1);
  METRIS_ASSERT(ivar < 0 || (ivar >= 0 && ivar < tdim + 1));

  // ------------------------------------------------------------
  // Geometry value: canonical-reference Jacobian.
  // jmat is tdim x gdim.
  // ------------------------------------------------------------
  double coopr[gdim];
  double jmat[tdim*gdim];

  if constexpr(tdim == 2){
    eval2<gdim,1>(msh.coord, ent2pol, msh.getBasis(),
                  DifVar::Bary, DifVar::None,
                  bary, coopr, jmat, NULL);
  }else{
    eval3<gdim,1>(msh.coord, ent2pol, msh.getBasis(),
                  DifVar::Bary, DifVar::None,
                  bary, coopr, jmat, NULL);
  }

  // ------------------------------------------------------------
  // Jreg_T = Jreg^T = J0^{-T} Jcanonical^T.
  // Shape: tdim x gdim.
  // ------------------------------------------------------------
  double Jreg_T[tdim*gdim];

  for(int i = 0; i < tdim; i++){
    for(int a = 0; a < gdim; a++){
      Jreg_T[i*gdim+a] = 0.0;
      for(int k = 0; k < tdim; k++){
        Jreg_T[i*gdim+a] +=
          Constants::invtJ_0[hana::type_c<double>][tdim][i*tdim+k]
          *jmat[k*gdim+a];
      }
    }
  }

  // ------------------------------------------------------------
  // Frozen metric.
  // Usually passed by metqua. Fallback retained for standalone calls.
  // ------------------------------------------------------------
  double met_local[nnmet];
  const double* met = met_;

  if(met == NULL){
    msh.met.getMetFullinfo(asdmet, DifVar::None, MetSpace::Exp,
                           ent2pol, tdim, bary, coopr,
                           met_local, NULL);
    met = met_local;
  }

  // ------------------------------------------------------------
  // Value-only path.
  // ------------------------------------------------------------
  if(ivar < 0){
    ftype Jreg_T_f[tdim*gdim];
    ftype met_f[nnmet];

    for(int i = 0; i < tdim*gdim; i++) Jreg_T_f[i] = (ftype)Jreg_T[i];
    for(int i = 0; i < nnmet; i++)     met_f[i]    = (ftype)met[i];

    ftype phi;
    eval_phi_grad_impl<gdim,tdim,ftype>(
        Jreg_T_f, met_f, NULL,
        msh.param->step_distance_p,
        msh.param->step_distance_regularization,
        &phi, NULL);

    return phi;
  }

  // ------------------------------------------------------------
  // P1 gradient of active shape function in regular reference.
  // ------------------------------------------------------------
  double gradN_d[tdim];
  get_gradN_P1_regular<gdim,tdim>(ivar, gradN_d);

  // ------------------------------------------------------------
  // First derivative: analytic expression with double/ftype.
  // ------------------------------------------------------------
  ftype Jreg_T_f[tdim*gdim];
  ftype met_f[nnmet];
  ftype gradN_f[tdim];

  for(int i = 0; i < tdim*gdim; i++) Jreg_T_f[i] = (ftype)Jreg_T[i];
  for(int i = 0; i < nnmet; i++)     met_f[i]    = (ftype)met[i];
  for(int i = 0; i < tdim; i++)      gradN_f[i]  = (ftype)gradN_d[i];

  ftype phi;
  ftype dphi[gdim];

  eval_phi_grad_impl<gdim,tdim,ftype>(
      Jreg_T_f, met_f, gradN_f,
      msh.param->step_distance_p,
      msh.param->step_distance_regularization,
      &phi, dphi);

  for(int i = 0; i < gdim; i++){
    dquael[i] = dphi[i];
  }

  // ------------------------------------------------------------
  // Hessian: AD of analytic gradient.
  //
  // Only J varies. Metric is frozen.
  // d(Jreg_T[i,a])/d(X_q)_k = gradN[i] delta_{ak}.
  // ------------------------------------------------------------
  if(hquael != NULL){

    using S = SANS::SurrealS<gdim,double>;

    S Jreg_TS[tdim*gdim];
    S metS[nnmet];
    S gradNS[tdim];

    for(int i = 0; i < tdim; i++){
      for(int a = 0; a < gdim; a++){
        const int idx = i*gdim + a;

        Jreg_TS[idx].value() = Jreg_T[idx];

        for(int k = 0; k < gdim; k++){
          Jreg_TS[idx].deriv(k) = (a == k) ? gradN_d[i] : 0.0;
        }
      }
    }

    for(int im = 0; im < nnmet; im++){
      metS[im].value() = met[im];
      for(int k = 0; k < gdim; k++){
        metS[im].deriv(k) = 0.0;
      }
    }

    for(int i = 0; i < tdim; i++){
      gradNS[i].value() = gradN_d[i];
      for(int k = 0; k < gdim; k++){
        gradNS[i].deriv(k) = 0.0;
      }
    }

    S phiS;
    S dphiS[gdim];

    eval_phi_grad_impl<gdim,tdim,S>(
        Jreg_TS, metS, gradNS,
        msh.param->step_distance_p,
        msh.param->step_distance_regularization,
        &phiS, dphiS);

    for(int i = 0; i < gdim; i++){
      for(int j = i; j < gdim; j++){
        hquael[sym2idx(i,j)] = (ftype)dphiS[i].deriv(j);
      }
    }
  }

  return phi;
}

// While cumbersome, this replaces a bunch of manual instantiations, about to
// be made worse the day we add tdimn as a template argument.
#define EXPAND_TEMPLATE(z,SEQ) \
                  INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ),\
                              BOOST_PP_SEQ_ELEM(1, SEQ))

#define INSTANTIATE(MFT_VAL,FTYPE)\
template FTYPE d_quafun_stepDistance< MFT_VAL , 2, 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int* ent2pol, \
                   const double*__restrict__ bary,\
                   const double*__restrict__ met_,\
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, \
                   FTYPE*__restrict__ hquael);\
template FTYPE d_quafun_stepDistance< MFT_VAL , 3, 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int* ent2pol, \
                   const double*__restrict__ bary,\
                   const double*__restrict__ met_,\
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, \
                   FTYPE*__restrict__ hquael);\
template FTYPE d_quafun_stepDistance< MFT_VAL , 3, 3, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int* ent2pol, \
                   const double*__restrict__ bary,\
                   const double*__restrict__ met_,\
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, \
                   FTYPE*__restrict__ hquael);
BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TEMPLATE,\
                              (MFT_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE
#undef EXPAND_TEMPLATE

#undef MFT_SEQ
#undef ASDEG_SEQ


}
