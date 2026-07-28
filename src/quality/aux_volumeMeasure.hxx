// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __AUX_VOLUMEMEASURE__
#define __AUX_VOLUMEMEASURE__

#ifdef STEPDISTANCE

#include "../linalg/symidx.hxx"
#include "../libs/SANS/Surreal/SurrealS.h"

namespace Metris{

namespace VolumeMeasureHelpers{

template<typename T>
inline double get_val(const T& x){
  return (double)x;
}

template<int nvar>
inline double get_val(const SANS::SurrealS<nvar,double>& x){
  return x.value();
}

template<int n, typename T>
T det_full(const T* A){
  if constexpr(n == 2){
    return A[0]*A[3] - A[1]*A[2];
  }else{
    return
        A[0]*(A[4]*A[8] - A[5]*A[7])
      - A[1]*(A[3]*A[8] - A[5]*A[6])
      + A[2]*(A[3]*A[7] - A[4]*A[6]);
  }
}

template<int n, typename T>
void inv_full(const T* A, T* Ai){
  T detA = det_full<n,T>(A);
  METRIS_ENFORCE(std::abs(get_val(detA)) > 1.0e-14);

  if constexpr(n == 2){
    Ai[0] =  A[3]/detA;
    Ai[1] = -A[1]/detA;
    Ai[2] = -A[2]/detA;
    Ai[3] =  A[0]/detA;
  }else{
    Ai[0] =  (A[4]*A[8] - A[5]*A[7])/detA;
    Ai[1] = -(A[1]*A[8] - A[2]*A[7])/detA;
    Ai[2] =  (A[1]*A[5] - A[2]*A[4])/detA;

    Ai[3] = -(A[3]*A[8] - A[5]*A[6])/detA;
    Ai[4] =  (A[0]*A[8] - A[2]*A[6])/detA;
    Ai[5] = -(A[0]*A[5] - A[2]*A[3])/detA;

    Ai[6] =  (A[3]*A[7] - A[4]*A[6])/detA;
    Ai[7] = -(A[0]*A[7] - A[1]*A[6])/detA;
    Ai[8] =  (A[0]*A[4] - A[1]*A[3])/detA;
  }
}

template<int gdim, typename T>
void sym_packed_to_full(const T* met, T* M){
  for(int ii = 0; ii < gdim*gdim; ii++) M[ii] = T(0);

  if constexpr(gdim == 2){
    M[0] = met[0];
    M[1] = met[1];
    M[2] = met[1];
    M[3] = met[2];
  }else{
    M[0] = met[0]; M[1] = met[1]; M[2] = met[3];
    M[3] = met[1]; M[4] = met[2]; M[5] = met[4];
    M[6] = met[3]; M[7] = met[4]; M[8] = met[5];
  }
}


// -----------------------------------------------------------------------------
// Frozen-metric theta:
//
// theta = sqrt(det(J^T J)) sqrt(det(M))
//
// Only J varies. M is treated as constant.
// Jreg_T is tdim x gdim, i.e. J^T, from ideal to physical
// -----------------------------------------------------------------------------
template<int gdim, int tdim, typename T>
void eval_theta_fixed_metric_grad(
    const T*__restrict__ Jreg_T,
    const T*__restrict__ met,
    const T*__restrict__ gradN,
    T*__restrict__ theta,
    T*__restrict__ dtheta)
{
  // G = J^T J = Jreg_T * Jreg_T^T.
  T G[tdim*tdim];
  for(int i = 0; i < tdim; i++){
    for(int j = 0; j < tdim; j++){
      G[i*tdim+j] = T(0);
      for(int a = 0; a < gdim; a++){
        G[i*tdim+j] += Jreg_T[i*gdim+a]*Jreg_T[j*gdim+a];
      }
    }
  }

  T M[gdim*gdim];
  sym_packed_to_full<gdim,T>(met, M);

  T detG = det_full<tdim,T>(G);
  T detM = det_full<gdim,T>(M);

  METRIS_ENFORCE(get_val(detG) > 0.0);
  METRIS_ENFORCE(get_val(detM) > 0.0);

  *theta = sqrt(detG)*sqrt(detM);

  if(dtheta == NULL) return;

  T Ginv[tdim*tdim];
  inv_full<tdim,T>(G, Ginv);

  T Ginv_gradN[tdim];
  for(int i = 0; i < tdim; i++){
    Ginv_gradN[i] = T(0);
    for(int j = 0; j < tdim; j++){
      Ginv_gradN[i] += Ginv[i*tdim+j]*gradN[j];
    }
  }

  // J G^{-1} gradN.
  // Since Jreg_T = J^T, component a is sum_i Jreg_T[i,a]*(Ginv_gradN)_i.
  for(int a = 0; a < gdim; a++){
    dtheta[a] = T(0);
    for(int i = 0; i < tdim; i++){
      dtheta[a] += Jreg_T[i*gdim+a]*Ginv_gradN[i];
    }
    dtheta[a] *= (*theta);
  }
}

// -----------------------------------------------------------------------------
// Build Surreal inputs for frozen-metric theta Hessian.
//
// Only Jreg_T varies:
//
// d(Jreg_T[i,a])/d(X_q)_k = gradN[i] delta_{a k}
//
// met is constant.
// gradN is constant.
// -----------------------------------------------------------------------------
template<int gdim, int tdim>
void eval_theta_fixed_metric_hess_by_surreal(
    const double*__restrict__ Jreg_T,
    const double*__restrict__ met,
    const double*__restrict__ gradN,
    double*__restrict__ htheta)
{
  constexpr int nnmet = (gdim*(gdim+1))/2;
  using S = SANS::SurrealS<gdim,double>;

  S Jreg_TS[tdim*gdim];
  S metS[nnmet];
  S gradNS[tdim];

  for(int i = 0; i < tdim; i++){
    for(int a = 0; a < gdim; a++){
      const int idx = i*gdim + a;

      Jreg_TS[idx].value() = Jreg_T[idx];

      for(int k = 0; k < gdim; k++){
        Jreg_TS[idx].deriv(k) = (a == k) ? gradN[i] : 0.0;
      }
    }
  }

  for(int im = 0; im < nnmet; im++){
    metS[im].value() = met[im];
    for(int k = 0; k < gdim; k++) metS[im].deriv(k) = 0.0;
  }

  for(int i = 0; i < tdim; i++){
    gradNS[i].value() = gradN[i];
    for(int k = 0; k < gdim; k++) gradNS[i].deriv(k) = 0.0;
  }

  S thetaS;
  S dthetaS[gdim];

  eval_theta_fixed_metric_grad<gdim,tdim,S>(
      Jreg_TS, metS, gradNS,
      &thetaS, dthetaS);

  for(int i = 0; i < gdim; i++){
    for(int j = i; j < gdim; j++){
      htheta[sym2idx(i,j)] = dthetaS[i].deriv(j);
    }
  }
}

// -----------------------------------------------------------------------------
// Metric-volume collapse barrier.
//
// rho = sqrt(det(J^T M J)) is the element volume in metric space, normalized
// by the regular reference simplex. The barrier is deliberately not multiplied
// by theta: beta * max(0,log(rho0/rho))^4 must diverge as rho -> 0.
//
// Only J varies here. The metric is frozen consistently with StepDistance.
// -----------------------------------------------------------------------------
template<int gdim, int tdim, typename T>
void eval_metric_volume_barrier_fixed_metric_grad(
    const T*__restrict__ Jreg_T,
    const T*__restrict__ met,
    const T*__restrict__ gradN,
    double rho0,
    double beta,
    T*__restrict__ rho,
    T*__restrict__ barrier,
    T*__restrict__ dbarrier)
{
  T M[gdim*gdim];
  sym_packed_to_full<gdim,T>(met,M);

  T A[tdim*tdim];
  for(int i = 0; i < tdim; i++){
    for(int j = 0; j < tdim; j++){
      A[i*tdim+j] = T(0);
      for(int a = 0; a < gdim; a++){
        for(int b = 0; b < gdim; b++){
          A[i*tdim+j] +=
              Jreg_T[i*gdim+a]*M[a*gdim+b]*Jreg_T[j*gdim+b];
        }
      }
    }
  }

  const T detA = det_full<tdim,T>(A);
  METRIS_ENFORCE(get_val(detA) > 0.0);
  *rho = sqrt(detA);
  *barrier = T(0);

  if(dbarrier != NULL){
    for(int a = 0; a < gdim; a++) dbarrier[a] = T(0);
  }
  if(beta <= 0.0 || rho0 <= 0.0 || get_val(*rho) >= rho0) return;

  const T h = log(T(rho0)/(*rho));
  const T h2 = h*h;
  *barrier = T(beta)*h2*h2;

  if(dbarrier == NULL) return;

  T Ainv[tdim*tdim];
  inv_full<tdim,T>(A,Ainv);

  // d log(rho) / d X_q = M J A^{-1} gradN.
  T Ainv_gradN[tdim];
  for(int i = 0; i < tdim; i++){
    Ainv_gradN[i] = T(0);
    for(int j = 0; j < tdim; j++){
      Ainv_gradN[i] += Ainv[i*tdim+j]*gradN[j];
    }
  }

  T J_Ainv_gradN[gdim];
  for(int a = 0; a < gdim; a++){
    J_Ainv_gradN[a] = T(0);
    for(int i = 0; i < tdim; i++){
      J_Ainv_gradN[a] += Jreg_T[i*gdim+a]*Ainv_gradN[i];
    }
  }

  for(int a = 0; a < gdim; a++){
    T dlogrho = T(0);
    for(int b = 0; b < gdim; b++){
      dlogrho += M[a*gdim+b]*J_Ainv_gradN[b];
    }
    dbarrier[a] = T(-4.0*beta)*h2*h*dlogrho;
  }
}

template<int gdim, int tdim>
void eval_metric_volume_barrier_fixed_metric_hess_by_surreal(
    const double*__restrict__ Jreg_T,
    const double*__restrict__ met,
    const double*__restrict__ gradN,
    double rho0,
    double beta,
    double*__restrict__ hbarrier)
{
  constexpr int nnmet = (gdim*(gdim+1))/2;
  using S = SANS::SurrealS<gdim,double>;

  S Jreg_TS[tdim*gdim];
  S metS[nnmet];
  S gradNS[tdim];

  for(int i = 0; i < tdim; i++){
    for(int a = 0; a < gdim; a++){
      const int idx = i*gdim+a;
      Jreg_TS[idx].value() = Jreg_T[idx];
      for(int k = 0; k < gdim; k++){
        Jreg_TS[idx].deriv(k) = (a == k) ? gradN[i] : 0.0;
      }
    }
  }
  for(int im = 0; im < nnmet; im++){
    metS[im].value() = met[im];
    for(int k = 0; k < gdim; k++) metS[im].deriv(k) = 0.0;
  }
  for(int i = 0; i < tdim; i++){
    gradNS[i].value() = gradN[i];
    for(int k = 0; k < gdim; k++) gradNS[i].deriv(k) = 0.0;
  }

  S rhoS;
  S barrierS;
  S dbarrierS[gdim];
  eval_metric_volume_barrier_fixed_metric_grad<gdim,tdim,S>(
      Jreg_TS,metS,gradNS,rho0,beta,&rhoS,&barrierS,dbarrierS);

  for(int i = 0; i < gdim; i++){
    for(int j = i; j < gdim; j++){
      hbarrier[sym2idx(i,j)] = dbarrierS[i].deriv(j);
    }
  }
}

} //namespace VolumeMeasureHelpers

} // namespace Metris

#endif // STEPDISTANCE

#endif // __AUX_VOLUMEMEASURE__
