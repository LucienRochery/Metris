
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

#include <algorithm>
#include <cmath>

namespace Metris{

namespace{

// ShapeVolume is deliberately singular as det(A) -> 0.  Do not pass a
// numerically singular trial element to the spectral value/derivative code:
// cavity operations must see a large, finite quality and reject the trial in
// their normal rollback path.
constexpr double shape_volume_minimum_determinant = 1.0e-12;
constexpr double shape_volume_minimum_scaled_determinant = 1.0e-12;

template<typename T>
double step_distance_primal_value(const T& value){
  return static_cast<double>(value);
}

template<int nvar>
double step_distance_primal_value(
    const SANS::SurrealS<nvar,double>& value){
  return value.value();
}

template<int tdim, typename T>
bool shape_volume_matrix_is_numerically_admissible(
    const T*__restrict__ A){
  long double matrix[tdim*tdim];
  for(int i = 0; i < tdim*tdim; i++){
    matrix[i] = step_distance_primal_value(A[i]);
    if(!std::isfinite(matrix[i])) return false;
  }

  long double trace = 0.;
  for(int i = 0; i < tdim; i++) trace += matrix[i*tdim+i];
  if(!(trace > 0.) || !std::isfinite(trace)) return false;

  long double determinant;
  if constexpr(tdim == 2){
    determinant = matrix[0]*matrix[3] - matrix[1]*matrix[2];
  }else{
    determinant =
          matrix[0]*(matrix[4]*matrix[8] - matrix[5]*matrix[7])
        - matrix[1]*(matrix[3]*matrix[8] - matrix[5]*matrix[6])
        + matrix[2]*(matrix[3]*matrix[7] - matrix[4]*matrix[6]);
  }
  if(!(determinant >= shape_volume_minimum_determinant)
     || !std::isfinite(determinant)) return false;

  const long double mean_eigenvalue = trace/tdim;
  const long double scaled_determinant =
      determinant/std::pow(mean_eigenvalue,tdim);
  return std::isfinite(scaled_determinant)
      && scaled_determinant >= shape_volume_minimum_scaled_determinant;
}

} // anonymous namespace

bool objective_strictly_improves(double value_new,
                                 double value_old,
                                 double relative_improvement){
  METRIS_ENFORCE(relative_improvement >= 0.);
  return value_new <= (1. - relative_improvement)*value_old;
}

bool cavity_target_average_global_filter_accepts(
    double local_old,
    double local_new,
    double global_current,
    double global_trial,
    double global_best,
    double old_region_unit_weight,
    double new_region_unit_weight,
    double global_unit_weight,
    double global_relative_tolerance,
    double global_gain_fraction){
  (void)local_old;
  (void)local_new;
  (void)global_best;
  (void)old_region_unit_weight;
  (void)new_region_unit_weight;
  (void)global_unit_weight;
  (void)global_relative_tolerance;
  (void)global_gain_fraction;
  return objective_strictly_improves(global_trial,global_current);
}

double StepDistanceObjectiveState::value() const{
  return step_distance_region_objective(
      numerator,element_count,true);
}

double StepDistanceObjectiveState::region_value(
    double region_numerator,
    int region_element_count,
    double region_unit_weight) const{
  (void)region_unit_weight;
  return step_distance_region_objective(
      region_numerator,region_element_count,true);
}

double StepDistanceObjectiveState::replaced_value(
    double old_region_numerator,
    int old_region_element_count,
    double old_region_unit_weight,
    double new_region_numerator,
    int new_region_element_count,
    double new_region_unit_weight) const{
  (void)old_region_unit_weight;
  (void)new_region_unit_weight;
  return step_distance_replaced_region_objective(
      numerator,old_region_numerator,new_region_numerator,
      element_count,old_region_element_count,new_region_element_count);
}

bool StepDistanceObjectiveState::accepts_replacement(
    double old_region_numerator,
    int old_region_element_count,
    double old_region_unit_weight,
    double new_region_numerator,
    int new_region_element_count,
    double new_region_unit_weight,
    double relative_improvement) const{
  const double global_old = value();
  const double global_new = replaced_value(
      old_region_numerator,old_region_element_count,old_region_unit_weight,
      new_region_numerator,new_region_element_count,new_region_unit_weight);
  return objective_strictly_improves(
      global_new,global_old,relative_improvement);
}

void StepDistanceObjectiveState::replace(
    double old_region_numerator,
    int old_region_element_count,
    double old_region_unit_weight,
    double new_region_numerator,
    int new_region_element_count,
    double new_region_unit_weight){
  (void)old_region_unit_weight;
  (void)new_region_unit_weight;
  numerator += new_region_numerator - old_region_numerator;
  element_count += new_region_element_count - old_region_element_count;
  METRIS_ENFORCE(element_count > 0);
  target_weight = element_count;
  best_objective = std::min(best_objective,value());
  if(best_objective_storage != nullptr){
    *best_objective_storage = std::min(*best_objective_storage,best_objective);
  }
}

double step_distance_region_objective(double elemental_sum,
                                      double element_count,
                                      bool cavity_target_average){
  if(cavity_target_average){
    METRIS_ENFORCE_MSG(element_count > 0,
                       "Nonpositive StepDistance element count");
    return elemental_sum/element_count;
  }
  return elemental_sum;
}

double step_distance_region_contribution(double element_value,
                                         double unit_weight,
                                         bool cavity_target_average){
  (void)unit_weight;
  (void)cavity_target_average;
  return element_value;
}

double step_distance_replaced_region_objective(
    double global_elemental_sum,
    double old_region_elemental_sum,
    double new_region_elemental_sum,
    double global_element_count,
    double old_region_element_count,
    double new_region_element_count){
  const double new_global_element_count =
      global_element_count - old_region_element_count
      + new_region_element_count;
  METRIS_ENFORCE_MSG(new_global_element_count > 0,
                     "Nonpositive replaced StepDistance element count");
  const double new_global_sum =
      global_elemental_sum - old_region_elemental_sum
      + new_region_elemental_sum;
  return new_global_sum/new_global_element_count;
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

  if(msh.param->step_distance_shape_volume
     && !shape_volume_matrix_is_numerically_admissible<tdim>(A)){
    return ftype(step_distance_shape_volume_rejection_quality);
  }

  ftype eigval[tdim];
  ftype eigvec[tdim*tdim];

  geteigsym<tdim,ftype>(Ap, eigval, eigvec);

  // ------------------------------------------------------------
  // Step Distance Shape Volume:
  //
  //   Ahat = A / det(A)^(1/tdim)
  //   v(A) = det(A) - 1/det(A)
  //   d^2  = ||log(Ahat)||_F^2 + v(A)^2/(4*tdim),
  //   psi  = (d^2 + eps^2)^(p/2) - eps^p, p >= 1.
  //
  // This is the squared distance associated with the injective Euclidean
  // embedding A -> (log(Ahat),v(A)/(2*sqrt(tdim))).  It agrees to second
  // order with ||log(A)||_F^2 at A = I and diverges algebraically when A
  // approaches the SPD boundary.
  // ------------------------------------------------------------
  if(msh.param->step_distance_shape_volume){
    ftype detA = ftype(1);
    ftype mean_log = ftype(0);
    ftype loglam[tdim];

    for(int i = 0; i < tdim; i++){
      METRIS_ENFORCE(eigval[i] > ftype(0));
      detA *= eigval[i];
      loglam[i] = log(eigval[i]);
      mean_log += loglam[i]/ftype(tdim);
    }

    ftype distance2 = ftype(0);
    for(int i = 0; i < tdim; i++){
      const ftype centered_log = loglam[i] - mean_log;
      distance2 += centered_log*centered_log;
    }

    const ftype volume_coordinate = detA - ftype(1)/detA;
    distance2 +=
        volume_coordinate*volume_coordinate/ftype(4*tdim);
    const double power = msh.param->objective_p;
    const double regularization
        = msh.param->step_distance_regularization;
    const ftype regularization_squared
        = ftype(regularization*regularization);
    return pow(distance2 + regularization_squared,power/2.0)
         - ftype(std::pow(regularization,power));
  }

  // ------------------------------------------------------------
  // Original StepDistance:
  //   psi = (||log(A)||_F^2 + eps^2)^(p/2) - eps^p.
  // A small smooth norm regularization is used so p < 2 remains
  // differentiable at A = I.
  // ------------------------------------------------------------
  ftype distance2 = ftype(0);

  for(int i = 0; i < tdim; i++){
    METRIS_ENFORCE(eigval[i] > ftype(0));
    ftype li = log(eigval[i]);
    distance2 += li*li;
  }

  const double power = msh.param->objective_p;
  const double regularization = msh.param->step_distance_regularization;
  const ftype regularization_squared
      = ftype(regularization*regularization);
  return pow(distance2 + regularization_squared,power/2.0)
       - ftype(std::pow(regularization,power));
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
void eval_step_distance_psi_gradient(const T* Jreg_T,
                                     const T* met,
                                     const T* gradN,
                                     double power,
                                     double regularization,
                                     T* psi,
                                     T* dpsi){

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

  const T regularization_squared = T(regularization*regularization);
  *psi = pow(distance2 + regularization_squared,power/2.0)
       - T(std::pow(regularization,power));

  if(dpsi == NULL) return;

  const T power_scale = T(power/2.0)
      *pow(distance2 + regularization_squared,power/2.0 - 1.0);

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

  // dpsi = 4 power_scale M u.
  for(int a = 0; a < gdim; a++){
    dpsi[a] = T(0);
    for(int b = 0; b < gdim; b++){
      const int imet = symidx_met<gdim,T>(a,b);
      dpsi[a] += met[imet]*u[b];
    }
    // The final factor is (p/2)(distance2 + eps^2)^(p/2 - 1)
    // times the p=2 derivative of distance2.
    dpsi[a] *= T(4)*power_scale;
  }
}

template<int gdim, int tdim, typename T>
void eval_shape_volume_psi_gradient(const T* Jreg_T,
                                    const T* met,
                                    const T* gradN,
                                    double power,
                                    double regularization,
                                    T* psi,
                                    T* dpsi){

  constexpr int nnmet_tdim = (tdim*(tdim+1))/2;

  // A = J^T M J = Jreg_T M Jreg_T^T.
  T A[tdim*tdim];
  for(int i = 0; i < tdim; i++){
    for(int j = 0; j < tdim; j++){
      A[i*tdim+j] = T(0);
      for(int a = 0; a < gdim; a++){
        for(int b = 0; b < gdim; b++){
          const int imet = symidx_met<gdim,T>(a,b);
          A[i*tdim+j] +=
              Jreg_T[i*gdim+a]*met[imet]*Jreg_T[j*gdim+b];
        }
      }
    }
  }

  if(!shape_volume_matrix_is_numerically_admissible<tdim>(A)){
    *psi = T(step_distance_shape_volume_rejection_quality);
    if(dpsi != NULL){
      for(int a = 0; a < gdim; a++) dpsi[a] = T(0);
    }
    return;
  }

  T Ap[nnmet_tdim];
  sym_full_to_packed_tdim<tdim,T>(A,Ap);

  T eigval[tdim];
  T eigvec[tdim*tdim];
  geteigsym<tdim,T>(Ap,eigval,eigvec);

  T detA = T(1);
  T mean_log = T(0);
  T loglam[tdim];
  for(int i = 0; i < tdim; i++){
    METRIS_ENFORCE(eigval[i] > T(0));
    detA *= eigval[i];
    loglam[i] = log(eigval[i]);
    mean_log += loglam[i]/T(tdim);
  }

  T distance2 = T(0);
  T centered_log[tdim];
  for(int i = 0; i < tdim; i++){
    centered_log[i] = loglam[i] - mean_log;
    distance2 += centered_log[i]*centered_log[i];
  }

  const T invdetA = T(1)/detA;
  const T volume_coordinate = detA - invdetA;
  distance2 += volume_coordinate*volume_coordinate/T(4*tdim);
  const T regularization_squared = T(regularization*regularization);
  *psi = pow(distance2 + regularization_squared,power/2.0)
       - T(std::pow(regularization,power));

  if(dpsi == NULL) return;

  // For the unpowered squared distance s(A) above,
  //
  //   d s / d A = 2 A^{-1} log(Ahat)
  //               + (det(A)^2-det(A)^(-2))/(2*tdim) A^{-1}.
  //
  // The power_scale below is d psi/ds.
  //
  // If one P1 vertex moves, dJ = dx gradN^T and therefore
  // grad_x psi = 2 M J (d psi/dA) gradN.
  const T volume_scale =
      (detA*detA - invdetA*invdetA)/T(2*tdim);
  const T power_scale = T(power/2.0)
      *pow(distance2 + regularization_squared,power/2.0 - 1.0);
  T gradA_eig[tdim];
  for(int i = 0; i < tdim; i++){
    gradA_eig[i] = power_scale
        *(T(2)*centered_log[i] + volume_scale)/eigval[i];
  }

  T gradA[tdim*tdim];
  spectral_matrix_from_eigs<tdim,T>(eigval,eigvec,gradA_eig,gradA);

  T v[tdim];
  for(int i = 0; i < tdim; i++){
    v[i] = T(0);
    for(int j = 0; j < tdim; j++){
      v[i] += gradA[i*tdim+j]*gradN[j];
    }
  }

  T Jv[gdim];
  for(int a = 0; a < gdim; a++){
    Jv[a] = T(0);
    for(int i = 0; i < tdim; i++){
      Jv[a] += Jreg_T[i*gdim+a]*v[i];
    }
  }

  for(int a = 0; a < gdim; a++){
    dpsi[a] = T(0);
    for(int b = 0; b < gdim; b++){
      const int imet = symidx_met<gdim,T>(a,b);
      dpsi[a] += met[imet]*Jv[b];
    }
    dpsi[a] *= T(2);
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

    ftype psi;
    if(msh.param->step_distance_shape_volume){
      eval_shape_volume_psi_gradient<gdim,tdim,ftype>(
          Jreg_T_f,met_f,NULL,
          msh.param->objective_p,
          msh.param->step_distance_regularization,
          &psi,NULL);
    }else{
      eval_step_distance_psi_gradient<gdim,tdim,ftype>(
          Jreg_T_f, met_f, NULL,
          msh.param->objective_p,
          msh.param->step_distance_regularization,
          &psi, NULL);
    }

    return psi;
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

  ftype psi;
  ftype dpsi[gdim];

  if(msh.param->step_distance_shape_volume){
    eval_shape_volume_psi_gradient<gdim,tdim,ftype>(
        Jreg_T_f,met_f,gradN_f,
        msh.param->objective_p,
        msh.param->step_distance_regularization,
        &psi,dpsi);
  }else{
    eval_step_distance_psi_gradient<gdim,tdim,ftype>(
        Jreg_T_f, met_f, gradN_f,
        msh.param->objective_p,
        msh.param->step_distance_regularization,
        &psi, dpsi);
  }

  for(int i = 0; i < gdim; i++){
    dquael[i] = dpsi[i];
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

    S psiS;
    S dpsiS[gdim];

    if(msh.param->step_distance_shape_volume){
      eval_shape_volume_psi_gradient<gdim,tdim,S>(
          Jreg_TS,metS,gradNS,
          msh.param->objective_p,
          msh.param->step_distance_regularization,
          &psiS,dpsiS);
    }else{
      eval_step_distance_psi_gradient<gdim,tdim,S>(
          Jreg_TS, metS, gradNS,
          msh.param->objective_p,
          msh.param->step_distance_regularization,
          &psiS, dpsiS);
    }

    for(int i = 0; i < gdim; i++){
      for(int j = i; j < gdim; j++){
        hquael[sym2idx(i,j)] = (ftype)dpsiS[i].deriv(j);
      }
    }
  }

  return psi;
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
