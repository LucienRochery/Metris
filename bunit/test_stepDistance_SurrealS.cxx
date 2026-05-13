//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_stepDistance_Surreal

#include "common_setup.hxx"

#include "libs/SANS/Surreal/SurrealS.h"
#include "linalg/eigen.hxx"

using namespace Metris;


// -----------------------------------------------------------------------------
// Helper for SurrealS/double value extraction
// -----------------------------------------------------------------------------
double get_value(const double& x)
{
  return x;
}

template<int nvar>
double get_value(const SANS::SurrealS<nvar,double>& x)
{
  return x.value();
}

// -----------------------------------------------------------------------------
// Regular reference simplex inverse edge matrix.
// E0 = [X1hat-X0hat, X2hat-X0hat, ...]
// -----------------------------------------------------------------------------
template<int dim, typename T>
void get_regular_E0inv(T* E0inv)
{
  for(int i = 0; i < dim*dim; i++) E0inv[i] = T(0.0);

  if constexpr(dim == 2){
    const T invsqrt3 = T(1.0 / std::sqrt(3.0));

    E0inv[0] = T(1.0);
    E0inv[1] = -invsqrt3;

    E0inv[2] = T(0.0);
    E0inv[3] = T(2.0) * invsqrt3;
  }

  if constexpr(dim == 3){
    const T invsqrt3 = T(1.0 / std::sqrt(3.0));
    const T invsqrt6 = T(1.0 / std::sqrt(6.0));
    const T sqrt32   = T(std::sqrt(3.0 / 2.0));

    E0inv[0] = T(1.0);
    E0inv[1] = -invsqrt3;
    E0inv[2] = -invsqrt6;

    E0inv[3] = T(0.0);
    E0inv[4] = T(2.0) * invsqrt3;
    E0inv[5] = -invsqrt6;

    E0inv[6] = T(0.0);
    E0inv[7] = T(0.0);
    E0inv[8] = sqrt32;
  }
}

// call analytical metric field
template<int dim>
void eval_test_anamet(const double* x, double scale, int idif1,
                      double* met, double* dmet)
{
  if constexpr(dim == 2){
    anamet2D_7(nullptr, x, scale, idif1, met, dmet); // corner singularity
  }else{
    anamet3D_12(nullptr, x, scale, idif1, met, dmet); // Yano BL function
  }
}

// -----------------------------------------------------------------------------
// Convert symmetric packed metric to full matrix.
//
// Current assumed packing:
// 2D: [m00, m01, m11]
// 3D: [m00, m01, m11, m02, m12, m22]
// This matches the examples you showed.
// -----------------------------------------------------------------------------
template<int dim, typename T>
void sym_packed_to_full(const T* met, T* M)
{
  for(int i = 0; i < dim*dim; i++) M[i] = T(0.0);

  if constexpr(dim == 2){
    M[0] = met[0];
    M[1] = met[1];
    M[2] = met[1];
    M[3] = met[2];
  }

  if constexpr(dim == 3){
    M[0] = met[0];
    M[1] = met[1];
    M[2] = met[3];

    M[3] = met[1];
    M[4] = met[2];
    M[5] = met[4];

    M[6] = met[3];
    M[7] = met[4];
    M[8] = met[5];
  }
}

template<int dim, typename T>
void sym_full_to_packed(const T* A, T* Ap)
{
  if constexpr(dim == 2){
    // packed: [a00, a01, a11]
    Ap[0] = A[0]; // a00
    Ap[1] = A[1]; // a01
    Ap[2] = A[3]; // a11
  }

  if constexpr(dim == 3){
    // packed: [a00, a01, a11, a02, a12, a22]
    Ap[0] = A[0]; // a00
    Ap[1] = A[1]; // a01
    Ap[2] = A[4]; // a11
    Ap[3] = A[2]; // a02
    Ap[4] = A[5]; // a12
    Ap[5] = A[8]; // a22
  }
}

// -----------------------------------------------------------------------------
// Lift double metric + spatial metric derivatives into SurrealS.
//
// metS.deriv(ivar) = sum_k dmet/dx_k * xS[k].deriv(ivar)
//
// -----------------------------------------------------------------------------
template<int dim, int nvar>
void lift_metric_to_surreal(const SANS::SurrealS<nvar,double>* xS,
                            const double* met,
                            const double* dmet,
                            SANS::SurrealS<nvar,double>* metS)
{
  constexpr int nnmet = dim*(dim+1)/2;

  for(int im = 0; im < nnmet; im++){
    metS[im].value() = met[im];

    for(int ivar = 0; ivar < nvar; ivar++){
      metS[im].deriv(ivar) = 0.0;

      for(int k = 0; k < dim; k++){
        metS[im].deriv(ivar) += dmet[k*nnmet + im] * xS[k].deriv(ivar);
      }
    }
  }
}

// -----------------------------------------------------------------------------
// Compute phi = ||log(A)||_F^2 at the barycenter.
//
// X is an array of dim+1 vertices, each in R^dim.
// X[active_vertex] may contain SurrealS entries.
// Metric is evaluated at x = average vertices, with dmet lifted into Surreal.
// -----------------------------------------------------------------------------
template<int dim, typename T>
T compute_phi(const T X[dim+1][dim])
{
  constexpr int nnmet = dim*(dim+1)/2;
  constexpr double scale = 1.0;

  // E_K = [X1-X0, X2-X0, ...], shape dim x dim.
  T EK[dim*dim];

  for(int col = 0; col < dim; col++){
    const int ivert = col + 1;
    for(int row = 0; row < dim; row++){
      EK[dim*row + col] = X[ivert][row] - X[0][row];
    }
  }

  // J = E_K E_0^{-1}.
  T E0inv[dim*dim];
  get_regular_E0inv<dim,T>(E0inv);

  T J[dim*dim];
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      J[dim*i + j] = T(0.0);
      for(int k = 0; k < dim; k++){
        J[dim*i + j] += EK[dim*i + k] * E0inv[dim*k + j];
      }
    }
  }

  // Physical point at barycenter.
  T x[dim];
  for(int k = 0; k < dim; k++){
    x[k] = T(0.0);
    for(int a = 0; a < dim+1; a++){
      x[k] += X[a][k] / T(dim + 1);
    }
  }

  // Evaluate analytical metric with double coordinates.
  double xd[dim];
  for(int k = 0; k < dim; k++){
    xd[k] = get_value(x[k]);
  }

  double met[nnmet];
  double dmet[dim*nnmet];

  eval_test_anamet<dim>(xd, scale, 1, met, dmet);

  // Build metric in type T.
  T metT[nnmet];

  if constexpr(std::is_same<T,double>::value){
    for(int im = 0; im < nnmet; im++){
      metT[im] = met[im];
    }
  }else{
    // This branch assumes T = SurrealS<dim,double>.
    lift_metric_to_surreal<dim,dim>(x, met, dmet, metT);
  }

  // Full matrix M.
  T M[dim*dim];
  sym_packed_to_full<dim,T>(metT, M);

  // A = J^T M J.
  T A[dim*dim];
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      A[dim*i + j] = T(0.0);
      for(int a = 0; a < dim; a++){
        for(int b = 0; b < dim; b++){
          A[dim*i + j] += J[dim*a + i] * M[dim*a + b] * J[dim*b + j];
        }
      }
    }
  }

  T Ap[dim*(dim+1)/2];
  sym_full_to_packed<dim,T>(A, Ap);

  T eigval[dim];
  T eigvec[dim*dim];

  geteigsym<dim,T>(Ap, eigval, eigvec);

  T phi = T(0.0);
  for(int i = 0; i < dim; i++){
    T li = log(eigval[i]);
    phi += li*li;
  }

  return phi;
}

// -----------------------------------------------------------------------------
// Helper to compute determinant
// -----------------------------------------------------------------------------
template<int dim, typename T>
T det_full(const T* A)
{
  if constexpr(dim == 2){
    return A[0]*A[3] - A[1]*A[2];
  }

  if constexpr(dim == 3){
    return
        A[0]*(A[4]*A[8] - A[5]*A[7])
      - A[1]*(A[3]*A[8] - A[5]*A[6])
      + A[2]*(A[3]*A[7] - A[4]*A[6]);
  }
}

// -----------------------------------------------------------------------------
// Compute theta = sqrt( det(J^T J) ) * sqrt( det(M) )
//
// X is an array of dim+1 vertices, each in R^dim.
// X[active_vertex] may contain SurrealS entries.
// Metric is evaluated at x = average vertices, with dmet lifted into Surreal.
// -----------------------------------------------------------------------------
template<int dim, typename T>
T compute_theta(const T X[dim+1][dim])
{
  constexpr int nnmet = dim*(dim+1)/2;
  constexpr double scale = 1.0;

  // E_K = [X1-X0, X2-X0, ...], shape dim x dim.
  T EK[dim*dim];

  for(int col = 0; col < dim; col++){
    const int ivert = col + 1;
    for(int row = 0; row < dim; row++){
      EK[dim*row + col] = X[ivert][row] - X[0][row];
    }
  }

  // J = E_K E_0^{-1}.
  T E0inv[dim*dim];
  get_regular_E0inv<dim,T>(E0inv);

  T J[dim*dim];
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      J[dim*i + j] = T(0.0);
      for(int k = 0; k < dim; k++){
        J[dim*i + j] += EK[dim*i + k] * E0inv[dim*k + j];
      }
    }
  }

  // Physical point at barycenter.
  T x[dim];
  for(int k = 0; k < dim; k++){
    x[k] = T(0.0);
    for(int a = 0; a < dim+1; a++){
      x[k] += X[a][k] / T(dim + 1);
    }
  }

  // Evaluate analytical metric with double coordinates.
  double xd[dim];
  for(int k = 0; k < dim; k++){
    xd[k] = get_value(x[k]);
  }

  double met[nnmet];
  double dmet[dim*nnmet];

  eval_test_anamet<dim>(xd, scale, 1, met, dmet);

  // Build metric in type T.
  T metT[nnmet];

  if constexpr(std::is_same<T,double>::value){
    for(int im = 0; im < nnmet; im++){
      metT[im] = met[im];
    }
  }else{
    // This branch assumes T = SurrealS<dim,double>.
    lift_metric_to_surreal<dim,dim>(x, met, dmet, metT);
  }

  // Full matrix M.
  T M[dim*dim];
  sym_packed_to_full<dim,T>(metT, M);

  // Ggeo = J^T J.
  T Ggeo[dim*dim];
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      Ggeo[dim*i + j] = T(0.0);
      for(int a = 0; a < dim; a++){
        Ggeo[dim*i + j] += J[dim*a + i] * J[dim*a + j];
      }
    }
  }

  // theta = sqrt(det(J^T J)) * sqrt(det(M*)).
  T theta = sqrt(det_full<dim,T>(Ggeo)) * sqrt(det_full<dim,T>(M));

  return theta;
}

// -----------------------------------------------------------------------------
// Copy double vertices into T vertices.
// -----------------------------------------------------------------------------
template<int dim, typename T>
void copy_vertices_to_T(const double Xd[dim+1][dim], T X[dim+1][dim])
{
  for(int a = 0; a < dim+1; a++){
    for(int k = 0; k < dim; k++){
      X[a][k] = T(Xd[a][k]);
    }
  }
}

// -----------------------------------------------------------------------------
// Run one AD-vs-FD test on phi only.
// active_vertex can be 0, ..., dim.
// -----------------------------------------------------------------------------
template<int dim>
void run_stepdist_surreal_test(
    const double Xd[dim+1][dim],
    int active_vertex)
{
  using S = SANS::SurrealS<dim,double>;

  METRIS_ENFORCE(active_vertex >= 0 && active_vertex < dim+1);

  S Xs[dim+1][dim];
  copy_vertices_to_T<dim,S>(Xd, Xs);

  // Initialize values and clear derivatives.
  for(int a = 0; a < dim+1; a++){
    for(int k = 0; k < dim; k++){
      Xs[a][k].value() = Xd[a][k];

      for(int j = 0; j < dim; j++){
        Xs[a][k].deriv(j) = 0.0;
      }
    }
  }

  // Seed active vertex derivatives.
  for(int k = 0; k < dim; k++){
    Xs[active_vertex][k].deriv(k) = 1.0;
  }

  S phi_ad = compute_phi<dim,S>(Xs);

  const double eps = 1.0e-6;

  fmt::print("\n-- Testing {}D, active vertex {}\n", dim, active_vertex);
  fmt::print("phi = {:.16e}\n", phi_ad.value());

  BOOST_CHECK(std::isfinite(phi_ad.value()));

  for(int k = 0; k < dim; k++){
    double Xp[dim+1][dim];
    double Xm[dim+1][dim];

    for(int a = 0; a < dim+1; a++){
      for(int j = 0; j < dim; j++){
        Xp[a][j] = Xd[a][j];
        Xm[a][j] = Xd[a][j];
      }
    }

    Xp[active_vertex][k] += eps;
    Xm[active_vertex][k] -= eps;

    const double phi_p = compute_phi<dim,double>(Xp);
    const double phi_m = compute_phi<dim,double>(Xm);

    const double fd = (phi_p - phi_m) / (2.0 * eps);
    const double ad = phi_ad.deriv(k);

    fmt::print("dphi/dX[{}][{}] AD = {:.16e}, FD = {:.16e}, err = {:.3e}\n",
               active_vertex, k, ad, fd, std::abs(ad - fd));

    BOOST_CHECK(std::isfinite(ad));
    BOOST_CHECK(std::isfinite(fd));
    BOOST_CHECK_CLOSE(ad, fd, 1.0e-5);
  }
}

// -----------------------------------------------------------------------------
// Run one AD-vs-FD test on theta only.
// active_vertex can be 0, ..., dim.
// -----------------------------------------------------------------------------
template<int dim>
void run_geometry_surreal_test(
    const double Xd[dim+1][dim],
    int active_vertex)
{
  using S = SANS::SurrealS<dim,double>;

  METRIS_ENFORCE(active_vertex >= 0 && active_vertex < dim+1);

  S Xs[dim+1][dim];
  copy_vertices_to_T<dim,S>(Xd, Xs);

  // Initialize values and clear derivatives.
  for(int a = 0; a < dim+1; a++){
    for(int k = 0; k < dim; k++){
      Xs[a][k].value() = Xd[a][k];

      for(int j = 0; j < dim; j++){
        Xs[a][k].deriv(j) = 0.0;
      }
    }
  }

  // Seed active vertex derivatives.
  for(int k = 0; k < dim; k++){
    Xs[active_vertex][k].deriv(k) = 1.0;
  }

  S theta_ad = compute_theta<dim,S>(Xs);

  const double eps = 1.0e-6;

  fmt::print("\n-- Testing {}D, active vertex {}\n", dim, active_vertex);
  fmt::print("theta = {:.16e}\n", theta_ad.value());

  BOOST_CHECK(std::isfinite(theta_ad.value()));

  for(int k = 0; k < dim; k++){
    double Xp[dim+1][dim];
    double Xm[dim+1][dim];

    for(int a = 0; a < dim+1; a++){
      for(int j = 0; j < dim; j++){
        Xp[a][j] = Xd[a][j];
        Xm[a][j] = Xd[a][j];
      }
    }

    Xp[active_vertex][k] += eps;
    Xm[active_vertex][k] -= eps;

    const double theta_p = compute_theta<dim,double>(Xp);
    const double theta_m = compute_theta<dim,double>(Xm);

    const double fd = (theta_p - theta_m) / (2.0 * eps);
    const double ad = theta_ad.deriv(k);

    fmt::print("dtheta/dX[{}][{}] AD = {:.16e}, FD = {:.16e}, err = {:.3e}\n",
               active_vertex, k, ad, fd, std::abs(ad - fd));

    BOOST_CHECK(std::isfinite(ad));
    BOOST_CHECK(std::isfinite(fd));
    BOOST_CHECK_CLOSE(ad, fd, 1.0e-5);
  }
}

// -----------------------------------------------------------------------------
// Inverse of full dim x dim matrix.
// -----------------------------------------------------------------------------
template<int dim>
void inv_full(const double* A, double* Ai)
{
  const double detA = det_full<dim>(A);
  METRIS_ENFORCE(std::abs(detA) > 1.0e-14);

  if constexpr(dim == 2){
    Ai[0] =  A[3]/detA;
    Ai[1] = -A[1]/detA;
    Ai[2] = -A[2]/detA;
    Ai[3] =  A[0]/detA;
  }

  if constexpr(dim == 3){
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


// -----------------------------------------------------------------------------
// C = A * B, full row-major.
// -----------------------------------------------------------------------------
template<int dim>
void matmul_full(const double* A, const double* B, double* C)
{
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      C[dim*i+j] = 0.0;
      for(int k = 0; k < dim; k++){
        C[dim*i+j] += A[dim*i+k]*B[dim*k+j];
      }
    }
  }
}


// -----------------------------------------------------------------------------
// y = A x, full row-major.
// -----------------------------------------------------------------------------
template<int dim>
void matvec_full(const double* A, const double* x, double* y)
{
  for(int i = 0; i < dim; i++){
    y[i] = 0.0;
    for(int j = 0; j < dim; j++){
      y[i] += A[dim*i+j]*x[j];
    }
  }
}

// -----------------------------------------------------------------------------
// Build spectral matrix F(A) = sum_i f(lambda_i) q_i q_i^T.
// Assumes eigvec stores eigenvectors by rows: q_i[k] = eigvec[dim*i+k].
// This matches the swapping logic in dsyevq.
// -----------------------------------------------------------------------------
template<int dim>
void spectral_matrix_from_eigs(const double* eigval,
                               const double* eigvec,
                               const double* fval,
                               double* F)
{
  for(int i = 0; i < dim*dim; i++) F[i] = 0.0;

  for(int a = 0; a < dim; a++){
    for(int i = 0; i < dim; i++){
      for(int j = 0; j < dim; j++){
        F[dim*i+j] += fval[a] * eigvec[dim*a+i] * eigvec[dim*a+j];
      }
    }
  }
}

// -----------------------------------------------------------------------------
// Gradient of regular reference shape function N_q wrt reference coordinates y.
// -----------------------------------------------------------------------------
template<int dim>
void get_regular_gradN(int q, double* gradN)
{
  double E0inv[dim*dim];
  get_regular_E0inv<dim,double>(E0inv);

  for(int i = 0; i < dim; i++) gradN[i] = 0.0;

  if(q == 0){
    for(int a = 0; a < dim; a++){
      for(int i = 0; i < dim; i++){
        gradN[i] -= E0inv[dim*a+i]; // row a
      }
    }
  }else{
    const int row = q - 1;
    for(int i = 0; i < dim; i++){
      gradN[i] = E0inv[dim*row+i];
    }
  }
}

// -----------------------------------------------------------------------------
// Convert packed dmet component dM/dx_k to full matrix.
// -----------------------------------------------------------------------------
template<int dim>
void dmet_packed_component_to_full(const double* dmet, int k, double* dMdxk)
{
  constexpr int nnmet = dim*(dim+1)/2;
  sym_packed_to_full<dim,double>(&dmet[k*nnmet], dMdxk);
}


// -----------------------------------------------------------------------------
// Compute analytic gradient of one barycenter contribution
//
// I_r = phi * theta
//
// with g = M*.
// Weight w_r omitted here
// -----------------------------------------------------------------------------
template<int dim>
void compute_I_analytic_grad(const double X[dim+1][dim],
                                         int active_vertex,
                                         double* grad)
{
  constexpr int nnmet = dim*(dim+1)/2;
  constexpr double scale = 1.0;

  METRIS_ENFORCE(active_vertex >= 0 && active_vertex < dim+1);

  // ---- EK
  double EK[dim*dim];
  for(int col = 0; col < dim; col++){
    const int ivert = col + 1;
    for(int row = 0; row < dim; row++){
      EK[dim*row+col] = X[ivert][row] - X[0][row];
    }
  }

  // ---- J = EK E0inv
  double E0inv[dim*dim];
  get_regular_E0inv<dim,double>(E0inv);

  double J[dim*dim];
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      J[dim*i+j] = 0.0;
      for(int k = 0; k < dim; k++){
        J[dim*i+j] += EK[dim*i+k]*E0inv[dim*k+j];
      }
    }
  }

  // ---- barycenter x and N_q
  double x[dim];
  for(int k = 0; k < dim; k++){
    x[k] = 0.0;
    for(int a = 0; a < dim+1; a++){
      x[k] += X[a][k] / double(dim+1);
    }
  }

  const double Nq = 1.0 / double(dim+1);

  double gradN[dim];
  get_regular_gradN<dim>(active_vertex, gradN);

  // ---- metric M and dM
  double met[nnmet];
  double dmet[dim*nnmet];

  eval_test_anamet<dim>(x, scale, 1, met, dmet);

  double M[dim*dim];
  sym_packed_to_full<dim,double>(met, M);

  // ---- A = J^T M J
  double A[dim*dim];
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      A[dim*i+j] = 0.0;
      for(int a = 0; a < dim; a++){
        for(int b = 0; b < dim; b++){
          A[dim*i+j] += J[dim*a+i]*M[dim*a+b]*J[dim*b+j];
        }
      }
    }
  }

  // ---- eig(A), log(A), A^{-1}log(A)
  double Ap[nnmet];
  sym_full_to_packed<dim,double>(A, Ap);

  double eigval[dim];
  double eigvec[dim*dim];
  geteigsym<dim,double>(Ap, eigval, eigvec);

  double loglam[dim];
  double invloglam[dim];

  double phi = 0.0;
  for(int i = 0; i < dim; i++){
    loglam[i] = std::log(eigval[i]);
    invloglam[i] = loglam[i] / eigval[i];
    phi += loglam[i]*loglam[i];
  }

  double logA[dim*dim];
  double AinvlogA[dim*dim];

  spectral_matrix_from_eigs<dim>(eigval, eigvec, loglam, logA);
  spectral_matrix_from_eigs<dim>(eigval, eigvec, invloglam, AinvlogA);

  // ---- C = J^T J, theta = sqrt(det(C)) sqrt(det(M))
  double C[dim*dim];
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      C[dim*i+j] = 0.0;
      for(int a = 0; a < dim; a++){
        C[dim*i+j] += J[dim*a+i]*J[dim*a+j];
      }
    }
  }

  const double detC = det_full<dim>(C);
  const double detM = det_full<dim>(M);
  METRIS_ENFORCE(detC > 0.0);
  METRIS_ENFORCE(detM > 0.0);

  const double theta = std::sqrt(detC) * std::sqrt(detM);

  // ---- m^phi_i = 2 tr(A^{-1}log(A) J^T dM/dx_i J)
  double mphi[dim];

  for(int k = 0; k < dim; k++){
    double dMdxk[dim*dim];
    dmet_packed_component_to_full<dim>(dmet, k, dMdxk);

    double B[dim*dim]; // B = J^T dM J
    for(int i = 0; i < dim; i++){
      for(int j = 0; j < dim; j++){
        B[dim*i+j] = 0.0;
        for(int a = 0; a < dim; a++){
          for(int b = 0; b < dim; b++){
            B[dim*i+j] += J[dim*a+i]*dMdxk[dim*a+b]*J[dim*b+j];
          }
        }
      }
    }

    mphi[k] = 0.0;
    for(int i = 0; i < dim; i++){
      for(int j = 0; j < dim; j++){
        mphi[k] += AinvlogA[dim*i+j] * B[dim*j+i];
      }
    }
    mphi[k] *= 2.0;
  }

  // ---- m^theta_i = tr(M^{-1} dM/dx_i), because g = M*
  double Minv[dim*dim];
  inv_full<dim>(M, Minv);

  double mtheta[dim];

  for(int k = 0; k < dim; k++){
    double dMdxk[dim*dim];
    dmet_packed_component_to_full<dim>(dmet, k, dMdxk);

    mtheta[k] = 0.0;
    for(int i = 0; i < dim; i++){
      for(int j = 0; j < dim; j++){
        mtheta[k] += Minv[dim*i+j] * dMdxk[dim*j+i];
      }
    }
  }

  // ---- term_phi = 4 M J A^{-1}log(A) gradN + Nq mphi
  double MJ[dim*dim];
  matmul_full<dim>(M, J, MJ);

  double MJAinvlog[dim*dim];
  matmul_full<dim>(MJ, AinvlogA, MJAinvlog);

  double term_phi[dim];
  matvec_full<dim>(MJAinvlog, gradN, term_phi);

  for(int k = 0; k < dim; k++){
    term_phi[k] = 4.0*term_phi[k] + Nq*mphi[k];
  }

  // ---- term_theta = J (J^T J)^{-1} gradN + Nq/2 mtheta
  double Cinv[dim*dim];
  inv_full<dim>(C, Cinv);

  double JCinv[dim*dim];
  matmul_full<dim>(J, Cinv, JCinv);

  double term_theta[dim];
  matvec_full<dim>(JCinv, gradN, term_theta);

  for(int k = 0; k < dim; k++){
    term_theta[k] += 0.5*Nq*mtheta[k];
  }

  // ---- grad I = theta term_phi + phi theta term_theta
  for(int k = 0; k < dim; k++){
    grad[k] = theta*term_phi[k] + phi*theta*term_theta[k];
  }
}


// -----------------------------------------------------------------------------
// Run full comparison test:
//
//  - AD
//  - finite difference
//  - analytic
//
// for:
//  - phi
//  - theta
//  - I = phi * theta
//
// active_vertex can be 0, ..., dim.
// -----------------------------------------------------------------------------
template<int dim>
void run_full_stepdist_test(
    const double Xd[dim+1][dim],
    int active_vertex)
{
  using S = SANS::SurrealS<dim,double>;

  METRIS_ENFORCE(active_vertex >= 0 && active_vertex < dim+1);

  // ---------------------------------------------------------------------------
  // Build Surreal vertices
  // ---------------------------------------------------------------------------
  S Xs[dim+1][dim];
  copy_vertices_to_T<dim,S>(Xd, Xs);

  for(int a = 0; a < dim+1; a++){
    for(int k = 0; k < dim; k++){

      Xs[a][k].value() = Xd[a][k];

      for(int j = 0; j < dim; j++){
        Xs[a][k].deriv(j) = 0.0;
      }
    }
  }

  // Seed active vertex
  for(int k = 0; k < dim; k++){
    Xs[active_vertex][k].deriv(k) = 1.0;
  }

  // ---------------------------------------------------------------------------
  // AD quantities
  // ---------------------------------------------------------------------------
  S phi_ad   = compute_phi<dim,S>(Xs);
  S theta_ad = compute_theta<dim,S>(Xs);
  S I_ad     = phi_ad * theta_ad;

  // ---------------------------------------------------------------------------
  // Analytic gradient for I
  // ---------------------------------------------------------------------------
  double gradI_an[dim];
  compute_I_analytic_grad<dim>(Xd, active_vertex, gradI_an);

  // ---------------------------------------------------------------------------
  // Finite difference setup
  // ---------------------------------------------------------------------------
  const double eps = 1.0e-6;

  fmt::print("\n=====================================================\n");
  fmt::print("FULL STEP-DIST TEST {}D active_vertex {}\n",
             dim, active_vertex);
  fmt::print("=====================================================\n");

  fmt::print("phi   = {:.16e}\n", phi_ad.value());
  fmt::print("theta = {:.16e}\n", theta_ad.value());
  fmt::print("I     = {:.16e}\n", I_ad.value());

  BOOST_CHECK(std::isfinite(phi_ad.value()));
  BOOST_CHECK(std::isfinite(theta_ad.value()));
  BOOST_CHECK(std::isfinite(I_ad.value()));

  // ---------------------------------------------------------------------------
  // Loop over active coordinate direction
  // ---------------------------------------------------------------------------
  for(int k = 0; k < dim; k++){

    // -------------------------------------------------------------------------
    // Perturbed geometries
    // -------------------------------------------------------------------------
    double Xp[dim+1][dim];
    double Xm[dim+1][dim];

    for(int a = 0; a < dim+1; a++){
      for(int j = 0; j < dim; j++){
        Xp[a][j] = Xd[a][j];
        Xm[a][j] = Xd[a][j];
      }
    }

    Xp[active_vertex][k] += eps;
    Xm[active_vertex][k] -= eps;

    // -------------------------------------------------------------------------
    // FD quantities
    // -------------------------------------------------------------------------
    const double phi_p   = compute_phi<dim,double>(Xp);
    const double phi_m   = compute_phi<dim,double>(Xm);

    const double theta_p = compute_theta<dim,double>(Xp);
    const double theta_m = compute_theta<dim,double>(Xm);

    const double I_p = phi_p * theta_p;
    const double I_m = phi_m * theta_m;

    const double dphi_fd   = (phi_p   - phi_m  ) / (2.0*eps);
    const double dtheta_fd = (theta_p - theta_m) / (2.0*eps);
    const double dI_fd     = (I_p     - I_m    ) / (2.0*eps);

    // -------------------------------------------------------------------------
    // AD quantities
    // -------------------------------------------------------------------------
    const double dphi_ad   = phi_ad.deriv(k);
    const double dtheta_ad = theta_ad.deriv(k);
    const double dI_ad     = I_ad.deriv(k);

    // -------------------------------------------------------------------------
    // Analytic quantities
    // -------------------------------------------------------------------------
    const double dI_an = gradI_an[k];

    // -------------------------------------------------------------------------
    // Print
    // -------------------------------------------------------------------------
    fmt::print("\n--------------------------------------------------\n");
    fmt::print("Coordinate direction k = {}\n", k);
    fmt::print("--------------------------------------------------\n");

    fmt::print("dphi/dX:\n");
    fmt::print("  AD = {:.16e}\n", dphi_ad);
    fmt::print("  FD = {:.16e}\n", dphi_fd);
    fmt::print("  |AD-FD| = {:.3e}\n",
               std::abs(dphi_ad - dphi_fd));

    fmt::print("\ndtheta/dX:\n");
    fmt::print("  AD = {:.16e}\n", dtheta_ad);
    fmt::print("  FD = {:.16e}\n", dtheta_fd);
    fmt::print("  |AD-FD| = {:.3e}\n",
               std::abs(dtheta_ad - dtheta_fd));

    fmt::print("\ndI/dX:\n");
    fmt::print("  AD = {:.16e}\n", dI_ad);
    fmt::print("  FD = {:.16e}\n", dI_fd);
    fmt::print("  AN = {:.16e}\n", dI_an);

    fmt::print("  |AD-FD| = {:.3e}\n",
               std::abs(dI_ad - dI_fd));

    fmt::print("  |AD-AN| = {:.3e}\n",
               std::abs(dI_ad - dI_an));

    fmt::print("  |FD-AN| = {:.3e}\n",
               std::abs(dI_fd - dI_an));

    // -------------------------------------------------------------------------
    // Checks
    // -------------------------------------------------------------------------
    BOOST_CHECK(std::isfinite(dphi_ad));
    BOOST_CHECK(std::isfinite(dtheta_ad));
    BOOST_CHECK(std::isfinite(dI_ad));

    BOOST_CHECK(std::isfinite(dphi_fd));
    BOOST_CHECK(std::isfinite(dtheta_fd));
    BOOST_CHECK(std::isfinite(dI_fd));

    BOOST_CHECK(std::isfinite(dI_an));

    BOOST_CHECK_CLOSE(dphi_ad, dphi_fd, 1.0e-5);
    BOOST_CHECK_CLOSE(dtheta_ad, dtheta_fd, 1.0e-5);

    BOOST_CHECK_CLOSE(dI_ad, dI_fd, 1.0e-5);
    BOOST_CHECK_CLOSE(dI_ad, dI_an, 1.0e-5);
  }
}

//==============================================


BOOST_AUTO_TEST_CASE(test_stepdist_surreal_phi_gradient)
{
  {
    // Keep triangle away from origin because anamet2D_7 has a corner singularity
    double X[3][2] = {
      {0.50, 0.40},
      {1.60, 0.55},
      {0.95, 1.35}
    };

    for(int active_vertex = 0; active_vertex < 3; active_vertex++){
      run_stepdist_surreal_test<2>(X, active_vertex);
    }
  }

  {
    // 3D BL metric varies strongly with x[0], so avoid extreme coordinates.
    double X[4][3] = {
      {0.10, 0.00,  0.00},
      {1.10, 0.10, -0.05},
      {0.45, 0.92,  0.05},
      {0.55, 0.25,  0.85}
    };

    for(int active_vertex = 0; active_vertex < 4; active_vertex++){
      run_stepdist_surreal_test<3>(X, active_vertex);
    }
  }
}

BOOST_AUTO_TEST_CASE(test_surreal_theta_gradient)
{
  {
    // Keep triangle away from origin because anamet2D_7 has a corner singularity
    double X[3][2] = {
      {0.50, 0.40},
      {1.60, 0.55},
      {0.95, 1.35}
    };

    for(int active_vertex = 0; active_vertex < 3; active_vertex++){
      run_geometry_surreal_test<2>(X, active_vertex);
    }
  }

  {
    // 3D BL metric varies strongly with x[0], so avoid extreme coordinates.
    double X[4][3] = {
      {0.10, 0.00,  0.00},
      {1.10, 0.10, -0.05},
      {0.45, 0.92,  0.05},
      {0.55, 0.25,  0.85}
    };

    for(int active_vertex = 0; active_vertex < 4; active_vertex++){
      run_geometry_surreal_test<3>(X, active_vertex);
    }
  }
}

BOOST_AUTO_TEST_CASE(test_full_stepdist_gradient)
{
  {
    double X[3][2] = {
      {0.50, 0.40},
      {1.60, 0.55},
      {0.95, 1.35}
    };

    for(int active_vertex = 0; active_vertex < 3; active_vertex++){
      run_full_stepdist_test<2>(X, active_vertex);
    }
  }

  {
    double X[4][3] = {
      {0.10, 0.00,  0.00},
      {1.10, 0.10, -0.05},
      {0.45, 0.92,  0.05},
      {0.55, 0.25,  0.85}
    };

    for(int active_vertex = 0; active_vertex < 4; active_vertex++){
      run_full_stepdist_test<3>(X, active_vertex);
    }
  }
}
