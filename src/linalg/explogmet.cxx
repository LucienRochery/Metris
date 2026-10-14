//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "../linalg/explogmet.hxx"
#include "../linalg/eigen.hxx"
#include "../aux_exceptions.hxx"

#include "../low_geo/misc.hxx"
#include "../linalg/matprods.hxx"
#include "../linalg/utils.hxx"

#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"

#include <cmath>
#include <algorithm>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>
#include "fmt/format.h"

namespace Metris{

template<int ndim>
void getexpmet_frechet(const double *logm, const double *dlogm, int nderiv,
                      double *met, double *dmet){
  static_assert(ndim == 2 || ndim == 3);
  constexpr int nm = ndim*(ndim+1)/2;
  METRIS_ENFORCE(logm != nullptr && met != nullptr && nderiv >= 0);
  METRIS_ENFORCE(nderiv == 0 || (dlogm != nullptr && dmet != nullptr));
  using Matrix = Eigen::Matrix<double,ndim,ndim>;
  Matrix L;
  for(int i = 0; i < ndim; i++) for(int j = 0; j < ndim; j++)
    L(i,j) = logm[sym2idx(i,j)];
  METRIS_ENFORCE_MSG(L.allFinite(), "Nonfinite log metric in matrix exponential");
  const Eigen::SelfAdjointEigenSolver<Matrix> eigen(L);
  METRIS_ENFORCE(eigen.info() == Eigen::Success);
  const Matrix Q = eigen.eigenvectors();
  const auto lambda = eigen.eigenvalues();
  const auto expLambda = lambda.array().exp().eval();
  METRIS_ENFORCE_MSG(expLambda.allFinite() && (expLambda > 0.).all(),
                    "Metric exponential is outside the finite SPD range");
  const Matrix M = Q*expLambda.matrix().asDiagonal()*Q.transpose();
  for(int i = 0; i < ndim; i++) for(int j = 0; j <= i; j++)
    met[sym2idx(i,j)] = M(i,j);
  if(nderiv == 0) return;

  // Divided differences use their continuous limit at repeated eigenvalues.
  Matrix dividedDifference;
  for(int i = 0; i < ndim; i++) for(int j = 0; j < ndim; j++){
    const double gap = std::abs(lambda[i]-lambda[j]);
    dividedDifference(i,j) = std::exp(std::max(lambda[i],lambda[j]))
        *(gap == 0. ? 1. : -std::expm1(-gap)/gap);
  }
  for(int r = 0; r < nderiv; r++){
    Matrix derivative;
    for(int i = 0; i < ndim; i++) for(int j = 0; j < ndim; j++)
      derivative(i,j) = dlogm[r*nm+sym2idx(i,j)];
    METRIS_ENFORCE(derivative.allFinite());
    const Matrix rotated = Q.transpose()*derivative*Q;
    const Matrix DM = Q*dividedDifference.cwiseProduct(rotated)*Q.transpose();
    METRIS_ENFORCE(DM.allFinite());
    for(int i = 0; i < ndim; i++) for(int j = 0; j <= i; j++)
      dmet[r*nm+sym2idx(i,j)] = DM(i,j);
  }
}

template void getexpmet_frechet<2>(const double*, const double*, int, double*, double*);
template void getexpmet_frechet<3>(const double*, const double*, int, double*, double*);

// -----------------------------------------------------------------------------
template<int ndim, typename T>
void getlogmet_inp(T *met){
	static_assert(ndim == 2 || ndim == 3);

	T eigval[ndim], eigvec[ndim*ndim];

	geteigsym<ndim,T>(met,eigval,eigvec);
  bool iok = false;
  if constexpr(std::is_same<T,double>::value){
    iok = !(std::isnan(log(eigval[0])) || std::isinf(log(eigval[0])));
  }else{
    iok = eigval[0].value() > 0;
    iok = iok && !std::isnan(log(eigval[0].value()));
  }
  
	if(!iok){

    #ifdef METRIS_USE_LAPACK
    if constexpr(std::is_same<T,double>::value){
      printf("## INVALID METRIC ! eigvals = ");
      dblAr1(ndim,eigval).print();
      printf(" TRY CORRECT WITH LAPACK\n");
      const int nwork = 20;
      double rwork[nwork];
      geteigsym_LAPACK<ndim>(met, nwork, rwork, eigval, eigvec);
      iok = !(std::isnan(log(eigval[0])) || std::isinf(log(eigval[0])));
      if(iok) goto fixed;
    }
    #endif

    int nnmet = (ndim*(ndim+1))/2;
    fmt::print(stderr,"Invalid metric: {}\n", MeshArray1D<T, int>(nnmet, met));
    fmt::print(stderr,"eigvals: {}\n",MeshArray1D<T, int>(ndim,eigval));

    METRIS_THROW_MSG("Negative eigenvalues");
  }

#ifdef METRIS_USE_LAPACK
fixed:
#endif

	for(int ii = 0; ii < ndim ; ii++) eigval[ii] = log(eigval[ii]);

	eig2met<ndim,T>(eigval,eigvec,met);
}


// inp and out can be the same pointer
template <int ndim, typename T>
void getexpmet_Eigen(T* inp, T* out){
  typedef Eigen::Matrix<T,ndim,ndim> MatrixN;
  MatrixN met_Eigen;
  for(int jj = 0; jj < ndim; jj++) 
    for(int ii = 0; ii < ndim; ii++) 
      met_Eigen(ii, jj) = inp[sym2idx(ii,jj)];
  met_Eigen.exp().evalTo(met_Eigen);

  for(int jj = 0; jj < ndim; jj++) 
    for(int ii = jj; ii < ndim; ii++) 
      out[sym2idx(ii,jj)] = met_Eigen(ii, jj);
}

template <int ndim, typename T>
void getexpmet_dsyevq(T* met){
  T eigval[ndim], eigvec[ndim*ndim];

  geteigsym<ndim,T>(met,eigval,eigvec);

  for(int ii = 0; ii < ndim ; ii++) eigval[ii] = exp(eigval[ii]);

  eig2met<ndim,T>(eigval,eigvec,met);
}
template void getexpmet_dsyevq<2,double>(double*);
template void getexpmet_dsyevq<2,SANS::SurrealS<2,double>>(SANS::SurrealS<2,double>*);
template void getexpmet_dsyevq<3,double>(double*);
template void getexpmet_dsyevq<3,SANS::SurrealS<3,double>>(SANS::SurrealS<3,double>*);


template <int ndim, typename T>
void getexpmet_inp(T* met){
  // In dim 3, using Eigen's exp is faster than using the eigendecomposition
  if constexpr(ndim == 3 && std::is_same<T,double>::value){
    getexpmet_Eigen<ndim, T>(met, met);
  }else{
    getexpmet_dsyevq<ndim,T>(met);
  }
}


template void getlogmet_inp<2,double>(double*);
template void getlogmet_inp<2,SANS::SurrealS<2,double>>(SANS::SurrealS<2,double>*);
template void getlogmet_inp<3,double>(double*);
template void getlogmet_inp<3,SANS::SurrealS<3,double>>(SANS::SurrealS<3,double>*);

template void getexpmet_inp<2,double>(double*);
template void getexpmet_inp<2,SANS::SurrealS<2,double>>(SANS::SurrealS<2,double>*);
template void getexpmet_inp<3,double>(double*);
template void getexpmet_inp<3,SANS::SurrealS<3,double>>(SANS::SurrealS<3,double>*);


template<>
void getexpmet_cpy<2>(const double* met ,double*  __restrict__ expm){
  for(int ii = 0; ii < 3; ii++) expm[ii] = met[ii];
  getexpmet_inp<2, double>(expm);
}
template<>
void getexpmet_cpy<3>(const double* met ,double*  __restrict__ expm){
  for(int ii = 0; ii < 6; ii++) expm[ii] = met[ii];
  getexpmet_inp<3, double>(expm);
}


} // End namespace
