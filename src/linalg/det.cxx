//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "det.hxx"
#include <Eigen/Dense>
#include "../aux_exceptions.hxx"

namespace Metris{

template<int ndim,typename ftype>
ftype detsym_Eigen_LLT(const ftype *met){
  typedef Eigen::Matrix<ftype,ndim,ndim> MatrixN;

  MatrixN met_Eigen;
  for(int ii = 0; ii < ndim; ii++) 
    for(int jj = 0; jj < ndim; jj++) 
      met_Eigen(ii, jj) = met[sym2idx(ii,jj)];
  Eigen::LLT<MatrixN> llt(met_Eigen);
  if(llt.info() != Eigen::Success) METRIS_THROW(GeomExcept());

  MatrixN LL = llt.matrixL();

  ftype detE = LL(0,0);
  for(int ii = 1; ii < ndim; ii++) detE *= LL(ii,ii);
  return detE*detE;
}

template  double detsym_Eigen_LLT<2>(const double *met);
template  double detsym_Eigen_LLT<3>(const double *met);

#ifdef METRIS_USE_LAPACK
template<int ndim>
double detsym_LAPACK(const double* met){
  double A[ndim][ndim]; 
  for(int ii = 0; ii < ndim ; ii++){
    for(int jj = 0; jj < ndim ;jj++){
      A[ii][jj] = met[sym2idx(ii,jj)];
    }
  }
  int ipiv[ndim] = {-1};
  int info;
  int n = ndim;
  dgetrf_(&n,&n,A[0],&n,ipiv,&info);
  // It's actually not a permutation vector but lists permutations per iteration...
  int nn = (ipiv[0] != 1);
  double det = A[0][0];
  for(int ii = 1; ii < ndim; ii++){
    nn += ipiv[ii] != (ii+1);
    det*= A[ii][ii];
  }
  return nn % 2 == 0 ? det : -det; 
}

template  double detsym_LAPACK<2>(const double *met);
template  double detsym_LAPACK<3>(const double *met);
#endif


}