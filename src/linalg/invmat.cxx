//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "../linalg/invmat.hxx"

#include "../linalg/det.hxx"
#include "../linalg/symidx.hxx"

#include "../aux_exceptions.hxx"
#include "../low_geo/misc.hxx"
#include "../metris_constants.hxx"

#include "../libs/SANS/Surreal/SurrealS.h"

#include <Eigen/Dense>

namespace Metris{

// -----------------------------------------------------------------------------
// met must be positive definite. Otherwise use inv3sym.
#ifdef METRIS_USE_LAPACK
int invspd_LAPACK(int nmat, double met[]){
	char c = 'U';
	int info;

	dpptrf_(&c,&nmat,met,&info);
  if(info != 0) return(abs(info));

	dpptri_(&c,&nmat,met,&info);
  if(info != 0) return(abs(info));

  return 0;
}
#endif

// inp and out can be the same
template<int ndim, typename T>
int invspd_Eigen(T *inp, T *out){
  typedef Eigen::Matrix<T, ndim, ndim> MatrixN;

  MatrixN met_Eigen;
  for(int jj = 0; jj < ndim; jj++)
    for(int ii = jj; ii < ndim; ii++)
      met_Eigen(ii,jj) = inp[sym2idx(ii,jj)];

  Eigen::LLT<MatrixN> llt(met_Eigen);
  if(llt.info() != Eigen::Success) return 1;

  MatrixN invmet_Eigen = llt.solve(Eigen::Matrix<double,ndim,ndim>::Identity());

  for(int jj = 0; jj < ndim; jj++)
    for(int ii = jj; ii < ndim; ii++)
      out[sym2idx(ii,jj)] = invmet_Eigen(ii,jj);

  return 0;
}
template int invspd_Eigen<2,double>(double *inp, double *out);
template int invspd_Eigen<3,double>(double *inp, double *out);

template<int ndim, typename T>
int invspd(T* met){
  // LAPACK never worthwhile here
  //#ifdef METRIS_USE_LAPACK
  //return invspd_LAPACK(nmat, met);
  //#else
  return invspd_Eigen<ndim>(met, met);
  //#endif
}
template int invspd<2,double>(double *met);
template int invspd<3,double>(double *met);



// Only call if matrix is known to be invertible!
// mat and invmat can be the same
// using PartialPivLU (faster but less robust);
template<int ndim, typename T>
int invmat_EigenLUPP(T* mat, T* inv){
  typedef Eigen::Matrix<T, ndim, ndim> MatrixN;
  MatrixN mat_Eigen;
  // Store mat^T in mat_Eigen
  for(int ii = 0; ii < ndim; ii++)
    for(int jj = 0; jj < ndim; jj++)
      mat_Eigen(jj,ii) = mat[ndim*ii+jj];

  Eigen::PartialPivLU<MatrixN> lu(mat_Eigen);

  MatrixN invmat_Eigen = lu.inverse();

  // Copy back invmat_Eigen^T into mat
  for(int ii = 0; ii < ndim; ii++)
    for(int jj = 0; jj < ndim; jj++)
      inv[ndim*ii+jj] = invmat_Eigen(jj,ii);

  return 0;
}
template int invmat_EigenLUPP<2,double>(double *mat, double* inv);
template int invmat_EigenLUPP<3,double>(double *mat, double* inv);

// mat and invmat can be the same
// using FullPivLU (faster but less robust)
template<int ndim, typename T>
int invmat_EigenLUFP(T* mat, T* inv){
  typedef Eigen::Matrix<T, ndim, ndim> MatrixN;
  MatrixN mat_Eigen;
  // Store mat^T in mat_Eigen
  for(int ii = 0; ii < ndim; ii++)
    for(int jj = 0; jj < ndim; jj++)
      mat_Eigen(jj,ii) = mat[ndim*ii+jj];

  Eigen::FullPivLU<MatrixN> lu(mat_Eigen);

  if(!lu.isInvertible()) return 1;

  MatrixN invmat_Eigen = lu.inverse();

  // Copy back invmat_Eigen^T into mat
  for(int ii = 0; ii < ndim; ii++)
    for(int jj = 0; jj < ndim; jj++)
      inv[ndim*ii+jj] = invmat_Eigen(jj,ii);

  return 0;
}
template int invmat_EigenLUFP<2,double>(double *mat, double* inv);
template int invmat_EigenLUFP<3,double>(double *mat, double* inv);

#ifdef METRIS_USE_LAPACK
  // Matrix stored line first in C fashion
  int invmat_LAPACK(int n, double mat[]){
  	METRIS_ENFORCE_MSG(n <= 3, "invmat expecting n <= 3");
  	int ipiv[3];
  	constexpr int nwork = 20;
    int nwork_ = nwork;
  	double rwork[nwork];
  	int info;

  	dgetrf_(&n,&n,mat,&n,ipiv,&info);
    if(info != 0) return(abs(info));

  	dgetri_(&n,mat,&n,ipiv,rwork,&nwork_,&info);
    if(info != 0) return(abs(info));

    return 0;
  }
#endif

template<>
int invmat_naive<3>(double *mat){
  // Calculate determinant
  double det = detmat<3>(mat);
  if(abs(det) < Constants::detTol) return 1;

  double invdet = 1.0 / det;

  // Calculate cofactor matrix and transpose to get adjugate
  double inv[9];
  inv[0] =  (mat[4]*mat[8] - mat[5]*mat[7]) * invdet;
  inv[1] = -(mat[1]*mat[8] - mat[2]*mat[7]) * invdet;
  inv[2] =  (mat[1]*mat[5] - mat[2]*mat[4]) * invdet;
  inv[3] = -(mat[3]*mat[8] - mat[5]*mat[6]) * invdet;
  inv[4] =  (mat[0]*mat[8] - mat[2]*mat[6]) * invdet;
  inv[5] = -(mat[0]*mat[5] - mat[2]*mat[3]) * invdet;
  inv[6] =  (mat[3]*mat[7] - mat[4]*mat[6]) * invdet;
  inv[7] = -(mat[0]*mat[7] - mat[1]*mat[6]) * invdet;
  inv[8] =  (mat[0]*mat[4] - mat[1]*mat[3]) * invdet;

  // Copy result back to input matrix
  for(int i = 0; i < 9; i++){
    mat[i] = inv[i];
  }

  return 0;
}

template<>
int invmat_naive<2>(double *mat){
  double det = mat[0]*mat[3] - mat[1]*mat[2];
  if(abs(det) < Constants::detTol) return 1;

  double tmp = mat[0];
  mat[0] = mat[3] / det;
  mat[3] = tmp / det;
  mat[1] = -mat[1] / det;
  mat[2] = -mat[2] / det;
  return 0;
}

template<>
int invmat_naive<1>(double *mat){
  if(abs(*mat) < Constants::detTol) return 1;
  *mat = 1.0 / (*mat);
  return 0;
}


template<int ndim>
int invmat(double *mat){
  // naive formula no less stable in dim 2 and much faster
  if constexpr(ndim <= 2) return invmat_naive<ndim>(mat);
  else{
    #ifdef METRIS_USE_LAPACK
    return invmat_LAPACK(ndim, mat);
    #else
    return invmat_EigenLUFP<ndim>(mat,mat);
    #endif
  }
}
template int invmat<1>(double *mat);
template int invmat<2>(double *mat);
template int invmat<3>(double *mat);

} // End namespace
