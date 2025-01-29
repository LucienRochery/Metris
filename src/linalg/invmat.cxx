//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "../linalg/invmat.hxx"

#include "../linalg/det.hxx"
#include "../linalg/symidx.hxx"

#include "../aux_exceptions.hxx"
#include "../metris_constants.hxx"

#include "../SANS/Surreal/SurrealS.h"
#include <lapacke.h>


namespace Metris{

/*
	LAPACK WRAPPERS & AUX
*/
// -----------------------------------------------------------------------------
// met must be positive definite. Otherwise use inv3sym. 
int invspd(int nmat, double met[]){
	char c = 'U';
	int info;

	LAPACK_dpptrf(&c,&nmat,met,&info);
  if(info != 0) return(abs(info));

	LAPACK_dpptri(&c,&nmat,met,&info);
  if(info != 0) return(abs(info));

  return 0;
}

// 
template <int n>
int invsym(double* met){
	char c = 'U';
	int info;
	int ipiv[n];
	double work[n];
  int nmat = n;

	LAPACK_dsptrf(&c,&nmat,met,ipiv,&info);
  if(info != 0) return(abs(info));

	LAPACK_dsptri(&c,&nmat,met,ipiv,work,&info);
  if(info != 0) return(abs(info));

  return 0;
}
template int invsym<2>(double *met);
template int invsym<3>(double *met);

template<typename T>
int inv3sym(T *met, T *inv){
  T det = detsym<3,T>(met);
  if(abs(det) < Constants::detTol) return 1;
  inv[sym2idx(0,0)] =   (met[sym2idx(1,1)]*met[sym2idx(2,2)] - met[sym2idx(1,2)]*met[sym2idx(2,1)]) / det;
  inv[sym2idx(1,0)] = - (met[sym2idx(1,0)]*met[sym2idx(2,2)] - met[sym2idx(1,2)]*met[sym2idx(2,0)]) / det;
  inv[sym2idx(2,0)] =   (met[sym2idx(1,0)]*met[sym2idx(2,1)] - met[sym2idx(1,1)]*met[sym2idx(2,0)]) / det;

  inv[sym2idx(1,1)] =   (met[sym2idx(0,0)]*met[sym2idx(2,2)] - met[sym2idx(0,2)]*met[sym2idx(2,0)]) / det;
  inv[sym2idx(2,1)] = - (met[sym2idx(0,0)]*met[sym2idx(2,1)] - met[sym2idx(2,0)]*met[sym2idx(0,1)]) / det;

  inv[sym2idx(2,2)] =   (met[sym2idx(0,0)]*met[sym2idx(1,1)] - met[sym2idx(0,1)]*met[sym2idx(1,0)]) / det;
  return 0;
}
template int inv3sym<double>(double *met, double *inv);
template int inv3sym<SANS::SurrealS<3,double>>(SANS::SurrealS<3,double> *met, SANS::SurrealS<3,double> *inv);



// Matrix stored line first in C fashion
int invmat(int n, double mat[]){
	METRIS_ENFORCE_MSG(n <= 3, "invmat expecting n <= 3");
	int ipiv[3];
	constexpr int nwork = 20;
	double rwork[nwork];
	int info;

	LAPACK_dgetrf(&n,&n,mat,&n,ipiv,&info);
  if(info != 0) return(abs(info)); 

	LAPACK_dgetri(&n,mat,&n,ipiv,rwork,&nwork,&info);
  if(info != 0) return(abs(info)); 

  return 0;
}


template<>
int invmat<2>(double *mat){
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
int invmat<3>(double *mat){
  METRIS_THROW_MSG(TODOExcept(),"invmat 3");
}

template<>
int invmat<1>(double *mat){
  if(abs(*mat) < Constants::detTol) return 1; 
  *mat = 1.0 / (*mat);
  return 0;
}

} // End namespace
