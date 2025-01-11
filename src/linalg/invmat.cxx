//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "../linalg/invmat.hxx"

#include "../linalg/det.hxx"
#include "../linalg/sym3idx.hxx"

#include "../aux_exceptions.hxx"

#include "../SANS/Surreal/SurrealS.h"
#include <lapacke.h>


namespace Metris{

/*
	LAPACK WRAPPERS & AUX
*/
// -----------------------------------------------------------------------------
// met must be positive definite. Otherwise use inv3sym. 
void invspd(int nmat, double met[]){
	char c = 'U';
	int info;

	LAPACK_dpptrf(&c,&nmat,met,&info);
  if(info != 0) METRIS_THROW_MSG(AlgoExcept(),
  	"LAPACK_dpptrf FAILED INFO = "<<info<<"\n" 
    <<"mat = "<<met[0]<<" "<<met[1]<<" "<<met[2]<<"\n");

	LAPACK_dpptri(&c,&nmat,met,&info);
  if(info != 0) METRIS_THROW_MSG(AlgoExcept(),
  	"LAPACK_dpptri FAILED INFO = "<<info<<"\n");
}
// -----------------------------------------------------------------------------
// met must be positive definite. Otherwise use inv3sym. 
void inv2spd(double met[]){
  invspd(2,met);
	//double det = detsym<2>(met); 
	//double nrm0 = (abs(met[0]) + abs(met[1]) + abs(met[2]))/3;
	//nrm0 = 1 > nrm0 ? 1 : nrm0;
	//if(det < 1.0e-16*nrm0){
	//	printf("## SPD MATRIX WITH NEGATIVE DETERMINANT %15.7e ?? MATRIX = \n",det);
	//	printf("%15.7e %15.7e %15.7e\n",met[0],met[1],met[2]);
	//	exit(1);
	//}

	//double tmp = met[0];
	//met[0] =  met[2] / det ;
	//met[1] = -met[1] / det;
	//met[2] =  tmp    / det;

	//#ifndef NDEBUG
	//if(isnan(met[0]) || isnan(met[1]) || isnan(met[2])){
	//	printf("## METRIC IS NAN! %15.7e  %15.7e  %15.7e  DET = %15.7e\n",met[0],met[1],met[2],det);
	//	exit(1);
	//}
	//#endif
}

// 
template <int n>
void invsym(double* met){
	char c = 'U';
	int info;
	int ipiv[n];
	double work[n];
  int nmat = n;

	LAPACK_dsptrf(&c,&nmat,met,ipiv,&info);
  if(info != 0) METRIS_THROW_MSG(AlgoExcept(),
  	"LAPACK_dsytrf FAILED INFO = "<<info<< "det = "<<detsym<n>(met));

	LAPACK_dsptri(&c,&nmat,met,ipiv,work,&info);
  if(info != 0) METRIS_THROW_MSG(AlgoExcept(),
  	"LAPACK_dsytri FAILED INFO = "<<info);
}
template void invsym<2>(double *met);
template void invsym<3>(double *met);

template<typename T>
void inv3sym(T *met, T *inv){
  T det = detsym<3,T>(met);
  inv[sym3idx(0,0)] =   (met[sym3idx(1,1)]*met[sym3idx(2,2)] - met[sym3idx(1,2)]*met[sym3idx(2,1)]) / det;
  inv[sym3idx(1,0)] = - (met[sym3idx(1,0)]*met[sym3idx(2,2)] - met[sym3idx(1,2)]*met[sym3idx(2,0)]) / det;
  inv[sym3idx(2,0)] =   (met[sym3idx(1,0)]*met[sym3idx(2,1)] - met[sym3idx(1,1)]*met[sym3idx(2,0)]) / det;

  inv[sym3idx(1,1)] =   (met[sym3idx(0,0)]*met[sym3idx(2,2)] - met[sym3idx(0,2)]*met[sym3idx(2,0)]) / det;
  inv[sym3idx(2,1)] = - (met[sym3idx(0,0)]*met[sym3idx(2,1)] - met[sym3idx(2,0)]*met[sym3idx(0,1)]) / det;

  inv[sym3idx(2,2)] =   (met[sym3idx(0,0)]*met[sym3idx(1,1)] - met[sym3idx(0,1)]*met[sym3idx(1,0)]) / det;
}
template void inv3sym<double>(double *met, double *inv);
template void inv3sym<SANS::SurrealS<3,double>>(SANS::SurrealS<3,double> *met, SANS::SurrealS<3,double> *inv);



// Matrix stored line first in C fashion
void invmat(int n, double mat[]){
	METRIS_ENFORCE_MSG(n <= 3, "invmat expecting n <= 3");
	int ipiv[3];
	constexpr int nwork = 20;
	double rwork[nwork];
	int info;

	LAPACK_dgetrf(&n,&n,mat,&n,ipiv,&info);
	if(n == 3){
	  if(info != 0) METRIS_THROW_MSG(AlgoExcept(),
	  	"LAPACK_dpptrf FAILED INFO = "<<info<<" matrix = "
	  	<<mat[0]<<" "<<mat[1]<<" "<<mat[2]<<" "
	  	<<mat[3]<<" "<<mat[4]<<" "<<mat[5]<<" "
	  	<<mat[6]<<" "<<mat[7]<<" "<<mat[8]<<" det = "<<detmat<3>(mat));
	}else{
	  if(info != 0) METRIS_THROW_MSG(AlgoExcept(),
	  	"LAPACK_dpptrf FAILED INFO = "<<info<<" matrix = "
	  	<<mat[0]<<" "<<mat[1]<<" "<<mat[2]<<" "
	  	<<mat[3]<<" "<<mat[4]<<" det = "<<detmat<2>(mat));
	}

	LAPACK_dgetri(&n,mat,&n,ipiv,rwork,&nwork,&info);
  if(info != 0) METRIS_THROW_MSG(AlgoExcept(),
  	"LAPACK_dgetri FAILED INFO = "<<info);
}


template<>
void invmat<2>(double *mat){
  double det = mat[0]*mat[3] - mat[1]*mat[2];
  if(abs(det) < 1.0e-30){
    printf("## SINGULAR MATRIX DETERMINANT\n");
    for(int ii = 0; ii < 2; ii++){
      for(int jj = 0; jj < 2; jj++){
        printf("%d,%d: %24.16e \n",ii,jj,mat[ii*2+jj]);
      }
    }
    printf("det = %24.16e\n",det);
    METRIS_THROW(GeomExcept());
  } 
  double tmp = mat[0];
  mat[0] = mat[3] / det;
  mat[3] = tmp / det;
  mat[1] = -mat[1] / det;
  mat[2] = -mat[2] / det;
}

template<>
void invmat<3>(double *mat){
  METRIS_THROW_MSG(TODOExcept(),"invmat 3");
}

template<>
void invmat<1>(double *mat){
  if(abs(*mat) < 1.0e-30) METRIS_THROW_MSG(GeomExcept(), "Singular matrix det = "<<*mat);
  *mat = 1.0 / (*mat);
}

} // End namespace
