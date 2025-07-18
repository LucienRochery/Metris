//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "../linalg/explogmet.hxx"
#include "../linalg/eigen.hxx"
#include "../aux_exceptions.hxx"

#include "../low_geo.hxx"
#include "../linalg/matprods.hxx"
#include "../linalg/utils.hxx"

#include <cmath>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>

namespace Metris{

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

    printf("Invalid metric: \n");
    int nnmet = (ndim*(ndim+1))/2;
    if constexpr (std::is_same<T, double>::value){
      for(int ii = 0 ; ii < nnmet; ii++) printf(" %23.15e ",met[ii]); 
    }else{
      for(int ii = 0 ; ii < nnmet; ii++) std::cout<<met[ii]<<" ";
    }
    std::cout<<"\n";

    if constexpr(std::is_same<T,double>::value){
      std::cout<<"eigvals:";
      dblAr1(ndim,eigval).print();
    }

    METRIS_THROW_MSG(RealExcept(),"Negative eigenvalues");
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

