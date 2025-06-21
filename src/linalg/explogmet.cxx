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

    #ifdef USE_LAPACK
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

fixed:

	for(int ii = 0; ii < ndim ; ii++) eigval[ii] = log(eigval[ii]);

	eig2met<ndim,T>(eigval,eigvec,met);
}


template <int ndim, typename T>
void getexpmet_inp(T* met){
	T eigval[ndim], eigvec[ndim*ndim];

	geteigsym<ndim,T>(met,eigval,eigvec);

	for(int ii = 0; ii < ndim ; ii++) eigval[ii] = exp(eigval[ii]);

	eig2met<ndim,T>(eigval,eigvec,met);
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

