//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../linalg/eigen.hxx"
#include "../linalg/dsyevq.hxx"

#include "../aux_exceptions.hxx"
#include "../SANS/Surreal/SurrealS.h"
#include <lapacke.h>


namespace Metris{

// -----------------------------------------------------------------------------

template<>
void geteigsym<3,double>(const double* met,int nwork,double* rwork,double* eigval,double* eigvec){
  //double eigve2[9];
  eigvec[3*0+0] = met[0];
  eigvec[3*1+0] = met[1];
  eigvec[3*1+1] = met[2];
  eigvec[3*2+0] = met[3];
  eigvec[3*2+1] = met[4];
  eigvec[3*2+2] = met[5];

  char c1 = 'V', c2 = 'U';
  int three = 3, info;
  LAPACK_dsyev(&c1,&c2,&three,eigvec,&three,eigval,rwork,&nwork,&info);
  if(info != 0) METRIS_THROW_MSG(AlgoExcept(),
  	"LAPACK_dsyev FAILED INFO = "<<info<<"\n");


}

//template<>
//void geteigsym(const SANS::SurrealS<3,double>* __restrict__ met,
//             SANS::SurrealS<3,double>* __restrict__ eigval,
//             SANS::SurrealS<3,double>* __restrict__ eigvec){
//  int ierro = dsyevq3<SANS::SurrealS<3,double>>(met,eigvec,eigval);
//  if(ierro != 0)METRIS_THROW_MSG(AlgoExcept(),
//  	"dsyevq3 FAILED INFO = "<<ierro<<"\n");
//}

template<int ndimn, typename T>
void geteigsym(const T* __restrict__ met,
                     T* __restrict__ eigval,
                     T* __restrict__ eigvec){

  //if constexpr(ndimn == 2){
//
  //}else{
    int ierro = dsyevq<ndimn,T>(met,eigvec,eigval);
    if(ierro != 0)METRIS_THROW_MSG(AlgoExcept(),
  	 "dsyevq3 FAILED INFO = "<<ierro<<"inp ="<<met[0]<<" "<<met[1]<<" "<<met[2]<<"\n"
      );
  //}
}

template void geteigsym<2,double>(const double* __restrict__ met,
                                        double* __restrict__ eigval,
                                        double* __restrict__ eigvec);
template void geteigsym<3,double>(const double* __restrict__ met,
                                        double* __restrict__ eigval,
                                        double* __restrict__ eigvec);

template void geteigsym<2,SANS::SurrealS<2,double>>(const SANS::SurrealS<2,double>* __restrict__ met,
                                                          SANS::SurrealS<2,double>* __restrict__ eigval,
                                                          SANS::SurrealS<2,double>* __restrict__ eigvec);
template void geteigsym<3,SANS::SurrealS<3,double>>(const SANS::SurrealS<3,double>* __restrict__ met,
                                                          SANS::SurrealS<3,double>* __restrict__ eigval,
                                                          SANS::SurrealS<3,double>* __restrict__ eigvec);

//template<> void geteigsym<SANS::SurrealS<3,double>>(const SANS::SurrealS<3,double>* met, 
//              int nwork, double* rwork, SANS::SurrealS<3,double>* eigval, SANS::SurrealS<3,double>* eigvec);



} // End namespace
