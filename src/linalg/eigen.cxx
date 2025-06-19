//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "eigen.hxx"
#include "dsyevq.hxx"
#include "Metris_LAPACK.hxx"

#include "../aux_exceptions.hxx"
#include "../SANS/Surreal/SurrealS.h"


namespace Metris{

// -----------------------------------------------------------------------------

#ifdef USE_LAPACK
template<>
void geteigsym_LAPACK<2>(const double* met,int nwork,double* rwork,double* eigval,double* eigvec){
  //double eigve2[9];
  eigvec[2*0+0] = met[0];
  eigvec[2*1+0] = met[1];
  eigvec[2*1+1] = met[2];

  char c1 = 'V', c2 = 'U';
  int two = 2, info;
  // The LAPACKE version does memory allocation...
  dsyev_(&c1,&c2,&two,eigvec,&two,eigval,rwork,&nwork,&info);
  if(info != 0)METRIS_THROW_MSG(AlgoExcept(),
   "1 dsyev failed last ierro = "<<info<<"\n");

}
template<>
void geteigsym_LAPACK<3>(const double* met,int nwork,double* rwork,double* eigval,double* eigvec){
  //double eigve2[9];
  eigvec[3*0+0] = met[0];
  eigvec[3*1+0] = met[1];
  eigvec[3*1+1] = met[2];
  eigvec[3*2+0] = met[3];
  eigvec[3*2+1] = met[4];
  eigvec[3*2+2] = met[5];

  char c1 = 'V', c2 = 'U';
  int three = 3, info;
  dsyev_(&c1,&c2,&three,eigvec,&three,eigval,rwork,&nwork,&info);
  if(info != 0) METRIS_THROW_MSG(AlgoExcept(),
   "2 dsyev failed last ierro = "<<info<<"\n");

}
#endif

//template<>
//void geteigsym(const SANS::SurrealS<3,double>* __restrict__ met,
//             SANS::SurrealS<3,double>* __restrict__ eigval,
//             SANS::SurrealS<3,double>* __restrict__ eigvec){
//  int ierro = dsyevq3<SANS::SurrealS<3,double>>(met,eigvec,eigval);
//  if(ierro != 0)METRIS_THROW_MSG(AlgoExcept(),
//   "dsyevq3 FAILED INFO = "<<ierro<<"\n");
//}

// This function can take SANS::SurrealS as input.
template<int ndimn, typename T>
void geteigsym(const T* __restrict__ met,
                     T* __restrict__ eigval,
                     T* __restrict__ eigvec){

  int ierro = dsyevq<ndimn,T>(met,eigvec,eigval);
  if(ierro != 0){
    //if constexpr(std::is_same<T,double>::value){
    //  const int nwork = 20;
    //  double rwork[nwork];
    //  geteigsym<ndimn>(met, nwork, rwork, eigval, eigvec);
    //}else{
      METRIS_THROW_MSG(AlgoExcept(),
        "dsyevq3 FAILED INFO = "<<ierro<<"inp ="<<met[0]<<" "<<met[1]<<" "<<met[2]<<"\n");
    //}
  }
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
