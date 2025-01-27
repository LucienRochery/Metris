//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../Arrays/aux_msharrays.hxx"
#include "../SurrealS_inc.hxx"
#include <egads.h>
#include "../metris_constants.hxx"

namespace Metris{



template<typename T,typename INT1>
MeshArray1D<T,INT1>::MeshArray1D(){
  n1 = m1 = 0;
  array   = NULL;
  ialloc  = 0;
}

template<typename T,typename INT1>
MeshArray1D<T,INT1>::MeshArray1D(INT1 m){
  if(m < 0)METRIS_THROW_MSG(WArgExcept(),"NEGATIVE SIZE PASSED TO MeshArray1D");
  m1 = m;
  n1 = 0; 
  array   = new T[m];
  ialloc  = 1;
}

//template<typename T,typename INT1>
//MeshArray1D<T,INT1>::MeshArray1D(const MeshArray1D &x){
//  ialloc = 0;
//  array  = x.array;
//  m1 = x.m1;
//  n1 = x.n1; 
//}

template<typename T,typename INT1>
MeshArray1D<T,INT1>::MeshArray1D(INT1 n, T *a){
  array  = a;
  ialloc = 0;
  n1 = m1 = n;
}

template<typename T,typename INT1>
MeshArray1D<T,INT1>::MeshArray1D(INT1 n, const T *a){
  ialloc  = 1;
  n1 = m1 = n;
  array = new T[n];
  for(INT1 i = 0; i < n; i++) array[i] = a[i];
}

template<typename T,typename INT1>
MeshArray1D<T,INT1>::MeshArray1D(const std::initializer_list<T> & list){
  n1 = m1 = list.size();
  array = new T[m1];
  std::copy_n(list.begin(),list.size(),array);
}

template<typename T,typename INT1>
MeshArray1D<T,INT1>::MeshArray1D( MeshArray1D &&cpy ){
  ialloc     = 1;
  cpy.ialloc = 0;
  array = cpy.array;
  m1 = cpy.m1;
  n1 = cpy.n1;
}

template<typename T,typename INT1>
void MeshArray1D<T,INT1>::set_buffer(INT1 m, T *a){
  METRIS_ASSERT(m >= 0);
  array   = a;
  ialloc  = 0;
  n1 = 0;
  m1 = m;
}

template<typename T,typename INT1>
MeshArray1D<T,INT1>& MeshArray1D<T,INT1>::operator=(const MeshArray1D &cpy){
  ialloc = 0;
  array  = cpy.array;
  m1 = cpy.m1; 
  n1 = cpy.n1; 
  return *this;
}

template<typename T,typename INT1>
MeshArray1D<T,INT1>& MeshArray1D<T,INT1>::operator=(MeshArray1D &&cpy){
  METRIS_ENFORCE(ialloc == 0 || cpy.array==NULL);
  ialloc     = 1;
  cpy.ialloc = 0;
  array  = cpy.array;
  m1 = cpy.m1;
  n1 = cpy.n1;
  return *this;
}


template<typename T,typename INT1>
bool MeshArray1D<T,INT1>::allocate(INT1 m){

  if(m <= m1) return false;
  
  METRIS_ASSERT(m >= m1 || n1 <= 0);


  bool newarr = false;
  T* array_old = array; 
  bool ialloc0 = (ialloc > 0);
  if(array_old == NULL || m > m1){
    array  = new T[m];
    ialloc = 1;
    newarr = true;
  }


  if(array_old != NULL && array != array_old){
    for(INT1 ii = 0; ii < n1; ii++){
      array[ii] = array_old[ii]; 
    }
    if(ialloc0) delete[] array_old;
  }

  m1 = MAX(m,m1);
  return newarr;
  //n1 = n;

  //if(realloc && array != NULL){
  //  T* array_old = array;
  //  array = new T[m];
  //  if(array == NULL) METRIS_THROW(DMemExcept());
  //  for(INT1 ii = 0; ii < m1; ii++) array[ii] = array_old[ii];
  //  m1 = m;
  //  if(ialloc > 0) delete[] array_old;
  //  ialloc = 1;
  //}else{
  //  m1 = m;
  //  if(ialloc > 0) delete[] array;
  //  ialloc  = 1;
  //  array   = new T[m];
  //  if(array == NULL) METRIS_THROW(DMemExcept());
  //}
}

template<typename T,typename INT1>
void MeshArray1D<T,INT1>::free(){
  if(array != NULL && ialloc == 1) delete[] array;
  ialloc  = 0;
  n1 = m1 = 0;
  array   = NULL;
}


template<typename T,typename INT1>
void MeshArray1D<T,INT1>::fill(INT1 m, T x){
  for(INT1 i = 0; i < m; i++) array[i] = x;
}

template<typename T,typename INT1>
void MeshArray1D<T,INT1>::fill(T x){
  for(INT1 i=0; i < m1; i++)
      array[i] = x;
}

template<typename T,typename INT1>
void MeshArray1D<T,INT1>::copyTo(MeshArray1D<T,INT1> &out, INT1 ncopy) const{
  if(ncopy < 0) ncopy = m1;
  if(out.size() < ncopy) METRIS_THROW_MSG(DMemExcept(),"Increase out size or decrease ncopy");
  memcpy(&out[0],array,ncopy*sizeof(T));
}

template<typename T,typename INT1>
MeshArray1D<T,INT1>::~MeshArray1D(){
  free();
}


template<typename T,typename INT1>
void MeshArray1D<T,INT1>::print(INT1 n) const{
  INT1 m = n < m1 ? n : m1;
  for(INT1 i = 0; i < m; i++){
    std::cout<<array[i]<<" ";
  }
  std::cout<<"\n";
}
template<typename T,typename INT1>
void MeshArray1D<T,INT1>::print() const{
  this->print(n1);
}


template<typename T,typename INT1>
void MeshArray1D<T,INT1>::stack(T val){
  if(n1 >= m1){
    m1 = MAX(n1+1,m1 * Defaults::mem_growfac); 
    T *arr_new = new T[m1]; 
    for(INT1 ii = 0; ii < n1; ii++) arr_new[ii] = array[ii];
    if(ialloc) delete[] array;
    ialloc = 1;
    array = arr_new;
  }
  array[n1] = val;
  n1++;
}


template<typename T,typename INT1>
T MeshArray1D<T,INT1>::pop(){
  METRIS_ASSERT(n1 > 0); 
  n1--;
  return array[n1]; 
}


#include "../aux_pp_inc.hxx"
#define T_SEQ (bool)(ego)(int)(double)(float4)(float8)(SurrealS2)(SurrealS3)
#define INT1_SEQ (int32_t)(int64_t)
#define INSTANTIATE(T,INT1) template class MeshArray1D<T,INT1>;
#define EXPAND_TUP(r,SEQ) INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ),BOOST_PP_SEQ_ELEM(1, SEQ))
BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TUP,(T_SEQ)(INT1_SEQ))

} // End namespace 
