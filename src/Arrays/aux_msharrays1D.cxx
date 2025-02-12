//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../Arrays/aux_msharrays.hxx"
#include "../SurrealS_inc.hxx"
#include <egads.h>
#include "../metris_constants.hxx"
#include <memory>


namespace Metris{


template<typename T,typename INT1>
MeshArray1D<T,INT1>::MeshArray1D(){
  n1 = m1 = 0;
  array    = NULL;
  array_ro = NULL;
  array_sp = NULL;
}

template<typename T,typename INT1>
MeshArray1D<T,INT1>::MeshArray1D(INT1 m){
  METRIS_ENFORCE_MSG(m >= 0, "Negative size passed to MeshArray1D");
  m1 = m;
  n1 = 0;
  // Unfortunately, make_shared for array types is only a c++20+ feature
  array_sp = cpp17_make_shared<T[]>(m);
  array    = array_sp.get();
  array_ro = array; 
}

template<typename T,typename INT1>
MeshArray1D<T,INT1>::MeshArray1D(INT1 n, T *a){
  array_sp = NULL;
  array    = a;
  array_ro = array;
  n1 = m1 = n;
}

template<typename T,typename INT1>
MeshArray1D<T,INT1>::MeshArray1D(INT1 n, const T *a){
  array_sp = NULL;
  array    = NULL;
  array_ro = a;
  n1 = m1 = n;
}

#if 0
template<typename T,typename INT1>
template<typename INT1_2, typename INT2_2>
MeshArray1D<T,INT1>::MeshArray1D(MeshArray2D<T,INT1_2,INT2_2> &arr2){
  using INTL = typename MeshArray2D<T,INT1_2,INT2_2>::INTL;
  // Check our integer type is large enough.
  using IOK = typename std::conditional< ( std::numeric_limits<INTL>::max()
                                         <=std::numeric_limits<INT1>::max()),
                                          std::true_type, std::false_type>::type;
  static_assert(std::is_same<IOK,std::true_type>::value);

  m1 = arr2.size();
  n1 = ((INT1) arr2.get_n()) * ((INT1) arr2.get_stride());

  array_sp = arr2.get_sp();
  array    = array_sp.get();
  array_ro = array; 
}
#endif


template<typename T,typename INT1>
MeshArray1D<T,INT1>::MeshArray1D(const std::initializer_list<T>& list){
  n1 = m1 = list.size();
  array_sp = cpp17_make_shared<T[]>(m1);
  array    = array_sp.get();
  array_ro = array; 
  std::copy_n(list.begin(),list.size(),array);
}

template<typename T,typename INT1>
MeshArray1D<T,INT1>::MeshArray1D(MeshArray1D &&cpy){
  *this = std::move(cpy);
}

//template<typename T,typename INT1>
//void MeshArray1D<T,INT1>::set_buffer(INT1 m, T *a){
//  METRIS_ASSERT(m >= 0);
//  this->free();
//  array   = a;
//  iowner  = false;
//  n1 = 0;
//  m1 = m;
//}

template<typename T,typename INT1>
void MeshArray1D<T,INT1>::set_sp(INT1 n, std::shared_ptr<T[]> a){
  METRIS_ASSERT(n >= 0);
  this->free();
  array_sp = a;
  array    = array_sp.get();
  m1 = n;
  n1 = 0;
}


template<typename T,typename INT1>
MeshArray1D<T,INT1>& MeshArray1D<T,INT1>::operator=(const MeshArray1D &cpy){
  this->free();
  array_sp = cpy.array_sp;
  array    = array_sp.get();
  array_ro = array; 
  m1 = cpy.m1;
  n1 = cpy.n1;
  return *this;
}

template<typename T,typename INT1>
MeshArray1D<T,INT1>& MeshArray1D<T,INT1>::operator=(MeshArray1D &&cpy){
  this->free();
  array_sp = std::move(cpy.array_sp);
  array    = array_sp.get();
  array_ro = array; 
  m1 = cpy.m1;
  n1 = cpy.n1;
  cpy.free();
  return *this;
}


template<typename T,typename INT1>
bool MeshArray1D<T,INT1>::allocate(INT1 m){
  METRIS_ASSERT(m >= 0);

  if(m <= m1 || m == 0) return false;

  std::shared_ptr<T[]> new_array_sp = cpp17_make_shared<T[]>(m);
  T* new_array = new_array_sp.get();

  METRIS_ASSERT(array != NULL || n1 <= 0);
    
  for(int ii = 0; ii < n1; ii++) new_array[ii] = array[ii];
  
  array_sp = new_array_sp;
  array    = new_array;
  array_ro = array; 
  
  m1 = m;

  return true;
}

template<typename T,typename INT1>
void MeshArray1D<T,INT1>::free(){
  array_sp.reset();
  array    = NULL;
  array_ro = NULL; 
  n1 = m1 = 0;
}


//template<typename T,typename INT1>
//void MeshArray1D<T,INT1>::fill(INT1 m, T x){
//  for(INT1 ii = 0; ii < m; ii++) array[ii] = x;
//}

template<typename T,typename INT1>
void MeshArray1D<T,INT1>::fill(T x){
  for(INT1 ii=0; ii < n1; ii++) array[ii] = x;
}

template<typename T,typename INT1>
void MeshArray1D<T,INT1>::copyTo(MeshArray1D<T,INT1> &out, INT1 ncopy) const{
  if(ncopy < 0) ncopy = m1;
  if(out.size() < ncopy) 
    METRIS_THROW_MSG(DMemExcept(),"Increase out size or decrease ncopy");
  memcpy(&out[0],array,ncopy*sizeof(T));
}

template<typename T,typename INT1>
MeshArray1D<T,INT1>::~MeshArray1D(){
  this->free();
}


template<typename T,typename INT1>
void MeshArray1D<T,INT1>::print(INT1 n) const{
  INT1 m = n < m1 ? n : m1;
  for(INT1 i = 0; i < m; i++){
    std::cout<<array_ro[i]<<" ";
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
    // allocate already copies from old to new
    this->allocate(MAX(n1 + 1, m1*Defaults::mem_growfac)); 
  }
  array[n1] = val;
  n1++;
}


template<typename T,typename INT1>
T MeshArray1D<T,INT1>::pop(){
  METRIS_ENFORCE(n1 > 0); 
  n1--;
  return array_ro[n1]; 
}


#include "../aux_pp_inc.hxx"
#define T_SEQ (bool)(ego)(int)(double)(float4)(float8)(SurrealS2)(SurrealS3)
#define INT1_SEQ (int32_t)(int64_t)
#define INSTANTIATE(T,INT1) template class MeshArray1D<T,INT1>;
#define EXPAND_TUP(r,SEQ) INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ),BOOST_PP_SEQ_ELEM(1, SEQ))
BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TUP,(T_SEQ)(INT1_SEQ))
#undef INSTANTIATE

//template MeshArray1D<bool,int32_t>::MeshArray1D<bool,int32_t,int32_t>(MeshArray2D<bool,int32_t,int32_t> &arr2);
//
//#define INSTANTIATE(T,INT1,INT1_2,INT2_2)\
//template MeshArray1D<T,INT1>::MeshArray1D<T,INT1_2,INT2_2>(MeshArray2D<T,INT1_2,INT2_2> &arr2);
//#define INSTANTIATE_r(r,data,T) INSTANTIATE(T,int32_t,int32_t,int32_t)
//BOOST_PP_SEQ_FOR_EACH(INSTANTIATE_r,~,T_SEQ)
//#define INSTANTIATE_r(r,data,T) INSTANTIATE(T,int64_t,int32_t,int64_t)
//BOOST_PP_SEQ_FOR_EACH(INSTANTIATE_r,~,T_SEQ)
//#define INSTANTIATE_r(r,data,T) INSTANTIATE(T,int64_t,int64_t,int32_t)
//BOOST_PP_SEQ_FOR_EACH(INSTANTIATE_r,~,T_SEQ)
//#undef INSTANTIATE
//#undef INSTANTIATE_r


} // End namespace 
