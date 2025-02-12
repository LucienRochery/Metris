//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../Arrays/aux_msharrays.hxx"
#include "../SANS/Surreal/SurrealS.h"
#include "../metris_constants.hxx"


namespace Metris{




template<typename T, typename INT1, typename INT2>
MeshArray2D<T,INT1,INT2>::MeshArray2D(){
  nmemalc = m1 = n1 = 0;
  stride  = 0;
  array_sp = NULL;
  array    = NULL;
  array_ro = NULL; 
}
template<typename T, typename INT1, typename INT2>
MeshArray2D<T,INT1,INT2>::MeshArray2D(INT1 m, INT2 s){
  METRIS_ENFORCE_MSG(m >= 0 && s > 0, "MeshArray2D initialized with m < 0 or s <= 0")
  stride  = s;
  m1      = m;
  n1      = 0;
  nmemalc = ((INTL)m1)*((INTL)stride);
  array_sp = cpp17_make_shared<T[]>(nmemalc);
  array    = array_sp.get();
  array_ro = array; 
}

#if 0
template<typename T, typename INT1, typename INT2>
MeshArray2D<T,INT1,INT2>::MeshArray2D(INT2 s, const std::initializer_list<T> & list){
  METRIS_ASSERT(s > 0);
  METRIS_ASSERT(list.size()%s == 0);

  stride  = s;
  m1 = n1 = list.size()%s;
  nmemalc = list.size();
  array   = new T[nmemalc];

  std::copy_n(list.begin(),list.size(),array);
  iowner  = true;
}
#endif

// Dangerous: unamanaged memory. The caller is responsible. 
template<typename T, typename INT1, typename INT2>
MeshArray2D<T,INT1,INT2>::MeshArray2D(INT1 n, INT2 s, const T* ar){
  m1 = n1 = n; 
  stride   = s;
  array_sp = NULL;
  array    = NULL;
  array_ro = ar; 
  nmemalc  = ((INTL)s)*((INTL)n);
}

// Dangerous: unamanaged memory. The caller is responsible. 
template<typename T, typename INT1, typename INT2>
MeshArray2D<T,INT1,INT2>::MeshArray2D(INT1 n, INT2 s, T* ar){
  m1 = n1 = n; 
  stride   = s;
  array_sp = NULL;
  array    = ar;
  array_ro = array; 
  nmemalc  = ((INTL)s)*((INTL)n);
}

//template<typename T, typename INT1, typename INT2>
//MeshArray2D<T,INT1,INT2>::MeshArray2D(const MeshArray2D &cpy){
//  dbgid   = 0; 
//  iowner = 0;
//  array  = cpy.array;
//  nmemalc= cpy.nmemalc;
//  stride = cpy.stride;
//  #ifdef DEBUG_ARRAYS_FULL
//    narray = NULL;
//  #endif
//}
template<typename T, typename INT1, typename INT2>
MeshArray2D<T,INT1,INT2>::MeshArray2D(MeshArray2D &&cpy){
  *this = std::move(cpy);
}






// Reallocates to different major size and stride, copying old info. 
// If the new stride is smaller, then the old data is truncated 
template<typename T, typename INT1, typename INT2>
bool MeshArray2D<T,INT1,INT2>::allocate(INT1 m, INT2 s){

  METRIS_ASSERT_MSG(m >= n1," Trying to allocate size "<<m<<" < n1 = "<<n1); 
  METRIS_ASSERT(s >= 0); 

  // No need to reallocate nor copy if the stride hasn't changed and we have enough room.
  if(m <= m1 && s == stride) return false;

  INTL nmemalc_new = ((INTL)m)*((INTL)s);
  // No need to reallocate, only to copy, if the overall size hasn't increased
  bool irealloc = nmemalc_new > nmemalc;
  std::shared_ptr<T[]> new_array_sp = irealloc ? cpp17_make_shared<T[]>(nmemalc_new)
                                               : array_sp;
  T* new_array = irealloc ? new_array_sp.get()
                          : array;

  // Now we copy. 
  METRIS_ASSERT(array != NULL || n1 <= 0);
  // If the new stride is smaller, we copy from left to right. 
  // If is is greater, we copy from right to left, otherwise we overwrite data
  // to copy.
  // Note, this is only necessary if we haven't reallocated. But no harm otherwise.
  if(s <= stride){
    INT2 scopy = s;
    for(INT1 ii = 0; ii < n1; ii++){
      for(INT2 jj = 0; jj < scopy; jj++){
        new_array[ii*s + jj] = array[ii*stride + jj];
      }
    }
  }else{
    INT2 scopy = stride;
    for(INT1 ii = n1-1; ii >= 0; ii--){
      for(INT2 jj = 0; jj < scopy; jj++){
        new_array[ii*s + jj] = array[ii*stride + jj];
      }
    }
  }

  m1       = m;
  stride   = s;
  nmemalc  = nmemalc_new;
  array_sp = new_array_sp;
  array    = new_array;
  array_ro = array; 

  return irealloc;
}

template<typename T, typename INT1, typename INT2>
void MeshArray2D<T,INT1,INT2>::free(){
  array_sp.reset();
  array    = NULL;
  array_ro = NULL; 
  m1 = n1 = stride = nmemalc = 0;
}


template<typename T, typename INT1, typename INT2>
void MeshArray2D<T,INT1,INT2>::inc_n(){
  if(n1 >= m1){
    INTL m1_new = MAX(MAX(n1+1,n1 * Defaults::mem_growfac), 
                      m1 * Defaults::mem_growfac); 
    this->allocate(m1_new,stride);
  }
  n1++;
}

template<typename T, typename INT1, typename INT2>
void MeshArray2D<T,INT1,INT2>::fill(INT1 n, INT2 s, T x){
  METRIS_ASSERT(n <= m1 || s == 0); 
  METRIS_ASSERT(s <= stride || n == 0); 
  for(INT1 ii = 0; ii < n; ii++) 
      for(INT2 jj = 0; jj < s; jj++)
        (*this)(ii,jj) = x;
}
template<typename T, typename INT1, typename INT2>
void MeshArray2D<T,INT1,INT2>::fill(T x){
  for(INTL i = 0;i < nmemalc; i++) array[i] = x;
}


template<typename T, typename INT1, typename INT2>
void MeshArray2D<T,INT1,INT2>::copyTo(MeshArray2D<T,INT1,INT2> &out, INT1 ncopy) const{
  METRIS_ENFORCE(ncopy > 0);
  
  if(ncopy < 0) ncopy = nmemalc/stride;
  if(out.get_stride() < stride) METRIS_THROW_MSG(WArgExcept(), 
                       "Out stride = " << out.get_stride()<<" this = "<<stride);
  if(out.size()/out.get_stride() < ncopy) METRIS_THROW_MSG(DMemExcept(),
                       "Increase out size or decrease ncopy");
  
  for(INT1 ii = 0; ii < ncopy; ii++){
    for(INT2 jj = 0; jj < stride; jj++){
      out[ii][jj] = (*this)[ii][jj];
    }
  }
}

template<typename T, typename INT1, typename INT2>
void MeshArray2D<T,INT1,INT2>::print(INT1 n) const{
  INT1 m = n < (nmemalc/stride) ? n : (nmemalc/stride);
  for(INT1 ii = 0; ii < m; ii++){
      std::cout<<ii<<":";
      for(INT2 jj = 0; jj < stride; jj++){
          std::cout<<" "<<array_ro[ii*stride+jj]<<" ";
      }
      std::cout<<"\n";
  }
}
template<typename T, typename INT1, typename INT2>
void MeshArray2D<T,INT1,INT2>::print() const{
  this->print(nmemalc); 
}



template<typename T, typename INT1, typename INT2>
MeshArray2D<T,INT1,INT2>::~MeshArray2D(){
  free();
}

//template<typename T, typename INT1, typename INT2>
//MeshArray2D<T,INT1,INT2>& MeshArray2D<T,INT1,INT2>::operator=(const std::initializer_list<T> & list){
//  if(list.size()%stride != 0) METRIS_THROW_MSG(WArgExcept(),
//                               "INITIALIZER LIST NOT A MULTIPLE OF STRIDE");
//  if(list.size() > nmemalc)   METRIS_THROW_MSG(WArgExcept(),
//                               "INITIALIZER LIST TOO LARGE");
//  std::copy_n(list.begin(),list.size(),array);
//  nmemalc = list.size();
//  return *this;
//}

//template<typename T, typename INT1, typename INT2>
//MeshArray2D<T,INT1,INT2>& MeshArray2D<T,INT1,INT2>::operator=(const MeshArray2D &cpy){
//  iowner = 0;
//  array  = cpy.array;
//  nmemalc= cpy.nmemalc;
//  stride = cpy.stride;
//  return *this;
//}


template<typename T, typename INT1, typename INT2>
MeshArray2D<T,INT1,INT2>& MeshArray2D<T,INT1,INT2>::operator=(MeshArray2D &&cpy){
  this->free();
  array_sp = std::move(cpy.array_sp);
  array    = array_sp.get();
  array_ro = array; 
  m1       = cpy.m1;
  n1       = cpy.n1;
  stride   = cpy.stride;
  nmemalc  = cpy.nmemalc;
  cpy.free();
  return *this;
}

// ent2poi: large n, small s. ent2tag: small n, large s
template class MeshArray2D<int,int32_t,int32_t>;
template class MeshArray2D<int,int64_t,int32_t>;
template class MeshArray2D<int,int32_t,int64_t>;
// Coordinates, metric, Lag2Bez work arrays : large n but small n
template class MeshArray2D<double,int32_t,int32_t>;
template class MeshArray2D<double,int64_t,int32_t>;
template class MeshArray2D<SANS::SurrealS<2,double>,int32_t,int32_t>;
template class MeshArray2D<SANS::SurrealS<3,double>,int32_t,int32_t>;



}// End namespace
