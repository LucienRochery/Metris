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
  array   = NULL;
  ialloc  = 0;
}
template<typename T, typename INT1, typename INT2>
MeshArray2D<T,INT1,INT2>::MeshArray2D(INT1 m, INT2 s){
  if(m < 0 || s <= 0)METRIS_THROW(WArgExcept());

  stride  = s;
  m1      = m;
  n1      = 0;
  nmemalc = m1*stride;
  array   = new T[nmemalc]; 
  ialloc  = 1;
}

template<typename T, typename INT1, typename INT2>
MeshArray2D<T,INT1,INT2>::MeshArray2D(INT2 s, const std::initializer_list<T> & list){
  METRIS_ASSERT(s > 0);
  METRIS_ASSERT(list.size()%s == 0);

  stride  = s;
  m1 = n1 = list.size()%s;
  nmemalc = list.size();
  array   = new T[nmemalc];

  std::copy_n(list.begin(),list.size(),array);
  ialloc  = 1;
}

template<typename T, typename INT1, typename INT2>
MeshArray2D<T,INT1,INT2>::MeshArray2D(INT1 n, INT2 s, const T* ar){
  METRIS_ASSERT(n >= 0 && s > 0); 
  stride = s;
  m1 = n1 = n; 
  nmemalc = n*s;
  array = new T[n*s];
  ialloc = 1;
  for(INT1 ii = 0; ii < nmemalc; ii++) array[ii] = ar[ii];
}

template<typename T, typename INT1, typename INT2>
MeshArray2D<T,INT1,INT2>::MeshArray2D(INT1 n, INT2 s, T* ar){
  ialloc = 0;
  m1 = n1 = n; 
  stride = s;
  array  = ar;
  nmemalc = s*n;
}
//template<typename T, typename INT1, typename INT2>
//MeshArray2D<T,INT1,INT2>::MeshArray2D(const MeshArray2D &cpy){
//  dbgid   = 0; 
//  ialloc = 0;
//  array  = cpy.array;
//  nmemalc= cpy.nmemalc;
//  stride = cpy.stride;
//  #ifdef DEBUG_ARRAYS_FULL
//    narray = NULL;
//  #endif
//}
template<typename T, typename INT1, typename INT2>
MeshArray2D<T,INT1,INT2>::MeshArray2D( MeshArray2D &&cpy){
  METRIS_ASSERT(cpy.ialloc == 1); 
  ialloc     = 1;
  cpy.ialloc = 0;
  array   = cpy.array;
  nmemalc = cpy.nmemalc;
  m1      = cpy.m1;
  n1      = cpy.n1;
  stride  = cpy.stride;
}


template<typename T, typename INT1, typename INT2>
bool MeshArray2D<T,INT1,INT2>::allocate(INT1 m, INT2 s){

  METRIS_ASSERT(m >= 0); 
  METRIS_ASSERT(s >= 0); 
  
  //INT1 m_old = m1;
  INT2 s_old = stride;

  m1      = m;
  stride  = s;

  INTL nmemalc_new = ((INTL)m)*((INTL)s);
  if(nmemalc_new <= nmemalc && stride == s) return false;

  
  bool newarr = false;

  T* array_old = array; 
  bool ialloc0 = (ialloc > 0);
  if(array_old == NULL || nmemalc_new > nmemalc){
    nmemalc = nmemalc_new;
    array  = new T[nmemalc]; 
    ialloc = 1;
    newarr = true;
  }

  //std::cout<<"m_old = "<<m_old<<" m = "<<m<<" n = "<<n1<<" s_old = "<<s_old<<" s = "<<stride<<"\n";
  if(array_old != NULL && array != array_old){
    for(INT1 ii = 0; ii < n1; ii++){
      for(INT2 jj = 0; jj < s_old; jj++){
        //printf("ii = %d jj = %d \n",ii,jj);fflush(stdout);
        array[ii*stride + jj] = array_old[ii*s_old + jj];
      }
    }
    if(ialloc0) delete[] array_old;
  }

  return newarr;

  //INT2 s_old = stride;
  //INT1 m_old = nmemalc / (s_old > 0 ? s_old : 1); // Only 0 if non init, then not used anyways

  //METRIS_ASSERT(s >= s_old);

  //stride  = s;

  //// If the array was owned by this object
  //// and it was already of sufficient size, then exit.
  //if(((int64_t)m)*((int64_t)s) < nmemalc && ialloc) return;
  //// In both other cases, we'll want to (re)allocate the array

  //if(realloc && array != NULL){ 
  //  T *array_old = array; 
  //  INT1 nmem_old = nmemalc; 
  //  nmemalc  = ((int64_t)m)*((int64_t)stride);
  //  METRIS_ENFORCE(nmemalc >= nmem_old);
  //  array = new T[nmemalc];
  //  if(array == NULL) METRIS_THROW_MSG(DMemExcept(), 
  //      "MeshArray2D failed new returns NULL for size "<<nmemalc);
  //  for(INT1 ii = 0; ii < m_old; ii++){
  //    for(INT2 jj = 0; jj < s_old; jj++){
  //      array[ii*stride + jj] = array_old[ii*s_old + jj];
  //    }
  //  }
  //  //for(INT1 ii = 0; ii < nmem_old; ii++) array[ii] = array_old[ii];
  //  if(ialloc > 0)  delete[] array_old; 
  //}else{ 
  //  nmemalc = ((int64_t)m)*((int64_t)stride);
  //  if(ialloc > 0)   delete[] array;
  //  if(nmemalc == 0) return;
  //  ialloc = 1;
  //  array = new T[nmemalc];
  //  if(array == NULL) METRIS_THROW_MSG(DMemExcept(),
  //      "MeshArray2D failed new returns NULL for size "<<nmemalc);
  //}
}

template<typename T, typename INT1, typename INT2>
void MeshArray2D<T,INT1,INT2>::free(){
  m1 = n1 = stride = nmemalc = 0;
  if(array != NULL && ialloc > 0) delete[] array;
  ialloc = 0;
  array  = NULL;
}


template<typename T, typename INT1, typename INT2>
void MeshArray2D<T,INT1,INT2>::inc_n(){
  //std::cout<<"inc_n pre nmemalc = "<<nmemalc<<" m = "<<m1<<" n = "<<n1<<" stride = "<<stride<<"\n";
  if(n1 >= m1){
    m1 = MAX(MAX(n1+1,n1 * Defaults::mem_growfac), 
                      m1 * Defaults::mem_growfac); 
    nmemalc = ((INTL)m1) * ((INTL)stride);
    T *arr_new = new T[nmemalc]; 
    for(INT1 ii = 0; ii < n1; ii++){
      for(INT2 jj = 0; jj < stride; jj++){
        arr_new[ii*stride + jj] = array[ii*stride + jj];
      }
    }
    if(ialloc) delete[] array;
    ialloc = 1;
    array = arr_new;
  }
  n1++;
  //std::cout<<"inc_n pst nmemalc = "<<nmemalc<<" m = "<<m1<<" n = "<<n1<<" stride = "<<stride<<"\n";
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
  //memcpy(out[0],array,ncopy*stride*sizeof(T));
}

template<typename T, typename INT1, typename INT2>
void MeshArray2D<T,INT1,INT2>::print(INT1 n) const{
  INT1 m = n < (nmemalc/stride) ? n : (nmemalc/stride);
  for(INT1 ii = 0; ii < m; ii++){
      std::cout<<ii<<":";
      for(INT2 jj = 0; jj < stride; jj++){
          std::cout<<" "<<array[ii*stride+jj]<<" ";
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
//  ialloc = 0;
//  array  = cpy.array;
//  nmemalc= cpy.nmemalc;
//  stride = cpy.stride;
//  return *this;
//}


template<typename T, typename INT1, typename INT2>
MeshArray2D<T,INT1,INT2>& MeshArray2D<T,INT1,INT2>::operator=(MeshArray2D &&cpy){
  METRIS_ENFORCE(ialloc == 0); 
  METRIS_ENFORCE(cpy.ialloc == 1 || cpy.array == NULL); 
  ialloc = 1;
  cpy.ialloc = 0;
  array  = cpy.array;
  nmemalc= cpy.nmemalc;
  m1     = cpy.m1;
  n1     = cpy.n1;
  stride = cpy.stride;
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
