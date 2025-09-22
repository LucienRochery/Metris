//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "tuple_hashtable.hxx"

#include "aux_misc.hxx"
#include "../aux_exceptions.hxx"
#include "../metris_defaults.hxx"

#include "../../libs/rapidhash/rapidhash.h"

#include <boost/math/special_functions/prime.hpp>

namespace Metris{




template<int n, typename T>
TupleHashTable<n,T>::TupleHashTable(){
  nhead = -1;
  nreserve = 0;
}


template<int n, typename T>
TupleHashTable<n,T>::TupleHashTable(int nreserve_){
  METRIS_ASSERT(nreserve_ >= 1);
  nreserve = nreserve_;

  nhead = MAX(nreserve / Defaults::TupleHashTable_ncoll, 1);
  uint32_t next_prime = nextPrime(nhead);
  nhead = (uint32_t) next_prime;

  head.allocate(nhead);
  head.set_n(nhead);
  head.fill(-1);

  lkeys.allocate(nreserve, 1+n);
  lkeys.set_n(0);

  lvals.allocate(nreserve);
  lvals.set_n(0);
}


template<int n, typename T>
void TupleHashTable<n,T>::reserve(int nreserve_){
  METRIS_ASSERT(nreserve_ >= 1);

  if(nreserve_ < nreserve) return;
  nreserve = nreserve_;

  lkeys.allocate(nreserve, 1+n);
  lvals.allocate(nreserve);

  if(nhead == -1){
    nhead = MAX(nreserve / Defaults::TupleHashTable_ncoll, 1);
    uint32_t next_prime = nextPrime(nhead);
    nhead = (uint32_t) next_prime;

    head.allocate(nhead);
    head.set_n(nhead);
    head.fill(-1);
    lkeys.set_n(0);
    lvals.set_n(0);
  }

}




template<int n, typename T>
int TupleHashTable<n,T>::find(const uint32_t* key) const{

  uint32_t hash = TupleHashFunction<n>((uint32_t *) key) % nhead;

  int ientry = head[hash];


  while(ientry >= 0){
    if(TupleKeyCompare<n>(key, lkeys[ientry])) return ientry;
    ientry = lkeys(ientry, n);
  }

  // Here if initially -1 or only at the end, in all cases:
  return -1-hash;


}

template<int n, typename T>
bool TupleHashTable<n,T>::insert(const uint32_t* key, T val){

  uint32_t hash = TupleHashFunction<n>((uint32_t *) key) % nhead;

  int ientry = head[hash];

  // First key of this hash
  if(ientry < 0){
    int ilist = lkeys.get_n();
    METRIS_ASSERT(ilist == lvals.get_n());
    lkeys.inc_n();
    lvals.inc_n();

    lvals[ilist] = val;
    for(int ii = 0; ii < n; ii++) lkeys(ilist, ii) = key[ii];
    lkeys(ilist, n) = -1;

    head[hash] = ilist;

    return true;
  }

  // Key hash already exists
  int ientrp = ientry;
  while(ientry >= 0){
    if(TupleKeyCompare<n>(key, lkeys[ientry])) return false;
    ientrp = ientry;
    ientry = lkeys(ientry, n);
  }

  // ientrp is last 
  int ilist = lkeys.get_n();
  METRIS_ASSERT(ilist == lvals.get_n());
  lkeys.inc_n();
  lvals.inc_n();

  lvals[ilist] = val;
  for(int ii = 0; ii < n; ii++) lkeys(ilist, ii) = key[ii];
  lkeys(ilist, n) = -1;

  lkeys(ientrp, n) = ilist;

  return true;
}


template<int n, typename T>
T TupleHashTable<n,T>::operator[](int ilist){
  return lvals[ilist];
}

template<int n, typename T>
void TupleHashTable<n,T>::stat(int* ncol_min, int* ncol_max, double* ncol_avg, int *nempty){
  *nempty = 0;
  *ncol_avg = 0;
  *ncol_min = lkeys.get_n()+1;
  *ncol_max = -1;
  for(int ihash = 0; ihash < nhead; ihash++){
    int ilist = head[ihash];
    if(ilist < 0){
      (*nempty)++;
      continue;
    }
    int ncol = 0;
    while(ilist >= 0){
      ncol++;
      ilist = lkeys(ilist, n);
    }
    (*ncol_min) = MIN(*ncol_min, ncol);
    (*ncol_max) = MAX(*ncol_max, ncol);
    (*ncol_avg) += ncol;
  }
  if(*nempty < nhead) (*ncol_avg) /= (nhead - *nempty);
  return;
}


template class TupleHashTable<2, int>;
template class TupleHashTable<2, double>;
template class TupleHashTable<3, int>;
template class TupleHashTable<4, double>;



/*

Helpers

*/

uint32_t nextPrime(int n){

  int ipow = log2(n);
  return pow(2, ipow+1);

  int n1 = 0, n2 = 10000;

  uint32_t prime;
  while(n1 < n2 - 1){
    int navg = (n1 + n2) / 2;
    prime = boost::math::prime(navg);
    printf("navg = {} n1 = {} n2 = {} prime = {}\n",navg,n1,n2,prime);
    if(prime >= n){
      n2 = navg;
    }else{
      n1 = navg;
    }
  }

  // prime was evaluated at n1.
  if(prime < n) prime = boost::math::prime(n1+1);

  return prime;
}


// Implements FNV hash. 
// https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
// has a good summary of parameter values.
// This is for 32 bit hashing:
static const uint32_t FNV_prime = 16777619;
static const uint32_t FNV_offset_basis = 2166136261;

template<int n>
uint32_t TupleHashFunction(const uint32_t* key){

  return rapidhash(key, n*sizeof(uint32_t));

  uint32_t hash_val = FNV_offset_basis;

  return TupleHashFunction0<n>((uint8_t*) key, hash_val);
}

template<int n>
uint32_t TupleHashFunction0(const uint8_t *byte, uint32_t hash){

  hash ^= *byte;
  hash *= FNV_prime;

  if constexpr (n == 1){
    return hash;
  }else{
    return TupleHashFunction0<n-1>(byte++, hash);
  }
}

template uint32_t TupleHashFunction0<1>(const uint8_t* byte, uint32_t hash);
template uint32_t TupleHashFunction0<2>(const uint8_t* byte, uint32_t hash);
template uint32_t TupleHashFunction0<3>(const uint8_t* byte, uint32_t hash);
template uint32_t TupleHashFunction0<4>(const uint8_t* byte, uint32_t hash);

template uint32_t TupleHashFunction<1>(const uint32_t* key);
template uint32_t TupleHashFunction<2>(const uint32_t* key);
template uint32_t TupleHashFunction<3>(const uint32_t* key);
template uint32_t TupleHashFunction<4>(const uint32_t* key);



template<int n>
bool TupleKeyCompare(const uint32_t* key1,const uint32_t* key2){
  if constexpr (n == 1){
    return *key1 == *key2;
  }else if (n==2){
    return key1[1] == key2[1] && key1[0] == key2[0];
  }else{
    return key1[0] == key2[0] && TupleKeyCompare<(n-1)>(key1 + 1, key2 + 1);
  }
}

//template<>
//bool TupleKeyCompare<1>(const uint32_t* key1,const uint32_t* key2){
//  printf("debug n = 1\n");
//  fflush(stdout);
//  return key1[0] == key2[0];
//}
template bool TupleKeyCompare<1>(const uint32_t* key1,const uint32_t* key2);
template bool TupleKeyCompare<2>(const uint32_t* key1,const uint32_t* key2);
template bool TupleKeyCompare<3>(const uint32_t* key1,const uint32_t* key2);
template bool TupleKeyCompare<4>(const uint32_t* key1,const uint32_t* key2);


}// namespace