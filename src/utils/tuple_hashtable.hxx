//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_TUPLE_HASHTABLE__
#define __METRIS_TUPLE_HASHTABLE__

#include "../types_arrays.hxx"


namespace Metris{

// Keys are n integers, value is one T
template<int n, typename T>
class TupleHashTable{

public:
  TupleHashTable();
  TupleHashTable(int nreserve);

  void reserve(int nreserve);

  // Returns <0 if not found, ival >= 0 otherwise
  int find(const uint32_t* key) const;

  // Feed ival to get value at found entry.
  T operator[](int ival);

  // Insert value at key. Returns 1 if new, 0 if already in.
  bool insert(const uint32_t* key, T val);
  // After a ihead = -find() call
  bool insert(int ihead, T val);

  void stat(int* ncol_min, int* ncol_max, double* ncol_avg, int *nempty);

  int get_nhash() const {return head.get_n();}

private:
  int nhead, nreserve;
  intAr1 head;
  MeshArray2D<uint32_t> lkeys; // 0..n-1 store keys, n store next
  MeshArray1D<T> lvals;
};


template<int n>
uint32_t TupleHashFunction(const uint32_t* key);
template<int n>
uint32_t TupleHashFunction0(const uint8_t *byte, uint32_t hash);


template<int n>
bool TupleKeyCompare(const uint32_t* key1,const uint32_t* key2);


uint32_t nextPrime(int n);

}// namespace
#endif