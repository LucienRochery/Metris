//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __SRC__MESH_TYPES__
#define __SRC__MESH_TYPES__

#ifdef USE_ABSL
  #include <absl/container/flat_hash_map.h>
#else
  #include "aux_hashtab.hxx"
  #include <unordered_map>
#endif
#include <boost/multiprecision/cpp_bin_float.hpp>

#ifdef USE_METRIS_HASH
#include "utils/tuple_hashtable.hxx"
#endif

#include "types_arrays.hxx"


namespace Metris{


#ifdef USE_ABSL
typedef absl::flat_hash_map<std::tuple<int, int>, int>      HshTab_I2I;
typedef absl::flat_hash_map<std::tuple<int, int>, std::pair<int,int>> HshTab_I2I2;
typedef absl::flat_hash_map<std::tuple<int, int>, double>   HshTab_I2R;
typedef absl::flat_hash_map<std::tuple<int, int, int>, int> HshTab_I3I;
// For quality hashing in cavity. Only in 3D
typedef absl::flat_hash_map<std::tuple<int, int, int, int>, double> HshTab_I4R;
#elif USE_METRIS_HASH
typedef TupleHashTable<2,int> HshTab_I2I;
typedef std::unordered_map<std::tuple<int, int>, std::pair<int,int>, tup2_hash::hash> HshTab_I2I2;
typedef TupleHashTable<2,double> HshTab_I2R;
typedef TupleHashTable<3,int> HshTab_I3I;
//typedef std::unordered_map<std::tuple<int, int, int>, std::pair<int,int>, tup3_hash::hash> HshTab_I3Ix2;
typedef TupleHashTable<4,double> HshTab_I4R;
#else
typedef std::unordered_map<std::tuple<int, int>, int,      tup2_hash::hash> HshTab_I2I;
typedef std::unordered_map<std::tuple<int, int>, std::pair<int,int>, tup2_hash::hash> HshTab_I2I2;
typedef std::unordered_map<std::tuple<int, int>, double,   tup2_hash::hash> HshTab_I2R;
typedef std::unordered_map<std::tuple<int, int, int>, int, tup3_hash::hash> HshTab_I3I;
//typedef std::unordered_map<std::tuple<int, int, int>, std::pair<int,int>, tup3_hash::hash> HshTab_I3Ix2;
typedef std::unordered_map<std::tuple<int, int, int, int>, double, tup4_hash::hash> HshTab_I4R;
#endif

} // End namespace


#endif