//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __SRC__MESH_TYPES__
#define __SRC__MESH_TYPES__

#include "Arrays/aux_msharrays.hxx"
#ifdef USE_ABSL
  #include <absl/container/flat_hash_map.h>
#else
  #include "aux_hashtab.hxx"
  #include <unordered_map>
#endif
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <egads.h>

namespace Metris{


#ifdef METRIS_LARGE_MESH

#define METRIS_INT1 int64_t
#define METRIS_INT2 int32_t

#else

#define METRIS_INT1 int32_t
#define METRIS_INT2 int32_t

#endif

using intAr1 = MeshArray1D<int   ,METRIS_INT1>;
using bolAr1 = MeshArray1D<int   ,METRIS_INT1>;
using egoAr1 = MeshArray1D<ego   ,METRIS_INT1>;
using dblAr1 = MeshArray1D<double,METRIS_INT1>;

using intAr2  = MeshArray2D<int   ,METRIS_INT1,METRIS_INT2>;
using intAr2r = MeshArray2D<int   ,METRIS_INT2,METRIS_INT1>;
using dblAr2  = MeshArray2D<double,METRIS_INT1,METRIS_INT2>;

using intLoop = Loop<int, METRIS_INT1>;

typedef MeshArray3D<int   > intAr3;
typedef MeshArray3D<double> dblAr3;

#ifdef USE_ABSL
typedef absl::flat_hash_map<std::tuple<int, int>, int>      HshTab_I2I;
typedef absl::flat_hash_map<std::tuple<int, int>, std::pair<int,int>> HshTab_I2I2;
typedef absl::flat_hash_map<std::tuple<int, int>, double>   HshTab_I2R;
typedef absl::flat_hash_map<std::tuple<int, int, int>, int> HshTab_I3I;
// For quality hashing in cavity. Only in 3D
typedef absl::flat_hash_map<std::tuple<int, int, int, int>, double> HshTab_I4R;
//typedef absl::flat_hash_map<std::tuple<int, int, int>, std::pair<int,int>> HshTab_I3Ix2;
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