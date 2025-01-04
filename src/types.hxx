//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __SRC__MESH_TYPES__
#define __SRC__MESH_TYPES__

#include "Arrays/aux_msharrays.hxx"
#include <absl/container/flat_hash_map.h>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <egads.h>

namespace Metris{


#ifdef METRIS_LARGE_MESH

#define INT1 int64_t
#define INT2 int32_t

#else

#define INT1 int32_t
#define INT2 int32_t

#endif

using intAr1 = MeshArray1D<int   ,INT1>;
using bolAr1 = MeshArray1D<int   ,INT1>;
using egoAr1 = MeshArray1D<ego   ,INT1>;
using dblAr1 = MeshArray1D<double,INT1>;

using intAr2  = MeshArray2D<int   ,INT1,INT2>;
using intAr2r = MeshArray2D<int   ,INT2,INT1>;
using dblAr2  = MeshArray2D<double,INT1,INT2>;

//template<typename INT1 = int, typename INT2 = int>
//using intAr3 = MeshArray3D<int   ,INT1,INT2>;
//template<typename INT1 = int, typename INT2 = int>
//using dblAr3 = MeshArray3D<double,INT1,INT2>;

#undef INT1
#undef INT2



//template<typename INT1 = int>
//using intAr1 = MeshArray1D<int   ,INT1>;
//template<typename INT1 = int>
//using egoAr1 = MeshArray1D<ego   ,INT1>;
//template<typename INT1 = int>
//using dblAr1 = MeshArray1D<double,INT1>;

//template<typename INT1 = int, typename INT2 = int>
//using intAr2 = MeshArray2D<int   ,INT1,INT2>;
//template<typename INT1 = int, typename INT2 = int>
//using dblAr2 = MeshArray2D<double,INT1,INT2>;

//template<typename INT1 = int, typename INT2 = int>
//using intAr3 = MeshArray3D<int   ,INT1,INT2>;
//template<typename INT1 = int, typename INT2 = int>
//using dblAr3 = MeshArray3D<double,INT1,INT2>;




//typedef MeshArray1D<int   > intAr1;
//typedef MeshArray1D<ego   > egoAr1;
//typedef MeshArray1D<double> dblAr1;

//typedef MeshArray2D<int   > intAr2;
//typedef MeshArray2D<double> dblAr2;

typedef MeshArray3D<int   > intAr3;
typedef MeshArray3D<double> dblAr3;

typedef absl::flat_hash_map<std::tuple<int, int>, int>      HshTabInt2;
typedef absl::flat_hash_map<std::tuple<int, int>, double>   HshTabDbl2;
typedef absl::flat_hash_map<std::tuple<int, int, int>, int> HshTabInt3;


} // End namespace


#endif