//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_TYPES_ARRAYS__
#define __METRIS_TYPES_ARRAYS__

#include "Arrays/aux_msharrays.hxx"
#include <egads.h>

#ifdef METRIS_LARGE_MESH

#define METRIS_INT1 int64_t
#define METRIS_INT2 int32_t

#else

#define METRIS_INT1 int32_t
#define METRIS_INT2 int32_t

#endif

namespace Metris{


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


template<typename T> class WorkArray1D;
using intWrkAr1 = WorkArray1D<int   >;
using dblWrkAr1 = WorkArray1D<double>;



}

#endif