//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __MAIN_INCLUDES__
#define __MAIN_INCLUDES__

// Concatenate function name and _ for fortran routines

#ifndef ALWAYS_INLINE

// ALWAYS_INLINE is a macro to further encourage the compiler to inline a function

#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__clang__)
#define ALWAYS_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define ALWAYS_INLINE __forceinline
#else
#warning Not forcing inline with this compiler... (Please add this compiler to tools/always_inline.h)
#define ALWAYS_INLINE inline
#endif
#endif


#include "Mesh/Mesh.hxx"
#include "metris_constants.hxx"
#include "Arrays/aux_msharrays.hxx"
#include "ho_constants.hxx"
#include "types.hxx"
#include "msh_structs.hxx"
#include "Mesh/Mesh.hxx"

#include "aux_timer.hxx"
#include "aux_exceptions.hxx"
#include "../SANS/tools/minmax.h"


//#include "aux_utils.hxx"
//#include "low_eval.hxx"
//#include "low_topo.hxx"
//
#include "low_geo.hxx"
#include "linalg/matprods.hxx"
//#include "aux_topo.hxx"
//
////#include "aux_hashtab.hxx"
////#include "fortran_headers.hxx"
////
////#include "../libs/libmeshb7.h"
//#include <absl/container/flat_hash_map.h>
//
//#include <iostream>
//#include <string>
//#include <cstdio>
//#include <cmath>
//#include <cstdio>
//#include <cassert>
//#include <tuple>
//
//#include "../SANS/Surreal/SurrealS.h"



//#include <boost/program_options.hpp>
//#include <boost/exception/all.hpp>


//#define METRIS_MAX_DEG 2 defined in mod_hoconstants.f
//const int edgnpps[1+10] = {1 , 2 , 3  , 4  , 5  , 6  , 7  , 8   , 9   , 10  , 11};
//const int facnpps[1+10] = {1 , 3 , 6  , 10 , 15 , 21 , 28 , 36  , 45  , 55  , 66};
//const int tetnpps[1+10] = {1 , 4 , 10 , 20 , 35 , 56 , 84 , 120 , 165 , 220 , 286};

//// Edge to Vertex for FAces
//// Order : 110 011 101 ; this is a bit stupid but it is consistent with HO numbering
//extern int lnoed2___DISREGARD[6];
//extern PtrSizePair2D<int> lnoed2;
//
//// Face to Vertex ELements
//// Likewise. This means the 1st face is not opposed to the first vertex... 
//extern int lnofa3___DISREGARD[12];
//extern PtrSizePair2D<int> lnofa3;
//
//
//// EDge to Vertex for ELements
//// Order : 1100 0110 1010 1001 0101 0011
//extern int lnoed3___DISREGARD[12];
//extern PtrSizePair2D<int> lnoed3;
//

#endif