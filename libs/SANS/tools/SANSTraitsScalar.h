// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef METRIS_TRAITS_SCALAR_H
#define METRIS_TRAITS_SCALAR_H

#include "../LinearAlgebra/DenseLinAlg/StaticSize/MatrixS_Type.h"
//#include "types_scalar.hxx"


namespace Metris
{

//Used to extract the scalar associated with a type. May not be POD, i.e. could be Surreal
template<class T>
struct Scalar { typedef T type; };


// C++ 20 required
//template<typename T>
//concept real_type = std::is_floating_point_v<T> 
//                  || requires(T a){
//  std::is_same<T,float4>::value == true || std::is_same<T,float8>::value == true;
//};

}  // namespace Metris

#endif  // METRIS_TRAITS_SCALAR_H
