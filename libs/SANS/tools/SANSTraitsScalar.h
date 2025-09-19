// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SANSTRAITSSCALAR_H
#define SANSTRAITSSCALAR_H

#include "../LinearAlgebra/DenseLinAlg/StaticSize/MatrixS_Type.h"
//#include "types_scalar.hxx"


namespace SANS
{

//Used to extract the scalar associated with a type. May not be POD, i.e. could be Surreal
template<class T>
struct Scalar { typedef T type; };


//template<class T>
//struct Scalar< DLA::MatrixD<T> >
//{
//  typedef typename Scalar<T>::type type;
//};
//
//template<class T>
//struct Scalar< DLA::VectorD<T> >
//{
//  typedef typename Scalar<T>::type type;
//};

template<int M, int N, class T>
struct Scalar< DLA::MatrixS<M,N,T> >
{
  typedef typename Scalar<T>::type type;
};

template<int M, class T>
struct Scalar< DLA::MatrixSymS<M,T> >
{
  typedef typename Scalar<T>::type type;
};

template<int M, class T>
struct Scalar< DLA::VectorS<M,T> >
{
  typedef typename Scalar<T>::type type;
};

// C++ 20 required
//template<typename T>
//concept real_type = std::is_floating_point_v<T> 
//                  || requires(T a){
//  std::is_same<T,float4>::value == true || std::is_same<T,float8>::value == true;
//};

}  // namespace SANS

#endif  // SANSTRAITSSCALAR_H
