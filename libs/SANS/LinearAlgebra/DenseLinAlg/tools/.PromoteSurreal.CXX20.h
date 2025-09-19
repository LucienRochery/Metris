// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef DLA_PROMOTE_SURREAL_H
#define DLA_PROMOTE_SURREAL_H

#include "SANS//Surreal/PromoteSurreal.h"
#include "SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS_Type.h"

namespace SANS
{

// It is important here that if the right type is a Matrix/Vector, then the resulting type must also be
// the same sized Matrix/Vector

template<std::integral Int, int M, class T>
struct promote_Surreal< Int, DLA::VectorS<M,T> > { typedef DLA::VectorS<M,T> type; };

template<std::integral Int, int M, int N, class T>
struct promote_Surreal< Int, DLA::MatrixS<M,N,T> > { typedef DLA::MatrixS<M,N,T> type; };


template<SANS::real_type RealType, int M, class T>
struct promote_Surreal< RealType, DLA::VectorS<M,T> > { typedef DLA::VectorS<M,T> type; };

template<SANS::real_type RealType,int M, int N, class T>
struct promote_Surreal< RealType, DLA::MatrixS<M,N,T> > { typedef DLA::MatrixS<M,N,T> type; };

template<SANS::real_type RealType,int M, class T>
struct promote_Surreal< RealType, DLA::MatrixSymS<M,T> > { typedef DLA::MatrixSymS<M,T> type; };


template<SANS::real_type RealType,int M, class T>
struct promote_Surreal< DLA::VectorS<M,T>, RealType > { typedef DLA::VectorS<M,T> type; };

template<SANS::real_type RealType,int M, int N, class T>
struct promote_Surreal< DLA::MatrixS<M,N,T>, RealType> { typedef DLA::MatrixS<M,N,T> type; };

template<SANS::real_type RealType,int M, class T>
struct promote_Surreal< DLA::MatrixSymS<M,T>, RealType > { typedef DLA::MatrixSymS<M,T> type; };


template<int SN, int M, class T0, class T1>
struct promote_Surreal< SurrealS<SN,T0>, DLA::VectorS<M,T1> >
{
  typedef typename promote_Surreal<SurrealS<SN,T0>, T1>::type T;
  typedef DLA::VectorS<M,T> type;
};

template<int SN, int M, int N, class T0, class T1>
struct promote_Surreal< SurrealS<SN,T0>, DLA::MatrixS<M,N,T1> >
{
  typedef typename promote_Surreal<SurrealS<SN,T0>, T1>::type T;
  typedef DLA::MatrixS<M,N,T> type;
};


template<int M, class T1, class T2>
struct promote_Surreal< DLA::VectorS<M,T1>, DLA::VectorS<M,T2> >
{
  typedef typename promote_Surreal<T1, T2>::type T;
  typedef DLA::VectorS<M,T> type;
};

template<int M, int N, class T1, class T2>
struct promote_Surreal< DLA::MatrixS<M,N,T1>, DLA::VectorS<N,T2> >
{
  typedef typename promote_Surreal<T1, T2>::type T;
  typedef DLA::VectorS<N,T> type;
};

template<int M, int N, int K, class T1, class T2>
struct promote_Surreal< DLA::MatrixS<M,N,T1>, DLA::MatrixS<N,K,T2> >
{
  typedef typename promote_Surreal<T1, T2>::type T;
  typedef DLA::MatrixS<N,K,T> type;
};

} //namespace SANS

#endif //DLA_PROMOTE_SURREAL_H
