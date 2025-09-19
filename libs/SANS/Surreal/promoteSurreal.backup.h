// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef PROMOTE_SURREAL_H
#define PROMOTE_SURREAL_H

#include <cstddef> // std::size_t

#include "SANS/tools/SANSTraitsScalar.h"
#include "SANS//Surreal/SurrealS_fwd.h"

// Forward declare
class SurrealD;

namespace SANS
{

//=============================================================================
// Allows to use promote Surreal for any number of template arguments
//
// typedef Real Arg0;
// typedef SurrealS<1> Arg1;
// typedef Real Arg2;
// typedef typename promote_SurrealN<Arg0,Arg1,Arg2>::type T;
// static_assert( std::is_same< SurrealS<1>, T>::value, "T should be a Surreal" );
//
template<class... Args>
struct promote_Surreal;

//Recursively pick off the Surreal class
template<class L, class R, class... Args>
struct promote_Surreal<L, R, Args...>
{
  typedef typename promote_Surreal< typename promote_Surreal< L, R >::type, Args...>::type type;
};

// Given the choice if two types, choose the surreal type
template<class T1, class T2>
struct promote_Surreal<T1, T2> { typedef typename promote_Surreal< typename Scalar<T1>::type, typename Scalar<T2>::type >::type type; };


template<>
struct promote_Surreal<int,int> { typedef int type; };

template<>
struct promote_Surreal<Real,Real> { typedef Real type; };

template<>
struct promote_Surreal<Real,int> { typedef Real type; };

template<>
struct promote_Surreal<int,Real> { typedef Real type; };

template<>
struct promote_Surreal<Real,std::size_t> { typedef Real type; };

template<>
struct promote_Surreal<std::size_t,Real> { typedef Real type; };


template<int N, class T>
struct promote_Surreal< int, SurrealS<N,T> > { typedef SurrealS<N,T> type; };

template<int N, class T>
struct promote_Surreal< Real, SurrealS<N,T> > { typedef SurrealS<N,T> type; };

template<int N, class T>
struct promote_Surreal< SurrealS<N,T>, int > { typedef SurrealS<N,T> type; };

template<int N, class T>
struct promote_Surreal< SurrealS<N,T>, Real > { typedef SurrealS<N,T> type; };

template<int N, class T>
struct promote_Surreal< SurrealS<N,T>, SurrealS<N,T> > { typedef SurrealS<N,T> type; };

template<>
struct promote_Surreal< int, SurrealD > { typedef SurrealD type; };

template<>
struct promote_Surreal< Real, SurrealD > { typedef SurrealD type; };

template<>
struct promote_Surreal< SurrealD, int > { typedef SurrealD type; };

template<>
struct promote_Surreal< SurrealD, Real > { typedef SurrealD type; };

template<>
struct promote_Surreal< SurrealD, SurrealD > { typedef SurrealD type; };



} // namespace SANS

#endif //PROMOTE_SURREAL_H
