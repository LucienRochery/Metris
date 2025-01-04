// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_MUL_H
#define MATRIXS_MUL_H

#include <type_traits>

// Use boost static assert to show the integers in the compiler error messages.
// C++11 static_assert lacks this ability
#include <boost/mpl/assert.hpp>

// Surreals are registered with boost::is_arithmetic.
// not sure how to register them with std::is_arithmetic
#include <boost/type_traits/is_arithmetic.hpp>

#include "SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS_Type.h"
#include "SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS_MatMul_Native.h"
#include "SANS/LinearAlgebra/DenseLinAlg/tools/PromoteSurreal.h"
#include "SANS//Surreal/SurrealS_Type.h"

namespace SANS
{
namespace DLA
{

//Multiplication with MatrixS

//-----------------------------------------------------------------------------
// OpMulS_impl implements the actual multiplication algorithm
template<class TL, class TR, class S, class T>
class OpMulS_impl;


template<int ML, int K, class TL, int NR, class TR, class S, class T>
class OpMulS_impl<MatrixS<ML,K,TL>, MatrixS<K,NR,TR>, S, MatrixS<ML,NR,T>>
{
public:

  static inline void value(const MatrixS<ML,K,TL>& ml, const MatrixS<K,NR,TR>& mr, const S& sgn, MatrixS<ML,NR,T>& res )
  {
    MatrixS_MatMul_Native<TL,TR,S,T>::value(ml, mr, sgn, res);
  }

  static inline void plus(const MatrixS<ML,K,TL>& ml, const MatrixS<K,NR,TR>& mr, const S& sgn, MatrixS<ML,NR,T>& res )
  {
    MatrixS_MatMul_Native<TL,TR,S,T>::plus(ml, mr, sgn, res);
  }
};

template<int K, class TL, class TR, class S, class T>
class OpMulS_impl<MatrixS<1,K,TL>, MatrixS<K,1,TR>, S, MatrixS<1,1,T>>
{
public:
  static inline void value(const MatrixS<1,K,TL>& ml, const MatrixS<K,1,TR>& mr, const S& sgn, MatrixS<1,1,T>& res )
  {
    res = 0;
    plus(ml,mr,sgn,res);
  }

  static inline void plus(const MatrixS<1,K,TL>& ml, const MatrixS<K,1,TR>& mr, const S& sgn, MatrixS<1,1,T>& res )
  {
    for (int k = 0; k < K; k++)
      res(0,0) += sgn*(ml(0,k)*mr(k,0));
  }
};

template<int K, class TL, class TR, class S, class T>
class OpMulS_impl<MatrixS<1,K,TL>, MatrixS<K,1,TR>, S, T>
{
public:
  static inline void value(const MatrixS<1,K,TL>& ml, const MatrixS<K,1,TR>& mr, const S& sgn, T& res )
  {
    res = 0;
    plus(ml,mr,sgn,res);
  }

  static inline void plus(const MatrixS<1,K,TL>& ml, const MatrixS<K,1,TR>& mr, const S& sgn, T& res )
  {
    for (int k = 0; k < K; k++)
      res += sgn*(ml(0,k)*mr(k,0));
  }
};



#if 0
//-----------------------------------------------------------------------------
template<int M, int N>
class OpMulS_impl<Real>
{
public:
  template<int ML, int NL, int MR, int NR>
  static inline void value(const MatrixS<ML,NL,T>& ml, const MatrixS<MR,NR,T>& mr, const T& sgn, MatrixS<ML,NR,T>& res )
  {
    if (M*N > 1) //Could add a lower limit for the size to call BLAS here
      MatMul_BLAS<Real>::value(ml, mr, sgn, res);
    else
      MatMul_Native<Real>::value(ml, mr, sgn, res);
  }

  template<int ML, int NL, int MR, int NR>
  static inline void plus(const MatrixS<ML,NL,T>& ml, const MatrixS<MR,NR,T>& mr, const T& sgn, MatrixS<ML,NR,T>& res )
  {
    if (M*N > 1) //Could add a lower limit for the size to call BLAS here
      MatMul_BLAS<Real>::plus(ml, mr, sgn, res);
    else
      MatMul_Native<Real>::plus(ml, mr, sgn, res);
  }
};
#endif

template<class T>
struct dummy_scalar_conversion { explicit dummy_scalar_conversion(T&&) {} };

//-----------------------------------------------------------------------------
// This is a general multiplication of two expressions, i.e.
// M = (A + B)*(C + D)
// The expressions must be evaluated into temporary variables in order to evaluate
// the multiplication
template<class ExprL, class ExprR>
class OpMulS : public MatrixSType< OpMulS<ExprL, ExprR>, true, true >
{
public:
  typedef typename ExprL::Ttype TL;
  typedef typename ExprR::Ttype TR;
  typedef typename promote_Surreal<TL,TR>::type Ttype;

  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::M );

  typedef MatrixS<ExprL::M, ExprL::N, TL> MatrixSL;
  typedef MatrixS<ExprR::M, ExprR::N, TR> MatrixSR;

  static const int M = ExprL::M;
  static const int N = ExprR::N;

  OpMulS(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR) {}

  template<class Scalar, class T>
  inline void value(const Scalar& sgn, T& res) const
  {
    MatrixSL ml( eL );
    MatrixSR mr( eR );

    OpMulS_impl<MatrixSL,MatrixSR,Scalar,T>::value( ml, mr, sgn, res );
  }
  template<class Scalar, class T>
  inline void plus(const Scalar& sgn, T& res) const
  {
    MatrixSL ml( eL );
    MatrixSR mr( eR );

    OpMulS_impl<MatrixSL,MatrixSR,Scalar,T>::plus( ml, mr, sgn, res );
  }
/*
  template< class MatrixL >
  inline void value(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<T> ML( eL );
    MatrixS<T> MR( eR );
    MatrixS<T> res( tuple.m(), tuple.n() );

    OpMulS_impl<T>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixL >
  inline void plus(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<T> ML( eL );
    MatrixS<T> MR( eR );
    MatrixS<T> res( tuple );

    OpMulS_impl<T>::plus( ML, MR, sgn, res );
    tuple = res;
  }
*/

  typedef typename std::conditional<M==1 && N==1, Ttype, dummy_scalar_conversion<Ttype>>::type scalar_type;
  operator scalar_type() const
  {
    MatrixS<1,1,Ttype> A(*this);
    return A(0,0);
  }

  //Means to access the left and right entries
  const ExprL& left()  const { return eL; }
  const ExprR& right() const { return eR; }

  inline const OpMulS&
  operator+() const { return *this; }
private:
  const ExprL& eL;
  const ExprR& eR;
};

//-----------------------------------------------------------------------------
// This specialization is for multiplication between a matrix and an expression, i.e.
// M = A*(B + C)
// only the expression needs to be evaluated in a temporary variable
template<int M_, int N_, class TL, class Expr>
class OpMulS< MatrixS<M_,N_,TL>, Expr > : public MatrixSType< OpMulS< MatrixS<M_,N_,TL>, Expr >, true, true >
{
public:
  typedef typename Expr::Ttype TR;
  typedef typename promote_Surreal<TL,TR>::type Ttype;

  BOOST_MPL_ASSERT_RELATION( N_, ==, Expr::M );

  typedef MatrixS<M_     , N_     , TL> MatrixSL;
  typedef MatrixS<Expr::M, Expr::N, TR> MatrixSR;

  static const int M = M_;
  static const int N = Expr::N;

  OpMulS(const MatrixSL& ml, const Expr& e) :  ml(ml), e(e) {}

  template<class Scalar, class T>
  inline void value(const Scalar& sgn, T& res) const
  {
    MatrixSR mr( e );
    OpMulS_impl<MatrixSL,MatrixSR,Scalar,T>::value( ml, mr, sgn, res );
  }
  template<class Scalar, class T>
  inline void plus(const Scalar& sgn, T& res) const
  {
    MatrixSR mr( e );
    OpMulS_impl<MatrixSL,MatrixSR,Scalar,T>::plus( ml, mr, sgn, res );
  }
/*
  template< class MatrixL >
  inline void value(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<T> MR( e );
    MatrixS<T> res( tuple.m(), tuple.n() );

    OpMulS_impl<T>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixL >
  inline void plus(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<T> MR( e );
    MatrixS<T> res( tuple );

    OpMulS_impl<T>::plus( ML, MR, sgn, res );
    tuple = res;
  }
*/

  typedef typename std::conditional<M==1 && N==1, Ttype, dummy_scalar_conversion<Ttype>>::type scalar_type;
  operator scalar_type() const
  {
    MatrixS<1,1,Ttype> A(*this);
    return A(0,0);
  }

  //Means to access the left and right entries
  const MatrixSL& left()  const { return ml; }
  const Expr&     right() const { return e; }

  inline const OpMulS&
  operator+() const { return *this; }
private:
  const MatrixSL& ml;
  const Expr& e;
};

//-----------------------------------------------------------------------------
// This specialization is for multiplication between an expression and a matrix, i.e.
// M = (A + B)*C
// only the expression needs to be evaluated in a temporary variable
template<class Expr, int M_, int N_, class TR>
class OpMulS< Expr, MatrixS<M_,N_,TR> > : public MatrixSType< OpMulS< Expr, MatrixS<M_,N_,TR> >, true, true >
{
public:
  typedef typename Expr::Ttype TL;
  typedef typename promote_Surreal<TL,TR>::type Ttype;

  BOOST_MPL_ASSERT_RELATION( Expr::N, ==, M_ );

  typedef MatrixS<Expr::M, Expr::N, TL> MatrixSL;
  typedef MatrixS<M_     , N_     , TR> MatrixSR;

  static const int M = Expr::M;
  static const int N = N_;

  OpMulS(const Expr& e, const MatrixSR& mr) : e(e), mr(mr) {}

  template<class Scalar, class T>
  inline void value(const Scalar& sgn, T& res) const
  {
    MatrixSL ml( e );
    OpMulS_impl<MatrixSL,MatrixSR,Scalar,T>::value( ml, mr, sgn, res );
  }
  template<class Scalar, class T>
  inline void plus(const Scalar& sgn, T& res) const
  {
    MatrixSL ml( e );
    OpMulS_impl<MatrixSL,MatrixSR,Scalar,T>::plus( ml, mr, sgn, res );
  }
/*
  template< class MatrixL >
  inline void value(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<T> ML( e );
    MatrixS<T> res( tuple.m(), tuple.n() );

    OpMulS_impl<T>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixL >
  inline void plus(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<T> ML( e );
    MatrixS<T> res( tuple );

    OpMulS_impl<T>::plus( ML, MR, sgn, res );
    tuple = res;
  }
*/

  typedef typename std::conditional<M==1 && N==1, Ttype, dummy_scalar_conversion<Ttype>>::type scalar_type;
  operator scalar_type() const
  {
    MatrixS<1,1,Ttype> A(*this);
    return A(0,0);
  }

  //Means to access the left and right entries
  const Expr&     left()  const { return e; }
  const MatrixSR& right() const { return mr; }

  inline const OpMulS&
  operator+() const { return *this; }
private:
  const Expr& e;
  const MatrixSR& mr;
};

//-----------------------------------------------------------------------------
// This specialization is for multiplication between two matrices, i.e.
// M = A*B
// no temporary variables are required
template<int ML, int NL, class TL, int MR, int NR, class TR>
class OpMulS< MatrixS<ML,NL,TL>, MatrixS<MR,NR,TR> > : public MatrixSType< OpMulS< MatrixS<ML,NL,TL>, MatrixS<MR,NR,TR> >, true, true >
{
public:
  typedef typename promote_Surreal<TL,TR>::type Ttype;
  BOOST_MPL_ASSERT_RELATION( NL, ==, MR );

  typedef MatrixS<ML,NL,TL> MatrixSL;
  typedef MatrixS<MR,NR,TR> MatrixSR;

  static const int M = ML;
  static const int N = NR;

  OpMulS(const MatrixSL& ml, const MatrixSR& mr) : ml(ml), mr(mr) {}

  template<class Scalar, class T>
  inline void value(const Scalar& sgn, T& res) const
  {
    OpMulS_impl<MatrixSL,MatrixSR,Scalar,T>::value( ml, mr, sgn, res );
  }
  template<class Scalar, class T>
  inline void plus(const Scalar& sgn, T& res) const
  {
    OpMulS_impl<MatrixSL,MatrixSR,Scalar,T>::plus( ml, mr, sgn, res );
  }
/*
  //A temporary variable cannot be avoided when res is a MatrixTuple
  template< class MatrixL >
  inline void value(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<M, N, T> res( tuple.m(), tuple.n() );
    OpMulS_impl<M, N, T>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixL >
  inline void plus(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<M, N, T> res( tuple );
    OpMulS_impl<M, N, T>::plus( ML, MR, sgn, res );
    tuple = res;
  }
*/

  typedef typename std::conditional<M==1 && N==1, Ttype, dummy_scalar_conversion<Ttype>>::type scalar_type;
  operator scalar_type() const
  {
    MatrixS<1,1,Ttype> A(*this);
    return A(0,0);
  }

  //Means to access the left and right entries
  const MatrixSL& left()  const { return ml; }
  const MatrixSR& right() const { return mr; }

  inline const OpMulS&
  operator+() const { return *this; }
private:
  const MatrixSL& ml;
  const MatrixSR& mr;
};

//-----------------------------------------------------------------------------
// This specialization is for multiplication between two 'scalar' symmetric matrices, i.e.
// M = A*B
// no temporary variables are required
template<class TL, class TR>
class OpMulS< MatrixSymS<1,TL>, MatrixSymS<1,TR> > : public MatrixSType< OpMulS< MatrixSymS<1,TL>, MatrixSymS<1,TR> >, true, true >
{
public:
  typedef typename promote_Surreal<TL,TR>::type Ttype;

  static const int M = 1;
  static const int N = 1;

  OpMulS(const MatrixSymS<1,TL>& ml, const MatrixSymS<1,TR>& mr) : ml(ml), mr(mr) {}

  template<class Scalar, class T>
  inline void value(const Scalar& sgn, MatrixS<1,1,T>& res) const
  {
    res(0,0) = sgn*ml(0,0)*mr(0,0);
  }
  template<class Scalar, class T>
  inline void plus(const Scalar& sgn, MatrixS<1,1,T>& res) const
  {
    res(0,0) += sgn*ml(0,0)*mr(0,0);
  }
/*
  //A temporary variable cannot be avoided when res is a MatrixTuple
  template< class MatrixL >
  inline void value(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<M, N, T> res( tuple.m(), tuple.n() );
    OpMulS_impl<M, N, T>::value( ML, MR, sgn, res );
    tuple = res;
  }
  template< class MatrixL >
  inline void plus(const T& sgn, MatrixSTuple<MatrixL>& tuple) const
  {
    MatrixS<M, N, T> res( tuple );
    OpMulS_impl<M, N, T>::plus( ML, MR, sgn, res );
    tuple = res;
  }
*/

  typedef typename std::conditional<M==1 && N==1, Ttype, dummy_scalar_conversion<Ttype>>::type scalar_type;
  operator scalar_type() const
  {
    return ml(0,0)*mr(0,0);
  }

  //Means to access the left and right entries
  const MatrixSymS<1,TL>& left()  const { return ml; }
  const MatrixSymS<1,TR>& right() const { return mr; }

  inline const OpMulS&
  operator+() const { return *this; }
private:
  const MatrixSymS<1,TL>& ml;
  const MatrixSymS<1,TR>& mr;
};

//-----------------------------------------------------------------------------
// Operator for generating an OpMulS representation of matrix multiplication including expressions
template<class ExprL, bool useRFL, bool MatrixFullL, class ExprR, bool useRFR, bool MatrixFullR>
inline OpMulS<ExprL, ExprR>
operator*(const MatrixSType<ExprL, useRFL, MatrixFullL>& eL, const MatrixSType<ExprR, useRFR, MatrixFullR>& eR)
{
  return OpMulS<ExprL, ExprR>( eL.cast(), eR.cast() );
}


// Multiplication between a MatrixS expression and a scalar quantity


//-----------------------------------------------------------------------------
// This represents multiplication between a matrix or matrix expression and a scalar, i.e.
// M = s*(A + B)
// or
// M = (A + B)*s
// where s is not POD, but is a scalar quantity. The class also represents division with s by storing 1/s.
template<class Expr, class S, bool useRF, bool MatrixFull>
class OpMulSScalar : public MatrixSType< OpMulSScalar<Expr, S, useRF, MatrixFull>, useRF, MatrixFull >
{
public:
  typedef typename promote_Surreal< S, typename Expr::Ttype >::type Ttype;
  static const int M = Expr::M;
  static const int N = Expr::N;

  OpMulSScalar(const Expr& e, const S s) : e(e), s(s) {}

  template< class Scalar, class T>
  inline void value(const Scalar& sgn, T& res) const
  {
    e.value((S)(s*sgn), res);
  }

  template<class Scalar, class T>
  inline void plus(const Scalar& sgn, T& res) const
  {
    e.plus((S)(s*sgn), res);
  }

  //Element-wise expression
  inline Ttype value(const int& i) const
  {
    return s*e.value(i);
  }

  typedef typename std::conditional<M==1 && N==1, Ttype, dummy_scalar_conversion<Ttype>>::type scalar_type;
  operator scalar_type() const
  {
#if defined(__GNUC__) && !defined(__clang__)
    // suppress -Wmaybe-uninitialized
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
    MatrixS<1,1,Ttype> A(*this);
    return A(0,0);
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
  }

  inline const OpMulSScalar&
  operator+() const { return *this; }

  const Expr& e;
  const S s;
};


template<class T> struct Real_or_T { typedef T type; };
template<> struct Real_or_T<int> { typedef Real type; };

//=============================================================================
// Overloaded operators to represent multiplication between a scalar and a matrix expression
template<class Expr, bool useRF, bool MatrixFull, typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, OpMulSScalar<Expr, typename Real_or_T<T>::type, useRF, MatrixFull > >::type
operator*(const MatrixSType<Expr, useRF, MatrixFull>& e, const T& s)
{
  return OpMulSScalar<Expr, typename Real_or_T<T>::type, useRF, MatrixFull >( e.cast(), s );
}

template<class Expr, bool useRF, bool MatrixFull, typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, OpMulSScalar<Expr, typename Real_or_T<T>::type, useRF, MatrixFull > >::type
operator*(const T& s, const MatrixSType<Expr, useRF, MatrixFull>& e)
{
  return OpMulSScalar<Expr, typename Real_or_T<T>::type, useRF, MatrixFull>( e.cast(), s );
}

template<class Expr, bool useRF, bool MatrixFull, typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, OpMulSScalar<Expr, typename Real_or_T<T>::type, useRF, MatrixFull > >::type
operator/(const MatrixSType<Expr, useRF, MatrixFull>& e, const T& s)
{
  return OpMulSScalar<Expr, typename Real_or_T<T>::type, useRF, MatrixFull>( e.cast(), Real(1)/s );
}

#ifdef SURREAL_LAZY

//=============================================================================
// Overloaded operators specific to working with SurrealS
template<class ExprL, class ExprR, bool useRF, bool MatrixFull, typename T>
inline OpMulSScalar<ExprL, SurrealS<ExprR::N,T>, useRF, MatrixFull >
operator*(const MatrixSType<ExprL, useRF, MatrixFull>& e, const SurrealS<ExprR::N,T>& s)
{
  return OpMulSScalar<ExprL, SurrealS<ExprR::N,T>, useRF, MatrixFull >( e.cast(), s );
}

template<class ExprL, class ExprR, bool useRF, bool MatrixFull, typename T>
inline OpMulSScalar<ExprR, SurrealS<ExprL::N,T>, useRF, MatrixFull >
operator*(const SurrealSType<ExprL, T>& s, const MatrixSType<ExprR, useRF, MatrixFull>& e)
{
  return OpMulSScalar<ExprR, SurrealS<ExprL::N,T>, useRF, MatrixFull>( e.cast(), s );
}

#endif

//=============================================================================
// This is a special case when multiplies scalars are multiplying from two sides, i.e. B = 2*A*2;
// This reduces the complexity of the expression tree and hence reduces code bloat
template<class Expr, typename T, bool useRF, bool MatrixFull, class S>
inline typename std::enable_if< boost::is_arithmetic<S>::value, OpMulSScalar<Expr, T, useRF, MatrixFull> >::type
operator*(const OpMulSScalar<Expr, T, useRF, MatrixFull>& MulScal, const S& s)
{
  return OpMulSScalar<Expr, T, useRF, MatrixFull>( MulScal.e, MulScal.s*s );
}

template<class Expr, typename T, bool useRF, bool MatrixFull, class S>
inline typename std::enable_if< boost::is_arithmetic<S>::value, OpMulSScalar<Expr, T, useRF, MatrixFull> >::type
operator/(const OpMulSScalar<Expr, T, useRF, MatrixFull>& MulScal, const S& s)
{
  return OpMulSScalar<Expr, T, useRF, MatrixFull>( MulScal.e, MulScal.s/s );
}

template<class Expr, typename T, bool useRF, bool MatrixFull, class S>
inline typename std::enable_if< boost::is_arithmetic<S>::value, OpMulSScalar<Expr, T, useRF, MatrixFull> >::type
operator*(const S& s, const OpMulSScalar<Expr, T, useRF, MatrixFull>& MulScal)
{
  return OpMulSScalar<Expr, T, useRF, MatrixFull>( MulScal.e, MulScal.s*s );
}


//-----------------------------------------------------------------------------
// This specialization is for multiplication that involve a scalar and translates
// (s*A)*x into s*(A*x)
// no temporary variables are required
template<class ExprL, class S, bool useRFL, bool MatrixFullL, class ExprR, bool MatrixFull>
class OpMulSFactor< OpMulSScalar<ExprL, S, useRFL, MatrixFullL>, ExprR, MatrixFull >
  : public MatrixSType< OpMulSFactor< OpMulSScalar<ExprL, S, useRFL, MatrixFullL>, ExprR, MatrixFull >, true, MatrixFull >
{
public:
  typedef typename ExprR::Ttype Ttype;
  static const int M = ExprL::M;
  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::M );
  static const int N = ExprR::N;

  OpMulSFactor(const OpMulSScalar<ExprL,S,useRFL,MatrixFullL>& eL, const ExprR& eR) : eL(eL), eR(eR) {}

  template< class Scalar, class T>
  inline void value(const Scalar& sgn, T& res) const
  {
    S tmp = eL.s*sgn;
    res = tmp*(eL.e*eR);
  }
  template< class Scalar, class T>
  inline void plus(const Scalar& sgn, T& res) const
  {
    S tmp = eL.s*sgn;
    res += tmp*(eL.e*eR);
  }

  typedef typename std::conditional<M==1 && N==1, Ttype, dummy_scalar_conversion<Ttype>>::type scalar_type;
  operator scalar_type() const
  {
    MatrixS<1,1,Ttype> A(*this);
    return A(0,0);
  }

  inline const OpMulSFactor&
  operator+() const { return *this; }
private:
  const OpMulSScalar<ExprL, S, useRFL, MatrixFullL>& eL;
  const ExprR& eR;
};

//-----------------------------------------------------------------------------
// This specialization is for multiplication that involve a scalar and translates
// A*(x*s) into (A*x)*s
// no temporary variables are required
template<class ExprL, class ExprR, class S, bool useRFR, bool MatrixFullR, bool MatrixFull>
class OpMulSFactor< ExprL, OpMulSScalar<ExprR, S, useRFR, MatrixFullR>, MatrixFull >
  : public MatrixSType< OpMulSFactor< ExprL, OpMulSScalar<ExprR, S, useRFR, MatrixFullR>, MatrixFull >, true, MatrixFull >
{
public:
  typedef typename ExprL::Ttype Ttype;
  static const int M = ExprL::M;
  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::M );
  static const int N = ExprR::N;

  OpMulSFactor(const ExprL& eL, const OpMulSScalar<ExprR, S, useRFR, MatrixFullR>& eR) : eL(eL), eR(eR) {}

  template<class Scalar, class T>
  inline void value(const Scalar& sgn, T& res) const
  {
    S tmp = eR.s*sgn;
    res = tmp*(eL*eR.e);
  }
  template<class Scalar, class T>
  inline void plus(const Scalar& sgn, T& res) const
  {
    S tmp = eR.s*sgn;
    res += tmp*(eL*eR.e);
  }

  typedef typename std::conditional<M==1 && N==1, Ttype, dummy_scalar_conversion<Ttype>>::type scalar_type;
  operator scalar_type() const
  {
    MatrixS<1,1,Ttype> A(*this);
    return A(0,0);
  }

  inline const OpMulSFactor&
  operator+() const { return *this; }
private:
  const ExprL& eL;
  const OpMulSScalar<ExprR, S, useRFR, MatrixFullR>& eR;
};


//-----------------------------------------------------------------------------
// Operator for generating an OpMulS representation of matrix multiplication including expressions
template<class ExprL, class S, bool useRFL, bool MatrixFullL, class ExprR, bool useRFR, bool MatrixFullR>
inline OpMulSFactor< OpMulSScalar<ExprL, S, useRFL, MatrixFullL>, ExprR, MatrixFullL || MatrixFullR>
operator*(const OpMulSScalar<ExprL, S, useRFL, MatrixFullL>& eL, const MatrixSType<ExprR, useRFR, MatrixFullR>& eR)
{
  return OpMulSFactor< OpMulSScalar<ExprL, S, useRFL, MatrixFullL>, ExprR, MatrixFullL || MatrixFullR>( eL, eR.cast() );
}

template<class ExprL, bool useRFL, bool MatrixFullL, class ExprR, class S, bool useRFR, bool MatrixFullR>
inline OpMulSFactor<ExprL, OpMulSScalar<ExprR, S, useRFR, MatrixFullR>, MatrixFullL || MatrixFullR >
operator*(const MatrixSType<ExprL, useRFL, MatrixFullL>& eL, const OpMulSScalar<ExprR, S, useRFR, MatrixFullR>& eR)
{
  return OpMulSFactor<ExprL, OpMulSScalar<ExprR, S, useRFR, MatrixFullR>, MatrixFullL || MatrixFullR >( eL.cast(), eR );
}


} //namespace DLA
} //namespace SANS


#endif //MATRIXS_MUL_H
