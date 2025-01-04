// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_ADD_H
#define MATRIXS_ADD_H

// Use boost static assert to show the integers in the compiler error messages.
// C++11 static_assert lacks this ability
#include <boost/mpl/assert.hpp>

#include "SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS_Type.h"
#include "SANS/LinearAlgebra/DenseLinAlg/tools/PromoteSurreal.h"

namespace SANS
{
namespace DLA
{

// Lazy expressions based on recursive function calls

// Addition and Subtraction
// No temporary variables are needed for these operations

template<class ExprL, class ExprR, bool useRF, bool MatrixFull>
class OpAddS : public MatrixSType< OpAddS<ExprL, ExprR, useRF, MatrixFull>, useRF, MatrixFull >
{
public:
  typedef typename ExprL::Ttype TL;
  typedef typename ExprR::Ttype TR;
  typedef typename promote_Surreal<TL,TR>::type Ttype;
  BOOST_MPL_ASSERT_RELATION( ExprL::M, ==, ExprR::M );
  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::N );
  static const int M = ExprL::M;
  static const int N = ExprL::N;

  OpAddS(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR) {}

  template< class Scalar, class T >
  inline void value(const Scalar& sgn, MatrixS<M, N, T>& res) const
  {
    eL.value(sgn, res);
    eR.plus(sgn, res);
  }
  template< class Scalar, class T >
  inline void plus(const Scalar& sgn, MatrixS<M, N, T>& res) const
  {
    eL.plus(sgn, res);
    eR.plus(sgn, res);
  }

  //Element-wise expression
  inline Ttype value(const int& i) const
  {
    return eL.value(i) + eR.value(i);
  }

/*
  template< class MatrixL >
  inline void value(const T& sgn, MatrixSTuple<MatrixL>& res) const
  {
    eL.value(sgn, res);
    eR.plus(sgn, res);
  }
  template< class MatrixL >
  inline void plus(const T& sgn, MatrixSTuple<MatrixL>& res) const
  {
    eL.plus(sgn, res);
    eR.plus(sgn, res);
  }
*/
  inline const OpAddS&
  operator+() const { return *this; }
private:
  const ExprL& eL;
  const ExprR& eR;
};


template<class ExprL, class ExprR, bool useRF, bool MatrixFull >
class OpSubS : public MatrixSType< OpSubS<ExprL, ExprR, useRF, MatrixFull>, useRF, MatrixFull >
{
public:
  typedef typename ExprL::Ttype TL;
  typedef typename ExprR::Ttype TR;
  typedef typename promote_Surreal<TL,TR>::type Ttype;
  BOOST_MPL_ASSERT_RELATION( ExprL::M, ==, ExprR::M );
  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::N );
  static const int M = ExprL::M;
  static const int N = ExprL::N;

  OpSubS(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR) {}

  template< class Scalar, class T >
  inline void value(const Scalar& sgn, MatrixS<M, N, T>& res) const
  {
    Scalar nsgn = -sgn;
    eL.value(sgn, res);
    eR.plus(nsgn, res);
  }
  template< class Scalar, class T >
  inline void plus(const Scalar& sgn, MatrixS<M, N, T>& res) const
  {
    Scalar nsgn = -sgn;
    eL.plus(sgn, res);
    eR.plus(nsgn, res);
  }

  //Element-wise expression
  inline Ttype value(const int& i) const
  {
    return eL.value(i) - eR.value(i);
  }

/*
  template< class MatrixL >
  inline void value(const T& sgn, MatrixSTuple<MatrixL>& res) const
  {
    eL.value(sgn, res);
    eR.plus(-sgn, res);
  }
  template< class MatrixL >
  inline void plus(const T& sgn, MatrixSTuple<MatrixL>& res) const
  {
    eL.plus(sgn, res);
    eR.plus(-sgn, res);
  }
*/
  inline const OpSubS&
  operator+() const { return *this; }
private:
  const ExprL& eL;
  const ExprR& eR;
};


//Generator functions to generate an addition/subtraction operation
template<class ExprL, bool useRFL, bool MatrixFullL, class ExprR, bool useRFR, bool MatrixFullR>
inline OpAddS<ExprL, ExprR, useRFL || useRFR, MatrixFullL || MatrixFullR>
operator+(const MatrixSType<ExprL, useRFL, MatrixFullL>& L, const MatrixSType<ExprR, useRFR, MatrixFullR>& R)
{
  return OpAddS<ExprL, ExprR, useRFL || useRFR, MatrixFullL || MatrixFullR>( L.cast(), R.cast() );
}

template<class ExprL, bool useRFL, bool MatrixFullL, class ExprR, bool useRFR, bool MatrixFullR>
inline OpSubS<ExprL, ExprR, useRFL || useRFR, MatrixFullL || MatrixFullR>
operator-(const MatrixSType<ExprL, useRFL, MatrixFullL>& L, const MatrixSType<ExprR, useRFR, MatrixFullR>& R)
{
  return OpSubS<ExprL, ExprR, useRFL || useRFR, MatrixFullL || MatrixFullR>( L.cast(), R.cast() );
}

} //namespace DLA
} //namespace SANS


/* This is dissabled until we decide that it's a good idea to add/subtract scalar quantities with a matrix
namespace SANS
{
namespace DLA
{

//Addition and Subtraction with scalar quantities

template<class Expr>
class OpScalar : public MatrixSType< OpScalar<Expr> >
{
public:
  typedef typename Expr::node_type node_type;

  OpScalar(const Expr& e, const Real esgn, const Real s) : e(e), esgn(esgn), s(s) {}

  template< class T >
  inline void value(const T& sgn, MatrixSView<T>& res) const
  {
    e.value(esgn*sgn, res);
    res.value() += sgn*s;
  }
  template< class T >
  inline void plus(const T& sgn, MatrixSView<T>& res) const
  {
    e.plus(esgn*sgn, res);
    res.value() += sgn*s;
  }

  inline const OpScalar&
  operator+() const { return *this; }
  int m() const { return e.m(); }
  int n() const { return e.n(); }
  int size() const { return m()*n(); }
private:
  const Expr& e;
  const Real esgn;
  const Real s;
};

}
}


template<class Expr>
inline SurrealExp::OpScalar<Expr>
operator+(const MatrixSType<Expr>& e, const Real& s)
{
  return SurrealExp::OpScalar<Expr>( e.cast(), 1, s );
}

template<class Expr>
inline SurrealExp::OpScalar<Expr>
operator+(const Real& s, const MatrixSType<Expr>& e)
{
  return SurrealExp::OpScalar<Expr>( e.cast(), 1, s );
}

template<class Expr>
inline SurrealExp::OpScalar<Expr>
operator-(const MatrixSType<Expr>& e, const Real& s)
{
  return SurrealExp::OpScalar<Expr>( e.cast(), 1, -s );
}

template<class Expr>
inline SurrealExp::OpScalar<Expr>
operator-(const Real& s, const MatrixSType<Expr>& e)
{
  return SurrealExp::OpScalar<Expr>( e.cast(), -1, s );
}
*/

#endif //MATRIXS_ADD_H
