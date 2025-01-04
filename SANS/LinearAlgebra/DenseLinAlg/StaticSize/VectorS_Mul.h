// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef VECTORS_MUL_H
#define VECTORS_MUL_H

#include <type_traits>

// Use boost static assert to show the integers in the compiler error messages.
// C++11 static_assert lacks this ability
#include <boost/mpl/assert.hpp>

// Surreals are registered with boost::is_arithmetic.
// not sure how to register them with std::is_arithmetic
#include <boost/type_traits/is_arithmetic.hpp>

#include "SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS_Type.h"
#include "SANS/LinearAlgebra/DenseLinAlg/tools/PromoteSurreal.h"

namespace SANS
{
namespace DLA
{

// Multiplication with VectorS

//-----------------------------------------------------------------------------
// This specialization is for multiplication between two VectorS, i.e.
// M = A*B
// Where M is VectorS<M, Vector<N, T2>>, A is VectorS<M, T0>, and B is VectorS<N, T1>
// A and B are not commutative
// no temporary variables are required
//
template<class, class >
class OpMulVecS;

template<int ML, class TL, class Expr>
class OpMulVecS< VectorS<ML,TL>, Expr > : public MatrixSType< OpMulVecS< VectorS<ML,TL>, Expr >, true, true >
{
public:
  typedef typename Expr::Ttype TR;
  static const int MR = Expr::M;
  BOOST_MPL_ASSERT_RELATION( Expr::N, ==, 1 );

  typedef typename promote_Surreal<TL,TR>::type Ttype;

  typedef VectorS<ML,TL> VectorSL;

  static const int M = ML;
  static const int N = 1;

  OpMulVecS(const VectorSL& vl, const Expr& er) : vl_(vl), er_(er) {}

  template<class Scalar>
  inline void value(const Scalar& sgn, MatrixS<ML, 1, VectorS<MR, Ttype>>& res) const
  {
    for (int i = 0; i < ML; i++)
      for (int j = 0; j < MR; j++)
        res(i,0)[j] = sgn*vl_[i]*er_.value(j);
  }
  template<class Scalar>
  inline void plus(const Scalar& sgn, MatrixS<ML, 1, VectorS<MR, Ttype>>& res) const
  {
    for (int i = 0; i < ML; i++)
      for (int j = 0; j < MR; j++)
        res(i,0)[j] += sgn*vl_[i]*er_.value(j);
  }

  //Means to access the left and right entries
  const VectorSL& left()  const { return vl_; }
  const Expr& right() const { return er_; }

  inline const OpMulVecS&
  operator+() const { return *this; }
private:
  const VectorSL& vl_;
  const Expr& er_;
};

//-----------------------------------------------------------------------------
// Operator for generating an OpMulS representation of VectorS<ML, VectorS<MR, T>> = VectorS<ML,TL> * VectorS<MR,TR>
template<int ML, class TL, int MR, class TR>
inline typename std::enable_if< boost::is_arithmetic<TL>::value && boost::is_arithmetic<TR>::value &&
                                (ML > 1 || MR > 1),
                                OpMulVecS<VectorS<ML, TL>, VectorS<MR, TR>> >::type
operator*(const VectorS<ML, TL>& vL, const VectorS<MR, TR>& vR)
{
  return OpMulVecS<VectorS<ML, TL>, VectorS<MR, TR>>( vL, vR );
}

//-----------------------------------------------------------------------------
// Operator for generating an OpMulS representation of VectorS<ML, VectorS<MR, T>> = VectorS<ML,TL> * VectorS<MR,TR>
template<int ML, class TL, class Expr>
inline typename std::enable_if< boost::is_arithmetic<TL>::value && boost::is_arithmetic<typename Expr::Ttype>::value &&
                                (ML > 1 && Expr::M > 1 && Expr::N == 1),
                                OpMulVecS<VectorS<ML, TL>, Expr>
                              >::type
operator*(const VectorS<ML, TL>& vL, const MatrixSType<Expr, false, true>& er)
{
  return OpMulVecS<VectorS<ML, TL>, Expr>( vL, er.cast() );
}


} //namespace DLA
} //namespace SANS


#endif //VECTORS_MUL_H
