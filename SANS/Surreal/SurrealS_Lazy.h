// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SANS_SURREALS_LAZY_H
#define SANS_SURREALS_LAZY_H

//  overloaded derivative operator
//  ref: derivify.h (Google it)

#include <cmath>
#include <algorithm> // min/max
#include <iostream>
#include <string>

#include <type_traits>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/preprocessor/cat.hpp>

#include "SANS/tools/SANSnumerics.h"     // Real
#include "SANS/tools/SANSException.h"
#include "SANS/tools/SANSTraitsPOD.h"
#include "SANS/tools/CacheLineSize.h"
#include "SANS/tools/AlignMem.h"


#include "SANS//Surreal/SurrealS_Type.h"

//#define SURREALS_LOOP_UNROLL

//Some compilers define min/max macros...
#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

namespace SANS
{

namespace SurrealSExpr
{
template<class L, class R, class T > class OpMul;
class OpMul_impl;
}


//----------------------------------------------------------------------------//
// SurrealS:  value, N derivatives
//
// Operators with Lazy Expressions
//
// statically defined derivative array
//----------------------------------------------------------------------------//
template<int N_, class T>
class SurrealS : public SurrealSType< SurrealS<N_, T>, T >
{
public:
  static const int N = N_;

  //The default constructor is intentionally empty here. This means Surreal is
  //not initialized when declared, which is consistent with regular numbers. This also
  //improves performance.
  ALWAYS_INLINE SurrealS() {}
  ALWAYS_INLINE SurrealS( const SurrealS& z );
  ALWAYS_INLINE SurrealS( const int v0 );
  ALWAYS_INLINE SurrealS( const T& v0 );
  ALWAYS_INLINE SurrealS( const typename SANS::POD<T>::type& v0 );
  ALWAYS_INLINE SurrealS( const Real v0, const int d0[], int n );
  ALWAYS_INLINE SurrealS( const Real v0, const Real d0[], int n );
  template<class Expr>
  ALWAYS_INLINE SurrealS( const SurrealSType<Expr, T>& r ) : v_(0) { operator=(r); }
  ALWAYS_INLINE ~SurrealS() {}

  ALWAYS_INLINE int size() const { return N; }

  // value accessor operators
  ALWAYS_INLINE       T& value()       { return v_; }
  ALWAYS_INLINE const T& value() const { return v_; }

  // derivative accessor operators
  ALWAYS_INLINE       T& deriv( const int i=0 )       { return d_[i]; }
  ALWAYS_INLINE const T& deriv( const int i=0 ) const { return d_[i]; }

  // assignment
  ALWAYS_INLINE SurrealS& operator=( const SurrealS& );
  ALWAYS_INLINE SurrealS& operator=( const int& );
  ALWAYS_INLINE SurrealS& operator=( const Real& );

  template<class Expr> ALWAYS_INLINE SurrealS& operator= ( const SurrealSType<Expr, T>& );
  template<class Expr> ALWAYS_INLINE SurrealS& operator+=( const SurrealSType<Expr, T>& );
  template<class Expr> ALWAYS_INLINE SurrealS& operator-=( const SurrealSType<Expr, T>& );

  // unary operators; no side effects
  ALWAYS_INLINE const SurrealS& operator+() const;

  // binary accumulation operators
  ALWAYS_INLINE SurrealS& operator+=( const Real& );
  ALWAYS_INLINE SurrealS& operator-=( const Real& );
  ALWAYS_INLINE SurrealS& operator*=( const SurrealS& );
  ALWAYS_INLINE SurrealS& operator*=( const Real& );
  ALWAYS_INLINE SurrealS& operator/=( const SurrealS& );
  ALWAYS_INLINE SurrealS& operator/=( const Real& );

#if 0
  // classification functions <cmath>
  friend bool isfinite( const SurrealS& );
  friend bool isinf( const SurrealS& );
  friend bool isnan( const SurrealS& );
#endif

  // input/output
  template<int M> friend std::istream& operator>>( std::istream&, SurrealS<M, T>& );

protected:
  ALIGN_MEM T d_[N];   // derivative array
  ALIGN_MEM T v_;      // value
};

//Constructors

template<int N, class T>
ALWAYS_INLINE
SurrealS<N,T>::SurrealS( const SurrealS& z )
{
  v_ = z.v_;
PRAGMA_IVDEP
  for (int i = 0; i < N; i++)
    d_[i] = z.d_[i];
}
template<int N, class T>
ALWAYS_INLINE
SurrealS<N,T>::SurrealS( const int v0 )
{
  v_ = v0;
PRAGMA_IVDEP
  for (int i = 0; i < N; i++)
    d_[i] = 0;
}
template<int N, class T>
ALWAYS_INLINE
SurrealS<N,T>::SurrealS( const T& v0 )
{
  v_ = v0;
PRAGMA_IVDEP
  for (int i = 0; i < N; i++)
    d_[i] = 0;
}
template<int N, class T>
ALWAYS_INLINE
SurrealS<N,T>::SurrealS( const typename SANS::POD<T>::type& v0 )
{
  v_ = v0;
PRAGMA_IVDEP
  for (int i = 0; i < N; i++)
    d_[i] = 0;
}
template<int N, class T>
ALWAYS_INLINE
SurrealS<N,T>::SurrealS( const Real v0, const int d0[], int n )
{
  SANS_ASSERT( n == N );
  v_ = v0;
PRAGMA_IVDEP
  for (int i = 0; i < N; i++)
    d_[i] = d0[i];
}
template<int N, class T>
ALWAYS_INLINE
SurrealS<N,T>::SurrealS( const Real v0, const Real d0[], int n )
{
  SANS_ASSERT( n == N );
  v_ = v0;
PRAGMA_IVDEP
  for (int i = 0; i < N; i++)
    d_[i] = d0[i];
}


namespace SurrealSExpr
{

// Lazy expressions


// Addition and Subtraction

template<class L, class R, class T>
class OpAdd : public SurrealSType< OpAdd<L,R,T>, T >
{
public:
  static const int N = L::N;
  BOOST_MPL_ASSERT_RELATION( L::N, ==, R::N );

  ALWAYS_INLINE
  OpAdd(const L& Ll, const R& Rr) : Ll(Ll), Rr(Rr) {}

  ALWAYS_INLINE T value() const { return Ll.value() + Rr.value(); }
  ALWAYS_INLINE T deriv(const int& i) const { return Ll.deriv(i) + Rr.deriv(i); }

  ALWAYS_INLINE const OpAdd&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return Ll.size(); }
private:
  const L& Ll;
  const R& Rr;
};

}

template<class L, class R, class T>
ALWAYS_INLINE SurrealSExpr::OpAdd<L,R,T>
operator+(const SurrealSType<L,T>& Ll, const SurrealSType<R,T>& Rr)
{
  return SurrealSExpr::OpAdd<L,R,T>( Ll.cast(), Rr.cast() );
}

namespace SurrealSExpr
{

template<class L, class R, class T>
class OpSub : public SurrealSType< OpSub<L,R,T>, T >
{
public:
  static const int N = L::N;
  BOOST_MPL_ASSERT_RELATION( L::N, ==, R::N );

  ALWAYS_INLINE
  OpSub(const L& Ll, const R& Rr) : Ll(Ll), Rr(Rr) {}

  ALWAYS_INLINE T value() const { return Ll.value() - Rr.value(); }
  ALWAYS_INLINE T deriv(const int& i) const { return Ll.deriv(i) - Rr.deriv(i); }

  ALWAYS_INLINE const OpSub&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return Ll.size(); }
private:
  const L& Ll;
  const R& Rr;
};

}

template<class L, class R, class T>
ALWAYS_INLINE SurrealSExpr::OpSub<L,R,T>
operator-(const SurrealSType<L,T>& Ll, const SurrealSType<R,T>& Rr)
{
  return SurrealSExpr::OpSub<L,R,T>( Ll.cast(), Rr.cast() );
}

//Addition and Subtraction with scalar quantities

namespace SurrealSExpr
{

template<class Expr, class T>
class OpScalar : public SurrealSType< OpScalar<Expr, T>, T >
{
public:
  static const int N = Expr::N;

  ALWAYS_INLINE
  OpScalar(const Expr& e, const Real esgn, const Real s) : e(e), esgn(esgn), s(s) {}

  ALWAYS_INLINE T value() const { return esgn*e.value() + s; }
  ALWAYS_INLINE T deriv(const int& i) const { return esgn*e.deriv(i); }

  ALWAYS_INLINE const OpScalar&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }
private:
  const Expr& e;
  const Real esgn;
  const Real s;
};

}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr,T>
operator+(const SurrealSType<Expr,T>& e, const Real& s)
{
  return SurrealSExpr::OpScalar<Expr,T>( e.cast(), 1, s );
}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr,T>
operator+(const Real& s, const SurrealSType<Expr,T>& e)
{
  return SurrealSExpr::OpScalar<Expr,T>( e.cast(), 1, s );
}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr,T>
operator-(const SurrealSType<Expr,T>& e, const Real& s)
{
  return SurrealSExpr::OpScalar<Expr,T>( e.cast(), 1, -s );
}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr,T>
operator-(const Real& s, const SurrealSType<Expr,T>& e)
{
  return SurrealSExpr::OpScalar<Expr,T>( e.cast(), -1, s );
}


//Multiplication with SurrealSs

namespace SurrealSExpr
{

template<class ExprL, class ExprR, class T>
class OpMul : public SurrealSType< OpMul<ExprL, ExprR, T>, T >
{
public:
  static const int N = ExprL::N;
  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::N );

  ALWAYS_INLINE
  OpMul(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), eL_val(eL.value()), eR_val(eR.value()) {}

  ALWAYS_INLINE T value() const { return eL_val*eR_val; }
  ALWAYS_INLINE T deriv(const int& i) const { return eL_val*eR.deriv(i) + eL.deriv(i)*eR_val; }

  ALWAYS_INLINE const OpMul&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const T eL_val, eR_val;
};

template<class ExprL, class T>
class OpMul<ExprL, Real, T> : public SurrealSType< OpMul<ExprL, Real, T>, T >
{
public:
  static const int N = ExprL::N;

  ALWAYS_INLINE
  OpMul(const ExprL& e, const Real s) : e(e), s(s) {}

  ALWAYS_INLINE T value() const { return e.value()*s; }
  ALWAYS_INLINE T deriv(const int& i) const { return e.deriv(i)*s; }

  ALWAYS_INLINE const OpMul&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }

  const ExprL& e;
  const Real s;
};
}

//=============================================================================
template<class ExprL, class ExprR, class T>
ALWAYS_INLINE SurrealSExpr::OpMul<ExprL, ExprR, T>
operator*(const SurrealSType<ExprL, T>& z1, const SurrealSType<ExprR, T>& z2)
{
  return SurrealSExpr::OpMul<ExprL, ExprR, T>( z1.cast(), z2.cast() );
}

template<class Expr, class T, typename S>
ALWAYS_INLINE typename std::enable_if< is_arithmetic_not_SurrealS<S>::value,
                                         SurrealSExpr::OpMul<Expr, Real, T> >::type
operator*(const SurrealSType<Expr, T>& e, const S& s)
{
  return SurrealSExpr::OpMul<Expr, Real, T>( e.cast(), s );
}

template<class Expr, class T, typename S>
ALWAYS_INLINE typename std::enable_if< is_arithmetic_not_SurrealS<S>::value,
                                         SurrealSExpr::OpMul<Expr, Real, T> >::type
operator/(const SurrealSType<Expr, T>& e, const S& s)
{
  return SurrealSExpr::OpMul<Expr, Real, T>( e.cast(), Real(1)/s );
}

template<class Expr, class T, typename S>
ALWAYS_INLINE typename std::enable_if< is_arithmetic_not_SurrealS<S>::value,
                                         SurrealSExpr::OpMul<Expr, Real, T> >::type
operator*(const S& s, const SurrealSType<Expr, T>& e)
{
  return SurrealSExpr::OpMul<Expr, Real, T>( e.cast(), s );
}

//=============================================================================
// This is a special case when multiplies scalars are multiplying from two sides, i.e. B = 2*A*2;
// This reduces the complexity of the expression tree and hence reduces code bloat
template<class Expr, class T, typename S>
ALWAYS_INLINE typename std::enable_if< is_arithmetic_not_SurrealS<S>::value,
                                         SurrealSExpr::OpMul<Expr, Real, T> >::type
operator*(const SurrealSExpr::OpMul<Expr, Real, T>& MulScal, const S& s)
{
  return SurrealSExpr::OpMul<Expr, Real, T>( MulScal.e, MulScal.s*s );
}

template<class Expr, class T, typename S>
ALWAYS_INLINE typename std::enable_if< is_arithmetic_not_SurrealS<S>::value,
                                         SurrealSExpr::OpMul<Expr, Real, T> >::type
operator/(const SurrealSExpr::OpMul<Expr, Real, T>& MulScal, const S& s)
{
  return SurrealSExpr::OpMul<Expr, Real, T>( MulScal.e, MulScal.s/s );
}

template<class Expr, class T, typename S>
ALWAYS_INLINE typename std::enable_if< is_arithmetic_not_SurrealS<S>::value,
                                         SurrealSExpr::OpMul<Expr, Real, T> >::type
operator*(const S& s, const SurrealSExpr::OpMul<Expr, Real, T>& MulScal)
{
  return SurrealSExpr::OpMul<Expr, Real, T>( MulScal.e, MulScal.s*s );
}


//=============================================================================
//Division with SurrealSs

namespace SurrealSExpr
{

template<class ExprL, class ExprR, class T>
class OpDiv : public SurrealSType< OpDiv<ExprL, ExprR, T>, T >
{
public:
  static const int N = ExprL::N;
  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::N );

  ALWAYS_INLINE
  OpDiv(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), eL_val(eL.value()), eR_val(eR.value())
                                          , vali(1/(eR_val*eR_val)) {}

  ALWAYS_INLINE T value() const { return eL_val/eR_val; }
  ALWAYS_INLINE T deriv(const int& i) const { return (eR_val*eL.deriv(i) - eR.deriv(i)*eL_val)*vali; }

  ALWAYS_INLINE const OpDiv&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const T eL_val, eR_val, vali;
};

}

template< class ExprL, class ExprR, class T >
ALWAYS_INLINE SurrealSExpr::OpDiv<ExprL, ExprR, T>
operator/(const SurrealSType<ExprL, T>& eL, const SurrealSType<ExprR, T>& eR)
{
  return SurrealSExpr::OpDiv<ExprL, ExprR, T>( eL.cast(), eR.cast() );
}


namespace SurrealSExpr
{

template<class Expr, class T>
class OpDivScalarNumerator : public SurrealSType< OpDivScalarNumerator<Expr, T>, T >
{
public:
  static const int N = Expr::N;

  ALWAYS_INLINE
  OpDivScalarNumerator(const Expr& e, const Real& s) : e(e), s(s), e_val(e.value()), se_val2i(s/(e_val*e_val)) {}

  ALWAYS_INLINE T value() const { return s/e_val; }
  ALWAYS_INLINE T deriv(const int& i) const { return -se_val2i*e.deriv(i); }

  ALWAYS_INLINE const OpDivScalarNumerator&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }
private:
  const Expr& e;
  const Real s;
  const T e_val, se_val2i;
};

}

template<class Expr, class T, typename S>
ALWAYS_INLINE typename std::enable_if< is_arithmetic_not_SurrealS<S>::value,
                                         SurrealSExpr::OpDivScalarNumerator<Expr,T> >::type
operator/(const S& s, const SurrealSType<Expr,T>& e)
{
  return SurrealSExpr::OpDivScalarNumerator<Expr,T>( e.cast(), s );
}


// assignment

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator=( const SurrealS& z )
{
  //Do nothing if assigning self to self
  if ( &z == this ) return *this;

  v_ = z.v_;
PRAGMA_IVDEP
  for (int i = 0; i < N; i++)
    d_[i] = z.d_[i];
  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator=( const int& r )
{
  v_ = r;
PRAGMA_IVDEP
  for (int i = 0; i < N; i++)
    d_[i] = 0;
  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator=( const Real& r )
{
  v_ = r;
PRAGMA_IVDEP
  for (int i = 0; i < N; i++)
    d_[i] = 0;
  return *this;
}

template<int N, class T>
template< class Expr >
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator=( const SurrealSType<Expr,T>& r )
{
  const Expr& Tree = r.cast();

  BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

PRAGMA_IVDEP
  for (int i = 0; i < N; ++i)
    d_[i] = Tree.deriv(i);

  //Value must be set last as it might be used in the derivative calculation
  v_ = Tree.value();

  return *this;
}

template<int N, class T>
template< class Expr >
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator+=( const SurrealSType<Expr,T>& r )
{
  const Expr& Tree = r.cast();

  BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

PRAGMA_IVDEP
  for (int i = 0; i < N; ++i)
    d_[i] += Tree.deriv(i);

  //Value must be set last as it might be used in the derivative calculation
  v_ += Tree.value();

  return *this;
}

template<int N, class T>
template< class Expr >
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator-=( const SurrealSType<Expr,T>& r )
{
  const Expr& Tree = r.cast();

  BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

PRAGMA_IVDEP
  for (int i = 0; i < N; ++i)
    d_[i] -= Tree.deriv(i);

  //Value must be set last as it might be used in the derivative calculation
  v_ -= Tree.value();

  return *this;
}


// unary operators; no side effects

template<int N, class T>
ALWAYS_INLINE const SurrealS<N,T>&
SurrealS<N,T>::operator+() const
{
  return *this;
}

template< class Expr, class T >
ALWAYS_INLINE const SurrealSExpr::OpMul<Expr, Real, T>
operator-(SurrealSType<Expr,T> const& e)
{
  return SurrealSExpr::OpMul<Expr, Real, T>( e.cast(), -1 );
}

template< class Expr, class T >
ALWAYS_INLINE const SurrealSExpr::OpMul<Expr, Real, T>
operator-(SurrealSExpr::OpMul<Expr, Real, T> const& Mul)
{
  return SurrealSExpr::OpMul<Expr, Real, T>( Mul.e, -1*Mul.s );
}

// binary accumulation operators


template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator+=( const Real& r )
{
  v_ += r;
  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator-=( const Real& r )
{
  v_ -= r;
  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator*=( const SurrealS& z )
{
PRAGMA_IVDEP
  for (int i = 0; i < N; i++)
    d_[i] = v_*z.d_[i] + d_[i]*z.v_;
  v_ *= z.v_;
  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator*=( const Real& r )
{
PRAGMA_IVDEP
  for (int i = 0; i < N; i++)
    d_[i] *= r;
  v_ *= r;
  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator/=( const SurrealS& z)
{
  Real tmp = 1./(z.v_*z.v_);
PRAGMA_IVDEP
  for (int i = 0; i < N; i++)
    d_[i] = (z.v_*d_[i] - v_*z.d_[i])*tmp;
  v_ /= z.v_;
  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator/=( const Real& r )
{
  Real tmp = 1./r;
PRAGMA_IVDEP
  for (int i = 0; i < N; i++)
    d_[i] *= tmp;
  v_ *= tmp;
  return *this;
}

// relational operators

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE bool
operator==( const SurrealSType<ExprL, T>& lhs, const SurrealSType<ExprR, T>& rhs )
{
  return lhs.value() == rhs.value();
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator==( const SurrealSType<Expr, T>& lhs, const Real& rhs )
{
  return lhs.value() == rhs;
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator==( const Real& lhs, const SurrealSType<Expr, T>& rhs )
{
  return lhs == rhs.value();
}

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE bool
operator!=( const SurrealSType<ExprL, T>& lhs, const SurrealSType<ExprR, T>& rhs )
{
  return lhs.value() != rhs.value();
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator!=( const SurrealSType<Expr, T>& lhs, const Real& rhs )
{
  return lhs.value() != rhs;
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator!=( const Real& lhs, const SurrealSType<Expr, T>& rhs )
{
  return lhs != rhs.value();
}

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE bool
operator>( const SurrealSType<ExprL, T>& lhs, const SurrealSType<ExprR, T>& rhs )
{
  return lhs.value() > rhs.value();
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator>( const SurrealSType<Expr, T>& lhs, const Real& rhs )
{
  return lhs.value() > rhs;
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator>( const Real& lhs, const SurrealSType<Expr, T>& rhs )
{
  return lhs > rhs.value();
}

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE bool
operator<( const SurrealSType<ExprL, T>& lhs, const SurrealSType<ExprR, T>& rhs )
{
  return lhs.value() < rhs.value();
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator<( const SurrealSType<Expr, T>& lhs, const Real& rhs )
{
  return lhs.value() < rhs;
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator<( const Real& lhs, const SurrealSType<Expr, T>& rhs )
{
  return lhs < rhs.value();
}

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE bool
operator>=( const SurrealSType<ExprL, T>& lhs, const SurrealSType<ExprR, T>& rhs )
{
  return lhs.value() >= rhs.value();
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator>=( const SurrealSType<Expr, T>& lhs, const Real& rhs )
{
  return lhs.value() >= rhs;
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator>=( const Real& lhs, const SurrealSType<Expr, T>& rhs )
{
  return lhs >= rhs.value();
}

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE bool
operator<=( const SurrealSType<ExprL, T>& lhs, const SurrealSType<ExprR, T>& rhs )
{
  return lhs.value() <= rhs.value();
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator<=( const SurrealSType<Expr, T>& lhs, const Real& rhs )
{
  return lhs.value() <= rhs;
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator<=( const Real& lhs, const SurrealSType<Expr, T>& rhs )
{
  return lhs <= rhs.value();
}


//Functions for SurrealSs
#define SURREALS_FUNC1( NAME, FUNC, DERIV ) \
using ::NAME; \
namespace SurrealSExpr \
{  \
template<class Expr, class T> \
class BOOST_PP_CAT(SurrealS_, NAME) : public SurrealSType< BOOST_PP_CAT(SurrealS_, NAME)<Expr, T>, T > \
{ /*This is for functions when the argument is an expression*/ \
public: \
  static const int N = Expr::N; \
  \
  ALWAYS_INLINE \
  BOOST_PP_CAT(SurrealS_, NAME)(const Expr& e) : e(e), z(e.value()), der(DERIV) {} \
  \
  ALWAYS_INLINE T value() const { return FUNC; } \
  ALWAYS_INLINE T deriv(const int& i) const { return der*e.deriv(i); } \
  \
  ALWAYS_INLINE const BOOST_PP_CAT(SurrealS_, NAME)& \
  operator+() const { return *this; } \
  ALWAYS_INLINE int size() const { return e.size(); } \
private: \
  const Expr& e; \
  const T z, der; \
}; \
} \
template<class Expr, class T> \
ALWAYS_INLINE SurrealSExpr::BOOST_PP_CAT(SurrealS_, NAME)<Expr, T> \
NAME(const SurrealSType<Expr, T>& z) { return SurrealSExpr::BOOST_PP_CAT(SurrealS_, NAME)<Expr, T>( z.cast() ); }


#define SURREALS_FUNC2( NAME, FUNC, DERIV ) \
using ::NAME; \
namespace SurrealSExpr \
{  \
template<class ExprL, class ExprR, class T> \
class BOOST_PP_CAT(SurrealS_, NAME) : public SurrealSType< BOOST_PP_CAT(SurrealS_, NAME)<ExprL, ExprR, T>, T > \
{ /*This is for functions when the argument is an expression*/ \
public: \
  static const int N = ExprL::N; \
  BOOST_MPL_ASSERT_RELATION(ExprL::N, ==, ExprR::N ); \
  \
  ALWAYS_INLINE \
  BOOST_PP_CAT(SurrealS_, NAME)(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), z1(eL.value()), z2(eR.value()), \
                                                                 der(DERIV) {} \
  \
  ALWAYS_INLINE T value() const { return FUNC; } \
  ALWAYS_INLINE T deriv(const int& i) const { return der*(z2*eL.deriv(i) - z1*eR.deriv(i)); } \
  \
  ALWAYS_INLINE const BOOST_PP_CAT(SurrealS_, NAME)& \
  operator+() const { return *this; } \
  ALWAYS_INLINE int size() const { return eL.size(); } \
private: \
  const ExprL& eL; \
  const ExprR& eR; \
  const T z1, z2, der; \
}; \
  \
} \
template<class ExprL, class ExprR, class T> \
ALWAYS_INLINE SurrealSExpr::BOOST_PP_CAT(SurrealS_, NAME)<ExprL, ExprR, T> \
NAME(const SurrealSType<ExprL, T>& z1, const SurrealSType<ExprR, T>& z2) \
{ return SurrealSExpr::BOOST_PP_CAT(SurrealS_, NAME)<ExprL, ExprR, T>( z1.cast(), z2.cast() ); }

// trig functions <cmath>

SURREALS_FUNC1( cos, cos(z), -sin(z) )
SURREALS_FUNC1( sin, sin(z),  cos(z) )
SURREALS_FUNC1( tan, tan(z),  Real(1)/(cos(z)*cos(z)) )
SURREALS_FUNC1( acos, acos(z), -Real(1)/sqrt(1 - z*z) )
SURREALS_FUNC1( asin, asin(z),  Real(1)/sqrt(1 - z*z) )
SURREALS_FUNC1( atan, atan(z),  Real(1)/(1 + z*z) )

SURREALS_FUNC2( atan2, atan2(z1, z2),  Real(1)/(z1*z1 + z2*z2) )

// hyperbolic functions <cmath>

SURREALS_FUNC1( cosh, cosh(z), sinh(z) )
SURREALS_FUNC1( sinh, sinh(z), cosh(z) )
SURREALS_FUNC1( tanh, tanh(z), Real(1)/(cosh(z)*cosh(z)) )

// exp and log functions <cmath>

SURREALS_FUNC1( exp, exp(z), exp(z) )
SURREALS_FUNC1( expm1, expm1(z), exp(z) )
SURREALS_FUNC1( log, log(z), Real(1)/z )
SURREALS_FUNC1( log10, log10(z), Real(1)/(z*log(10.)) )
SURREALS_FUNC1( log1p, log1p(z), Real(1)/( 1 + z ) )

// error-functions <cmath>

SURREALS_FUNC1( erf , erf(z) ,  Real(2./sqrt(M_PI))*exp(-(z*z)) )
SURREALS_FUNC1( erfc, erfc(z), -Real(2./sqrt(M_PI))*exp(-(z*z)) )

// power functions <cmath>

using ::pow;

namespace SurrealSExpr
{

template<class ExprL, class ExprR, class T>
class SurrealS_pow : public SurrealSType< SurrealS_pow<ExprL, ExprR, T>, T >
{ /*This is for functions when the argument is an expression*/
public:
  static const int N = ExprL::N;
  BOOST_MPL_ASSERT_RELATION(ExprL::N, ==, ExprR::N);

  ALWAYS_INLINE
  SurrealS_pow(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), a(eL.value()), b(eR.value()),
                                                   powab(pow(a,b)),
                                                   tmp1( (a == T(0)) ? ((b == T(1)) ? T(1) : T(0)) : b*pow(a, b - 1) ),
                                                   tmp2( (a == T(0)) ? T(0) : powab*log(a) ) {}

  ALWAYS_INLINE T value() const { return powab; }
  ALWAYS_INLINE T deriv(const int& i) const { return tmp1*eL.deriv(i) + tmp2*eR.deriv(i); }

  ALWAYS_INLINE const SurrealS_pow&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const T a, b, powab, tmp1, tmp2;
};

template<class ExprL, class T>
class SurrealS_pow<ExprL, Real, T> : public SurrealSType< SurrealS_pow<ExprL, Real, T>, T >
{ /*This is optimized when the argument is SurrealS and Real*/
public:
  static const int N = ExprL::N;

  ALWAYS_INLINE
  SurrealS_pow(const ExprL& eL, const Real& b) : eL(eL), a(eL.value()),
                                                 powab(pow(a,b)),
                                                 tmp1( (a == T(0)) ? ((b == 1) ? T(1) : T(0)) : b*pow(a, b - 1) ) {}

  ALWAYS_INLINE T value() const { return powab; }
  ALWAYS_INLINE T deriv(const int& i) const { return tmp1*eL.deriv(i); }

  ALWAYS_INLINE const SurrealS_pow
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const T a, powab, tmp1;
};


template<class ExprR, class T>
class SurrealS_pow<Real, ExprR, T> : public SurrealSType< SurrealS_pow<Real, ExprR, T>, T >
{ /*This is optimized when the argument is a Real and SurrealS*/
public:
  static const int N = ExprR::N;

  ALWAYS_INLINE
  SurrealS_pow(const Real& a, const ExprR& eR) : eR(eR), b(eR.value()),
                                                 powab( (b == T(0)) ? T(1) : pow(a,b) ),
                                                 tmp2( (a == 0) ? T(0) : powab*log(a) ) {}

  ALWAYS_INLINE T value() const { return powab; }
  ALWAYS_INLINE T deriv(const int& i) const { return tmp2*eR.deriv(i); }

  template<int I>
  ALWAYS_INLINE T    get_deriv() const { return tmp2*eR.template get_deriv<I>(); }

  ALWAYS_INLINE const SurrealS_pow
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eR.size(); }
private:
  const ExprR& eR;
  const T b, powab, tmp2;

};

} // namespace SurrealSExpr

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE SurrealSExpr::SurrealS_pow<ExprL, ExprR, T>
pow(const SurrealSType<ExprL, T>& a, const SurrealSType<ExprR, T>& b)
{
  return SurrealSExpr::SurrealS_pow<ExprL, ExprR, T>( a.cast(), b.cast() );
}

template<class Expr, class T, typename S>
ALWAYS_INLINE typename std::enable_if< is_arithmetic_not_SurrealS<S>::value,
                                  SurrealSExpr::SurrealS_pow<Expr, Real, T> >::type
pow(const SurrealSType<Expr, T>& a, const S& b )
{
  return SurrealSExpr::SurrealS_pow<Expr, Real, T>( a.cast(), b );
}

template<class Expr, class T, typename S>
ALWAYS_INLINE typename std::enable_if< is_arithmetic_not_SurrealS<S>::value,
                                  SurrealSExpr::SurrealS_pow<Real, Expr, T> >::type
pow(const S& a, const SurrealSType<Expr, T>& b)
{
  return SurrealSExpr::SurrealS_pow<Real, Expr, T>( a, b.cast() );
}

using ::sqrt;

namespace SurrealSExpr
{

template<class Expr, class T>
class SurrealS_sqrt : public SurrealSType< SurrealS_sqrt<Expr, T>, T >
{ /*This is optimized when the argument is an Expression*/
public:
  static const int N = Expr::N;

  ALWAYS_INLINE
  SurrealS_sqrt(const Expr& e) : e(e), sqrtv( sqrt(e.value()) ), tmp( sqrtv == 0 ? sqrtv : 0.5/sqrtv ) {}

  ALWAYS_INLINE T value() const { return sqrtv; }
  ALWAYS_INLINE T deriv(const int& i) const { return tmp*e.deriv(i); }

  template<int I>
  ALWAYS_INLINE T get_deriv() const { return tmp*e.template get_deriv<I>(); }

  ALWAYS_INLINE const SurrealS_sqrt
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }
private:
  const Expr& e;
  const T sqrtv, tmp;
};
} // namespace SurrealExpr

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::SurrealS_sqrt<Expr, T>
sqrt(const SurrealSType<Expr, T>& z)
{
  return SurrealSExpr::SurrealS_sqrt<Expr, T>( z.cast() );
}


// rounding functions <cmath>

SURREALS_FUNC1( ceil, ceil(z), 0 )
SURREALS_FUNC1( floor, floor(z), 0 )

// misc functions <cmath>

using ::abs;

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, Real, T>
abs( const SurrealSType<Expr, T>& z )
{
  return (z.value() < 0) ?
         SurrealSExpr::OpMul<Expr, Real, T>( z.cast(), -1 ) :
         SurrealSExpr::OpMul<Expr, Real, T>( z.cast(),  1 );
}

using ::fabs;

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, Real, T>
fabs( const SurrealSType<Expr, T>& z )
{
  return (z.value() < 0) ?
         SurrealSExpr::OpMul<Expr, Real, T>( z.cast(), -1 ) :
         SurrealSExpr::OpMul<Expr, Real, T>( z.cast(),  1 );
}

using std::max;

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE SurrealS<ExprL::N,T>
max( const SurrealSType<ExprL, T>& a, const SurrealSType<ExprR, T>& b )
{
  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::N );
  return a.cast().value() > b.cast().value() ? a : b;
}

template<class ExprR, class T>
ALWAYS_INLINE SurrealS<ExprR::N,T>
max( const Real& a, const SurrealSType<ExprR, T>& b )
{
  if ( a > b.cast().value() )
    return a;
  else
    return b;
}

template<class ExprL, class T>
ALWAYS_INLINE SurrealS<ExprL::N,T>
max( const SurrealSType<ExprL, T>& a, const Real& b )
{
  if ( a.cast().value() > b )
    return a;
  else
    return b;
}

using std::min;

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE SurrealS<ExprL::N,T>
min( const SurrealSType<ExprL, T>& a, const SurrealSType<ExprR, T>& b )
{
  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::N );
  return a.cast().value() < b.cast().value() ? a : b;
}

template<class ExprR, class T>
ALWAYS_INLINE SurrealS<ExprR::N,T>
min( const Real& a, const SurrealSType<ExprR, T>& b )
{
  if ( a < b.cast().value() )
    return a;
  else
    return b;
}

template<class ExprL, class T>
ALWAYS_INLINE SurrealS<ExprL::N,T>
min( const SurrealSType<ExprL, T>& a, const Real& b )
{
  if ( a.cast().value() < b )
    return a;
  else
    return b;
}

// I/O

template<int N, class T>
std::istream&
operator>>( std::istream& is, SurrealS<N,T>& z )
{
  Real v = 0;
  Real d[10] = {0};
  char c = 0;
  int n = 0;

  is >> c;
  if (c == '(')
  {
    is >> v;

    is >> c;
    bool done = false;
    while (! done)
    {
      if (c != ')') is.clear(std::ios::badbit);
      if (c == ',')
      {
        is >> d[n]; n++;
      }
      else if (c == ')')
      {
        done = true;
      }
    }
  }
  else
  {
    is.putback(c);
    is >> v;
  }

  if (is) z = SurrealS<N,T>(v, d, n);
  return is;
}


template<class Expr, class T>
std::ostream&
operator<<( std::ostream& os, const SurrealSType<Expr, T>& ztype )
{
  const Expr& z = ztype.cast();
  os << '(' << z.value() << ';';
  for (int i = 0; i < Expr::N - 1; i++)
    os << z.deriv(i) << ',';
  os << z.deriv(Expr::N - 1) << ')';
  return os;
}

} // namespace SANS

//Created specialized version of fpt_abs to work with the boost unit testing framework

namespace boost
{
namespace test_tools
{
namespace tt_detail
{

template<typename FPT>
FPT
fpt_abs( FPT fpv );

template<typename FPT>
FPT
safe_fpt_division( FPT f1, FPT f2 );


template<int N, class T>
inline Real
fpt_abs( const SANS::SurrealSExpr::OpMul<SANS::SurrealS<N,T>, Real, T>& fpv )
{
  Real val = fpv.value();
  return fpt_abs( val );
}

template<int N, class T>
inline Real
fpt_abs( const SANS::SurrealSExpr::OpSub<SANS::SurrealS<N,T>, SANS::SurrealS<N,T>, T>& fpv )
{
  Real val = fpv.value();
  return fpt_abs( val );
}

// both f1 and f2 are unsigned here
template<int N, class T>
inline Real
safe_fpt_division( const SANS::SurrealS<N,T>& f1, const SANS::SurrealS<N,T>& f2 )
{
  Real val1 = f1.value();
  Real val2 = f2.value();
  return safe_fpt_division( val1, val2 );
}

} // namespace tt_detail
} // namespace test_tools
} // namespace boost


//Clean up macro definitions
#undef SURREALS_FUNC1
#undef SURREALS_FUNC2

#endif // SANS_SURREALS_LAZY_H
