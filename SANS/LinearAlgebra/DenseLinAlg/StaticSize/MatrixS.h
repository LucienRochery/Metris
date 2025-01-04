// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_CLASS_H
#define MATRIXS_CLASS_H

// Matrix class with compile-time size

// NOTES:
// - indexing: zero-based
// - C-style matrix storage (row first)

#include <array>
#include <iostream>
#include <stdint.h> // uintptr_t
#include <initializer_list>

// Surreals are registered with boost::is_arithmetic.
// not sure how to register them with std::is_arithmetic
#include <boost/type_traits/is_arithmetic.hpp>

#include "SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS_Type.h"
#include "SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS_Add.h"
#include "SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS_Mul.h"

#include "SANS/LinearAlgebra/DenseLinAlg/tools/Matrix_Util.h"
#include "SANS/LinearAlgebra/DenseLinAlg/tools/Identity.h"

#include "SANS/tools/SANSTraitsPOD.h"
#include "SANS/tools/SANSException.h"
#include "SANS/tools/SANSTraitsInitListAssign.h"

namespace SANS
{
#ifdef __INTEL_COMPILER
namespace DLA
{

//Forward declaration
template< int M,  int N, class T >
class MatrixS;

}

//Create a specialization so to allow for the syntax
//   VectorD< VectorS<2,Real> >
//      v = { {3,3}, {3,2} };
//
// This is completely unnessary if the intel compiler could use templated initializer_list functions....
//
template<int M, int N, class T>
struct initializer_list_assign< DLA::MatrixS< M, N, T > >
{
  template<class U>
  initializer_list_assign(DLA::MatrixS< M, N, T >& val, const std::initializer_list< std::initializer_list<U> >& s) { val = s; }
};
#endif

namespace DLA
{

//----------------------------------------------------------------------------//
// MatrixS<M,N,T>:  matrix w/ MxN data elements, each of type T
//
// template parameters:
//   MxN          matrix dimension
//   T            matrix data type (e.g. double)
//
// Notes:
//  -implemented as 1-D array with C-style (row-first) matrix storage
//----------------------------------------------------------------------------//

template <int M_, int N_, class T>
class MatrixS : public MatrixSType< MatrixS<M_, N_, T>, !boost::is_arithmetic<T>::value, true >
{
public:
  typedef T Ttype;
  static constexpr int M = M_;
  static constexpr int N = N_;

  //Default constructor does not initialize the data. This makes it behave more like POD
  //and valgrind can catch any use of uninitialized data
  MatrixS() = default;
  MatrixS(MatrixS const&) = default;
  MatrixS(MatrixS&&) = default;
  MatrixS& operator=(MatrixS const&) = default;
  MatrixS& operator=(MatrixS&&) = default;

  template<class Z>
  MatrixS( const MatrixSymS<M,Z>& m ) { operator=(m); }
  MatrixS( const Identity& I ) { (*this) = I; };
  MatrixS( const T s );  // missing 'explicit' allows MatrixS<M,N,Real> q = 0
  MatrixS( const T s[], int n );
  MatrixS( const std::initializer_list< const std::initializer_list<T> >& s ) { operator=(s); }
  MatrixS( const typename SANS::POD<T>::type& s );


  //Sub matrix
  MatrixS<1,N_,T> row( const int i ) const;
  MatrixS<M_,1,T> col( const int j ) const;

  //This allows a matrix to be assigned an expression as it is constructed, i.e.
  //MatrixS<M,N,T> C = A + B;
  template<class Expr, bool useRF, bool MatrixFull>
  MatrixS( const MatrixSType<Expr, useRF, MatrixFull>& r ) { *this = r; }

  // accessor operators
        T& operator()( int i, int j )       { return data[N*i + j]; }
  const T& operator()( int i, int j ) const { return data[N*i + j]; }

  // assignment
  template<class Z>
  MatrixS& operator=( const MatrixSymS<M,Z>& m );
  MatrixS& operator=( const T& s );
  MatrixS& operator=( const typename SANS::POD<T>::type& s );
  MatrixS& operator=( const Identity& I );
  MatrixS& operator=( const std::initializer_list< const std::initializer_list<T> >& s );

  // unary operators; no side effects
  const MatrixS& operator+() const;

  // binary accumulation operators
  MatrixS& operator+=( const MatrixS& m );
  MatrixS& operator-=( const MatrixS& m );
  MatrixS& operator+=( const Identity& I );
  MatrixS& operator-=( const Identity& I );
  MatrixS& operator+=( const T& s );
  MatrixS& operator-=( const T& s );
  MatrixS& operator*=( const T& s );
  MatrixS& operator/=( const T& s );
  MatrixS& operator+=( const typename SANS::POD<T>::type& s );
  MatrixS& operator-=( const typename SANS::POD<T>::type& s );
  MatrixS& operator*=( const typename SANS::POD<T>::type& s );
  MatrixS& operator/=( const typename SANS::POD<T>::type& s );

  // Lazy expression with recursive functions assignment and binary accumulation
  template<class Expr, bool Full> MatrixS& operator= ( const MatrixSType<Expr, true, Full>& );
  template<class Expr, bool Full> MatrixS& operator+=( const MatrixSType<Expr, true, Full>& );
  template<class Expr, bool Full> MatrixS& operator-=( const MatrixSType<Expr, true, Full>& );
  template<class Expr, bool useRF, bool Full> MatrixS& operator*=( const MatrixSType<Expr, useRF, Full>& );

  // Lazy expression element-wise assignment and binary accumulation
  template<class Expr> MatrixS& operator= ( const MatrixSType<Expr, false, true>& );
  template<class Expr> MatrixS& operator+=( const MatrixSType<Expr, false, true>& );
  template<class Expr> MatrixS& operator-=( const MatrixSType<Expr, false, true>& );

  void dump( int indentSize=0, std::ostream& out = std::cout ) const;

  //Column manipulations
  inline void swap_cols(const int i, const int j );

  //Row manipulations
  inline void swap_rows(const int i, const int j );
  template< class Ta >
  inline void scale_row(const int i, const Ta& a, const int start = 0, int end = 0);
  template< class Ta >
  inline void axpy_rows(const int i, const int j, const Ta& a, const int start = 0, int end = 0);
  inline int max_row_in_col(const int j, const int start = 0) const;

  //A unique ID for the matrix so the matrix can be uniquely identified
  uintptr_t ID() const { return (uintptr_t)this; }

  //Lazy Expression Operators
  // Assign the value of this matrix to res
  template<class Scalar, class Tres>
  inline void value(const Scalar& sgn, MatrixS<M,N,Tres>& res) const
  {
    if ( ID() == res.ID() && sgn == 1 ) return; //Assigning this to this, nothing to do.

    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N; ++j)
        res(i,j) = sgn*(*this)(i,j);
  }

  // Add the value of this matrix to res
  template<class Scalar, class Tres>
  inline void plus(const Scalar& sgn, MatrixS<M,N,Tres>& res) const
  {
    SANS_ASSERT( ID() != res.ID() ); //The variable on the left of the assignment cannot also appear on the right

    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N; ++j)
        res(i,j) += sgn*(*this)(i,j);
  }

  //Element-wise expression
  inline const T& value(const int& i) const
  {
    return data[i];
  }
protected:
  std::array<T, M*N> data;
};


// constructors
template <int M, int N, class T>
inline
MatrixS<M,N,T>::MatrixS( const T s[], int n )
{
  SANS_ASSERT(n == M*N);
  for (n = 0; n < M*N; n++)
    data[n] = s[n];
}

template <int M, int N, class T>
inline
MatrixS<M,N,T>::MatrixS( const T s )
{
  for (int n = 0; n < M*N; n++)
    data[n] = s;
}

// needed for MatrixS<M,N, Surreal<Z>>(Real)
template <int M, int N, class T>
inline
MatrixS<M,N,T>::MatrixS( const typename SANS::POD<T>::type& s )
{
  for (int n = 0; n < M*N; n++)
    data[n] = s;
}

// sub-matrix

template <int M, int N, class T>
inline MatrixS<1,N,T>
MatrixS<M,N,T>::row( const int i ) const
{
  MatrixS<1,N,T> m;
  for (int j = 0; j < N; j++)
    m(0,j) = (*this)(i,j);
  return m;
}

template <int M, int N, class T>
inline MatrixS<M,1,T>
MatrixS<M,N,T>::col( const int j ) const
{
  MatrixS<M,1,T> m;
  for (int i = 0; i < M; i++)
    m(i,0) = (*this)(i,j);
  return m;
}

// assignment
template <int M, int N, class T>
template <class Z>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator=( const MatrixSymS<M,Z>& m )
{
  // can only assign a symmetric matrix to a square matrix
  BOOST_MPL_ASSERT_RELATION(M, ==, N);

  for (int i = 0; i < M; i++)
    for (int j = 0; j < N; j++)
      (*this)(i,j) = m(i,j);

  return *this;
}

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator=( const T& s )
{
  for (int i = 0; i < M*N; i++)
    data[i] = s;
  return *this;
}

// needed for  MatrixS<M,N, Surreal> q; q = 0;
template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator=( const typename SANS::POD<T>::type& s )
{
  for (int i = 0; i < M*N; i++)
    data[i] = s;
  return *this;
}

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator=( const Identity& I )
{
  *this = 0;
  if ( M == N )
    for (int i = 0; i < M; i++)
      (*this)(i,i) = I;
  else
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        if ( i == j )
          (*this)(i,j) = I;

  return *this;
}

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator=( const std::initializer_list< const std::initializer_list<T> >& s )
{
  SANS_ASSERT(s.size() == M);
  size_t n = 0;
  for (const auto& row: s)
  {
    SANS_ASSERT(row.size() == N);
    for (const auto& col: row)
      data.at(n++) = col;
  }
  SANS_ASSERT(n == data.size());
  return *this;
}

// unary operators; no side effects

template <int M, int N, class T>
inline const MatrixS<M,N,T>&
MatrixS<M,N,T>::operator+() const
{
  return *this;
}


// binary accumulation operators

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator+=( const MatrixS& m )
{
  for (int i = 0; i < M*N; i++)
    data[i] += m.data[i];
  return *this;
}

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator-=( const MatrixS& m )
{
  for (int i = 0; i < M*N; i++)
    data[i] -= m.data[i];
  return *this;
}

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator+=( const Identity& I )
{
  if ( M == N )
    for (int i = 0; i < M; i++)
      (*this)(i,i) += I;
  else
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        if ( i == j )
          (*this)(i,j) += I;

  return *this;
}

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator-=( const Identity& I )
{
  if ( M == N )
    for (int i = 0; i < M; i++)
      (*this)(i,i) -= I;
  else
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++)
        if ( i == j )
          (*this)(i,j) -= I;

  return *this;
}

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator+=( const T& s )
{
  for (int i = 0; i < M*N; i++)
    data[i] += s;
  return *this;
}

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator-=( const T& s )
{
  for (int i = 0; i < M*N; i++)
    data[i] -= s;
  return *this;
}

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator*=( const T& s )
{
  for (int i = 0; i < M*N; i++)
    data[i] *= s;
  return *this;
}

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator/=( const T& s )
{
  for (int i = 0; i < M*N; i++)
    data[i] /= s;
  return *this;
}

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator+=( const typename SANS::POD<T>::type& s )
{
  for (int i = 0; i < M*N; i++)
    data[i] += s;
  return *this;
}

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator-=( const typename SANS::POD<T>::type& s )
{
  for (int i = 0; i < M*N; i++)
    data[i] -= s;
  return *this;
}

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator*=( const typename SANS::POD<T>::type& s )
{
  for (int i = 0; i < M*N; i++)
    data[i] *= s;
  return *this;
}

template <int M, int N, class T>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator/=( const typename SANS::POD<T>::type& s )
{
  for (int i = 0; i < M*N; i++)
    data[i] /= s;
  return *this;
}

// Lazy expression with recursive functions assignment and binary accumulation

// Generically assume that expressions do not reduce to scalars
template<class Expr, class Matrix>
struct cast_to_scalar
{
  static void value(const Expr& Tree, Matrix& A)
  {
    BOOST_MPL_ASSERT_RELATION( Matrix::M, ==, Expr::M );
    BOOST_MPL_ASSERT_RELATION( Matrix::N, ==, Expr::N );

    Tree.value(1., A);
  }

  static void plus(const Real sgn, const Expr& Tree, Matrix& A)
  {
    BOOST_MPL_ASSERT_RELATION( Matrix::M, ==, Expr::M );
    BOOST_MPL_ASSERT_RELATION( Matrix::N, ==, Expr::N );

    Tree.plus(sgn, A);
  }
};

// Default for multiplication expression that does not reduce to a scalar
template<int M, int N, class ExprT, class Matrix>
struct cast_to_scalar_mul
{
  template<class Expr>
  static void value(const Expr& Tree, Matrix& A)
  {
    BOOST_MPL_ASSERT_RELATION( Matrix::M, ==, Expr::M );
    BOOST_MPL_ASSERT_RELATION( Matrix::N, ==, Expr::N );

    Tree.value(1., A);
  }

  template<class Expr>
  static void plus(const Real sgn, const Expr& Tree, Matrix& A)
  {
    BOOST_MPL_ASSERT_RELATION( Matrix::M, ==, Expr::M );
    BOOST_MPL_ASSERT_RELATION( Matrix::N, ==, Expr::N );

    Tree.plus(sgn, A);
  }
};

// Used with MatrixS = MatrixS<MatrixS> * MatrixS<MatrixS>, i.e. the multiplication reduces to the 'scalar' MatrixS
template<int M, int N, class T>
struct cast_to_scalar_mul<1,1,MatrixS<M,N,T>,MatrixS<M,N,T>>
{
  template<class Expr>
  static void value(const Expr& Tree, MatrixS<M,N,T>& A)
  {
    Tree.value(1., A);
  }

  template<class Expr>
  static void plus(const Real sgn, const Expr& Tree, MatrixS<M,N,T>& A)
  {
    Tree.plus(sgn, A);
  }
};

template<int M, class T>
struct cast_to_scalar_mul<1,1,VectorS<M,T>,MatrixS<M,1,T>>
{
  template<class Expr>
  static void value(const Expr& Tree, MatrixS<M,1,T>& A)
  {
    Tree.value(1., A);
  }

  template<class Expr>
  static void plus(const Real sgn, const Expr& Tree, MatrixS<M,1,T>& A)
  {
    Tree.plus(sgn, A);
  }
};


// Multiplication could reduce to the 'scalar' in MatrixS<MatrixS>
template<class ExprL, class ExprR, int M, int N, class T>
struct cast_to_scalar<OpMulS<ExprL,ExprR>, MatrixS<M,N,T>>
{
  typedef OpMulS<ExprL,ExprR> Expr;
  static void value(const Expr& Tree, MatrixS<M,N,T>& A)
  {
    cast_to_scalar_mul<Expr::M, Expr::N, typename Expr::Ttype, MatrixS<M,N,T>>::value(Tree,A);
  }

  static void plus(const Real sgn, const Expr& Tree, MatrixS<M,N,T>& A)
  {
    cast_to_scalar_mul<Expr::M, Expr::N, typename Expr::Ttype, MatrixS<M,N,T>>::plus(sgn, Tree,A);
  }
};

template<class ExprL, class ExprR, bool Full, int M, int N, class T>
struct cast_to_scalar<OpMulSFactor<ExprL,ExprR,Full>, MatrixS<M,N,T>>
{
  typedef OpMulSFactor<ExprL,ExprR,Full> Expr;
  static void value(const Expr& Tree, MatrixS<M,N,T>& A)
  {
    cast_to_scalar_mul<Expr::M, Expr::N, typename Expr::Ttype, MatrixS<M,N,T>>::value(Tree,A);
  }

  static void plus(const Real sgn, const Expr& Tree, MatrixS<M,N,T>& A)
  {
    cast_to_scalar_mul<Expr::M, Expr::N, typename Expr::Ttype, MatrixS<M,N,T>>::plus(sgn, Tree, A);
  }
};

template<class ExprS, class S, bool useRF, bool Full, int M, int N, class T>
struct cast_to_scalar<OpMulSScalar<ExprS,S,useRF,Full>, MatrixS<M,N,T>>
{
  typedef OpMulSScalar<ExprS,S,useRF,Full> Expr;
  static void value(const Expr& Tree, MatrixS<M,N,T>& A)
  {
    cast_to_scalar_mul<Expr::M, Expr::N, typename Expr::Ttype, MatrixS<M,N,T>>::value(Tree,A);
  }

  static void plus(const Real sgn, const Expr& Tree, MatrixS<M,N,T>& A)
  {
    cast_to_scalar_mul<Expr::M, Expr::N, typename Expr::Ttype, MatrixS<M,N,T>>::plus(sgn, Tree,A);
  }
};


template <int M, int N, class T>
template< class Expr, bool Full >
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator=( const MatrixSType<Expr, true, Full>& r )
{
  // Just to suppress some strange clang analyzer warnings
#ifdef __clang_analyzer__
  for ( int i = 0; i < M*N; i++ )
    data[i] = 0;
#endif

  // Used with MatrixS = MatrixS<MatrixS> * MatrixS<MatrixS>
  // only needed for Full which can be a matrix multipliction
  cast_to_scalar<Expr, MatrixS<M,N,T>>::value(r.cast(),*this);

  return *this;
}

template <int M, int N, class T>
template< class Expr, bool Full>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator+=( const MatrixSType<Expr, true, Full>& r )
{
  // Used with MatrixS = MatrixS<MatrixS> * MatrixS<MatrixS>
  // only needed for Full which can be a matrix multipliction
  cast_to_scalar<Expr, MatrixS<M,N,T>>::plus(1., r.cast(),*this);

  return *this;
}

template <int M, int N, class T>
template< class Expr, bool Full>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator-=( const MatrixSType<Expr, true, Full>& r )
{
  // Used with MatrixS = MatrixS<MatrixS> * MatrixS<MatrixS>
  // only needed for Full which can be a matrix multipliction
  cast_to_scalar<Expr, MatrixS<M,N,T>>::plus(-1., r.cast(),*this);

  return *this;
}

template <int M, int N, class T>
template< class Expr, bool useRF, bool Full >
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator*=( const MatrixSType<Expr, useRF, Full>& r )
{
  MatrixS<M,N,T> tmp(*this);
  *this = tmp*r;
  return *this;
}

// Element wise lazy expression assignment and binary accumulation

template <int M, int N, class T>
template< class Expr>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator=( const MatrixSType<Expr, false, true>& r )
{
  const Expr& Tree = r.cast();

  BOOST_MPL_ASSERT_RELATION( M, ==, Expr::M );
  BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

  // Just to suppress some strange clang analyzer warnings
#ifdef __clang_analyzer__
  for ( int i = 0; i < M*N; i++ )
    data[i] = 0;
#endif

  for ( int i = 0; i < M*N; i++ )
    data[i] = Tree.value(i);

  return *this;
}

template <int M, int N, class T>
template< class Expr>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator+=( const MatrixSType<Expr, false, true>& r )
{
  const Expr& Tree = r.cast();

  BOOST_MPL_ASSERT_RELATION( M, ==, Expr::M );
  BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

  for ( int i = 0; i < M*N; i++ )
    data[i] += Tree.value(i);

  return *this;
}

template <int M, int N, class T>
template< class Expr>
inline MatrixS<M,N,T>&
MatrixS<M,N,T>::operator-=( const MatrixSType<Expr, false, true>& r )
{
  const Expr& Tree = r.cast();

  BOOST_MPL_ASSERT_RELATION( M, ==, Expr::M );
  BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

  for ( int i = 0; i < M*N; i++ )
    data[i] -= Tree.value(i);

  return *this;
}

// debug dump of private data
template <int M, int N, class T>
void
MatrixS<M,N,T>::dump( int indentSize, std::ostream& out ) const
{
  std::string indent(indentSize, ' ');
  out << indent << "MatrixS<" << M << "," << N << ",T>:" << std::endl;
#if 1
  out << indent << "  data = ";
  for (int n = 0; n < M*N; n++)
    out << data[n] << " ";
  out << std::endl;
#else     // only works for class T with member function dump()
  for (int n = 0; n < M*N; n++)
  {
    out << indent << "  data[" << n << "] = ";
    data[n].dump(indentSize);
  }
#endif
}

//-----------------------------------------------------------------------------
//Column operations

//Swaps columns i and j
template< int M, int N, class T >
inline void
MatrixS<M,N,T>::swap_cols(const int i, const int j )
{
  if ( i == j ) return; //Nothing to do...
  //Swap two columns
  T* __restrict x = &(*this)(0,i);
  T* __restrict y = &(*this)(0,j);
  for (int k = 0; k < M; ++k)
  {
    T tmp = x[N*k];
    x[N*k]  = y[N*k];
    y[N*k]  = tmp;
  }
}

//-----------------------------------------------------------------------------
//Row operations

//Swaps rows i and j
template< int M, int N, class T >
inline void
MatrixS<M,N,T>::swap_rows(const int i, const int j )
{
  if ( i == j ) return; //Nothing to do...
  //Swap two rows
  T* __restrict x = &(*this)(i,0);
  T* __restrict y = &(*this)(j,0);
  for (int k = 0; k < N; ++k)
  {
    T tmp = x[k];
    x[k]  = y[k];
    y[k]  = tmp;
  }
}

//Scales the row i
template< int M, int N, class T >
template< class Ta >
inline void
MatrixS<M,N,T>::scale_row(const int i, const Ta& a, const int start, int end)
{
  if (end == 0) end = N;
  MatrixUtil_Native<T,Ta>::scal(&(*this)(i,0) + start, a, end - start);
}

//Performs y = a*x + y where x = row(i) and y = row(j)
template< int M, int N, class T >
template< class Ta >
inline void
MatrixS<M,N,T>::axpy_rows(const int i, const int j, const Ta& a, const int start, int end)
{
  if (end == 0) end = N;
  MatrixUtil_Native<T,Ta>::axpy(&(*this)(i,0) + start, &(*this)(j,0) + start, a, end - start);
}

//Finds the row with maximum value in a column
template< int M, int N, class T >
inline int
MatrixS<M,N,T>::max_row_in_col(const int j, const int start) const
{
  return MatrixUtil_Native<T,T>::max_row_in_col(&(*this)(0,j), M, N, start);
}

//Unary negation

template< class Expr, bool useRF, bool Full >
inline const OpMulSScalar<Expr, Real, useRF, Full>
operator-(MatrixSType<Expr, useRF, Full> const& e)
{
  return OpMulSScalar<Expr, Real, useRF, Full>( e.cast(), -1 );
}

// I/O

template <int M, int N, class T>
std::ostream&
operator<<( std::ostream& out, const MatrixS<M,N,T>& m )
{
//#define USE_MATHEMATICA_STYLE
#ifdef USE_MATHEMATICA_STYLE
  out << "{";
  for (int i = 0; i < M; i++)
  {
    out << "{" << m(i,0);
    for (int j = 1; j < N; j++)
      out << ", " << m(i,j);
    if (i < M-1) out << "}, ";
  }
  out << "}}";
#else
  for (int i = 0; i < M; i++)
  {
    out << "{";
    for (int j = 0; j < N; j++)
    {
      out << m(i,j);
      if ( j == N-1 ) out << "}";
      else out << ", ";
    }
    if ( i < M-1 )out << std::endl;
  }
#endif
  return out;
}

template <int M, int N, class T>
std::istream&
operator>>( std::istream &in, MatrixS<M,N,T>& m )
{
  for (int i = 0; i < M; i++)
    for (int j = 0; j < N; j++)
      in >> m(i,j);
  return in;
}

} //namespace DLA
} //namespace SANS

#endif // MatrixS_H
