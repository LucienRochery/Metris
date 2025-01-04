// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef VECTORS_CLASS_H
#define VECTORS_CLASS_H

// Vector class with compile-time size

// NOTES:
// - indexing: zero-based


#include <iostream>
#include <vector>

#include "SANS/tools/SANSTraitsPOD.h"
#include "SANS/tools/SANSException.h"
#include "SANS/tools/SANSTraitsInitListAssign.h"

#include "SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"
#include "SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS_Mul.h"

namespace SANS
{
#ifdef __INTEL_COMPILER
namespace DLA
{

//Forward declaration
template< int M, class T >
class VectorS;

}

//Create a specialization so to allow for the syntax
//   VectorD< VectorS<2,Real> >
//      v = { {3,3}, {3,2} };
//
// This is completely unnessary if the intel compiler could use templated initializer_list functions....
//
template<int M, class T>
struct initializer_list_assign< DLA::VectorS< M, T > >
{
  template<class U>
  initializer_list_assign(DLA::VectorS< M, T >& val, const std::initializer_list<U>& s) { val = s; }
};
#endif

namespace DLA
{

//----------------------------------------------------------------------------//
// VectorS<M,T>:  vector w/ M data elements, each of type T
//
// template parameters:
//   M            array dimension
//   T            array data type (e.g. double)
//----------------------------------------------------------------------------//

template <int M_, class T>
class VectorS : public MatrixS<M_, 1, T>
{
public:
  typedef T Ttype;
  typedef MatrixS<M_, 1, T> BaseType;
  using BaseType::M;

  //Default constructor does not initialize the data. This makes it behave more like POD
  //and valgrind can catch any use of uninitialized data
  VectorS() = default;
  VectorS(VectorS const&) = default;
  VectorS(VectorS &&) = default;
  VectorS& operator=(VectorS const&) = default;
  VectorS& operator=(VectorS &&) = default;

  VectorS( const T& s );   // missing 'explicit' allows VectorS<M,Real> q = 0
  VectorS( const T s[], int nsize );
  VectorS( const std::initializer_list<T>& s ) { operator=(s); }
  VectorS( const std::vector<T>& s ) { operator=(s); }
  VectorS( const typename SANS::POD<T>::type& s );

  //This allows a vector to be assigned an expression as it is constructed, i.e.
  //VectorS<M,T> C = A + B;
  template<class Expr, bool useRF, bool Full>
  VectorS( const MatrixSType<Expr, useRF, Full>& r ) { *this = r; }

  // accessor operators
        T& operator[]( int n )       { return data[n]; }
  const T& operator[]( int n ) const { return data[n]; }
        T& operator()( int n )       { return data[n]; }
  const T& operator()( int n ) const { return data[n]; }

  // assignment
  VectorS& operator=( const T& s );
  VectorS& operator=( const typename SANS::POD<T>::type& s );
  VectorS& operator=( const std::initializer_list<T>& s );
  VectorS& operator=( const std::vector<T>& s );
  template<class T2>
  VectorS& operator=( const std::initializer_list<T2>& s );
  template<class Expr, bool Full> VectorS& operator= ( const MatrixSType<Expr, true, Full>& );
  template<class Expr> VectorS& operator= ( const MatrixSType<Expr, false, true>& );

  // compound assignment
  using BaseType::operator+=;
  VectorS& operator+=( const std::initializer_list<T>& s );

  void dump( int indentSize=0, std::ostream& out = std::cout ) const;

private:
  using MatrixS<M_,1,T>::data;
};


// constructors

template <int M, class T>
inline
VectorS<M,T>::VectorS( const T& s )
{
  for (int i = 0; i < M; i++)
    data[i] = s;
}

template <int M, class T>
inline
VectorS<M,T>::VectorS( const T s[], int nsize )
{
  SANS_ASSERT(nsize == M);
  for (int n = 0; n < M; n++)
    data[n] = s[n];
}

// needed for VectorS<M, Surreal<M> >(Real)
template <int M, class T>
inline
VectorS<M,T>::VectorS( const typename SANS::POD<T>::type& s )
{
  for (int i = 0; i < M; i++)
    data[i] = s;
}

// assignment

template <int M, class T>
inline VectorS<M,T>&
VectorS<M,T>::operator=( const T& s )
{
  for (int i = 0; i < M; i++)
    data[i] = s;
  return *this;
}

// needed for  VectorS<M, Surreal> q; q = 0;
template <int M, class T>
inline VectorS<M,T>&
VectorS<M,T>::operator=( const typename SANS::POD<T>::type& s )
{
  for (int i = 0; i < M; i++)
    data[i] = s;
  return *this;
}

template <int M, class T>
inline VectorS<M,T>&
VectorS<M,T>::operator=( const std::initializer_list<T>& s )
{
  SANS_ASSERT_MSG(s.size() == (std::size_t)M, "s.size() = %d, M = %d", s.size(), M);
  int n = 0;
  auto row = s.begin();
  for (std::size_t i = 0; i < M; ++i, row++)
      data[n++] = *row;

  return *this;
}

template <int M, class T>
inline VectorS<M,T>&
VectorS<M,T>::operator=( const std::vector<T>& s )
{
  SANS_ASSERT_MSG(s.size() == (std::size_t)M, "s.size() = %d, M = %d", s.size(), M);
  int n = 0;
  auto row = s.begin();
  for (std::size_t i = 0; i < M; ++i, row++)
      data[n++] = *row;

  return *this;
}

template <int M, class T>
template<class T2>
inline VectorS<M,T>&
VectorS<M,T>::operator=( const std::initializer_list<T2>& s )
{
  SANS_ASSERT_MSG(s.size() == (std::size_t)M, "s.size() = %d, M = %d", s.size(), M);
  int n = 0;
  auto row = s.begin();
  for (std::size_t i = 0; i < M; ++i, row++)
    data[n++] = *row;

  return *this;
}

template <int M, class T>
template< class Expr, bool Full >
inline VectorS<M,T>&
VectorS<M,T>::operator=( const MatrixSType<Expr, true, Full>& r )
{
  const Expr& Tree = r.cast();

  BOOST_MPL_ASSERT_RELATION( M, ==, Expr::M );
  BOOST_MPL_ASSERT_RELATION( 1, ==, Expr::N );

  // casting removes unneeded instantiations with VectorS
  Tree.value(1., static_cast<MatrixS<M,1,T>&>(*this));

  return *this;
}

template <int M, class T>
template< class Expr >
inline VectorS<M,T>&
VectorS<M,T>::operator=( const MatrixSType<Expr, false, true>& r )
{
  const Expr& Tree = r.cast();

  BOOST_MPL_ASSERT_RELATION( M, ==, Expr::M );
  BOOST_MPL_ASSERT_RELATION( 1, ==, Expr::N );

  for ( int i = 0; i < M; i++ )
    data[i] = Tree.value(i);

  return *this;
}

template <int M, class T>
inline VectorS<M,T>&
VectorS<M,T>::operator+=( const std::initializer_list<T>& s )
{
  SANS_ASSERT_MSG(s.size() == (std::size_t)M, "s.size() = %d, M = %d", s.size(), M);
  int n = 0;
  auto row = s.begin();
  for (std::size_t i = 0; i < M; ++i, row++)
      data[n++] += *row;

  return *this;
}


// debug dump of private data
template <int M, class T>
void
VectorS<M,T>::dump( int indentSize, std::ostream& out ) const
{
  std::string indent(indentSize, ' ');
  out << indent << "VectorS<" << M << ",T>:" << std::endl;
#if 1
  out << indent << "  data = ";
  for (int n = 0; n < M; n++)
    out << data[n] << " ";
  out << std::endl;
#else     // only works for class T with member function dump()
  for (int n = 0; n < M; n++)
  {
    out << indent << "  data[" << n << "] = ";
    data[n].dump(indentSize, out);
  }
#endif
}

template <int k, int M, class T>
const VectorS<M,T>&
get(const VectorS<M,T>& v)
{
  static_assert( k == 0 || k == -1, "k should be zero or -1 if the argument to get is VectorS<M,T>");
  return v;
}

// I/O

template <int M, class T>
std::ostream&
operator<<( std::ostream& out, const VectorS<M,T>& v )
{
  for (int i = 0; i < M-1; i++)
    out << v[i] << ", ";
  out << v[M-1];
  return out;
}


} //namespace DLA
} //namespace SANS

#endif // VECTORS_CLASS_H
