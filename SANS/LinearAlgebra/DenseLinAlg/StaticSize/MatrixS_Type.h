// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_TYPE_H
#define MATRIXS_TYPE_H


namespace SANS
{
namespace DLA
{

//Forward declarations
template <int M, int N, class T>
class MatrixS;

template <int M, class T>
class MatrixSymS;

template <int M, class T>
class VectorS;

template<class L, class R, bool useRF, bool MatrixFull>
class OpAddS;

template<class L, class R, bool useRF, bool MatrixFull>
class OpSubS;

template<class ExprL, class ExprR>
class OpMulS;

template<class Expr, class S, bool useRF, bool MatrixFull>
class OpMulSScalar;

template<class ExprL, class ExprR, bool MatrixFull>
class OpMulSFactor;

template<int M, class T>
class MatrixSDiag;

template<int M_, int N_, class T>
class MatrixSTranspose;

//Base class to distinguish what types belong to the static size matrix family
//useRF is a flag to indicate the expression includes a matmul
//MatrixFull indicates if a full matrix is part of the expression
template< class Derived, bool useRF, bool MatrixFull >
struct MatrixSType
{
  //A convenient method for casting to the derived type
  inline       Derived& cast()       { return static_cast<      Derived&>(*this); }
  inline const Derived& cast() const { return static_cast<const Derived&>(*this); }
};

}
}

#endif //MATRIXS_TYPE_H
