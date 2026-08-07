// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef METRIS_DLA_MATRIXS_MATMUL_NATIVE_H
#define METRIS_DLA_MATRIXS_MATMUL_NATIVE_H

#include "tools/SANSnumerics.h"     // Real
#include "tools/minmax.h"

#include "tools/CacheLineSize.h"
#include "LinearAlgebra/DenseLinAlg/StaticSize/MatrixS_Type.h"

namespace Metris
{
namespace DLA
{

template<class TL, class TR, class S, class T>
class MatrixS_MatMul_Native
{
  static const int CacheItemsL = MAX( 1, CACHE_LINE_SIZE / sizeof(TL) );
  static const int CacheItemsR = MAX( 1, CACHE_LINE_SIZE / sizeof(TR) );

public:
//-----------------------------------------------------------------------------
  template<int ML, int NL, int MR, int NR>
  static inline void value(const MatrixS<ML,NL,TL>& ml, const MatrixS<MR,NR,TR>& mr, const S& sgn, MatrixS<ML,NR,T>& res )
  {
    res = 0;
    plus(ml, mr, sgn, res);
  }

//-----------------------------------------------------------------------------
  template<int ML, int NL, int MR, int NR>
  static void plus(const MatrixS<ML,NL,TL>& ml, const MatrixS<MR,NR,TR>& mr, const S& sgn, MatrixS<ML,NR,T>& res );

};

} //namespace DLA
} //namespace Metris

#endif // METRIS_DLA_MATRIXS_MATMUL_NATIVE_H
