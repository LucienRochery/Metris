// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MATRIXS_MUL_NATIVE_H
#define MATRIXS_MUL_NATIVE_H

#include "SANS/tools/SANSnumerics.h"     // Real
#include "SANS/tools/minmax.h"

#include "SANS/tools/CacheLineSize.h"
#include "SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS_Type.h"

namespace SANS
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
} //namespace SANS

#endif //MATRIXS_MUL_NATIVE_H
