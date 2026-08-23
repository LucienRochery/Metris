// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef METRIS_DLA_MATRIXS_MATMUL_NATIVE_IMPL_H
#define METRIS_DLA_MATRIXS_MATMUL_NATIVE_IMPL_H

// This file is only for explicit instantiation.

#if !defined(METRIS_MATRIXS_MATMUL_NATIVE_INSTANTIATE) && !defined(METRIS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "MatrixS_MatMul_Native.h"
#include "../../../../tools/SANSException.h"
#include "../MatrixS.h"

#include <boost/mpl/assert.hpp>

namespace Metris
{
namespace DLA
{

template<class TL, class TR, class S, class T>
template<int ML, int NL, int MR, int NR>
void MatrixS_MatMul_Native<TL,TR,S,T>::plus(const MatrixS<ML,NL,TL>& ml, const MatrixS<MR,NR,TR>& mr,
                                            const S& sgn, MatrixS<ML,NR,T>& res )
{
  BOOST_MPL_ASSERT_RELATION(NL, ==, MR);

  METRIS_SUPPORT_ASSERT(res.ID() != ml.ID());
  METRIS_SUPPORT_ASSERT(res.ID() != mr.ID());

  static const int Rstride = NR;

  int i1, j1, k1, j2, k2;
  T *__restrict rC;

  for (i1 = 0; i1 < ML; ++i1)
  {
    T *C = &res(i1,0);
    const TL *L1 = &ml(i1,0);
    const TR *R1 = &mr(0,0);

//---------------
    for (j1 = 0; j1 < (NL/CacheItemsL)*CacheItemsL; j1 += CacheItemsL, R1 += CacheItemsL*Rstride)
    {
      const TL *__restrict rL1 = &L1[j1];
      for (k1 = 0; k1 < (NR / CacheItemsR)*CacheItemsR; k1 += CacheItemsR)
      {
        const TR *__restrict rR1 = &R1[k1];
        rC = &C[k1];
        for (j2 = 0; j2 < CacheItemsL; ++j2, rR1 += Rstride)
        {
          const TL& Lv = rL1[j2];
          for (k2 = 0; k2 < CacheItemsR; ++k2)
            rC[k2] += sgn*(Lv*rR1[k2]);
        }
      }

      if ( NR % CacheItemsR )
      {
        k1 = (NR / CacheItemsR)*CacheItemsR;
        const TR *__restrict rR1 = &R1[k1];
        rC = &C[k1];
        for (j2 = 0; j2 < CacheItemsL; ++j2, rR1 += Rstride)
        {
          const TL& Lv = rL1[j2];
          for (k2 = 0; k2 < NR % CacheItemsR; ++k2)
            rC[k2] += sgn*(Lv*rR1[k2]);
        }
      }
    }
//---------------

    if ( NL % CacheItemsL )
    {
      for (j1 = (NL/CacheItemsL)*CacheItemsL; j1 < NL; ++j1, R1 += Rstride)
      {
        const TL& Lv = L1[j1];
        for (k1 = 0; k1 < (NR / CacheItemsR)*CacheItemsR; k1 += CacheItemsR)
        {
          const TR *__restrict rR1 = &R1[k1];
          rC = &C[k1];
          for (k2 = 0; k2 < CacheItemsR; ++k2)
            rC[k2] += sgn*(Lv*rR1[k2]);
        }

        if ( NR % CacheItemsR )
        {
          k1 = (NR / CacheItemsR)*CacheItemsR;
          const TR *__restrict rR1 = &R1[k1];
          rC = &C[k1];
          for (k2 = 0; k2 < NR % CacheItemsR; ++k2)
            rC[k2] += sgn*(Lv*rR1[k2]);
        }
      }
    }
  }
}

} //namespace DLA
} // namespace Metris

#define METRIS_MATRIXS_MATMUL_NATIVE(ML, NL,  MR, NR,  TL, TR, S, T) \
  template void MatrixS_MatMul_Native<TL, TR, S, T>::plus(const MatrixS< ML, NL, TL >& ml, const MatrixS< MR, NR, TR >& mr, \
                                                          const S& sgn, MatrixS< ML, NR, T >& res )

#endif // METRIS_DLA_MATRIXS_MATMUL_NATIVE_IMPL_H
