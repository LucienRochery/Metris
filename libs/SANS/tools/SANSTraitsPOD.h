// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef METRIS_TRAITS_POD_H
#define METRIS_TRAITS_POD_H

#include "SANSnumerics.h"   // Metris::Real
#include "SANSException.h"

namespace Metris
{

// Traits class to allow templated classes to allow Metris::Real arguments in arithmetic
// operations; examples of classes that use this are DenseMatrix<MatrixS<N,Metris::Real>>
// and VectorS<N, SurrealS<M>>

struct MetrisDummyType
{
  operator int() const { METRIS_SUPPORT_ASSERT(false); return 0; }
};

template< class T >
struct POD { typedef Metris::Real type; };

template<>
struct POD<Metris::Real> { typedef MetrisDummyType type; };

template<>
struct POD<int> { typedef MetrisDummyType type; };

template<>
struct POD<unsigned int> { typedef MetrisDummyType type; };

}  // namespace Metris

#endif  // METRIS_TRAITS_POD_H
