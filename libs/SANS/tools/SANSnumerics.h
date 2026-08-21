// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef METRIS_NUMERICS_H
#define METRIS_NUMERICS_H

#ifndef METRIS_ALWAYS_INLINE
// METRIS_ALWAYS_INLINE is a macro to further encourage the compiler to inline a function
#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__clang__)
#define METRIS_ALWAYS_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define METRIS_ALWAYS_INLINE __forceinline
#else
#warning Not forcing inline with this compiler... (Please add this compiler to tools/always_inline.h)
#define METRIS_ALWAYS_INLINE inline
#endif
#endif


namespace Metris
{

typedef double Real;

typedef long LongInt;

// math constants
static const Real PI = 3.14159265358979323846;

} // namespace Metris

#endif  // METRIS_NUMERICS_H
