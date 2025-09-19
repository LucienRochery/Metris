// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef ALWAYS_INLINE
// ALWAYS_INLINE is a macro to further encourage the compiler to inline a function
#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__clang__)
#define ALWAYS_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define ALWAYS_INLINE __forceinline
#else
#warning Not forcing inline with this compiler... (Please add this compiler to tools/always_inline.h)
#define ALWAYS_INLINE inline
#endif
#endif


#ifndef SANSNUMERICS_H
#define SANSNUMERICS_H



typedef double Real;

typedef long LongInt;

// math constants
static const Real PI = 3.14159265358979323846;

#endif  // SANSNUMERICS_H
