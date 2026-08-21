// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef METRIS_MINMAX_H
#define METRIS_MINMAX_H

#include <algorithm>

//Some compilers make min/max macros...
#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#ifndef METRIS_MAX
#define METRIS_MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef METRIS_MIN
#define METRIS_MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

// Put the standard min/max functions in the Metris namespace.
namespace Metris
{
using std::min;
using std::max;
}

#if 0
#include <type_traits>
#include <boost/type_traits/is_arithmetic.hpp>

namespace Metris
{
template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type&
min(const T& a, const T& b) { return METRIS_MIN(a,b); }

template<typename T>
inline typename std::enable_if< boost::is_arithmetic<T>::value, T>::type&
max(const T& a, const T& b) { return METRIS_MAX(a,b); }
}
#endif

#endif //METRIS_MINMAX_H
