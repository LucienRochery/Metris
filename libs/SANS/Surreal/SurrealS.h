// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef METRIS_SURREALS_H
#define METRIS_SURREALS_H

#include <boost/type_traits/is_arithmetic.hpp>

#include "SurrealS_Lazy.h"

//Let Surreal be treated as a arithmetic type
namespace boost
{
template<int N, class T>
struct is_arithmetic< Metris::SurrealS<N, T> > : boost::true_type {};
}


#endif // METRIS_SURREALS_H
