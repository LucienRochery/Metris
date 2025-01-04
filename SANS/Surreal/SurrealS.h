// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SANS_SURREALS_H
#define SANS_SURREALS_H

#include <boost/type_traits/is_arithmetic.hpp>

////#ifdef SURREAL_RVO
////#include "SurrealS_RVO.h"
//#if defined(SURREAL_TRAD)
//#include "SurrealS_Trad.h"
//#elif defined(SURREAL_LAZY)
#include "SANS//Surreal/SurrealS_Lazy.h"
//#elif defined(SURREAL_REVERSE)
//#include "SurrealS_Reverse.h"
//#else
//#error "Please define SURREAL_TRAD, SURREAL_LAZY or SURREAL_REVERSE"
//#endif

//Let Surreal be treated as a arithmetic type
namespace boost
{
template<int N, class T>
struct is_arithmetic< SANS::SurrealS<N, T> > : boost::true_type {};
}


#endif // SANS_SURREALS_H
