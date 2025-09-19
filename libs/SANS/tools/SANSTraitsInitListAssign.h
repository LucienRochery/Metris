// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SANSTRAITSINITLISTASSING_H
#define SANSTRAITSINITLISTASSING_H

#ifdef __INTEL_COMPILER
//Maybe someday the intel compiler will get fixed and we won't need any of this mess...

#include "tools/SANSnumerics.h"   // Real
#include <initializer_list>

namespace SANS
{

//Used to assign list initializer to special types that can take the list initalizer. Otherwise it does nothing
template<class T>
struct initializer_list_assign
{
  typedef T type;
  template< class U >
  initializer_list_assign(T& val, const std::initializer_list<U>& s) {}
};

}  // namespace SANS

#endif

#endif  // SANSTRAITSINITLISTASSING_H
