// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef METRIS_TRAITS_INIT_LIST_ASSIGN_H
#define METRIS_TRAITS_INIT_LIST_ASSIGN_H

#ifdef __INTEL_COMPILER
//Maybe someday the intel compiler will get fixed and we won't need any of this mess...

#include "SANSnumerics.h"   // Metris::Real
#include <initializer_list>

namespace Metris
{

//Used to assign list initializer to special types that can take the list initalizer. Otherwise it does nothing
template<class T>
struct initializer_list_assign
{
  typedef T type;
  template< class U >
  initializer_list_assign(T& val, const std::initializer_list<U>& s) {}
};

}  // namespace Metris

#endif

#endif  // METRIS_TRAITS_INIT_LIST_ASSIGN_H
