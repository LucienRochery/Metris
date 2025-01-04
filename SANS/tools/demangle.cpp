// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "demangle.h"

#ifdef __GNUG__

#include <cstdlib>
#include <memory>
#include <cxxabi.h>

namespace SANS
{

std::string demangle(const char* name)
{
  int status = -4; // some arbitrary value to eliminate the compiler warning

  std::unique_ptr<char, void(*)(void*)> res(
    abi::__cxa_demangle(name, NULL, NULL, &status), std::free );

  return (status==0) ? res.get() : name;
}

}// namespace SANS

#else

namespace SANS
{

// does nothing if not g++
std::string demangle(const char* name)
{
  return name;
}

}// namespace SANS

#endif
