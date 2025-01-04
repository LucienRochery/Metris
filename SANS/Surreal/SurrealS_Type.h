// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SURREALS_TYPE_H
#define SURREALS_TYPE_H

#include  <type_traits>

#include "SANS/tools/SANSnumerics.h"
//#include "SANS//Surreal/always_inline.h"

#include "SANS//Surreal/SurrealS_fwd.h"


namespace SANS
{

// Forward declaration of SurrealS types

class SurrealSTypeBase {};

template< class Derived, class T >
struct SurrealSType : SurrealSTypeBase
{
  //A convenient method for casting to the derived type
  ALWAYS_INLINE const Derived& cast() const { return static_cast<const Derived&>(*this); }

  //A simple way to call value without having to case first
  ALWAYS_INLINE T value() const { return cast().value(); }
};

template<class T>
struct is_arithmetic_not_SurrealS
{
  static const bool value = 
   (std::is_arithmetic<T>::value) //|| 
   	//boost::is_arithmetic<T>::value ||
   	//boost::multiprecision::detail::is_arithmetic<T>::value)
   && !std::is_base_of<SurrealSTypeBase, T>::value;
};

} // namespace SANS

#endif //SURREALS_TYPE_H
