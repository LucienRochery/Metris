//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_CT_LOOP__
#define __METRIS_CT_LOOP__

#include <boost/hana.hpp> 
#include "inc_hana.hxx"

/*
Use these for compile time looping. Included/excluded bounds. Example: 
CT_FOR0_INC(1,nnode,ivar){
  // ... use ivar as a compile-time constant
}CT_FOR1(ivar);

Guidelines: 
 - Minimize usage. If a routine needs a given variable several times as a 
 compile-time constant, only one loop should be used or the variable should be 
 passed as a template param.
 - Minimize loop size. Code is pasted around; place as close to where needed. 
*/
#define CT_FOR0_INC(i1,i2,var) hana::while_(hana::less_equal.than(hana::int_c<i2>), hana::int_c<i1>, [&](auto c_##var){\
  constexpr int var = c_##var;
#define CT_FOR0_EXC(i1,i2,var) hana::while_(hana::less.than(hana::int_c<i2>), hana::int_c<i1>, [&](auto c_##var){\
  constexpr int var = c_##var;
#define CT_FOR1(var) return c_##var+1_c;})
#define CT_CONTINUE(var) return c_##var+1_c;


#endif
