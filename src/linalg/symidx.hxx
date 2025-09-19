//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_SYMIDX__
#define __METRIS_SYMIDX__

#include "utils/aux_pp_inc.hxx"

namespace Metris{



// Constexpr here just allows these functions to be used in a constexpr context
inline constexpr int sym2idx(int i, int j){
  return i > j ? (i*(i+1))/2 + j : (j*(j+1))/2 + i;
}


inline constexpr int sym3idx(int i, int j, int k){
  if(i >= j && i >= k){
    return (i*(i+1)*(i+2))/6 + sym2idx(j,k);
  }else if(j >= i && j>= k){
    return sym3idx(j,k,i);
  }else if(k >= i && k >= j){
    return sym3idx(k,i,j);
  }
  return -1;
}

// The rest follows recursively as needed. 
// We might even be able to generate using macros up to MAXDEG + 1
// Or else generate in a codegen script (preferably)



} // End namespace
#endif