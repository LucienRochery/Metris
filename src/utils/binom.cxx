//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include <cstdio>
#include <array>
#include "aux_misc.hxx"
#include "../ho_constants.hxx"

namespace Metris{


// Binomial coefficient k among n = n!/(k!(n-k)!)
// Occurs in Bernstein computations (e.g. prod) and much less risk of overflow
// this way than computing factorials. 
int binom(int k, int n){
  static std::array<std::array<int, METRIS_MAX_DEG_ORDERING>,
                                    METRIS_MAX_DEG_ORDERING > cache({{}}); 
  if(cache[k][n] == 0){
    if(k == 0 || k == n) cache[n-k][n] = cache[k][n] = 1;
    else                 cache[n-k][n] = cache[k][n] = binom(k-1,n-1) + binom(k,n-1);
  }
  return cache[k][n];
}


// (i,j,k) among i+j+k
int multinom(int i, int j, int k){
  static std::array<std::array<int, getnnod2(METRIS_MAX_DEG_ORDERING)>,
                                    METRIS_MAX_DEG_ORDERING > cache({{}});
  int inode = mul2nod(i,j,k);
  if(cache[i+j+k][inode] == 0){
         if(i == 0) cache[i+j+k][inode] = binom(j,j+k);
    else if(j == 0) cache[i+j+k][inode] = binom(k,k+i);
    else if(k == 0) cache[i+j+k][inode] = binom(i,i+j);
    else            cache[i+j+k][inode] = multinom(i-1,j,k)  
                                        + multinom(i,j-1,k) 
                                        + multinom(i,j,k-1);
  }
  return cache[i+j+k][inode];
}



}