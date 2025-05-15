//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_BINOM__
#define __METRIS_BINOM__

namespace Metris{


//// Binomial coefficient n among m = m!/(n!(m-n)!)
//// Occurs in Bernstein computations (e.g. prod) and much less risk of overflow
//// this way than computing factorials. 
//template<int n, int m>
//inline constexpr int binom(){
//  //B(n, k) = B(n − 1, k − 1) + B(n − 1, k)
//  if constexpr(n == m){
//    return 1;
//  }else if(n == 0){
//    return 1;
//  }else{
//    return binom<n-1,m-1>() + binom<n,m-1>();
//  }
//}
//
//template<int m>
//inline constexpr int binom<0,m>(){
//  return 1;
//}


struct ct_Binom{
  constexpr ct_Binom(): get(){
    for(int ii = 0; ii <= METRIS_MAX_DEG_ORDERING; ii++)
      for(int jj = 0; jj <= METRIS_MAX_DEG_ORDERING; jj++)
        get[ii][jj] = -1;

    get[0][0] = 1;
    for(int ii = 1; ii <= METRIS_MAX_DEG_ORDERING; ii++){
      get[ii][0] = 1;
      for(int jj = 1; jj < ii; jj++)
        get[ii][jj] = get[ii-1][jj] + get[ii-1][jj-1];
      get[ii][ii] = 1;
    }
  }
  // [n][k] -> n!/k!(n-k)!
  int get[1+METRIS_MAX_DEG_ORDERING][METRIS_MAX_DEG_ORDERING+1];
};

constexpr const ct_Binom ct_binom;

int binom(int m, int n);


}
#endif