//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php




#ifndef __METRIS_SYM3IDX__
#define __METRIS_SYM3IDX__

//#include <array>
namespace Metris{

//constexpr std::array<std::array<int,3>,3> sym2idx{[]() constexpr{
//  std::array<std::array<int,3>,3> ret{};
//  ret[0][0] = 0;
//  ret[0][1] = ret[1][0] = 1;
//  ret[1][1] = 2;
//  ret[0][2] = ret[2][0] = 3;
//  ret[1][2] = ret[2][1] = 4;
//  ret[2][2] = 5;
//  return ret;
//}()};


// Constexpr here just allows these functions to be used in a constexpr context
inline constexpr int sym2idx(int i, int j){
  return i > j ? (i*(i+1))/2 + j : (j*(j+1))/2 + i;
}


inline constexpr int sym3idx(int i, int j, int k){
  if(i >= j && i >= k){
    if(j >= k){
      return (i*(i+1)*(i+2))/6 + (j*(j+1))/2 + k;
    }else{
      return (i*(i+1)*(i+2))/6 + (k*(k+1))/2 + j;
    }
  }else if(j >= i && j>= k){
    return sym3idx(j,k,i);
  }else if(k >= i && k >= j){
    return sym3idx(k,i,j);
  }
  return -1;
}

//// We could have used [] but C++23 needed
//class sym2idx_C{
//public:
//   int operator()  (int i, int j) const{
//  // N(i) + j
//  // N(i) is the total number of indices before (i,0)
//  // this is simply sum_1^i k  = (i(i+1))/2
//  // if j >= i, invert
//    return i > j ? (i*(i+1))/2 + j : (j*(j+1))/2 + i;
//  }
//};
//
//constexpr const sym2idx_C sym2idx;

} // End namespace
#endif