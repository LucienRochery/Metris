//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php




#ifndef __METRIS_SYM3IDX__
#define __METRIS_SYM3IDX__

//#include <array>
namespace Metris{

//constexpr std::array<std::array<int,3>,3> sym3idx{[]() constexpr{
//  std::array<std::array<int,3>,3> ret{};
//  ret[0][0] = 0;
//  ret[0][1] = ret[1][0] = 1;
//  ret[1][1] = 2;
//  ret[0][2] = ret[2][0] = 3;
//  ret[1][2] = ret[2][1] = 4;
//  ret[2][2] = 5;
//  return ret;
//}()};

inline int sym3idx(int i, int j){
  return i > j ? (i*(i+1))/2 + j : (j*(j+1))/2 + i;
}

//// We could have used [] but C++23 needed
//class sym3idx_C{
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
//constexpr const sym3idx_C sym3idx;

} // End namespace
#endif