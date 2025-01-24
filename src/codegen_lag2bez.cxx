//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "codegen_lag2bez.hxx"

namespace Metris{

template<> void lag2bez1<0,1>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){}

template<> void lag2bez1<1,1>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  rfld1[lfld[  0]][0] = rfld0[lfld[  0]][0];
  rfld1[lfld[  1]][0] = rfld0[lfld[  1]][0];
}

template<> void lag2bez1<1,2>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 2; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
  }
}

template<> void lag2bez1<1,3>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 3; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
  }
}

template<> void lag2bez1<1,6>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 6; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
  }
}

template<> void lag2bez1<2,1>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  rfld1[lfld[  0]][0] = rfld0[lfld[  0]][0];
  rfld1[lfld[  1]][0] = rfld0[lfld[  1]][0];
  rfld1[lfld[  2]][0] = -1*rfld0[lfld[  0]][0]/2 + -1*rfld0[lfld[  1]][0]/2 + 4*rfld0[lfld[  2]][0]/2;
}

template<> void lag2bez1<2,2>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 2; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = -1*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + 4*rfld0[lfld[  2]][i]/2;
  }
}

template<> void lag2bez1<2,3>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 3; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = -1*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + 4*rfld0[lfld[  2]][i]/2;
  }
}

template<> void lag2bez1<2,6>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 6; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = -1*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + 4*rfld0[lfld[  2]][i]/2;
  }
}

template<> void lag2bez1<3,1>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  rfld1[lfld[  0]][0] = rfld0[lfld[  0]][0];
  rfld1[lfld[  1]][0] = rfld0[lfld[  1]][0];
  rfld1[lfld[  2]][0] = -5*rfld0[lfld[  0]][0]/6 + 1*rfld0[lfld[  1]][0]/3 + 6*rfld0[lfld[  2]][0]/2 + -3*rfld0[lfld[  3]][0]/2;
  rfld1[lfld[  3]][0] = 1*rfld0[lfld[  0]][0]/3 + -5*rfld0[lfld[  1]][0]/6 + -3*rfld0[lfld[  2]][0]/2 + 6*rfld0[lfld[  3]][0]/2;
}

template<> void lag2bez1<3,2>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 2; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = -5*rfld0[lfld[  0]][i]/6 + 1*rfld0[lfld[  1]][i]/3 + 6*rfld0[lfld[  2]][i]/2 + -3*rfld0[lfld[  3]][i]/2;
    rfld1[lfld[  3]][i] = 1*rfld0[lfld[  0]][i]/3 + -5*rfld0[lfld[  1]][i]/6 + -3*rfld0[lfld[  2]][i]/2 + 6*rfld0[lfld[  3]][i]/2;
  }
}

template<> void lag2bez1<3,3>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 3; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = -5*rfld0[lfld[  0]][i]/6 + 1*rfld0[lfld[  1]][i]/3 + 6*rfld0[lfld[  2]][i]/2 + -3*rfld0[lfld[  3]][i]/2;
    rfld1[lfld[  3]][i] = 1*rfld0[lfld[  0]][i]/3 + -5*rfld0[lfld[  1]][i]/6 + -3*rfld0[lfld[  2]][i]/2 + 6*rfld0[lfld[  3]][i]/2;
  }
}

template<> void lag2bez1<3,6>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 6; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = -5*rfld0[lfld[  0]][i]/6 + 1*rfld0[lfld[  1]][i]/3 + 6*rfld0[lfld[  2]][i]/2 + -3*rfld0[lfld[  3]][i]/2;
    rfld1[lfld[  3]][i] = 1*rfld0[lfld[  0]][i]/3 + -5*rfld0[lfld[  1]][i]/6 + -3*rfld0[lfld[  2]][i]/2 + 6*rfld0[lfld[  3]][i]/2;
  }
}

template<> void lag2bez2<0,1>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){}

template<> void lag2bez2<1,1>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  rfld1[lfld[  0]][0] = rfld0[lfld[  0]][0];
  rfld1[lfld[  1]][0] = rfld0[lfld[  1]][0];
  rfld1[lfld[  2]][0] = rfld0[lfld[  2]][0];
}

template<> void lag2bez2<1,2>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 2; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
  }
}

template<> void lag2bez2<1,3>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 3; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
  }
}

template<> void lag2bez2<1,6>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 6; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
  }
}

template<> void lag2bez2<2,1>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  rfld1[lfld[  0]][0] = rfld0[lfld[  0]][0];
  rfld1[lfld[  1]][0] = rfld0[lfld[  1]][0];
  rfld1[lfld[  2]][0] = rfld0[lfld[  2]][0];
  rfld1[lfld[  3]][0] = 0*rfld0[lfld[  0]][0]/2 + -1*rfld0[lfld[  1]][0]/2 + -1*rfld0[lfld[  2]][0]/2 + 4*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2;
  rfld1[lfld[  4]][0] = -1*rfld0[lfld[  0]][0]/2 + 0*rfld0[lfld[  1]][0]/2 + -1*rfld0[lfld[  2]][0]/2 + 0*rfld0[lfld[  3]][0]/2 + 4*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2;
  rfld1[lfld[  5]][0] = -1*rfld0[lfld[  0]][0]/2 + -1*rfld0[lfld[  1]][0]/2 + 0*rfld0[lfld[  2]][0]/2 + 0*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + 4*rfld0[lfld[  5]][0]/2;
}

template<> void lag2bez2<2,2>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 2; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = 0*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + 4*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2;
    rfld1[lfld[  4]][i] = -1*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 4*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2;
    rfld1[lfld[  5]][i] = -1*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 4*rfld0[lfld[  5]][i]/2;
  }
}

template<> void lag2bez2<2,3>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 3; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = 0*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + 4*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2;
    rfld1[lfld[  4]][i] = -1*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 4*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2;
    rfld1[lfld[  5]][i] = -1*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 4*rfld0[lfld[  5]][i]/2;
  }
}

template<> void lag2bez2<2,6>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 6; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = 0*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + 4*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2;
    rfld1[lfld[  4]][i] = -1*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 4*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2;
    rfld1[lfld[  5]][i] = -1*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 4*rfld0[lfld[  5]][i]/2;
  }
}

template<> void lag2bez2<3,1>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  rfld1[lfld[  0]][0] = rfld0[lfld[  0]][0];
  rfld1[lfld[  1]][0] = rfld0[lfld[  1]][0];
  rfld1[lfld[  2]][0] = rfld0[lfld[  2]][0];
  rfld1[lfld[  3]][0] = 0*rfld0[lfld[  0]][0]/2 + -5*rfld0[lfld[  1]][0]/6 + 1*rfld0[lfld[  2]][0]/3 + 6*rfld0[lfld[  3]][0]/2 + -3*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2;
  rfld1[lfld[  4]][0] = 0*rfld0[lfld[  0]][0]/2 + 1*rfld0[lfld[  1]][0]/3 + -5*rfld0[lfld[  2]][0]/6 + -3*rfld0[lfld[  3]][0]/2 + 6*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2;
  rfld1[lfld[  5]][0] = 1*rfld0[lfld[  0]][0]/3 + 0*rfld0[lfld[  1]][0]/2 + -5*rfld0[lfld[  2]][0]/6 + 0*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + 6*rfld0[lfld[  5]][0]/2 + -3*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2;
  rfld1[lfld[  6]][0] = -5*rfld0[lfld[  0]][0]/6 + 0*rfld0[lfld[  1]][0]/2 + 1*rfld0[lfld[  2]][0]/3 + 0*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + -3*rfld0[lfld[  5]][0]/2 + 6*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2;
  rfld1[lfld[  7]][0] = -5*rfld0[lfld[  0]][0]/6 + 1*rfld0[lfld[  1]][0]/3 + 0*rfld0[lfld[  2]][0]/2 + 0*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 6*rfld0[lfld[  7]][0]/2 + -3*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2;
  rfld1[lfld[  8]][0] = 1*rfld0[lfld[  0]][0]/3 + -5*rfld0[lfld[  1]][0]/6 + 0*rfld0[lfld[  2]][0]/2 + 0*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + -3*rfld0[lfld[  7]][0]/2 + 6*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2;
  rfld1[lfld[  9]][0] = 1*rfld0[lfld[  0]][0]/3 + 1*rfld0[lfld[  1]][0]/3 + 1*rfld0[lfld[  2]][0]/3 + -3*rfld0[lfld[  3]][0]/4 + -3*rfld0[lfld[  4]][0]/4 + -3*rfld0[lfld[  5]][0]/4 + -3*rfld0[lfld[  6]][0]/4 + -3*rfld0[lfld[  7]][0]/4 + -3*rfld0[lfld[  8]][0]/4 + 9*rfld0[lfld[  9]][0]/2;
}

template<> void lag2bez2<3,2>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 2; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = 0*rfld0[lfld[  0]][i]/2 + -5*rfld0[lfld[  1]][i]/6 + 1*rfld0[lfld[  2]][i]/3 + 6*rfld0[lfld[  3]][i]/2 + -3*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  4]][i] = 0*rfld0[lfld[  0]][i]/2 + 1*rfld0[lfld[  1]][i]/3 + -5*rfld0[lfld[  2]][i]/6 + -3*rfld0[lfld[  3]][i]/2 + 6*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  5]][i] = 1*rfld0[lfld[  0]][i]/3 + 0*rfld0[lfld[  1]][i]/2 + -5*rfld0[lfld[  2]][i]/6 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 6*rfld0[lfld[  5]][i]/2 + -3*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  6]][i] = -5*rfld0[lfld[  0]][i]/6 + 0*rfld0[lfld[  1]][i]/2 + 1*rfld0[lfld[  2]][i]/3 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + -3*rfld0[lfld[  5]][i]/2 + 6*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  7]][i] = -5*rfld0[lfld[  0]][i]/6 + 1*rfld0[lfld[  1]][i]/3 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 6*rfld0[lfld[  7]][i]/2 + -3*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  8]][i] = 1*rfld0[lfld[  0]][i]/3 + -5*rfld0[lfld[  1]][i]/6 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + -3*rfld0[lfld[  7]][i]/2 + 6*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  9]][i] = 1*rfld0[lfld[  0]][i]/3 + 1*rfld0[lfld[  1]][i]/3 + 1*rfld0[lfld[  2]][i]/3 + -3*rfld0[lfld[  3]][i]/4 + -3*rfld0[lfld[  4]][i]/4 + -3*rfld0[lfld[  5]][i]/4 + -3*rfld0[lfld[  6]][i]/4 + -3*rfld0[lfld[  7]][i]/4 + -3*rfld0[lfld[  8]][i]/4 + 9*rfld0[lfld[  9]][i]/2;
  }
}

template<> void lag2bez2<3,3>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 3; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = 0*rfld0[lfld[  0]][i]/2 + -5*rfld0[lfld[  1]][i]/6 + 1*rfld0[lfld[  2]][i]/3 + 6*rfld0[lfld[  3]][i]/2 + -3*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  4]][i] = 0*rfld0[lfld[  0]][i]/2 + 1*rfld0[lfld[  1]][i]/3 + -5*rfld0[lfld[  2]][i]/6 + -3*rfld0[lfld[  3]][i]/2 + 6*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  5]][i] = 1*rfld0[lfld[  0]][i]/3 + 0*rfld0[lfld[  1]][i]/2 + -5*rfld0[lfld[  2]][i]/6 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 6*rfld0[lfld[  5]][i]/2 + -3*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  6]][i] = -5*rfld0[lfld[  0]][i]/6 + 0*rfld0[lfld[  1]][i]/2 + 1*rfld0[lfld[  2]][i]/3 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + -3*rfld0[lfld[  5]][i]/2 + 6*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  7]][i] = -5*rfld0[lfld[  0]][i]/6 + 1*rfld0[lfld[  1]][i]/3 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 6*rfld0[lfld[  7]][i]/2 + -3*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  8]][i] = 1*rfld0[lfld[  0]][i]/3 + -5*rfld0[lfld[  1]][i]/6 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + -3*rfld0[lfld[  7]][i]/2 + 6*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  9]][i] = 1*rfld0[lfld[  0]][i]/3 + 1*rfld0[lfld[  1]][i]/3 + 1*rfld0[lfld[  2]][i]/3 + -3*rfld0[lfld[  3]][i]/4 + -3*rfld0[lfld[  4]][i]/4 + -3*rfld0[lfld[  5]][i]/4 + -3*rfld0[lfld[  6]][i]/4 + -3*rfld0[lfld[  7]][i]/4 + -3*rfld0[lfld[  8]][i]/4 + 9*rfld0[lfld[  9]][i]/2;
  }
}

template<> void lag2bez2<3,6>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 6; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = 0*rfld0[lfld[  0]][i]/2 + -5*rfld0[lfld[  1]][i]/6 + 1*rfld0[lfld[  2]][i]/3 + 6*rfld0[lfld[  3]][i]/2 + -3*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  4]][i] = 0*rfld0[lfld[  0]][i]/2 + 1*rfld0[lfld[  1]][i]/3 + -5*rfld0[lfld[  2]][i]/6 + -3*rfld0[lfld[  3]][i]/2 + 6*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  5]][i] = 1*rfld0[lfld[  0]][i]/3 + 0*rfld0[lfld[  1]][i]/2 + -5*rfld0[lfld[  2]][i]/6 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 6*rfld0[lfld[  5]][i]/2 + -3*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  6]][i] = -5*rfld0[lfld[  0]][i]/6 + 0*rfld0[lfld[  1]][i]/2 + 1*rfld0[lfld[  2]][i]/3 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + -3*rfld0[lfld[  5]][i]/2 + 6*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  7]][i] = -5*rfld0[lfld[  0]][i]/6 + 1*rfld0[lfld[  1]][i]/3 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 6*rfld0[lfld[  7]][i]/2 + -3*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  8]][i] = 1*rfld0[lfld[  0]][i]/3 + -5*rfld0[lfld[  1]][i]/6 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + -3*rfld0[lfld[  7]][i]/2 + 6*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  9]][i] = 1*rfld0[lfld[  0]][i]/3 + 1*rfld0[lfld[  1]][i]/3 + 1*rfld0[lfld[  2]][i]/3 + -3*rfld0[lfld[  3]][i]/4 + -3*rfld0[lfld[  4]][i]/4 + -3*rfld0[lfld[  5]][i]/4 + -3*rfld0[lfld[  6]][i]/4 + -3*rfld0[lfld[  7]][i]/4 + -3*rfld0[lfld[  8]][i]/4 + 9*rfld0[lfld[  9]][i]/2;
  }
}

template<> void lag2bez3<0,1>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){}

template<> void lag2bez3<1,1>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  rfld1[lfld[  0]][0] = rfld0[lfld[  0]][0];
  rfld1[lfld[  1]][0] = rfld0[lfld[  1]][0];
  rfld1[lfld[  2]][0] = rfld0[lfld[  2]][0];
  rfld1[lfld[  3]][0] = rfld0[lfld[  3]][0];
}

template<> void lag2bez3<1,2>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 2; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = rfld0[lfld[  3]][i];
  }
}

template<> void lag2bez3<1,3>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 3; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = rfld0[lfld[  3]][i];
  }
}

template<> void lag2bez3<1,6>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 6; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = rfld0[lfld[  3]][i];
  }
}

template<> void lag2bez3<2,1>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  rfld1[lfld[  0]][0] = rfld0[lfld[  0]][0];
  rfld1[lfld[  1]][0] = rfld0[lfld[  1]][0];
  rfld1[lfld[  2]][0] = rfld0[lfld[  2]][0];
  rfld1[lfld[  3]][0] = rfld0[lfld[  3]][0];
  rfld1[lfld[  4]][0] = -1*rfld0[lfld[  0]][0]/2 + -1*rfld0[lfld[  1]][0]/2 + 0*rfld0[lfld[  2]][0]/2 + 0*rfld0[lfld[  3]][0]/2 + 4*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2;
  rfld1[lfld[  5]][0] = 0*rfld0[lfld[  0]][0]/2 + -1*rfld0[lfld[  1]][0]/2 + -1*rfld0[lfld[  2]][0]/2 + 0*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + 4*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2;
  rfld1[lfld[  6]][0] = -1*rfld0[lfld[  0]][0]/2 + 0*rfld0[lfld[  1]][0]/2 + -1*rfld0[lfld[  2]][0]/2 + 0*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 4*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2;
  rfld1[lfld[  7]][0] = -1*rfld0[lfld[  0]][0]/2 + 0*rfld0[lfld[  1]][0]/2 + 0*rfld0[lfld[  2]][0]/2 + -1*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 4*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2;
  rfld1[lfld[  8]][0] = 0*rfld0[lfld[  0]][0]/2 + -1*rfld0[lfld[  1]][0]/2 + 0*rfld0[lfld[  2]][0]/2 + -1*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 4*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2;
  rfld1[lfld[  9]][0] = 0*rfld0[lfld[  0]][0]/2 + 0*rfld0[lfld[  1]][0]/2 + -1*rfld0[lfld[  2]][0]/2 + -1*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 4*rfld0[lfld[  9]][0]/2;
}

template<> void lag2bez3<2,2>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 2; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = rfld0[lfld[  3]][i];
    rfld1[lfld[  4]][i] = -1*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 4*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  5]][i] = 0*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 4*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  6]][i] = -1*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 4*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  7]][i] = -1*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + -1*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 4*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  8]][i] = 0*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + -1*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 4*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  9]][i] = 0*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + -1*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 4*rfld0[lfld[  9]][i]/2;
  }
}

template<> void lag2bez3<2,3>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 3; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = rfld0[lfld[  3]][i];
    rfld1[lfld[  4]][i] = -1*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 4*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  5]][i] = 0*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 4*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  6]][i] = -1*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 4*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  7]][i] = -1*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + -1*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 4*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  8]][i] = 0*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + -1*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 4*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  9]][i] = 0*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + -1*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 4*rfld0[lfld[  9]][i]/2;
  }
}

template<> void lag2bez3<2,6>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 6; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = rfld0[lfld[  3]][i];
    rfld1[lfld[  4]][i] = -1*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 4*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  5]][i] = 0*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 4*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  6]][i] = -1*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 4*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  7]][i] = -1*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + -1*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 4*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  8]][i] = 0*rfld0[lfld[  0]][i]/2 + -1*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + -1*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 4*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2;
    rfld1[lfld[  9]][i] = 0*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + -1*rfld0[lfld[  2]][i]/2 + -1*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 4*rfld0[lfld[  9]][i]/2;
  }
}

template<> void lag2bez3<3,1>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  rfld1[lfld[  0]][0] = rfld0[lfld[  0]][0];
  rfld1[lfld[  1]][0] = rfld0[lfld[  1]][0];
  rfld1[lfld[  2]][0] = rfld0[lfld[  2]][0];
  rfld1[lfld[  3]][0] = rfld0[lfld[  3]][0];
  rfld1[lfld[  4]][0] = -5*rfld0[lfld[  0]][0]/6 + 1*rfld0[lfld[  1]][0]/3 + 0*rfld0[lfld[  2]][0]/2 + 0*rfld0[lfld[  3]][0]/2 + 6*rfld0[lfld[  4]][0]/2 + -3*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2 + 0*rfld0[lfld[ 10]][0]/2 + 0*rfld0[lfld[ 11]][0]/2 + 0*rfld0[lfld[ 12]][0]/2 + 0*rfld0[lfld[ 13]][0]/2 + 0*rfld0[lfld[ 14]][0]/2 + 0*rfld0[lfld[ 15]][0]/2 + 0*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[  5]][0] = 1*rfld0[lfld[  0]][0]/3 + -5*rfld0[lfld[  1]][0]/6 + 0*rfld0[lfld[  2]][0]/2 + 0*rfld0[lfld[  3]][0]/2 + -3*rfld0[lfld[  4]][0]/2 + 6*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2 + 0*rfld0[lfld[ 10]][0]/2 + 0*rfld0[lfld[ 11]][0]/2 + 0*rfld0[lfld[ 12]][0]/2 + 0*rfld0[lfld[ 13]][0]/2 + 0*rfld0[lfld[ 14]][0]/2 + 0*rfld0[lfld[ 15]][0]/2 + 0*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[  6]][0] = 0*rfld0[lfld[  0]][0]/2 + -5*rfld0[lfld[  1]][0]/6 + 1*rfld0[lfld[  2]][0]/3 + 0*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 6*rfld0[lfld[  6]][0]/2 + -3*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2 + 0*rfld0[lfld[ 10]][0]/2 + 0*rfld0[lfld[ 11]][0]/2 + 0*rfld0[lfld[ 12]][0]/2 + 0*rfld0[lfld[ 13]][0]/2 + 0*rfld0[lfld[ 14]][0]/2 + 0*rfld0[lfld[ 15]][0]/2 + 0*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[  7]][0] = 0*rfld0[lfld[  0]][0]/2 + 1*rfld0[lfld[  1]][0]/3 + -5*rfld0[lfld[  2]][0]/6 + 0*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + -3*rfld0[lfld[  6]][0]/2 + 6*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2 + 0*rfld0[lfld[ 10]][0]/2 + 0*rfld0[lfld[ 11]][0]/2 + 0*rfld0[lfld[ 12]][0]/2 + 0*rfld0[lfld[ 13]][0]/2 + 0*rfld0[lfld[ 14]][0]/2 + 0*rfld0[lfld[ 15]][0]/2 + 0*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[  8]][0] = 1*rfld0[lfld[  0]][0]/3 + 0*rfld0[lfld[  1]][0]/2 + -5*rfld0[lfld[  2]][0]/6 + 0*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 6*rfld0[lfld[  8]][0]/2 + -3*rfld0[lfld[  9]][0]/2 + 0*rfld0[lfld[ 10]][0]/2 + 0*rfld0[lfld[ 11]][0]/2 + 0*rfld0[lfld[ 12]][0]/2 + 0*rfld0[lfld[ 13]][0]/2 + 0*rfld0[lfld[ 14]][0]/2 + 0*rfld0[lfld[ 15]][0]/2 + 0*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[  9]][0] = -5*rfld0[lfld[  0]][0]/6 + 0*rfld0[lfld[  1]][0]/2 + 1*rfld0[lfld[  2]][0]/3 + 0*rfld0[lfld[  3]][0]/2 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + -3*rfld0[lfld[  8]][0]/2 + 6*rfld0[lfld[  9]][0]/2 + 0*rfld0[lfld[ 10]][0]/2 + 0*rfld0[lfld[ 11]][0]/2 + 0*rfld0[lfld[ 12]][0]/2 + 0*rfld0[lfld[ 13]][0]/2 + 0*rfld0[lfld[ 14]][0]/2 + 0*rfld0[lfld[ 15]][0]/2 + 0*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[ 10]][0] = -5*rfld0[lfld[  0]][0]/6 + 0*rfld0[lfld[  1]][0]/2 + 0*rfld0[lfld[  2]][0]/2 + 1*rfld0[lfld[  3]][0]/3 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2 + 6*rfld0[lfld[ 10]][0]/2 + -3*rfld0[lfld[ 11]][0]/2 + 0*rfld0[lfld[ 12]][0]/2 + 0*rfld0[lfld[ 13]][0]/2 + 0*rfld0[lfld[ 14]][0]/2 + 0*rfld0[lfld[ 15]][0]/2 + 0*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[ 11]][0] = 1*rfld0[lfld[  0]][0]/3 + 0*rfld0[lfld[  1]][0]/2 + 0*rfld0[lfld[  2]][0]/2 + -5*rfld0[lfld[  3]][0]/6 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2 + -3*rfld0[lfld[ 10]][0]/2 + 6*rfld0[lfld[ 11]][0]/2 + 0*rfld0[lfld[ 12]][0]/2 + 0*rfld0[lfld[ 13]][0]/2 + 0*rfld0[lfld[ 14]][0]/2 + 0*rfld0[lfld[ 15]][0]/2 + 0*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[ 12]][0] = 0*rfld0[lfld[  0]][0]/2 + -5*rfld0[lfld[  1]][0]/6 + 0*rfld0[lfld[  2]][0]/2 + 1*rfld0[lfld[  3]][0]/3 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2 + 0*rfld0[lfld[ 10]][0]/2 + 0*rfld0[lfld[ 11]][0]/2 + 6*rfld0[lfld[ 12]][0]/2 + -3*rfld0[lfld[ 13]][0]/2 + 0*rfld0[lfld[ 14]][0]/2 + 0*rfld0[lfld[ 15]][0]/2 + 0*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[ 13]][0] = 0*rfld0[lfld[  0]][0]/2 + 1*rfld0[lfld[  1]][0]/3 + 0*rfld0[lfld[  2]][0]/2 + -5*rfld0[lfld[  3]][0]/6 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2 + 0*rfld0[lfld[ 10]][0]/2 + 0*rfld0[lfld[ 11]][0]/2 + -3*rfld0[lfld[ 12]][0]/2 + 6*rfld0[lfld[ 13]][0]/2 + 0*rfld0[lfld[ 14]][0]/2 + 0*rfld0[lfld[ 15]][0]/2 + 0*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[ 14]][0] = 0*rfld0[lfld[  0]][0]/2 + 0*rfld0[lfld[  1]][0]/2 + -5*rfld0[lfld[  2]][0]/6 + 1*rfld0[lfld[  3]][0]/3 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2 + 0*rfld0[lfld[ 10]][0]/2 + 0*rfld0[lfld[ 11]][0]/2 + 0*rfld0[lfld[ 12]][0]/2 + 0*rfld0[lfld[ 13]][0]/2 + 6*rfld0[lfld[ 14]][0]/2 + -3*rfld0[lfld[ 15]][0]/2 + 0*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[ 15]][0] = 0*rfld0[lfld[  0]][0]/2 + 0*rfld0[lfld[  1]][0]/2 + 1*rfld0[lfld[  2]][0]/3 + -5*rfld0[lfld[  3]][0]/6 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2 + 0*rfld0[lfld[ 10]][0]/2 + 0*rfld0[lfld[ 11]][0]/2 + 0*rfld0[lfld[ 12]][0]/2 + 0*rfld0[lfld[ 13]][0]/2 + -3*rfld0[lfld[ 14]][0]/2 + 6*rfld0[lfld[ 15]][0]/2 + 0*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[ 16]][0] = 0*rfld0[lfld[  0]][0]/2 + 1*rfld0[lfld[  1]][0]/3 + 1*rfld0[lfld[  2]][0]/3 + 1*rfld0[lfld[  3]][0]/3 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + -3*rfld0[lfld[  6]][0]/4 + -3*rfld0[lfld[  7]][0]/4 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2 + 0*rfld0[lfld[ 10]][0]/2 + 0*rfld0[lfld[ 11]][0]/2 + -3*rfld0[lfld[ 12]][0]/4 + -3*rfld0[lfld[ 13]][0]/4 + -3*rfld0[lfld[ 14]][0]/4 + -3*rfld0[lfld[ 15]][0]/4 + 9*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[ 17]][0] = 1*rfld0[lfld[  0]][0]/3 + 0*rfld0[lfld[  1]][0]/2 + 1*rfld0[lfld[  2]][0]/3 + 1*rfld0[lfld[  3]][0]/3 + 0*rfld0[lfld[  4]][0]/2 + 0*rfld0[lfld[  5]][0]/2 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + -3*rfld0[lfld[  8]][0]/4 + -3*rfld0[lfld[  9]][0]/4 + -3*rfld0[lfld[ 10]][0]/4 + -3*rfld0[lfld[ 11]][0]/4 + 0*rfld0[lfld[ 12]][0]/2 + 0*rfld0[lfld[ 13]][0]/2 + -3*rfld0[lfld[ 14]][0]/4 + -3*rfld0[lfld[ 15]][0]/4 + 0*rfld0[lfld[ 16]][0]/2 + 9*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[ 18]][0] = 1*rfld0[lfld[  0]][0]/3 + 1*rfld0[lfld[  1]][0]/3 + 0*rfld0[lfld[  2]][0]/2 + 1*rfld0[lfld[  3]][0]/3 + -3*rfld0[lfld[  4]][0]/4 + -3*rfld0[lfld[  5]][0]/4 + 0*rfld0[lfld[  6]][0]/2 + 0*rfld0[lfld[  7]][0]/2 + 0*rfld0[lfld[  8]][0]/2 + 0*rfld0[lfld[  9]][0]/2 + -3*rfld0[lfld[ 10]][0]/4 + -3*rfld0[lfld[ 11]][0]/4 + -3*rfld0[lfld[ 12]][0]/4 + -3*rfld0[lfld[ 13]][0]/4 + 0*rfld0[lfld[ 14]][0]/2 + 0*rfld0[lfld[ 15]][0]/2 + 0*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 9*rfld0[lfld[ 18]][0]/2 + 0*rfld0[lfld[ 19]][0]/2;
  rfld1[lfld[ 19]][0] = 1*rfld0[lfld[  0]][0]/3 + 1*rfld0[lfld[  1]][0]/3 + 1*rfld0[lfld[  2]][0]/3 + 0*rfld0[lfld[  3]][0]/2 + -3*rfld0[lfld[  4]][0]/4 + -3*rfld0[lfld[  5]][0]/4 + -3*rfld0[lfld[  6]][0]/4 + -3*rfld0[lfld[  7]][0]/4 + -3*rfld0[lfld[  8]][0]/4 + -3*rfld0[lfld[  9]][0]/4 + 0*rfld0[lfld[ 10]][0]/2 + 0*rfld0[lfld[ 11]][0]/2 + 0*rfld0[lfld[ 12]][0]/2 + 0*rfld0[lfld[ 13]][0]/2 + 0*rfld0[lfld[ 14]][0]/2 + 0*rfld0[lfld[ 15]][0]/2 + 0*rfld0[lfld[ 16]][0]/2 + 0*rfld0[lfld[ 17]][0]/2 + 0*rfld0[lfld[ 18]][0]/2 + 9*rfld0[lfld[ 19]][0]/2;
}

template<> void lag2bez3<3,2>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 2; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = rfld0[lfld[  3]][i];
    rfld1[lfld[  4]][i] = -5*rfld0[lfld[  0]][i]/6 + 1*rfld0[lfld[  1]][i]/3 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 6*rfld0[lfld[  4]][i]/2 + -3*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  5]][i] = 1*rfld0[lfld[  0]][i]/3 + -5*rfld0[lfld[  1]][i]/6 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + -3*rfld0[lfld[  4]][i]/2 + 6*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  6]][i] = 0*rfld0[lfld[  0]][i]/2 + -5*rfld0[lfld[  1]][i]/6 + 1*rfld0[lfld[  2]][i]/3 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 6*rfld0[lfld[  6]][i]/2 + -3*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  7]][i] = 0*rfld0[lfld[  0]][i]/2 + 1*rfld0[lfld[  1]][i]/3 + -5*rfld0[lfld[  2]][i]/6 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + -3*rfld0[lfld[  6]][i]/2 + 6*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  8]][i] = 1*rfld0[lfld[  0]][i]/3 + 0*rfld0[lfld[  1]][i]/2 + -5*rfld0[lfld[  2]][i]/6 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 6*rfld0[lfld[  8]][i]/2 + -3*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  9]][i] = -5*rfld0[lfld[  0]][i]/6 + 0*rfld0[lfld[  1]][i]/2 + 1*rfld0[lfld[  2]][i]/3 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + -3*rfld0[lfld[  8]][i]/2 + 6*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 10]][i] = -5*rfld0[lfld[  0]][i]/6 + 0*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 6*rfld0[lfld[ 10]][i]/2 + -3*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 11]][i] = 1*rfld0[lfld[  0]][i]/3 + 0*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + -5*rfld0[lfld[  3]][i]/6 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + -3*rfld0[lfld[ 10]][i]/2 + 6*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 12]][i] = 0*rfld0[lfld[  0]][i]/2 + -5*rfld0[lfld[  1]][i]/6 + 0*rfld0[lfld[  2]][i]/2 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 6*rfld0[lfld[ 12]][i]/2 + -3*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 13]][i] = 0*rfld0[lfld[  0]][i]/2 + 1*rfld0[lfld[  1]][i]/3 + 0*rfld0[lfld[  2]][i]/2 + -5*rfld0[lfld[  3]][i]/6 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + -3*rfld0[lfld[ 12]][i]/2 + 6*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 14]][i] = 0*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + -5*rfld0[lfld[  2]][i]/6 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 6*rfld0[lfld[ 14]][i]/2 + -3*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 15]][i] = 0*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + 1*rfld0[lfld[  2]][i]/3 + -5*rfld0[lfld[  3]][i]/6 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + -3*rfld0[lfld[ 14]][i]/2 + 6*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 16]][i] = 0*rfld0[lfld[  0]][i]/2 + 1*rfld0[lfld[  1]][i]/3 + 1*rfld0[lfld[  2]][i]/3 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + -3*rfld0[lfld[  6]][i]/4 + -3*rfld0[lfld[  7]][i]/4 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + -3*rfld0[lfld[ 12]][i]/4 + -3*rfld0[lfld[ 13]][i]/4 + -3*rfld0[lfld[ 14]][i]/4 + -3*rfld0[lfld[ 15]][i]/4 + 9*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 17]][i] = 1*rfld0[lfld[  0]][i]/3 + 0*rfld0[lfld[  1]][i]/2 + 1*rfld0[lfld[  2]][i]/3 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + -3*rfld0[lfld[  8]][i]/4 + -3*rfld0[lfld[  9]][i]/4 + -3*rfld0[lfld[ 10]][i]/4 + -3*rfld0[lfld[ 11]][i]/4 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + -3*rfld0[lfld[ 14]][i]/4 + -3*rfld0[lfld[ 15]][i]/4 + 0*rfld0[lfld[ 16]][i]/2 + 9*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 18]][i] = 1*rfld0[lfld[  0]][i]/3 + 1*rfld0[lfld[  1]][i]/3 + 0*rfld0[lfld[  2]][i]/2 + 1*rfld0[lfld[  3]][i]/3 + -3*rfld0[lfld[  4]][i]/4 + -3*rfld0[lfld[  5]][i]/4 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + -3*rfld0[lfld[ 10]][i]/4 + -3*rfld0[lfld[ 11]][i]/4 + -3*rfld0[lfld[ 12]][i]/4 + -3*rfld0[lfld[ 13]][i]/4 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 9*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 19]][i] = 1*rfld0[lfld[  0]][i]/3 + 1*rfld0[lfld[  1]][i]/3 + 1*rfld0[lfld[  2]][i]/3 + 0*rfld0[lfld[  3]][i]/2 + -3*rfld0[lfld[  4]][i]/4 + -3*rfld0[lfld[  5]][i]/4 + -3*rfld0[lfld[  6]][i]/4 + -3*rfld0[lfld[  7]][i]/4 + -3*rfld0[lfld[  8]][i]/4 + -3*rfld0[lfld[  9]][i]/4 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 9*rfld0[lfld[ 19]][i]/2;
  }
}

template<> void lag2bez3<3,3>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 3; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = rfld0[lfld[  3]][i];
    rfld1[lfld[  4]][i] = -5*rfld0[lfld[  0]][i]/6 + 1*rfld0[lfld[  1]][i]/3 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 6*rfld0[lfld[  4]][i]/2 + -3*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  5]][i] = 1*rfld0[lfld[  0]][i]/3 + -5*rfld0[lfld[  1]][i]/6 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + -3*rfld0[lfld[  4]][i]/2 + 6*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  6]][i] = 0*rfld0[lfld[  0]][i]/2 + -5*rfld0[lfld[  1]][i]/6 + 1*rfld0[lfld[  2]][i]/3 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 6*rfld0[lfld[  6]][i]/2 + -3*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  7]][i] = 0*rfld0[lfld[  0]][i]/2 + 1*rfld0[lfld[  1]][i]/3 + -5*rfld0[lfld[  2]][i]/6 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + -3*rfld0[lfld[  6]][i]/2 + 6*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  8]][i] = 1*rfld0[lfld[  0]][i]/3 + 0*rfld0[lfld[  1]][i]/2 + -5*rfld0[lfld[  2]][i]/6 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 6*rfld0[lfld[  8]][i]/2 + -3*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  9]][i] = -5*rfld0[lfld[  0]][i]/6 + 0*rfld0[lfld[  1]][i]/2 + 1*rfld0[lfld[  2]][i]/3 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + -3*rfld0[lfld[  8]][i]/2 + 6*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 10]][i] = -5*rfld0[lfld[  0]][i]/6 + 0*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 6*rfld0[lfld[ 10]][i]/2 + -3*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 11]][i] = 1*rfld0[lfld[  0]][i]/3 + 0*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + -5*rfld0[lfld[  3]][i]/6 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + -3*rfld0[lfld[ 10]][i]/2 + 6*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 12]][i] = 0*rfld0[lfld[  0]][i]/2 + -5*rfld0[lfld[  1]][i]/6 + 0*rfld0[lfld[  2]][i]/2 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 6*rfld0[lfld[ 12]][i]/2 + -3*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 13]][i] = 0*rfld0[lfld[  0]][i]/2 + 1*rfld0[lfld[  1]][i]/3 + 0*rfld0[lfld[  2]][i]/2 + -5*rfld0[lfld[  3]][i]/6 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + -3*rfld0[lfld[ 12]][i]/2 + 6*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 14]][i] = 0*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + -5*rfld0[lfld[  2]][i]/6 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 6*rfld0[lfld[ 14]][i]/2 + -3*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 15]][i] = 0*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + 1*rfld0[lfld[  2]][i]/3 + -5*rfld0[lfld[  3]][i]/6 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + -3*rfld0[lfld[ 14]][i]/2 + 6*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 16]][i] = 0*rfld0[lfld[  0]][i]/2 + 1*rfld0[lfld[  1]][i]/3 + 1*rfld0[lfld[  2]][i]/3 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + -3*rfld0[lfld[  6]][i]/4 + -3*rfld0[lfld[  7]][i]/4 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + -3*rfld0[lfld[ 12]][i]/4 + -3*rfld0[lfld[ 13]][i]/4 + -3*rfld0[lfld[ 14]][i]/4 + -3*rfld0[lfld[ 15]][i]/4 + 9*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 17]][i] = 1*rfld0[lfld[  0]][i]/3 + 0*rfld0[lfld[  1]][i]/2 + 1*rfld0[lfld[  2]][i]/3 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + -3*rfld0[lfld[  8]][i]/4 + -3*rfld0[lfld[  9]][i]/4 + -3*rfld0[lfld[ 10]][i]/4 + -3*rfld0[lfld[ 11]][i]/4 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + -3*rfld0[lfld[ 14]][i]/4 + -3*rfld0[lfld[ 15]][i]/4 + 0*rfld0[lfld[ 16]][i]/2 + 9*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 18]][i] = 1*rfld0[lfld[  0]][i]/3 + 1*rfld0[lfld[  1]][i]/3 + 0*rfld0[lfld[  2]][i]/2 + 1*rfld0[lfld[  3]][i]/3 + -3*rfld0[lfld[  4]][i]/4 + -3*rfld0[lfld[  5]][i]/4 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + -3*rfld0[lfld[ 10]][i]/4 + -3*rfld0[lfld[ 11]][i]/4 + -3*rfld0[lfld[ 12]][i]/4 + -3*rfld0[lfld[ 13]][i]/4 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 9*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 19]][i] = 1*rfld0[lfld[  0]][i]/3 + 1*rfld0[lfld[  1]][i]/3 + 1*rfld0[lfld[  2]][i]/3 + 0*rfld0[lfld[  3]][i]/2 + -3*rfld0[lfld[  4]][i]/4 + -3*rfld0[lfld[  5]][i]/4 + -3*rfld0[lfld[  6]][i]/4 + -3*rfld0[lfld[  7]][i]/4 + -3*rfld0[lfld[  8]][i]/4 + -3*rfld0[lfld[  9]][i]/4 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 9*rfld0[lfld[ 19]][i]/2;
  }
}

template<> void lag2bez3<3,6>(const int* __restrict__ lfld,
                                const dblAr2& __restrict__ rfld0,
                                dblAr2& __restrict__ rfld1){
  for(int i=0; i< 6; i++){
    rfld1[lfld[  0]][i] = rfld0[lfld[  0]][i];
    rfld1[lfld[  1]][i] = rfld0[lfld[  1]][i];
    rfld1[lfld[  2]][i] = rfld0[lfld[  2]][i];
    rfld1[lfld[  3]][i] = rfld0[lfld[  3]][i];
    rfld1[lfld[  4]][i] = -5*rfld0[lfld[  0]][i]/6 + 1*rfld0[lfld[  1]][i]/3 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + 6*rfld0[lfld[  4]][i]/2 + -3*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  5]][i] = 1*rfld0[lfld[  0]][i]/3 + -5*rfld0[lfld[  1]][i]/6 + 0*rfld0[lfld[  2]][i]/2 + 0*rfld0[lfld[  3]][i]/2 + -3*rfld0[lfld[  4]][i]/2 + 6*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  6]][i] = 0*rfld0[lfld[  0]][i]/2 + -5*rfld0[lfld[  1]][i]/6 + 1*rfld0[lfld[  2]][i]/3 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 6*rfld0[lfld[  6]][i]/2 + -3*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  7]][i] = 0*rfld0[lfld[  0]][i]/2 + 1*rfld0[lfld[  1]][i]/3 + -5*rfld0[lfld[  2]][i]/6 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + -3*rfld0[lfld[  6]][i]/2 + 6*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  8]][i] = 1*rfld0[lfld[  0]][i]/3 + 0*rfld0[lfld[  1]][i]/2 + -5*rfld0[lfld[  2]][i]/6 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 6*rfld0[lfld[  8]][i]/2 + -3*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[  9]][i] = -5*rfld0[lfld[  0]][i]/6 + 0*rfld0[lfld[  1]][i]/2 + 1*rfld0[lfld[  2]][i]/3 + 0*rfld0[lfld[  3]][i]/2 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + -3*rfld0[lfld[  8]][i]/2 + 6*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 10]][i] = -5*rfld0[lfld[  0]][i]/6 + 0*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 6*rfld0[lfld[ 10]][i]/2 + -3*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 11]][i] = 1*rfld0[lfld[  0]][i]/3 + 0*rfld0[lfld[  1]][i]/2 + 0*rfld0[lfld[  2]][i]/2 + -5*rfld0[lfld[  3]][i]/6 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + -3*rfld0[lfld[ 10]][i]/2 + 6*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 12]][i] = 0*rfld0[lfld[  0]][i]/2 + -5*rfld0[lfld[  1]][i]/6 + 0*rfld0[lfld[  2]][i]/2 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 6*rfld0[lfld[ 12]][i]/2 + -3*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 13]][i] = 0*rfld0[lfld[  0]][i]/2 + 1*rfld0[lfld[  1]][i]/3 + 0*rfld0[lfld[  2]][i]/2 + -5*rfld0[lfld[  3]][i]/6 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + -3*rfld0[lfld[ 12]][i]/2 + 6*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 14]][i] = 0*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + -5*rfld0[lfld[  2]][i]/6 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 6*rfld0[lfld[ 14]][i]/2 + -3*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 15]][i] = 0*rfld0[lfld[  0]][i]/2 + 0*rfld0[lfld[  1]][i]/2 + 1*rfld0[lfld[  2]][i]/3 + -5*rfld0[lfld[  3]][i]/6 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + -3*rfld0[lfld[ 14]][i]/2 + 6*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 16]][i] = 0*rfld0[lfld[  0]][i]/2 + 1*rfld0[lfld[  1]][i]/3 + 1*rfld0[lfld[  2]][i]/3 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + -3*rfld0[lfld[  6]][i]/4 + -3*rfld0[lfld[  7]][i]/4 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + -3*rfld0[lfld[ 12]][i]/4 + -3*rfld0[lfld[ 13]][i]/4 + -3*rfld0[lfld[ 14]][i]/4 + -3*rfld0[lfld[ 15]][i]/4 + 9*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 17]][i] = 1*rfld0[lfld[  0]][i]/3 + 0*rfld0[lfld[  1]][i]/2 + 1*rfld0[lfld[  2]][i]/3 + 1*rfld0[lfld[  3]][i]/3 + 0*rfld0[lfld[  4]][i]/2 + 0*rfld0[lfld[  5]][i]/2 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + -3*rfld0[lfld[  8]][i]/4 + -3*rfld0[lfld[  9]][i]/4 + -3*rfld0[lfld[ 10]][i]/4 + -3*rfld0[lfld[ 11]][i]/4 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + -3*rfld0[lfld[ 14]][i]/4 + -3*rfld0[lfld[ 15]][i]/4 + 0*rfld0[lfld[ 16]][i]/2 + 9*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 18]][i] = 1*rfld0[lfld[  0]][i]/3 + 1*rfld0[lfld[  1]][i]/3 + 0*rfld0[lfld[  2]][i]/2 + 1*rfld0[lfld[  3]][i]/3 + -3*rfld0[lfld[  4]][i]/4 + -3*rfld0[lfld[  5]][i]/4 + 0*rfld0[lfld[  6]][i]/2 + 0*rfld0[lfld[  7]][i]/2 + 0*rfld0[lfld[  8]][i]/2 + 0*rfld0[lfld[  9]][i]/2 + -3*rfld0[lfld[ 10]][i]/4 + -3*rfld0[lfld[ 11]][i]/4 + -3*rfld0[lfld[ 12]][i]/4 + -3*rfld0[lfld[ 13]][i]/4 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 9*rfld0[lfld[ 18]][i]/2 + 0*rfld0[lfld[ 19]][i]/2;
    rfld1[lfld[ 19]][i] = 1*rfld0[lfld[  0]][i]/3 + 1*rfld0[lfld[  1]][i]/3 + 1*rfld0[lfld[  2]][i]/3 + 0*rfld0[lfld[  3]][i]/2 + -3*rfld0[lfld[  4]][i]/4 + -3*rfld0[lfld[  5]][i]/4 + -3*rfld0[lfld[  6]][i]/4 + -3*rfld0[lfld[  7]][i]/4 + -3*rfld0[lfld[  8]][i]/4 + -3*rfld0[lfld[  9]][i]/4 + 0*rfld0[lfld[ 10]][i]/2 + 0*rfld0[lfld[ 11]][i]/2 + 0*rfld0[lfld[ 12]][i]/2 + 0*rfld0[lfld[ 13]][i]/2 + 0*rfld0[lfld[ 14]][i]/2 + 0*rfld0[lfld[ 15]][i]/2 + 0*rfld0[lfld[ 16]][i]/2 + 0*rfld0[lfld[ 17]][i]/2 + 0*rfld0[lfld[ 18]][i]/2 + 9*rfld0[lfld[ 19]][i]/2;
  }
}


} // End namespace

