//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_GEN_BARY__
#define __METRIS_GEN_BARY__

#include "types_arrays.hxx"

#include <random>

namespace Metris{
static void genBary(int nsamp, int tdim, dblAr2 &bary){
  bary.allocate(nsamp, tdim+1);
  bary.set_n(nsamp);
  bary.set_stride(tdim+1);
  std::uniform_real_distribution<double> unif(0.0,1.0);
  std::default_random_engine rng(0);
  for(int isamp = 0; isamp < nsamp; isamp++){
    double sum = 0;
    do{
      for(int jj = 0; jj < tdim+1; jj++){
        bary(isamp,jj) = unif(rng);
        sum += bary(isamp,jj);
      }
    }while(abs(sum) < 1.0e-16);

    for(int jj = 0; jj < tdim+1; jj++){
      bary(isamp,jj) /= sum;
    }
  }
}


}
#endif