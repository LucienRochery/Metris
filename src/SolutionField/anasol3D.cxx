//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "anasol.hxx"

#include "../../SANS/Surreal/SurrealS.h"
#include "../aux_exceptions.hxx"
#include "../linalg/symidx.hxx"

#include <cmath>


namespace Metris{


double anasol3D_1([[maybe_unused]] void *ctx, 
                  [[maybe_unused]] const double*__restrict__ crd, 
                  std::initializer_list<double*> dfun){

  int ndiff = dfun.size();

  if(ndiff > 0){

    auto it = dfun.begin();


    if(ndiff >= 1){
      double *d1fun = *it;
      METRIS_ASSERT(d1fun != NULL);
      for(int ii = 0; ii < 3; ii++) d1fun[ii] = 0;
    }

    if(ndiff >= 2){
      it++;
      double *d2fun = *it;
      METRIS_ASSERT(d2fun != NULL);
      for(int ii = 0; ii < 3; ii++) 
        for(int jj = ii; jj < 3; jj++) 
          d2fun[sym2idx(ii,jj)] = 0;
    }
    if(ndiff >= 3){
      it++;
      double *d3fun = *it;
      METRIS_ASSERT(d3fun != NULL);
      for(int ii = 0; ii < 3; ii++) 
        for(int jj = ii; jj < 3; jj++) 
          for(int kk = jj; kk < 3; kk++) 
            d3fun[sym3idx(ii,jj,kk)] = 0;
    }
    if(ndiff >= 4){
      METRIS_THROW_MSG(TODOExcept(), "ndiff >= 4 not implemented. Might be more pleasant to use codegen.")
    }

  }


  return 1.0;
}


} // End namespace
