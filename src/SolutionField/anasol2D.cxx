//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "anasol.hxx"

#include "../../SANS/Surreal/SurrealS.h"
#include "../aux_exceptions.hxx"
#include "../linalg/symidx.hxx"

#include <cmath>
#include <cstdarg>


namespace Metris{


// Constant function = 1
double anasol2D_1([[maybe_unused]] void *ctx, 
                  [[maybe_unused]] const double*__restrict__ crd, 
                  std::initializer_list<double*> dfun){

  int ndiff = dfun.size();

  if(ndiff > 0){

    auto it = dfun.begin();


    if(ndiff >= 1){
      double *d1fun = *it;
      METRIS_ASSERT(d1fun != NULL);
      for(int ii = 0; ii < 2; ii++) d1fun[ii] = 0;
    }

    if(ndiff >= 2){
      it++;
      double *d2fun = *it;
      METRIS_ASSERT(d2fun != NULL);
      for(int ii = 0; ii < 2; ii++) 
        for(int jj = ii; jj < 2; jj++) 
          d2fun[sym2idx(ii,jj)] = 0;
    }
    if(ndiff >= 3){
      it++;
      double *d3fun = *it;
      METRIS_ASSERT(d3fun != NULL);
      for(int ii = 0; ii < 2; ii++) 
        for(int jj = ii; jj < 2; jj++) 
          for(int kk = jj; kk < 2; kk++) 
            d3fun[sym3idx(ii,jj,kk)] = 0;
    }
    if(ndiff >= 4){
      METRIS_THROW_MSG(TODOExcept(), "ndiff >= 4 not implemented. Might be more pleasant to use codegen.")
    }

  }


  return 1.0;
}


// Linear function = x + y + 1 if ctx == NULL, otherwise interpret ctx as double*
// and get C, d_x and d_y in this order; not implemented
double anasol2D_2([[maybe_unused]] void *ctx, 
                  [[maybe_unused]] const double*__restrict__ crd, 
                  std::initializer_list<double*> dfun){

  int ndiff = dfun.size();

  if(ndiff > 0){

    auto it = dfun.begin();

    if(ndiff >= 1){
      double *d1fun = *it;
      METRIS_ASSERT(d1fun != NULL);
      for(int ii = 0; ii < 2; ii++) d1fun[ii] = 1;
    }

    if(ndiff >= 2){
      it++;
      double *d2fun = *it;
      METRIS_ASSERT(d2fun != NULL);
      for(int ii = 0; ii < 2; ii++) 
        for(int jj = ii; jj < 2; jj++) 
          d2fun[sym2idx(ii,jj)] = 0;
    }
    if(ndiff >= 3){
      it++;
      double *d3fun = *it;
      METRIS_ASSERT(d3fun != NULL);
      for(int ii = 0; ii < 2; ii++) 
        for(int jj = ii; jj < 2; jj++) 
          for(int kk = jj; kk < 2; kk++) 
            d3fun[sym3idx(ii,jj,kk)] = 0;
    }
    if(ndiff >= 4){
      METRIS_THROW_MSG(TODOExcept(), "ndiff >= 4 not implemented. Might be more pleasant to use codegen.")
    }

  }


  return 1.0 + crd[0] + crd[1];
}


// Quadratic x^2 + 2xy + y^2 + x + y + 1
// TODO: read parameters from ctx. 
double anasol2D_3([[maybe_unused]] void *ctx, 
                  [[maybe_unused]] const double*__restrict__ crd, 
                  std::initializer_list<double*> dfun){
  int ndiff = dfun.size();

  if(ndiff > 0){

    auto it = dfun.begin();

    if(ndiff >= 1){
      double *d1fun = *it;
      METRIS_ASSERT(d1fun != NULL);
      d1fun[1] = d1fun[0] = 1 + 2*(crd[0]+crd[1]);
    }

    if(ndiff >= 2){
      it++;
      double *d2fun = *it;
      METRIS_ASSERT(d2fun != NULL);
      for(int ii = 0; ii < 2; ii++) 
        for(int jj = ii; jj < 2; jj++) 
          d2fun[sym2idx(ii,jj)] = 2;
    }
    if(ndiff >= 3){
      it++;
      double *d3fun = *it;
      METRIS_ASSERT(d3fun != NULL);
      for(int ii = 0; ii < 2; ii++) 
        for(int jj = ii; jj < 2; jj++) 
          for(int kk = jj; kk < 2; kk++) 
            d3fun[sym3idx(ii,jj,kk)] = 0;
    }
    if(ndiff >= 4){
      METRIS_THROW_MSG(TODOExcept(), "ndiff >= 4 not implemented. Might be more pleasant to use codegen.")
    }

  }


  return 1.0 + crd[0] + crd[1] + crd[0]*crd[0] + 2*crd[0]*crd[1] + crd[1]*crd[1];
}


// Cubic x^2 + y^2 + x^3 + y^3 
// TODO: read parameters from ctx. 
double anasol2D_4([[maybe_unused]] void *ctx, 
                  [[maybe_unused]] const double*__restrict__ crd, 
                  std::initializer_list<double*> dfun){
  int ndiff = dfun.size();

  if(ndiff > 0){

    auto it = dfun.begin();

    if(ndiff >= 1){
      double *d1fun = *it;
      METRIS_ASSERT(d1fun != NULL);
      d1fun[0] = 2*crd[0] + 3*crd[0]*crd[0];
      d1fun[1] = 2*crd[1] + 3*crd[1]*crd[1];
    }

    if(ndiff >= 2){
      it++;
      double *d2fun = *it;
      METRIS_ASSERT(d2fun != NULL);
      for(int ii = 0; ii < 2; ii++) 
        for(int jj = ii; jj < 2; jj++) 
          d2fun[sym2idx(ii,jj)] = 0;
      d2fun[sym2idx(0,0)] = 2 + 6*crd[0];
      d2fun[sym2idx(1,1)] = 2 + 6*crd[1];
    }
    if(ndiff >= 3){
      it++;
      double *d3fun = *it;
      METRIS_ASSERT(d3fun != NULL);
      for(int ii = 0; ii < 2; ii++) 
        for(int jj = ii; jj < 2; jj++) 
          for(int kk = jj; kk < 2; kk++) 
            d3fun[sym3idx(ii,jj,kk)] = 0;
      d3fun[sym3idx(0,0,0)] = 6;
      d3fun[sym3idx(1,1,1)] = 6;
    }
    if(ndiff >= 4){
      METRIS_THROW_MSG(TODOExcept(), "ndiff >= 4 not implemented. Might be more pleasant to use codegen.")
    }

  }

  return crd[0]*crd[0] + crd[1]*crd[1] +  crd[0]*crd[0]*crd[0] + crd[1]*crd[1]*crd[1];
}




// sin(2pi*x) sin(2pi*y) ; cf Arthur Bawin's thesis function f_1 pp28
double anasol2D_5([[maybe_unused]] void *ctx, 
                  [[maybe_unused]] const double*__restrict__ crd, 
                  std::initializer_list<double*> dfun){
  int ndiff = dfun.size();

  double pi = 3.14159265358979323846;


  if(ndiff > 0){

    auto it = dfun.begin();

    if(ndiff >= 1){
      double *d1fun = *it;
      METRIS_ASSERT(d1fun != NULL);
      d1fun[0] =  2*pi*cos(2*pi*crd[0])*cos(2*pi*crd[1]);
      d1fun[1] = -2*pi*sin(2*pi*crd[0])*sin(2*pi*crd[1]);
    }

    if(ndiff >= 2){
      it++;
      double *d2fun = *it;
      METRIS_ASSERT(d2fun != NULL);
      d2fun[sym2idx(0,0)] = -4*pi*pi*sin(2*pi*crd[0])*cos(2*pi*crd[1]);
      d2fun[sym2idx(0,1)] = -4*pi*pi*cos(2*pi*crd[0])*sin(2*pi*crd[1]);
      d2fun[sym2idx(1,1)] = -4*pi*pi*sin(2*pi*crd[0])*cos(2*pi*crd[1]);
    }
    if(ndiff >= 3){
      it++;
      double *d3fun = *it;
      METRIS_ASSERT(d3fun != NULL);
      d3fun[sym3idx(0,0,0)] = -8*pi*pi*pi*cos(2*pi*crd[0])*cos(2*pi*crd[1]);
      d3fun[sym3idx(0,0,1)] =  8*pi*pi*pi*sin(2*pi*crd[0])*sin(2*pi*crd[1]);
      d3fun[sym3idx(0,1,1)] = -8*pi*pi*pi*cos(2*pi*crd[0])*cos(2*pi*crd[1]);
      d3fun[sym3idx(1,1,1)] =  8*pi*pi*pi*sin(2*pi*crd[0])*sin(2*pi*crd[1]);
    }
    if(ndiff >= 4){
      METRIS_THROW_MSG(TODOExcept(), "ndiff >= 4 not implemented. Might be more pleasant to use codegen.")
    }

  }

  return sin(2*pi*crd[0])*cos(2*pi*crd[1]);
}

} // End namespace
