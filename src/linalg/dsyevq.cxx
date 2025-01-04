//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



// ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
// ----------------------------------------------------------------------------
#include <cmath>
#include "../linalg/dsytrd.hxx"
#include "../linalg/dsyevq.hxx"
#include "../SANS/Surreal/SurrealS.h"
#include "../aux_utils.hxx"

// Macros
#define SQR(x)      ((x)*(x))                        // x^2 


namespace Metris{


// ----------------------------------------------------------------------------
template<int ndim, typename T>
int dsyevq(const T* __restrict__ A, T* __restrict__  Q, T* __restrict__  w)
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using the QL algorithm with implicit shifts, preceded by a
// Householder reduction to tridiagonal form.
// The function accesses only the diagonal and upper triangular parts of A.
// The access is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
// Dependencies:
//   dsytrd3()
// ----------------------------------------------------------------------------
{
  T e[ndim];                   // The third element is used only as temporary workspace
  T g, r, p, f, b, s, c, t; // Intermediate storage
  int nIter;
  int m;

  // Transform A to real tridiagonal form by the Householder method
  dsytrd<T,ndim>(A, Q, w, e);
  
  // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
  // with the QL method
  //
  // Loop over all off-diagonal elements
  for (int l=0; l < ndim-1; l++)
  {
    nIter = 0;
    while (1)
    {
      // Check for convergence and exit iteration loop if off-diagonal
      // element e(l) is zero
      for (m=l; m <= ndim-2; m++)
      {
        g = abs(w[m])+abs(w[m+1]);
        if (abs(e[m]) + g == g)
          break;
      }
      if (m == l)
        break;
      
      if (nIter++ >= 30)
        return -1;

      // Calculate g = d_m - k
      g = (w[l+1] - w[l]) / (e[l] + e[l]);
      r = sqrt(SQR(g) + 1.0);
      if (g > 0)
        g = w[m] - w[l] + e[l]/(g + r);
      else
        g = w[m] - w[l] + e[l]/(g - r);

      s = c = 1.0;
      p = 0.0;
      for (int i=m-1; i >= l; i--)
      {
        f = s * e[i];
        b = c * e[i];
        if (abs(f) > abs(g))
        {
          c      = g / f;
          r      = sqrt(SQR(c) + 1.0);
          e[i+1] = f * r;
          c     *= (s = 1.0/r);
        }
        else
        {
          s      = f / g;
          r      = sqrt(SQR(s) + 1.0);
          e[i+1] = g * r;
          s     *= (c = 1.0/r);
        }
        
        g = w[i+1] - p;
        r = (w[i] - g)*s + 2.0*c*b;
        p = s * r;
        w[i+1] = g + p;
        g = c*r - b;

        // Form eigenvectors
        for (int k=0; k < ndim; k++)
        {
          t = Q[ndim*(i+1) + k];
          Q[ndim*(i+1) + k] = s*Q[ndim*i + k] + c*t;
          Q[ndim*i     + k] = c*Q[ndim*i + k] - s*t;
        }
      }
      w[l] -= p;
      e[l]  = g;
      e[m]  = 0.0;
    }
  }


  //int idx[3];

  if constexpr(ndim == 3){
    if(w[0] < w[1] && w[0] < w[2]){// 0 ..
      if(w[1] < w[2]){ // 0 1 2 already sorted
        return 0;
      }else{ // 0 2 1 
        swi(w[1],w[2]);
        for(int i = 0; i < ndim; i++) swi(Q[ndim*1+i], Q[ndim*2+i]);
      }
    }else{ // ...
      if(w[1] < w[0] && w[1] < w[2]){  // 1 ..
        if(w[0] < w[2]){ // 1 0 2
          swi(w[1],w[0]);
          for(int i = 0; i < ndim; i++) swi(Q[ndim*1+i], Q[ndim*0+i]);
        }else{ // 1 2 0 
          swi(w[0],w[2]);
          for(int i = 0; i < ndim; i++) swi(Q[ndim*2+i], Q[ndim*0+i]);

          swi(w[0],w[1]);
          for(int i = 0; i < ndim; i++) swi(Q[ndim*0+i], Q[ndim*1+i]);
        }
      }else{// 2 ..
        if(w[0] < w[1]){ // 2 0 1
          swi(w[2],w[0]);
          for(int i = 0; i < ndim; i++) swi(Q[ndim*2+i], Q[ndim*0+i]);

          swi(w[2],w[1]);
          for(int i = 0; i < ndim; i++) swi(Q[ndim*2+i], Q[ndim*1+i]);
        }else{ // 2 1 0
          swi(w[2],w[0]);
          for(int i = 0; i < ndim; i++) swi(Q[ndim*2+i], Q[ndim*0+i]);
        }
      }
    }
  }else{
    if(w[1] < w[0]){
      swi(w[1],w[0]);
      for(int i = 0; i < ndim; i++) swi(Q[ndim*1+i], Q[ndim*0+i]);
    }
  }

  return 0;
}


template int dsyevq<2,double>(const double* __restrict__ A, double* __restrict__  Q, double* __restrict__  w);
template int dsyevq<3,double>(const double* __restrict__ A, double* __restrict__  Q, double* __restrict__  w);
template int dsyevq<2,SANS::SurrealS<2,double>>(const SANS::SurrealS<2,double>* __restrict__ A, SANS::SurrealS<2,double>* __restrict__  Q, SANS::SurrealS<2,double>* __restrict__  w);
template int dsyevq<3,SANS::SurrealS<2,double>>(const SANS::SurrealS<2,double>* __restrict__ A, SANS::SurrealS<2,double>* __restrict__  Q, SANS::SurrealS<2,double>* __restrict__  w);
template int dsyevq<3,SANS::SurrealS<3,double>>(const SANS::SurrealS<3,double>* __restrict__ A, SANS::SurrealS<3,double>* __restrict__  Q, SANS::SurrealS<3,double>* __restrict__  w);

#undef SQR

} // End namespace

