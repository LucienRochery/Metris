//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
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
#include "../utils/aux_misc.hxx"

// Macros
#define SQR(x)      ((x)*(x))                        // x^2 


namespace Metris{


// ----------------------------------------------------------------------------
template<int ndim, typename T>
int dsyevq(const T* __restrict__ mat, T* __restrict__ eigvec, T* __restrict__  eigval)
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
  dsytrd<T,ndim>(mat, eigvec, eigval, e);
  
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
        g = abs(eigval[m])+abs(eigval[m+1]);
        if (abs(e[m]) + g == g)
          break;
      }
      if (m == l)
        break;
      
      if (nIter++ >= 30)
        return -1;

      // Calculate g = d_m - k
      g = (eigval[l+1] - eigval[l]) / (e[l] + e[l]);
      r = sqrt(SQR(g) + 1.0);
      if (g > 0)
        g = eigval[m] - eigval[l] + e[l]/(g + r);
      else
        g = eigval[m] - eigval[l] + e[l]/(g - r);

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
        
        g = eigval[i+1] - p;
        r = (eigval[i] - g)*s + 2.0*c*b;
        p = s * r;
        eigval[i+1] = g + p;
        g = c*r - b;

        // Form eigenvectors
        for (int k=0; k < ndim; k++)
        {
          t = eigvec[ndim*(i+1) + k];
          eigvec[ndim*(i+1) + k] = s*eigvec[ndim*i + k] + c*t;
          eigvec[ndim*i     + k] = c*eigvec[ndim*i + k] - s*t;
        }
      }
      eigval[l] -= p;
      e[l]  = g;
      e[m]  = 0.0;
    }
  }


  //int idx[3];

  #if 0
  int idx[ndim];
  if constexpr(ndim == 3){
    if(eigval[0] < eigval[1] && eigval[0] < eigval[2]){// 0 ..
      if(eigval[1] < eigval[2]){ // 0 1 2 already sorted
        return 0;
      }else{ // 0 2 1 
        swi(eigval[1],eigval[2]);
        for(int i = 0; i < ndim; i++) swi(eigvec[ndim*1+i], eigvec[ndim*2+i]);
      }
    }else{ // ...
      if(eigval[1] < eigval[0] && eigval[1] < eigval[2]){  // 1 ..
        if(eigval[0] < eigval[2]){ // 1 0 2
          swi(eigval[1],eigval[0]);
          for(int i = 0; i < ndim; i++) swi(eigvec[ndim*1+i], eigvec[ndim*0+i]);
        }else{ // 1 2 0 
          swi(eigval[0],eigval[2]);
          for(int i = 0; i < ndim; i++) swi(eigvec[ndim*2+i], eigvec[ndim*0+i]);

          swi(eigval[0],eigval[1]);
          for(int i = 0; i < ndim; i++) swi(eigvec[ndim*0+i], eigvec[ndim*1+i]);
        }
      }else{// 2 ..
        if(eigval[0] < eigval[1]){ // 2 0 1
          swi(eigval[2],eigval[0]);
          for(int i = 0; i < ndim; i++) swi(eigvec[ndim*2+i], eigvec[ndim*0+i]);

          swi(eigval[2],eigval[1]);
          for(int i = 0; i < ndim; i++) swi(eigvec[ndim*2+i], eigvec[ndim*1+i]);
        }else{ // 2 1 0
          swi(eigval[2],eigval[0]);
          for(int i = 0; i < ndim; i++) swi(eigvec[ndim*2+i], eigvec[ndim*0+i]);
        }
      }
    }
  }else{
    if(eigval[1] < eigval[0]){
      swi(eigval[1],eigval[0]);
      for(int i = 0; i < ndim; i++) swi(eigvec[ndim*1+i], eigvec[ndim*0+i]);
    }
  }
  #endif

  return 0;
}


template int dsyevq<2,double>(const double* __restrict__ A, double* __restrict__  Q, double* __restrict__  w);
template int dsyevq<3,double>(const double* __restrict__ A, double* __restrict__  Q, double* __restrict__  w);
template int dsyevq<2,SANS::SurrealS<2,double>>(const SANS::SurrealS<2,double>* __restrict__ A, SANS::SurrealS<2,double>* __restrict__  Q, SANS::SurrealS<2,double>* __restrict__  w);
template int dsyevq<3,SANS::SurrealS<2,double>>(const SANS::SurrealS<2,double>* __restrict__ A, SANS::SurrealS<2,double>* __restrict__  Q, SANS::SurrealS<2,double>* __restrict__  w);
template int dsyevq<3,SANS::SurrealS<3,double>>(const SANS::SurrealS<3,double>* __restrict__ A, SANS::SurrealS<3,double>* __restrict__  Q, SANS::SurrealS<3,double>* __restrict__  w);

#undef SQR

} // End namespace

