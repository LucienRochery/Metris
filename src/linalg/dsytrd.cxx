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
#include "../linalg/sym3idx.hxx"
#include "../SANS/Surreal/SurrealS.h"


// Macros
#define SQR(x) ((x)*(x))                        // x^2 

namespace Metris{

// ----------------------------------------------------------------------------
template<typename T>
void dsytrd3(const T* __restrict__ A, T* __restrict__  Q, T* __restrict__  d, T* __restrict__  e)
// ----------------------------------------------------------------------------
// Reduces a symmetric 3x3 matrix to tridiagonal form by applying
// (unitary) Householder transformations:
//            [ d[0]  e[0]       ]
//    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
//            [       e[1]  d[2] ]
// The function accesses only the diagonal and upper triangular parts of
// A. The access is read-only.
// ---------------------------------------------------------------------------
{
  const int n = 3;
  T u[n], q[n];
  T omega, f;
  T K, h, g;
  //int idx[3];
  
  // Initialize Q to the identitity matrix
#ifndef EVALS_ONLY
  for (int i=0; i < n; i++)
  {
    Q[3*i+i] = 1.0;
    for (int j=0; j < i; j++)
      Q[3*i+j] = Q[3*j+i] = 0.0;
  }
#endif

  // Bring first row and column to the desired form 
  h = SQR(A[sym3idx(0,1)]) + SQR(A[sym3idx(0,2)]);
  if (A[sym3idx(0,1)] > 0)
    g = -sqrt(h);
  else
    g = sqrt(h);
  e[0] = g;
  f    = g * A[sym3idx(0,1)];
  u[1] = A[sym3idx(0,1)] - g;
  u[2] = A[sym3idx(0,2)];
  
  omega = h - f;
  if (omega > 0.0)
  {
    omega = 1.0 / omega;
    K     = 0.0;
    for (int i=1; i < n; i++)
    {
      f    = A[sym3idx(1,i)] * u[1] + A[sym3idx(i,2)] * u[2];
      q[i] = omega * f;                  // p
      K   += u[i] * f;                   // u* A u
    }
    K *= 0.5 * SQR(omega);

    for (int i=1; i < n; i++)
      q[i] = q[i] - K * u[i];
    
    d[0] = A[sym3idx(0,0)];
    d[1] = A[sym3idx(1,1)] - 2.0*q[1]*u[1];
    d[2] = A[sym3idx(2,2)] - 2.0*q[2]*u[2];

    //if(d[0] < d[1] && d[0] < d[2]){
    //  idx[0] = 0;
    //    if(d[1] < d[2]){
    //      idx[1] = 1;
    //      idx[2] = 2;
    //    }else{
    //      idx[1] = 2;
    //      idx[2] = 1;
    //    }
    //}else{
    //  if(d[1] < d[0] && d[1] < d[2]){
    //    idx[0] = 1;
    //    if(d[0] < d[2]){
    //      idx[1] = 0;
    //      idx[2] = 2;
    //    }else{
    //      idx[1] = 2;
    //      idx[2] = 0;
    //    }
    //  }else{
    //    idx[0] = 2;
    //    if(d[0] < d[1]){
    //      idx[1] = 0;
    //      idx[2] = 1;
    //    }else{
    //      idx[1] = 1;
    //      idx[2] = 0;
    //    }
    //  }
    //}

    //sortupto8_inc(d,3);
    
    // Store inverse Householder transformation in Q
    for (int j=1; j < n; j++)
    {
      f = omega * u[j];
      for (int i=1; i < n; i++)
        //Q[3*idx[i] + j] = Q[3*idx[i] + j] - f*u[idx[i]];
        Q[3*i + j] = Q[3*i + j] - f*u[i];
    }

    // Calculate updated A[sym3idx[1][+ ]2] and store it in e[1]
    e[1] = A[sym3idx(1,2)] - q[1]*u[2] - u[1]*q[2];
  }
  else
  {
    for (int i=0; i < n; i++)
      d[i] = A[sym3idx(i,i)];
    e[1] = A[sym3idx(1,2)];
  }
}


// ----------------------------------------------------------------------------
template<typename T>
void dsytrd2(const T* __restrict__ A, T* __restrict__  Q, T* __restrict__  d, T* __restrict__  e)
// ----------------------------------------------------------------------------
// Reduces a symmetric 3x3 matrix to tridiagonal form by applying
// (unitary) Householder transformations:
//            [ d[0]  e[0]       ]
//    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
//            [       e[1]  d[2] ]
// The function accesses only the diagonal and upper triangular parts of
// A. The access is read-only.
// ---------------------------------------------------------------------------
{
  const int n = 2;
  T u[n], q[n];
  T omega, f;
  T K, h, g;
  
  // Initialize Q to the identitity matrix
#ifndef EVALS_ONLY
  for (int i=0; i < n; i++)
  {
    Q[2*i+i] = 1.0;
    for (int j=0; j < i; j++)
      Q[2*i+j] = Q[2*j+i] = 0.0;
  }
#endif

  // Bring first row and column to the desired form 
  h = SQR(A[sym3idx(0,1)]);
  if (A[sym3idx(0,1)] > 0)
    g = -sqrt(h);
  else
    g = sqrt(h);
  e[0] = g;
  f    = g * A[sym3idx(0,1)];
  u[1] = A[sym3idx(0,1)] - g;
  
  omega = h - f;
  if (omega > 0.0)
  {
    omega = 1.0 / omega;
    K     = 0.0;
    f    = A[sym3idx(1,1)] * u[1];
    q[1] = omega * f;                  // p
    K   += u[1] * f;                   // u* A u
    K *= 0.5 * SQR(omega);

    q[1] = q[1] - K * u[1];
    
    d[0] = A[sym3idx(0,0)];
    d[1] = A[sym3idx(1,1)] - 2.0*q[1]*u[1];
    
    // Store inverse Householder transformation in Q
    f = omega * u[1];
    Q[2*1 + 1] = Q[2*1 + 1] - f*u[1];

    // Calculate updated A[sym3idx[1][+ ]2] and store it in e[1]
    e[1] = 0; 
  }
  else
  {
    for (int i=0; i < n; i++)
      d[i] = A[sym3idx(i,i)];
    e[1] = 0;
  }
}

template void dsytrd3<double>(const double* __restrict__ A, double* __restrict__  Q, double* __restrict__  d, double* __restrict__  e);
template void dsytrd3<SANS::SurrealS<3,double>>(const SANS::SurrealS<3,double>* __restrict__ A, SANS::SurrealS<3,double>* __restrict__  Q, 
                                                      SANS::SurrealS<3,double>* __restrict__  d, SANS::SurrealS<3,double>* __restrict__  e);
template void dsytrd3<SANS::SurrealS<2,double>>(const SANS::SurrealS<2,double>* __restrict__ A, SANS::SurrealS<2,double>* __restrict__  Q, 
                                                      SANS::SurrealS<2,double>* __restrict__  d, SANS::SurrealS<2,double>* __restrict__  e);

template void dsytrd2<double>(const double* __restrict__ A, double* __restrict__  Q, double* __restrict__  d, double* __restrict__  e);
template void dsytrd2<SANS::SurrealS<2,double>>(const SANS::SurrealS<2,double>* __restrict__ A, SANS::SurrealS<2,double>* __restrict__  Q, 
                                                  SANS::SurrealS<2,double>* __restrict__  d, SANS::SurrealS<2,double>* __restrict__  e);

} // End namespace
#undef SQR
