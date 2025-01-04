// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "../Eigen.h"

#include "../MatrixSymS.h"
#include "../VectorS.h"

#include <cmath> // sqrt
#include <sstream>
#include <iomanip> // std::setprecision

// See Kopp_2008_Efficient_numerical_diagonalization_of_hermitian_3x3_matrices.pdf
// in the PXLibrary for references on eigen value calculations
// See also: Numerical diagonalization of 3x3 matrices
// (https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/)

// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using the Jacobi algorithm. Based on dsyevj3 in 3x3-C

namespace SANS
{
namespace DLA
{

template< int M, class T >
void
EigenSystem_Jacobi(const MatrixSymS<M,T>& A, VectorS<M,T>& L, MatrixS<M,M,T>& E )
{
  MatrixSymS<M,T> Atmp = A; //create copy of A

  T so;              // Sums of off-diagonal elements
  T s, c, t;         // sin(phi), cos(phi), tan(phi) and temporary storage
  T g, h, z, theta;  // More temporary storage
  T thresh;

  // Initialize E to the identity matrix
  E = Identity();

  // Initialize L to diag(A)
  for (int i = 0; i < M; i++)
    L[i] = Atmp(i,i);

  // Main iteration loop
  for (int nIter = 0; nIter < 50; nIter++)
  {
    // Test for convergence
    so = 0.0;
    for (int p = 0; p < M; p++)
      for (int q = p+1; q < M; q++)
        so += fabs(Atmp(p,q));

    if (so == 0.0)
      return;

    if (nIter < 4)
      thresh = 0.2 * so / ((T) M*M);
    else
      thresh = 0.0;

    // Do sweep
    for (int p = 0; p < M; p++)
      for (int q = p+1; q < M; q++)
      {
        g = 100.0 * fabs(Atmp(p,q));
        if (nIter >= 4
            &&  fabs(L[p]) + g == fabs(L[p])
            &&  fabs(L[q]) + g == fabs(L[q]))
        {
          Atmp(p,q) = 0.0;
        }
        else if (fabs(Atmp(p,q)) > thresh)
        {
          // Calculate Jacobi transformation
          h = L[q] - L[p];
          if (fabs(h) + g == fabs(h))
          {
            t = Atmp(p,q) / h;
          }
          else
          {
            theta = 0.5 * h / Atmp(p,q);
            if (theta < 0.0)
              t = -1.0 / (sqrt(1.0 + theta*theta) - theta);
            else
              t = 1.0 / (sqrt(1.0 + theta*theta) + theta);
          }
          c = 1.0/sqrt(1.0 + t*t);
          s = t * c;
          z = t * Atmp(p,q);

          // Apply Jacobi transformation
          Atmp(p,q) = 0.0;
          L[p] -= z;
          L[q] += z;
          for (int r = 0; r < p; r++)
          {
            t = Atmp(r,p);
            Atmp(r,p) = c*t - s*Atmp(r,q);
            Atmp(r,q) = s*t + c*Atmp(r,q);
          }
          for (int r = p+1; r < q; r++)
          {
            t = Atmp(p,r);
            Atmp(p,r) = c*t - s*Atmp(r,q);
            Atmp(r,q) = s*t + c*Atmp(r,q);
          }
          for (int r = q+1; r < M; r++)
          {
            t = Atmp(p,r);
            Atmp(p,r) = c*t - s*Atmp(q,r);
            Atmp(q,r) = s*t + c*Atmp(q,r);
          }

          // Update eigenvectors
          for (int r = 0; r < M; r++)
          {
            t = E(r,p);
            E(r,p) = c*t - s*E(r,q);
            E(r,q) = s*t + c*E(r,q);
          }
        }
      }
  } //main loop

  std::stringstream ss;

  ss << "EigenSystem<M,T> - Jacobi algorithm did not converge!" << std::endl;
  ss << std::endl;
  ss << "Failing matrix: " << std::endl;
  ss << std::setprecision(16) << A << std::endl;

  SANS_DEVELOPER_EXCEPTION(ss.str());
}

}
}
