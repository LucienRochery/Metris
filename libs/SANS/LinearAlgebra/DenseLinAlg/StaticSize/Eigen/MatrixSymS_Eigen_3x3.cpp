// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "MatrixSymS_Eigen_impl.h"

#include "../MatrixS_Trace.h"
#include "../MatrixS_Det.h"

#include "../../../../Surreal/SurrealS.h"

#include "../../../../tools/SANSnumerics.h"

#include "../../tools/dot.h"

#include <cmath> // sqrt
#include <limits>

namespace Metris
{

namespace DLA
{

template< int M, class T >
void
EigenValues(const MatrixSymS<M,T>& A, VectorS<M,T>& L )
{
  MatrixS<M,M,T> E;
  EigenSystem( A, L, E );
}

template< int M, class T >
void
EigenVectors(const MatrixSymS<M,T>& A, MatrixS<M,M,T>& E )
{
  VectorS<M,T> L;
  EigenSystem( A, L, E );
}



// Constants
#define METRIS_EIGEN_SQRT3    1.73205080756887729352744634151   // sqrt(3)

// Macros
#define METRIS_EIGEN_SQR(x)      ((x)*(x))                        // x^2


// ----------------------------------------------------------------------------
template<class T>
void dsyevc3(const MatrixSymS<3,T>& A, VectorS<3,T>& L)
// ----------------------------------------------------------------------------
// Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
// analytical algorithm.
// Only the diagonal and upper triangular parts of A are accessed. The access
// is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   L: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
{
  T m, c1, c0;

  // Determine coefficients of characteristic poynomial. We write
  //       | a   d   f  |
  //  A =  | d*  b   e  |
  //       | f*  e*  c  |

  const T& a = A(0,0);
  const T& b = A(1,1);
  const T& c = A(2,2);
  const T& d = A(1,0);
  const T& e = A(2,1);
  const T& f = A(2,0);

  T de = d * e;                                   // d * e
  T dd = METRIS_EIGEN_SQR(d);                                  // d^2
  T ee = METRIS_EIGEN_SQR(e);                                  // e^2
  T ff = METRIS_EIGEN_SQR(f);                                  // f^2
  m  = a + b + c;
  c1 = (a*b + a*c + b*c) - (dd + ee + ff);        // a*b + a*c + b*c - d^2 - e^2 - f^2

  c0 = c*dd + a*ee + b*ff - a*b*c - 2.0 * f*de;    // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e

  T p, sqrt_p, q, cn, sn, phi;
  p = METRIS_EIGEN_SQR(m) - 3.0*c1;
  q = m*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
  sqrt_p = sqrt(fabs(p));

  phi = 27.0 * ( 0.25*METRIS_EIGEN_SQR(c1)*(p - c1) + c0*(q + 27.0/4.0*c0));
  if (unlikely(phi == 0 && q == 0)) // avoid divide by zero
    phi = 0;
  else
    phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);

  cn = sqrt_p*cos(phi);
  sn = (1.0/METRIS_EIGEN_SQRT3)*sqrt_p*sin(phi);

  L[1]  = (1.0/3.0)*(m - cn);
  L[2]  = L[1] + sn;
  L[0]  = L[1] + cn;
  L[1] -= sn;
}

template< int M, class T >
void
EigenSystem(const MatrixSymS<M,T>& A, VectorS<M,T>& L, MatrixS<M,M,T>& E )
{
  T norm;          // Squared norm or inverse norm of current eigenvector
  T error;         // Estimated maximum roundoff error
  T t, u;          // Intermediate storage

  // Calculate eigenvalues
  dsyevc3(A, L);

  t = fabs(L[0]);
  if ((u=fabs(L[1])) > t)
    t = u;
  if ((u=fabs(L[2])) > t)
    t = u;
  if (t < 1.0)
    u = t;
  else
    u = METRIS_EIGEN_SQR(t);
  error = 256.0 * std::numeric_limits<Real>::epsilon() * METRIS_EIGEN_SQR(u);

  //       | a   d   f  |
  //  A =  | d*  b   e  |
  //       | f*  e*  c  |

  const T& a = A(0,0);
  const T& b = A(1,1);
  const T& d = A(1,0);
  const T& e = A(2,1);
  const T& f = A(2,0);

  E(0,1) = d*e - f*b;
  E(1,1) = f*d - e*a;
  E(2,1) = METRIS_EIGEN_SQR(d);

  // Calculate first eigenvector by the formula
  //   v[0] = (A - L[0]).e1 x (A - L[0]).e2
  E(0,0) = E(0,1) + f*L[0];
  E(1,0) = E(1,1) + e*L[0];
  E(2,0) = (a - L[0]) * (b - L[0]) - E(2,1);
  norm   = METRIS_EIGEN_SQR(E(0,0)) + METRIS_EIGEN_SQR(E(1,0)) + METRIS_EIGEN_SQR(E(2,0));

  // If vectors are nearly linearly dependent, or if there might have
  // been large cancellations in the calculation of A[i][i] - w[0], fall
  // back to Jacobi algorithm
  // Note that this simultaneously ensures that multiple eigenvalues do
  // not cause problems: If L[0] = L[1], then A - L[0] * I has rank 1,
  // i.e. all columns of A - L[0] * I are linearly dependent.
  if (norm <= error)
  {
    EigenSystem_Jacobi(A, L, E);
    return;
  }
  else                      // This is the standard branch
  {
    norm = sqrt(1.0 / norm);
    for (int j = 0; j < M; j++)
      E(j,0) = E(j,0) * norm;
  }

  // Calculate second eigenvector by the formula
  //   v[1] = (A - L[1]).e1 x (A - L[1]).e2
  E(0,1)  = E(0,1) + f*L[1];
  E(1,1)  = E(1,1) + e*L[1];
  E(2,1)  = (a - L[1]) * (b - L[1]) - E(2,1);
  norm    = METRIS_EIGEN_SQR(E(0,1)) + METRIS_EIGEN_SQR(E(1,1)) + METRIS_EIGEN_SQR(E(2,1));
  if (norm <= error)
  {
    EigenSystem_Jacobi(A, L, E);
    return;
  }
  else
  {
    norm = sqrt(1.0 / norm);
    for (int j = 0; j < M; j++)
      E(j,1) = E(j,1) * norm;
  }

  // Calculate third eigenvector according to
  //   v[2] = v[0] x v[1]
  E(0,2) = E(1,0)*E(2,1) - E(2,0)*E(1,1);
  E(1,2) = E(2,0)*E(0,1) - E(0,0)*E(2,1);
  E(2,2) = E(0,0)*E(1,1) - E(1,0)*E(0,1);
}

#define METRIS_INSTANTIATE_EIGEN(T) \
template void EigenValues<3,T>(const MatrixSymS<3,T>& A, VectorS<3,T>& L ); \
template void EigenVectors<3,T>(const MatrixSymS<3,T>& A, MatrixS<3,3,T>& E ); \
template void EigenSystem<3,T>(const MatrixSymS<3,T>& A, VectorS<3,T>& L, MatrixS<3,3,T>& E );

METRIS_INSTANTIATE_EIGEN(Real)
METRIS_INSTANTIATE_EIGEN(SurrealS<1>)
METRIS_INSTANTIATE_EIGEN(SurrealS<6>)
METRIS_INSTANTIATE_EIGEN(SurrealS<12>)
METRIS_INSTANTIATE_EIGEN(SurrealS<24>)

}
}
