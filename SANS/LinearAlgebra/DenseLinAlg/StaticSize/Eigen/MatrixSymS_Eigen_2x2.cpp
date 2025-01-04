// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "../Eigen.h"

#include "../MatrixSymS.h"
#include "../VectorS.h"

#include "Surreal/SurrealS.h"

#include "tools/SANSnumerics.h"

#include <cmath>  // sqrt

namespace SANS
{
namespace DLA
{

#define DELTA 1e-3

template< int M, class T >
void
EigenValues(const MatrixSymS<M,T>& A, VectorS<M,T>& L )
{
  const T& a = A(0,0);
  const T& b = A(1,0); const T& c = A(1,1);

  // Check if the off diagonal is small enough to be considered zero
  if ( DELTA*b + a == a &&
       DELTA*b + c == c )
  {
    // A is diagonal.
    L[0] = a;
    L[1] = c;
  }
  else
  {
    T sm = a + c;
    T df = a - c;
    T rt = sqrt(df*df + 4.0*b*b);
    T t;

    if (sm > 0.0)
    {
      L[0] = 0.5 * (sm + rt);
      t    = 1.0/L[0];
      L[1] = (a*t)*c - (b*t)*b;
    }
    else if (sm < 0.0)
    {
      L[1] = 0.5 * (sm - rt);
      t    = 1.0/L[1];
      L[0] = (a*t)*c - (b*t)*b;
    }
    else       // This case needs to be treated separately to avoid div by 0
    {
      L[0] =  0.5 * rt;
      L[1] = -0.5 * rt;
    }

#if 0
    T t = a + c; //trace
    T u = a - c;
    T v = u*u + 4.0*b*b;

    SANS_ASSERT( v >= 0 );

    T w = sqrt(v);

    L[0] = (t - w)/2.0;
    L[1] = (t + w)/2.0;
#endif
  }
}

template< int M, class T >
void
EigenVectors(const MatrixSymS<M,T>& A, MatrixS<M,M,T>& E )
{
  VectorS<M,T> L;
  EigenSystem( A, L, E );
}

template< int M, class T >
void
EigenSystem(const MatrixSymS<M,T>& A, VectorS<M,T>& L, MatrixS<M,M,T>& E )
{
  const T& a = A(0,0);
  const T& b = A(1,0); const T& c = A(1,1);

  // Check if the off diagonal is small enough to be considered zero
  if ( DELTA*b + a == a &&
       DELTA*b + c == c )
  {
    // A is diagonal.
    L[0] = a;
    L[1] = c;

    E(0,0) = 1;
    E(1,0) = 0;

    E(0,1) = 0;
    E(1,1) = 1;
  }
  else
  {
    // ----------------------------------------------------------------------------
    // Calculates the eigensystem of a real symmetric 2x2 matrix
    //    [ A  B ]
    //    [ B  C ]
    // in the form
    //    [ A  B ]  =  [ cs  -sn ] [ rt1   0  ] [  cs  sn ]
    //    [ B  C ]     [ sn   cs ] [  0   rt2 ] [ -sn  cs ]
    // where rt1 >= rt2. Note that this convention is different from the one used
    // in the LAPACK routine DLAEV2, where |rt1| >= |rt2|.
    // ----------------------------------------------------------------------------

    T sm = a + c;
    T df = a - c;
    T rt = sqrt(df*df + 4.0*b*b);
    T t, cs, sn;

    if (sm > 0.0)
    {
      L[0] = 0.5 * (sm + rt);
      t    = 1.0/L[0];
      L[1] = (a*t)*c - (b*t)*b;
    }
    else if (sm < 0.0)
    {
      L[1] = 0.5 * (sm - rt);
      t    = 1.0/L[1];
      L[0] = (a*t)*c - (b*t)*b;
    }
    else       // This case needs to be treated separately to avoid div by 0
    {
      L[0] =  0.5 * rt;
      L[1] = -0.5 * rt;
    }

    // Calculate eigenvectors
    if (df > 0.0)
      cs = df + rt;
    else
      cs = df - rt;

    // '>=' is required in the case of df == 0 and b << 1
    // using '>' results in a 1/b for any derivatives for a and c
    if (fabs(cs) >= 2.0*fabs(b))
    {
      t   = -2.0 * b / cs;
      sn = 1.0 / sqrt(1.0 + t*t);
      cs = t * sn;
    }
    else
    {
      t  = -0.5 * cs / b;
      cs = 1.0 / sqrt(1.0 + t*t);
      sn = t * cs;
    }

    if (df > 0.0)
    {
      t  = cs;
      cs = -sn;
      sn = t;
    }

    E(0,0) = cs; E(0,1) = -sn;
    E(1,0) = sn; E(1,1) =  cs;

#if 0
    T t = a + c; //trace
    T u = a - c;
    T v = u*u + 4.0*b*b;

    SANS_ASSERT( v >= 0 );

    T w = sqrt(v);

    //Eigenvalues
    L[0] = (t - w)/2.0;
    L[1] = (t + w)/2.0;

    //Eigenvectors
    T v0[2], v1[2];

    if (u < 0)
    {
      v0[0] = (u - w)/2.0;
      v0[1] = b;

      v1[0] = b;
      v1[1] = (-u + w)/2.0;
    }
    else
    {
      v0[0] = b;
      v0[1] = -(u + w)/2.0;

      v1[0] = -(u + w)/2.0;
      v1[1] = -b;
    }

    T norm0 = sqrt( v0[0]*v0[0] + v0[1]*v0[1] );
    T norm1 = sqrt( v1[0]*v1[0] + v1[1]*v1[1] );

    E(0,0) = v0[0]/norm0;
    E(1,0) = v0[1]/norm0;

    E(0,1) = v1[0]/norm1;
    E(1,1) = v1[1]/norm1;
#endif
  }
}

#define INSTANTIATE_EIGEN(T) \
template void EigenValues<2,T>(const MatrixSymS<2,T>& A, VectorS<2,T>& L ); \
template void EigenVectors<2,T>(const MatrixSymS<2,T>& A, MatrixS<2,2,T>& E ); \
template void EigenSystem<2,T>(const MatrixSymS<2,T>& A, VectorS<2,T>& L, MatrixS<2,2,T>& E );

INSTANTIATE_EIGEN(Real)
INSTANTIATE_EIGEN(SurrealS<1>)
INSTANTIATE_EIGEN(SurrealS<2>)
INSTANTIATE_EIGEN(SurrealS<3>)
INSTANTIATE_EIGEN(SurrealS<6>)
INSTANTIATE_EIGEN(SurrealS<9>)

}
}
