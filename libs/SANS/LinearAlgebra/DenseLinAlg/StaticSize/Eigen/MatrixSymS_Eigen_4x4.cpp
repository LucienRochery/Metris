// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "MatrixSymS_Eigen_impl.h"

#include "Surreal/SurrealS.h"

#include "tools/SANSnumerics.h"

namespace SANS
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

template< int M, class T >
void
EigenSystem(const MatrixSymS<M,T>& A, VectorS<M,T>& L, MatrixS<M,M,T>& E )
{
  // iterate with Jacobi to find the solution
  EigenSystem_Jacobi(A, L, E);
}

#define INSTANTIATE_EIGEN(T) \
template void EigenValues<4,T>(const MatrixSymS<4,T>& A, VectorS<4,T>& L ); \
template void EigenVectors<4,T>(const MatrixSymS<4,T>& A, MatrixS<4,4,T>& E ); \
template void EigenSystem<4,T>(const MatrixSymS<4,T>& A, VectorS<4,T>& L, MatrixS<4,4,T>& E );

INSTANTIATE_EIGEN(Real)
INSTANTIATE_EIGEN(SurrealS<1>)
INSTANTIATE_EIGEN(SurrealS<10>)
INSTANTIATE_EIGEN(SurrealS<20>)
INSTANTIATE_EIGEN(SurrealS<50>)

}
}
