// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "../Eigen.h"

#include "../MatrixSymS.h"
#include "../VectorS.h"

#include "../MatrixS_Det.h"
#include "../MatrixS_Trace.h"

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
  //Trivial
  L[0] = A(0,0);
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

  EigenValues( A, L );

  //Trivial
  E = 1;
}

#define INSTANTIATE_EIGEN(T) \
template void EigenValues<1,T>(const MatrixSymS<1,T>& A, VectorS<1,T>& L ); \
template void EigenVectors<1,T>(const MatrixSymS<1,T>& A, MatrixS<1,1,T>& E ); \
template void EigenSystem<1,T>(const MatrixSymS<1,T>& A, VectorS<1,T>& L, MatrixS<1,1,T>& E );

INSTANTIATE_EIGEN(Real)
INSTANTIATE_EIGEN(SurrealS<1>)
INSTANTIATE_EIGEN(SurrealS<2>)

}
}
