// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXS_MATMUL_NATIVE_INSTANTIATE
#include "MatrixS_MatMul_Native_impl.h"

#include "SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"
#include "SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"

namespace Metris
{
namespace DLA
{
//-----------------------------------------------------------------------------
// Instantiates the expressions
//
// C = s * A x B
//
// where s is a scalar, and A, B, C are matrices. The template arguments are:
//
// MATRIXS_MATMUL_NATIVE( Am,An, Bm,Bn, A-type, B-type, s-type, C-type );
//
// Am, An are the sizes of the A matrix (similarly for B). The C matrix is deduced from A and B sizes.

typedef VectorS<1,Real> VectorS1;
typedef MatrixS<1,1,Real> MatrixS11;

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, VectorS1, Real, VectorS1 );

MATRIXS_MATMUL_NATIVE( 1,2, 2,1, VectorS1, VectorS1, Real, Real );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, VectorS1, VectorS1, Real, VectorS1 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, VectorS1, VectorS1, Real, VectorS1 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, VectorS1, VectorS1, Real, VectorS1 );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS11, VectorS1, Real, VectorS1 );

MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS11, MatrixS11, Real, Real );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS11, Real,      Real, MatrixS11 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS11, MatrixS11, Real, MatrixS11 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS11, MatrixS11, Real, MatrixS11 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS11, MatrixS11, Real, MatrixS11 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS11, MatrixS11, Real, MatrixS11 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, MatrixS11, MatrixS11, Real, MatrixS11 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS11, MatrixS11, Real, MatrixS11 );

typedef VectorS<2,Real> VectorS2;
typedef MatrixS<2,2,Real> MatrixS22;

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, VectorS2, Real, VectorS2 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, VectorS2, Real, VectorS2 );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, MatrixS11, Real, MatrixS11 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, MatrixS22, Real, MatrixS22 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, MatrixS22, Real, MatrixS22 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS22, MatrixS22, Real, MatrixS22 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS22, MatrixS22, Real, MatrixS22 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS22, MatrixS22, Real, MatrixS22 );
MATRIXS_MATMUL_NATIVE( 1,3, 3,1, MatrixS22, MatrixS22, Real, MatrixS22 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS22, MatrixS22, Real, MatrixS22 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, MatrixS22, MatrixS22, Real, MatrixS22 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS22, MatrixS22, Real, MatrixS22 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS22, VectorS2, Real, VectorS2 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS22, VectorS2, Real, VectorS2 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS22, VectorS2, Real, VectorS2 );
MATRIXS_MATMUL_NATIVE( 2,1, 1,1, MatrixS22, VectorS2, Real, VectorS2 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS22, VectorS2, Real, VectorS2 );
MATRIXS_MATMUL_NATIVE( 3,1, 1,1, MatrixS22, VectorS2, Real, VectorS2 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS22, VectorS2, Real, VectorS2 );

typedef VectorS<3,Real> VectorS3;
typedef MatrixS<3,3,Real> MatrixS33;

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, Real, VectorS3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, VectorS3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, VectorS3, Real, VectorS3 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, Real, MatrixS33, Real, MatrixS33 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, MatrixS33, Real, MatrixS33 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, MatrixS33, Real, MatrixS33 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS33, Real, Real, MatrixS33 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS33, Real, Real, MatrixS33 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, VectorS3, MatrixS11, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 1,1, 1,1, VectorS3, Real, Real, VectorS3 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS33, MatrixS33, Real, MatrixS33 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS33, MatrixS33, Real, MatrixS33 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS33, MatrixS33, Real, MatrixS33 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS33, MatrixS33, Real, MatrixS33 );
MATRIXS_MATMUL_NATIVE( 1,3, 3,1, MatrixS33, MatrixS33, Real, MatrixS33 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, MatrixS33, MatrixS33, Real, MatrixS33 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS33, MatrixS33, Real, MatrixS33 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS33, VectorS3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 2,1, 1,1, MatrixS33, VectorS3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS33, VectorS3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS33, VectorS3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS33, VectorS3, Real, VectorS3 );
MATRIXS_MATMUL_NATIVE( 1,3, 3,1, MatrixS33, VectorS3, Real, VectorS3 );

typedef VectorS<4,Real> VectorS4;
typedef MatrixS<4,4,Real> MatrixS44;

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, VectorS4, Real, VectorS4 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, VectorS4, Real, VectorS4 );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, MatrixS44, Real, MatrixS44 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, MatrixS44, Real, MatrixS44 );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44, Real, Real, MatrixS44 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS44, Real, Real, MatrixS44 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS44, MatrixS44, Real, MatrixS44 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS44, MatrixS44, Real, MatrixS44 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS44, MatrixS44, Real, MatrixS44 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44, MatrixS44, Real, MatrixS44 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, MatrixS44, MatrixS44, Real, MatrixS44 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS44, MatrixS44, Real, MatrixS44 );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS44, VectorS4, Real, VectorS4 );
MATRIXS_MATMUL_NATIVE( 2,1, 1,1, MatrixS44, VectorS4, Real, VectorS4 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44, VectorS4, Real, VectorS4 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS44, VectorS4, Real, VectorS4 );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS44, VectorS4, Real, VectorS4 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS44, VectorS4, Real, VectorS4 );


typedef VectorS<5,Real> VectorS5;
typedef MatrixS<5,5,Real> MatrixS55;

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, VectorS5, Real, VectorS5 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, VectorS5, Real, VectorS5 );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, MatrixS55, Real, MatrixS55 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, MatrixS55, MatrixS55, Real, MatrixS55 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, MatrixS55, Real, MatrixS55 );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS55, Real, Real, MatrixS55 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS55, Real, Real, MatrixS55 );

MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS55, VectorS5, Real, VectorS5 );
MATRIXS_MATMUL_NATIVE( 2,1, 1,1, MatrixS55, VectorS5, Real, VectorS5 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS55, VectorS5, Real, VectorS5 );

MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS55, MatrixS55, Real, MatrixS55 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS55, MatrixS55, Real, MatrixS55 );

MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS55, VectorS5, Real, VectorS5 );

typedef VectorS<6,Real> VectorS6;
typedef MatrixS<6,6,Real> MatrixS66;

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, VectorS6, Real, VectorS6 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, VectorS6, Real, VectorS6 );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, MatrixS66, Real, MatrixS66 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, MatrixS66, Real, MatrixS66 );

MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS66, Real, Real, MatrixS66 );

MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS66, VectorS6, Real, VectorS6 );
MATRIXS_MATMUL_NATIVE( 2,1, 1,1, MatrixS66, VectorS6, Real, VectorS6 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS66, VectorS6, Real, VectorS6 );

MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS66, MatrixS66, Real, MatrixS66 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS66, MatrixS66, Real, MatrixS66 );

MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS66, VectorS6, Real, VectorS6 );

typedef VectorS<7,Real> VectorS7;
typedef MatrixS<7,7,Real> MatrixS77;

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, VectorS7, Real, VectorS7 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, MatrixS77, Real, MatrixS77 );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS77, VectorS7, Real, VectorS7 );

MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, VectorS7, Real, VectorS7 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, MatrixS77, Real, MatrixS77 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS77, VectorS7, Real, VectorS7 );

typedef VectorS<8,Real> VectorS8;
typedef MatrixS<8,8,Real> MatrixS88;

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS88, VectorS8, Real, VectorS8 );
MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS88, VectorS8, Real, VectorS8 );

MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS88, MatrixS88, Real, MatrixS88 );

} //namespace DLA
} //namespace Metris
