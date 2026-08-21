// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define METRIS_MATRIXS_MATMUL_NATIVE_INSTANTIATE
#include "MatrixS_MatMul_Native_impl.h"

#include "../MatrixS.h"
#include "../VectorS.h"

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
// METRIS_MATRIXS_MATMUL_NATIVE( Am,An, Bm,Bn, A-type, B-type, s-type, C-type );
//
// Am, An are the sizes of the A matrix (similarly for B). The C matrix is deduced from A and B sizes.

typedef VectorS<1,Metris::Real> VectorS1;
typedef MatrixS<1,1,Metris::Real> MatrixS11;

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, VectorS1, Metris::Real, VectorS1 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, VectorS1, VectorS1, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, VectorS1, VectorS1, Metris::Real, VectorS1 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, VectorS1, VectorS1, Metris::Real, VectorS1 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, VectorS1, VectorS1, Metris::Real, VectorS1 );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS11, VectorS1, Metris::Real, VectorS1 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS11, MatrixS11, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS11, Metris::Real,      Metris::Real, MatrixS11 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS11, MatrixS11, Metris::Real, MatrixS11 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS11, MatrixS11, Metris::Real, MatrixS11 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS11, MatrixS11, Metris::Real, MatrixS11 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS11, MatrixS11, Metris::Real, MatrixS11 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, MatrixS11, MatrixS11, Metris::Real, MatrixS11 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS11, MatrixS11, Metris::Real, MatrixS11 );

typedef VectorS<2,Metris::Real> VectorS2;
typedef MatrixS<2,2,Metris::Real> MatrixS22;

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, VectorS2, Metris::Real, VectorS2 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, VectorS2, Metris::Real, VectorS2 );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, MatrixS11, Metris::Real, MatrixS11 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, MatrixS22, Metris::Real, MatrixS22 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, MatrixS22, Metris::Real, MatrixS22 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS22, MatrixS22, Metris::Real, MatrixS22 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS22, MatrixS22, Metris::Real, MatrixS22 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS22, MatrixS22, Metris::Real, MatrixS22 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,1, MatrixS22, MatrixS22, Metris::Real, MatrixS22 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS22, MatrixS22, Metris::Real, MatrixS22 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, MatrixS22, MatrixS22, Metris::Real, MatrixS22 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS22, MatrixS22, Metris::Real, MatrixS22 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS22, VectorS2, Metris::Real, VectorS2 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS22, VectorS2, Metris::Real, VectorS2 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS22, VectorS2, Metris::Real, VectorS2 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,1, MatrixS22, VectorS2, Metris::Real, VectorS2 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS22, VectorS2, Metris::Real, VectorS2 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,1, MatrixS22, VectorS2, Metris::Real, VectorS2 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS22, VectorS2, Metris::Real, VectorS2 );

typedef VectorS<3,Metris::Real> VectorS3;
typedef MatrixS<3,3,Metris::Real> MatrixS33;

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, Metris::Real, VectorS3, Metris::Real, VectorS3 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, VectorS3, Metris::Real, VectorS3 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, VectorS3, Metris::Real, VectorS3 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, Metris::Real, MatrixS33, Metris::Real, MatrixS33 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, MatrixS33, Metris::Real, MatrixS33 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, MatrixS33, Metris::Real, MatrixS33 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS33, Metris::Real, Metris::Real, MatrixS33 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS33, Metris::Real, Metris::Real, MatrixS33 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, VectorS3, MatrixS11, Metris::Real, VectorS3 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, VectorS3, Metris::Real, Metris::Real, VectorS3 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS33, MatrixS33, Metris::Real, MatrixS33 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS33, MatrixS33, Metris::Real, MatrixS33 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS33, MatrixS33, Metris::Real, MatrixS33 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS33, MatrixS33, Metris::Real, MatrixS33 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,1, MatrixS33, MatrixS33, Metris::Real, MatrixS33 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, MatrixS33, MatrixS33, Metris::Real, MatrixS33 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS33, MatrixS33, Metris::Real, MatrixS33 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS33, VectorS3, Metris::Real, VectorS3 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,1, MatrixS33, VectorS3, Metris::Real, VectorS3 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS33, VectorS3, Metris::Real, VectorS3 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS33, VectorS3, Metris::Real, VectorS3 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS33, VectorS3, Metris::Real, VectorS3 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,1, MatrixS33, VectorS3, Metris::Real, VectorS3 );

typedef VectorS<4,Metris::Real> VectorS4;
typedef MatrixS<4,4,Metris::Real> MatrixS44;

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, VectorS4, Metris::Real, VectorS4 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, VectorS4, Metris::Real, VectorS4 );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, MatrixS44, Metris::Real, MatrixS44 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, MatrixS44, Metris::Real, MatrixS44 );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44, Metris::Real, Metris::Real, MatrixS44 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS44, Metris::Real, Metris::Real, MatrixS44 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS44, MatrixS44, Metris::Real, MatrixS44 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS44, MatrixS44, Metris::Real, MatrixS44 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS44, MatrixS44, Metris::Real, MatrixS44 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44, MatrixS44, Metris::Real, MatrixS44 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, MatrixS44, MatrixS44, Metris::Real, MatrixS44 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS44, MatrixS44, Metris::Real, MatrixS44 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS44, VectorS4, Metris::Real, VectorS4 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,1, MatrixS44, VectorS4, Metris::Real, VectorS4 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44, VectorS4, Metris::Real, VectorS4 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS44, VectorS4, Metris::Real, VectorS4 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS44, VectorS4, Metris::Real, VectorS4 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS44, VectorS4, Metris::Real, VectorS4 );


typedef VectorS<5,Metris::Real> VectorS5;
typedef MatrixS<5,5,Metris::Real> MatrixS55;

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, VectorS5, Metris::Real, VectorS5 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, VectorS5, Metris::Real, VectorS5 );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, MatrixS55, Metris::Real, MatrixS55 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, MatrixS55, MatrixS55, Metris::Real, MatrixS55 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, MatrixS55, Metris::Real, MatrixS55 );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS55, Metris::Real, Metris::Real, MatrixS55 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS55, Metris::Real, Metris::Real, MatrixS55 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS55, VectorS5, Metris::Real, VectorS5 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,1, MatrixS55, VectorS5, Metris::Real, VectorS5 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS55, VectorS5, Metris::Real, VectorS5 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS55, MatrixS55, Metris::Real, MatrixS55 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS55, MatrixS55, Metris::Real, MatrixS55 );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS55, VectorS5, Metris::Real, VectorS5 );

typedef VectorS<6,Metris::Real> VectorS6;
typedef MatrixS<6,6,Metris::Real> MatrixS66;

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, VectorS6, Metris::Real, VectorS6 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, VectorS6, Metris::Real, VectorS6 );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, MatrixS66, Metris::Real, MatrixS66 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, MatrixS66, Metris::Real, MatrixS66 );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS66, Metris::Real, Metris::Real, MatrixS66 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS66, VectorS6, Metris::Real, VectorS6 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,1, MatrixS66, VectorS6, Metris::Real, VectorS6 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS66, VectorS6, Metris::Real, VectorS6 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS66, MatrixS66, Metris::Real, MatrixS66 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS66, MatrixS66, Metris::Real, MatrixS66 );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS66, VectorS6, Metris::Real, VectorS6 );

typedef VectorS<7,Metris::Real> VectorS7;
typedef MatrixS<7,7,Metris::Real> MatrixS77;

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, VectorS7, Metris::Real, VectorS7 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, MatrixS77, Metris::Real, MatrixS77 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS77, VectorS7, Metris::Real, VectorS7 );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, VectorS7, Metris::Real, VectorS7 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, MatrixS77, Metris::Real, MatrixS77 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS77, VectorS7, Metris::Real, VectorS7 );

typedef VectorS<8,Metris::Real> VectorS8;
typedef MatrixS<8,8,Metris::Real> MatrixS88;

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS88, VectorS8, Metris::Real, VectorS8 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS88, VectorS8, Metris::Real, VectorS8 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS88, MatrixS88, Metris::Real, MatrixS88 );

} //namespace DLA
} //namespace Metris
