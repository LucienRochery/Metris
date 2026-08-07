// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXS_MATMUL_NATIVE_INSTANTIATE
#include "MatrixS_MatMul_Native_impl.h"

namespace Metris
{
namespace DLA
{

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, int, int, Real, Real );

MATRIXS_MATMUL_NATIVE( 1,1, 1,1, int, int, Real, int );
MATRIXS_MATMUL_NATIVE( 1,2, 2,1, int, int, Real, int );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, int, int, Real, int );
MATRIXS_MATMUL_NATIVE( 2,2, 2,1, int, int, Real, int );
MATRIXS_MATMUL_NATIVE( 2,2, 2,2, int, int, Real, int );

MATRIXS_MATMUL_NATIVE( 2,3, 3,2, int, int, Real, int );
MATRIXS_MATMUL_NATIVE( 2,3, 3,3, int, int, Real, int );

MATRIXS_MATMUL_NATIVE( 3,2, 2,2, int, int, Real, int );
MATRIXS_MATMUL_NATIVE( 3,2, 2,3, int, int, Real, int );

} //namespace DLA
} //namespace Metris
