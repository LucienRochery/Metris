// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define MATRIXS_MATMUL_NATIVE_INSTANTIATE
#include "MatrixS_MatMul_Native_impl.h"

#include <boost/preprocessor/repetition/repeat_from_to.hpp>

namespace SANS
{
namespace DLA
{

MATRIXS_MATMUL_NATIVE( 1,1, 1,2, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,1, 1,3, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,1, 1,4, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,1, 1,5, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,1, 1,6, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,1, 1,7, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 1,2, 2,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,2, 2,2, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,2, 2,3, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,2, 2,4, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,2, 2,5, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,2, 2,6, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,2, 2,7, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 1,3, 3,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,3, 3,2, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,3, 3,3, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,3, 3,4, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,3, 3,5, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,3, 3,6, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,3, 3,7, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 1,4, 4,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,4, 4,3, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,4, 4,4, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,4, 4,5, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,4, 4,6, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,4, 4,7, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 1,5, 5,4, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,5, 5,5, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,5, 5,6, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,5, 5,7, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 1,6, 6,5, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,6, 6,6, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,6, 6,7, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 1,7, 7,6, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 1,7, 7,7, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 1,8, 8,1, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 2,1, 1,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 2,1, 1,3, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 2,1, 1,4, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 2,1, 1,5, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 2,1, 1,6, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 2,1, 1,7, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 2,2, 2,3, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 2,3, 3,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 2,3, 3,2, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 2,3, 3,3, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 2,4, 4,1, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 2,5, 5,1, Real, Real, Real, Real ); // FP-IBL3 coupling
MATRIXS_MATMUL_NATIVE( 2,6, 6,1, Real, Real, Real, Real ); // FP-IBL3 coupling
MATRIXS_MATMUL_NATIVE( 2,8, 8,1, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 3,1, 1,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 3,1, 1,2, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 3,1, 1,3, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 3,1, 1,4, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 3,1, 1,5, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 3,1, 1,6, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 3,1, 1,7, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 3,2, 2,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 3,2, 2,2, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 3,2, 2,3, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 3,3, 3,2, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 3,8, 8,1, Real, Real, Real, Real );

// begin pcaplan
MATRIXS_MATMUL_NATIVE( 3,4, 4,4, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 3,3, 3,4, Real, Real, Real, Real );
// end pcaplan
MATRIXS_MATMUL_NATIVE( 3,4, 4,1, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 4,1, 1,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 4,1, 1,3, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 4,1, 1,4, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 4,1, 1,5, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 4,1, 1,6, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 4,1, 1,7, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 4,2, 2,1, Real, Real, Real, Real ); // for neural net with input dim 2 and width 4
MATRIXS_MATMUL_NATIVE( 4,3, 3,1, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 4,4, 4,1, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 4,8, 8,1, Real, Real, Real, Real ); // for 2D IBLAmpLagSplit

MATRIXS_MATMUL_NATIVE( 5,1, 1,4, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 5,1, 1,5, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 5,1, 1,6, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 5,1, 1,7, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 5,2, 2,1, Real, Real, Real, Real ); // for neural net with input dim 2 and width 5
MATRIXS_MATMUL_NATIVE( 5,3, 3,1, Real, Real, Real, Real ); // for neural net with input dim 2 and width 5

MATRIXS_MATMUL_NATIVE( 5,5, 5,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 5,6, 6,1, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 6,1, 1,5, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 6,1, 1,6, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 6,1, 1,7, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 6,2, 2,1, Real, Real, Real, Real ); // FP-IBL3 coupling
MATRIXS_MATMUL_NATIVE( 6,5, 5,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 6,6, 6,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 6,8, 8,1, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 7,1, 1,6, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 7,1, 1,7, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 7,7, 7,1, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 8,2, 2,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 8,3, 3,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 8,4, 4,1, Real, Real, Real, Real );
MATRIXS_MATMUL_NATIVE( 8,6, 6,1, Real, Real, Real, Real );

MATRIXS_MATMUL_NATIVE( 8,8, 8,1, Real, Real, Real, Real );

//Explicit instantiation of square matrix multiplications
#define DECL(z, n, text) \
MATRIXS_MATMUL_NATIVE( n,n, n,n, Real, Real, Real, Real );

BOOST_PP_REPEAT_FROM_TO(1, 16, DECL, )

} //namespace DLA
} //namespace SANS
