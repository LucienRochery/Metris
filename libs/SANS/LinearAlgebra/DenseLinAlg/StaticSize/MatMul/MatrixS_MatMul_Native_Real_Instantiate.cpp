// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define METRIS_MATRIXS_MATMUL_NATIVE_INSTANTIATE
#include "MatrixS_MatMul_Native_impl.h"

#include <boost/preprocessor/repetition/repeat_from_to.hpp>

namespace Metris
{
namespace DLA
{

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,2, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,3, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,4, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,5, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,6, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,7, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,3, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,4, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,5, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,6, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,7, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,2, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,3, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,4, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,5, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,6, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,7, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 1,4, 4,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,4, 4,3, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,4, 4,4, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,4, 4,5, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,4, 4,6, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,4, 4,7, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 1,5, 5,4, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,5, 5,5, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,5, 5,6, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,5, 5,7, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 1,6, 6,5, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,6, 6,6, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,6, 6,7, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 1,7, 7,6, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 1,7, 7,7, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 1,8, 8,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,3, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,4, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,5, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,6, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,7, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,3, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 2,3, 3,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 2,3, 3,2, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 2,3, 3,3, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 2,4, 4,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 2,5, 5,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real ); // FP-IBL3 coupling
METRIS_MATRIXS_MATMUL_NATIVE( 2,6, 6,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real ); // FP-IBL3 coupling
METRIS_MATRIXS_MATMUL_NATIVE( 2,8, 8,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 3,5, 5,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 4,5, 5,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 5,4, 4,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );


METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,2, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,3, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,4, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,5, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,6, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,7, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 3,2, 2,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 3,2, 2,2, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 3,2, 2,3, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,2, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 3,8, 8,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

// begin pcaplan
METRIS_MATRIXS_MATMUL_NATIVE( 3,4, 4,4, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,4, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
// end pcaplan
METRIS_MATRIXS_MATMUL_NATIVE( 3,4, 4,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 4,1, 1,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 4,1, 1,3, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 4,1, 1,4, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 4,1, 1,5, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 4,1, 1,6, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 4,1, 1,7, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 4,2, 2,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real ); // for neural net with input dim 2 and width 4
METRIS_MATRIXS_MATMUL_NATIVE( 4,3, 3,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 4,8, 8,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real ); // for 2D IBLAmpLagSplit

METRIS_MATRIXS_MATMUL_NATIVE( 5,1, 1,4, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 5,1, 1,5, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 5,1, 1,6, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 5,1, 1,7, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 5,2, 2,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real ); // for neural net with input dim 2 and width 5
METRIS_MATRIXS_MATMUL_NATIVE( 5,3, 3,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real ); // for neural net with input dim 2 and width 5

METRIS_MATRIXS_MATMUL_NATIVE( 5,5, 5,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 5,6, 6,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 6,1, 1,5, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 6,1, 1,6, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 6,1, 1,7, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 6,2, 2,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real ); // FP-IBL3 coupling
METRIS_MATRIXS_MATMUL_NATIVE( 6,5, 5,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 6,6, 6,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 6,8, 8,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 7,1, 1,6, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 7,1, 1,7, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 7,7, 7,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 8,2, 2,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 8,3, 3,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 8,4, 4,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 8,6, 6,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 8,8, 8,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

METRIS_MATRIXS_MATMUL_NATIVE( 20,11, 11,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real ); // for Prt NN with 11 inputs and 20 neurons
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real);
METRIS_MATRIXS_MATMUL_NATIVE( 1,20, 20,1, Metris::Real, Metris::Real, Metris::Real, Metris::Real );
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,20, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

//Explicit instantiation of square matrix multiplications
#define METRIS_MATMUL_DECL(z, n, text) \
METRIS_MATRIXS_MATMUL_NATIVE( n,n, n,n, Metris::Real, Metris::Real, Metris::Real, Metris::Real );

BOOST_PP_REPEAT_FROM_TO(1, 16, METRIS_MATMUL_DECL, )

} //namespace DLA
} //namespace Metris
