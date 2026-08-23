// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define METRIS_MATRIXS_MATMUL_NATIVE_INSTANTIATE
#include "MatrixS_MatMul_Native_impl.h"

#include "../../../../Surreal/SurrealS.h"
#include "../VectorS.h"

//#include "metris_constants.hxx"

namespace Metris
{
using Metris::SurrealS;

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


METRIS_MATRIXS_MATMUL_NATIVE( 2, 2, 2, 2, SurrealS<2>, SurrealS<2>, double, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 2, 2, 2, 2,      double, SurrealS<2>, double, SurrealS<2> );

METRIS_MATRIXS_MATMUL_NATIVE( 2, 2, 2, 3, SurrealS<2>, SurrealS<2>, double, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 2, 2, 2, 3,      double, SurrealS<2>, double, SurrealS<2> );

METRIS_MATRIXS_MATMUL_NATIVE( 2, 2, 2, 3, SurrealS<3>, SurrealS<3>, double, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 2, 2, 2, 3,      double, SurrealS<3>, double, SurrealS<3> );

METRIS_MATRIXS_MATMUL_NATIVE( 3, 3, 3, 3, SurrealS<3>, SurrealS<3>, double, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 3, 3, 3, 3,      double, SurrealS<3>, double, SurrealS<3> );


#if 0 

//---- SurrealS = Metris::Real * Metris::Real x Metris::Real ----------------------------------------
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,4, Metris::Real, Metris::Real, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 5,5, 5,5, Metris::Real, Metris::Real, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 6,6, 6,6, Metris::Real, Metris::Real, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 7,7, 7,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<1> );

METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Metris::Real, Metris::Real, Metris::Real, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Metris::Real, Metris::Real, Metris::Real, SurrealS<2> );

METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Metris::Real, Metris::Real, Metris::Real, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Metris::Real, Metris::Real, Metris::Real, SurrealS<3> );

METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Metris::Real, Metris::Real, Metris::Real, SurrealS<4> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<4> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,4, Metris::Real, Metris::Real, Metris::Real, SurrealS<4> );

METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,3, Metris::Real, Metris::Real, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 5,5, 5,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 5,5, 5,5, Metris::Real, Metris::Real, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,11, 11,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<5> ); // for Prt NN with 11 inputs and 20 neurons
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<5>);
METRIS_MATRIXS_MATMUL_NATIVE( 1,20, 20,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,20, Metris::Real, Metris::Real, Metris::Real, SurrealS<5> );

METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Metris::Real, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Metris::Real, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,3, Metris::Real, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 6,6, 6,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 6,6, 6,6, Metris::Real, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,11, 11,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<6> ); // for Prt NN with 11 inputs and 20 neurons
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<6>);
METRIS_MATRIXS_MATMUL_NATIVE( 1,20, 20,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,20, Metris::Real, Metris::Real, Metris::Real, SurrealS<6> );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<7> );
METRIS_MATRIXS_MATMUL_NATIVE( 7,7, 7,1, Metris::Real, Metris::Real, Metris::Real, SurrealS<7> );


//---- SurrealS = Metris::Real * Metris::Real x SurrealS ----------------------------------------
METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, Metris::Real, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, Metris::Real, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Metris::Real, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 7,7, 7,1, Metris::Real, SurrealS<1>, Metris::Real, SurrealS<1> );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, Metris::Real, SurrealS<2>, Metris::Real, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, SurrealS<2>, Metris::Real, SurrealS<2> );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, Metris::Real, SurrealS<3>, Metris::Real, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, SurrealS<3>, Metris::Real, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Metris::Real, SurrealS<3>, Metris::Real, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, SurrealS<3>, Metris::Real, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Metris::Real, SurrealS<3>, Metris::Real, SurrealS<3> );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, SurrealS<4>, Metris::Real, SurrealS<4> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, SurrealS<4>, Metris::Real, SurrealS<4> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,1, Metris::Real, SurrealS<4>, Metris::Real, SurrealS<4> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,4, Metris::Real, SurrealS<4>, Metris::Real, SurrealS<4> );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, SurrealS<5>, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,3, Metris::Real, SurrealS<5>, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, SurrealS<5>, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,3, 3,1, Metris::Real, SurrealS<5>, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,1, Metris::Real, SurrealS<5>, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 5,3, 3,1, Metris::Real, SurrealS<5>, Metris::Real, SurrealS<5> ); // for neural net with input dim 2 and width 5
METRIS_MATRIXS_MATMUL_NATIVE( 5,5, 5,1, Metris::Real, SurrealS<5>, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,11, 11,1, Metris::Real, SurrealS<5>, Metris::Real, SurrealS<5> ); // for Prt NN with 11 inputs and 20 neurons
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,1, Metris::Real, SurrealS<5>, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,20, 20,1, Metris::Real, SurrealS<5>, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,20, Metris::Real, SurrealS<5>, Metris::Real, SurrealS<5> );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,1, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,3, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,3, 3,1, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,1, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 5,3, 3,1, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> ); // for neural net with input dim 2 and width 5
METRIS_MATRIXS_MATMUL_NATIVE( 5,5, 5,1, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 6,6, 6,1, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 6,6, 6,6, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,11, 11,1, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> ); // for Prt NN with 11 inputs and 20 neurons
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,1, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,20, 20,1, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,20, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
// METRIS_MATRIXS_MATMUL_NATIVE( 11,20, 1,11, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> ); // for Prt NN with 11 inputs and 20 neurons
// METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 1,20, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );
// METRIS_MATRIXS_MATMUL_NATIVE( 20,1, 1,20, Metris::Real, SurrealS<6>, Metris::Real, SurrealS<6> );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, SurrealS<7>, Metris::Real, SurrealS<7> );
METRIS_MATRIXS_MATMUL_NATIVE( 7,7, 7,1, Metris::Real, SurrealS<7>, Metris::Real, SurrealS<7> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,11, 11,1, Metris::Real, SurrealS<7>, Metris::Real, SurrealS<7> ); // for Prt NN with 11 inputs and 20 neurons
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,1, Metris::Real, SurrealS<7>, Metris::Real, SurrealS<7>);
METRIS_MATRIXS_MATMUL_NATIVE( 1,20, 20,1, Metris::Real, SurrealS<7>, Metris::Real, SurrealS<7> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,20, Metris::Real, SurrealS<7>, Metris::Real, SurrealS<7> );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, SurrealS<8>, Metris::Real, SurrealS<8> );
METRIS_MATRIXS_MATMUL_NATIVE( 8,8, 8,1, Metris::Real, SurrealS<8>, Metris::Real, SurrealS<8> );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Metris::Real, SurrealS<9>, Metris::Real, SurrealS<9> );

METRIS_MATRIXS_MATMUL_NATIVE( 20,11, 11,1, Metris::Real, SurrealS<10>, Metris::Real, SurrealS<10> ); // for Prt NN with 11 inputs and 20 neurons
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,1, Metris::Real, SurrealS<10>, Metris::Real, SurrealS<10> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,20, 20,1, Metris::Real, SurrealS<10>, Metris::Real, SurrealS<10> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,20, Metris::Real, SurrealS<10>, Metris::Real, SurrealS<10> );

METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,3, Metris::Real, SurrealS<12>, Metris::Real, SurrealS<12> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Metris::Real, SurrealS<12>, Metris::Real, SurrealS<12> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Metris::Real, SurrealS<12>, Metris::Real, SurrealS<12> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,11, 11,1, Metris::Real, SurrealS<12>, Metris::Real, SurrealS<12> ); // for Prt NN with 11 inputs and 20 neurons
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,1, Metris::Real, SurrealS<12>, Metris::Real, SurrealS<12>);
METRIS_MATRIXS_MATMUL_NATIVE( 1,20, 20,1, Metris::Real, SurrealS<12>, Metris::Real, SurrealS<12> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,20, Metris::Real, SurrealS<12>, Metris::Real, SurrealS<12> );

METRIS_MATRIXS_MATMUL_NATIVE( 20,11, 11,1, Metris::Real, SurrealS<14>, Metris::Real, SurrealS<14> ); // for Prt NN with 11 inputs and 20 neurons
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,1, Metris::Real, SurrealS<14>, Metris::Real, SurrealS<14>);
METRIS_MATRIXS_MATMUL_NATIVE( 1,20, 20,1, Metris::Real, SurrealS<14>, Metris::Real, SurrealS<14> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,20, Metris::Real, SurrealS<14>, Metris::Real, SurrealS<14> );

METRIS_MATRIXS_MATMUL_NATIVE( 1,4, 4,4, Metris::Real, SurrealS<20>, Metris::Real, SurrealS<20> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,4, Metris::Real, SurrealS<20>, Metris::Real, SurrealS<20> );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Metris::Real, SurrealS<24>, Metris::Real, SurrealS<24> );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Metris::Real, SurrealS<48>, Metris::Real, SurrealS<48> );

typedef SurrealS<1,SurrealS<1>> SurrealS1_1;
typedef SurrealS<3,SurrealS<3>> SurrealS3_3;
typedef SurrealS<4,SurrealS<4>> SurrealS4_4;
typedef SurrealS<5,SurrealS<5>> SurrealS5_5;
typedef SurrealS<9,SurrealS<9>> SurrealS9_9;
typedef SurrealS<12,SurrealS<12>> SurrealS12_12;
typedef SurrealS<24,SurrealS<24>> SurrealS24_24;
typedef SurrealS<48,SurrealS<48>> SurrealS48_48;

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Metris::Real, SurrealS3_3  , Metris::Real, SurrealS3_3   );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Metris::Real, SurrealS9_9  , Metris::Real, SurrealS9_9   );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, Metris::Real, SurrealS12_12, Metris::Real, SurrealS12_12 );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, Metris::Real, SurrealS4_4  , Metris::Real, SurrealS4_4   );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Metris::Real, SurrealS24_24, Metris::Real, SurrealS24_24 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, Metris::Real, SurrealS48_48, Metris::Real, SurrealS48_48 );

METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,1, Metris::Real, SurrealS5_5  , Metris::Real, SurrealS5_5   );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS1_1, Metris::Real,  Metris::Real, SurrealS1_1 );

//---- SurrealS = SurrealS * Metris::Real x Metris::Real ----------------------------------------
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,3, Metris::Real, Metris::Real, SurrealS<1>, SurrealS<1> );

//---- SurrealS = Metris::Real * int x SurrealS -----------------------------------------
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<1>, int, Metris::Real, SurrealS<1> );

//---- SurrealS = Metris::Real * SurrealS x Metris::Real ----------------------------------------
METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, SurrealS<1>, Metris::Real, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, SurrealS<1>, Metris::Real, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<1>, Metris::Real, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<1>, Metris::Real, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<1>, Metris::Real, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,1, SurrealS<1>, Metris::Real, Metris::Real, SurrealS<1> );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, SurrealS<2>, Metris::Real, Metris::Real, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<2>, Metris::Real, Metris::Real, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<2>, Metris::Real, Metris::Real, SurrealS<2> );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, SurrealS<3>, Metris::Real, Metris::Real, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<3>, Metris::Real, Metris::Real, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<3>, Metris::Real, Metris::Real, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<3>, Metris::Real, Metris::Real, SurrealS<3> );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, SurrealS<4>, Metris::Real, Metris::Real, SurrealS<4> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, SurrealS<4>, Metris::Real, Metris::Real, SurrealS<4> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<4>, Metris::Real, Metris::Real, SurrealS<4> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<4>, Metris::Real, Metris::Real, SurrealS<4> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,1, SurrealS<4>, Metris::Real, Metris::Real, SurrealS<4> );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, SurrealS<5>, Metris::Real, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<5>, Metris::Real, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<5>, Metris::Real, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 5,5, 5,1, SurrealS<5>, Metris::Real, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,11, 11,1, SurrealS<5>, Metris::Real, Metris::Real, SurrealS<5> ); // for Prt NN with 11 inputs and 20 neurons
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,1, SurrealS<5>, Metris::Real, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,20, 20,1, SurrealS<5>, Metris::Real, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,20, SurrealS<5>, Metris::Real, Metris::Real, SurrealS<5> );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, SurrealS<6>, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, SurrealS<6>, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<6>, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<6>, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,3, SurrealS<6>, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS<6>, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<6>, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 6,6, 6,1, SurrealS<6>, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 6,6, 6,6, SurrealS<6>, Metris::Real, Metris::Real, SurrealS<6> );
// METRIS_MATRIXS_MATMUL_NATIVE( 6,1, 1,1, SurrealS<6>, Metris::Real, Metris::Real, SurrealS<6> ); //dPrt mult
METRIS_MATRIXS_MATMUL_NATIVE( 20,11, 11,1, SurrealS<6>, Metris::Real, Metris::Real, SurrealS<6> ); // for Prt NN with 11 inputs and 20 neurons
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,1, SurrealS<6>, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,20, 20,1, SurrealS<6>, Metris::Real, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 20,20, 20,20, SurrealS<6>, Metris::Real, Metris::Real, SurrealS<6> );

METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,3, SurrealS<7>, Metris::Real, Metris::Real, SurrealS<7> );
METRIS_MATRIXS_MATMUL_NATIVE( 7,7, 7,1, SurrealS<7>, Metris::Real, Metris::Real, SurrealS<7> );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<8>, Metris::Real, Metris::Real, SurrealS<8> );

METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,1, SurrealS<12>, Metris::Real, Metris::Real, SurrealS<12> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS<12>, Metris::Real, Metris::Real, SurrealS<12> );

METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,4, SurrealS<20>, Metris::Real, Metris::Real, SurrealS<20> );

//---- SurrealS = Metris::Real * SurrealS x SurrealS -----------------------------------
METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, SurrealS<1>, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, SurrealS<1>, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, SurrealS<1>, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<1>, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<1>, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<1>, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,1, SurrealS<1>, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,2, SurrealS<1>, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,3, SurrealS<1>, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,4, SurrealS<1>, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 5,5, 5,1, SurrealS<1>, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 6,6, 6,1, SurrealS<1>, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 7,7, 7,1, SurrealS<1>, SurrealS<1>, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 8,8, 8,1, SurrealS<1>, SurrealS<1>, Metris::Real, SurrealS<1> );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, SurrealS<2>, SurrealS<2>, Metris::Real, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<2>, SurrealS<2>, Metris::Real, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<2>, SurrealS<2>, Metris::Real, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<2>, SurrealS<2>, Metris::Real, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<9>, SurrealS<9>, Metris::Real, SurrealS<9> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<12>, SurrealS<12>, Metris::Real, SurrealS<12> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<2>, SurrealS<2>, Metris::Real, SurrealS<2> );

METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,1, SurrealS<2>, SurrealS<2>, Metris::Real, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,2, SurrealS<2>, SurrealS<2>, Metris::Real, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,3, SurrealS<2>, SurrealS<2>, Metris::Real, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,4, SurrealS<2>, SurrealS<2>, Metris::Real, SurrealS<2> );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<3>, SurrealS<3>, Metris::Real, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<3>, SurrealS<3>, Metris::Real, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<3>, SurrealS<3>, Metris::Real, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS<3>, SurrealS<3>, Metris::Real, SurrealS<3> );

METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,1, SurrealS<4>, SurrealS<4>, Metris::Real, SurrealS<4> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,2, SurrealS<4>, SurrealS<4>, Metris::Real, SurrealS<4> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,3, SurrealS<4>, SurrealS<4>, Metris::Real, SurrealS<4> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,4, SurrealS<4>, SurrealS<4>, Metris::Real, SurrealS<4> );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, SurrealS<5>, SurrealS<5>, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,3, SurrealS<5>, SurrealS<5>, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 5,5, 5,1, SurrealS<5>, SurrealS<5>, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 5,5, 5,5, SurrealS<5>, SurrealS<5>, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 5,5, 5,5,        Metris::Real, SurrealS<5>, Metris::Real, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 5,5, 5,5, SurrealS<5>,        Metris::Real, Metris::Real, SurrealS<5> );

METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<6>, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, SurrealS<6>, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, SurrealS<6>, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<6>, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS<6>, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS<6>, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 6,6, 6,1, SurrealS<6>, SurrealS<6>, Metris::Real, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 6,6, 6,6, SurrealS<6>, SurrealS<6>, Metris::Real, SurrealS<6> );

METRIS_MATRIXS_MATMUL_NATIVE( 6,6, 6,1,        Metris::Real,        Metris::Real, Metris::Real, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 6,6, 6,1,        Metris::Real, SurrealS<1>, Metris::Real, SurrealS<1> );

METRIS_MATRIXS_MATMUL_NATIVE( 7,7, 7,1, SurrealS<7>, SurrealS<7>, Metris::Real, SurrealS<7> );
METRIS_MATRIXS_MATMUL_NATIVE( 7,7, 7,7, SurrealS<7>, SurrealS<7>, Metris::Real, SurrealS<7> );

METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<8>, SurrealS<8>, Metris::Real, SurrealS<8> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, SurrealS<8>, SurrealS<8>, Metris::Real, SurrealS<8> );
METRIS_MATRIXS_MATMUL_NATIVE( 8,8, 8,1, SurrealS<8>, SurrealS<8>, Metris::Real, SurrealS<8> );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, SurrealS<9>, SurrealS<9>, Metris::Real, SurrealS<9> );

METRIS_MATRIXS_MATMUL_NATIVE( 5,5, 5,1, SurrealS<10>, SurrealS<10>, Metris::Real, SurrealS<10> );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS<12>, SurrealS<12>, Metris::Real, SurrealS<12> );

METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,4, SurrealS<20>, SurrealS<20>, Metris::Real, SurrealS<20> );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS<24>, SurrealS<24>, Metris::Real, SurrealS<24> );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS<48>, SurrealS<48>, Metris::Real, SurrealS<48> );


METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS9_9  , SurrealS9_9  , Metris::Real, SurrealS9_9   );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,2, SurrealS12_12, SurrealS12_12, Metris::Real, SurrealS12_12 );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS24_24, SurrealS24_24, Metris::Real, SurrealS24_24 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, SurrealS48_48, SurrealS48_48, Metris::Real, SurrealS48_48 );


//---- SurrealS = SurrealS * Metris::Real x Metris::Real ------------------------------------
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Metris::Real, Metris::Real, SurrealS<2>, SurrealS<2> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Metris::Real, Metris::Real, SurrealS<3>, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Metris::Real, Metris::Real, SurrealS<4>, SurrealS<4> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Metris::Real, Metris::Real, SurrealS<5>, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,3, Metris::Real, Metris::Real, SurrealS<5>, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, Metris::Real, Metris::Real, SurrealS<6>, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,3, Metris::Real, Metris::Real, SurrealS<6>, SurrealS<6> );

//---- SurrealS = SurrealS * SurrealS x SurrealS ------------------------------------
METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, SurrealS<1>, SurrealS<1>, SurrealS<1>, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<1>, SurrealS<1>, SurrealS<1>, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,3, SurrealS<1>, SurrealS<1>, SurrealS<1>, SurrealS<1> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<3>, SurrealS<3>, SurrealS<3>, SurrealS<3> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<5>, SurrealS<5>, SurrealS<5>, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 3,1, 1,3, SurrealS<5>, SurrealS<5>, SurrealS<5>, SurrealS<5> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<6>, SurrealS<6>, SurrealS<6>, SurrealS<6> );
METRIS_MATRIXS_MATMUL_NATIVE( 2,1, 1,2, SurrealS<8>, SurrealS<8>, SurrealS<8>, SurrealS<8> );

//-----------------------------------------------------------------------------------
// Multiplications of block matrices
//-----------------------------------------------------------------------------------
typedef VectorS<1,Metris::Real> Vector1;
typedef VectorS<1,SurrealS<1>> VectorS1;
typedef MatrixS<1,1,SurrealS<1>> MatrixS11;
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS11, Vector1,  Metris::Real, VectorS1 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS11, VectorS1, Metris::Real, VectorS1 );

typedef VectorS<2,Metris::Real> Vector2;
typedef VectorS<2,SurrealS<2>> VectorS2;
typedef MatrixS<2,2,SurrealS<2>> MatrixS22;
METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS22, Vector2, Metris::Real, VectorS2 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS22, Vector2, Metris::Real, VectorS2 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS22, Vector2, Metris::Real, VectorS2 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS22, VectorS2, Metris::Real, VectorS2 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS22, VectorS2, Metris::Real, VectorS2 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS22, VectorS2, Metris::Real, VectorS2 );

typedef VectorS<3,Metris::Real> Vector3;
typedef VectorS<3,SurrealS<1>> Vector3_S1;
typedef MatrixS<3,3,SurrealS<1>> Matrix33_S1;
typedef VectorS<3,SurrealS<3>> VectorS_S3;
typedef MatrixS<3,3,SurrealS<3>> MatrixS33_S3;
METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS33_S3, MatrixS33_S3, Metris::Real, MatrixS33_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS33_S3, MatrixS33_S3, Metris::Real, MatrixS33_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS33_S3, MatrixS33_S3, Metris::Real, MatrixS33_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS33_S3, MatrixS33_S3, Metris::Real, MatrixS33_S3 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS33_S3, VectorS_S3, Metris::Real, VectorS_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS33_S3, VectorS_S3, Metris::Real, VectorS_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS33_S3, VectorS_S3, Metris::Real, VectorS_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS33_S3, VectorS_S3, Metris::Real, VectorS_S3 );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS33_S3, VectorS_S3, Metris::Real, VectorS_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,3, MatrixS33_S3, VectorS_S3, Metris::Real, VectorS_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,1, MatrixS33_S3, VectorS_S3, Metris::Real, VectorS_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS33_S3, VectorS_S3, Metris::Real, VectorS_S3 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS33_S3, Vector3, Metris::Real, VectorS_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS33_S3, Vector3, Metris::Real, VectorS_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS33_S3, Vector3, Metris::Real, VectorS_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS33_S3, Vector3, Metris::Real, VectorS_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS33_S3, Vector3, Metris::Real, VectorS_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,3, MatrixS33_S3, Vector3, Metris::Real, VectorS_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,3, 3,1, MatrixS33_S3, Vector3, Metris::Real, VectorS_S3 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,3, MatrixS33_S3, Vector3, Metris::Real, VectorS_S3 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, Matrix33_S1, Vector3, Metris::Real, Vector3_S1 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, Matrix33_S1, Vector3, Metris::Real, Vector3_S1 );

typedef VectorS<4,Metris::Real> VectorS4_R;
typedef VectorS<4,SurrealS<1>> VectorS4_S1;
typedef MatrixS<4,4,SurrealS<1>> MatrixS44_S1;
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44_S1, VectorS4_S1, Metris::Real, VectorS4_S1 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44_S1, VectorS4_R,  Metris::Real, VectorS4_S1 );

typedef VectorS<4,Metris::Real> Vector4;
typedef VectorS<4,SurrealS<4>> VectorS4;
typedef MatrixS<4,4,SurrealS<4>> MatrixS44;
METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS44, MatrixS44, Metris::Real, MatrixS44 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44, MatrixS44, Metris::Real, MatrixS44 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS44, MatrixS44, Metris::Real, MatrixS44 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS44, MatrixS44, Metris::Real, MatrixS44 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,1, 1,1, MatrixS44, VectorS4, Metris::Real, VectorS4 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44, VectorS4, Metris::Real, VectorS4 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS44, VectorS4, Metris::Real, VectorS4 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,1, MatrixS44, VectorS4, Metris::Real, VectorS4 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS44, VectorS4, Metris::Real, VectorS4 );

METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,1, MatrixS44, VectorS4, Metris::Real, VectorS4 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,4, 4,4, MatrixS44, VectorS4, Metris::Real, VectorS4 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,4, 4,1, MatrixS44, VectorS4, Metris::Real, VectorS4 );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4, 4,4, MatrixS44, VectorS4, Metris::Real, VectorS4 );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44, Vector4, Metris::Real, VectorS4 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS44, Vector4, Metris::Real, VectorS4 );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS44, Metris::Real, Metris::Real, MatrixS44 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS44, Metris::Real, Metris::Real, MatrixS44 );


typedef VectorS<5,Metris::Real> Vector5;
typedef VectorS<5,SurrealS<5>> VectorS5;
typedef MatrixS<5,5,SurrealS<5>> MatrixS55;

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS55, Metris::Real, Metris::Real, MatrixS55 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS55, Metris::Real, Metris::Real, MatrixS55 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS55, Vector5, Metris::Real, VectorS5 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS55, Vector5, Metris::Real, VectorS5 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS55, Vector5, Metris::Real, VectorS5 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS55, VectorS5, Metris::Real, VectorS5 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS55, VectorS5, Metris::Real, VectorS5 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS55, VectorS5, Metris::Real, VectorS5 );
METRIS_MATRIXS_MATMUL_NATIVE( 5,5, 5,1, MatrixS55, VectorS5, Metris::Real, VectorS5 );

typedef VectorS<5,Metris::Real> VectorS5_R;
typedef VectorS<5,SurrealS<1>> VectorS5_S1;
typedef MatrixS<5,5,SurrealS<1>> MatrixS55_S1;
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS55_S1, VectorS5_R, Metris::Real, VectorS5_S1 );

typedef VectorS<5,SurrealS<2>> VectorS5_S2;
typedef MatrixS<5,5,SurrealS<2>> MatrixS55_S2;
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS55_S2, VectorS5_R, Metris::Real, VectorS5_S2 );

typedef VectorS<5,SurrealS<3>> VectorS5_S3;
typedef MatrixS<5,5,SurrealS<3>> MatrixS55_S3;
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS55_S3, VectorS5_R, Metris::Real, VectorS5_S3 );

typedef VectorS<5,SurrealS<4>> VectorS5_S4;
typedef MatrixS<5,5,SurrealS<4>> MatrixS55_S4;
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS55_S4, VectorS5_R, Metris::Real, VectorS5_S4 );

typedef VectorS<5,SurrealS<10>> VectorS10_5;
typedef MatrixS<5,5,SurrealS<10>> MatrixS10_55;

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS10_55, Vector5, Metris::Real, VectorS10_5 );

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS10_55, VectorS10_5, Metris::Real, VectorS10_5 );

typedef VectorS<6,Metris::Real> Vector6;
typedef VectorS<6,SurrealS<6>> VectorS6;
typedef MatrixS<6,6,SurrealS<6>> MatrixS66;

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS66, Metris::Real, Metris::Real, MatrixS66 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS66, Vector6, Metris::Real, VectorS6);
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS66, Vector6, Metris::Real, VectorS6 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS66, Vector6, Metris::Real, VectorS6 );

METRIS_MATRIXS_MATMUL_NATIVE( 1,2, 2,2, MatrixS66, VectorS6, Metris::Real, VectorS6 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS66, VectorS6, Metris::Real, VectorS6 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS66, VectorS6, Metris::Real, VectorS6 );

typedef VectorS<6,SurrealS<1>> VectorS61;
typedef MatrixS<6,6,SurrealS<1>> MatrixS661;

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS661, Vector6, Metris::Real, VectorS61 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS661, Vector6, Metris::Real, VectorS61 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS661, VectorS61, Metris::Real, VectorS61 );

typedef VectorS<7,Metris::Real> Vector7;
typedef VectorS<7,SurrealS<7>> VectorS7;
typedef MatrixS<7,7,SurrealS<7>> MatrixS77;

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS77, Vector7, Metris::Real, VectorS7 );
METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS77, VectorS7, Metris::Real, VectorS7 );

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS77, Vector7, Metris::Real, VectorS7 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS77, VectorS7, Metris::Real, VectorS7 );

typedef VectorS<7,SurrealS<1>> VectorS71;
typedef MatrixS<7,7,SurrealS<1>> MatrixS771;

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS771, Vector7, Metris::Real, VectorS71 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS771, VectorS71, Metris::Real, VectorS71 );

typedef VectorS<8,Metris::Real> Vector8;
typedef VectorS<8,SurrealS<8>> VectorS8;
typedef MatrixS<8,8,SurrealS<8>> MatrixS88;

METRIS_MATRIXS_MATMUL_NATIVE( 2,2, 2,1, MatrixS88, VectorS8, Metris::Real, VectorS8 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS88, VectorS8, Metris::Real, VectorS8 );

typedef VectorS<8,SurrealS<1>> VectorS81;
typedef MatrixS<8,8,SurrealS<1>> MatrixS881;

METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS881, Vector8, Metris::Real, VectorS81 );
METRIS_MATRIXS_MATMUL_NATIVE( 3,3, 3,1, MatrixS881, VectorS81, Metris::Real, VectorS81 );

typedef MatrixS<4,4,Metris::Real> MatrixR44;
typedef MatrixS<4,4,SurrealS<10,Metris::Real>> MatrixS44_S10;
typedef MatrixS<1,4,Metris::Real> MatrixR14;
typedef MatrixS<1,4,SurrealS<10,Metris::Real>> MatrixS14_S10;
METRIS_MATRIXS_MATMUL_NATIVE( 4,4,4,4, MatrixS44_S10, MatrixR44 , Metris::Real, MatrixS44_S10 );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4,4,4, MatrixS44_S10, MatrixS44_S10 , Metris::Real, MatrixS44_S10 );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4,4,4, MatrixR44, MatrixS44_S10 , Metris::Real, MatrixS44_S10 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,4,4,4, MatrixR14, MatrixS44_S10, Metris::Real, MatrixS14_S10 );
METRIS_MATRIXS_MATMUL_NATIVE( 1,4,4,4, SurrealS<10>, SurrealS<10>, Metris::Real, SurrealS<10> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4,4,4, SurrealS<10>, SurrealS<10>, Metris::Real, SurrealS<10> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4,4,4, SurrealS<10>, Metris::Real, Metris::Real, SurrealS<10> );
METRIS_MATRIXS_MATMUL_NATIVE( 4,4,4,4, Metris::Real, SurrealS<10>, Metris::Real, SurrealS<10> );
METRIS_MATRIXS_MATMUL_NATIVE( 1,4,4,4, Metris::Real, SurrealS<10>, Metris::Real, SurrealS<10> );

METRIS_MATRIXS_MATMUL_NATIVE( 4,4,4,4, Metris::Real,SurrealS<50>,Metris::Real,SurrealS<50>);
METRIS_MATRIXS_MATMUL_NATIVE( 4,4,4,4, SurrealS<50>,SurrealS<50>,Metris::Real,SurrealS<50>);
#endif

} //namespace DLA
} //namespace Metris
