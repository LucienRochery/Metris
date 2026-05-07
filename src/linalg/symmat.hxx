//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_SYMMAT__
#define __METRIS_SYMMAT__

#include "symidx.hxx"
#include "../libs/SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"


namespace SANS::DLA{

//template<int M, class T>
//using MatSymS = VectorS<(M*(M+1))/2, T>;

template<int M, class T>
class MatSymS : public VectorS<(M*(M+1))/2, T>{
};

}


#endif