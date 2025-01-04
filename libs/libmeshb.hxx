//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LIBMESHB_INTERFACE__
#define __METRIS_LIBMESHB_INTERFACE__

#include "../libs/libmeshb7.h"

#define __MAX_LIBMESHB_DEG__ 4

namespace Metris{
namespace libmeshb{


const int edgeKwds[1+4] = {0 , GmfEdges      , GmfEdgesP2      , GmfEdgesP3      , GmfEdgesP4};
const int faceKwds[1+4] = {0 , GmfTriangles  , GmfTrianglesP2  , GmfTrianglesP3  , GmfTrianglesP4};
const int elemKwds[1+4] = {0 , GmfTetrahedra , GmfTetrahedraP2 , GmfTetrahedraP3 , GmfTetrahedraP4};

const int edgeOrdKwds[1+4] = {0 , 0 , GmfEdgesP2Ordering      , GmfEdgesP3Ordering      , GmfEdgesP4Ordering};
const int faceOrdKwds[1+4] = {0 , 0 , GmfTrianglesP2Ordering  , GmfTrianglesP3Ordering  , GmfTrianglesP4Ordering};
const int elemOrdKwds[1+4] = {0 , 0 , GmfTetrahedraP2Ordering , GmfTetrahedraP3Ordering , GmfTetrahedraP4Ordering};

}
}


#endif
