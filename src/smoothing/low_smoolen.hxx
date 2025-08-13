//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php
#ifndef __METRIS_LOW_SMOOLEN__
#define __METRIS_LOW_SMOOLEN__

#include "../Mesh/MeshFwd.hxx"
#include "../types_arrays.hxx"

namespace Metris{

class MshCavity;
struct EdgeSeed;

template<class MFT>
int movePointCavLen(Mesh<MFT>& msh, const MshCavity &cav, int tdim, int iseed, int miter, int ithrd1);
template<class MFT>
int movePointCavLen(Mesh<MFT>& msh, const MshCavity &cav, const EdgeSeed &insertionSeed, int miter, int ithrd1);

template<class MFT>
int smoopoilen(Mesh<MFT>& msh, int ipmov, const intAr1 &lpoin, int miter, int tdimp = -1, int iseed = -1);

}//namespace Metris
#endif