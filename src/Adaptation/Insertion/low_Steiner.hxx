//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LOW_STEINER__
#define __METRIS_LOW_STEINER__

#include "../../Mesh/MeshFwd.hxx"
#include "../../types_arrays.hxx"

namespace Metris{

class MshCavity;
struct EdgeSeed;
struct CavWrkArrs;

// Return 0 if done nothing, 1 if error, -1 if inserted new point
template<class MFT>
int insertSteiner(Mesh<MFT>& msh, 
                  const EdgeSeed &insertionSeed,
                  MshCavity &cav, CavWrkArrs &work, 
                  intAr1 &lcaverr, int ithrd1, int ithrd2);



} // end namespace

#endif