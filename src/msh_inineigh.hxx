//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __SRC_MSH_ININEIGH__
#define __SRC_MSH_ININEIGH__


namespace Metris{

class MeshBase;

// Call this one
template<int ideg> void iniMeshNeighbours(MeshBase &msh);

template<int ideg> void iniMeshNeighbours2D(MeshBase &msh);
template<int ideg> void iniMeshNeighbours3D(MeshBase &msh);


} // End namespace

#endif
