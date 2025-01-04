//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_LINEFORCE__
#define __METRIS_MSH_LINEFORCE__


#include "../Mesh/MeshFwd.hxx"

namespace Metris{

// Evaluate points on lines and forcibly reinsert them in the mesh, breaking 
// everything on the way 
template<class MFT>
void reinsertLines(Mesh<MFT> &msh, int ithrd1 = 0, int ithrd2 = 1);

}// end Namespace

#endif