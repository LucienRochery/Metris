//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_REINSERT_FLAT__
#define __METRIS_MSH_REINSERT_FLAT__

#include "../Mesh/MeshFwd.hxx"

namespace Metris{



// Reinsert vertices that create almost flat elements. 
template<class MFT, int gdim, int ideg>
int reinsertFlat(Mesh<MFT> &msh, bool allow_collapse, int ithrd1 = 0);


}//namespace Metris

#endif