//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __SRC_LOWCCOEF_D__
#define __SRC_LOWCCOEF_D__

#include "types.hxx"

namespace Metris{

class MeshBase;

/* 
Linear components (derivatives) of Bézier coefficients of the Jacobian determinant
*/


void getccoef_map(const MeshBase &msh, int ientt, double *ccoef, dblAr2 &d_ccoef);
void rev_getccoef_map(const MeshBase &msh, int ientt, double *ccoef, dblAr2 &d_ccoef);


}

#endif