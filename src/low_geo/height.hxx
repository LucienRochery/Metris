//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __LOW_HEIGHT__
#define __LOW_HEIGHT__

#include "../Mesh/MeshFwd.hxx"

namespace Metris{


template<class MFT, int tdim>
void getheightentP1_aniso(const Mesh<MFT> &msh, int ientt, double *height);


} // End namespace
#endif