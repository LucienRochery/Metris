//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSHCHECKTOPO__
#define __METRIS_MSHCHECKTOPO__


namespace Metris{

class MeshBase;

void check_topo(MeshBase &msh,
                int nbpoi, int npoin, int nedge, int nface, int nelem, int ithread = 0); 
void check_topo(MeshBase &msh, int ithread = 0); 

}// end namespace

#endif