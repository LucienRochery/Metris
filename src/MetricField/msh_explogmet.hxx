//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __MSH_METRIC__
#define __MSH_METRIC__

#include "../types.hxx"


namespace Metris{

class MeshBase;

template<int ndimn>
void setLogMetMesh0(const MeshBase &msh,  dblAr2 &metfld);

template<int ndimn>
void setExpMetMesh0(const MeshBase &msh,  dblAr2 &metfld);


} // End namespace


#endif