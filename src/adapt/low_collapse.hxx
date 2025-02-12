//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LOW_COLLAPSE__
#define __METRIS_LOW_COLLAPSE__


#include "../Mesh/MeshFwd.hxx"
#include "../types.hxx"


namespace Metris{

class MshCavity;
struct CavWrkArrs;


// Collapse edge iedl of triangle iface
template<class MFT>
int colledgsurf(Mesh<MFT>& msh, int iface, int iedl, double qmax_suf, 
                MshCavity &cav, CavWrkArrs &work, 
                intAr1 &lerro, int ithrd1 = 0, int ithrd2 = 1);


template<class MFT>
int collversurf(Mesh<MFT>& msh, int iface, int iver, double qmax_suf, 
                MshCavity &cav, CavWrkArrs &work, 
                intAr1 &lerro, int ithrd1 = 0, int ithrd2 = 1);


} // end namespace

#endif