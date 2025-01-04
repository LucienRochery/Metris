//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_MSH_COLLAPSE__
#define __METRIS_MSH_COLLAPSE__

#include "../Mesh/MeshFwd.hxx"


namespace Metris{

struct swapOptions;

// Collapse edge iedl of triangle iface
template<class MFT,int gdim,int ideg>
int swapface(Mesh<MFT>& msh, int iface, swapOptions opt,
             double *qumx0, double *qumx1, int ithread = 0);


} // end namespace

#endif