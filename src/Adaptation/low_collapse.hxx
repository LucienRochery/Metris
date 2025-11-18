//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LOW_COLLAPSE__
#define __METRIS_LOW_COLLAPSE__


#include "../Mesh/MeshFwd.hxx"
#include "../types.hxx"
#include "../aux_badEntHandler.hxx"


namespace Metris{

class MshCavity;
struct CavWrkArrs;


// Collapse edge iedl of triangle iface
template<class MFT>
int collapseEdge(Mesh<MFT>& msh, int tdim, int ientt, int iedl, double qmax_suf,
                 MshCavity &cav, CavWrkArrs &work,
                 intAr1 &lerro,
                 #ifdef TESTQUALITYALGO
                 BadEntHandler& handler,
                 #endif
                 int ithrd1, int ithrd2, int ithrd3);

template<class MFT>
int collapseEdge2(Mesh<MFT>& msh, int tdim, int ientt, int iedl, double qmax_suf,
                  MshCavity &cav, CavWrkArrs &work,
                  intAr1 &lerro, int ithrd1, int ithrd2, int ithrd3);

template<class MFT>
int collapseVertex(Mesh<MFT>& msh, int ipcol, double qmax_suf,
                   MshCavity &cav, CavWrkArrs &work,
                   intAr1 &lerro, int ithrd1, int ithrd2);



} // end namespace

#endif