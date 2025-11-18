//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_LOW_INSERT__
#define __METRIS_LOW_INSERT__

#include "../../Mesh/MeshFwd.hxx"
#include "../../types.hxx"

#include "insert_errors.hxx"
#include "../../aux_badEntHandler.hxx"

namespace Metris{

class MshCavity;
struct CavWrkArrs;
struct EdgeSeed;

// Collapse edge iedl of triangle iface
// bar1 is t along the edge with 1 if lnoed[iedl][0]
template<class MFT>
int insertEdge(Mesh<MFT>& msh,
               const EdgeSeed &insertionSeed,
               double lenqua_short_max,
               bool icollapse,
               MshCavity &cav, CavWrkArrs &work,
               intAr1 &lerro,
               #ifdef TESTQUALITYALGO
               BadEntHandler& handler,
               #endif
               int ithrd1, int ithrd2);


} // end namespace

#endif