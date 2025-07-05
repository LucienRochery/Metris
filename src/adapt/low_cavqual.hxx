//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LOW_CAVQUAL__
#define __METRIS_LOW_CAVQUAL__


#include "../Mesh/MeshFwd.hxx"
#include "../types.hxx"

#include <unordered_set>

namespace Metris{

class MshCavity;
struct CavWrkArrs;

// Reject proposed cavity based on edge length score same as in swaps. 
// If filter_long/short is set, then the output edges that are long/short
// are ignored. In the case of an insertion, new long edges can be ignored. 
// Set grow_check to true to inspect elements outside the cavity. Only useful
// for insertions. 
// nocomp is work array
template<class MFT>
int collrejcav_lenqua(Mesh<MFT>& msh, MshCavity &cav, 
                      bool filter_long, bool filter_short, bool grow_check,
                      double lenqua_short_max,
                      std::unordered_set<std::tuple<int,int>,tup2_hash::hash> &nocomp,
                      int ithrd1);


template<class MFT>
int collrejcav_dens(Mesh<MFT>& msh, MshCavity &cav, int ithrd1, int ithrd2);


//// Reject proposed cavity based on edge length:
//// if more long edges are created than short edges are destroyed, reject
//template<class MFT>
//int collrejcav_len(Mesh<MFT>& msh, MshCavity &cav, int ithrd1);




}// namespace
#endif