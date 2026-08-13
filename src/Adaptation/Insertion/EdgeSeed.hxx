//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_INSERTION_EDGE_SEED__
#define __METRIS_INSERTION_EDGE_SEED__

#include <egads.h>

namespace Metris{

class MshCavity;
class MeshBase;

struct EdgeSeed{

  EdgeSeed() = delete;
  // icollapse rejects edges a collapse cannot merge, then seeds the balls of both
  // ends, a collapse removing them. tdimp is read from the seeded cavity, so it
  // must be complete here. The caller must check ierro before using the seed.
  EdgeSeed(MeshBase& msh, MshCavity& cav, int tdim_adp, int tdim_ent, int ientt, int iedl,
           bool icollapse = false, int ithrd = 0);

  int ierro; // INS2D_ERR_*, 0 if the seed is usable
  int tdim_adp; // Context of insertion
  int tdimp; // Topo dim of a point on this edge
  int iseed; // Lowest dimensional entity that contains edge
  int iref;  // Ref of iseed  
  int ipedg[2]; // Edge points
  ego obj;   // If msh.CAD(), the corresponding CAD object.
  MshCavity& cav; // Shell of edge
};



}// namespace Metris

#endif