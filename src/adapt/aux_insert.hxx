//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_AUX_INSERT__
#define __METRIS_AUX_INSERT__

#include "../Mesh/MeshFwd.hxx"
#include "../types.hxx"

namespace Metris{

class MshCavity;
struct CavWrkArrs;

template<class MFT>
int aux_bisecPointLen(Mesh<MFT> &msh, 
                      int tdim, int ientt, int iedl,
                      int ibins, int tdimp, int iseed, int iref,
                      bool icollapse,
                      const MshCavity &cav);

// Correct point location in case of cavity construction error (e.g. short edge)
template<class MFT>
int aux_movePointCav(Mesh<MFT>& msh, MshCavity &cav, 
                     int tdimp, int iseed, int iref, double *algnd);

template<class MFT>
int aux_findCloseConstrained(Mesh<MFT>& msh, MshCavity &cav, 
                             int ithrd1, int ithrd2);

}
#endif