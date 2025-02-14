//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LOW_INCREASECAV__
#define __METRIS_LOW_INCREASECAV__

#include "../Mesh/MeshFwd.hxx"


namespace Metris{

class MshCavity;
struct CavOprOpt;

// Increase cavity based on validity only 
int increase_cavity2D(MeshBase &msh, MshCavity &cav, int ithread);


// Increase cavity for Delaunay criterion on ipoin 
template<class MFT>
int increase_cavity_Delaunay(MeshMetric<MFT> &msh, MshCavity &cav, int ithread);

// Increase cavity to avoid short edges (add pts to collapse)
// return nprem ++points to collapse
template<class MFT>
int increase_cavity_lenedg(MeshMetric<MFT> &msh, MshCavity &cav, CavOprOpt &opts, 
                              int ipins, int ithrd1, int ithrd2);
template<class MFT, int gdim>
int increase_cavity_lenedg0(MeshMetric<MFT> &msh, MshCavity &cav, 
                              int ipins, int ithrd1, int ithrd2);

// Tag cav.ipins's surface references if any. Used to filter in the other routines
void aux_taginsrefs(MeshBase &msh, MshCavity &cav, int ithread);
  
} // end namespace





#endif