//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LOW_INCREASECAV__
#define __METRIS_LOW_INCREASECAV__

#include "../Mesh/MeshFwd.hxx"
#include "../types_arrays.hxx"
#include "../aux_hashtab.hxx"
#include <unordered_set>

/*
Cavity extension routines:
- increase_cavity increases for validity and Delaunay
- increase_cavity_lenedg increases/rejects for short edges

The point ipins must be constructed by newpoitopo followed by newbpotopo 
of the appropriate dimension.
*/


namespace Metris{

class MshCavity;
struct CavOprOpt;
struct EdgeSeed;


template<class MFT>
int movePointCavLen(Mesh<MFT>& msh, const MshCavity &cav, int tdim, int iseed, int miter, int ithrd1);
template<class MFT>
int movePointCavLen(Mesh<MFT>& msh, const MshCavity &cav, const EdgeSeed &insertionSeed, int miter, int ithrd1);

template<class MFT>
int setCavityInsertion(Mesh<MFT>& msh, MshCavity &cav, const CavOprOpt &opts, 
                       const EdgeSeed &insertionSeed, int mgrow, double lenqua_short_max, 
                       std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp,
                       int ithrd1, int ithrd2);

template<class MFT>
int setCavityInsertion2(Mesh<MFT>& msh, MshCavity &cav, 
                        const EdgeSeed &insertionSeed,
                        int miter, int ithrd1, int ithrd2);
                        
template<class MFT>
int setCavityInsertion2(Mesh<MFT>& msh, MshCavity &cav, const CavOprOpt &opts, 
                       int mgrow, double lenqua_short_max, 
                       std::unordered_set<std::tuple<int,int>,tup2_hash::hash> nocomp,
                       int ithrd1, int ithrd2);


// Check if any removed points; only those > 1/sqrt(2) from ipins if chklen
// This can possibly be reworked to be faster, for now we check everything every
// time, even though this is called in iterative cavity building.
template<class MFT>
void check_cavity_rempoint(MeshMetric<MFT> &msh, MshCavity &cav, const CavOprOpt &opts,
                           intAr1 &lrempoi, bool chklen, int ithrd1);

// Increase for validity and Delaunay (if idelaunay == true) both. 
template<class MFT>
int increase_cavity(MeshMetric<MFT> &msh, MshCavity &cav, 
                    bool idelaunay, int ithrd1, int ithrd2);

// Increase cavity based on validity only 
int increase_cavity_validity(MeshBase &msh, MshCavity &cav, int ithread);

// Increase cavity for Delaunay criterion on ipoin 
// normal is only necessary if dimension 3 and cavity has faces
template<class MFT>
int increase_cavity_Delaunay(MeshMetric<MFT> &msh, MshCavity &cav, 
                             int tdim, int ngrow, int ithread);

// Increase cavity to avoid short edges (add pts to collapse)
// return nprem ++points to collapse
// If an edge ipins-ipoin is short:
// - if dim(ipoin) < dim(ipins), error (always)
// - if dim(ipoin) == dim(ipins), error iff !opts.allow_remove_points
// - if dim(ipoin) > dim(ipins), error iff !opts.allow_remove_points_superdim
// I.e. set opts.allow_remove_points = false 
//      and opts.allow_remove_points_superdim = true
// to force a boundary operation except if it would conflict with boundary points.
template<class MFT>
int increase_cavity_lenedg(MeshMetric<MFT> &msh, MshCavity &cav, const CavOprOpt &opts, 
                              int ipins, int ithrd1, int ithrd2);
template<class MFT, int gdim>
int increase_cavity_lenedg0(MeshMetric<MFT> &msh, MshCavity &cav, const CavOprOpt &opts, 
                              int ipins, int ithrd1, int ithrd2);

// Tag cav.ipins's surface references if any. Used to filter in the other routines
void aux_taginsrefs(MeshBase &msh, MshCavity &cav, int ithread);
  
} // end namespace





#endif