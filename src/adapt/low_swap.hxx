//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifndef __METRIS_LOW_SWAP__
#define __METRIS_LOW_SWAP__

#include "../Mesh/MeshFwd.hxx"


namespace Metris{

class MshCavity;
struct CavWrkArrs;
struct CavOprOpt;
struct swapOptions;

// Swap edge between two triangles (including surface w/ tets)
// looking for best swap involving iface (low quality triangle). 
// Return 0 if nothing done, 1 if error, -1 if swap done
// Compute using norm specified in opt: if 0, take max.
// If norm is -1, use edge length instead.
// low_swapface.cxx
template<class MFT, int gdim, int ideg>
int swapface(Mesh<MFT>& msh, int iface, swapOptions opt,
             MshCavity &cav, CavWrkArrs &work, 
             double *qumx0, double *qumx1, int ithread);


// low_swaptetra.cxx
template<class MFT,int ideg>
int swaptetra(Mesh<MFT>& msh, int itetr, swapOptions opt,
              MshCavity &cav, CavWrkArrs &work, 
              double *qumx0, double *qumx1, int ithrd1, int ithrd2);

template<class MFT,int ideg>
int aux_swaptetface(Mesh<MFT>& msh, swapOptions opt, int itetr, int ifacl, double quae1,
                    MshCavity &cav, CavOprOpt &opts, CavWrkArrs &work,
                    double *qnrm0_, double *qnrm1_,
                    int ithread);

template<class MFT,int ideg>
int aux_swaptetedge(Mesh<MFT>& msh, swapOptions opt, int itetr, int ifacl, double quae1,
                    MshCavity &cav, CavOprOpt &opts, CavWrkArrs &work,
                    double *qnrm0_, double *qnrm1_,
                    int ithrd1, int ithrd2);

} // end namespace

#endif