//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_LINEADAPT__
#define __METRIS_MSH_LINEADAPT__


#include "../Mesh/MeshFwd.hxx"
#include "../types.hxx"




namespace Metris{

class MshCavity;



template<class MFT>
void adaptGeoLines(Mesh<MFT> &msh, int ithrd1 = 0, int ithrd2 = 1);
template<class MFT>
void getCADCurveLengths(Mesh<MFT> &msh, double tol, dblAr1 &crv_len);


/* ---------------------------------------------
// Functions auxiliary to msh_lineadapt.hxx 
// Functions that don't have a use outside of breaking up adaptGeoLines
// and making it more readable. 
   --------------------------------------------- */

// From ipoistart, iedgstart, keep walking along the curve until adjusted_tarlen
// has been met. 
// Return edg2pol with first and last point seen, ibpend the ibpo of 
// last pt, ifaceend the face adjacent to iedgend, in turn the last seen edge.
template<class MFT>
static int aux_walk_line(Mesh<MFT> &msh, MshCavity& cav, ego obj, 
                         int ipoistart, int iedgstart, double adjusted_tarlen, 
                         int *edg2pol, int *ibpend, int *iedgend, int *ifacend, 
                         double *lenend, double *szend, bool *ifin, 
                         double *erredg, int *npterr, int *nEGrro,
                         int ithrd1);

// In case an insertion is needed: generate new point using bisection.

template<class MFT>
static int gen_newp_line(Mesh<MFT> &msh, MshCavity& cav, ego obj,
                         int ibpo0, int ibpnw,
                         int iedgseed, int ifacseed, const int *edg2pol, 
                         double *sz, int miter_bisection,
                         double tarlen, double lentolabs,
                         double *adjusted_tarlen,  // in/out
                         double *lenend, int *nEGrro);


}// end Namespace


#endif