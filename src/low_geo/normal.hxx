//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __LOW_NORMAL__
#define __LOW_NORMAL__

#include "../types.hxx"
#include "../Mesh/MeshFwd.hxx"

namespace Metris{

enum class AsDeg;


// ------------------------------------------------------------------------
// Compute normals
// ------------------------------------------------------------------------
void getnorfacP1(const int *fac2pol, const dblAr2 &coord, double *nrmal);

// This should be preferred if the face is in the mesh:
void getnorfac(const MeshBase&__restrict__ msh, int iface, 
               const double*__restrict__ bary, AsDeg asdmsh, double*__restrict__ nrmal);
// This is only for temporary faces:
void getnorfac(const MeshBase&__restrict__ msh, const int*__restrict__ fac2pol, 
               const double*__restrict__ bary, AsDeg asdmsh, double*__restrict__ nrmal);

// Average normals at the vertices
// To compute the CAD normal, the safest is to average the vertex normals.
// This is because taking the average of the (u,v)'s can send us just about
// anywhere.
int getnorfacCAD(const MeshBase &msh, int iface, double *nrmal);

// Return outgoing normal of edge (2D only). 
int getnorpoiCAD1(const MeshBase &msh, int ipoin, std::map<ego,int> &edgorient, 
                  double *norpoi);

int getnorpoiCAD2(const MeshBase &msh, int ibpoi, double *norpoi);

// iref >= 0 filters, iref < 0 ignored
template <int ideg>
void getnorballref(MeshBase &msh, const intAr1 &lball, int iref, double* norpoi);

void getnorpoiref(const MeshBase &msh, int ipoin, int iref, double* norpoi);

// Specify either iedg0 or iref0 (or both). Handles both CAD and discrete case. 
// - if iedg0 >= 0 && iref0 < 0: get tangent from iedg0 side 
//   (in CAD case, no difference if not corner)
int gettanpoiref(const MeshBase &msh, int ipoin, int iref, double* tanpoi);



// ------------------------------------------------------------------------
// Compute normal deviation
// ------------------------------------------------------------------------
// Normalized dotprod difference to 1 of CAD/elt normals accumulated over the nodes
// - norfac can be provided if already computed by caller, only works in P1 and not necessary.
template<int ideg>
double getnordev(const MeshBase&__restrict__ msh, int iface, const double*__restrict__ norfac = NULL);

// This is only for temporary faces:
// nod2bpo holds the ibpois for the element nodes
template<int ideg>
double getnordev(const MeshBase&__restrict__ msh, const int*__restrict__ fac2pol,
                 const int*__restrict__ nod2bpo, const double*__restrict__ norfac);

}// namespace

#endif