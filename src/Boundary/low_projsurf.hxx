//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LOW_PROJSURF__
#define __METRIS_LOW_PROJSURF__


namespace Metris{

class MeshBase;

// Project on P1 edge: no CAD used. 
template <int gdim, int ideg>
int projptedg(MeshBase &msh, const double*__restrict__ coop, 
              int iedge, 
              double*__restrict__ bary,
              double*__restrict__ coopr);

template <int gdim, int ideg>
int projptedg(MeshBase &msh, const double*__restrict__ coop, 
              const int*__restrict__ edg2pol, 
              double*__restrict__ bary,
              double*__restrict__ coopr);


// Project point coop on edge iedge, return barycentric of proj and CAD param
// Requires msh to be linked to CAD object. 
template <int gdim>
int projptedg(MeshBase &msh, const double*__restrict__ coop, 
              double tol, int iedge, 
              double*__restrict__ param, double*__restrict__ bary,
              double*__restrict__ coopr);

template <int gdim>
int projptedg(MeshBase &msh, const double*__restrict__ coop, 
              double tol, int iref,
              const int*__restrict__ edg2pol, 
              double*__restrict__ param, double*__restrict__ bary,
              double*__restrict__ coopr);


}// namespace Metris

#endif