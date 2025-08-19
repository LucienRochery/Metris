//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __LOW_TOPO__
#define __LOW_TOPO__


#include "Mesh/MeshFwd.hxx"
#include "types.hxx"


namespace Metris{

// Gathers all entities
// lbedg, lbfac and lbtet can be size 0 (allocated to 0)
// in that case, they will not be filled.
int ball(MeshBase& msh, int ipoin,
         intAr1 &lbedg, intAr1 &lbfac, intAr1 &lbtet,
         int *iopen, bool append, int ithrd);


// Get the points surrounding a point, among dimension tdim elements.
// Optionally fills lbent the ball of ipoin. 
int poi2poi(MeshBase& msh, int ipoin, int tdim, intAr1 &lpoin, intAr1 *lbent, int ithrd1);

// Get the points surrounding a point, among dimension tdim elements.
// Caller provides ball lbent.
void poi2poi(MeshBase& msh, int ipoin, int tdim, const intAr1 &lbent, intAr1 &lpoin, int ithrd1);

int ball3(MeshBase& __restrict__ msh,
          int ipoin  ,int iele0, 
          intAr1 &lball,
          int* __restrict__ iopen,
          int ithrd = 0);

void ball3_nm(MeshBase& __restrict__ msh,
              int ipoin  ,int iele0, 
              int* __restrict__ nball_ ,
              int* __restrict__ nbfac_ ,
              intAr1 &lball,
              intAr1 &lbfac, 
              int* __restrict__ iopen,
              int ithread);

void ball3_full(MeshBase& __restrict__ msh,
               int ipoin  ,int tdimn, int iseed, 
               int* __restrict__ nbtet,
               int* __restrict__ nbfac,
               int* __restrict__ nbedg,
               intAr1 &lbtet,
               intAr1 &lbfac, 
               intAr1 &lbedg,
               int ithread);

int ball2(MeshBase& __restrict__ msh,
          int ipoin     ,int ifac0, 
          intAr1&            lbfac,
          intAr1&            lbedg,
          int* __restrict__  iopen,
          bool* __restrict__ imani,
          int ithread);

// iopen = 0 if closed
//       = 1 if open 
void shell(const MeshBase& msh,
           int ipoi1, int ipoi2, 
           int tdim , int iele0, 
           intAr1& lsedg, intAr1& lsfac, intAr1& lstet, int* iopen);

// iopen = 0 if shell is closed
// = boundary shell element otherwise
void shell3(const MeshBase& msh,
	          int ipoi1, int ipoi2, int iele0, 
            intAr1& lshell, intAr1& lbdry,
            int* iopen);


// Gather triangles surrounding non-manifold edge. 
void shell2_nm(const MeshBase& msh,
               int iedge ,
               intAr1&           lshell);

void shell2_nm(const MeshBase& msh,
               int iface ,
               int iedl  ,
               intAr1&           lshell);

} // End namespace

#endif
