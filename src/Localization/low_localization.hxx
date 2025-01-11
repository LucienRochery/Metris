//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __LOW_LOCALIZATION__
#define __LOW_LOCALIZATION__

#include "../types.hxx"

namespace Metris{

class MeshBase;

// tol: absolute tolerance on 2-norm between coor0 and coopr
// Return value: 
//  0: in element
//  1: converged to outside
//  2: error or unconverged 
template<int gdim, int ideg> 
int inveval(const MeshBase &msh,
            const int* ent2pol,
            const dblAr2 &coord,
            const double* coor0, 
            double* __restrict__ coopr, double* __restrict__ bary,
            double tol);

template<int gdim, int ideg> 
int inveval(MeshBase &msh, int ientt, const double*__restrict__ coor0, 
            double*__restrict__ coopr, double*__restrict__ bary,
            double tol);


// Quick reject using bounding box to avoid exceptions in optim where points 
// are liable to end up outside. 
template<int gdim>
int locMeshQuick(MeshBase &msh, const double *coor0);


} // End namespace


#endif