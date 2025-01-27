//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __LOW_LOCALIZATION__
#define __LOW_LOCALIZATION__

#include "../types.hxx"

namespace Metris{

class MeshBase;

/* Call these ones in priority */
// tol: absolute tolerance on 2-norm between coor0 and coopr
// Return value: 
//  0: in element
//  1: converged to outside
template<int gdim, int ideg> 
int inveval_badNewton(MeshBase &msh, int ientt, const double*__restrict__ coor0, 
                      double*__restrict__ coopr, double*__restrict__ bary,
                      double tol);

// tol: absolute tolerance on 2-norm between coor0 and coopr
// Return value: 
//  0: in element
//  1: converged to outside
template<int gdim, int ideg> 
int inveval(MeshBase &msh, int ientt, const double*__restrict__ coor0, 
            double*__restrict__ coopr, double*__restrict__ bary,
            double tol);

// tol: absolute tolerance on 2-norm between coor0 and coopr
// Return value: 
//  0: in element
//  1: converged to outside
//  2: error or unconverged 
template<int gdim, int ideg> 
int inveval_badNewton0(const MeshBase &msh,
                       const int* ent2pol,
                       const dblAr2 &coord,
                       const double* coor0, 
                       double* __restrict__ coopr, double* __restrict__ bary,
                       double tol);

// tol: absolute tolerance on 2-norm between coor0 and coopr
// Return value: 
//  0: in element
//  1: converged to outside
template<int gdim, int ideg> 
int inveval0(const MeshBase &msh,
             const int* ent2pol,
             const dblAr2 &coord,
             const double* coor0, 
             double* __restrict__ coopr, 
             double* __restrict__ bary,
             dblAr1 &work,
             double tol);



// Quick reject using bounding box to avoid exceptions in optim where points 
// are liable to end up outside. 
template<int gdim>
int locMeshQuick(MeshBase &msh, const double *coor0);


// Passing in coord lets us apply preconditionning
struct invevalfun_data{
  const MeshBase &msh;
  const int* ent2pol;
  const dblAr2& coord;
  const double *coor0;
  double *coopr;
  invevalfun_data(const MeshBase &msh_, const int* ent2pol_, 
                  const dblAr2 &coord_,
                  const double* coor0_, double* coopr_)
  : msh(msh_), ent2pol(ent2pol_), coord(coord_), coor0(coor0_), coopr(coopr_){}
};


// For debugging the derivatives and nlopt
template<int gdim, int ideg> 
double invevalfun(const MeshBase &msh,
                  const int* ent2pol,
                  const dblAr2 &coord, // this can differ from msh if precond
                  const double* coor0, 
                  const double* bary ,
                  int ihess,
                  double *coopr,
                  double *grad,
                  double *hess);

template<int gdim,int ideg>
double invevalfun_nlointf(unsigned int nvar, const double *x, 
                          double *grad, void *f_data);


} // End namespace


#endif