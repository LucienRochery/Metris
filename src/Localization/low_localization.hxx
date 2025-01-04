//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __LOW_LOCALIZATION__
#define __LOW_LOCALIZATION__


namespace Metris{

class MeshBase;
//// -----------
//// Invert eval3 within single element 
//template<int ideg, int ilag> int inveval3_old(Mesh &msh, int ielem, const double* coor0, double* __restrict__ coopr, double* __restrict__ bary);
//template<int ideg, int ilag> double fg_inveval3_old(unsigned int n, const double* bar, double* grad, void* f_data);
////struct StructInveval3Funargs{double *rfld;int *lfld;};
//// -----------

int constrain_bary_desc(int gdim, bool *icstr, double* __restrict__ bary,
                        double* __restrict__ bar0, 
                        double* __restrict__ desc, int iprt = 0);

// Tess covering tessellations of ielem for whether the point is in there. 
// Return 1 if point inside tess, 0 otherwise. 
// Put bary guess in bary
template<int gdim, int ideg>
int inveval_tessreject(MeshBase &msh, int ielem, const double* coor0, double*__restrict__ bary, 
                       double tol = 1.0e-6, int iprt = 0);


template<int ideg> 
int inveval3_q2_Newton(MeshBase &msh, int ielem, const double* coor0, 
	            double* __restrict__ coopr, double* __restrict__ bary,
	            double tol, int *ncall, bool ilazy, int iprt = 0);

// ilazy = true: exit as soon as a barycentric is negative
// tol: absolute tolerance on 2-norm between coor0 and coopr
// Return value:
//  0: in tetra
//  1: error (non invertible Jacobian, didn't converge)
// -1: (if lazy) outside of tetra
template<int gdim, int ideg> 
int inveval(MeshBase &msh, int ielem, const double* coor0, 
	          double* __restrict__ coopr, double* __restrict__ bary,
	          double tol = 1.0e-6, bool ilazy = false, int iprt = 0);


template<int gdim, int ideg> 
int inveval_linesearch(MeshBase &msh, int ielem, const double* coor0, 
	          double* __restrict__ coopr, double* __restrict__ bary,
	          double tol, int *ncall, bool ilazy = false, int iprt = 0);



// Quick reject using bounding box to avoid exceptions in optim where points 
// are liable to end up outside. 
template<int gdim>
int locMeshQuick(MeshBase &msh, const double *coor0);


template<int ideg> 
int inveval3_old(MeshBase &msh, int ielem, const double* coor0, 
	           double* __restrict__ coopr, double* __restrict__ bary,
	           double tol = 1.0e-6);

template<int ideg>
double fg_inveval3_old(unsigned int n, const double* bar, double* grad, void* f_data);

} // End namespace


#endif