//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LOW_MEAS_HXX__
#define __METRIS_LOW_MEAS_HXX__

#include "../types_arrays.hxx"
#include "../Mesh/MeshFwd.hxx"


namespace Metris{
struct MetrisParameters;

// Isotropic length, area or volume. 
template <int gdim>
double getmeasentP1(const int *ent2pol,const dblAr2& coord);

// Anisotropic
template <class MFT, int gdim, int ideg>
double getmeasent_aniso(const MeshMetric<MFT> &msh, int ientt);

// Checks whether element valid under the tolerances, also uses normal deviation checks for surface elements.
// - inp norref (can be NULL) is a reference normal, e.g. from CAD
// - output meas (can be NULL) is element isotropic measure
// - inpute nordev_tol is maximum normal deviation, -1 for no check
template <int gdim, int tdim>
bool isvalideltP1(const MeshBase&__restrict__ msh, int ientt, const double*__restrict__ norref = NULL, 
                  double*__restrict__ meas = NULL, double nordev_tol = -1);
template <int gdim, int tdim>
bool isvalideltP1(const MeshBase&__restrict__ msh, int ientt, const int*__restrict__ nod2bpo, 
                  const double*__restrict__ norref, double*__restrict__ meas, double nordev_tol);
template <int gdim, int tdim>
bool isvalideltP1(const MeshBase&__restrict__ msh, const int*__restrict__ ent2pol, const int*__restrict__ nod2bpo, 
                  const double*__restrict__ norref, double*__restrict__ meas, double nordev_tol);



// --- Others
// Equivalent to isvalideltP1 with non-NULL meas. In practice, we very rarely care about the volume, only the validity.
template <int gdim, int tdim>
double getmeasentP1(const MeshBase&__restrict__ msh, int ientt, 
                    const double*__restrict__ norref, bool*__restrict__ iflat, double nordev_tol = -1);

// This is only for unfinished elements:
template<int gdim, int tdim>
double getmeasentP1(const MeshBase&__restrict__ msh, const int*__restrict__ ent2pol, 
                    const int*__restrict__ nod2bpo,
                    const double*__restrict__ norref, bool*__restrict__ iflat, double nordev_tol);



template <int gdim>
void getmeasentP1grad(const int *ent2pol, const dblAr2& coord, int idof, double *grad);

}// namespace Metris
#endif