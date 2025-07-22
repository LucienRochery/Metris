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
template <int gdim, int tdim>
bool isvalidelt(const MeshBase& msh, int ientt, const double* norref, double *meas);



// --- Others
// Equivalent to isvalidelt with non-NULL meas. In practice, we very rarely care about the volume, only the validity.
template <int gdim, int tdim>
double getmeasentP1(const MeshBase &msh, const int* ent2pol, 
                    const double *norref, bool *iflat);

template <int gdim, int tdim>
double getmeasentP1(const MetrisParameters *param, 
                    const int* ent2pol, 
                    const dblAr2 &coord,
                    const double* norref, bool* iflat);



template <int gdim>
void getmeasentP1grad(const int *ent2pol, const dblAr2& coord, int idof, double *grad);

}// namespace Metris
#endif