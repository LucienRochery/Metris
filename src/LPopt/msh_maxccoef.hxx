//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_MAXCCOEF__
#define __METRIS_MSH_MAXCCOEF__

#include "LPsolver.hxx"
#include "../Mesh/MeshFwd.hxx"


namespace Metris{


enum class OptDoF{Full, HO};

/* Jacobian Correcting Linear Program - global coordinate seperation */
template <int gdim, int tdim, int ideg>
double maximizeCcoef(MeshBase &msh, OptDoF idofs, LPMethod method, LPLib lib);

/* Incorporates the MF into LP */
template <class MFT, int gdim, int tdim, int ideg>
double maximizeMetCcoef(Mesh<MFT> &msh, OptDoF idofs, LPMethod method, LPLib lib, 
                        const bool ccoefPosFlag = true);


/* Attempt to enforce displacements while keeping validity */
template<int gdim, int tdim, int ideg>
double bezGapsLP(MeshBase &msh, const intAr1 &idx_point, const dblAr2 &pos_ctrlp,
                 LPMethod method,  LPLib lib);


/* Helper functions */ 
template<int gdim, int ideg>
double getminccoef(MeshBase &msh);


}// end namespace 



#endif