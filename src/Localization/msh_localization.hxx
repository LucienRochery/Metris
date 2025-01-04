//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_LOCALIZATION__
#define __METRIS_MSH_LOCALIZATION__


#include "../Mesh/MeshFwd.hxx"


namespace Metris{

enum Localization_Errors {
  LOC_ERR_NOERR = 0,
  LOC_ERR_OUTBB = 1,
  LOC_ERR_ALLPOS = 2,
  LOC_ERR_FAILP1 = 3,
  LOC_ERR_PROJ = 4, // Failed the projection 
  LOC_WARN_PROJ = 5 // Warning: a projection has occurred. 
};


template<class MetricFieldType, int bdeg>
void interpFrontBack(Mesh<MetricFieldType> &msh, MeshBack &bak, int ipoi0 = 0);

/*
Localize point in mesh 
- *ielem (in/out): first guess and final in which loc
- coop (in): point to loc
- pdim: topo dim of point (> 0)
- uvsrf: what's found in bpo2rbi: t (if pdim = 1) or (u,v) (if pdim = 2). Can be 
  NULL if pdim == msh.get_tdim(). Otherwise, mandatory.
- iref: CAD ref line (if pdim = 1) or face (if pdim = 2)
- algnd: surf tangent (if pdim = 1) or normal (if pdim = 2) guess. Only necessary
  if pdim > tdim.
- coopr (out): final projection
- bary (out): bary in *ielem
- tol (in): HO only
- iexpensive: stack all prospective, very expensive, never used 
Note:
coop must remain loose as it may not pertain to a vertex (even virtual) in the 
mesh to localize in. Typically front -> back.
*/
template<int gdim, int tdim, int ideg>
int locMesh(MeshBase &msh, int *ielem, 
            const double* coop, int pdim, const double* uvsrf, 
            int iref, const double* algnd_, 
	          double* coopr, double* bary,
	          double tol = 1.0e-6, int ithrd = 0, bool iexpensive = false);



} // End namespace

#endif
