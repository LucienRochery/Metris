//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_BEZ2GAPS__
#define __METRIS_MSH_BEZ2GAPS__


#include "Mesh/MeshFwd.hxx"


namespace Metris{

/*
Prefer calling one of these two routines
*/
template <int ideg>
void bez2gaps(MeshBase &msh);

template <int ideg>
void gaps2bez(MeshBase &msh);



/*
These can be called to translate a single element. For debugging purposes only (and aux to above). 
*/
template <int ndimn, int ideg>
void eltbez2gaps(MeshBase &msh, int ielem);

template <int ndimn, int ideg>
void eltgaps2bez(MeshBase &msh, int ielem);


/*
No reason to ever call these directly. 
*/
/*
-1 corresponds to bez -> gaps
+1                gaps-> bez
Can be used to transform a single tetra, by first incrementing tag[0]. 
*/
template <int ndimn, int ideg>
void eltgapsbezconv(MeshBase &msh, int ielem, int isig, bool confirm = false);



} // End namespace

#endif