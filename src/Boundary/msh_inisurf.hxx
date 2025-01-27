//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __SRC_INI_SURF__
#define __SRC_INI_SURF__


#include "../types.hxx"


namespace Metris{

class MeshBase;
 

/*
Points from rank 0 to nbpo0 excluded have been read from file. Their (u,v)s are set. 
Those are only verified. 
Points from nbpo0 included to nbpoi excluded have been re-created. Their (u,v)s will be projected. 
typ_proj = 1: project
           0: evaluate at uv guesses
*/
void prjMeshPoints(MeshBase &msh, int nbpo0, bool onlyproj=false,  bool updtX=false);

/*
	Helper functions for iniMeshNeighbours. The main routine already initializes certain surface entities
in the main loops, such as between different reference (interior) elements. 
*/

template <int ideg> void iniMeshBdryTriangles(MeshBase &msh, HshTabInt3 &intfHshTab);
template <int ideg> void iniMeshBdryEdges(MeshBase &msh);
void iniMeshBdryCorners(MeshBase &msh);
template <int ideg>
int iniMeshBdryPoints(MeshBase &msh);

///*
//Generate corner information: 
//	- lpoin size npoin, lpoin[i] > 0 iff corner
//	- lpoin[i] <= ncorn number of corners
//	- 
//*/
//int genCornerList(MeshBase &msh, int mcorn, int *lpoin, int *lbpoi, int offs = 0);
//// For compatibility with Vizir
//void genCornerIdx(MeshBase &msh, int *lbpoi);

/*
Generate lists for VerticesOnGeometric(Vertices|Edges|Faces)
incre = increment ( +1 for file IO )
*/
void genOnGeometricEntLists(MeshBase &msh, intAr1& lcorn, intAr1& lpoic,
                                           intAr2& lgpoe, dblAr2& rgpoe,
                                           intAr2& lgpof, dblAr2& rgpof,
                            int incre = 0);

int getNumCorners(MeshBase &msh); 


} // End namespace

#endif
