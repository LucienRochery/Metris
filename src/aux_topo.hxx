//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __SRC_AUX_TOPO__
#define __SRC_AUX_TOPO__

#include "types.hxx"



namespace Metris{

class MeshBase;
	

// Get a tetrahedron point belongs to. 
int getpoitet(const MeshBase &msh, int ipoin, int iprt=0);
// Get a triangle point belongs to even if tet attached. 
int getpoifac(const MeshBase &msh, int ipoin);
// Get an edge point belongs to even if tet or tri attached. 
int getpoiedg(const MeshBase &msh, int ipoin);
// Get entity of topo dimn tdimn 
int getpoient(const MeshBase &msh, int ipoin, int tdimn);


// Get global edges or faces (hash tables)
// Note: morally const, but hashtab internals 
// Return -1 on not found
int getedgglo(const MeshBase &msh, int i1, int i2);
int getfacglo(const MeshBase &msh, int i1, int i2, int i3);


// Get various sub-entities
int getedgfacOpp(const MeshBase &msh, int iface, int i1, int i2);
int getedgfac(const MeshBase &msh, int iface, int i1, int i2);
int getfactetOpp(const MeshBase &msh, int ielem, int i1, int i2, int i3);
int getfactet(const MeshBase &msh, int ielem, int i1, int i2, int i3);
int getneitet(const MeshBase &msh, int iele1, int iele2);
int getedgtet(const MeshBase &msh, int ielem, int i1, int i2);

// In case this is changed in the future; e.g. -1 but what if unsigned? 
inline bool isdeadent(int ient, const intAr2 &lent){
	return (lent[ient][0] == lent[ient][1] || lent[ient][0] < 0);
}
inline void killent(int ient, intAr2 &lent){
	lent[ient][0] = -1;
	//lent[ient][0] = lent[ient][1];
}


// Is i1, i2 the ied-th edge of ielem?
int isedgtet(const intAr2 &tet2poi,int ielem, int ied,int i1,int i2);



// Debug prints
void print_bpolist(MeshBase &msh, int ibpoi);


// gather the ibpois of the points on ientt of topological 
// dim tdim == 1 (edge) or == 2 (triangle)
template<int ideg,int tdim>
void getbpois(const MeshBase &msh, int ientt, int *lbpoi);


/*
Cycle non manifold triangles around edge.
iface is only used to return true/false for convenience. 
ifac2 is the one used. 

 Example:
	ifac2 = iface;
	i1, i2 (or ied), iface
	ied   = getedgfac(msh,iface,i1,i2)
  while(getnextfacnm(msh,iface,&ied,&ifac2)){
		// All other faces than iface in NM triangle shell
  }
*/
bool getnextfacnm(const MeshBase &msh, int iface, int i1, int i2,
	                int* ifac2, int* ied);
bool getnextedgnm(const MeshBase &msh, int iedge, int ipoin, 
                        int* iedg2, int* inei);




/*
Helpers: copy. Vertices copied or not, no guarantee. Mesh should be conforming in the first place.
*/
// Copy global edge iedge onto iedfa-th edge of iface
template<int ideg> 
void cpy_gloedg2facedg(MeshBase &msh, int iedge, int iface, int iedfa);
// Copy edge from ifac1 to ifac2
template<int ideg> 
void cpy_facedg2facedg(MeshBase &msh, int ifac1, int iedf1, int ifac2, int iedf2);

// Copy global edge edge onto iedel-th edge of ielem
template<int ideg> 
void cpy_gloedg2tetedg(MeshBase &msh, int iedge, int ielem, int iedel);
// Copy edge from iele1 to iele2
template<int ideg> 
void cpy_tetedg2tetedg(MeshBase &msh, int iele1, int iede1, int iele2, int iede2);
// Copy edge from iface to ielem
template<int ideg> 
void cpy_facedg2tetedg(MeshBase &msh, int iface, int iedfa, int ielem, int iedel);

// Copy global face iface onto ifael-th face of ielem
template<int ideg> 
void cpy_glofac2tetfac(MeshBase &msh, int iface, int ielem, int ifael);
// Copy face from iele1 to iele2
template<int ideg> 
void cpy_tetfac2tetfac(MeshBase &msh, int iele1, int ifae1, int iele2, int ifae2);





} // End namespace

#endif