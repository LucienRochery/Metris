//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "aux_topo.hxx"

#include "aux_utils.hxx"
#include "ho_constants.hxx"
#include "aux_exceptions.hxx"
#include "metris_constants.hxx"
#include "Mesh/MeshBase.hxx"

#include "types.hxx"
#include "mprintf.hxx"

#include <boost/preprocessor/iteration/local.hpp>

#include <tuple>
#include <assert.h>
#include <stdio.h>


namespace Metris{
	

int getpoitet(const MeshBase &msh, int ipoin, int iprt){
	//printf("Debug ipoin = %d \n",ipoin);
	assert(ipoin >= 0 && ipoin < msh.npoin);

  if(msh.nelem <= 0) return -1;

  if(msh.poi2ent(ipoin,1) == 3) return msh.poi2ent(ipoin,0); 
  
  if(iprt > 0) printf("Start getpoitet ipoin = %d poi2bpo %d \n",ipoin,msh.poi2bpo[ipoin]);
  
	// Note: poi2bpo is not unique for ipoin. 
	int ibpoi = msh.poi2bpo[ipoin];
  METRIS_ASSERT(ibpoi >= 0); 
bdry:
	int itype = msh.bpo2ibi(ibpoi,1);
  if(iprt > 0) printf("(re)start ibpoi itype %d %d \n",ibpoi,itype);

	if(itype == 2){
		int iface = msh.bpo2ibi(ibpoi,2);
    if(iprt > 0) printf("Type 2 iface = %d \n",iface);
		METRIS_ASSERT_MSG(iface >= 0 && iface < msh.nface,
			"face out of bounds nface = "<<msh.nface<<" iface = "<<iface<<" due to ipoin "<<ipoin);
		return msh.fac2tet(iface,0);
	}
	if(itype == 1){
		int iedge = msh.bpo2ibi(ibpoi,2);
    if(iprt > 0) printf("Type 1 iedge = %d \n",iedge);
		METRIS_ASSERT_MSG(iedge >= 0 && iedge < msh.nedge,
			"Edge out of bounds nedge = "<<msh.nedge<<" iedge = "<<iedge<<" due to ipoin "<<ipoin);
		int iface = msh.edg2fac[iedge];
		// There is no triangle, thus no tetrahedron, attached to this edge
		if(iface < 0) return -1;
    METRIS_ASSERT(iface >= 0 && iface < msh.nface);
		return msh.fac2tet(iface,0);
	}

	// The point is a corner. We'll have to go through bpo2ibi to find
	// a boundary entity. Then, we'll proceed as previously. 
	// We only need to jump once as there can only be one corner
	// attached to the point

	ibpoi = msh.bpo2ibi(ibpoi,3);
  if(iprt > 0) printf("Type 0 next = %d ?\n",ibpoi);
//	printf("Debug link to next -> %d \n",ibpoi);
	// There is no entity. 
	if(ibpoi == -1) return -1;
goto bdry; 

	return -1;
}


int getpoifac(const MeshBase &msh, int ipoin){
	assert(ipoin >= 0 && ipoin < msh.npoin);
  if(msh.nface <= 0) return -1;

  if(msh.poi2ent(ipoin,1) == 2) return msh.poi2ent(ipoin,0);

  // If the point is not tagged as boundary, it cannot then be attached to a triangle.
  if(msh.isboundary_faces()) return -1;

	int ibpoi = msh.poi2bpo[ipoin];
  METRIS_ASSERT(ibpoi >= 0);

nocor:
	int itype = msh.bpo2ibi(ibpoi,1);

	// There can be several ibpois per ipoin, linked by bpo2ibi(ibpoi,3).
	// The first entry (stored in poi2bpo) is always a lowest-dim :

	// 2. Only triangles attached, and only one at that. 
	if(itype == 2 && msh.idim == 2) METRIS_THROW_MSG(TopoExcept(), "Triangle flagged as bdry element in 2D");
	if(itype == 2) return msh.bpo2ibi(ibpoi,2);

	// 1. Edge case is simple because we can use edg2fac:
	if(itype == 1){
		int iedge = msh.bpo2ibi(ibpoi,2);
		assert(iedge >= 0 && iedge < msh.nedge);
		return msh.edg2fac[iedge];
	}

	// 0. Finally, corner case. We don't store anything related to corners directly. 
	// No need to cycle through all of bpo2ibi, it suffice to get one higher-dim
	// (cannot all be corners) and jump back. 
	
	ibpoi = msh.bpo2ibi(ibpoi,3);
	if(ibpoi == -1) return -1; // Case where detached corner (shouldn't exist though)
goto nocor;

	return -1;
}


int getpoiedg(const MeshBase &msh, int ipoin){
	assert(ipoin >= 0 && ipoin < msh.npoin);

  if(msh.poi2ent(ipoin,1) == 1) return msh.poi2ent(ipoin,0); 

  if(msh.isboundary_edges()) return -1;

  int ibpoi = msh.poi2bpo[ipoin];
  METRIS_ASSERT(ibpoi >= 0);

nocor:
	int itype = msh.bpo2ibi(ibpoi,1);

	// There can be several ibpois per ipoin, linked by bpo2ibi(ibpoi,3).
	// The first entry (stored in poi2bpo) is always a lowest-dim :

	// 2: Lowest dim is a triangle, this is surface interior, no edges. 
	if(itype == 2) return -1;

	// 1. Edge 
	if(itype == 1){
		int iedge = msh.bpo2ibi(ibpoi,2);
		assert(iedge >= 0 && iedge < msh.nedge);
		return iedge;
	}

	// 0. Finally, corner case. We don't store anything related to corners directly. 
	// No need to cycle through all of bpo2ibi, it suffice to get one higher-dim
	// (cannot all be corners) and jump back. 
	ibpoi = msh.bpo2ibi(ibpoi,3);
	if(ibpoi == -1) return -1; // Case where detached corner (shouldn't exist though)
goto nocor;

	return -1;
}

int getpoient(const MeshBase &msh, int ipoin, int tdimn){
  if(tdimn == 1)     return getpoiedg(msh,ipoin);
  else if(tdimn == 2)return getpoifac(msh,ipoin); 
  else if(tdimn == 3)return getpoitet(msh,ipoin); 
  METRIS_THROW_MSG(WArgExcept(),"tdimn out of bounds");
}



int getedgglo(const MeshBase &msh, int i1, int i2){
  auto key = stup2(i1,i2);
  auto t = msh.edgHshTab.find(key);
  if(t == msh.edgHshTab.end()) return -1;
  return t->second;
}

int getfacglo(const MeshBase &msh, int i1, int i2, int i3){
  auto key = stup3(i1,i2,i3);
  auto t = msh.facHshTab.find(key);
  if(t == msh.facHshTab.end()) return -1;
  return t->second;
}



int getedgfac(const MeshBase &msh, int iface, int i1, int i2){
  METRIS_ASSERT(iface >= 0 && iface < msh.nface);
	for(int i=0;i<3;i++){
		int j1 = msh.fac2poi(iface,lnoed2[i][0]);
		int j2 = msh.fac2poi(iface,lnoed2[i][1]);
		if((i1 == j1 && i2 == j2) || (i1 == j2 && i2 == j1)) return i;
	}
	METRIS_THROW_MSG(TopoExcept(),
    "EDGE NOT IN TRIANGLE iface = "<<iface<<" vertices "<<
    msh.fac2poi(iface,0)<<" "<<msh.fac2poi(iface,1)<<" "<<
    msh.fac2poi(iface,2)<<" "<<
    " i1 = "<<i1<<" i2 = "<<i2);
}

int getedgfacOpp(const MeshBase &msh, int iface, int i1, int i2){
  METRIS_ASSERT(iface >= 0 && iface < msh.nface);
	for(int i=0;i<3;i++){
		int j1 = msh.fac2poi(iface,lnoed2[i][0]);
		int j2 = msh.fac2poi(iface,lnoed2[i][1]);
		if(i1 == j2 && i2 == j1) return i;
	}
	METRIS_THROW_MSG(TopoExcept(),"EDGE NOT IN TRIANGLE");
}

int getfactetOpp(const MeshBase &msh, int ielem, int i1, int i2, int i3){
  METRIS_ASSERT(ielem >= 0 && ielem < msh.nelem);
	for(int i=0;i<4;i++){
		int j1 = msh.tet2poi(ielem,lnofa3[i][0]);
		int j2 = msh.tet2poi(ielem,lnofa3[i][1]);
		int j3 = msh.tet2poi(ielem,lnofa3[i][2]);
		if((i1 == j1 && i2 == j3 && i3 == j2) 
		|| (i1 == j2 && i2 == j1 && i3 == j3) 
		|| (i1 == j3 && i2 == j2 && i3 == j1)) return i;
	}
	METRIS_THROW_MSG(TopoExcept(),"FACE NOT BACKWARDS (OR NOT AT ALL) IN TET");
}

int getfactet(const MeshBase &msh, int ielem, int i1, int i2, int i3){
  METRIS_ASSERT(ielem >= 0 && ielem < msh.nelem);
	for(int i=0;i<4;i++){
		int j1 = msh.tet2poi(ielem,lnofa3[i][0]);
		int j2 = msh.tet2poi(ielem,lnofa3[i][1]);
		int j3 = msh.tet2poi(ielem,lnofa3[i][2]);
		if((i1 == j1 && i2 == j3 && i3 == j2) 
		|| (i1 == j1 && i2 == j2 && i3 == j3)
		|| (i1 == j2 && i2 == j1 && i3 == j3) 
		|| (i1 == j2 && i2 == j3 && i3 == j1) 
		|| (i1 == j3 && i2 == j2 && i3 == j1)
		|| (i1 == j3 && i2 == j1 && i3 == j2) ) return i;
	}
	METRIS_THROW_MSG(TopoExcept(),"FACE NOT IN TET");
}

// Get the face ifa2 such that iele2's ifa-th neighbour is iele1
int getneitet(const MeshBase &msh, int iele1, int iele2){
	for(int i = 0; i < 4; i++){
		if(msh.tet2tet(iele2,i) == iele1) return i;
	}
	METRIS_THROW_MSG(TopoExcept(),"TET NEIGHBOUR DOES NOT EXIST");
}

int getedgtet(const MeshBase &msh, int ielem,  int i1, int i2){
	for(int i=0;i<6;i++){
		int j1 = msh.tet2poi(ielem,lnoed3[i][0]);
		int j2 = msh.tet2poi(ielem,lnoed3[i][1]);
		if((i1 == j1 && i2 == j2) || (i1 == j2 && i2 == j1)) return i;
	}
	METRIS_THROW_MSG(TopoExcept(),"EDGE NOT IN TET");
}




int isedgtet(const intAr2 &tet2poi,int ielem, int ied,int i1,int i2){
  int j1 = tet2poi(ielem,lnoed3[ied][0]);
  int j2 = tet2poi(ielem,lnoed3[ied][1]);

  return ( (j1 == i1 && j2 == i2)
			   ||(j1 == i2 && j2 == i1));
}


void print_bpolist(MeshBase &msh, int ibpoi){
	if(ibpoi < 0) return;

  GETVDEPTH(msh);
	MPRINTF("-- START printing list for ibpoi = %d \n",ibpoi);

	int ibpo2 = ibpoi;
	int nlist = 0;
	do{
		nlist++;
		if(nlist > 100){
			printf("## LIST > 100\n");
			return;
		}
		MPRINTF("   - %d: %d = (%d %d %d %d)\n",nlist,ibpo2,msh.bpo2ibi(ibpo2,0)
			      ,msh.bpo2ibi(ibpo2,1),msh.bpo2ibi(ibpo2,2),msh.bpo2ibi(ibpo2,3));
		ibpo2 = msh.bpo2ibi(ibpo2,3);
	}while(ibpo2 >= 0 && ibpo2 != ibpoi);
}



// gather the ibpois of the points on ientt of topological dim tdim == 1 (edge) or == 2 (triangle)
template<int ideg,int tdim>
void getbpois(const MeshBase &msh, int ientt, int *lbpoi){

	if constexpr(tdim == 1){
		for(int irnk = 0; irnk < edgnpps[ideg]; irnk++){
			int ipoin = msh.edg2poi(ientt,irnk);
			assert(ipoin >= 0 && ipoin < msh.npoin);
			int ibpo0 = msh.poi2bpo[ipoin];
			assert(ibpo0 >= 0 && ibpo0 < msh.nbpoi);
			// If the first pointer is already of type 1, then there is no lower-dim entity
			// In that case, there is no guarantee there will be an entry for ientt specifically
			// as there is no ambiguity. 
			if(msh.bpo2ibi(ibpo0,1) == 1){
				lbpoi[irnk] = ibpo0;
				continue;
			}
			int ibpoi = ibpo0;
			int nloop = 0;
			while((msh.bpo2ibi(ibpoi,2) != ientt || msh.bpo2ibi(ibpoi,1) != 1) && msh.bpo2ibi(ibpoi,3) != ibpo0){
				ibpoi = msh.bpo2ibi(ibpoi,3);
				nloop++;
				if(nloop > 100) METRIS_THROW_MSG(TopoExcept(),
					"100 BOUNDARY POINTS FOR ONE POINT? INFINITE LOOP");
			}
			assert(msh.bpo2ibi(ibpoi,2) == ientt);
			lbpoi[irnk] = ibpoi;
		}
	}else{
		for(int irnk = 0; irnk < facnpps[ideg]; irnk++){
			int ipoin = msh.fac2poi(ientt,irnk);
			assert(ipoin >= 0 && ipoin < msh.npoin);
			int ibpo0 = msh.poi2bpo[ipoin];
			assert(ibpo0 >= 0 && ibpo0 < msh.nbpoi);
			// If the first pointer is already of type 2, then there is no lower-dim entity
			// In that case, there is no guarantee there will be an entry for ientt specifically
			// as there is no ambiguity. 
			if(msh.bpo2ibi(ibpo0,1) == 2){
				lbpoi[irnk] = ibpo0;
				continue;
			}
			int ibpoi = ibpo0;
			int nloop = 0;
			while((msh.bpo2ibi(ibpoi,2) != ientt || msh.bpo2ibi(ibpoi,1) != 2) 
				  && msh.bpo2ibi(ibpoi,3) != ibpo0){
				ibpoi = msh.bpo2ibi(ibpoi,3);
				assert(ibpoi < msh.nbpoi);
				nloop++;
				if(nloop > 100)METRIS_THROW_MSG(TopoExcept(),
					"100 BOUNDARY POINTS FOR ONE POINT? INFINITE LOOP");
			}
			assert(msh.bpo2ibi(ibpoi,2) == ientt);
			lbpoi[irnk] = ibpoi;
		}
	}
}
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template void getbpois< n ,1>(const MeshBase &msh, int ientt, int *lbpoi);\
template void getbpois< n ,2>(const MeshBase &msh, int ientt, int *lbpoi);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()




bool getnextfacnm(const MeshBase &msh, int iface, int i1, int i2,
	                      int* ifac2, int* ied){
  *ifac2 = - msh.fac2fac(*ifac2,*ied) - 2;
  METRIS_ASSERT_MSG(*ifac2 >= 0 && *ifac2 < msh.nface, "Invalid next face ifac2 = "<<*ifac2);
  *ied  = getedgfac(msh,*ifac2,i1,i2);
	return *ifac2 != iface;
}


bool getnextedgnm(const MeshBase &msh, int iedg0, int ipoin, 
                        int* iedg2, int* inei){
  *iedg2 = - msh.edg2edg(*iedg2,*inei) - 2;
  METRIS_ASSERT_MSG(*iedg2 >= 0 && *iedg2 < msh.nedge, "Invalid next nm edge = "<<*iedg2);
  *inei = -1;
  for(int ii = 0; ii < 2; ii++){
    if(msh.edg2poi(*iedg2,1-ii) == ipoin) *inei = ii;
  }
  return *iedg2 != iedg0;
}



//template<>
//void newbpotopo<0>(MeshBase &msh, int ipoin, int ientt, bool lnkent){
//	// Create new corner if does not exist
//	// Do not update edg2bpo yet because we only have partial information
//	// Namely, some corners (boundary) will be created after this loop. 
//	int ibpoi = msh.poi2bpo[ipoin]; 
//	if(ibpoi < 0) ibpoi = -1;
//
//	if(msh.nbpoi >= msh.mbpoi)METRIS_THROW_MSG(DMemExcept(),
//		"INCREASE MBPO (iniMeshNeighbours)");
//
//	// Corner always the lowest-dim
//	msh.poi2bpo[ipoin]        = msh.nbpoi;
//	msh.bpo2ibi[msh.nbpoi][0] = ipoin; 
//	msh.bpo2ibi[msh.nbpoi][1] = 0;  // Corner type
//	msh.bpo2ibi[msh.nbpoi][2] = ientt; // Ref
//	msh.bpo2ibi[msh.nbpoi][3] = ibpoi; // Link to next
//
//	msh.poi2ent[ipoin] = ientt;
//
//	if(ibpoi >= 0){
//		// If already a bpoi for this ipoin, insert corner in linked list
//		int tmp = msh.bpo2ibi(ibpoi,3);
//		msh.bpo2ibi(ibpoi,3) = msh.nbpoi;
//		msh.bpo2ibi[msh.nbpoi][3] = tmp;
//		if(tmp < 0) msh.bpo2ibi[msh.nbpoi][3] = ibpoi;
//	}
//
//	msh.nbpoi++;
//	return;
//}
//
//
//template<>
//void newbpotopo<1>(MeshBase &msh, int ipoin, int ientt, bool lnkent){
//
//	int ibpoi = msh.poi2bpo[ipoin]; 
//	if(ibpoi < 0) ibpoi = -1;
//	if(msh.nbpoi >= msh.mbpoi)METRIS_THROW_MSG(DMemExcept(),
//		"INCREASE MBPO (iniMeshNeighbours)");
//
//
//	// Only update if new ibpoi is the lowest-dimensional (or unknown yet)
//	if(ibpoi < 0 || msh.bpo2ibi(ibpoi,1) > 1) msh.poi2bpo[ipoin] = msh.nbpoi;
//	msh.bpo2ibi[msh.nbpoi][0] = ipoin; 
//	msh.bpo2ibi[msh.nbpoi][1] = 1;  // Edge type
//	msh.bpo2ibi[msh.nbpoi][2] = ientt; // Ref
//	msh.bpo2ibi[msh.nbpoi][3] = ibpoi; // Link to next
//
//	msh.bpo2rbi[msh.nbpoi][0] = 0;
//	msh.bpo2rbi[msh.nbpoi][1] = 0;
////	printf("Debug newbpotopo ipoin %d ientt %d ibpoi %d \n",ipoin,ientt,ibpoi);
//
//	if(ibpoi >= 0){
//		// There is no particular order here, just put it at the start
//		int tmp = msh.bpo2ibi(ibpoi,3);
//		msh.bpo2ibi(ibpoi,3) = msh.nbpoi;
//		msh.bpo2ibi[msh.nbpoi][3] = tmp;
//		if(tmp < 0) msh.bpo2ibi[msh.nbpoi][3] = ibpoi;
//	}
//	
////		// Check whether link should be carried out.
////		// Only if this is the lowest dim link (or only link) so far
////		if( lnkent && ientt >= 0 && (ibpoi < 0 || ibpoi >= 0 && msh.bpo2ibi(ibpoi,1) < 1) ){
////			int iver = getveredg<ideg>(ientt,msh.edg2poi,ipoin);
////	#ifndef NDEBUG
////	//		printf("Dbug link new edge ipoin %d ientt %d \n",ipoin,ientt);
////			if(iver < 0){
////				printf("## IVER = %d < 0\n",iver);
////				printf("iedge = %d = (%d,%d) ipoin = %d\n",ientt,msh.edg2poi(ientt,0),msh.edg2poi(ientt,1),ipoin);
////			}
////	#endif
////			assert(iver >= 0);
////	
////			std::tuple<int,int> key = {ientt,iver};
////			auto t = msh.edg2bpo.find(key); 
////			if(t != msh.edg2bpo.end()){
////				printf("## THIS EDGE POINT COMBINATION IS NOT UNIQUE\n");
////				exit(1);
////			}
////			msh.edg2bpo[key] = msh.nbpoi;
////		}
//	
//	msh.nbpoi++;
//}
//
//// Todo: get rid of lnkent.
//template<> void newbpotopo<2>(MeshBase &msh, int ipoin, int ientt, bool lnkent){
//	int ibpoi = msh.poi2bpo[ipoin]; 
//	if(ibpoi < 0) ibpoi = -1;
//	if(msh.nbpoi >= msh.mbpoi)METRIS_THROW_MSG(DMemExcept(),
//		"INCREASE MBPO (iniMeshNeighbours)");
//
//	
//	// Only update if new ibpoi is unknown yet (2 -> always highest dim)
//	if(ibpoi < 0) msh.poi2bpo[ipoin] = msh.nbpoi;
//	msh.bpo2ibi[msh.nbpoi][0] = ipoin; 
//	msh.bpo2ibi[msh.nbpoi][1] = 2;  // Face type
//	msh.bpo2ibi[msh.nbpoi][2] = ientt; // Ref
//	msh.bpo2ibi[msh.nbpoi][3] = ibpoi; // Link to next
//	
//	msh.bpo2rbi[msh.nbpoi][0] = 0;
//	msh.bpo2rbi[msh.nbpoi][1] = 0;
//	
//	if(ibpoi >= 0){
//		// There is no particular order here, just put it at the start
//		int tmp = msh.bpo2ibi(ibpoi,3);
//		msh.bpo2ibi(ibpoi,3) = msh.nbpoi;
//		msh.bpo2ibi[msh.nbpoi][3] = tmp;
//		if(tmp < 0) msh.bpo2ibi[msh.nbpoi][3] = ibpoi;
//	}
//	
//
//// Getting rid of fac2bpo. 
////		// If ityp == 2 (triangle) this ibpoi is NEVER the min type, unless ibpoi < 0
////		// In both cases, no need to walk linked list to find current min type
////		// No need to update poi2bpo either as already done or never to do
////		// Update fac2bpo
////		if( lnkent && ientt >= 0 && (ibpoi < 0 || ibpoi >= 0 && msh.bpo2ibi(ibpoi,1) < 2) ){
////			int iver = getverfac<ideg>(ientt,msh.fac2poi,ipoin);
////			#ifndef NDEBUG
////			if(iver < 0){
////				printf("## FAILED TO FIND VERTEX %d IN FACE %d =",ipoin,ientt);
////				for(int ii = 0; ii < facnpps[ideg]; ii++){printf(" %d ",msh.fac2poi(ientt,ii));}
////				printf("\n");
////				printf("ideg = %d facnpps %d \n",ideg,facnpps[ideg]);
////			}
////			#endif
////			assert(iver >= 0);
////	
////			std::tuple<int,int> key = {ientt,iver};
////			auto t = msh.fac2bpo.find(key); 
////			if(t != msh.fac2bpo.end()){
////				printf("## THIS TRIANGLE POINT COMBINATION IS NOT UNIQUE\n");
////				printf("Inserting ipoin %d initial ibpoi %d ientt = %d \n",
////					ipoin,ibpoi,ientt);
////				printf("ientt iver %d %d \n",ientt,iver);
////				print_bpolist(msh,ibpoi);
////				print_fac2bpo(msh,100);
////				exit(1);
////			}
////			msh.fac2bpo[key] = msh.nbpoi;
////		}
//	msh.nbpoi++;
//	return;
//}








/* 
	NODE COPY HELPER FUNCTIONS 
*/
// Copy global edge iedge onto iedfa-th edge of iface
template<int ideg> 
void cpy_gloedg2facedg(MeshBase &msh, int iedge, int iface, int iedfa){
	if constexpr(ideg == 1) return;

	int ip1 = msh.fac2poi(iface,lnoed2[iedfa][0]);
	int ip2 = msh.fac2poi(iface,lnoed2[iedfa][1]);
	int irn0 = 3 + iedfa*(edgnpps[ideg] - 2);
	assert((ip1 == msh.edg2poi(iedge,0) && ip2 == msh.edg2poi(iedge,1)) 
		||   (ip2 == msh.edg2poi(iedge,0) && ip1 == msh.edg2poi(iedge,1)) );
	if(msh.edg2poi(iedge,0) == ip1){
		for(int ii=0; ii < edgnpps[ideg]-2; ii++){
			msh.fac2poi[iface][irn0 + ii] = msh.edg2poi[iedge][2+ii];
		}
	}else{
		for(int ii=0; ii < edgnpps[ideg]-2; ii++){
			msh.fac2poi[iface][irn0 + ii] = msh.edg2poi[iedge][edgnpps[ideg]-1-ii];
		}
	}
}


// Copy edge from ifac1 to ifac2
template<int ideg> 
void cpy_facedg2facedg(MeshBase &msh, int ifac1, int iedf1, int ifac2, int iedf2){
	if constexpr(ideg == 1) return;
	
	int ip1 = msh.fac2poi(ifac1,lnoed2[iedf1][0]);
	int ip2 = msh.fac2poi(ifac1,lnoed2[iedf1][1]);

	int jp1 = msh.fac2poi(ifac2,lnoed2[iedf2][0]);
	int jp2 = msh.fac2poi(ifac2,lnoed2[iedf2][1]);

	assert((ip1 == jp1 && ip2 == jp2) 
		||   (ip2 == jp1 && ip1 == jp2) );

	int irn2 = 3 + iedf2*(edgnpps[ideg] - 2);
	if(ip1 == jp1){
		int irn1 = 3 + iedf1*(edgnpps[ideg] - 2);
		for(int ii=0; ii < edgnpps[ideg]-2; ii++){
			msh.fac2poi[ifac2][irn2 + ii] = msh.fac2poi[ifac1][irn1 + ii];
		}
	}else{
		int irn1 = 3 + (iedf1+1)*(edgnpps[ideg] - 2) ;
		for(int ii=0; ii < edgnpps[ideg]-2; ii++){
			msh.fac2poi[ifac2][irn2 + ii] = msh.fac2poi[ifac1][irn1 - 1 - ii];
		}
	}
}





// Copy global edge iedge onto iedel-th edge of ielem
template<int ideg> 
void cpy_gloedg2tetedg(MeshBase &msh, int iedge, int ielem, int iedel){
	if constexpr(ideg == 1) return;

	int ip1 = msh.tet2poi(ielem,lnoed3[iedel][0]);
	int ip2 = msh.tet2poi(ielem,lnoed3[iedel][1]);
	int irn0 = 4 + iedel*(edgnpps[ideg] - 2);
	assert((ip1 == msh.edg2poi(iedge,0) && ip2 == msh.edg2poi(iedge,1)) 
		||   (ip2 == msh.edg2poi(iedge,0) && ip1 == msh.edg2poi(iedge,1)) );
	if(msh.edg2poi(iedge,0) == ip1){
		for(int ii=0; ii < edgnpps[ideg]-2; ii++){
			msh.tet2poi[ielem][irn0 + ii] = msh.edg2poi[iedge][2+ii];
		}
	}else{
		for(int ii=0; ii < edgnpps[ideg]-2; ii++){
			msh.tet2poi[ielem][irn0 + ii] = msh.edg2poi[iedge][edgnpps[ideg]-1-ii];
		}
	}
}

// Copy edge from iele2 to iele2
template<int ideg> 
void cpy_facedg2tetedg(MeshBase &msh, int iface, int iedfa, int ielem, int iedel){
	if constexpr(ideg == 1) return;

	int ip1 = msh.fac2poi(iface,lnoed2[iedfa][0]);
	int ip2 = msh.fac2poi(iface,lnoed2[iedfa][1]);

	int jp1 = msh.tet2poi(ielem,lnoed3[iedel][0]);
	int jp2 = msh.tet2poi(ielem,lnoed3[iedel][1]);

	assert((ip1 == jp1 && ip2 == jp2) 
		||   (ip2 == jp1 && ip1 == jp2) );

	int irn2 = 4 + iedel*(edgnpps[ideg] - 2);
	if(ip1 == jp1){
		int irn1 = 3 + iedfa*(edgnpps[ideg] - 2);
		for(int ii=0; ii < edgnpps[ideg]-2; ii++){
			msh.tet2poi[ielem][irn2 + ii] = msh.fac2poi[iface][irn1 + ii];
		}
	}else{
		int irn1 = 3 + (iedfa+1)*(edgnpps[ideg] - 2) ;
		for(int ii=0; ii < edgnpps[ideg]-2; ii++){
			msh.tet2poi[ielem][irn2 + ii] = msh.fac2poi[iface][irn1 - 1 - ii];
		}
	}
}





// Copy edge from iele2 to iele2
template<int ideg> 
void cpy_tetedg2tetedg(MeshBase &msh, int iele1, int iede1, int iele2, int iede2){
	if constexpr(ideg == 1) return;

	int ip1 = msh.tet2poi(iele1,lnoed3[iede1][0]);
	int ip2 = msh.tet2poi(iele1,lnoed3[iede1][1]);

	int jp1 = msh.tet2poi(iele2,lnoed3[iede2][0]);
	int jp2 = msh.tet2poi(iele2,lnoed3[iede2][1]);

	assert((ip1 == jp1 && ip2 == jp2) 
		||   (ip2 == jp1 && ip1 == jp2) );

	int irn2 = 4 + iede2*(edgnpps[ideg] - 2);
	if(ip1 == jp1){
		int irn1 = 4 + iede1*(edgnpps[ideg] - 2);
		for(int ii=0; ii < edgnpps[ideg]-2; ii++){
			msh.tet2poi[iele2][irn2 + ii] = msh.tet2poi[iele1][irn1 + ii];
		}
	}else{
		int irn1 = 4 + (iede1+1)*(edgnpps[ideg] - 2) ;
		for(int ii=0; ii < edgnpps[ideg]-2; ii++){
			msh.tet2poi[iele2][irn2 + ii] = msh.tet2poi[iele1][irn1 - 1 - ii];
		}
	}
}





// Copy global face iface onto ifael-th face of ielem
template<int ideg> 
void cpy_glofac2tetfac(MeshBase &msh, int iface, int ielem, int ifael){
	if constexpr(ideg == 1) return;

	int ip1 = msh.tet2poi(ielem,lnofa3[ifael][0]);
	int ip2 = msh.tet2poi(ielem,lnofa3[ifael][1]);

	int ii = 2;
			 if(ip1 == msh.fac2poi(iface,0)) ii = 0;
	else if(ip1 == msh.fac2poi(iface,1)) ii = 1;

	int perm[4];
	perm[ii] = 0;

	// Now check the order: is the next our second or our 3rd?
	int jj = (ii+1)%3;
	if(msh.fac2poi(iface,jj) == ip2){
		// Same order
		perm[jj] = 1;
		perm[(jj+1)%3] = 2;
	}else{
		// Reverse order
		perm[jj] = 2;
		perm[(jj+1)%3] = 1;
	}

	#ifndef NDEBUG
	if((msh.tet2poi[ielem][lnofa3[ifael][perm[0]]] != msh.fac2poi(iface,0))
	 ||(msh.tet2poi[ielem][lnofa3[ifael][perm[1]]] != msh.fac2poi(iface,1))
	 ||(msh.tet2poi[ielem][lnofa3[ifael][perm[2]]] != msh.fac2poi(iface,2))){
		printf("## FACES OR PERMUTATION DO NOT CORRESPOND\n");
		printf("perm = %d %d %d seed ii = %d \n",perm[0],perm[1],perm[2],ii);
		printf("Glofac vertices (natural order) %d %d %d \n",msh.fac2poi(iface,0),msh.fac2poi(iface,1),msh.fac2poi(iface,2));
		printf("Tetra vertices (nat ord %d %d %d \n"
			,msh.tet2poi(ielem,lnofa3[ifael][0])
			,msh.tet2poi(ielem,lnofa3[ifael][1])
			,msh.tet2poi(ielem,lnofa3[ifael][2]));
		printf("Tetra vertices (per ord %d %d %d \n"
			,msh.tet2poi[ielem][lnofa3[ifael][perm[0]]]
			,msh.tet2poi[ielem][lnofa3[ifael][perm[1]]]
			,msh.tet2poi[ielem][lnofa3[ifael][perm[2]]]);
	}

	assert(msh.tet2poi[ielem][lnofa3[ifael][perm[0]]] == msh.fac2poi(iface,0));
	assert(msh.tet2poi[ielem][lnofa3[ifael][perm[1]]] == msh.fac2poi(iface,1));
	assert(msh.tet2poi[ielem][lnofa3[ifael][perm[2]]] == msh.fac2poi(iface,2));
	#endif

	int idx_tet[4];
  for(int irnk1 = 0;irnk1 < facnpps[ideg]; irnk1++){
   idx_tet[ifael] = 0; // The opposite vertex gets 0
   idx_tet[lnofa3[ifael][perm[0]]] = ordfac.s[ideg][irnk1][0];
   idx_tet[lnofa3[ifael][perm[1]]] = ordfac.s[ideg][irnk1][1];
   idx_tet[lnofa3[ifael][perm[2]]] = ordfac.s[ideg][irnk1][2];
   int irnkt = mul2nod(idx_tet[0],idx_tet[1],idx_tet[2],idx_tet[3]);
   msh.tet2poi(ielem,irnkt) = msh.fac2poi(iface,irnk1);
  }
}

// Copy face from iele1 to iele2
template<int ideg> 
void cpy_tetfac2tetfac(MeshBase &msh, int iele1, int ifae1, int iele2, int ifae2){
	if constexpr(ideg == 1) return;

	assert(ifae1 >= 0 && ifae1 < 4);
	assert(ifae2 >= 0 && ifae2 < 4);

	// What we copy from
	int jp1 = msh.tet2poi(iele1,lnofa3[ifae1][0]);
	int jp2 = msh.tet2poi(iele1,lnofa3[ifae1][1]);

	// What we copy into
	int ip1 = msh.tet2poi(iele2,lnofa3[ifae2][0]);
	int ip2 = msh.tet2poi(iele2,lnofa3[ifae2][1]);

	int ii = 2;
			 if(ip1 == jp1) ii = 0;
	else if(ip1 == jp2) ii = 1;

	int perm[4];
	perm[ii] = 0;

	// Now check the order: is the next our second or our 3rd?
	int jj = (ii+1)%3;
	if(msh.tet2poi(iele1,lnofa3[ifae1][jj]) == ip2){
		// Same order
		perm[jj] = 1;
		perm[(jj+1)%3] = 2;
	}else{
		// Reverse order
		perm[jj] = 2;
		perm[(jj+1)%3] = 1;
	}

	assert(msh.tet2poi[iele2][lnofa3[ifae2][perm[0]]] == msh.tet2poi(iele1,lnofa3[ifae1][0]));
	assert(msh.tet2poi[iele2][lnofa3[ifae2][perm[1]]] == msh.tet2poi(iele1,lnofa3[ifae1][1]));
	assert(msh.tet2poi[iele2][lnofa3[ifae2][perm[2]]] == msh.tet2poi(iele1,lnofa3[ifae1][2]));

	int idx_tet1[4];
	int idx_tet2[4];
  for(int irnk1 = 0;irnk1 < facnpps[ideg]; irnk1++){
   idx_tet1[ifae1] = 0; // The opposite vertex gets 0
   idx_tet2[ifae2] = 0; // The opposite vertex gets 0
   for(int i = 0; i<3 ;i++){
	   idx_tet1[lnofa3[ifae1][     i ]] = ordfac.s[ideg][irnk1][i];
	   idx_tet2[lnofa3[ifae2][perm[i]]] = ordfac.s[ideg][irnk1][i];
   }
   assert(idx_tet1[0] + idx_tet1[1] + idx_tet1[2] + idx_tet1[3] == ideg);
   assert(idx_tet2[0] + idx_tet2[1] + idx_tet2[2] + idx_tet2[3] == ideg);

   int irnkt1 = mul2nod(idx_tet1[0],idx_tet1[1],idx_tet1[2],idx_tet1[3]);
   int irnkt2 = mul2nod(idx_tet2[0],idx_tet2[1],idx_tet2[2],idx_tet2[3]);
   msh.tet2poi(iele2,irnkt2) = msh.tet2poi(iele1,irnkt1);
  }

}


// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template void cpy_gloedg2facedg< n >(MeshBase &, int, int, int);\
template void cpy_facedg2facedg< n >(MeshBase &, int, int, int, int);\
template void cpy_glofac2tetfac< n >(MeshBase &, int, int, int);\
template void cpy_tetfac2tetfac< n >(MeshBase &, int, int, int, int);\
template void cpy_gloedg2tetedg< n >(MeshBase &, int, int, int);\
template void cpy_facedg2tetedg< n >(MeshBase &, int, int, int, int);\
template void cpy_tetedg2tetedg< n >(MeshBase &, int, int, int, int);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()







} // End namespace

