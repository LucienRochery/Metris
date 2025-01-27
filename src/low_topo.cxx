//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "low_topo.hxx"

#include "Mesh/Mesh.hxx"

#include "aux_exceptions.hxx"
#include "aux_topo.hxx"

#include <assert.h>
#include <tuple>

namespace Metris{

std::tuple<int,int,int> stup3(int i1,int i2,int i3);


// Basic manifold ball, disregarding internal faces. 
// lball(:,0) stores elements
// lball(:,1) stores rank of point in elements
int ball3(MeshBase& __restrict__ msh,
           int ipoin  ,int iele0, 
           intAr1&           lball,
           int* __restrict__ iopen,
           int ithrd){
	int iball,ielem,i,iele2;

	assert("Input iele0 within bounds " && iele0 >= 0 && iele0 < msh.nelem);
	assert("Sought point is inside element " && (ipoin == msh.tet2poi(iele0,0)
						                                 ||ipoin == msh.tet2poi(iele0,1)
						                                 ||ipoin == msh.tet2poi(iele0,2)
						                                 ||ipoin == msh.tet2poi(iele0,3)));
	assert("Input ielem alive" && !isdeadent(iele0,msh.tet2poi));

	msh.tag[ithrd] += 1;
	*iopen = 0;

	ielem = iele0;

	for(i = 0;i < 4; i++){
		if(msh.tet2poi(iele0,i) == ipoin) break;
	}

  lball.set_n(0); 
  lball.stack(iele0); 
	msh.tet2tag(ithrd,iele0) = msh.tag[ithrd];

	iball = 0;
	while(iball < lball.get_n()){
		ielem = lball[iball];
    int iver = getvertet<1>(ielem, msh.tet2poi, ipoin); 

		assert("Ball element is alive " && !isdeadent(ielem,msh.tet2poi));

//  Loop over neighbours skipping those already in ball
		for(i=0;i<4;i++){
			// This face is opposite the vertex: ball boundary
			if(iver == i) continue;

			iele2 = msh.tet2tet(ielem,i);
			assert("Neighbour table correct " && iele2 >= -1);

//    No neighbour here means the ball is open
//    Flag and move on
			if(iele2 == -1){
				*iopen = 1;
				continue;
			}

//    This element is already in the ball
			if(msh.tet2tag(ithrd,iele2) >= msh.tag[ithrd]) continue;
      msh.tet2tag(ithrd,iele2) = msh.tag[ithrd];

//    Add to stack 
      lball.stack(iele2); 
		}


		iball++;
	}

  return 0;
}


// This version also gathers internal faces to the ball (including boundary if open)
// Does not gather edges. 
void ball3_nm(MeshBase& __restrict__ msh,
           	  int ipoin  ,int iele0, 
           	  int* __restrict__ nball_,
           	  int* __restrict__ nbfac_,
           	  intAr1&           lball,
           	  intAr1&           lbfac, 
           	  int* __restrict__ iopen,
              int ithread){
	int iball,ielem,i,iele2,nball,nbfac;

	assert("Input iele0 within bounds " && iele0 >= 0 && iele0 < msh.nelem);
	assert("Sought point is inside element " && (ipoin == msh.tet2poi(iele0,0)
						                                 ||ipoin == msh.tet2poi(iele0,1)
						                                 ||ipoin == msh.tet2poi(iele0,2)
						                                 ||ipoin == msh.tet2poi(iele0,3)));
	assert("Input ielem alive" && !isdeadent(iele0,msh.tet2poi));

	msh.tag[ithread] += 1;
	*iopen = 0;

	ielem = iele0;

	for(i = 0;i < 4; i++){
		if(msh.tet2poi(iele0,i) == ipoin) break;
	}

	lball[0] = iele0;
	nbfac = 0;
	nball = 1;
	msh.tet2tag(ithread,iele0) = msh.tag[ithread];

	iball = 0;
	while(iball < nball){
		ielem = lball[iball];
    int iver = getvertet<1>(ielem, msh.tet2poi, ipoin); 

		assert("Ball element is alive " && !isdeadent(ielem,msh.tet2poi));

//  Loop over neighbours skipping those already in ball
		for(i=0;i<4;i++){
			// This face is opposite the vertex: ball boundary
			if(iver == i) continue;

			// The array was devised to avoid unnecessary hash table lookups
			if(msh.tet2ftg[ielem]){
				int iface = msh.tetfac2glo(ielem,i);
				if(iface < 0 || iface >= msh.nface) 
					METRIS_THROW_MSG(TopoExcept(),"Face missing or invalid in hash tab "<<iface<<"\n");
      	nbfac++; 
				if(nbfac > lbfac.size1()) METRIS_THROW(SMemExcept());
				lbfac[nbfac-1] = iface; 
			}

			iele2 = msh.tet2tet(ielem,i);
			assert("Neighbour table correct " && iele2 >= -1);

//    No neighbour here means the ball is open
//    Flag and move on
			if(iele2 == -1){
				*iopen = 1;
				continue;
			}

//    This element is already in the ball
			if(msh.tet2tag(ithread,iele2) >= msh.tag[ithread]) continue;

//    Add to stack and get vertex
			nball++;
			if(nball >= lball.size1()) METRIS_THROW(SMemExcept());

			msh.tet2tag(ithread,iele2) = msh.tag[ithread];
			lball[nball-1] = iele2;
		}


		iball++;
	}

	*nball_ = nball;
	*nbfac_ = nbfac;
}



// Gather everything: tetras, triangles, edges
// Seed can be any: edge, triangle or tetra
void ball3_full(MeshBase& __restrict__ msh,
                int ipoin  ,int tdimn, int iseed, 
                int* __restrict__ nbtet,
                int* __restrict__ nbfac,
                int* __restrict__ nbedg,
                intAr1&           lbtet,
                intAr1&           lbfac, 
                intAr1&           lbedg,
                int ithread){

  METRIS_THROW_MSG(TODOExcept(),"All wrong here, reimplement")

  /*
  METRIS_ASSERT(tdimn >= 1 && tdimn <= 3);
  METRIS_ASSERT(lbedg.n > 0 && lbfac.n > 0 && (lbtet.n > 0 || msh.nelem <= 0));

  *nbedg = 0;
  *nbfac = 0;
  *nbtet = 0;


  METRIS_ASSERT(ipoin >= 0 && ipoin < msh.npoin);
  int ibpoi = msh.poi2bpo[ipoin];

  if(ibpoi >= 0){
    int itype = msh.bpo2ibi(ibpoi,1);
    if(itype == 0){ // Corner case: we can potentially get all boundary entities

    }
  }

  int iedg0=-1, ifac0=-1, iele0=-1;
  if(tdimn == 1){
    METRIS_THROW_MSG(TODOExcept(), "edge ball3 not implemented");
    if(msh.isboundary_edges()) METRIS_ASSERT(msh.poi2bpo[ipoin] >= 0);
    iedg0 = iseed;
  }else if(tdimn == 2){
    ifac0 = iseed;
    int ibpoi = msh.poi2bpo[ipoin];
    if(msh.isboundary_faces()) METRIS_ASSERT(ibpoi >= 0);
    // This will probably always hold but just in case, this is an assumption
    METRIS_ASSERT(msh.isboundary_edges());
    // Unroll poi2bpo, we get the edges for free. 
    #ifndef NDEBUG 
      // Just in case, we'll check for duplicates. 
      msh.tag[ithread]++;
    #endif
    int ibpo2 = ibpoi;
    do{
      int itype = msh.bpo2ibi(ibpo2,1); 
      int ientt = msh.bpo2ibi(ibpo2,2);
      METRIS_ASSERT(ientt >= 0);
      if(itype == 1){
        if(*nbedg >= lbedg.n) METRIS_THROW_MSG(DMemExcept(),"Increase lbedg.ne")
        #ifndef NDEBUG 
          METRIS_ASSERT(msh.edg2tag(ithread,ientt) < msh.tag[ithread]);
          msh.edg2tag(ithread,ientt) = msh.tag[ithread];
        #endif
        lbedg[*nbedg] = ientt;
        (*nbedg)++;
      }else if(itype == 2){
        if(*nbfac >= lbfac.n) METRIS_THROW_MSG(DMemExcept(),"Increase lbfac.ne")
        #ifndef NDEBUG 
          METRIS_ASSERT(msh.fac2tag(ithread,ientt) < msh.tag[ithread]);
          msh.fac2tag(ithread,ientt) = msh.tag[ithread];
        #endif
        lbfac[*nbfac] = ientt;
        (*nbfac)++;
      }
      ibpo2 = msh.bpo2ibi(ibpo2,3);
    }while(ibpo2 != ibpoi && ibpo2 >= 0);

    // Call basic ball3 for tetrahedra
    if(msh.nelem >= 0){
      iele0 = getpoitet(msh,ipoin);
      int iopen;
      ball3(msh, ipoin ,iele0, lbtet.n, nbtet, lbtet, &iopen, ithread);
    }

    return;
  }else{
    METRIS_THROW_MSG(TODOExcept(), "eleme ball3 not implemented")
    iele0 = iseed;    
  }
  */
}



// Same as bfac3 but using triangle as seed. Can handle non manifold meshes. 
// imani returns wheter manifold (true) or non-manifold (false)
// 
// lbfac only stores iface, not the index. Just call getverfac if needed. 
// 
// If lbedg.n > 0, gather internal edges + bdry if iopen 
int ball2(MeshBase& __restrict__ msh,
          int ipoin  ,int ifac0, 
          //int* __restrict__ nbfac_,
          intAr1&           lbfac,
          //int* __restrict__ nbedg_,
          intAr1&           lbedg,
          int* __restrict__ iopen,
          bool* __restrict__ imani,
          int ithread){
  int ibfac,iface,ifac2;

  assert("Input ifac0 within bounds " && ifac0 >= 0 && ifac0 < msh.nface);
  assert("Sought point is inside element " && (ipoin == msh.fac2poi(ifac0,0)
                                             ||ipoin == msh.fac2poi(ifac0,1)
                                             ||ipoin == msh.fac2poi(ifac0,2)));
  assert("Input iface alive" && !isdeadent(ifac0,msh.fac2poi));

  *imani = true; 

  lbfac.set_n(0);
  lbedg.set_n(0);

  if(msh.isboundary_faces()){
    // If the point is attached to an edge or corner, bpo2ibi stores 
    // all the triangles it belongs to. We can exploit this to go faster.
    int ibpoi = msh.poi2bpo[ipoin];
    assert(ibpoi >= 0 && ibpoi < msh.nbpoi);
    int itype = msh.bpo2ibi(ibpoi,1);
    //nbfac = 0;
    //nbedg = 0;
    // If the point is only edge, and we want edges, we won't get them this way
    if(itype == 1 && lbedg.size1() == 0 
    || itype == 0){
      int ibpo2 = msh.bpo2ibi(ibpoi,3);
      while(ibpo2 != -1){
        if(msh.bpo2ibi(ibpo2,1) == 2){

          int iface = msh.bpo2ibi(ibpo2,2);
          METRIS_ASSERT(iface >= 0 && iface < msh.nface);

          lbfac.stack(iface); 

        }if(lbedg.size1() > 0 && msh.bpo2ibi(ibpo2,1) == 1){

          int iedge = msh.bpo2ibi(ibpo2,2);
          METRIS_ASSERT(iedge >= 0 && iedge < msh.nedge);

          lbedg.stack(iedge); 

        }
        ibpo2 = msh.bpo2ibi(ibpo2,3);
      }

      //*nbfac_ = nbfac;
      //*nbedg_ = nbedg;
      return 0;
    }
  }
  // If the point is only attached to triangles, we have to do it 
  // the usual way. 
  msh.tag[ithread] += 1;
  *iopen = 0;

  iface = ifac0;

  int i;
  for(i = 0;i < 3; i++){
    if(msh.fac2poi(ifac0,i) == ipoin) break;
  }

  lbfac.stack(ifac0);
  msh.fac2tag(ithread,ifac0) = msh.tag[ithread];

  ibfac = 0;
  while(ibfac < lbfac.get_n()){
    iface = lbfac[ibfac];

    METRIS_ASSERT("bfac element is alive " && !isdeadent(iface,msh.fac2poi));

//  Loop over neighbours skipping those already in bfac
    for(int ied = 0; ied < 3; ied++){

      // Edge opposite to ipoin does not contain ipoin !
      if(msh.fac2poi(iface,ied) == ipoin) continue; 

      ifac2 = msh.fac2fac(iface,ied);


      int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
      int ip2 = msh.fac2poi(iface,lnoed2[ied][1]);

      if(ifac2 < -1){
        *imani = false;
        // Non-manifold case
        // This is where we'd glean the edges if we wanted them. 
        int ied2 = ied;
        int ifac3 = iface;
        while(getnextfacnm(msh,iface,ip1,ip2,&ifac3,&ied2)){
          if(msh.fac2tag(ithread,ifac3) >= msh.tag[ithread]){
            printf("## DEBUG THIS CASE SHOULD NEVER HAPPEN\n");
            printf("Indeed if a nm face has been seen, all should have been.\n");
            exit(1);
          }
          lbfac.stack(ifac3);
        }

        if(lbedg.size1() > 0){
          int iedge = getedgglo(msh,ip1,ip2);
          METRIS_ASSERT(iedge >= 0);
          if(msh.edg2tag(ithread,iedge) < msh.tag[ithread]){
            msh.edg2tag(ithread,iedge) = msh.tag[ithread];
            lbedg.stack(iedge);
          }
        }

        continue;
      }

//    No neighbour here means the bfac is open
//    Flag and move on
      if(ifac2 == -1){
        *iopen = 1;
        if(lbedg.size1() > 0){
          int iedge = getedgglo(msh,ip1,ip2);
          METRIS_ASSERT(iedge >= 0);
          if(msh.edg2tag(ithread,iedge) < msh.tag[ithread]){
            msh.edg2tag(ithread,iedge) = msh.tag[ithread];
            lbedg.stack(iedge); 
          }
        }
        continue;
      }


      // There can be an edge between two triangles ! 
      // Why is this not after the next continue? Because while the next face 
      // may not be new, the edge may be first time seen !
      if(lbedg.size1() > 0){
        int iedge = getedgglo(msh,ip1,ip2);
        if(iedge >= 0){
          if(msh.edg2tag(ithread,iedge) < msh.tag[ithread]){
            msh.edg2tag(ithread,iedge) = msh.tag[ithread];
            lbedg.stack(iedge); 
          }
        }
      }

//    This element is already in the bfac
      if(msh.fac2tag(ithread,ifac2) >= msh.tag[ithread]) continue;


//    Add to stack and get vertex
      lbfac.stack(ifac2);
      msh.fac2tag(ithread,ifac2) = msh.tag[ithread];
    }

    ibfac++;
  }

//  *nbfac_ = nbfac;
//  if(lbedg.size1() > 0) *nbedg_ = nbedg;

  return 0;
}






///*
//Shell3 was taking over half the time in degree elevate by this point. 
//Runtime on medium.meshb: 8.9S
//Adding __restrict__ to every variable here dropped that to 7.9s !
//Unfortunately, GCC doesn't seem to support generalized restrict as 
//a compiler option. 
//ICC has -fno-alias we can use. 
//
//Note: we don't need to keep the edge that has ip1, ip2 because we're never looping on edges.
//Reduce storage and only compute if needed (outside). 
//*/
//
//int shell3(int npoin  ,int nelem  ,
//            const intAr1 & __restrict__ tet2poi,const intAr1 & __restrict__ tet2tet,
//            int * __restrict__ tag_,intAr1 &__restrict__ ttag,
//            int ipoi1  ,int ipoi2  ,int iele0  ,
//            int mshell ,int *__restrict__ nshell_,int *__restrict__ lshell,int *__restrict__ iopen){
//
//	int ishell,iface,ielem,i,j,iele2;
//	int nshell;
//
//	assert("Mshell initialized" && mshell >= 1);
//	assert("Input iele0 within bounds" && iele0 >= 0 && iele0 < nelem);
//
//	int tag = *tag_;
//	tag += 1;
//	*iopen = 0;
//
//	assert("Edge inside provided element " && getedgtet(iele0,tet2poi,ipoi1,ipoi2) > -1);
//
//	lshell[0] = iele0;
//	ttag(iele0) = tag;
//	nshell = 1;
//
//	ishell = 0;
//	while(ishell < nshell){
//		ielem = lshell[ishell];
//		assert("Shell element is alive " && !isdeadent(ielem,tet2poi));
//
//
//
//
//		for(i = 0;i<4;i++){
//			// If the face is opposite either vertex, it cannot contain the edge
//			int ip = tet2poi(ielem,i);
//			if(ip == ipoi1 || ip == ipoi2) continue;
//
//			iele2 = tet2tet(ielem,i);
//			assert("Neighbour table correct " && iele2 >= -1 && iele2 < nelem);
//			if(iele2 == -1){
//				*iopen = 1;
//				continue;
//			}
//
//			if(ttag(iele2) >= tag) continue;
//
//			if(nshell > mshell) return 1;
//
//			ttag(iele2) = tag;
//			lshell[nshell] = iele2;
//			nshell++;
//
//			assert("Neighbour has the edge" && getedgtet(iele2,tet2poi,ipoi1,ipoi2) > -1);
//
//		}
//
//		ishell ++;
//	}
//
//	*nshell_ = nshell;
//	*tag_ = tag;
//	return 0;
//}


// iopen = -1 if shell is closed
// = boundary shell element otherwise
void shell3(const MeshBase& msh,
	          int ipoi1, int ipoi2, int iele0, 
       		  int* __restrict__ nshell,
       		  intAr1&           lshell,
       		  int* __restrict__ iopen){

	assert("Input iele0 within bounds" && iele0 >= 0 && iele0 < msh.nelem);

	*iopen = -1;

	#ifndef NDEBUG
	// Will throw error if edge not in the tetra.
	getedgtet(msh,iele0,ipoi1,ipoi2);
	#endif

	lshell[0] = iele0;
	*nshell = 1;

//	printf("Ok so far \n");fflush(stdout);

	for(int ifa = 0; ifa < 4; ifa++){
		int ip = msh.tet2poi(iele0,ifa);
		if(ip == ipoi1 || ip == ipoi2) continue;
//		printf("ifa = %d \n",ifa);fflush(stdout);

		int iele2 = msh.tet2tet(iele0,ifa);
		int iele1 = iele0;
//		printf("iele1 = %d iele2 = %d \n",iele1,iele2);fflush(stdout);
		// Only two faces remain after this check. 
		// Unroll shell from either side untill hitting iele0 or boundary.
		// If iele0 hit, shell finished. Otherwise, repeat with second face. 
		while(iele2 > -1 && iele2 != iele0){
			if(*nshell >= lshell.size1()) METRIS_THROW(SMemExcept());
			lshell[*nshell] = iele2;
			(*nshell)++; 

//			printf("stack %d at %d \n",iele2,(*nshell)-1);fflush(stdout);
			for(int ifa2 = 0; ifa2 < 4; ifa2++){
				int ip2 = msh.tet2poi(iele2,ifa2); 
				if(ip2 == ipoi1 || ip2 == ipoi2) continue;

				int itmp = msh.tet2tet(iele2,ifa2); 
//				printf("  ifa2 = %d,inei = %d \n",ifa2,itmp);fflush(stdout);
				if(itmp == iele1) continue;

				iele1 = iele2;
				iele2 = itmp; 
				break;
			}
		}
		if(iele2 == -1) *iopen = iele1;
		if(iele2 == iele0) break;
	}
}



// Gather triangles surrounding edge. 
// There is no need to merge this with shell3 as we'll be operating
// completely differently. We'll be using non-manifold neighbours directly,
// rather than asking the hash table if there's a triangle stuffed
// between two tets. 
// shell2 cannot be open or closed. 
void shell2_nm(const MeshBase& msh,
	          	 int iedge, 
       		  	 intAr1&           lshell){
  METRIS_ASSERT(iedge >= 0 && iedge < msh.nedge);


  METRIS_ASSERT(!isdeadent(iedge,msh.edg2poi));

  lshell.set_n(0);

  int ip1 = msh.edg2poi(iedge,0);
  int ip2 = msh.edg2poi(iedge,1);

	int iface = msh.edg2fac[iedge];
	assert(iface >= 0 && iface < msh.nface);
  lshell.stack(iface); 

	int ied   = getedgfac(msh,iface,ip1,ip2);
	assert(ied >= 0);

	int ifac2 = iface;
	while(getnextfacnm(msh,iface,ip1,ip2,&ifac2,&ied)){
    lshell.stack(ifac2); 
	}
}



// Same as previous but start from a triangle and local edge index iedl
void shell2_nm(const MeshBase& msh,
               int iface, 
               int iedl,
               intAr1&           lshell){
  METRIS_ASSERT(iface >= 0 && iface < msh.nface);
  METRIS_ASSERT(iedl >= 0 && iedl < 3);
  METRIS_ASSERT(!isdeadent(iface,msh.fac2poi));

  int ip1 = msh.fac2poi(iface,lnoed2[iedl][0]);
  int ip2 = msh.fac2poi(iface,lnoed2[iedl][1]);

  lshell.set_n(0);
  lshell.stack(iface); 

  int ifac2 = iface;
  while(getnextfacnm(msh,iface,ip1,ip2,&ifac2,&iedl)){
    lshell.stack(ifac2); 
  }
}



} // End namespace


