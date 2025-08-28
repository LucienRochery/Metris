//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "low_topo.hxx"

#include "Mesh/Mesh.hxx"

#include "aux_exceptions.hxx"
#include "aux_topo.hxx"
#include "utils/mprintf.hxx"
#include "utils/fmt_formatters.hxx"

#include <assert.h>
#include <tuple>

namespace Metris{

std::tuple<int,int,int> stup3(int i1,int i2,int i3);

// Get the points surrounding a point, among dimension tdim elements.
// Optionally fills lbent the ball of ipoin. 
int poi2poi(MeshBase& msh, int ipoin, int tdim, intAr1 &lpoin, intAr1 *lbent, int ithrd1){
  GETVDEPTH(msh.param);

  if(tdim <= 0) tdim = msh.get_tdim();

  int medge = 0, mface = 0, mtetr = 0;
  if(tdim == 1) medge = 10;
  if(tdim == 2) mface = 30;
  if(tdim == 3) mtetr = 100;

  intWrkAr1 lbedg_ = msh.get_iwork(medge);
  intWrkAr1 lbfac_ = msh.get_iwork(mface);
  intWrkAr1 lbtet_ = msh.get_iwork(mtetr);

  intAr1 dum;
  intAr1 &lbedg = tdim == 1 ? (lbent == NULL ? lbedg_.get_array() : *lbent) : dum;
  intAr1 &lbfac = tdim == 2 ? (lbent == NULL ? lbfac_.get_array() : *lbent) : dum;
  intAr1 &lbtet = tdim == 3 ? (lbent == NULL ? lbtet_.get_array() : *lbent) : dum;


  int iopen;
  int ierro = ball(msh, ipoin, lbedg, lbfac, lbtet, &iopen, false, ithrd1);
  if(ierro != 0) return ierro;

  CPRINTF1(" - ball gathered edge {} face {} tetra {}\n",
           lbedg.get_n(), lbfac.get_n(), lbtet.get_n());

  intAr1 *lbent_ = lbent;
  if(lbent == NULL){
    lbent_ = tdim == 1 ? &lbedg :
             tdim == 2 ? &lbfac : &lbtet; 
  }else{
    CPRINTF1(" - using provided lbent with size {}\n",
             lbent_->get_n());
  }
  poi2poi(msh, ipoin, tdim, *lbent_, lpoin, ithrd1);

  return 0;
}


// Get the points surrounding a point, among dimension tdim elements.
// Caller provides ball lbent.
void poi2poi(MeshBase& msh, int ipoin, int tdim, const intAr1 &lbent, intAr1 &lpoin, int ithrd1){
  GETVDEPTH(msh.param);
  const intAr2 &ent2poi = msh.ent2poi(tdim);

  lpoin.allocate(lbent.get_n());
  lpoin.set_n(0);

  const int nnode = getnnode(tdim,msh.curdeg);

  
  CPRINTF1("-- START poi2poi lbent.n = {}\n",lbent.get_n());

  msh.tag[ithrd1]++;
  for(int ientt : lbent){
    INCVDEPTH(msh.param);
    for(int inode = 0; inode < nnode; inode++){
      int ipoi2 = ent2poi(ientt,inode);
      if(ipoi2 == ipoin) continue;
      if(msh.poi2tag(ithrd1,ipoi2) >= msh.tag[ithrd1]) continue;
      msh.poi2tag(ithrd1,ipoi2) = msh.tag[ithrd1];
      lpoin.stack(ipoi2);
      CPRINTF1(" - point {} shares element {} with {}\n",ipoi2,ientt,ipoin);
    }
  }

}

// lbedg, lbfac and lbtet can be size 0 (allocated to 0)
// in that case, they will not be filled. 
int ball(MeshBase& msh, int ipoin,
         intAr1 &lbedg, intAr1 &lbfac, intAr1 &lbtet,
         int *iopen, bool append, int ithrd){

  GETVDEPTH(msh.param);

  int pdim = msh.getpoitdim(ipoin);

  bool doedg = lbedg.size1() > 0 && msh.nedge > 0 && pdim <= 1;
  bool dofac = lbfac.size1() > 0 && msh.nface > 0 && pdim <= 2;
  bool dotet = lbtet.size1() > 0 && msh.nelem > 0 && pdim <= 3;

  METRIS_ASSERT(doedg || dofac || dotet);

  CPRINTF1("-- START ball ipoin {} pdim {} gather edge {} face {} tetra {} \n",
           ipoin,pdim,doedg,dofac,dotet);

  if(lbedg.size1() && !append) lbedg.set_n(0);
  if(lbfac.size1() && !append) lbfac.set_n(0);
  if(lbtet.size1() && !append) lbtet.set_n(0);

  msh.tag[ithrd]++;

  if(append){
    CPRINTF1(" - append mode with: {} edges, {} faces, {} tetras\n",
             lbedg.get_n(),lbfac.get_n(),lbtet.get_n());
    for(int ientt : lbedg) msh.edg2tag(ithrd,ientt) = msh.tag[ithrd];
    for(int ientt : lbfac) msh.fac2tag(ithrd,ientt) = msh.tag[ithrd];
    for(int ientt : lbtet) msh.tet2tag(ithrd,ientt) = msh.tag[ithrd];
  }

  // In this case, we have boundary info that we can use.
  if(pdim < msh.idim && (  (doedg && msh.isboundary_edges()) 
                        || (dofac && msh.isboundary_faces())) ){
    int iedg0;
    for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
      int tdim = msh.bpo2ibi(ibpoi,1);
      CPRINTF1(" - ibpoi {} tdim {} ientt {} \n",ibpoi,tdim,msh.bpo2ibi(ibpoi,2));
      if(tdim == 0) continue;
      if(tdim == 1 && !doedg) continue;
      int ientt = msh.bpo2ibi(ibpoi,2);
      if(tdim == 1){
        iedg0 = ientt;
        CPRINTF1(" - edge {} tagged {} <? {} = tag\n",
                 ientt, msh.edg2tag(ithrd,ientt), msh.tag[ithrd]);
        if(msh.edg2tag(ithrd, ientt) < msh.tag[ithrd]){
          lbedg.stack(ientt);
          msh.edg2tag(ithrd, ientt) = msh.tag[ithrd];
        }
      }
      else{
        if(msh.fac2tag(ithrd, ientt) < msh.tag[ithrd]){
          lbfac.stack(ientt);
          msh.fac2tag(ithrd, ientt) = msh.tag[ithrd];
        }
      }
    }

    // If the point is dim 1, we still need the edges, but we're seeded
    if(pdim == 1 && doedg){
      METRIS_ASSERT(lbedg.get_n() > 0);
      METRIS_ASSERT(iedg0 >= 0);
      METRIS_ASSERT(!isdeadent(iedg0,msh.edg2poi));
      int iver = msh.getveredg<1>(iedg0, ipoin);
      #ifndef NDEBUG
      if(iver < 0){
        PRINTF("iedg0 {} ipoin {} \n",iedg0, ipoin);
        PRINTF("ipoin full bpo list:\n");
        for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
          PRINTF("{} : {}\n",ibpoi, intAr1(nibi,msh.bpo2ibi[ibpoi]));
        }
        PRINTF("lbedg = {}\n",lbedg);
        for(int ii : lbedg){
          PRINTF("{} : {}\n",ii,intAr1(2,msh.edg2poi[ii]));
        }
      }
      METRIS_ASSERT_MSG(iver >= 0,"got iver < 0 with iedg0 = {} ipoin = {}",iedg0, ipoin);
      #endif
      int iedg1 = msh.edg2edg(iedg0,1-iver);
      CPRINTF1(" - case pdim = 1 w/ edge, seed iedg0 = {} neighbour {}\n",iedg0,iedg1);
      if(iedg1 >= 0){
        if(msh.edg2tag(ithrd,iedg1) < msh.tag[ithrd]){
          msh.edg2tag(ithrd,iedg1) = msh.tag[ithrd];
          lbedg.stack(iedg1);
        }
      }else if(iedg1 < -1){
        int inei;
        iedg1 = -iedg1 - 2;
        if(msh.edg2tag(ithrd,iedg1) < msh.tag[ithrd]){
          msh.edg2tag(ithrd,iedg1) = msh.tag[ithrd];
          lbedg.stack(iedg1);
        }
        while(getnextedgnm(msh,iedg0,ipoin,&iedg1,&inei)){
          if(msh.edg2tag(ithrd,iedg1) < msh.tag[ithrd]){
            msh.edg2tag(ithrd,iedg1) = msh.tag[ithrd];
            lbedg.stack(iedg1);
          }
        }
      }
    }

    CPRINTF1(" - ball gathered boundary: {} edges {} faces\n",
             lbedg.get_n(), lbfac.get_n());
  }

  // Edges will always have been filled by the previous case. 
  // Faces are filled if the point is dim 1 and the mesh is dim 3
  if(!dotet && pdim <= 1 && msh.idim == 3) return 0;


  // Seed ball
  int iface = -1;
  if(pdim <= 1){
    int iedge = msh.poi2ent(ipoin,0);
    iface = msh.edg2fac[iedge];
  }else if(pdim == 2){
    iface = msh.poi2ent(ipoin,0);
  }

  if(dofac){
    METRIS_ASSERT(iface >= 0);
    if(msh.fac2tag(ithrd,iface) < msh.tag[ithrd]){
      msh.fac2tag(ithrd,iface) = msh.tag[ithrd];
      lbfac.stack(iface);
    }
  }

  if(dotet){
    int itetr = -1;
    if(pdim <= 2){
      METRIS_ASSERT(iface >= 0);
      itetr = msh.fac2tet(iface,0);
      if(itetr < 0) itetr = msh.fac2tet(iface,1);
    }else{
      itetr = msh.poi2ent(ipoin,0);
    }
    METRIS_ASSERT(itetr >= 0);
    if(msh.tet2tag(ithrd,itetr) < msh.tag[ithrd]){
      msh.tet2tag(ithrd,itetr) = msh.tag[ithrd];
      lbtet.stack(itetr);
    }
  }

  if(!dofac && !dotet) return 0;

  // Minimum tdim. 1 is done. 2 iff pdim > 1 or idim = 2
  int tdi0 = ((pdim > 1 || msh.idim == 2) && dofac) ? 2 : 3;
  int tdi1 = msh.get_tdim();
  if(tdi1 == 3 && !dotet) tdi1 = 2;

  for(int tdim = tdi0; tdim <= tdi1; tdim++){

    const intAr2 &ent2poi = msh.ent2poi(tdim);
    const intAr2 &ent2ent = msh.ent2ent(tdim);
          intAr2 &ent2tag = msh.ent2tag(tdim);

    intAr1 &lbent = tdim == 2 ? lbfac : lbtet;
    METRIS_ASSERT_MSG(lbent.get_n() > 0, "lbent is empty with "
      "tdim = {} do_edge = {} do_face = {} do_tetra = {}", 
      tdim, doedg, dofac, dotet);

    ent2tag(ithrd, lbent[0]) = msh.tag[ithrd];

    int ibent = 0;
    while(ibent < lbent.get_n()){
      int ientt = lbent[ibent];

      ibent++;

      // If appending, check the element containts the point.
      if(append){
        int iver = msh.getverent(ientt, tdim, ipoin);
        if(iver < 0) continue;
      }

      for(int inei = 0; inei < tdim + 1; inei++){
        // facet opposite ipoin does not contain ipoin
        if(ent2poi(ientt,inei) == ipoin) continue;

        int ient2 = ent2ent(ientt,inei);
        if(ient2 >= 0 && ent2tag(ithrd,ient2) >= msh.tag[ithrd]) continue;

        if(ient2 >= 0){
          ent2tag(ithrd,ient2) = msh.tag[ithrd];
          lbent.stack(ient2);
          continue;
        }

        if(ient2 == -1){
          *iopen = 1;
          continue;
        }

        // Lastly ient2 < -1 only if tdim == 2, non manifold case
        METRIS_ASSERT(tdim == 2);

        // If one has been seen, all nm faces have been seen. 
        if(msh.fac2tag(ithrd,-ient2-2) >= msh.tag[ithrd]) continue;

        int ip1 = msh.fac2poi(ientt,lnoed2[inei][0]);
        int ip2 = msh.fac2poi(ientt,lnoed2[inei][1]);
        int ied2 = inei;
        int ifac3 = ientt;
        while(getnextfacnm(msh,ientt,ip1,ip2,&ifac3,&ied2)){
          // If a nm face has been seen, all should have been seen. 
          METRIS_ASSERT(msh.fac2tag(ithrd,ifac3) < msh.tag[ithrd]);
          lbfac.stack(ifac3);
          msh.fac2tag(ithrd,ifac3) = msh.tag[ithrd];
        }// getnextfacnm


      }// for inei
    }// while ibent
    CPRINTF1(" - ball gathered {} dim {} entities\n",lbent.get_n(),tdim);

  }

  return 0;
}



#if 0
int shell(MeshBase& msh, int ipoi1, int ipoi2,
          intAr1 &lbedg, intAr1 &lbfac, intAr1 &lbtet,
          int *iopen, int ithrd){

  GETVDEPTH(msh.param);

  bool doedg = lbedg.size1() > 0 && msh.nedge > 0;
  bool dofac = lbfac.size1() > 0 && msh.nface > 0;
  bool dotet = lbtet.size1() > 0 && msh.nelem > 0;

  METRIS_ASSERT(doedg || dofac || dotet);

  CPRINTF1("-- START shell ipoi1 {} ipoi2 {} gather edge {} face {} tetra {} \n",
           ipoi1,ipoi2,doedg,dofac,dotet);

  if(doedg) lbedg.set_n(0);
  if(dofac) lbfac.set_n(0);
  if(dotet) lbtet.set_n(0);

  if(doedg){

  }

  // In this case, we have boundary info that we can use.
  if(pdim < msh.idim && (  doedg && msh.isboundary_edges() 
                        || dofac && msh.isboundary_faces()) ){
    for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
      int tdim = msh.bpo2ibi(ibpoi,1);
      if(tdim == 0) continue;
      if(tdim == 1 && !doedg) continue;
      int ientt = msh.bpo2ibi(ibpoi,2);
      if(tdim == 1) lbedg.stack(ientt);
      else          lbfac.stack(ientt);
    }
    CPRINTF1(" - ball gathered boundary: {} edges {} faces\n",
             lbedg.get_n(), lbfac.get_n());
  }

}
#endif







// Basic manifold ball, disregarding internal faces. 
// lball(:,0) stores elements
// lball(:,1) stores rank of point in elements
int ball3(MeshBase& __restrict__ msh,
           int ipoin    , int iele0, 
           intAr1& lball, 
           int* __restrict__ iopen,
           int ithrd){
	int iball,ielem,i,iele2;

	METRIS_ASSERT_MSG(iele0 >= 0 && iele0 < msh.nelem,"Input iele0 within bounds");
	METRIS_ASSERT_MSG(( ipoin == msh.tet2poi(iele0,0)
                    ||ipoin == msh.tet2poi(iele0,1)
                    ||ipoin == msh.tet2poi(iele0,2)
                    ||ipoin == msh.tet2poi(iele0,3)),"Sought point is inside element");
	METRIS_ASSERT_MSG(!isdeadent(iele0,msh.tet2poi),"Input ielem alive");

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
    int iver = msh.getvertet<1>(ielem, ipoin); 

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
              [[maybe_unused]] intAr1&           lbfac, 
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
    int iver = msh.getvertet<1>(ielem, ipoin); 

		assert("Ball element is alive " && !isdeadent(ielem,msh.tet2poi));

//  Loop over neighbours skipping those already in ball
		for(i=0;i<4;i++){
			// This face is opposite the vertex: ball boundary
			if(iver == i) continue;

			//// The array was devised to avoid unnecessary hash table lookups
			//if(msh.tet2ftg[ielem]){
			//	int iface = msh.tetfac2glo(ielem,i);
			//	if(iface < 0 || iface >= msh.nface) 
			//		METRIS_THROW_MSG("Face missing or invalid in hash tab "<<iface<<"\n");
      //  nbfac++; 
			//	METRIS_ASSERT(nbfac <= lbfac.size1());
			//	lbfac[nbfac-1] = iface; 
			//}

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
      METRIS_ENFORCE(nball < lball.size1());

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
void ball3_full([[maybe_unused]] MeshBase& __restrict__ msh,
                [[maybe_unused]] int ipoin  ,
                [[maybe_unused]] int tdimn, 
                [[maybe_unused]] int iseed, 
                [[maybe_unused]] int* __restrict__ nbtet,
                [[maybe_unused]] int* __restrict__ nbfac,
                [[maybe_unused]] int* __restrict__ nbedg,
                [[maybe_unused]] intAr1&           lbtet,
                [[maybe_unused]] intAr1&           lbfac, 
                [[maybe_unused]] intAr1&           lbedg,
                [[maybe_unused]] int ithread){

  METRIS_THROW_MSG("TODO: All wrong here, reimplement")

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
    METRIS_THROW_MSG("TODO: edge ball3 not implemented");
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
        if(*nbedg >= lbedg.n) METRIS_THROW_MSG("Increase lbedg.ne")
        #ifndef NDEBUG 
          METRIS_ASSERT(msh.edg2tag(ithread,ientt) < msh.tag[ithread]);
          msh.edg2tag(ithread,ientt) = msh.tag[ithread];
        #endif
        lbedg[*nbedg] = ientt;
        (*nbedg)++;
      }else if(itype == 2){
        if(*nbfac >= lbfac.n) METRIS_THROW_MSG("Increase lbfac.ne")
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
    METRIS_THROW_MSG("TODO: eleme ball3 not implemented")
    iele0 = iseed;    
  }
  */
}



// Can handle non manifold meshes. 
// imani returns wheter manifold (true) or non-manifold (false)
// 
// lbfac only stores iface, not the index. Just call getverfac if needed. 
// 
// If lbedg.n > 0, gather internal edges + bdry if iopen 
int ball2(MeshBase& __restrict__ msh,
          int ipoin  ,int ifac0, 
          intAr1&           lbfac,
          intAr1&           lbedg,
          int* __restrict__ iopen,
          bool* __restrict__ imani,
          int ithread){
  int ibfac,iface,ifac2;
  GETVDEPTH(msh.param);

  METRIS_ASSERT_MSG(ifac0 >= 0 && ifac0 < msh.nface,"Input ifac0 within bounds");
  METRIS_ASSERT_MSG(ipoin == msh.fac2poi(ifac0,0)
                  ||ipoin == msh.fac2poi(ifac0,1)
                  ||ipoin == msh.fac2poi(ifac0,2),"Sought point is inside element");
  METRIS_ASSERT_MSG(!isdeadent(ifac0,msh.fac2poi),"Input iface alive");

  *imani = true; 

  lbfac.set_n(0);
  lbedg.set_n(0);

  int pdim = msh.getpoitdim(ipoin);
  if(pdim == 3){
    CPRINTF1("-- ball2: no faces attached to dim 3 point {} \n",ipoin);
    return 0;
  }

  if(msh.isboundary_faces()){
    CPRINTF1("-- ball2 using boundary info to gather triangles and edges\n");
    for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
      int tdim = msh.bpo2ibi(ibpoi,1);
      if(tdim == 0) continue;
      if(tdim == 1 && lbedg.size1() == 0) continue;
      int ientt = msh.bpo2ibi(ibpoi,2);
      if(tdim == 1) lbedg.stack(ientt);
      else          lbfac.stack(ientt);
    }
    return 0;
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
            PRINTF("## DEBUG THIS CASE SHOULD NEVER HAPPEN\n");
            PRINTF("Indeed if a nm face has been seen, all should have been.\n");
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





// = boundary shell element otherwise
// Store boundary elements in lbdry if provided with non-zero allocation
// does not support non-manifold -> only faces are on the boundary
void shell(const MeshBase& msh,
           int ipoi1, int ipoi2, 
           int tdim , int iele0, 
           intAr1& lsedg, intAr1& lsfac, intAr1& lstet, int* iopen){

  METRIS_ASSERT(tdim >= 1 && tdim <= 3);
  METRIS_ASSERT_MSG(iele0 >= 0 && iele0 < msh.nentt(tdim), "Input iele0 within bounds");
  METRIS_ASSERT(getedgent(msh,tdim,iele0,ipoi1,ipoi2) >= 0);

  GETVDEPTH(msh.param);

  bool ibdry;
  *iopen = 0;

  int iedg0 = -1, ifac0 = -1, itet0 = -1;
       if(tdim == 1) iedg0 = iele0;
  else if(tdim == 2) ifac0 = iele0;
  else               itet0 = iele0;

  bool doedg = lsedg.size1() > 0 && msh.nedge > 0;
  bool dofac = lsfac.size1() > 0 && msh.nface > 0;
  bool dotet = lstet.size1() > 0 && msh.nelem > 0;

  if(doedg) lsedg.set_n(0);
  if(dofac) lsfac.set_n(0);
  if(dotet) lstet.set_n(0);


  CPRINTF1("-- START shell seed {} dim {} edge {} {} doedg {} dofac {} dotet {}\n",
           iele0,tdim,ipoi1,ipoi2, doedg, dofac, dotet);

  // Skip straight to shell3-like; only if tet seeded and edge not requested
  if(tdim == 3 && dotet && !doedg){
    CPRINTF1(" - using tet seed, skip to regular shell as no edges required\n");
    goto doshell3;
  }
  // Seed tet if edge not needed and go to shell3
  if(tdim == 2 && dotet && !doedg){
    itet0 = msh.fac2tet(ifac0,0);
    if(itet0 < 0) itet0 = msh.fac2tet(ifac0,1);
    METRIS_ASSERT(itet0 >= 0);
    CPRINTF1(" - using face seed, get attached tet and skip to regular shell as no edges required\n");
    goto doshell3;
  }


  ibdry = msh.poi2bpo[ipoi1] >= 0 && msh.poi2bpo[ipoi2] >= 0;


  // Seed edge. We need this only if doedg is true. 
  // Both in cases 2D and 3D, this can be an edge only if 
  // ipoi1 and ipoi2 are both boundary.
  if(doedg && ibdry){
    iedg0 = getedgglo(msh,ipoi1,ipoi2);
    CPRINTF1(" - looking for edge {} {} -> got {}\n",ipoi1,ipoi2,iedg0);
  }
  // Either seeded or provided as argument, the only possible shell edge: 
  if(iedg0 >= 0) lsedg.stack(iedg0);

  // Seed face. We need this only if:
  //  - dofac is true
  //  - dotet is true and tdim == 1 (if it is tdim 2, we already have the face)
  if(ifac0 < 0 && (dofac || (dotet && tdim == 1))){
    if(iedg0 >= 0){
      ifac0 = msh.edg2fac[iedg0];
    }
    // if iedg0 < 0, we are necessarily in case tdim == 3
    // we can't discover the face by using just the seed tetrahedron, it could
    // be interior to an open ball. Hence we delay. 
  }

  // Seed tetra. At this stage, if tdim was originally 1, we have ifac0. 
  // If it was originally 2, we also have ifac0. 
  if(itet0 < 0 && dotet){
    itet0 = msh.fac2tet(ifac0,0);
    if(itet0 < 0) itet0 = msh.fac2tet(ifac0,1);
    METRIS_ASSERT(itet0 >= 0);
  }

  // Fill triangle shell
  // If tdim == 3 and iedg0 < 0, we do not have a face seed yet. 
  // If dotet, do nothing as we'll gather the faces naturally with the tet shell.
  // Otherwise, try filling using boundary info. 
  // If not ibdry, there is no shell in case tdim == 3 which assumes idim == 3. 
  if(ifac0 < 0 && dofac && !dotet && ibdry){
    METRIS_ASSERT(tdim == 3);
    METRIS_THROW_MSG("TODO: Implement nm / surf shell in case tet seed and no tet shell")
  }else if(dofac && !dotet){
    METRIS_ASSERT(msh.idim == 2 || msh.nelem == 0);
    METRIS_ASSERT(ifac0 >= 0);
    // 2D case 
    lsfac.stack(ifac0);

    int ied = getedgfac(msh,ifac0,ipoi1,ipoi2);
    METRIS_ASSERT(ied >= 0);
    int ifac1 = msh.fac2fac(ifac0,ied);

    if(ifac1 >= 0){
      lsfac.stack(ifac1);
      goto endfun;
    }

    if(ifac1 == -1){
      *iopen = 1;
      goto endfun;
    }

    ifac1 = ifac0;
    while(getnextfacnm(msh,ifac0,ipoi1,ipoi2,&ifac1,&ied)){
      lsfac.stack(ifac1);
    }
    
    goto endfun;
  }




  if(!dotet) goto endfun;

  doshell3:

  METRIS_ASSERT(itet0 >= 0);
  lstet.stack(itet0);
  

  for(int ifa = 0; ifa < 4; ifa++){
    int ip = msh.tet2poi(itet0,ifa);
    if(ip == ipoi1 || ip == ipoi2) continue;

    int iele2 = msh.tet2tet(itet0,ifa);
    int iele1 = itet0;
    int ifa1  = ifa;
    // Only two faces remain after this check. 
    // Unroll shell from either side untill hitting itet0 or boundary.
    // If itet0 hit, shell finished. Otherwise, repeat with second face. 
    while(iele2 >= 0 && iele2 != itet0){
      INCVDEPTH(msh.param);
      lstet.stack(iele2);
      CPRINTF1(" - check iele2 = {} \n",iele2);

      for(int ifa2 = 0; ifa2 < 4; ifa2++){
        int ipoin = msh.tet2poi(iele2,ifa2); 
        if(ipoin == ipoi1 || ipoin == ipoi2) continue;

        int itmp = msh.tet2tet(iele2,ifa2);
        if(itmp == iele1) continue;

        ifa1  = ifa2;
        iele1 = iele2;
        iele2 = itmp; 
        break;
      }
    }
    if(iele2 < 0){
      if(dofac){
        int iface = msh.tetfac2glo(iele1,ifa1);
        CPRINTF1(" - iele2 < 0 iface = {} is face {} of {}\n",
                 iface,ifa,iele1);
        METRIS_ASSERT(iface >= 0);
        lsfac.stack(iface);
      }
      *iopen = 1;
    }
    if(iele2 == itet0) break;
  }

  endfun:
  CPRINTF1("-- END shell nedge {} nface {} ntetr {}\n",
           lsedg.get_n(), lsfac.get_n(), lstet.get_n());
}







// iopen = -1 if shell is closed
// = boundary shell element otherwise
// Store boundary elements in lbdry if provided with non-zero allocation
// does not support non-manifold -> only faces are on the boundary
void shell3(const MeshBase& msh,
	          int ipoi1, int ipoi2, int iele0, 
        	  intAr1& lshell, intAr1& lbdry, int* iopen){

	METRIS_ASSERT_MSG(iele0 >= 0 && iele0 < msh.nelem, "Input iele0 within bounds");
  GETVDEPTH(msh.param);

	*iopen = -1;

	#ifndef NDEBUG
	// Will throw error if edge not in the tetra.
	METRIS_ASSERT(getedgtet(msh,iele0,ipoi1,ipoi2) >= 0);
	#endif

  bool dofac = lbdry.size1() > 0;
  if(dofac) lbdry.set_n(0);

  lshell.allocate(10);
  lshell.set_n(0);
  lshell.stack(iele0);

  CPRINTF1("-- START shell3 seed {} edge {} {}\n",iele0,ipoi1,ipoi2);

	for(int ifa = 0; ifa < 4; ifa++){
		int ip = msh.tet2poi(iele0,ifa);
		if(ip == ipoi1 || ip == ipoi2) continue;

		int iele2 = msh.tet2tet(iele0,ifa);
		int iele1 = iele0;
    int ifa1  = ifa;
		// Only two faces remain after this check. 
		// Unroll shell from either side untill hitting iele0 or boundary.
		// If iele0 hit, shell finished. Otherwise, repeat with second face. 
		while(iele2 >= 0 && iele2 != iele0){
      INCVDEPTH(msh.param);
      lshell.stack(iele2);
      CPRINTF1(" - check iele2 = {} \n",iele2);

			for(int ifa2 = 0; ifa2 < 4; ifa2++){
				int ipoin = msh.tet2poi(iele2,ifa2); 
				if(ipoin == ipoi1 || ipoin == ipoi2) continue;

				int itmp = msh.tet2tet(iele2,ifa2);
				if(itmp == iele1) continue;

        ifa1  = ifa2;
				iele1 = iele2;
				iele2 = itmp; 
				break;
			}
		}
		if(iele2 < 0){
      if(dofac){
        int iface = msh.tetfac2glo(iele1,ifa1);
        CPRINTF1(" - iele2 < 0 iface = {} is face {} of {}\n",
                 iface,ifa,iele1);
        METRIS_ASSERT(iface >= 0);
        lbdry.stack(iface);
      }
      *iopen = iele1;
    }
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


