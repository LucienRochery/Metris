//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../cavity/msh_cavity.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../aux_topo.hxx"
#include "../ho_constants.hxx"
#include "../low_geo.hxx"
#include "../quality/low_metqua.hxx"
#include "../mprintf.hxx"


namespace Metris{

// The boundary here is a set of edges. They will be reconnected to 
// ipins. 
// Note: if the initial cavity includes edges, these have been reconnected
// and the result is in lnewed
template <class MetricFieldType, int ideg>
int reconnect_faccav(Mesh<MetricFieldType> &msh, const MshCavity& cav, 
                     CavOprOpt &opts, CavWrkArrs &work,
                     int nedg0, double *qmax, int ithread){
  GETVDEPTH(msh);

  const int ncfac = cav.lcfac.get_n();
  const int ncedg = cav.lcedg.get_n();

  int ierro = CAV_NOERR;


	if(ncfac <= 0) return 0;

  *qmax = -1.0e30;

  bool check_qua = (opts.qmax_nec > 0 && msh.get_tdim() == 2)
                || (opts.qmax_suf > 0 && msh.get_tdim() == 2)
                || (opts.qmax_iff > 0 && msh.get_tdim() == 2);
  bool check_val = opts.fast_reject 
                || opts.max_increase_cav_geo <= 0
                || check_qua; 

  CPRINTF1(" - start reconnect_faccav check_qua = %d \n",check_qua);

  if(check_qua) *qmax = -1;

	// Edges from ipins to boundary cavity ; EDge EXists
	// Can be a glo edge or a triangle edge: typ = 1 (glo edg) or 2 (tri)
	// ent stores the entity. 
	// Initially, this is only from the new edges, nedg0 to nedge. 
	// Later as triangles are created, and new ho points, we store refs here as well. 
	//const int medex = 100;
  intAr1 &edtyp = work.edtyp;
  intAr1 &edent = work.edent;
  edtyp.set_n(0);
  edent.set_n(0);

  intAr1 &lfcco = work.lfcco;
  lfcco.set_n(0); 

  dblAr2 &lnorf = work.lnorf;
  lnorf.set_n(0);

  int nfac0 = msh.nface;

  int iref0 = -1;
  //bool multiref = false; // Used later, tag number of face refs
  
	msh.tag[ithread]++;
  // Store this so we can still substract the tag even after reconnect_tetcav
  work.tagf0 = msh.tag[ithread];

	// Tag cavity edges and triangles to determine cavity boundary. 
	for(int iedge : cav.lcedg){
		METRIS_ASSERT(iedge >= 0);
    // avoid duplicates in cavity 
    METRIS_ASSERT(msh.edg2tag(ithread,iedge) < msh.tag[ithread]);
		msh.edg2tag(ithread,iedge) = msh.tag[ithread];
	}

  // We need to colour the faces per connex components. 
  // Also compute normals per connex component and store into lfcco. 
	for(int iface : cav.lcfac){
		METRIS_ASSERT(iface >= 0);
    // avoid duplicates in cavity 
    METRIS_ASSERT(msh.fac2tag(ithread,iface) < msh.tag[ithread]);
		msh.fac2tag(ithread,iface) = msh.tag[ithread];
	}
  int ncoco = 0;
  for(int iface : cav.lcfac){
    if(msh.fac2tag(ithread,iface) > msh.tag[ithread]) continue;

    ncoco++;
    lnorf.set_n(ncoco);

    if(msh.idim == 3){
      for(int ii = 0; ii < 3; ii++) lnorf(ncoco-1,ii) = 0;
    }

    lfcco.stack(iface);
    msh.fac2tag(ithread,iface) = msh.tag[ithread] + ncoco;

    // The algo is add any old face to lfcco, then stack as work array
    // the work section is only from ncoco:-
    // At the end, seeds to the cocos are left. 
    // We still need to allow popping this guy or we're not getting started.
    while(lfcco.get_n() >= ncoco){
      int ifacs = lfcco.pop();

      if(msh.idim == 3){
        double nrmal[3];
        if(msh.CAD()){
          getnorfacCAD(msh, ifacs, nrmal);
        }else{
          getnorfacP1(msh.fac2poi[ifacs], msh.coord, nrmal);
        }
        for(int ii = 0; ii < 3; ii++) lnorf(ncoco-1,ii) += nrmal[ii];
      }

      for(int ied = 0; ied < 3; ied++){
        int ifac2 = msh.fac2fac(ifacs,ied);
        // nm -> edge -> other coco
        if(ifac2 < 0) continue;

        if(msh.fac2tag(ithread,ifac2) < msh.tag[ithread]) continue;

        //METRIS_ASSERT(msh.fac2tag(ithread,ifac2) == msh.tag[ithread]
        //      || msh.fac2tag(ithread,ifac2) == msh.tag[ithread] + ncoco);

        // Already stacked (tag on stack)
        if(msh.fac2tag(ithread,ifac2) > msh.tag[ithread]) continue;


        int itmp = msh.facedg2glo(ifacs,ied);
        if(itmp >= 0) continue;

        CPRINTF2(" - add %d to coco %d \n",ifacs,ncoco-1);
        // tag & stack
        lfcco.stack(ifac2);
        msh.fac2tag(ithread,ifac2) = msh.tag[ithread] + ncoco;
      }
    }

    lfcco.stack(iface);

  }

  METRIS_ASSERT(lfcco.get_n() == ncoco);
  if(msh.param->dbgfull){
    for(int ii = 0; ii < ncoco; ii++){
      int iface = lfcco[ii];
      METRIS_ASSERT(iface >= 0);
      int icof = msh.fac2tag(ithread,iface) - msh.tag[ithread] - 1;
      METRIS_ASSERT(icof == ii);
    }
  }

  CPRINTF1(" - %d connex component(s) \n",ncoco);

	for(int iedge = nedg0; iedge < msh.nedge; iedge++){
		METRIS_ASSERT(!isdeadent(iedge,msh.edg2poi));
		int ipoin;
		if(msh.edg2poi(iedge,0) == cav.ipins){
			ipoin = msh.edg2poi(iedge,1);
		}else if(msh.edg2poi(iedge,1) == cav.ipins){
			ipoin = msh.edg2poi(iedge,0);
		}else{
			METRIS_THROW_MSG(TopoExcept(), "Either nedg0 wrong or corrupted new edges (no ipins)");
		}
    int iedex = edtyp.get_n(); 
		msh.poi2tag(ithread,ipoin) = msh.tag[ithread] + iedex;
		edtyp.stack(1);
		edent.stack(iedge);
	}

  // Other category of iedex, if ipins on the boundary:
  // triangle edges that contain ipins, against the cavity bdry
  // edtyp = -2 as they are triangles but not internal edges
  for(int ifacl = 0; ifacl < ncfac; ifacl++){
    int iface = cav.lcfac[ifacl];
    for(int ied = 0; ied < 3; ied++){
      int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
      int ip2 = msh.fac2poi(iface,lnoed2[ied][1]);
      int ipoin;
      if(ip1 == cav.ipins){
        ipoin = ip2;
      }else if(ip2 == cav.ipins){
        ipoin = ip1;
      }else{
        continue;
      }

      if(msh.poi2tag(ithread,ipoin) >= msh.tag[ithread]) continue;

      int ifac2 = msh.fac2fac(iface,ied);
      if(ifac2 >= 0){
        METRIS_ASSERT(!isdeadent(ifac2,msh.fac2poi));
        if(msh.fac2tag(ithread,ifac2) >= msh.tag[ithread]) continue;
        // Standard neighbour, exterior face
  
        int iedex = edtyp.get_n(); 
        msh.poi2tag(ithread,ipoin) = msh.tag[ithread] + iedex;

        edtyp.stack(-2);
        edent.stack(ifac2);

      }else if(ifac2 == -1){ 
        // In this case, it's a actually an edge
        // But it may be an edge prior to the cavity call (not in the cavity)
        // In this case, we don't want to update the edge just yet.

        int iedge = getedgglo(msh,ip1,ip2);
        if(iedge < 0) continue;

        // This one no longer exists (cavity edge)
        if(msh.edg2tag(ithread,iedge) >= msh.tag[ithread]) continue;
      
        int iedex = edtyp.get_n();
        msh.poi2tag(ithread,ipoin) = msh.tag[ithread] + iedex;

        edtyp.stack(-1);
        edent.stack(iedge);

      }else{ // Non manifold
        // Some of the nm faces could be in the cavity, others not. 
        // We'll leave this case be for now 

        METRIS_THROW_MSG(TODOExcept(), "Non manifold when ipins on bdry todo")
      }
    }
  }  

  if(DOPRINTS2()){
    CPRINTF2(" - nedex = %d :\n",edtyp.get_n());
    for(int ii = 0; ii < edtyp.get_n(); ii++){
      CPRINTF2(" - %d : typ %d ent %d \n",ii,edtyp[ii],edent[ii]);
    }
  }


  // This is more of an assert type situation
  // Ensure internal edges had been added to the cavity 
  #ifndef NDEBUG
    for(int ifacl = 0; ifacl < ncfac; ifacl++){
      int iface = cav.lcfac[ifacl];
      for(int ied = 0; ied < 3 ;ied++){
        int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
        int ip2 = msh.fac2poi(iface,lnoed2[ied][1]);
        int ifac2 = msh.fac2fac(iface,ied);
        if(ifac2 == -1) continue;

        if(ifac2 >= 0){
          if(msh.fac2tag(ithread,ifac2) < msh.tag[ithread]) continue;
        }else if(ifac2 < -1){// Non manifold
          // Check any of the other nm faces are in the cavity. 
          // If yes, the edge needs to be as well. 


          int ifac3 = iface;
          bool iint = false;
          int ied2 = ied;
          while(getnextfacnm(msh,iface,ip1,ip2,&ifac3,&ied2)){
            if(msh.fac2tag(ithread,ifac3) >= msh.tag[ithread]){
              iint = true;
              break;
            }
          }

          if(!iint) continue;
        }

        int iedge = getedgglo(msh,ip1,ip2);

        if(iedge < 0 && ifac2 < 0) METRIS_THROW_MSG(TopoExcept(),"Non manifold and no edge");

        if(iedge < 0) continue;

        if(msh.edg2tag(ithread,iedge) < msh.tag[ithread]){
          CPRINTF1("## edge %d is internal but was not in cavity\n",iedge);
          // This is not always an error in the sense of an assert.
          // The assert has proved useful to spot legitimate bugs but let's downgrade it now
          ierro = CAV_ERR_INTEDG;
          goto cleanup;
//          METRIS_THROW_MSG(TopoExcept(),   "Internal edge not in cavity")
        }
      }
    }
  #endif



	for(int iface : cav.lcfac){
		// If there is a new triangle, its ref will be this triangle's ref.
		int iref  = msh.fac2ref[iface]; 
    if(iref0 < 0) iref0 = iref;
    //if(iref0 != iref) multiref = true;
    #ifndef NDEBUG
      if(msh.fac2tag(ithread,iface) <= msh.tag[ithread]){
        printf("## TAG mismatch %d <= %d \n",msh.fac2tag(ithread,iface),msh.tag[ithread]);
        printf("iface = %d \n",iface);
      }
      METRIS_ASSERT(msh.fac2tag(ithread,iface) > msh.tag[ithread]);
    #endif

		for(int ied = 0; ied < 3; ied++){

      //int ifac2 = msh.fac2fac(iface,ied);
      ierro = crenewfa<MetricFieldType,ideg>(msh,cav,opts,work,
                                 iface,ied,iref,
                                 check_val,check_qua,
                                 nedg0,nfac0,qmax,ithread);
      if(ierro > 0) goto cleanup;

		} // for ied
	}



  cleanup:
  msh.tag[ithread] += ncoco;

	//untag
	while(edtyp.get_n() > 0){
		int ityp  = edtyp.pop();
		int ientt = edent.pop();
    METRIS_ASSERT(ientt >= 0);
		if(abs(ityp) == 1){
			for(int ii = 0; ii < abs(ityp) + 1; ii++){
				int ipoin = msh.edg2poi(ientt,ii);
				msh.poi2tag(ithread,ipoin) = 0;
			}
		}else if(abs(ityp) == 2){
			for(int ii = 0; ii < abs(ityp) + 1; ii++){
				int ipoin = msh.fac2poi(ientt,ii);
				msh.poi2tag(ithread,ipoin) = 0;
			}
		}else{
			METRIS_THROW(TopoExcept());
		}
	}


	return ierro;
}

// ---------- Forward declarations
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int reconnect_faccav<MetricFieldAnalytical, n >(Mesh<MetricFieldAnalytical> &msh, const MshCavity& cav, CavOprOpt &opts,\
                        CavWrkArrs&, int nedg0, double* qmax, int ithread);\
template int reconnect_faccav<MetricFieldFE        , n >(Mesh<MetricFieldFE        > &msh, const MshCavity& cav, CavOprOpt &opts,\
                        CavWrkArrs&, int nedg0, double* qmax, int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


/*
New triangle ifacn created by cavity triangle ifac0 
Point ip, update bpoi -> this is only for points that already existed !

if ip is ipins, we may not have the info ! that is if ip is a new pt 
the driver must do the update ! we only do it topo

- If ip is face, we already knew
- If ip is edge, we know all possible (u,v)s there. 
  - We need one face in same connex component (coco) as ifac0 
    -> Assume by contradiction 2 faces i1 i2 in same coco have different (u,v)
       Then some path from i1 to i2 collapses to a point -> a squeeze, changes geometry topo: CQFD
  - Connect components correspond to same tag. 
*/
template<class MetricFieldType>
static void aux_bpo_update_fac(Mesh<MetricFieldType> &msh, 
    int nedg0, int ip, int ifacn, int ifac0, int ipins, int ithread){

  int ib = msh.poi2bpo[ip];
  METRIS_ASSERT(ib >= 0);

  // Easy case: same dim 
  if(msh.bpo2ibi(ib,1) == 2){
    int ifa = msh.bpo2ibi(ib,2);
    if(msh.fac2tag(ithread,ifa) >= msh.tag[ithread] // if face about to get deleted
    && msh.bpo2ibi(ib,3) == -1){ // and have not yet done what we're about to do  || note unlike edges there is no higher tdim
      // then create new bpo and copy the old uvs here
      int ibn = msh.template newbpotopo<2>(ip,ifacn);
      for(int jj = 0; jj < nrbi; jj++) msh.bpo2rbi(ibn,jj) = msh.bpo2rbi(ib,jj);
    }
    return;
  }

  // Edge or corner 
  // Depends whether new, or old edge. 
  int ibn = msh.template newbpotopo<2>(ip,ifacn);

  for(int jj = 0; jj < nrbi; jj++) msh.bpo2rbi(ibn,jj) = 0.0;

  // Find the entry relating to the mother face ifac0
  // Note, this may not exist, if we're not currently looking at the edge ifac0 gave ifacn. 
  ib = msh.poi2ebp(ip,2,ifac0,-1);
  //ib = getent2bpo(msh, ib, ifac0, 2);

  if(ib >= 0){
    for(int jj = 0; jj < nrbi; jj++) msh.bpo2rbi(ibn,jj) = msh.bpo2rbi(ib,jj);
    return;
  }

  // !! driver beware !!
  if(ip == ipins) return; 

  // Get a face in same connex component as ifac0, which has ip. 
  ib = msh.poi2bpo[ip];
  // We already know first is not face. 
  ib = msh.bpo2ibi(ib,3);
  int ib2 = ib;
  int tag0 = msh.fac2tag(ithread,ifac0);
  do{
    if(msh.bpo2ibi(ib2,1) == 2){
      int iface = msh.bpo2ibi(ib2,2);
      if(msh.fac2tag(ithread,iface) == tag0){
        // Found ! 
        // Note this could be a virtual face, as it might in the 
        // case of an insertion
        for(int kk = 0; kk < nrbi; kk++) msh.bpo2rbi(ibn,kk) = msh.bpo2rbi(ib2,kk);
        return;
      }
    }
    ib2 = msh.bpo2ibi(ib2,3);
  }while(ib2 != ib && ib2 > 0);

  METRIS_THROW_MSG(TopoExcept(),"Failed to find ib in old faces");
}

//int iedge = msh.bpo2ibi(ib,2);
//// Dead edge is virtual 
//if(!isdeadent(iedge,msh.edg2poi) && iedge < nedg0){ // Old edge, relevant information available.
//  ib = getent2bpo(msh, ib, ifac0, 2);
//  METRIS_ASSERT(ib >= 0);
//  for(int jj = 0; jj < nrbi; jj++) msh.bpo2rbi(ibn,jj) = msh.bpo2rbi(ib,jj);
//}

template void aux_bpo_update_fac<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, 
    int nedg0, int ip, int ifacn, int ifac0, int ipins, int ithread);
template void aux_bpo_update_fac<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, 
    int nedg0, int ip, int ifacn, int ifac0, int ipins, int ithread);




// ifac1 spawns new element from its edge ied
template<class MetricFieldType, int ideg>
int crenewfa(Mesh<MetricFieldType> &msh, const MshCavity& cav, 
             CavOprOpt &opts, CavWrkArrs &work,
             int ifac1, int ied, int iref ,
             bool check_val, bool check_qua,
             int nedg0, int nfac0, double* qmax, int ithread){
  GETVDEPTH(msh);
  METRIS_ASSERT(ithread >= 0 && ithread < METRIS_MAXTAGS);

  constexpr int nnode = facnpps[ideg];

  int ibins = msh.poi2bpo[cav.ipins];
  int tdimi = -1;
  if(msh.isboundary_faces()){
    METRIS_ASSERT(ibins >= 0);
    tdimi = msh.bpo2ibi(ibins,1);
  }

  // Neighbour/attached entity
  int ifac2 = -1;
  ifac2 = msh.fac2fac(ifac1,ied);
  if(ifac2 >= 0){ // Manifold neighbour
    METRIS_ASSERT(!isdeadent(ifac2,msh.fac2poi));
    // Belongs to cavity: continue
    if(msh.fac2tag(ithread,ifac2) >= msh.tag[ithread]) return 0;
    // Outside edge: this will create a triangle. 
  }else{ // Edge neighbour (non manifold or 2D boundary)
    int iedge = msh.facedg2glo(ifac1,ied);
    if(iedge < 0) METRIS_THROW_MSG(TopoExcept(),"No edge where face neighbourless or nm");
    // If edge in edge cavity, dead
    if(msh.edg2tag(ithread,iedge) >= msh.tag[ithread]) return 0;
    // Outside edge or non manifold but not inside cavity! this should have been caught by 
    // check cavity topo. So we can assume this is a proper cavity boundary edge.
  }



  // Whether 2D bdry or cavity bdry, this edge is not shared nor in edge cavity.
  // A new triangle will be created
  // We need to check if the other edges already exist or not. 
  int ip1 = msh.fac2poi(ifac1,lnoed2[ied][0]);
  if(ip1 == cav.ipins) return 0;
  int ip2 = msh.fac2poi(ifac1,lnoed2[ied][1]);
  if(ip2 == cav.ipins) return 0;

  int ifacn = msh.nface;

  msh.set_nface(msh.nface+1);
  msh.fac2poi(ifacn,0) = cav.ipins;
  msh.fac2poi(ifacn,1) = ip1;
  msh.fac2poi(ifacn,2) = ip2;

  for(int ii = 0; ii < METRIS_MAXTAGS; ii++) msh.fac2tag(ii,ifacn) = 0;

  // Store connex component for (u,v) updates down the line. We also use here to 
  // retrieve normal for getmeasentP1
  int icoco = msh.fac2tag(ithread,ifac1) - msh.tag[ithread] - 1;
  METRIS_ASSERT_MSG(icoco >= 0 && icoco < work.lnorf.get_n(), 
                    "icoco = "<<icoco<< " n "<<work.lnorf.get_n());
  // We don't use msh.tag[ithread] because the next functions are free
  // to touch that. 
  msh.fac2tag(ithread,ifacn) = icoco;

  if(check_val){
    bool iflat;
    double meas0;
    if(msh.idim == 2){
      meas0 = getmeasentP1<2,2>(msh, msh.fac2poi[ifacn], NULL, &iflat);
    }else if(msh.idim == 3){
      //double nrmal[3]; 
      //if(tdimi <= 1){ // NO depends on side for periodic
      //  getnorpoiref<1>(msh,cav.ipins,iref,cav.lcfac,nrmal);
      //  //msh.getpoinormal(cav.ipins,iref,nrmal);
      //}else{
      //  for(int ii = 0; ii < 3; ii++) nrmal[ii] = cav.nrmal[ii];
      //}

      //use ifac1 to compute the normal
      //double normal[3];
      //if(msh.CAD()){
      //  getnorfacCAD(msh, ifac1, normal);
      //}else{
      //getnorfacP1(msh.fac2poi[ifac1], msh.coord, normal);
      //}
      if(DOPRINTS2()){
        CPRINTF2(" - test using nrmal = ");
        dblAr1(3,work.lnorf[icoco]).print();
      }
      meas0 = getmeasentP1<3,2>(msh, msh.fac2poi[ifacn], work.lnorf[icoco], &iflat);
    }else{
      METRIS_THROW_MSG(TopoExcept(),"reconncting faces in 1D mesh...");
    }
    if(iflat){
      CPRINTF1(" - iflat ! return ip1 ip2 ip3 = %d %d %d meas = %15.7e \n", 
               cav.ipins,ip1,ip2,meas0); 
      if(DOPRINTS1()){
        if(msh.idim == 3){
          CPRINTF1(" normal = ");
          dblAr1(msh.idim,cav.nrmal).print();
        }
      }
      return CAV_ERR_FLATFAC;
    }     
    if(meas0 < 0){
      CPRINTF1(" - meas0 < 0\n");
      return CAV_ERR_NEGFAC;
    }
  }

  msh.fac2ref[ifacn]   = iref;
  msh.fac2fac(ifacn,0) = ifac2; // This neighbour is free (opposite ipins)
  if(ifac2 >= 0) METRIS_ASSERT(!isdeadent(ifac2,msh.fac2poi));
  if(ifac2 < -1) METRIS_ASSERT(!isdeadent(-ifac2-2,msh.fac2poi));

  // Set the other neighbours to zero for now. 
  msh.fac2fac(ifacn,1) = -1;
  msh.fac2fac(ifacn,2) = -1;
  if(msh.nelem > 0){
    msh.fac2tet(ifacn,0) = -1;
    msh.fac2tet(ifacn,1) = -1;
  }


  if(msh.isboundary_faces()){ // Create bpois and get uvs
    // Is it really 3 ? why no HO nodes?
    for(int ii = 0; ii < nnode; ii++){
      int ip = msh.fac2poi(ifacn,ii);
      aux_bpo_update_fac(msh,nedg0,ip,ifacn,ifac1,cav.ipins,ithread);
    }
  }
 
  // Next up, edge interior nodes and neighbours. 

  // Not very pretty but let's not assume 1, 2 is edge 0. (though it is and should remain)
  for(int iedn = 0; iedn < 3; iedn++){
    int jp1 = msh.fac2poi(ifacn,lnoed2[iedn][0]);
    METRIS_ASSERT(jp1 >= 0 && jp1 < msh.npoin);
    int jp2 = msh.fac2poi(ifacn,lnoed2[iedn][1]);
    METRIS_ASSERT(jp2 >= 0 && jp2 < msh.npoin);


    //if(jp1 != cav.ipins && jp2 != cav.ipins) continue;
    int iedex = -1;
    if(jp1 == ip1 && jp2 == ip2){ // This edge is inherited from ifac1
      cpy_facedg2facedg<ideg>(msh,ifac1,ied,ifacn,0);
      // The following two edges will use iedex

      if(msh.isboundary_faces()){ // Don't duplicate the work done if edge
        int idx0 = 3 + iedn * (edgnpps[ideg] - 2);
        for(int ii = 0; ii < edgnpps[ideg] - 2; ii++){
          int ip = msh.fac2poi[ifacn][idx0 + ii];
          aux_bpo_update_fac(msh,nedg0,ip,ifacn,ifac1,cav.ipins,ithread);
        }
      }

      continue;
    }else if(jp1 == ip2 && jp2 == cav.ipins){
      iedex = msh.poi2tag(ithread,ip2) - msh.tag[ithread];
    }else if(jp1 == cav.ipins && jp2 == ip1){
      iedex = msh.poi2tag(ithread,ip1) - msh.tag[ithread];
    }else{
      METRIS_THROW(TopoExcept());
    }

    // From now on, either an interior edge, or an exterior edge that ifac1 does
    // not know about ! 
    CPRINTF2(" - iedn = %d iedex = %d \n",iedn,iedex);


    if(iedex >= 0){ // This edge exists! Copy and update neighbours when possible
      int ityp = work.edtyp[iedex];
      METRIS_ASSERT(abs(ityp) == 1 || abs(ityp) == 2);
      CPRINTF2(" - iedex = %d ityp = %d \n",iedex,ityp);
      // The edge can either be attached to an edge, or a triangle. 
      // Copy from respective entity 
      // Then update internal neighbours
      // Can be non manifold in both cases. 
      // In case of edge, we might not know a second face yet. 
      int ifac2 = -1;
      if(abs(ityp) == 1){
        int iedge = work.edent[iedex];
        cpy_gloedg2facedg<ideg>(msh,iedge,ifacn,iedn);
        // Next up we update or exploit edg2fac. 
        if(ityp > 0){ // A new edge
          ifac2 = msh.edg2fac[iedge];
          CPRINTF2(" -> ifac2 from edg2fac = %d \n",ifac2);
          // Empty, initialize
          if(ifac2 == -1) msh.edg2fac[iedge] = ifacn;
          else if(isdeadent(ifac2,msh.fac2poi)) msh.edg2fac[iedge] = ifacn;
        }

        if(msh.isboundary_faces()){
          int idx0 = 3 + iedn * (edgnpps[ideg] - 2);
          if(ityp < 0){ // Old edge
            for(int ii = 0; ii < edgnpps[ideg] - 2; ii++){
              int ip = msh.fac2poi[ifacn][idx0 + ii];
              aux_bpo_update_fac(msh,nedg0,ip,ifacn,ifac1,cav.ipins,ithread);
            }
          }else{
            METRIS_ASSERT(cav.nrmal != NULL || !msh.CAD());
            // We don't do it now. But later, we will need to regenerate (u,v) from t. 
            // For now just create ibpoi. 
            for(int ii = 0; ii < edgnpps[ideg] - 2; ii++){
              int ip = msh.fac2poi[ifacn][idx0 + ii];
              msh.template newbpotopo<2>(ip,ifacn);
            }
          }
        }
      }else{
        ifac2 = work.edent[iedex];
        CPRINTF2(" -> ifac2 from edent = %d \n",ifac2);
        METRIS_ASSERT(ifac2 >= 0);
      }

      // Either edge is second time seen, or it's a geometric edge (ityp == 1)
      if(ifac2 >= 0){
        METRIS_ASSERT(ifac2 < msh.nface);
        METRIS_ASSERT(ifac2 != ifacn); // Never careful enough

        int ied2 = getedgfac(msh,ifac2,jp1,jp2);
        METRIS_ASSERT(ied2 >= 0);

        // Note this could be outside the if(ifac2 >= 0) as that is always verified
        if(abs(ityp) == 2)cpy_facedg2facedg<ideg>(msh,ifac2,ied2,ifacn,iedn);

        // Don't duplicate the work done if edge
        if(msh.isboundary_faces() && abs(ityp) == 2){ 
          int idx0 = 3 + iedn * (edgnpps[ideg] - 2);
          for(int ii = 0; ii < edgnpps[ideg] - 2; ii++){
            int ip = msh.fac2poi[ifacn][idx0 + ii];
            aux_bpo_update_fac(msh,nedg0,ip,ifacn,ifac1,cav.ipins,ithread);
          }
        }

        // Update neighbours. 
        if(ityp > 0){ // A new face

          // Find out what ifac2 knows already 
          int ifac3 = msh.fac2fac(ifac2,ied2);
          // Note this is a cavity interior edge, there can be no outside triangles here
          METRIS_ASSERT(ifac3 < 0 || ifac3 >= nfac0);

          // This cannot happen. 
          // Indeed, only if ifacn was the edg2fac link could someone have made us their neighbour
          // but then we would not be here as we would have had ifac2 == ifacn. 
          METRIS_ASSERT(ifac3 != ifacn);

          if(ifac3 == -1){ // easiest: both are brand new, update
            msh.fac2fac(ifacn,iedn) = ifac2;
            msh.fac2fac(ifac2,ied2) = ifacn;
          }else if(ifac3 < -1){ // Non manifold but already established as such
            // Since we are not the edg2fac, we cannot be in this nm neighbourhood yet, no-one knows about us
            ifac3 = - ifac3 - 2;
            METRIS_ASSERT(ifac3 != ifac2 && ifac3 != ifacn); // Two things that cannot happen for different reasons
            CPRINTF2(" - update ifacn(ied) nm nei = %d(%d) - %d(?) - %d(%d)\n",
              ifacn,iedn,ifac3,ifac2,ied2);
            // Right now we have ifacn to insert, ifac2 and ifac3 two neighbours. 
            // Let's sandwich ifacn between the two. 
            msh.fac2fac(ifacn,iedn) = - ifac3 - 2; // we point to ifac3..
            msh.fac2fac(ifac2,ied2) = - ifacn - 2; // and ifac2 switches from ifac3 to us
          }else{ // New non manifold
            // ifac3 and ifac2 are regular neighbours but, actually, there is a third or more triangle here
            //ifac3 = - ifac3 - 2;
            int ied3 = getedgfac(msh,ifac3,jp1,jp2);
            METRIS_ASSERT(ied3 >= 0);
            CPRINTF2(" - update ifacn(ied) nm nei = %d(%d) - %d(%d) - %d(%d)\n",
              ifacn,iedn,ifac3,ied3,ifac2,ied2);
            msh.fac2fac(ifacn,iedn) = - ifac3 - 2;
            msh.fac2fac(ifac3,ied3) = - ifac2 - 2;
            msh.fac2fac(ifac2,ied2) = - ifacn - 2;
          }
        } else{ // An old face
          if(ifac2 < 0){
            METRIS_THROW_MSG(TODOExcept(), "Handle non manifold ext neighbour");
          }else{
            msh.fac2fac(ifacn,iedn) = ifac2;
            // Only difference is we don't update ifac2 back
          }
        }// end if ityp > 0
      }
    }else{ // New edge. Create and add to edex.
      int jb1, jb2;
      if(msh.isboundary_faces()){
        // Note if the edge did not exist (as is the case here), 
        // then the new points cannot be edge points but only face. 
        jb1 = msh.poi2ebp(jp1,2,-1,iref);
        METRIS_ASSERT(jb1 >= 0);

        jb2 = msh.poi2ebp(jp2,2,-1,iref);
        METRIS_ASSERT(jb2 >= 0);
      }
      // Create new ho pts
      int idx[3] = {0};
      for(int ii = 0; ii < ideg - 1; ii++){
        idx[lnoed2[iedn][0]] = 1+ii;
        idx[lnoed2[iedn][1]] = ideg - 1 - ii;
        double u1 = idx[lnoed2[iedn][0]] / (double) ideg;
        double u2 = idx[lnoed2[iedn][1]] / (double) ideg;

        int ipnew = msh.newpoitopo(2,ifacn);
        

        for(int jj = 0; jj < msh.idim; jj++){
          msh.coord(ipnew,jj) = u1*msh.coord(jp1,jj) + u2*msh.coord(jp2,jj);
        }
        if(msh.isboundary_faces()){
          int ibnew = msh.template newbpotopo<2>(ipnew,ifacn);
          for(int jj = 0; jj < nrbi; jj++){
            msh.bpo2rbi(ibnew,jj) = u1*msh.bpo2rbi(jb1,jj) + u2*msh.bpo2rbi(jb2,jj);
          }
        }
        msh.fac2poi[ifacn][mul2nod(idx[0],idx[1],idx[2])] = ipnew;
      }

      // Update edex 
      int iedex = work.edtyp.get_n();

      work.edtyp.stack(2);
      work.edent.stack(ifacn);
      if(jp1 == cav.ipins) msh.poi2tag(ithread,jp2) = msh.tag[ithread] + iedex;
      else                 msh.poi2tag(ithread,jp1) = msh.tag[ithread] + iedex;
    } // End new edge
  }



  // This is more of a retroactive check on the edges.
  // See that we did not create an edge that is on the bdry of the cavity. 
  // This could be recovered but it becomes messy when CAD and HO is involved so just reject... 
  // To recover, "simply" disable (but no such logic exists) CAD projection for the HO nodes
  // Also, copy the HO nodes from whatever informs this face... 
  if(msh.nedge - nedg0 > 20) printf("\n\n ## IMPLEMENT HASHTAB FOR HTIS\n");
  for (int ii = 0; ii < 3 ;ii++){
    // This is only bad if there is already a triangle on the other side (not boundary !)
    int ifac2 = msh.fac2fac(ifacn,ii);
    if(ifac2 < 0) continue; 
    if(msh.fac2tag(ithread,ifac2) < msh.tag[ithread]) continue;
    int kp1 = msh.fac2poi(ifacn,lnoed2[ii][0]);
    int kp2 = msh.fac2poi(ifacn,lnoed2[ii][1]);
    for(int iedge = nedg0; iedge < msh.nedge; iedge++){
      int jp1 = msh.edg2poi(iedge,0);
      int jp2 = msh.edg2poi(iedge,1);
      if(jp1 == kp1 && jp2 == kp2 || jp1 == kp2 && jp2 == kp1){
        CPRINTF1(" previously created edge coincides with this edge");
        return CAV_ERR_DUPEDG2;
      }
    }
  }



  if constexpr(ideg >= 3){

    // Face interior nodes.
    int irnk0 = 3 + 3*(ideg-1); 
    int ipn1 = msh.fac2poi(ifacn,0);
    int ipn2 = msh.fac2poi(ifacn,1);
    int ipn3 = msh.fac2poi(ifacn,2);
    METRIS_ASSERT(ipn1 >= 0 && ipn2 >= 0 && ipn3 >= 0);
    int ibn1,ibn2,ibn3;
    if(msh.isboundary_faces()){
      ibn1 = msh.poi2bpo[ipn1];
      ibn2 = msh.poi2bpo[ipn2];
      ibn3 = msh.poi2bpo[ipn3];
      METRIS_ASSERT(ibn1 >= 0 && ibn2 >= 0 && ibn3 >= 0);
    }
    for(int irnk = irnk0; irnk < nnode; irnk++){
      double bary[3] = {
        ((double)ordfac.s[ideg][irnk][0])/ideg,
        ((double)ordfac.s[ideg][irnk][2])/ideg,
        ((double)ordfac.s[ideg][irnk][2])/ideg
      };

      int ipnew = msh.newpoitopo(2,ifacn);
      for(int jj = 0; jj < msh.idim; jj++){
        msh.coord(ipnew,jj) = bary[0]*msh.coord(ipn1,jj) 
                                 + bary[1]*msh.coord(ipn2,jj)
                                 + bary[2]*msh.coord(ipn3,jj);
      }
      if(msh.isboundary_faces()){
        msh.template newbpotopo<2>(ipnew,ifacn);
        for(int jj = 0; jj < nrbi; jj++){
          msh.bpo2rbi[msh.nbpoi - 1][jj] = bary[0]*msh.bpo2rbi(ibn1,jj) 
                                         + bary[1]*msh.bpo2rbi(ibn2,jj)
                                         + bary[2]*msh.bpo2rbi(ibn3,jj);
        }
      }
      msh.fac2poi(ifacn,irnk) = ipnew;
    }
  }




  if(check_qua){
    double quael;
    if(msh.idim == 2)
      quael = metqua<MetricFieldType,2,2>(msh,AsDeg::Pk,AsDeg::P1,
                                          ifacn,opts.qpower,
                                          opts.qpnorm,1.0);
    else
      quael = metqua<MetricFieldType,3,2>(msh,AsDeg::Pk,AsDeg::P1,
                                          ifacn,opts.qpower,
                                          opts.qpnorm,1.0);
    CPRINTF1(" - new triangle %d = %d %d %d from %d nei = %d %d %d conf error = %f \n",ifacn,
       msh.fac2poi(ifacn,0), msh.fac2poi(ifacn,1), msh.fac2poi(ifacn,2),ifac1,
       msh.fac2fac(ifacn,0), msh.fac2fac(ifacn,1), msh.fac2fac(ifacn,2), quael);
    *qmax = quael > *qmax ? quael : *qmax;
    if(quael > opts.qmax_nec && opts.qmax_nec > 0.0) return CAV_ERR_QMAXNEC; // Run rejected
    if(quael > opts.qmax_iff && opts.qmax_iff > 0.0) return CAV_ERR_QMAXIFF; // Run rejected
    if(quael < 0) return CAV_ERR_QFACNEG;  // Run rejected
  }

  return 0;
}


// ---------- Forward declarations
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int crenewfa<MetricFieldAnalytical, n>(Mesh<MetricFieldAnalytical> &msh, \
              const MshCavity& cav, CavOprOpt &opts, CavWrkArrs &work, \
              int ifac1, int ied, int iref,\
              bool check_val, bool check_qua,\
              int nedg0, int nfac0, double* qmax, int ithread);\
template int crenewfa<MetricFieldFE        , n>(Mesh<MetricFieldFE        > &msh, \
              const MshCavity& cav, CavOprOpt &opts, CavWrkArrs &work, \
              int ifac1, int ied, int iref,\
              bool check_val, bool check_qua,\
              int nedg0, int nfac0, double* qmax, int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


} // End namespace
