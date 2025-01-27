//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../cavity/msh_cavity.hxx"

#include "../aux_topo.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../Mesh/Mesh.hxx"


namespace Metris{

template<class MFT>
int check_cavity_topo(Mesh<MFT> &msh, MshCavity &cav, 
                      CavOprOpt &opts, //RoutineWorkMemory<int> &iwrk, 
                      int ithread){
  int ierro = 0;
  int iverb = msh.param->iverb;

  // Check no rem pts
  if(!opts.allow_remove_points){
    
    msh.tag[ithread]++;
    for(int iedge : cav.lcedg) msh.edg2tag(ithread,iedge) = msh.tag[ithread];
    for(int iface : cav.lcfac) msh.fac2tag(ithread,iface) = msh.tag[ithread];
    for(int ielem : cav.lctet) msh.tet2tag(ithread,ielem) = msh.tag[ithread];

    // Points to be removed are those that are surrounded by only cavity elements.
    // Hence, loop over cavity elements and tag any points that belong to a 
    // non-cavity neighbour. 
    // Lastly, count untagged vertices. 

    // If it belongs to any lower dim elements, that should be in the cavity. 
    // It suffice there is one, as if it doesnt belong to all, it would be tagged. 

    int tdimn = cav.lctet.get_n() > 0 ? 3 
              : cav.lcfac.get_n() > 0 ? 2 
                                      : 1;
    const intAr1& lcent = cav.lcent(tdimn);
    const intAr2& ent2ent = msh.ent2ent(tdimn);
    const intAr2& ent2poi = msh.ent2poi(tdimn);
    const intAr2& ent2tag = msh.ent2tag(tdimn);


    for(int ientt : lcent){
      for(int ii = 0; ii < tdimn + 1; ii++){
        int ipoin = ent2poi(ientt,ii);
        // Cycle neighbours that have ii (i.e. all but ii-th neighbour)
        for(int jj = 0; jj < tdimn + 1; jj++){
          if(jj == ii) continue;
          int ient2 = ent2ent(ientt,jj);
          if(ient2 < 0) continue;
          // Tag point if the adjacent element is not in the cavity 
          // This point is not set to be deleted. 
          if(ent2tag(ithread,ient2) < msh.tag[ithread]){
            msh.poi2tag(ithread,ipoin) = msh.tag[ithread];
            if(iverb >= METRIS_CAV_PRTLEV + 1) printf("  - not rem point %d \n", ipoin);

          }

        }
      }
    }

    // Go over elements, counting vertices that have not been tagged.
    for(int ientt : lcent){
      for(int ii = 0; ii < tdimn + 1; ii++){
        int ipoin = ent2poi(ientt,ii);
        if(ipoin == cav.ipins) continue;
        if(msh.poi2tag(ithread,ipoin) >= msh.tag[ithread]) continue;
        if(iverb >= METRIS_CAV_PRTLEV + 1) printf("  - rem pt ? %d \n", ipoin);

        // point going to be deleted, but only if any existing lower dim entities
        // are also in the cavity. 
        if(tdimn == 3){
          // If there is a face attached, check it is in the cavity.
          int iface = getpoifac(msh, ipoin);
          if(iface >= 0){
            // If not, this point won't be removed. Continue. 
            if(msh.fac2tag(ithread,iface) < msh.tag[ithread]) continue;
          }

        } 

        if(tdimn >= 2){
          // If there is an edge attached, check it is in the cavity.
          int iedge = getpoiedg(msh,ipoin);
          if(iedge >= 0){
            // If not, this point won't be removed. Continue. 
            if(msh.edg2tag(ithread,iedge) < msh.tag[ithread]) continue;
          }
        }

        // If we're here, that means that there are either no attached lower dim
        // or there are and they are all in the cavity; indeed, assume there exist
        // at least one, and at least one not in the cav. Then the point is not
        // tagged. Then we wouldn't be here. 

        if(iverb >= METRIS_CAV_PRTLEV) 
          printf(" ## norempts and point %d will be removed \n",ipoin);
        
        return 1;


      }
    }

  } // if(!opts.allow_remove_points)


  return ierro;
}







#if 0


// This proceeds in two main phases:
// First, from highest to lowest topo dim, add to tdim-1 cavity
// those elements stuffed in tdim cavity. Example: triangle between
// two tetrahedra. 
// Second, from lowest to highest topo dim, add higher tdim cavity 
// elements. This has some overlap but also new elements. 

// Example:
// Initial cavity contains two triangles. Between these is edge iedge
// but it does not belong to the cavity. This cannot be as 
// the edge between triangles would be collapsed without iedge being
// accounted for. 
// Thus the edge is added to the cavity. 
// On phase two, this edge's triangle shell is computed. 
// This includes the original two triangles, but it can contain
// many more. These cannot be ignored either, as they share an edge
// with triangles to be collapsed. 
// Idem tetrahedra. 

// For each cavity edge, all adjacent triangles must be in triangle cavity
// For each cavity triangle, all adjacent tetrahedra must be in tet cavity
// If allow_remove_points == false, no points can be interior to the cavity
// If allow_remove_corners == false, no corners can be interior to the cavity
// Else, only one corner can be removed, and cav.ipins must be a corner.
// If cav.ipins is a corner, one corner must be removed. 

// Todo: add check rempts for interior triangle, tetra
// BENCHMARK THIS

// Previously a RoutineWorkMemory allocated here but "fast" in boost_fast_allocator is a "suggestion" only (72.8% time spent here)
// lball size 2 mball, lbfac size mball
template<class MFT>
int check_cavity_topo(Mesh<MFT> &msh, MshCavity &cav, CavOprOpt &opts, RoutineWorkMemory<int> &iwrk){
//                      int mball, int *lball, int *lbfac){


	static_assert(METRIS_MAXTAGS >= 2);
	// ball3 uses tags
	const int ITAG = 1;
  const int ITAG2 = ITAG + 1;
  if(ITAG2 > METRIS_MAXTAGS) METRIS_THROW_MSG(TopoExcept(), "Increase METRIS_MAXTAGS by 1") 

	int ierro = 0;


  // As some points will be tagged multiple times, we must ensure to properly increase tag in the end.
	cav.nrempts = 0;

  #ifndef NDEBUG
    if(iverb >= METRIS_CAV_PRTLEV){
      printf("## ADD CHECKS ON LINE LOOPS?: -- i -- o\n");
      printf("                                \\-/ \n");
      printf("## Add 'top to bottom' checks? tets adding triangles and edges etc.\n");
    }
  #endif

	//RoutineWorkMemory<int> iwrk(msh.iwrkmem);
	int mball = 100;
  intAr1 lball(mball,  iwrk.allocate(mball));
  intAr1 lbfac(mball, iwrk.allocate(mball));
	//int *lball = iwrk.allocate(2*mball);
	//int *lbfac = iwrk.allocate(  mball);

	// We are going to populate the cavity, and some checks should 
	// only concern the original elements. 
	const int nced0 = cav.lcedg.size1();
	const int ncfa0 = cav.lcfac.size1();
	const int ncte0 = cav.lctet.size1();

	// Start by tagging all cavity elements. 
	msh.tag[ITAG]++;
	int cavtag = msh.tag[ITAG];


/* 
	I. Edge checks:
	 a) Check attached to existing cavity triangles and tetrahedra. 
	 b) Gather rempoints and rem corner. 
*/
	for(int ifacl = 0; ifacl < cav.ncfac; ifacl++){
		int iface = cav.lcfac[ifacl]; 
		METRIS_ASSERT(iface >= 0 && iface < msh.nface);
		METRIS_ASSERT_MSG(!isdeadent(iface,msh.fac2poi),"Face "<<iface<<" dead\n");
		msh.fac2tag(ITAG,iface) = cavtag;
	}
	for(int ielel = 0; ielel < cav.nctet; ielel++){
		int ielem = cav.lctet[ielel];
		METRIS_ASSERT(ielem >= 0 && ielem < msh.nelem);
		METRIS_ASSERT(!isdeadent(ielem,msh.tet2poi));
		msh.tet2tag(ITAG,ielem) = cavtag;
	}

	/*
	Edges STEP 1: gather collapsed points and check balls in cav
	After this, only edges with no collapsed points will need shell3
	called on them. 
	We may also exit here if corners are removed or no points should be 
	collapsed.
	Points will be tagged: 
		< tag - 1: never seen
    = tag - 1: seen once
    = tag    : seen twice (regular collapse) or more (corner collapse)
	*/
	msh.tag[ITAG] += 2;
	int ncorn = 0;
	for(int iedgl = 0; iedgl < cav.ncedg; iedgl++){
		int iedge = cav.lcedg[iedgl];
		for(int ip = 0; ip < 2; ip++){
			int ipoin = msh.edg2poi(iedge,ip); 
			if(msh.poi2tag(ITAG,ipoin) < msh.tag[ITAG] - 1){
				msh.poi2tag(ITAG,ipoin) = msh.tag[ITAG] - 1;
			}else if(msh.poi2tag(ITAG,ipoin) == msh.tag[ITAG] - 1){
				msh.poi2tag(ITAG,ipoin) = msh.tag[ITAG]; 
				// This point will be removed
				if(!opts.allow_remove_points && ipoin != cav.ipins) return 5;
				cav.nrempts ++;

				if(iverb >=3) printf(" -> collapse %d nrempt = %d\n",ipoin,cav.nrempts);

				// Check all higher-dim (tri and tet) entities having ipoin are in the cavity
				// If point attached to tet, get both tets and faces this way. 
				int iele0 = getpoitet(msh,ipoin);
				if(iverb >= METRIS_CAV_PRTLEV) printf(" -> found (possibly not a) tetrahedron %d\n",iele0);
				if(iele0 > 0){
					int nball,nbfac,iopen;
					ball3_nm(msh, ipoin, iele0, &nball, &nbfac, lball, lbfac, &iopen, ITAG2);
					if(iopen) printf("## WARNING OPEN BALL3\n");
	
					for(int ii = 0; ii < nball; ii++){
						int ielem = lball[ii];
						if(msh.tet2tag(ITAG,ielem) == cavtag) continue;
						// An element in the ball does not belong to the cavity
						if(!opts.allow_topological_correction) return 3;
						if(cav.nctet >= cav.lctet.size()) return 4;
						cav.lctet[cav.nctet] = ielem;
						cav.nctet++;
						msh.tet2tag(ITAG,ielem) = cavtag;
					}
	
					for(int ii = 0; ii < nbfac; ii++){
						int iface = lbfac[ii];
						if(msh.fac2tag(ITAG,iface) == cavtag) continue;
						// An element in the ball does not belong to the cavity
						if(!opts.allow_topological_correction) return 3;
						if(cav.ncfac >= cav.lcfac.size()) return 4;
						cav.lcfac[cav.ncfac] = iface;
						cav.ncfac++;
						msh.fac2tag(ITAG,iface) = cavtag;
					}
				}else{
					int ifac0 = getpoifac(msh,ipoin);
          METRIS_ASSERT_MSG(!isdeadent(ifac0,msh.fac2poi),"getpoifac returns invalid face for ipoin = "
            <<ipoin<<" get iface (dead) = "<<ifac0);
					if(iverb >= METRIS_CAV_PRTLEV)printf(" -> found (possibly not a) triangle %d\n",ifac0);
					if(ifac0 > 0){
						int iopen, nball;
            bool imani;
            intAr1 dum(0);
						ierro = ball2(msh,ipoin,ifac0,&nball,lbfac,
                                 NULL, dum, 
                                 &iopen,&imani,ITAG2);
            METRIS_ASSERT(ierro == 0);
						for(int ii = 0; ii < nball; ii++){
							int iface = lbfac[ii];
							if(msh.fac2tag(ITAG,iface) == cavtag) continue;
							// An element in the ball does not belong to the cavity
							if(!opts.allow_topological_correction) return 3;
							if(cav.ncfac >= cav.lcfac.size()) return 4;
							cav.lcfac[cav.ncfac] = iface;
							cav.ncfac++;
							msh.fac2tag(ITAG,iface) = cavtag;
						}
					}
				}

				// Is this point a corner? If so: check option. 
				int ibpoi = msh.poi2bpo[ipoin];
				if(ibpoi < 0) continue;
				METRIS_ASSERT("Valid ibpoi " && ibpoi < msh.nbpoi); 
				METRIS_ASSERT("ibpoi pts back to ipoin" && msh.bpo2ibi(ibpoi,0) == ipoin);
				METRIS_ASSERT("Correct poi type" && msh.bpo2ibi(ibpoi,1) <= 1); // Edge or corner
				// If corner:
				if(msh.bpo2ibi(ibpoi,1) == 0){
					ncorn ++;
					// If this is the second corner we see, return error. 
					if(ncorn > 1){
						if(iverb >= 1){
							printf("## ERROR SECOND CORNER SEEN, PREVIOUS = %d CURRENT %d\n",cav.iremcor,ipoin);
						}
						return 8;
					}
					// If corner and corners not allowed to be removed, return. 
					if(!opts.allow_remove_corners && ipoin != cav.ipins) return 6; 
					// If corner and point to insert isn't one, return error. 
					// This would alter the topology of the surface
					// This is safe as we checked in the beginning >= 0
					cav.iremcor = ipoin;
				}
			}
		}
	}

	/*
	Edges STEP 2: check shells in the cavity. 
	Note that point balls include the shells of neighbouring edges. 
	As such, we only check those edges with no collapsed points. 
	*/
	for(int iedgl = 0; iedgl < cav.ncedg; iedgl++){
		int iedge = cav.lcedg[iedgl];
		METRIS_ASSERT(iedge >= 0 && iedge < msh.nedge);
		int ip1 = msh.edg2poi(iedge,0);
		if(msh.poi2tag(ITAG,ip1) == msh.tag[ITAG]) continue;
		int ip2 = msh.edg2poi(iedge,1);
		if(msh.poi2tag(ITAG,ip2) == msh.tag[ITAG]) continue;

		// This edge has non-collapsed points. We need both tet and triangle
		// shell. 
		// We can't use point-bound information because it might send us 
		// away from the edge. Instead, we use directly edg2fac and 
		// fac2tet. 


		// 1. Fetch triangles 
		int iface = msh.edg2fac[iedge];
		// Case of an edge not attached to triangles. 
		// a fortiori not to tetrahedra either, so just skip. 
		if(iface < 0) continue;

		if(msh.fac2tag(ITAG,iface) != cavtag){
			// Face not in cavity. 
			if(!opts.allow_topological_correction) return 3;
			cav.lcfac[cav.ncfac] = iface;
			cav.ncfac++;
			msh.fac2tag(ITAG,iface)   = cavtag;
		}
		int ied = getedgfac(msh,iface,ip1,ip2);
		assert(ied >= 0);

    int ifac2 = msh.fac2fac(iface,ied);
    if(ifac2 >= 0){
      if(msh.fac2tag(ITAG,ifac2) == cavtag) continue;
      if(!opts.allow_topological_correction) return 3;
      cav.lcfac[cav.ncfac] = ifac2;
      cav.ncfac++;
      msh.fac2tag(ITAG,ifac2)   = cavtag;
    }else if(ifac2 < -1){
      ifac2 = iface;
      while(getnextfacnm(msh,iface,ip1,ip2,&ifac2,&ied)){
        if(msh.fac2tag(ITAG,ifac2) == cavtag) continue;
        // Face not in cavity. 
        if(!opts.allow_topological_correction) return 3;
        cav.lcfac[cav.ncfac] = ifac2;
        cav.ncfac++;
        msh.fac2tag(ITAG,ifac2)   = cavtag;
      }
    }

		if(msh.nelem <= 0) continue; 

		// 2. Fetch tetrahedra
		int ielem = msh.fac2tet(iface,0);
		// Triangle is not attached to any tets (entry 0 filled first)
		if(ielem < 0) continue;

		int iopen, nball;
		shell3(msh, ip1, ip2, ielem, &nball, lbfac, &iopen);

		for(int ii = 0; ii < nball; ii++){
			// Shell does not store the edge rank.
			int iele1 = lbfac[ii];
			if(msh.tet2tag(ITAG,iele1) == cavtag) continue;
			// Tet not in cavity. 
			if(!opts.allow_topological_correction) return 3;
			cav.lctet[cav.nctet] = iele1;
			cav.nctet++;
			msh.tet2tag(ITAG,iele1)   = cavtag;
		}
	}





/* 
	II. Triangle checks 
	All tetrahedra attached to cavity triangles must belong to the cavity.
	A number of these are already accounted for from the previous step:
	 triangles attached to a cavity edge can be skipped. 	
  
*/
	for(int ielel = 0; ielel < cav.nctet; ielel++){
		int ielem = cav.lctet[ielel]; 
		if(ielem < 0 || ielem > msh.nelem) return 1;
		if(isdeadent(ielem,msh.tet2poi)){
			if(iverb >= 1)
				printf("Dead tetrahedron in cavity %d nctet = %d msh.nelem = %d \n",
								ielem,cav.nctet,msh.nelem);
			return 2;
		}
		msh.tet2tag(ITAG,ielem) = msh.tag[ITAG];
	}

	for(int ifacl = 0; ifacl < cav.ncfac; ifacl++){
		int iface = cav.lcfac[ifacl]; 
		int ip1 = msh.fac2poi(iface,0);
		int ip2 = msh.fac2poi(iface,1);
		int ip3 = msh.fac2poi(iface,2);

		if(msh.isboundary_faces())METRIS_ASSERT(
			msh.poi2bpo[ip1] >= 0 && msh.poi2bpo[ip2] >= 0 && msh.poi2bpo[ip3] >= 0);

		if(msh.nelem <= 0) continue;

		for(int i=0;i<2;i++){
			int ielem = msh.fac2tet(iface,i);
			if(ielem < 0) continue;
			if(msh.tet2tag(ITAG,ielem) < msh.tag[ITAG]){
				// This face is not in the cavity. 
				if(!opts.allow_topological_correction) return 3;
				if(cav.nctet >= cav.lctet.size())  	   return 4;
				cav.lctet[cav.nctet] = ielem;
				msh.tet2tag(ITAG,ielem) = msh.tag[ITAG];
				cav.nctet++;
			}
		}
	}

	if(iverb >= METRIS_CAV_PRTLEV) printf("-- Cavity to remove %d points \n",cav.nrempts);

	return 0;
}
#endif

template int 
check_cavity_topo<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, 
 MshCavity &cav, CavOprOpt &opts,int ithread);
template int 
check_cavity_topo<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, 
 MshCavity &cav, CavOprOpt &opts,int ithread);


} // End namespace
