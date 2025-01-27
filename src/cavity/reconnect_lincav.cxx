//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../cavity/msh_cavity.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../aux_topo.hxx"
#include "../low_geo.hxx"
#include "../mprintf.hxx"


namespace Metris{

template <class MetricFieldType, int ideg>
int reconnect_lincav(Mesh<MetricFieldType> &msh, const MshCavity& cav, CavOprOpt &opts, 
                     double *qumin, int ithread){
  GETVDEPTH(msh);

  const int ncedg = cav.lcedg.get_n();
  
  if(ncedg <= 0) return 0;
  // if(msh.get_tdim() == 1) METRIS_THROW_MSG(TODOExcept(), "Add quality checks in reconnect lincav for tdim = 1 meshes");
  int ierro = CAV_NOERR;

  msh.tag[ithread]++;

  // Tag cavity edges to determine cavity boundary. 
  for(int iedge : cav.lcedg) msh.edg2tag(ithread,iedge) = msh.tag[ithread];




  // Let's count boundary points first, this will inform whether 2 or more new edges -> easy neighbours
  int necbp = 0;
  int inei2 = -1; // Only needed if necbp == 2 
  // This is the neighbour against ipins which we wouldn't see otherwise 

  for(int iedgl = 0; iedgl < ncedg; iedgl++){
    int iedge = cav.lcedg[iedgl];
    for(int ifa = 0; ifa < 2; ifa++){
      int iedg2 = msh.edg2edg(iedge,ifa);
      if(iedg2 >= 0){// Manifold neighbour
        // Belongs to cavity: continue
        if(msh.edg2tag(ithread,iedg2) >= msh.tag[ithread]) continue;
      }
      int ipoin = msh.edg2poi[iedge][1-ifa]; 
      if(ipoin == cav.ipins){
        inei2 = iedg2;
        continue;
      }

      necbp++;
    }
  }

  //if(necbp == 1) {
  //  if(ncedg != 2) METRIS_THROW_MSG(TODOExcept(),"Interpret this case (lazy assert)")
  //}

  if(necbp == 0) METRIS_THROW_MSG(TODOExcept(),"Interpret this case 2 (lazy assert)")

 
  // Tag possible refs for reconnection
  if(msh.isboundary_edges()){
    int ibins = msh.poi2bpo[cav.ipins];
    METRIS_ASSERT(ibins >= 0);

    // Tag possible edge references 
    int ityp = msh.bpo2ibi(ibins,1);

    if(ityp > 1){
      METRIS_THROW_MSG(TopoExcept(), "Non line new pt on line not caught before ??");
    }


    if(ityp == 0){ // Corner may have several
      // Moreover it is necessarily a reinsertion
      // Thus only the inciding edges in the cavity are considered !
      int ibpo2 = ibins;
      int nn = 0;
      while(ibpo2 != -1){
        METRIS_ENFORCE(nn++ <= METRIS_MAX_WHILE);
        int ity2 = msh.bpo2ibi(ibpo2,1);
        if(ity2 == 1){
          int iedge = msh.bpo2ibi(ibpo2,2);
          int iref  = msh.edg2ref[iedge];
          METRIS_ASSERT_MSG(iref >= 0 && iref < msh.CAD.ncaded,
                                                   " lincav init iref = "<<iref);
          // Only if cavity element. 
          if(msh.edg2tag(ithread,iedge) >= msh.tag[ithread]) 
                                  msh.ced2tag(ithread,iref) = msh.tag[ithread];
        }
        ibpo2 = msh.bpo2ibi(ibpo2,3);
      }
    }else if(ityp == 1){
      int iedge = msh.bpo2ibi(ibins,2);
      METRIS_ASSERT(iedge >= 0);
      int iref  = msh.edg2ref[iedge];
      METRIS_ASSERT(iref >= 0 && iref < msh.CAD.ncaded);
      msh.ced2tag(ithread,iref) = msh.tag[ithread]; 
    }
  }


  int nedg0 = msh.nedge;

  for(int iedgl = 0; iedgl < ncedg; iedgl++){
    INCVDEPTH(msh);
    int iedge = cav.lcedg[iedgl];
    int iref = msh.edg2ref[iedge];
    if(msh.isboundary_edges()){
      METRIS_ASSERT(iref >= 0 && iref < msh.CAD.ncaded);
      if(msh.ced2tag(ithread,iref) < msh.tag[ithread]) return CAV_ERR_LINETOPO; 
      // This error corresponds to changing the CAD topology by joining a line 
      // to another when they didn't meet (or not at this node, if one)
      // It can happen VERY easily:
      // Set points 1, 2, 3 forming an L. Point 2 is the corner. 
      // Edge 1: 1-2 is ref 1; edge 2: 2-3 is ref 2
      // Cavity is edge 1. ipins is point 3. 
      // Cavity boundary is points 1, 2 seen from edge 1 of ref 1 -> new edges 
      // of ref 1 are 1-3 and 2-3
      // First, edge 2-3 already existed (and had ref 2) -> this is yet another 
      // error, but we do not check it as it requires accessing the hash table. 
      // I think this "pinching" may perhaps always be caught checking instead
      // that the edge has the wrong ref. 
    }
    double nrm2 = -1;
    for(int ifa = 0; ifa < 2; ifa++){ // It's a facet (point) 
      int iedg2 = msh.edg2edg(iedge,ifa);
      if(iedg2 >= 0){// Manifold neighbour
        // Belongs to cavity: continue
        if(msh.edg2tag(ithread,iedg2) >= msh.tag[ithread]) continue;
      }

      int ipseed = msh.edg2poi(iedge,1-ifa);
      if(ipseed == cav.ipins) continue;

      // Check edge does not already exist
      // This can happen when collapsing one point in 2 + 1 edge
      {
        int itmp = getedgglo(msh,cav.ipins,ipseed);
        if(itmp >= 0){
          if(msh.edg2tag(ithread,itmp) < msh.tag[ithread]){
            ierro = CAV_ERR_DUPEDG;
            goto cleanup; 
          }
        }
      }


      int iedgn = msh.nedge;
      msh.set_nedge(msh.nedge + 1);
      msh.edg2poi(iedgn,0) = cav.ipins;
      msh.edg2poi(iedgn,1) = ipseed;

      CPRINTF1(" - new edge %d = %d %d    \n",iedgn,
                msh.edg2poi(iedgn,0), msh.edg2poi(iedgn,1) );

      if(msh.isboundary_edges()){
        for(int ii = 0; ii < 2; ii++){
          int ip = msh.edg2poi(iedgn,ii);
          int ib = msh.poi2bpo[ip];
          METRIS_ASSERT(ib >= 0);
          if(msh.bpo2ibi(ib,1) == 1){
            int ied = msh.bpo2ibi(ib,2);
            // if edge about to get deleted:
            if(msh.edg2tag(ithread,ied) >= msh.tag[ithread]){ 
              if(msh.bpo2ibi(ib,3) != -1){
                int ib2 = msh.bpo2ibi(ib,3);
                // This would only happen if we just added another edge here:
                if(msh.bpo2ibi(ib2,1) != 1){ 
                  // then create new bpo and copy the old uvs here
                  int ibn = msh.template newbpotopo<1>(ip,iedgn);
                  for(int jj = 0; jj < nrbi; jj++) 
                    msh.bpo2rbi(ibn,jj) = msh.bpo2rbi(ib,jj);
                  CPRINTF2(" - (1) newbpo ip = %d ibn = %d from ib = %d, t = %f\n",
                    ip,ibn,ib,msh.bpo2rbi(ibn,0));
                }
              }else{
                // then create new bpo and copy the old uvs here
                int ibn = msh.template newbpotopo<1>(ip,iedgn);
                for(int jj = 0; jj < nrbi; jj++) 
                  msh.bpo2rbi(ibn,jj) = msh.bpo2rbi(ib,jj);
                CPRINTF2(" - (2) newbpo ip = %d ibn = %d from ib = %d, t = %f\n",
                  ip,ibn,ib,msh.bpo2rbi(ibn,0));
              }
            }
          }else{ // Corner, if face something went wrong 
            METRIS_ASSERT(msh.bpo2ibi(ib,1) == 0);

            //// We're not in the business of adding nodes
            //METRIS_ASSERT(ip == cav.ipins);

            int ibn = msh.template newbpotopo<1>(ip,iedgn);
            CPRINTF2(" - (3) newbpo ip = %d ibn = %d CORNER case\n",ip,ibn);

            // Let's assume most likely, this is not a loop. 
            // Thus count same ref and, if nn == 1, let's go. 
            // Otherwise, loop again (but happens rarely)
            // and do something more fancy
            int nn = 0; // count edges in the mesh of same ref this corner in
            int ib3 = -1;
            for(int ib2 = ib; ib2 >= 0; ib2 = msh.bpo2ibi(ib2,3)){
              if(msh.bpo2ibi(ib2,1) != 1) continue;
              // Skip the one we just created 
              if(msh.bpo2ibi(ib2,2) == iedgn) continue;

              int iref2 = msh.edg2ref[msh.bpo2ibi(ib2,2)];
              METRIS_ASSERT_MSG(msh.bpo2ibi(ib2,2) >= 0 
                             && msh.bpo2ibi(ib2,2) < msh.nedge - 1,
                "ib2 pointing to invalid edge? msh.nedge = "<<msh.nedge
                <<" pts to "<<msh.bpo2ibi(ib2,2)<<" iedgn = "<<iedgn<<"\n"
                <<" ip = "<<ip<<" ib = "<<ib<<" ib2 = "<<ib2);
              METRIS_ASSERT(iref2 >= 0 && iref2 < msh.CAD.ncaded); // Valgrind mostly
              if(iref2 == iref){
                nn++;
                ib3 = ib2;
              }

            }

            CPRINTF2(" - (3) CORNER step 1 ib3 = %d nn = %d\n",ib3, nn);

            if(nn == 1){
              for(int jj = 0; jj < nrbi; jj++) 
                msh.bpo2rbi(ibn,jj) = msh.bpo2rbi(ib3,jj);
              CPRINTF2(" - (3) CORNER (1) update from ib3 %d ip = %d t = %f\n",
                       ib3,msh.bpo2ibi(ib3,0),msh.bpo2rbi(ib3,0));

            }else{ 
              // In this case, there are several ibpois for the same edge ref
              // For each, we start a walk from the pointed entity (edge)
              // and continue walking away from ipins for as long as we 
              // remain in the cavity. 
              // If the last point is ipseed - the other point in this new edge -
              // then that starting ibpoi had the correct t to update 
              // the ibpoi associated to this new edge. 

              bool ifndg = false;
              for(int ib2 = ib; ib2 >= 0; ib2 = msh.bpo2ibi(ib2,3)){
                // Loop only over edge bpos
                if(msh.bpo2ibi(ib2,1) != 1) continue;

                int iedg3 = msh.bpo2ibi(ib2,2);
                int iref2 = msh.edg2ref[iedg3];
                //Only edges on same CAD edge
                if(iref2 != iref) continue;

                // Now the walking
                bool ifnd = false; // have we found ipseed ?
                int ipoip = cav.ipins; // The point we walk away from

                do{
                  // Stop as soon as cavity left (including first step)
                  if(msh.edg2tag(ithread,iedg3) < msh.tag[ithread]) break;

                  int ip1 = msh.edg2poi(iedg3,0);
                  int ip2 = msh.edg2poi(iedg3,1);
                  if(ip1 == ipseed || ip2 == ipseed){ 
                    // we found the other edge point!
                    ifnd = true;
                    break;
                  }

                  if(ip1 == ipoip){ // Run away from ipoip -> neighbour opp to it
                    ipoip = ip2;
                    iedg3 = msh.edg2edg(iedg3,0);
                  }else if(ip2 == ipoip){
                    ipoip = ip1;
                    iedg3 = msh.edg2edg(iedg3,1);
                  }else{
                    METRIS_THROW_MSG(TopoExcept(),
                      "Neighbours didnt share the expected vertex");
                  }

                // THANKFULLY we need only move as long as manifold (otherwise corner)
                }while(iedg3 >= 0); 

                // Great, we found our connex component. 
                // Then this ib2 is the correct one to copy uvs from
                if(ifnd){
                  for(int ii = 0; ii < nrbi; ii++) 
                    msh.bpo2rbi(ibn,ii) = msh.bpo2rbi(ib2,ii);
                  ifndg = true;
                  CPRINTF2(" - (3) CORNER (2) update from ib2 = %d t = %f \n",
                           ib2, msh.bpo2rbi(ib2,0));
                  break;
                }

              }
              if(msh.param->dbgfull) METRIS_ASSERT(ifndg);

            }
          }
        }
      }


      // Check geometric approximation and validity 
      {
        double
        nrm1 = msh.idim == 2 ? geterrl2<2>(msh.coord[cav.ipins], msh.coord[ipseed])
                             : geterrl2<3>(msh.coord[cav.ipins], msh.coord[ipseed]);

        if(nrm1 < Defaults::ltol * Defaults::ltol){
          CPRINTF1(" - small edge len = %f ! return ip1 ip2 = %d %d \n",
                   sqrt(nrm1),cav.ipins,ipseed); 
          return CAV_ERR_FLATEDG;
        }
        nrm1 = 1.0/sqrt(nrm1);
        int ip1 = msh.edg2poi(iedge,ifa);
        int ip2 = msh.edg2poi(iedge,1-ifa);
        if(nrm2 < 0){ 
          nrm2 = msh.idim == 2 ? geterrl2<2>(msh.coord[ip1], msh.coord[ip2])
                               : geterrl2<3>(msh.coord[ip1], msh.coord[ip2]);
          if(nrm2 < Defaults::ltol * Defaults::ltol) 
            METRIS_THROW_MSG(GeomExcept(), "small nrm2 = "<<nrm2);
          nrm2 = 1.0/sqrt(nrm2);
        }

        // Compute dot product and ensure above geotol
        // (ipins - ipseed) . (ip1 - ip2)
        double 
        dtprd = msh.idim == 2 ? getprdl2<2>(msh.coord[cav.ipins], msh.coord[ipseed],
                                            msh.coord[ip1      ], msh.coord[ip2  ])
                              : getprdl2<3>(msh.coord[cav.ipins], msh.coord[ipseed],
                                            msh.coord[ip1      ], msh.coord[ip2  ]);
        CPRINTF2(" - iedge old %d new %d dtprd %f \n",iedge,iedgn,dtprd*nrm1*nrm2);

        if(dtprd*nrm1*nrm2 < 1 - opts.geodev1){
          CPRINTF1(" ## GEODEV REJECT dtprd = %f ipins %d ipseed %d ip1 %d ip2 %d \n",
                   dtprd,cav.ipins, ipseed, ip1, ip2);
          return CAV_ERR_GEODEVLIN;
        }
      }


      msh.edg2ref[iedgn]   = iref;
      msh.edg2edg(iedgn,0) = iedg2; // Neighbour opposite ipins is free (previous against same vertex)

      if(necbp == 1){
        msh.edg2edg(iedgn,1) = inei2;
      }else if(necbp == 2){
        // Two boundary points creating edges looking at each other.
        if(iedgn == nedg0) msh.edg2edg(iedgn,1) = iedgn+1;
        else               msh.edg2edg(iedgn,1) = iedgn-1;
      }else{
        // Non-manifold edge neighbourhood. We'll have them pointing to each other
        // in same sequence as this loop. Don't forget last one points back to first.
        if(iedgn < nedg0 + necbp - 1)
          msh.edg2edg(iedgn,1) = - (iedgn + 1) - 2;
        else
          msh.edg2edg(iedgn,1) = - nedg0 - 2;
      }

      if(msh.nface >= 0) msh.edg2fac[iedgn] = -1;

      for(int ii = 0; ii < METRIS_MAXTAGS; ii++) msh.edg2tag(ii,iedgn) = 0;


      // Next, we create new HO points if needed
      if constexpr(ideg == 1) continue;


      // Now create control points with coords and t CAD coord. 
      // We need to find the ibpoi corresponding to ipins with the expected ref
      // Otherwise we won't have the (right) t coord in the CAD edge. 
      // We don't have to worry about loops so any edge with the correct ref 
      // is good as there can only be one !
      int ibins, ibseed;
      if(msh.isboundary_edges()){
        // If this fails, seed to iedge 
        ibseed = msh.poi2ebp(ipseed,1,iedgn,-1);
        METRIS_ASSERT(ibseed >= 0);

        // ipins is 0-th vertex of iedgn.
        ibins = msh.poi2ebp(cav.ipins,1,iedgn,-1);
        METRIS_ASSERT(ibins >= 0); 
      }// endif(msh.isboundary_edges())
      
      // Create new control points
      for(int j = 0; j < ideg - 1; j++){
        int ipnew = msh.npoin;
        msh.newpoitopo(1,iedgn);
        msh.edg2poi(iedgn,2+j) = ipnew;
        msh.poi2ent(ipnew,0) = iedgn;
        msh.poi2ent(ipnew,1) = 1;
        msh.poi2bpo[ipnew] = -1;
        double t = (1.0 + j)/ideg; 
        for(int k = 0; k < msh.idim; k++){
          msh.coord(ipnew,k) = (1.0-t)*msh.coord(cav.ipins,k)
                             +      t *msh.coord(ipseed   ,k);
        }

        if(msh.isboundary_edges()){
          // New bpo attached to edge. 
          msh.template newbpotopo<1>(ipnew,iedgn);
          for(int ll = 0; ll < nrbi; ll++){
            msh.bpo2rbi(msh.nbpoi-1,ll) = (1 - t) * msh.bpo2rbi(ibins ,ll)
                                             + t  * msh.bpo2rbi(ibseed,ll);
          }
        }
      }
      


    }
  }


  cleanup:
  // actually nothing to do

  return ierro;
}




/*
// Create new edges directly in mesh.
// Fetch edge cavity boundary points. 
// Skip ipins if one of those. 
// If ipins is a corner, arbitrarily many end points possible. 
// Otherwise, only 2 (or 1 if ipins included)

// TODO: for corners, get refs of incident edges and replicate on new edges
template <class MetricFieldType, int ideg>
int reconnect_lincav(Mesh<MetricFieldType> &msh, const MshCavity& cav, CavOprOpt &opts, double *qumin, int ithread){
	if(ncedg <= 0) return 0;

  if(msh.get_tdim() == 1) METRIS_THROW_MSG(TODOExcept(), "Add quality checks in reconnect lincav for tdim = 1 meshes");

	int msh.param->iverb = opts.msh.param->iverb - 1;

	// n edge cavity boundary pts
	int necbp = 0;
	const int mecbp = 10;
	// List edge cavity boundary points
	// List edge cavity boundary neighbours
	// List edge cavity boundary refs: for corners case
	int lecbp[mecbp], lecbn[mecbp], lecbr[mecbp];

	// Fetch boundary of input cavity
	// Tag whether seen once or twice. 
	msh.tag[ithread]+=2;
	for(int iedgl = 0; iedgl < ncedg; iedgl++){
		int iedge = cav.lcedg[iedgl];
		msh.edg2tag(ithread,iedge) = msh.tag[ithread];
		for(int i=0; i < 2;i++){
			int ipoin = msh.edg2poi(iedge,i);
			if(msh.poi2tag(ithread,ipoin) == msh.tag[ithread]-1){
				// Becomes seen twice or more
				msh.poi2tag(ithread,ipoin) = msh.tag[ithread];
			}else if(msh.poi2tag(ithread,ipoin) < msh.tag[ithread]-1){
				// Point hadn't been seen yet. 
				msh.poi2tag(ithread,ipoin) = msh.tag[ithread]-1;
			}
		}
	}

	if(msh.param->iverb >= 1){
		printf("Print line cavity with points. Tagged == tag have been seen twice.\n");
		for(int iedgl = 0; iedgl < ncedg; iedgl++){
			int iedge = cav.lcedg[iedgl];
			for(int i=0;i<2;i++){
				int ipoin = msh.edg2poi(iedge,i);
				printf("Debug iedgl = %d (%d) ip = %d ptag %d tag %d \n",
					iedgl,iedge,ipoin,msh.poi2tag(ithread,ipoin),msh.tag[ithread]);
			}
		}
	}

	// If ipins already belongs, store one edge having it. 
	// If seen by two edges, it isn't on the boundary, set back to -1
//	int ipinsbdry = -1;
	for(int iedgl = 0; iedgl < ncedg; iedgl++){
		int iedge = cav.lcedg[iedgl];
		msh.edg2tag(ithread,iedge) = msh.tag[ithread];
		for(int i=0; i < 2;i++){
			int ipoin = msh.edg2poi(iedge,i);
			if(msh.poi2tag(ithread,ipoin) == msh.tag[ithread]-1){
				// Point was seen only once. This is always a boundary point.
				// Skip if ipins 
				if(ipoin == cav.ipins){
					//if(ipinsbdry == -1)
					//  ipinsbdry = iedge;
					//else
					//	ipinsbdry = -2;
				}else{
					// If the point is not ipins, it is a boundary point we will need to 
					// connect against. 
					// Because ipins may have several adjacent edges of different refs, 
					// we need to keep in memory 
					if(necbp >= mecbp) METRIS_THROW_MSG(SMemExcept(),"Increase mecbp");
					
					if(msh.param->iverb > 1)
						printf("    -- Add edge cav bpo (1) = %d neigh = %d ref = %d \n",
							ipoin,msh.edg2edg[iedge][1-i],msh.edg2ref[iedge]);
					
					lecbp[necbp] = ipoin;
					lecbn[necbp] = msh.edg2edg[iedge][1-i]; // Neighbour adjacent to the pt
					lecbr[necbp] = msh.edg2ref[iedge];
					necbp++;
				}
			}else if(msh.bpo2ibi[msh.poi2bpo[ipoin]][1] == 0 && ipoin != cav.iremcor){
				// Corner, not to be removed (reinserted)
				if(ipoin == cav.ipins){
					//if(ipinsbdry == -1)
					//  ipinsbdry = iedge;
					//else
					//	ipinsbdry = -2;
				}else{
					if(necbp >= mecbp) METRIS_THROW_MSG(SMemExcept(), "Increase mecbp");
					if(msh.param->iverb > 1) 
						printf("    -- Add edge cav bpo (2) = %d neigh = %d ref = %d \n",
							ipoin,msh.edg2edg[iedge][1-i],msh.edg2ref[iedge]);
					lecbp[necbp] = ipoin;
					lecbn[necbp] = msh.edg2edg[iedge][1-i]; // Neighbour adjacent to the pt
					lecbr[necbp] = msh.edg2ref[iedge];
					necbp++;
				}
			}
		}
	}
	// At this stage: 
	// ipins = -1 if never seen 
	//       = -2 if seen twice or more
	//       = iedge if seen once by edge iedge

	// If > 2 boundary points, ipins is a >= triple point thus a corner
	// Idem if 2 bdry points + ipinsbdry seen once
//  if(necbp > 2 || necbp == 2 && ipinsbdry >= 0){
	if(necbp > 2){
		if(msh.param->iverb >= 1) printf("Point to insert is a corner\n");
    int ibpoi = msh.poi2bpo[cav.ipins];
    if(ibpoi >= 0){ // Point already existed, just check the entry is there
      if(msh.bpo2ibi(ibpoi,1) != 0) return 2;
    }else{ // Point did not exist, create new bpoi for it 
      msh.template newbpotopo<0>(cav.ipins, cav.ipins);
    }
	}

	// If only 1 bdry point and ipinsbdry never seen
	//if(necbp == 1 && ipinsbdry == -1){
  //  METRIS_THROW_MSG(TopoExcept(),"## INSUFFICIENT EDGE CAVITY BDRY PTS ! "<<necbp)
	//}

	if(msh.param->iverb >= 1){
		printf("  -- reconnect_lincav nbdry pts = %d :",necbp);
		for(int i=0; i<necbp; i++) printf(" %d ",lecbp[i]);
		printf("\n");
		//if(ipinsbdry>=0){
		//	printf("  + ipins (%d) is an edge cavity bdry point belonging to %d \n",cav.ipins,ipinsbdry);
		//}
	}

	//if(ipinsbdry >= 0 && necbp != 1){
	//	printf("## CASE NOT HANDLED YET !\n");
	//	printf("Ipins on the boundary as edge %d and necbp = %d \n",ipinsbdry,necbp);
	//	// In this case, the point is a corner. Neighbours must be non manifold, see below on creation. 
	//	exit(1);
	//}

 // METRIS_ENFORCE_MSG(necbp == 1,"## Not sure how to handle this yet.\n");

	//int ipinsnei = -2;
	//if(ipinsbdry >= 0){
	//	if(msh.edg2poi(ipinsbdry,0) == cav.ipins){
	//		ipinsnei = msh.edg2edg(ipinsbdry,1);
	//	}else{
	//		ipinsnei = msh.edg2edg(ipinsbdry,0);
	//	}
	//}


	// Each boundary point creates an edge with ipins. 
	int ncree = 0;
	int iedg0 = msh.nedge - 1;
	for(int iecbp = 0; iecbp < necbp; iecbp++){
		if(msh.nedge >= msh.medge){
			return 2;
		}

		msh.edg2poi[msh.nedge][0] = cav.ipins;
		msh.edg2poi[msh.nedge][1] = lecbp[iecbp];
    // Actually edges can't be flat. 
    //if(opts.fast_reject){
    //  bool iflat;
    //  double meas0 = getmeasentP1<gdim>(msh.edg2poi[msh.nedge], msh.coord, opts.vtol, &iflat);
    //}
		msh.edg2ref[msh.nedge]    = lecbr[iecbp];
    msh.edg2fac[msh.nedge]    = -1;
		for(int kk = 0; kk < METRIS_MAXTAGS; kk++) 
			msh.edg2tag[kk][msh.nedge] = 0; 


		// There are two cases for the neighbours:
		// - if necbp == 2, very straightforward. Two new edges and they're each other's neighbour.
		// - Otherwise, ipins will be at the center of a non-manifold edge neighbourhood.
		// Recall neighbours opposite i-th entity

		// This one is independent of the above
		msh.edg2edg[msh.nedge][0] = lecbn[iecbp];

		if(necbp <= 2){
			// Two boundary points creating edges looking at each other.
			if(iecbp == 0){
				msh.edg2edg[msh.nedge][1] = msh.nedge+1;
			}else{
				msh.edg2edg[msh.nedge][1] = msh.nedge-1;
			}
		}else{
			// Non-manifold edge neighbourhood. We'll have them pointing to each other
			// in same sequence as this loop. Don't forget last one points back to first.
			if(iecbp < necbp - 1)
				msh.edg2edg[msh.nedge][1] = - (msh.nedge + 1) - 2;
			else
				msh.edg2edg[msh.nedge][1] = - (    iedg0 + 1) - 2;
		}

		if constexpr(ideg > 1){

			// Now create control points with coords and t CAD coord. 
			// We need to find the ibpoi corresponding to ipins with the expected ref
			// Otherwise we won't have the (right) t coord in the CAD edge. 
			// We don't have to worry about loops so any edge with the correct ref 
			// is good as there can only be one !
			int ibins, ibseed;
			if(msh.isboundary_edges()){
				ibins = getpoiref2edgbpo(msh, lecbr[iecbp], cav.ipins);
				METRIS_ASSERT(ibins >= 0);
				ibseed = getpoiref2edgbpo(msh, lecbr[iecbp], lecbp[iecbp]);
				METRIS_ASSERT(ibseed >= 0);
			}
			
			// Create new control points
			for(int j=0; j < ideg - 1; j++){
				int ipnew = msh.npoin;
				msh.newpoitopo(1,msh.nedge);
				msh.edg2poi[msh.nedge][2+j] = ipnew;
				msh.poi2ent[ipnew] = msh.nedge;
				msh.poi2bpo[ipnew] = -1;
				double t = (1.0 + j)/ideg; 
				for(int k = 0; k < msh.idim; k++){
					msh.coord(ipnew,k) = (1.0-t)*msh.coord[cav.ipins][k]
						                      +      t *msh.coord[lecbp[iecbp]][k];
				}

				if(msh.isboundary_edges()){
					// New bpo attached to edge. 
					msh.template newbpotopo<1>(ipnew,msh.nedge);
					for(int ll = 0; ll < 2; ll++){
						msh.bpo2rbi[msh.nbpoi-1][ll] = (1 - t) * msh.bpo2rbi(ibins,ll)
							                                + t  * msh.bpo2rbi(ibseed,ll);
					}
				}
	
			}

		}
	}
	return 0;
}	*/
// ---------- Forward declarations
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int reconnect_lincav<MetricFieldAnalytical, n >(Mesh<MetricFieldAnalytical> &msh,\
        const MshCavity& cav, CavOprOpt &opts, double * qumin, int ithread);\
template int reconnect_lincav<MetricFieldFE        , n >(Mesh<MetricFieldFE        > &msh,\
        const MshCavity& cav, CavOprOpt &opts, double * qumin, int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


}// End namespace

