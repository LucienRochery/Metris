//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "Boundary/msh_inisurf.hxx"

#include "../Mesh/MeshBase.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../ho_constants.hxx"
#include "../aux_topo.hxx"
#include "../low_topo.hxx"
#include "../low_geo.hxx"
#include "../mprintf.hxx"
#include "../CT_loop.hxx"
#include "../io_libmeshb.hxx"

#include <tuple>


namespace Metris{

int getNumCorners(MeshBase &msh){
  int ret = 0;
  for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
    ret += msh.bpo2ibi(ibpoi,1) == 0;
  }
  return ret;
}



/*
Points from rank 0 to nbpo0 excluded have been read from file. Their (u,v)s are set. 
Those are only verified. 
Points from nbpo0 included to nbpoi excluded have been re-created. Their (u,v)s will be projected. 
*/
void prjMeshPoints(MeshBase &msh, int nbpo0, bool onlyproj, bool updtX){
  GETVDEPTH(msh);
	if(!msh.CAD()) METRIS_THROW_MSG(TopoExcept(),
		"EMPTY EGADS CONTEXT");

	if(msh.CAD.EGADS_model == NULL) METRIS_THROW_MSG(TopoExcept(),
		"EMPTY EGADS MODEL !");

  const int ithrd = 0;

  if(onlyproj){
    CPRINTF1("-- prjMeshPoints start: project %d\n",msh.nbpoi-nbpo0);
  }else{
    CPRINTF1("-- prjMeshPoints start: verify %d, project %d\n",nbpo0,msh.nbpoi-nbpo0);
  }


	double errl2[3] = {0}; 
	double errli[3] = {-1.0};
  int imax[3] = {0};
	int nent[3] = {0};
  int nerr[3] = {0};
  int ierro;
	ego obj;
	double result[18];


	if(onlyproj) goto doproj;


	for(int ibpoi = 0; ibpoi < nbpo0; ibpoi++){
    INCVDEPTH(msh);
    int ipoin = msh.bpo2ibi(ibpoi,0);
		int ientt = msh.bpo2ibi(ibpoi,2);
    METRIS_ASSERT_MSG(ientt >= 0,"ipoin = "<<ipoin<<" ientt "<<ientt
      <<" ibpoi "<<ibpoi);
		int bdim = msh.bpo2ibi(ibpoi,1);
		METRIS_ASSERT(bdim >= 0 && bdim <= 2 && "bdim within bounds");

		int iref = bdim == 0 ? ientt :
               bdim == 1 ? msh.edg2ref[ientt]
                         : msh.fac2ref[ientt];

    obj = bdim == 0 ? msh.CAD.cad2nod[iref] :
          bdim == 1 ? msh.CAD.cad2edg[iref]
                    : msh.CAD.cad2fac[iref];
		METRIS_ASSERT(obj!=NULL);
		
		if(EG_evaluate(obj, msh.bpo2rbi[ibpoi], result) != 0){
			nerr[bdim]++;
			continue;
		}

    double err;
    if(msh.idim == 2){
      err = geterrl2<2>(msh.coord[ipoin],result);
    }else{
      err = geterrl2<3>(msh.coord[ipoin],result);
    }
		if(err > errli[bdim]){
			imax[bdim]  = ipoin;
			errli[bdim] = err;
		}
		errl2[bdim] += err;
		nent[bdim]++;
		
	}
	if(nbpo0 > 0){
		if(nent[0] > 0)errl2[0] = sqrt(errl2[0])/nent[0];
		if(nent[1] > 0)errl2[1] = sqrt(errl2[1])/nent[1];
		if(nent[2] > 0)errl2[2] = sqrt(errl2[2])/nent[2];
    CPRINTF2(" - Tested %d points. Errors:\n",nbpo0);
    if(nent[0] > 0)CPRINTF2("%8d Corners: %8.3e max (%d), %8.3e avg L2, %d errs\n",
                            nent[0],errli[0],imax[0],errl2[0],nerr[0]);
    if(nent[1] > 0)CPRINTF2("%8d Edges  :   %8.3e max (%d), %8.3e avg L2, %d errs\n",
                            nent[1],errli[1],imax[1],errl2[1],nerr[1]);
    if(nent[2] > 0)CPRINTF2("%8d Faces  :   %8.3e max (%d), %8.3e avg L2, %d errs\n",
                            nent[2],errli[2],imax[2],errl2[2],nerr[2]);
	}


doproj:

	errl2[0]=0.0; errli[0]=-1.0;imax[0]=0;
	errl2[1]=0.0; errli[1]=-1.0;imax[1]=0;
	errl2[2]=0.0; errli[2]=-1.0;imax[2]=0;
	nent[0]=0; nent[1]=0; nent[2]=0;
	nerr[1]=0;nerr[2]=0;

  intAr1 lbad(10);
  lbad.set_n(0);
  msh.tag[ithrd]++;
  int btag = msh.tag[ithrd]; // So we can move tag later for faces
  // Corners and such may be delayed once 
  // Double full loop is overkill but simplest to write and few involved. 
  int ndelay = 0;
  for(int irep = 0; irep < 2; irep++){
  	for(int ibpoi = nbpo0; ibpoi < msh.nbpoi; ibpoi++){
      INCVDEPTH(msh);
      int ipoin = msh.bpo2ibi(ibpoi,0);
      if(msh.poi2ent(ipoin,0) < 0) continue;
      if(msh.bpo2tag(ithrd,ibpoi) >= btag) continue;

  		int ientt = msh.bpo2ibi(ibpoi,2);
  		int bdim  = msh.bpo2ibi(ibpoi,1);
      // Corners don't need projecting 
      if(bdim == 0) continue;
  		METRIS_ASSERT(bdim >= 0 && bdim <= 2);

      int ibpo0 = msh.poi2bpo[ipoin]; 
      int pdim  = msh.bpo2ibi(ibpo0,1);

      //if(pdim == 2 && msh.idim == 3 || pdim == 1 && msh.idim == 2){
      //  // one and done points 
      //  msh.poi2tag(ithrd,ipoin) = msh.tag[ithrd];
      //}

      int iref = bdim == 1 ? msh.edg2ref[ientt] : msh.fac2ref[ientt];
      METRIS_ASSERT(iref >= 0);

      // In this case, we produce a guess to the t or (u,v)
      if(pdim < bdim){
        const int nnode = msh.nnode(bdim);
        const intAr2& ent2poi = msh.ent2poi(bdim);

        // Trying the other vertices, look for one that can provide a (u,v).
        bool delayp = true;
        for(int iver = 0; iver < nnode; iver++){
          int ipoi2 = ent2poi(ientt,iver);
          if(ipoi2 == ipoin) continue;

          int ibpo2 = msh.poi2ebp(ipoi2,bdim,ientt,-1);
          METRIS_ASSERT(ibpo2 >= 0);

          if(DOPRINTS1()){
            CPRINTF1(" - %d -> %d guess ipoin %d ibpo0 %d ientt %d ibpoi %d ipoi2 %d ibpo2 %d\n",
                   pdim,bdim,ipoin,ibpo0,msh.bpo2ibi(ibpoi,2),ibpoi,ipoi2,ibpo2);
            int pdim2 = msh.bpo2ibi(msh.poi2bpo[ipoi2],1);
            if(ibpo2 >= ibpoi && irep == 0){
              CPRINTF1(" -> delay \n");
            }else{
              CPRINTF1(" - ibpo2 dim %d (u,v) = %f %f \n",pdim2,
                        msh.bpo2rbi(ibpo2,0),msh.bpo2rbi(ibpo2,1));
            }
          }

          if(msh.bpo2tag(ithrd,ibpo2) < btag){
            // This point can be revisited on second loop integration, or the 
            // other non-ipoin vertex might give us what we need. 
            continue;
          }

          delayp = false;
          for(int ii = 0; ii < nrbi; ii++) 
            msh.bpo2rbi(ibpoi,ii) = msh.bpo2rbi(ibpo2,ii);
          
        }

        if(delayp){
          CPRINTF1(" -> delay \n");
          if(irep == 1) lbad.stack(ibpoi);
          ndelay++;
          continue;
        }
      }

      msh.bpo2tag(ithrd,ibpoi) = btag;
      
  		double err; 

			obj = bdim == 1 ? msh.CAD.cad2edg[iref] : msh.CAD.cad2fac[iref];

      if(pdim == bdim){
        ierro = EG_invEvaluate(obj, msh.coord[ipoin], msh.bpo2rbi[ibpoi], result);
      }else{
        ierro = EG_invEvaluateGuess(obj, msh.coord[ipoin], msh.bpo2rbi[ibpoi], result);
      }

			if(ierro != 0){
				nerr[bdim]++;
				continue;
			}

      if(msh.idim == 2){
        err = geterrl2<2>(msh.coord[ipoin],result);
      }else{
        err = geterrl2<3>(msh.coord[ipoin],result);
      }


			if(err > errli[bdim]){
				imax[bdim] = ipoin;
				errli[bdim] = err;
			}
			errl2[bdim] += err;
			nent[bdim]++;


      if(updtX){
        for(int ii = 0; ii < msh.idim && updtX; ii++) 
          msh.coord(ipoin,ii) = result[ii];
      }

  	}

    if(ndelay == 0){
      CPRINTF2(" - no delayed points -> break\n");
      break;
    }

  }


  if(lbad.get_n() > 0) 
    CPRINTF1("-- End main loop %d bad points to fix \n",lbad.get_n());

  int niter = 0;
  intAr1 lfacl(100);
  const int nnod2 = msh.nnode(2);
  while(lbad.get_n() > 0){
    INCVDEPTH(msh);
    if(niter++ > 100) METRIS_THROW_MSG(GeomExcept(), 
                         "Could not fix "<<lbad.get_n()<<" points in CAD proj")

    int ibpoi = lbad.pop();
    METRIS_ASSERT(msh.bpo2tag(ithrd,ibpoi) < btag);
    int ipoin = msh.bpo2ibi(ibpoi,0);
    METRIS_ASSERT(msh.poi2ent(ipoin,0) >= 0)

    int ifac0 = msh.bpo2ibi(ibpoi,2);
    int bdim  = msh.bpo2ibi(ibpoi,1);
    // Corners can't end up here. Edges neither for now (we could change this, but doing irep above)
    METRIS_ASSERT_MSG(bdim == 2,"Corners cannot end up here");

    int iref = bdim == 1 ? msh.edg2ref[ifac0] : msh.fac2ref[ifac0];
    METRIS_ASSERT(iref >= 0);
    obj = msh.CAD.cad2fac[iref];

    // Seeing as this is face only case, we're going to go by the point's ball
    // seeded by the ifac0 and that is not allowed to cross edges even if same
    // ref (e.g. periodic)
    lfacl.set_n(1);
    lfacl[0] = ifac0;
    msh.tag[ithrd]++;
    for(int ii = 0; ii < lfacl.get_n(); ii++){
      int iface = lfacl[ii];

      int iver = -1;
      bool solvedpt = false;
      for(int inode = 0; inode < nnod2; inode++){
        int ipoi2 = msh.fac2poi(iface,inode);
        if(ipoi2 == ipoin){
          iver = inode;
          continue;
        }

        int ibpo2 = msh.poi2ebp(ipoi2,2,iface,-1);
        METRIS_ASSERT(ibpo2 >= 0);

        if(msh.bpo2tag(ithrd,ibpo2) < btag) continue;
        // Found one ! 
        CPRINTF1(" - ibpoi %d ipoin %d ifac0 %d using guess face %d ibpo2 %d ipoi2 %d (u,v) = %f %f \n",
                     ibpoi,ipoin,ifac0,iface,ibpo2,ipoi2,msh.bpo2rbi(ibpo2,0),msh.bpo2rbi(ibpo2,1));
        
        for(int ii = 0; ii < nrbi; ii++) 
          msh.bpo2rbi(ibpoi,ii) = msh.bpo2rbi(ibpo2,ii);

        ierro = EG_invEvaluateGuess(obj, msh.coord[ipoin], msh.bpo2rbi[ibpoi], result);

        if(ierro == 0){
          solvedpt = true;
          msh.bpo2tag(ithrd,ibpoi) = btag;
          if(updtX){
            for(int ii = 0; ii < msh.idim && updtX; ii++) 
              msh.coord(ipoin,ii) = result[ii];
          }
          CPRINTF1(" -> got (u,v) = %f %f \n",msh.bpo2rbi(ibpoi,0),msh.bpo2rbi(ibpoi,1));
          break;
        }else{
          CPRINTF1("## EG_invEvaluateGuess error %d \n",ierro);
        }
      }

      if(solvedpt) break;

      METRIS_ASSERT(iver >= 0);

      // We failed to find a good guess. Add neighbours that share ipoin 

      for(int ied = 0; ied < 3; ied++){
        // Only neighbours that contain ipoin
        if(ied == iver) continue;
        int ineil = msh.fac2fac(iface,ied);

        // non-manifold
        if(ineil < 0) continue; 

        // Already seen
        if(msh.fac2tag(ithrd,ineil) >= msh.tag[ithrd]) continue;
        msh.fac2tag(ithrd,ineil) = msh.tag[ithrd];

        // Other surf
        int iref2 = msh.fac2ref[ineil];
        if(iref2 != iref) continue;

        // Sandwiched edge -> periodic case
        int iedgl = msh.fac2edg(iface,ied);
        if(iedgl >= 0) continue;

        // Now ineil is same ref, same left/right neighb of edge if periodic
        lfacl.stack(ineil);
      } // for int ied

    }


  }






	if(msh.nbpoi - nbpo0 > 0){
		if(nent[1] > 0)errl2[1] = sqrt(errl2[1])/nent[1];
		if(nent[2] > 0)errl2[2] = sqrt(errl2[2])/nent[2];
		CPRINTF1("Projected %d points. Errors:\n",msh.nbpoi - nbpo0);
		if(nent[1] > 0)CPRINTF1("%8d Edges:   %8.3e max (%d), %8.3e avg L2, %d errs\n",
			                     nent[1],errli[1],imax[1],errl2[1],nerr[1]);
		if(nent[2] > 0)CPRINTF1("%8d Faces:   %8.3e max (%d), %8.3e avg L2, %d errs\n",
			                     nent[2],errli[2],imax[2],errl2[2],nerr[2]);
	}



}



/*   Helper function for iniMeshNeighbours
If the file supplies no triangles, we have to rebuild:
	- triangles next to tetrahedra without neighbours (boundary)
	- triangles between tetrahedra with different refs
The latter have already been rebuilt in iniMesshNeighbours
If the file supplies triangles, same goes but skip those already in msh.facHshTab. 
*/
template<int ideg>
void iniMeshBdryTriangles(MeshBase &msh, HshTabInt3 &intfHshTab){
	if(msh.idim == 2) METRIS_THROW_MSG(TopoExcept(), "Calling iniMeshBdryTriangles on 2D meshes: NO!");

	int ncref = 0;
	// In this case, we don't need to check that triangles already exist: save some time.
  // Populate with intfHshTab and update fac2tet. 
  if(msh.nface == 0){
    printf("~~ No file-supplied faces; using ininei table\n");
    printf("   -- New nface = %d \n",msh.nface);
    msh.facHshTab.merge(intfHshTab);
    // This cannot recover interior faces: each face is only seen once.
    for(auto t : msh.facHshTab){
      int ielem = t.second;
      msh.facHshTab[t.first] = msh.nface;
      printf("Debug ielem  %d \n",ielem);
      msh.fac2tet[msh.nface][0] = ielem;
      // Correct orienttation by fetching from element. 
      int i1 = std::get<0>(t.first);
      int i2 = std::get<1>(t.first);
      int i3 = std::get<2>(t.first);
      int ifa = getfactet(msh,ielem,i1,i2,i3);
      // Invert from tetrahedron
      msh.newfactopo<ideg>(ielem, ifa, 0);
      ncref ++;
    }
#ifndef NDEBUG 
    printf("Debug print triangles\n");
    msh.fac2poi.print(msh.nface < 10 ? msh.nface : 10);
#endif
  }

  // This is the case where either part of the faces were provided or, more likely, interior faces were
  // not provided but have been reconstructed already. 

  for(auto t : intfHshTab){
  	// Check if already exists, then skip
  	auto s = msh.facHshTab.find(t.first);
  	if(s != msh.facHshTab.end()) continue;

    int ielem = t.second;
    msh.fac2tet[msh.nface][0] = ielem;

    msh.facHshTab[t.first] = msh.nface;

    int i1 = std::get<0>(t.first);
    int i2 = std::get<1>(t.first);
    int i3 = std::get<2>(t.first);
    int ifa = getfactet(msh,ielem,i1,i2,i3);
    msh.newfactopo<ideg>(ielem, ifa, 0);
    ncref ++;
  }

	if(ncref > 0)printf(" - %d new faces \n",ncref);
}
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template void iniMeshBdryTriangles< n >(MeshBase &msh, HshTabInt3 &intfHshTab);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


/*
The main routine has already rebuilt: 
	- edges between different ref triangles
	- non-manifold (triple+) edges 
We now need to rebuild boundary edges for open surfaces (2D).
To handle the nm case, we have not eliminated edges from the hash table as 
new neighbours have been created. 
As such, we have no other choice but to loop over all triangles and create where
neighbour = -1. We don't even need the hash table created by iniMshNeighbours. 
*/
template <int ideg>
void iniMeshBdryEdges(MeshBase &msh){
	// If the mesh is manifold and single ref, we have not created any edges yet. 
	// If this routine becomes a bottleneck in that case, specialize and remove 
	// checks "is in hash table" as for previous. 
	// But I suspect this won't be the case.

	int ncree = 0;
	for(int iface = 0; iface < msh.nface; iface++){
    if(isdeadent(iface,msh.fac2poi)) continue;
		for(int ied = 0; ied < 3; ied++){
			int ifac2 = msh.fac2fac(iface,ied); 
			if(ifac2 >= 0) continue; // Has neighbour, nothing to add
			if(ifac2 < -1) continue; // Non-manifold, already added
			int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
			int ip2 = msh.fac2poi(iface,lnoed2[ied][1]);

			if(getedgglo(msh,ip1,ip2) >= 0) continue; // Already exists
			printf("Debug creating edge from iface %d, ied %d neigh = %d ip1, ip2 = %d %d \n",iface,ied,ifac2,ip1,ip2);
      if(ncree == 0 && msh.param->iverb >= 2){
        printf(" -> first creation, output mesh \n");
        writeMesh("debug_iniMeshBdryEdges",msh);
      }
			// Create edge
			msh.newedgtopo<ideg>(iface,ied,0);
			ncree ++;
		}
	}

	if(ncree > 0)printf(" - %d new edges \n",ncree);
}
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template void iniMeshBdryEdges< n >(MeshBase &msh);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


void iniMeshBdryCorners(MeshBase &msh){
	int ncrec = 0;
	for(int iedge = 0; iedge < msh.nedge; iedge++){
    if(isdeadent(iedge,msh.edg2poi)) continue;
		for(int ive = 0; ive < 2; ive++){
			int iedg2 = msh.edg2edg[iedge][1-ive]; 
			if(iedg2 >= 0) continue; // Has neighbour, nothing to add
			if(iedg2 < -1) continue; // Non-manifold, already added
			int ipoin = msh.edg2poi(iedge,ive);
			int ibpoi = msh.poi2bpo[ipoin];
			if(ibpoi >= 0){
				// Does this already exist as a corner? 
				int ibpo2 = ibpoi;
    		int minty = 3;
    		int nloop = 0;
    		do{
    		  minty = minty < msh.bpo2ibi(ibpo2,1) ? minty : msh.bpo2ibi(ibpo2,1);
    		  ibpo2 = msh.bpo2ibi(ibpo2,3);
    		  if(nloop > 100){
    		  	printf("100 times duplicated boundary point = fishy !\n");
    		  	printf("cf iniMeshBdryCorners\n");
    		  	exit(1);
    		  }
    		}while(ibpo2 >= 0 && ibpo2 != ibpoi);
    		if(minty == 0) continue;
			}
			// Create corner
			msh.newbpotopo<0>(ipoin);
			ncrec ++;
		}
	}

	if(ncrec > 0)printf(" - %d new edges \n",ncrec);
}

/*
	This is called after iniMeshNeighbour
	It creates all the remaining ibpois 
	The file may have supplied most of them through VerticesOnGeometricEdges, etc. 
	But if we're missing that information, we can reconstruct at least the topological link. 
	This MUST NOT create a single corner as that has presumably already been handled. 

  VerticesOnGeometricX entries can be negative if we're in refine convention, 
  in which case they point to CAD entities directly. This is also converted to 
  Metris format ibpois. 

	We have two things to update. Some boundary points do not exist at all. 
	In that case, we create a boundary point and update (edg|fac)2bpo. 

	Some boundary points exist, i.e. poi2bpo[ipoin] >= 0, but the link to the geometric
	edge or face is not set. This is only for corners (edges) or edge points (faces). 

	It is also possible that VerticesOnGeometric(Edges|Triangles) were supplied. 
	In that case, the entries (edg|fac)2bpo have already been initialized. 
	These must be skipped. 
*/
int iniMeshBdryPoints(MeshBase &msh, int ithread){
  GETVDEPTH(msh);

  const int ideg = msh.curdeg;

  if(msh.isboundary_faces() && msh.param->refineConventionsInp)
    METRIS_THROW_MSG(TODOExcept(), "Surface bpois not handled iniMeshBdryPoints "
      "with refineConventionsInp == true.");

  intAr1 lrbpo(10);
  int ncor0 = 0, ncor1 = 0;

  intAr1 dum, lbfac(20);
  intAr1 lcof1(20), lcof2(20);

  // If in refine convention, begin by propagating negative refs (CAD entity)
  if(!msh.param->refineConventionsInp) goto noRefine;


  if(msh.isboundary_faces()) 
    MPRINTF("\n## WARNING: Faces untested in iniMeshBdryPoints w/ refineConventionsInp\n");


  for(int tdim = 1; tdim <= 2; tdim++){
    if(tdim == 2 && !msh.isboundary_faces()) break;
    int nentt = msh.nentt(tdim);
    const intAr2& ent2poi = msh.ent2poi(tdim);
    const intAr1& ent2ref = msh.ent2ref(tdim);
    const int nnode = tdim == 1 ? edgnpps[ideg] : facnpps[ideg];

    for(int ientt = 0; ientt < nentt; ientt++){
      INCVDEPTH(msh);
      if(isdeadent(ientt,ent2poi)) continue;
      
      int iref1 = ent2ref[ientt];

      for(int iver = 0; iver < nnode; iver++){
        int ipoin = ent2poi(ientt,iver);
        int ibpo0 = msh.poi2bpo[ipoin];

        METRIS_ENFORCE_MSG(ibpo0 >= 0, "Missing t or (u,v) entries");

        int pdim = msh.bpo2ibi(ibpo0,1);

        if(ipoin == 2){
          printf("## DEBUG ipoin 2\n");
          for(int ibpoi = ibpo0; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
            printf("ibpoi %d :",ibpoi);
            intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
          }// for ibpoi
        }

        lrbpo.set_n(0);
        for(int ibpoi = ibpo0; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
          int itype = msh.bpo2ibi(ibpoi,1);
          if(itype != tdim) continue;
          int ientt = msh.bpo2ibi(ibpoi,2);
          if(ientt >= 0) continue;
          // In refine convention, the onGeometricEdges entry stores the ref
          // we put here - the entry. 
          int iref2 = - ientt - 1;
          if(iref2 != iref1) continue;
          // The ref is correct, but this could still be a loop (one ref, two t's)
          // Simply stack and deal with later
          lrbpo.stack(ibpoi);
        }// for ibpoi



        if(lrbpo.get_n() == 0) continue;

        CPRINTF1(" - tdim %d ientt %d ipoin %d nrbpo = %d \n",
                 tdim,ientt,ipoin,lrbpo.get_n());

        if(lrbpo.get_n() == 1){
          int ibpoi = lrbpo.pop();
          msh.bpo2ibi(ibpoi,2) = ientt;
          CPRINTF2(" - create link ipoin %d ibpoi %d -> dim %d ent %d\n",
                   ipoin,ibpoi,tdim,ientt);
          ncor1 ++;
          continue;
        }

        // case several for one point
        METRIS_ENFORCE_MSG(pdim < tdim, "Point interior to dim "<<tdim<<
          " geom but given several t/(u,v) coordinates");

        // Both of these are edge only
        if(tdim == 1){
          // Never too safe
          METRIS_ENFORCE_MSG(iver < 2, "Edge HO node given several t coordinates.");

          METRIS_ENFORCE_MSG(lrbpo.get_n() == 2,
                  "At most 2 t coords can be given per edge point but "
                  <<lrbpo.get_n()<<" provided");
        }


        if(tdim == 1){
          // Get the other edge, find which gets which.
          int iedg2 = msh.edg2edg(ientt,1-iver);
          if(iedg2 < 0){
            int inei;
            bool ifnd = false;
            while(getnextedgnm(msh,ientt,ipoin,&iedg2,&inei)){
              int iref2 = ent2ref[iedg2];
              if(iref2 == iref1){
                ifnd = true;
                break;
              }
            }
            METRIS_ENFORCE_MSG(ifnd, "Failed to find other same ref edge in 2 t coord case");
          }

          // Now ientt has another vertex, as does iedg2, each with a t coordinate
          int ipoi1 = ent2poi(ientt,1-iver);
          int iver2 = getveredg<1>(iedg2,ent2poi,ipoin);
          METRIS_ASSERT_MSG(iver2 >= 0, "Neighbour does not share vertex");
          int ipoi2 = ent2poi(iedg2,1-iver2);

          // We can't handle these two being corners yet. To handle this, we'll need
          // several passes, using the hopefully non corner other t neighbours. 

          // We might not have fixed these points negative refs yet. 
          int ibpos[2];
          for(int iipoi = 0; iipoi < 2; iipoi++){
            int ipoi3 = iipoi == 0 ? ipoi1 : ipoi2;

            int ibpo1 = msh.poi2bpo[ipoi3];
            METRIS_ENFORCE_MSG(ibpo1 > 0,
              "Only some ts are given but not all: fix input VerticesOnGeometricEdges")
            // For now, don't handle even one being corner
            METRIS_ENFORCE_MSG(msh.bpo2ibi(ibpo1,1) == 1,
              "TODO: handle CAD edges with no interior nodes in iniMeshBdryPoints");

            CPRINTF2(" - using ipoi3 = %d ibpo1 = %d\n",ipoi3,ibpo1);

            int ibpoi; 
            bool ifnd = false;
            for(ibpoi = ibpo1; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
              int itype = msh.bpo2ibi(ibpoi,1);
              if(itype != 1) continue;
              int ientt = msh.bpo2ibi(ibpoi,2);
              if(ientt >= 0){
                METRIS_ASSERT_MSG(ipoi3 == ipoi1 && ientt == ientt
                               || ipoi3 == ipoi2 && ientt == iedg2,
                               "Already fixed ibpoi does not point to correct edge")
                ifnd = true;
                break;
              }
              int iref2 = - ientt - 1;
              if(iref2 != iref1) continue;
              // fix the ref and break
              if(ipoi3 == ipoi1){
                msh.bpo2ibi(ibpoi,2) = ientt;
              }else{
                msh.bpo2ibi(ibpoi,2) = iedg2;
              }
              ifnd = true;
              break;
            }// for ibpoi

            // The ibpoi here points to the correct entry for ipoi1 or ipoi2
            ibpos[iipoi] = ibpoi;
          }

          CPRINTF2(" - neighbour ibpois: %d, t = %f ; %d, t = %f\n",
                   ibpos[0],msh.bpo2rbi(ibpos[0],0),ibpos[1],msh.bpo2rbi(ibpos[0],0));

          // Now we have the ibpois for ipoi1, ipoi2, corrected.
          // Start by updating for ientt:
          int iused = -1;
          for(int ied = 0; ied < 2; ied++){
            double dst1 = abs(msh.bpo2rbi(lrbpo[0],0) - msh.bpo2rbi(ibpos[ied],0));
            double dst2 = abs(msh.bpo2rbi(lrbpo[1],0) - msh.bpo2rbi(ibpos[ied],0));
            if(dst1 < dst2){
              CPRINTF2(" - t coordinate distances %10.3e < %10.3e -> update %d\n",
                       dst1,dst2,lrbpo[0]);
              METRIS_ENFORCE_MSG(iused != 0, "t coordinte already used, distances too close?")
              iused = 0;
            }else{
              CPRINTF2(" - t coordinate distances %10.3e > %10.3e -> update %d\n",
                       dst1,dst2,lrbpo[1]);
              METRIS_ENFORCE_MSG(iused != 1, "t coordinte already used, distances too close?")
              iused = 1;
            }
            int ibpoi = lrbpo[iused];
            msh.bpo2ibi(ibpoi,2) = ied == 0 ? ientt : iedg2;
          }

          ncor0++;
        }else if(tdim == 2){
          // Face case. Similar but over ball connex components. 
          int iopen;
          bool imani;
          int ierro = ball2(msh, ipoin, ientt, lbfac, dum, &iopen, &imani, ithread);
          METRIS_ENFORCE_MSG(ierro == 0, "ball2 failed in inisurf")
          METRIS_ENFORCE_MSG(lbfac.get_n() > 0, "empty ball2");

          // Split the ball by connex components. Also get one interior ibpoi 
          // per connex component for proximity check. 
          lcof1.set_n(0);
          lcof2.set_n(0);

          msh.tag[ithread]++;
          int ntagfa = 0;
          for(int iface : lbfac){
            msh.fac2tag(ithread,iface) = msh.tag[ithread];
            ntagfa++;
          }

          int icoco = 0;
          int ibpoc[2] = {-1,-1};
          while(ntagfa > 0){
            METRIS_ENFORCE_MSG(icoco <= 1, "More than 2 connex components in ball2");
            intAr1 &lcofa = icoco == 0 ? lcof1 : lcof2;

            lcofa.stack(lbfac[0]);
            msh.fac2tag(ithread,lbfac[0]) = msh.tag[ithread] - 1; // Untag
            ntagfa--;

            int nadded = 0;
            int istack = 0;
            do{
              int iface = lcofa[istack];
              istack++;

              // Seize the opportunity to get an interior ibpoi per coco
              if(ibpoc[icoco] < 0){
                for(int ii = 0; ii < 3; ii++){
                  int ipoi2 = msh.fac2poi(iface,ii);
                  int ibpoi = msh.poi2bpo[ipoi2];
                  METRIS_ASSERT(ibpoi >= 0);
                  if(msh.bpo2ibi(ibpoi,1) == 2){
                    ibpoc[icoco] = ibpoi;
                    break;
                  }
                }
              }//if ibpoc 

              // Add people to connex component
              for(int ii = 0; ii < 3; ii++){
                int ifnei = msh.fac2fac(iface,ii);
                if(ifnei < 0) continue; // Non-manifold : other connex component

                // Not in ball
                if(msh.fac2tag(ithread,ifnei) < msh.tag[ithread]) continue;

                int iedge = msh.facedg2glo(iface,ii);
                if(iedge >= 0) continue; // Interior edge: other connex component

                // In ball and same connex component
                lcofa.stack(ifnei);
                msh.fac2tag(ithread,ifnei) = msh.tag[ithread] - 1;
                nadded++;
              }

            }while(nadded > 0);

            CPRINTF1(" - connex component %d has %d faces \n",icoco, lcofa.get_n());
            icoco++;
          } // while ntagfa

          METRIS_ENFORCE(ibpoc[0] >= 0 && (ibpoc[1] >= 0 || lcof2.get_n() == 0))

          int iused = -1;
          for(int icoco = 0; icoco <= 1; icoco++){
            intAr1 &lcofa = icoco == 0 ? lcof1 : lcof2;
            if(lcofa.get_n() == 0) continue;

            int ibpoi = ibpoc[icoco];

            double dst1 = abs(msh.bpo2rbi(lrbpo[0],0) - msh.bpo2rbi(ibpoi,0));
            double dst2 = abs(msh.bpo2rbi(lrbpo[1],0) - msh.bpo2rbi(ibpoi,0));
            if(dst1 < dst2){
              CPRINTF2(" - (u,v) coordinate distances %10.3e < %10.3e -> update %d",
                       dst1,dst2,lrbpo[0]);
              METRIS_ENFORCE_MSG(iused != 0, "t coordinte already used, distances too close?")
              iused = 0;
            }else{
              CPRINTF2(" - (u,v) coordinate distances %10.3e > %10.3e -> update %d",
                       dst1,dst2,lrbpo[1]);
              METRIS_ENFORCE_MSG(iused != 1, "t coordinte already used, distances too close?")
              iused = 1;
            }

            // Update the negative ibpoi with the seed face 
            msh.bpo2ibi(lrbpo[iused],1) = ientt;

            // Create new ibpois for the other faces
            int nneg = 0;
            for(int iface : lcofa){
              if(iface == ientt) continue;
              msh.newbpotopo<2>(ipoin, iface);
            }// for iface 

          }// for icoco

        }// if tdim == 1

      } // for iver
    }// for ientt

    if(tdim == 1){
      CPRINTF1(" - VerticesOnGeometricEdges w/ refine convention: " 
               "translated %d open, %d loop bpois \n",ncor1,ncor0);
    }else{
      CPRINTF1(" - VerticesOnGeometricFaces w/ refine convention: " 
               "translated %d bpois \n",ncor1);
    }
  }// for tdim







  noRefine:

	// Start with edges. Corners are all initialized already. 
  int ncre1 = 0;
	for(int iedge = 0; iedge < msh.nedge; iedge++){
    INCVDEPTH(msh);
		if(isdeadent(iedge,msh.edg2poi)) continue;
		for(int iver = 0; iver < edgnpps[ideg]; iver++){
			int ipoin = msh.edg2poi(iedge,iver);
			METRIS_ASSERT(ipoin >= 0 && ipoin < msh.npoin);
			int ibpoi = msh.poi2bpo[ipoin];


			// Create new bpo link either if point new bdry or if corner
			if(ibpoi < 0 || msh.bpo2ibi(ibpoi,1) == 0){
				msh.newbpotopo<1>(ipoin,iedge);
				ncre1++;
        CPRINTF2(" - new edge bpo ipoin = %d iedge = %d ncre1 = %d\n",ipoin,iedge,ncre1);
			}
		}
	}

	// Triangles are only boundary entities in dimension 3+
  int ncre2 = 0;
	if(msh.idim >= 3){
		// We can now do faces as we needed to know about edge points. 
		for(int iface = 0; iface < msh.nface; iface++){
      INCVDEPTH(msh);
			if(isdeadent(iface,msh.fac2poi)) continue;
			for(int iver = 0; iver < facnpps[ideg]; iver++){
				int ipoin = msh.fac2poi(iface,iver);
				METRIS_ASSERT(ipoin >= 0 && ipoin < msh.npoin);
				int ibpoi = msh.poi2bpo[ipoin];


				// New bpo link if either new or (edge or corner) point. 
				if(ibpoi < 0 || msh.bpo2ibi(ibpoi,1) < 2){
					msh.newbpotopo<2>(ipoin,iface);
					ncre2++;
          CPRINTF2(" - new face bpo ipoin = %d iface = %d ncre2 = %d\n",ipoin,iface,ncre2);
					continue;
				}
			}
		}
	}

  CPRINTF1("-- Created %d edge, %d face bpois \n",ncre1,ncre2);

	return ncre1 + ncre2;
}

///*
//Helper function for writeMesh
//Returns lpoin and lbpoi to be written as VerticesOnGeometricVertices
//*/
//int genCornerList(Mesh &msh, int mcorn, int *lpoin, int *lbpoi, int offs){
//	// Corners are unique because they're always the lowest-dim representation of a point.
//	// -> no additional checks needed
//	int ncorn = 0;
//
//	for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
//		if(msh.bpo2ibi(ibpoi,1) > 0) continue;
//		if(ncorn >= mcorn){
//			printf("## INCREASE MCORN genCornerList\n");
//			exit(1);
//		}
//		lbpoi[ncorn] = ibpoi + offs;
//		lpoin[ncorn] = msh.bpo2ibi(ibpoi,0) + offs;
//		ncorn++;
//	}
//
//	return ncorn;
//}
//
///*
//Helper function for writeMesh
//For compatibility with Vizir with Corners keyword
//Fills lbpoi of size npoin with 0 where not corner, ref of corner otherwise
//*/
//void genCornerIdx(Mesh &msh, int *lbpoi){
//	for(int ipoin=0;ipoin<msh.npoin;ipoin++){
//		lbpoin[ipoin] = 0;
//		if(msh.poi2bpo[ipoin] < 0) continue;
//		int ibpoi = msh.poi2bpo[ipoin];
//		if(msh.bpo2ibi(ibpoi,1) > 0) continue;
//		lbpoin[ipoin] = msh.bpo2ibi(ibpoi,2);
//	}
//}


/*
Generate lists to write VerticesOnGeometricVertices/Edges/Triangles cf libmeshb
- corn: corners/geometric nodes
- gpoe: geometric points on edges
- gpof: geometric points on faces
*/
void genOnGeometricEntLists(const MeshBase &msh, intAr1& lcorn, intAr1& lpoic,
	                                               intAr2& lgpoe, dblAr2& rgpoe,
	                                               intAr2& lgpof, dblAr2& rgpof,
                                                 int incre){

  GETVDEPTH(msh);

  METRIS_ASSERT(lgpoe.get_stride() == 2);
  METRIS_ASSERT(lgpof.get_stride() == 2);//not a typo

  METRIS_ASSERT(rgpoe.get_stride() == 2);
  METRIS_ASSERT(rgpof.get_stride() == 3);

  bool do_lpoic = lpoic.size1() > 0;

  if(do_lpoic) lpoic.set_n(msh.npoin);
  lcorn.set_n(0); 
  lgpoe.set_n(0); 
  rgpoe.set_n(0); 
  lgpof.set_n(0); 
  rgpof.set_n(0); 
  

  CPRINTF2("-- START genOnGeometricEntLists nbpoi = %d\n",msh.nbpoi);

  for(int ipoin = 0; ipoin < msh.npoin ;ipoin++){
    INCVDEPTH(msh);
    if(msh.poi2ent(ipoin,0) < 0){
      if(do_lpoic) lpoic[ipoin] = 0;
      continue;
    }
  	if(do_lpoic) lpoic[ipoin] = 0;
  	int ibpoi = msh.poi2bpo[ipoin];
  	if(ibpoi < 0) continue;
    METRIS_ASSERT_MSG(msh.bpo2ibi(ibpoi,0) == ipoin, 
      "ibpoi mismatch? ipoin = "<<ipoin<<" ibpoi = "<<ibpoi<<" = "<<
      msh.bpo2ibi(ibpoi,0)<<" "<<msh.bpo2ibi(ibpoi,1)<<" "<<
      msh.bpo2ibi(ibpoi,2)<<" "<<msh.bpo2ibi(ibpoi,3)<<" "
      <<"poi2ent = "<<msh.poi2ent(ipoin,0)<<" "<<msh.poi2ent(ipoin,1));
  	METRIS_ASSERT(msh.bpo2ibi(ibpoi,0) == ipoin);
  	if(ibpoi >= 0 && msh.bpo2ibi(ibpoi,1) == 0){
  		// Rank of the corner (may not be ncorn)
  		if(do_lpoic) lpoic[ipoin] = msh.bpo2ibi(ibpoi,2) + incre; 
      lcorn.stack(ipoin + incre); 
  	}
  }

  CPRINTF2(" - nbpoi = %d ncorn = %d\n",msh.nbpoi,lcorn.get_n());
  if(msh.param->refineConventionsOut) goto refineConvention;

  for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
    INCVDEPTH(msh);
  	int ipoin = msh.bpo2ibi(ibpoi,0);
  	if(ipoin < 0) continue;
    if(msh.poi2ent(ipoin,0) < 0) continue;
    if(ipoin >= msh.npoin){
      printf("ipoin = %d >= npoin = %d \n",ipoin,msh.npoin);
      printf("poi2ent = %d tdim %d \n",msh.poi2ent(ipoin,0),msh.poi2ent(ipoin,1));
      printf("ibpoi = %d :\n",ibpoi);
      intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
      printf("ipoin = %d \n",ipoin);
      METRIS_THROW(TopoExcept());
    }
    METRIS_ASSERT(ipoin < msh.npoin);
  	int itype = msh.bpo2ibi(ibpoi,1);
  	if(itype == 2){
      //face
      int ngpof = lgpof.get_n(); 
      lgpof.inc_n();
      rgpof.inc_n();

  		lgpof[ngpof][0] = ipoin + incre;
      lgpof[ngpof][1] = msh.bpo2ibi(ibpoi,2) + incre;

  		rgpof[ngpof][0] = msh.bpo2rbi(ibpoi,0);
  		rgpof[ngpof][1] = msh.bpo2rbi(ibpoi,1);
  		rgpof[ngpof][2] = 0.0; // Placeholder: should be distance to ent
  	}else if(itype == 1){
      //edge
      int ngpoe = lgpoe.get_n(); 
      lgpoe.inc_n();
      rgpoe.inc_n();
      
  		lgpoe[ngpoe][0] = ipoin + incre;
      lgpoe[ngpoe][1] = msh.bpo2ibi(ibpoi,2) + incre;

  		rgpoe[ngpoe][0] = msh.bpo2rbi(ibpoi,0);
  		rgpoe[ngpoe][1] = 0.0; // Placeholder: should be distance to ent
  	}
  }
  return;

  // With refine convention, we need to eliminate duplicates. 
  refineConvention:
  intAr1 lebpo(10), leref(10);
  intAr1 lfbpo(10), lfref(10);
  // We'll be looping over the chained list for each point
  for(int ibpo0 = 0; ibpo0 < msh.nbpoi; ibpo0++){ 
    INCVDEPTH(msh);
    int ipoin = msh.bpo2ibi(ibpo0,0);
    if(ipoin < 0) continue;

    if(msh.poi2bpo[ipoin] != ibpo0) continue;

    lebpo.set_n(0);
    leref.set_n(0);

    lfbpo.set_n(0);
    lfref.set_n(0);
    for(int ibpoi = ibpo0; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
      int itype = msh.bpo2ibi(ibpoi,1);
      if(itype < 1) continue;
      int ientt = msh.bpo2ibi(ibpoi,2);
      int iref = itype == 1 ? msh.edg2ref[ientt] : msh.fac2ref[ientt];
      intAr1 &lbpo = itype == 1 ? lebpo : lfbpo;
      intAr1 &lref = itype == 1 ? leref : lfref;

      // If one of same ref and close t/(u,v), skip
      bool iskip = false;
      for(int ii = 0; ii < lref.get_n(); ii++){
        int iref2 = lref[ii];
        // If found of same ref, check not same t. 
        if(iref2 == iref){
          int ibpo2 = lbpo[ii];
          double dist = abs(msh.bpo2rbi(ibpo2,0) - msh.bpo2rbi(ibpoi,0));
          if(itype == 2) dist += abs(msh.bpo2rbi(ibpo2,1) - msh.bpo2rbi(ibpoi,1));
          if(dist < Constants::CADparamTol){
            CPRINTF2(" - skip ibpoi %d t = %f\n",ibpo2,msh.bpo2rbi(ibpo2,1));
            iskip = true;
            break;
          }
        }
      }

      if(iskip) continue;

      CPRINTF2(" - stack ibpoi %d iref %d \n",ibpoi,iref);
      lref.stack(iref);
      lbpo.stack(ibpoi);

    }// for ibpoi


    // Created stacks with unique bpos: add to the global lists
    for(int ii = 0; ii < lebpo.get_n(); ii++){
      int ibpoi = lebpo[ii];
      int iref  = leref[ii];

      METRIS_ASSERT_MSG(ibpoi >= 0 && ibpoi < msh.nbpoi,
        "invalid ibpoi in refineConventionsOut genOnGeometricEntLists edge case");

      int ngpoe = lgpoe.get_n(); 
      lgpoe.inc_n();
      rgpoe.inc_n();

      lgpoe(ngpoe,0) = ipoin + incre;
      lgpoe(ngpoe,1) = iref; // irefs should not be incremented as they've been in the writer

      rgpoe(ngpoe,0) = msh.bpo2rbi(ibpoi,0);
      rgpoe(ngpoe,1) = 0.0; // Placeholder: should be distance to ent
    }


    for(int ii = 0; ii < lfbpo.get_n(); ii++){
      int ibpoi = lfbpo[ii];
      int iref  = lfref[ii];

      METRIS_ASSERT_MSG(ibpoi >= 0 && ibpoi < msh.nbpoi,
        "invalid ibpoi in refineConventionsOut genOnGeometricEntLists face case");

      int ngpof = lgpof.get_n(); 
      lgpof.inc_n();
      rgpof.inc_n();
      

      lgpof(ngpof,0) = ipoin + incre;
      lgpof(ngpof,1) = iref; // irefs should not be incremented as they've been in the writer

      rgpof(ngpof,0) = msh.bpo2rbi(ibpoi,0);
      rgpof(ngpof,1) = msh.bpo2rbi(ibpoi,1);
      rgpof(ngpof,2) = 0.0; // Placeholder: should be distance to ent
    }
  }


}







} // End namespace












