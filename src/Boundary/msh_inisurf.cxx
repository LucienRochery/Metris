//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "Boundary/msh_inisurf.hxx"

#include "../Mesh/MeshBase.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../types.hxx"
#include "../ho_constants.hxx"
#include "../aux_topo.hxx"
#include "../low_topo.hxx"
#include "../low_geo/misc.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"
#include "../utils/CT_loop.hxx"
#include "../io_libmeshb.hxx"

#include <tuple>


namespace Metris{

int getNumCorners(const MeshBase &msh){
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
  GETVDEPTH(msh.param);
	if(!msh.CAD()) METRIS_THROW_MSG(
		"EMPTY EGADS CONTEXT");

	if(msh.CAD.EGADS_model == NULL) METRIS_THROW_MSG(
		"EMPTY EGADS MODEL !");

  const int ithrd = 0;

  if(onlyproj){
    CPRINTF2("-- prjMeshPoints start: project {}\n",msh.nbpoi-nbpo0);
  }else{
    CPRINTF2("-- prjMeshPoints start: verify {}, project {}\n",nbpo0,msh.nbpoi-nbpo0);
  }


	double errl2[3] = {0}; 
	double errli[3] = {-1.0};
  int imax[3] = {0};
	int nent[3] = {0};
  int nerr[3] = {0};
  int ierro;
	ego obj;
	double result[18];

  msh.tag[ithrd]++;
  int btag = msh.tag[ithrd]; // So we can move tag later for faces

	if(onlyproj) goto doproj;


	for(int ibpoi = 0; ibpoi < nbpo0; ibpoi++){
    INCVDEPTH(msh.param);
    msh.bpo2tag(ithrd,ibpoi) = btag;
    int ipoin = msh.bpo2ibi(ibpoi,0);
		int ientt = msh.bpo2ibi(ibpoi,2);
    METRIS_ASSERT_MSG(ientt >= 0," ientt < 0: ipoin = {} ientt {} ibpoi {}",
      ipoin,ientt,ibpoi);
		int bdim = msh.bpo2ibi(ibpoi,1);
		METRIS_ASSERT(bdim >= 0 && bdim <= 2 && "bdim within bounds");

		int iref = bdim == 0 ? ientt :
               bdim == 1 ? msh.edg2ref[ientt]
                         : msh.fac2ref[ientt];

    obj = bdim == 0 ? msh.CAD.cad2nod[iref] :
          bdim == 1 ? msh.CAD.cad2edg[iref]
                    : msh.CAD.cad2fac[iref];
    if(obj == NULL){
      MPRINTF("## NULL OBJ bdim = {} iref = {}  \n",bdim,iref);
      egoAr1& cad2ego =  bdim == 0 ? msh.CAD.cad2nod :
                         bdim == 1 ? msh.CAD.cad2edg
                                   : msh.CAD.cad2fac;
      int ncad = cad2ego.get_n();
      MPRINTF("ncad = {}, cad2edg = {}\n",ncad,cad2ego);
    }
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
    CPRINTF2(" - Tested {} points. Errors:\n",nbpo0);
    if(nent[0] > 0)CPRINTF2("{:8} Corners: {:8.3e} max ({}), {:8.3e} avg L2, {} errs\n",
                            nent[0],errli[0],imax[0],errl2[0],nerr[0]);
    if(nent[1] > 0)CPRINTF2("{:8} Edges  :   {:8.3e} max ({}), {:8.3e} avg L2, {} errs\n",
                            nent[1],errli[1],imax[1],errl2[1],nerr[1]);
    if(nent[2] > 0)CPRINTF2("{:8} Faces  :   {:8.3e} max ({}), {:8.3e} avg L2, {} errs\n",
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
  // Corners and such may be delayed once 
  // Double full loop is overkill but simplest to write and few involved. 
  int ndelay = 0;
  for(int irep = 0; irep < 2; irep++){
    INCVDEPTH(msh.param);
    CPRINTF3(" - prjMeshPoints outer rep {}/2\n",irep+1);
    for(int ibpoi = nbpo0; ibpoi < msh.nbpoi; ibpoi++){
      INCVDEPTH(msh.param);
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

      CPRINTF3(" - try project point {} ibpoi {} ibpo0 {} bdim {}\n",
              ipoin,ibpoi,ibpo0,bdim);

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

          if(DOPRINTS3()){
            CPRINTF3(" - {} -> {} guess ipoin {} ibpo0 {} ientt {} ibpoi {} ipoi2 {} ibpo2 {}\n",
                   pdim,bdim,ipoin,ibpo0,msh.bpo2ibi(ibpoi,2),ibpoi,ipoi2,ibpo2);
            int pdim2 = msh.bpo2ibi(msh.poi2bpo[ipoi2],1);
            if(ibpo2 >= ibpoi && irep == 0){
              CPRINTF3(" -> delay \n");
            }else{
              CPRINTF3(" - ibpo2 dim {} (u,v) = {} {} \n",pdim2,
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
          CPRINTF3(" -> delay");
          if(irep == 1){
            lbad.stack(ibpoi);
            if(DOPRINTS3()) PRINTF(" and stack to lbad\n");
          }else if(DOPRINTS3()){
            PRINTF("\n");
          }
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

      CPRINTF3(" - proj ipoin {} pdim {} bdim {} iref {} ibpoi {} ierro {} coord (u,v)/t = {} {}\n",
               ipoin,pdim,bdim,iref,ibpoi,ierro,msh.bpo2rbi(ibpoi,0),msh.bpo2rbi(ibpoi,1));

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
    CPRINTF2("-- End main loop {} bad points to fix \n",lbad.get_n());

  intAr1 lfacl(100);
  const int nnod2 = msh.nnode(2);
  for(int niter = 0; niter < 100; niter++){
    INCVDEPTH(msh.param);
    
    int ncorr = 0;
    int nbad = lbad.get_n();
    for(int ibad = 0; ibad < nbad; ibad++){
      INCVDEPTH(msh.param);
      int ibpoi = lbad[ibad];

      METRIS_ASSERT(msh.bpo2tag(ithrd,ibpoi) < btag);
      int ipoin = msh.bpo2ibi(ibpoi,0);
      METRIS_ASSERT(msh.poi2ent(ipoin,0) >= 0)

      int ifac0 = msh.bpo2ibi(ibpoi,2);
      int bdim  = msh.bpo2ibi(ibpoi,1);
      // Corners can't end up here. Edges neither for now (we could change this, but doing irep above)
      METRIS_ASSERT_MSG(bdim == 2,"Neither corners nor edge points can end up here");

      int iref = msh.fac2ref[ifac0];
      METRIS_ASSERT(iref >= 0);
      obj = msh.CAD.cad2fac[iref];

      CPRINTF3(" - try correct ibpoi {} ipoin {} ifac0 = {}\n",ibpoi,ipoin,ifac0);

      // Seeing as this is face only case, we're going to go by the point's ball
      // seeded by the ifac0 and that is not allowed to cross edges even if same
      // ref (e.g. periodic)
      lfacl.set_n(1);
      lfacl[0] = ifac0;
      msh.tag[ithrd]++;
      bool solvedpt = false;
      for(int ii = 0; ii < lfacl.get_n(); ii++){
        INCVDEPTH(msh.param);
        int iface = lfacl[ii];

        int iver = -1;
        for(int inode = 0; inode < nnod2; inode++){
          int ipoi2 = msh.fac2poi(iface,inode);
          if(ipoi2 == ipoin){
            iver = inode;
            continue;
          }


          int ibpo2 = msh.poi2ebp(ipoi2,2,iface,-1);
          METRIS_ASSERT(ibpo2 >= 0);
          CPRINTF3(" - try point seed {} ipoi2 with ibpo2 = {} tag {} <? {}\n",
                  ipoi2,ibpo2,msh.bpo2tag(ithrd,ibpo2),btag);

          if(msh.bpo2tag(ithrd,ibpo2) < btag) continue;
          // Found one ! 
          CPRINTF3(" - ibpoi {} ipoin {} ifac0 {} using guess face {} ibpo2 {} ipoi2 {} (u,v) = {} {} \n",
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
            CPRINTF3(" -> got (u,v) = {} {} \n",msh.bpo2rbi(ibpoi,0),msh.bpo2rbi(ibpoi,1));
            break;
          }else{
            CPRINTF2("## EG_invEvaluateGuess error {} \n",ierro);
          }
        }// for inode

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

      if(solvedpt){
        lbad[ibad] = -1;
        ncorr++;
      }

    }// for ibad

    int ibad = 0;
    while(ibad < lbad.get_n()){
      int ibpoi = lbad[ibad];
      if(ibpoi >= 0){
        ibad++;
        continue;
      }

      int ibpo1 = lbad.pop();
      // edge case of array ends up empty
      if(ibpo1 >= 0) lbad[ibad] = ibpo1;
    }

    CPRINTF1(" - corr loop {:3} ncorr {}, {} remaining\n",niter,ncorr,lbad.get_n());
    if(lbad.get_n() == 0) break;
  }// for niter

  #if 0
  while(lbad.get_n() > 0){
    INCVDEPTH(msh.param);
    if(niter++ > 100) METRIS_THROW_MSG( 
                         "## Could not fix "<<lbad.get_n()<<" points in CAD proj")

    int ibpoi = lbad.pop();
    METRIS_ASSERT(msh.bpo2tag(ithrd,ibpoi) < btag);
    int ipoin = msh.bpo2ibi(ibpoi,0);
    METRIS_ASSERT(msh.poi2ent(ipoin,0) >= 0)

    int ifac0 = msh.bpo2ibi(ibpoi,2);
    int bdim  = msh.bpo2ibi(ibpoi,1);
    // Corners can't end up here. Edges neither for now (we could change this, but doing irep above)
    METRIS_ASSERT_MSG(bdim == 2,"Neither corners nor edge points can end up here");

    int iref = msh.fac2ref[ifac0];
    METRIS_ASSERT(iref >= 0);
    obj = msh.CAD.cad2fac[iref];

    CPRINTF3(" - try correct ibpoi {} ipoin {} ifac0 = {}\n",ibpoi,ipoin,ifac0);

    // Seeing as this is face only case, we're going to go by the point's ball
    // seeded by the ifac0 and that is not allowed to cross edges even if same
    // ref (e.g. periodic)
    lfacl.set_n(1);
    lfacl[0] = ifac0;
    msh.tag[ithrd]++;
    for(int ii = 0; ii < lfacl.get_n(); ii++){
      INCVDEPTH(msh.param);
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
        CPRINTF3(" - try point seed {} ipoi2 with ibpo2 = {} tag {} <? {}\n",
                ipoi2,ibpo2,msh.bpo2tag(ithrd,ibpo2),btag);

        if(msh.bpo2tag(ithrd,ibpo2) < btag) continue;
        // Found one ! 
        CPRINTF3(" - ibpoi {} ipoin {} ifac0 {} using guess face {} ibpo2 {} ipoi2 {} (u,v) = {} {} \n",
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
          CPRINTF3(" -> got (u,v) = {} {} \n",msh.bpo2rbi(ibpoi,0),msh.bpo2rbi(ibpoi,1));
          break;
        }else{
          CPRINTF2("## EG_invEvaluateGuess error {} \n",ierro);
        }
      }// for inode

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

    if(!solvedpt){

    }


  }
  #endif






	if(msh.nbpoi - nbpo0 > 0){
		if(nent[1] > 0)errl2[1] = sqrt(errl2[1])/nent[1];
		if(nent[2] > 0)errl2[2] = sqrt(errl2[2])/nent[2];
		CPRINTF2("Projected {} points. Errors:\n",msh.nbpoi - nbpo0);
		if(nent[1] > 0)CPRINTF2("{:8} Edges:   {:8.3e} max ({}), {:8.3e} avg L2, {} errs\n",
			                     nent[1],errli[1],imax[1],errl2[1],nerr[1]);
		if(nent[2] > 0)CPRINTF2("{:8} Faces:   {:8.3e} max ({}), {:8.3e} avg L2, {} errs\n",
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
void iniMeshBdryTriangles(MeshBase &msh, HshTab_I3I &intfHshTab){
	if(msh.idim == 2) METRIS_THROW_MSG( "Calling iniMeshBdryTriangles on 2D meshes: NO!");
  GETVDEPTH(msh.param);
	int ncref = 0;
	// In this case, we don't need to check that triangles already exist: save some time.
  // Populate with intfHshTab and update fac2tet. 
  if(msh.nface == 0){
    CPRINTF1("~~ No file-supplied faces; using ininei table\n");
    CPRINTF1("-- New nface = {} \n",msh.nface);
    msh.facHshTab.merge(intfHshTab);
    // This cannot recover interior faces: each face is only seen once.
    for(auto t : msh.facHshTab){
      INCVDEPTH(msh.param);
      int ielem = t.second;
      msh.facHshTab[t.first] = msh.nface;
      CPRINTF1(" - debug ielem {} \n",ielem);
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

	if(ncref > 0) CPRINTF1(" - {} new faces \n",ncref);
}
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template void iniMeshBdryTriangles< n >(MeshBase &msh, HshTab_I3I &intfHshTab);
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
  GETVDEPTH(msh.param);
	int ncree = 0;
	for(int iface = 0; iface < msh.nface; iface++){
    INCVDEPTH(msh.param);
    if(isdeadent(iface,msh.fac2poi)) continue;
		for(int ied = 0; ied < 3; ied++){
			int ifac2 = msh.fac2fac(iface,ied); 
			if(ifac2 >= 0) continue; // Has neighbour, nothing to add
			if(ifac2 < -1) continue; // Non-manifold, already added
			int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
			int ip2 = msh.fac2poi(iface,lnoed2[ied][1]);

			if(getedgglo(msh,ip1,ip2) >= 0) continue; // Already exists
			CPRINTF1("Debug creating edge from iface {}, ied {} neigh = {} ip1, ip2 = {} {} \n",iface,ied,ifac2,ip1,ip2);
      if(ncree == 0 && DOPRINTS2()){
        CPRINTF2(" -> first creation, output mesh \n");
        writeMesh("debug_iniMeshBdryEdges",msh);
      }
			// Create edge
			msh.newedgtopo<ideg>(iface,ied,0);
			ncree ++;
		}
	}

	if(ncree > 0) CPRINTF1(" - {} new edges \n",ncree);
}
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template void iniMeshBdryEdges< n >(MeshBase &msh);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


void iniMeshBdryCorners(MeshBase &msh){
  GETVDEPTH(msh.param);
	int ncrec = 0;
	for(int iedge = 0; iedge < msh.nedge; iedge++){
    INCVDEPTH(msh.param);
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
        METRIS_ENFORCE_MSG(nloop <= 100, "100 times duplicated boundary point = fishy ! cf iniMeshBdryCorners\n");
     	}while(ibpo2 >= 0 && ibpo2 != ibpoi);
     	if(minty == 0) continue;
			}
			// Create corner
			msh.newbpotopo(ipoin,0);
			ncrec ++;
		}
	}

	if(ncrec > 0) CPRINTF1(" - {} new edges \n",ncrec);
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
int iniMeshBdryPoints(MeshBase &msh, int *nbpo0, int ithread){
  GETVDEPTH(msh.param);

  const int ideg = msh.curdeg;

  //if(msh.isboundary_faces() && msh.param->refineConventionsInp)
  //  METRIS_THROW_MSG("TODO: Surface bpois not handled iniMeshBdryPoints "
  //    "with refineConventionsInp == true.");

  intAr1 lrbpo(10);
  int ncor0 = 0, ncor1 = 0, ncrebp = 0;

  intAr1 dum, lbfac(20);
  //MeshArray1D<MeshArray1D<int, METRIS_INT1>, METRIS_INT1> lcofa;
  //lcofa.allocate(4);
  //for(int ii = 0; ii < 4; ii++) lcofa[ii].allocate(20);
  //intAr1 lcof1, lcof2;

  // If in refine convention, begin by propagating negative refs (CAD entity)
  if(!msh.param->refineConventionsInp) goto noRefine;


  if(msh.isboundary_faces()) 
    MPRINTF("\n## WARNING: Faces untested in iniMeshBdryPoints w/ refineConventionsInp\n");

  CPRINTF2(" - Refine convention bpoi intialization\n");
  #ifndef NDEBUG
  // Check all bpois point to refs, not entities. 
  for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
    if(msh.bpo2ibi(ibpoi,1) < 1)continue;
    METRIS_ASSERT(msh.bpo2ibi(ibpoi,2) < 0);
  }
  #endif

  for(int tdim = 1; tdim <= 2; tdim++){
    INCVDEPTH(msh.param);
    if(tdim == 2 && !msh.isboundary_faces()) break;
    int nentt = msh.nentt(tdim);
    const intAr2& ent2poi = msh.ent2poi(tdim);
    const intAr1& ent2ref = msh.ent2ref(tdim);
    const int nnode = getnnode(tdim,ideg);

    for(int ientt = 0; ientt < nentt; ientt++){
      INCVDEPTH(msh.param);
      if(isdeadent(ientt,ent2poi)) continue;
      
      int iref1 = ent2ref[ientt];

      for(int iver = 0; iver < nnode; iver++){
        int ipoin = ent2poi(ientt,iver);
        int ibpo0 = msh.poi2bpo[ipoin];

        METRIS_ENFORCE_MSG(ibpo0 >= 0, "Missing t or (u,v) entries");

        int pdim = msh.bpo2ibi(ibpo0,1);

        lrbpo.set_n(0);
        for(int ibpoi = ibpo0; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
          INCVDEPTH(msh.param);
          int itype = msh.bpo2ibi(ibpoi,1);
          CPRINTF2(" - ibpoi {} bpo tdim {}\n",ibpoi,itype);
          if(itype != tdim) continue;
          int ientt = msh.bpo2ibi(ibpoi,2);
          CPRINTF2(" - ientt {}\n",ientt);
          if(ientt >= 0) continue;
          // In refine convention, the onGeometricEdges entry stores the ref
          // we put here - the entry. 
          int iref2 = - ientt - 1;
          CPRINTF2(" - iref2 {} iref1 {} \n",iref2, iref1);
          if(iref2 != iref1) continue;
          // The ref is correct, but this could still be a loop (one ref, two t's)
          // Simply stack and deal with later
          lrbpo.stack(ibpoi);
        }// for ibpoi

        CPRINTF2(" - tdim {} ientt {} iref {} ipoin {} nrbpo = {} \n",
                 tdim,ientt,iref1,ipoin,lrbpo.get_n());

        if(lrbpo.get_n() == 0) continue;


        if(lrbpo.get_n() == 1){
          int ibpoi = lrbpo[0];
          msh.bpo2ibi(ibpoi,2) = ientt;
          CPRINTF2(" - create link ipoin {} ibpoi {} -> dim {} ent {}\n",
                   ipoin,ibpoi,tdim,ientt);
          if(tdim <= 1) ncor1++;
          //continue;
        }else{
          // case several for one point
          METRIS_ENFORCE_MSG(pdim < tdim, 
            "Point interior to dim {} geom but given several t/(u,v) coordinates",tdim);
        }


        if(tdim == 1){
          // Never too safe
          METRIS_ENFORCE_MSG(iver < 2, "Edge HO node given several t coordinates.");
          METRIS_ENFORCE_MSG(lrbpo.get_n() <= 2,
                  "At most 2 t coords can be given per edge point but {} provided",
                  lrbpo.get_n());
          if(msh.getpoitdim(ipoin) >= 2){
            METRIS_THROW_MSG( "Face point was given edge t coordinate\n");
          }else if(msh.getpoitdim(ipoin) == 1){
            CPRINTF2(" - edge interior point has been updated, continue\n");
            continue;
          }
        }

        if(tdim == 2){
          if(msh.getpoitdim(ipoin) > 2){
            METRIS_THROW_MSG(
            "Interior point was given face (u,v) coordinate\n");
          }else if(msh.getpoitdim(ipoin) == 2){
            CPRINTF2(" - face interior point has been updated, continue\n");
            continue;
          }
        }



        if(tdim == 1){
          // Get the other edge, find which gets which.
          int iedg2 = msh.edg2edg(ientt,1-iver);
          CPRINTF2(" - check neighbour {}\n",iedg2);
          if(iedg2 == -1){
            CPRINTF2(" -> no neighbour, skip\n");
            continue;
          }
          else if(iedg2 < -1){
            bool ifnd = false;
            iedg2 = ientt;
            int inei = 1-iver;
            while(getnextedgnm(msh,ientt,ipoin,&iedg2,&inei)){
              CPRINTF2(" - next edge: {}\n",iedg2);
              METRIS_ASSERT(iedg2 >= 0);
              int iref2 = ent2ref[iedg2];
              if(iref2 == iref1){
                CPRINTF2(" - found edge {} of same ref {}\n",iedg2, iref1);
                ifnd = true;
                break;
              }
            }
            if(!ifnd){
              CPRINTF2(" - node at 3+ node but ref {} is no loop: finish\n",iref1);
              continue;
            }
          }else{
            int iref2 = msh.edg2ref[iedg2];
            CPRINTF2(" - neighbour {} ref {}\n",iedg2,iref2);
            if(iref2 != iref1){
              CPRINTF2(" - node with 2 edges of different refs {} {}: finish\n",
                       iref1,iref2);
              continue;
            }
          }

          // Now ientt has another vertex, as does iedg2, each with a t coordinate
          int ipoi1 = ent2poi(ientt,1-iver);
          int iver2 = msh.getveredg<1>(iedg2,ipoin);
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

            CPRINTF2(" - using ipoi3 = {} ibpo1 = {}\n",ipoi3,ibpo1);

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

            METRIS_ASSERT_MSG(ifnd, "Failed to find ibpoi")
            // The ibpoi here points to the correct entry for ipoi1 or ipoi2
            ibpos[iipoi] = ibpoi;
          }

          CPRINTF2(" - neighbour ibpois: {}, t = {} ; {}, t = {}\n",
                   ibpos[0],msh.bpo2rbi(ibpos[0],0),ibpos[1],msh.bpo2rbi(ibpos[0],0));

          // Now we have the ibpois for ipoi1, ipoi2, corrected.
          // Start by updating for ientt:
          int iused = -1;
          for(int ied = 0; ied < 2; ied++){
            double dst1 = abs(msh.bpo2rbi(lrbpo[0],0) - msh.bpo2rbi(ibpos[ied],0));
            double dst2 = abs(msh.bpo2rbi(lrbpo[1],0) - msh.bpo2rbi(ibpos[ied],0));
            if(dst1 < dst2){
              CPRINTF2(" - t coordinate distances {:10.3e} < {:10.3e} -> update {}\n",
                       dst1,dst2,lrbpo[0]);
              METRIS_ENFORCE_MSG(iused != 0, "t coordinte already used, distances too close?")
              iused = 0;
            }else{
              CPRINTF2(" - t coordinate distances {:10.3e} > {:10.3e} -> update {}\n",
                       dst1,dst2,lrbpo[1]);
              METRIS_ENFORCE_MSG(iused != 1, "t coordinte already used, distances too close?")
              iused = 1;
            }
            int ibpoi = lrbpo[iused];
            msh.bpo2ibi(ibpoi,2) = ied == 0 ? ientt : iedg2;
          }

          ncor0++;

        }else if(tdim == 2){

          METRIS_ASSERT_MSG(lrbpo.get_n() <= 2, "Face meets itself {} at same point", lrbpo.get_n());

          // There can be up to 2 connex components to the point's ball, that
          // have same ref as ientt. One contains ientt. The other, we don't know,
          // and it perhaps doesn't exist. As we fill up ientt's connex component,
          // we will find an element to seed the 2nd connex component iff the 
          // latter exists. 
          int lseed[2] = {ientt ,-1};
          msh.tag[ithread]++;
          for(int icoco = 0; icoco < 2; icoco++){
            int iseed = lseed[icoco];
            if(iseed < 0) break;

            lbfac.set_n(0);
            lbfac.stack(iseed);
            msh.fac2tag(ithread,iseed) = msh.tag[ithread];

            CPRINTF2(" - start new connex component seeded from {}\n",iseed);

            // Do not replace with range-based for loop as lbfac changes.
            for(int ibfac = 0; ibfac < lbfac.get_n(); ibfac++){
              int iface = lbfac[ibfac];
              for(int inei = 0; inei < 3; inei++){
                if(msh.fac2poi(iface,inei) == ipoin) continue;
                int ifac2 = msh.fac2fac(iface, inei);
                // < 0 even if nm same ref is not to share the same (u,v)
                if(ifac2 < 0) continue;
                int iref2 = msh.fac2ref[ifac2];
                if(iref2 != iref1) continue;
                // Cross a connex component
                if(msh.facedg2glo(iface,inei) >= 0){
                  if(iseed == 0 && lseed[1] < 0){
                    CPRINTF2("   - found seed for second ball component {}\n",iface);
                    lseed[1] = ifac2;
                  }
                  continue;
                }
                // Now we are in the same connex component.
                if(msh.fac2tag(ithread,ifac2) >= msh.tag[ithread]) continue;
                msh.fac2tag(ithread,ifac2) = msh.tag[ithread];
                // Same ref face in same connex component, add to ball and update
                //int ibpon = msh.newbpotopo(ipoin,2,ifac2);
                lbfac.stack(ifac2);
              }// for inei
            }// for ibfac

            // See if we can make the update without looking for closest (u,v).
            int irbpo = -1;
            if(icoco == 0){
              if(lseed[1] == -1){
                METRIS_ENFORCE_MSG(lrbpo.get_n() == 1,
                  "{} (u,v)s provided but only one ball connex component found for ipoin {}.",
                  lrbpo.get_n(),ipoin);
                // No other choice
                irbpo = 0;
                CPRINTF2(" - only one connex component -> update using {} = {} {} \n",
                     lrbpo[0], msh.bpo2rbi(lrbpo[0],0), msh.bpo2rbi(lrbpo[0],1));
                goto update_bpo_face;
              }else{
                METRIS_ENFORCE_MSG(lrbpo.get_n() == 2,
                  "Second connex component found but only one (u,v) provided.");
              }
            }else{
              if(lrbpo[0] < 0) irbpo = 1;
              else             irbpo = 0;
              CPRINTF2(" - second connex component -> update using {} = {} {} \n",
                     lrbpo[irbpo], msh.bpo2rbi(lrbpo[irbpo],0), msh.bpo2rbi(lrbpo[irbpo],1));
              goto update_bpo_face;
            }

            // Create new ibpois. First, get average (u,v) from neighbouring 
            // points. This will determine which lrbpo is used for which. 
            {// namespace (goto)
            double uvavg[2] = {0,0};
            int navg = 0;
            double dstmin = 1.0e30, dstmax = -1.0e30;
            for(int iface : lbfac){
              for(int iver = 0; iver < nnode; iver++){
                int ipoi2 = msh.fac2poi(iface,iver);
                int ibpo2 = msh.poi2ebp(ipoi2, 2, iface, iref1);
                if(ibpo2 < 0) continue;
                if(msh.bpo2ibi(ibpo2,2) < 0) continue;
                CPRINTF2("   - point {} attached to face {} has (u,v) = {} {}\n",
                         ipoi2,iface,msh.bpo2rbi(ibpo2,0),msh.bpo2rbi(ibpo2,1));
                uvavg[0] += msh.bpo2rbi(ibpo2,0);
                uvavg[1] += msh.bpo2rbi(ibpo2,1);
                navg++;
              }// for iver
            }// for iface

            METRIS_ENFORCE_MSG(navg > 0, "No nearby (u,v)s were found\n"
              "This error could be handled by delaying a connex component"
              " and proceding by elimination...");
            uvavg[0] /= navg; 
            uvavg[1] /= navg;

            for(int ii : lrbpo){
              int ibpoi = lrbpo[ii];
              METRIS_ASSERT_MSG(ibpoi >= 0, "How are we here with an eliminated rbpoi?");
              double dist = geterrl2<2>(msh.bpo2rbi[ibpoi], uvavg);
              if(dist < dstmin) irbpo = ii;
              dstmin = MIN(dstmin, dist);
              dstmax = MAX(dstmax, dist);
              CPRINTF2(" - lrbpo[{}] = {} (u,v) = {} {} dist = {}\n",
                 ii, ibpoi, msh.bpo2rbi(ibpoi,0), msh.bpo2rbi(ibpoi,1), dist);
            }
            METRIS_ENFORCE_MSG(1 - dstmin/dstmax < 0.1, 
              "Very low (u,v) 'contrast'... max = {} min = {}", dstmax, dstmin);

            }// namespace (goto)


            update_bpo_face:
            METRIS_ASSERT(irbpo >= 0);
            int ibpo_file = lrbpo[irbpo];
            lrbpo[irbpo] = -1;

            // Now we can update the ball connex component bpois. 
            for(int iface : lbfac){
              METRIS_ASSERT(iface >= 0);
              int ibpoi = msh.poi2ebp(ipoin, 2, iface, -1);
              if(ibpoi >= 0){
                METRIS_ASSERT(iface == ientt);
              }else{
                ibpoi = msh.newbpotopo(ipoin, 2, iface);
              }

              for(int ii = 0; ii < nrbi; ii++) 
                msh.bpo2rbi(ibpoi,ii) = msh.bpo2rbi(ibpo_file, ii);
            }

          }// for icoco

        }// if tdim == 1

      } // for iver
    }// for ientt

    if(tdim == 1){
      CPRINTF2(" - VerticesOnGeometricEdges w/ refine convention: " 
               "translated {} open, {} loop bpois \n",ncor1,ncor0);
    }else{
      CPRINTF2(" - VerticesOnGeometricFaces w/ refine convention: " 
               "translated {} bpois created {}\n",ncor1,ncrebp);
    }
  }// for tdim







  noRefine:
  *nbpo0 = msh.nbpoi;

	// Start with edges. Corners are all initialized already. 
  int ntry1 = 0;
	for(int iedge = 0; iedge < msh.nedge; iedge++){
    INCVDEPTH(msh.param);
		if(isdeadent(iedge,msh.edg2poi)) continue;
		for(int iver = 0; iver < getnnod1(ideg); iver++){
      INCVDEPTH(msh.param);
			int ipoin = msh.edg2poi(iedge,iver);
			METRIS_ASSERT(ipoin >= 0 && ipoin < msh.npoin);
			int ibpoi = msh.poi2bpo[ipoin];


			// Create new bpo link either if point new bdry or if corner
			if(ibpoi < 0 || msh.bpo2ibi(ibpoi,1) == 0){
				msh.newbpotopo(ipoin,1,iedge);
				ntry1++;
        CPRINTF3(" - new edge bpo ipoin = {} iedge = {} ntry1 = {}\n",ipoin,iedge,ntry1);
			}
		}
	}
  int ncre1 = msh.nbpoi - *nbpo0;


	// Triangles are only boundary entities in dimension 3+
  int ntry2 = 0;
  int nbpo1 = msh.nbpoi;
	if(msh.idim >= 3){
		// We can now do faces as we needed to know about edge points. 
		for(int iface = 0; iface < msh.nface; iface++){
      INCVDEPTH(msh.param);
			if(isdeadent(iface,msh.fac2poi)) continue;
			for(int iver = 0; iver < getnnod2(ideg); iver++){
        INCVDEPTH(msh.param);
				int ipoin = msh.fac2poi(iface,iver);
				METRIS_ASSERT(ipoin >= 0 && ipoin < msh.npoin);
				int ibpoi = msh.poi2bpo[ipoin];


				// New bpo link if either new or (edge or corner) point. 
				if(ibpoi < 0 || msh.bpo2ibi(ibpoi,1) < 2){
					msh.newbpotopo(ipoin,2,iface);
					ntry2++;
          CPRINTF3(" - new face bpo ipoin = {} iface = {} ntry2 = {}\n",ipoin,iface,ntry2);
					continue;
				}
			}
		}
	}
  int ncre2 = msh.nbpoi - nbpo1;

  CPRINTF2("-- Created {}/{} edge, {}/{} face bpois \n",ncre1,ntry1,ncre2,ntry2);

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

  GETVDEPTH(msh.param);

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
  

  CPRINTF3("-- START genOnGeometricEntLists nbpoi = {}\n",msh.nbpoi);

  for(int ipoin = 0; ipoin < msh.npoin ;ipoin++){
    INCVDEPTH(msh.param);
    if(msh.poi2ent(ipoin,0) < 0){
      if(do_lpoic) lpoic[ipoin] = 0;
      continue;
    }
   if(do_lpoic) lpoic[ipoin] = 0;
   int ibpoi = msh.poi2bpo[ipoin];
   if(ibpoi < 0) continue;
    METRIS_ASSERT_MSG(msh.bpo2ibi(ibpoi,0) == ipoin, 
      "ibpoi mismatch? ipoin = {} ibpoi = {} = {} {} {} {} {}",
      ipoin,ibpoi,
      msh.bpo2ibi(ibpoi,0),msh.bpo2ibi(ibpoi,1),
      msh.bpo2ibi(ibpoi,2),msh.bpo2ibi(ibpoi,3),
      msh.poi2ent(ipoin,0),msh.poi2ent(ipoin,1));
   METRIS_ASSERT(msh.bpo2ibi(ibpoi,0) == ipoin);
   if(ibpoi >= 0 && msh.bpo2ibi(ibpoi,1) == 0){
   	// Rank of the corner (may not be ncorn)
   	if(do_lpoic) lpoic[ipoin] = msh.bpo2ibi(ibpoi,2) + incre; 
      lcorn.stack(ipoin + incre); 
   }
  }

  CPRINTF3(" - nbpoi = {} ncorn = {}\n",msh.nbpoi,lcorn.get_n());
  if(msh.param->refineConventionsOut) goto refineConvention;

  for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
    INCVDEPTH(msh.param);
    int ipoin = msh.bpo2ibi(ibpoi,0);
    if(ipoin < 0) continue;
    if(msh.poi2ent(ipoin,0) < 0) continue;
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
    INCVDEPTH(msh.param);
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
            CPRINTF3(" - skip ibpoi {} t = {}\n",ibpo2,msh.bpo2rbi(ibpo2,1));
            iskip = true;
            break;
          }
        }
      }

      if(iskip) continue;

      CPRINTF3(" - stack ibpoi {} iref {} \n",ibpoi,iref);
      lref.stack(iref+1);
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












