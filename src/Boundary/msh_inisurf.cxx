//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "Boundary/msh_inisurf.hxx"

#include "../Mesh/MeshBase.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../ho_constants.hxx"
#include "../aux_topo.hxx"
#include "../low_geo.hxx"
#include "../io_libmeshb.hxx"

#include <tuple>


namespace Metris{

int getNumCorners(MeshBase &msh){
  int ret = 0;
  for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
    ret += msh.bpo2ibi[ibpoi][1] == 0;
  }
  return ret;
}



/*
Points from rank 0 to nbpo0 excluded have been read from file. Their (u,v)s are set. 
Those are only verified. 
Points from nbpo0 included to nbpoi excluded have been re-created. Their (u,v)s will be projected. 
*/
void prjMeshPoints(MeshBase &msh, int nbpo0, bool onlyproj, int typ_proj, bool updtX){
	if(!msh.CAD()) METRIS_THROW_MSG(TopoExcept(),
		"EMPTY EGADS CONTEXT");

	if(msh.CAD.EGADS_model == NULL) METRIS_THROW_MSG(TopoExcept(),
		"EMPTY EGADS MODEL !");

  const int iverb = msh.param->iverb;
  const int ithrd = 0;

  if(iverb >= 1){
    if(onlyproj){
      printf("  prjMeshPoints start: project %d\n",msh.nbpoi-nbpo0);
    }else{
      printf("  prjMeshPoints start: verify %d, project %d\n",nbpo0,msh.nbpoi-nbpo0);
    }
  }


	double errl2_cor=0.0, errli_cor=-1.0;int imax_cor=0;
	double errl2_edg=0.0, errli_edg=-1.0;int imax_edg=0;
	double errl2_fac=0.0, errli_fac=-1.0;int imax_fac=0;

	int n_cor=0, n_edg=0, n_fac=0;
	int ierro;
	int nerr_cor=0,nerr_edg=0,nerr_fac=0;
	ego obj;
	double result[18];


	if(onlyproj) goto doproj;



	for(int ibpoi = 0; ibpoi < nbpo0; ibpoi++){
		int ientt = msh.bpo2ibi[ibpoi][2];
		int ityp = msh.bpo2ibi[ibpoi][1];
		assert(ityp >= 0 && ityp <= 2 && "ityp within bounds");
		int ipoin = msh.bpo2ibi[ibpoi][0];

		if(ityp == 0){
			// Type Node
			int iref = ientt;
			if(iref >= msh.CAD.ncadno) METRIS_THROW_MSG(TopoExcept(),
				"INVALID NODE REFERENCE ! "<<iref<<" >= "<<msh.CAD.ncadno);
			obj = msh.CAD.cad2nod[iref];
			assert(obj!=NULL);
			
			ierro = EG_evaluate(obj, NULL, result);
			if(ierro != 0){
				nerr_cor++;
				continue;
			}
			double err = geterrl2<3>(msh.coord[ipoin],result);
			if(err > errli_cor){
				imax_cor = ipoin;
				errli_cor = err;
			}
			errl2_cor += err;
			n_cor++;
		}else if(ityp == 1){
			// Type Edge
			int iref = msh.edg2ref[ientt];
			if(iref >= msh.CAD.ncaded)METRIS_THROW_MSG(TopoExcept(),
				"INVALID EDGE REFERENCE "<<iref<<" >= "<<msh.CAD.ncaded);

			obj = msh.CAD.cad2edg[iref];
			assert(obj!=NULL);
			
			double t[2];
			t[0] = msh.bpo2rbi[ibpoi][0];
			ierro = EG_evaluate(obj, t, result);
			if(ierro != 0){
				nerr_edg++;
				continue;
			}
      double err;
      if(msh.idim == 2){
        err = geterrl2<2>(msh.coord[ipoin],result);
      }else{
        err = geterrl2<3>(msh.coord[ipoin],result);
      }
			if(err > errli_edg){
				imax_edg = ipoin;
				errli_edg = err;
			}
			errl2_edg += err;
			n_edg++;
		}else if(ityp == 2 && msh.idim >= 3){
			int iref = msh.fac2ref[ientt];
			// Type Face
			if(iref >= msh.CAD.ncadfa)METRIS_THROW_MSG(TopoExcept(),
				"INVALID FACE REFERENCE ! "<<iref<<" >= "<<msh.CAD.ncadfa);

			obj = msh.CAD.cad2fac[iref];
			assert(obj!=NULL);
			
			double uv[2] = {msh.bpo2rbi[ibpoi][0],msh.bpo2rbi[ibpoi][1]};
			ierro = EG_evaluate(obj, uv, result);
			if(ierro != 0){
				nerr_fac++;
				continue;
			}
      double err;
      if(msh.idim == 2){
        err = geterrl2<2>(msh.coord[ipoin],result);
      }else{
        err = geterrl2<3>(msh.coord[ipoin],result);
      }
			if(err > errli_fac){
				imax_fac = ipoin;
				errli_fac = err;
			}
			errl2_fac += err;
			n_fac++;
		}
		//EG_evaluate()
	}
	if(nbpo0 > 0){
		if(n_cor > 0)errl2_cor = sqrt(errl2_cor)/n_cor;
		if(n_edg > 0)errl2_edg = sqrt(errl2_edg)/n_edg;
		if(n_fac > 0)errl2_fac = sqrt(errl2_fac)/n_fac;
    if(iverb >= 1){
      printf("Tested %d points. Errors:\n",nbpo0);
      if(n_cor > 0)printf("%8d Corners: %8.3e max (%d), %8.3e avg L2, %d errs\n",
                             n_cor,errli_cor,imax_cor,errl2_cor,nerr_cor);
      if(n_edg > 0)printf("%8d Edges:   %8.3e max (%d), %8.3e avg L2, %d errs\n",
                             n_edg,errli_edg,imax_edg,errl2_edg,nerr_edg);
      if(n_fac > 0)printf("%8d Faces:   %8.3e max (%d), %8.3e avg L2, %d errs\n",
                             n_fac,errli_fac,imax_fac,errl2_fac,nerr_fac);
    }
	}


doproj:

	errl2_cor=0.0; errli_cor=-1.0;imax_cor=0;
	errl2_edg=0.0; errli_edg=-1.0;imax_edg=0;
	errl2_fac=0.0; errli_fac=-1.0;imax_fac=0;
	n_cor=0; n_edg=0; n_fac=0;
	nerr_cor=0;nerr_edg=0;nerr_fac=0;


  double coop[3];

  msh.tag[ithrd]++;
  // Corners and such may be delayed once 
  // Double full loop is overkill but simplest to write and few involved. 
  int ndelay = 0;
  for(int irep = 0; irep < 2; irep++){

  	for(int ibpoi = nbpo0; ibpoi < msh.nbpoi; ibpoi++){
      int ipoin = msh.bpo2ibi[ibpoi][0];
      if(msh.poi2ent(ipoin,0) < 0) continue;
      if(msh.poi2tag(ithrd,ipoin) >= msh.tag[ithrd]) continue;

  		int ientt = msh.bpo2ibi[ibpoi][2];
  		int ityp  = msh.bpo2ibi[ibpoi][1];
  		METRIS_ASSERT(ityp >= 0 && ityp <= 2);

      int ibpo0 = msh.poi2bpo[ipoin]; 
      int ityp0 = msh.bpo2ibi[ibpo0][1];

      int iref = ityp == 0 ? ientt 
               : ityp == 1 ? msh.edg2ref[ientt] : msh.fac2ref[ientt];

      if(ityp0 < ityp){
        if(ityp == 1){
          // Guess using ientt. 
          int ipoi2 = msh.edg2poi(ientt,0) == ipoin ? 
                      msh.edg2poi(ientt,1) : msh.edg2poi(ientt,0);
          int ibpo2 = msh.poi2bpo[ipoi2]; 
          METRIS_ASSERT(ibpo2 >= 0);

          if(iverb >= 3){
            printf(" - Corner guess ipoin %d ibpo0 %d ibpoi %d ipoi2 %d ibpo2 %d\n",
                   ipoin,ibpo0,ibpoi,ipoi2,ibpo2);
            if(ibpo2 >= ibpoi && irep == 0) printf(" -> delay \n");
            else printf(" - ibpo2 t = %f\n",msh.bpo2rbi(ibpo2,0));
          }

          if(ibpo2 >= ibpoi && irep == 0){
            // This point will need to be revisited on second loop iteration 
            ndelay++;
            continue;
          }

          msh.bpo2rbi(ibpoi,0) = msh.bpo2rbi(ibpo2,0);

        }else{
          METRIS_THROW_MSG(TODOExcept(), "Implement face (u,v) guess");
        }
      }else{
        // one and done points 
        msh.poi2tag(ithrd,ipoin) = msh.tag[ithrd];
      }


      // Only carry out projection if this is the lowest tdim representation of ipoin 
      // if(ityp2 != ityp) continue;
   

      for(int ii = 0; ii < msh.idim; ii++) coop[ii] = msh.coord[ipoin][ii];
      for(int ii = msh.idim; ii < 3; ii++) coop[ii] = 0;
      


  		double err; 
  		if(ityp == 0){
  			// Type Node
        continue;
  		}else if(ityp == 1){
  			// Type Edge
  			if(iref < 0) continue;

  			if(iref >= msh.CAD.ncaded){
  				printf("## INVALID EDGE INDEX ! %d >= %d\n",iref,msh.CAD.ncaded);
  				exit(1);
  			}
  			obj = msh.CAD.cad2edg[iref];
  			
  			//if(nprt1 < 10){
  			//	printf("Debug 10, type 1 bpo2rbi %f \n", msh.bpo2rbi[ibpoi][0]);
  			//} 

  			if(typ_proj == 1){
          if(ityp0 == ityp){
            ierro = EG_invEvaluate(obj, coop, msh.bpo2rbi[ibpoi], result);
          }else{
            if(ipoin == 1){
              printf("## DEBUG IPOIN 1 projection irep %d guess %15.7e\n",irep,
                msh.bpo2rbi(ibpoi,0));
            }
            ierro = EG_invEvaluateGuess(obj, coop, msh.bpo2rbi[ibpoi], result);
            if(ipoin == 1){
              printf("## DEBUG IPOIN 1 projection t after %15.7e\n",
                msh.bpo2rbi(ibpoi,0));
            }
          }
  				//ierro = EG_invEvaluateGuess(obj, msh.coord[ipoin], msh.bpo2rbi[ibpoi], result);
  			}else{
  				ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
  			}
  			if(ierro != 0){
  				nerr_edg++;
  				continue;
  			}
        if(msh.idim == 2){
          err = geterrl2<2>(msh.coord[ipoin],result);
        }else{
          err = geterrl2<3>(msh.coord[ipoin],result);
        }


  			//if(nprt1 < 10){
  			//	printf("Debug 10, type 1 err %10.3e \n", err);
  			//	nprt1++;
  			//}
  			if(err > errli_edg){
  				imax_edg = ipoin;
  				errli_edg = err;
  			}
  			errl2_edg += err;
  			n_edg++;

        //if(ityp2 == ityp){ // Use lowest tdim rep for coord update. We still needed the uvs though.
        //  for(int i = 0; i < msh.idim; i++)
        //    msh.coord[ipoin][i] = result[i];
        //}
  		}else if(ityp == 2 && msh.idim >= 3){
  			// Type Face
  			//if(nprt2 < 10){
  			//	printf("Debug 10, type 2 bpo2rbi %f %f \n", msh.bpo2rbi[ibpoi][0], msh.bpo2rbi[ibpoi][1]);
  			//}
  			if(iref < 0) continue;

  			if(iref >= msh.CAD.ncadfa)METRIS_THROW_MSG(TopoExcept(),
  				"INVALID FACE INDEX ! "<<iref<<" >= "<<msh.CAD.ncadfa);
  				
  			obj = msh.CAD.cad2fac[iref];

  			if(typ_proj == 1){ 
          if(irep == 0){
            ierro = EG_invEvaluate(obj, coop, msh.bpo2rbi[ibpoi], result);
          }else{
            ierro = EG_invEvaluateGuess(obj, coop, msh.bpo2rbi[ibpoi], result);
          }
  				//ierro = EG_invEvaluateGuess(obj, msh.coord[ipoin], msh.bpo2rbi[ibpoi], result);
  			}else{
  				ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
  			}
  			if(ierro != 0){
  				nerr_fac++;
  				continue;
  			}
        if(msh.idim == 2){
          err = geterrl2<2>(msh.coord[ipoin],result);
        }else{
          err = geterrl2<3>(msh.coord[ipoin],result);
        } 
  			//if(nprt2 < 10){
  			//	printf("Debug 10, type 2 err %10.3e \n", err);
  			//	nprt2++;
  			//}
  			if(err > errli_fac){
  				imax_fac = ipoin;
  				errli_fac = err;
  			}
  			errl2_fac += err;
  			n_fac++;

        //if(ityp2 == ityp){
    		//	for(int i = 0; i < msh.idim; i++)
  			//	  msh.coord[ipoin][i] = result[i];
        //}
  		}

      if(updtX){
        for(int i = 0; i < msh.idim; i++) msh.coord[ipoin][i] = result[i];
      }

  	}

    if(ndelay == 0){
      if(iverb >= 1) printf(" - no delayed points -> break\n");
      break;
    }

  }

	if(msh.nbpoi - nbpo0 > 0){
		if(n_edg > 0)errl2_edg = sqrt(errl2_edg)/n_edg;
		if(n_fac > 0)errl2_fac = sqrt(errl2_fac)/n_fac;
    if(iverb >= 1){
  		if(typ_proj == 1) printf("Projected ");
  		else              printf("Evaluated ");
  		printf("%d points. Errors:\n",msh.nbpoi - nbpo0);
  		if(n_edg > 0)printf("%8d Edges:   %8.3e max (%d), %8.3e avg L2, %d errs\n",
  			                     n_edg,errli_edg,imax_edg,errl2_edg,nerr_edg);
  		if(n_fac > 0)printf("%8d Faces:   %8.3e max (%d), %8.3e avg L2, %d errs\n",
  			                     n_fac,errli_fac,imax_fac,errl2_fac,nerr_fac);
    }
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
    		  minty = minty < msh.bpo2ibi[ibpo2][1] ? minty : msh.bpo2ibi[ibpo2][1];
    		  ibpo2 = msh.bpo2ibi[ibpo2][3];
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

	We have two things to update. Some boundary points do not exist at all. 
	In that case, we create a boundary point and update (edg|fac)2bpo. 

	Some boundary points exist, i.e. poi2bpo[ipoin] >= 0, but the link to the geometric
	edge or face is not set. This is only for corners (edges) or edge points (faces). 

	It is also possible that VerticesOnGeometric(Edges|Triangles) were supplied. 
	In that case, the entries (edg|fac)2bpo have already been initialized. 
	These must be skipped. 
*/
template <int ideg>
int iniMeshBdryPoints(MeshBase &msh){

	int ncrea = 0;
  int iverb = msh.param->iverb; 


	// Start with edges. Corners are all initialized already. 
	for(int iedge = 0; iedge < msh.nedge; iedge++){
		if(isdeadent(iedge,msh.edg2poi)) continue;
		for(int iver = 0; iver < edgnpps[ideg]; iver++){
			int ipoin = msh.edg2poi(iedge,iver);
			assert(ipoin >= 0 && ipoin < msh.npoin);
			int ibpoi = msh.poi2bpo[ipoin];


			// Create new bpo link either if point new bdry or if corner
			if(ibpoi < 0 || msh.bpo2ibi[ibpoi][1] == 0){
				msh.newbpotopo<1>(ipoin,iedge);
				ncrea++;
        if(iverb >= 3) printf(" - new edge bpo ipoin = %d iedge = %d ncrea = %d\n",ipoin,iedge,ncrea);
			}
		}
	}

	// Triangles are only boundary entities in dimension 3+
	if(msh.idim >= 3){
		// We can now do faces as we needed to know about edge points. 
		for(int iface = 0; iface < msh.nface; iface++){
			if(isdeadent(iface,msh.fac2poi)) continue;
			for(int iver = 0; iver < facnpps[ideg]; iver++){
				int ipoin = msh.fac2poi(iface,iver);
				assert(ipoin >= 0 && ipoin < msh.npoin);
				int ibpoi = msh.poi2bpo[ipoin];


				// New bpo link if either now or (edge or corner) point. 
				if(ibpoi < 0 || msh.bpo2ibi[ibpoi][1] < 2){
					msh.newbpotopo<2>(ipoin,iface);
					ncrea++;
          if(iverb >= 3) printf(" - new face bpo ipoin = %d iface = %d ncrea = %d\n",ipoin,iface,ncrea);
					continue;
				}
			}
		}
	}

//	print_bpolist(msh,0);

	return ncrea;
}
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int iniMeshBdryPoints< n >(MeshBase &msh);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


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
//		if(msh.bpo2ibi[ibpoi][1] > 0) continue;
//		if(ncorn >= mcorn){
//			printf("## INCREASE MCORN genCornerList\n");
//			exit(1);
//		}
//		lbpoi[ncorn] = ibpoi + offs;
//		lpoin[ncorn] = msh.bpo2ibi[ibpoi][0] + offs;
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
//		if(msh.bpo2ibi[ibpoi][1] > 0) continue;
//		lbpoin[ipoin] = msh.bpo2ibi[ibpoi][2];
//	}
//}


/*
Generate lists to write VerticesOnGeometricVertices/Edges/Triangles cf libmeshb
- corn: corners/geometric nodes
- gpoe: geometric points on edges
- gpof: geometric points on faces
*/
void genOnGeometricEntLists(MeshBase &msh, intAr1& lcorn, intAr1& lpoic,
	                                         intAr2& lgpoe, dblAr2& rgpoe,
	                                         intAr2& lgpof, dblAr2& rgpof,
                                           int incre){

  const int iverb = msh.param->iverb;

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
  
  if(iverb >= 2) printf(" - genOnGeometricEntListsstart nbpoi = %d\n",msh.nbpoi);

  for(int ipoin = 0; ipoin < msh.npoin ;ipoin++){
    if(msh.poi2ent(ipoin,0) < 0){
      if(do_lpoic) lpoic[ipoin] = 0;
      continue;
    }
  	if(do_lpoic) lpoic[ipoin] = 0;
  	int ibpoi = msh.poi2bpo[ipoin];
  	if(ibpoi < 0) continue;
    METRIS_ASSERT_MSG(msh.bpo2ibi[ibpoi][0] == ipoin, 
      "ibpoi mismatch? ipoin = "<<ipoin<<" ibpoi = "<<ibpoi<<" = "<<
      msh.bpo2ibi[ibpoi][0]<<" "<<msh.bpo2ibi[ibpoi][1]<<" "<<
      msh.bpo2ibi[ibpoi][2]<<" "<<msh.bpo2ibi[ibpoi][3]<<" "
      <<"poi2ent = "<<msh.poi2ent(ipoin,0)<<" "<<msh.poi2ent(ipoin,1));
  	METRIS_ASSERT(msh.bpo2ibi[ibpoi][0] == ipoin);
  	if(ibpoi >= 0 && msh.bpo2ibi[ibpoi][1] == 0){
  		// Rank of the corner (may not be ncorn)
  		if(do_lpoic) lpoic[ipoin] = msh.bpo2ibi[ibpoi][2] + incre; 
      lcorn.stack(ipoin + incre); 
  	}
  }

  if(iverb >= 2) printf(" - nbpoi = %d ncorn = %d\n",msh.nbpoi,lcorn.get_n());

  for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
  	int ipoin = msh.bpo2ibi[ibpoi][0];
  	if(ipoin < 0) continue;
    if(msh.poi2ent(ipoin,0) < 0) continue;
    if(ipoin >= msh.npoin){
      printf("ipoin = %d >= npoin = %d \n",ipoin,msh.npoin);
      printf("poi2ent = %d tdim %d \n",msh.poi2ent(ipoin,0),msh.poi2ent(ipoin,1));
      printf("ibpoi = %d :\n",ibpoi);
      intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
      printf("iopin = %d \n",ipoin);
      exit(1);
    }
    METRIS_ASSERT(ipoin < msh.npoin);
  	int itype = msh.bpo2ibi[ibpoi][1];
  	if(itype == 2){
      //face
      int ngpof = lgpof.get_n(); 
      lgpof.inc_n();
      rgpof.inc_n();

  		lgpof[ngpof][0] = ipoin + incre;
  		lgpof[ngpof][1] = msh.bpo2ibi[ibpoi][2] + incre;

  		rgpof[ngpof][0] = msh.bpo2rbi[ibpoi][0];
  		rgpof[ngpof][1] = msh.bpo2rbi[ibpoi][1];
  		rgpof[ngpof][2] = 0.0; // Placeholder: should be distance to ent
  	}else if(itype == 1){
      //edge
      int ngpoe = lgpoe.get_n(); 
      lgpoe.inc_n();
      rgpoe.inc_n();
      
  		lgpoe[ngpoe][0] = ipoin + incre;
  		lgpoe[ngpoe][1] = msh.bpo2ibi[ibpoi][2] + incre;
  		rgpoe[ngpoe][0] = msh.bpo2rbi[ibpoi][0];
  		rgpoe[ngpoe][1] = 0.0; // Placeholder: should be distance to ent
  	}
  }
}







} // End namespace












