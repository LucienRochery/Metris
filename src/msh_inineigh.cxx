//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



//#include "msh_inineigh.hxx"
//#include "Arrays/aux_msharrays.hxx"
//#include "aux_topo.hxx"
//#include "utils/aux_misc.hxx"
//#include "utils/aux_timer.hxx"
//#include "aux_hashtab.hxx"
//
//#include "Boundary/msh_inisurf.hxx"
//
//#include "common_includes.hxx"
//#include <absl/container/flat_hash_map.h>


#include "msh_inineigh.hxx"
#ifdef USE_ABSL
#include <absl/container/flat_hash_map.h>          // for FlatHashMapPolicy
#include <absl/hash/hash.h>                        // for Hash
#endif
#include <assert.h>                                // for assert
#include <stdio.h>                                 // for printf, fflush
#include <boost/preprocessor/iteration/local.hpp>  // for BOOST_PP_LOCAL_ITE...
#include <iostream>                                // for basic_ostream, ope...
#include <memory>                                  // for allocator
#include "aux_topo.hxx"                        // for getedgfac, isdeadent
#include "utils/aux_misc.hxx"                       // for stup2, stup3
#include "utils/mprintf.hxx"                       // for stup2, stup3
#include "Boundary/msh_inisurf.hxx"                     // for iniMeshBdryCorners
#include <sstream>                                 // for basic_stringstream
#include <tuple>                                   // for tuple
#include "aux_exceptions.hxx"                  // for METRIS_THROW_MSG
#include "Mesh/Mesh.hxx"                    
#include "types.hxx"                           // for HshTab_I2I, HshTab...
/*
Build topo links, neighbours. Also surface reconstruction. 
Since we're building face, edge hash tables, compare to those created on reading file. 
If those are empty, populate with the ones built here. 
This is free to do because we're already insert-and-deleting tet faces, triangle edges, 
so the only ones left in the end are the bdry ones. 
This must be done in the order tet -> face -> edge, as each is liable to populate the next.
*/


namespace Metris{



template <int ideg>
void iniMeshNeighbours(MeshBase &msh){

  if(msh.idim == 2){
    iniMeshNeighbours2D<ideg>(msh);
  }else if(msh.idim == 3){
    iniMeshNeighbours3D<ideg>(msh);
  }else METRIS_THROW_MSG(WArgExcept(),"Unsupported dimension "<<msh.idim);
}



// See 3D version (original) for more comments
template <int ideg>
void iniMeshNeighbours2D(MeshBase &msh){
  GETVDEPTH(msh.param);

  // Ref to give when no refs
  int iref_dum = 0;

  msh.edg2edg.fill(msh.nedge,2,-1);
  msh.fac2fac.fill(msh.nface,3,-1);

  // Note the absence of a "point hash table". We'll handle edge neighbours 
  // using msh.poi2tag[0] (the index is an implicit value)
  HshTab_I2I edgeHshTab(1.5*msh.nface + msh.nedge); // index by edge, fetch iface

  //  Algo:
  //  For each edge taken from triangle, check if in hashtab. 
  //       If not: insert, move along. Add triplet (ip1,ip2,iface)
  //       If  so: get ifac2 from inserted triplet: this is the last face that saw the edge
  //               if ifac2 > 0
  //                  Check ifac2 neighbour for this edge. 
  //                  If none (0):
  //                     add iface to ifac2 neighbours, ifac2 to iface neighbours, as positive entry (manifold)
  //                     replace ifac2 by iface in hashtab triplet. ?
  //                  if manifold (>0):
  //                     there is already a neighbour ifac3, 
  //                     ifac2 becomes the seed, it points to ifac3 (or rather, -ifac3) 
  //                     ifac3 points to -iface
  //                     iface points to -ifac2
  //                     Replace ifac2 in hashtab triplet by -iface, -ifac2 or -ifac3
  //               if ifac2 < 0: (non manifold)
  //                  This is a non manifold edge 
  //                  Replace ifac2 neighbour -ifac3 through this edge with -iface
  //                  set iface neighbour through this edge to -ifac3 
  //                  Done
  // Routine tested using nm.meshb
  //  int npp = faceNpps[msh.strdeg]; 
  //  int ifnd;
  
  // Also seize opportunity to hash mshfachsh with supplied bfaces.
  // This will be used to set up the fac -> tet and tet -> fac links in the next loop on tets.

  // In this case we can't do the find and remove trick (extract)
  // We'll have to loop over all entries in the end. 

  // + Add non manifold edges when they don't exist yet
  // + Add edges between refs 
  // That's all if the surface is closed. 
  int ncree = 0;
  for(int iface = 0; iface < msh.nface; iface++){
    if(isdeadent(iface,msh.fac2poi)) continue;
    for(int ied1 = 0; ied1 < 3; ied1++){

      int i1 = msh.fac2poi(iface,lnoed2[ied1][0]);
      int i2 = msh.fac2poi(iface,lnoed2[ied1][1]);

      auto key = stup2(i1,i2);

      //if(iprt > 0) std::cout<<"Tri "<<iface<<" edg "<<ied1<<" points "<<i1<<" "<<i2<<std::endl;

      int iinte = 0, ifac2 = -1;
      auto s = edgeHshTab.extract(key);
      if(s.empty()){
        edgeHshTab.insert({key,iface});
      }else{
        ifac2 = s.mapped();
        //if(iprt > 0) std::cout<<" - Found ifac2 = "<<ifac2<<std::endl;
        // Previous face who saw this edge is manifold, may or may not already have neighbour
        // Fetch neighbour to this edge
        int ied2; 
        try{
          ied2 = getedgfacOpp(msh,ifac2,i1,i2);
        }catch(const MetrisExcept &e){
          std::cout<<" (1) ## SOME PROBLEM WITH THIS HASH TABLE33 !"<<std::endl;
          std::cout<<"Did not find edge in ifac2 of verts "<<
          msh.fac2poi(ifac2,0)<<" "<<
          msh.fac2poi(ifac2,1)<<" "<<
          msh.fac2poi(ifac2,2)<<" "<<std::endl;
          std::cout<<"Edge = "<<i1<<" "<<i2<<std::endl;
          throw e;
        }
        msh.fac2fac(ifac2,ied2) = iface ;
        msh.fac2fac(iface,ied1) = ifac2;

        // If the edge corresponds to 2 triangles of different refs, we will need to create it
        if(msh.fac2ref[iface] != msh.fac2ref[ifac2]) iinte = 1;
      }

      // Now we look for the edge in the boundary edge hash table
      auto t = msh.edgHshTab.find(key);

      if(t!=msh.edgHshTab.end()){
        int iedge = t->second;
        METRIS_ASSERT_MSG(iedge >= 0 && iedge < msh.nedge,"Glo edg hsh tab valid");
        msh.edg2fac[iedge] = iface;
      }else if(iinte == 1){
        // If the edge was sandwiched between facs of diff refs, we must create it
        // Note that this only applies if the previous failed, i.e. face does not exist
        msh.newedgtopo<ideg>(iface, ied1, iref_dum);
        msh.edg2fac[msh.nedge-1] = iface;
        ncree++;
      }
    }
    //if(iprt > 0) std::cout<<"iter end  "<<iface<<" print all nei"<<std::endl;
  }
  if(ncree > 0)CPRINTF1(" - created %d edges \n",ncree);
  //if(iprt > 0)printf("-- Faces finished\n");


  iniMeshBdryEdges<ideg>(msh);


  // Do the same for edges
  // For this, we simply use msh.poi2tag[0] to store the "value";
  // edge facets are simply vertices !
  // Add corners at triple+ points and diff refs
  for(int iedge = 0; iedge < msh.nedge; iedge++){
    if(isdeadent(iedge,msh.edg2poi)) continue;
    msh.poi2tag[0][msh.edg2poi(iedge,0)] = -1;
    msh.poi2tag[0][msh.edg2poi(iedge,1)] = -1;
  }

  for(int iedge = 0; iedge < msh.nedge; iedge++){
    if(isdeadent(iedge,msh.edg2poi)) continue;
    //if(iprt > 0) printf("Start iedge = %d vertices = %d %d \n",iedge,msh.edg2poi(iedge,0),msh.edg2poi(iedge,1));
    for(int ive1 = 0; ive1 < 2; ive1++){
      int ip = msh.edg2poi(iedge,ive1);
      //if(iprt>0)printf("  ive1 = %d ip = %d \n",ive1,ip);
      int iedge2 = msh.poi2tag(0,ip); 
      //if(iprt>0)printf("  iedge2 = %d \n",iedge2);

      // Case no neighbour yet
      if(iedge2 == -1){
        msh.poi2tag(0,ip) = iedge; 
        continue;
      }
      //if(iprt > 0) std::cout<<" - Found iedge2 = "<<iedge2<<std::endl;

      int inewc = 0;
      // If iedge2 >= 0, that point has already seen another edge (iedge2). 
      if(iedge2 >= 0){
        // Previous point who saw this edge is manifold, may or may not already have neighbour
        // Fetch neighbour to this edge
        int                           ive2 = 0;
        if(msh.edg2poi(iedge2,1) == ip) ive2 = 1;
        //if(iprt > 0)printf("  msh.edg2poi iedge2: %d %d \n",msh.edg2poi(iedge2,0),msh.edg2poi(iedge2,1));
        assert(("Hash table (1)" && msh.edg2poi(iedge2,ive2) == ip));
        int iedge3 = msh.edg2edg[iedge2][1-ive2];

        if(iedge3 == -1){ 
          // There is a neighbour but the point is still manifold
          //if(iprt>0)printf("  Initialize neighbour (%d,%d) -> %d and (%d,%d) -> %d \n",
          //  iedge2,1-ive2,iedge,iedge,1-ive1,iedge2);
          msh.edg2edg[iedge2][1-ive2] = iedge  ;
          msh.edg2edg[iedge ][1-ive1] = iedge2 ;
          // Tag point to be added as corner if refs differ
          if(msh.edg2ref[iedge] != msh.edg2ref[iedge2]) inewc = 1;
        }else{ 
          // The point is about to become nm, as iedge2 already had a (mfold) neigh iedge3.
          // Tag point to be added as corner
          inewc = 1;
          assert(("Skipped msh.poi2tag[0] update" && iedge3 > -1));

          //if(iprt>0)printf("  Found neighbour %d \n",iedge3);

          assert(("Array corruption" && iedge3 != iedge)); 

          // Non manifold edges
          // Find by which vertex iedge3 is being non manifold here
          int                           ive3 = 0;  
          if(msh.edg2poi(iedge3,1) == ip) ive3 = 1;

          assert(("Hash table (2)" && msh.edg2poi(iedge3,ive3) == ip));

          //if(iprt>0)printf("  Update using 1,2,3 = (%d,%d) (%d,%d) (%d,%d)\n",
          //                              iedge,1-ive1,iedge2,1-ive2,iedge3,1-ive3);
          msh.edg2edg[iedge2][1-ive2] = -(iedge3+2) ;
          msh.edg2edg[iedge3][1-ive3] = -(iedge +2) ;
          msh.edg2edg[iedge ][1-ive1] = -(iedge2+2) ;

          // We must update the value in the hash table to reflect this. 
          // This should certainly not lead to an insertion.
          //if(iprt>0)printf("  Update ip = %d tag to %d \n",ip,-iedge3-1);
          msh.poi2tag(0,ip) =  -(iedge3+2);
        }
      }else{
        // We already know the face is non manifold
        iedge2 = - (iedge2 + 2) ; 
        assert(("Array corruption" && iedge2 != iedge));

        //if(iprt>0)printf("Already non manifold edge, get iedge2 = %d (For)\n",iedge2+1);

        int ive2 = -1;
        if(msh.edg2poi(iedge2,0) == ip) ive2 = 0;
        if(msh.edg2poi(iedge2,1) == ip) ive2 = 1;
        assert(("Hash table (3)" && ive2 >= 0));

        int iedge3  = msh.edg2edg[iedge2][1-ive2];
        assert(("Array corruption" && iedge3 < -1));

        iedge3 = - (iedge3 + 2);

        //if(iprt>0)printf("  Got iedge3 = %d \n",iedge3);

        msh.edg2edg[iedge2][1-ive2] = -(iedge  + 2);
        msh.edg2edg[iedge ][1-ive1] = -(iedge3 + 2);
      }
      if(inewc > 0){
        // Create new corner if does not exist
        // Do not update edg2bpo yet because we only have partial information
        // Namely, some corners (boundary) will be created after this loop. 
        // Points always point to the corner if it exists
        int ibpoi = msh.poi2bpo[ip];
        if(ibpoi < 0 || (
          ibpoi > 0 && msh.bpo2ibi(ibpoi,1) > 0)){
          msh.newbpotopo(ip,0);
        }
      }
    }


    //if(iprt > 1) std::cout<<"iter end  "<<iedge<<" print all nei"<<std::endl;
    //if(iprt > 1) msh.edg2edg.print(msh.nedge); 
  }

  // Reset tag to 0, it's few points
  for(int iedge = 0; iedge < msh.nedge; iedge++){
    if(isdeadent(iedge,msh.edg2poi)) continue;
    msh.poi2tag[0][msh.edg2poi(iedge,0)] = 0;
    msh.poi2tag[0][msh.edg2poi(iedge,1)] = 0;
  }
  iniMeshBdryCorners(msh);
}



template <int ideg>
void iniMeshNeighbours3D(MeshBase &msh){
  GETVDEPTH(msh.param);

  int iref_dum = 0;

  msh.edg2edg.fill(msh.nedge,2,-1);
  msh.fac2fac.fill(msh.nface,3,-1);
  msh.tet2tet.fill(msh.nelem,4,-1);

  // Note the absence of a "point hash table". We'll handle edge neighbours 
  // using msh.poi2tag[0] (the index is an implicit value)
  HshTab_I2I edgeHshTab(1.5*msh.nface + msh.nedge); // index by edge, fetch iface
  HshTab_I3I intfHshTab(2*msh.nelem + msh.nface);

  // Input: empty intfHshTab
  // Output: intfHshTab contains faces seen only once (bdry)
  //         tet2tet filled
  //         fac2tet filled for those faces already in facHshTab 
  //         faces ween twice with different tet refs have been added to glo hash tab
  //           and nface updated if does not exist already
  //for(int ielem = 0; ielem < msh.nelem; ielem++){
  //  if(isdeadent(ielem,msh.tet2poi))continue;
  //  int ip1 = 124, ip2 = 316, ip3 = 71;
  //  int ifa = getfactet(msh,ielem,ip1,ip2,ip3);
  //  if(ifa >= 0){
  //    printf("## DEBUG FOUND IELEM %d ifa = %d \n",ielem,ifa);
  //  }
  //}
  //wait();
  CPRINTF2("-- Start tetra loop nelem = %d \n",msh.nelem);
  for(int ielem = 0; ielem < msh.nelem; ielem++){
    INCVDEPTH(msh.param);
    if(isdeadent(ielem,msh.tet2poi)) continue;
    for(int ifa1 = 0; ifa1 < 4; ifa1++){
      INCVDEPTH(msh.param);

      int i1 = msh.tet2poi(ielem,lnofa3[ifa1][0]);
      int i2 = msh.tet2poi(ielem,lnofa3[ifa1][1]);
      int i3 = msh.tet2poi(ielem,lnofa3[ifa1][2]);

      auto key = stup3(i1,i2,i3);

      CPRINTF3(" - tetra %d face %d points %d %d %d \n",ielem,ifa1,i1,i2,i3);

      int iintf = 0, ielem2 = -1;
      auto s = intfHshTab.extract(key);
      if(s.empty()){
        intfHshTab.insert({key,ielem});
      }else{
        ielem2 = s.mapped();
        // Previous face who saw this edge is manifold, may or may not already have neighbour
        // Fetch neighbour to this edge
        int ifa2 = getfactetOpp(msh,ielem2,i1,i2,i3);
        CPRINTF3(" - Found ielem2 = %d ifa2 = %d \n",ielem2,ifa2);
        if(ifa2 < 0){
          std::cout<<" (2) ## SOME PROBLEM WITH THIS HASH TABLE33 !"<<std::endl;
          std::cout<<"Did not find face in ielem2 of verts "<<
          msh.tet2poi(ielem2,0)<<" "<<
          msh.tet2poi(ielem2,1)<<" "<<
          msh.tet2poi(ielem2,2)<<" "<<
          msh.tet2poi(ielem2,3)<<" "<<std::endl;
          std::cout<<"Face = "<<i1<<" "<<i2<<" "<<i3<<std::endl;
          METRIS_THROW(TopoExcept())
        }
        msh.tet2tet(ielem2,ifa2) = ielem ;
        msh.tet2tet(ielem ,ifa1) = ielem2;

        // If the face corresponds to 2 tetrahedra of different refs, we will need to create it
        if(msh.tet2ref[ielem] != msh.tet2ref[ielem2]) iintf = 1;
      }

      // Now we look for the face in the boundary face hash table
      auto t = msh.facHshTab.find(key);

      if(t!=msh.facHshTab.end()){
        // Found this face in the global face hash table: update msh.fac2tet
        int iface = t->second;
        if(msh.fac2tet(iface,0) > -1){
          msh.fac2tet(iface,1) = ielem;
        }else{
          msh.fac2tet(iface,0) = ielem;
        }
      }else if(iintf == 1){
        // If the face was sandwiched between tets of diff refs, we must create it
        // Note that this only applies if the previous failed, i.e. face does not exist
        msh.newfactopo<ideg>(ielem, ifa1, iref_dum, ielem2);
      }
    }
  }
  // We are now going to use the previous information to finish creating all 
  // missing triangles. As we already created interior faces, the only missing ones 
  // are boundary triangles contained in intfHshTab.
  iniMeshBdryTriangles<ideg>(msh,  intfHshTab);

  //  Algo:
  //  For each edge taken from triangle, check if in hashtab. 
  //       If not: insert, move along. Add triplet (ip1,ip2,iface)
  //       If  so: get ifac2 from inserted triplet: this is the last face that saw the edge
  //               if ifac2 > 0
  //                  Check ifac2 neighbour for this edge. 
  //                  If none (0):
  //                     add iface to ifac2 neighbours, ifac2 to iface neighbours, as positive entry (manifold)
  //                     replace ifac2 by iface in hashtab triplet. ?
  //                  if manifold (>0):
  //                     there is already a neighbour ifac3, 
  //                     ifac2 becomes the seed, it points to ifac3 (or rather, -ifac3) 
  //                     ifac3 points to -iface
  //                     iface points to -ifac2
  //                     Replace ifac2 in hashtab triplet by -iface, -ifac2 or -ifac3
  //               if ifac2 < 0: (non manifold)
  //                  This is a non manifold edge 
  //                  Replace ifac2 neighbour -ifac3 through this edge with -iface
  //                  set iface neighbour through this edge to -ifac3 
  //                  Done
      // Routine tested using nm.meshb
  //  int npp = faceNpps[msh.strdeg]; 
  //  int ifnd;
  
  // Also seize opportunity to hash mshfachsh with supplied bfaces.
  // This will be used to set up the fac -> tet and tet -> fac links in the next loop on tets.

  // In this case we can't do the find and remove trick (extract)
  // We'll have to loop over all entries in the end. 

  // + Add non manifold edges when they don't exist yet
  // + Add edges between refs 
  // That's all if the surface is closed. 
  int ncree = 0;
  for(int iface = 0; iface < msh.nface; iface++){
    if(isdeadent(iface,msh.fac2poi)) continue;
    INCVDEPTH(msh.param);

    for(int ied1 = 0; ied1 < 3; ied1++){
      INCVDEPTH(msh.param);
      int i1 = msh.fac2poi(iface,lnoed2[ied1][0]);
      int i2 = msh.fac2poi(iface,lnoed2[ied1][1]);
      auto key = stup2(i1,i2);
      CPRINTF3(" - face %d edge %d points %d %d \n",iface,ied1,i1,i2);
      int inewe = 0;

      auto t = edgeHshTab.find(key);
      if(t == edgeHshTab.end()){
        CPRINTF3(" new edge, add to hash tab\n");
        edgeHshTab.insert({key,iface}); 
      }else{
        int ifac2 = t->second;
        CPRINTF3(" - found ifac2 = %d\n",ifac2);
        if(ifac2 >= 0){

          // Previous face who saw this edge is manifold, may or may not already have neighbour
          // Fetch neighbour to this edge
          int ied2 = getedgfac(msh,ifac2,i1,i2);
//          assert("Triangle has expected edge " && ied2 > -1);

          int ifac3  = msh.fac2fac(ifac2,ied2);
          CPRINTF3(" - found ifac3 = %d\n",ifac3);
          if(ifac3 == -1){
            // In this case, the previous edge did not yet have a neighbour. 
            msh.fac2fac(ifac2,ied2) = iface ;
            msh.fac2fac[iface ][ied1] = ifac2;
            // Tag edge to be added on grounds of different refs 
            if(msh.fac2ref[ifac2] != msh.fac2ref[iface]) inewe = 1;            
          }else{
            // If this case is verified, the previous case ifac3 == -1 has already been 
            // seen with this edge before. Still, we need to set inewe = 1 because 
            // the edge is now non-manifold. 
            assert("Manifold neighbour doesn't have non manifold neighbour" && ifac3 >= -1);
            assert("Neighbour doesn't already know us as neighbour" && ifac3!=iface);

            inewe = 1;

            // Faces are, in fact, non manifold !
            // Find by which edge ifac3 is being non manifold here
            int ied3 = getedgfac(msh,ifac3,i1,i2);         
//            assert("Neighbour's neighbour has edge" && ied3 >=0);

            msh.fac2fac(ifac2,ied2) = -(ifac3 + 2);
            msh.fac2fac(ifac3,ied3) = -(iface  + 2);
            msh.fac2fac[iface ][ied1] = -(ifac2 + 2);

            // We must update the value in the hash table to reflect this. 
            // This should certainly not lead to an insertion.
            edgeHshTab[key] =  -(ifac3 + 2);
          }
        }else{
          // We already knew this edge was non manifold: it has already been added
          // to the glo hash table if needed. 
          if(ifac2 == -1) METRIS_THROW_MSG(TopoExcept(),
            "SOME PROB (4) WITH THIS HASH TABLE")
          ifac2 = -ifac2 -2;
          // We already know the face is non manifold
          int ied2 = getedgfac(msh,ifac2,i1,i2);
          int ifac3  = msh.fac2fac(ifac2,ied2);
          if(ifac3 >= 0) METRIS_THROW_MSG(TopoExcept(),
            "SOME PROB (7) WITH THIS HASH TABLE");

          msh.fac2fac(ifac2,ied2) = -(iface  + 2);
          msh.fac2fac[iface ][ied1] = ifac3;

        }
      }

      // Lastly, update edg2fac and create edge if inewe = 1
      auto tglo = msh.edgHshTab.find(key);
      if(tglo != msh.edgHshTab.end()){
        int iedge = tglo -> second;
        assert("Glo edg hsh tab valid" && iedge >= 0 && iedge < msh.nedge);
        msh.edg2fac[iedge] = iface;
      }else if(inewe == 1){
        ncree++;
        msh.newedgtopo<ideg>(iface,ied1,iref_dum);
        msh.edg2fac[msh.nedge-1] = iface;
//        #ifndef NDEBUG
//          printf("Debug created edge ideg = %d iface = %d ied1 = %d \n",
//                   ideg,iface,ied1);
//        #endif
      }
    }
  }
  if(ncree > 0) CPRINTF2("-- Created %d edges \n",ncree);
  CPRINTF3("-- Faces finished\n");

  iniMeshBdryEdges<ideg>(msh);


  // Do the same for edges
  // For this, we simply use msh.poi2tag[0] to store the "value";
  // edge facets are simply vertices !
  // Add corners at triple+ points and diff refs
  for(int iedge = 0; iedge < msh.nedge; iedge++){
    if(isdeadent(iedge,msh.edg2poi)) continue;
    msh.poi2tag[0][msh.edg2poi(iedge,0)] = -1;
    msh.poi2tag[0][msh.edg2poi(iedge,1)] = -1;
  }
 

  for(int iedge = 0; iedge < msh.nedge; iedge++){
    if(isdeadent(iedge,msh.edg2poi)) continue;
    INCVDEPTH(msh.param);

    for(int ive1 = 0; ive1 < 2; ive1++){
      INCVDEPTH(msh.param);
      int ip = msh.edg2poi(iedge,ive1);
      CPRINTF3(" - iedge %d ive1 = %d ip = %d \n",iedge,ive1,ip);

      int iedge2 = msh.poi2tag(0,ip); 
      CPRINTF3("  iedge2 = %d \n",iedge2);

      // Case no neighbour yet
      if(iedge2 == -1){
        msh.poi2tag(0,ip) = iedge; 
        continue;
      }
      CPRINTF3(" - found iedge2 = %d \n",iedge2);

      int inewc = 0;
      // If iedge2 >= 0, that point has already seen another edge (iedge2). 
      if(iedge2 >= 0){
        // Previous point who saw this edge is manifold, may or may not already have neighbour
        // Fetch neighbour to this edge
        int                           ive2 = 0;
        if(msh.edg2poi(iedge2,1) == ip) ive2 = 1;
        CPRINTF3("  msh.edg2poi iedge2: %d %d \n",msh.edg2poi(iedge2,0),msh.edg2poi(iedge2,1));
        assert(("Hash table (1)" && msh.edg2poi(iedge2,ive2) == ip));
        int iedge3 = msh.edg2edg[iedge2][1-ive2];

        if(iedge3 == -1){ 
          // There is a neighbour but the point is still manifold
          CPRINTF3("  Initialize neighbour (%d,%d) -> %d and (%d,%d) -> %d \n",
                   iedge2,1-ive2,iedge,iedge,1-ive1,iedge2);
          msh.edg2edg[iedge2][1-ive2] = iedge  ;
          msh.edg2edg[iedge ][1-ive1] = iedge2 ;
          // Tag point to be added as corner if refs differ
          if(msh.edg2ref[iedge] != msh.edg2ref[iedge2]) inewc = 1;
        }else{ 
          // The point is about to become nm, as iedge2 already had a (mfold) neigh iedge3.
          // Tag point to be added as corner
          inewc = 1;
          assert(("Skipped msh.poi2tag[0] update" && iedge3 > -1));

          CPRINTF3("  Found neighbour %d \n",iedge3);

          assert(("Array corruption" && iedge3 != iedge)); 

          // Non manifold edges
          // Find by which vertex iedge3 is being non manifold here
          int                           ive3 = 0;  
          if(msh.edg2poi(iedge3,1) == ip) ive3 = 1;

          assert(("Hash table (2)" && msh.edg2poi(iedge3,ive3) == ip));

          CPRINTF3("  Update using 1,2,3 = (%d,%d) (%d,%d) (%d,%d)\n",
                                        iedge,1-ive1,iedge2,1-ive2,iedge3,1-ive3);
                                        msh.edg2edg[iedge2][1-ive2] = -(iedge3+2) ;
          msh.edg2edg[iedge3][1-ive3] = -(iedge +2) ;
          msh.edg2edg[iedge ][1-ive1] = -(iedge2+2) ;

          // We must update the value in the hash table to reflect this. 
          // This should certainly not lead to an insertion.
          CPRINTF3("  Update ip = %d tag to %d \n",ip,-iedge3-1);
          msh.poi2tag(0,ip) =  -(iedge3+2);
        }
      }else{
        // We already know the face is non manifold
        iedge2 = - (iedge2 + 2) ; 
        assert(("Array corruption" && iedge2 != iedge));

        CPRINTF3("Already non manifold edge, get iedge2 = %d (For)\n",iedge2+1);

        int ive2 = -1;
        if(msh.edg2poi(iedge2,0) == ip) ive2 = 0;
        if(msh.edg2poi(iedge2,1) == ip) ive2 = 1;
        assert(("Hash table (3)" && ive2 >= 0));

        int iedge3  = msh.edg2edg[iedge2][1-ive2];
        assert(("Array corruption" && iedge3 < -1));

        iedge3 = - (iedge3 + 2);

        CPRINTF3("  Got iedge3 = %d \n",iedge3);

        msh.edg2edg[iedge2][1-ive2] = -(iedge  + 2);
        msh.edg2edg[iedge ][1-ive1] = -(iedge3 + 2);
      }
      if(inewc > 0){
        // Create new corner if does not exist
        // Do not update edg2bpo yet because we only have partial information
        // Namely, some corners (boundary) will be created after this loop. 
        // Points always point to the corner if it exists
        int ibpoi = msh.poi2bpo[ip];
        if(ibpoi < 0 || (ibpoi > 0 && msh.bpo2ibi(ibpoi,1) > 0)){
          msh.newbpotopo(ip,0);
        }
      }
    }

  }

  // Reset tag to 0, it's few points
  for(int iedge = 0; iedge < msh.nedge; iedge++){
    if(isdeadent(iedge,msh.edg2poi)) continue;
    msh.poi2tag[0][msh.edg2poi(iedge,0)] = 0;
    msh.poi2tag[0][msh.edg2poi(iedge,1)] = 0;
  }
  iniMeshBdryCorners(msh);
}



// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template void iniMeshNeighbours< n >(MeshBase &msh);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


#if 0
// See 3D version (original) for more comments
template <int ideg>
void iniMeshNeighbours0(MeshBase &msh){
  GETVDEPTH(msh.param);

  msh.edg2edg.fill(msh.nedge,2,-1);
  msh.fac2fac.fill(msh.nface,3,-1);
  msh.tet2tet.fill(msh.nelem,4,-1);

  // Note the absence of a "point hash table". We'll handle edge neighbours 
  // using msh.poi2tag[0] (the index is an implicit value)
  HshTab_I2I edgeHshTab(1.5*msh.nface + msh.nedge); // index by edge, fetch iface
  HshTab_I3I intfHshTab(msh.nelem > 0 ? 2*msh.nelem + msh.nface : 0); // not used if 2D

  // Go from highest to lowest tdim, building missing lower dim entities as needed
  // Use CT LOOP to use tdim as time (e.g. for stupn)
  CT_FOR0_INC(0,2,tdim0){
    constexpr int tdim = 3 - tdim0;
    constexpr int nneil = tdim + 1;
    const int nentt = msh.nentt(tdim);
    const intAr2 &ent2poi = msh.ent2poi(tdim);
          intAr2 &ent2ent = msh.ent2ent(tdim);

    constexpr auto lnosub = getlnosub<tdim>();
    auto subHshTab = tdim == 3 ? intfHshTab : edgeHshTab;

    auto subHshTabMesh = tdim == 3 ? msh.facHshTab : msh.edgHshTab;

    for(int ientt = 0; ientt < nentt; ientt++){
      for(int ineil = 0; ineil < nneil; ineil++){

        if constexpr(tdim >= 2){

          int idxsub[tdim]; 
          for(int ii = 0; ii < tdim; ii++) 
            idxsub[ii] = ent2poi(ientt,lnosub[ineil][ii]);
          auto key = stupn<tdim>(idxsub);
          auto sfnd = subHshTab.extract(key);
          int ient2 = -1;
          if(sfnd.empty()){
            subHshTab.insert({key,ientt});
          }else{
            ient2 = sfnd.mapped();
            int inei2 = -1;
            bool iints = false;
            if constexpr(tdim == 2){
              inei2 = getedgfacOpp(msh,ient2,idxsub[0],idxsub[1]);
            }else{
              inei2 = getfactetOpp(msh,ient2,idxsub[0],idxsub[1],idxsub[2]);
            }
            CPRINTF3(" - tdim %d found ient2 = %d inei2 = %d \n",ielem2,inei2);
            METRIS_ASSERT(inei2 >= 0);

            ent2ent(ient2, inei2) = ientt;
            ent2ent(ientt, ineil) = ient2;

            // If the subent corresponds to 2 entts of different refs, we will 
            // need to create it if it doesn't exist.
            if(msh.ent2ref(tdim)[ientt] != msh.ent2ref(tdim)[ient2]) iints = true;
          }// if snfd.empty()

          // Now we look for the facet (subent) in the boundary face hash table
          auto tfnd = subHshTabMesh.find(key);

          if(tfnd != msh.facHshTab.end()){
            // Found this subent in the global subent hash table: update msh.fac2tet
            // or msh.edg2fac
            int isube = tfnd->second;
            if constexpr(tdim == 3){
              if(msh.fac2tet(isube,0) > -1){
                msh.fac2tet(isube,1) = ientt;
              }else{
                msh.fac2tet(isube,0) = ientt;
              }
            }else{
              msh.edg2fac[isube] = ientt;
            }
          }else if(iintf == 1){
            // If the face was sandwiched between tets of diff refs, we must create it
            // Note that this only applies if the previous failed, i.e. face does not exist
            if constexpr(tdim == 3){
              msh.newfactopo<ideg>(ientt, ineil, 0, ient2);
            }else{
              msh.newedgtopo<ideg>(ientt, ineil, 0);
              msh.edg2fac[msh.nedge-1] = ientt;
            }
          }// if tfnd != msh.facHshTab.end()

        }// if tdim >= 2

      }// for ineil
    }// for ientt

    if constexpr(tdim == 3){
      iniMeshBdryTriangles<ideg>(msh,  intfHshTab);
    }else if(tdim == 2){

    }

  }CT_FOR1(tdim0);

}

// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template void iniMeshNeighbours0< n >(MeshBase &msh);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()
#endif



} // End namespace


