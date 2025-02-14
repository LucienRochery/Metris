//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../aux_utils.hxx"
#include "../ho_constants.hxx"
#include "../io_libmeshb.hxx"
#include "../CT_loop.hxx"
#include "../mprintf.hxx"
#include "../linalg/det.hxx"
#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../low_geo.hxx"

namespace Metris{





// redge and rface hold the uvs as provided by reconnect_ routines when new faces are created
// this is because we lose the link between which cavity element created which element. 
// they are dimensioned for nface - nfac0, nedge - nedg0.
// redge has info for the point that is not ipins
// rface is dblAr3 as 2 nodes, each 2 parameters
template<class MFT, int ideg>
int update_cavity(Mesh<MFT> &msh, const MshCavity &cav, const CavWrkArrs &work,
                  int npoi0, int nedg0, int nfac0, int nele0, 
                  int ithread){

  GETVDEPTH(msh);

  const int ncedg = cav.lcedg.get_n();
  const int ncfac = cav.lcfac.get_n();
  const int nctet = cav.lctet.get_n();

  if(DOPRINTS2()){
    MshCavity cav2(msh.nelem-nele0,msh.nedge-nedg0,msh.nface-nfac0);
    for(int ii = nele0; ii < msh.nelem; ii++) cav2.lctet.stack(ii);
    for(int ii = nfac0; ii < msh.nface; ii++) cav2.lcfac.stack(ii);
    for(int ii = nedg0; ii < msh.nedge; ii++) cav2.lcedg.stack(ii);
    cav2.ipins = cav.ipins;
    writeMeshCavity("cavity1",msh,cav2);
  }

  // Tag cavity entities 
  msh.tag[ithread]++;
  CT_FOR0_INC(1,3,tdimn){
    for(int ientt : cav.lcent<tdimn>()){
      msh.template ent2tag<tdimn>()[ithread][ientt] = msh.tag[ithread];
    }
  }CT_FOR1(tdimn);


  // -- 1 Manage bpois

  // -- 1.0  Get (u,v)s in case of a pdim < 2 point. If point is tdim 2, (u,v) is 
  // computed by the driver. We need to do this before deleting bpois.
  int pdim = -1;
  int ibins = msh.poi2bpo[cav.ipins];
  if(ibins >= 0) pdim = msh.bpo2ibi(ibins,1);
  if(pdim < 2 && msh.isboundary_faces()){ 
    // We could be clever in the case where ipins already belonged to the mesh
    // and recycle known (u,v)s. But we're not doing that yet, let's keep it 
    // generic. 
    // What info do we have? Surface connex components are stored in lfcco, 
    // which points to seeds in initial cavity. 
    // We'll use those to do an EG_invEvaluateGuess.
    // Similarly, we could use EG_getEdgeUV, but then we'd have to deal with 
    // isens, and with the corner case. So let's remain sane and use EG_invEvaluateGuess. 
    // This could still be done dumbly by using the seed and computing distance
    // to obtained (u,v)s, keeping the closest as the correct one. 
    double result[18];
    for(ibins = msh.poi2bpo[cav.ipins]; ibins >= 0; ibins = msh.bpo2ibi(ibins,3)){
      int bdim = msh.bpo2ibi(ibins,1);
      if(bdim != 2) continue;

      // This can be either an old or a new face as we've got bpois for both 
      // at this stage. 
      int iface = msh.bpo2ibi(ibins,2);
      METRIS_ASSERT(iface >= 0 && iface < msh.nface);

      // If it's an old one, skip.
      if(iface < nfac0) continue;

      // We can get its connex component:
      int icoco = msh.fac2tag(ithread,iface); 
      METRIS_ASSERT(icoco >= 0 && icoco < work.lfcco.get_n());

      // Get a face from before, from this connex component
      int ifaco = work.lfcco[icoco];

      for(int ii = 0; ii < nrbi; ii++) msh.bpo2rbi(ibins,ii) = 0;

      // Use only the vertices, that'll be enough. 
      for(int ii = 0; ii < 3; ii++){
        int ipoin = msh.fac2poi(ifaco, ii);
        METRIS_ASSERT(ipoin >= 0);
        int ibpoi = msh.poi2ebp(ipoin, 2, ifaco, -1);
        METRIS_ASSERT(ibpoi >= 0);
        for(int ii = 0; ii < nrbi; ii++) 
          msh.bpo2rbi(ibins,ii) += msh.bpo2rbi(ibpoi,ii) / 3.0;
      }

      int iref = msh.fac2ref[iface];
      METRIS_ASSERT(iref >= 0);
      ego obj = msh.CAD.cad2fac[iref];

      int ierro = EG_invEvaluateGuess(obj, msh.coord[cav.ipins], msh.bpo2rbi[ibins], result);
      if(ierro != 0){
        CPRINTF1("## EG_invEvaluateGuess ERROR %d \n",ierro);
        goto cleanup;
      }


    }
  }


  // -- 1.1 remove old ibpois

  //bool dowait = false;
  CT_FOR0_INC(1,2,tdimn){
    if (msh.isboundary_tdimn(tdimn)){
    int nnode = msh.nnode(tdimn);
    for(int ientt : cav.lcent<tdimn>()){
      INCVDEPTH(msh);
      for(int ii = 0; ii < nnode; ii++){
        int ip = msh.template ent2poi<tdimn>()[ientt][ii];
        if(msh.poi2tag(ithread,ip) >= msh.tag[ithread]) continue;
        msh.poi2tag(ithread,ip) = msh.tag[ithread];
        if(DOPRINTS1()){
          CPRINTF1(" - ip = %d clean bpo pre:\n",ip);
          print_bpolist(msh,msh.poi2bpo[ip]);
        }
        msh.rembpotag(ip,ithread);
        if(DOPRINTS1()){
          CPRINTF1(" - bpo post:\n");
          print_bpolist(msh,msh.poi2bpo[ip]);
        }
      }
    }
    }
  }CT_FOR1(tdimn);

  //if(dowait){
  //  printf("## WAITING BECAUSE OF IP 2 REMBPO\n");
  //  printf("New bpolist:\n");
  //  for(int ibpoi = msh.poi2bpo[2]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
  //    printf("%d t = %f : ",ibpoi,msh.bpo2rbi(ibpoi,0));
  //    intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
  //    if(abs(msh.bpo2rbi(ibpoi,0) - 1) < 1.0e-6) dowait = false;
  //  }
  //  if(dowait) wait();
  //}


  //if(msh.isboundary_edges()){
  //  int nnode = edgnpps[msh.curdeg];
  //  for(int iedgl = 0; iedgl < ncedg; iedgl++){
  //    int iedge = cav.lcedg[iedgl];
  //    for(int ii = 0; ii < nnode; ii++){
  //      int ip = msh.edg2poi(iedge,ii);
  //      if(msh.poi2tag(ithread,ip) >= msh.tag[ithread]) continue;
  //      printf("debug ip = %d \n",ip);
  //      msh.poi2tag(ithread,ip) = msh.tag[ithread];
  //      if(iverb >= METRIS_CAV_PRTLEV + 1){
  //        printf("   - ip = %d clean bpo pre:\n",ip);
  //        print_bpolist(msh,msh.poi2bpo[ip]);
  //      }
  //      msh.rembpotag(ip,ithread);
  //      if(iverb >= METRIS_CAV_PRTLEV + 1){
  //        printf("   - bpo post:\n");
  //        print_bpolist(msh,msh.poi2bpo[ip]);
  //      }
  //    }
  //  }
  //}
  //if(msh.isboundary_faces()){
  //  int nnodf = facnpps[msh.curdeg];
  //  for(int ifacl = 0; ifacl < ncfac; ifacl++){
  //    int iface = cav.lcfac[ifacl];
  //    msh.fac2tag(ithread,iface) = msh.tag[ithread];
  //    for(int ii = 0; ii < nnodf; ii++){
  //      int ip = msh.fac2poi(iface,ii);
  //      if(msh.poi2tag(ithread,ip) >= msh.tag[ithread]) continue;
  //      printf("debug ip = %d \n",ip);
  //      msh.poi2tag(ithread,ip) = msh.tag[ithread];
  //      if(iverb >= METRIS_CAV_PRTLEV + 1){
  //        printf("   - ip = %d clean bpo pre:\n",ip);
  //        print_bpolist(msh,msh.poi2bpo[ip]);
  //      }
  //      msh.rembpotag(ip,ithread);
  //      if(iverb >= METRIS_CAV_PRTLEV + 1){
  //        printf("   - bpo post:\n");
  //        print_bpolist(msh,msh.poi2bpo[ip]);
  //      }
  //    }
  //  }
  //}

  // ipins may be on a "virtual" dead edge/face 
  if(msh.poi2bpo[cav.ipins] >= 0) msh.rembpotag(cav.ipins,ithread);


  // During reconnect_faccav and reconnect_lincav, we called newbpotopo on each
  // new entity. Now, we need to clean up. This was necessary as not stacking any
  // new bpos might have made a point unlinked after old elt deletion. 
  for(int tdimn : {1,2}){
    // Edges are probably always boundary but still.
    if(!msh.isboundary_tdimn(tdimn)) continue;
    int nnode = msh.nnode(tdimn);
    int nentt = msh.nentt(tdimn);
    int nent0 = tdimn == 1 ? nedg0 : nfac0;
    intAr2 &ent2poi = msh.ent2poi(tdimn);
    for(int ientt = nent0; ientt < nentt; ientt++){
      for(int ii = 0; ii < nnode; ii++){
        int ip = ent2poi(ientt,ii);
        if(msh.poi2tag(ithread,ip) >= msh.tag[ithread]+1) continue;
        msh.poi2tag(ithread,ip) = msh.tag[ithread]+1;

        // First link and minimum topo dimn
        int ibpo0 = msh.poi2bpo[ip];
        METRIS_ASSERT(ibpo0 >= 0); 
        int mindi = msh.bpo2ibi(ibpo0,1);

        // Init straight at next and current dimn
        int ibpoi = msh.bpo2ibi(ibpo0,3);
        if(ibpoi >= 0){
          int curdi = msh.bpo2ibi(ibpoi,1);
          while(curdi <= mindi){
            int ibpon = msh.bpo2ibi(ibpoi,3);
            msh.bpo2ibi(ibpo0,3) = ibpon;
            ibpoi = ibpon;
            if(ibpoi < 0) break; 
            curdi = msh.bpo2ibi(ibpoi,1);
          }
        }//endif(ibpoi >= 0)

      }
    }
  }



  // Deleting entities in the initial cavity. 
  CT_FOR0_INC(1,2,tdime){
    for(int ientt : cav.lcent<tdime>()){
      auto key = stupn<tdime+1>(msh.template ent2poi<tdime>()[ientt]);
      // This removes if exists
      msh.template hshTab<tdime>().extract(key);
    }
  }CT_FOR1(tdime);
  //for(int iedgl = 0; iedgl < ncedg; iedgl++){
  //  int iedge = cav.lcedg[iedgl];
  //  msh.edg2tag(ithread,iedge) = msh.tag[ithread];

  //  int ip1 = msh.edg2poi(iedge,0);
  //  int ip2 = msh.edg2poi(iedge,1);
  //  auto key = stup2(ip1,ip2);
  //  // This removes if exists
  //  msh.edgHshTab.extract(key); 
  //}
  //for(int ifacl = 0; ifacl < ncfac; ifacl++){
  //  int iface = cav.lcfac[ifacl];
  //  msh.fac2tag(ithread,iface) = msh.tag[ithread];

  //  int ip1 = msh.fac2poi(iface,0);
  //  int ip2 = msh.fac2poi(iface,1);
  //  int ip3 = msh.fac2poi(iface,2);
  //  auto key = stup3(ip1,ip2,ip3);
  //  // This removes if exists
  //  msh.facHshTab.extract(key); 
  //}
  //for(int ielel = 0; ielel < nctet; ielel++){
  //  int ielem = cav.lctet[ielel];
  //  msh.tet2tag(ithread,ielem) = msh.tag[ithread];
  //}

  // remove poi2ents 
  for(int tdime = 1; tdime <= 3; tdime++){
    const intAr1 &lcent = tdime == 1 ? cav.lcedg : 
                          tdime == 2 ? cav.lcfac : cav.lctet;
    int ncent = lcent.get_n();
    intAr2 &ent2poi = msh.ent2poi(tdime);
    intAr2 &ent2tag = msh.ent2tag(tdime);
    // The degree is not a template argument -> non constexpr 
    int nnode = msh.nnode(tdime); 

    for(int ientl = 0; ientl < ncent; ientl++){
      int ientt = lcent[ientl];
      for(int ii = 0; ii < nnode; ii++){
        int ipoin = ent2poi(ientt,ii);
        int ientt = msh.poi2ent(ipoin,0);
        int tdimp = msh.poi2ent(ipoin,1); 
        if(tdimp != tdime) continue;
        if(ent2tag(ithread,ientt) < msh.tag[ithread]) continue;
        msh.poi2ent(ipoin,0) = -1;
        msh.poi2ent(ipoin,1) = -1;
      }
    }
  }

  //// 1.2 Add new bpois 
  //if(msh.isboundary_edges()){
  //  int nnode = edgnpps[msh.curdeg];
  //  for(int iedge = nedg0 ; iedge < msh.nedge; iedge++){
  //    for(int ii = 0; ii < nnode ;ii ++){
  //      int ip = msh.edg2poi(iedge,ii);
  //      METRIS_ASSERT(ip >= 0 && ip < msh.npoin);
  //      int ib = msh.poi2bpo[ip];
  //      // Normal rules don't apply as we just removed a bunch of bpois
  //      // If ib == -1, or if ib >= 0 and itype == 0, then insert this edge. 
  //      // NOTE: if the point is a corner, this info has not been lost! thus there is no
  //      // risk of missing to add an edge of a reference not represented. 
  //      if(ib >= 0){
  //        if(msh.bpo2ibi(ib,0) == -1) msh.poi2bpo[ip] = -1;
  //        if(msh.bpo2ibi(ib,1) == 0 || msh.bpo2ibi(ib,1) == 2 || msh.bpo2ibi(ib,0) == -1) ib = -1; 
  //      }
  //      if(ib == -1){
  //        ib = msh.newbpotopo(ip,1,iedge);
  //        if(ii < 2 && ip != cav.ipins){
  //          msh.bpo2rbi(ib,0) = redge[iedge-nedg0];
  //        }
  //      }
  //    }
  //  } 
  //}

  //if(msh.isboundary_faces()){
  //  int nnode = facnpps[msh.curdeg];
  //  for(int iface = nfac0 ; iface < msh.nface; iface++){
  //    int jj = 0;
  //    for(int ii = 0; ii < nnode ;ii ++){
  //      int ip = msh.fac2poi(iface,ii);
  //      METRIS_ASSERT(ip >= 0 && ip < msh.npoin);
  //      int ib = msh.poi2bpo[ip];
  //      // Same as previous, triangle goes in only if nothing yet (ib == -1)
  //      // or corner / edge (itype <= 1)
  //      if(ib >= 0){
  //        if(msh.bpo2ibi(ib,0) == -1) msh.poi2bpo[ip] = -1;
  //        if(msh.bpo2ibi(ib,1) <= 1 || msh.bpo2ibi(ib,0) == -1) ib = -1; 
  //      }
  //      if(ib == -1){
  //        ib = msh.newbpotopo(ip,2,iface);
  //        //msh.bpo2rbi(ib,0) = redge[iedge-nedg0];
  //        if(ii < 3 && ip != cav.ipins){
  //          msh.bpo2rbi(ib,0) = rface(iface-nfac0,jj,0);
  //          msh.bpo2rbi(ib,1) = rface(iface-nfac0,jj,1);
  //          jj++; // node in rface index
  //        }
  //      }
  //    }
  //  } 
  //}
  


  // Remove dead non manifold neighbours 
  for(int ii = 0; ii < ncedg; ii++){
    int iedge = cav.lcedg[ii];
    for(int jj = 0; jj < 2; jj++){
      int iedg2 = msh.edg2edg(iedge,jj);
      // >= 0 valid, -1 no neighbour
      if(iedg2 >= -1) continue; 
      // We're in a linked list, but we don't know previous. 
      // Get iedg1 = next, iedg0 = previous. 
      int iedg0 = -1, iedg1 = iedg2; 
      int ipoin = msh.edg2poi[iedge][1-jj];
      int iedgc = iedge, ineic = jj;
      int inei0;
      //int nnei = 0;
      while(getnextedgnm(msh,iedge,ipoin,&iedgc,&ineic)){
        iedg0 = iedgc;
        inei0 = ineic;
        //nnei++;
      }
      //METRIS_ASSERT(iedg0 != -1);
      if(iedg0 != -1){
        msh.edg2edg(iedg0,inei0) = iedg1; 
      }else{ // Case where we collapsed all the edges surrounding a corner
        // Meaning only iedg1 is left: simply set it to -1
        iedg1 = - iedg1 - 2; // This was never set positive
        if(msh.edg2poi(iedg1,0) == ipoin){
          msh.edg2edg(iedg1,1) = -1;
        }else{
          msh.edg2edg(iedg1,0) = -1;
        }
      }
      //// iedg0, iedg1 are the resp previous, next edges 
      //inei1 = msh.edg2poi[iedg1][1-0] == ipoin ? 0 : 
      //        msh.edg2poi[iedg1][1-1] == ipoin ? 1 : -1;
      //METRIS_ASSERT(inei1 != -1);
      //if(nnei > 2){
      //}else{ // ((((((No longer non manifold)))))) WRONG !! WRONG !! Topology does not change!!
      //  iedg1 = - iedg1 - 2; // This was never set positive
      //  msh.edg2edg(iedg0,inei0) = iedg1;
      //  if(msh.edg2poi(iedg1,0) == ipoin){
      //    msh.edg2edg(iedg1,1) = iedg0;
      //  }else{
      //    msh.edg2edg(iedg1,0) = iedg0;
      //  }
      //}
    }
  }


  // Update edge hash table and poi2ent 
  for(int iedge = nedg0; iedge < msh.nedge; iedge++){
    int ip[2]; 
    for(int ii = 0; ii < 2; ii++){
      ip[ii] = msh.edg2poi(iedge,ii);
      if(msh.poi2ent[ip[ii]][1] >= 1 || msh.poi2ent[ip[ii]][1] <= 0){
        msh.poi2ent[ip[ii]][0] = iedge;
        msh.poi2ent[ip[ii]][1] = 1;
      }
    }
    constexpr int nnode = edgnpps[ideg];
    for(int ii = 2; ii < nnode; ii++){
      int ipoin = msh.edg2poi(iedge,ii);
      if(msh.poi2ent[ipoin][1] >= 1 || msh.poi2ent[ipoin][1] <= 0){
        msh.poi2ent[ipoin][0] = iedge;
        msh.poi2ent[ipoin][1] = 1;
      }
    }

    auto key = stup2(ip[0],ip[1]);
    msh.edgHshTab.insert({key,iedge});
  }


  // Internal neighbours are updated directly in reconnect_lincav. 
  // External neighbours
  for(int iedge = nedg0; iedge < msh.nedge; iedge++){
    for(int inei = 0; inei < 2; inei++){
      int iedg2 = msh.edg2edg(iedge,inei); 
      if(iedg2 >= nedg0) continue;
      if(iedg2 == -1) continue;
      if(iedg2 < -1 && - iedg2 - 2 >= nedg0) continue;

      int ipoin = msh.edg2poi[iedge][1-inei];

      if(iedg2 >= 0){// Simple manifold neighbour.
        int ine2 = msh.edg2poi[iedg2][1-0] == ipoin ? 0 : 
                   msh.edg2poi[iedg2][1-1] == ipoin ? 1 : -1;
        msh.edg2edg(iedg2,ine2) = iedge;
      }else if(iedg2 < -1){ // Non manifold
        // What happens here is that iedge points to some loop, but is not part of that loop
        // Simply read two elements and sandwich it there
        int iedg0 = - iedg2 - 2;
        int ine0 = msh.edg2poi[iedg0][1-0] == ipoin ? 0 : 
                   msh.edg2poi[iedg0][1-1] == ipoin ? 1 : -1;
        int iedg1 = msh.edg2edg(iedg0,ine0);
        METRIS_ASSERT(iedg1 < -1);
        // Now we make iedg0 point to iedge and iedge to iedg1
        msh.edg2edg(iedg0,ine0) = - iedge - 2;
        msh.edg2edg(iedge,inei) = iedg1; // Note hasn't been flipped
      }
    }
  }

  //msh.tag[ithread]++;
  //// Internal neighbours. 
  //for(int iedge = nedg0; iedge < msh.nedge; iedge++){
  //  for(int inei = 0; inei < 2; inei++){
  //    int ip = msh.edg2poi[iedge][1-nei];
  //    if(msh.poi2tag(ithread,ip) > msh.tag[ithread]){
  //      // This is an internal neighbour
  //      int iedg2 = msh.poi2tag(ithread,ip) - msh.tag[ithread] - 1;
  //      // Mostly a check that poi2tag[ithread] is valid
  //      assert("Internal neighbour is a new edge" && iedg2 >= nedg0);
  //      #ifndef NDEBUG
  //        printf("    -- int nei / %d, iedg2 = %d (%d , %d)\n",ip,iedg2,msh.edg2poi(iedg2,0),msh.edg2poi(iedg2,1));
  //        fflush(stdout);
  //      #endif
  //      msh.edg2edg(iedge,inei) = iedg2;
  //      if(msh.edg2poi(iedg2,0) == ip){
  //        msh.edg2edg(iedg2,1) = iedge;
  //      }else if(msh.edg2poi(iedg2,1) == ip){
  //        msh.edg2edg(iedg2,0) = iedge;
  //      }else{
  //        printf("## ERROR EDGE NEIGHBOUR DOES NOT HAVE VERTEX \n");
  //        exit(1);
  //      }
  //      msh.poi2tag(ithread,ip) = msh.tag[ithread];
  //      continue;
  //    }
  //    if(msh.poi2tag(ithread,ip) == msh.tag[ithread]){
  //      // This means the point was seen twice already. 
  //      // This cannot happen before we have implemented non manifold point insertion (corner)
  //      printf("## CORNER CASE! TEST THOROUGHLY: IMPLEMENT REFS ! \n");
  //    }
  //    msh.poi2tag(ithread,ip) = msh.tag[ithread] + 1 + iedge;
  //    
  //    int ienei = msh.edg2edg(iedge,inei);
  //    if(ienei < 0) continue;

  //    #ifndef NDEBUG
  //      printf("    -- ext nei / %d, ienei = %d (%d , %d)\n",ip,ienei,msh.edg2poi(ienei,0),msh.edg2poi(ienei,1));
  //    #endif

  //    if(msh.edg2poi(ienei,0) == ip){
  //      msh.edg2edg(ienei,1) = iedge;
  //    }else if(msh.edg2poi(ienei,1) == ip){
  //      msh.edg2edg(ienei,0) = iedge;
  //    }else{
  //      printf("## ERROR EDGE NEIGHBOUR DOES NOT HAVE VERTEX\n");
  //      exit(1);
  //    }
  //  }
  //}


  //for(int iedge = nedg0; iedge < msh.nedge; iedge++){
  //  for(int i=0; i<2;i++){
  //    msh.poi2tag[ithread][msh.edg2poi(iedge,i)] = 0;
  //  }
  //}

  // Tag cavity faces
  for(int iface : cav.lcfac) msh.fac2tag(ithread,iface) = msh.tag[ithread];

 
  // Inform cavity neighbours that their neighbours are defunct. 
  // This is necessry because some edges (or facets) of the cavity boundary need 
  // not have created new elements !
  // Also update edg2fac when the fac in question is outside the cavity 
  for(int ifacl = 0; ifacl < ncfac; ifacl++){
    int iface = cav.lcfac[ifacl];
    for(int ied = 0; ied < 3; ied++){
      int ifac2 = msh.fac2fac(iface,ied);
      if(ifac2 == -1) continue;
      if(ifac2 >= 0 && msh.fac2tag(ithread,ifac2) >= msh.tag[ithread]) continue;
      if(ifac2 < -1 && msh.fac2tag[ithread][-ifac2-2] >= msh.tag[ithread]) continue;

      // Even if the edge was not originally boundary, it may be now!! 
      // This is precisely why we do this. 
      int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
      int ip2 = msh.fac2poi(iface,lnoed2[ied][1]);
      int iedge = getedgglo(msh,ip1,ip2);
      if(iedge >= 0){
        if(ifac2 >= 0) msh.edg2fac[iedge] = ifac2;
        else           msh.edg2fac[iedge] = - ifac2 - 2;
      }
    }
  }

  if(msh.nelem > 0) METRIS_THROW_MSG(TODOExcept(),
    "See dead triangle neighbour update -> tets nelem = "<<msh.nelem)




  // New face updates, namely neighbours 
  for(int ifanw = nfac0; ifanw < msh.nface; ifanw++){
    INCVDEPTH(msh);
    // poi2ent update
    int ip[3]; 
    for(int ii = 0; ii < 3; ii++){
      ip[ii] = msh.fac2poi(ifanw,ii);
      if(msh.poi2ent[ip[ii]][1] >= 2  || msh.poi2ent[ip[ii]][1] <= 0){
        msh.poi2ent[ip[ii]][0] = ifanw;
        msh.poi2ent[ip[ii]][1] = 2;
      }
    }
    constexpr int nnode = facnpps[ideg];
    for(int ii = 3; ii < nnode; ii++){
      int ipoin = msh.fac2poi(ifanw,ii);
      if(msh.poi2ent[ipoin][1] >= 2  || msh.poi2ent[ipoin][1] <= 0){
        msh.poi2ent[ipoin][0] = ifanw;
        msh.poi2ent[ipoin][1] = 2;
      }
    }
    // insert in hashtab 
    if(msh.get_tdim() == 3){
      auto key = stup3(ip[0],ip[1],ip[2]);
      msh.facHshTab.insert({key,ifanw});
    }

    bool ok = false;
    for(int ii = 0; ii < 3; ii++){
      // Not all cavity boundary edges are opposite ipins ! 
      // We need to fetch neighbour (which should always be init)
      // and see if new (update ok) or old
      int ifaca = msh.fac2fac(ifanw,ii);


      if(ifaca >= nfac0) METRIS_ASSERT(!isdeadent(ifaca,msh.fac2poi));
      if(ifaca >= nfac0) continue; // Internal neighbours are known.
      //if(msh.fac2poi(ifanw,ii) != cav.ipins && ifaca >= nfac0) continue;
      //if(msh.fac2poi(ifanw,ii) != cav.ipins && ifaca < 0) continue;


      ok = true;
      int jp1   = msh.fac2poi(ifanw,lnoed2[ii][0]);
      int jp2   = msh.fac2poi(ifanw,lnoed2[ii][1]);

      CPRINTF1(" - ifanw = %d ii = %d ineigh = %d jp1 %d jp2 %d \n",
               ifanw,ii,ifaca,jp1,jp2);

      if(ifaca < 0){ // Either surface boundary or non manifold
        int iedge = msh.facedg2glo(ifanw, ii); 

        CPRINTF1(" - edge neighbour = %d <? %d = nedg0 \n",iedge,nedg0);

        // Check valid AND not a new edge otherwise what happened??
        METRIS_ENFORCE_MSG(iedge >= 0,"ifanw = "<<ifanw<<" ifaca = "<<ifaca<<" iedg = "<<ii);

        if(iedge >= nedg0) continue; // New edge -> already handled

        METRIS_ASSERT(!isdeadent(iedge,msh.edg2poi));

        // Figure out if surface bdry or non manifold. 
        // Start by getting face attached to edge. 
        int ifaed = msh.edg2fac[iedge];
        METRIS_ASSERT(ifaed >= 0);
        CPRINTF1(" - already attached ifaed = %d dead ? %d \n",
                 ifaed,isdeadent(ifaed,msh.fac2poi));
        if(!isdeadent(ifaed,msh.fac2poi) && DOPRINTS1()){
          CPRINTF1(" - ifaed vertices: ");
          intAr1(facnpps[msh.curdeg],msh.fac2poi[ifaed]).print();
        }

        // If an old cavity element, or a new one already updated
        if(isdeadent(ifaed,msh.fac2poi)
        || msh.fac2tag(ithread,ifaed) >= msh.tag[ithread]){
          CPRINTF1(" - edg2fac link update iedge = %d : %d <- %d (new <- old)\n",
                    iedge,ifanw,msh.edg2fac[iedge]);
          msh.edg2fac[iedge] = ifanw;
          continue;
        }

        // In this case, ifaed is external to the cavity. 
        int ieed;
        try{
          ieed = getedgfac(msh,ifaed,jp1,jp2);
        }catch(const MetrisExcept &e){
          printf("### FATAL ERROR \n");
          printf("Initial ifanw = %d nodes = %d %d %d ; current edge = %d %d \n",
            ifanw,msh.fac2poi(ifanw,0),msh.fac2poi(ifanw,1),msh.fac2poi(ifanw,2),
            jp1,jp2);
          if(ifaca < -1){
            printf("First neighbour (ifaca) = %d nodes = %d %d %d \n",ifaca,
              msh.fac2poi[-ifaca-2][0],msh.fac2poi[-ifaca-2][1],msh.fac2poi[-ifaca-2][2]);
          }else{
            printf("Boundary -> no first neighbour (ifaca = -1)\n");
          }
          printf("glo edge = %d nodes %d %d \n",iedge,msh.edg2poi(iedge,0)
            ,msh.edg2poi(iedge,1));
          printf("Edge points to ifaed = %d nodes = %d %d %d \n",ifaed
            ,msh.fac2poi(ifaed,0),msh.fac2poi(ifaed,1),msh.fac2poi(ifaed,2));
          throw(e);
        }
        METRIS_ASSERT(ieed >= 0);
        // Neighbour to the edge of the face originally attached to edge
        int ifanm = msh.fac2fac(ifaed,ieed);
        CPRINTF1(" - original edge -> fac neighbour ifanm = %d \n",ifanm);

        // If there was no neighbour, then edge was pointed to by single face. 
        // This is the easiest case. 
        if(ifanm == -1){
          if(msh.fac2tag(ithread,ifaed) >= msh.tag[ithread]){
            msh.fac2fac(ifanw,ii) = -1; // Actually nothing to do here, for clarity and robustness. 
            msh.edg2fac[iedge] = ifanw;  // Update edge to face link. 
            CPRINTF1(" - edg2fac link update iedge = %d : %d <- %d (new <- old)\n",
                     iedge,ifanw,msh.edg2fac[iedge]);
          }else{
            // Edge was pointed to by a single face, but this face was (and still is!!) outside the cavity... 
            // the topology has changed, no no
            METRIS_THROW_MSG(TopoExcept(),"Surface topology changed!!")
          }
          // Only one to link to, no further info needed from previous face, do update now
          msh.edg2fac[iedge] = ifanw;
          CPRINTF1(" - edg2fac link update iedge = %d : %d <- %d (new <- old)\n",
                   iedge,ifanw,msh.edg2fac[iedge]);
        }else{ // Non manifold case, leave to future self
          METRIS_THROW_MSG(TODOExcept(),"Implement non manifold neighbour update cavity")
        }

      }else{ // Nicely manifold. 
        int ied2 = getedgfac(msh,ifaca,jp1,jp2);
        METRIS_ASSERT(ied2 >= 0 && ied2 < 3);
        METRIS_ASSERT(!isdeadent(ifanw,msh.fac2poi));
        METRIS_ASSERT(!isdeadent(ifaca,msh.fac2poi));
        msh.fac2fac(ifaca,ied2) = ifanw;
        // Note, even if glo edge sandwiched, no need to update as neighbour is exterior, edg2fac remains valid. 
      }
    }
    METRIS_ASSERT(ok == true);
  }


  // Edge case... quite litterally: 2 faces one = the other (flipped) become 1 edge
  if(msh.nface == nfac0 && ncfac > 0){
    if(msh.nedge != nedg0 + 1){
      printf("## ERROR cavity ipins = %d \n",cav.ipins);
      if(cav.lcedg.get_n() > 0){
        printf("## Edge cavity (%d): ",cav.lcedg.get_n());
        for(int ii = 0; ii < cav.lcedg.get_n(); ii++){
          printf("%d : %d = ",ii,cav.lcedg[ii]);
          intAr1(edgnpps[ideg],msh.edg2poi[ii]).print();
        }
      }
      if(cav.lcfac.get_n() > 0){
        printf("## Face cavity (%d): ",cav.lcfac.get_n());
        for(int ii = 0; ii < cav.lcfac.get_n(); ii++){
          printf("%d : %d = ",ii,cav.lcfac[ii]);
          intAr1(facnpps[ideg],msh.fac2poi[ii]).print();
        }
      }
      if(cav.lctet.get_n() > 0){
        printf("## tet  cavity (%d): ",cav.lctet.get_n());
        for(int ii = 0; ii < cav.lctet.get_n(); ii++){
          printf("%d : %d = ",ii,cav.lctet[ii]);
          intAr1(tetnpps[ideg],msh.tet2poi[ii]).print();
        }
      }
      if(msh.nedge > nedg0){
        printf("## nedg0 %d nedge %d\n",nedg0,msh.nedge);
        for(int ii = nedg0; ii < msh.nedge; ii++){
          printf("%d = ",ii);
          intAr1(edgnpps[ideg],msh.edg2poi[ii]).print();
        }
      }
      if(msh.nface > nfac0){
        printf("## nfac0 %d nface %d\n",nfac0,msh.nface);
        for(int ii = nfac0; ii < msh.nface; ii++){
          printf("%d = ",ii);
          intAr1(facnpps[ideg],msh.fac2poi[ii]).print();
        }
      }
      if(msh.nelem > nele0){
        printf("## nele0 %d nelem %d\n",nele0,msh.nelem);
        for(int ii = nele0; ii < msh.nelem; ii++){
          printf("%d = ",ii);
          intAr1(tetnpps[ideg],msh.tet2poi[ii]).print();
        }
      }
      writeMesh("fatal",msh);
      writeMeshCavity("fatal.cav",msh,cav,0);
      METRIS_THROW_MSG(TODOExcept(), 
        "Not the 2 face -> 1 edge case, inspect. msh.nedge = "
        <<msh.nedge<<" nedg0 = "<<nedg0);
    }
    // Basically these are the -2 typ entries in edent/edtyp reconnect_faccav
    // Let's just recompute it
    int ifac0 = -1;
    int nedex = 0;
    for(int ifacl = 0; ifacl < ncfac; ifacl++){
      int iface = cav.lcfac[ifacl];
      for(int ied = 0; ied < 3; ied++){
        int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
        int ip2 = msh.fac2poi(iface,lnoed2[ied][1]);

        //int ipoin;
        //if(ip1 == cav.ipins){
        //  ipoin = ip2;
        //}else if(ip2 == cav.ipins){
        //  ipoin = ip1;
        //}else{
        //  continue;
        //}
        if(ip1 != cav.ipins && ip2 != cav.ipins) continue;


        int ifac2 = msh.fac2fac(iface,ied);
        if(ifac2 >= 0){
          METRIS_ASSERT(!isdeadent(ifac2,msh.fac2poi));
          if(msh.fac2tag(ithread,ifac2) >= msh.tag[ithread]) continue;
          // Standard neighbour, exterior face

          if(nedex >= 2) METRIS_THROW(SMemExcept());

          if(nedex == 0){
            nedex++;
            ifac0 = ifac2;
          }else{
            nedex++;
            int ifac1 = ifac2; 
            int ied1 = getedgfac(msh,ifac1,ip1,ip2);
            int ied0 = getedgfac(msh,ifac0,ip1,ip2);
            METRIS_ASSERT(ied0 >= 0);
            msh.fac2fac(ifac0,ied0) = ifac1;
            msh.fac2fac(ifac1,ied1) = ifac0;
          }
        }

      }
    }
  }



 
  // Inform cavity neighbours that their neighbours are defunct. 
  // This is necessry because some edges (or facets) of the cavity boundary need not have created new elements !
  // Also update edg2fac when the fac in question is outside the cavity 
  for(int ifacl = 0; ifacl < ncfac; ifacl++){
    int iface = cav.lcfac[ifacl];
    for(int ied = 0; ied < 3; ied++){
      int ifac2 = msh.fac2fac(iface,ied);
      if(ifac2 == -1) continue;

      // Even if the edge was not originally boundary, it may be now!! This is precisely why we do this. 
      int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
      int ip2 = msh.fac2poi(iface,lnoed2[ied][1]); 

      if(ifac2 >= 0){
        if(msh.fac2tag(ithread,ifac2) >= msh.tag[ithread]) continue;

        int ied2 = getedgfac(msh,ifac2,ip1,ip2);
        METRIS_ASSERT(ied2 >= 0);

        if(msh.fac2fac(ifac2,ied2) == iface) msh.fac2fac(ifac2,ied2) = -1;

      }else{
        METRIS_THROW_MSG(TODOExcept(),"Investigate this case");
      }
    }
  }

  if(msh.nelem > 0) METRIS_THROW_MSG(TODOExcept(),
        "See dead triangle neighbour update -> tets nelem = "<<msh.nelem)




  // poi2ent updates. Order matters (lowest tdim prevails)
  // First reset link for cavity vertices 

//  for(int ielem = nele0; ielem < msh.nelem; ielem++){
//    for(int ii = 0; ii < tetnpps[msh.curdeg]; ii++){
//      int ipoin = msh.tet2poi(ielem,ii);
//      msh.poi2ent[ipoin] = ielem;
//    }
//  }
//  for(int iface = nfac0; iface < msh.nface; iface++){
//    for(int ii = 0; ii < facnpps[msh.curdeg]; ii++){
//      int ipoin = msh.fac2poi(iface,ii);
//      msh.poi2ent[ipoin] = iface;
//    }
//  }
//  for(int iedge = nedg0; iedge < msh.nedge; iedge++){
//    for(int ii = 0; ii < edgnpps[msh.curdeg]; ii++){
//      int ipoin = msh.edg2poi(iedge,ii);
//      msh.poi2ent[ipoin] = iedge;
//    }
//  }


  if(msh.isboundary_edges()){
    for(int iedgl = 0; iedgl < ncedg; iedgl++){
      int iedge = cav.lcedg[iedgl];
      int nnode = edgnpps[msh.curdeg]; 
      for(int ii = 0; ii < nnode; ii++){
        int ipoin = msh.edg2poi(iedge,ii);
        if(msh.poi2ent(ipoin,0) != -1) continue;
        int ibpoi = msh.poi2bpo[ipoin];
        while(ibpoi >= 0){
          msh.bpo2ibi(ibpoi,0) = -1;
          ibpoi = msh.bpo2ibi(ibpoi,3);
        }
        msh.poi2bpo[ipoin] = -1;
      }
    }
  }
  if(msh.isboundary_faces()){
    for(int ifacl = 0; ifacl < ncfac; ifacl++){
      int iface = cav.lcfac[ifacl];
      int nnode = facnpps[msh.curdeg]; 
      for(int ii = 0; ii < nnode; ii++){
        int ipoin = msh.fac2poi(iface,ii);
        if(msh.poi2ent(ipoin,0) != -1) continue;
        int ibpoi = msh.poi2bpo[ipoin];
        while(ibpoi >= 0){
          msh.bpo2ibi(ibpoi,0) = -1;
          ibpoi = msh.bpo2ibi(ibpoi,3);
        }
        msh.poi2bpo[ipoin] = -1;
      }
    }
  }

//  { // new HO point updates: metric and CAD
//    // To avoid polluting the scope
//  int nentt = msh.nelem > 0 ? msh.nelem : 
//              msh.nface > 0 ? msh.nface : 
//                              msh.nedge;
//  int nent0 = msh.nelem > 0 ? nele0 : 
//              msh.nface > 0 ? nfac0 : 
//                              nedg0;
//  int nnode = msh.nelem > 0 ? tetnpps[msh.curdeg] : 
//              msh.nface > 0 ? facnpps[msh.curdeg] : 
//                              edgnpps[msh.curdeg];
//  intAr2& ent2poi = msh.nelem > 0 ? msh.tet2poi : 
//                    msh.nface > 0 ? msh.fac2poi : 
//                                    msh.edg2poi;
//  for(int ientt = nent0; ientt < nentt; ientt++){
//    for(int ii = 0; ii < nnode; ii++){
//      int ipoin = ent2poi(ientt,ii);
//      // Old points need this update because elements changed
//      msh.poi2ent[ipoin] = ientt;
//
//      // However only new pts need a metric / CAD update
//      if(ipoin < npoi0) continue;
//    }
//  }
//
//
//
//  } // poi2ent update

  if(msh.isboundary_faces() && msh.curdeg > 1 && msh.CAD()){ 
    // No updates necessary otherwise
    METRIS_ASSERT(cav.nrmal != NULL);

    for(int ipoin = npoi0; ipoin < msh.npoin; ipoin++){
      // regenerate (u,v)s from t
      int ib = msh.poi2bpo[ipoin];
      int ityp = msh.bpo2ibi(ib,1);
      if(ityp >= 2) continue;

      // Cannot have created a corner. 
      METRIS_ASSERT(ityp == 1);

      // Material is to the left of the curve when positively walked
      // Meaning if the normal is going towards us, edge is going up, 
      // then triangles to the left that share the edge, have 
      // the height (or the other edges) with positive scalar product 
      // to the vector product between normal and edge... 

      // edge = e 
      // normal = n
      // h = project 3rd triangle vertex (not on edge) -> 3rd triangle vertex
      // compute vec(n,e) = l
      // if (l,h) > 0 triangle is left -> positive
      //          < 0 triangle is right -> negative 

      double t = msh.bpo2rbi(ib,0);

      int iedge = msh.bpo2ibi(ib,2);
      METRIS_ASSERT(iedge >= 0);

      int irefe = msh.edg2ref[iedge];
      ego EG_edg = msh.CAD.cad2edg[irefe];
      METRIS_ASSERT(EG_edg != NULL);


      double vec[3];

      if constexpr(ideg >= 2){
        int jj = -1;
        for(int ii = 2; ii < edgnpps[ideg]; ii++){
          int ip = msh.edg2poi(iedge,ii);
          if(ip == ipoin){
            jj = ii;
            break;
          }
        }
        METRIS_ASSERT(jj>=0);
    
        double bar1[2];
        bar1[0] = ordedg.s[ideg][jj][0] / (ideg * 1.0);
        bar1[1] = ordedg.s[ideg][jj][1] / (ideg * 1.0);
      
        double dum[3], tang[3];
        eval1<3,ideg>(msh.coord,msh.edg2poi[iedge],msh.getBasis(),
                      DifVar::Bary,DifVar::None,bar1,dum,tang,NULL);

        // This vector is the definition of left
        vecprod(cav.nrmal,tang,vec);
      }



      ib = msh.bpo2ibi(ib,3);
      METRIS_ASSERT(ib >= 0);
      do{
        if(msh.bpo2ibi(ib,1) == 2){
          // Found a triangle, now update. 
          int iface = msh.bpo2ibi(ib,2);
          int ireff = msh.fac2ref[iface];
          ego EG_fac  = msh.CAD.cad2fac[ireff];

          METRIS_ASSERT(EG_fac != NULL);

          // Material is left if edge + 
          // Check to other triangle vertices. 
          // Note triangle may or may not share a full edge so just take max dtprd
          // For robustness take both min and max and keep highest abs value
          double dtprm =  1.0;
          double dtprM = -1.0;
          double du[3];
          for(int ii = 0; ii < 3; ii++){
            int ip = msh.fac2poi(iface,ii);
            if(ip == ipoin) continue;
            for(int jj = 0; jj < 3; jj++) du[jj] = 
              msh.coord(ip,jj) - msh.coord(ipoin,jj);
            double dtprd = getprdl2<3>(du,vec);
            if(dtprd > dtprM) dtprM = dtprd;
            if(dtprd < dtprm) dtprm = dtprd;
          }

          // If both dtprM and dtprm > 0, then |dtprM| > |dtprm|
          // If both dtprM and dtprm < 0, then |dtprm| > |dtprM|

          // If dtprM > 0 and dtprm < 0, we assume this is roff error 
          // of the smallest in absolute value. 
          // In this case, if |dtprM| > |dtprm|, then dtprM > 0 and dtprm > 0 (in reality)

          // In summary, the test is |dtprM| > |dtprm| <=> dtprM > 0 && dtprm > 0
          int sg = 1;
          if(abs(dtprM) < abs(dtprm)){
            sg = -1;
          }

          int icode = EG_getEdgeUV(EG_fac, EG_edg, sg, t, msh.bpo2rbi[ib]);

          if(icode != 0) METRIS_THROW_MSG(GeomExcept(),"EG_getEdgeUV failed !");

        }
        ib = msh.bpo2ibi(ib,3);
      }while(ib >= 0);




    } 
  }


  // Deleting entities in the initial cavity. 
  for(int iedgl = 0; iedgl < ncedg; iedgl++){
    int iedge = cav.lcedg[iedgl];
    //    msh.edg2tag(ithread,iedge) = msh.tag[ithread];

    //int ip1 = msh.edg2poi(iedge,0);
    //int ip2 = msh.edg2poi(iedge,1);
    //auto key = stup2(ip1,ip2);
    // This removes if exists
    killent(iedge, msh.edg2poi);
  }
  for(int ifacl = 0; ifacl < ncfac; ifacl++){
    int iface = cav.lcfac[ifacl];
    //    msh.fac2tag(ithread,iface) = msh.tag[ithread];

    //int ip1 = msh.fac2poi(iface,0);
    //int ip2 = msh.fac2poi(iface,1);
    //int ip3 = msh.fac2poi(iface,2);
    //auto key = stup3(ip1,ip2,ip3);
    // This removes if exists
    killent(iface, msh.fac2poi);
  }
  for(int ielel = 0; ielel < nctet; ielel++){
    int ielem = cav.lctet[ielel];
    //    msh.tet2tag(ithread,ielem) = msh.tag[ithread];
    killent(ielem, msh.tet2poi);
  }


  cleanup:
  msh.tag[ithread]++; // Because points were tag+1
  return 0;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template int update_cavity<MetricFieldAnalytical,n>(Mesh<MetricFieldAnalytical> &msh,\
                            const MshCavity &cav, const CavWrkArrs &work, int npoi0, int nedg0, \
                            int nfac0, int nele0, \
                            int ithread);\
template int update_cavity<MetricFieldFE        ,n>(Mesh<MetricFieldFE        > &msh, \
                            const MshCavity &cav, const CavWrkArrs &work, int npoi0, int nedg0, \
                            int nfac0, int nele0, \
                            int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()









} // end namespace