//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../utils/aux_misc.hxx"
#include "../ho_constants.hxx"
#include "../io_libmeshb.hxx"
#include "../utils/CT_loop.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"
#include "../linalg/det.hxx"
#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../low_geo/misc.hxx"

namespace Metris{





// redge and rface hold the uvs as provided by reconnect_ routines when new faces are created
// this is because we lose the link between which cavity element created which element. 
// they are dimensioned for nface - nfac0, nedge - nedg0.
// redge has info for the point that is not ipins
// rface is dblAr3 as 2 nodes, each 2 parameters
template<class MFT, int ideg>
int update_cavity(Mesh<MFT> &msh, MshCavity &cav, [[maybe_unused]] const CavWrkArrs &work,
                  [[maybe_unused]] int npoi0, int nedg0, int nfac0, int nele0, 
                  int ithread){

  int ierro = 0;
  GETVDEPTH(msh.param);

  if(msh.param->dbgfull){
    bool ok[3] = {true, true, true};
    for(int tdim = 1; tdim <= 3; tdim++){
      for(int ientt : cav.lcent(tdim)){
        if(msh.ent2tag(tdim)(ithread,ientt) >= msh.tag[ithread]) continue;
        ok[tdim-1] = false;
        break;
      }
    }
    if(!ok[0] || !ok[1] || !ok[2]){
      PRINTF("## Untagged cavity entities tag = {} \n",msh.tag[ithread]);
      for(int tdim = 1; tdim <= 3; tdim++){
        for(int ientt : cav.lcent(tdim)){
          PRINTF("tdim {} ientt {} tag {} \n",tdim,ientt,msh.ent2tag(tdim)(ithread,ientt));
        }
      }
      METRIS_THROW();
    }
  }

  const int ncedg = cav.lcedg.get_n();
  const int ncfac = cav.lcfac.get_n();
  const int nctet = cav.lctet.get_n();
  HshTab_I3I facHsh; 

  if(DOPRINTS2()){
    MshCavity cav2(msh.nelem-nele0,msh.nedge-nedg0,msh.nface-nfac0);
    for(int ii = nele0; ii < msh.nelem; ii++) cav2.lctet.stack(ii);
    for(int ii = nfac0; ii < msh.nface; ii++) cav2.lcfac.stack(ii);
    for(int ii = nedg0; ii < msh.nedge; ii++) cav2.lcedg.stack(ii);
    cav2.ipins = cav.ipins;
    writeMeshCavity("cavity2",msh,cav2);
    cav2.print(msh);
  }

  // -- 1 Manage bpois


  if(msh.param->dbgfull){
    bool ok[3] = {true, true, true};
    for(int tdim = 1; tdim <= 3; tdim++){
      for(int ientt : cav.lcent(tdim)){
        if(msh.ent2tag(tdim)(ithread,ientt) >= msh.tag[ithread]) continue;
        ok[tdim-1] = false;
        break;
      }
    }
    if(!ok[0] || !ok[1] || !ok[2]){
      PRINTF("## 1 Untagged cavity entities tag = {} \n",msh.tag[ithread]);
      for(int tdim = 1; tdim <= 3; tdim++){
        for(int ientt : cav.lcent(tdim)){
          PRINTF("tdim {} ientt {} tag {} \n",tdim,ientt,msh.ent2tag(tdim)(ithread,ientt));
        }
      }
      METRIS_THROW();
    }
  }
  // -- 1.1 remove old ibpois

  //bool dowait = false;
  int ptag0 = msh.tag[ithread] + 1;
  cav.maxtag = MAX(cav.maxtag,ptag0);
  CT_FOR0_INC(1,2,tdimn){
    if (msh.isboundary_tdim(tdimn)){
    int nnode = msh.nnode(tdimn);
    for(int ientt : cav.lcent<tdimn>()){
      INCVDEPTH(msh.param);
      for(int ii = 0; ii < nnode; ii++){
        int ip = msh.template ent2poi<tdimn>()(ientt,ii);
        if(msh.poi2tag(ithread,ip) >= ptag0) continue;
        msh.poi2tag(ithread,ip) = ptag0;
        if(DOPRINTS3()){
          CPRINTF3(" - ip = {} clean bpo pre:\n",ip);
          print_bpolist(msh,msh.poi2bpo[ip]);
        }
        msh.rembpotag(ip,ithread);
        if(DOPRINTS3()){
          CPRINTF3(" - bpo post:\n");
          print_bpolist(msh,msh.poi2bpo[ip]);
        }
      }
    }
    }
  }CT_FOR1(tdimn);

  // ipins may be on a "virtual" dead edge/face 
  if(msh.poi2bpo[cav.ipins] >= 0) msh.rembpotag(cav.ipins,ithread);

  if(msh.param->dbgfull){
    bool ok[3] = {true, true, true};
    for(int tdim = 1; tdim <= 3; tdim++){
      for(int ientt : cav.lcent(tdim)){
        if(msh.ent2tag(tdim)(ithread,ientt) >= msh.tag[ithread]) continue;
        ok[tdim-1] = false;
        break;
      }
    }
    if(!ok[0] || !ok[1] || !ok[2]){
      PRINTF("## 2 Untagged cavity entities tag = {} \n",msh.tag[ithread]);
      for(int tdim = 1; tdim <= 3; tdim++){
        for(int ientt : cav.lcent(tdim)){
          PRINTF("tdim {} ientt {} tag {} \n",tdim,ientt,msh.ent2tag(tdim)(ithread,ientt));
        }
      }
      METRIS_THROW();
    }
  }

  // During reconnect_faccav and reconnect_lincav, we called newbpotopo on each
  // new entity. Now, we need to clean up. This was necessary as not stacking any
  // new bpos might have made a point unlinked after old elt deletion. 
  ptag0++;
  cav.maxtag = MAX(cav.maxtag,ptag0);
  for(int tdim : {1,2}){
    // Edges are probably always boundary but still.
    if(!msh.isboundary_tdim(tdim)) continue;
    int nnode = msh.nnode(tdim);
    int nentt = msh.nentt(tdim);
    int nent0 = tdim == 1 ? nedg0 : nfac0;
    intAr2 &ent2poi = msh.ent2poi(tdim);
    for(int ientt = nent0; ientt < nentt; ientt++){
      INCVDEPTH(msh.param);
      for(int ii = 0; ii < nnode; ii++){
        int ip = ent2poi(ientt,ii);
        if(msh.poi2tag(ithread,ip) >= ptag0) continue;
        msh.poi2tag(ithread,ip) = ptag0;

        CPRINTF2(" - update bpo tdim {} ientt {} ipoin {}\n",tdim,ientt,ip);
        // First link and minimum topo dimn
        int ibpo0 = msh.poi2bpo[ip];
        METRIS_ASSERT(ibpo0 >= 0); 
        int mindi = msh.bpo2ibi(ibpo0,1);

        // Init straight at next and current dimn
        int ibpoi = msh.bpo2ibi(ibpo0,3);
        if(ibpoi >= 0){
          int curdi = msh.bpo2ibi(ibpoi,1);
          while(curdi <= mindi){
            METRIS_ASSERT(curdi == mindi);
            int ibpon = msh.bpo2ibi(ibpoi,3);
            msh.bpo2ibi(ibpo0,3) = ibpon;
            msh.bpo2ibi(ibpoi,0) = -1; // Kill the ibpoi left hanging.
            ibpoi = ibpon;
            if(ibpoi < 0) break; 
            curdi = msh.bpo2ibi(ibpoi,1);
          }
        }//endif(ibpoi >= 0)

      }
    }
  }


  if(msh.param->dbgfull){
    bool ok[3] = {true, true, true};
    for(int tdim = 1; tdim <= 3; tdim++){
      for(int ientt : cav.lcent(tdim)){
        if(msh.ent2tag(tdim)(ithread,ientt) >= msh.tag[ithread]) continue;
        ok[tdim-1] = false;
        break;
      }
    }
    if(!ok[0] || !ok[1] || !ok[2]){
      PRINTF("## 3 Untagged cavity entities tag = {} \n",msh.tag[ithread]);
      for(int tdim = 1; tdim <= 3; tdim++){
        for(int ientt : cav.lcent(tdim)){
          PRINTF("tdim {} ientt {} tag {} \n",tdim,ientt,msh.ent2tag(tdim)(ithread,ientt));
        }
      }
      METRIS_THROW();
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
      INCVDEPTH(msh.param);
      int ientt = lcent[ientl];
      for(int ii = 0; ii < nnode; ii++){
        int ipoin = ent2poi(ientt,ii);
        int ientt = msh.poi2ent(ipoin,0);
        if(ientt < -1) ientt = - ientt - 2; // control point
        int tdimp = msh.poi2ent(ipoin,1); 
        if(tdimp != tdime) continue;
        if(ent2tag(ithread,ientt) < msh.tag[ithread]) continue;
        msh.set_poi2ent(Vertex{ipoin}, -1, -1); // type doesn't matter for -1
      }
    }
  }


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
      int ipoin = msh.edg2poi(iedge,1-jj);
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
      if(msh.poi2ent(ip[ii],1) >= 1 || msh.poi2ent(ip[ii],1) <= 0){
        msh.set_poi2ent(Vertex{ip[ii]}, 1, iedge);
      }
    }
    constexpr int nnode = getnnod1(ideg);
    for(int ii = 2; ii < nnode; ii++){
      int ipoin = msh.edg2poi(iedge,ii);
      if(msh.poi2ent(ipoin,1) >= 1 || msh.poi2ent(ipoin,1) <= 0){
        msh.set_poi2ent(CtrlPt{ipoin}, 1, iedge);
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

      int ipoin = msh.edg2poi(iedge,1-inei);

      if(iedg2 >= 0){// Simple manifold neighbour.
        int ine2 = msh.edg2poi(iedg2,1-0) == ipoin ? 0 : 
                   msh.edg2poi(iedg2,1-1) == ipoin ? 1 : -1;
        msh.edg2edg(iedg2,ine2) = iedge;
      }else if(iedg2 < -1){ // Non manifold
        // What happens here is that iedge points to some loop, but is not part of that loop
        // Simply read two elements and sandwich it there
        int iedg0 = - iedg2 - 2;
        int ine0 = msh.edg2poi(iedg0,1-0) == ipoin ? 0 : 
                   msh.edg2poi(iedg0,1-1) == ipoin ? 1 : -1;
        int iedg1 = msh.edg2edg(iedg0,ine0);
        METRIS_ASSERT(iedg1 < -1);
        // Now we make iedg0 point to iedge and iedge to iedg1
        msh.edg2edg(iedg0,ine0) = - iedge - 2;
        msh.edg2edg(iedge,inei) = iedg1; // Note hasn't been flipped
      }
    }
  }
 
  // Inform cavity neighbours that their neighbours are defunct. 
  // This is necessry because some edges (or facets) of the cavity boundary may 
  // not have created new elements !
  // Also update edg2fac when the fac in question is outside the cavity 
  for(int iface : cav.lcfac){
    for(int ied = 0; ied < 3; ied++){
      int ifac2 = msh.fac2fac(iface,ied);
      if(ifac2 == -1) continue;
      if(ifac2 >= 0 && msh.fac2tag(ithread,ifac2) >= msh.tag[ithread]) continue;
      if(ifac2 < -1 && msh.fac2tag(ithread,-ifac2-2) >= msh.tag[ithread]) continue;

      // ifac2 is a non-cavity element, this edge is on cavity bdry

      // Even if the edge was not originally boundary, it may be now!
      int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
      int ip2 = msh.fac2poi(iface,lnoed2[ied][1]);
      int iedge = getedgglo(msh,ip1,ip2);
      if(iedge >= 0){
        if(ifac2 >= 0) msh.edg2fac[iedge] = ifac2;
        else           msh.edg2fac[iedge] = - ifac2 - 2;
      }
    }
  }



  // Remove old faces from hashtable
  for(int iface : cav.lcfac){
    int ip1 = msh.fac2poi(iface,0);
    int ip2 = msh.fac2poi(iface,1);
    int ip3 = msh.fac2poi(iface,2);
    auto key = stup3(ip1,ip2,ip3);
    msh.facHshTab.erase(key);
  }

  // New face updates, namely neighbours 
  for(int ifanw = nfac0; ifanw < msh.nface; ifanw++){
    INCVDEPTH(msh.param);
    // poi2ent update
    int ip[3]; 
    for(int ii = 0; ii < 3; ii++){
      ip[ii] = msh.fac2poi(ifanw,ii);
      if(msh.poi2ent(ip[ii],1) >= 2  || msh.poi2ent(ip[ii],1) <= 0){
        msh.set_poi2ent(Vertex{ip[ii]}, 2, ifanw);
      }
    }
    constexpr int nnode = getnnod2(ideg);
    for(int ii = 3; ii < nnode; ii++){
      int ipoin = msh.fac2poi(ifanw,ii);
      if(msh.poi2ent(ipoin,1) >= 2  || msh.poi2ent(ipoin,1) <= 0){
        msh.set_poi2ent(CtrlPt{ipoin}, 2, ifanw);
      }
    }
    // insert in hashtab 
    if(msh.idim >= 3){
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

      CPRINTF1(" - ifanw = {} ii = {} ineigh = {} jp1 {} jp2 {} \n",
               ifanw,ii,ifaca,jp1,jp2);

      if(ifaca < 0){ // Either surface boundary or non manifold
        int iedge = msh.facedg2glo(ifanw, ii); 

        CPRINTF1(" - edge neighbour = {} <? {} = nedg0 \n",iedge,nedg0);

        // Check valid AND not a new edge otherwise what happened??
        METRIS_ENFORCE_MSG(iedge >= 0,"ifanw = {} ifaca = {} iedg = {}", ifanw, ifaca, ii);

        if(iedge >= nedg0) continue; // New edge -> already handled

        METRIS_ASSERT(!isdeadent(iedge,msh.edg2poi));

        // Figure out if surface bdry or non manifold. 
        // Start by getting face attached to edge. 
        int ifaed = msh.edg2fac[iedge];
        METRIS_ASSERT(ifaed >= 0);
        CPRINTF1(" - already attached ifaed = {} dead ? {} \n",
                 ifaed,isdeadent(ifaed,msh.fac2poi));
        if(!isdeadent(ifaed,msh.fac2poi) && DOPRINTS1()){
          CPRINTF1(" - ifaed vertices: {}\n",intAr1(getnnod2(ideg),msh.fac2poi[ifaed]));
        }

        // If an old cavity element, or a new one already updated
        if(isdeadent(ifaed,msh.fac2poi)
        || msh.fac2tag(ithread,ifaed) >= msh.tag[ithread]){
          CPRINTF1(" - edg2fac link update iedge = {} : {} <- {} (new <- old)\n",
                    iedge,ifanw,msh.edg2fac[iedge]);
          msh.edg2fac[iedge] = ifanw;
          continue;
        }

        // In this case, ifaed is external to the cavity. 
        int ieed;
        try{
          ieed = getedgfac(msh,ifaed,jp1,jp2);
        }catch(const MetrisExcept &e){
          PRINTF("### FATAL ERROR \n");
          PRINTF("Initial ifanw = {} nodes = {} {} {} ; current edge = {} {} \n",
            ifanw,msh.fac2poi(ifanw,0),msh.fac2poi(ifanw,1),msh.fac2poi(ifanw,2),
            jp1,jp2);
          if(ifaca < -1){
            PRINTF("First neighbour (ifaca) = {} nodes = {} {} {} \n",ifaca,
              msh.fac2poi(-ifaca-2,0),msh.fac2poi(-ifaca-2,1),msh.fac2poi(-ifaca-2,2));
          }else{
            PRINTF("Boundary -> no first neighbour (ifaca = -1)\n");
          }
          PRINTF("glo edge = {} nodes {} {} \n",iedge,msh.edg2poi(iedge,0)
            ,msh.edg2poi(iedge,1));
          PRINTF("Edge points to ifaed = {} nodes = {} {} {} \n",ifaed
            ,msh.fac2poi(ifaed,0),msh.fac2poi(ifaed,1),msh.fac2poi(ifaed,2));
          throw(e);
        }
        METRIS_ASSERT(ieed >= 0);
        // Neighbour to the edge of the face originally attached to edge
        int ifanm = msh.fac2fac(ifaed,ieed);
        CPRINTF1(" - original edge -> fac neighbour ifanm = {} \n",ifanm);

        // If there was no neighbour, then edge was pointed to by single face. 
        // This is the easiest case. 
        if(ifanm == -1){
          if(msh.fac2tag(ithread,ifaed) >= msh.tag[ithread]){
            msh.fac2fac(ifanw,ii) = -1; // Actually nothing to do here, for clarity and robustness. 
            msh.edg2fac[iedge] = ifanw;  // Update edge to face link. 
            CPRINTF1(" - edg2fac link update iedge = {} : {} <- {} (new <- old)\n",
                     iedge,ifanw,msh.edg2fac[iedge]);
          }else{
            // Edge was pointed to by a single face, but this face was (and still is!!) outside the cavity... 
            // the topology has changed, no no
            METRIS_THROW_MSG("Surface topology changed!!")
          }
          // Only one to link to, no further info needed from previous face, do update now
          msh.edg2fac[iedge] = ifanw;
          CPRINTF1(" - edg2fac link update iedge = {} : {} <- {} (new <- old)\n",
                   iedge,ifanw,msh.edg2fac[iedge]);
        }else{ // Non manifold case, leave to future self
          METRIS_THROW_MSG("TODO: Implement non manifold neighbour update cavity")
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
      PRINTF("## ERROR cavity ipins = {} \n",cav.ipins);
      PRINTF("Cavity:\n");
      cav.print(msh,10);
      if(msh.nedge > nedg0){
        PRINTF("## nedg0 {} nedge {}\n",nedg0,msh.nedge);
        for(int ii = nedg0; ii < msh.nedge; ii++){
          PRINTF("{} = {}\n",ii,intAr1(getnnod1(ideg),msh.edg2poi[ii]));
        }
      }
      if(msh.nface > nfac0){
        PRINTF("## nfac0 {} nface {}\n",nfac0,msh.nface);
        for(int ii = nfac0; ii < msh.nface; ii++){
          PRINTF("{} = {}\n",ii,intAr1(getnnod2(ideg),msh.fac2poi[ii]));
        }
      }
      if(msh.nelem > nele0){
        PRINTF("## nele0 {} nelem {}\n",nele0,msh.nelem);
        for(int ii = nele0; ii < msh.nelem; ii++){
          PRINTF("{} = {}\n",ii,intAr1(getnnod3(ideg),msh.tet2poi[ii]));
        }
      }
      writeMesh("fatal",msh);
      writeMeshCavity("fatal.cav",msh,cav);
      METRIS_THROW_MSG( 
        "Not the 2 face -> 1 edge case, inspect. msh.nedge = {}, nedg0 = {}",
        msh.nedge, nedg0);
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

          METRIS_ENFORCE(nedex == 2);

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
  for(int iface : cav.lcfac){
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
        METRIS_THROW_MSG("TODO: Investigate this case");
      }
    }
  }

  //if(msh.nelem > 0) METRIS_THROW_MSG(
  //      "See dead triangle neighbour update -> tets nelem = "<<msh.nelem)

  // Next two steps are to update exterior neighbours and fac2tet. 
  // Tetra internal neighbour relationships have already been set in reconnect_tetcav.
  // We already have one-way links from new elements to exterior elements, but
  // not the other way around. 
  // We do not have these links when the faces contain ipins as these did not 
  // generate new elements, nor are they interior to the cavity. 

  // 1. Using cavity elements, a) hash the faces that contain ipins to update
  // exterior neighbours. 
  // b) look for any boundary faces and update their fac2tet. Most obvious is 
  // setting to -1. But in the case of a collapse, some volume can be carved out
  // exposing previously mesh-interior faces as now boundary faces. 
  // Hence we must check if there is a neighbour that is outside the cavity 
  // and set that as fac2tet instead of -1. These faces are invisible to new
  // cavity elements as there are none there.
  for(int ielem : cav.lctet){
    INCVDEPTH(msh.param);

    // If there is an outside element, kill self as neighbour. 
    // These faces will not be seen from new cavity elements in some cases 
    // (collapse), hence this is the only chance we have of correcting the outside
    // tet neighbour to none. 
    // Additionally, see if there is a boundary face here. Then set its fac2tet
    // to the neighbour. (see point b) above)
    for(int ifa = 0; ifa < 4; ifa++){
      int iele2 = msh.tet2tet(ielem,ifa);
      if(iele2 < 0) continue;
      if(msh.tet2tag(ithread,iele2) >= msh.tag[ithread]) continue;
      bool ifnd = false;
      for(int if2 = 0; if2 < 4; if2++){
        if(msh.tet2tet(iele2,if2) == ielem){
          msh.tet2tet(iele2,if2) = -1;
          ifnd = true;
          break;
        }
      }
      METRIS_ASSERT(ifnd);

      int iface = msh.tetfac2glo(ielem, ifa);
      if(iface < 0) continue;

      CPRINTF1(" - tetra {} interior face {} = {} {} {} becomes boundary iface = {} \n",
        ielem, ifa, 
        msh.fac2poi(iface,0),msh.fac2poi(iface,1),msh.fac2poi(iface,2),iface);
      if(msh.fac2tet(iface,0) < 0 
        || msh.tet2tag(ithread,msh.fac2tet(iface,0)) >= msh.tag[ithread]){
        // Nothing here yet
        msh.fac2tet(iface,0) = iele2;
      }else{
        // Points to dead element.
        msh.fac2tet(iface,1) = iele2;
      }

    }
    int iver = msh.template getvertet<1>(ielem, cav.ipins);
    if(iver < 0) continue;
    // Go over the other faces and see if they are cavity bdry 
    for(int ifa = 0; ifa < 4; ifa++){
      //CPRINTF1("Debug consider ielem {} ifa {}\n",ielem,ifa);
      if(ifa == iver) continue;
      int iele2 = msh.tet2tet(ielem,ifa);
      if(iele2 < 0) continue; // this is handled next
      if(iele2 >= 0 && msh.tet2tag(ithread,iele2) >= msh.tag[ithread]) continue;
      int ip1 = msh.tet2poi(ielem,lnofa3[ifa][0]);
      int ip2 = msh.tet2poi(ielem,lnofa3[ifa][1]);
      int ip3 = msh.tet2poi(ielem,lnofa3[ifa][2]);
      CPRINTF1(" - add {} to neighb hashtable {} {} {} from {} face {} \n",iele2,
               ip1,ip2,ip3,ielem,ifa);
      auto key = stup3(ip1,ip2,ip3);
      facHsh[key] = iele2;
    }
  }


  for(int ielem = nele0; ielem < msh.nelem; ielem++){
    INCVDEPTH(msh.param)
    for(int ifa = 0; ifa < 4; ifa++){
      int iele2 = msh.tet2tet(ielem,ifa);

      if(iele2 >= 0){

        // Lastly, we can be pointing to dead or cav elements. This is only if ipins 
        // was boundary: removing bdry elements exposed tetrahedra and made some 
        // inernal faces boundary. 
        if(isdeadent(iele2,msh.tet2poi) || msh.tet2tag(ithread,iele2) >= msh.tag[ithread]){
          CPRINTF1(" - found dead tetra neighbour -> remove\n");
          msh.tet2tet(ielem,ifa) = -1;
        }else{
          //if(msh.tet2tag(ithread,iele2) >= msh.tag[ithread]){
          //  METRIS_ASSERT(msh.tet2tet(iele2,0) == ielem 
          //              ||msh.tet2tet(iele2,1) == ielem 
          //              ||msh.tet2tet(iele2,2) == ielem 
          //              ||msh.tet2tet(iele2,3) == ielem);
          //  continue;
          //}

          int if2 = getfactet(msh,iele2,msh.tet2poi(ielem,lnofa3[ifa][0])
                                       ,msh.tet2poi(ielem,lnofa3[ifa][2])
                                       ,msh.tet2poi(ielem,lnofa3[ifa][1]));
          METRIS_ASSERT(if2 >= 0);

          msh.tet2tet(iele2, if2) = ielem;
          CPRINTF1(" - tetra neighbour update {} -> {}\n",iele2,ielem);
        }
      }

      int iface = msh.tetfac2glo(ielem, ifa);
      CPRINTF1(" - ielem {} ifa {} glo face? {} \n",ielem,ifa,iface);
      if(iface >= 0 && msh.fac2tet(iface,0) != ielem 
                    && msh.fac2tet(iface,1) != ielem){
        CPRINTF1(" - current fac2tet = {} {} \n",msh.fac2tet(iface,0),msh.fac2tet(iface,1));
        if(msh.fac2tet(iface,0) < 0 
        || msh.tet2tag(ithread,msh.fac2tet(iface,0)) >= msh.tag[ithread]){
          // Nothing here yet
          msh.fac2tet(iface,0) = ielem;
        }else{
          // Points to dead element.
          msh.fac2tet(iface,1) = ielem;
        }
      }

      if(msh.tet2poi(ielem, ifa) != cav.ipins){
        // If not face opposite ipins, check facHsh
        int ip1 = msh.tet2poi(ielem,lnofa3[ifa][0]);
        int ip2 = msh.tet2poi(ielem,lnofa3[ifa][1]);
        int ip3 = msh.tet2poi(ielem,lnofa3[ifa][2]);
        auto key = stup3(ip1,ip2,ip3);
        auto tt = facHsh.find(key);
        if(tt != facHsh.end()){
          int iele2 = tt->second;
          CPRINTF1(" - found outside neighbour {} in hash tab\n",iele2);
          int ifa2 = getfactet(msh,iele2,ip1,ip2,ip3);
          METRIS_ASSERT(ifa2 >= 0);
          msh.tet2tet(ielem,ifa)  = iele2;
          msh.tet2tet(iele2,ifa2) = ielem;
        }
      }
    }

    int nnode = getnnod3(ideg);
    for(int ii = 0; ii < nnode; ii++){
      int ipoin = msh.tet2poi(ielem, ii);
      METRIS_ASSERT(ipoin >= 0);
      if(msh.poi2ent(ipoin,1) <= 2 && msh.poi2ent(ipoin,1) >= 1) continue;
      if(ii < 4) msh.set_poi2ent(Vertex{ipoin}, 3, ielem);
      else       msh.set_poi2ent(CtrlPt{ipoin}, 3, ielem);
    }
  }



  if(msh.isboundary_edges()){
    for(int iedgl = 0; iedgl < ncedg; iedgl++){
      int iedge = cav.lcedg[iedgl];
      int nnode = getnnod1(ideg);
      for(int ii = 0; ii < nnode; ii++){
        int ipoin = msh.edg2poi(iedge,ii);
        if((msh.poi2ent(ipoin,0) >= 0  && ii <  2)
         ||(msh.poi2ent(ipoin,0) <= -2 && ii >= 2)) continue;
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
      int nnode = getnnod2(ideg);
      for(int ii = 0; ii < nnode; ii++){
        int ipoin = msh.fac2poi(iface,ii);
        if((msh.poi2ent(ipoin,0) >= 0  && ii <  3)
         ||(msh.poi2ent(ipoin,0) <= -2 && ii >= 3)) continue;
        int ibpoi = msh.poi2bpo[ipoin];
        while(ibpoi >= 0){
          msh.bpo2ibi(ibpoi,0) = -1;
          ibpoi = msh.bpo2ibi(ibpoi,3);
        }
        msh.poi2bpo[ipoin] = -1;
      }
    }
  }

  #if 0
  // Kill points left hanging (poi2ent still -1)
  const int tdimc = cav.get_tdim();
  {// scope
  const intAr1 &lcent = cav.lcent(tdimc);
  intAr2 &ent2poi = msh.ent2poi(tdimc);
  // The degree is not a template argument -> non constexpr 
  int nnode = msh.nnode(tdimc); 

  for(int ientt : lcent){
    INCVDEPTH(msh.param);
    for(int ii = 0; ii < nnode; ii++){
      int ipoin = ent2poi(ientt,ii);
      if((msh.poi2ent(ipoin,0) >= 0  && ii <  tdimc)
       ||(msh.poi2ent(ipoin,0) <= -2 && ii >= tdimc)) continue;
      msh.killpoint(ipoin);
    }
  }
  }// scope
  #endif

  // Deleting entities in the initial cavity. 
  for(int iedgl = 0; iedgl < ncedg; iedgl++){
    int iedge = cav.lcedg[iedgl];
    killent(iedge, msh.edg2poi);
  }
  for(int ifacl = 0; ifacl < ncfac; ifacl++){
    int iface = cav.lcfac[ifacl];
    killent(iface, msh.fac2poi);
  }
  for(int ielel = 0; ielel < nctet; ielel++){
    int ielem = cav.lctet[ielel];
    killent(ielem, msh.tet2poi);
  }


  return ierro;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template int update_cavity<MetricFieldAnalytical,n>(Mesh<MetricFieldAnalytical> &msh,\
                            MshCavity &cav, const CavWrkArrs &work, int npoi0, int nedg0, \
                            int nfac0, int nele0, \
                            int ithread);\
template int update_cavity<MetricFieldFE        ,n>(Mesh<MetricFieldFE        > &msh, \
                            MshCavity &cav, const CavWrkArrs &work, int npoi0, int nedg0, \
                            int nfac0, int nele0, \
                            int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()









} // end namespace