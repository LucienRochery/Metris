//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../msh_checktopo.hxx"
#include "../aux_utils.hxx"
#include "../mprintf.hxx"


namespace Metris{

template<class MetricFieldType>
void Mesh<MetricFieldType>::cleanup(int ithread){

  GETVDEPTH((*this));

  intAr1& lentt = this->iwork;

  intAr1 lpoin(this->npoin);
  lpoin.set_n(this->npoin);



 
  // ----- Geometric links 
  int nbpon = 0;
  lentt.set_n(this->nbpoi);
  for(int ibpoi = 0; ibpoi < this->nbpoi; ibpoi++){
    int ipoin = this->bpo2ibi(ibpoi,0); 
    if(ipoin < 0) continue;

    int ibpon = nbpon; 
    lentt[ibpoi] = ibpon; 
    nbpon++; 

    for(int ii = 0; ii < nibi; ii++) 
      this->bpo2ibi(ibpon,ii) = this->bpo2ibi(ibpoi,ii);

    for(int ii = 0; ii < nrbi; ii++) 
      this->bpo2rbi(ibpon,ii) = this->bpo2rbi(ibpoi,ii);
    
    if(this->poi2bpo[ipoin] == ibpoi) this->poi2bpo[ipoin] = ibpon; 
  }

  if(nbpon != this->nbpoi){

    CPRINTF2(" - cleanup bdry links %d -> %d \n",this->nbpoi,nbpon);

    // Update link to next 
    for(int ibpoi = 0; ibpoi < nbpon; ibpoi++){
      int ibnxt = this->bpo2ibi(ibpoi,3); 
      if(ibnxt < 0) continue;
      this->bpo2ibi(ibpoi,3) = lentt[ibnxt]; 
    }
  }

  // ----- Vertices 
  // Concerned arrays:
  // coord
  // met 
  // poi2bpo 
  // poi2bak 
  // poi2ent -> recompute 
  const int idim = this->idim;
  const int nnmet = (idim*(idim+1))/2;
  int nponn = 0;
  for(int ipoin = 0; ipoin < this->npoin; ipoin++){
    if(this->poi2ent(ipoin,0) < 0) continue;

    int iponn = nponn;
    nponn++;

    lpoin[ipoin] = iponn; 

    if(iponn == ipoin) continue;

    for(int ii = 0; ii < idim; ii++) 
      this->coord(iponn,ii) = this->coord(ipoin,ii);
    for(int ii = 0; ii < nnmet; ii++) 
      this->met(iponn,ii) = this->met(ipoin,ii);
    
    this->poi2bpo[iponn] = this->poi2bpo[ipoin]; 
    for(int ii = 1; ii <= this->get_tdim(); ii++)
      this->poi2bak(iponn,ii-1) = this->poi2bak(ipoin,ii-1); 
  }

  if(this->npoin == nponn) goto update_tetras; 

  CPRINTF2(" - cleanup vertices %d -> %d \n",this->npoin,nponn);

  for(int ibpoi = 0; ibpoi < nbpon; ibpoi++){
    int ipoin = this->bpo2ibi(ibpoi,0); 
    if(ipoin < 0) continue;
    int ityp = this->bpo2ibi(ibpoi,1); 

    this->bpo2ibi(ibpoi,0) = lpoin[ipoin]; 
    if(ityp == 0) this->bpo2ibi(ibpoi,2) = lpoin[ipoin]; 
  }



  // ----- Tetras
  update_tetras:
  int nelen = 0; 
  lentt.set_n(this->nelem);
  for(int ielem = 0; ielem < this->nelem; ielem++){
    if(isdeadent(ielem,this->tet2poi)) continue;

    int ielen    = nelen; 
    lentt[ielem] = ielen; 
    nelen++; 

    int nnode = tetnpps[this->curdeg]; 
    // Not just copy but also translation using lpoin 
    for(int ii = 0; ii < nnode; ii++) 
      this->tet2poi(ielen,ii) = lpoin[this->tet2poi(ielem,ii)]; 

    // The next updates are only needed if the index changed. 
    if(ielen == ielem) continue; 

    this->tet2ref[ielen] = this->tet2ref[ielem]; 

    for(int ii = 0; ii < 4; ii++) 
      this->tet2tet(ielen,ii) = this->tet2tet(ielem,ii);
  }
  
  if(nelen == this->nelem) goto update_faces;

  CPRINTF2(" - cleanup tetras %d -> %d \n",this->nelem,nelen);

  for(int iface = 0; iface < this->nface; iface++){
    for(int ii = 0; ii < 2; ii++){
      int ielem = this->fac2tet(iface,ii);
      if(ielem < 0) continue;
      this->fac2tet(iface,ii) = lentt[ielem]; 
    }
  }


  for(int ielem = 0; ielem < nelen; ielem++){
    for(int ii = 0; ii < 4; ii++){
      int ineig = this->tet2tet(ielem,ii); 
      if(ineig < 0) continue;
      this->tet2tet(ielem,ii) = lentt[ineig]; 
    }
  }


  // ----- Faces
  update_faces:
  // Concerned arrays: 
  // fac2poi 
  // fac2ref 
  // fac2fac 
  // fac2tet 
  // poi2bpo when ityp == 2
  int nfacn = 0; 
  lentt.set_n(this->nface);
  this->facHshTab.erase(this->facHshTab.begin(), this->facHshTab.end());
  for(int iface = 0; iface < this->nface; iface++){
    if(isdeadent(iface,this->fac2poi)) continue;

    int ifacn    = nfacn; 
    lentt[iface] = ifacn; 
    nfacn++; 

    int nnode = facnpps[this->curdeg]; 
    // Not just copy but also translation using lpoin 
    for(int ii = 0; ii < nnode; ii++) 
      this->fac2poi(ifacn,ii) = lpoin[this->fac2poi(iface,ii)]; 

    //// Necessary to update both iface -> ifacn but also lpoin 
    //auto key_old = stup3(this->fac2poi(iface,0), 
    //                     this->fac2poi(iface,1),
    //                     this->fac2poi(iface,2));
    //this->facHshTab.erase(key_old); 

    auto key_new = stup3(this->fac2poi(ifacn,0), 
                         this->fac2poi(ifacn,1),
                         this->fac2poi(ifacn,2));
    this->facHshTab.insert({key_new,ifacn}); 

    // The next updates are only needed if the index changed. 
    if(ifacn == iface) continue; 

    this->fac2ref[ifacn] = this->fac2ref[iface]; 
    if(idim >= 3){
      this->fac2tet(ifacn,0) = this->fac2tet(iface,0); 
      this->fac2tet(ifacn,1) = this->fac2tet(iface,1); 
    }

    for(int ii = 0; ii < 3; ii++) 
      this->fac2fac(ifacn,ii) = this->fac2fac(iface,ii);

  }

  if(nfacn == this->nface) goto update_edges;

  CPRINTF2(" - cleanup faces %d -> %d \n",this->nface,nfacn);

  for(int ibpoi = 0; ibpoi < nbpon; ibpoi++){
    int ipoin = this->bpo2ibi(ibpoi,0); 
    if(ipoin < 0) continue;
    int ityp = this->bpo2ibi(ibpoi,1); 

    if(ityp == 2){
      this->bpo2ibi(ibpoi,2) = lentt[this->bpo2ibi(ibpoi,2)]; 
    }
  }

  for(int iface = 0; iface < nfacn; iface++){
    for(int ii = 0; ii < 3; ii++){
      int ineig = this->fac2fac(iface,ii); 
      if(ineig == -1) continue;
      if(ineig < 0){
        ineig = - ineig - 2;
        ineig = lentt[ineig]; 
        this->fac2fac(iface,ii) = - ineig - 2;
      }else{
        this->fac2fac(iface,ii) = lentt[ineig]; 
      }
    }
  }

  for(int iedge = 0; iedge < this->nedge; iedge++){
    int iface = this->edg2fac[iedge];
    if(iface < 0) continue;
    this->edg2fac[iedge] = lentt[iface]; 
  }


  // ----- Edges
  update_edges:
  // Concerned arrays: 
  // edg2poi 
  // edg2ref 
  // edg2fac 
  // poi2bpo when ityp == 1
  // edg2edg 
  int nedgn = 0; 
  lentt.set_n(this->nedge);
  // Outright erase edgHshTab. 
  // This avoids shrinking the buffers. 
  this->edgHshTab.erase(this->edgHshTab.begin(), this->edgHshTab.end());
  for(int iedge = 0; iedge < this->nedge; iedge++){
    if(isdeadent(iedge,this->edg2poi)) continue;

    int iedgn    = nedgn; 
    lentt[iedge] = iedgn; 
    nedgn++; 

    int nnode = edgnpps[this->curdeg]; 
    // Not just copy but also translation using lpoin 
    for(int ii = 0; ii < nnode; ii++){
      METRIS_ASSERT(lpoin[this->edg2poi(iedge,ii)] >= 0 
         && lpoin[this->edg2poi(iedge,ii)] < nponn);
      this->edg2poi(iedgn,ii) = lpoin[this->edg2poi(iedge,ii)]; 
    }

    //// Necessary to update both iedge -> iedgn but also lpoin 
    //auto key_old = stup2(this->edg2poi(iedge,0), this->edg2poi(iedge,1));
    //// ATTENTION! key_old may, per chance, be the nodes of a NEW edge !
    //// Only erase if this is pointing to iedge. 
    //if(this->edgHshTab[key_old] == iedge) this->edgHshTab.erase(key_old); 
    

    auto key_new = stup2(this->edg2poi(iedgn,0), this->edg2poi(iedgn,1));
    this->edgHshTab.insert({key_new,iedgn}); 

    // The next updates are only needed if the index changed. 
    if(iedgn == iedge) continue; 

    this->edg2ref[iedgn] = this->edg2ref[iedge]; 
    this->edg2fac[iedgn] = this->edg2fac[iedge]; 


    for(int ii = 0; ii < 2; ii++) 
      this->edg2edg(iedgn,ii) = this->edg2edg(iedge,ii);

  }


  if(nedgn == this->nedge) goto update_final; 


  CPRINTF2(" - cleanup edges %d -> %d \n",this->nedge,nedgn);

  for(int ibpoi = 0; ibpoi < nbpon; ibpoi++){
    int ipoin = this->bpo2ibi(ibpoi,0); 
    if(ipoin < 0) continue;
    METRIS_ASSERT(ipoin < this->npoin);
    int ityp = this->bpo2ibi(ibpoi,1); 

    // Now edge update proper. 
    if(ityp == 1){
      this->bpo2ibi(ibpoi,2) = lentt[this->bpo2ibi(ibpoi,2)]; 
    }
  }

  for(int iedge = 0; iedge < nedgn; iedge++){
    for(int ii = 0; ii < 2; ii++){
      int ineig = this->edg2edg(iedge,ii); 
      if(ineig == -1) continue;
      if(ineig < 0){
        ineig = - ineig - 2;
        ineig = lentt[ineig]; 
        this->edg2edg(iedge,ii) = - ineig - 2;
      }else{
        this->edg2edg(iedge,ii) = lentt[ineig]; 
      }
    }
  }

  // ----- Final updates
  update_final:
  this->set_nedge(nedgn);
  this->set_nface(nfacn);
  this->set_nelem(nelen);
  this->set_npoin(nponn);
  this->set_nbpoi(nbpon);

  // Recompute poi2ent 
  int tdimm = this->get_tdim();
  for(int tdimn = tdimm; tdimn >= 1; tdimn--){
    intAr2 &ent2poi = this->ent2poi(tdimn);
    int nnode = this->nnode(tdimn);
    int nentt = this->nentt(tdimn); 
    for(int ientt = 0; ientt < nentt; ientt++){
      for(int ii = 0; ii < nnode; ii++){
        int ipoin = ent2poi(ientt,ii);
        this->poi2ent(ipoin,0) = ientt;
        this->poi2ent(ipoin,1) = tdimn;
      }
    }
  }

  if(this->param->dbgfull) check_topo(*this);


}
template void Mesh<MetricFieldFE>::cleanup(int ithread);
template void Mesh<MetricFieldAnalytical>::cleanup(int ithread);


} // End namespace
