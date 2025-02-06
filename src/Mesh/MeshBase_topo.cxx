//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "../Mesh/MeshBase.hxx"
#include "../ho_constants.hxx"
#include "../CT_loop.hxx"
#include "../aux_utils.hxx"
#include "../aux_topo.hxx"
#include "../mprintf.hxx"


namespace Metris{

int MeshBase::fac2edg(int iface, int iedl){
  int ipoi1 = fac2poi(iface,lnoed2[iedl][0]);
  int ipoi2 = fac2poi(iface,lnoed2[iedl][1]);
  return getedgglo(*this, ipoi1, ipoi2);
}

int MeshBase::tet2fac(int ielem, int ifal){
  int ipoi1 = tet2poi(ielem,lnofa3[ifal][0]);
  int ipoi2 = tet2poi(ielem,lnofa3[ifal][1]);
  int ipoi3 = tet2poi(ielem,lnofa3[ifal][2]);
  return getfacglo(*this, ipoi1, ipoi2, ipoi3);
}

int MeshBase::poi2ebp(int ipoin, int tdim, int ientt, int iref) const {
  METRIS_ASSERT(ipoin >= 0);
  METRIS_ASSERT(tdim >= 1 && tdim <= 3);
  
  int ibpoi = poi2bpo[ipoin];
  if(ibpoi < 0) return -1;

  int pdim = bpo2ibi(ibpoi,1);

  for(int ibpo2 = ibpoi; ibpo2 >= 0; ibpo2 = bpo2ibi(ibpo2,3)){
    int itype = bpo2ibi(ibpo2,1);
    if(itype < tdim) continue;
    if(itype > tdim) return -1;

    // This is the right dimension, now 

    // If nothing specified, return this 
    if(ientt < 0 && iref < 0) return ibpo2;

    // If the point is same dim, only one entry exists and this one will do.
    if(pdim == tdim) return ibpo2;

    // Otherwise either we match the entity:
    int ient2 = bpo2ibi(ibpo2,2);
    if(ient2 == ientt) return ibpo2;

    // or the reference:
    int iref2 = tdim == 1 ? edg2ref[ient2] : fac2ref[ient2];
    if(iref == iref2) return ibpo2;

  }

  return -1;
}


int MeshBase::newpoitopo(int tdimn, int ientt){
  set_npoin(npoin+1);
  for(int ii = 0; ii < METRIS_MAXTAGS; ii++) poi2tag[ii][npoin-1] = 0;
  poi2ent[npoin-1][0] = ientt;
  poi2ent[npoin-1][1] = tdimn;
  poi2bpo[npoin-1]    = -1;
  return npoin-1;
}

void MeshBase::killpoint(int ipoin){
  for(int ibpoi = poi2bpo[ipoin]; ibpoi >= 0; ibpoi = bpo2ibi(ibpoi,3)){
    bpo2ibi(ibpoi,0) = -1;
  }
  poi2bpo[ipoin] = -1;
  poi2ent(ipoin,0) = -1;
  poi2ent(ipoin,1) = -1;
}

// Create new face by copying from tetrahedron
template <int ideg>
void MeshBase::newfactopo(int ielem, int ifael, int iref, int iele2){

  int ifacn = nface;
  set_nface(nface+1);
  //if(nface >= mface)METRIS_THROW_MSG(DMemExcept(),
  //  "INCREASE MFACE (iniMeshNeighbours)");
  METRIS_ASSERT(ifael >= 0 && ifael < 4);
  METRIS_ASSERT(!isdeadent(ielem, tet2poi));

  int ip1 = tet2poi(ielem,lnofa3[ifael][0]);
  int ip2 = tet2poi(ielem,lnofa3[ifael][1]);
  int ip3 = tet2poi(ielem,lnofa3[ifael][2]);

  int ib1 = ip1;
  int ib2 = ip2; 
  int ib3 = ip3; 

  // Don't overwrite edge/corner links
  if(ib1 >= 0 && bpo2ibi(ib1,1) > 1){
    poi2ent(ip1,0) = ifacn;
    poi2ent(ip1,1) = 2;
  }
  if(ib2 >= 0 && bpo2ibi(ib2,1) > 1){
    poi2ent(ip2,0) = ifacn;
    poi2ent(ip2,1) = 2;
  }
  if(ib3 >= 0 && bpo2ibi(ib3,1) > 1){
    poi2ent(ip3,0) = ifacn;
    poi2ent(ip3,1) = 2;
  }

  fac2poi(ifacn,0) = ip1;
  fac2poi(ifacn,1) = ip2;
  fac2poi(ifacn,2) = ip3;


  if constexpr(ideg > 1){
    int nppf = facnpps[ideg] - 3;
    int nppe = edgnpps[ideg] - 2;
    int idx0 = 4 + 6*nppe + ifael*nppf;
    for(int i=0; i < nppf; i++){
      int ipoin = tet2poi[ielem][idx0 + i];
      fac2poi[ifacn][3+i] = ipoin;
      // Tet face interior points can be neither corners nor edge points:
      poi2ent(ipoin,0) = ifacn;
      poi2ent(ipoin,1) = 2;
    }
  }
  fac2ref[ifacn]    = iref;
  for(int ii = 0; ii < METRIS_MAXTAGS; ii++) fac2tag(ii,ifacn) = 0; 

  fac2fac(ifacn,0) = -1;
  fac2fac(ifacn,1) = -1;
  fac2fac(ifacn,2) = -1;

  fac2tet(ifacn,0) = ielem;
  fac2tet(ifacn,1) = iele2;

  auto key = stup3(fac2poi(ifacn,0),fac2poi(ifacn,1),fac2poi(ifacn,2));
  facHshTab[key] = ifacn;

}


int MeshBase::newfacvirtual(int iref){
  int iface = nface;
  set_nface(nface+1);
  killent(iface,fac2poi);
  fac2ref[iface] = iref;
  for(int ii = 0; ii < METRIS_MAXTAGS; ii++) fac2tag(ii,iface) = tag[ii];
  return iface;
}

int MeshBase::newedgvirtual(int iref){
  int iedge = nedge;
  set_nedge(nedge+1);
  killent(iedge,edg2poi);
  edg2ref[iedge] = iref;
  for(int ii = 0; ii < METRIS_MAXTAGS; ii++) edg2tag(ii,iedge) = tag[ii];
  return iedge;
}



template <int ideg>
void MeshBase::newedgtopo(int iface, int iedfa, int iref){
  int iedgn = nedge;
  set_nedge(nedge+1);
  assert(iedfa >= 0 && iedfa < 3);
  assert(!isdeadent(iface, fac2poi));

  int ip1 = fac2poi(iface,lnoed2[iedfa][0]);
  int ip2 = fac2poi(iface,lnoed2[iedfa][1]); 
  edg2poi(iedgn,0) = ip1;
  edg2poi(iedgn,1) = ip2;

  int ib1 = poi2bpo[ip1];
  int ib2 = poi2bpo[ip2]; 

  // End vertices can be corners. Make sure not to overwrite poi2ent. 
  if(ib1 > 0 && bpo2ibi(ib1,1) > 0){
    poi2ent(ip1,0) = iedgn;
    poi2ent(ip1,1) = 1;
  } 
  if(ib2 > 0 && bpo2ibi(ib2,1) > 0){
    poi2ent(ip2,0) = iedgn;
    poi2ent(ip2,1) = 1;
  } 


  if constexpr(ideg > 1){
    int nppe = edgnpps[ideg] - 2;
    int idx0 = 3 + iedfa*nppe;
    for(int i = 0; i < nppe; i++){
      int ipoin = fac2poi[iface][idx0 + i];
      edg2poi[iedgn][2+i] = ipoin;
      // Edge interior vertices cannot be corners. 
      poi2ent(ipoin,0) = iedgn;
      poi2ent(ipoin,1) = iedgn;
    }
  }

  edg2ref[iedgn]    = iref;
  for(int ii = 0; ii < METRIS_MAXTAGS; ii++) edg2tag(ii,iedgn) = 0; 

  edg2edg(iedgn,0) = -1;
  edg2edg(iedgn,1) = -1;

  edg2fac[iedgn] = iface;

  auto key = stup2(edg2poi(iedgn,0),edg2poi(iedgn,1));
  edgHshTab[key] = iedgn;

}

#define BOOST_PP_LOCAL_MACRO(n)\
template void MeshBase::newfactopo< n >(int ielem, int ifael, int iref, int iele2);\
template void MeshBase::newedgtopo< n >(int iface, int iedfa, int iref);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


// Note: this is a stack more than a linked list, really.  
template<int tdim>
int MeshBase::newbpotopo(int ipoin, int ientt){
  static_assert(tdim >= 0 && tdim < 3);
  
  if(tdim == 2 && !isboundary_faces()) return -1;
  if(tdim == 1 && !isboundary_edges()) return -1;

  METRIS_ASSERT(ipoin < npoin);

  int ibpon = nbpoi;

  // Get any ibpoi already attached 
  int ibpoi = poi2bpo[ipoin]; 
  if(ibpoi < 0) ibpoi = -1;

  // Only update poi2bpo if new ibpon is the lowest-dimensional, or none yet
  if(ibpoi < 0 || bpo2ibi(ibpoi,1) > tdim){

    // As this is lowest-dim, put at start. 
    poi2bpo[ipoin] = ibpon;
    if(ientt >= 0){
      poi2ent(ipoin,0) = ientt;
      poi2ent(ipoin,1) = tdim;
    }

    // Create new ibpoi
    set_nbpoi(nbpoi+1);
    bpo2ibi(ibpon,3) = ibpoi; // Link to next
    
  }else{
    // If not then after all tdim <= current passed and check ientt not already
    // present 
    int ibprv = ibpoi;
    for(;ibpoi >= 0 && bpo2ibi(ibpoi,1) <= tdim; ibpoi = bpo2ibi(ibpoi,3)){
      ibprv = ibpoi;
      int itype = bpo2ibi(ibpoi,1);
      if(itype < tdim) continue;
      // This case same dim 
      int ient0 = bpo2ibi(ibpoi,2); 
      // If same entity, return this ibpoi and do not create a new one. 
      // Negative entry is specific to initialization to store a ref. These 
      // can be duplicated. 
      if(ient0 == ientt && ient0 >= 0) return ibpoi;
    }

    //int ibpo2 = ibpoi;
    //int ibprv = -1;
    //do{
    //  ibprv = ibpo2;
    //  ibpo2 = bpo2ibi(ibpo2,3);
    //  if(ibpo2 < 0 || bpo2ibi(ibpo2,1) > tdim) break;
    //}while(ibpo2 >= 0 && ibpo2 != ibpoi);

    // Create new ibpoi
    set_nbpoi(nbpoi+1);
    int tmp = bpo2ibi(ibprv,3);
    bpo2ibi(ibprv,3) = ibpon;
    bpo2ibi(ibpon,3) = tmp;

  }
  for(int ii = 0; ii < METRIS_MAXTAGS; ii++) bpo2tag(ii,ibpon) = 0;
  bpo2ibi(ibpon,0) = ipoin; 
  bpo2ibi(ibpon,1) = tdim;  // Type
  bpo2ibi(ibpon,2) = ientt; // Ref
  //bpo2ibi(ibpon,3) = ibpoi; // Link to next

  for(int i = 0; i < nrbi ;i++) bpo2rbi(ibpon,i) = 0;
  

  //if(ibpoi >= 0){
  //  // There is no particular order here, just put it at the start
  //  int tmp = bpo2ibi(ibpoi,3);
  //  bpo2ibi(ibpoi,3) = ibpon;
  //  bpo2ibi(ibpon,3) = tmp;
  //  //if(tmp < 0) bpo2ibi(ibpon,3) = ibpoi; // Not a loop
  //}

  return ibpon;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template int MeshBase::newbpotopo< n >(int ipoin, int ientt);
#define BOOST_PP_LOCAL_LIMITS     (0,2)
#include BOOST_PP_LOCAL_ITERATE()


void MeshBase::rembpotag(int ipoin, int ithread){
  METRIS_ASSERT(ithread >= 0 && ithread < METRIS_MAXTAGS)
  METRIS_ASSERT(ipoin >= 0 && ipoin < npoin);
  int ibpoi = poi2bpo[ipoin];
  METRIS_ASSERT(ibpoi >= 0 && ibpoi < nbpoi);

  int itag  = tag[ithread];
  int ibpoc = ibpoi; // current 
  int ibpop = -1;    // previous
  int nn = 0;
  do{
    nn++;
    METRIS_ASSERT(nn <= METRIS_MAX_WHILE);
    METRIS_ASSERT(ibpoc >= 0 && ibpoc < nbpoi);
    int tdim  = bpo2ibi(ibpoc,1);
    int ientt = bpo2ibi(ibpoc,2);

    bool rement = false;

    if(ientt < 0) rement = true;
    else{
      int ietag = tdim == 1 ? edg2tag(ithread,ientt) :
                  tdim == 2 ? fac2tag(ithread,ientt) : itag-1;
      if(ietag >= itag) rement = true;
    }
    int ibpon = bpo2ibi(ibpoc,3); // next
    
    if(rement){ // Remove this entry 
      // Update previous link. If ibpop == -1, this is still poi2bpo. 
      if(ibpop == -1){
        poi2bpo[ipoin] = ibpon;
      }else{
        bpo2ibi(ibpop,3) = ibpon;
      }
      // Invalidate current ibpoi. Note we could tuck index behind our ear for later
      for(int ii = 0; ii < nibi; ii++) bpo2ibi(ibpoc,ii) = -1;
    }else{
      ibpop = ibpoc;
    }
    ibpoc = ibpon;
  }while(ibpoc >= 0 && ibpoc != ibpoi);

}



int MeshBase::getpoitdim(int ipoin) const{
  int ibpoi = poi2bpo[ipoin]; 
  if(ibpoi < 0){
    if(idim >= 3) return 3; 
    else if(!isboundary_faces()) return 2;
    else if(!isboundary_edges()) return 1;
    else return 0;
  }

  return bpo2ibi(ibpoi,1);
}


int MeshBase::facedg2glo(int iface, int iedl) const{
  METRIS_ASSERT(iface >= 0 && iface < nface);
  METRIS_ASSERT(iedl >= 0 && iedl < 3);
  int ip1 = fac2poi(iface,lnoed2[iedl][0]);
  int ip2 = fac2poi(iface,lnoed2[iedl][1]);

  return getedgglo(*this,ip1,ip2);
}


int MeshBase::tetfac2glo(int ielem, int ifal) const{
  METRIS_ASSERT(ielem >= 0 && ielem < nelem);
  METRIS_ASSERT(ifal >= 0 && ifal < 4);
  int ip1 = tet2poi(ielem,lnofa3[ifal][0]);
  int ip2 = tet2poi(ielem,lnofa3[ifal][1]);
  int ip3 = tet2poi(ielem,lnofa3[ifal][2]);

  return getfacglo(*this,ip1,ip2,ip3);
}

int MeshBase::tetedg2glo(int ielem, int iedl) const{
  METRIS_ASSERT(ielem >= 0 && ielem < nelem);
  METRIS_ASSERT(iedl >= 0 && iedl < 6);
  int ip1 = tet2poi(ielem,lnoed3[iedl][0]);
  int ip2 = tet2poi(ielem,lnoed3[iedl][1]);

  return getedgglo(*this,ip1,ip2);
}



int MeshBase::getverent(int ientt, int tdimn, int ipoin){
  int iver;
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == curdeg){
    if(tdimn == 1){
      iver = getveredg<ideg>(ientt,edg2poi,ipoin);
    }else if(tdimn == 2){
      iver = getverfac<ideg>(ientt,fac2poi,ipoin);
    }else{
      iver = getvertet<ideg>(ientt,tet2poi,ipoin);
    }
  }}CT_FOR1(ideg);
  return iver;
}

template <int ideg>
int MeshBase::getverent(int ientt, int tdimn, int ipoin){
  if(tdimn == 1){
    return getveredg<ideg>(ientt,edg2poi,ipoin);
  }else if(tdimn == 2){
    return getverfac<ideg>(ientt,fac2poi,ipoin);
  }else{
    return getvertet<ideg>(ientt,tet2poi,ipoin);
  }
  return -2;
}


#define BOOST_PP_LOCAL_MACRO(n)\
template int MeshBase::getverent<n>(int ientt, int tdimn, int ipoin);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()





}// end namespace