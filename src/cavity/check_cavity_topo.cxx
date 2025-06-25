//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../cavity/msh_cavity.hxx"

#include "../aux_topo.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../Mesh/Mesh.hxx"
#include "../utils/mprintf.hxx"


namespace Metris{

// Note: tag is incremented outside. 
int check_cavity_topo(MeshBase &msh, MshCavity &cav, 
                      CavOprOpt &opts, //RoutineWorkMemory<int> &iwrk, 
                      int ithread){

  cav.maxtag = ++msh.tag[ithread];

  GETVDEPTH(msh.param);

  //// Check tetrahedra supported by faces are removed if those faces are removed
  //bool icheck1 = cav.lcfac.get_n() > 0 && msh.nelem > 0;
  //// Check triangles supported by edges are removed if those edges are removed
  //bool icheck2 = cav.lcedg.get_n() > 0 && msh.nface > 0;
  //// Check point removal rules are respected.
  //bool icheck3 = !opts.allow_remove_points;


  // Tag while checking no duplicates (hard error).
  for(int tdim = 1; tdim <= msh.get_tdim(); tdim++){
    const intAr1& lcent = cav.lcent(tdim);
    intAr2& ent2tag = msh.ent2tag(tdim);
    const int nnode = msh.nnode(tdim);
    const intAr2 &ent2poi = msh.ent2poi(tdim);
    for(int ientt : lcent){
      METRIS_ASSERT(ent2tag(ithread,ientt) < msh.tag[ithread]);
      ent2tag(ithread,ientt) = msh.tag[ithread];
      if(!cav.inewp) continue;
      // Check if ipins is present in the cavity elements.
      for(int inode = 0; inode < nnode; inode++){
        int ipoin = ent2poi(ientt, inode);
        if(ipoin == cav.ipins) cav.inewp = false;
      }
    }
  }

  if(opts.skip_topo_checks){
    CPRINTF1(" - flag skip_topo_checks set: skipping topological checks.\n");
    return 0;
  }


  // Check that supported triangles are collapsed together with edges
  if(msh.nface > 0){
    for(int iedge : cav.lcedg){
      int iface = msh.edg2fac[iedge];
      if(iface < 0) continue;
      if(msh.fac2tag(ithread,iface) < msh.tag[ithread]){
        CPRINTF1(" # cavity edge %d supports face %d not in cavity\n",iedge,iface);
        return CAV_ERR_NOEDGFAC;
      }
      int ip1 = msh.edg2poi(iedge,0);
      int ip2 = msh.edg2poi(iedge,1);
      int ied = getedgfac(msh, iface, ip1, ip2);
      METRIS_ASSERT(ied >= 0);
      int ifac2 = msh.fac2fac(iface, ied);
      if(ifac2 == -1) continue;
      if(ifac2 >= 0){
        if(msh.fac2tag(ithread,ifac2) < msh.tag[ithread]){
          CPRINTF1(" # cavity edge %d supports n-m face %d not in cavity\n",iedge,ifac2);
          return CAV_ERR_NOEDGFAC;
        }
        continue;
      }
      ifac2 = iface;
      while(getnextfacnm(msh, iface, ip1, ip2, &ifac2, &ied)){
        CPRINTF1(" - check face %d nm face nei %d\n",iface,ifac2);
        METRIS_ASSERT(ifac2 >= 0 && ifac2 < msh.nface);
        METRIS_ASSERT(!isdeadent(ifac2, msh.fac2poi));
        METRIS_ASSERT(ied >= 0);
        if(msh.fac2tag(ithread,ifac2) < msh.tag[ithread]){
          CPRINTF1(" # cavity edge %d supports n-m face %d not in cavity\n",iedge,ifac2);
          return CAV_ERR_NOEDGFAC;
        }
      }
    }
  }
  // Check that supported tetrahedra are collapsed together with faces
  if(msh.nelem > 0){
    for(int iface : cav.lcfac){
      for(int ii = 0; ii < 2; ii++){
        int ielem = msh.fac2tet(iface,ii);
        if(ielem < 0) continue;
        if(msh.tet2tag(ithread,ielem) < msh.tag[ithread]){
          CPRINTF1(" # cavity face %d supports tet %d not in cavity\n",iface,ielem);
          return CAV_ERR_NOFACTET;
        }
      }
    }
  }

  // Ensure internal edges had been added to the cavity 
  for(int iface : cav.lcfac){
    for(int ied = 0; ied < 3 ;ied++){
      int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
      int ip2 = msh.fac2poi(iface,lnoed2[ied][1]);
      int ifac2 = msh.fac2fac(iface,ied);
      if(ifac2 == -1) continue;

      if(ifac2 >= 0){
        if(msh.fac2tag(ithread,ifac2) < msh.tag[ithread]) continue;
      }else if(ifac2 < -1){// Non manifold
        // Check any of the other nm faces are in the cavity. 
        // If yes, the edge needs to be as well. 


        int ifac3 = iface;
        bool iint = false;
        int ied2 = ied;
        while(getnextfacnm(msh,iface,ip1,ip2,&ifac3,&ied2)){
          if(msh.fac2tag(ithread,ifac3) >= msh.tag[ithread]){
            iint = true;
            break;
          }
        }

        if(!iint) continue;
      }

      int iedge = getedgglo(msh,ip1,ip2);

      if(iedge < 0 && ifac2 < 0) METRIS_THROW_MSG(TopoExcept(),"Non manifold and no edge");

      if(iedge < 0) continue;

      if(msh.edg2tag(ithread,iedge) < msh.tag[ithread]){
        CPRINTF1("## edge %d is internal but was not in cavity\n",iedge);
        // This is not always an error in the sense of an assert.
        // The assert has proved useful to spot legitimate bugs but let's downgrade it now
        return CAV_ERR_INTEDG;
      }
    }
  }


  for(int itetr : cav.lctet){
    for(int ifa = 0; ifa < 4; ifa++){
      int itet2 = msh.tet2tet(itetr, ifa);
      if(itet2 < 0) continue;
      if(msh.tet2tag(ithread,itet2) < msh.tag[ithread]) continue;
      // Don't use refs here so we can support sheet bodies next (surface internal to single tet domain)
      int ip1 = msh.tet2poi(itetr, lnofa3[ifa][0]);
      int ip2 = msh.tet2poi(itetr, lnofa3[ifa][1]);
      int ip3 = msh.tet2poi(itetr, lnofa3[ifa][2]);
      int iface = getfacglo(msh, ip1, ip2, ip3);
      if(iface < 0) continue;
      if(msh.fac2tag(ithread,iface) < msh.tag[ithread]){
        CPRINTF1("## face %d is internal but was not in cavity\n",iface);
        return CAV_ERR_INTFAC;
      }
    }
  }



  // Check no rem pts
  if(!opts.allow_remove_points){
    
    for(int iedge : cav.lcedg) msh.edg2tag(ithread,iedge) = msh.tag[ithread];

    // Points to be removed are those that are surrounded by only cavity elements.
    // Hence, loop over cavity elements and tag any points that belong to a 
    // non-cavity neighbour. 
    // Lastly, count untagged vertices. 

    // If it belongs to any lower dim elements, that should be in the cavity. 
    // It suffice there is one, as if it doesnt belong to all, it would be tagged. 

    int tdimn = cav.lctet.get_n() > 0 ? 3 
              : cav.lcfac.get_n() > 0 ? 2 
                                      : 1;
    const intAr1&  lcent = cav.lcent(tdimn);
    const intAr2&  ent2ent = msh.ent2ent(tdimn);
    const intAr2&  ent2poi = msh.ent2poi(tdimn);
    const intAr2r& ent2tag = msh.ent2tag(tdimn);

    // ipins should always be seeded with a newbpotopo if it is going to be bdry
    const int pdim_ipins = msh.getpoitdim(cav.ipins);
    METRIS_ASSERT_MSG(pdim_ipins >= 0 && pdim_ipins <= msh.get_tdim(),
      "pdim_ipins = "<<pdim_ipins);

    // Tag points that won't be deleted: there is at least one elt outside
    // the cavity that has the point. 
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
            CPRINTF2("  - not rem point %d \n", ipoin);
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
        CPRINTF2("  - rem pt ? %d \n", ipoin);

        // Check the point dimension wrt to option allow_remove_points_superdim
        int pdim = msh.getpoitdim(ipoin);
        if(pdim > pdim_ipins && opts.allow_remove_points_superdim){
          CPRINTF1(" - point dim %d > %d = dim(ipins) " 
                   "with allow_remove_points_superdim, skip check\n",
                   pdim, pdim_ipins);
          continue;
        }

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

        CPRINTF1(" ## norempts and point %d will be removed \n",ipoin);
        
        return CAV_ERR_REMPT;

      }
    }

  } // if(!opts.allow_remove_points)


  return 0;
}


//template int 
//check_cavity_topo<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, 
// MshCavity &cav, CavOprOpt &opts,int ithread);
//template int 
//check_cavity_topo<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, 
// MshCavity &cav, CavOprOpt &opts,int ithread);


} // End namespace
