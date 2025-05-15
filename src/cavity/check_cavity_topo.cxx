//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
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

  bool icheck1 = cav.lcfac.get_n() > 0 && msh.nelem > 0;
  bool icheck2 = !opts.allow_remove_points;

  if(icheck1 || icheck2){
    for(int ielem : cav.lctet) msh.tet2tag(ithread,ielem) = msh.tag[ithread];
  }

  // Check that supported tetrahedra are collapsed together with faces
  if(cav.lcfac.get_n() > 0 && msh.nelem > 0){
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

  // Check no rem pts
  if(!opts.allow_remove_points){
    
    for(int iedge : cav.lcedg) msh.edg2tag(ithread,iedge) = msh.tag[ithread];
    for(int iface : cav.lcfac) msh.fac2tag(ithread,iface) = msh.tag[ithread];

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
