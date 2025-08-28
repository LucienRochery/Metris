//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../cavity/msh_cavity.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../aux_topo.hxx"
#include "../ho_constants.hxx"
#include "../low_geo/normal.hxx"
#include "../low_geo/misc.hxx"
#include "../quality/low_metqua.hxx"
#include "../io_libmeshb.hxx"

#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"
#include "../utils/aux_misc.hxx"

namespace Metris{


int update_bpois_newp(MeshBase &msh, const MshCavity &cav, CavWrkArrs &work,
                      int npoi0, int nfac0, int ithread){

  GETVDEPTH(msh.param);

  CPRINTF1("-- START update_bpois_newp {} <= ipoin < {} + ipins = {}\n",
           npoi0, msh.npoin, cav.ipins);

  if(!msh.isboundary_faces()){
    CPRINTF1("-- END update_bpois_newp, mesh dim {} <= 2\n",msh.idim);
    return 0;
  }

  if(!msh.CAD()){
    CPRINTF1("-- END update_bpois_newp, no CAD link\n");
    return 0;
  }

  bool ipins_uv = msh.getpoitdim(cav.ipins) < 2 && cav.inewp;
  if(ipins_uv){
    // In this case, we only have ipins as a dim 1 point. We'll use EG_getEdgeUV
    METRIS_ASSERT(msh.getpoitdim(cav.ipins) == 1);

    // Face connex component bounding edge information
    const auto& edcco = work.edcco;

    // To each connex component, a sign for each edge that bounds it is associated. 
    for(int iface = nfac0; iface < msh.nface; iface++){
      INCVDEPTH(msh.param);

      int icoco = -(msh.fac2tag(ithread,iface) - msh.tag[ithread]) - 1;
      METRIS_ASSERT_MSG(icoco >= 0,
        "face {} icoco {}", iface, icoco);
      METRIS_ASSERT(icoco < work.lfcco.get_n());
      METRIS_ASSERT_MSG(edcco[icoco].get_n() > 0,"edcco[icoco] empty in update_bpois_newp");
      int ireff = msh.fac2ref[iface];
      ego EG_fac = msh.CAD.cad2fac[ireff];

      CPRINTF1(" - iface {} iref {} update ipins bpoi\n",iface,ireff);

      // Find an edge ibpoi for an edge that bounds this connex component
      bool found_update = false;
      for(int ibpoe = msh.poi2bpo[cav.ipins]; ibpoe >= 0; ibpoe = msh.bpo2ibi(ibpoe,3)){
        int tdim = msh.bpo2ibi(ibpoe,1);
        if(tdim != 1) continue;
        int iedge = msh.bpo2ibi(ibpoe,2);
        int irefe = msh.edg2ref[iedge];
        METRIS_ASSERT(irefe >= 0);
        CPRINTF1("   - iedge {} iref {} bpo link potential coco bound t {}\n",
                 iedge,irefe,msh.bpo2rbi(ibpoe,0));
        for(auto ipair : edcco[icoco]){
          CPRINTF1("     - pair {} {} \n",ipair.first,ipair.second);
          if(irefe != ipair.first) continue;
          found_update = true;
          int isens = ipair.second;

          // Use this ibpoe to update the point bpo.
          ego EG_edg = msh.CAD.cad2edg[irefe];
          double tedg = msh.bpo2rbi(ibpoe,0);

          // fetch or create face ibpoi
          int ibpoi = msh.poi2ebp(cav.ipins, 2, iface, -1);
          CPRINTF1("   - existing ibpoi ? {}\n",ibpoi);
          if(ibpoi < 0) ibpoi = msh.newbpotopo(cav.ipins,2,iface);
          else CPRINTF1("   - old (u,v) = {} {} \n"
                        ,msh.bpo2rbi(ibpoi,0),msh.bpo2rbi(ibpoi,1));

          int icode = EG_getEdgeUV(EG_fac, EG_edg, isens, 
                                   tedg, msh.bpo2rbi[ibpoi]);
          CPRINTF1("   - new (u,v) = {} {}\n",
                   msh.bpo2rbi(ibpoi,0),msh.bpo2rbi(ibpoi,1));
          METRIS_ENFORCE_MSG(icode == 0,"EG_getEdgeUV error {}", icode);

          break;
        }
        if(found_update) break;
      }
      METRIS_ASSERT(found_update);
    }// for iface

  }else{
    CPRINTF1(" - skip ipins update, point dim = {}\n",msh.getpoitdim(cav.ipins));
  }// if ipins_uv



  if(msh.curdeg == 1){
    CPRINTF1("-- END update_bpois_newp, no HO nodes\n");
    return 0;
  };

  // What folllows is only for the HO nodes.

  // For line points, regenerate (u,v) from t. 
  if(msh.isboundary_faces() && msh.CAD()){ 
    // No updates necessary otherwise
    METRIS_ASSERT(msh.idim == 3);

    // The update makes some assumptions when it comes to HO points, namely 
    // if they belong to both an edge and a face, then the edge belongs to the 
    // face.
    for(int ipoin = npoi0; ipoin < msh.npoin; ipoin++){

      int ibpoi = msh.poi2bpo[ipoin];
      if(ibpoi < 0) continue;
      int ityp = msh.bpo2ibi(ibpoi,1);
      if(ityp >= 2) continue;

      // Cannot have created a corner. 
      METRIS_ASSERT(ityp == 1);

      // To use EG_getEdgeUV, we need to provide the sign (direction) of the curve 
      // in the face. In EGADS convention, material is to the left of the curve
      // when positively walked.
      // Meaning if the normal is going towards us, edge is going "up", 
      // then triangles to the left that share the edge, have 
      // the height (or the other edges) with positive scalar product 
      // to the vector product between normal and edge... 

      // edge tangent = e 
      // normal = n
      // h = project 3rd triangle vertex (not on edge) -> 3rd triangle vertex
      // compute l = n x e (vector product)
      // if (l,h) > 0 triangle is left -> positive
      //          < 0 triangle is right -> negative 

      // We can simplify all of this and base it only on a topological criterion:
      // is the triangle edge ordered the same as the edge edge?
      // If so, then the edge going up, triangle is trigo, normal is going towards
      // us, this yields positive edge. So positive = same ordering. 


      int iedge = msh.bpo2ibi(ibpoi,2);
      METRIS_ASSERT(iedge >= 0);
      int ip1 = msh.edg2poi(iedge,0);
      int ip2 = msh.edg2poi(iedge,1);


      // Get EG object for this edge
      double tt = msh.bpo2rbi(ibpoi,0);
      int irefe = msh.edg2ref[iedge];
      ego EG_edg = msh.CAD.cad2edg[irefe];
      METRIS_ASSERT(EG_edg != NULL);

      for(ibpoi = msh.bpo2ibi(ibpoi,3); ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
        if(msh.bpo2ibi(ibpoi,1) != 2) continue;

        // Found a triangle, now update.
        int iface = msh.bpo2ibi(ibpoi,2);
        int ireff = msh.fac2ref[iface];
        ego EG_fac  = msh.CAD.cad2fac[ireff];

        METRIS_ASSERT(EG_fac != NULL);

        // Base on ordering
        int ied = getedgfac(msh,iface,ip1,ip2);
        METRIS_ASSERT(ied >= 0);
        int sg = 1;
        if(ip1 != msh.fac2poi(iface,lnoed2[ied][0])) sg = -1;

        if(ipoin == 71925){
          fmt::print(" debug 71925 previous (u,v) = {} {}\n",
            msh.bpo2rbi(ibpoi,0),msh.bpo2rbi(ibpoi,1));
        }
        int icode = EG_getEdgeUV(EG_fac, EG_edg, sg, tt, msh.bpo2rbi[ibpoi]);
        if(ipoin == 71925){
          fmt::print(" debug using sg = {} tt = {} new (u,v) = {} {}\n",
            sg,tt,msh.bpo2rbi(ibpoi,0),msh.bpo2rbi(ibpoi,1));
          double uvdbg[4];
          icode = EG_getEdgeUV(EG_fac, EG_edg, -sg, tt, uvdbg);
          fmt::print(" debug if sg = -sg icode {} : {} {} \n",icode,uvdbg[0],uvdbg[1]);
        }

        static int nwarnprt1 = 0;
        if(nwarnprt1++ < 10) fmt::print("## DISABLE THIS DEBUG\n");

        if(icode != 0) METRIS_THROW_MSG("EG_getEdgeUV failed !");

      }

    } 
  }

  return 0;
}




}// namespace Metris