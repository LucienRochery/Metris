//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "msh_degelev.hxx"

#include "aux_topo.hxx"
#include "low_topo.hxx"
#include "aux_exceptions.hxx"
#include "aux_utils.hxx"

#include "Mesh/Mesh.hxx"


namespace Metris{

/*
Routine reworked by mid-June after surface topological handling was rewritten. 
It now seems more straightforward to proceed from edges, to triangles, to tetrahedra rather than the opposite. 
Let us also save on evaluations by first checking the neighbours and only then generating as needed. 
This implements topological and (u,v) guess init of boundary points.
*/
template<class MFT, int mshdeg, int tardeg>
void deg_elevate(Mesh<MFT> &msh){
  if(msh.idim != 2 && msh.idim != 3) METRIS_THROW_MSG(TODOExcept(), "Implement dim "<<msh.idim<<" degelev.");

  // This actually works, so why?
  //if constexpr(mshdeg != 1) METRIS_THROW_MSG(TODOExcept(), "Write into temporary arrays and update at the end.")

  int tagedl[6];
  int nshell = 0;
  intAr1 lshell(256);

  msh.edg2poi.set_stride(edgnpps[tardeg]);
  msh.fac2poi.set_stride(facnpps[tardeg]);
  msh.tet2poi.set_stride(tetnpps[tardeg]);

  //intAr2 newe2poi(msh.nedge,edgnpps[tardeg-2]);
  //intAr2 newf2poi(msh.nface,facnpps[tardeg-3]);
  //intAr2 newt2poi(msh.nelem,tetnpps[tardeg-4]);

  //dblAr2 metnew(msh.mpoin,msh.met.getnnmet());

  //int nvert = msh.idim == 2 ? msh.npoin * verppoi2[msh.curdeg] : 
  //                            msh.npoin * verppoi3[msh.curdeg] ;
  //int nedgi = msh.idim == 2 ? nvert * edgpver2 :
  //                            nvert * edgpver3 ;
  //int nfaci = msh.idim == 2 ? msh.nface :
  //                            nvert * edgpver3 ;


  //const int c1=(tardeg-1)*(tardeg>1), // pure edge 
  //          c2=(tardeg-1)*(tardeg-2)*(tardeg>2), // pure triangle 
  //          c3=(tardeg-1)*(tardeg-2)*(tardeg-3)*(tardeg>3); // pure tetra 
  //int mpnew = c1 * nedgi + c2 * nfaci + c3 * msh.nelem;
  //mpnew *= 1.5;

  //dblAr2 metnew(mpnew,msh.met.getnnmet());


  msh.met.setSpace(MetSpace::Log);
  FEBasis metbas0 = msh.met.getBasis();
  bool hasmet = (metbas0 != FEBasis::Undefined);
  msh.met.setBasis(FEBasis::Lagrange);

  //msh.setBasis(FEBasis::Bezier);

  msh.tag[0]++;

  CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
    
    double bary[4];
    int lbpoi[facnpps[mshdeg]];
    double newpt[tetnpps[tardeg]][gdim];
    double newuv[tetnpps[tardeg]][2];
    constexpr int nnmet = (gdim*(gdim+1))/2;
    double newmtl[tetnpps[tardeg]][nnmet];
    

    int lbpon[facnpps[tardeg]];
    
    constexpr int nppf0 = (facnpps[tardeg] - 3 - 3*(edgnpps[tardeg]-2));
    
    /* PHASE 1: EDGES*/
    for(int iedge = 0; iedge < msh.nedge; iedge++){
      if(isdeadent(iedge,msh.edg2poi)){
        for(int ii = 2; ii < edgnpps[tardeg]; ii++) msh.edg2poi(iedge,ii) = -1;
        continue;
      }
    
      // Always need to generate the points first because we'll overwrite nodes
      // We could get away with in place interior points but only in Lagrange
      gen_nodes_degelev<mshdeg,tardeg,1,gdim>(msh,iedge,newpt[0]);

      if(msh.isboundary_edges()){
        // Similarly, we need the lbpois now 
        getbpois<mshdeg,1>(msh,iedge,lbpoi);
        // Evaluate CAD coordinate. 
        gen_uvs_degelev<mshdeg,tardeg,1,2>(msh,lbpoi,newuv[0]);
      }


      // Evaluate metrics
      if(hasmet) msh.met.template getMetNodes<gdim,1,mshdeg,tardeg>(iedge,newmtl[0]);

      int irn0 = 2; 
      for(int irnk = irn0; irnk < edgnpps[tardeg]; irnk++){

        int ipnew;
        int ibnew;
        if(irnk < edgnpps[mshdeg]){
          // Inherits poi2ent
          ipnew = msh.edg2poi(iedge,irnk);
          // Inherits attachment to this edge, triangles, etc. Same topological make. 
          ibnew = msh.poi2bpo[ipnew];
        }else{
          ipnew = msh.newpoitopo(1,iedge);
          ibnew = msh.newbpotopo(ipnew,1,iedge);
        }

        if(ipnew == 1517) printf("## DEBUG 1517 CREATED BY EDGE %d \n",iedge);

        for(int ii = 0; ii < gdim; ii++) msh.coord(ipnew,ii) = newpt[irnk][ii];


        msh.edg2poi(iedge,irnk) = ipnew;
        //newe2poi[iedge][irnk-2] = ipnew;

        msh.bpo2rbi(ibnew,0) = newuv[irnk][0];
        msh.bpo2rbi(ibnew,1) = newuv[irnk][1];

        // Do not convert to expmet, we will expmet everyone at the end
        if(hasmet) for(int ii = 0; ii < nnmet; ii++) msh.met(ipnew,ii) = newmtl[irnk][ii];
      }

    }
    
    
    
    /* PHASE 2: TRIANGLES */
    for(int iface = 0; iface < msh.nface; iface++){
      if(isdeadent(iface,msh.fac2poi)){
        for(int ii = 3; ii < facnpps[tardeg]; ii++) msh.fac2poi(iface,ii) = -1;
        continue;
      } 
      msh.fac2tag(0,iface) = msh.tag[0];
    
      // Some of the boundary points are not correlated to a new point. 
      // Let's simply store them in a facnpps sized array
      for(int irnk = 0; irnk < facnpps[tardeg]; irnk++)
        lbpon[irnk] = -1;
    
      // Always need to generate the points first because we'll overwrite nodes
      // We could get away with in place interior points but only in Lagrange
      gen_nodes_degelev<mshdeg,tardeg,2,gdim>(msh,iface,newpt[0]);



      if(msh.isboundary_faces()){
        // Similarly, we need the lbpois now 
        getbpois<mshdeg,2>(msh,iface,lbpoi);
        // Evaluate CAD coordinate. 
        gen_uvs_degelev<mshdeg,tardeg,2,2>(msh,lbpoi,newuv[0]);
      }

      // Evaluate metrics
      // Figure out this .template mistery, compiler error (clang) otherwise:
      // error: missing 'template' keyword prior to dependent template name
      if(hasmet) msh.met.template getMetNodes<gdim,2,mshdeg,tardeg>(iface,newmtl[0]);
      
      // Start with interior points if any.
      if constexpr(tardeg >= 3){
        int irn0     = 3 + 3*(edgnpps[tardeg]-2);
        int irn0_prv = 3 + 3*(edgnpps[mshdeg]-2);

        for(int irnk = irn0; irnk < facnpps[tardeg]; irnk++){

          int ipnew, ibnew;
          if(irn0_prv + irnk - irn0 < facnpps[mshdeg]){
            METRIS_ASSERT(tardeg > 3);
            ipnew = msh.fac2poi[iface][irn0_prv + irnk - irn0];
            ibnew = msh.poi2bpo[ipnew];
          }else{
            ipnew = msh.newpoitopo(2,iface);
            ibnew = msh.newbpotopo(ipnew,2,iface);
          }

          if(ipnew == 1517) printf("## DEBUG 1517 CREATED BY FACE %d \n",iface);

          for(int ii = 0; ii < gdim; ii++) msh.coord(ipnew,ii) = newpt[irnk][ii];

          msh.fac2poi(iface,irnk) = ipnew;

          lbpon[irnk] = ibnew;
    
          // Now the metric (if back == front)
          // Do not convert to expmet, we will expmet everyone at the end
          //if(hasmet) for(int ii = 0; ii < nnmet; ii++) metnew[ipnew][ii] = newmtl[irnk][ii];// msh.met(ipnew,ii) = newmtl[irnk][ii];
          if(hasmet) for(int ii = 0; ii < nnmet; ii++) msh.met(ipnew,ii) = newmtl[irnk][ii];
        }
      }
    
    
      // Edge points. First check if neighbour init. Then check if geo edge. 
      // ATTENTION !!!!!!!!!!!!!!!
      // Walk backwards so as not to overwrite indices of previous degree points !
      // It might be safer to copy over into an array before any modifs, but this warning 
      // will do for now. 
      // ATTENTION !!!!!!!!!!!!!!!
      for(int ied = 2; ied >= 0; ied--){
        int ifac2 = msh.fac2fac(iface,ied);
        // We must start with checking if an edge is here because of bpoi links
        int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
        int ip2 = msh.fac2poi(iface,lnoed2[ied][1]);

        int irn0 = 3 +  ied   *(edgnpps[tardeg]-2);
        int irn1 = 3 + (ied+1)*(edgnpps[tardeg]-2);
        int irn0_prv = 3 +  ied   *(edgnpps[mshdeg]-2);
        int irn1_prv = 3 + (ied+1)*(edgnpps[mshdeg]-2);

        int iedge = msh.facedg2glo(iface,ied);
        if(iedge >= 0){
          cpy_gloedg2facedg<tardeg>(msh,iedge,iface,ied);

          if(msh.idim > 2){
            for(int irnk = irn0; irnk < irn1; irnk++){
              if(irn0_prv + irnk - irn0 < irn1_prv){
                lbpon[irnk] = lbpoi[irn0_prv + irnk - irn0];
                msh.bpo2ibi[lbpon[irnk]][0] = msh.fac2poi(iface,irnk);
              }else{
                msh.newbpotopo(msh.fac2poi(iface,irnk),2,iface);
                lbpon[irnk] = msh.nbpoi - 1;
              }
            }
          }
        }else if(ifac2 < 0){
          METRIS_THROW_MSG(TopoExcept(),"ERROR NO EDGE BUT NON MANIFOLD OR SURF BDRY")

        // VV From here on, manifold neighbour
        }else if(msh.fac2tag(0,ifac2) >= msh.tag[0]){
          // Copy from neighbour
          int ied2 = getedgfac(msh,ifac2,ip1,ip2);
          cpy_facedg2facedg<tardeg>(msh,ifac2,ied2,iface,ied);
        }else{
          // Create new
          for(int irnk = irn0; irnk < irn1; irnk++){

            int ipnew, ibnew;
            if(irn0_prv + irnk - irn0 < irn1_prv){
              ipnew = msh.fac2poi[iface][irn0_prv + irnk - irn0];
              ibnew = msh.poi2bpo[ipnew];
            }else{
              ipnew = msh.newpoitopo(2,iface);
              ibnew = msh.newbpotopo(ipnew,2,iface);
            }

            for(int ii = 0; ii < gdim; ii++) msh.coord(ipnew,ii) = newpt[irnk][ii];


            msh.fac2poi(iface,irnk) = ipnew;
            //newf2poi[iface][irnk-3] = ipnew;
            // Edge points not on an edge don't need to be linked
            lbpon[irnk] = ibnew;
            // Now the metric (if back == front)
            // Do not convert to expmet, we will expmet everyone at the end
            //if(hasmet) for(int ii = 0; ii < nnmet; ii++) metnew[ipnew][ii] = newmtl[irnk][ii];//msh.met(ipnew,ii) = newmtl[irnk][ii];
            if(hasmet) for(int ii = 0; ii < nnmet; ii++) msh.met(ipnew,ii) = newmtl[irnk][ii];
          }
        }
      }
    
      if(msh.isboundary_faces()){
        // Update rbpois 
        for(int irnk = 3; irnk < facnpps[tardeg]; irnk++){
          //int ipoin = msh.fac2poi(iface,irnk);
          //int ipoin = newf2poi[iface][irnk-3];
          // Boundary points were created in the same order as points
          int ibpoi = lbpon[irnk];
          if(ibpoi < 0) continue;
          METRIS_ASSERT(ibpoi >= 0 && ibpoi < msh.nbpoi);
          //METRIS_ASSERT_MSG(msh.bpo2ibi(ibpoi,0) == ipoin, "\nError at irnk = "<<irnk<<" ipoin = "<<ipoin<<" ibpoi "<<ibpoi<<" bpo2ibi = "<<msh.bpo2ibi(ibpoi,0));
          msh.bpo2rbi(ibpoi,0) = newuv[irnk][0];
          msh.bpo2rbi(ibpoi,1) = newuv[irnk][1];
        }
      }

    }
    
    




    /* PHASE 3: TETRAHEDRA */
    if constexpr (gdim >= 3){ // Note: constexpr avoids compilation of invalid template param combinations such as tdim = 3, gdim = 2
      for(int ielem = 0; ielem < msh.nelem; ielem++){
        if(isdeadent(ielem,msh.tet2poi)){
          for(int ii = 4; ii < tetnpps[tardeg]; ii++) msh.tet2poi(ielem,ii) = -1;
          continue;
        } 
        msh.tet2tag(0,ielem) = msh.tag[0];
        gen_nodes_degelev<mshdeg,tardeg,3,gdim>(msh,ielem,newpt[0]);

        // Evaluate metrics
        if(hasmet) msh.met.template getMetNodes<gdim,3,mshdeg,tardeg>(ielem,newmtl[0]);
    
        if constexpr(tardeg >= 4){
          // Interior points
          int irnk0 = 4 + 6 * (edgnpps[tardeg] - 2)
                        + 4 * nppf0;

          for(int irnk = irnk0; irnk < tetnpps[tardeg]; irnk++){
            int ipnew = msh.newpoitopo(3,ielem);

            for(int ii = 0; ii < gdim; ii++) msh.coord(ipnew,ii) = newpt[irnk][ii];
            

            //newt2poi[ielem][irnk-4] = ipnew;
            msh.tet2poi(ielem,irnk) = ipnew;
            bary[0] = ordtet.s[tardeg][irnk][0]/(tardeg*1.0);
            bary[1] = ordtet.s[tardeg][irnk][1]/(tardeg*1.0);
            bary[2] = ordtet.s[tardeg][irnk][2]/(tardeg*1.0);
            bary[3] = ordtet.s[tardeg][irnk][3]/(tardeg*1.0);
            // Now the metric (if back == front)
            // Do not convert to expmet, we will expmet everyone at the end
            //if(hasmet) for(int ii = 0; ii < nnmet; ii++) metnew[ipnew][ii] = newmtl[irnk][ii];//msh.met(ipnew,ii) = newmtl[irnk][ii];
            if(hasmet) for(int ii = 0; ii < nnmet; ii++) msh.met(ipnew,ii) = newmtl[irnk][ii];
          }
        }
    
        // Start with faces. If neighbour or boundary face, copy. 
        // Tag copied face edges as done. 
        // Otherwise, create face interior points, but not edge. 
        // After faces, loop over edges, check edgHshTab, check shell. 
        // If found, copy, otherwise create edge nodes. 
    
        // In the P2 case, this saves 3 shell calls if a neighbour is HO
    
        for(int ii = 0;ii<6;ii++) 
          tagedl[ii] = 0;
    
        for(int ifa = 0; ifa < 4; ifa++){
          int iele2 = msh.tet2tet(ielem,ifa);
          if(iele2 < 0){
            int iface = msh.tetfac2glo(ielem,ifa);
            if(iface < 0) METRIS_THROW_MSG(TopoExcept(),
              "FACE MISSING FROM HASH TABLE")

            cpy_glofac2tetfac<tardeg>(msh,iface,ielem,ifa);
            continue;
          }
          
          // In this case, we have a neighbour. Is it already HO?
          if(msh.tet2tag(0,iele2) >= msh.tag[0]){
            // Yes, copy over
            int ifa2 = getneitet(msh,ielem,iele2);
            METRIS_ASSERT(ifa2>=0);
            cpy_tetfac2tetfac<tardeg>(msh,iele2,ifa2,ielem,ifa);
            // Tag edges as already done with
            for(int i = 0; i < 3 ; i++)
              tagedl[ledfa3[ifa][i]] = 1;
          }else{
            // No, create interior points
            int irnk0 = 4 + 6*(edgnpps[tardeg]-2) + ifa*nppf0;
            int irnk1 = 4 + 6*(edgnpps[tardeg]-2) + (ifa+1)*nppf0; 

            for(int irnk = irnk0; irnk < irnk1; irnk++){
              int ipnew = msh.newpoitopo(3,ielem);
              for(int ii = 0; ii < gdim; ii++) msh.coord(ipnew,ii) = newpt[irnk][ii];

              msh.tet2poi(ielem,irnk) = ipnew;
              // Now the metric (if back == front)
              // Do not convert to expmet, we will expmet everyone at the end
              //if(hasmet) for(int ii = 0; ii < nnmet; ii++) metnew[ipnew][ii] = newmtl[irnk][ii];//msh.met(ipnew,ii) = newmtl[irnk][ii];
              if(hasmet) for(int ii = 0; ii < nnmet; ii++) msh.met(ipnew,ii) = newmtl[irnk][ii];
            }
          }
        } 
    
        for(int ied = 0; ied < 6; ied++){
          if(tagedl[ied] == 1) continue;
          int ip1 = msh.tet2poi(ielem,lnoed3[ied][0]);
          int ip2 = msh.tet2poi(ielem,lnoed3[ied][1]);
          
          // First check if in global hash table
          int iedge = msh.tetedg2glo(ielem,ied);
          if(iedge >= 0){
            // This edge has already been initialized as a global edge
            cpy_gloedg2tetedg<tardeg>(msh,iedge,ielem,ied);
            continue;
          }
    
    
          int iopen;
          shell3(msh, ip1 ,ip2 ,ielem ,&nshell ,lshell, &iopen); // i = dum
    
    
          // If the shell is open, that means the edge is on the boundary. 
          // We already know it's not a geometric edge. 
          // We must fetch one of the boundary triangles closing the shell.
          if(iopen >= 0){
            int iele2 = iopen;
            //int ied2 = getedgtet(msh,iele2,ip1,ip2);
            int ifnd = 0;
            for(int ifa2 = 0; ifa2 < 4; ifa2++){
              int ip = msh.tet2poi(iele2,ifa2);
              if(ip == ip1 || ip == ip2) continue;
              if(msh.tet2tet(iele2,ifa2) >= 0) continue;
              // This neighbour is a face that contains the edge
              ifnd = 1;
              int iface = msh.tetfac2glo(iele2,ifa2);
              if(iface < 0) METRIS_THROW_MSG(TopoExcept(),
                                                 "FACE MISSING FROM HASH TABLE")

              int iedf = getedgfac(msh,iface,ip1,ip2);
              cpy_facedg2tetedg<tardeg>(msh,iface,iedf,ielem,ied);
            }
            if(ifnd == 0)METRIS_THROW_MSG(TopoExcept(),
                                "FAILED TO FIND BOUNDARY FACE WITH THE EDGE")
          }else{
            int ifnd = 0;
            for(int ishell = 0;ishell < nshell;ishell++){
              int iele2 = lshell[ishell];
              if(iele2 == ielem) continue;
              if(msh.tet2tag(0,iele2) < msh.tag[0]) continue;
    
              int ied2 = getedgtet(msh,iele2,ip1,ip2);
              METRIS_ASSERT(ied2 >= 0);
              cpy_tetedg2tetedg<tardeg>(msh,iele2,ied2,ielem,ied);
              METRIS_ASSERT(msh.tet2poi[ielem][4+ied*(tardeg-1)] == 
                     msh.tet2poi[iele2][4+ied2*(tardeg-1)]
                    || msh.tet2poi[ielem][4+ied*(tardeg-1)]
                    == msh.tet2poi[iele2][4+ied2*(tardeg-1) + tardeg-2]); 
              ifnd = 1;
              break; 
            } // End shell loop
    
    
            if(ifnd == 0){ // No shell element could help us
              int irn0 = 4 + (ied  )*(edgnpps[tardeg]-2);
              int irn1 = 4 + (ied+1)*(edgnpps[tardeg]-2);
              METRIS_ASSERT(irn1 <= tetnpps[tardeg]);
    
              for(int irnk = irn0; irnk < irn1;irnk++){
                int ipnew = msh.newpoitopo(3,ielem);
                msh.tet2poi(ielem,irnk) = ipnew;
                //newt2poi[ielem][irnk-4] = ipnew;
                for(int i = 0; i < gdim; i++) msh.coord(ipnew,i) = newpt[irnk][i];

                // Now the metric (if back == front)
                // Do not convert to expmet, we will expmet everyone at the end
                if(hasmet) for(int ii = 0; ii < nnmet; ii++) msh.met(ipnew,ii) = newmtl[irnk][ii];
              }
            }
          }
    
        }
    
      }
    }

  }}CT_FOR1(gdim);

  msh.curdeg = tardeg;
  if(hasmet){
    msh.met.setSpace(MetSpace::Exp);
    msh.met.setBasis(metbas0);
  }

  //if(msh.param->iverb >= 1){
  //  FEBasis ibas0 = msh.getBasis();
  //  msh.setBasis(FEBasis::Lagrange);
  //  writeMesh("degelev",msh);
  //  msh.setBasis(ibas0);
  //  if(msh.idim == 2){
  //    maximizeMetCcoef<MFT,2,2,tardeg>(msh, OptDoF::HO, LPMethod::IPM, LPLib::alglib);
  //  }else{
  //    maximizeMetCcoef<MFT,3,3,tardeg>(msh, OptDoF::HO, LPMethod::IPM, LPLib::alglib);
  //  }
  //}

}
// Explicit instantiation.
#include <src/msh_degelev.ixx>


} // End namespace
