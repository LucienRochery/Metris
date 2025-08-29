//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../cavity/msh_cavity.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../aux_topo.hxx"
#include "../ho_constants.hxx"
#include "../low_geo/measure.hxx"
#include "../quality/low_metqua.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/aux_misc.hxx"
#include "../io_libmeshb.hxx"

namespace Metris{


template <class MFT, int ideg>
int reconnect_tetcav(Mesh<MFT> &msh, 
                     MshCavity& cav, 
                     CavOprOpt  &opts, 
                     CavOprInfo &info,
                     int nfac0, double *qmax, int ithread){

  const int nctet = cav.lctet.get_n();
  if(nctet <= 0) return 0;

  bool check_qua = (opts.qmax_nec > 0 && msh.get_tdim() == 3)
                || (opts.qmax_suf > 0 && msh.get_tdim() == 3)
                || (opts.qmax_iff > 0 && msh.get_tdim() == 3);

  //// When would we not though?
  //bool check_val = opts.fast_reject 
  //              || opts.max_increase_cav_geo <= 0
  //              || check_qua;

  const int qpnorm = msh.param->opt_pnorm;
  *qmax = -1;
  info.qcav3 = qpnorm == 0 ? -1 : 0;

  // We need to increment this independently per element, avoid conflict with
  // ithread
  int cfatag = msh.tag[ithread];


  GETVDEPTH(msh.param);
  CPRINTF1(" - START reconnect_tetcav check_qua = {} qpnorm {} qpower {}\n",
            check_qua,qpnorm,msh.param->opt_power);

  const int nele0 = msh.nelem;

  // In the future, handle connex components here. For now, unique ref. 
  // If the cavity has two+ tet connex components, then ipins is necessarily
  // a tdim 2 or 1 point. It is on the interface beetween the connex cpts. 
  // In that case, the faces nfac0 through nface do not generate tetrahedra. 
  // Hence the only case where tets are created with new faces is if point 
  // is interior. In that case it is quite simple, the tetrahedra are oriented
  // same as face with ipins in first position, and we are always the first fac2tet.
  for(int ielem : cav.lctet){
    msh.tet2tag(ithread, ielem) = msh.tag[ithread];
  }

  // Keep track of faces that contain ipins. These are either new or prior
  // cavity boundary faces, or new cavity interior faces. 
  // As they all share ipins, it is enough to hash the remaining two vertices. 
  // Such faces hit at most 2 elements; we need only keep track of one, as 
  // when we are the second, we update info and can trash the face. 
  HshTab_I2I2 facHsh; 
  for(int iface = nfac0; iface < msh.nface; iface++){
    int iver = msh.template getverfac<1>(iface, cav.ipins);
    METRIS_ASSERT(iver >= 0);
    int ip1 = msh.fac2poi(iface,lnoed2[iver][0]);
    int ip2 = msh.fac2poi(iface,lnoed2[iver][1]);
    auto key = stup2(ip1,ip2);
    facHsh[key] = {2, iface};
  }
  //// Cavity boundary need not be comprised of new triangles:
  //for(int iele0 : cav.lctet){
  //  for(int ifa0 = 0; ifa0 < 4; ifa0++){
  //    int ienei = msh.tet2tet(iele0,ifa0);
  //    if(ienei >= 0 && msh.tet2tag(ithread,ienei) >= msh.tag[ithread]) continue;

  //    bool iskip = true;
  //    for(int ii = 0; ii < 3; ii++){
  //      if(msh.tet2poi(iele0, lnoed3[ifa0][ii]) != cav.ipins) continue;
  //      iskip = false;
  //      break;
  //    }
  //    if(iskip) continue;
  //  }// for ifa0
  //}// for iele0

  // To check if all points are on the boundary. 
  const int pdim_ipins = msh.getpoitdim(cav.ipins);

  // Next loop over cavity boundary faces to generate new tets. If there are 
  // new faces, they do not generate tetrahedra, as those would have ipins twice.
  for(int iele0 : cav.lctet){
    INCVDEPTH(msh.param);
    for(int ifa0 = 0; ifa0 < 4; ifa0++){
      int ienei = msh.tet2tet(iele0,ifa0);
      if(ienei >= 0 && msh.tet2tag(ithread,ienei) >= msh.tag[ithread]){
        CPRINTF1(" - skip tetra creation from {} face {} due to cavity neighbour {} \n",
                 iele0, ifa0, ienei);
        continue;
      }

      // Also skip if face contains ipins already. 
      bool iskip = false;
      for(int ii = 0; ii < 3; ii++){
        if(msh.tet2poi(iele0, lnofa3[ifa0][ii]) != cav.ipins) continue;
        iskip = true;
        break;
      }
      if(iskip){
        CPRINTF1(" - skip tetra creation from {} face {} which contains ipins\n",
                 iele0, ifa0);
        continue;
      }

      // Additionally, skip any faces that coincide with triangles which were
      // in the cavity. These have disappeared from the mesh. 
      {// namespace
        auto key = stup3(msh.tet2poi(iele0, lnofa3[ifa0][0]),
                         msh.tet2poi(iele0, lnofa3[ifa0][1]),
                         msh.tet2poi(iele0, lnofa3[ifa0][2]));
        auto tt = msh.facHshTab.find(key);
        if(tt != msh.facHshTab.end()){
          int iface = tt->second;
          if(msh.fac2tag(ithread,iface) >= msh.tag[ithread]){
            CPRINTF1(" - skip tetra creation from {} face {} is glo face {} in cavity\n",
                     iele0, ifa0 ,iface);
            continue;
          }
        }
      }


      // Boundary face, create tet. 
      int ielen = msh.nelem;
      msh.set_nelem(msh.nelem+1);

      msh.tet2tet(ielen, 0) = -1;
      msh.tet2tet(ielen, 1) = -1;
      msh.tet2tet(ielen, 2) = -1;
      msh.tet2tet(ielen, 3) = -1;
      msh.tet2tet(ielen, ifa0) = ienei;
      msh.tet2poi(ielen, ifa0) = cav.ipins;
      msh.tet2poi(ielen, lnofa3[ifa0][0]) = msh.tet2poi(iele0, lnofa3[ifa0][0]);
      msh.tet2poi(ielen, lnofa3[ifa0][1]) = msh.tet2poi(iele0, lnofa3[ifa0][1]);
      msh.tet2poi(ielen, lnofa3[ifa0][2]) = msh.tet2poi(iele0, lnofa3[ifa0][2]);
      
      CPRINTF1(" - new tetra = {} from iele0 = {} ifa = {} vertices: {} {} {} {}\n",
               ielen,iele0,ifa0,msh.tet2poi(ielen,0),msh.tet2poi(ielen,1),
               msh.tet2poi(ielen,2),msh.tet2poi(ielen,3));

      double meas0;
      if(!isvalideltP1<3,3>(msh, ielen,NULL,&meas0)){
        CPRINTF1(" - iflat ! return ip1 ip2 ip3 ip4 = {} {} {} {} meas = {:15.7e} \n", 
                 msh.tet2poi(ielen,0),msh.tet2poi(ielen,1),msh.tet2poi(ielen,2),
                 msh.tet2poi(ielen,3),meas0); 
        return CAV_ERR_FLATTET;
      }
      CPRINTF2(" - new tetra volume {} \n",meas0);


      // If the tet has all vertices on boundary, check faces orientations.
      if(pdim_ipins < 3){

        for(int ii = 0; ii < 3; ii++){
          int pdim = msh.getpoitdim(msh.tet2poi(ielen, lnofa3[ifa0][ii]));
          if(pdim < 3) continue;
          CPRINTF1(" face point {} has dim 3 -> skip\n",msh.tet2poi(ielen, lnofa3[ifa0][ii]));
          goto check_bdry_done;
        }

        int itmp1 = iverb__;
        int itmp2 = ivdepth__;
        bool ipause = false;
        //if(msh.tet2poi(ielen,0) == 30 && msh.tet2poi(ielen,1) == 478
        //  && msh.tet2poi(ielen,0) == 29 && msh.tet2poi(ielen,1) == 118
        //  || ielen == 4353){
        //  printf("\n\n ## DEBUG 4 BDRY VERTEX TET CASE MAX PRINTS ielen {} vertices ",ielen);
        //  intAr1(4,msh.tet2poi[ielen]).print();
        //  writeMesh("debug",msh);
        //  ipause = true;
        //  iverb__ = 5;
        //  ivdepth__ = 10;
        //  wait();
        //}

        // The previous test was more lenient: it allowed tetrahedra to have
        // all four vertices on the boundary, provided that the tet faces 
        // were reverse oriented to the boundary faces they supported. 
        // This meant the tet was inside the domain. 
        // One issue is a tetrahedron could be created with only one (perhaps even none)
        // face on the boundary, and later be "uncovered". 
        // New criterion is more conservative but does not depend on future topology
        // we simply check if no four vertices are on the same face reference.

        // tag[ithread] has been incremented for other uses already, and no
        // conflict with cfa2tag:
        cfatag++;
        cav.maxtag = MAX(cav.maxtag,cfatag + 4);
        for(int iver = 0; iver < 4; iver++){
          //INCVDEPTH(msh.param);
          int ipoin = msh.tet2poi(ielen, iver);
          for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
            int tdimb = msh.bpo2ibi(ibpoi,1);
            if(tdimb != 2) continue;
            int iface = msh.bpo2ibi(ibpoi,2);
            int ireff = msh.fac2ref[iface];
            int nseenm1 = msh.cfa2tag(ithread,ireff) - cfatag;
            CPRINTF1(" - check ver {} ipoin {} face {} ref {} already seen {}x\n",
                     iver,ipoin,iface,ireff,nseenm1+1);
            // Same ref can be seen twice by same vertex, meaning it has already 
            // been updated.
            METRIS_ASSERT_MSG(nseenm1 <= iver, "nseenm1 = {} iver = {}", nseenm1, iver);

            // We must be careful to skip refs that have been skipped by at 
            // least one prior point, because our method of counting is 
            // setting tag to tag + iver. Otherwise, if only the last point sees
            // a face ref, we'd think all did.
            if(nseenm1 + 1 < iver && iver != 0){
              CPRINTF1(" - inactive ref -> skip\n");
              continue;
            }

            // If has already been seen 3 times before,
            if(nseenm1 + 1 == 3 && iver == 3){
              CPRINTF1(" # error 4 points on same face ref\n");
              return CAV_ERR_BDRYTET2;
            }
            // An edge point can see the same face ref several times. We mustn't
            // increment.
            msh.cfa2tag(ithread,ireff) = cfatag + iver;
          }
        }
        cfatag += 4;

        #if 0
        // All vertices are on the boundary. Go over faces and compare orientation
        bool nogood = false;
        for(int ifal = 0; ifal < 4; ifal++){
          int iface = msh.tet2fac(ielen, ifal);
          CPRINTF1(" - tet face {} -> glo face {}\n",ifal, iface);
          if(iface < 0) continue;
          if(msh.fac2tet(iface,0) >= 0 && msh.fac2tet(iface,1) >= 0){
            if(msh.param->dbgfull) 
              METRIS_THROW_MSG("TODO: Handle internal surface case here");
            CPRINTF1(" # REJECT: tet face {} matches glo face {} which is internal\n",
              ifal, iface);
            return CAV_ERR_BDRYTET2;
          }
          int ip1 = msh.tet2poi(ielen, lnofa3[ifal][0]);
          int ip2 = msh.tet2poi(ielen, lnofa3[ifal][1]);
          int ip3 = msh.tet2poi(ielen, lnofa3[ifal][2]);
          int jp1 = msh.fac2poi(iface, 0);
          int jp2 = msh.fac2poi(iface, 1);
          int jp3 = msh.fac2poi(iface, 2);
          CPRINTF1(" - glo face {} = {} {} {}: check opposite orientation\n",
                  iface, jp1, jp2, jp3);
          if(ip1 == jp1){
            if(ip2 == jp2) nogood = true;
          }else if(ip1 == jp2){
            if(ip2 == jp3) nogood = true;
          }else if(ip1 == jp3){
            if(ip2 == jp1) nogood = true;
          }else{
            METRIS_THROW_MSG( "Tet local <-> global face mismatch");
          }
          if(nogood) break;
        }
        if(nogood){
          CPRINTF1(" # REJECT: tet face matches glo face orientation\n");
          return CAV_ERR_BDRYTET;
        }
      #endif
        if(ipause){
          iverb__ = itmp1;
          ivdepth__ = itmp2;
          printf("## DEBUG CASE ACCEPTED \n\n");
          wait();
        }
      }
      check_bdry_done:


      if(check_qua){
        // Regardless of degree, verify underlying P1 element is decent enough
        double quael;
        if(opts.cache_tetra_quality){
          auto key = stup4(msh.tet2poi(ielen,0),msh.tet2poi(ielen,1),
                           msh.tet2poi(ielen,2),msh.tet2poi(ielen,3));
          auto tt = cav.qtetr.find(key);
          if(tt != cav.qtetr.end()){
            quael = tt->second;
            CPRINTF2(" - found cached quality\n");
          }else{
            quael = metqua<MFT,3,3>(msh,AsDeg::P1,AsDeg::P1,ielen,1.0);
            cav.qtetr[key] = quael;
          }
        }else{
          quael = metqua<MFT,3,3>(msh,AsDeg::P1,AsDeg::P1,ielen,1.0);
        }
        CPRINTF1(" - new tetra {} = {} {} {} {} from {} conf error = {} \n",
           ielen,
           msh.tet2poi(ielen,0), msh.tet2poi(ielen,1), msh.tet2poi(ielen,2),
           msh.tet2poi(ielen,3) ,iele0,quael);
        *qmax = MAX(quael, *qmax);
        if(qpnorm == 0){
          info.qcav3 = MAX(info.qcav3, quael);
        }else{
          info.qcav3 += pow(abs(quael), qpnorm);
        }
        if(quael > opts.qmax_nec && opts.qmax_nec > 0.0){
          CPRINTF1(" # quael = {} > {} = qmax_nec -> reject\n",quael, opts.qmax_nec);
          return CAV_ERR_QMAXNEC; // Run rejected
        }
        if(quael > opts.qmax_iff && opts.qmax_iff > 0.0) return CAV_ERR_QMAXIFF; // Run rejected
        if(quael < 0) return CAV_ERR_QFACNEG;  // Run rejected
      }

      msh.tet2ref[ielen] = msh.tet2ref[iele0];


      // Copy ctrl pts; capping at ideg-1 and preventing ii == jj == 0 is enough
      // to skip the vertices. 
      for(int ii = 0; ii <= ideg-1; ii++){
        for(int jj = (int)(ii == 0); jj <= ideg-1; jj++){
          int kk = ideg - ii - jj;

          int inode = ifa0 == 0 ? mul2nod(0,ii,jj,kk)
                    : ifa0 == 1 ? mul2nod(kk,0,ii,jj)
                    : ifa0 == 2 ? mul2nod(jj,kk,0,ii)
                                : mul2nod(ii,jj,kk,0);

          msh.tet2poi(ielen,inode) = msh.tet2poi(iele0, inode);
        }
      }

      // For the three faces impinging on ipins, we need to check if they exist
      // in the hash table. Otherwise create and then add them. 
      // Loop over opposite face edges. 
      for(int ied = 0; ied < 3; ied++){
        // Let's break this down. ifa is the node we want to avoid (ipins). 
        // lnoed[ied][0] and lnoed[ied][1] will span all combinations of 
        // 2 among 0,1,2
        // to this we add 1, and go over ifa + 1, ifa + 2 and ifa + 3 
        // modulo 4, which spans all pairs of vertices not ifa. 
        int ip1 = msh.tet2poi(ielen, (ifa0 + lnoed2[ied][0]+1)%4); 
        int ip2 = msh.tet2poi(ielen, (ifa0 + lnoed2[ied][1]+1)%4); 

        // This makes the face opposite ifa + ied. 
        int ifan = (ifa0 + ied + 1)%4;

        CPRINTF1(" - with ifa0 {} have ifan {} \n",ifa0, ifan);
        //CPRINTF1(" debug face ifa0 {} {} {} \n",msh.tet2poi(ielen,lnofa3[ifa0][0])
        //      ,msh.tet2poi(ielen,lnofa3[ifa0][1]),msh.tet2poi(ielen,lnofa3[ifa0][2]));
        //CPRINTF1(" debug face ifan {} {} {} \n",msh.tet2poi(ielen,lnofa3[ifan][0])
        //      ,msh.tet2poi(ielen,lnofa3[ifan][1]),msh.tet2poi(ielen,lnofa3[ifan][2]));


        auto key = stup2(ip1, ip2);
        auto tt = facHsh.find(key);
        if(tt != facHsh.end()){
          std::pair<int,int> hshfac = tt->second;
          int tdimf = hshfac.first;
          METRIS_ASSERT(tdimf == 2 || tdimf == 3);
          int ielef = hshfac.second;
          METRIS_ASSERT(ielef >= 0 && ielef < msh.nentt(tdimf));

          CPRINTF1(" - found internal or new boundary face tdim {} entt {} \n",
                   tdimf,ielef);

          // Copy the nodes from the entity. Could be dimension 2 or 3. 
          if(tdimf == 2){

            CPRINTF1(" - bdry face nodes {} {} {} \n",msh.fac2poi(ielef,0)
                     ,msh.fac2poi(ielef,1),msh.fac2poi(ielef,2)); 
            CPRINTF1("   match with {} {} {} \n",msh.tet2poi(ielen,lnofa3[ifan][0])
              ,msh.tet2poi(ielen,lnofa3[ifan][1]),msh.tet2poi(ielen,lnofa3[ifan][2]));


            cpy_glofac2tetfac<ideg>(msh, ielef, ielen, ifan);

            // This is a new face. Update the fac2tet
            METRIS_ASSERT(ielef >= nfac0);
            if(msh.fac2tet(ielef,0) < 0){
              // Nothing here yet
              msh.fac2tet(ielef,0) = ielen;
            }else if(msh.fac2tet(ielef,0) >= nele0){
              // 0 points to a new one, put in 1
              msh.fac2tet(ielef,1) = ielen;
            }else{
              METRIS_THROW_MSG( "New face points to dead element")
            }

          }else{
            METRIS_ASSERT(ielef >= nele0);
            // Copy from a tetrahedron. For this we need to know the ifa for 
            // the tet. 
            int ifaf = getfactet(msh, ielef, ip1, ip2, cav.ipins);
            METRIS_ASSERT(ifaf >= 0);

            cpy_tetfac2tetfac<ideg>(msh, ielef, ifaf, ielen, ifan);
            
            // Also update the neighbour. This is a cavity interior face and both are
            // new elements. 
            msh.tet2tet(ielen,ifan) = ielef;
            msh.tet2tet(ielef,ifaf) = ielen;
          }

        }else{// tt != facHsh.end()
          // Create control points and add this interior face to hash table. 
          // Loop over face indices all but vertices
          int idx[4] = {};
          for(int ii = 0; ii <= ideg-1; ii++){
            for(int jj = (int)(ii == 0); jj <= ideg-1; jj++){
              int kk = ideg - ii - jj;
              idx[lnofa3[ifan][0]] = ii;
              idx[lnofa3[ifan][1]] = jj;
              idx[lnofa3[ifan][2]] = kk;
              int inod3 = mul2nod(4,idx);
              int ipnew = msh.newpoitopo(3, ielen);
              msh.tet2poi(ielen,inod3) = ipnew;
              for(int mm = 0; mm < 3; mm++){
                msh.coord(ipnew,mm) = (1.0+ii)/ideg*msh.coord(msh.tet2poi(ielen,lnofa3[ifan][0]),mm)
                                    + (1.0+jj)/ideg*msh.coord(msh.tet2poi(ielen,lnofa3[ifan][1]),mm)
                                    + (1.0+kk)/ideg*msh.coord(msh.tet2poi(ielen,lnofa3[ifan][2]),mm);
              }
            }
          }

          // Additionally, add the new face to hash table. 
          facHsh[key] = {3, ielen};

        }// tt != facHsh.end()

      }// for ied


      // Create interior control points. These are all those of the form (1,1,1,1) + (i,j,k,l)
      for(int ii = 0; ii < ideg - 4; ii++){
        for(int jj = 0; jj < ideg - 4; jj++){
          for(int kk = 0; kk < ideg - 4; kk++){
            int ll = 1 - ii - jj - kk;
            int inod3 = mul2nod(1+ii,1+jj,1+kk,1+ll);
            int ipnew = msh.newpoitopo(3, ielen);
            msh.tet2poi(ielen,inod3) = ipnew;

            for(int mm = 0; mm < 3; mm++){
              msh.coord(ipnew,mm) = (1.0+ii)/ideg*msh.coord(msh.tet2poi(ielen,0),mm)
                                  + (1.0+jj)/ideg*msh.coord(msh.tet2poi(ielen,1),mm)
                                  + (1.0+kk)/ideg*msh.coord(msh.tet2poi(ielen,2),mm)
                                  + (1.0+ll)/ideg*msh.coord(msh.tet2poi(ielen,3),mm);
            }
          }
        }
      }// for ii 

      for(int ii = 0; ii < METRIS_MAXTAGS; ii++) msh.tet2tag(ii,ielen) = 0;


    }// for iver
  }// for iele0

  if(check_qua && qpnorm > 0){
    info.qcav3 = pow(info.qcav3, 1.0 / qpnorm);
    CPRINTF1(" - Final tetra cavity quality = {:.3f}\n",info.qcav3);
  }


  return 0;
}


// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template int reconnect_tetcav<MetricFieldAnalytical, n >(Mesh<MetricFieldAnalytical> &msh,\
 MshCavity& cav, CavOprOpt &opts, CavOprInfo &info, int nfac0, double *qmax, int ithread);\
template int reconnect_tetcav<MetricFieldFE        , n >(Mesh<MetricFieldFE        > &msh,\
 MshCavity& cav, CavOprOpt &opts, CavOprInfo &info, int nfac0, double *qmax, int ithread); 
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



}