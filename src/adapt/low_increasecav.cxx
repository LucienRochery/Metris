//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../MetrisRunner/MetrisParameters.hxx"
#include "../adapt/low_increasecav.hxx"
#include "../adapt/low_delaunay.hxx"
#include "../low_geo/normal.hxx"
#include "../low_geo/measure.hxx"
#include "../low_geo/lenedg.hxx"
#include "../aux_topo.hxx"
#include "../low_topo.hxx"
#include "../utils/mprintf.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../Mesh/Mesh.hxx"


//#define NODELSURF

namespace Metris{


// Check if any removed points; only those > 1/sqrt(2) from ipins if chklen
// This can possibly be reworked to be faster, for now we check everything every
// time, even though this is called in iterative cavity building.
template<class MFT>
void check_cavity_rempoint(MeshMetric<MFT> &msh, MshCavity &cav, const CavOprOpt &opts,
                           intAr1 &lrempoi, bool chklen, int ithrd1){
  GETVDEPTH(msh.param);

  lrempoi.set_n(0);

  for(int ientt : cav.lcedg) msh.edg2tag(ithrd1,ientt) = msh.tag[ithrd1];
  for(int ientt : cav.lcfac) msh.fac2tag(ithrd1,ientt) = msh.tag[ithrd1];
  for(int ientt : cav.lctet) msh.tet2tag(ithrd1,ientt) = msh.tag[ithrd1];

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
        if(ent2tag(ithrd1,ient2) < msh.tag[ithrd1]){
          msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1];
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
      if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
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
        // If not, this point won't be removed. Continue. 
        if(iface >= 0 && msh.fac2tag(ithrd1,iface) < msh.tag[ithrd1]) continue;
      } 

      if(tdimn >= 2){
        // If there is an edge attached, check it is in the cavity.
        int iedge = getpoiedg(msh,ipoin);
        // If not, this point won't be removed. Continue. 
        if(iedge >= 0 && msh.edg2tag(ithrd1,iedge) < msh.tag[ithrd1]) continue;
      }

      // If we're here, that means that there are either no attached lower dim
      // or there are and they are all in the cavity; indeed, assume there exist
      // at least one, and at least one not in the cav. Then the point is not
      // tagged. Then we wouldn't be here. 

      CPRINTF1(" ## point %d will be removed \n",ipoin);
      // tag point so we don't check for it again
      msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1];
      if(chklen){
        int edg2pol[2] = {cav.ipins, ipoin};
        double sz[2];
        double len = msh.idim == 2 ? getlenedg_geosz<MFT,2,1>(msh,edg2pol,sz)
                                   : getlenedg_geosz<MFT,3,1>(msh,edg2pol,sz);
        CPRINTF1(" -> found len = %e >? 1/sqrt(2): %d\n",len,len*sqrt(2) > 1);
        if(len > 1.0/sqrt(2)) lrempoi.stack(ipoin);
      }else{
        lrempoi.stack(ipoin);
      }
    }
  }

  return;
}

template void check_cavity_rempoint<MetricFieldAnalytical>
  (MeshMetric<MetricFieldAnalytical> &msh, MshCavity &cav, const CavOprOpt &opts,
   intAr1 &lrempoi, bool chklen, int ithrd1);
template void check_cavity_rempoint<MetricFieldFE        >
  (MeshMetric<MetricFieldFE        > &msh, MshCavity &cav, const CavOprOpt &opts,
   intAr1 &lrempoi, bool chklen, int ithrd1);


// Increase for validity and Delaunay (if idelaunay == true) both. 
template<class MFT>
int increase_cavity(MeshMetric<MFT> &msh, MshCavity &cav, 
                    bool idelaunay, int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);
  METRIS_ASSERT(ithrd1 != ithrd2);
  METRIS_ASSERT(cav.ipins >= 0 && cav.ipins < msh.npoin);


  //#ifdef NODELSURF
  //static int nwarn = 0;
  //// Disable surf
  //if(msh.get_tdim() < msh.idim && msh.param->iflag1 == 0 && idelaunay){
  //  if(nwarn++ < 10) MPRINTF("## WARNING DELAUNAY SURFACE DISABLED\n");
  //  idelaunay = false;
  //}
  //#endif

  msh.tag[ithrd1]++;
  if(idelaunay) msh.tag[ithrd2]++;

  // Tag entities and references
  for(int tdim = 1; tdim <= 3; tdim++){
    const intAr1& lcent = cav.lcent(tdim);
    for(int ientt : lcent){
      METRIS_ASSERT(ientt >= 0 && ientt < msh.nentt(tdim));
      METRIS_ASSERT(!isdeadent(ientt,msh.ent2poi(tdim)));
      msh.ent2tag(tdim)(ithrd1,ientt) = msh.tag[ithrd1];

      int iref = msh.ent2ref(tdim)[ientt];
      METRIS_ASSERT(iref >= 0);
      if(msh.ref2tag(tdim)(ithrd1,iref) < msh.tag[ithrd1]){
        CPRINTF1(" - ipins has edge ref %d \n",iref);
      }
      msh.ref2tag(tdim)(ithrd1,iref) = msh.tag[ithrd1];
    }
  }


  int ibins = msh.poi2bpo[cav.ipins];
  int pdim  = msh.get_tdim();
  if(ibins >= 0) pdim = msh.bpo2ibi(ibins,1);

  CPRINTF1("-- START increase_cavity ipins %d dim %d list initial cavity:\n", cav.ipins, pdim);
  if(DOPRINTS1()){
    if(cav.lcedg.get_n() > 0){
      CPRINTF1(" - Edge cavity: ");
      cav.lcedg.print();
    }
    if(cav.lcfac.get_n() > 0){
      CPRINTF1(" - Face cavity: ");
      cav.lcfac.print();
    }
    if(cav.lctet.get_n() > 0){
      CPRINTF1(" - Tetra cavity: ");
      cav.lctet.print();
    }
  }
  if(DOPRINTS2()){
    for(int tdim = 1; tdim <= 3; tdim++){
      if(cav.lcent(tdim).get_n() <= 0) continue;
      if(tdim == 1){
        CPRINTF2(" - Edge cavity: \n");
      }else if(tdim == 2){
        CPRINTF2(" - Face cavity: \n");
      }else{
        CPRINTF2(" - Tetra cavity: \n");
      }
      for(int ientt : cav.lcent(tdim)){
        CPRINTF2("%d : ",ientt);
        for(int ii = 0; ii < msh.nnode(tdim); ii++){
          printf(" %d ",msh.ent2poi(tdim)(ientt,ii));
        }
        printf("\n");
      }
    }
  }

  int ent2pol[4];
  ent2pol[0] = cav.ipins;

  bool restart;
  int niter = 0;
  int ient0[2] = {0,0};


  int nnmet = (msh.idim * (msh.idim + 1)) / 2;
  double metl[6], lmet[6];
  double *metl_p; 
  if(idelaunay){
    if(msh.met.getSpace() == MetSpace::Log){
      for(int ii = 0; ii < nnmet; ii++) lmet[ii] = msh.met(cav.ipins,ii);
      if(msh.idim == 2){
        getexpmet_cpy<2>(lmet, metl);
      }else{
        getexpmet_cpy<3>(lmet, metl);
      }
      metl_p = metl;
    }else{
      metl_p = msh.met[cav.ipins];
    }
  }

  do{

    restart = false;
    if(niter++ > 100){
      #ifndef NDEBUG
      MPRINTF("# Possibly infinite cavity ext iterations 100\n");
      printf("## WAIT\n");
      wait();
      #endif
      return 1;
    }

    CPRINTF1(" - inccav iter %d ifac0 %d itetr0 %d \n",niter,ient0[0],ient0[1]);


    int ient01[2] = {cav.lcfac.get_n(), cav.lctet.get_n()};

    for(int tdim = 2; tdim <= msh.get_tdim(); tdim++){

      intAr1 &lcent = cav.lcent(tdim);
      const intAr2 &ent2poi = msh.ent2poi(tdim);
      const intAr2 &ent2ent = msh.ent2ent(tdim);
      const intAr1 &ent2ref = msh.ent2ref(tdim);
      const intAr2 &ref2tag = msh.ref2tag(tdim);
      intAr2 &ent2tag = msh.ent2tag(tdim);
      const intAr2 &sub2tag = msh.ent2tag(tdim-1);


      // Note the bound is reeval'd, can't use range based
      for(int ientl = ient0[tdim-2]; ientl < lcent.get_n(); ientl++){
        INCVDEPTH(msh.param)
        int ientt = lcent[ientl];
        if(tdim == 2){
          CPRINTF1(" - inccav try %d / %d = %d (%d,%d,%d) \n",
                   ientl,lcent.get_n(),ientt,ent2poi(ientt,0),
                   ent2poi(ientt,1),ent2poi(ientt,2));
        }else{
          CPRINTF1(" - inccav try %d / %d = %d (%d,%d,%d,%d) \n",
                   ientl,lcent.get_n(),ientt,ent2poi(ientt,0),
                   ent2poi(ientt,1),ent2poi(ientt,2),ent2poi(ientt,3));
        }

        double norCAD[3];
        if(msh.idim == 3 && tdim == 2){
          // If dimension 3 topo dim 2, get a normal for this face. 
          if(msh.CAD()){
            getnorfacCAD(msh,ientt,norCAD);
          }else{
            getnorfacP1(ent2poi[ientt],msh.coord,norCAD);
          }
        }

        for(int inei = 0; inei < tdim + 1; inei++){

          bool iskip = false;
          for(int ii = 0; ii < tdim + 1; ii++){
            int ip = ent2poi(ientt, ii);
            if(ip != cav.ipins) continue;
            iskip = true;
            break;
          }
          if(iskip) continue;

          int ienei = ent2ent(ientt,inei);

          CPRINTF1("   - inei %d ienei = %d\n", inei, ienei);

          if(ienei >= 0){
            if(ent2tag(ithrd1,ienei) >= msh.tag[ithrd1]){
              CPRINTF1("   - ienei = %d is tagged %d >= %d\n",
                                 ienei,ent2tag(ithrd1,ienei),msh.tag[ithrd1]);
              continue;
            }
          }

          // tdim 2: if there's an edge here and it's in the cavity, then it will 
          // be split and we'll get no face from it. 
          // tdim 3: idem, faces. 
          int isube = -1;
          if(tdim == 2){
            isube = msh.fac2edg(ientt,inei);
          }else{
            isube = msh.tet2fac(ientt,inei);
          }
          if(isube >= 0 && sub2tag(ithrd1,isube) >= msh.tag[ithrd1]){
            CPRINTF1("   - ientt %d -> isube %d is tagged, skip\n",ientt,isube);
            continue;
          }

          // New face is ipins, ip1, ip2  
          if(tdim == 2){
            ent2pol[1] = ent2poi(ientt,lnoed2[inei][0]);
            ent2pol[2] = ent2poi(ientt,lnoed2[inei][1]);
          }else{
            ent2pol[lnofa3[0][0]] = ent2poi(ientt,lnofa3[inei][0]);
            ent2pol[lnofa3[0][1]] = ent2poi(ientt,lnofa3[inei][1]);
            ent2pol[lnofa3[0][2]] = ent2poi(ientt,lnofa3[inei][2]);
          }

          // First, check if this is a sliver
          bool iflat;
          double meas0;
          meas0 = msh.idim == 2 ? getmeasentP1<2,2>(msh, ent2pol, norCAD, &iflat)
                :     tdim == 2 ? getmeasentP1<3,2>(msh, ent2pol, norCAD, &iflat)
                                : getmeasentP1<3,3>(msh, ent2pol, norCAD, &iflat);
          
          CPRINTF1("   - inccav pdim %d tdim %d ent %d = ",pdim,tdim,ientt);
          if(DOPRINTS1()){
            if(tdim == 2){
              printf(" %d %d %d ",ent2pol[0],ent2pol[1],ent2pol[2]);
            }else{
              printf(" %d %d %d %d ",ent2pol[0],ent2pol[1],ent2pol[2],ent2pol[3]);
            }
          }
          CPRINTF1(" w/ vtol = %f got iflat = %d meas0 = %15.7e neighbour = %d\n",
                   msh.param->vtol,iflat,meas0,ienei);

          #if 0
          // Next check geodev 
          // Actually not because adding more faces will only damage the cavity further
          // Do this in the future as pre reject, possibly. 
          // Also depends on Pk etc. Probably best to leave in cav.
          double nrmal[3]; 
          if(msh.idim == 3 && pdim < 2){
            // Get the normal in the case we're on an edge in 3D, and get only 
            // the correct side.
            int iref = msh.fac2ref[ientt];
            getnorballref<1>(msh,lcent,iref,nrmal);
          }
          #endif

          // if element created with this facet is negative, add the neighbour
          // to cavity. 
          if(iflat || meas0 < 0){
            if(ienei >= 0){
              if(ref2tag(ithrd1,ent2ref[ienei]) < msh.tag[ithrd1]){
                CPRINTF1("   - ienei = %d is wrong ref %d -> cannot correct\n",
                         ienei,ent2ref[ienei]);
                return 1;
              }
              lcent.stack(ienei);
              ent2tag(ithrd1,ienei) = msh.tag[ithrd1];
              CPRINTF1("   - inccav added entt %d to stack \n", ienei);
              // If this is a face, we must also add the supported tets
              if(tdim == 2 && msh.nelem > 0){
                for(int ii = 0; ii < 2; ii++){
                  int ielem = msh.fac2tet(ienei, ii);
                  if(ielem < 0) continue;
                  if(msh.tet2tag(ithrd1,ielem) >= msh.tag[ithrd1]) continue;
                  int iref = msh.tet2ref[ielem];
                  if(msh.dom2tag(ithrd1,iref) < msh.tag[ithrd1]){
                    CPRINTF1("   - iface %d -> itetr %d is wrong ref %d\n",ienei,ielem,iref);
                    return 1;
                  }
                  cav.lctet.stack(ielem);
                  msh.tet2tag(ithrd1,ielem) = msh.tag[ithrd1];
                  CPRINTF1("   - inccav added tet %d to stack \n", ielem);
                }
              }
            }

            // There are two cases:
            // - ienei >= 0, then this entity is sandwiched and needs to be added
            // - ienei < 0, then the only hope of correction is adding this entity
            // Hence, in any case, if there is a subdim entity here, add it. 
            if(isube >= 0){
              // If the point tdim is greater than this element's dim, it cannot
              // be added.
              if(pdim > tdim-1){
                CPRINTF1("   - ientt %d dim %d < ipins dim %d\n",isube,tdim-1,pdim);
                return 2;
              }
              // Add the boundary entity, but only if in allowed refs. 
              int iref = msh.ent2ref(tdim-1)[isube];
              if(msh.ref2tag(tdim-1)(ithrd1,iref) < msh.tag[ithrd1]){
                CPRINTF1("   - ientt %d -> isube %d is wrong ref %d\n",ienei,isube,iref);
                return 1;
              }
              cav.lcent(tdim-1).stack(isube);
              msh.ent2tag(tdim-1)(ithrd1,isube) = msh.tag[ithrd1];
              CPRINTF1("   - inccav added dim %d ent %d to stack \n", tdim-1, 
                       isube);
              // We added a lower dim entity, hence restart will be required. 
              restart = true;
            }


            // If added due to validity, skip delaunay
            continue;
          }// if iflat

          // Only apply Delaunay on highest tdim elements.
          if(idelaunay && ienei >= 0 && tdim == msh.get_tdim()){
            if(ent2tag(ithrd2,ienei) >= msh.tag[ithrd2]){
              CPRINTF1("   - ienei = %d has already been checked for delaunay -> skip\n",
                ienei);
              continue;
            }
            ent2tag(ithrd2,ienei) = msh.tag[ithrd2];

            if(ref2tag(ithrd1,ent2ref[ienei]) < msh.tag[ithrd1]){
              CPRINTF1("   - ienei = %d is wrong ref %d -> skip Delaunay\n",
                       ienei,ent2ref[ienei]);
              continue;
            }


            // Check if Delaunay 
            bool isinsph;
            if(tdim == 2){
              if(msh.idim == 2){
                isinsph = indelsphere<2,2>(msh, msh.coord[cav.ipins], metl_p, 
                                           ent2poi[ienei]);
              }else{
                isinsph = indelsphere<3,2>(msh, msh.coord[cav.ipins], metl_p, 
                                           ent2poi[ienei]);
              }
            }else{
              isinsph = indelsphere<3,3>(msh, msh.coord[cav.ipins], metl_p, 
                                         ent2poi[ienei]);
            }

            if(isinsph){
              lcent.stack(ienei);
              ent2tag(ithrd1,ienei) = msh.tag[ithrd1];

              if(isube >= 0){
                // Add the boundary entity, but only if in allowed refs. 
                int iref = msh.ent2ref(tdim-1)[isube];
                if(msh.ref2tag(tdim-1)(ithrd1,iref) < msh.tag[ithrd1]){
                  CPRINTF1("   - ientt %d -> isube %d is wrong ref %d\n",ienei,isube,iref);
                  return 1;
                }
                cav.lcent(tdim-1).stack(isube);
                msh.ent2tag(tdim-1)(ithrd1,isube) = msh.tag[ithrd1];
                CPRINTF1("   - inccav added dim %d ent %d to stack \n", tdim-1, 
                         isube);
                // We added a lower dim entity, hence restart will be required. 
                restart = true;
              }
            }

          }// if idelaunay
      
        } // for int inei

      } // for int ientl
    } // for int tdim

    ient0[0] = ient01[0];
    ient0[1] = ient01[1];

  }while(restart);

  return 0;
}


template int increase_cavity(MeshMetric<MetricFieldAnalytical> &msh, MshCavity &cav, 
                    bool idelaunay, int ithrd1, int ithrd2);
template int increase_cavity(MeshMetric<MetricFieldFE        > &msh, MshCavity &cav, 
                    bool idelaunay, int ithrd1, int ithrd2);








// Increase for validity. Only allow same refs as ipins already has. 
int increase_cavity_validity(MeshBase &msh, MshCavity &cav, int ithread){
  GETVDEPTH(msh.param);

  METRIS_ASSERT(cav.ipins >= 0 && cav.ipins < msh.npoin);

  msh.tag[ithread]++;

  // Tag point's surface references if any. Filter entities
  aux_taginsrefs(msh,cav,ithread);

  for(int ielem : cav.lctet){
    METRIS_ASSERT(ielem >= 0 && ielem < msh.nelem);
    METRIS_ASSERT(!isdeadent(ielem,msh.tet2poi));
    msh.tet2tag(ithread,ielem) = msh.tag[ithread];
  }

  for(int iface : cav.lcfac){
    METRIS_ASSERT(iface >= 0 && iface < msh.nface);
    METRIS_ASSERT(!isdeadent(iface,msh.fac2poi));
    msh.fac2tag(ithread,iface) = msh.tag[ithread];
    if(!msh.isboundary_faces()) continue;

    int iref = msh.fac2ref[iface];
    METRIS_ASSERT(iref >= 0);
    METRIS_ASSERT(msh.cfa2tag(ithread,iref) <= msh.tag[ithread]);
    if(msh.cfa2tag(ithread,iref) < msh.tag[ithread]){
      CPRINTF1("## ERROR increase_cavity_validity: cavity face ref %d is not a ipins bdry ref\n",iref);
      return 2;
    }
  }

  for(int iedge : cav.lcedg){
    METRIS_ASSERT(iedge >= 0 && iedge < msh.nedge);
    METRIS_ASSERT(!isdeadent(iedge,msh.edg2poi));
    msh.edg2tag(ithread,iedge) = msh.tag[ithread];
    if(!msh.isboundary_edges()) continue;

    int iref = msh.edg2ref[iedge];
    METRIS_ASSERT(msh.ced2tag(ithread,iref) <= msh.tag[ithread]);
    if(msh.ced2tag(ithread,iref) < msh.tag[ithread]){
      CPRINTF1("## ERROR increase_cavity_validity: cavity edge is not a ipins bdry ref\n");
      return 2;
    }
  }

  CPRINTF1("-- START increase_cavity_validity ipins %d list initial cavity:\n", cav.ipins);
  if(DOPRINTS1()){
    if(cav.lcedg.get_n() > 0){
      CPRINTF1(" - Edge cavity: ");
      cav.lcedg.print();
    }
    if(cav.lcfac.get_n() > 0){
      CPRINTF1(" - Face cavity: ");
      cav.lcfac.print();
    }
    if(cav.lctet.get_n() > 0){
      CPRINTF1(" - Tetra cavity: ");
      cav.lctet.print();
    }
  }
  if(DOPRINTS2()){
    for(int tdim = 1; tdim <= 3; tdim++){
      intAr1 &lcent = cav.lcent(tdim);
      int ncent = lcent.get_n();
      if(ncent <= 0) continue;
      intAr2 &ent2poi = msh.ent2poi(tdim);

      if(tdim == 1){
        CPRINTF2(" - Edge cavity: \n");
      }else if(tdim == 2){
        CPRINTF2(" - Face cavity: \n");
      }else{
        CPRINTF2(" - Tetra cavity: \n");
      }
      int nnode = msh.nnode(tdim);
      for(int ientt : lcent){
        CPRINTF2("%d : ",ientt);
        for(int ii = 0; ii < nnode; ii++){
          printf(" %d ",ent2poi(ientt,ii));
        }
        printf("\n");
      }
    }
  }

  int ibins = msh.poi2bpo[cav.ipins];
  int pdim  = msh.get_tdim();
  if(ibins >= 0) pdim = msh.bpo2ibi(ibins,1);

  int ent2pol[4];
  ent2pol[0] = cav.ipins;

  bool restart;
  int niter = 0;
  int ient0[2] = {0,0};

  do{

    restart = false;
    if(niter++ > 100){
      #ifndef NDEBUG
      MPRINTF("# Possibly infinite cavity ext iterations 100\n");
      printf("## WAIT\n");
      wait();
      #endif
      return 1;
    }

    CPRINTF1(" - inccav iter %d ifac0 %d itetr0 %d \n",niter,
             ient0[0],ient0[1]);


    int ient01[2] = {cav.lcfac.get_n(), cav.lctet.get_n()};

    // Note the bound is reeval'd, can't use range based
    for(int tdim = 2; tdim <= 3; tdim++){

      intAr1 &lcent = cav.lcent(tdim);
      const intAr2 &ent2poi = msh.ent2poi(tdim);
      const intAr2 &ent2ent = msh.ent2ent(tdim);
      intAr2 &ent2tag = msh.ent2tag(tdim);


      for(int ientl = ient0[tdim-2]; ientl < lcent.get_n(); ientl++){
        INCVDEPTH(msh.param)
        int ientt = lcent[ientl];
        if(tdim == 2){
          CPRINTF1(" - inccav try %d / %d = %d (%d,%d,%d) \n",
                   ientl,lcent.get_n(),ientt,ent2poi(ientt,0),
                   ent2poi(ientt,1),ent2poi(ientt,2));
        }else{
          CPRINTF1(" - inccav try %d / %d = %d (%d,%d,%d,%d) \n",
                   ientl,lcent.get_n(),ientt,ent2poi(ientt,0),
                   ent2poi(ientt,1),ent2poi(ientt,2),ent2poi(ientt,3));
        }

        double norCAD[3];
        if(msh.idim == 3 && tdim == 2){
          // If dimension 3 topo dim 2, get a normal for this face. 
          if(msh.CAD()){
            getnorfacCAD(msh,ientt,norCAD);
          }else{
            getnorfacP1(ent2poi[ientt],msh.coord,norCAD);
          }
        }

        for(int inei = 0; inei < tdim + 1; inei++){

          bool iskip = false;
          for(int ii = 0; ii < tdim + 1; ii++){
            int ip = ent2poi(ientt, ii);
            if(ip != cav.ipins) continue;
            iskip = true;
            break;
          }
          if(iskip) continue;

          int ienei = ent2ent(ientt,inei);

          CPRINTF1("   - inei %d ienei = %d\n", inei, ienei);

          if(ienei >= 0){
            if(ent2tag(ithread,ienei) >= msh.tag[ithread]){
              CPRINTF1("   - ienei = %d is tagged %d >= %d\n",
                                 ienei,ent2tag(ithread,ienei),msh.tag[ithread]);
              continue;
            }
            if(tdim == 2){
              int iref = msh.fac2ref[ienei];
              if(msh.cfa2tag(ithread,iref) < msh.tag[ithread] && msh.isboundary_faces()){
                CPRINTF1("   - ienei = %d is wrong bdry ref %d\n",ienei,iref);
                continue;
              }
            }else{
              if(msh.tet2ref[ienei] != msh.tet2ref[ientt]){
                CPRINTF1("   - ienei %d ref = %d != ientt %d ref %d -> skip\n",
                ienei,msh.tet2ref[ienei],ientt,msh.tet2ref[ientt]);
                continue;
              }
            }
          }

          // tdim 2: if there's an edge here and it's in the cavity, then it will 
          // be split and we'll get no face from it. 
          // tdim 3: idem, faces. 
          int iedge = -1, iface = -1;
          if(tdim == 2){
            iedge = msh.fac2edg(ientt,inei);
            if(iedge >= 0){
              if(msh.edg2tag(ithread,iedge) >= msh.tag[ithread]){
                CPRINTF1("   - iface %d -> iedge %d is tagged, skip\n",ientt,iedge);
                continue;
              }
              //int iref = msh.edg2ref[iedge];
              //if(msh.ced2tag(ithread,iref) < msh.tag[ithread] && msh.isboundary_edges()){
              //  CPRINTF1("   - iface %d -> iedge %d is wrong bdry ref %d\n",ienei,iedge,iref);
              //  continue;
              //}
            }
          }else{
            iface = msh.tet2fac(ientt,inei);
            if(iface >= 0){
              if(msh.fac2tag(ithread,iface) >= msh.tag[ithread]){
                CPRINTF1("   - itetr %d -> iface %d is tagged, skip\n",ientt,iface);
                continue;
              }
              //int iref = msh.fac2ref[iface];
              //if(msh.cfa2tag(ithread,iref) < msh.tag[ithread]){
              //  CPRINTF1("   - itetr %d -> iface %d is wrong bdry ref %d\n",ientt,iface,iref);
              //  continue;
              //}
            }
          }

          // New face is ipins, ip1, ip2  
          if(tdim == 2){
            ent2pol[1] = ent2poi(ientt,lnoed2[inei][0]);
            ent2pol[2] = ent2poi(ientt,lnoed2[inei][1]);
          }else{
            ent2pol[lnofa3[0][0]] = ent2poi(ientt,lnofa3[inei][0]);
            ent2pol[lnofa3[0][1]] = ent2poi(ientt,lnofa3[inei][1]);
            ent2pol[lnofa3[0][2]] = ent2poi(ientt,lnofa3[inei][2]);
          }

          // First, check if this is a sliver
          bool iflat;
          double meas0;
          meas0 = msh.idim == 2 ? getmeasentP1<2,2>(msh, ent2pol, norCAD, &iflat)
                :     tdim == 2 ? getmeasentP1<3,2>(msh, ent2pol, norCAD, &iflat)
                                : getmeasentP1<3,3>(msh, ent2pol, norCAD, &iflat);
          
          CPRINTF1("   - inccav pdim %d tdim %d ent %d = ",pdim,tdim,ientt);
          if(DOPRINTS1()){
            if(tdim == 2){
              printf(" %d %d %d ",ent2pol[0],ent2pol[1],ent2pol[2]);
            }else{
              printf(" %d %d %d %d ",ent2pol[0],ent2pol[1],ent2pol[2],ent2pol[3]);
            }
          }
          CPRINTF1(" w/ vtol = %f got iflat = %d meas0 = %15.7e neighbour = %d\n",
                   msh.param->vtol,iflat,meas0,ienei);

          #if 0
          // Next check geodev 
          // Actually not because adding more faces will only damage the cavity further
          // Do this in the future as pre reject, possibly. 
          // Also depends on Pk etc. Probably best to leave in cav.
          double nrmal[3]; 
          if(msh.idim == 3 && pdim < 2){
            // Get the normal in the case we're on an edge in 3D, and get only 
            // the correct side.
            int iref = msh.fac2ref[ientt];
            getnorballref<1>(msh,lcent,iref,nrmal);
          }
          #endif
          // ignore ienei < 0 as it could be bdry -> edge remeshing
          if((iflat || meas0 < 0)){
            //if(ienei == -1) return 1;
            //// Cannot be corrected 
            //if(ienei < 0){
            //  METRIS_ASSERT(iedge >= 0 && tdim == 2 || iface >= 0 && tdim ==3);
            //  CPRINTF1(" # abort flat no neighbour: meas %23.15e\n", meas0);
            //  return 1;
            //}

            if(ienei >= 0){
              lcent.stack(ienei);
              ent2tag(ithread,ienei) = msh.tag[ithread];
              CPRINTF1("   - inccav added entt %d to stack \n", ienei);
            }else{
              // Add the boundary entity, but only if in allowed refs. 
              if(tdim == 2){
                int iref = msh.edg2ref[iedge];
                if(msh.ced2tag(ithread,iref) < msh.tag[ithread] && msh.isboundary_edges()){
                  CPRINTF1("   - iface %d -> iedge %d is wrong bdry ref %d\n",ienei,iedge,iref);
                  return 1;
                }
              }else{
                int iref = msh.fac2ref[iface];
                if(msh.cfa2tag(ithread,iref) < msh.tag[ithread]){
                  CPRINTF1("   - itetr %d -> iface %d is wrong bdry ref %d\n",ienei,iedge,iref);
                  return 1;
                }
              }
              restart = true;
            }

            // If a subdim entity was sandwiched here, we need to add it
            // Also true if no neighbour -> add bdry entity. 
            if((tdim == 2 && iedge >= 0) || (tdim == 3 && iface >= 0)){
              if(tdim == 2){
                cav.lcedg.stack(iedge);
                msh.edg2tag(ithread,iedge) = msh.tag[ithread];
              }else{
                cav.lcfac.stack(iface);
                msh.fac2tag(ithread,iface) = msh.tag[ithread];
              }
              CPRINTF1("   - inccav added dim %d ent %d to stack \n", tdim - 1, 
                       tdim == 2 ? iedge : iface);
            }

            // If this is a face, we must also add the supported tets
            if(tdim == 2 && msh.nelem > 0){
              for(int ii = 0; ii < 2; ii++){
                int ielem = msh.fac2tet(ientt, ii);
                if(ielem < 0) continue;
                if(msh.tet2tag(ithread,ielem) >= msh.tag[ithread]) continue;
                msh.tet2tag(ithread,ielem) = msh.tag[ithread];
                cav.lctet.stack(ielem);
                CPRINTF1("   - inccav added tet %d to stack \n", ielem);
              }
            }

          }
      
        } // for int inei

      } // for int ientl
    } // for int tdim

    ient0[0] = ient01[0];
    ient0[1] = ient01[1];

  }while(restart);

  return 0;
}




// Increase cavity for Delaunay criterion on ipoin 
// Normal only needed in 3D case if cavity has faces
template<class MFT>
int increase_cavity_Delaunay(MeshMetric<MFT> &msh, MshCavity &cav, 
                             int ngrow, int ithread){

  GETVDEPTH(msh.param);

  //#ifdef NODELSURF
  //static int nwarn = 0;

  //// Disable surf
  //if(msh.get_tdim() < msh.idim && msh.param->iflag1 == 0){
  //  if(nwarn++ < 10) MPRINTF("## WARNING DELAUNAY SURFACE DISABLED\n");
  //  return 0;
  //}
  //#endif


  //if(msh.get_tdim() == 3) 
  //  METRIS_THROW_MSG(TODOExcept(), "Unit test this for n = 3. Implement gettetfac instead of getfacedg");
  // Simply disable surface Delaunay for now 

  int nnmet = (msh.idim * (msh.idim + 1)) / 2;

  int tdim = cav.lctet.get_n() > 0 ? 3 : 2;

  msh.tag[ithread]++;

  // Tag point's surface references if any. Filter entities
  aux_taginsrefs(msh,cav,ithread);

  for(int ielem : cav.lctet){
    METRIS_ASSERT(ielem >= 0 && ielem < msh.nelem);
    METRIS_ASSERT(!isdeadent(ielem,msh.tet2poi));
    msh.tet2tag(ithread,ielem) = msh.tag[ithread];
  }

  for(int iface : cav.lcfac){
    METRIS_ASSERT(iface >= 0 && iface < msh.nface);
    METRIS_ASSERT(!isdeadent(iface,msh.fac2poi));
    msh.fac2tag(ithread,iface) = msh.tag[ithread];
  }

  for(int iedge : cav.lcedg){
    METRIS_ASSERT(iedge >= 0 && iedge < msh.nedge);
    METRIS_ASSERT(!isdeadent(iedge,msh.edg2poi));
    msh.edg2tag(ithread,iedge) = msh.tag[ithread];
  }

  // Actually do only the one dimension, keep this in for the future, maybe.
  //for(int tdim = msh.get_tdim(); tdim <= msh.get_tdim(); tdim++){
  intAr1 &lcent = cav.lcent(tdim);
  //if(lcent.get_n() == 0) continue;
  intAr1 &lcsub = cav.lcent(tdim-1);

  CPRINTF1(" - Delaunay dim %d\n",tdim);
  const intAr2&  ent2ent = msh.ent2ent(tdim);
  const intAr2&  ent2poi = msh.ent2poi(tdim);
        intAr2r& ent2tag = msh.ent2tag(tdim);
        intAr2r& sub2tag = msh.ent2tag(tdim-1);


  double metl[6], lmet[6];
  double *metl_p; 
  if(msh.met.getSpace() == MetSpace::Log){
    for(int jj = 0; jj < nnmet; jj++) lmet[jj] = msh.met(cav.ipins,jj);
    if(msh.idim == 2){
      getexpmet_cpy<2>(lmet, metl);
    }else{
      getexpmet_cpy<3>(lmet, metl);
    }
    metl_p = metl;
  }else{
    metl_p = msh.met[cav.ipins];
  }


  int icen0 = 0, icen1 = lcent.get_n();
  for(int igrow = 0; igrow < ngrow || ngrow < 0; igrow++){

    for(int icent = icen0; icent < icen1; icent++){
      int ientt = lcent[icent];
      for(int jj = 0; jj < tdim + 1; jj++){
        int ienei = ent2ent(ientt,jj);
        if(ienei < 0) continue; // Non manifold skip

        if(ent2tag(ithread,ienei) >= msh.tag[ithread]){
          CPRINTF1("   - ienei = %d is tagged %d >= %d\n",
                   ienei,ent2tag(ithread,ienei),msh.tag[ithread]);
          continue;
        }

        int isube = -1;
        if(tdim == 2){
          int iref2 = msh.fac2ref[ienei];
          if(msh.cfa2tag(ithread,iref2) < msh.tag[ithread] && msh.isboundary_faces()){
            CPRINTF1("   - ienei = %d is wrong bdry ref %d\n",ienei,iref2);
            continue;
          }
          int isube = msh.fac2edg(ientt,jj);
          if(isube >= 0){
            if(msh.edg2tag(ithread,isube) >= msh.tag[ithread]){
              CPRINTF1("   - iface %d -> iedge %d is tagged, skip\n",ientt,isube);
              continue;
            }
            int iref1 = msh.edg2ref[isube];
            if(msh.ced2tag(ithread,iref1) < msh.tag[ithread] && msh.isboundary_edges()){
              CPRINTF1("   - iface %d -> iedge %d is wrong bdry ref %d\n",ienei,isube,iref1);
              continue;
            }
          }
        }else{
          if(msh.tet2ref[ienei] != msh.tet2ref[ientt]){
            CPRINTF1("   - ienei %d ref = %d != ientt %d ref %d -> skip\n",
                     ienei,msh.tet2ref[ienei],ientt,msh.tet2ref[ientt]);
          }
          int isube = msh.tet2fac(ientt,jj);
          if(isube >= 0){
            if(msh.fac2tag(ithread,isube) >= msh.tag[ithread]){
              CPRINTF1("   - itetr %d -> iface %d is tagged, skip\n",ientt,isube);
              continue;
            }
            int iref1 = msh.fac2ref[isube];
            if(msh.ced2tag(ithread,iref1) < msh.tag[ithread]){
              CPRINTF1("   - itetr %d -> iface %d is wrong bdry ref %d\n",ienei,isube,iref1);
              continue;
            }
          }
        }

        ent2tag(ithread,ienei) = msh.tag[ithread];

        bool isinsph;
        if(tdim == 2){
          if(msh.idim == 2){
            isinsph = indelsphere<2,2>(msh, msh.coord[cav.ipins], metl_p, 
                                       ent2poi[ienei]);
          }else{
            isinsph = indelsphere<3,2>(msh, msh.coord[cav.ipins], metl_p, 
                                       ent2poi[ienei]);
          }
        }else{
          isinsph = indelsphere<3,3>(msh, msh.coord[cav.ipins], metl_p, 
                                     ent2poi[ienei]);
        }
        if(isinsph){
          lcent.stack(ienei);
          if(isube >= 0){
            sub2tag(ithread,isube) = msh.tag[ithread];
            lcsub.stack(isube);
          }
        }
        
      }// for j = 0,tdim
    }// for icent


    icen0 = icen1;
    icen1 = lcent.get_n();
    CPRINTF1(" - del grow %d / %d + %d ent\n",igrow,ngrow,icen1-icen0);
    if(icen1 == icen0) break;
  }// for igrow
  //}// for tdim

  return 0;
}

template int increase_cavity_Delaunay(MeshMetric<MetricFieldAnalytical> &msh, 
                                      MshCavity &cav, int ngrow, int ithread);
template int increase_cavity_Delaunay(MeshMetric<MetricFieldFE        > &msh, 
                                      MshCavity &cav, int ngrow, int ithread);





template<class MFT>
int increase_cavity_lenedg(MeshMetric<MFT> &msh, MshCavity &cav, 
                           CavOprOpt &opts,
                           int ipins,int ithrd1, int ithrd2){
  int nprem = 0;
//  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    CT_FOR0_INC(2,3,gdim){if(msh.idim == gdim){
      nprem = increase_cavity_lenedg0<MFT,gdim>(msh,cav,opts,ipins,ithrd1,ithrd2);
    }}CT_FOR1(gdim);
//  }}CT_FOR1(ideg);
  return nprem;
}

template<class MFT, int gdim>
int increase_cavity_lenedg0(MeshMetric<MFT> &msh, MshCavity &cav, 
                            CavOprOpt &opts,
                            int ipins, int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);


  // Note ipins must be seeded with newbpotopo
  const int pdim_ipins = msh.getpoitdim(ipins);

  //const intAr2 &ent2ent = msh.ent2ent(tdim);
  msh.tag[ithrd1]++;
  for(int tdim = 1; tdim <= msh.get_tdim(); tdim++){
    for(int ientt : cav.lcent(tdim)){
      msh.ent2tag(tdim)(ithrd1,ientt) = msh.tag[ithrd1];
    }
  }

  // Tag point's surface references if any. Filter entities
  aux_taginsrefs(msh,cav,ithrd1);



  intAr1 lbtet(20), lbfac(20), lbedg(20);
  int iopen;
  bool imani;

  int nprem = 0;

  int edg2pol[2];
  edg2pol[0] = ipins;
  double sz[2];

  //int ncomp = 0;
  //int ncav0 = lcent.get_n();

  // NB: loop bounds MUST be reevaluated ! don't range-for this 
  int cdim = 0;
       if(cav.lctet.get_n() > 0) cdim = 3;
  else if(cav.lcfac.get_n() > 0) cdim = 2;
  else if(cav.lcedg.get_n() > 0) cdim = 1;
  const int nedgl = (cdim*(cdim+1))/2;
  const intAr2 lnoed(nedgl,2,cdim == 2 ? lnoed2[0] : lnoed3[0]);

  intAr1 &lcent = cav.lcent(cdim); 
  for(int ii = 0; ii < lcent.get_n(); ii++){
    INCVDEPTH(msh.param);
    int ientt = lcent[ii];
    METRIS_ASSERT_MSG(!isdeadent(ientt, msh.ent2poi(cdim)),
      "entity "<<ientt<<" tdim "<<cdim<<" is dead");


    #if 0
    for(int ifa = 0; ifa < tdim + 1; ifa++){
      int ientn = ent2ent(ientt,ifa);
      if(ientn >= 0){
        if(ent2tag(ithrd1,ientn) >= msh.tag[ithrd1]) continue;
      }
      // Cavity boundary 
      // Loop over face nodes 
      int kk = -1;
      for(int ii = 0; ii < tdim; ii++){
        // Increment and skip when == to ifa (= not on facet)
        kk += 1 + ((kk + 1) == ifa);
        int ipoin = ent2poi(ientt,kk);
        if(ipoin == ipins) continue;
        if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
        msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1];

        edg2pol[1] = ipoin;
        double len = getlenedg_geosz<MFT,gdim,1>(msh,edg2pol,sz);
        ncomp++;
        if(len > 1.0/sqrt(2)) continue;


        // Short edge

        if(!opts.allow_remove_points) return -1; 
        if constexpr (tdim == 2){
          ball2(msh,ipoin,ientt,lbfac,dum,&iopen,&imani,ithrd2);
        }else{
          ball3(msh,ipoin,ientt,lbfac,&iopen,ithrd2);
        }
        nprem++;
        for(int ient2 : lbfac){
          if(ent2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          ent2tag(ithrd1,ient2) = msh.tag[ithrd1];
          lcent.stack(ient2);
        }
      }
    }
    #else
    for(int inode = 0; inode < cdim + 1; inode++){
      int ipoin = msh.ent2poi(cdim)(ientt,inode);
      if(ipoin == ipins) continue;
      if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
      msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1];

      edg2pol[1] = ipoin;
      double len = getlenedg_geosz<MFT,gdim,1>(msh,edg2pol,sz);


      CPRINTF1(" - check len ipoin %d len = %f <? 1/sqrt(2) %d\n",
                ipoin,len,len <= 1.0/sqrt(2));

      if(len <= 1.0/sqrt(2)){
        int pdim = msh.getpoitdim(ipoin);

        if(pdim < pdim_ipins){
          CPRINTF1(" - short edge and other end has dim %d < %d = dim ipins -> reject\n",
            pdim, pdim_ipins);
          return -1;
        }

        if(pdim == pdim_ipins && !opts.allow_remove_points){
          CPRINTF1(" - short edge and other end has dim %d = %d = dim ipins "
                  "w/ opts.allow_remove_points == false -> reject\n",
                 pdim, pdim_ipins);
          return -1;
        }

        if(pdim > pdim_ipins && !opts.allow_remove_points_superdim){
          CPRINTF1(" - short edge and other end has dim %d > %d = dim ipins "
                  "w/ opts.allow_remove_points_superdim == false -> reject\n",
                 pdim, pdim_ipins);
          return -1;
        }

        lbedg.set_n(0);
        lbfac.set_n(0);
        lbtet.set_n(0);
        if(cdim == 2){
          ball2(msh,ipoin,ientt,lbfac,lbedg,&iopen,&imani,ithrd2);
        }else{
          ball3(msh,ipoin,ientt,lbtet,&iopen,ithrd2);
          if(pdim <= 2){
            // Also get ball2 of point
            int iface = -1;
            if(pdim == 1){
              int iedge = msh.poi2ent(ipoin,0);
              iface = msh.edg2fac[iedge];
            }else{
              iface = msh.poi2ent(ipoin,0);
            }
            METRIS_ASSERT(iface >= 0 && iface < msh.nface);
            ball2(msh,ipoin,iface,lbfac,lbedg,&iopen,&imani,ithrd2);
          }
        }
        int ncel0 = cav.lctet.get_n();
        int ncfa0 = cav.lcfac.get_n();
        int nced0 = cav.lcedg.get_n();

        bool ifail = false;
        for(int ient2 : lbedg){
          if(msh.edg2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          int iref = msh.edg2ref[ient2];
          if(msh.ced2tag(ithrd1,iref) < msh.tag[ithrd1]){
            ifail = true;
            goto failed;
          }
        }
        for(int ient2 : lbfac){
          if(msh.fac2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          int iref = msh.fac2ref[ient2];
          if(msh.cfa2tag(ithrd1,iref) < msh.tag[ithrd1]){
            ifail = true;
            goto failed;
          }
        }
        for(int ient2 : lbtet){
          if(msh.tet2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          int iref = msh.tet2ref[ient2];
          if(msh.dom2tag(ithrd1,iref) < msh.tag[ithrd1]){
            ifail = true;
            goto failed;
          }
        }

        for(int ient2 : lbedg){
          if(msh.edg2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          msh.edg2tag(ithrd1,ient2) = msh.tag[ithrd1];
          cav.lcedg.stack(ient2);
        }
        for(int ient2 : lbfac){
          if(msh.fac2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          msh.fac2tag(ithrd1,ient2) = msh.tag[ithrd1];
          cav.lcfac.stack(ient2);
        }
        for(int ient2 : lbtet){
          if(msh.tet2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          msh.tet2tag(ithrd1,ient2) = msh.tag[ithrd1];
          cav.lctet.stack(ient2);
        }

        nprem++;

        failed:
        if(ifail){
          CPRINTF1(" - Failed to add point %d to collapse\n",ipoin);
          cav.lcedg.set_n(nced0);
          cav.lcfac.set_n(ncfa0);
          cav.lctet.set_n(ncel0);
        }
      }
    }
    #endif

    //// Control height, only in dimension 2d.
    //if(tdim == 2){

    //}else{
    //  METRIS_THROW_MSG(TODOExcept(), 
    //    "Implement height control in increase_cavity_lenedg 3D");
    //}
  }

  //printf("Debug ncavity init = %d final = %d ncomp = %d \n",ncav0,lcent.get_n(),ncomp);

  return nprem;
}

template int increase_cavity_lenedg(MeshMetric<MetricFieldAnalytical> &msh, 
           MshCavity &cav, CavOprOpt &opts, int ipins, int ithrd1, int ithrd2);
template int increase_cavity_lenedg(MeshMetric<MetricFieldFE        > &msh, 
           MshCavity &cav, CavOprOpt &opts, int ipins, int ithrd1, int ithrd2);


template int increase_cavity_lenedg0<MetricFieldAnalytical,2>(
                            MeshMetric<MetricFieldAnalytical> &msh, 
           MshCavity &cav, CavOprOpt &opts, int ipins, int ithrd1, int ithrd2);
template int increase_cavity_lenedg0<MetricFieldFE        ,2>(
                            MeshMetric<MetricFieldFE        > &msh, 
           MshCavity &cav, CavOprOpt &opts, int ipins, int ithrd1, int ithrd2);
template int increase_cavity_lenedg0<MetricFieldAnalytical,3>(
                            MeshMetric<MetricFieldAnalytical> &msh, 
           MshCavity &cav, CavOprOpt &opts, int ipins, int ithrd1, int ithrd2);
template int increase_cavity_lenedg0<MetricFieldFE        ,3>(
                            MeshMetric<MetricFieldFE        > &msh, 
           MshCavity &cav, CavOprOpt &opts, int ipins, int ithrd1, int ithrd2);





void aux_taginsrefs(MeshBase &msh, MshCavity &cav, int ithread){
  GETVDEPTH(msh.param);
  for(int iedge : cav.lcedg){
    int iref = msh.edg2ref[iedge];
    METRIS_ASSERT(iref >= 0);
    if(msh.ced2tag(ithread,iref) < msh.tag[ithread]){
      CPRINTF1(" - ipins has edge ref %d \n",iref);
    }
    msh.ced2tag(ithread,iref) = msh.tag[ithread];
  }
  for(int iface : cav.lcfac){
    int iref = msh.fac2ref[iface];
    METRIS_ASSERT(iref >= 0);
    if(msh.ced2tag(ithread,iref) < msh.tag[ithread]){
      CPRINTF1(" - ipins has face ref %d \n",iref);
    }
    msh.cfa2tag(ithread,iref) = msh.tag[ithread];
  }
  for(int ielem : cav.lctet){
    int iref = msh.tet2ref[ielem];
    METRIS_ASSERT_MSG(iref >= 0, "ielem = "<<ielem<<" invalid iref = "<<iref);
    if(msh.dom2tag(ithread,iref) < msh.tag[ithread]){
      CPRINTF1(" - ipins has tetra ref %d \n",iref);
    }
    msh.dom2tag(ithread,iref) = msh.tag[ithread];
  }
  #if 0
  for(int ibpoi = msh.poi2bpo[cav.ipins]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
    int bdim = msh.bpo2ibi(ibpoi,1);
    if(bdim == 0) continue;
    int ientt = msh.bpo2ibi(ibpoi,2);
    if(bdim == 1){
      int iref = msh.edg2ref[ientt];
      METRIS_ASSERT(iref >= 0);
      if(msh.ced2tag(ithread,iref) < msh.tag[ithread]){
        CPRINTF1(" - ipins has edge ref %d \n",iref);
      }
      msh.ced2tag(ithread,iref) = msh.tag[ithread];
    }else{
      int iref = msh.fac2ref[ientt];
      METRIS_ASSERT(iref >= 0);
      if(msh.cfa2tag(ithread,iref) < msh.tag[ithread]){
        CPRINTF1(" - ipins has face ref %d \n",iref);
      }
      msh.cfa2tag(ithread,iref) = msh.tag[ithread];
    }
  }
  #endif
}



} // end namespace
