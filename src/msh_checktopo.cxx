//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"

#include "msh_checktopo.hxx"
#include "aux_exceptions.hxx"
#include "aux_topo.hxx"
#include "low_geo.hxx"
#include "aux_utils.hxx"
#include "CT_loop.hxx"
#include "io_libmeshb.hxx"
#include "low_ccoef.hxx"
#include "mprintf.hxx"

namespace Metris{



void check_topo(MeshBase &msh, int ithread){
  check_topo(msh,msh.nbpoi,msh.npoin,msh.nedge,msh.nface,msh.nelem,ithread);
}


void check_topo(MeshBase &msh, 
                int nbpoi, int npoin, int nedge, int nface, int nelem, int ithread){
  GETVDEPTH(msh);

  try{

    if(DOPRINTS2()) printf("-- check_topo start \n");
    // Absolute difference CAD -> point
    const double geotol = 1.0e+6;

    int tdim = msh.idim;
    const int jdeg = tdim * (msh.curdeg - 1);
    const int ncoef = tdim == 2 ? facnpps[jdeg]
                                : tetnpps[jdeg];
    dblAr1 ccoef(ncoef);
    ccoef.set_n(ncoef);

    if(tdim > msh.get_tdim()) goto aftervalidity;


    CT_FOR0_INC(2,3,idim){if(idim == tdim){
      int nentt = msh.nentt(tdim);
      intAr2 &ent2poi = msh.ent2poi(tdim);
      for(int ientt = 0; ientt < nentt; ientt++){
        if(isdeadent(ientt,ent2poi)) continue;
        bool iflat;
        double meas = getmeasentP1<idim,idim>(msh,ent2poi[ientt],NULL,&iflat);
        if(iflat || meas <= 0){
          printf("## FLAT ELEMENT %d \n",ientt);
          writeMesh("flat"+std::to_string(ientt),msh);
          METRIS_THROW(GeomExcept());
        }
        if(msh.curdeg > 1){
          CT_FOR0_INC(2,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
            getsclccoef<idim,idim,ideg>(msh,ientt,NULL,&ccoef[0],&iflat);
          }}CT_FOR1(ideg);
          if(iflat){
            printf("## NEGATIVE JACOBIAN %d \n ",ientt);
            writeMesh("inva"+std::to_string(ientt),msh);
            METRIS_THROW(GeomExcept());
          }
        }

      }
    }}CT_FOR1(idim);

    aftervalidity:

    // Check poi2ent 
    for(int tdim = 1; tdim <= msh.get_tdim(); tdim++){
      int nentt = msh.nentt(tdim);
      const intAr2& ent2poi = msh.ent2poi(tdim);
      int nnode = msh.nnode(tdim);
      // Check all points owned by an element do in fact have a poi2ent of this 
      // dimension
      for(int ientt = 0; ientt < nentt; ientt++){
        if(isdeadent(ientt,ent2poi)) continue;
        for(int ii = 0; ii < nnode; ii++){
          int ipoin = ent2poi(ientt,ii);
          int ient2 = msh.poi2ent(ipoin,0);
          METRIS_ENFORCE(ient2 >= 0 && ient2 < nentt);
          int tdim2 = msh.poi2ent(ipoin,1);
          METRIS_ENFORCE(tdim2 <= tdim);
          int ibpoi = msh.poi2bpo[ipoin];
          if(ibpoi < 0) continue;
          METRIS_ENFORCE(tdim2 < msh.idim);
        }
      }
    }
    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      int ientt = msh.poi2ent(ipoin,0);
      if(ientt < 0) continue;
      int tdim2 = msh.poi2ent(ipoin,1);
      METRIS_ENFORCE(0 <= tdim2 && tdim2 <= msh.get_tdim());
      if(tdim2 >= 1){
        int iver = msh.getverent(ientt,tdim2,ipoin);
        METRIS_ENFORCE_MSG(iver >= 0, "ipoin = "<<ipoin<<" poi2ent "<<
          ientt<<" tdim "<< tdim2<<" not found");
      }else{
        int ibpoi = msh.poi2bpo[ipoin];
        METRIS_ENFORCE(ibpoi >= 0);
        int itype = msh.bpo2ibi(ibpoi,1);
        METRIS_ENFORCE(itype == 0);
        METRIS_ENFORCE(ientt == ipoin);
      }
    }


    if(msh.meshClass() == MeshClass::Mesh && msh.metricClass() == MetricClass::MetricFieldFE){
      const intAr2 &poi2bak = msh.metricClass() == MetricClass::MetricFieldFE ? 
       ((Mesh<MetricFieldFE> *)(&msh))->poi2bak 
      :((Mesh<MetricFieldAnalytical> *)(&msh))->poi2bak;

      for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
        if(msh.poi2ent(ipoin,msh.get_tdim()-1) < 0) continue;

        int pdim = msh.getpoitdim(ipoin);
        if(pdim < 0 || pdim > msh.get_tdim()){
          printf(" ## INVLIAD pdim %d \n",pdim);
          printf("iopin = %d \n",ipoin);
          writeMesh("debug_checktopo",msh);
          METRIS_THROW(TopoExcept());
        }
        if(pdim == 0) continue;
        int iebak = poi2bak(ipoin,pdim-1);
        if(iebak < 0){
          printf("ipoin %d pdim %d poi2bak %d \n",ipoin,pdim,iebak);
          printf("dump all:\n");
          for(int tdim = 1; tdim <= msh.get_tdim(); tdim++){
            printf("dim %d : %d \n",tdim, poi2bak(ipoin,tdim-1));
          }
        }
        METRIS_ENFORCE(iebak >= 0);
        const intAr2& ent2pob = msh.metricClass() == MetricClass::MetricFieldFE ?  
        ((Mesh<MetricFieldFE> *)(&msh))->bak->ent2poi(pdim):
        ((Mesh<MetricFieldAnalytical> *)(&msh))->bak->ent2poi(pdim);
        if(isdeadent(iebak,ent2pob)){
          printf("Point %d tdim %d has dead back seed %d \n",ipoin,pdim,iebak);
        }
        METRIS_ENFORCE_MSG(!isdeadent(iebak,ent2pob),
          "Point "<<ipoin<<" tdim "<<pdim<<" has back seed "<<iebak<<" which is dead");

        //METRIS_ENFORCE(iebak < msh_->bak->nentt(pdim));
      }
    }


    // Check tags correct 
    for(int itag = 0; itag < METRIS_MAXTAGS; itag++){
      for(int ielem = 0; ielem < nelem; ielem++){
        if(isdeadent(ielem,msh.tet2poi)) continue;
        METRIS_ENFORCE(msh.tet2tag(itag,ielem) <= msh.tag[itag] && msh.tet2tag(itag,ielem) >= 0);
      }
      for(int iface = 0; iface < nface; iface++){
        if(isdeadent(iface,msh.fac2poi)) continue;
        METRIS_ENFORCE_MSG(msh.fac2tag(itag,iface) <= msh.tag[itag] && msh.fac2tag(itag,iface) >= 0,
          "failed itag = "<<itag<<" iface = "<<iface<<" face tag = "<<msh.fac2tag(itag,iface)
          <<" tag bound = "<<msh.tag[itag]);
      }
      for(int iedge = 0; iedge < nedge; iedge++){
        if(isdeadent(iedge,msh.edg2poi)) continue;
        METRIS_ENFORCE_MSG(msh.edg2tag(itag,iedge) <= msh.tag[itag] && msh.edg2tag(itag,iedge) >= 0,"failed itag = "<<itag<<" iedge = "<<iedge
          <<" edgtag = "<<msh.edg2tag(itag,iedge)<<" tag = "<<msh.tag[itag]
          <<" = "<<msh.edg2poi(iedge,0)<<" "<<msh.edg2poi(iedge,1));
      }
      for(int ipoin = 0; ipoin < npoin; ipoin++){
        METRIS_ENFORCE(msh.poi2tag(itag,ipoin) <= msh.tag[itag] && msh.poi2tag(itag,ipoin) >= 0);
      }
    }

    for(auto t : msh.edgHshTab){
      std::tuple<int,int> k = t.first;
      int iedge = t.second;
      METRIS_ENFORCE(!isdeadent(iedge,msh.edg2poi));

      int ip1 = std::get<0>(k);
      int ip2 = std::get<1>(k);

      int jp1 = msh.edg2poi(iedge,0);
      int jp2 = msh.edg2poi(iedge,1);
      METRIS_ENFORCE_MSG(ip1 == jp1 && ip2 == jp2 ||
                         ip1 == jp2 && ip2 == jp1 ,
      "global hashtable has ip1 = "<<ip1<< " ip2 = "<<ip2<<
      " and iedge = "<<iedge<<" but iedge has nodes jp1 = "<<jp1<<" jp2 = "<<jp2<<"\n");

    }


    // check edg2fac
    for(int iedge = 0; iedge < nedge; iedge++){
      if(isdeadent(iedge,msh.edg2poi)) continue;
      int iface = msh.edg2fac[iedge];
      // No detached edges, and in bounds 
      METRIS_ENFORCE(msh.get_tdim() == 1 || (iface >= 0 && iface < msh.nface));
      int ip1 = msh.edg2poi(iedge,0);
      int ip2 = msh.edg2poi(iedge,1);

      int iedl; 
      try{
        iedl = getedgfac(msh,iface,ip1,ip2);
        METRIS_ENFORCE(iedl >= 0);
      }catch(const MetrisExcept &e){
        printf("edge not in triangle (edg2fac link)\n");
        throw(e);
      }
    }

    // check getedgglo
    for(int iedge = 0; iedge < nedge; iedge++){
      if(isdeadent(iedge,msh.edg2poi)) continue;
      int ip1 = msh.edg2poi(iedge,0);
      int ip2 = msh.edg2poi(iedge,1);

      int iedg2 = getedgglo(msh, ip1, ip2);
      METRIS_ENFORCE(iedg2 == iedge); 
    }



    for(int ipoin = 0; ipoin < npoin; ipoin++){
      int ientt = msh.poi2ent(ipoin,0);
      if(ientt < 0) continue;

      METRIS_ENFORCE_MSG(msh.poi2ent(ipoin,1) >= 1 
                     &&  msh.poi2ent(ipoin,1) <= 3 , "Wrong poi2ent ipoin = "<<ipoin
                 <<" value at 1 = "<<msh.poi2ent(ipoin,1)); 

      int ibpoi = msh.poi2bpo[ipoin];
      int tdimn = msh.get_tdim();
      if(ibpoi >= 0) tdimn = msh.bpo2ibi(ibpoi,1);

      METRIS_ENFORCE_MSG(tdimn == msh.poi2ent(ipoin,1) || tdimn == 0,
        "Wrong tdimn? ipoin = "<<ipoin<<" tdimn using bpo = "<<tdimn<<
        " in poi2ent = "<<msh.poi2ent(ipoin,1)); 

      // Corners -> edges
      if(tdimn == 0) tdimn = 1;

      intAr2 &ent2poi = msh.ent2poi(tdimn);
      int nentt = msh.nentt(tdimn);
      if(ientt >= nentt){
        printf("poi2ent link defective? ipoin = %d ibpoi = %d ientt = %d  tdimn = %d \n",ipoin,ibpoi,ientt,tdimn);
        printf("nedge = %d nface = %d nelem = %d \n",msh.nedge,msh.nface,msh.nelem);

        printf("face nodes: \n");
        int nnod2 = facnpps[msh.curdeg];
        intAr1(nnod2,msh.fac2poi[ientt]).print();

        int ibpo2 = ibpoi;
        while(ibpo2 >= 0){
          printf("ibpoi %d : ",ibpo2);
          intAr1(nibi,msh.bpo2ibi[ibpo2]).print();
          ibpo2 = msh.bpo2ibi(ibpo2,3);
        }

        if(ibpoi>=0){
          printf("bpo2ibi = ");
          intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
        }
        printf("ientt = %d nodes = ",ientt);
        int nnode = msh.nnode(tdimn);
        intAr1(nnode,ent2poi[ientt]).print();
        for(int tdim2 = 1; tdim2 <= msh.get_tdim(); tdim2++){
          int nnod2 = msh.nnode(tdim2);
          intAr2 &ent2po2 = msh.ent2poi(tdim2);
          printf("idim %d entt nodes ",tdim2);
          intAr1(nnod2,ent2po2[ientt]).print();
        }
      }
      METRIS_ENFORCE(ientt < nentt);

      bool ifnd = false;
      int nnode = msh.nnode(tdimn);
      for(int ii = 0; ii < nnode; ii++){
        if(ent2poi(ientt,ii) == ipoin){
          ifnd = true;
          break;
        }
      }
      if(!ifnd){
        printf("2 poi2ent link defective? ipoin = %d ibpoi = %d tdimn = %d \n",
               ipoin,ibpoi,tdimn);
        if(ibpoi>=0){
          printf("bpo2ibi = ");
          intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
        }
        printf("ientt = %d nodes = ",ientt);
        intAr1(nnode,ent2poi[ientt]).print();
        for(int tdim2 = 1; tdim2 <= msh.get_tdim(); tdim2++){
          int nnod2 = msh.nnode(tdim2);
          intAr2 &ent2po2 = msh.ent2poi(tdim2);
          printf("idim %d entt nodes ",tdim2);
          intAr1(nnod2,ent2po2[ientt]).print();
        }
      }
      METRIS_ENFORCE_MSG(ifnd, "poi2ent not fnd, ipoin = "<<ipoin
                         <<" poi2ent = "<<ientt);
    }


    for(int ielem = 0; ielem < nelem; ielem++){
      if(isdeadent(ielem,msh.tet2poi)) continue;
      for(int ifa = 0; ifa < 4 ; ifa++){
        int i1 = msh.tet2poi(ielem,lnofa3[ifa][0]);
        int i2 = msh.tet2poi(ielem,lnofa3[ifa][1]);
        int i3 = msh.tet2poi(ielem,lnofa3[ifa][2]);

        METRIS_ENFORCE_MSG(i1 >= 0 && i1 < npoin && i2 >= 0 && i2 < npoin && i3 >= 0 && i3 < npoin,
         "Face vertices within bounds 0 <= i < "<<npoin<<"observed: "<<i1<<" "<<i2<<" "<<i3<<"\n");

        int iele2 = msh.tet2tet(ielem,ifa); 

        METRIS_ENFORCE_MSG(iele2 >= -1 && iele2 < nelem,"Neighbour either inexistent or valid index");

        if(iele2 == -1 || msh.tet2ref[iele2] != msh.tet2ref[ielem]){
          int iface = msh.tetfac2glo(ielem,ifa);
          METRIS_ENFORCE_MSG(iface >= 0 && iface < nface,"No neighbour means neighbour is a bface");
          int j1 = msh.fac2poi(iface,0);
          int j2 = msh.fac2poi(iface,1);
          int j3 = msh.fac2poi(iface,2);
          int ifal = getfactet(msh,ielem,j1,j2,j3);
          METRIS_ENFORCE_MSG(ifal == ifa,"Face neighbour has same vertices as local face ");
          if(iele2 == -1){
            METRIS_ENFORCE_MSG(msh.fac2tet(iface,0) == ielem,
              "ielem = "<<ielem<<" fac2tet = "<<msh.fac2tet(iface,0)<<","<<msh.fac2tet(iface,1));
          }else{
            METRIS_ENFORCE_MSG(msh.fac2tet(iface,0) == ielem || msh.fac2tet(iface,1) == ielem, 
              "ielem = "<<ielem<<" fac2tet = "<<msh.fac2tet(iface,0)<<","<<msh.fac2tet(iface,1));
          }
          continue;
        }

        int ifa2 = getfactetOpp(msh,iele2,i1,i2,i3);
        METRIS_ENFORCE_MSG(ifa2 >= 0 && ifa2 < 4,"Returned face within 0-3 range");
        METRIS_ENFORCE_MSG(msh.tet2tet(iele2,ifa2) == ielem,"Neighbour has ielem as neighbour at common face");
      }
    }

    HshTabInt2 hshTab2(1.5*msh.nface + msh.nedge);

    msh.tag[ithread]++;
    msh.ced2tag.fill(0);
    for(int iface = 0; iface < nface; iface++){
      if(isdeadent(iface,msh.fac2poi)) continue;

      // Check not all three vertices on same ref edge. 
      if(msh.ced2tag.get_n() > 0){
        // This should only fail if routine called too early (before CAD init)

        int iedg[3] = {0,0,0};
        for(int ii = 0; ii < 3; ii++){
          int ipoin = msh.fac2poi(iface,ii);

          // Corners can tag an edge twice
          for(int jj = 0; jj < msh.ced2tag.get_stride(); jj++) msh.ced2tag(1,jj) = 0;

          for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
            int bdim = msh.bpo2ibi(ibpoi,1);
            if(bdim != 1) continue;
            int ientt = msh.bpo2ibi(ibpoi,2);
            METRIS_ASSERT(ientt >= 0);
            int iref = msh.edg2ref[ientt];
            METRIS_ASSERT(iref >= 0);
            if(msh.ced2tag(1,iref) == 1) continue;
            msh.ced2tag(1,iref) = 1;
            msh.ced2tag(0,iref)++;
            iedg[ii] = 1;
          }
        }
        if(iedg[0] + iedg[1] + iedg[2] == 3){
          for(int jj = 0; jj < msh.ced2tag.get_stride(); jj++){
            if(msh.ced2tag(0,jj) == 3){
              printf("## All triangle vertices on same edge ref %d \n",jj);
              printf("Triangle %d vertices ",iface);
              intAr1(3,msh.fac2poi[iface]).print();

              for(int ii = 0; ii < 3; ii++){
                int ipoin = msh.fac2poi(iface,ii);
                printf(" ipoin %d : \n",ipoin);
                for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
                  int bdim = msh.bpo2ibi(ibpoi,1);
                  if(bdim != 1) continue;
                  int ientt = msh.bpo2ibi(ibpoi,2);
                  int iref  = msh.edg2ref[ientt];
                  printf(" edge bpoi %d ientt %d iref %d \n",ibpoi,ientt,iref);
                }
              }
              METRIS_THROW(TopoExcept());
            }
          }
        }
        msh.ced2tag.fill(0);
      }


      // Need more context or to call EGADS primitives. Leave this for after geo refactoring
      //bool iflat;
      //double meas0;
      //if(msh.idim == 2){
      //  meas0 = getmeasentP1<2,2>(msh.fac2poi[iface], msh.coord, Defaults::vtol, NULL, &iflat);
      //}else if(msh.idim == 3){
      //  double nrmal[3];
      //  meas0 = getmeasentP1<3,2>(msh.fac2poi[iface], msh.coord, Defaults::vtol, nrmal, &iflat);
      //}

      //if(iflat){
      //  METRIS_THROW_MSG(GeomExcept(), "Flat triangle vol = "<<meas0);
      //}
      //if(meas0 < 0){
      //  METRIS_THROW_MSG(GeomExcept(), "Negative triangle vol = "<<meas0);
      //}


      for(int ied = 0; ied < 3; ied++){
        int i1 = msh.fac2poi(iface,lnoed2[ied][0]);
        int i2 = msh.fac2poi(iface,lnoed2[ied][1]);

        auto key = stup2(i1,i2);
        auto  t = hshTab2.find(key);
        if(t != hshTab2.end()){
          int ifac2 = t->second;
          if(msh.fac2fac(iface,ied) != ifac2){
            printf("iface = %d i1 %d i2 %d found ifac2 = %d but not in mesh neighbour table\n",
              iface,i1,i2,ifac2);
          }
          METRIS_ASSERT(msh.fac2fac(iface,ied) == ifac2);
        }else{
          hshTab2[key] = iface;
        }

        METRIS_ENFORCE_MSG(i1 >= 0 && i1 < npoin 
                        && i2 >= 0 && i2 < npoin,"Face vertices within bounds ");

        int ifac2 = msh.fac2fac(iface,ied); 

        METRIS_ENFORCE_MSG(ifac2 > -nface-2 && ifac2 < nface,"Neighbour either inexistent or valid index"); 
        if(ifac2 >= 0) METRIS_ENFORCE_MSG(!isdeadent(ifac2,msh.fac2poi),
          "Valid neighbour iface = "<<iface<<" ied = "<<ied<<" ifac2 = "<<ifac2);
        if(ifac2 < -1) METRIS_ENFORCE_MSG(!isdeadent(-ifac2-2,msh.fac2poi),
          "Valid neighbour iface = "<<iface<<" ied = "<<ied<<" ifac2 = "<<ifac2);

        if(msh.idim == 2){
          METRIS_ENFORCE_MSG(ifac2 >= -1,"Manifold in dimension 2");
        }

        if(ifac2 <= -1 || msh.fac2ref[ifac2] != msh.fac2ref[iface]){
          int iedge = msh.facedg2glo(iface,ied);
          METRIS_ENFORCE_MSG(iedge >= 0 && iedge < msh.nedge, "Edge invalid = "<<iedge
            <<" iface "<<iface<<" ied "<<ied<<" ifac2 "<<ifac2);

          METRIS_ENFORCE_MSG(!isdeadent(iedge,msh.edg2poi),"Triangle points to valid edge");

          if(msh.edg2tag(ithread,iedge) == msh.tag[ithread]){
            // This edge has been seen by another triangle. They must be neighbours !
            METRIS_ENFORCE_MSG(ifac2 != -1, "Edge seen twice but not neighbours !");
          }
          msh.edg2tag(ithread,iedge) = msh.tag[ithread];

          METRIS_ENFORCE_MSG(iedge >= 0 && iedge < nedge, "No face neighbour -> is bdry edge")
          int j1 = msh.edg2poi(iedge,0);
          int j2 = msh.edg2poi(iedge,1);
          int iedl;
          METRIS_TRY0(iedl = getedgfac(msh,iface,j1,j2))
          METRIS_ENFORCE_MSG(iedl == ied,"iface = "<<iface<<" ied = "<<ied<<" verts i1 = "<<i1<<" i2 = "<<i2 << 
            " neighbour to iedge "<<iedge<<" verts j1 = "<<j1<<" j2 = "<<j2<<" get iedl = "<<iedl);
  //        METRIS_ENFORCE_MSG(iedl == ied,"Shared edge vertices differ iface = "<<iface<<" ied = "<<ied<<
  //          " verts i1 = "<<i1<<" i2 = "<<i2 << 
  //          " ifac2 = "<<ifac2);
          continue;
        }

        if(ifac2 < -1){
          ifac2 = -ifac2-2;
          METRIS_ENFORCE(ifac2 < nface);
                  // Check this edge is a geometric edge
          auto t = msh.edgHshTab.find(stup2(i1,i2));
          METRIS_ENFORCE_MSG(t!=msh.edgHshTab.end(),"Non manifold edge accounted for in edge hash table");
                  // Unroll neighbours and check we arrive back at iface.
          msh.tag[ithread] ++;
          msh.fac2tag(ithread,iface) = msh.tag[ithread];
          while(ifac2 != iface){
            METRIS_ENFORCE_MSG(msh.fac2tag(ithread,ifac2) < msh.tag[ithread],"First time seeing this face");
            if(msh.fac2tag(ithread,ifac2) >= msh.tag[ithread]){
              printf("Failure with face %d nface %d fac2tag[ithread] = %d tag = %d \n",ifac2,nface,msh.fac2tag(ithread,ifac2),msh.tag[ithread]);
              printf("Debug print non-manifold neighbours\n");
              int ifacd0=iface;
              do{
                int ied0;
                METRIS_TRY0(ied0 = getedgfac(msh,ifacd0,i1,i2);)
                printf("Iface %d neighbour %d \n",ifacd0,msh.fac2fac(ifacd0,ied0));
                ifacd0 = - msh.fac2fac(ifacd0,ied0) - 2;
              }while(ifacd0 != iface);
            }
            msh.fac2tag(ithread,ifac2) = msh.tag[ithread];
            int ied2;
            METRIS_TRY0(ied2 = getedgfac(msh,ifac2,i1,i2);)
            METRIS_ENFORCE_MSG(ied2 >= 0 && ied2 < 3,"Returned edge within 0-2 range" );
            ifac2 = -msh.fac2fac(ifac2,ied2) - 2;
            METRIS_ENFORCE(ifac2 < nface);
            METRIS_ENFORCE_MSG(ifac2 >= 0 && ifac2 < nface,"Next face in line is valid ");
          }
          continue;
        }

        int ied2;
        METRIS_TRY0(ied2 = getedgfac(msh,ifac2,i1,i2))
        METRIS_ENFORCE_MSG(ied2 >= 0 && ied2 < 3,"Returned edge within 0-2 range");
        METRIS_ENFORCE_MSG(msh.fac2fac(ifac2,ied2) == iface,"Neighbour has iface as neighbour at common face");

      }

      // Only in dimension >=3 should faces be part of a hash table
      if(msh.get_tdim() == 3){
        int i1 = msh.fac2poi(iface,0);
        int i2 = msh.fac2poi(iface,1);
        int i3 = msh.fac2poi(iface,2);
        auto key = stup3(i1,i2,i3);
        auto t = msh.facHshTab.find(key);
        METRIS_ENFORCE_MSG(t != msh.facHshTab.end(),"Face is in hashtable ");
        METRIS_ENFORCE_MSG(iface == t->second,"Hashtab entry gives face ");
      }



      // -- Check all points have bpois if faces are boundary
      if(msh.isboundary_faces()){
        for(int ii = 0; ii < facnpps[msh.curdeg]; ii++){
          int ip = msh.fac2poi(iface,ii);
          METRIS_ENFORCE(ip >= 0 && ip < msh.npoin);
          int ib = msh.poi2bpo[ip];
          METRIS_ENFORCE_MSG(ib >= 0 && ib < msh.nbpoi,"ip = "<<ip<<" has ib = "<<ib);

          int itype = msh.bpo2ibi(ib,1);
          if(itype == 2){
            // If this is a face point, we need only check the given face has the same ref
            int ifac2 = msh.bpo2ibi(ib,2);
            METRIS_ENFORCE(ifac2 >= 0 && ifac2 < msh.nface);
            if(msh.fac2ref[ifac2] != msh.fac2ref[iface]){
              printf("Missing bpo entries ? for ip = %d ib = %d \n",ip,ib);
              print_bpolist(msh,ib);
            }
            METRIS_ENFORCE(msh.fac2ref[ifac2] == msh.fac2ref[iface]);
          }else{
            // Otherwise we certainly find the triangle in the bpoi list
            int ib2 = ib; 
            int nn = 0;
            bool ifnd = false;
            do{
              nn++;
              METRIS_ENFORCE(nn <= METRIS_MAX_WHILE);
              int ityp = msh.bpo2ibi(ib2,1);
              if(ityp == 2){
                int ifac2 = msh.bpo2ibi(ib2,2);
                if(ifac2 == iface){
                  ifnd = true;
                  break;
                }
              }
              ib2 = msh.bpo2ibi(ib2,3);
            }while(ib2 >= 0 && ib2 != ib);
            if(!ifnd){
              printf("Failed to find face in bpoi. face = %d ip = %d ib = %d \n",iface,ip,ib);
              print_bpolist(msh,ib);
            }
            METRIS_ENFORCE_MSG(ifnd == true, "face "<<iface<<" not found in ipoin "<<ip<<" bpois");
          }

        }
      }
    }


    if(nface > 0){
      for(int iedge = 0; iedge < nedge; iedge++){
        if(isdeadent(iedge,msh.edg2poi)) continue;
        int i1 = msh.edg2poi(iedge,0);
        int i2 = msh.edg2poi(iedge,1);
        auto key = stup2(i1,i2);
        auto  t = msh.edgHshTab.find(key);
        METRIS_ENFORCE_MSG(t!=msh.edgHshTab.end(),"Edge is in hashtable iedge = "<<iedge<<" i1= "<<i1<<" i2= "<<i2<<"\n");
        METRIS_ENFORCE_MSG(iedge==t->second,"Hashtab points to some other edge ? should be "<<iedge<<" is "<<t->second);
      }
    }


    for(int iedge = 0; iedge < msh.nedge; iedge++){
      if(isdeadent(iedge,msh.edg2poi)) continue;

      // Check neighbours
      for(int inei = 0; inei < 2; inei ++){
        int ipoin = msh.edg2poi[iedge][1-inei];
        METRIS_ENFORCE(ipoin >= 0 && ipoin < msh.npoin);
        int iedg2 = msh.edg2edg(iedge,inei);

        if(iedg2 >= 0){
          for(int jj = 0; jj < 2 ; jj++){
            if(msh.edg2poi(iedg2,jj) == ipoin){
              METRIS_ENFORCE(msh.edg2edg[iedg2][1-jj] == iedge);
              goto ifnde;
            }
          }
          METRIS_THROW_MSG(TopoExcept(),"Did not find vertex in neighbour");
          ifnde:
          // Check point is corner in case differing refs
          if(msh.edg2ref[iedge] != msh.edg2ref[iedg2]){
            int ib = msh.poi2bpo[ipoin];
            METRIS_ENFORCE(ib >= 0 && ib < msh.nbpoi);
            int itype = msh.bpo2ibi(ib,1);
            METRIS_ENFORCE_MSG(itype == 0, "iedge = "<<iedge<<" iedg2 = "<<iedg2<<
              " ip = "<<ipoin<<" ib = "<<ib<<" bdry t dim = "<<itype<<" but differing refs should be corner!");
          }
        }else{// Non manifold or boundary 
          // Check point is corner
          int ib = msh.poi2bpo[ipoin];
          METRIS_ENFORCE(ib >= 0 && ib < msh.nbpoi);
          int itype = msh.bpo2ibi(ib,1);
          METRIS_ENFORCE_MSG(itype == 0, "iedge = "<<iedge<<
            " ip = "<<ipoin<<" ib = "<<ib<<" bdry t dim = "<<itype
            <<" but neighbour = "<<iedg2<<" ! should be corner !");

          if(iedg2 == -1) continue;

          int inein = inei;
          int iedgn = iedge;
          int nn = 0;
          while(getnextedgnm(msh,iedge,ipoin,&iedgn,&inein)){
            nn++;
            METRIS_ENFORCE_MSG(nn <= METRIS_MAX_WHILE,"Infinite nm edge start from "<<iedge<<" ip = "<<ipoin);
            METRIS_ENFORCE_MSG(msh.edg2poi[iedgn][1-inein] == ipoin,
              "non manifold (edge) point"<<ipoin<<" initial edge "<<iedge<<" current "<<iedgn
              <<" neigh supposed "<<inein<<" vertices ip1 = "<<msh.edg2poi(iedgn,0)<<" ip2 = "
              <<msh.edg2poi(iedgn,0)<<"\n neighbours 1 : "<<msh.edg2edg(iedgn,0)<<" 2: "
              <<msh.edg2edg(iedgn,1)); // redundant but for clarity
          }
          // Loop over edges, check neighbourhood ok. 
        }
      }



      // -- Check all points have bpois if edges are boundary
      if(msh.isboundary_edges()){
        for(int ii = 0; ii < edgnpps[msh.curdeg]; ii++){
          int ip = msh.edg2poi(iedge,ii);
          METRIS_ENFORCE(ip >= 0 && ip < msh.npoin);
          int ib = msh.poi2bpo[ip];
          METRIS_ENFORCE_MSG(ib >= 0 && ib < msh.nbpoi, "out of bounds ib = "<<ib
            <<" nbpoi = "<<msh.nbpoi<<" ipoin = "<<ip);

          int itype = msh.bpo2ibi(ib,1);
          METRIS_ENFORCE_MSG(itype < 2,"boundary ip = "<<ip<<" ib = "<<ib<<" on edge but itype = 2");

          if(itype == 1){
            // If this is an edge point, we need only check the given edge has the same ref
            int iedg2 = msh.bpo2ibi(ib,2);
            METRIS_ENFORCE(iedg2 >= 0 && iedg2 < msh.nedge);
            METRIS_ENFORCE_MSG(msh.edg2ref[iedg2] == msh.edg2ref[iedge], "Edge point "<<ip<<" belongs to edges of ref1 = "<<
              msh.edg2ref[iedge]<<" ref2 = "<<msh.edg2ref[iedg2]);
            // this should throw, but just in case of future changes
            CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
              try{
                int iver = msh.getveredg<ideg>(iedg2, ip);
                METRIS_ENFORCE(iver >= 0);
              }catch(const MetrisExcept& e){
                printf("ip = %d ib = %d not in edge %d = (%d,%d) from iedge = %d \n",
                  ip,ib,iedg2,
                  msh.edg2poi(iedg2,0),msh.edg2poi(iedg2,1),iedge);
                throw(e);
              }
            }}CT_FOR1(ideg);
          }else{
            // Otherwise corner, then this edge is in the bpo list
            int ib2 = ib; 
            int nn = 0;
            bool ifnd = false;
            do{
              nn++;
              METRIS_ENFORCE(nn <= METRIS_MAX_WHILE);
              int ityp = msh.bpo2ibi(ib2,1);
              if(ityp == 1){
                int iedg2 = msh.bpo2ibi(ib2,2);
                if(iedg2 == iedge){
                  ifnd = true;
                  break;
                }
              }
              ib2 = msh.bpo2ibi(ib2,3);
            }while(ib2 >= 0 && ib2 != ib);
            if(!ifnd){
              printf("Failed to find edge in bpoi. edge = %d ip = %d ib = %d \n",iedge,ip,ib);
              print_bpolist(msh,ib);
            }
            METRIS_ENFORCE_MSG(ifnd == true, "Found edge in point bpois");
          }
        }
      }
    }





    // Check non-duplication of lowest-dimensionality ibpoi
    // + check link back to ipoin correct
    for(int ibpoi = 0; ibpoi < nbpoi; ibpoi++){
      int ipoin = msh.bpo2ibi(ibpoi,0);
      if(ipoin < 0) continue;
      METRIS_ENFORCE_MSG(ipoin < npoin,"ibpoi = "<<ibpoi<<" points to ipoin = "<<ipoin<<" but npoin = "<<npoin);
      METRIS_ENFORCE(msh.poi2bpo[ipoin] >= 0);
      METRIS_ENFORCE(msh.bpo2ibi(msh.poi2bpo[ipoin],0) == ipoin);
      METRIS_ENFORCE(msh.poi2ent(ipoin,0) >= 0); 

      int bdim = msh.bpo2ibi(ibpoi,1);
      METRIS_ENFORCE(bdim == 0 || bdim == 1 || bdim == 2);
      if(bdim > 0){
        int ientt = msh.bpo2ibi(ibpoi,2);
        METRIS_ENFORCE(ientt >= 0);
        int iver = msh.getverent(ientt,bdim,ipoin);
        METRIS_ENFORCE(iver>=0);
      }

      int ibpo2 = msh.poi2bpo[ipoin];
      int minty = 3;
      int nloop = 0;
      do{
        minty = minty < msh.bpo2ibi(ibpo2,1) ? minty : msh.bpo2ibi(ibpo2,1);
        ibpo2 = msh.bpo2ibi(ibpo2,3);
  //      printf("In loop ibpo2 = %d ibpoin = %d \n",ibpo2,ibpoi);
        nloop++;
        if(nloop > METRIS_MAX_WHILE){
          METRIS_THROW_MSG(TopoExcept(),"LOOP TOO LONG ipoin = "<<ipoin);
          break;
        }
      }while(ibpo2 >= 0 && ibpo2 != ibpoi);
      if(nloop > METRIS_MAX_WHILE) break;

      ibpo2 = msh.poi2bpo[ipoin];
      int ninty = 0;
      do{
        ninty += (msh.bpo2ibi(ibpo2,1) == minty);
        ibpo2 = msh.bpo2ibi(ibpo2,3);
      }while(ibpo2 >= 0 && ibpo2 != ibpoi);
      METRIS_ENFORCE_MSG(ninty==1,"ninty = "<<ninty<<" with minty = "<<minty<<
        " ipoin = "<<msh.bpo2ibi(ibpoi,0)<<" ibpoi = "<<ibpoi<<"\n");

      if(msh.CAD()){
        double result[18];
        ego obj;
        int tdim  = msh.bpo2ibi(ibpoi,1);
        int ientt = msh.bpo2ibi(ibpoi,2);
        if(tdim == 1){
          METRIS_ENFORCE_MSG(!isdeadent(ientt,msh.edg2poi),
             "ipoin = "<<ipoin<<" ibpoi "<<ibpoi<<" points to dead edge "<<ientt);
          int iref = msh.edg2ref[ientt];
          obj = msh.CAD.cad2edg[iref];
        }else if(tdim == 2){
          METRIS_ENFORCE_MSG(!isdeadent(ientt,msh.fac2poi),"ibpoi points to non-dead entity");
          int iref = msh.fac2ref[ientt];
          obj = msh.CAD.cad2fac[iref];
        }
        if(tdim == 2 || tdim == 1){
          int ierro = EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
          METRIS_ENFORCE_MSG(ierro==0, "EG_evaluate failed ierro = "<<ierro);
          double nrm;
          if(msh.idim == 3){
            nrm = geterrl2<3>(result,msh.coord[ipoin]);
          }else{
            nrm = geterrl2<2>(result,msh.coord[ipoin]);
          }
          METRIS_ENFORCE_MSG(nrm < geotol*geotol,"large nrm = "<<nrm<<" point = "
            <<ipoin<<" ibpoi = "<<ibpoi<<" type = "<<tdim<<
            " uv = "<<msh.bpo2rbi(ibpoi,0)<<" "<<msh.bpo2rbi(ibpoi,1));
        }
      }

    }


  }catch(const MetrisExcept& e){
    printf("Check_topo failed, dumping mesh:\n");
    writeMesh("check_topo_fail.meshb",msh);
    printf("DONE \n");
    throw(e);
  }
  //if(stopend){
  //  auto test_id = boost::unit_test::framework::current_test_case().p_id;
  //  METRIS_ENFORCE(boost::unit_test::results_collector.results(test_id).passed());
  //}

}








}// end namespace

