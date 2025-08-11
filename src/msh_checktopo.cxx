//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"

#include "msh_checktopo.hxx"
#include "aux_exceptions.hxx"
#include "aux_topo.hxx"
#include "low_geo/measure.hxx"
#include "io_libmeshb.hxx"
#include "low_geo/ccoef.hxx"
#include "low_geo/nrml2.hxx"
#include "low_geo/normal.hxx"

#include "utils/mprintf.hxx"
#include "utils/fmt_formatters.hxx"
#include "utils/aux_misc.hxx"
#include "utils/CT_loop.hxx"

namespace Metris{



void check_topo(MeshBase &msh, int ithread){
  check_topo(msh,msh.nbpoi,msh.npoin,msh.nedge,msh.nface,msh.nelem,ithread);
}


void check_topo(MeshBase &msh, 
                int nbpoi, int npoin, int nedge, int nface, int nelem, int ithread){
  INCVDEPTH(msh.param);

  static int ncall_this = 0;
  ncall_this++;

  if(msh.param->iflag3 > 0 && ncall_this%msh.param->iflag3 != 0) return;
    
  try{

    CPRINTF2("-- START check_topo\n");

    const int jdeg = msh.idim * (msh.curdeg - 1);
    dblAr1 ccoef(getnnode(msh.idim,jdeg));


    // 1 check not nan coords
    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      for(int ii = 0; ii < msh.idim; ii++){
        if(std::isnan(msh.coord(ipoin,ii))){
          MPRINTF("## NAN COORDINATE FOR ipoin {} ",ipoin);
          dblAr1(msh.idim,msh.coord[ipoin]).print();
        }
        METRIS_ENFORCE(!std::isnan(msh.coord(ipoin,ii)))
      }
    }


    CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
    CT_FOR0_INC(2,gdim,tdim){if(tdim <= msh.get_tdim()){
      int nentt = msh.nentt(tdim);
      const intAr2 &ent2poi = msh.ent2poi(tdim);
      for(int ientt = 0; ientt < nentt; ientt++){
        if(isdeadent(ientt,ent2poi)) continue;

        bool iflat = !isvalideltP1<gdim,tdim>(msh, ientt, NULL, NULL, -1 );
        METRIS_ENFORCE(!iflat);
        if(msh.curdeg > 1){
          CT_FOR0_INC(2,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
            getsclccoef<gdim,tdim,ideg>(msh,ientt,NULL,&ccoef[0],&iflat);
          }}CT_FOR1(ideg);
          if(iflat){
            MPRINTF("## NEGATIVE JACOBIAN {} \n ",ientt);
            writeMesh("inva"+std::to_string(ientt),msh);
            METRIS_THROW(GeomExcept());
          }
        }
      }
    }}CT_FOR1(tdim);
    }}CT_FOR1(gdim);


    //aftervalidity:

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
          METRIS_ENFORCE_MSG(ient2 >= 0,
            "Invalid poin2ent? tdim = "<<tdim<<" ientt "<<ientt<<" ipoin "<<ipoin
            <<" poi2ent = "<<ient2<<" anounced dim "<<msh.poi2ent(ipoin,1));
          int tdim2 = msh.poi2ent(ipoin,1);
          METRIS_ENFORCE(tdim2 <= tdim && tdim2 >= 1);
          int ibpoi = msh.poi2bpo[ipoin];
          if(ibpoi < 0) continue;
          METRIS_ENFORCE(tdim2 < msh.idim);
        }
      }
    }// for tdim

    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      int ientt = msh.poi2ent(ipoin,0);
      if(ientt < 0) continue;
      int tdim = msh.poi2ent(ipoin,1);
      METRIS_ENFORCE_MSG(1 <= tdim && tdim <= msh.get_tdim(), 
        "ipoin "<<ipoin<<" poi2ent = "<<ientt<<" "<<tdim);
      if(tdim >= 1){
        int iver = msh.getverent(ientt,tdim,ipoin);
        if(iver < 0){
          MPRINTF("failed to find ipoin in its poi2ent. ipoin = {} ientt {} tdim {} ",
                   ipoin, ientt, tdim);
          const intAr2 &ent2poi = msh.ent2poi(tdim);
          MPRINTF("vertices: ");
          intAr1(tdim+1,ent2poi[ientt]).print();
        }
        METRIS_ENFORCE_MSG(iver >= 0, "did not find ipoin = "<<ipoin<<" poi2ent "<<
          ientt<<" tdim "<< tdim<<" not found");
      }else{
        int ibpoi = msh.poi2bpo[ipoin];
        METRIS_ENFORCE(ibpoi >= 0);
        int itype = msh.bpo2ibi(ibpoi,1);
        METRIS_ENFORCE(itype == 0);
        METRIS_ENFORCE(ientt == ipoin);
      }
    }

    // Check poi2bak in case MetricFieldFE
    if(msh.meshClass() == MeshClass::Mesh && msh.metricClass() == MetricClass::MetricFieldFE){
      const intAr1 &poi2bak = msh.metricClass() == MetricClass::MetricFieldFE ? 
       ((Mesh<MetricFieldFE> *)(&msh))->poi2bak 
      :((Mesh<MetricFieldAnalytical> *)(&msh))->poi2bak;

      for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
        if(msh.poi2ent(ipoin,0) < 0) continue;


        int pdim = msh.getpoitdim(ipoin);
        if(pdim < 0 || pdim > msh.get_tdim()){
          MPRINTF(" ## INVALID pdim {} \n",pdim);
          MPRINTF("iopin = {} \n",ipoin);
          writeMesh("debug_checktopo",msh);
          METRIS_THROW(TopoExcept());
        }
        if(pdim == 0) continue;
        int iebak = poi2bak[ipoin];
        if(iebak < 0){
          MPRINTF("ipoin {} pdim {} poi2bak {} negative\n",ipoin,pdim,iebak);
        }
        METRIS_ENFORCE(iebak >= 0);
        int nentb = msh.metricClass() == MetricClass::MetricFieldFE ?  
          ((Mesh<MetricFieldFE> *)(&msh))->bak->nentt(pdim):
          ((Mesh<MetricFieldAnalytical> *)(&msh))->bak->nentt(pdim);
        if(iebak >= nentb){
          MPRINTF("ipoin {} pdim {} poi2bak {} >= nentb {}\n",ipoin,pdim,iebak,nentb);
        }
        METRIS_ENFORCE(iebak < nentb);
        const intAr2& ent2pob = msh.metricClass() == MetricClass::MetricFieldFE ?  
          ((Mesh<MetricFieldFE> *)(&msh))->bak->ent2poi(pdim):
          ((Mesh<MetricFieldAnalytical> *)(&msh))->bak->ent2poi(pdim);
        if(isdeadent(iebak,ent2pob)){
          MPRINTF("Point {} tdim {} has dead back seed {} \n",ipoin,pdim,iebak);
        }
        METRIS_ENFORCE_MSG(!isdeadent(iebak,ent2pob),
          "Point "<<ipoin<<" tdim "<<pdim<<" has back seed "<<iebak<<" which is dead");

        //METRIS_ENFORCE(iebak < msh_->bak->nentt(pdim));
      }
    }

    // Check metric in case analytical
    
    if(msh.meshClass() == MeshClass::Mesh && msh.metricClass() == MetricClass::MetricFieldAnalytical){
      Mesh<MetricFieldAnalytical>& msh_a = (Mesh<MetricFieldAnalytical>&)msh;
      for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
        if(msh.poi2ent(ipoin,0) < 0) continue;
        double metl[6];
        msh_a.met.getMetPhys(DifVar::None, msh_a.met.getSpace(),
                             msh.coord[ipoin], metl, NULL);
        double errmet = msh.idim == 2 ? geterrl2<3>(metl, msh_a.met[ipoin]) :
                                        geterrl2<6>(metl, msh_a.met[ipoin]);
        double nrm1 = msh.idim == 2 ? getnrml2<3>(msh_a.met[ipoin]) :
                                      getnrml2<6>(msh_a.met[ipoin]);
        METRIS_ENFORCE_MSG(errmet < 1.0E-15*nrm1, "Large metric error\n");
      }
    }


    // Check tags correct 
    for(int itag = 0; itag < METRIS_MAXTAGS; itag++){
      for(int ielem = 0; ielem < nelem; ielem++){
        if(isdeadent(ielem,msh.tet2poi)) continue;
        METRIS_ENFORCE_MSG(msh.tet2tag(itag,ielem) <= msh.tag[itag],
          "Tetra "<<ielem<<" has tag = "<<msh.tet2tag(itag,ielem)<<" > tag = "
          <<msh.tag[itag]<<" for itag = "<<itag
          <<" neighbours: "<<msh.tet2poi(ielem,0)<<" "<<msh.tet2poi(ielem,1)
                      <<" "<<msh.tet2poi(ielem,2)<<" "<<msh.tet2poi(ielem,3));
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
        METRIS_ENFORCE_MSG(msh.poi2tag(itag,ipoin) <= msh.tag[itag] 
                        && msh.poi2tag(itag,ipoin) >= 0,
                        "Wrong point tag ithread "<<itag<<" ipoin "
                        << ipoin <<" has "<<msh.poi2tag(itag,ipoin)<<
                        " tag is "<<msh.tag[itag]);
      }

      for(int ii = 0; ii < msh.cfa2tag.get_stride(); ii++){
        METRIS_ENFORCE(msh.cfa2tag(itag, ii) <= msh.tag[itag]);
        msh.cfa2tag(ithread, ii) = 0;
      }
      for(int ii = 0; ii < msh.ced2tag.get_stride(); ii++){
        METRIS_ENFORCE(msh.ced2tag(itag, ii) <= msh.tag[itag]);
        msh.ced2tag(ithread, ii) = 0;
      }
      for(int ii = 0; ii < msh.cno2tag.get_stride(); ii++){
        METRIS_ENFORCE(msh.cno2tag(itag, ii) <= msh.tag[itag]);
        msh.cno2tag(ithread, ii) = 0;
      }
      for(int ii = 0; ii < msh.dom2tag.get_stride(); ii++){
        METRIS_ENFORCE(msh.dom2tag(itag, ii) <= msh.tag[itag]);
        msh.dom2tag(ithread, ii) = 0;
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
      METRIS_ENFORCE_MSG(!isdeadent(iface,msh.fac2poi),
        "iedge "<<iedge<<" points to face "<<iface<<" which is dead"
        <<"\nedge vertices: "<<msh.edg2poi(iedge,0)<<" "<<msh.edg2poi(iedge,1));
      int ip1 = msh.edg2poi(iedge,0);
      int ip2 = msh.edg2poi(iedge,1);

      int iedl; 
      try{
        iedl = getedgfac(msh,iface,ip1,ip2);
        METRIS_ENFORCE(iedl >= 0);
      }catch(const MetrisExcept &e){
        MPRINTF("edge not in triangle (edg2fac link)\n");
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

    // check getfacglo -> exists and no dupes
    for(int iface = 0; iface < nface && msh.idim >= 3; iface++){
      if(isdeadent(iface,msh.fac2poi)) continue;
      int ip1 = msh.fac2poi(iface,0);
      int ip2 = msh.fac2poi(iface,1);
      int ip3 = msh.fac2poi(iface,2);

      int ifac2 = getfacglo(msh, ip1, ip2, ip3);
      METRIS_ENFORCE_MSG(ifac2 == iface, 
        "Face either doesnt exist or is duplicated iface = "<<iface<<" ifac2 = "<<ifac2); 
    }


    for(int itetr = 0; itetr < nelem; itetr++){
      if(isdeadent(itetr,msh.tet2poi)) continue;
      int idom1 = msh.tet2ref[itetr];
      for(int ifa = 0; ifa < 4; ifa++){
        int ip1 = msh.tet2poi(itetr,lnofa3[ifa][0]);
        int ip2 = msh.tet2poi(itetr,lnofa3[ifa][1]);
        int ip3 = msh.tet2poi(itetr,lnofa3[ifa][2]);

        int iface = getfacglo(msh, ip1, ip2, ip3);

        int itet2 = msh.tet2tet(itetr, ifa);
        if(itet2 < 0){
          METRIS_ENFORCE_MSG(iface >= 0, "No face on domain boundary");
        }else{
          int idom2 = msh.tet2ref[itetr];
          if(idom1 == idom2){
            METRIS_ENFORCE_MSG(iface < 0, "Face found between same ref tetras")
          }else{
            METRIS_ENFORCE_MSG(iface >= 0, "No face between two ref tetras")
          }
        }
      }
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

      if(!(tdimn == msh.poi2ent(ipoin,1) || tdimn == 0)){
        MPRINTF("ipoin = {} poi2ent = {} {}",ipoin,msh.poi2ent(ipoin,0)
               ,msh.poi2ent(ipoin,1));
        for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
          MPRINTF("ibpoi {} : ",ibpoi);
          intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
        }
      }
      METRIS_ENFORCE_MSG(tdimn == msh.poi2ent(ipoin,1) || tdimn == 0,
        "Wrong tdimn? ipoin = "<<ipoin<<" tdimn using bpo = "<<tdimn<<
        " in poi2ent = "<<msh.poi2ent(ipoin,1)); 

      // Corners -> edges
      if(tdimn == 0) tdimn = 1;

      intAr2 &ent2poi = msh.ent2poi(tdimn);
      int nentt = msh.nentt(tdimn);
      if(ientt >= nentt){
        MPRINTF("poi2ent link defective? ipoin = {} ibpoi = {} ientt = {}  tdimn = {} \n",ipoin,ibpoi,ientt,tdimn);
        MPRINTF("nedge = {} nface = {} nelem = {} \n",msh.nedge,msh.nface,msh.nelem);

        MPRINTF("face nodes: \n");
        int nnod2 = getnnod2(msh.curdeg);
        intAr1(nnod2,msh.fac2poi[ientt]).print();

        int ibpo2 = ibpoi;
        while(ibpo2 >= 0){
          MPRINTF("ibpoi {} : ",ibpo2);
          intAr1(nibi,msh.bpo2ibi[ibpo2]).print();
          ibpo2 = msh.bpo2ibi(ibpo2,3);
        }

        if(ibpoi>=0){
          MPRINTF("bpo2ibi = ");
          intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
        }
        MPRINTF("ientt = {} nodes = ",ientt);
        int nnode = msh.nnode(tdimn);
        intAr1(nnode,ent2poi[ientt]).print();
        for(int tdim2 = 1; tdim2 <= msh.get_tdim(); tdim2++){
          int nnod2 = msh.nnode(tdim2);
          intAr2 &ent2po2 = msh.ent2poi(tdim2);
          MPRINTF("idim {} entt nodes ",tdim2);
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
        MPRINTF("2 poi2ent link defective? ipoin = {} ibpoi = {} tdimn = {} \n",
               ipoin,ibpoi,tdimn);
        if(ibpoi>=0){
          MPRINTF("bpo2ibi = ");
          intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
        }
        MPRINTF("ientt = {} nodes = ",ientt);
        intAr1(nnode,ent2poi[ientt]).print();
        for(int tdim2 = 1; tdim2 <= msh.get_tdim(); tdim2++){
          int nnod2 = msh.nnode(tdim2);
          intAr2 &ent2po2 = msh.ent2poi(tdim2);
          MPRINTF("idim {} entt nodes ",tdim2);
          intAr1(nnod2,ent2po2[ientt]).print();
        }
      }
      METRIS_ENFORCE_MSG(ifnd, "poi2ent not fnd, ipoin = "<<ipoin
                         <<" poi2ent = "<<ientt);
    }

    for(int iface = 0; iface < nface && nelem > 0; iface++){
      if(isdeadent(iface,msh.fac2poi)) continue;
      METRIS_ENFORCE_MSG(msh.fac2tet(iface,0) >= 0,
                         "iface "<<iface<<" vertices "<<msh.fac2poi(iface,0)
                         <<" "<<msh.fac2poi(iface,1)<<" "<<msh.fac2poi(iface,2)
                         <<" has fac2tet = "<<msh.fac2tet(iface,0)<<" "
                         <<msh.fac2tet(iface,1));
      METRIS_ENFORCE_MSG(msh.fac2tet(iface,0) != msh.fac2tet(iface,1),
                         "iface "<<iface<<" vertices "<<msh.fac2poi(iface,0)
                         <<" "<<msh.fac2poi(iface,1)<<" "<<msh.fac2poi(iface,2)
                         <<" has repeated fac2tet");

      for(int ii = 0; ii < 2; ii++){
        int ielem = msh.fac2tet(iface,ii);
        if(ielem < 0) continue;
        if(isdeadent(ielem,msh.tet2poi)){
          MPRINTF("## iface {} points to dead tetra = {}\n",iface,ielem);
          MPRINTF("Nodes: ");
          intAr1(3,msh.fac2poi[iface]).print();
          MPRINTF("fac2tet {} {} \n",msh.fac2tet(iface,0),msh.fac2tet(iface,1));
        }
        METRIS_ENFORCE(!isdeadent(ielem,msh.tet2poi));
      }
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
        if(iele2 >= 0)
        METRIS_ENFORCE_MSG(!isdeadent(iele2,msh.tet2poi),
          "Neighbour invalid "<<iele2<<" of tet "
          <<ielem<<" against iface "<<ifa<<" vertices "<<i1<<" "<<i2<<" "<<i3);

        if(iele2 == -1 || msh.tet2ref[iele2] != msh.tet2ref[ielem]){
          int iface = msh.tetfac2glo(ielem,ifa);
          if(iface < 0 || iface >= nface){
            MPRINTF("No tet neighbour and no bface. ielem = {} ifa = {} iface {}\n",
                   ielem,ifa,iface);
          }
          METRIS_ENFORCE_MSG(iface >= 0 && iface < nface,"No neighbour means neighbour is a bface");
          int j1 = msh.fac2poi(iface,0);
          int j2 = msh.fac2poi(iface,1);
          int j3 = msh.fac2poi(iface,2);
          int ifal = getfactet(msh,ielem,j1,j2,j3);
          METRIS_ENFORCE_MSG(ifal == ifa,"Face neighbour has same vertices as local face ");
          if(iele2 == -1){
            if(msh.fac2tet(iface,0) != ielem && msh.fac2tet(iface,1) != ielem){
              MPRINTF(" fac2tet does not match. iface {} ielem {} \n",iface,ielem);
              MPRINTF("face: ");
              intAr1(3,msh.fac2poi[iface]).print();
              MPRINTF("tetra: ");
              intAr1(4,msh.tet2poi[ielem]).print();
              MPRINTF("fac2tet : {} {} \n",msh.fac2tet(iface,0),msh.fac2tet(iface,1));
              for(int ii = 0; ii < 2; ii++){
                iele2 = msh.fac2tet(iface,ii);
                if(iele2 < 0) continue;
                MPRINTF("iele {} : ",iele2);
                intAr1(4,msh.tet2poi[iele2]).print();
              }
            }
            METRIS_ENFORCE_MSG(msh.fac2tet(iface,0) == ielem || msh.fac2tet(iface,1) == ielem,
              "iface "<<iface<<" ielem = "<<ielem<<" fac2tet = "<<msh.fac2tet(iface,0)<<","<<msh.fac2tet(iface,1));
          }else{
            METRIS_ENFORCE_MSG(msh.fac2tet(iface,0) == ielem || msh.fac2tet(iface,1) == ielem, 
              "ielem = "<<ielem<<" fac2tet = "<<msh.fac2tet(iface,0)<<","<<msh.fac2tet(iface,1));
          }
          continue;
        }

        int ifa2 = getfactetOpp(msh,iele2,i1,i2,i3);
        if(ifa2 < 0){
          MPRINTF("## Did not find face {} {} {} in tet.",i1,i2,i3);
          MPRINTF(" ielem = {} vertices:",iele2);
          intAr1(4,msh.tet2poi[iele2]).print();
        }
        METRIS_ENFORCE_MSG(ifa2 >= 0 && ifa2 < 4,"Returned face within 0-3 range");
        METRIS_ENFORCE_MSG(msh.tet2tet(iele2,ifa2) == ielem,"Neighbour has ielem as neighbour at common face");
      }
    }

    HshTab_I2I hshTab2(1.5*msh.nface + msh.nedge);

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
            METRIS_ENFORCE(ientt >= 0);
            int iref = msh.edg2ref[ientt];
            METRIS_ENFORCE(iref >= 0);
            if(msh.ced2tag(1,iref) == 1) continue;
            msh.ced2tag(1,iref) = 1;
            msh.ced2tag(0,iref)++;
            iedg[ii] = 1;
          }
        }
        if(iedg[0] + iedg[1] + iedg[2] == 3){
          for(int jj = 0; jj < msh.ced2tag.get_stride(); jj++){
            if(msh.ced2tag(0,jj) == 3){
              MPRINTF("## All triangle vertices on same edge ref {} \n",jj);
              MPRINTF("Triangle {} vertices ",iface);
              intAr1(3,msh.fac2poi[iface]).print();

              for(int ii = 0; ii < 3; ii++){
                int ipoin = msh.fac2poi(iface,ii);
                MPRINTF(" ipoin {} : \n",ipoin);
                for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
                  int bdim = msh.bpo2ibi(ibpoi,1);
                  if(bdim != 1) continue;
                  int ientt = msh.bpo2ibi(ibpoi,2);
                  int iref  = msh.edg2ref[ientt];
                  MPRINTF(" edge bpoi {} ientt {} iref {} \n",ibpoi,ientt,iref);
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
            MPRINTF("iface = {} i1 {} i2 {} found ifac2 = {} but not in mesh neighbour table\n",
              iface,i1,i2,ifac2);
          }
          METRIS_ENFORCE(msh.fac2fac(iface,ied) == ifac2);
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
              MPRINTF("Failure with face {} nface {} fac2tag[ithread] = {} tag = {} \n",ifac2,nface,msh.fac2tag(ithread,ifac2),msh.tag[ithread]);
              MPRINTF("Debug print non-manifold neighbours\n");
              int ifacd0=iface;
              do{
                int ied0;
                METRIS_TRY0(ied0 = getedgfac(msh,ifacd0,i1,i2);)
                MPRINTF("Iface {} neighbour {} \n",ifacd0,msh.fac2fac(ifacd0,ied0));
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
        for(int ii = 0; ii < getnnod2(msh.curdeg); ii++){
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
              MPRINTF("Missing bpo entries ? for ip = {} ib = {} \n",ip,ib);
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
              MPRINTF("Failed to find face in bpoi. face = {} ip = {} ib = {} \n",iface,ip,ib);
              print_bpolist(msh,ib);
            }
            METRIS_ENFORCE_MSG(ifnd == true, "face "<<iface<<" not found in ipoin "<<ip<<" bpois");
          }

        }
      }
    }

    
    // Check closed surface in 3D and manifold
    if(!msh.is_nonmanifold() && msh.idim == 3){
      for(int iface = 0; iface < nface; iface++){
        for(int ied = 0; ied < 3; ied++){
          int ifnei = msh.fac2fac(iface,ied);
          if(ifnei >= 0) continue;
          METRIS_ENFORCE(ifnei == -1);
          int iedge = msh.fac2edg(iface, ied);
          METRIS_ENFORCE(iedge >= 0 && iedge < nedge);
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
        for(int ii = 0; ii < getnnod1(msh.curdeg); ii++){
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
            METRIS_ENFORCE_MSG(!isdeadent(iedg2,msh.edg2poi), "Point "<<ip<<
              " has bpo "<<ib<<" tdim 1 pointing to dead edge "<<iedg2);
            METRIS_ENFORCE_MSG(msh.edg2ref[iedg2] == msh.edg2ref[iedge], 
              "Edge point "<<ip<<" belongs to edges of ref1 = "<<
              msh.edg2ref[iedge]<<" ref2 = "<<msh.edg2ref[iedg2]);
            // this should throw, but just in case of future changes
            CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
              try{
                int iver = msh.getveredg<ideg>(iedg2, ip);
                METRIS_ENFORCE(iver >= 0);
              }catch(const MetrisExcept& e){
                MPRINTF("ip = {} ib = {} not in edge {} = ({},{}) from iedge = {} \n",
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
              MPRINTF("Failed to find edge in bpoi. edge = {} ip = {} ib = {} \n",iedge,ip,ib);
              print_bpolist(msh,ib);
            }
            METRIS_ENFORCE_MSG(ifnd == true, "Found edge in point bpois");
          }
        }
      }
    }


    // Check bpois link to points and back
    for(int ibpoi = 0; ibpoi < nbpoi; ibpoi++){
      int ipoin = msh.bpo2ibi(ibpoi,0);
      if(ipoin < 0) continue;
      METRIS_ENFORCE_MSG(ipoin < npoin,"ibpoi = "<<ibpoi<<" points to ipoin = "<<ipoin<<" but npoin = "<<npoin);
      if(msh.poi2bpo[ipoin] < 0){
        MPRINTF("ibpoi {} : ",ibpoi);
        intAr1(nibi,msh.bpo2ibi[ibpoi]).print();
        MPRINTF("ipoin = {} poi2ent {} {} \n",ipoin,msh.poi2ent(ipoin,0),msh.poi2ent(ipoin,1));
      }
      METRIS_ENFORCE_MSG(msh.poi2bpo[ipoin] >= 0,"ibpoi = "<<ibpoi<<
        " points to ipoin = "<<ipoin<<" but poi2bpo[ipoin] = "<<msh.poi2bpo[ipoin]);

      bool ifnd = false;
      // Check ibpoi is found in ipoin's list
      for(int ibpo0 = msh.poi2bpo[ipoin]; ibpo0 >= 0; ibpo0 = msh.bpo2ibi(ibpo0,3)){
        if(ibpo0 == ibpoi) ifnd = true;
        // Check ibpoi in question is pointing to ipoin
        METRIS_ENFORCE(msh.bpo2ibi(ibpo0,0) == ipoin);
      }
      if(!ifnd){
        MPRINTF("ibpoi = {} links to ipoin {} but not found in its linked list\n",
               ibpoi,ipoin);
        MPRINTF("Full bpo list starting from ibpoi:\n");
        print_bpolist(msh,ibpoi);
        MPRINTF("Full bpo list starting from poi2bpo[ipoin] = {}:\n",msh.poi2bpo[ipoin]);
        print_bpolist(msh,msh.poi2bpo[ipoin]);

      }
      METRIS_ENFORCE(ifnd);
    }

    // Check non-duplication of lowest-dimensionality ibpoi
    // + check link back to ipoin correct
    for(int ibpoi = 0; ibpoi < nbpoi; ibpoi++){
      int ipoin = msh.bpo2ibi(ibpoi,0);
      if(ipoin < 0) continue;
      METRIS_ENFORCE(msh.bpo2ibi(msh.poi2bpo[ipoin],0) == ipoin);
      METRIS_ENFORCE(msh.poi2ent(ipoin,0) >= 0); 

      int bdim = msh.bpo2ibi(ibpoi,1);
      METRIS_ENFORCE(bdim == 0 || bdim == 1 || bdim == 2);
      if(bdim > 0){
        int ientt = msh.bpo2ibi(ibpoi,2);
        METRIS_ENFORCE(!(ientt < 0 || ientt >= msh.nentt(bdim)));
        int iver = msh.getverent(ientt,bdim,ipoin);
        if(iver < 0){
          MPRINTF("ibpoi = {} ipoin = {} bdim = {} ientt = {} iref = {}\n",
                 ibpoi,ipoin,bdim,ientt,msh.ent2ref(bdim)[ientt]);
          MPRINTF("nbpoi = {} npoin = {}\n",msh.nbpoi,msh.npoin);
          MPRINTF("poi2ent = {} {} \n",msh.poi2ent(ipoin,0),msh.poi2ent(ipoin,1));
          MPRINTF("Nodes: ");
          intAr1(bdim+1,msh.ent2poi(bdim)[ientt]).print();
          MPRINTF("Vertex not found\n");
        }
        METRIS_ENFORCE(iver >= 0);
      }

      if(ibpoi == msh.poi2bpo[ipoin]){
        // Check the first entry is the lowest dimensional
        // Also check no other entries of the same dimension
        // In other worst, check all other entries are dimension > bdim.
        for(int ibpo2 = msh.bpo2ibi(ibpoi,3); ibpo2 >= 0; ibpo2 = msh.bpo2ibi(ibpo2,3)){
          if(msh.bpo2ibi(ibpo2,1) <= bdim){
            MPRINTF("ibpoi = {} bdim = {} ibpo2 = {} dim = {}\n",
                   ibpoi,bdim,ibpo2,msh.bpo2ibi(ibpo2,1));
            MPRINTF("Full bpo list:\n");
            print_bpolist(msh,ibpoi);
          }
          METRIS_ENFORCE(msh.bpo2ibi(ibpo2,1) > bdim);
        }
      }
    }


    // Check no duplicate entries in ibpoi.
    for(int itmp = 0; itmp < msh.nbpoi; itmp++){
      int ipoin = msh.bpo2ibi(itmp,0);
      if(ipoin < 0) continue;

      msh.tag[ithread]++;

      bool icor = false;
      for(int ibpoi = msh.poi2bpo[ipoin]; ibpoi >= 0; ibpoi = msh.bpo2ibi(ibpoi,3)){
        int tdim  = msh.bpo2ibi(ibpoi,1);
        int ientt = msh.bpo2ibi(ibpoi,2);
        if(tdim == 0){
          METRIS_ENFORCE_MSG(!icor, "Duplicate corner entries");
          icor = true;
          continue;
        }
        METRIS_ENFORCE_MSG(msh.ent2tag(tdim)(ithread, ientt) < msh.tag[ithread],
          "Entity found twice in point ibpois.");
        msh.ent2tag(tdim)(ithread, ientt) = msh.tag[ithread];
      }
    }


    // Check (u,v)s and ts are correct.
    if(msh.CAD()){
      // Compute a minimum distance between each point and their neighbour. 
      // This will be used as an epsilon.
      // Loop over highest dim boundary entities.
      const double EG_tol = 1.0e-10;
      double result[18];
      int tdims = 1;
      if(msh.isboundary_faces() && msh.nface > 0) tdims = 2;
      msh.tag[ithread]++;
      for(int ientt = 0; ientt < msh.nentt(tdims); ientt++){
        if(isdeadent(ientt, msh.ent2poi(tdims))) continue;
        for(int iver = 0; iver < tdims+1; iver++){
          int ipoin = msh.ent2poi(tdims)(ientt, iver);
          msh.poi2tag(ithread,ipoin) = msh.tag[ithread];
        }
      }

      for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
        int ipoin = msh.bpo2ibi(ibpoi,0);
        if(ipoin < 0) continue;
        if(msh.poi2tag(ithread,ipoin) < msh.tag[ithread]) continue;
        int tdim = msh.bpo2ibi(ibpoi,1);
        if(tdim == 0) continue;
        int ientt = msh.bpo2ibi(ibpoi,2);
        int iref = msh.ent2ref(tdim)[ientt];
        ego obj = tdim == 1 ? msh.CAD.cad2edg[iref] : msh.CAD.cad2fac[iref];
        EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
        double err = msh.idim == 2 ? geterrl2<2>(result, msh.coord[ipoin])
                                   : geterrl2<3>(result, msh.coord[ipoin]);
        if(err >= EG_tol){
          MPRINTF("## geo error {} EG_tol {} PRINT all bpoi for ipoin {} \n",
                 err,EG_tol,ipoin);
          for(int ibpo0 = msh.poi2bpo[ipoin]; ibpo0 >= 0; ibpo0 = msh.bpo2ibi(ibpo0,3)){
            int tdim1 = msh.bpo2ibi(ibpo0,1);
            int irefdbg = -1;
            double err1 = 0;
            if(tdim1 > 0){
              irefdbg = msh.ent2ref(tdim1)[msh.bpo2ibi(ibpo0,2)];
              ego obj1 = tdim1 == 1 ? msh.CAD.cad2edg[irefdbg] : msh.CAD.cad2fac[irefdbg];
              EG_evaluate(obj1, msh.bpo2rbi[ibpo0], result);
              err1 = msh.idim == 2 ? geterrl2<2>(result, msh.coord[ipoin])
                                   : geterrl2<3>(result, msh.coord[ipoin]);
            }

            MPRINTF("{} : {} {} {} ref {}; {} {}; err = {}\n",ibpo0, msh.bpo2ibi(ibpo0,0), 
              msh.bpo2ibi(ibpo0,1), msh.bpo2ibi(ibpo0,2), irefdbg,
              msh.bpo2rbi(ibpo0,0),msh.bpo2rbi(ibpo0,1), err1);
          }
        }
        METRIS_ENFORCE(err <= EG_tol);
      } 

      #if 0 
      dblAr1 rbpoi(msh.npoin);
      rbpoi.fill(1.0e30);
      for(int ientt = 0; ientt < msh.nentt(tdims); ientt++){
        if(isdeadent(ientt, msh.ent2poi(tdims))) continue;
        for(int iver = 0; iver < msh.nnode(tdims); iver++){
          int ipoi1 = msh.ent2poi(tdims)(ientt, iver);
          for(int ive2 = iver+1; ive2 < msh.nnode(tdims); ive2++){
            int ipoi2 = msh.ent2poi(tdims)(ientt, ive2);
            double dst = msh.idim == 2 ? geterrl2<2>(msh.coord[ipoi1], msh.coord[ipoi2])
                                       : geterrl2<3>(msh.coord[ipoi1], msh.coord[ipoi2]);
            if(dst < 1.0e-16){
              MPRINTF(" ## DEBUG < 1.0e-16 DISt = {} ipoi1 = {} ipoi2 = {}\n",dst,
                     ipoi1,ipoi2);
            }
            rbpoi[ipoi1] = MIN(rbpoi[ipoi1], dst);
            rbpoi[ipoi2] = MIN(rbpoi[ipoi2], dst);
          }
        }
      }

      bool ifail = false;
      for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
        int ipoin = msh.bpo2ibi(ibpoi,0);
        if(ipoin < 0) continue;
        int tdim = msh.bpo2ibi(ibpoi,1);
        if(tdim == 0) continue;
        int ientt = msh.bpo2ibi(ibpoi,2);
        int iref = msh.ent2ref(tdim)[ientt];
        ego obj = tdim == 1 ? msh.CAD.cad2edg[iref] : msh.CAD.cad2fac[iref];
        EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
        double err = msh.idim == 2 ? geterrl2<2>(result, msh.coord[ipoin])
                                   : geterrl2<3>(result, msh.coord[ipoin]);

        if(err >= rbpoi[ipoin]*EG_tol){
          MPRINTF("## geo error {} rbpoi {} EG_tol {} PRINT all bpoi for ipoin {} \n",
                 err,rbpoi[ipoin],EG_tol,ipoin);
          for(int ibpo0 = msh.poi2bpo[ipoin]; ibpo0 >= 0; ibpo0 = msh.bpo2ibi(ibpo0,3)){
            int tdim1 = msh.bpo2ibi(ibpo0,1);
            int irefdbg = -1;
            double err1 = 0;
            if(tdim1 > 0){
              irefdbg = msh.ent2ref(tdim1)[msh.bpo2ibi(ibpo0,2)];
              ego obj1 = tdim1 == 1 ? msh.CAD.cad2edg[irefdbg] : msh.CAD.cad2fac[irefdbg];
              EG_evaluate(obj1, msh.bpo2rbi[ibpo0], result);
              err1 = msh.idim == 2 ? geterrl2<2>(result, msh.coord[ipoin])
                                   : geterrl2<3>(result, msh.coord[ipoin]);
            }

            MPRINTF("{} : {} {} {} ref {}; {} {}; err = {}\n",ibpo0, msh.bpo2ibi(ibpo0,0), 
              msh.bpo2ibi(ibpo0,1), msh.bpo2ibi(ibpo0,2), irefdbg,
              msh.bpo2rbi(ibpo0,0),msh.bpo2rbi(ibpo0,1), err1);
          }
        }
        if(err >= rbpoi[ipoin]*EG_tol){
          std::cout<< "High point surface error "<<
            err<<" with eps "<<rbpoi[ipoin]<<" bpoi "<<ibpoi<<" is dim "<<tdim<<" ientt "<<ientt
            <<" iref "<<iref<<" rbi "<<msh.bpo2rbi(ibpoi,0)<<" "<<msh.bpo2rbi(ibpoi,1)
            <<" ipoin = "<<ipoin<<"\n";
          ifail = true;
        }
      }

      if(ifail){
        int npoi0 = msh.npoin;
        for(int ibpoi = 0; ibpoi < msh.nbpoi; ibpoi++){
          int ipoin = msh.bpo2ibi(ibpoi,0);
          if(ipoin < 0) continue;
          int tdim = msh.bpo2ibi(ibpoi,1);
          if(tdim == 0) continue;
          int ientt = msh.bpo2ibi(ibpoi,2);
          int iref = msh.ent2ref(tdim)[ientt];
          ego obj = tdim == 1 ? msh.CAD.cad2edg[iref] : msh.CAD.cad2fac[iref];
          EG_evaluate(obj, msh.bpo2rbi[ibpoi], result);
          int ipnew = msh.newpoitopo(tdim, -1);
          for(int ii = 0; ii < msh.idim; ii++) msh.coord(ipnew,ii) = result[ii];
          msh.newbpotopo(ipnew, 0, ipnew);
        }
        writeMesh("EGEvals",msh);
        for(int ipoin = npoi0; ipoin < msh.npoin; ipoin++) msh.killpoint(ipoin);
      }

      METRIS_ENFORCE_MSG(!ifail,"High surface error (wrong ts/(u,v)s)");
    #endif

    }

  }catch(const MetrisExcept& e){
    MPRINTF("Check_topo failed, dumping mesh:\n");
    writeMesh("check_topo_fail.meshb",msh);
    MPRINTF("DONE \n");
    throw(e);
  }
  //if(stopend){
  //  auto test_id = boost::unit_test::framework::current_test_case().p_id;
  //  METRIS_ENFORCE(boost::unit_test::results_collector.results(test_id).passed());
  //}

  CPRINTF2("-- END check_topo\n");
}








}// end namespace

