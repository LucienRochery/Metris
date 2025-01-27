//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../MetrisRunner/MetrisParameters.hxx"
#include "../adapt/low_increasecav.hxx"
#include "../adapt/low_delaunay.hxx"
#include "../low_geo.hxx"
#include "../aux_topo.hxx"
#include "../low_lenedg.hxx"
#include "../low_topo.hxx"
#include "../mprintf.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../Mesh/Mesh.hxx"


namespace Metris{

int increase_cavity2D(MeshBase &msh, double *coop, 
                      const CavOprOpt &opts, MshCavity &cav, int ithread){
  GETVDEPTH(msh);

  METRIS_ASSERT(cav.ipins >= 0 && cav.ipins < msh.npoin);


  msh.tag[ithread]++;
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

  CPRINTF1("-- START increase_cavity2D ipins %d list initial cavity:\n", cav.ipins);
  if(DOPRINTS1()){
    if(cav.lcedg.get_n() > 0){
      CPRINTF1(" - Edge cavity: ");
      cav.lcedg.print(cav.lcedg.get_n());
    }
    if(cav.lcfac.get_n() > 0){
      CPRINTF1(" - Face cavity: ");
      cav.lcfac.print(cav.lcfac.get_n());
    }
    if(cav.lctet.get_n() > 0){
      CPRINTF1(" - Tetra cavity: ");
      cav.lctet.print(cav.lctet.get_n());
    }
  }
  if(DOPRINTS2()){
    for(int tdimn = 1; tdimn <= 3; tdimn++){
      intAr1 &lcent = cav.lcent(tdimn);
      int ncent = lcent.get_n();
      if(ncent <= 0) continue;
      intAr2 &ent2poi = msh.ent2poi(tdimn);

      if(tdimn == 1){
        CPRINTF2(" - Edge cavity: \n");
      }else if(tdimn == 2){
        CPRINTF2(" - Face cavity: \n");
      }else{
        CPRINTF2(" - Tetra cavity: \n");
      }
      int nnode = msh.nnode(tdimn);
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

  int fac2pol[3];
  fac2pol[0] = cav.ipins;

  // Note the bound is reeval'd, can't use range based
  for(int ifacl = 0; ifacl < cav.lcfac.get_n(); ifacl++){
    INCVDEPTH(msh)
    int iface = cav.lcfac[ifacl];
    CPRINTF1(" - inccav try %d / %d = %d (%d,%d,%d) \n",
                          ifacl,cav.lcfac.get_n(),iface,msh.fac2poi(iface,0)
                          ,msh.fac2poi(iface,1),msh.fac2poi(iface,2));

    // If dimension 3, get a normal for this face. 
    double norCAD[3];
    if(msh.idim == 3){
      getnorfacCAD(msh,iface,norCAD);
    }

    for(int iedl = 0; iedl < 3; iedl ++){

      int ip1 = msh.fac2poi(iface,lnoed2[iedl][0]);
      if(ip1 == cav.ipins) continue;
      int ip2 = msh.fac2poi(iface,lnoed2[iedl][1]);
      if(ip2 == cav.ipins) continue;

      int ifnei = msh.fac2fac(iface,iedl);

      CPRINTF1("   - iedl %d ifnei = %d\n", iedl, ifnei);

      if(ifnei >= 0){
        if(msh.fac2tag(ithread,ifnei) >= msh.tag[ithread]){
          CPRINTF1("   - ifnei = %d is tagged %d >= %d\n",
                             ifnei,msh.fac2tag(ithread,ifnei),msh.tag[ithread]);
          continue;
        }
      }

      // If there's an edge here and it's in the cavity, then it will be split 
      // and we'll get no face from it. 
      int iedge = msh.fac2edg(iface,iedl);
      if(iedge >= 0){
        if(msh.edg2tag(ithread,iedge) >= msh.tag[ithread]){
          CPRINTF1("   - iface %d -> iedge %d is tagged, skip\n",
                                iface,iedge);
          continue;
        }
      }

      // New face is ipins, ip1, ip2
      fac2pol[1] = ip1;
      fac2pol[2] = ip2;

      // First, check if this is a sliver
      bool iflat;
      double meas0;
      meas0 = msh.idim == 2 ? getmeasentP1<2,2>(msh, fac2pol, norCAD, &iflat)
                            : getmeasentP1<3,2>(msh, fac2pol, norCAD, &iflat);
      
      CPRINTF1("   - inccav pdim %d tested fac %d = %d %d %d w/ vtol = %f "
               "got iflat = %d meas0 = %15.7e neighbour = %d\n",pdim,
               iface,ip1,ip2,cav.ipins,msh.param->vtol,iflat,meas0,ifnei);

      #if 0
      // Next check geodev 
      // Actually not because adding more faces will only damage the cavity further
      // Do this in the future as pre reject, possibly. 
      // Also depends on Pk etc. Probably best to leave in cav.
      double nrmal[3]; 
      if(msh.idim == 3 && pdim < 2){
        // Get the normal in the case we're on an edge in 3D, and get only 
        // the correct side.
        int iref = msh.fac2ref[iface];
        getnorballref<1>(msh,cav.lcfac,iref,nrmal);
      }
      #endif
      // ignore ifnei < 0 as it could be bdry -> edge remeshing
      if((iflat || meas0 < 0)){
        //if(ifnei == -1) return 1;
        // Cannot be corrected 
        if(ifnei < 0){
          CPRINTF1(" # abort iflat %d meas0 %f \n",iflat,meas0);
          return 1;
        }
        cav.lcfac.stack(ifnei);
        msh.fac2tag(ithread,ifnei) = msh.tag[ithread];
        CPRINTF1("   - inccav added %d to stack \n", ifnei);
      }
  
    } // for int iedl

  } // for int ifacl

  return 0;
}


// Increase cavity for Delaunay criterion on ipoin 
template<class MFT>
void increase_cavity_Delaunay(MeshMetric<MFT> &msh, MshCavity &cav, 
                              int ipins,int ithread){


  int tdimn = msh.get_tdim();
  if(tdimn == 3) METRIS_THROW_MSG(TODOExcept(), "Unit test this for n = 3");
  int nnmet = (tdimn * (tdimn + 1)) / 2;


  const intAr2 &ent2ent = msh.ent2ent(tdimn);
  const intAr2 &ent2poi = msh.ent2poi(tdimn);
        intAr2 &ent2tag = msh.ent2tag(tdimn);
  intAr1 &lcent = cav.lcent(tdimn); 
  // NB: loop bounds MUST be reevaluated ! don't range-for this 
  msh.tag[ithread]++;
  for(int ientt : lcent){
    ent2tag(ithread,ientt) = msh.tag[ithread];
  }

  double metl[6], lmet[6];
  double *metl_p; 
  // If the metric field is log 
  if(msh.met.getSpace() == MetSpace::Log) metl_p = metl;

  for(int ii = 0; ii < lcent.get_n(); ii++){
    int ientt = lcent[ii];
    if(msh.met.getSpace() == MetSpace::Log){
      for(int ii = 0; ii < nnmet; ii++) lmet[ii] = msh.met(ipins,ii);
      if(tdimn == 2){
        getexpmet_cpy<2>(lmet, metl);
      }else{
        getexpmet_cpy<3>(lmet, metl);
      }
    }else{
      metl_p = msh.met[ipins];
    }
    for(int jj = 0; jj < tdimn + 1; jj++){
      int ienei = ent2ent(ientt,jj);
      if(ienei < 0) continue; // Non manifold also excluded
      if(ent2tag(ithread,ienei) >= msh.tag[ithread]) continue;
      ent2tag(ithread,ienei) = msh.tag[ithread];

      // Check if Delaunay 
      bool isinsph;
      if(tdimn == 2){
        isinsph = indelsphere<2>(msh.coord[ipins], metl_p, 
                                 msh.coord, ent2poi[ienei]);
      }else{
        isinsph = indelsphere<3>(msh.coord[ipins], metl_p, 
                                 msh.coord, ent2poi[ienei]);
      }
      if(isinsph) lcent.stack(ienei);
      
    }
  }

}

template void increase_cavity_Delaunay(MeshMetric<MetricFieldAnalytical> &msh, 
                            MshCavity &cav, int ipins,  int ithread);
template void increase_cavity_Delaunay(MeshMetric<MetricFieldFE        > &msh, 
                            MshCavity &cav, int ipins,  int ithread);





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
  //const int iverb = msh.param->iverb;
  constexpr int tdimn = gdim;
  METRIS_ASSERT(tdimn == msh.get_tdim());

  //const intAr2 &ent2ent = msh.ent2ent(tdimn);
  const intAr2 &ent2poi = msh.ent2poi(tdimn);
        intAr2 &ent2tag = msh.ent2tag(tdimn);
  intAr1 &lcent = cav.lcent(tdimn); 
  // NB: loop bounds MUST be reevaluated ! don't range-for this 
  msh.tag[ithrd1]++;
  for(int ientt : lcent){
    ent2tag(ithrd1,ientt) = msh.tag[ithrd1];
  }


  const int nedgl = (tdimn*(tdimn+1))/2;
  const intAr2 lnoed(nedgl,2,tdimn == 2 ? lnoed2[0] : lnoed3[0]);

  intAr1 lball(20);
  int iopen;
  bool imani;
  intAr1 dum;

  int nprem = 0;

  int edg2pol[2];
  edg2pol[0] = ipins;
  double sz[2];

  //int ncomp = 0;
  //int ncav0 = lcent.get_n();

  for(int ii = 0; ii < lcent.get_n(); ii++){
    int ientt = lcent[ii];


    #if 0
    for(int ifa = 0; ifa < tdimn + 1; ifa++){
      int ientn = ent2ent(ientt,ifa);
      if(ientn >= 0){
        if(ent2tag(ithrd1,ientn) >= msh.tag[ithrd1]) continue;
      }
      // Cavity boundary 
      // Loop over face nodes 
      int kk = -1;
      for(int ii = 0; ii < tdimn; ii++){
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
        if constexpr (tdimn == 2){
          ball2(msh,ipoin,ientt,lball,dum,&iopen,&imani,ithrd2);
        }else{
          ball3(msh,ipoin,ientt,lball,&iopen,ithrd2);
        }
        nprem++;
        for(int ient2 : lball){
          if(ent2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          ent2tag(ithrd1,ient2) = msh.tag[ithrd1];
          lcent.stack(ient2);
        }
      }
    }
    #else
    for(int inode = 0; inode < tdimn + 1; inode++){
      int ipoin = ent2poi(ientt,inode);
      if(ipoin == ipins) continue;
      if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
      msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1];

      edg2pol[1] = ipoin;
      double len = getlenedg_geosz<MFT,gdim,1>(msh,edg2pol,sz);
      //ncomp++;

      if(len <= 1.0/sqrt(2)){
        if(!opts.allow_remove_points) return -1; 
        if constexpr (tdimn == 2){
          ball2(msh,ipoin,ientt,lball,dum,&iopen,&imani,ithrd2);
        }else{
          ball3(msh,ipoin,ientt,lball,&iopen,ithrd2);
        }
        nprem++;
        for(int ient2 : lball){
          if(ent2tag(ithrd1,ient2) >= msh.tag[ithrd1]) continue;
          ent2tag(ithrd1,ient2) = msh.tag[ithrd1];
          lcent.stack(ient2);
        }
      }
    }
    #endif

    // Control height, only in dimension 2d.
    if constexpr(tdimn == 2){

    }else{
      METRIS_THROW_MSG(TODOExcept(), 
        "Implement height control in increase_cavity_lenedg 3D");
    }
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


} // end namespace
