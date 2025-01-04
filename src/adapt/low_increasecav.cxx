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
#include "../cavity/msh_cavity.hxx"
#include "../Mesh/Mesh.hxx"


namespace Metris{

int increase_cavity2D(MeshBase &msh, int iref, double *coop, 
                      const CavOprOpt &opts, MshCavity &cav, int ithread){

  const int iverb = msh.param->iverb;

  METRIS_ASSERT(cav.ipins >= 0 && cav.ipins < msh.npoin);


  msh.tag[ithread]++;
  for(int iface : cav.lcfac){
    METRIS_ASSERT(iface >= 0 && iface < msh.nface);
    METRIS_ASSERT(!isdeadent(iface,msh.fac2poi));

    msh.fac2tag(ithread,iface) = msh.tag[ithread];
  }

  for(int iedge : cav.lcedg){
    msh.edg2tag(ithread,iedge) = msh.tag[ithread];
  }

  if(iverb >= 4){
    printf("  -- start increase_cavity2D ipins %d list initial cavity:\n",
          cav.ipins);
    for(int iedge : cav.lcedg)
      printf("   - %d = %d %d \n",iedge, msh.edg2poi(iedge,0)
        , msh.edg2poi(iedge,1));
    for(int iface : cav.lcfac)
      printf("   - %d = %d %d %d \n",iface, msh.fac2poi(iface,0)
        , msh.fac2poi(iface,1), msh.fac2poi(iface,2));
  }

  int ibins = msh.poi2bpo[cav.ipins];
  int tdimi = -1;
  if(msh.isboundary_faces()){
    METRIS_ASSERT(ibins >= 0);
    tdimi = msh.bpo2ibi[ibins][1];
  }

  int fac2pol[3];
  fac2pol[0] = cav.ipins;

  // Note the bound is reeval'd
  for(int ifacl = 0; ifacl < cav.lcfac.get_n(); ifacl++){
    int iface = cav.lcfac[ifacl];
    if(iverb >= 4) printf("   - inccav try %d / %d = %d (%d,%d,%d) \n",
      ifacl,cav.lcfac.get_n(),iface,msh.fac2poi(iface,0)
      ,msh.fac2poi(iface,1),msh.fac2poi(iface,2));

    for(int iedl = 0; iedl < 3; iedl ++){
      int ifnei = msh.fac2fac(iface,iedl);

      int ip1 = msh.fac2poi(iface,lnoed2[iedl][0]);
      if(ip1 == cav.ipins) continue;
      int ip2 = msh.fac2poi(iface,lnoed2[iedl][1]);
      if(ip2 == cav.ipins) continue;


      if(iverb >= 4) printf("   - iedl %d ifnei = %d\n", iedl, ifnei);

      if(ifnei < 0){
        int iedge = getedgglo(msh,ip1,ip2);
        // cavity edge will be split -> no face from this 
        if(msh.edg2tag(ithread,iedge) >= msh.tag[ithread]) continue;
      }


      if(ifnei >= 0){
        if(msh.fac2tag(ithread,ifnei) >= msh.tag[ithread]){
          if(iverb >= 4) printf("   - ifnei = %d is tagged %d >= %d\n",
                             ifnei,msh.fac2tag(ithread,ifnei),msh.tag[ithread]);
          continue;
        }
      }


      fac2pol[1] = ip1;
      fac2pol[2] = ip2;
      
      bool iflat;
      double meas0;
      if(msh.idim == 2){
        meas0 = getmeasentP1<2,2>(fac2pol, msh.coord, msh.param->vtol, cav.nrmal, 
                                  &iflat,iverb-1);
        if(iverb >= 4)
          printf("    - inccav tested fac %d = %d %d %d w/ vtol = %f"
                 "got iflat = %d neighbour = %d\n", 
                 iface,ip1,ip2,cav.ipins,msh.param->vtol,iflat,ifnei);
      }else{
        METRIS_ASSERT(msh.idim == 3); // suppress iflat warning
        double nrmal[3]; 
        if(tdimi <= 1){ // NO depends on side for periodic
          getnorpoiref<1>(msh,cav.ipins,iref,cav.lcfac,nrmal);
          //msh.getpoinormal(cav.ipins,iref,nrmal);
        }else{
          for(int ii = 0; ii < 3; ii++) nrmal[ii] = cav.nrmal[ii];
        }
        meas0 = getmeasentP1<3,2>(fac2pol, msh.coord, msh.param->vtol, nrmal, 
                                  &iflat,iverb-1);
      }

      // ignore ifnei < 0 as it could be bdry -> edge remeshing
      if((iflat || meas0 < 0)){
        //if(ifnei == -1) return 1;
        // Cannot be corrected 
        if(ifnei < 0) return 1;
        cav.lcfac.stack(ifnei);
        msh.fac2tag(ithread,ifnei) = msh.tag[ithread];
        if(iverb >= 4) printf("   - inccav added %d to stack \n", ifnei);
      }
  
    }

  }

  return 0;
}


// Increase cavity for Delaunay criterion on ipoin 
template<class MFT>
void increase_cavity_Delaunay(MeshMetric<MFT> &msh, MshCavity &cav, 
                              int ipins,int ithread){

  //const int iverb = msh.param->iverb;
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
      for(int ii = 0; ii < nnmet; ii++) lmet[ii] = msh.met[ipins][ii];
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
        if(msh.poi2tag[ithrd1][ipoin] >= msh.tag[ithrd1]) continue;
        msh.poi2tag[ithrd1][ipoin] = msh.tag[ithrd1];

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
      if(msh.poi2tag[ithrd1][ipoin] >= msh.tag[ithrd1]) continue;
      msh.poi2tag[ithrd1][ipoin] = msh.tag[ithrd1];

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
