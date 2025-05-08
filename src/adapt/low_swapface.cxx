//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_swap.hxx"
#include "msh_swap.hxx" // for swapOptions

#include "../Mesh/Mesh.hxx"

#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../low_topo.hxx"
#include "../low_geo.hxx"
#include "../utils/aux_misc.hxx"
#include "../low_lenedg.hxx"
#include "../utils/mprintf.hxx"
#include "../low_normal.hxx"
#include "../io_libmeshb.hxx"
#include "../quality/low_metqua.hxx"

namespace Metris{


// Swap edge between two triangles (including surface w/ tets)
// Return 0 if nothing done, 1 if error, -1 if swap done
// Compute using norm specified in opt: if 0, take max.
// If norm is -1, use edge length instead.
template<class MFT, int gdim, int ideg>
int swapface(Mesh<MFT>& msh, int iface, swapOptions opt, 
             MshCavity &cav, CavWrkArrs &work, 
             double *qnrm0_, double *qnrm1_, int ithread){
  INCVDEPTH(msh.param);

  double &qnrm0 = *qnrm0_;
  double &qnrm1 = *qnrm1_;
  const int spnorm = opt.swap_norm; // swap norm -> < 0 is length, >= 0 is Lp over the 2

  if(spnorm > 0) METRIS_ASSERT(spnorm == msh.param->opt_pnorm);


  constexpr int tdim = 2;
  constexpr AsDeg asdmet = AsDeg::P1;


  if(isdeadent(iface,msh.fac2poi)) return 0; 

  CavOprOpt opts;
  CavOprInfo info;
  opts.allow_topological_correction = false;
  opts.skip_topo_checks = true;
  opts.allow_remove_points = false;
  opts.dryrun = false;
  cav.lcedg.set_n(0);
  cav.lcfac.set_n(0);
  cav.lctet.set_n(0);


  int iele0 = -1;
  if(msh.get_tdim() == 3){
    iele0 = msh.fac2tet(iface,0);
    if(iele0 < 0) iele0 = msh.fac2tet(iface,1);
    METRIS_ASSERT(iele0 >= 0);
  }


  double quae1;

  if(spnorm >= 0){
    quae1 = metqua<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,iface,1.0);
    METRIS_ASSERT_MSG(quae1 > -1.0e-16, "Negative quae1 "<<quae1<<" iface "<<iface);
  }


  CPRINTF1("-- START swapface iface = %d",iface);
  if(DOPRINTS1() && spnorm >= 0){
    printf(" initial quality = %f \n",quae1);
  }else if(DOPRINTS1()){
    printf("\n");
  }

  int iref = msh.fac2ref[iface];

  // Old qualities associated to each possible swap. 
  // If pnorm >= 0, this is the p-norm of quality accross 
  // Otherwise it is the length of the edge. 
  double quaol[3];
  for(int ied = 0; ied < 3; ied++){
    quaol[ied] = -1; // In bounds 0, 1, -1 is disregarded 

    int ifac2 = msh.fac2fac(iface,ied);
    if(ifac2 < 0) continue; // Can't swap across nm edge or bdry 

    // Different refs, skip 
    if(iref != msh.fac2ref[ifac2]) continue;

    // Note: manifold but edge in-between is ineig >= 0
    int iedge = msh.fac2edg(iface, ied);
    if(iedge >= 0) continue;


    if(spnorm >= 0){
      quaol[ied] = metqua<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,ifac2,1.0);
    }else{  
      double sz[2], len;
      len = getlenedg_geosz<MFT,gdim,ideg>(msh, iface, 2, ied, sz);
      // Attribute quality between 0 and 1, multiplicatively symmetric: 
      // i.e. q(sqrt2) = q(1/sqrt2). 
      quaol[ied] = len < 1.0 ? 1.0 - len 
                             : 1.0 - 1.0 / len;
      CPRINTF1(" - edge %d length %f quality %15.7e\n",ied,len, quaol[ied]);
    }
    CPRINTF1(" - candidate neighbour %d has qual %e\n",ifac2,quaol[ied]);
  }
  

  int idx[3] = {0,1,2};
  sortupto8_dec<double>(quaol,idx,3);

  // Simulate swaps as P1 
  // Improve when curvature added to cavity 
  //intAr2 fac2pol(2,3);
  int edg2pol[2]; // only for length based 
  //fac2pol.set_n(2);
  cav.lcfac.set_n(2); 
  cav.lcfac[0] = iface;
  for(int iix = 0; iix < 3; iix++){
    int ied = idx[iix];
    double quae2 = quaol[ied]; 
    if(quae2 < 0) continue;
    CPRINTF1(" - consider swap qface = %f qneigh = %f \n", quae1,quae2);

    // Quality of previous configuration 
    if(spnorm == 0){
      qnrm0 = MAX(quae1,quae2);
    }else if (spnorm > 0){
      qnrm0 = pow(pow(quae1,spnorm) + pow(quae2,spnorm), 1.0/spnorm);
    }else{
      qnrm0 = quae2;
    }

    int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
    int ip2 = msh.fac2poi(iface,lnoed2[ied][1]);

    int ifac2 = msh.fac2fac(iface,ied);
    METRIS_ASSERT(ifac2 >= 0);
    int ie2 = getedgfac(msh,ifac2,ip1,ip2);

    int nfac0 = msh.nface;
    msh.set_nface(msh.nface+2);

    msh.fac2ref[nfac0+0] = iref;
    msh.fac2ref[nfac0+1] = iref;

    msh.fac2poi(nfac0+0,0) = msh.fac2poi(iface,ied);
    msh.fac2poi(nfac0+0,1) = ip1;
    msh.fac2poi(nfac0+0,2) = msh.fac2poi(ifac2,ie2);

    msh.fac2poi(nfac0+1,0) = msh.fac2poi(iface,ied);
    msh.fac2poi(nfac0+1,1) = msh.fac2poi(ifac2,ie2);
    msh.fac2poi(nfac0+1,2) = ip2;

    bool iflat = false;
    for(int ii = 0; ii <= 1; ii++){
      getmeasentP1<gdim,2>(msh, msh.fac2poi[nfac0+ii], NULL, &iflat);
      if(iflat){
        CPRINTF1(" - new face %d: %d %d %d would be flat",ii+1,msh.fac2poi(nfac0+ii,0),
          msh.fac2poi(nfac0+ii,1),msh.fac2poi(nfac0+ii,2));
        break;
      }
    }

    if(iflat){
      msh.set_nface(nfac0);
      continue;
    }

    // Need to create new ibpoi entries as well for quality function nordev.
    if constexpr (gdim >= 3){

      int ibpon, ibpoo;
      // ip1
      ibpoo = msh.poi2ebp(ip1, 2, iface, iref);
      ibpon = msh.newbpotopo(ip1, 2, nfac0+0);
      for(int ii = 0; ii < nrbi; ii++) 
        msh.bpo2rbi(ibpon, ii) = msh.bpo2rbi(ibpoo, ii);
      // ip2
      ibpoo = msh.poi2ebp(ip2, 2, iface, iref);
      ibpon = msh.newbpotopo(ip2, 2, nfac0+1);
      for(int ii = 0; ii < nrbi; ii++) 
        msh.bpo2rbi(ibpon, ii) = msh.bpo2rbi(ibpoo, ii);
      // ip3 (msh.fac2poi(iface,ied))
      ibpoo = msh.poi2ebp(msh.fac2poi(iface,ied), 2, iface, iref);
      ibpon = msh.newbpotopo(msh.fac2poi(iface,ied), 2, nfac0+1);
      for(int ii = 0; ii < nrbi; ii++) 
        msh.bpo2rbi(ibpon, ii) = msh.bpo2rbi(ibpoo, ii);
      // ip4 (msh.fac2poi(ifac2,ie2))
      ibpoo = msh.poi2ebp(msh.fac2poi(ifac2,ie2), 2, ifac2, iref);
      ibpon = msh.newbpotopo(msh.fac2poi(ifac2,ie2), 2, nfac0+1);
      for(int ii = 0; ii < nrbi; ii++) 
        msh.bpo2rbi(ibpon, ii) = msh.bpo2rbi(ibpoo, ii);

      // For rembpotags. 
      msh.tag[ithread]++;
      msh.fac2tag(ithread,nfac0+0) = msh.tag[ithread];
      msh.fac2tag(ithread,nfac0+1) = msh.tag[ithread];

    }



    double qunw1, qunw2;
    if(spnorm >= 0){
      //qunw1 = metqua0<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,fac2pol[0],1.0);
      qunw1 = metqua<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,nfac0+0,1.0);
    }else{
      edg2pol[0] = msh.fac2poi(iface,ied);
      edg2pol[1] = msh.fac2poi(ifac2,ie2);
      double sz[2], len;
      len = getlenedg_geosz<MFT,gdim,ideg>(msh, edg2pol, sz);
      qunw1 = len < 1.0 ? 1.0 - len 
                        : 1.0 - 1.0 / len;
      CPRINTF1(" - new w/ edge %d length %f quality %15.7e\n",ied,len, qunw1);
    }
    CPRINTF1(" - new face quality = %f \n",qunw1);
    // Can skip already if using max
    if(spnorm == 0 && qunw1 + opt.swap_thres > qnrm0){
      msh.set_nface(nfac0);
      if constexpr (gdim >= 3){
        msh.rembpotag(ip1, ithread);
        msh.rembpotag(ip2, ithread);
        msh.rembpotag(msh.fac2poi(iface,ied), ithread);
        msh.rembpotag(msh.fac2poi(ifac2,ie2), ithread);
      }
      continue;
    }

    // If edge length, only one "quality" to consider. If worse, skip already. 
    if(spnorm  < 0 && qunw1 + opt.swap_thres > qnrm0){
      msh.set_nface(nfac0);
      if constexpr (gdim >= 3){
        msh.rembpotag(ip1, ithread);
        msh.rembpotag(ip2, ithread);
        msh.rembpotag(msh.fac2poi(iface,ied), ithread);
        msh.rembpotag(msh.fac2poi(ifac2,ie2), ithread);
      }
      continue;
    }

    if(spnorm >= 0){
      //qunw2 = metqua0<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,fac2pol[1],1.0);
      qunw2 = metqua<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,nfac0+1,1.0);
    }

    CPRINTF1(" - new face quality = %f \n",qunw2);
    
    msh.set_nface(nfac0);
    if constexpr (gdim >= 3){
      msh.rembpotag(ip1, ithread);
      msh.rembpotag(ip2, ithread);
      msh.rembpotag(msh.fac2poi(iface,ied), ithread);
      msh.rembpotag(msh.fac2poi(ifac2,ie2), ithread);
    }

    // Quality of new configuration 
    if(spnorm == 0){
      qnrm1 = MAX(qunw1,qunw2);
    }else if (spnorm > 0){
      qnrm1 = pow(pow(qunw1,spnorm) + pow(qunw2,spnorm), 1.0/spnorm);
    }else{
      qnrm1 = qunw1;
    }
    if(spnorm >= 0 && qnrm1 + opt.swap_thres > qnrm0) continue; 

    cav.lcfac[1] = ifac2;
    cav.ipins = msh.fac2poi(iface,ied);

    cav.lctet.set_n(0);
    if(msh.get_tdim() == 3){
      // Fill tet cavity with edge shell
      intAr1 dum;
      int iopen;
      shell(msh, ip1, ip2, 3, iele0, dum, dum, cav.lctet, &iopen);
    }

    CPRINTF1(" - enact swap ||(%f,%f)|| = %f -> ||(%f,%f)|| = %f \n ",
                                           quae1,quae2,qnrm0,qunw1,qunw2,qnrm1);
    if(tdim == 3) CPRINTF1(" - cavity ntetr %d \n",cav.lctet.get_n());

    int ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithread);
  
    if(info.done && ierro == 0){
      CPRINTF1("-- END swapface did %d - %d -> %d - %d \n",iface,
                                                 ifac2,msh.nface-2,msh.nface-1);
      return -1; // Return did op
    }
  }

  return 0;
}


#define BOOST_PP_LOCAL_MACRO(n)\
template int swapface<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical>& msh, \
                                    int iface, swapOptions opt, \
                                    MshCavity &cav, CavWrkArrs &work, \
                                    double *qumx0, double *qnrm1, int ithread);\
template int swapface<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        >& msh, \
                                    int iface, swapOptions opt, \
                                    MshCavity &cav, CavWrkArrs &work, \
                                    double *qumx0, double *qnrm1, int ithread);\
template int swapface<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical>& msh, \
                                    int iface, swapOptions opt, \
                                    MshCavity &cav, CavWrkArrs &work, \
                                    double *qumx0, double *qnrm1,int ithread);\
template int swapface<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        >& msh, \
                                    int iface, swapOptions opt, \
                                    MshCavity &cav, CavWrkArrs &work, \
                                    double *qumx0, double *qnrm1, int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

} // end namespace