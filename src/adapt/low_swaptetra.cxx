//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_swap.hxx"
#include "msh_swap.hxx" // for swapOptions

#include "../Mesh/Mesh.hxx"

#include "../linalg/det.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../low_topo.hxx"
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
int swaptetra(Mesh<MFT>& msh, int itetr, swapOptions opt, 
             MshCavity &cav, CavWrkArrs &work, 
             double *qnrm0_, double *qnrm1_, int ithread){
  INCVDEPTH(msh.param);


  // If an edge is long, look for a edge to face swap (3->2, 4->4)
  // If a face is large, look for a face to edge swap (2->3). 

  // Edges can be on the boundary, only way to know is compute shell and check closed.
  // Thus start with faces

  // Compute metric determinants at the vertices: both for face area and edge len
  // computations. 
  double dets[4];
  for(int iver = 0; iver < 4; iver++){
    int ipoin = msh.tet2poi(itetr,iver);
    if(msh.met.getSpace() == MetSpace::Exp){
      dets[iver] = detsym<3,double>(msh.met[ipoin]);
    }else{
      dets[iver] = exp(msh.met(ipoin,0) + msh.met(ipoin,2) + msh.met(ipoin,5));
    }
  }
  for(int ifa = 0; ifa < 4; ifa++){
    double det = 0;
//    double bary[4] = 
  }

  bool irej[6] = {false, false, false, false, false, false};
  bool ifnd = false;

  int ilmax = -1;
  while(!ifnd){
    int imax = -1;
    double lmax = -1;
    double sz[2];
    for(int ied = 0; ied < 6; ied++){
      if(irej[ied]) continue;
      const int edg2pol[2] = {msh.tet2poi(itetr, lnoed3[ied][0]),
                              msh.tet2poi(itetr, lnoed3[ied][1])};
      double sz[2];
      double len = getlenedg_geosz<MFT,gdim,1>(msh,edg2pol,sz);
      if(len > lmax){
        lmax = len;
        imax = ied;
      }
    }
    if(lmax > 1){

    }
  }

  double &qnrm0 = *qnrm0_;
  double &qnrm1 = *qnrm1_;
  const int pnorm = opt.swap_norm;

  constexpr int tdim = 2;
  constexpr AsDeg asdmet = AsDeg::P1;


  if(isdeadent(itetr,msh.fac2poi)) return 0; 

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
    iele0 = msh.fac2tet(itetr,0);
    if(iele0 < 0) iele0 = msh.fac2tet(itetr,1);
    METRIS_ASSERT(iele0 >= 0);
  }


  double quae1;

  if(pnorm >= 0){
    quae1 = metqua<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,itetr,opts.qpower,
                                  opts.qpnorm,1.0);
    METRIS_ASSERT_MSG(quae1 > -1.0e-16, "Negative quae1 "<<quae1<<" itetr "<<itetr);
  }


  CPRINTF1("-- START swaptetra itetr = %d",itetr);
  if(DOPRINTS1() && pnorm >= 0){
    printf(" initial quality = %f \n",quae1);
  }else if(DOPRINTS1()){
    printf("\n");
  }

  // Old qualities associated to each possible swap. 
  // If pnorm >= 0, this is the p-norm of quality accross 
  // Otherwise it is the length of the edge. 
  double quaol[3];
  for(int ied = 0; ied < 3; ied++){
    quaol[ied] = -1; // In bounds 0, 1, -1 is disregarded 

    int ifac2 = msh.fac2fac(itetr,ied);
    if(ifac2 < 0) continue; // Can't swap across nm edge or bdry 

    // Note: manifold but edge in-between is ineig >= 0
    int iedge = msh.fac2edg(itetr, ied);
    if(iedge >= 0) continue;

    if(pnorm >= 0){
      quaol[ied] = metqua<MFT,gdim,tdim>(msh,AsDeg::P1, asdmet, 
                                         ifac2,opts.qpower,opts.qpnorm,1.0);
    }else{  
      double sz[2], len;
      len = getlenedg_geosz<MFT,gdim,ideg>(msh, itetr, 2, ied, sz);
      // Attribute quality between 0 and 1, multiplicatively symmetric: 
      // i.e. q(sqrt2) = q(1/sqrt2). 
      quaol[ied] = len < 1.0 ? 1.0 - len 
                             : 1.0 - 1.0 / len;
      CPRINTF1(" - edge %d length %f quality %15.7e\n",ied,len, quaol[ied]);
    }
  }
  

  int idx[3] = {0,1,2};
  sortupto8_dec<double>(quaol,idx,3);

  // Simulate swaps as P1 
  // Improve when curvature added to cavity 
  intAr2 fac2pol(2,3);
  int edg2pol[2]; // only for length based 
  fac2pol.set_n(2);
  cav.lcfac.set_n(2); 
  cav.lcfac[0] = itetr;
  for(int iix = 0; iix < 3; iix++){
    int ied = idx[iix];
    double quae2 = quaol[ied]; 
    if(quae2 < 0) continue;
    CPRINTF1(" - consider swap qface = %f qneigh = %f \n", quae1,quae2);

    // Quality of previous configuration 
    if(pnorm == 0){
      qnrm0 = MAX(quae1,quae2);
    }else if (pnorm > 0){
      qnrm0 = pow(pow(quae1,pnorm) + pow(quae2,pnorm), 1.0/pnorm);
    }else{
      qnrm0 = quae2;
    }

    int ip1 = msh.fac2poi(itetr,lnoed2[ied][0]);
    int ip2 = msh.fac2poi(itetr,lnoed2[ied][1]);

    int ifac2 = msh.fac2fac(itetr,ied);
    METRIS_ASSERT(ifac2 >= 0);
    int ie2 = getedgfac(msh,ifac2,ip1,ip2);

    fac2pol(0,0) = msh.fac2poi(itetr,ied);
    fac2pol(0,1) = ip1;
    fac2pol(0,2) = msh.fac2poi(ifac2,ie2);

    fac2pol(1,0) = msh.fac2poi(itetr,ied);
    fac2pol(1,1) = msh.fac2poi(ifac2,ie2);
    fac2pol(1,2) = ip2;

    double qunw1, qunw2; 
    if(pnorm >= 0){
      qunw1 = metqua0<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,fac2pol[0],
                                     opts.qpower, opts.qpnorm,1.0);
    }else{
      edg2pol[0] = msh.fac2poi(itetr,ied);
      edg2pol[1] = msh.fac2poi(ifac2,ie2);
      double sz[2], len;
      len = getlenedg_geosz<MFT,gdim,ideg>(msh, edg2pol, sz);
      qunw1 = len < 1.0 ? 1.0 - len 
                        : 1.0 - 1.0 / len;
      CPRINTF1(" - new w/ edge %d length %f quality %15.7e\n",ied,len, qunw1);
    }
    CPRINTF1(" - new face quality = %f \n",qunw1);
    // Can skip already if using max
    if(pnorm == 0 && qunw1 + opt.swap_thres > qnrm0) continue; 

    // If edge length, only one "quality" to consider. If worse, skip already. 
    if(pnorm  < 0 && qunw1 + opt.swap_thres > qnrm0) continue;

    if(pnorm >= 0){
      #ifndef NDEBUG
      try{
      #endif
      qunw2 = metqua0<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,fac2pol[1],
                                     opts.qpower, opts.qpnorm,1.0);
      #ifndef NDEBUG
      }catch(const MetrisExcept &e){
        printf(" ## METQUA FAILED DUE TO FAC2POL ? = \n");
        fac2pol.print();
        writeMesh("debugExcept",msh);
        METRIS_THROW(e);
      }
      #endif
    }

    CPRINTF1(" - new face quality = %f \n",qunw2);

    // Quality of new configuration 
    if(pnorm == 0){
      qnrm1 = MAX(qunw1,qunw2);
    }else if (pnorm > 0){
      qnrm1 = pow(pow(qunw1,pnorm) + pow(qunw2,pnorm), 1.0/pnorm);
    }else{
      qnrm1 = qunw1;
    }
    if(pnorm >= 0 && qnrm1 + opt.swap_thres > qnrm0) continue; 

    cav.lcfac[1] = ifac2;
    cav.ipins = msh.fac2poi(itetr,ied);

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
      CPRINTF1("-- END swaptetra did %d - %d -> %d - %d \n",itetr,
                                                 ifac2,msh.nface-2,msh.nface-1);
      //#ifndef NDEBUG
      //  if(iverb >= 4){
      //    writeMesh("debug_swap1.meshb",msh);
      //    for(int ifanw = msh.nface-2; ifanw < msh.nface; ifanw++){
      //      qunw1 = metqua<MFT,gdim,tdim,ideg,asdmet,double>(msh,ifanw,opts.qpower,
      //                                                          opts.qpnorm,1.0);
      //      printf(" - debug after cavity new face qua%d = %f \n",ifanw-msh.nface-1,qunw1);
      //    }
      //  } 
      //#endif
      return -1; // Return did op
    }
  }

  return 0;
}


#define BOOST_PP_LOCAL_MACRO(n)\
template int swaptetra<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical>& msh, \
                                    int itetr, swapOptions opt, \
                                    MshCavity &cav, CavWrkArrs &work, \
                                    double *qumx0, double *qnrm1,int ithread);\
template int swaptetra<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        >& msh, \
                                    int itetr, swapOptions opt, \
                                    MshCavity &cav, CavWrkArrs &work, \
                                    double *qumx0, double *qnrm1, int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

} // end namespace