//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "low_swap2D.hxx"
#include "msh_swap2D.hxx" // for swapOptions

#include "../Mesh/Mesh.hxx"

#include "../cavity/msh_cavity.hxx"
#include "../aux_topo.hxx"
#include "../aux_utils.hxx"
#include "../low_lenedg.hxx"
#include "../mprintf.hxx"
#include "../low_geo.hxx"
#include "../io_libmeshb.hxx"
#include "../quality/low_metqua.hxx"

namespace Metris{



// Return 0 if done nothing, 1 if error, -1 if done swap
// Compute using norm specified in opt: if 0, take max. 
// If norm is -1, use edge length instead. 
template<class MFT, int gdim, int ideg>
int swapface(Mesh<MFT>& msh, int iface, swapOptions opt, 
             MshCavity &cav, CavWrkArrs &work, 
             double *qnrm0_, double *qnrm1_, int ithread){
  INCVDEPTH(msh);

  double &qnrm0 = *qnrm0_;
  double &qnrm1 = *qnrm1_;
  const int pnorm = opt.swap_norm;

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

  double quae1;

  if(pnorm >= 0){
    quae1 = metqua<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,iface,opts.qpower,
                                  opts.qpnorm,1.0);
    METRIS_ASSERT(quae1 > 0);
  }


  CPRINTF1("-- START swap2D iface = %d",iface);
  if(DOPRINTS1() && pnorm >= 0){
    printf(" initial quality = %f \n",quae1);
  }else if(DOPRINTS1()){
    printf("\n");
  }

  // Accept any swap that gives conformity error norm lower than this element's
  double nrmal[3];
  // Surface case, never tested
  if constexpr (gdim == 3){
    getnorfacP1(msh.fac2poi[iface],msh.coord,nrmal);
    cav.nrmal = nrmal;
  }

  // Old qualities associated to each possible swap. 
  // If pnorm >= 0, this is the p-norm of quality accross 
  // Otherwise it is the length of the edge. 
  double quaol[3];
  for(int ied = 0; ied < 3; ied++){
    quaol[ied] = -1; // In bounds 0, 1, -1 is disregarded 

    int ifac2 = msh.fac2fac(iface,ied);
    if(ifac2 < 0) continue; // Can't swap across nm edge or bdry 

    // Note: manifold but edge in-between is ineig >= 0
    int iedge = msh.fac2edg(iface, ied);
    if(iedge >= 0) continue;

    if(pnorm >= 0){
      quaol[ied] = metqua<MFT,gdim,tdim>(msh,AsDeg::P1, asdmet, 
                                         ifac2,opts.qpower,opts.qpnorm,1.0);
    }else{  
      double sz[2], len;
      len = getlenedg_geosz<MFT,gdim,ideg>(msh, iface, 2, ied, sz);
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
  cav.lcfac[0] = iface;
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

    int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
    int ip2 = msh.fac2poi(iface,lnoed2[ied][1]);

    int ifac2 = msh.fac2fac(iface,ied);
    METRIS_ASSERT(ifac2 >= 0);
    int ie2 = getedgfac(msh,ifac2,ip1,ip2);

    fac2pol(0,0) = msh.fac2poi(iface,ied);
    fac2pol(0,1) = ip1;
    fac2pol(0,2) = msh.fac2poi(ifac2,ie2);

    fac2pol(1,0) = msh.fac2poi(iface,ied);
    fac2pol(1,1) = msh.fac2poi(ifac2,ie2);
    fac2pol(1,2) = ip2;

    double qunw1, qunw2; 
    if(pnorm >= 0){
      qunw1 = metqua0<MFT,gdim,tdim>(msh,AsDeg::P1,asdmet,fac2pol[0],
                                     opts.qpower, opts.qpnorm,1.0);
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
    cav.ipins = msh.fac2poi(iface,ied);

    CPRINTF1(" - enact swap ||(%f,%f)|| = %f -> ||(%f,%f)|| = %f \n ",
                                           quae1,quae2,qnrm0,qunw1,qunw2,qnrm1);

    int ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithread);
  
    if(info.done && ierro == 0){
      CPRINTF1("-- END swap2D did %d - %d -> %d - %d \n",iface,
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


#if 0


// Return 0 if done nothing, 1 if error, -1 if done swap
template<class MFT>
int swapface(Mesh<MFT>& msh, int iface, int iverb, int ithread){
  int iret = 0;

  if(isdeadent(iface,msh.fac2poi)) return 0; 

  int mcfac = 2; // more than 2 is a corner collapse: no!
  const int nwork = 2;
  int iwork[nwork];
  MshCavity cav(0,mcfac,0,nwork,iwork);
  CavOprOpt opts;
  CavOprInfo info;
  CavWrkArrs work;
  opts.allow_topological_correction = false;
  opts.skip_topo_checks = true;
  opts.allow_remove_points = false;
  opts.dryrun = true;

  const int tdim = 2;

  double quael;
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    if(msh.idim == 2)
      metqua<MFT,2,tdim,ideg,AsDeg::Pk,double>(msh,iface,opts.qpower,
                                                         &quael,opts.qpnorm,1.0);
    else
      metqua<MFT,3,tdim,ideg,AsDeg::Pk,double>(msh,iface,opts.qpower,
                                                         &quael,opts.qpnorm,1.0);
  }}CT_FOR1(ideg);

  METRIS_ASSERT(quael > 0);

  CPRINTF1("-- START swap2D iface = %d initial quality = %f \n",iface,quael);
  #ifndef NDEBUG
    if(iverb >= 4) writeMesh("debug_swap0.meshb",msh);
  #endif

  // Accept any swap that gives maximum conformity error lower than this element's
  opts.qmax_suf = quael;


  cav.lcfac.set_n(2); 
  cav.lcfac[0] = iface;

  // In this case let's even just take the face normal
  double nrmal[3];
  if(msh.idim == 3){
    getnorfacP1(msh.fac2poi[iface],msh.coord,nrmal);
    cav.nrmal = nrmal;
  }


  for(int ied = 0; ied < 3; ied++){
    int ifac2 = msh.fac2fac(iface,ied);
    if(ifac2 < 0) continue; // Can't swap across geo edge

    METRIS_ASSERT(!isdeadent(ifac2,msh.fac2poi));

    cav.lcfac[1] = ifac2;

    cav.ipins = msh.fac2poi(iface,ied);

    int ip1 = msh.fac2poi(iface,lnoed2[ied][0]);
    int ip2 = msh.fac2poi(iface,lnoed2[ied][1]);
    int iedge = getedgglo(msh,ip1,ip2);

    if(iedge >= 0){
      continue;
    }


    double quae2;
    if(iverb >= 4){
      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
        if(msh.idim == 2)
          metqua<MFT,2,tdim,ideg,AsDeg::Pk,double>(msh,ifac2,
                                          opts.qpower,&quae2,opts.qpnorm,1.0);
        else
          metqua<MFT,3,tdim,ideg,AsDeg::Pk,double>(msh,ifac2,
                                          opts.qpower,&quae2,opts.qpnorm,1.0);
      }}CT_FOR1(ideg);

      printf(" - consider swap %d - %d quaels = %f , %f \n",iface,ifac2,quael,quae2);
    }

    int ierro;
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
      ierro = cavity_operator<MFT,ideg>(msh,cav,opts,work,info,ithread);
    }}CT_FOR1(ideg);


    if(info.done){
      CPRINTF1("-- END swap2D did %d - %d -> %d - %d \n",iface,
                                                  ifac2,msh.nface-2,msh.nface-1);
      #ifndef NDEBUG
        if(iverb >= 4) writeMesh("debug_swap1.meshb",msh);
      #endif
      CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
        if(msh.idim == 2){
          metqua<MFT,2,tdim,ideg,AsDeg::Pk,double>(msh,msh.nface-2,
                                            opts.qpower,&quael,opts.qpnorm,1.0);
          metqua<MFT,2,tdim,ideg,AsDeg::Pk,double>(msh,msh.nface-1,
                                            opts.qpower,&quae2,opts.qpnorm,1.0);
        }else{
          metqua<MFT,3,tdim,ideg,AsDeg::Pk,double>(msh,msh.nface-2,
                                            opts.qpower,&quael,opts.qpnorm,1.0);
          metqua<MFT,3,tdim,ideg,AsDeg::Pk,double>(msh,msh.nface-1,
                                            opts.qpower,&quae2,opts.qpnorm,1.0);
        }
      }}CT_FOR1(ideg);

      if(iverb >= 4){
        printf(" - new quaels = %f , %f \n",quael,quae2);
        //wait();
      }

      return -1; // Return did op
    }

  }


  return iret;
}


#endif

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