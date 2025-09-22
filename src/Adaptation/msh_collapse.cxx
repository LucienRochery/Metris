//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_collapse.hxx"
#include "low_collapse.hxx"
#include "msh_swap.hxx"
#include "Insertion/low_insert.hxx"

#include "../Mesh/Mesh.hxx"

#include "../low_geo/lenedg.hxx"
#include "../low_geo/misc.hxx"
#include "../low_geo/height.hxx"
#include "../aux_topo.hxx"
#include "../io_libmeshb.hxx"
#include "../utils/aux_timer.hxx"
#include "../utils/mprintf.hxx"
#include "../linalg/det.hxx"
#include "../cavity/msh_cavity.hxx"
#include "../aux_histogram.hxx"
#include "../msh_lenedg.hxx"


namespace Metris{


// qmax_suf: quality threshold to accept collapses
// Exterior to this as may depend on swaps, inserts etc. 
// Prints level 1 routine 
// ithrdcstr tracks constrained points
template<class MFT, int gdim, int ideg>
double collapseShortEdges(Mesh<MFT> &msh, int tdim, double qmax_suf, int *ncoll,
                          int ithrd1, int ithrd2, int ithrd3, int ithrd4){

  GETVDEPTH(msh.param);

  METRIS_ASSERT(ithrd1 >= 0 && ithrd1 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd2 >= 0 && ithrd2 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd3 >= 0 && ithrd3 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd4 >= 0 && ithrd4 < METRIS_MAXTAGS);
  METRIS_ASSERT(ithrd1 != ithrd2);
  METRIS_ASSERT(ithrd1 != ithrd3);
  METRIS_ASSERT(ithrd2 != ithrd3);
  METRIS_ASSERT(ithrd4 != ithrd1);
  METRIS_ASSERT(ithrd4 != ithrd2);
  METRIS_ASSERT(ithrd4 != ithrd3);

  const int nedgl = (tdim*(tdim+1))/2;
  const intAr2 lnoed(nedgl,2,tdim == 1 ? lnoed1[0] :
                             tdim == 2 ? lnoed2[0] : lnoed3[0]);

  double stat = 0; 
  int ierro; 

  const bool ctrl_height = true && tdim >= 2;
  //const double isvolsmall = sqrt(3)/2 / 10;

  const int mcaverr = CAV_ERR_NERROR;
  intAr1 lcaver1(mcaverr), lcaver2(mcaverr);
  const int mcolerr = INS2D_ERR_NERROR;
  intAr1 lcoler1(mcolerr), lcoler2(mcolerr);

  msh.met.setSpace(MetSpace::Exp);


  const int miter = 10;
  int niter = 0;

  CPRINTF2("-- START collapseShortEdges miter = {} \n",miter);

  int ncoll1 = 0, ncoll2 = 0;
  *ncoll = 0;

  msh.tag[ithrd1]++;

  // More than 2 edges is a corner collapse
  MshCavity cav(100,100,2);
  CavWrkArrs work;

  intAr2& ent2tag = msh.ent2tag(tdim);
  intAr2& ent2poi = msh.ent2poi(tdim);
  intAr2& ent2ent = msh.ent2ent(tdim);


  // 2D: Ideal height is sqrt(3)/2
  // 3D: Ideal height is sqrt(2/3)
  // h0 / sqrt(2) is the smallest admissible height
  const double htmin = tdim == 2 ? sqrt(3.0/2)/2 : sqrt(1.0/3);


  do{
    INCVDEPTH(msh.param);

    int nerro1 = 0, nerro2 = 0;
    int nedgt = 0;
    ncoll1 = ncoll2 = 0;
    int nent0 = msh.nentt(tdim);
    lcoler1.fill(0);
    lcaver1.fill(0);
    lcoler2.fill(0);
    lcaver2.fill(0);

    double minl = 1.0e30;
    double maxl = -1.0;

    double t0 = get_cpu_time();



    double minht = 1.0e30;
    // Try collapsing for small height against bdry edge now
    // We don't have height control for surface yet (implement getheightentP1_aniso<gdim,tdim>)
    for(int ientt = 0; ientt < nent0 && ctrl_height && tdim == gdim; ientt++){
      INCVDEPTH(msh.param);
      if(ent2tag(ithrd1,ientt) >= msh.tag[ithrd1]) continue;
      if(isdeadent(ientt,ent2poi)) continue;

      CPRINTF1("-- height check ientt {} dim {}\n",ientt,tdim);


      // Determine which pts are opposite bdry facets, if any
      bool skipelt = true;
      bool ibdry[4];
      if(tdim == 3){
        if(msh.is_nonmanifold()){
          for(int ifa = 0; ifa < 4; ifa++){
            ibdry[ifa] = false;
            if(msh.tet2fac(ientt,ifa) < 0) continue;
            ibdry[ifa] = true;
            skipelt = false;
          }
        }else{
          for(int ifa = 0; ifa < 4; ifa++){
            ibdry[ifa] = false;
            if(msh.tet2tet(ientt,ifa) >= 0) continue;
            ibdry[ifa] = true;
            skipelt = false;
          }
        }
      }else{
        for(int ied = 0; ied < 3; ied++){
          ibdry[ied] = false;
          if(msh.fac2edg(ientt,ied) < 0) continue;
          ibdry[ied] = true;
          skipelt = false;
        }
      }

      if(skipelt) continue;

      // Get heights
      double height[4];
      getheightentP1_aniso<MFT,gdim>(msh, ientt, height);

      for(int iver = 0; iver < tdim + 1; iver++){
        INCVDEPTH(msh.param);
        if(!ibdry[iver]) continue;

        CPRINTF1(" - vertex {} height {} >=? h0 = {}\n",iver,height[iver],
                 height[iver] >= htmin);
        minht = MIN(minht, height[iver]);

        if(height[iver] >= htmin) continue;

        int ipcol = ent2poi(ientt,iver);
        int pdim = msh.getpoitdim(ipcol);
        if(pdim == 0) continue;

        CPRINTF1(" - call collapseVertex ipcol {} pdim {}\n",ipcol,pdim);

        ierro = collapseVertex(msh, ipcol, qmax_suf, cav, work, lcaver2, ithrd2, ithrd3);
        if(ierro > 0){
          nerro2++;
        }else{
          ncoll2++;
        }
      }

    }// for ientt

    if(ctrl_height && tdim == gdim) CPRINTF1(" - min bdry height = {}\n",minht);



    // Collapse short edges 
    for(int ientt = 0; ientt < nent0; ientt++){
      INCVDEPTH(msh.param);
      CPRINTF2(" - debug ientt {} tag {} <? {} \n",ientt,ent2tag(ithrd1,ientt),msh.tag[ithrd1]);
      if(ent2tag(ithrd1,ientt) >= msh.tag[ithrd1]) continue;
      if(isdeadent(ientt,ent2poi)) continue;


      // If an operation goes through, the element goes away, then this does nothing
      // Otherwise, an operation does not happen, thus the element is inert.
      ent2tag(ithrd1,ientt) = msh.tag[ithrd1];


      for(int ied = 0; ied < nedgl; ied++){

        // Consistency with insertion. Very close
        double len = getlenedg_geosz<MFT,gdim,ideg>(msh,ientt,tdim,ied);
        minl = len < minl ? len : minl;
        maxl = len > maxl ? len : maxl;

        //int ipoi0 = ent2poi(ientt, lnoed2[ied][0]);
        //int ipoi1 = ent2poi(ientt, lnoed2[ied][1]);
        //int ityp0 = 2;
        //int ityp1 = 2;
        //int tmp = msh.poi2bpo[ipoi0];
        //if(tmp >= 0) ityp0 = msh.bpo2ibi(tmp,1);
        //tmp = msh.poi2bpo[ipoi1];
        //if(tmp >= 0) ityp1 = msh.bpo2ibi(tmp,1);


        nedgt++;

        if(len >= 1.0/sqrt(2)) continue;  

        if(tdim == 2){
          CPRINTF2(" - found short edge {} {} len = {} \n",
            ent2poi(ientt,lnoed(ied,0)),ent2poi(ientt,lnoed(ied,1)),len);
        }else{
          CPRINTF2(" - found short edge {} {} len = {} \n",
            ent2poi(ientt,lnoed(ied,0)),ent2poi(ientt,lnoed(ied,1)),len);
        }

        int nent00 = msh.nentt(tdim);

        int icorner = -1;
        if(msh.getpoitdim(ent2poi(ientt,lnoed(ied,0))) == 0) icorner = 0;
        if(msh.getpoitdim(ent2poi(ientt,lnoed(ied,1))) == 0) icorner = 1;
        
        if(icorner < 0){
          ierro = collapseEdge(msh, tdim, ientt, ied, qmax_suf, cav, work, lcaver1, ithrd2, ithrd3, ithrd4);
        }else{
          CPRINTF1(" # collapseEdge vertex {} = {} is corner, call collapseVertex\n",
                   icorner,ent2poi(ientt,lnoed(ied,icorner)));
          int ipcol = ent2poi(ientt,lnoed(ied,1-icorner));
          ierro = collapseVertex(msh, ipcol, qmax_suf, cav, work, lcaver1, ithrd2, ithrd3);
        }


        if(ierro > 0){
          lcoler1[ierro - 1]++;
          nerro1++;
          //fmt::print("# DEBUG IERRO = {}", ierro);
          //wait();
          continue;
        }else{
          ncoll1++;
          for(int ientn = nent00; ientn < msh.nentt(tdim); ientn++){
            for(int ii = 0; ii < tdim + 1 ; ii++){
              int ineig = ent2ent(ientn,ii);
              if(ineig < 0) continue;
              METRIS_ASSERT(!isdeadent(ineig,ent2poi));
              ent2tag(ithrd1,ineig) = msh.tag[ithrd1] - 1;
            }
          }
        }

        break;
      }// for ied
    }// for ientt



    double t1 = get_cpu_time();
    int ncallps = 1000*(int)(((ncoll1+ncoll2) / (t1-t0)) / 1000);
    CPRINTF2(" - Loop end t = {:.2e} ncoll1 = {} ncoll2 = {} tot =  {} /s; nerro1 {} nerro2 {}\n",
      t1-t0,ncoll1,ncoll2,ncallps,nerro1,nerro2);
    CPRINTF2(" {} < len < {} \n",minl,maxl);
    if(DOPRINTS2()){
      if(nerro1 > 0){
        CPRINTF2(" - short edge cavity ierro list:\n");
        for(int ii = 0; ii < mcaverr; ii++){
          if(lcaver1[ii] == 0) continue;
          CPRINTF2("   ierro = {} : {} \n",ii+1,lcaver1[ii]);
        }
        CPRINTF2(" - short edge inspoi ierro list:\n");
        for(int ii = 0; ii < mcolerr; ii++){
          if(lcoler1[ii] == 0) continue;
          CPRINTF2("   ierro = {} : {} \n",ii+1,lcoler1[ii]);
        }
      }
      if(nerro2 > 0){
        CPRINTF2(" - low height cavity ierro list:\n");
        for(int ii = 0; ii < mcaverr; ii++){
          if(lcaver2[ii] == 0) continue;
          CPRINTF2("   ierro = {} : {} \n",ii+1,lcaver2[ii]);
        }
        CPRINTF2(" - low height inspoi ierro list:\n");
        for(int ii = 0; ii < mcolerr; ii++){
          if(lcoler2[ii] == 0) continue;
          CPRINTF2("   ierro = {} : {} \n",ii+1,lcoler2[ii]);
        }
      }
    }

    if(nedgt == 0) stat = MAX(stat, 0);
    else           stat = MAX(stat, (double)(ncoll1 + ncoll2) / (double)nedgt);

    *ncoll += ncoll1 + ncoll2;

  }while(ncoll1 + ncoll2 > 0 && niter++ < miter);

  //if(DOPRINTS2()){
  //  intAr2 ilned;
  //  dblAr1 rlned;
  //  dblAr1 lenbds = {1.0/sqrt(2), sqrt(2)};
  //  double pct_unit = getLengthEdges<MFT>(msh,ilned,rlned);
  //  print_histogram(msh,rlned,IntrpTyp::Linear,lenbds,"l","Edge length (collapse)");
  //}

  return stat;

}


// ---------- Forward declarations
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double collapseShortEdges<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh,int tdim,\
                           double qmax_suf, int* ncoll, int ithrd1, int ithrd2, int ithrd3, int ithrd4);\
template double collapseShortEdges<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh,int tdim,\
                           double qmax_suf, int* ncoll, int ithrd1, int ithrd2, int ithrd3, int ithrd4);\
template double collapseShortEdges<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh,int tdim,\
                           double qmax_suf, int* ncoll, int ithrd1, int ithrd2, int ithrd3, int ithrd4);\
template double collapseShortEdges<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh,int tdim,\
                           double qmax_suf, int* ncoll, int ithrd1, int ithrd2, int ithrd3, int ithrd4);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()







} // end namespace