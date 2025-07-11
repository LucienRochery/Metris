//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../adapt/msh_collapse.hxx"
#include "../adapt/low_collapse.hxx"
#include "../adapt/msh_swap.hxx"

#include "../Mesh/Mesh.hxx"

#include "../low_lenedg.hxx"
#include "../low_geo.hxx"
#include "../low_height.hxx"
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
template<class MFT, int gdim, int ideg>
double collapseShortEdges(Mesh<MFT> &msh, double qmax_suf, int *ncoll,
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

  const int tdim = msh.get_tdim();
  const int nedgl = (tdim*(tdim+1))/2;

  double stat = 0; 
  int ierro; 

  const bool ctrl_height = true;
  //const double isvolsmall = sqrt(3)/2 / 10;

  const int merror = CAV_ERR_NERROR;
  intAr1 lerro1(merror), lerro2(merror);

  msh.met.setSpace(MetSpace::Exp);


  const int miter = 10;
  int niter = 0;

  CPRINTF2("-- START collapseShortEdges miter = %d \n",miter);

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
    lerro1.fill(0);
    lerro2.fill(0);

    double minl = 1.0e30;
    double maxl = -1.0;

    double t0 = get_wall_time();



    double minht = 1.0e30;
    // Try collapsing for small height against bdry edge now
    // We don't have height control for surface yet (implement getheightentP1_aniso<gdim,tdim>)
    for(int ientt = 0; ientt < nent0 && ctrl_height && tdim == gdim; ientt++){
      INCVDEPTH(msh.param);
      if(ent2tag(ithrd1,ientt) >= msh.tag[ithrd1]) continue;
      if(isdeadent(ientt,ent2poi)) continue;


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
        if(!ibdry[iver]) continue;

        CPRINTF1(" - vertex %d height %e >=? h0 = %d\n",iver,height[iver],
                 height[iver] >= htmin);
        minht = MIN(minht, height[iver]);

        if(height[iver] >= htmin) continue;

        int ipcol = ent2poi(ientt,iver);
        int pdim = msh.getpoitdim(ipcol);
        if(pdim == 0) continue;

        ierro = collapseVertex(msh, ipcol, qmax_suf, cav, work, lerro2, ithrd2, ithrd3);
        if(ierro > 0){
          nerro2 ++;
        }else{
          ncoll2 ++;
        }
      }

    }// for ientt

    if(ctrl_height && tdim == gdim) CPRINTF1(" - min bdry height = %e\n",minht);



    // Collapse short edges 
    for(int ientt = 0; ientt < nent0; ientt++){
      INCVDEPTH(msh.param);
      CPRINTF2(" - debug ientt %d tag %d <? %d \n",ientt,ent2tag(ithrd1,ientt),msh.tag[ithrd1]);
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
          CPRINTF2(" - found short edge %d %d len = %f \n",
            ent2poi(ientt,lnoed2[ied][0]),ent2poi(ientt,lnoed2[ied][1]),len);
        }else{
          CPRINTF2(" - found short edge %d %d len = %f \n",
            ent2poi(ientt,lnoed3[ied][0]),ent2poi(ientt,lnoed3[ied][1]),len);
        }

        int nent00 = msh.nentt(tdim);
        ierro = collapseEdge(msh, tdim, ientt, ied, qmax_suf, cav, work, lerro1, ithrd2, ithrd3, ithrd4);


        if(ierro > 0){
          nerro1 ++;
          continue;
        }else{
          ncoll1 ++;
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



    double t1 = get_wall_time();
    int ncallps = 1000*(int)(((ncoll1+ncoll2) / (t1-t0)) / 1000);
    CPRINTF2(" - Loop end t = %f ncoll1 = %d ncoll2 = %d tot =  %d /s; nerro1 %d nerro2 %d\n",
      t1-t0,ncoll1,ncoll2,ncallps,nerro1,nerro2);
    CPRINTF2(" %f < len < %f \n",minl,maxl);
    if(DOPRINTS2()){
      if(nerro1 > 0){
        CPRINTF2(" - ierro list short edge:\n");
        for(int ii = 0; ii < merror; ii++){
          if(lerro1[ii] == 0) continue;
          CPRINTF2(" ierro = %d : %d \n",ii+1,lerro1[ii]);
        }
      }
      if(nerro2 > 0){
        CPRINTF2(" - ierro list low height bdry:\n");
        for(int ii = 0; ii < merror; ii++){
          if(lerro2[ii] == 0) continue;
          CPRINTF2(" ierro = %d : %d \n",ii+1,lerro2[ii]);
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
template double collapseShortEdges<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh,\
                           double qmax_suf, int* ncoll, int ithrd1, int ithrd2, int ithrd3, int ithrd4);\
template double collapseShortEdges<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh,\
                           double qmax_suf, int* ncoll, int ithrd1, int ithrd2, int ithrd3, int ithrd4);\
template double collapseShortEdges<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh,\
                           double qmax_suf, int* ncoll, int ithrd1, int ithrd2, int ithrd3, int ithrd4);\
template double collapseShortEdges<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh,\
                           double qmax_suf, int* ncoll, int ithrd1, int ithrd2, int ithrd3, int ithrd4);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()







} // end namespace