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

  constexpr int nnmet = (gdim*(gdim+1))/2;
  const int tdim = msh.get_tdim();
  const int nedgl = (tdim*(tdim+1))/2;

  double stat = 0; 
  int ierro; 

  const bool ctrl_height = false;
  // Bad:
  const bool ctrl_small_bdry = false;
  const double isvolsmall = sqrt(3)/2 / 10;

  const int merror = CAV_ERR_NERROR;
  intAr1 lerror(merror);

  msh.met.setSpace(MetSpace::Exp);

  //msh.met.setSpace(MetSpace::Log);
  if(ctrl_height || ctrl_small_bdry){
    METRIS_ENFORCE_MSG(msh.met.getSpace() == MetSpace::Log,
      "Front metric in log space: is it faster in collapse with options?");
  }

  const int miter = 10;
  int niter = 0;

  CPRINTF2("-- START collapseShortEdges miter = %d \n",miter);

  int ncoll1 = 0, ncoll2 = 0, ncoll3 = 0;
  *ncoll = 0;

  msh.tag[ithrd1]++;

  // More than 2 edges is a corner collapse
  MshCavity cav(100,100,2);
  CavWrkArrs work;

  intAr2& ent2tag = msh.ent2tag(tdim);
  intAr2& ent2poi = msh.ent2poi(tdim);
  intAr2& ent2ent = msh.ent2ent(tdim);



  do{
    INCVDEPTH(msh.param);


    int nerro1 = 0, nerro2 = 0, nerro3 = 0;
    int nedgt = 0;
    ncoll1 = ncoll2 = ncoll3 = 0;
    int nent0 = msh.nentt(tdim);
    lerror.fill(0);

    double minl = 1.0e30;
    double maxl = -1.0;

    double t0 = get_wall_time();


    // Don't do this for tdim == 3 (implement later)
    static int warn1 = 0;
    if(warn1++ < 10 && ctrl_height) 
      printf("## Disabled ctrl_height for tdim == 3\n");

    for(int ientt = 0; ientt < nent0 && ctrl_height && tdim == 2; ientt++){
      INCVDEPTH(msh.param);
      if(ent2tag(ithrd1,ientt) >= msh.tag[ithrd1]) continue;
      if(isdeadent(ientt,ent2poi)) continue;

      // If an operation goes through, the element goes away, then this does nothing
      // Otherwise, an operation does not happen, thus the element is inert.
      ent2tag(ithrd1,ientt) = msh.tag[ithrd1];


      // Try collapsing for small height against bdry edge now
      for(int ied = 0; ied < tdim + 1; ied++){

        int ipcol = ent2poi(ientt,ied);
        if(msh.poi2bpo[ipcol] >= 0){
          // Skip corners.
          if(msh.bpo2ibi[msh.poi2bpo[ipcol]][1] == 0) continue;
        }

        int ipoi1 = ent2poi(ientt,lnoed2[ied][0]);
        if(msh.poi2bpo[ipoi1] < 0) continue;
        int ipoi2 = ent2poi(ientt,lnoed2[ied][1]);
        if(msh.poi2bpo[ipoi2] < 0) continue;

        // Only if geometric edge 
        if(getedgglo(msh,ipoi1,ipoi2) < 0) continue;


        int ipoie = msh.newpoitopo(2,-1);

        double bary[3] = {};
        bary[lnoed2[ied][0]] = 0.5;
        bary[lnoed2[ied][1]] = 0.5;

        if(tdim == 2){
          eval2<gdim,ideg>(msh.coord,ent2poi[ientt],msh.getBasis(),
                           DifVar::None, DifVar::None, 
                           bary, msh.coord[ipoie], NULL, NULL);
        }else{
          if constexpr(gdim == 3){ // linker
            eval3<gdim,ideg>(msh.coord,ent2poi[ientt],msh.getBasis(),
                             DifVar::None, DifVar::None, 
                             bary, msh.coord[ipoie], NULL, NULL);
          }
        }


        msh.met.getMetBary(AsDeg::P1,DifVar::None, 
                           msh.met.getSpace(), ent2poi[ientt], 
                           tdim, bary, 
                           msh.met[ipoie], NULL) ;

        int edg2pol[2] = {ipoie, ent2poi(ientt,ied)};

        double sz[2];
        double len = getlenedg_geosz<MFT,gdim,ideg>(msh,edg2pol,sz);

        msh.set_npoin(msh.npoin-1);

        // Ideal height is sqrt(3)/2
        // sqrt(3)/2 / sqrt(2) is the smallest admissible height
        if(len > sqrt(3)/(2*sqrt(2))) continue;


        CPRINTF2(" - collapse flat %d height = %f \n",ientt,len);
        for(int ied2 = 0; ied2 < 3; ied2++){
          int edg2po2[2] = {ent2poi(ientt,lnoed2[ied][0]), 
                            ent2poi(ientt,lnoed2[ied][1])};
          double dd2s = getlenedg_geosz<MFT,gdim,ideg>(msh,edg2po2,sz);
          CPRINTF2(" - DEBUG ied = %d %d len = %f \n",ent2poi(ientt,lnoed2[ied][0]),
                    ent2poi(ientt,lnoed2[ied][1]), dd2s);
        }
        try{
          ierro = collversurf(msh, ientt, ied, qmax_suf, cav, work, lerror, ithrd2, ithrd3);
        }catch(const MetrisExcept &e){
          printf("## FATAL ERROR IN MSH_COLLAPSE\n");
          writeMesh("error_collapse.meshb",msh);
          throw(e);
        }
        if(ierro > 0){
          nerro2 ++;
        }else{
          ncoll2 ++;
        }

        break;
      }// for ied

    }// for ientt



    // Collapse small triangles (bad idea)
    static int warn2 = 0;
    if(warn2++ < 10 && ctrl_small_bdry) 
      printf("## Disabled ctrl_small_bdry for tdim == 3\n");
    for(int ientt = 0; ientt < nent0 && ctrl_small_bdry && tdim == 2; ientt++){
      INCVDEPTH(msh.param);
      if(ent2tag(ithrd1,ientt) >= msh.tag[ithrd1]) continue;
      if(isdeadent(ientt,ent2poi)) continue;

      // If an operation goes through, the element goes away, then this does nothing
      // Otherwise, an operation does not happen, thus the element is inert.
      ent2tag(ithrd1,ientt) = msh.tag[ithrd1];


      // Check at least one non boundary to collapse 
      bool iskip = true;
      int icol;
      for(int ied = 0; ied < 3; ied++){
        if(ent2ent(ientt,ied) >= 0) continue;
        int ipoin = ent2poi(ientt,ied);
        if(msh.poi2bpo[ipoin] >= 0) continue;

        iskip = false;
        icol  = ied;
      }
      if(iskip) continue;

      bool iflat;
      double bary[4];
      for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0 / (tdim + 1);
      double metl[nnmet];
      double detm, volM;
      double meas0 = tdim == 2 ? getmeasentP1<gdim,2>(msh, ent2poi[ientt], NULL, &iflat)
                               : getmeasentP1<gdim,3>(msh, ent2poi[ientt], NULL, &iflat);
      if(iflat) goto do_collapse;

      // If not flat, compute volume 

      msh.met.getMetBary(AsDeg::P1,DifVar::None, 
                         msh.met.getSpace(), ent2poi[ientt], 
                         tdim, bary, 
                         metl, NULL) ;
      
      detm = detsym<gdim>(metl);
      volM = meas0*sqrt(detm)/2/(tdim == 2 ? 1 : 3); // factorial tdim
      
      if(volM >= isvolsmall) continue;
      CPRINTF2(" - collapse small triangle %d vol = %f \n",ientt,volM);
      CPRINTF2(" - meas0 = %f detm = %f \n",meas0,detm);

      do_collapse:

      int nent00 = msh.nentt(tdim);
      try{
        ierro = collversurf(msh, ientt, icol, qmax_suf, cav, work, lerror, ithrd2, ithrd3);
      }catch(const MetrisExcept &e){
        printf("## FATAL ERROR IN MSH_COLLAPSE\n");
        writeMesh("error_collapse.meshb",msh);
        throw(e);
      }
      if(ierro > 0){
        nerro3 ++;
      }else{
        ncoll3 ++;
        for(int ientn = nent00; ientn < msh.nentt(tdim); ientn++){
          for(int ii = 0; ii < tdim + 1 ; ii++){
            int ineig = ent2ent(ientn,ii);
            if(ineig < 0) continue;
            METRIS_ASSERT(!isdeadent(ineig,ent2poi));
            ent2tag(ithrd1,ineig) = msh.tag[ithrd1] - 1;
          }
        }
      }
    }// for ientt



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
        try{
          ierro = colledgsurf(msh, tdim, ientt, ied, qmax_suf, cav, work, lerror, ithrd2, ithrd3, ithrd4);
        }catch(const MetrisExcept &e){
          printf("## FATAL ERROR IN MSH_COLLAPSE\n");
          writeMesh("error_collapse.meshb",msh);
          throw(e);
        }


        if(ierro > 0){
          nerro1 ++;
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
    CPRINTF2(" - Loop end t = %f ncoll1 = %d ncoll2 = %d ncoll3 = %d tot =  %d /s; nerro1 %d nerro2 %d nerro3 %d\n",
      t1-t0,ncoll1,ncoll2,ncoll3,ncallps,nerro1,nerro2,nerro3);
    CPRINTF2(" %f < len < %f \n",minl,maxl);
    if(DOPRINTS2()){
      if(nerro1 + nerro2 + nerro3 > 0){
        CPRINTF2(" - ierro list:\n");
        for(int ii = 0; ii < merror; ii++){
          if(lerror[ii] == 0) continue;
          CPRINTF2(" ierro = %d : %d \n",ii+1,lerror[ii]);
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