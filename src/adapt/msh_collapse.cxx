//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "../adapt/msh_collapse.hxx"
#include "../adapt/low_collapse.hxx"
#include "../adapt/msh_swap2D.hxx"

#include "../Mesh/Mesh.hxx"

#include "../low_lenedg.hxx"
#include "../low_geo.hxx"
#include "../aux_topo.hxx"
#include "../io_libmeshb.hxx"
#include "../aux_timer.hxx"
#include "../linalg/det.hxx"
#include "../cavity/msh_cavity.hxx"


namespace Metris{


// qmax_suf: quality threshold to accept collapses
// Exterior to this as may depend on swaps, inserts etc. 
// Prints level 1 routine 
template<class MFT, int gdim, int ideg>
double  collapseShortEdges(Mesh<MFT> &msh, double qmax_suf, int ithread ){
  constexpr int nnmet = (gdim*(gdim+1))/2;
  constexpr int tdim = 2;
  if(msh.get_tdim() != 2) METRIS_THROW_MSG(TODOExcept(), 
    "Implement collapseShortEdges on tdim != 2, got tdim = "<<msh.get_tdim());

  int iverb = msh.param->iverb;
  double stat = 0; 
  int ierro; 

  const bool ctrl_height = false;
  // Bad:
  const bool ctrl_small_bdry = false;
  const double isvolsmall = sqrt(3)/2 / 10;

  const int merror = CAV_ERR_NERROR;
  intAr1 lerror(merror);
  lerror.set_n(merror);


  msh.met.setSpace(MetSpace::Log);
  // Absolutely cretin first prototype:

  const int miter = 10;
  int niter = 0;

  if(iverb >= 1) printf("  -- START collapseShortEdges miter = %d \n",miter);

  int ncoll1 = 0, ncoll2 = 0, ncoll3 = 0;
  do{

    double t0 = get_wall_time();

    double stat_swap = swap2D<MFT,gdim,ideg>(msh, Defaults::swapOptAdapt);

    double t2 = get_wall_time();
    if(iverb >= 2){
      printf("   - Post collapse swap time %f stat = %f \n",t2-t0,stat_swap);
    }


    int nerro1 = 0, nerro2 = 0, nerro3 = 0;
    int nedgt = 0;
    ncoll1 = ncoll2 = ncoll3 = 0;
    int nfac0 = msh.nface;
    lerror.fill(0);

    double minl = 1.0e30;
    double maxl = -1.0;

    t0 = get_wall_time();



    for(int iface = 0; iface < nfac0 && ctrl_height; iface++){
      if(isdeadent(iface,msh.fac2poi)) continue;

      // Try collapsing for small height against bdry edge now
      for(int ied = 0; ied < 3; ied++){

        int ipcol = msh.fac2poi(iface,ied);
        if(msh.poi2bpo[ipcol] >= 0){
          // Skip corners.
          if(msh.bpo2ibi[msh.poi2bpo[ipcol]][1] == 0) continue;
        }

        int ipoi1 = msh.fac2poi(iface,lnoed2[ied][0]);
        if(msh.poi2bpo[ipoi1] < 0) continue;
        int ipoi2 = msh.fac2poi(iface,lnoed2[ied][1]);
        if(msh.poi2bpo[ipoi2] < 0) continue;

        // Only if geometric edge 
        if(getedgglo(msh,ipoi1,ipoi2) < 0) continue;


        int ipoie = msh.newpoitopo(2,-1);

        double bary[3] = {0};
        bary[lnoed2[ied][0]] = 0.5;
        bary[lnoed2[ied][1]] = 0.5;

        eval2<gdim,ideg>(msh.coord,msh.fac2poi[iface],msh.getBasis(),
                         DifVar::None, DifVar::None, 
                         bary, msh.coord[ipoie], NULL, NULL);


        msh.met.getMetBary(AsDeg::P1,DifVar::None, 
                           msh.met.getSpace(), msh.fac2poi[iface], 
                           2,  bary, 
                           msh.met[ipoie], NULL) ;

        int edg2pol[2] = {ipoie, msh.fac2poi(iface,ied)};

        double sz[2];
        double len = getlenedg_geosz<MFT,gdim,ideg>(msh,edg2pol,sz);

        msh.set_npoin(msh.npoin-1);

        // Ideal height is sqrt(3)/2
        // sqrt(3)/2 / sqrt(2) is the smallest admissible height
        if(len > sqrt(3)/(2*sqrt(2))) continue;


        if(iverb >= 3) printf("   - collapse flat %d height = %f \n",iface,len);
        for(int ied2 = 0; ied2 < 3; ied2++){
          int edg2po2[2] = {msh.fac2poi(iface,lnoed2[ied][0]), 
                            msh.fac2poi(iface,lnoed2[ied][1])};
          double dd2s = getlenedg_geosz<MFT,gdim,ideg>(msh,edg2po2,sz);
          printf(" DEBUG ied = %d %d len = %f \n",msh.fac2poi(iface,lnoed2[ied][0]),
            msh.fac2poi(iface,lnoed2[ied][1]), dd2s);
        }
        try{
          ierro = collversurf(msh, iface, ied, qmax_suf, lerror, ithread);
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
      }

    }



    // Collapse small triangles (bad idea)
    for(int iface = 0; iface < nfac0 && ctrl_small_bdry; iface++){
      if(isdeadent(iface,msh.fac2poi)) continue;

      // Check at least one non boundary to collapse 
      bool iskip = true;
      int icol;
      for(int ied = 0; ied < 3; ied++){
        if(msh.fac2fac(iface,ied) >= 0) continue;
        int ipoin = msh.fac2poi(iface,ied);
        if(msh.poi2bpo[ipoin] >= 0) continue;

        iskip = false;
        icol  = ied;
      }
      if(iskip) continue;

      bool iflat;
      double bary[3] = {1.0/3.0}; 
      double metl[nnmet];
      double detm, volM;
      double meas0 = getmeasentP1<gdim,tdim>(msh.fac2poi[iface],msh.coord,Defaults::vtol,
                                                                NULL,&iflat,0);
      if(iflat) goto do_collapse;

      // If not flat, compute volume 

      msh.met.getMetBary(AsDeg::P1,DifVar::None, 
                         msh.met.getSpace(), msh.fac2poi[iface], 
                         2,  bary, 
                         metl, NULL) ;
      
      detm = detsym<gdim>(metl);
      volM = meas0*sqrt(detm)/2; // factorial tdim
      
      if(volM >= isvolsmall) continue;
      if(iverb >= 3) printf("   - collapse small triangle %d vol = %f \n",iface,volM);
      printf("meas0 = %f detm = %f \n",meas0,detm);

      do_collapse:

      try{
        ierro = collversurf(msh, iface, icol, qmax_suf, lerror, ithread);
      }catch(const MetrisExcept &e){
        printf("## FATAL ERROR IN MSH_COLLAPSE\n");
        writeMesh("error_collapse.meshb",msh);
        throw(e);
      }
      if(ierro > 0){
        nerro3 ++;
      }else{
        ncoll3 ++;
      }


    }

    // Collapse short edges 
    for(int iface = 0; iface < nfac0; iface++){
      if(isdeadent(iface,msh.fac2poi)) continue;

      for(int ied = 0; ied < 3; ied++){

        // Consistency with insertion. Very close
        double len = getlenedg_geosz<MFT,gdim,ideg>(msh,iface,2,ied);
        minl = len < minl ? len : minl;
        maxl = len > maxl ? len : maxl;

        //int ipoi0 = msh.fac2poi(iface, lnoed2[ied][0]);
        //int ipoi1 = msh.fac2poi(iface, lnoed2[ied][1]);
        //int ityp0 = 2;
        //int ityp1 = 2;
        //int tmp = msh.poi2bpo[ipoi0];
        //if(tmp >= 0) ityp0 = msh.bpo2ibi(tmp,1);
        //tmp = msh.poi2bpo[ipoi1];
        //if(tmp >= 0) ityp1 = msh.bpo2ibi(tmp,1);


        nedgt++;

        if(len >= 1.0/sqrt(2)) continue;

        if(iverb >= 3) printf("   - found short edge %d %d len = %f \n",
          msh.fac2poi(iface,lnoed2[ied][0]),msh.fac2poi(iface,lnoed2[ied][1]),len);
        try{
          ierro = colledgsurf(msh, iface, ied, qmax_suf, lerror, ithread);
        }catch(const MetrisExcept &e){
          printf("## FATAL ERROR IN MSH_COLLAPSE\n");
          writeMesh("error_collapse.meshb",msh);
          throw(e);
        }
        if(ierro > 0){
          nerro1 ++;
        }else{
          ncoll1 ++;
          //if(ityp0 == 0 && ityp1 < 2|| ityp1 == 0 && ityp0 < 2){
          //  printf("Debug collapsed a boundary point against corner %d \n",ipoi0,ipoi1);
          //  printf("Length was %f \n",len);
          //  writeMesh("debug_corner", msh);
          //  wait();
          //}
        }

        break;
      }

    }



    double t1 = get_wall_time();
    int ncallps = 1000*(int)(((ncoll1+ncoll2) / (t1-t0)) / 1000);
    if(iverb >= 1){
      printf("   - Loop end t = %f ncoll1 %d = ncoll2 = %d ncoll3 = %d tot =  %d /s; nerro1 %d nerro2 %d nerro3 %d\n",
        t1-t0,ncoll1,ncoll2,ncoll3,ncallps,nerro1,nerro2,nerro3);
      printf("  %f < len < %f \n",minl,maxl);
      if(iverb >= 2){
        if(nerro1 + nerro2 + nerro3 > 0){
          printf("   - ierro list:\n");
          for(int ii = 0; ii < merror; ii++){
            if(lerror[ii] == 0) continue;
            printf("     ierro = %d : %d \n",ii+1,lerror[ii]);
          }
        }
      }
    }

    stat = MAX(stat, (double)(ncoll1 + ncoll2) / (double)nedgt);


    //if(iverb >= 1) printf(" - Warning: disabled collapse looping\n");
    //break;
  }while(ncoll1 + ncoll2 > 0 && niter++ < miter);


  return stat;

}


// ---------- Forward declarations
// See https://www.boost.org/doc/libs/1_82_0/libs/preprocessor/doc/AppendixA-AnIntroductiontoPreprocessorMetaprogramming.html
// Section A.4.1.2 Vertical Repetition
#define BOOST_PP_LOCAL_MACRO(n)\
template double collapseShortEdges<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh, double qmax_suf,int ithread);\
template double collapseShortEdges<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh, double qmax_suf,int ithread);\
template double collapseShortEdges<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh, double qmax_suf,int ithread);\
template double collapseShortEdges<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh, double qmax_suf,int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()





} // end namespace