//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

/*
Routine for "direct" smoothing as P1. From each (facet, metric) pair, generate remaining vertex to be unit. Then average over ball.
Simplest possible approach.
*/

#include "../Mesh/Mesh.hxx" 
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../smoothing/msh_smooball.hxx"
#include "../smoothing/low_smooballdiff.hxx"

#include "../aux_topo.hxx"
#include "../aux_timer.hxx"
#include "../low_topo.hxx"
#include "../mprintf.hxx"
#include "../quality/low_metqua.hxx"
#include "../io_libmeshb.hxx"

//#include "../libs/lplib3.h"


namespace Metris{

template<class MFT>
double smoothInterior_Ball(Mesh<MFT> &msh, QuaFun iquaf, int ithrd1, int ithrd2){

  int tdimn = msh.get_tdim();

  if(tdimn == 1) METRIS_THROW(TODOExcept());

  // Geo and topo dimn must match otherwise surface specific 
  METRIS_ASSERT(tdimn == msh.idim);
  double noper;
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    if(tdimn == 2){
      noper = smoothInterior_Ball0<MFT,2,ideg>(msh,iquaf,ithrd1,ithrd2);
    }else{
      noper = smoothInterior_Ball0<MFT,3,ideg>(msh,iquaf,ithrd1,ithrd2);
    } 
  }}CT_FOR1(ideg);

  return noper; 
}

template double smoothInterior_Ball<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, 
                                             QuaFun iquaf, int ithrd1, int ithrd2);
template double smoothInterior_Ball<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, 
                                             QuaFun iquaf, int ithrd1, int ithrd2);




// idim: gdim = tdim
template<class MFT, int idim, int ideg>
double smoothInterior_Ball0(Mesh<MFT> &msh, QuaFun iquaf, 
                            int ithrd1, int ithrd2){
  GETVDEPTH(msh);

  constexpr int tdim = idim;
  //constexpr int gdim = idim;
  if(ideg > idim + 1){
    printf("## SMOOTHING DISABLED FOR DEGREE %d and dim %d \n",ideg,idim);
    return -1.0;
  }


  int nentt = msh.nentt(tdim);
  const intAr2 &ent2poi = msh.ent2poi(tdim); 
  const intAr2 &ent2ent = msh.ent2ent(tdim); 

  #undef USE_LPLIB_SMOOTHINTERIOR 

  #ifdef USE_LPLIB_SMOOTHINTERIOR
    // LPlib init
    int nproc = msh.param->nproc;
    int nthread = GetNumberOfCores();
    if(nthread <= 0){
      CPRINTF1("## WARNING: LPlib function GetNumberOfCores() returned " 
               "negative threads. Set to default %d.\n",METRIS_MAXTAGS);
      nthread = METRIS_MAXTAGS;
    }else{
      CPRINTF2("-- LPlib found ncore = %d \n",nthread);
      if(nthread > METRIS_MAXTAGS){
        CPRINTF1("## WARNING: must verify nthread <= METRIS_MAXTAGS = %d." 
                " Increase in metris_constants.hxx.\n", METRIS_MAXTAGS);
        nthread = METRIS_MAXTAGS;
      }
    }
    if(nproc > 0) nthread = MIN(nthread, nproc);
    int64_t LibIdx = InitParallel(nthread);
    int LP_elt = NewType(LibIdx, nentt);
    int LP_poi = NewType(LibIdx, msh.npoin);
    float LP_stat[2];
    BeginDependency(LibIdx, LP_elt, LP_poi);
    for(int ientt = 0; ientt < nentt; ientt++){
      for(int ii = 0; ii < tdim + 1; ii++) 
        AddDependency(LibIdx, ientt+1, ent2poi(ientt,ii)+1);
    }
    EndDependency(LibIdx, LP_stat);
    // END LPlib init

    int itag_shared = ithrd1;
    // ithrd2 can be used freely as it is for elements, whose collisions are
    // avoided by LPlib
  #endif

  msh.met.setSpace(MetSpace::Log);



  // Eventually move all constants to MetrisParameters 
  // L2 conformity error from 0 to 1
  const int qpower = msh.param->opt_power;
  const int qpnorm = msh.param->opt_pnorm; 
  const double difto = 1.0;
  const int miter = 10;
  //const double maxwt = 20.0;
  //const double qrthr = 2.0;
  const double tolavg = 0.005;
  const double tolmax = 0.005;

  dblAr1 work;
  if(msh.param->iflag2 != 0){
    MPRINTF("\n\n##WARNING EXPERIMENTAL SMOOTHING FUNCTION 2\n");
  }

  // 1 -> no maximum quality increase allowed 
  //const double maxinc_worst = 1.00;

  constexpr int nnmet = (idim*(idim+1))/2;

  METRIS_ENFORCE(qpower < 0); // Otherwise rework the mins / maxs
  // Otherwise not only edge nodes 
  METRIS_ENFORCE(ideg <= tdim + 1); 


  #ifndef USE_LPLIB_SMOOTHINTERIOR
  const int mball = 100;
  intAr1 lball(mball);
  dblAr1 qball(mball);
  #endif

  msh.tag[ithrd1]++;

  auto quafun = get_quafun<MFT,idim,idim>(iquaf);

  double noper = 0;
  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh);

    #if 0
    for(int ii = 0; ii < msh.npoin; ii++){
      rpoqe[ii] = 0.0;
      ipone[ii] = 0;
    }
    #endif

    double qmin = 1.0e30,qmax = -1.0e30, qnrm = 0.0;
    int imax = -1;
    int navg = 0;
    for(int ientt = 0; ientt < nentt; ientt++){

      if(isdeadent(ientt,ent2poi)) continue;

      navg++;

      double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,ientt,qpower,qpnorm,difto);

      qnrm += quael;
      qmin = MIN(qmin,quael);
      if(qmax < quael){
        imax = ientt;
        qmax = quael;
      }
    }

    qnrm /= navg;
    double t0 = get_wall_time();
    CPRINTF1(" - smoo iter %3d init %10.6e < q < %10.6e (at %d), avg = %10.6e " 
                   "(p = %d)\n",niter,qmin,qmax,imax,qnrm,qpnorm);
    //if(iverb >= 2 && qmax >= 1e10){
    //  printf("## HIGH QMAX mshdeg = %d \n",msh.curdeg);
    //  std::string fname = "qmax"+std::to_string(imax);
    //  writeMesh(fname,msh);
    //  //wait();
    //} 

    int nsucc = 0;
    int nmov  = 0;

    #if 0
    #ifdef USE_LPLIB_SMOOTHINTERIOR
    void (*inerloop_LPlib)(int,int,int,Mesh<MFT>*,int,int,double*,
          intAr1*,intAr1*,double,double)
    = [] (int ipoi0, int ipoi1, int ithread, Mesh<MFT> *msh,
          int itag_shared, int itag2, double *qmax, 
          intAr1 *nsuccthr, intAr1 *nmovthr, 
          double tolavg, double tolmax){

      const int mball = 100;
      constexpr int nnmet = (idim*(idim+1))/2;
      const int iverb = msh->param->iverb;
      intAr2 &ent2poi = tdim == 2 ? msh->fac2poi : msh->tet2poi;

      intAr1 lball(mball);
      dblAr1 qball(mball);

      for(int ipoin = ipoi0 - 1; ipoin < ipoi1; ipoin++){
        
        if(msh->poi2tag(itag_shared,ipoin) >= msh->tag[itag_shared]) continue;

        int ib = msh->poi2bpo[ipoin];
        if(ib >= 0) continue;

        int ientt = getpoient(*msh, ipoin, tdim);

        //double qpoin = rpoqe[ipoin] / ipone[ipoin];
        //if(qpoin > qnrm / qrthr){
          if(iverb >= 3){
            //printf("   - smoo pt %d seed elt %d quapt = %10.6e" 
            //  " qthrs = %10.6e qnrm = %10.6e\n",
            //  ipoin,ientt,qpoin,qrthr * qnrm,qnrm);
            printf("   - smoo pt %d seed elt %d \n", ipoin,ientt);
          }

          int iopen;
          bool imani = false;
          int ierro = 0,itmp = 0;
          if constexpr (idim == 2){
            intAr1 dum; 
            ierro = ball2(*msh,ipoin,ientt,lball,dum,&iopen,&imani,itag2);
          }else{
            ierro = ball3(*msh,ipoin,ientt,lball,&iopen,itag2);
          }
          METRIS_ASSERT(ierro == 0); 
          METRIS_ASSERT(iopen == 0);
          METRIS_ASSERT(imani == true);

          double coor0[idim];
          double met0[nnmet];
          for(int ii = 0; ii < idim; ii++) coor0[ii] = msh->coord(ipoin,ii);
          for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh->met(ipoin,ii);
          double qnrm0, qmax0, qnrm1, qmax1;
          try{
            //ierro = smooballdirect<MFT,idim,ideg>(msh,ipoin,lball,qball,
            //                       &qnrm0,&qmax0,&qnrm1,&qmax1,
            //                       qpower,qpnorm,difto,maxwt,inorm,iverb,itag2);
            ierro = smooballdiff<MFT,idim,ideg>(*msh,ipoin,lball,qball,
                                   &qnrm0,&qmax0,&qnrm1,&qmax1);//maxinc_worst,
            if(qmax1 > *qmax){
              if(iverb >= 2) printf("  - reject move, worst above last worst "
                " %15.7e > %15.7e\n", qmax1, *qmax);
              for(int ii = 0; ii < idim; ii++) 
                msh->coord(ipoin,ii) = coor0[ii];
              for(int ii = 0; ii < nnmet;ii++) 
                msh->met(ipoin,ii)   =  met0[ii];
              ierro = 1;
            }
          }catch(const MetrisExcept &e){
            printf("## FAILED  smooballdirect : METRIC INVALID\n");
            writeMesh("smooth_error.meshb",*msh);
            throw(e);
          }
          if(ierro == 0){
            (*nsuccthr)[ithread]++;
            if(iverb >= 3) printf("   - success smoothing %d q avg" 
                                     " %10.6e -> %10.6e max %10.6e -> %10.6e\n",
                                     ipoin,qnrm0,qnrm1,qmax0,qmax1);

            bool imov = false;
            // qnrm1 should be < qnrm0 for there to be progress 
            if(qnrm0 - qnrm1 > tolavg) imov = true;
            // idem qmax 
            if(qmax0 - qmax1 > tolmax) imov = true;
            if(imov){
              (*nmovthr)[ithread]++;
              for(int iele2 : lball){
                for(int ii = 0; ii < idim+1; ii++){
                  int ipoi2 = ent2poi(iele2,ii);
                  msh->poi2tag(itag_shared,ipoin) = msh->tag[itag_shared] - 1; // reactivate
                }
              }
            }else{
              msh->poi2tag(itag_shared,ipoin) = msh->tag[itag_shared]; // deactivate
            }
          }else{
            msh->poi2tag(itag_shared,ipoin) = msh->tag[itag_shared]; // deactivate
          }

        //} // if qpoin
      } // for ientt
      // Control sizes here if provided (hmin hmax)
    };

    float acc = LaunchParallelMultiArg(LibIdx, LP_elt, LP_poi, 
                                       (void*)inerloop_LPlib, 8, 
                                       &msh, itag_shared, ithrd2, &qmax, 
                                       &nsuccthr, &nmovthr, 
                                       tolavg, tolmax);
    CPRINTF1("Smoothing accel = %f \n",acc);
    for(int ii = 0; ii < nthread; ii++){
      nsucc += nsuccthr[ii];
      nmov  += nmovthr[ii];
    }

    #endif
    #endif


    #ifndef USE_LPLIB_SMOOTHINTERIOR
    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
      INCVDEPTH(msh);

      int ib = msh.poi2bpo[ipoin];
      if(ib >= 0) continue;

      int ientt = getpoient(msh, ipoin, tdim);  
      int iver = tdim == 2 ? getverfac<ideg>(ientt, ent2poi, ipoin)
                           : getvertet<ideg>(ientt, ent2poi, ipoin);

      METRIS_ASSERT(iver >= 0);


      CPRINTF1(" - smoo pt %d seed elt %d \n", ipoin,ientt);
      int ierro = 0;

      if(iver < tdim+1){
        bool imani = false;
        int iopen;
        // Vertex case 
        if constexpr (idim == 2){
          intAr1 dum; 
          ierro = ball2(msh,ipoin,ientt,lball,dum,&iopen,&imani,ithrd2);
        }else{
          ierro = ball3(msh,ipoin,ientt,lball,&iopen,ithrd2);
        }
        METRIS_ASSERT(iopen == 0);
        METRIS_ASSERT(imani == true);
      }else{
        // HO node
        if constexpr (tdim == 2){
          lball.set_n(0);
          lball.stack(ientt);
          int nppe = edgnpps[ideg] - 2;
          int ied = (iver - (tdim + 1)) / nppe;
          METRIS_ASSERT(ied < 4);

          int ifnei = ent2ent(ientt,ied);
          // nm not impl yet and bdry should never happen
          METRIS_ASSERT(ifnei >= 0);

          lball.stack(ifnei);
        }else{
          METRIS_THROW_MSG(TODOExcept(), "Call shell3");
        }

      }
      METRIS_ASSERT(ierro == 0); 

      double coor0[idim];
      double met0[nnmet];
      for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipoin,ii);
      for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh.met(ipoin,ii);
      double qnrm0, qmax0, qnrm1, qmax1;
      try{
        //ierro = smooballdirect<MFT,idim,ideg>(msh,ipoin,lball,qball,
        //                       &qnrm0,&qmax0,&qnrm1,&qmax1,
        //                       qpower,qpnorm,difto,maxwt,inorm,iverb,ithrd2);
        if(msh.param->iflag2 == 0){
          ierro = smooballdiff<MFT,idim,ideg>(msh,ipoin,lball,qball,
                                 &qnrm0,&qmax0,&qnrm1,&qmax1,iquaf);
        }else{
          ierro = smooballdiff_luksan<MFT,idim,ideg>(msh,ipoin,lball,qball,
                                     &qnrm0,&qmax0,&qnrm1,&qmax1,work,iquaf);
        }
        if(qmax1 > qmax){
          CPRINTF1(" - reject move, worst above last worst %15.7e > %15.7e\n", 
                   qmax1, qmax);
          for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii];
          for(int ii = 0; ii < nnmet;ii++) msh.met(ipoin,ii)   =  met0[ii];
          ierro = 1;
        }
      }catch(const MetrisExcept &e){
        printf("## FAILED  smooballdirect\n");
        writeMesh("smooth_error.meshb",msh);
        throw(e);
      }
      if(ierro == 0){
        nsucc++;
        CPRINTF1(" - success smoothing %d q avg %10.6e -> %10.6e " 
                 "max %10.6e -> %10.6e\n",ipoin,qnrm0,qnrm1,qmax0,qmax1);

        bool imov = false;
        // qnrm1 should be < qnrm0 for there to be progress 
        if(qnrm0 - qnrm1 > tolavg) imov = true;
        // idem qmax 
        if(qmax0 - qmax1 > tolmax) imov = true;
        if(imov){
          nmov ++;
          for(int iele2 : lball){
            for(int ii = 0; ii < idim+1; ii++){
              int ipoi2 = ent2poi(iele2,ii);
              msh.poi2tag(ithrd1,ipoi2) = msh.tag[ithrd1] - 1; // reactivate
            }
          }
        }else{
          msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1]; // deactivate
        }
      }else{
        msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1]; // deactivate
      }

    } // for ipoin // for ientt
    #endif

    double t1 = get_wall_time();
    CPRINTF1(" - Iteration end t = %f nsuccess = %d nmov = %d \n",
                          t1-t0,nsucc,nmov);
    noper += nmov;
    if(nmov == 0) break;
  } // end for niter

  return noper / (double) nentt;
}


#define BOOST_PP_LOCAL_MACRO(n)\
template double smoothInterior_Ball0<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh,\
                                        QuaFun iquaf, int ithrd1, int ithrd2);\
template double smoothInterior_Ball0<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh,\
                                        QuaFun iquaf, int ithrd1, int ithrd2);\
template double smoothInterior_Ball0<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh,\
                                        QuaFun iquaf, int ithrd1, int ithrd2);\
template double smoothInterior_Ball0<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh,\
                                        QuaFun iquaf, int ithrd1, int ithrd2);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


} // end namespace
