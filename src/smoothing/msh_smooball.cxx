//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
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
#include "../utils/aux_timer.hxx"
#include "../low_topo.hxx"
#include "../utils/mprintf.hxx"
#include "../quality/low_metqua.hxx"
#include "../io_libmeshb.hxx"

#include "lplib3/lplib3.h"

#undef USE_LPLIB_SMOOTHINTERIOR


namespace Metris{

template<class MFT>
double smoothInterior_Ball(Mesh<MFT> &msh, QuaFun iquaf, int ithrd1, int ithrd2){

  int tdimn = msh.get_tdim();

  METRIS_ASSERT_MSG(tdimn > 1, "TODO: edge smooth interior ball");

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


// =============================================================================================== //
// =============================================================================================== //

// idim: gdim = tdim
template<class MFT, int idim, int ideg>
double smoothInterior_Ball0(Mesh<MFT> &msh, QuaFun iquaf,
                            int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);

  constexpr int tdim = idim;
  //constexpr int gdim = idim;
  if(ideg > idim + 1){
    PRINTF("## SMOOTHING DISABLED FOR DEGREE {} and dim {} \n",ideg,idim);
    return -1.0;
  }


  int nentt = msh.nentt(tdim);
  const intAr2 &ent2poi = msh.ent2poi(tdim);
  const intAr2 &ent2ent = msh.ent2ent(tdim);


  #ifdef USE_LPLIB_SMOOTHINTERIOR
    // LPlib init
    int nproc = msh.param->nproc;
    int nthread = GetNumberOfCores();
    if(nthread <= 0){
      CPRINTF1("## WARNING: LPlib function GetNumberOfCores() returned "
               "negative threads. Set to default {}.\n",METRIS_MAXTAGS);
      nthread = METRIS_MAXTAGS;
    }else{
      CPRINTF2("-- LPlib found ncore = {} \n",nthread);
      if(nthread > METRIS_MAXTAGS){
        CPRINTF1("## WARNING: must verify nthread <= METRIS_MAXTAGS = {}."
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

  //msh.met.setSpace(MetSpace::Log);



  // Eventually move all constants to MetrisParameters
  // L2 conformity error from 0 to 1
  const double difto = 1.0;
  const int miter = msh.param->opt_smoo_niter;
  //const double maxwt = 20.0;
  //const double qrthr = 2.0;
  const double tolavg = msh.param->opt_smoo_tol;
  const double tolmax = msh.param->opt_smoo_tol;

  dblAr1 work;
  if(msh.param->iflag2 != 0){
    MPRINTF("\n\n##WARNING EXPERIMENTAL SMOOTHING FUNCTION 2\n");
  }

  // 1 -> no maximum quality increase allowed
  //const double maxinc_worst = 1.00;

  constexpr int nnmet = (idim*(idim+1))/2;

  METRIS_ASSERT(msh.param->opt_power == 1 || msh.param->opt_power == -1);
  // Otherwise not only edge nodes
  METRIS_ENFORCE(ideg <= tdim + 1);


  #ifndef USE_LPLIB_SMOOTHINTERIOR
  const int mball = 100;
  intAr1 lball(mball);
  #endif

  msh.tag[ithrd1]++;

  auto quafun = get_quafun<MFT,idim,idim>(iquaf);

  double noper = 0;
  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh.param);

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

      double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,ientt,difto);

      qnrm += quael;
      qmin = MIN(qmin,quael);
      if(qmax < quael){
        imax = ientt;
        qmax = quael;
      }
    }

    qnrm /= navg;
    double t0 = get_cpu_time();
    CPRINTF1(" - smoo iter {:3} init {:10.6e} < q < {:10.6e} (at {}), avg = {:10.6e} "
                   "(p = {})\n",niter,qmin,qmax,imax,qnrm,msh.param->opt_pnorm);
    //if(iverb >= 2 && qmax >= 1e10){
    //  printf("## HIGH QMAX mshdeg = {} \n",msh.curdeg);
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
            //printf("   - smoo pt {} seed elt {} quapt = {:10.6e}"
            //  " qthrs = {:10.6e} qnrm = {:10.6e}\n",
            //  ipoin,ientt,qpoin,qrthr * qnrm,qnrm);
            printf("   - smoo pt {} seed elt {} \n", ipoin,ientt);
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
                " {:15.7e} > {:15.7e}\n", qmax1, *qmax);
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
            if(iverb >= 3) printf("   - success smoothing {} q avg"
                                     " {:10.6e} -> {:10.6e} max {:10.6e} -> {:10.6e}\n",
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
    CPRINTF1("Smoothing accel = {} \n",acc);
    for(int ii = 0; ii < nthread; ii++){
      nsucc += nsuccthr[ii];
      nmov  += nmovthr[ii];
    }

    #endif
    #endif


    #ifndef USE_LPLIB_SMOOTHINTERIOR

    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){

      // if (ipoin != 15361) continue;

      // skip if dead point or tagged
      if(msh.isdeadpoint(ipoin)) continue;
      if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
      INCVDEPTH(msh.param);

      // check if boundary point
      int ibpoin = msh.poi2bpo[ipoin];
      bool pointOnEdge = false;
      bool pointOnFace = false;
      if (ibpoin >= 0){
        METRIS_ASSERT(msh.bpo2ibi(ibpoin,0) == ipoin);

        #if !defined(SMOOTHEDGES) && !defined(SMOOTFACES)
        continue;
        #endif

        if (msh.bpo2ibi(ibpoin,1) == 0) continue;           // point is a CAD corner, don't move

        #ifdef SMOOTHEDGES
        if (msh.bpo2ibi(ibpoin,1) == 1) pointOnEdge = true; // point is in CAD edge
        #else
        if (msh.bpo2ibi(ibpoin,1) == 1) continue;
        #endif

        #ifdef SMOOTHFACES
        if (msh.bpo2ibi(ibpoin,1) == 2) pointOnFace = true; // point is in CAD face
        #else
        if (msh.bpo2ibi(ibpoin,1) == 2) continue;
        #endif
      }
      const bool interiorPoint = !pointOnEdge && !pointOnFace;

      int ientt = getpoient(msh, ipoin, tdim);
      int iver = tdim == 2 ? msh.template getverfac<ideg>(ientt, ipoin)
                           : msh.template getvertet<ideg>(ientt, ipoin);

      METRIS_ASSERT(iver >= 0);


      CPRINTF1(" - smoo pt {} seed elt {} \n", ipoin,ientt);
      int ierro = 0;

      if(iver < tdim+1){
        int iopen;
        // Vertex case
        if constexpr (idim == 2){
          intAr1 dum;
          bool imani = false;
          ierro = ball2(msh,ipoin,ientt,lball,dum,&iopen,&imani,ithrd2);
          METRIS_ASSERT(imani == true);
        }else{
          ierro = ball3(msh,ipoin,ientt,lball,&iopen,ithrd2);
        }
        // METRIS_ASSERT(iopen == 0);
      }else{
        // HO node
        int nppe = getnnod1(ideg) - 2;
        int ied = (iver - (tdim + 1)) / nppe;
        if constexpr (tdim == 2){
          lball.set_n(0);
          lball.stack(ientt);
          METRIS_ASSERT(ied < 4);

          int ifnei = ent2ent(ientt,ied);
          // nm not impl yet and bdry should never happen
          METRIS_ASSERT(ifnei >= 0);

          lball.stack(ifnei);
        }else{
          int ip1 = msh.tet2poi(ientt, lnoed3[ied][0]);
          int ip2 = msh.tet2poi(ientt, lnoed3[ied][1]);
          static intAr1 dum;
          int iopen;
          shell3(msh, ip1, ip2, ientt, lball, dum, &iopen);
          // METRIS_ASSERT(iopen < 0);
        }

      }
      METRIS_ASSERT(ierro == 0);

      double coor0[idim];
      double met0[nnmet];
      for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipoin,ii);
      for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh.met(ipoin,ii);
      double tparam0 = -1e30;
      if (pointOnEdge) tparam0 = msh.bpo2rbi(ibpoin,0); // back up for t parameter
      double uparam0 = -1e30;
      double vparam0 = -1e30;
      if (pointOnFace){
        // backup u and v params
        uparam0 = msh.bpo2rbi(ibpoin,0);
        vparam0 = msh.bpo2rbi(ibpoin,1);
      }

      // std::cout << "printing qualities of elements in ball BEFORE smoobaldiff" << std::endl;

      // for (int ieleInBall : lball){

      //   double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,ieleInBall,difto);
      //   std::cout << "ele = " << ieleInBall << " q = " << quael << std::endl;
      // }

      double qnrm0, qmax0, qnrm1, qmax1;
      try{
        //ierro = smooballdirect<MFT,idim,ideg>(msh,ipoin,lball,qball,
        //                       &qnrm0,&qmax0,&qnrm1,&qmax1,
        //                       qpower,qpnorm,difto,maxwt,inorm,iverb,ithrd2);
        if(msh.param->iflag2 == 0){
          if (pointOnEdge)   ierro = smooballdiff_boundary<MFT,idim,ideg>(msh,ipoin,1,lball,&qnrm0,&qmax0,&qnrm1,&qmax1,iquaf);
          if (pointOnFace)   ierro = smooballdiff_boundary<MFT,idim,ideg>(msh,ipoin,2,lball,&qnrm0,&qmax0,&qnrm1,&qmax1,iquaf);
          if (interiorPoint) ierro = smooballdiff<MFT,idim,ideg>(msh,ipoin,lball,&qnrm0,&qmax0,&qnrm1,&qmax1,iquaf);
        }else{
          ierro = smooballdiff_luksan<MFT,idim,ideg>(msh,ipoin,lball,
                                     &qnrm0,&qmax0,&qnrm1,&qmax1,work,iquaf);
        }
        if(qmax1 > qmax){
          CPRINTF1(" - reject move, worst above last worst {:15.7e} > {:15.7e}\n",
                   qmax1, qmax);
          for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii];
          for(int ii = 0; ii < nnmet;ii++) msh.met(ipoin,ii)   =  met0[ii];
          if (pointOnEdge) msh.bpo2rbi(ibpoin,0) = tparam0;
          if (pointOnFace){
            msh.bpo2rbi(ibpoin,0) = uparam0;
            msh.bpo2rbi(ibpoin,1) = vparam0;
          }
          ierro = 1;
        }else if(pointOnEdge){

          // important: as we moved the point along the edge, we also need to update the (u,v) parameters
          // of the point as part of the faces incident to the edge
          double topt = msh.bpo2rbi(ibpoin,0);

          const int iedge = msh.bpo2ibi(ibpoin,2); // a mesh edge attached to the point
          METRIS_ASSERT(iedge >= 0 && iedge < msh.nedge);

          const int irefEdge = msh.edg2ref[iedge]; // get CAD edge reference
          ego cadEdge = msh.CAD.cad2edg[irefEdge]; // get actual CAD edge object

          int ibpoinRecord = ibpoin;
          while (ibpoinRecord != -1 && msh.bpo2ibi(ibpoinRecord,0) == ipoin){

            if (msh.bpo2ibi(ibpoinRecord,1) == 2){ // same point we moved but as living in a face

              // first get the CAD face object

              const int iface = msh.bpo2ibi(ibpoinRecord,2); // a mesh face attached to the point
              METRIS_ASSERT(iface >= 0 && iface < msh.nface);

              const int irefFace = msh.fac2ref[iface]; // get CAD face reference
              ego cadFace = msh.CAD.cad2fac[irefFace]; // get actual CAD face object

              // retrieve (u,v) of the point given the new topt location along edge
              double uv[2];
              int icode = EG_getEdgeUV(cadFace, cadEdge, 0, topt, uv);
              METRIS_ENFORCE_MSG(icode == EGADS_SUCCESS, "EG_getEdgeUV failed: {}", icode);

              msh.bpo2rbi(ibpoinRecord,0) = uv[0];
              msh.bpo2rbi(ibpoinRecord,1) = uv[1];

            } // finish updating (u,v) of a face incident to the edge

            // move on to next boundary point record in the linked list
            ibpoinRecord = msh.bpo2ibi(ibpoinRecord,3);

          }

        }
      }catch(const MetrisExcept &e){
        PRINTF("## FAILED  smoothing\n");
        MPRINTF("ipoin = {}", ipoin);
        writeMesh("smooth_error.meshb",msh);
        msh.met.writeMetricFile("smooth_error_metric.solb");
        for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii];
        for(int ii = 0; ii < nnmet;ii++) msh.met(ipoin,ii)   =  met0[ii];
        writeMesh("smooth_error_0.meshb",msh);
        msh.met.writeMetricFile("smooth_error_metric_0.solb");
        ierro = 1;
        // throw(e);
      }
      if(ierro == 0){
        nsucc++;
        CPRINTF1(" - success smoothing {} q avg {:10.6e} -> {:10.6e} "
                 "max {:10.6e} -> {:10.6e}\n",ipoin,qnrm0,qnrm1,qmax0,qmax1);


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

    double t1 = get_cpu_time();
    CPRINTF1(" - Iteration end time = {:.2e}s nsuccess = {} nmov = {} \n",
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

// =============================================================================================== //
// =============================================================================================== //

template<class MFT>
double smoothElement_Ball(Mesh<MFT> &msh, const int ientt, BadEntHandler& handler, QuaFun iquaf, int ithrd1, int ithrd2){

  int tdimn = msh.get_tdim();

  METRIS_ASSERT_MSG(tdimn > 1, "TODO: edge smooth interior ball");

  // Geo and topo dimn must match otherwise surface specific
  METRIS_ASSERT(tdimn == msh.idim);
  double noper;
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    if(tdimn == 2){
      noper = smoothElement_Ball0<MFT,2,ideg>(msh,ientt,handler,iquaf,ithrd1,ithrd2);
    }else{
      noper = smoothElement_Ball0<MFT,3,ideg>(msh,ientt,handler,iquaf,ithrd1,ithrd2);
    }
  }}CT_FOR1(ideg);

  return noper;
}

template double smoothElement_Ball<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, const int ientt, BadEntHandler& handler,
                                             QuaFun iquaf, int ithrd1, int ithrd2);
template double smoothElement_Ball<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, const int ientt, BadEntHandler& handler,
                                             QuaFun iquaf, int ithrd1, int ithrd2);

// =============================================================================================== //
// =============================================================================================== //


// idim: gdim = tdim
template<class MFT, int idim, int ideg>
double smoothElement_Ball0(Mesh<MFT> &msh, const int ientt, BadEntHandler& handler, QuaFun iquaf,
                           int ithrd1, int ithrd2){
  // TODO: for now just vertex smoothing. Add HO nodes in the future
  GETVDEPTH(msh.param);

  constexpr int tdim = idim;
  //constexpr int gdim = idim;
  if(ideg > idim + 1){
    PRINTF("## SMOOTHING DISABLED FOR DEGREE {} and dim {} \n",ideg,idim);
    return -1.0;
  }

  int nentt = msh.nentt(tdim);
  const intAr2 &ent2poi = msh.ent2poi(tdim);
  const intAr2 &ent2ent = msh.ent2ent(tdim);

  // Eventually move all constants to MetrisParameters
  // L2 conformity error from 0 to 1
  const double difto = 1.0;
  const int miter = msh.param->opt_smoo_niter;
  //const double maxwt = 20.0;
  //const double qrthr = 2.0;
  const double tolavg = msh.param->opt_smoo_tol;
  const double tolmax = msh.param->opt_smoo_tol;

  dblAr1 work;
  if(msh.param->iflag2 != 0){
    MPRINTF("\n\n##WARNING EXPERIMENTAL SMOOTHING FUNCTION 2\n");
  }

  constexpr int nnmet = (idim*(idim+1))/2;

  // METRIS_ENFORCE(msh.param->opt_power < 0); // Otherwise rework the mins / maxs
  // Otherwise not only edge nodes
  METRIS_ENFORCE(ideg <= tdim + 1);

  const int mball = 100;
  intAr1 lball(mball);

  msh.tag[ithrd1]++;

  auto quafun = get_quafun<MFT,idim,idim>(iquaf);

  double noper = 0;
  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh.param);

    #if 0
    for(int ii = 0; ii < msh.npoin; ii++){
      rpoqe[ii] = 0.0;
      ipone[ii] = 0;
    }
    #endif

    int nsucc = 0;
    int nmov  = 0;

    // loop over entity vertices
    for (int ii = 0; ii < tdim+1; ii++){

      const int ipoin = ent2poi(ientt, ii);

      // skip if dead point or tagged
      if(msh.isdeadpoint(ipoin)) continue;
      if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
      INCVDEPTH(msh.param);

     // check if boundary point
      int ibpoin = msh.poi2bpo[ipoin];
      bool pointOnEdge = false;
      bool pointOnFace = false;
      if (ibpoin >= 0){
        METRIS_ASSERT(msh.bpo2ibi(ibpoin,0) == ipoin);

        #if !defined(SMOOTHEDGES) && !defined(SMOOTFACES)
        continue;
        #endif

        if (msh.bpo2ibi(ibpoin,1) == 0) continue;           // point is a CAD corner, don't move

        #ifdef SMOOTHEDGES
        if (msh.bpo2ibi(ibpoin,1) == 1) pointOnEdge = true; // point is in CAD edge
        #else
        if (msh.bpo2ibi(ibpoin,1) == 1) continue;
        #endif

        #ifdef SMOOTHFACES
        if (msh.bpo2ibi(ibpoin,1) == 2) pointOnFace = true; // point is in CAD face
        #else
        if (msh.bpo2ibi(ibpoin,1) == 2) continue;
        #endif
      }
      const bool interiorPoint = !pointOnEdge && !pointOnFace;

      CPRINTF1(" - smoo pt {} seed elt {} \n", ipoin,ientt);
      int ierro = 0;

      int iopen;
      if constexpr (idim == 2){
        intAr1 dum;
        bool imani = false;
        ierro = ball2(msh,ipoin,ientt,lball,dum,&iopen,&imani,ithrd2);
        METRIS_ASSERT(imani == true);
      }else{
        ierro = ball3(msh,ipoin,ientt,lball,&iopen,ithrd2);
      }
      #ifndef NDEBUG
      // if (!pointOnEdge)  METRIS_ASSERT(iopen == 0);
      // else               METRIS_ASSERT(iopen == 1);
      METRIS_ASSERT(ierro == 0);
      #endif

      // we now have the ball for ipoin, get initial qualities

      double qmin = 1.0e30,qmax = -1.0e30, qnrm = 0.0;
      double qsum = 0.;
      int imax = -1;
      int navg = 0;
      for(auto ient2 : lball){

        if(isdeadent(ient2,ent2poi)) continue;

        navg++;

        double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,ient2,difto);

        qnrm += quael;
        qmin = MIN(qmin,quael);
        if(qmax < quael){
          imax = ient2;
          qmax = quael;
        }
      }

      qsum = qnrm;
      qnrm /= navg;
      double t0 = get_cpu_time();
      CPRINTF1(" - smoo iter {:3}, ipoin {} (ientt {}, point {}), ball init {:10.6e} < q < {:10.6e} (at {}), avg = {:10.6e} "
                     "(p = {})\n",niter,ipoin,ientt,ii,qmin,qmax,imax,qnrm,msh.param->opt_pnorm);

      double coor0[idim];
      double met0[nnmet];
      for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipoin,ii);
      for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh.met(ipoin,ii);
      double tparam0 = -1e30;
      double uparam0 = -1e30;
      double vparam0 = -1e30;
      if (pointOnEdge) tparam0 = msh.bpo2rbi(ibpoin,0); // back up for t parameter
      if (pointOnFace){
        // backup u and v params
        uparam0 = msh.bpo2rbi(ibpoin,0);
        vparam0 = msh.bpo2rbi(ibpoin,1);
      }
      double qnrm0, qmax0, qnrm1, qmax1;
      try{
        if(msh.param->iflag2 == 0){
          if (pointOnEdge)   ierro = smooballdiff_boundary<MFT,idim,ideg>(msh,ipoin,1,lball,&qnrm0,&qmax0,&qnrm1,&qmax1,iquaf);
          if (pointOnFace)   ierro = smooballdiff_boundary<MFT,idim,ideg>(msh,ipoin,2,lball,&qnrm0,&qmax0,&qnrm1,&qmax1,iquaf);
          if (interiorPoint) ierro = smooballdiff<MFT,idim,ideg>(msh,ipoin,lball,&qnrm0,&qmax0,&qnrm1,&qmax1,iquaf);
        }else{
          ierro = smooballdiff_luksan<MFT,idim,ideg>(msh,ipoin,lball,
                                     &qnrm0,&qmax0,&qnrm1,&qmax1,work,iquaf);
        }
        double qsumNew = 0.;
        double qmaxNew = 0.;
        for(auto ient2 : lball){

          if(isdeadent(ient2,ent2poi)) continue;

          double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,ient2,difto);
          qsumNew += quael;
          if (quael > qmaxNew) qmaxNew = quael;

        }
        bool improveQuaSum = handler.checkSuccess(qsumNew,qsum);
        bool improveQuaMax = true;
        #ifdef IMPROVEMAXQUAL
          improveQuaMax = qmaxNew < qmax;
        #endif
        if(!improveQuaSum || !improveQuaMax){
          CPRINTF1(" - reject move, quality error increased: {:15.7e} > {:15.7e}\n",
                   qsumNew, qsum);
          for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii];
          for(int ii = 0; ii < nnmet;ii++) msh.met(ipoin,ii)   =  met0[ii];
          #ifdef SMOOTHEDGES
          if (pointOnEdge) msh.bpo2rbi(ibpoin,0) = tparam0;
          #endif
          #ifdef SMOOTHFACES
          if (pointOnFace){
            msh.bpo2rbi(ibpoin,0) = uparam0;
            msh.bpo2rbi(ibpoin,1) = vparam0;
          }
          #endif

          ierro = 1;
        }else if(pointOnEdge){

          // important: as we moved the point along the edge, we also need to update the (u,v) parameters
          // of the point as part of the faces incident to the edge
          double topt = msh.bpo2rbi(ibpoin,0);

          const int iedge = msh.bpo2ibi(ibpoin,2); // a mesh edge attached to the point
          METRIS_ASSERT(iedge >= 0 && iedge < msh.nedge);

          const int irefEdge = msh.edg2ref[iedge]; // get CAD edge reference
          ego cadEdge = msh.CAD.cad2edg[irefEdge]; // get actual CAD edge object

          int ibpoinRecord = ibpoin;
          while (ibpoinRecord != -1 && msh.bpo2ibi(ibpoinRecord,0) == ipoin){

            if (msh.bpo2ibi(ibpoinRecord,1) == 2){ // same point we moved but as living in a face

              // first get the CAD face object

              const int iface = msh.bpo2ibi(ibpoinRecord,2); // a mesh face attached to the point
              METRIS_ASSERT(iface >= 0 && iface < msh.nface);

              const int irefFace = msh.fac2ref[iface]; // get CAD face reference
              ego cadFace = msh.CAD.cad2fac[irefFace]; // get actual CAD face object

              // retrieve (u,v) of the point given the new topt location along edge
              double uv[2];
              int icode = EG_getEdgeUV(cadFace, cadEdge, 0, topt, uv);
              METRIS_ENFORCE_MSG(icode == EGADS_SUCCESS, "EG_getEdgeUV failed: {}", icode);

              msh.bpo2rbi(ibpoinRecord,0) = uv[0];
              msh.bpo2rbi(ibpoinRecord,1) = uv[1];

            } // finish updating (u,v) of a face incident to the edge

            // move on to next boundary point record in the linked list
            ibpoinRecord = msh.bpo2ibi(ibpoinRecord,3);

          }

        }
      }catch(const MetrisExcept &e){
        PRINTF("## FAILED  smooballdirect\n");
        writeMesh("smooth_error.meshb",msh);
        throw(e);
      }
      if(ierro == 0){
        nsucc++;
        CPRINTF1(" - success smoothing {} q avg {:10.6e} -> {:10.6e} "
                 "max {:10.6e} -> {:10.6e}\n",ipoin,qnrm0,qnrm1,qmax0,qmax1);

        nmov++;
        for(int ient2 : lball){
          for(int ii = 0; ii < idim+1; ii++){
            int ipoi2 = ent2poi(ient2,ii);
            msh.poi2tag(ithrd1,ipoi2) = msh.tag[ithrd1] - 1; // reactivate
          }

          // inform affected entities: id and new quality
          handler.affectedEnttsAlive[ient2] = quafun(msh,AsDeg::Pk,AsDeg::Pk,ient2,difto);
        }
      }else{
        msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1]; // deactivate
      }

      double t1 = get_cpu_time();
      CPRINTF1(" - Iteration end time = {:.2e}s nsuccess = {} nmov = {} \n",
                            t1-t0,nsucc,nmov);
      noper += nmov;
      if(nmov == 0) break;
    } // ii (vertex in element)
  } // iter

  return noper / (double) nentt;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template double smoothElement_Ball0<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh,\
                                        const int ientt, BadEntHandler& handler, QuaFun iquaf, int ithrd1, int ithrd2);\
template double smoothElement_Ball0<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh,\
                                        const int ientt, BadEntHandler& handler, QuaFun iquaf, int ithrd1, int ithrd2);\
template double smoothElement_Ball0<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh,\
                                        const int ientt, BadEntHandler& handler, QuaFun iquaf, int ithrd1, int ithrd2);\
template double smoothElement_Ball0<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh,\
                                        const int ientt, BadEntHandler& handler, QuaFun iquaf, int ithrd1, int ithrd2);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

// =============================================================================================== //
// =============================================================================================== //

template<class MFT>
double smoothCavity(Mesh<MFT> &msh, MshCavity& cav, BadEntHandler& handler, QuaFun iquaf,
                    const double quaCav0, const double quaMaxCav0,
                    double& quaCav1, double& quaMaxCav1,
                    int ithrd1, int ithrd2){

  int tdimn = msh.get_tdim();

  METRIS_ASSERT_MSG(tdimn > 1, "TODO: edge smooth interior ball");

  // Geo and topo dimn must match otherwise surface specific
  METRIS_ASSERT(tdimn == msh.idim);
  double noper;
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    if(tdimn == 2){
      noper = smoothCavity0<MFT,2,ideg>(msh,cav,handler,iquaf,quaCav0,quaMaxCav0,quaCav1,quaMaxCav1,ithrd1,ithrd2);
    }else{
      noper = smoothCavity0<MFT,3,ideg>(msh,cav,handler,iquaf,quaCav0,quaMaxCav0,quaCav1,quaMaxCav1,ithrd1,ithrd2);
    }
  }}CT_FOR1(ideg);

  return noper;
}

template double smoothCavity<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, MshCavity& cav, BadEntHandler& handler, QuaFun iquaf,
                                                    const double quaCav0, const double quaMaxCav0,
                                                    double& quaCav1, double& quaMaxCav1,
                                                    int ithrd1, int ithrd2);
template double smoothCavity<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, MshCavity& cav, BadEntHandler& handler, QuaFun iquaf,
                                                    const double quaCav0, const double quaMaxCav0,
                                                    double& quaCav1, double& quaMaxCav1,
                                                    int ithrd1, int ithrd2);

// =============================================================================================== //
// =============================================================================================== //


// idim: gdim = tdim
template<class MFT, int idim, int ideg>
double smoothCavity0(Mesh<MFT> &msh, MshCavity& cav, BadEntHandler& handler, QuaFun iquaf, const double quaCav0, const double quaMaxCav0, double& quaCav1, double& quaMaxCav1,
                           int ithrd1, int ithrd2){
  // TODO: for now just vertex smoothing. Add HO nodes in the future
  GETVDEPTH(msh.param);

  constexpr int tdim = idim;
  //constexpr int gdim = idim;
  if(ideg > idim + 1){
    PRINTF("## SMOOTHING DISABLED FOR DEGREE {} and dim {} \n",ideg,idim);
    return -1.0;
  }

  int nentt = msh.nentt(tdim);
  const intAr2 &ent2poi = msh.ent2poi(tdim);
  const intAr2 &ent2ent = msh.ent2ent(tdim);

  // Eventually move all constants to MetrisParameters
  // L2 conformity error from 0 to 1
  const double difto = 1.0;
  // const int miter = msh.param->opt_smoo_niter;
  const int miter = 1;
  //const double maxwt = 20.0;
  //const double qrthr = 2.0;
  const double tolavg = msh.param->opt_smoo_tol;
  const double tolmax = msh.param->opt_smoo_tol;

  constexpr int nnmet = (idim*(idim+1))/2;

  METRIS_ENFORCE(ideg <= tdim + 1);

  msh.tag[ithrd1]++;

  double noper = 0;
  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh.param);

    int nsucc = 0;
    int nmov  = 0;

    const int ipins = cav.ipins;

    // TODO:
    // skip if boundary point
    // int ib = msh.poi2bpo[ipoin];
    // if(ib >= 0) continue;

    CPRINTF1(" - smoo cavity for insertion pt {} \n", ipins);
    int ierro = 0;

    double t0 = get_cpu_time();

    double coor0[idim];
    double met0[nnmet];
    for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipins,ii);
    for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh.met(ipins,ii);
    double qnrm0, qmax0, qnrm1, qmax1;
    try{

      ierro = smoocavdiff<MFT,idim,ideg>(msh,cav,quaCav1,quaMaxCav1,iquaf,ithrd1);

      // std::cout << "In smoothCavity0" << std::endl;
      // std::cout << "quaCav1 = " << quaCav1 << std::endl;
      // std::cout << "quaMaxCav1 = " << quaMaxCav1 << std::endl;
      // std::cout << "quaCav0 = " << quaCav0 << std::endl;
      // std::cout << "quaMaxCav0 = " << quaMaxCav0 << std::endl;

      if (ierro == 0){
        bool improveCavSum = handler.checkSuccess(quaCav1,quaCav0);
        bool improveCavMax = true;
        #ifdef IMPROVEMAXQUAL
        improveCavMax = quaMaxCav1 <= quaMaxCav0;
        #endif
        if(!improveCavSum || !improveCavMax){
          CPRINTF1(" - reject move, quality error increased ({:15.7e} > {:15.7e}) or max error increased ({:15.7e} > {:15.7e})\n",
                    quaCav1, quaCav0, quaMaxCav1, quaMaxCav0);
          for(int ii = 0; ii < idim; ii++) msh.coord(ipins,ii) = coor0[ii];
          for(int ii = 0; ii < nnmet;ii++) msh.met(ipins,ii)   =  met0[ii];

          quaCav1 = quaCav0;
          quaMaxCav1 = quaMaxCav0;
          ierro = 1;
        }
      }
      else{
        for(int ii = 0; ii < idim; ii++) msh.coord(ipins,ii) = coor0[ii];
        for(int ii = 0; ii < nnmet;ii++) msh.met(ipins,ii)   =  met0[ii];
      }
    }catch(const MetrisExcept &e){
      PRINTF("## FAILED  smoocavdiff\n");
      writeMesh("smooth_error.meshb",msh);
      throw(e);
    }
    if(ierro == 0){
      nsucc++;
      CPRINTF1(" - success smoothing {} q avg {:10.6e} -> {:10.6e} "
                "max {:10.6e} -> {:10.6e}\n",ipins,quaCav0,quaCav1,quaMaxCav0,quaMaxCav1);

    }

    double t1 = get_cpu_time();
    CPRINTF1(" - Iteration end time = {:.2e}s nsuccess = {} nmov = {} \n",
                          t1-t0,nsucc);
    noper += nsucc;
    if(nsucc == 0) break;
  } // iter

  return noper / (double) nentt;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template double smoothCavity0<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh,\
                                        MshCavity& cav, BadEntHandler& handler, QuaFun iquaf, const double quaCav0, const double quaMaxCav0, double& quaCav1, double& quaMaxCav1, \
                           int ithrd1, int ithrd2);\
template double smoothCavity0<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh,\
                                        MshCavity& cav, BadEntHandler& handler, QuaFun iquaf, const double quaCav0, const double quaMaxCav0, double& quaCav1, double& quaMaxCav1, \
                           int ithrd1, int ithrd2);\
template double smoothCavity0<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh,\
                                        MshCavity& cav, BadEntHandler& handler, QuaFun iquaf, const double quaCav0, const double quaMaxCav0, double& quaCav1, double& quaMaxCav1, \
                           int ithrd1, int ithrd2);\
template double smoothCavity0<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh,\
                                        MshCavity& cav, BadEntHandler& handler, QuaFun iquaf, const double quaCav0, const double quaMaxCav0, double& quaCav1, double& quaMaxCav1, \
                           int ithrd1, int ithrd2);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

} // end namespace
