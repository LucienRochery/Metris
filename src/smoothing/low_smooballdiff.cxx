//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

/*
Low level routine for "direct" P1 ball smoothing.
From each (facet, metric) pair, generate remaining vertex to be unit. Then average over ball.
Simplest possible approach.
*/


#include "../smoothing/low_smooballdiff.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"
#include "../quality/low_metqua.hxx"
#include "../aux_topo.hxx"
#include "../io_libmeshb.hxx"

#include "../Optimization/opt_generic.hxx"
#include "../quality/low_metqua_d.hxx"
#include "../low_geo/ccoef.hxx"
#include "../low_geo/measure.hxx"

#include "../utils/aux_timer.hxx"
#include "../utils/mprintf.hxx"
#include "../utils/fmt_formatters.hxx"

#include <cstdlib>
#include <fstream>

namespace Metris{

// inorm <= infi norm , p > 0 L^p norm (over ball)
template<class MFT, int idim, int ideg>
int smooballdiff(Mesh<MFT>& msh, int ipoin,
                   const intAr1 &lball,
                   double*__restrict__ qnrm0, double*__restrict__ qmax0,
                   double*__restrict__ qnrm1, double*__restrict__ qmax1,
                   QuaFun iquaf){

  GETVDEPTH(msh.param);

  CPRINTF1("-- START smooballdiff ipoin = {} nball {}\n",ipoin,lball.get_n());

  constexpr int tdim = idim;
  constexpr int gdim = idim;

  constexpr int nnmet = (idim*(idim+1))/2;
  constexpr int nhess = nnmet;

  const intAr2& ent2poi = msh.ent2poi(idim);

  const auto quafun = get_quafun<MFT,gdim,tdim>(iquaf);
  const auto d_quafun = get_d_quafun<MFT,gdim,tdim>(iquaf);

  int nball = lball.get_n();

  *qnrm0 = 0;
  *qmax0 = -1.0e30;

  // Optimization doesn't reinterpolate metric
  int miter1 = MAX(1,msh.param->iflag1);
  // int miter1 = 3;

  // Relative decrease tolerance
  const double ftol = 1.0e-2;
  // const double ftol = 1.0e-3;

  newton_drivertype_args<idim> nargs(msh.param);
  #ifdef TESTQUALITYALGO
  nargs.stpmin = 1.0e-6;
  #else
  nargs.stpmin = 1.0e-6;
  #endif
  nargs.wlfc1 = 0.1;
  nargs.wlfc2 = 10.0;
  nargs.ratnew= 0.5;
  nargs.maxit = 50;
  nargs.ftol  = ftol;

  int iflag = 0, ihess, ierro = 0;
  double  xcur[idim], coor0[idim], met0[nnmet], fcur;
  bool iinva;
  double d1qua[idim], d2qua[nhess];

  for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipoin,ii);
  for(int ii = 0; ii < idim; ii++) nargs.xopt[ii] = msh.coord(ipoin,ii); // a backup in case no updates
  for(int ii = 0; ii < nnmet;ii++) met0[ii]  = msh.met(ipoin,ii);

  for(int niter1 = 0; niter1 < miter1; niter1++){

    for(int ii = 0; ii < idim; ii++) xcur[ii]  = msh.coord(ipoin,ii);
    while(true){
      INCVDEPTH(msh.param);

      // std::cout << "BEFORE optim_newton_drivertype" << std::endl;
      // std::cout << "iflag = " << iflag << std::endl;
      // for(int ii = 0; ii < idim; ii++) std::cout << "xcur[" << ii << "]  = " << xcur[ii] << std::endl;

      ierro = optim_newton_drivertype(nargs, xcur, &fcur, d1qua, d2qua, &iflag, &ihess);

      // std::cout << "AFTER optim_newton_drivertype" << std::endl;
      // std::cout << "iflag = " << iflag << std::endl;
      // for(int ii = 0; ii < idim; ii++) std::cout << "xcur[" << ii << "]  = " << xcur[ii] << std::endl;

      //if(!fpreset){
      //  fpreset = true;
      //  fpre = fcur;
      //}else{
      //  if(abs(fpre - fcur) < ftol * abs(fpre)){
      //    CPRINTF1(" - Relative decrease {:15.7e} < {:15.7e} end\n",
      //                          abs(fpre - fcur) / abs(fpre), ftol);
      //    break;
      //  }
      //  fpre = fcur;
      //}
      if(ierro > 0){
        CPRINTF1(" # optim_newton_drivertype error {}\n",ierro);
        goto finish;
      }
      if(iflag <= 0) {
        CPRINTF1(" - iflag = 0 termination\n");
        break;
      }

      for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = xcur[ii];

      iinva = false;
      if constexpr (ideg == 1){
        for(int ientt : lball){
          iinva = !isvalideltP1<gdim,tdim>(msh,ientt);
          if(iinva) break;
        }
      }else{
        constexpr int jdeg = tdim*(ideg - 1);
        constexpr int ncoef = tdim == 2 ? getnnod2(jdeg)
                                        : getnnod3(jdeg);
        double ccoef[ncoef];
        for(int ientt : lball){
          getsclccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef,&iinva);
          if(iinva) break;
        }
      }

      if(iinva){
        fcur = 1.0e10;
        // radical solution for now
        CPRINTF1("# invalid config -> finish");
        goto finish;
      }


      fcur = 0;
      double targetWeightCurrent = 0.;
      double dqelt[idim], hqelt[nhess];
      for(int ii = 0; ii < idim; ii++) d1qua[ii] = 0;
      for(int ii = 0; ii < nhess;ii++) d2qua[ii] = 0;
      for(int iball = 0; iball < nball && !iinva; iball++){
        int ient2 = lball[iball];

        // std::cout << "ient2 = " << ient2;

        bool iflat = !isvalideltP1<idim,idim>(msh,ient2);
        if(iflat){
          fcur = 1.0e10;
          break;
        }

        int ivar  = msh.template getverent<ideg>(ient2,idim,ipoin);
        double quael;
        double quaelTrue;
        if(ihess){
          quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                           ient2,ivar,
                           msh.getBasis(),
                           DifVar::None,dqelt,hqelt,1);
          quaelTrue = quafun(msh,AsDeg::Pk,AsDeg::Pk,ient2,1);

        }else{
          quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                           ient2,ivar,
                           msh.getBasis(),
                           DifVar::None,dqelt,NULL,1);
        }
        // std::cout << " q = " << quaelTrue << std::endl;
        // std::cout << "d_quafun = " << quael << ", quafun = " << quaelTrue << std::endl;
        double regionWeight = 1.;
        if(iquaf == QuaFun::StepDistance
           && msh.param->step_distance_cavity_target_average){
          regionWeight =
              step_distance_element_target_weight<MFT,idim,idim>(
                  msh,AsDeg::Pk,ient2);
          targetWeightCurrent += regionWeight;
        }

        fcur += regionWeight*quael;
        for(int ii = 0; ii < idim; ii++){
          d1qua[ii] += regionWeight*dqelt[ii];
        }
        if(ihess)
          for(int ii = 0; ii < nhess;ii++){
            d2qua[ii] += regionWeight*hqelt[ii];
          }

        if(nargs.niter == 1 && niter1 == 0){
          *qnrm0 += regionWeight*quael;
          *qmax0  = MAX(quael,*qmax0);
        }
      }// for iball
      if(iquaf == QuaFun::StepDistance
         && msh.param->step_distance_cavity_target_average){
        METRIS_ENFORCE(targetWeightCurrent > 0.);
        fcur /= targetWeightCurrent;
        for(int ii = 0; ii < idim; ii++){
          d1qua[ii] /= targetWeightCurrent;
        }
        if(ihess){
          for(int ii = 0; ii < nhess; ii++){
            d2qua[ii] /= targetWeightCurrent;
          }
        }
        if(nargs.niter == 1 && niter1 == 0){
          *qnrm0 /= targetWeightCurrent;
        }
      }
      if(DOPRINTS1()){
        CPRINTF1(" - Newton iter {} fcur = {} xcur = {} {}",nargs.niter,fcur,xcur[0],xcur[1]);
        if(idim == 3){
          PRINTF(" {}\n",xcur[2]);
        }else{
          PRINTF("\n");
        }
        CPRINTF2(" - grad = {}\n",dblAr1(idim,d1qua));
      }

    } // end while true


    ierro = 0;

    finish:
    CPRINTF1(" -- END smooballdiff fopt = {} xopt = {}\n",
             nargs.fopt,dblAr1(idim,nargs.xopt));

    for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = nargs.xopt[ii];

    if(DOPRINTS2()) writeMesh("debug_smooth0.meshb",msh);

    ierro = msh.interpMetBack(ipoin);
    if(ierro > 0){
      CPRINTF1(" # smooballdiff interpMetBack failure ierro = {} \n",ierro);
      goto cleanup;
    }


    for(int iball = 0; iball < nball; iball++){
      int ient2 = lball[iball];
      bool iflat = !isvalideltP1<idim,idim>(msh,ient2);
      METRIS_ASSERT_MSG(!iflat,"Flat iball {} elt {}", iball, ient2);
    }

    *qnrm1 = 0;
    *qmax1 = -1.0e30;
    double targetWeightFinal = 0.;
    for(int iball = 0; iball < nball; iball++){
      int ient2 = lball[iball];
      double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,
                            ient2,1);

      double regionWeight = 1.;
      if(iquaf == QuaFun::StepDistance
         && msh.param->step_distance_cavity_target_average){
        regionWeight =
            step_distance_element_target_weight<MFT,idim,idim>(
                msh,AsDeg::Pk,ient2);
        targetWeightFinal += regionWeight;
      }
      *qnrm1 += regionWeight*quael;
      *qmax1  = MAX(quael,*qmax1);
    }
    if(iquaf == QuaFun::StepDistance
       && msh.param->step_distance_cavity_target_average){
      METRIS_ENFORCE(targetWeightFinal > 0.);
      *qnrm1 /= targetWeightFinal;
    }
    CPRINTF1(" - Newton update initial quality avg {:15.7e} "
                          "max {:15.7e} \n",*qnrm0,*qmax0);
    CPRINTF1(" -                 final quality avg {:15.7e} "
                          "max {:15.7e} \n",*qnrm1,*qmax1);
  }


  if(*qnrm1 > *qnrm0){
    ierro = 2;
    CPRINTF1(" # Local smoo reject: quality norm increase "
               "{} -> {} \n", *qnrm0, *qnrm1);
    goto cleanup;
  }

  if(msh.param->dbgfull){
    for(int ientt : lball){
      if constexpr (ideg == 1){
        METRIS_ENFORCE((isvalideltP1<idim,idim>(msh,ientt)));
      }else{
        constexpr int jdeg = tdim*(ideg - 1);
        constexpr int ncoef = tdim == 2 ? getnnod2(jdeg)
                                        : getnnod3(jdeg);
        double ccoef[ncoef];
        for(int ientt : lball){
          getsclccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef,&iinva);
          METRIS_ENFORCE(!iinva);
        }
      }
    }
  }

  return 0;

  cleanup:
  for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii];
  for(int ii = 0; ii < nnmet; ii++) msh.met(ipoin,ii) = met0[ii];
  *qnrm1 = *qnrm0;
  *qmax1 = *qmax0;

  if(msh.param->dbgfull){
    if constexpr (ideg >= 2){
      constexpr int jdeg = tdim*(ideg - 1);
      constexpr int ncoef = tdim == 2 ? getnnod2(jdeg)
                                      : getnnod3(jdeg);
      const double jtol = msh.param->jtol;
      double ccoef[ncoef];
      for(int ientt : lball){
        double vol = getmeasentP1<idim>(ent2poi[ientt], msh.coord);
        getccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef);
        for(int ii = 0; ii < ncoef; ii++){
          if(ccoef[ii] >= jtol * vol) continue;
          METRIS_THROW_MSG(" - 1 reject validity coef {:15.7e} scaled {:15.7e} \n",
                  ccoef[ii], ccoef[ii]/vol);
        }
      }
    }else{
      for(int ientt : lball){
        if(isvalideltP1<idim,idim>(msh,ientt)) continue;
        METRIS_THROW_MSG(" - 2 reject validity\n");
      }
    }
  }

  return ierro;
}




#define BOOST_PP_LOCAL_MACRO(n)\
template int smooballdiff<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        >& msh,\
 int ipoin, const intAr1 &lball,\
                   double*__restrict__ qnrm0, double*__restrict__ qmax0, \
                   double*__restrict__ qnrm1, double*__restrict__ qmax1,\
                   QuaFun iquaf);\
template int smooballdiff<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        >& msh,\
 int ipoin, const intAr1 &lball,\
                   double*__restrict__ qnrm0, double*__restrict__ qmax0, \
                   double*__restrict__ qnrm1, double*__restrict__ qmax1,\
                   QuaFun iquaf);\
template int smooballdiff<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical>& msh,\
 int ipoin, const intAr1 &lball,\
                   double*__restrict__ qnrm0, double*__restrict__ qmax0, \
                   double*__restrict__ qnrm1, double*__restrict__ qmax1,\
                   QuaFun iquaf);\
template int smooballdiff<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical>& msh,\
 int ipoin, const intAr1 &lball,\
                   double*__restrict__ qnrm0, double*__restrict__ qmax0, \
                   double*__restrict__ qnrm1, double*__restrict__ qmax1,\
                   QuaFun iquaf);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

// =============================================================================================== //
// =============================================================================================== //

template<class MFT, int idim, int ideg>
int smooballdiff_boundary(Mesh<MFT>& msh, int ipoin, const int cadDim,
                          const intAr1 &lball,
                          double*__restrict__ qnrm0, double*__restrict__ qmax0,
                          double*__restrict__ qnrm1, double*__restrict__ qmax1,
                          QuaFun iquaf){

  GETVDEPTH(msh.param);

  CPRINTF1("-- START smooballdiff_boundary ipoin = {} nball {}\n",ipoin,lball.get_n());

  constexpr int tdim = idim;
  constexpr int gdim = idim;

  constexpr int nnmet = (idim*(idim+1))/2;
  constexpr int nhess = nnmet;

  const intAr2& ent2poi = msh.ent2poi(idim);

  const auto quafun = get_quafun<MFT,gdim,tdim>(iquaf);
  const auto d_quafun = get_d_quafun<MFT,gdim,tdim>(iquaf);

  int nball = lball.get_n();

  *qnrm0 = 0;
  *qmax0 = -1.0e30;

  const int ibpoin = msh.poi2bpo[ipoin];

  // a few sanity checks
  METRIS_ASSERT_MSG(ibpoin >= 0, "Interior point detected! This is boundary smoothing. Should call smooballdiff instead.");
  METRIS_ASSERT(msh.bpo2ibi(ibpoin,0) == ipoin);
  METRIS_ASSERT_MSG(msh.bpo2ibi(ibpoin,1) != 0, "CAD corner detected! Should not be moved");
  METRIS_ASSERT(cadDim == 1 || cadDim == 2);
  #ifndef NDEBUG
  if      (cadDim == 1) {METRIS_ASSERT(msh.bpo2ibi(ibpoin,1) == 1);}
  else if (cadDim == 2) {METRIS_ASSERT(msh.bpo2ibi(ibpoin,1) == 2);}
  #endif

  if (cadDim == 1){

    const int iedge = msh.bpo2ibi(ibpoin,2); // a mesh edge attached to the point
    METRIS_ASSERT(iedge >= 0 && iedge < msh.nedge);

    const int iref = msh.edg2ref[iedge]; // get CAD edge reference
    ego cadEdge = msh.CAD.cad2edg[iref]; // get actual CAD edge object, needed to compute derivatives of the edge curve wrt parameter t

    double t0 = msh.bpo2rbi(ibpoin,0); // initial t parameter for current point location

    double coor0[idim], met0[nnmet];
    for(int ii=0; ii<idim;  ii++) coor0[ii] = msh.coord(ipoin,ii);
    for(int ii=0; ii<nnmet; ii++) met0[ii]  = msh.met(ipoin,ii);

    // Optimization doesn't reinterpolate metric
    int miter1 = MAX(1,msh.param->iflag1);
    // int miter1 = 3;

    // Relative decrease tolerance
    const double ftol = 1.0e-2;
    // const double ftol = 1.0e-3;

    constexpr bool inBoundary = true;
    constexpr int optDimension = 1; // optimization has just one variable: parameter t
    newton_drivertype_args<optDimension> nargs(msh.param);
    #ifdef TESTQUALITYALGO
    nargs.stpmin = 1.0e-6;
    #else
    nargs.stpmin = 1.0e-6;
    #endif
    nargs.wlfc1 = 0.1;
    nargs.wlfc2 = 10.0;
    nargs.ratnew= 0.5;
    nargs.maxit = 50;
    nargs.ftol  = ftol;

    int iflag = 0, ihess, ierro = 0;
    bool iinva = false;

    double tcur[1];
    tcur[0] = t0;
    nargs.xopt[0] = t0; // a backup in case no updated

    double fcur    = 0.;   // objective function value
    double d1t[1]  = {0.}; // first derivative wrt t
    double d2t[1]  = {0.}; // second derivative wrt t

    // derivatives wrt physical coord
    double gradX[idim], hessX[nhess];

    double Xt[3]  = {0.,0.,0.}; // edge tangent
    double Xtt[3] = {0.,0.,0};  // t derivative of tangent vector

    // CAD evaluation buffers
    double egParam[2] = {t0, 0.};
    double evalResult[18];

    for(int niter1 = 0; niter1 < miter1; niter1++){

      while(true){
        INCVDEPTH(msh.param);

        // std::cout << "Before calling Newton, niter1 = " << niter1 << ", nargs.niter = " << nargs.niter << std::endl;
        // std::cout << std::setprecision(16) << std::scientific;
        // std::cout << "tcur = " << tcur[0] << std::endl;
        // std::cout << std::setprecision(16) << std::scientific;
        // std::cout << "fcur = " << fcur << std::endl;
        // std::cout << std::setprecision(16) << std::scientific;
        // std::cout << "d1t = "  << d1t[0] << std::endl;
        // std::cout << std::setprecision(16) << std::scientific;
        // std::cout << "d2t = "  << d2t[0] << std::endl;
        // std::cout << "ihess = " << ihess << std::endl;

        ierro = optim_newton_drivertype<optDimension,inBoundary>(nargs, tcur, &fcur, d1t, d2t, &iflag, &ihess);

        if(ierro > 0){
          CPRINTF1(" # optim_newton_drivertype error {}\n",ierro);
          goto finishdim1;
        }
        if(iflag <= 0) {
          CPRINTF1(" - iflag = 0 termination\n");
          break;
        }

        // map parameter t to physical coordinates
        egParam[0] = tcur[0];
        int estat = EG_evaluate(cadEdge, egParam, evalResult);
        METRIS_ASSERT_MSG(estat == EGADS_SUCCESS, "EG_evaluate failed estat={}", estat);

        // set point coordinates from CAD curve evaluation
        for(int ii = 0; ii < idim; ii++) msh.coord(ipoin, ii) = evalResult[ii];

        // get curve t derivatives
        Xt[0]  = evalResult[3];  Xt[1]  = evalResult[4];  Xt[2]  = (idim==3 ? evalResult[5] : 0.0);
        Xtt[0] = evalResult[6];  Xtt[1] = evalResult[7];  Xtt[2] = (idim==3 ? evalResult[8] : 0.0);

        iinva = false;
        if constexpr (ideg == 1){
          for(int ientt : lball){
            iinva = !isvalideltP1<gdim,tdim>(msh,ientt);
            if(iinva) break;
          }
        }else{
          constexpr int jdeg = tdim*(ideg - 1);
          constexpr int ncoef = tdim == 2 ? getnnod2(jdeg)
                                          : getnnod3(jdeg);
          double ccoef[ncoef];
          for(int ientt : lball){
            getsclccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef,&iinva);
            if(iinva) break;
          }
        }

        if(iinva){
          fcur = 1.0e10;
          // radical solution for now
          CPRINTF1("# invalid config -> finish");
          goto finishdim1;
        }


        fcur = 0;
        double targetWeightCurrent = 0.;
        double dqelt[idim], hqelt[nhess];
        for(int ii = 0; ii < idim; ii++) gradX[ii] = 0;
        for(int ii = 0; ii < nhess;ii++) hessX[ii] = 0;
        for(int iball = 0; iball < nball && !iinva; iball++){
          int ient2 = lball[iball];

          // std::cout << "ient2 = " << ient2;

          bool iflat = !isvalideltP1<idim,idim>(msh,ient2);
          if(iflat){
            fcur = 1.0e10;
            break;
          }

          int ivar  = msh.template getverent<ideg>(ient2,idim,ipoin);
          double quael;
          double quael_True;
          if(ihess){
            quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                            ient2,ivar,
                            msh.getBasis(),
                            DifVar::None,dqelt,hqelt,1);
            quael_True = quafun(msh,AsDeg::Pk,AsDeg::Pk,
                              ient2,1);

          }else{
            quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                            ient2,ivar,
                            msh.getBasis(),
                            DifVar::None,dqelt,NULL,1);

            quael_True = quafun(msh,AsDeg::Pk,AsDeg::Pk,
                              ient2,1);
          }
          double regionWeight = 1.;
          if(iquaf == QuaFun::StepDistance
             && msh.param->step_distance_cavity_target_average){
            regionWeight =
                step_distance_element_target_weight<MFT,idim,idim>(
                    msh,AsDeg::Pk,ient2);
            targetWeightCurrent += regionWeight;
          }
          fcur += regionWeight*quael;
          for(int ii = 0; ii < idim; ii++){
            gradX[ii] += regionWeight*dqelt[ii];
          }
          if(ihess)
            for(int ii = 0; ii < nhess;ii++){
              hessX[ii] += regionWeight*hqelt[ii];
            }

          if(nargs.niter == 1 && niter1 == 0){
            *qnrm0 += regionWeight*quael;
            *qmax0  = MAX(quael,*qmax0);
          }
        }// for iball

        if(iquaf == QuaFun::StepDistance
           && msh.param->step_distance_cavity_target_average){
          METRIS_ENFORCE(targetWeightCurrent > 0.);
          fcur /= targetWeightCurrent;
          for(int ii = 0; ii < idim; ii++) gradX[ii] /= targetWeightCurrent;
          if(ihess){
            for(int ii = 0; ii < nhess; ii++) hessX[ii] /= targetWeightCurrent;
          }
          if(nargs.niter == 1 && niter1 == 0){
            *qnrm0 /= targetWeightCurrent;
          }
        }

        if(iinva){
          CPRINTF1("# invalid/flat config -> finish");
          goto finishdim1;
        }

        // now we have the gradient nad hessian of objective function with respect to vertex position
        // apply the chain rule to get t derivatives
        // first derivative: gradX * Xt
        double dfdt = 0.;
        for(int ii = 0; ii < idim; ii++) dfdt += gradX[ii] * Xt[ii];

        // second derivative: Xt^T * hessX * Xt + gradX * Xtt
        double XtHXt = 0.0;
        if (ihess){

          auto H = [&](int i, int j) -> double {
            if(i==j){
              return hessX[i]; // 0->H00, 1->H11, 2->H22 (when idim==3)
            }
            if(idim==2){
              // nhess==3: [H00,H11,H01]
              return hessX[2];
            }else{
              // nhess==6: [H00,H11,H22,H01,H02,H12]
              if((i==0 && j==1) || (i==1 && j==0)) return hessX[3];
              if((i==0 && j==2) || (i==2 && j==0)) return hessX[4];
              /* (1,2) */                             return hessX[5];
            }
          };

          for(int i = 0; i < idim; i++){
            for(int j = 0; j < idim; j++){
              XtHXt += Xt[i] * H(i,j) * Xt[j];
            }
          }
        }

        double gradXdotXtt = 0.0;
        for(int ii = 0; ii < idim; ii++) gradXdotXtt += gradX[ii] * Xtt[ii];

        double d2fdt2 = (ihess ? XtHXt : 0.0) + gradXdotXtt;

        // Feed optimizer with derivatives in parameter space
        d1t[0] = dfdt;
        if(ihess) d2t[0] = d2fdt2;


        if(DOPRINTS1()){
          CPRINTF1(" - Smoothing on edge: Newton iter {} fcur = {} tcur = {}",nargs.niter,fcur,tcur[0]);
          CPRINTF2(" - dquadt = {}\n",d1t[0]);
        }
      } // end while true


      ierro = 0;

      finishdim1:
      // CPRINTF1(" -- END smooballdiff fopt = {} xopt = {}\n",
      //         nargs.fopt,dblAr1(idim,nargs.xopt));

      // set final t
      double topt = nargs.xopt[0];
      msh.bpo2rbi(ibpoin,0) = topt;

      // map topt -> coords
      egParam[0] = topt;
      int estat = EG_evaluate(cadEdge, egParam, evalResult);
      METRIS_ASSERT_MSG(estat == EGADS_SUCCESS, "EG_evaluate failed estat={}", estat);

      for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = evalResult[ii];

      if(DOPRINTS2()) writeMesh("debug_smooth0.meshb",msh);

      ierro = msh.interpMetBack(ipoin);
      if(ierro > 0){
        CPRINTF1(" # smooballdiff interpMetBack failure ierro = {} \n",ierro);
        goto cleanupdim1;
      }


      for(int iball = 0; iball < nball; iball++){
        int ient2 = lball[iball];
        bool iflat = !isvalideltP1<idim,idim>(msh,ient2);
        METRIS_ASSERT_MSG(!iflat,"Flat iball {} elt {}", iball, ient2);
      }

      *qnrm1 = 0;
      *qmax1 = -1.0e30;
      double targetWeightFinal = 0.;
      for(int iball = 0; iball < nball; iball++){
        int ient2 = lball[iball];
        double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,
                              ient2,1);

        double regionWeight = 1.;
        if(iquaf == QuaFun::StepDistance
           && msh.param->step_distance_cavity_target_average){
          regionWeight =
              step_distance_element_target_weight<MFT,idim,idim>(
                  msh,AsDeg::Pk,ient2);
          targetWeightFinal += regionWeight;
        }
        *qnrm1 += regionWeight*quael;
        *qmax1  = MAX(quael,*qmax1);
      }
      if(iquaf == QuaFun::StepDistance
         && msh.param->step_distance_cavity_target_average){
        METRIS_ENFORCE(targetWeightFinal > 0.);
        *qnrm1 /= targetWeightFinal;
      }
      CPRINTF1(" - Newton update initial quality avg {:15.7e} "
                            "max {:15.7e} \n",*qnrm0,*qmax0);
      CPRINTF1(" -                 final quality avg {:15.7e} "
                            "max {:15.7e} \n",*qnrm1,*qmax1);
    }


    if(*qnrm1 > *qnrm0){
      ierro = 2;
      CPRINTF1(" # Local smoo reject: quality norm increase "
                "{} -> {} \n", *qnrm0, *qnrm1);
      goto cleanupdim1;
    }

    if(msh.param->dbgfull){
      for(int ientt : lball){
        if constexpr (ideg == 1){
          METRIS_ENFORCE((isvalideltP1<idim,idim>(msh,ientt)));
        }else{
          constexpr int jdeg = tdim*(ideg - 1);
          constexpr int ncoef = tdim == 2 ? getnnod2(jdeg)
                                          : getnnod3(jdeg);
          double ccoef[ncoef];
          for(int ientt : lball){
            getsclccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef,&iinva);
            METRIS_ENFORCE(!iinva);
          }
        }
      }
    }

    return 0;

    cleanupdim1:
    for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii];
    for(int ii = 0; ii < nnmet; ii++) msh.met(ipoin,ii) = met0[ii];
    msh.bpo2rbi(ibpoin,0) = t0;

    *qnrm1 = *qnrm0;
    *qmax1 = *qmax0;

    if(msh.param->dbgfull){
      if constexpr (ideg >= 2){
        constexpr int jdeg = tdim*(ideg - 1);
        constexpr int ncoef = tdim == 2 ? getnnod2(jdeg)
                                        : getnnod3(jdeg);
        const double jtol = msh.param->jtol;
        double ccoef[ncoef];
        for(int ientt : lball){
          double vol = getmeasentP1<idim>(ent2poi[ientt], msh.coord);
          getccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef);
          for(int ii = 0; ii < ncoef; ii++){
            if(ccoef[ii] >= jtol * vol) continue;
            METRIS_THROW_MSG(" - 1 reject validity coef {:15.7e} scaled {:15.7e} \n",
                    ccoef[ii], ccoef[ii]/vol);
          }
        }
      }else{
        for(int ientt : lball){
          if(isvalideltP1<idim,idim>(msh,ientt)) continue;
          METRIS_THROW_MSG(" - 2 reject validity\n");
        }
      }
    }

    return ierro;
  } // if cadDim == 1 (point on CAD edge)
  else if (cadDim == 2){

    const int iface = msh.bpo2ibi(ibpoin,2); // a mesh face attached to the point
    METRIS_ASSERT(iface >= 0 && iface < msh.nface);

    const int iref = msh.fac2ref[iface]; // get CAD face reference
    ego cadFace = msh.CAD.cad2fac[iref]; // get actual CAD face object, needed to compute derivatives of the face surface wrt parameters (u,v)

    double u0 = msh.bpo2rbi(ibpoin,0); // initial u parameter for current point location
    double v0 = msh.bpo2rbi(ibpoin,1); // initial v parameter for current point location

    double coor0[idim], met0[nnmet];
    for(int ii=0; ii<idim;  ii++) coor0[ii] = msh.coord(ipoin,ii);
    for(int ii=0; ii<nnmet; ii++) met0[ii]  = msh.met(ipoin,ii);

    // Optimization doesn't reinterpolate metric
    int miter1 = MAX(1,msh.param->iflag1);
    // int miter1 = 3;

    // Relative decrease tolerance
    const double ftol = 1.0e-2;
    // const double ftol = 1.0e-3;

    constexpr bool inBoundary = true;
    constexpr int optDimension = 2; // optimization has just two variables: parameters (u,v)
    newton_drivertype_args<optDimension> nargs(msh.param);
    #ifdef TESTQUALITYALGO
    nargs.stpmin = 1.0e-6;
    #else
    nargs.stpmin = 1.0e-6;
    #endif
    nargs.wlfc1 = 0.1;
    nargs.wlfc2 = 10.0;
    nargs.ratnew= 0.5;
    nargs.maxit = 50;
    nargs.ftol  = ftol;

    int iflag = 0, ihess, ierro = 0;
    bool iinva = false;

    double uvcur[2];
    uvcur[0] = u0;
    uvcur[1] = v0;
    nargs.xopt[0] = u0; // a backup in case no updated
    nargs.xopt[1] = v0; // a backup in case no updated

    double fcur    = 0.;   // objective function value
    double gradUV[2]  = {0., 0.};     // gradient wrt uv
    double hessUV[3]  = {0., 0., 0.}; // hessian wrt uv

    // derivatives wrt physical coord
    double gradX[idim], hessX[nhess];

    double Xu[3]  = {0.,0.,0.}; // surface tangent in u direction
    double Xuu[3] = {0.,0.,0};
    double Xv[3]  = {0.,0.,0.}; // surface tangent in v direction
    double Xvv[3] = {0.,0.,0};
    double Xuv[3] = {0.,0.,0};

    // CAD evaluation buffers
    double facParam[2] = {u0, v0};
    double evalResult[18];

    for(int niter1 = 0; niter1 < miter1; niter1++){

      while(true){
        INCVDEPTH(msh.param);

        // std::cout << "Before calling Newton, niter1 = " << niter1 << ", nargs.niter = " << nargs.niter << std::endl;
        // std::cout << std::setprecision(16) << std::scientific;
        // std::cout << "ucur = " << uvcur[0] << std::endl;
        // std::cout << std::setprecision(16) << std::scientific;
        // std::cout << "vcur = " << uvcur[1] << std::endl;
        // std::cout << std::setprecision(16) << std::scientific;
        // std::cout << "fcur = " << fcur << std::endl;
        // std::cout << std::setprecision(16) << std::scientific;
        // std::cout << "gradUV[0] = "  << gradUV[0] << std::endl;
        // std::cout << std::setprecision(16) << std::scientific;
        // std::cout << "gradUV[1] = "  << gradUV[1] << std::endl;
        // std::cout << std::setprecision(16) << std::scientific;
        // std::cout << "hessUV[00] = "  << hessUV[0] << std::endl;
        // std::cout << std::setprecision(16) << std::scientific;
        // std::cout << "hessUV[01] = "  << hessUV[1] << std::endl;
        // std::cout << std::setprecision(16) << std::scientific;
        // std::cout << "hessUV[11] = "  << hessUV[2] << std::endl;
        // std::cout << "ihess = " << ihess << std::endl;

        ierro = optim_newton_drivertype<optDimension,inBoundary>(nargs, uvcur, &fcur, gradUV, hessUV, &iflag, &ihess);

        if(ierro > 0){
          CPRINTF1(" # optim_newton_drivertype error {}\n",ierro);
          goto finishdim2;
        }
        if(iflag <= 0) {
          CPRINTF1(" - iflag = 0 termination\n");
          break;
        }

        // map parameter t to physical coordinates
        facParam[0] = uvcur[0];
        facParam[1] = uvcur[1];
        int estat = EG_evaluate(cadFace, facParam, evalResult);
        METRIS_ASSERT_MSG(estat == EGADS_SUCCESS, "EG_evaluate failed estat={}", estat);

        // set point coordinates from CAD surface evaluation
        for(int ii = 0; ii < idim; ii++) msh.coord(ipoin, ii) = evalResult[ii];

        // get surface uv derivatives
        Xu[0]  = evalResult[3];  Xu[1]  = evalResult[4];  Xu[2]  = (idim==3 ? evalResult[5] : 0.0);
        Xv[0]  = evalResult[6];  Xv[1]  = evalResult[7];  Xv[2]  = (idim==3 ? evalResult[8] : 0.0);

        // get surface uv second derivatives
        Xuu[0] = evalResult[9];   Xuu[1] = evalResult[10];  Xuu[2] = (idim==3 ? evalResult[11] : 0.0);
        Xuv[0] = evalResult[12];  Xuv[1] = evalResult[13];  Xuv[2] = (idim==3 ? evalResult[14] : 0.0);
        Xvv[0] = evalResult[15];  Xvv[1] = evalResult[16];  Xvv[2] = (idim==3 ? evalResult[17] : 0.0);


        iinva = false;
        if constexpr (ideg == 1){
          for(int ientt : lball){
            iinva = !isvalideltP1<gdim,tdim>(msh,ientt);
            if(iinva) break;
          }
        }else{
          constexpr int jdeg = tdim*(ideg - 1);
          constexpr int ncoef = tdim == 2 ? getnnod2(jdeg)
                                          : getnnod3(jdeg);
          double ccoef[ncoef];
          for(int ientt : lball){
            getsclccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef,&iinva);
            if(iinva) break;
          }
        }

        if(iinva){
          fcur = 1.0e10;
          // radical solution for now
          CPRINTF1("# invalid config -> finish");
          goto finishdim2;
        }


        fcur = 0;
        double targetWeightCurrent = 0.;
        double dqelt[idim], hqelt[nhess];
        for(int ii = 0; ii < idim; ii++) gradX[ii] = 0;
        for(int ii = 0; ii < nhess;ii++) hessX[ii] = 0;
        for(int iball = 0; iball < nball && !iinva; iball++){
          int ient2 = lball[iball];

          // std::cout << "ient2 = " << ient2;

          bool iflat = !isvalideltP1<idim,idim>(msh,ient2);
          if(iflat){
            fcur = 1.0e10;
            break;
          }

          int ivar  = msh.template getverent<ideg>(ient2,idim,ipoin);
          double quael;
          if(ihess){
            quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                            ient2,ivar,
                            msh.getBasis(),
                            DifVar::None,dqelt,hqelt,1);

          }else{
            quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                            ient2,ivar,
                            msh.getBasis(),
                            DifVar::None,dqelt,NULL,1);
          }
          double regionWeight = 1.;
          if(iquaf == QuaFun::StepDistance
             && msh.param->step_distance_cavity_target_average){
            regionWeight =
                step_distance_element_target_weight<MFT,idim,idim>(
                    msh,AsDeg::Pk,ient2);
            targetWeightCurrent += regionWeight;
          }
          fcur += regionWeight*quael;
          for(int ii = 0; ii < idim; ii++){
            gradX[ii] += regionWeight*dqelt[ii];
          }
          if(ihess)
            for(int ii = 0; ii < nhess;ii++){
              hessX[ii] += regionWeight*hqelt[ii];
            }

          if(nargs.niter == 1 && niter1 == 0){
            *qnrm0 += regionWeight*quael;
            *qmax0  = MAX(quael,*qmax0);
          }
        }// for iball

        if(iquaf == QuaFun::StepDistance
           && msh.param->step_distance_cavity_target_average){
          METRIS_ENFORCE(targetWeightCurrent > 0.);
          fcur /= targetWeightCurrent;
          for(int ii = 0; ii < idim; ii++) gradX[ii] /= targetWeightCurrent;
          if(ihess){
            for(int ii = 0; ii < nhess; ii++) hessX[ii] /= targetWeightCurrent;
          }
          if(nargs.niter == 1 && niter1 == 0){
            *qnrm0 /= targetWeightCurrent;
          }
        }

        if(iinva){
          CPRINTF1("# invalid/flat config -> finish");
          goto finishdim2;
        }

        // now we have the gradient nad hessian of objective function with respect to vertex position
        // apply the chain rule to get t derivatives
        // first derivative: gradX * Xt
        double dfdu = 0.;
        double dfdv = 0.;
        for(int ii = 0; ii < idim; ii++){

          dfdu += gradX[ii] * Xu[ii];
          dfdv += gradX[ii] * Xv[ii];
        }

        // second derivative wrt u: Xu^T * hessX * Xu + gradX * Xuu
        // second derivative wrt v: Xv^T * hessX * Xv + gradX * Xvv
        // mixed second derivative: Xu^T * hessX * Xv + gradX * Xuv
        double XuHXu = 0.0;
        double XvHXv = 0.0;
        double XuHXv = 0.0;
        if (ihess){

          auto H = [&](int i, int j) -> double {
            if(i==j){
              return hessX[i]; // 0->H00, 1->H11, 2->H22 (when idim==3)
            }
            if(idim==2){
              // nhess==3: [H00,H11,H01]
              return hessX[2];
            }else{
              // nhess==6: [H00,H11,H22,H01,H02,H12]
              if((i==0 && j==1) || (i==1 && j==0)) return hessX[3];
              if((i==0 && j==2) || (i==2 && j==0)) return hessX[4];
              /* (1,2) */                             return hessX[5];
            }
          };

          for(int i = 0; i < idim; i++){
            for(int j = 0; j < idim; j++){
              XuHXu += Xu[i] * H(i,j) * Xu[j];
              XvHXv += Xv[i] * H(i,j) * Xv[j];
              XuHXv += Xu[i] * H(i,j) * Xv[j];
            }
          }
        }

        double gradXdotXuu = 0.0;
        double gradXdotXvv = 0.0;
        double gradXdotXuv = 0.0;
        for(int ii = 0; ii < idim; ii++){

          gradXdotXuu += gradX[ii] * Xuu[ii];
          gradXdotXvv += gradX[ii] * Xvv[ii];
          gradXdotXuv += gradX[ii] * Xuv[ii];
        }

        double d2fdu2 = (ihess ? XuHXu : 0.0) + gradXdotXuu;
        double d2fdv2 = (ihess ? XvHXv : 0.0) + gradXdotXvv;
        double d2fduv = (ihess ? XuHXv : 0.0) + gradXdotXuv;

        // Feed optimizer with derivatives in parameter space
        gradUV[0] = dfdu;
        gradUV[1] = dfdv;
        if(ihess){

          hessUV[0] = d2fdu2;
          hessUV[1] = d2fduv;
          hessUV[2] = d2fdv2;
        }


        // if(DOPRINTS1()){
        //   CPRINTF1(" - Smoothing on edge: Newton iter {} fcur = {} tcur = {}",nargs.niter,fcur,tcur[0]);
        //   CPRINTF2(" - dquadt = {}\n",d1t[0]);
        // }
      } // end while true


      ierro = 0;

      finishdim2:
      // CPRINTF1(" -- END smooballdiff fopt = {} xopt = {}\n",
      //         nargs.fopt,dblAr1(idim,nargs.xopt));

      // set final (u,v)
      double uopt = nargs.xopt[0];
      double vopt = nargs.xopt[1];
      msh.bpo2rbi(ibpoin,0) = uopt;
      msh.bpo2rbi(ibpoin,1) = vopt;

      // map (u,v) opt -> coords
      facParam[0] = uopt;
      facParam[1] = vopt;
      int estat = EG_evaluate(cadFace, facParam, evalResult);
      METRIS_ASSERT_MSG(estat == EGADS_SUCCESS, "EG_evaluate failed estat={}", estat);

      for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = evalResult[ii];

      if(DOPRINTS2()) writeMesh("debug_smooth0.meshb",msh);

      ierro = msh.interpMetBack(ipoin);
      if(ierro > 0){
        CPRINTF1(" # smooballdiff interpMetBack failure ierro = {} \n",ierro);
        goto cleanupdim2;
      }

      for(int iball = 0; iball < nball; iball++){
        int ient2 = lball[iball];
        bool iflat = !isvalideltP1<idim,idim>(msh,ient2);
        METRIS_ASSERT_MSG(!iflat,"Flat iball {} elt {}", iball, ient2);
      }

      *qnrm1 = 0;
      *qmax1 = -1.0e30;
      double targetWeightFinal = 0.;
      for(int iball = 0; iball < nball; iball++){
        int ient2 = lball[iball];
        double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,
                              ient2,1);

        double regionWeight = 1.;
        if(iquaf == QuaFun::StepDistance
           && msh.param->step_distance_cavity_target_average){
          regionWeight =
              step_distance_element_target_weight<MFT,idim,idim>(
                  msh,AsDeg::Pk,ient2);
          targetWeightFinal += regionWeight;
        }
        *qnrm1 += regionWeight*quael;
        *qmax1  = MAX(quael,*qmax1);
      }
      if(iquaf == QuaFun::StepDistance
         && msh.param->step_distance_cavity_target_average){
        METRIS_ENFORCE(targetWeightFinal > 0.);
        *qnrm1 /= targetWeightFinal;
      }
      CPRINTF1(" - Newton update initial quality avg {:15.7e} "
                            "max {:15.7e} \n",*qnrm0,*qmax0);
      CPRINTF1(" -                 final quality avg {:15.7e} "
                            "max {:15.7e} \n",*qnrm1,*qmax1);
    }


    if(*qnrm1 > *qnrm0){
      ierro = 2;
      CPRINTF1(" # Local smoo reject: quality norm increase "
                "{} -> {} \n", *qnrm0, *qnrm1);
      goto cleanupdim2;
    }

    if(msh.param->dbgfull){
      for(int ientt : lball){
        if constexpr (ideg == 1){
          METRIS_ENFORCE((isvalideltP1<idim,idim>(msh,ientt)));
        }else{
          constexpr int jdeg = tdim*(ideg - 1);
          constexpr int ncoef = tdim == 2 ? getnnod2(jdeg)
                                          : getnnod3(jdeg);
          double ccoef[ncoef];
          for(int ientt : lball){
            getsclccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef,&iinva);
            METRIS_ENFORCE(!iinva);
          }
        }
      }
    }

    return 0;

    cleanupdim2:
    for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii];
    for(int ii = 0; ii < nnmet; ii++) msh.met(ipoin,ii) = met0[ii];
    msh.bpo2rbi(ibpoin,0) = u0;
    msh.bpo2rbi(ibpoin,1) = v0;

    *qnrm1 = *qnrm0;
    *qmax1 = *qmax0;

    if(msh.param->dbgfull){
      if constexpr (ideg >= 2){
        constexpr int jdeg = tdim*(ideg - 1);
        constexpr int ncoef = tdim == 2 ? getnnod2(jdeg)
                                        : getnnod3(jdeg);
        const double jtol = msh.param->jtol;
        double ccoef[ncoef];
        for(int ientt : lball){
          double vol = getmeasentP1<idim>(ent2poi[ientt], msh.coord);
          getccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef);
          for(int ii = 0; ii < ncoef; ii++){
            if(ccoef[ii] >= jtol * vol) continue;
            METRIS_THROW_MSG(" - 1 reject validity coef {:15.7e} scaled {:15.7e} \n",
                    ccoef[ii], ccoef[ii]/vol);
          }
        }
      }else{
        for(int ientt : lball){
          if(isvalideltP1<idim,idim>(msh,ientt)) continue;
          METRIS_THROW_MSG(" - 2 reject validity\n");
        }
      }
    }

    return ierro;
  }
  return 1;
}




#define BOOST_PP_LOCAL_MACRO(n)\
template int smooballdiff_boundary<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        >& msh,\
 int ipoin, const int cadDim, const intAr1 &lball,\
                   double*__restrict__ qnrm0, double*__restrict__ qmax0, \
                   double*__restrict__ qnrm1, double*__restrict__ qmax1,\
                   QuaFun iquaf);\
template int smooballdiff_boundary<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        >& msh,\
 int ipoin, const int cadDim,const intAr1 &lball,\
                   double*__restrict__ qnrm0, double*__restrict__ qmax0, \
                   double*__restrict__ qnrm1, double*__restrict__ qmax1,\
                   QuaFun iquaf);\
template int smooballdiff_boundary<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical>& msh,\
 int ipoin, const int cadDim, const intAr1 &lball,\
                   double*__restrict__ qnrm0, double*__restrict__ qmax0, \
                   double*__restrict__ qnrm1, double*__restrict__ qmax1,\
                   QuaFun iquaf);\
template int smooballdiff_boundary<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical>& msh,\
 int ipoin, const int cadDim, const intAr1 &lball,\
                   double*__restrict__ qnrm0, double*__restrict__ qmax0, \
                   double*__restrict__ qnrm1, double*__restrict__ qmax1,\
                   QuaFun iquaf);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

// =============================================================================================== //
// =============================================================================================== //

template<class MFT, int idim, int ideg>
int smoocavdiff(Mesh<MFT>& msh, MshCavity& cav,
                   double& quaCav1, double& quaMaxCav1,
                   double& targetWeightCav1,
                   QuaFun iquaf, int ithread){

  GETVDEPTH(msh.param);

  constexpr int tdim = idim;
  constexpr int gdim = idim;

  const intAr1& lcent = cav.lcent(tdim);
  const intAr2& ent2ent = msh.ent2ent(tdim);
  const int tmpEntt = msh.nentt(tdim)-1;
        intAr2& ent2poi = msh.ent2poi(tdim);
        intAr2& ent2tag = msh.ent2tag(tdim);

  const int ipins = cav.ipins;
  targetWeightCav1 = 0.;

  CPRINTF1("-- START smoocavdiff ipins = {} ncav {}\n",ipins,lcent.get_n());

  msh.tag[ithread]++;
  for (int ienttCav : lcent) ent2tag(ithread,ienttCav) = msh.tag[ithread];

  constexpr int nnmet = (idim*(idim+1))/2;
  constexpr int nhess = nnmet;

  const auto quafun = get_quafun<MFT,gdim,tdim>(iquaf);
  const auto d_quafun = get_d_quafun<MFT,gdim,tdim>(iquaf);

  int nenttCav = lcent.get_n();

  // Optimization doesn't reinterpolate metric
  int miter1 = MAX(1,msh.param->iflag1);
  // int miter1 = 3;

  // Relative decrease tolerance
  const double ftol = 1.0e-2;
  // const double ftol = 1.0e-3;

  newton_drivertype_args<idim> nargs(msh.param);
  #ifdef TESTQUALITYALGO
  // nargs.stpmin = 1.;
  nargs.stpmin = 1.0e-6;
  #else
  nargs.stpmin = 1.0e-6;
  #endif
  nargs.wlfc1 = 0.1;
  nargs.wlfc2 = 10.0;
  nargs.ratnew= 0.5;
  nargs.maxit = 50;
  nargs.ftol  = ftol;

  int iflag = 0, ihess, ierro = 0;
  double  xcur[idim], coor0[idim], met0[nnmet], fcur;
  bool iinva;
  double d1qua[idim], d2qua[nhess];

  for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipins,ii);
  for(int ii = 0; ii < idim; ii++) nargs.xopt[ii] = msh.coord(ipins,ii); // a backup in case no updates
  for(int ii = 0; ii < nnmet;ii++) met0[ii]  = msh.met(ipins,ii);

  for(int niter1 = 0; niter1 < miter1; niter1++){

    for(int ii = 0; ii < idim; ii++) xcur[ii]  = msh.coord(ipins,ii);
    while(true){
      INCVDEPTH(msh.param);

      ierro = optim_newton_drivertype(nargs, xcur, &fcur, d1qua, d2qua, &iflag, &ihess);

      if(ierro > 0){
        CPRINTF1(" # optim_newton_drivertype error {}\n",ierro);
        goto finish;
      }
      if(iflag <= 0) {
        CPRINTF1(" - iflag = 0 termination\n");
        break;
      }

      for(int ii = 0; ii < idim; ii++) msh.coord(ipins,ii) = xcur[ii];

      // reconnect on the fly, check validity and fill quality value and derivatives
      iinva = false;
      fcur = 0;
      double targetWeightCurrent = 0.;
      double dqelt[idim], hqelt[nhess];
      for(int ii = 0; ii < idim; ii++) d1qua[ii] = 0;
      for(int ii = 0; ii < nhess;ii++) d2qua[ii] = 0;
      if constexpr (ideg == 1){

        for (const int ienttCav : lcent){

          int ent2pol[4];
          for(int jj = 0; jj < tdim + 1; jj++){

            const int ienei = ent2ent(ienttCav,jj);

            // if neighbor tagged it is in cavity -> skip
            if (ienei >= 0 && ent2tag(ithread,ienei) >= msh.tag[ithread]) continue;

            // at this point, facet jj is at boundary of the cavity
            // need to probe quality of the reconnected element

            if constexpr (tdim == 2){

              // TODO: for now this function should not be called for boundary points
              // if (ipinsOnEdge){

              //   // check that facet jj, if also on boundary, is not in same boundary as the insertion edge. otherwise the new triangle would be flat
              //   int iedgeGlobal = msh.fac2edg(ienttCav,jj);
              //   if (iedgeGlobal >= 0){
              //     if (msh.edg2ref[iedgeGlobal] == msh.edg2ref[iedins]) continue;
              //   }
              // }

              // put together new triangle
              ent2pol[0] = ipins;
              ent2pol[lnoed2[0][0]] = ent2poi(ienttCav,lnoed2[jj][0]);
              ent2pol[lnoed2[0][1]] = ent2poi(ienttCav,lnoed2[jj][1]);

              if (ent2pol[1] == ipins || ent2pol[2] == ipins) continue;

              ent2poi(tmpEntt,0) = ent2pol[0];
              ent2poi(tmpEntt,1) = ent2pol[1];
              ent2poi(tmpEntt,2) = ent2pol[2];

              double meas;
              if (!isvalideltP1<2,2>(msh, tmpEntt, NULL, &meas)){

                iinva = true;
                break;
              }

              int ivar  = msh.template getverent<ideg>(tmpEntt,idim,ipins);
              double quael;
              if(ihess){
                quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                                tmpEntt,ivar,
                                msh.getBasis(),
                                DifVar::None,dqelt,hqelt,1.);
              }else{
                quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                                tmpEntt,ivar,
                                msh.getBasis(),
                                DifVar::None,dqelt,NULL,1.);
              }
              double regionWeight = 1.;
              if(iquaf == QuaFun::StepDistance
                 && msh.param->step_distance_cavity_target_average){
                regionWeight =
                    step_distance_element_target_weight<MFT,idim,idim>(
                        msh,AsDeg::Pk,tmpEntt);
                targetWeightCurrent += regionWeight;
              }
              fcur += regionWeight*quael;
              if (quael > quaMaxCav1) quaMaxCav1 = quael;

              for(int ii = 0; ii < idim; ii++){
                d1qua[ii] += regionWeight*dqelt[ii];
              }
              if(ihess)
                for(int ii = 0; ii < nhess;ii++){
                  d2qua[ii] += regionWeight*hqelt[ii];
                }

            }
            else{ // tdim == 3

              // TODO: for now, this function should not be called if ipins is on boundary
              // // if boundary face itself is in cavity (tagged), it will be split => no single cone tet
              // int ifaceGlobal = msh.tet2fac(ienttCav, jj);
              // if(ifaceGlobal >= 0 && msh.fac2tag(ithread, ifaceGlobal) >= msh.tag[ithread]) continue;

              int ent2pol[4];
              ent2pol[0] = ipins;
              ent2pol[lnofa3[0][0]] = ent2poi(ienttCav, lnofa3[jj][0]);
              ent2pol[lnofa3[0][1]] = ent2poi(ienttCav, lnofa3[jj][1]);
              ent2pol[lnofa3[0][2]] = ent2poi(ienttCav, lnofa3[jj][2]);

              if(ent2pol[1]==ipins || ent2pol[2]==ipins || ent2pol[3]==ipins) continue;

              // copy into tmpEntt for metqua
              ent2poi(tmpEntt,0)=ent2pol[0];
              ent2poi(tmpEntt,1)=ent2pol[1];
              ent2poi(tmpEntt,2)=ent2pol[2];
              ent2poi(tmpEntt,3)=ent2pol[3];

              // put this for debug build anyways
              double meas;
              if (!isvalideltP1<3,3>(msh, tmpEntt, NULL, &meas)){

                iinva = true;
                break;
              }

              int ivar  = msh.template getverent<ideg>(tmpEntt,idim,ipins);
              double quael;
              if(ihess){
                quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                                tmpEntt,ivar,
                                msh.getBasis(),
                                DifVar::None,dqelt,hqelt,1.);
              }else{
                quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                                tmpEntt,ivar,
                                msh.getBasis(),
                                DifVar::None,dqelt,NULL,1.);
              }
              double regionWeight = 1.;
              if(iquaf == QuaFun::StepDistance
                 && msh.param->step_distance_cavity_target_average){
                regionWeight =
                    step_distance_element_target_weight<MFT,idim,idim>(
                        msh,AsDeg::Pk,tmpEntt);
                targetWeightCurrent += regionWeight;
              }
              fcur += regionWeight*quael;
              if (quael > quaMaxCav1) quaMaxCav1 = quael;

              for(int ii = 0; ii < idim; ii++){
                d1qua[ii] += regionWeight*dqelt[ii];
              }
              if(ihess)
                for(int ii = 0; ii < nhess;ii++){
                  d2qua[ii] += regionWeight*hqelt[ii];
                }

            } // if tdim == 2 else 3
          } // for jj (bnd facets of ienttCav)
          if (iinva) break;
        } // for ienttCav
      }else{
        METRIS_THROW_MSG("Implement cavity smoothing for ideg > 1");
      }

      if(!iinva && iquaf == QuaFun::StepDistance
         && msh.param->step_distance_cavity_target_average){
        METRIS_ENFORCE(targetWeightCurrent > 0.);
        fcur /= targetWeightCurrent;
        for(int ii = 0; ii < idim; ii++){
          d1qua[ii] /= targetWeightCurrent;
        }
        if(ihess){
          for(int ii = 0; ii < nhess; ii++){
            d2qua[ii] /= targetWeightCurrent;
          }
        }
      }

      if(iinva){
        fcur = 1e15;
        quaCav1 = fcur;
        quaMaxCav1 = 1e15;
        // radical solution for now
        CPRINTF1("# invalid config -> finish");
        goto finish;
      }

      if(DOPRINTS1()){
        CPRINTF1(" - Newton iter {} fcur = {} xcur = {} {}",nargs.niter,fcur,xcur[0],xcur[1]);
        if(idim == 3){
          PRINTF(" {}\n",xcur[2]);
        }else{
          PRINTF("\n");
        }
        CPRINTF2(" - grad = {}\n",dblAr1(idim,d1qua));
      }

      quaCav1 = fcur;
    } // end while true

    ierro = 0;

    finish:
    CPRINTF1(" -- END smoocavdiff fopt = {} xopt = {}\n",
             nargs.fopt,dblAr1(idim,nargs.xopt));

    for(int ii = 0; ii < idim; ii++) msh.coord(ipins,ii) = nargs.xopt[ii];

    if(DOPRINTS2()) writeMesh("debug_smooth0.meshb",msh);

    ierro = msh.interpMetBack(ipins);
    if(ierro > 0){
      CPRINTF1(" # smoocavdiff interpMetBack failure ierro = {} \n",ierro);
      goto cleanup;
    }

    // Recompute after restoring xopt.
    quaCav1 = 0;
    quaMaxCav1 = -1.0e30;
    double targetWeightFinal = 0.;
    iinva = false;
    for(const int ienttCav : lcent){
      int ent2pol[4];
      for(int jj = 0; jj < tdim + 1; jj++){
        const int ienei = ent2ent(ienttCav,jj);
        if(ienei >= 0 && ent2tag(ithread,ienei) >= msh.tag[ithread]) continue;

        if constexpr (tdim == 2){
          ent2pol[0] = ipins;
          ent2pol[lnoed2[0][0]] = ent2poi(ienttCav,lnoed2[jj][0]);
          ent2pol[lnoed2[0][1]] = ent2poi(ienttCav,lnoed2[jj][1]);

          if(ent2pol[1] == ipins || ent2pol[2] == ipins) continue;

          ent2poi(tmpEntt,0) = ent2pol[0];
          ent2poi(tmpEntt,1) = ent2pol[1];
          ent2poi(tmpEntt,2) = ent2pol[2];

          double meas;
          if(!isvalideltP1<2,2>(msh, tmpEntt, NULL, &meas)){
            iinva = true;
            break;
          }
        }else{
          ent2pol[0] = ipins;
          ent2pol[lnofa3[0][0]] = ent2poi(ienttCav, lnofa3[jj][0]);
          ent2pol[lnofa3[0][1]] = ent2poi(ienttCav, lnofa3[jj][1]);
          ent2pol[lnofa3[0][2]] = ent2poi(ienttCav, lnofa3[jj][2]);

          if(ent2pol[1] == ipins || ent2pol[2] == ipins || ent2pol[3] == ipins) continue;

          ent2poi(tmpEntt,0) = ent2pol[0];
          ent2poi(tmpEntt,1) = ent2pol[1];
          ent2poi(tmpEntt,2) = ent2pol[2];
          ent2poi(tmpEntt,3) = ent2pol[3];

          double meas;
          if(!isvalideltP1<3,3>(msh, tmpEntt, NULL, &meas)){
            iinva = true;
            break;
          }
        }

        double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,tmpEntt,1.);
        double regionWeight = 1.;
        if(iquaf == QuaFun::StepDistance
           && msh.param->step_distance_cavity_target_average){
          regionWeight =
              step_distance_element_target_weight<MFT,idim,idim>(
                  msh,AsDeg::Pk,tmpEntt);
          targetWeightFinal += regionWeight;
        }
        quaCav1 += regionWeight*quael;
        quaMaxCav1 = MAX(quaMaxCav1, quael);
      }
      if(iinva) break;
    }
    if(iinva){
      CPRINTF1(" # smoocavdiff final xopt quality recompute invalid\n");
      goto cleanup;
    }
    if(iquaf == QuaFun::StepDistance
       && msh.param->step_distance_cavity_target_average){
      METRIS_ENFORCE(targetWeightFinal > 0.);
      targetWeightCav1 = targetWeightFinal;
    }

    // CPRINTF1(" - Newton update initial quality avg {:15.7e} "
    //                       "max {:15.7e} \n",*qnrm0,*qmax0);
    // CPRINTF1(" -                 final quality avg {:15.7e} "
    //                       "max {:15.7e} \n",*qnrm1,*qmax1);
  }


  // if(*qnrm1 > *qnrm0){
  //   ierro = 2;
  //   CPRINTF1(" # Local smoo reject: quality norm increase "
  //              "{} -> {} \n", *qnrm0, *qnrm1);
  //   goto cleanup;
  // }

  // if(msh.param->dbgfull){
  //   for(int ientt : lball){
  //     if constexpr (ideg == 1){
  //       METRIS_ENFORCE((isvalideltP1<idim,idim>(msh,ientt)));
  //     }else{
  //       constexpr int jdeg = tdim*(ideg - 1);
  //       constexpr int ncoef = tdim == 2 ? getnnod2(jdeg)
  //                                       : getnnod3(jdeg);
  //       double ccoef[ncoef];
  //       for(int ientt : lball){
  //         getsclccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef,&iinva);
  //         METRIS_ENFORCE(!iinva);
  //       }
  //     }
  //   }
  // }

  return 0;

  cleanup:
  for(int ii = 0; ii < idim; ii++) msh.coord(ipins,ii) = coor0[ii];
  for(int ii = 0; ii < nnmet; ii++) msh.met(ipins,ii) = met0[ii];
  // *qnrm1 = *qnrm0;
  // *qmax1 = *qmax0;

  // if(msh.param->dbgfull){
  //   if constexpr (ideg >= 2){
  //     constexpr int jdeg = tdim*(ideg - 1);
  //     constexpr int ncoef = tdim == 2 ? getnnod2(jdeg)
  //                                     : getnnod3(jdeg);
  //     const double jtol = msh.param->jtol;
  //     double ccoef[ncoef];
  //     for(int ientt : lball){
  //       double vol = getmeasentP1<idim>(ent2poi[ientt], msh.coord);
  //       getccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef);
  //       for(int ii = 0; ii < ncoef; ii++){
  //         if(ccoef[ii] >= jtol * vol) continue;
  //         METRIS_THROW_MSG(" - 1 reject validity coef {:15.7e} scaled {:15.7e} \n",
  //                 ccoef[ii], ccoef[ii]/vol);
  //       }
  //     }
  //   }else{
  //     for(int ientt : lball){
  //       if(isvalideltP1<idim,idim>(msh,ientt)) continue;
  //       METRIS_THROW_MSG(" - 2 reject validity\n");
  //     }
  //   }
  // }

  return ierro;
}


#define BOOST_PP_LOCAL_MACRO(n)\
template int smoocavdiff<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        >& msh,\
                  MshCavity& cav, \
                   double& quaCav1, double& quaMaxCav1, \
                   double& targetWeightCav1, \
                   QuaFun iquaf, int ithread);\
template int smoocavdiff<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        >& msh,\
                  MshCavity& cav, \
                   double& quaCav1, double& quaMaxCav1, \
                   double& targetWeightCav1, \
                   QuaFun iquaf, int ithread);\
template int smoocavdiff<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical>& msh,\
                  MshCavity& cav, \
                   double& quaCav1, double& quaMaxCav1, \
                   double& targetWeightCav1, \
                   QuaFun iquaf, int ithread);\
template int smoocavdiff<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical>& msh,\
                  MshCavity& cav, \
                   double& quaCav1, double& quaMaxCav1, \
                   double& targetWeightCav1, \
                   QuaFun iquaf, int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

// =============================================================================================== //
// =============================================================================================== //

template<class MFT, int idim, int ideg>
int smoocavdiff_boundary(Mesh<MFT>& msh, MshCavity& cav, const int cadDim,
                         double& quaCav1, double& quaMaxCav1,
                         double& targetWeightCav1,
                         QuaFun iquaf, int ithread){

  GETVDEPTH(msh.param);

  constexpr int tdim = idim;
  constexpr int gdim = idim;

  METRIS_ASSERT(cadDim == 1 || cadDim == 2);

  if(cadDim != 1){
    CPRINTF1(" # smoocavdiff_boundary only implemented for CAD edges for now\n");
    return 1;
  }

  const intAr1& lcent = cav.lcent(tdim);
  const intAr2& ent2ent = msh.ent2ent(tdim);
  const int tmpEntt = msh.nentt(tdim)-1;
        intAr2& ent2poi = msh.ent2poi(tdim);
        intAr2& ent2tag = msh.ent2tag(tdim);

  const int ipins = cav.ipins;
  targetWeightCav1 = 0.;

  CPRINTF1("-- START smoocavdiff_boundary ipins = {} ncav {}\n",ipins,lcent.get_n());

  msh.tag[ithread]++;
  for (int ienttCav : lcent) ent2tag(ithread,ienttCav) = msh.tag[ithread];

  constexpr int nnmet = (idim*(idim+1))/2;
  constexpr int nhess = nnmet;

  const auto quafun = get_quafun<MFT,gdim,tdim>(iquaf);
  const auto d_quafun = get_d_quafun<MFT,gdim,tdim>(iquaf);

  const int ibpoin = msh.poi2bpo[ipins];
  METRIS_ASSERT(ibpoin >= 0);
  METRIS_ASSERT(msh.bpo2ibi(ibpoin,0) == ipins);
  METRIS_ASSERT(msh.bpo2ibi(ibpoin,1) == 1);

  const int iedins = msh.bpo2ibi(ibpoin,2);
  METRIS_ASSERT(iedins >= 0 && iedins < msh.nedge);
  const int irefins = msh.edg2ref[iedins];
  ego cadEdge = msh.CAD.cad2edg[irefins];

  double range[2];
  int iperi;
  int ierro = 0;
  ierro = EG_getRange(cadEdge,range,&iperi);
  METRIS_ASSERT_MSG(ierro == EGADS_SUCCESS, "EG_getRange failed estat={}", ierro);
  if(range[0] > range[1]){
    double tmp = range[0];
    range[0] = range[1];
    range[1] = tmp;
  }

  int miter1 = MAX(1,msh.param->iflag1);

  const double ftol = 1.0e-2;

  constexpr bool inBoundary = true;
  constexpr int optDimension = 1;
  newton_drivertype_args<optDimension> nargs(msh.param);
  nargs.stpmin = 1.0e-6;
  nargs.wlfc1 = 0.1;
  nargs.wlfc2 = 10.0;
  nargs.ratnew= 0.5;
  nargs.maxit = 50;
  nargs.ftol  = ftol;

  int iflag = 0, ihess;
  bool iinva;
  double tcur[1], coor0[idim], met0[nnmet], fcur;
  double t0 = msh.bpo2rbi(ibpoin,0);
  double d1t[1] = {0.};
  double d2t[1] = {0.};

  for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipins,ii);
  for(int ii = 0; ii < nnmet;ii++) met0[ii]  = msh.met(ipins,ii);
  nargs.xopt[0] = t0;

  std::ofstream dbg;
  if(const char* dbgfile = std::getenv("METRIS_SMOOCAVDIFF_BOUNDARY_TRACE")){
    dbg.open(dbgfile,std::ios::app);
    dbg << "smoocavdiff_boundary start"
        << " ipins " << ipins
        << " cadDim " << cadDim
        << " iedins " << iedins
        << " irefins " << irefins
        << " ncav " << lcent.get_n()
        << " range [" << range[0] << " " << range[1] << "]"
        << " periodic " << iperi
        << " t0 " << t0
        << " coord0 [";
    for(int ii = 0; ii < idim; ii++){
      if(ii > 0) dbg << " ";
      dbg << msh.coord(ipins,ii);
    }
    dbg << "]\n";
  }

  double egParam[2] = {t0,0.};
  double evalResult[18];
  double Xt[3] = {0.,0.,0.};
  double Xtt[3] = {0.,0.,0.};

  for(int niter1 = 0; niter1 < miter1; niter1++){

    tcur[0] = msh.bpo2rbi(ibpoin,0);
    while(true){
      INCVDEPTH(msh.param);

      if(dbg.good()){
        dbg << "  driver_in"
            << " niter " << nargs.niter
            << " iflag " << iflag
            << " ihess " << ihess
            << " tcur " << tcur[0]
            << " fcur " << fcur
            << " d1t " << d1t[0]
            << " d2t " << d2t[0]
            << " fopt " << nargs.fopt
            << " xopt " << nargs.xopt[0]
            << "\n";
      }

      ierro = optim_newton_drivertype<optDimension,inBoundary>(nargs, tcur,
                                                               &fcur, d1t, d2t,
                                                               &iflag, &ihess);

      if(dbg.good()){
        dbg << "  driver_out"
            << " niter " << nargs.niter
            << " iflag " << iflag
            << " ihess " << ihess
            << " ierro " << ierro
            << " tcur " << tcur[0]
            << " fcur " << fcur
            << " d1t " << d1t[0]
            << " d2t " << d2t[0]
            << " fopt " << nargs.fopt
            << " xopt " << nargs.xopt[0]
            << "\n";
      }

      if(ierro > 0){
        CPRINTF1(" # optim_newton_drivertype error {}\n",ierro);
        goto finish;
      }
      if(iflag <= 0) {
        CPRINTF1(" - iflag = 0 termination\n");
        if(dbg.good()){
          dbg << "  driver_termination"
              << " iflag " << iflag
              << " ierro " << ierro
              << " tcur " << tcur[0]
              << " fopt " << nargs.fopt
              << " xopt " << nargs.xopt[0]
              << "\n";
        }
        break;
      }

      const double ttol = 1.0e-12*MAX(1.0,range[1] - range[0]);
      if(tcur[0] < range[0] - ttol || tcur[0] > range[1] + ttol){
        fcur = 1.0e15;
        quaCav1 = fcur;
        quaMaxCav1 = fcur;
        if(dbg.good()){
          dbg << "  trial_out_of_range"
              << " tcur " << tcur[0]
              << " range [" << range[0] << " " << range[1] << "]"
              << "\n";
        }
        continue;
      }
      if(tcur[0] < range[0]){
        tcur[0] = range[0];
      }else if(tcur[0] > range[1]){
        tcur[0] = range[1];
      }
      msh.bpo2rbi(ibpoin,0) = tcur[0];
      egParam[0] = tcur[0];
      int estat = EG_evaluate(cadEdge, egParam, evalResult);
      METRIS_ASSERT_MSG(estat == EGADS_SUCCESS, "EG_evaluate failed estat={}", estat);

      for(int ii = 0; ii < idim; ii++) msh.coord(ipins,ii) = evalResult[ii];
      Xt[0]  = evalResult[3]; Xt[1]  = evalResult[4]; Xt[2]  = idim == 3 ? evalResult[5] : 0.0;
      Xtt[0] = evalResult[6]; Xtt[1] = evalResult[7]; Xtt[2] = idim == 3 ? evalResult[8] : 0.0;

      iinva = false;
      fcur = 0;
      quaMaxCav1 = -1.0e30;
      double targetWeightCurrent = 0.;
      double gradX[idim], hessX[nhess];
      double dqelt[idim], hqelt[nhess];
      for(int ii = 0; ii < idim; ii++) gradX[ii] = 0;
      for(int ii = 0; ii < nhess;ii++) hessX[ii] = 0;

      if constexpr (ideg == 1){
        for (const int ienttCav : lcent){

          int ent2pol[4];
          for(int jj = 0; jj < tdim + 1; jj++){

            const int ienei = ent2ent(ienttCav,jj);

            if (ienei >= 0 && ent2tag(ithread,ienei) >= msh.tag[ithread]) continue;

            if constexpr (tdim == 2){

              int iedgeGlobal = msh.fac2edg(ienttCav,jj);
              if(iedgeGlobal >= 0 && msh.edg2ref[iedgeGlobal] == irefins) continue;

              ent2pol[0] = ipins;
              ent2pol[lnoed2[0][0]] = ent2poi(ienttCav,lnoed2[jj][0]);
              ent2pol[lnoed2[0][1]] = ent2poi(ienttCav,lnoed2[jj][1]);

              if (ent2pol[1] == ipins || ent2pol[2] == ipins) continue;

              ent2poi(tmpEntt,0) = ent2pol[0];
              ent2poi(tmpEntt,1) = ent2pol[1];
              ent2poi(tmpEntt,2) = ent2pol[2];

              double meas;
              if (!isvalideltP1<2,2>(msh, tmpEntt, NULL, &meas)){
                iinva = true;
                break;
              }

              int ivar  = msh.template getverent<ideg>(tmpEntt,idim,ipins);
              double quael;
              if(ihess){
                quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                                tmpEntt,ivar,
                                msh.getBasis(),
                                DifVar::None,dqelt,hqelt,1.);
              }else{
                quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                                tmpEntt,ivar,
                                msh.getBasis(),
                                DifVar::None,dqelt,NULL,1.);
              }
              double regionWeight = 1.;
              if(iquaf == QuaFun::StepDistance
                 && msh.param->step_distance_cavity_target_average){
                regionWeight =
                    step_distance_element_target_weight<MFT,idim,idim>(
                        msh,AsDeg::Pk,tmpEntt);
                targetWeightCurrent += regionWeight;
              }
              fcur += regionWeight*quael;
              if (quael > quaMaxCav1) quaMaxCav1 = quael;

              for(int ii = 0; ii < idim; ii++){
                gradX[ii] += regionWeight*dqelt[ii];
              }
              if(ihess)
                for(int ii = 0; ii < nhess;ii++){
                  hessX[ii] += regionWeight*hqelt[ii];
                }

            }else{
              METRIS_THROW_MSG("Implement CAD edge cavity smoothing in 3D");
            }
          }
          if (iinva) break;
        }
      }else{
        METRIS_THROW_MSG("Implement cavity smoothing for ideg > 1");
      }

      if(!iinva && iquaf == QuaFun::StepDistance
         && msh.param->step_distance_cavity_target_average){
        METRIS_ENFORCE(targetWeightCurrent > 0.);
        fcur /= targetWeightCurrent;
        for(int ii = 0; ii < idim; ii++){
          gradX[ii] /= targetWeightCurrent;
        }
        if(ihess){
          for(int ii = 0; ii < nhess; ii++){
            hessX[ii] /= targetWeightCurrent;
          }
        }
      }

      if(iinva){
        fcur = 1e15;
        quaCav1 = fcur;
        quaMaxCav1 = 1e15;
        CPRINTF1("# invalid config -> continue");
        if(dbg.good()){
          dbg << "  trial_invalid"
              << " tcur " << tcur[0]
              << "\n";
        }
        continue;
      }

      double dfdt = 0.;
      for(int ii = 0; ii < idim; ii++) dfdt += gradX[ii] * Xt[ii];

      double XtHXt = 0.;
      if(ihess){
        auto H = [&](int i, int j) -> double {
          if(i == j) return hessX[i];
          if constexpr (idim == 2){
            return hessX[2];
          }else{
            if((i == 0 && j == 1) || (i == 1 && j == 0)) return hessX[3];
            if((i == 0 && j == 2) || (i == 2 && j == 0)) return hessX[4];
            return hessX[5];
          }
        };
        for(int ii = 0; ii < idim; ii++)
          for(int jj = 0; jj < idim; jj++)
            XtHXt += Xt[ii] * H(ii,jj) * Xt[jj];
      }

      double gradXdotXtt = 0.;
      for(int ii = 0; ii < idim; ii++) gradXdotXtt += gradX[ii] * Xtt[ii];

      d1t[0] = dfdt;
      if(ihess) d2t[0] = XtHXt + gradXdotXtt;

      bool useGrad = false;
      double dtmaxGrad = 0.;
      if(ihess && std::abs(d2t[0]) > 0){
        double dtnew = -d1t[0] / d2t[0];
        if(dtnew * d1t[0] >= 0){
          dtmaxGrad = 0.1*(range[1] - range[0]);
          double distToBnd = d1t[0] > 0 ? tcur[0] - range[0]
                                        : range[1] - tcur[0];
          if(distToBnd > 0) dtmaxGrad = MIN(dtmaxGrad,0.5*distToBnd);
          dtmaxGrad = MAX(dtmaxGrad,1.0e-12);
          d2t[0] = std::abs(d1t[0]) / dtmaxGrad;
          useGrad = true;
        }
      }

      if(dbg.good()){
        dbg << "  objective"
            << " tcur " << tcur[0]
            << " coord [";
        for(int ii = 0; ii < idim; ii++){
          if(ii > 0) dbg << " ";
          dbg << msh.coord(ipins,ii);
        }
        dbg << "] fcur " << fcur
            << " ihess " << ihess
            << " gradX [";
        for(int ii = 0; ii < idim; ii++){
          if(ii > 0) dbg << " ";
          dbg << gradX[ii];
        }
        dbg << "] Xt [" << Xt[0] << " " << Xt[1] << " " << Xt[2] << "]"
            << " Xtt [" << Xtt[0] << " " << Xtt[1] << " " << Xtt[2] << "]"
            << " d1t " << d1t[0]
            << " d2t " << d2t[0];
        if(ihess && std::abs(d2t[0]) > 0){
          dbg << " newton_dt " << -d1t[0]/d2t[0];
        }
        if(useGrad) dbg << " use_grad_descent 1"
                         << " grad_dtmax " << dtmaxGrad;
        dbg << "\n";
      }

      if(DOPRINTS1()){
        CPRINTF1(" - Newton iter {} fcur = {} tcur = {}",nargs.niter,fcur,tcur[0]);
        PRINTF("\n");
        CPRINTF2(" - dquadt = {}\n",d1t[0]);
      }

      quaCav1 = fcur;
    }

    ierro = 0;

    finish:
    CPRINTF1(" -- END smoocavdiff_boundary fopt = {} xopt = {}\n",
             nargs.fopt,dblAr1(optDimension,nargs.xopt));
    if(dbg.good()){
      dbg << "  finish"
          << " fopt " << nargs.fopt
          << " xopt " << nargs.xopt[0]
          << " ierro " << ierro
          << "\n";
    }

    msh.bpo2rbi(ibpoin,0) = nargs.xopt[0];
    egParam[0] = nargs.xopt[0];
    int estat = EG_evaluate(cadEdge, egParam, evalResult);
    METRIS_ASSERT_MSG(estat == EGADS_SUCCESS, "EG_evaluate failed estat={}", estat);
    for(int ii = 0; ii < idim; ii++) msh.coord(ipins,ii) = evalResult[ii];

    if(DOPRINTS2()) writeMesh("debug_smooth0.meshb",msh);

    ierro = msh.interpMetBack(ipins);
    if(ierro > 0){
      CPRINTF1(" # smoocavdiff_boundary interpMetBack failure ierro = {} \n",ierro);
      goto cleanup;
    }

    quaCav1 = 0;
    quaMaxCav1 = -1.0e30;
    double targetWeightFinal = 0.;
    iinva = false;
    for(const int ienttCav : lcent){
      int ent2pol[4];
      for(int jj = 0; jj < tdim + 1; jj++){
        const int ienei = ent2ent(ienttCav,jj);
        if(ienei >= 0 && ent2tag(ithread,ienei) >= msh.tag[ithread]) continue;

        if constexpr (tdim == 2){
          int iedgeGlobal = msh.fac2edg(ienttCav,jj);
          if(iedgeGlobal >= 0 && msh.edg2ref[iedgeGlobal] == irefins) continue;

          ent2pol[0] = ipins;
          ent2pol[lnoed2[0][0]] = ent2poi(ienttCav,lnoed2[jj][0]);
          ent2pol[lnoed2[0][1]] = ent2poi(ienttCav,lnoed2[jj][1]);

          if(ent2pol[1] == ipins || ent2pol[2] == ipins) continue;

          ent2poi(tmpEntt,0) = ent2pol[0];
          ent2poi(tmpEntt,1) = ent2pol[1];
          ent2poi(tmpEntt,2) = ent2pol[2];

          double meas;
          if(!isvalideltP1<2,2>(msh, tmpEntt, NULL, &meas)){
            iinva = true;
            break;
          }
        }else{
          METRIS_THROW_MSG("Implement CAD edge cavity smoothing in 3D");
        }

        double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,tmpEntt,1.);
        double regionWeight = 1.;
        if(iquaf == QuaFun::StepDistance
           && msh.param->step_distance_cavity_target_average){
          regionWeight =
              step_distance_element_target_weight<MFT,idim,idim>(
                  msh,AsDeg::Pk,tmpEntt);
          targetWeightFinal += regionWeight;
        }
        quaCav1 += regionWeight*quael;
        quaMaxCav1 = MAX(quaMaxCav1, quael);
      }
      if(iinva) break;
    }
    if(iinva){
      CPRINTF1(" # smoocavdiff_boundary final xopt quality recompute invalid\n");
      goto cleanup;
    }
    if(iquaf == QuaFun::StepDistance
       && msh.param->step_distance_cavity_target_average){
      METRIS_ENFORCE(targetWeightFinal > 0.);
      targetWeightCav1 = targetWeightFinal;
    }
    if(dbg.good()){
      dbg << "  final_recompute"
          << " quaCav1 " << quaCav1
          << " quaMaxCav1 " << quaMaxCav1
          << " coord [";
      for(int ii = 0; ii < idim; ii++){
        if(ii > 0) dbg << " ";
        dbg << msh.coord(ipins,ii);
      }
      dbg << "]\n";
    }
  }

  if(dbg.good()) dbg << "smoocavdiff_boundary return 0\n";
  return 0;

  cleanup:
  for(int ii = 0; ii < idim; ii++) msh.coord(ipins,ii) = coor0[ii];
  for(int ii = 0; ii < nnmet; ii++) msh.met(ipins,ii) = met0[ii];
  msh.bpo2rbi(ibpoin,0) = t0;

  if(dbg.good()){
    dbg << "smoocavdiff_boundary cleanup"
        << " ierro " << ierro
        << " restore_t " << t0
        << "\n";
  }

  return ierro;
}


#define BOOST_PP_LOCAL_MACRO(n)\
template int smoocavdiff_boundary<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        >& msh,\
                  MshCavity& cav, const int cadDim, \
                   double& quaCav1, double& quaMaxCav1, \
                   double& targetWeightCav1, \
                   QuaFun iquaf, int ithread);\
template int smoocavdiff_boundary<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        >& msh,\
                   MshCavity& cav, const int cadDim, \
                   double& quaCav1, double& quaMaxCav1, \
                   double& targetWeightCav1, \
                   QuaFun iquaf, int ithread);\
template int smoocavdiff_boundary<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical>& msh,\
                   MshCavity& cav, const int cadDim, \
                   double& quaCav1, double& quaMaxCav1, \
                   double& targetWeightCav1, \
                   QuaFun iquaf, int ithread);\
template int smoocavdiff_boundary<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical>& msh,\
                   MshCavity& cav, const int cadDim, \
                   double& quaCav1, double& quaMaxCav1, \
                   double& targetWeightCav1, \
                   QuaFun iquaf, int ithread);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

// =============================================================================================== //
// =============================================================================================== //

template<class MFT, int idim, int ideg>
double smooballdiff_fun([[maybe_unused]] unsigned int nvar,
                        const double *xcur,
                        double *grad, void *f_data){
  constexpr int gdim = idim;
  constexpr int tdim = idim;

  smooballdiff_fun_data<MFT> *mydata = (smooballdiff_fun_data<MFT>*)(f_data);

  int ipoin = mydata->ipoin;
  Mesh<MFT> &msh = *(mydata->msh);
  const intAr1 &lball = *(mydata->lball);
  QuaFun iquaf = mydata->iquaf;

  double coor0[gdim];
  for(int ii = 0; ii < idim; ii++){
    coor0[ii] = msh.coord(ipoin,ii);
    msh.coord(ipoin,ii) = xcur[ii];
  }


  const auto quafun     = get_quafun<MFT,gdim,tdim>(iquaf);
  const auto d_quafun   = get_d_quafun<MFT,gdim,tdim>(iquaf);
  //const intAr2 &ent2poi = msh.ent2poi(idim);
  bool iinva = false;
  if constexpr (ideg == 1){
    for(int ientt : lball){
      iinva = !isvalideltP1<idim,idim>(msh,ientt);
      if(iinva) break;
    }
  }else{
    constexpr int jdeg = tdim*(ideg - 1);
    constexpr int ncoef = tdim == 2 ? getnnod2(jdeg)
                                    : getnnod3(jdeg);
    double ccoef[ncoef];
    for(int ientt : lball){
      getsclccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef,&iinva);
      if(iinva) break;
    }
  }

  if(iinva){
    // Set the gradient to move away from xcur.
    if(grad != NULL){
      for(int ii = 0; ii < idim; ii++){
        grad[ii] = 1.0e8 * (xcur[ii] - coor0[ii]);
      }
    }
    return 1.0e10;
  }

  double fcur = 0;
  for(int ii = 0; ii < idim && grad != NULL; ii++) grad[ii] = 0;
  for(int ientt : lball){

    bool iflat = !isvalideltP1<idim,idim>(msh,ientt);
    if(iflat) METRIS_THROW_MSG( "Flat after check??");

    int ivar = -1;
    if(grad != NULL) ivar = msh.template getverent<ideg>(ientt,idim,ipoin);
    double dqelt[idim];
    double quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                            ientt,ivar,
                            msh.getBasis(),
                            DifVar::None,dqelt,NULL,1);
    fcur += quael;
    for(int ii = 0; ii < idim && grad != NULL; ii++) grad[ii] += dqelt[ii];

    //if(!mydata->iqset){
    //  mydata->qnrm0 += quael;
    //  mydata->qmax0  = MAX(quael,mydata->qmax0);
    //}
  }

  //mydata->iqset = true;

  if(mydata->fopt > fcur){
    mydata->fopt = fcur;
    for(int ii = 0; ii < gdim; ii++) mydata->xopt[ii] = msh.coord(ipoin,ii);
  }

  return fcur;
}
#define BOOST_PP_LOCAL_MACRO(n)\
template double smooballdiff_fun<MetricFieldAnalytical,2,n>(unsigned int nvar, \
                          const double *x, double *grad, void *f_data);\
template double smooballdiff_fun<MetricFieldAnalytical,3,n>(unsigned int nvar, \
                          const double *x, double *grad, void *f_data);\
template double smooballdiff_fun<MetricFieldFE        ,2,n>(unsigned int nvar, \
                          const double *x, double *grad, void *f_data);\
template double smooballdiff_fun<MetricFieldFE        ,3,n>(unsigned int nvar, \
                          const double *x, double *grad, void *f_data);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


// inorm <= infi norm , p > 0 L^p norm (over ball)
template<class MFT, int idim, int ideg>
int smooballdiff_luksan(Mesh<MFT>& msh, int ipoin,
                        const intAr1 &lball,
                        double*__restrict__ qnrm0, double*__restrict__ qmax0,
                        double*__restrict__ qnrm1, double*__restrict__ qmax1,
                        dblAr1 &work,
                        QuaFun iquaf){

  GETVDEPTH(msh.param);
  constexpr int gdim = idim;
  constexpr int tdim = idim;
  constexpr int nnmet = (idim*(idim+1))/2;

  const int miter1 = MAX(msh.param->iflag1,1);
  const nlopt_algorithm algo = NLOPT_LD_TNEWTON_PRECOND_RESTART;


  const auto quafun = get_quafun<MFT,gdim,tdim>(iquaf);
  const intAr2& ent2poi = msh.ent2poi(idim);

  int ierro = 0;

  *qnrm0 = 0;
  *qmax0 = -1.0e30;
  //auto quafun
  //  = QuaFunList<MFT,gdim,tdim,ideg,AsDeg::Pk,AsDeg::Pk>{}.quafun(iquaf);
  //auto d_quafun
  //  = QuaFunList<MFT,gdim,tdim,ideg,AsDeg::Pk,AsDeg::Pk>{}.d_quafun(iquaf);

  double xopt[gdim];
  smooballdiff_fun_data<MFT> mydata(msh,lball,ipoin,iquaf,xopt);
  int nwork = luksan_pnet_worksize(gdim);
  work.allocate(nwork);
  work.set_n(nwork);
  double fstop = -1.0;
  double ftol_rel = 1.0e-9;
  double ftol_abs = -1e30;
  double lb[gdim], ub[gdim];

  double coor0[gdim], met0[nnmet];
  int poi2bak0;
  for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipoin,ii);
  for(int ii = 0; ii < nnmet;ii++) met0[ii]  = msh.met(ipoin,ii);
  poi2bak0 = msh.poi2bak[ipoin];

  for(int ii = 0; ii < gdim; ii++){
    lb[ii] = -HUGE_VAL;
    ub[ii] =  HUGE_VAL;
  }
  double xcur[gdim];
  for(int ii = 0; ii < idim; ii++) xcur[ii] = msh.coord(ipoin,ii);
  double fopt;

  *qnrm0 = 0;
  *qmax0 = -1.0e30;
  for(int ientt : lball){
    double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,ientt,1);

    *qnrm0 += quael;
    *qmax0  = MAX(quael,*qmax0);
  }

  double fpre = *qnrm0;

  for(int niter1 = 0; niter1 < miter1; niter1++){

    ierro = luksan_pnetS<gdim>(smooballdiff_fun<MFT,gdim,ideg>, &mydata,
                               lb, ub, /* bounds */
                               xcur, /* in: initial guess, out: minimizer */
                               &fopt,
                               //int mf, /* subspace dimension (0 for default) */
                               algo,
                               work,
                               fstop , ftol_rel, ftol_abs);

    CPRINTF1(" - end luksan_pnetS got ierro = {} \n",ierro);
    if(ierro == NLOPT_STOPVAL_REACHED
    || ierro == NLOPT_FTOL_REACHED
    || ierro == NLOPT_XTOL_REACHED) ierro = NLOPT_SUCCESS;

    // If algo failed but we caught a better iterate
    if(mydata.fopt < fpre && mydata.fopt < 1.0e10) ierro = NLOPT_SUCCESS;

    fpre = mydata.fopt;

    if(ierro != NLOPT_SUCCESS) goto cleanup;
    ierro = 0;

    for(int ii = 0; ii < gdim; ii++) msh.coord(ipoin,ii) = mydata.xopt[ii];

    if(msh.interpMetBack(ipoin, idim, lball[0], -1, NULL)) goto cleanup;

    *qnrm1 = 0;
    *qmax1 = -1.0e30;
    for(int ientt : lball){
      double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,ientt,1);

      *qnrm1 += quael;
      *qmax1  = MAX(quael,*qmax1);
    }

    if(*qnrm1 > *qnrm0){
      ierro = 2;
      CPRINTF1(" # Local smoo reject: quality norm increase "
                 "{} -> {} \n", *qnrm0, *qnrm1);
      goto cleanup;
    }

    CPRINTF1(" - local smoothing quality {} -> {} \n", *qnrm0, *qnrm1);


    for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipoin,ii);
    for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh.met(ipoin,ii);
    poi2bak0 = msh.poi2bak[ipoin];
  }

  return 0;


  cleanup:
  for(int ii = 0; ii < idim; ii++)  msh.coord(ipoin,ii) = coor0[ii];
  for(int ii = 0; ii < nnmet; ii++) msh.met(ipoin,ii) = met0[ii];
  msh.poi2bak[ipoin] = poi2bak0;

  *qnrm1 = *qnrm0;
  *qmax1 = *qmax0;

  if(msh.param->dbgfull){
    if constexpr (ideg >= 2){
      constexpr int jdeg = tdim*(ideg - 1);
      constexpr int ncoef = tdim == 2 ? getnnod2(jdeg)
                                      : getnnod3(jdeg);
      const double jtol = msh.param->jtol;
      double ccoef[ncoef];
      for(int ientt : lball){
        double vol = getmeasentP1<idim>(ent2poi[ientt], msh.coord);
        getccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef);
        for(int ii = 0; ii < ncoef; ii++){
          if(ccoef[ii] >= jtol * vol) continue;
          METRIS_THROW_MSG(" - 1 reject validity coef {:15.7e} scaled {:15.7e} \n",
                  ccoef[ii], ccoef[ii]/vol);
        }
      }
    }else{
      for(int ientt : lball){
        if(isvalideltP1<idim,idim>(msh,ientt)) continue;
        METRIS_THROW_MSG(" - 2 reject validity\n");
      }
    }
  }

  return ierro;
}




#define BOOST_PP_LOCAL_MACRO(n)\
template int smooballdiff_luksan<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        >& msh,\
 int ipoin, const intAr1 &lball,\
                   double*__restrict__ qnrm0, double*__restrict__ qmax0, \
                   double*__restrict__ qnrm1, double*__restrict__ qmax1,\
                   dblAr1 &work,\
                   QuaFun iquaf);\
template int smooballdiff_luksan<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        >& msh,\
 int ipoin, const intAr1 &lball,\
                   double*__restrict__ qnrm0, double*__restrict__ qmax0, \
                   double*__restrict__ qnrm1, double*__restrict__ qmax1,\
                   dblAr1 &work,\
                   QuaFun iquaf);\
template int smooballdiff_luksan<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical>& msh,\
 int ipoin, const intAr1 &lball,\
                   double*__restrict__ qnrm0, double*__restrict__ qmax0, \
                   double*__restrict__ qnrm1, double*__restrict__ qmax1,\
                   dblAr1 &work,\
                   QuaFun iquaf);\
template int smooballdiff_luksan<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical>& msh,\
 int ipoin, const intAr1 &lball,\
                   double*__restrict__ qnrm0, double*__restrict__ qmax0, \
                   double*__restrict__ qnrm1, double*__restrict__ qmax1,\
                   dblAr1 &work,\
                   QuaFun iquaf);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()





} // end namespace
