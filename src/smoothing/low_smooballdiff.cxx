//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
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
#include "../aux_timer.hxx"
#include "../io_libmeshb.hxx"
#include "../mprintf.hxx"

#include "../opt_generic.hxx"
#include "../quality/low_metqua_d.hxx"
#include "../low_ccoef.hxx"
#include "../low_geo.hxx"


namespace Metris{

// inorm <= infi norm , p > 0 L^p norm (over ball)
template<class MFT, int idim, int ideg>
int smooballdiff(Mesh<MFT>& msh, int ipoin, 
                   const intAr1 &lball, 
                   double*__restrict__ qnrm0, double*__restrict__ qmax0, 
                   double*__restrict__ qnrm1, double*__restrict__ qmax1,
                   QuaFun iquaf){

  GETVDEPTH(msh);

  constexpr int tdim = idim;
  constexpr int gdim = idim;

  constexpr int nnmet = (idim*(idim+1))/2;
  constexpr int nhess = nnmet;

  const intAr2& ent2poi = msh.ent2poi(idim); 
  const int qpower  = msh.param->opt_power;
  const int qpnorm  = msh.param->opt_pnorm;

  //auto quafun 
  //  = QuaFunList<MFT,gdim,tdim,ideg,AsDeg::Pk,AsDeg::Pk>{}.quafun(iquaf);
  //auto d_quafun 
  //  = QuaFunList<MFT,gdim,tdim,ideg,AsDeg::Pk,AsDeg::Pk>{}.d_quafun(iquaf);

  const auto quafun = get_quafun<MFT,gdim,tdim>(iquaf);
  const auto d_quafun = get_d_quafun<MFT,gdim,tdim>(iquaf);

  int nball = lball.get_n();

  *qnrm0 = 0;
  *qmax0 = -1.0e30;

  // Optimization doesn't reinterpolate metric 
  int miter1 = MAX(1,msh.param->iflag1);

  // Relative decrease tolerance 
  const double ftol = 1.0e-2;

  double xtol = -1,stpmin = 1.0e-6, 
         wlfc1 = 0.1 , wlfc2 = 10.0,ratnew = 0.5;
  int niter = 0, maxit = 50, iflag = 0, ihess, ierro = 0;
  const int nrwrk = 34, niwrk = 3;
  int iprt = 0;
  double rwrkN[nrwrk];
  int    iwrkN[niwrk];
  double xopt[idim], xcur[idim], coor0[idim], met0[nnmet], fopt = -1, fcur;
  double fpre; 
  bool fpreset = false;
  bool iinva;
  double d1qua[idim], d2qua[nhess];

  for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipoin,ii);
  for(int ii = 0; ii < nnmet;ii++) met0[ii]  = msh.met(ipoin,ii);

  for(int niter1 = 0; niter1 < miter1; niter1++){

    for(int ii = 0; ii < idim; ii++) xcur[ii]  = msh.coord(ipoin,ii);
    while(true){
      INCVDEPTH(msh);

      optim_newton_drivertype(idim  ,
                              xcur  , &fcur , d1qua  , d2qua ,
                              xtol  , stpmin, 1      ,
                              wlfc1 , wlfc2 , ratnew ,
                              &niter, maxit , iprt,
                              &iflag, &ihess,
                              nrwrk , rwrkN ,
                              niwrk , iwrkN ,
                              xopt  , &fopt , &ierro);
      if(!fpreset){
        fpreset = true;
        fpre = fcur; 
      }else{
        if(abs(fpre - fcur) < ftol * abs(fpre)){
          CPRINTF1(" - Relative decrease %15.7e < %15.7e end\n",
                                abs(fpre - fcur) / abs(fpre), ftol);
          break;
        }
        fpre = fcur;
      }
      if(ierro > 0){
        CPRINTF1(" # optim_newton_drivertype error %d\n",ierro);
        goto cleanup;
      }
      if(iflag <= 0) {
        CPRINTF1(" - iflag = 0 termination\n");
        break;
      }

      for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = xcur[ii];

      iinva = false;
      if constexpr (ideg == 1){
        for(int ientt : lball){
          getmeasentP1<idim,idim>(msh, ent2poi[ientt], NULL, &iinva);
          if(iinva) break;
        }
      }else{
        constexpr int jdeg = tdim*(ideg - 1);
        constexpr int ncoef = tdim == 2 ? facnpps[jdeg]
                                        : tetnpps[jdeg];
        double ccoef[ncoef];
        for(int ientt : lball){
          getsclccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef,&iinva); 
          if(iinva) break;
        }
      }

      if(iinva){
        fcur = 1.0e10;
        // radical solution for now 
        goto cleanup;
      }


      fcur = 0;
      double dqelt[idim], hqelt[nhess];
      for(int ii = 0; ii < idim; ii++) d1qua[ii] = 0;
      for(int ii = 0; ii < nhess;ii++) d2qua[ii] = 0;
      for(int iball = 0; iball < nball && !iinva; iball++){
        int ient2 = lball[iball];

        bool iflat;
        getmeasentP1<idim,idim>(msh, ent2poi[ient2], NULL, &iflat);
        if(iflat){
          fcur = 1.0e10;
          break;
        }

        int ivar  = msh.template getverent<ideg>(ient2,idim,ipoin);
        double quael;
        if(ihess){
          quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                           ient2,qpower,ivar,
                           msh.getBasis(),
                           DifVar::None,dqelt,hqelt,
                           qpnorm, 1);
        }else{
          quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                           ient2,qpower,ivar,
                           msh.getBasis(),
                           DifVar::None,dqelt,NULL,
                           qpnorm, 1);
        }
        fcur += quael;
        for(int ii = 0; ii < idim; ii++) d1qua[ii] += dqelt[ii];
        if(ihess)
          for(int ii = 0; ii < nhess;ii++) d2qua[ii] += hqelt[ii];

        if(niter == 1 && niter1 == 0){
          *qnrm0 += quael; 
          *qmax0  = MAX(quael,*qmax0);
        }
      }
    } // end while true

    ierro = 0;

    for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = xopt[ii];

    if(DOPRINTS2()) writeMesh("debug_smooth0.meshb",msh);

    ierro = msh.interpMetBack(ipoin, idim, lball[0], -1, NULL);
    if(ierro > 0){
      CPRINTF1(" # interpMetBack failure ierro = %d \n",ierro);
      //printf(" ipoin = %d \n",ipoin);
      //writeMesh("debug_smooth_interpMetBack.meshb",msh);
      //exit(1);
      goto cleanup;
    }


    for(int iball = 0; iball < nball; iball++){
      int ient2 = lball[iball];
      bool iflat;
      getmeasentP1<idim,idim>(msh, ent2poi[ient2], NULL, &iflat);
      if(iflat){
        CPRINTF1(" # Flat iball %d element %d \n",iball,ient2);
        goto cleanup;
      }
    }

    *qnrm1 = 0;
    *qmax1 = -1.0e30;
    for(int iball = 0; iball < nball; iball++){
      int ient2 = lball[iball];
      double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,
                            ient2,qpower,qpnorm,1);

      *qnrm1 += quael; 
      *qmax1  = MAX(quael,*qmax1);
    }
    CPRINTF1(" - Newton update initial quality avg %15.7e " 
                          "max %15.7e \n",*qnrm0,*qmax0);
    CPRINTF1(" -                 final quality avg %15.7e " 
                          "max %15.7e \n",*qnrm1,*qmax1);
  }


  if(*qnrm1 > *qnrm0){
    ierro = 2;
    CPRINTF1(" # Local smoo reject: quality norm increase "
               "%f -> %f \n", *qnrm0, *qnrm1);
    goto cleanup;
  }

  if(msh.param->dbgfull){
    for(int ientt : lball){
      if constexpr (ideg == 1){
        getmeasentP1<idim,idim>(msh, ent2poi[ientt], NULL, &iinva);
        METRIS_ENFORCE(!iinva);
      }else{
        constexpr int jdeg = tdim*(ideg - 1);
        constexpr int ncoef = tdim == 2 ? facnpps[jdeg]
                                        : tetnpps[jdeg];
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
      constexpr int ncoef = tdim == 2 ? facnpps[jdeg]
                                      : tetnpps[jdeg];
      const double jtol = msh.param->jtol;
      double ccoef[ncoef];
      for(int ientt : lball){
        double vol = getmeasentP1<idim>(ent2poi[ientt], msh.coord);
        getccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef); 
        for(int ii = 0; ii < ncoef; ii++){
          if(ccoef[ii] >= jtol * vol) continue;
          printf(" - 1 reject validity coef %15.7e scaled %15.7e \n",
                  ccoef[ii], ccoef[ii]/vol);
          METRIS_THROW(GeomExcept());
        }
      }
    }else{
      for(int ientt : lball){
        bool iflat;
        double meas0 = getmeasentP1<idim,idim>(msh, ent2poi[ientt], NULL, &iflat);
        if(!iflat && meas0 > 0) continue;
          printf(" - 2 reject validity\n");
        METRIS_THROW(GeomExcept());
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
  const intAr2 &ent2poi = msh.ent2poi(idim);
  bool iinva = false;
  if constexpr (ideg == 1){
    for(int ientt : lball){
      getmeasentP1<idim,idim>(msh, ent2poi[ientt], NULL, &iinva);
      if(iinva) break;
    }
  }else{
    constexpr int jdeg = tdim*(ideg - 1);
    constexpr int ncoef = tdim == 2 ? facnpps[jdeg]
                                    : tetnpps[jdeg];
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

  const int qpower  = msh.param->opt_power;
  const int qpnorm  = msh.param->opt_pnorm;
  double fcur = 0;
  for(int ii = 0; ii < idim && grad != NULL; ii++) grad[ii] = 0;
  for(int ientt : lball){

    bool iflat;
    getmeasentP1<idim,idim>(msh, ent2poi[ientt], NULL, &iflat);
    if(iflat) METRIS_THROW_MSG(GeomExcept(), "Flat after check??");

    int ivar = -1;
    if(grad != NULL) ivar = msh.template getverent<ideg>(ientt,idim,ipoin);
    double dqelt[idim];
    double quael = d_quafun(msh,AsDeg::Pk,AsDeg::Pk,
                            ientt,qpower,ivar,
                            msh.getBasis(),
                            DifVar::None,dqelt,NULL,
                            qpnorm, 1);
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

  GETVDEPTH(msh);
  constexpr int gdim = idim;
  constexpr int tdim = idim;
  constexpr int nnmet = (idim*(idim+1))/2;

  const int miter1 = MAX(msh.param->iflag1,1);
  const nlopt_algorithm algo = NLOPT_LD_TNEWTON_PRECOND_RESTART;


  const auto quafun = get_quafun<MFT,gdim,tdim>(iquaf);
  const int qpower  = msh.param->opt_power;
  const int qpnorm  = msh.param->opt_pnorm;
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
  int poi2bak0[3];
  for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipoin,ii);
  for(int ii = 0; ii < nnmet;ii++) met0[ii]  = msh.met(ipoin,ii);
  for(int ii = 0; ii < msh.get_tdim();ii++) poi2bak0[ii] = msh.poi2bak(ipoin,ii);

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
    double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,ientt,qpower,qpnorm,1);

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

    CPRINTF1(" - end luksan_pnetS got ierro = %d \n",ierro);
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
      double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,ientt,qpower,qpnorm,1);

      *qnrm1 += quael; 
      *qmax1  = MAX(quael,*qmax1);
    }

    if(*qnrm1 > *qnrm0){
      ierro = 2;
      CPRINTF1(" # Local smoo reject: quality norm increase "
                 "%f -> %f \n", *qnrm0, *qnrm1);
      goto cleanup;
    }

    CPRINTF1(" - local smoothing quality %f -> %f \n", *qnrm0, *qnrm1);


    for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipoin,ii);
    for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh.met(ipoin,ii);
    for(int ii = 0; ii < msh.get_tdim(); ii++) poi2bak0[ii] =  msh.poi2bak(ipoin,ii);
  }

  return 0;


  cleanup:
  for(int ii = 0; ii < idim; ii++)  msh.coord(ipoin,ii) = coor0[ii];
  for(int ii = 0; ii < nnmet; ii++) msh.met(ipoin,ii) = met0[ii];
  for(int ii = 0; ii < msh.get_tdim(); ii++) msh.poi2bak(ipoin,ii) = poi2bak0[ii];

  *qnrm1 = *qnrm0;
  *qmax1 = *qmax0;

  if(msh.param->dbgfull){
    if constexpr (ideg >= 2){
      constexpr int jdeg = tdim*(ideg - 1);
      constexpr int ncoef = tdim == 2 ? facnpps[jdeg]
                                      : tetnpps[jdeg];
      const double jtol = msh.param->jtol;
      double ccoef[ncoef];
      for(int ientt : lball){
        double vol = getmeasentP1<idim>(ent2poi[ientt], msh.coord);
        getccoef<gdim,tdim,ideg>(msh,ientt,NULL,ccoef); 
        for(int ii = 0; ii < ncoef; ii++){
          if(ccoef[ii] >= jtol * vol) continue;
          printf(" - 1 reject validity coef %15.7e scaled %15.7e \n",
                  ccoef[ii], ccoef[ii]/vol);
          METRIS_THROW(GeomExcept());
        }
      }
    }else{
      for(int ientt : lball){
        bool iflat;
        double meas0 = getmeasentP1<idim,idim>(msh, ent2poi[ientt], NULL, &iflat);
        if(!iflat && meas0 > 0) continue;
          printf(" - 2 reject validity\n");
        METRIS_THROW(GeomExcept());
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
