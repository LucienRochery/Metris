//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_newton_drivertype

#include <boost/test/included/unit_test.hpp>
#include "common_setup.hxx"

#include <boost/timer/progress_display.hpp>

#include "utils/aux_misc.hxx"
#include "quality/low_metqua.hxx"
#include "..libs/SANS/Surreal/SurrealS.h"
#include "..libs/SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"
//#include <src/msh_metric.hxx>
#include "Localization/low_localization.hxx"

#include "Optimization/opt_generic.hxx"

#include <boost/hana.hpp>
namespace hana = boost::hana;
using namespace hana::literals;

namespace utf = boost::unit_test;

using namespace Metris;


//double func(double *xcur, int ihess, double *gcur, double *hess){
//	double fcur = 0;
//	for(int i =0 ; i< 3; i++){
//		fcur += xcur[i]*xcur[i]/2;
//		gcur[i] = xcur[i];
//	}
//	if(ihess > 0){
//		hess[0] = 1;
//		hess[1] = 0;
//		hess[2] = 1;
//		hess[3] = 0;
//		hess[4] = 0;
//		hess[5] = 1;
//	}
//	return fcur;
//}


double func(double *xcur, int ihess, double *gcur, double *hess){
	double fcur = 0;
	for(int i = 0 ; i< 3; i++){
		fcur += xcur[i]*xcur[i]*xcur[i]*xcur[i]/24;
		gcur[i] = xcur[i]*xcur[i]*xcur[i]/6;
	}
	if(ihess > 0){
		hess[0] = xcur[0]*xcur[0]/2;
		hess[1] = 0;
		hess[2] = xcur[1]*xcur[1]/2;
		hess[3] = 0;
		hess[4] = 0;
		hess[5] = xcur[2]*xcur[2]/2;
	}
	return fcur;
}

BOOST_AUTO_TEST_CASE(test_newton_drivertype)
{//METRIS_MAX_DEG
	constexpr int nvar = 3;
	double fcur, xcur[3], gcur[3], hess[6];
	//double xtol = 1.0e-6,stpmin = 1.0e-12, wlfc1 = 1.0e-1, wlfc2 = -1;
	//int niter,maxit=50,iprt=0,iflag,ihess;
	//int nrwrk = 34,niwrk=3;
	//double rwork[nrwrk];
	//int    iwork[niwrk];
	//double xopt[3], fopt;
	//int ierro;

  MetrisParameters param;
  param.iverb = 0;

  newton_drivertype_args<nvar> args(&param);
  args.stpmin = 1.0e-12;
  args.ratnew = 0.5; // LS step decrease factor
  args.maxit = 500;
  args.wlfc1 = 1e-1;
  args.wlfc2 = 0.9;

	int ierro, ihess, iflag = 0;
  double fopt;

	for(int ii = 0; ii < 3; ii++) xcur[ii] = 1;

  fmt::print("-- Test basic Newton\n");
  bool istopped = false;
	do{
		//optim_newton_drivertype(nvar ,
		//                                  xcur ,&fcur  ,gcur   ,hess ,
		//                                  xtol ,stpmin,0, 1,
		//                                  wlfc1,wlfc2 ,0.1 ,
		//                                  &niter,maxit ,iprt   ,
		//                                  &iflag,&ihess ,
		//                                  nrwrk,rwork ,
		//                                  niwrk,iwork ,
		//                                  xopt ,&fopt ,&ierro);
    ierro = optim_newton_drivertype<nvar>(args,
                                       xcur, &fcur, gcur, hess, &iflag, &ihess);

    BOOST_CHECK(ierro <= 0);
    if(ierro > 0) break;

		if(iflag <= 0){
			//printf("Stop\n");
      istopped = true;
			break;
		}

		fcur = func(xcur,ihess,gcur,hess);
		//printf("-- Iteration %d flag %d value = %e opt = %f \n",args.niter,iflag,fcur,args.fopt);
		//printf("  - debug xcur %f %f %f grad %.2e %.2e %.2e \n",xcur[0],xcur[1],xcur[2]
		//	,gcur[0],gcur[1],gcur[2]);
	}while(iflag > 0);

  BOOST_CHECK(istopped);
  double err = getnrml2<nvar>(args.xopt);
  BOOST_CHECK_SMALL(err,1.0e-12);
  BOOST_CHECK_SMALL(args.fopt,1.0e-12);


  fmt::print("-- Test TNCG Newton\n");
  double buff[200];
  dblAr1 buf(200,buff);
  truncated_newton_work work(nvar,buf);

  args.niter = 0;
  //param.iverb = 3;
  //param.ivdepth = 5;

  iflag = 0;
  for(int ii = 0; ii < 3; ii++) xcur[ii] = 1;
  istopped = false;
  do{
    optim_newton_drivertype_TNCG<nvar>(args ,work, xcur, &fcur, gcur, hess, &iflag, &ihess);

    BOOST_CHECK(ierro <= 0);
    if(ierro > 0) break;

    if(iflag <= 0){
      istopped = true;
      break;
    }

    fcur = func(xcur,ihess,gcur,hess);
    //printf("-- Iteration %d flag %d value = %e opt = %f \n",args.niter,iflag,fcur,fopt);
    //printf("  - debug xcur %f %f %f grad %.2e %.2e %.2e \n",xcur[0],xcur[1],xcur[2]
    //  ,gcur[0],gcur[1],gcur[2]);
  }while(iflag > 0);

  BOOST_CHECK(istopped);
  double err_TNCG = getnrml2<nvar>(args.xopt);
  BOOST_CHECK_SMALL(err_TNCG,1.0e-5);
  BOOST_CHECK_SMALL(args.fopt,1.0e-5);

  fmt::print("-- END Newton tests final x err: basic {:.2e}, TNCG = {:.2e}\n", err, err_TNCG);

}









