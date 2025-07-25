//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 

#include <boost/test/included/unit_test.hpp> 
#include "common_setup.hxx"

#include <boost/timer/progress_display.hpp>

#include "../src/utils/aux_misc.hxx"
#include "../src/quality/low_metqua.hxx"
#include "../SANS/Surreal/SurrealS.h"
#include "../SANS/LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"
//#include <src/msh_metric.hxx>
#include "../src/Localization/low_localization.hxx"

#include "../src/Optimization/opt_generic.hxx"

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

BOOST_AUTO_TEST_CASE(test_eval3) 
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

  newton_drivertype_args<nvar> args(&param);
  args.stpmin = 1.0e-12;
  args.ratnew = 0.5; // LS step decrease factor 
  args.iprt = 0;
  args.maxit = 500;
  args.wlfc1 = 1e-1;
  args.wlfc2 = 0.9;

	int ierro, ihess, iflag = 0;
  double fopt;

	for(int ii = 0; ii < 3; ii++) xcur[ii] = 1;

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
		if(ierro > 0){
			printf("ERROR %d\n",ierro);
			exit(1);
		}
		if(iflag <= 0){
			printf("Stop\n");
			break;
		}

		fcur = func(xcur,ihess,gcur,hess);
		printf("-- Iteration %d flag %d value = %e opt = %f \n",args.niter,iflag,fcur,args.fopt);
		printf("  - debug xcur %f %f %f grad %f %f %f \n",xcur[0],xcur[1],xcur[2]
			,gcur[0],gcur[1],gcur[2]);
	}while(iflag > 0);


  printf("\n\n-- Newton TNCG\n");
  double buff[200];
  dblAr1 buf(200,buff);
  truncated_newton_work work(nvar,buf);

  args.niter = 0;
  args.iprt =5;

  iflag = 0;
  for(int ii = 0; ii < 3; ii++) xcur[ii] = 1;

  do{
    optim_newton_drivertype_TNCG<nvar>(args ,work, xcur, &fcur, gcur, hess, &iflag, &ihess);
    if(ierro > 0){
      printf("ERROR %d\n",ierro);
      exit(1);
    }
    if(iflag <= 0){
      printf("Stop\n");
      break;
    }

    fcur = func(xcur,ihess,gcur,hess);
    printf("-- Iteration %d flag %d value = %e opt = %f \n",args.niter,iflag,fcur,fopt);
    printf("  - debug xcur %f %f %f grad %f %f %f \n",xcur[0],xcur[1],xcur[2]
      ,gcur[0],gcur[1],gcur[2]);
  }while(iflag > 0);

}









