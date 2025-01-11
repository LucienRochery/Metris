//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

//#include "PETSc_helper.hxx"

#include "linalg/invmat.hxx"
#include "linalg/matprods.hxx"
#include "tools/minmax.h"
#include "types.hxx"
#include "opt_generic.hxx"


namespace Metris{


//    --- wlfc1 -> c1 , wlfc2 -> c2 for reference values ; assert c2 > c1
//    --- a greater c1 imposes a stricter decrease condition
//    --- a smaller c2 imposes a stricter conservation of curvature along iterations

//    -- Quick "see if this goes somewhere" params : large stepmin (ex. 1.0e-2)
//                                                   large c1 (ex. 1.0e-1)
//                                                   default c2
//    -- Set negative values to use defaults (1.0e-4 and 0.9d0)

//    rwork[0] = fpre
//    rwork[1] = step
//    rwork[2] = dot1 (gk,dk)
//    rwork(      3: nvar+3) = xpre
//    rwork( nvar+4:2nvar+3) = gpre
//    rwork(2nvar+4:3nvar+3) = desc
//    rwork(3nvar+4:7nvar+3) = previous descent directions
//    rwork(7nvar+4:10nvar+3)= previous iterates (at the end of LS)
//    rwork[10nvar+4-1] = first descent direction norm
//    nwork = 14 for nvar = 1
//            24 for nvar = 2
//            34 for nvar = 3

//    iwork[0] = index of oldest descent direction
//    iwork[1] = number of descent directions and iterates in memory
//    iwork[2] = n added since last correction

//    r = ratnew such that direction = r*(newton desc) + [1-r-1]*[-grad-1]
//    -> NOT TRUE (in 3d its the step ratio !! )

//    corrtrj works well for localization -> averages oscillating directions

// If isym > 0, the Hessian is given as :
// 1 2 4 
//   3 5
//     6

// Otherwise (e.g. Jacobian matrix), as [nvar][nvar]. If 0, hess is H, if -1, H^T. 

void optim_newton_drivertype(int nvar ,
                             double *xcur ,double *fcur  ,double *gcur   ,
                             double *hess ,
                             double xtol ,double stpmin,int isym ,
                             double wlfc1,double wlfc2 ,double ratnew ,
                             int *niter,int maxit ,int iprt   ,
                             int *iflag,int *ihess   , 
                             int nrwrk,double *rwork ,
                             int niwrk,int *iwork ,
                             double *xopt ,double *fopt ,int *ierro){

  *ierro = 0;
  
  const double alpha0 = 1.0;
  //const double c1 = 1.0e-4;
  const double c2 = 0.9;
  //const double stprat = 0.6;


  double epsilon = std::numeric_limits<double>::epsilon();

  double gnorm;

  int sg = 1; 


  if(*iflag == 0){
    if(nrwrk < 10*nvar+4) METRIS_THROW_MSG(DMemExcept(),"Increase nrwrk");

    *niter = 0;   
    *iflag = 1;
    *ierro = 0;
    *fopt  = 1.0e30;
    *ihess = 1;
    iwork[0] = 0;
    iwork[1] = 0;
    iwork[2] = 1;
    rwork[1] = alpha0;
    rwork[10*nvar+3] = -1.0;
    goto flag999;
  }

  if(iprt >= 3)printf(" - enter Newton niter %d fcur = %15.8f isym = %d\n",*niter,*fcur,isym);

  if(*fcur < *fopt){
    for(int ii = 0; ii < nvar; ii++) xopt[ii] = xcur[ii];
    *fopt = *fcur;
    //if(iprt >= 3) printf(" fopt update in newton algo %f \n",*fopt);
  }


  flag2000:
  if(*iflag == 1){

    //--- COMPUTE DESCENT DIRECTION
    if(nvar == 1){
      if(abs(hess[0]) > epsilon*(abs(gcur[0]))) {
        rwork[2*nvar+3]  = - gcur[0]/hess[0];
      }else{
        *ierro = 1;
        goto flag999;
      }
      gnorm = abs(rwork[2*nvar+3]);
    }else if(nvar == 2) {
      invsym<2>(hess);
      rwork[2*nvar+3] = -(hess[0]*gcur[0] + hess[1]*gcur[1]);
      rwork[2*nvar+4] = -(hess[1]*gcur[0] + hess[2]*gcur[1]);
      gnorm = MAX(abs(rwork[2*nvar+3]),
                  abs(rwork[2*nvar+4]));
    }else if(nvar == 3){
      if(isym > 0){
        invspd(3,hess);
        rwork[2*nvar+3] = -(hess[0]*gcur[0] + hess[1]*gcur[1] + hess[3]*gcur[2]);
        rwork[2*nvar+4] = -(hess[1]*gcur[0] + hess[2]*gcur[1] + hess[4]*gcur[2]);
        rwork[2*nvar+5] = -(hess[3]*gcur[0] + hess[4]*gcur[1] + hess[5]*gcur[2]);
      }else{
        // Actually a Jacobian, so not symmetric 
        invmat(3,hess);
        if(isym == 0){
          mat3vec(hess,gcur,&rwork[2*nvar+3]);
        }else{
          mat3vect(hess,gcur,&rwork[2*nvar+3]);
        }
        rwork[2*nvar+3] = -rwork[2*nvar+3];
        rwork[2*nvar+4] = -rwork[2*nvar+4];
        rwork[2*nvar+5] = -rwork[2*nvar+5];
      }
      gnorm = MAX(MAX(abs(rwork[2*nvar+3]),
                      abs(rwork[2*nvar+4])),
                      abs(rwork[2*nvar+5]));
    }

    if(*niter == 1){
      if(iprt >= 3)printf(" First dir norm %f \n",gnorm);
      rwork[10*nvar+4-1] = MAX(gnorm,1.0e-12);
    }else if(gnorm  <  xtol*rwork[10*nvar+4-1]) {
      if(iprt >= 3) printf(" debug gnorm termination \n");
      *iflag = 0;
      goto flag999;
    } 

  }

  if(*iflag == 1 || *iflag == 4){
    *ihess = 0;
    *iflag = 2;
    // fpre = fcur
    rwork[0] = *fcur;
    for(int ii = 0; ii < nvar; ii++) rwork[3+ii] = xcur[ii];

    // step = alpha0
    if(*niter == 1) {
      rwork[1] = alpha0;
    }else{
      rwork[1] = MIN(rwork[1] / ratnew, alpha0);
    }


    if(iprt >= 3){
      printf("-- start LS step = %15.7e dir = ",rwork[1]);
      dblAr1(nvar,&rwork[2*nvar+3]).print();
    } 

    rwork[2] = 0.0;
    for(int ii = 0; ii < nvar; ii++){
      rwork[2] += gcur[ii]*rwork[2*nvar+3+ii];
      xcur[ii] += sg*rwork[1]*rwork[2*nvar+3+ii];
    }
    goto flag999;

  }else if(*iflag == 2){

    double dot = rwork[2*nvar+3]*gcur[0];
    for(int ii = 1; ii < nvar; ii++){
      dot += rwork[2*nvar+3+ii]*gcur[ii];
    }

    double wc1 = wlfc1;
    double wc2 = wlfc2;
    //if(wlfc1  <  epsilon) wc1 = c1;
    if(wlfc2  <  epsilon) wc2 = c2;

    flag200:        
    if(wc1  >  wc2) {
      wc1 = MAX(wc1/2.0,wc2-1.0e-12);
      goto flag200;
    }

    if(iprt > 1) {
      printf(" (Newton) 1st cdt %15.7e <= %15.7e %d \n",
        *fcur, rwork[0] + wc1*rwork[1]*dot, ( *fcur      <=  (rwork[0] + wc1*rwork[1]*dot)) );

      //printf("dot = %f terms: rwork[2*nvar+3+:] = ",dot);
      //dblAr1(nvar,&rwork[2*nvar+3]).print();
      //printf("gcur = ");
      //dblAr1(nvar,gcur).print();
      //printf("\n");
      //printf("New eval %15.7e old %15.7e \n",*fcur,rwork[0]);
      //wait();

      printf(" (Newton) 2nd cdt %15.7e <= %15.7e %d \n",
        abs(dot),wc2*abs(rwork[2]),(abs(dot)  <=  wc2*abs(rwork[2])) );
      //printf(" (Newton) alt cdt %15.7e <= %15.7e %d \n",
      //  0.05,(rwork[0] - *fcur)/rwork[0], ((rwork[0] - *fcur)/rwork[0] > 0.05  ));
    }
    if(  ( *fcur      <=  (rwork[0] + wc1*rwork[1]*dot)
      &&   abs(dot)  <=  wc2*abs(rwork[2])             ) 
      //||  (rwork[0] - *fcur)/rwork[0] > 0.05  
      ) {
      if(iprt >= 3){
        printf(" ++ strong Wolfe conditions ok at xcur = ");
        dblAr1(nvar,xcur).print();
        printf("  relative decrease = %f\n",(*fcur - rwork[0])/rwork[0]);
      }
      *ierro = -1;
      if(*ihess <= 0) {
        *iflag = 1;
        *ihess = 1;
        goto flag999;
      }else{
        *iflag = 1;
        if(iprt >= 3) printf("++ since user set ihess to 1, it is assumed hessian is computed\n");
        if(iprt >= 3) printf("  -> thus going back to 1\n");
        goto flag2000;
      }
      goto flag999;
    }

    // --- Most times, the Wolfe condition is verified right away
    // -> this avoids having to compute f and grad twice to get hessian
    //    the second time around
    if(rwork[1] > sqrt(ratnew)*alpha0)  *ihess = 0;
    rwork[1] = rwork[1] * ratnew;
    if(rwork[1] < stpmin) {
      *iflag = 0;
      if(iprt >= 3) printf(" step < stepmin termination \n");
      goto flag999;
    }
    // xpre + step * desc
    //if(iprt >= 3 && nvar == 3) printf("New iterate from prev %f %f %f d = %f %f %f \n",rwork[3],rwork[4],rwork[5],
    //                                            rwork[3+2*nvar+0],rwork[3+2*nvar+1],rwork[3+2*nvar+2]);

    for(int ii = 0; ii < nvar; ii++){
      xcur[ii] = rwork[3+ii] + rwork[1]*rwork[3+2*nvar+ii];
    }

  }



  flag999:
  *niter += 1; 
  if(*niter  >  maxit) {
    if(iprt >= 3) printf(" newton max step exceeded %d \n",maxit);
    *iflag = 0;
  }

}






template <int nvar>
int optim_newton_drivertype(newton_drivertype_args<nvar> &args,
                            double *xcur ,double *fcur ,
                            double *gcur ,double *hess ,
                            int *iflag, int *ihess){

  int ierro = 0;

  const double alpha0 = 1.0;
  //const double c1 = 1.0e-4;
  const double c2 = 0.9;
  //const double stprat = 0.6;


  double epsilon = std::numeric_limits<double>::epsilon();

  double gnorm;

  int sg = 1; 


  if(*iflag == 0){
    args.niter = 0;   
    *iflag = 1;
    args.fopt  = 1.0e30;
    *ihess = 1;
    args.iwork[0] = 0;
    args.iwork[1] = 0;
    args.iwork[2] = 1;
    args.rwork[1] = alpha0;
    args.rwork[10*nvar+3] = -1.0;
    goto flag999;
  }

  if(args.iprt >= 3)printf(" - enter Newton niter %d fcur = %15.8f isym = %d\n",args.niter,*fcur,args.isym);

  if(*fcur < args.fopt){
    for(int ii = 0; ii < nvar; ii++) args.xopt[ii] = xcur[ii];
    args.fopt = *fcur;
    if(args.iprt >= 3) printf(" fopt update in newton algo %f \n",args.fopt);
  }


  flag2000:
  if(*iflag == 1){

    //--- COMPUTE DESCENT DIRECTION
    if(nvar == 1){
      if(abs(hess[0]) > epsilon*(abs(gcur[0]))) {
        args.rwork[2*nvar+3]  = - gcur[0]/hess[0];
      }else{
        ierro = 1;
        goto flag999;
      }
      gnorm = abs(args.rwork[2*nvar+3]);
    }else if(nvar == 2) {
      invsym<2>(hess);
      args.rwork[2*nvar+3] = -(hess[0]*gcur[0] + hess[1]*gcur[1]);
      args.rwork[2*nvar+4] = -(hess[1]*gcur[0] + hess[2]*gcur[1]);
      gnorm = MAX(abs(args.rwork[2*nvar+3]),
                  abs(args.rwork[2*nvar+4]));
    }else if(nvar == 3){
      if(args.isym > 0){
        invspd(3,hess);
        args.rwork[2*nvar+3] = -(hess[0]*gcur[0] + hess[1]*gcur[1] + hess[3]*gcur[2]);
        args.rwork[2*nvar+4] = -(hess[1]*gcur[0] + hess[2]*gcur[1] + hess[4]*gcur[2]);
        args.rwork[2*nvar+5] = -(hess[3]*gcur[0] + hess[4]*gcur[1] + hess[5]*gcur[2]);
      }else{
        // Actually a Jacobian, so not symmetric 
        invmat(3,hess);
        if(args.isym == 0){
          mat3vec(hess,gcur,&args.rwork[2*nvar+3]);
        }else{
          mat3vect(hess,gcur,&args.rwork[2*nvar+3]);
        }
        args.rwork[2*nvar+3] = -args.rwork[2*nvar+3];
        args.rwork[2*nvar+4] = -args.rwork[2*nvar+4];
        args.rwork[2*nvar+5] = -args.rwork[2*nvar+5];
      }
      gnorm = MAX(MAX(abs(args.rwork[2*nvar+3]),
                      abs(args.rwork[2*nvar+4])),
                      abs(args.rwork[2*nvar+5]));
    }

    if(args.niter == 1){
      if(args.iprt >= 3)printf(" First dir norm %f \n",gnorm);
      args.rwork[10*nvar+4-1] = MAX(gnorm,1.0e-12);
    }else if(gnorm < args.xtol*args.rwork[10*nvar+4-1]) {
      if(args.iprt >= 3) printf(" debug gnorm termination \n");
      *iflag = 0;
      goto flag999;
    } 

  }

  if(*iflag == 1 || *iflag == 4){
    *ihess = 0;
    *iflag = 2;
    // fpre = fcur
    args.rwork[0] = *fcur;
    for(int ii = 0; ii < nvar; ii++) args.rwork[3+ii] = xcur[ii];

    // step = alpha0
    if(args.niter == 1) {
      args.rwork[1] = alpha0;
    }else{
      args.rwork[1] = MIN(args.rwork[1] / args.ratnew, alpha0);
    }


    if(args.iprt >= 3){
      printf("-- start LS step = %15.7e dir = ",args.rwork[1]);
      dblAr1(nvar,&args.rwork[2*nvar+3]).print();
    } 

    args.rwork[2] = 0.0;
    for(int ii = 0; ii < nvar; ii++){
      args.rwork[2] += gcur[ii]*args.rwork[2*nvar+3+ii];
      xcur[ii] += sg*args.rwork[1]*args.rwork[2*nvar+3+ii];
    }
    goto flag999;

  }else if(*iflag == 2){

    double dot = args.rwork[2*nvar+3]*gcur[0];
    for(int ii = 1; ii < nvar; ii++){
      dot += args.rwork[2*nvar+3+ii]*gcur[ii];
    }

    double wc1 = args.wlfc1;
    double wc2 = args.wlfc2;
    //if(wlfc1  <  epsilon) wc1 = c1;
    if(args.wlfc2  <  epsilon) wc2 = c2;

    flag200:        
    if(wc1  >  wc2) {
      wc1 = MAX(wc1/2.0,wc2-1.0e-12);
      goto flag200;
    }

    if(args.iprt > 1) {
      printf(" (Newton) 1st cdt %15.7e <= %15.7e %d \n",
        *fcur, args.rwork[0] + wc1*args.rwork[1]*dot, ( *fcur      <=  (args.rwork[0] + wc1*args.rwork[1]*dot)) );

      //printf("dot = %f terms: args.rwork[2*nvar+3+:] = ",dot);
      //dblAr1(nvar,&args.rwork[2*nvar+3]).print();
      //printf("gcur = ");
      //dblAr1(nvar,gcur).print();
      //printf("\n");
      //printf("New eval %15.7e old %15.7e \n",*fcur,args.rwork[0]);
      //wait();

      printf(" (Newton) 2nd cdt %15.7e <= %15.7e %d \n",
        abs(dot),wc2*abs(args.rwork[2]),(abs(dot)  <=  wc2*abs(args.rwork[2])) );
      //printf(" (Newton) alt cdt %15.7e <= %15.7e %d \n",
      //  0.05,(args.rwork[0] - *fcur)/args.rwork[0], ((args.rwork[0] - *fcur)/args.rwork[0] > 0.05  ));
    }
    if(  ( *fcur      <=  (args.rwork[0] + wc1*args.rwork[1]*dot)
      &&   abs(dot)  <=  wc2*abs(args.rwork[2])             ) 
      //||  (args.rwork[0] - *fcur)/args.rwork[0] > 0.05  
      ){
      if(args.iprt >= 3){
        printf(" ++ strong Wolfe conditions ok at xcur = ");
        dblAr1(nvar,xcur).print();
        printf("  relative decrease = %f\n",(*fcur - args.rwork[0])/args.rwork[0]);
      }
      ierro = -1;
      if(*ihess <= 0) {
        *iflag = 1;
        *ihess = 1;
        goto flag999;
      }else{
        *iflag = 1;
        if(args.iprt >= 3) printf("++ since user set ihess to 1, it is assumed hessian is computed\n");
        if(args.iprt >= 3) printf("  -> thus going back to 1\n");
        goto flag2000;
      }
      goto flag999;
    }

    // --- Most times, the Wolfe condition is verified right away
    // -> this avoids having to compute f and grad twice to get hessian
    //    the second time around
    if(args.rwork[1] > sqrt(args.ratnew)*alpha0)  *ihess = 0;
    args.rwork[1] = args.rwork[1] * args.ratnew;
    if(args.rwork[1] < args.stpmin) {
      *iflag = 0;
      if(args.iprt >= 3) printf(" step < stepmin termination \n");
      goto flag999;
    }
    // xpre + step * desc
    //if(args.iprt >= 3 && nvar == 3) printf("New iterate from prev %f %f %f d = %f %f %f \n",args.rwork[3],args.rwork[4],args.rwork[5],
    //                                            args.rwork[3+2*nvar+0],args.rwork[3+2*nvar+1],args.rwork[3+2*nvar+2]);

    for(int ii = 0; ii < nvar; ii++)
      xcur[ii] = args.rwork[3+ii] + args.rwork[1]*args.rwork[3+2*nvar+ii];

  }


  flag999:
  args.niter++; 
  if(args.niter  >  args.maxit) {
    if(args.iprt >= 3) printf(" newton max step exceeded %d \n",args.maxit);
    *iflag = 0;
  }

  return ierro;
}
template
int optim_newton_drivertype<1>(newton_drivertype_args<1> &args,
                               double *xcur ,double *fcur ,
                               double *gcur ,double *hess ,
                               int *iflag, int *ihess);
template
int optim_newton_drivertype<2>(newton_drivertype_args<2> &args,
                               double *xcur ,double *fcur ,
                               double *gcur ,double *hess ,
                               int *iflag, int *ihess);
template
int optim_newton_drivertype<3>(newton_drivertype_args<3> &args,
                               double *xcur ,double *fcur ,
                               double *gcur ,double *hess ,
                               int *iflag, int *ihess);







#if 0

/*
rwork is reduced to 4 values:
 - rwork[0] : fcur at descent direction start
 - rwork[1] : current step 
 - rwork[2] : dot product between current gradient (L-S) and descent direction
 - rwork[3] : initial gradient norm for termination condition
Previous iterate is no longer stored, and the descent direction is stored explicitely.
*/
#ifdef USE_PETSC
int optim_newton_drivertype_PETSc(int nvar ,
                             Vec &XCUR   ,double *fcur  ,Vec &RHS, Mat &OJ, 
                             double xtol ,double stpmin, 
                             double wlfc1,double wlfc2 ,double ratnew ,
                             int *niter,int maxit ,int iprt   ,
                             int *iflag,int *ijaco , 
                             double *rwork,
                             Vec &DESC){

  const double alpha0 = 1.0;
  const double c1 = 1.0e-4;
  const double c2 = 0.9;
  const double stprat = 0.6;

  int ierro = 0;
  double epsilon = std::numeric_limits<double>::epsilon();
  const double minvp = epsilon;

  double eigval[3], eigvec[9];
  PetscReal gnorm;
  double rbuff[10];


  if(*iflag == 0){

    if(wlfc1 > wlfc2)
      METRIS_THROW_MSG(WArgExcept(),"Wolfe constants c1 < c2 not verified")


    *niter = 0;   
    *iflag = 1;
    ierro = 0;
    *ijaco = 1;
    rwork[1] = alpha0;
    rwork[3] = -1.0;
    goto flag999;
  }


  flag2000:
  if(*iflag == 1){
    // Compute descent direction. 
    // d = A^{-1}b 
    // -> ||dd|| <= |||dA^{-1}||| ||b|| 
    // We actually don't need the matrix norm of A^{-1}. 
    // Knowing rhs inf-norm, we can set |||dA^{-1}||| (abstol)
    // <= epsilon to guarantee ||dd|| 
    KSP ksp;
//    double ksp_tol = MAX(1.0e-5, rwork[1]/100);
    //double ksp_tol = MAX(1.0e-3, rwork[1]/100);
    double ksp_tol = MAX(1.0e-2, rwork[1]/100);
//    double ksp_tol = MAX(1.0e-5, rwork[1]/100);
    if(iprt >= 3) printf("Calling linear solver with tolerance %15.7e\n",ksp_tol);
    PetscKSPHelper KSOLVER(ksp, PETSC_COMM_SELF, OJ, PETSC_DEFAULT, ksp_tol);
    double t0 = get_wall_time();
    KSOLVER.solve(RHS,DESC);
    double t1 = get_wall_time();
    if(iprt >= 3) printf("KSP time %f \n",t1-t0);

    // Get current descent direction norm for termination condition
    PetscCall(VecNorm(DESC,NORM_2,&gnorm));
    if(iprt > 1) printf("Debug gnorm = %f\n",gnorm);

    if(*niter == 1) {
      if(iprt >= 3)printf(" First dir norm %f \n",gnorm);
      rwork[3] = gnorm;
    }else if(gnorm  <  xtol*rwork[3]) {
      if(iprt >= 3) printf(" debug gnorm termination \n");
      *iflag = 0;
      goto flag999;
    }
  }

  if(*iflag == 1 || (*iflag == 4)){

    if(iprt >= 3) printf("------> start line search \n");
    *iflag = 2;
    *ijaco = 0;
    // Store current function value
    rwork[0] = *fcur;
    printf("Debug set fcur to %f \n",*fcur);

    // Initialize line search. 
    // First is done using hard-coded alpha0. 
    // Subsequent line searches increase previous acceptable step by factor ratnew.
    if(*niter == 1) {
      rwork[1] = alpha0;
    }else{
      rwork[1] = MIN(rwork[1] / ratnew, alpha0);
    }
    // This dot product is used in the line search (curvature condition)
    PetscCall(VecDot(RHS,DESC,&rwork[2]));
    // Set xcur = xcur + rwork[1]*DESC
    PetscCall(VecAXPY(XCUR,rwork[1],DESC));
//    for(int ii = 0; ii < nvar; ii++){
//      xcur[ii] += sg*rwork[1]*dDESC[ii];
//    }
    goto flag999;


  }else if(*iflag == 2){

    double wc1 = wlfc1;
    double wc2 = wlfc2;
    if(wlfc1  <  epsilon) wc1 = c1;
    if(wlfc2  <  epsilon) wc2 = c2;

    //PetscVecArray dDESC(DESC);

    // This always fails the Wolfe condition (in case Jacobian failed)
    PetscScalar dot = *fcur/(wc1*rwork[1]) - rwork[0] + 1;
    if(*ijaco > 0){
      PetscCall(VecDot(RHS,DESC,&dot));
    }
    //double dot = dDESC[0]*RHS[0];;
    //for(int ii = 1; ii < nvar; ii++){
    //  dot += dDESC[ii]*gcur[ii];
    //}


    //while(wc1 > wc2){
    //  wc1 = MAX(wc1/2.0,wc2-1.0e-12);
    //}

    // Verify Wolfe conditions 
    if(iprt > 1) {
      printf(" (Newton) 1st cdt %15.7e <= %15.7e %d \n",
        *fcur, rwork[0] + wc1*rwork[1]*dot, ( *fcur      <=  (rwork[0] + wc1*rwork[1]*dot)) );
      printf(" (Newton) 2nd cdt %15.7e <= %15.7e %d \n",
        abs(dot),wc2*abs(rwork[2]),(abs(dot)  <=  wc2*abs(rwork[2])) );
      printf(" (Newton) alt cdt %15.7e <= %15.7e %d \n",
        0.15,(rwork[0] - *fcur)/rwork[0], ((rwork[0] - *fcur)/rwork[0] > 0.15  ));
    }
    if(  ( *fcur      <=  (rwork[0] + wc1*rwork[1]*dot)
      &&   abs(dot)  <=  wc2*abs(rwork[2])             ) 
      ||  (rwork[0] - *fcur)/rwork[0] > 0.15  ) {
      if(iprt >= 3) printf(" -> Step %f satisfies Wolfe conditions\n",rwork[1]);
      if(iprt >= 3) printf("  relative decrease = %f\n",(*fcur - rwork[0])/rwork[0]);
      ierro = -1;
      *iflag = 1;
      if(*ijaco <= 0){
        *ijaco = 1;
        goto flag999;
      }else{
        // If the Jacobian had been requested previously (e.g. first step success)
        // then restart immediately. 
        (*niter) += 1;
        goto flag2000;
      }
    }

    // Might have been reset by costfunc failure
    *ijaco = 1;

    // Wolfe conditions not verified. 
    // rwork[1] holds the current step. 
    if(rwork[1] > sqrt(ratnew)*alpha0)  *ijaco = 0;
    rwork[1] = rwork[1] * ratnew;
    if(rwork[1] < stpmin) {
      *iflag = 0;
      if(iprt >= 3) printf(" step < stepmin termination \n");
      goto flag999;
    }

    // Set new abscissa. 
    // First update is simply xcur + step * DESC
    // Afterwards, in order not to store xpre, we can simply do:
    // xcur = xcur + (stepnew - stepold) * DESC
    // This should be fine for steps not too close to machine zero. 
    PetscScalar stpeff = rwork[1]*(1 - 1/ratnew);
    // AXPY(y,a,x) : y <- y + a * x
    PetscCall(VecAXPY(XCUR, stpeff, DESC));
//    for(int ii = 0; ii < nvar; ii++){
//      xcur[ii] += stpeff * dDESC[ii];
//    }

  }



  flag999:
  *niter += 1; 
  if(*niter  >  maxit) {
    if(iprt >= 3) printf(" newton max step exceeded %d \n",maxit);
    *iflag = 0;
  }

  return ierro;
}

#endif
#endif //if 0

double brutal_samplingsearch(std::function<double (double)> func, double xmin, double xmax, int nsamp, int nrep){
  double xmi1 = xmin; 
  double xma1 = xmax;
  double xopt = xmin;
  for(int irep = 0; irep < nrep; irep++){
    double fopt = func(xmi1);
    xopt = xmi1;
    for(int isamp = 0; isamp < nsamp; isamp++){
      double xcur = (isamp*xma1 + (nsamp - 1 - isamp) * xmi1) / (nsamp - 1);
      double fcur = func(xcur);
      if(fcur < fopt){
        xopt = xcur;
        fopt = fcur;
      }
    }
    xmi1 = xopt - 1.0 / (nsamp - 1);
    xma1 = xopt + 1.0 / (nsamp - 1);
  }
  return xopt;
}



} // End namespace

