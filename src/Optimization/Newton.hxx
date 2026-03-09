//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_OPTIMIZATION_NEWTON__
#define __METRIS_OPTIMIZATION_NEWTON__

#ifdef METRIS_USE_PETSC
  #include <petscmat.h>
#endif

#include <functional>
#include <nlopt.h>

#include "../MetrisRunner/MetrisParameters.hxx"

namespace Metris{

template<int nvar>
struct newton_drivertype_args{
  double xtol, stpmin, wlfc1, wlfc2, ratnew, ftol;
  int niter, maxit, isym;
  int iwork[3];
  double rwork[4 + nvar*10];
  double xopt[nvar];
  double fopt;

  double fpre; // work


  newton_drivertype_args() = delete;

  newton_drivertype_args(MetrisParameters *param_) : param(param_){
    xtol   = -1;
    stpmin = 1.0e-3;
    wlfc1  = 0.1;
    wlfc2  = 10.0;
    ratnew = 0.5;
    ftol   = -1;

    niter = 0;
    maxit = 50;
    isym  = 1;
  }

  const MetrisParameters* param;
};
struct truncated_newton_work{
  // Using a dblAr1 should be safer (throw if undersized)
  truncated_newton_work(int ndim, dblAr1 &buf):
  pminor(ndim,&buf[0]), rminor(ndim,&buf[ndim]), dminor(ndim,&buf[2*ndim]),
  qminor(ndim,&buf[3*ndim]), pmajor(ndim, &buf[4*ndim])
  {
    METRIS_ASSERT(buf.get_n() >= 5*ndim);
  }
  dblAr1 pminor, rminor, dminor, qminor;
  dblAr1 pmajor;
};

void optim_newton_drivertype(const MetrisParameters &params,
                             int nvar ,
                             double *xcur ,double *fcur  ,double *gcur   ,double *hess ,
                             double xtol ,double stpmin,int isym,
                             double wlfc1,double wlfc2 ,double ratnew ,
                             int *niter,int maxit ,
                             int *iflag,int *ihess ,
                             int nrwrk,double *rwork ,
                             int niwrk,int *iwork ,
                             double *xopt ,double *fopt ,int *ierro);

template <int nvar,bool inBoundary=false>
int optim_newton_drivertype(newton_drivertype_args<nvar> &args,
                            double *xcur ,double *fcur ,
                            double *gcur ,double *hess ,
                            int *iflag, int *ihess);


template <int nvar>
int optim_newton_drivertype_TNCG(newton_drivertype_args<nvar> &args,
                                 truncated_newton_work &work,
                                 double *xcur ,double *fcur ,
                                 double *gcur ,double *hess ,
                                 int *iflag, int *ihess);

#ifdef METRIS_USE_PETSC
int optim_newton_drivertype_PETSc(int nvar ,
                             Vec &XCUR ,double *fcur  ,
                             Vec &RHS, Mat &OJ,
                             double xtol ,double stpmin,
                             double wlfc1,double wlfc2 ,double ratnew ,
                             int *niter,int maxit ,int iprt   ,
                             int *iflag,int *ijaco ,
                             double *rwork ,
                             Vec &DESC);
#endif


template<int ndim>
int truncated_newton_iteration(const MetrisParameters *params,
                               truncated_newton_work &work,
                               int outer_iter,
                               const double *gcur, const double *hcur,
                               double *desc);

double brutal_samplingsearch(std::function<double(double)> func, double xmin, double xmax, int nsamp, int nrep);


#define LUKSAN_PNET_MAXIT 100

// No mallocs in this one
template<int ndim>
nlopt_result luksan_pnetS(nlopt_func f, void *f_data,
                          const double *lb, const double *ub, /* bounds */
                          double *x, /* in: initial guess, out: minimizer */
                          double *fopt,
                          //int mf, /* subspace dimension (0 for default) */
                          nlopt_algorithm algorithm,
                          dblAr1 &lwork,
                          double fstop , double ftol_rel, double ftol_abs);

int luksan_pnet_worksize(int n);


} // End namespace


#endif