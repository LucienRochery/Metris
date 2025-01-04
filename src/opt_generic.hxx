//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_OPT_GENERIC__
#define __METRIS_OPT_GENERIC__

#ifdef USE_PETSC
  #include <petscmat.h>
#endif

#include <functional>

namespace Metris{


void optim_newton_drivertype(int nvar ,
                             double *xcur ,double *fcur  ,double *gcur   ,double *hess ,
                             double xtol ,double stpmin,int isym,
                             double wlfc1,double wlfc2 ,double ratnew ,
                             int *niter,int maxit ,int iprt   ,
                             int *iflag,int *ihess , 
                             int nrwrk,double *rwork ,
                             int niwrk,int *iwork ,
                             double *xopt ,double *fopt ,int *ierro);


#ifdef USE_PETSC
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

double brutal_samplingsearch(std::function<double(double)> func, double xmin, double xmax, int nsamp, int nrep);


} // End namespace


#endif