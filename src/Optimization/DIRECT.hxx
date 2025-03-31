//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_OPTIMIZATION_DIRECT__
#define __METRIS_OPTIMIZATION_DIRECT__

#include "../MetrisRunner/MetrisParameters.hxx"
#include "../types.hxx"

namespace Metris{


// Options and work
struct DIBLOB_args{
  // Options
  int miter;
  // Stop if abs(fmin) < ftol * abs(fscale)
  double ftol;
  // Stop if fmin_{n} - fmin_{n+1} < dftol * abs(fscale)
  double dftol;
  double fscale; // Scale of f, e.g. initial value
  // When a given region provides fmin nloc_switch iters in a row
  // the algo switches to local by discarding all other elements and only 
  // subdividing the lowest value provider in the future.
  int nloc_switch; 

  // Work
  int niter, iflag;
  intAr2 ent2pol;
  dblAr2 coorl;
  dblAr1 fuelt, rhull;
  intAr1 lhull;
  intAr1 fminhist; // Track number of iters in a row an element has provided fmin (inherited from parent)
  //bool ifirst;
  //double fmin0;
  double fmin_pre;
  bool iloc_mode;

  DIBLOB_args() = delete;
  DIBLOB_args(const MetrisParameters &param_) : param(&param_) {
    reset();
    fscale= 1;
    ftol  = 1.0e-3;
    dftol = -1;
    nloc_switch = 0;
  }

  void reset(){
    iflag = 0;
    iloc_mode = false;
  }

  const MetrisParameters *param;
};

// Optimize among nblob barycentric domains, returning:
// - ifmin: the bary domain of the optimum
// - barmin: the barycentric coordinate of the optimum
// - fmin: the minimum value
// Intermediate steps provide/ask for:
// - leval: list of evaluations to carry out, in bary dom leval(ii,1)
// - peval: corresponding barycentric coordinates
// - feval: resulting evaluations
void DIBLOB(DIBLOB_args &args,
            int idim, int nblob, 
            intAr2 &leval, dblAr2 &peval, dblAr1 &feval,
            int *ifmin, double *barmin ,double *fmin);
void aux_DIRECT_splittet(DIBLOB_args &args,int ielem,int ieglo,int ilev);
void aux_DIRECT_splittri(DIBLOB_args &args,int ielem,int ieglo,int ilev);
void aux_DIRECT_newevals(DIBLOB_args &args, int idim, int nele0, int ilev, int ieglo,
                         intAr2 &leval, dblAr2 &peval, dblAr1 &feval);



}// end namespace
#endif