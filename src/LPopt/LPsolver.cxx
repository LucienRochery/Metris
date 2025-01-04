//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "LPsolver.hxx"
#include "../Mesh/MeshBase.hxx"

namespace Metris{

// Only through the constructor can the lib be set (not to be changed after)
LPsolver::LPsolver(LPLib ilib, LPMethod method){
  // METRIS_THROW_MSG(TODOExcept(), "Implement LPsolver class for lib " << (int)ilib);
  this->ilib   = ilib;
  this->method = method;
  nrow = 0; 
  ncol = 0;
  #ifndef USE_CLP
    if(ilib == LPLib::clp) 
      METRIS_THROW_MSG(WArgExcept(), "Recompile with USE_CLP=True")
  #endif

  if(ilib == LPLib::alglib){
    INFINITY_m = alglib::fp_neginf;
    INFINITY_p = alglib::fp_posinf;
  }else{
    METRIS_THROW_MSG(TODOExcept(), "Set infinities for CLP");
  }
}

void LPsolver::allocate(int nrow, int ncol){

  METRIS_ENFORCE(this->nrow == 0 && this->ncol == 0);
  METRIS_ENFORCE(nrow != 0 && ncol != 0);

  this->nrow = nrow;
  this->ncol = ncol;

  using namespace alglib;
  minlpcreate(ncol, state);

  //using matrix = std::conditional_t<ilib == LPLib::alglib, 
  //                                  sparsematrix,
  //                                  CoinPackedMatrix>;

  //if constexpr(std::is_same<ilib,LPLib::alglib>::value){
  if(ilib == LPLib::alglib){

    sparsecreate(nrow, ncol, Acstr_ALG);
    LB_ALG.setlength(nrow);
    UB_ALG.setlength(nrow);
    obj_ALG.setlength(ncol);
    /**************************************************************************
      void minlpsetcost(minlpstate &state, const real_1d_array &c, 
       const xparams _xparams = alglib::xdefault);
      This function sets cost term for LP solver.
      By default, cost term is zero. ** In reality this does not happen :(
      INPUT PARAMETERS:
          State   -   structure which stores algorithm state
          C       -   cost term, array[N].
        -- ALGLIB --
          Copyright 19.07.2018 by Bochkanov Sergey
      *************************************************************************/

  }else if(ilib == LPLib::clp){

    #ifdef USE_CLP
    Acstr_CLP.setDimensions(nrow, ncol);
    CL.resize(ncol);
    CU.resize(ncol);
    LB_CLP.resize(nrow);
    UB_CLP.resize(nrow);
    //obj_CLP.resize(ncol, 0.0);
    obj_CLP.resize(ncol);
    #endif

  }

}

// To read/modify the objective: solver.obj(i) = x
double& LPsolver::obj(int i){
  METRIS_ASSERT(i >= 0 && i < ncol);
  if(ilib == LPLib::alglib){
    return obj_ALG[i];
  }else if(ilib == LPLib::clp){
    #ifdef USE_CLP
    return obj_CLP[i];
    #endif
  }
  METRIS_THROW_MSG(WArgExcept(),"Incorrect lib LPsolver");
  return obj_ALG[i];
}

void LPsolver::setConstraintMatrix(int i, int j, double x){
  METRIS_ASSERT(i >= 0 && i < nrow && j >= 0 && j < ncol);
  if(ilib == LPLib::alglib){
    sparseset(Acstr_ALG, i, j, x); 
  }else if(ilib == LPLib::clp){
    #ifdef USE_CLP
    Acstr_CLP.modifyCoefficient(i, j, x);
    #endif
  }
}

void LPsolver::setVarCstr(int i, double lb, double ub){
  if(ilib == LPLib::alglib){
    minlpsetbci(state, i, lb, ub);
  }else if(ilib == LPLib::clp){
    METRIS_THROW_MSG(TODOExcept(),"Implement setVarCstr in CLP case");
    #ifdef USE_CLP
    LB_CLP[i] = x;
    #endif
  }
}
void LPsolver::setConstraintLB(int i, double x){
  if(ilib == LPLib::alglib){
    LB_ALG[i] = x;
  }else if(ilib == LPLib::clp){
    #ifdef USE_CLP
    LB_CLP[i] = x;
    #endif
  }
}

void LPsolver::setConstraintUB(int i, double x){
  if(ilib == LPLib::alglib){
    UB_ALG[i] = x;
  }else if(ilib == LPLib::clp){
    #ifdef USE_CLP
    UB_CLP[i] = x;
    #endif
  }
}

void LPsolver::setConstraintLBinf(){
  if(ilib == LPLib::alglib){
    for(int ii = 0; ii < nrow; ii++) LB_ALG[ii] = INFINITY_m;
  }else if(ilib == LPLib::clp){
    #ifdef USE_CLP
    for(int ii = 0; ii < nrow; ii++) LB_CLP[ii] = -COIN_DBL_MAX;
    #endif
  }
}
void LPsolver::setConstraintUBinf(){
  if(ilib == LPLib::alglib){
    for(int ii = 0; ii < nrow; ii++) UB_ALG[ii] = INFINITY_p;
  }else if(ilib == LPLib::clp){
    #ifdef USE_CLP
    for(int ii = 0; ii < nrow; ii++) UB_CLP[ii] = COIN_DBL_MAX;
    #endif
  }
}
// For now, subject to change as class shall be abstacted from metric problem
// void LPsolver::setVariableConstraints(int noptim_points, int ncoefglob, bool metricOn){
//   if(ilib == LPLib::alglib){
//     if (metricOn){
//       minlpsetbci(state, ncol-1, 0, INFINITY);
//       for (int i = 0; i<noptim_points; i++) minlpsetbci(state, i, -INFINITY, INFINITY); 
//       for (int i = noptim_points; i<noptim_points+2*ncoefglob; i++) minlpsetbci(state, i, 0, INFINITY);
//     }else minlpsetbcall(state, -INFINITY, INFINITY);
//   }else if(ilib == LPLib::clp){
//     #ifdef USE_CLP
//     if(metricOn){
//       CL.resize(ncol, 0.0);
//       CU.resize(ncol, 0.0);
//       for (int i = 0; i<noptim_points; i++){
//         CL[i] = -COIN_DBL_MAX;
//         CU[i] =  COIN_DBL_MAX;
//       }
//       for (int i = noptim_points; i<noptim_points+2*ncoefglob+1; i++){
//         CL[i] = 0;
//         CU[i] = COIN_DBL_MAX;
//       }
//       LB_CLP.resize(nrow, 0.0);
//       UB_CLP.resize(nrow, 0.0);
//     }else{
//       CL.resize(ncol, -COIN_DBL_MAX);
//       CU.resize(ncol, COIN_DBL_MAX);
//       LB_CLP.resize(nrow, 0.0);
//       UB_CLP.resize(nrow, COIN_DBL_MAX);
//     }
//     #endif
//   }
// }

double LPsolver::optimize(){
  if(ilib == LPLib::alglib){
    using namespace alglib;
    minlpreport rep;
    minlpsetcost(state, obj_ALG);
    minlpsetlc2(state, Acstr_ALG, LB_ALG, UB_ALG, nrow);
    if(method == LPMethod::Simplex){
      alglib::trace_file("DSS.DETAILED,PREC.F15", "LPoutDSS.log");
      minlpsetalgodss(state, 1e-7);
    }else{
      alglib::trace_file("IPM.DETAILED,PREC.F15", "LPoutIPM.log");
      alglib::minlpsetalgoipm(state, 1e-5);
    }
    minlpoptimize(state);
    minlpresults(state, x_ALG, rep);
    // Won't be throwing later on, instead handle the error
    if(rep.terminationtype < 1) METRIS_THROW_MSG(AlgoExcept(),
                        "LP failed to be solved code = "<<rep.terminationtype)
    return rep.f;

  }else if(ilib == LPLib::clp){
    #ifdef USE_CLP
    model.loadProblem(Acstr_CLP, &CL[0], &CU[0], &obj_CLP[0], &LB_CLP[0], &UB_CLP[0]);
    // model.loadProblem(Acstr_CLP, NULL, NULL, &obj_CLP[0], &LB_CLP[0], &UB_CLP[0]);
    model.setLogLevel(0);
    model.setOptimizationDirection(1); // -1:maximize, +1:minimize:
    model.scaling(1);
    model.setMaximumSeconds(1.0);
    model.setPrimalTolerance(1e-2);
    model.setDualTolerance(1e-2);

    if(method == LPMethod::IPM){ 
      // TODO:: Implement Cholesky Factorization for performance
      //        Does not seem to be automatically handled by ClpInterior and Barrier Method
      model.setSolveType(ClpSolve::useBarrier);
      model.setMaximumBarrierIterations(100);

    }else if (method == LPMethod::Simplex){
      model.setSolveType(ClpSolve::useDual); 
    }
    model.primalDual();
    if (model.isProvenOptimal() == 0) METRIS_THROW_MSG(AlgoExcept(),
                "LP failed to be solved. Termination code = " <<model.status());
    #endif
  }
  return 0;
}


double LPsolver::updateCoord(MeshBase &msh, const intAr1 &idx_point, int icoor, double fac) const{
  const double* dx = ilib == LPLib::alglib ? &x_ALG[0] :
                                           #ifdef USE_CLP
                                           model.primalColumnSolution()
                                           #else
                                           NULL
                                           #endif
                                           ;
  double nrmdx = 0;
  for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
    int irank = idx_point[ipoin];
    if(irank < 0) continue;
    msh.coord[ipoin][icoor] += fac*dx[irank];
    nrmdx += dx[irank]*dx[irank];
    //printf("DEBUG dx[%d] = %15.6e\n",irank,dx[irank]);
  }
  return sqrt(nrmdx); 
}




} // namespace Metris
