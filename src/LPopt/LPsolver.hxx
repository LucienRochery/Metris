//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_LPSOLVER__
#define __METRIS_LPSOLVER__

#include "../libs/alglib-cpp/src/optimization.h"

#include "../types.hxx"

#ifdef USE_CLP
#include "ClpSimplex.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#include "ClpInterior.hpp"
#include "ClpSolve.hpp"
#include "ClpCholeskyBase.hpp"
#endif


namespace Metris{

class MeshBase;

enum class LPMethod{Simplex, IPM};
enum class LPLib{alglib, clp};


//template<LPLib ilib>
// Let's just use runtime checks. The vast majority of time will be spent solving
// anyways, and branch prediction will largely offset any ifs. 
class LPsolver{

public:
  // Only through the constructor can the lib be set (not to be changed after)
  LPsolver(LPLib ilib, LPMethod method);

  void allocate(int nrow, int ncol);

  // To read/modify the objective: solver.obj(i) = x
  double& obj(int i);

  void setConstraintMatrix(int i, int j, double x);
  
  // void setVariableConstraints(int noptim_points, int ncoefglob, bool metricOn);

  void setVarCstr(int i, double lb, double ub);

  void setConstraintLB(int i, double x);
  void setConstraintUB(int i, double x);

  void setConstraintLBinf();
  void setConstraintUBinf();
  // Returns the Value of the Objective Function
  double optimize();

  // Sets coordinates to what's provided by x_ALG or model (CLP) times fac. 
  // Only change coordinate icoor 
  double updateCoord(MeshBase &msh, const intAr1 &idx_point, int icoor, double fac) const;

  alglib::minlpstate state;


  double INFINITY_m;
  double INFINITY_p;
public:
  LPMethod method;
  LPLib    ilib;

  // Alglib
  alglib::sparsematrix Acstr_ALG;
  alglib::real_1d_array LB_ALG, UB_ALG;
  alglib::real_1d_array obj_ALG, x_ALG;

  // CLP 
  #ifdef USE_CLP
  CoinPackedMatrix Acstr_CLP;
  // ClpSimplex model;
  ClpInterior model;
  std::vector<double> LB_CLP, UB_CLP;
  std::vector<double> CL, CU, obj_CLP;
  #endif

  int nrow, ncol;
  // alglib::minlpstate state;

};




} // namespace Metris
#endif