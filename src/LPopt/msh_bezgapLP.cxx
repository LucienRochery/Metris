//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_maxccoef.hxx"
#include "LPsolver.hxx"

#include "../codegen_ccoef_d.hxx"
#include "../low_geo.hxx"
#include "../ho_constants.hxx"
#include "../io_libmeshb.hxx"

#include "../aux_timer.hxx"
#include "../aux_utils.hxx"

#include "../Mesh/MeshBase.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"


#ifdef USE_CLP
#include "ClpInterior.hpp"
#include "ClpSolve.hpp"
#include "ClpCholeskyBase.hpp"
#include "ClpSimplex.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#endif

using namespace alglib;


namespace Metris{

// pos_ctrlp is size 1:nvar, 1:gdim
// nvar is the number of >= 0 entries in idx_point 
// idx_point >= 0 entries point into pos_ctrlp 
// LP is:
// min t 
// dDi - dXi >= Pi - Pi* - |Pi - Pi*| (1)
// dDi + dXi >= Pi* - Pi - |Pi - Pi*| (2)
// dNj dX >= j0Vj - Nj                (3)
// t >= dDi - |Pi - Pi*|              (4)
template<int gdim, int tdim, int ideg>
double bezGapsLP(MeshBase &msh, const intAr1 &idx_point, 
                 const dblAr2 &pos_ctrlp, //const dblAr1 &weight,
                 LPMethod method,  LPLib lib){
  // METRIS_THROW_MSG(TODOExcept(), "maximizeCcoef not implemented for ideg = "<<ideg);

  const bool MAE = false;
  // We can convert the mesh, but doing that for the displacements would be a 
  // pain in the ass. Better the caller works in Bezier from the start. 
  METRIS_ENFORCE(msh.getBasis() == FEBasis::Bezier);

  constexpr int jdeg = tdim * (ideg - 1);
  constexpr int ncoef = tdim == 2 ? facnpps[jdeg] : tetnpps[jdeg];
  // Full node count 
  constexpr int nnode = tdim == 2 ? facnpps[ideg] : tetnpps[ideg];
  // Only P1 nodes 
  constexpr int vol0 = ifact<tdim>();

  const int iverb = msh.param->iverb;

  const intAr2 &ent2poi = msh.ent2poi(tdim); 
        dblAr2 &coord   = msh.coord;

  double jtol = msh.param->jtol;
  const double qupdt0 = 1.0;
  const double qupdt1 = 1.0;

  double obj_change[gdim] = {0.0, 0.0};
  double obj_prev[gdim] = {0.0, 0.0};


  double ccoef[ncoef];
  double min_ccoef = getminccoef<gdim,ideg>(msh);

  printf("-- Enter bezGapsLP with jtol = %e min ccoef = %e \n",jtol,min_ccoef);

  int nelems = msh.nentt(tdim);
  int ncoefglob = nelems * ncoef;
  
  int noptim_points = pos_ctrlp.get_n();
  
  printf("\n NUMBER OF OPTIM = %d \n", noptim_points);

  // DoFs are ordered as dX followed by the D_i, followed by virtual
  // variable t in the case of min max.
  int ncol =  2 * noptim_points + !MAE;
  int nrow = (2 + !MAE) * noptim_points + ncoefglob;
  
  LPsolver solver(lib, method);
  solver.allocate(nrow, ncol);


  if(MAE){
    // The L1 norm does not depend on the X_i:
    for(int ii = 0; ii < noptim_points; ii++) solver.obj(ii) = 0;
    // but on the dD_i
    for(int ii = noptim_points; ii < 2*noptim_points; ii++) solver.obj(ii) = 1.0; 
  }else{
    // If min-max, the cost function depends only on the virtual variable 
    for(int ii = 0; ii < ncol; ii++) solver.obj(ii) = 0;
    solver.obj(ncol - 1) = 1; // min t
  }

  for(int ii = 0; ii < ncol; ii++) 
    solver.setVarCstr(ii, solver.INFINITY_m, solver.INFINITY_p);
  //for(int ii = 0; ii < noptim_points; ii++){
  //  solver.setVarCstr(ii, solver.INFINITY_m, solver.INFINITY_p);
  //}
  //for(int ii = noptim_points; ii < ncol; ii++){
  //  solver.setVarCstr(ii, 0.0, solver.INFINITY_p);
  //}

  // These (dbl/int)Ar classes entail dynamic alloc so they should be initialized
  // the fewest amount of times (not in loops, ideally)
  dblAr2 d_ccoef(ncoef, nnode);
  d_ccoef.set_n(ncoef);

  constexpr int vfact = ifact<tdim>();


  const int miter = 30;
  const double obj_change_tol = 1.0e-5;
  double dxavg[2][gdim] = {{0.0}, {0.0}};
  double dxmax[2][gdim] = {{-1.0e30}, {-1.0e30}};
  int iwhich = 0;
  for(int niter = 0; niter < miter; niter++){

    double t = niter / (double) miter;


    const double qupdt = t * qupdt1 + (1.0 - t)*qupdt0; // this changes the beha
    // const double qupdt = qupdt0;
    // Get max norm along x or y displacement 
    for(int ii = 0; ii < gdim; ii++) dxavg[iwhich][ii] = 0.0;
    for(int ii = 0; ii < gdim; ii++) dxmax[iwhich][ii] = -1.0e30;
    for(int icoor = 0; icoor < gdim; icoor++){

      msh.tag[0]++;

      double t0 = get_wall_time();

      for(int ielem = 0; ielem < nelems; ielem++){
        if (isdeadent(ielem, ent2poi)) continue;
        double vol = getmeasentP1<gdim>(ent2poi[ielem], coord); 

        if constexpr(tdim==2){
          d_ccoef_genbez2<ideg>(ent2poi, coord, ielem, ccoef, icoor, d_ccoef);
        }else if(tdim==3){
          d_ccoef_genbez3<ideg>(ent2poi, coord, ielem, icoor, ccoef, d_ccoef);
        }else METRIS_THROW_MSG(TODOExcept(), "derivatives not implemented for gdim = "<<gdim);

        for(int jj = 0; jj < nnode; jj++){
          
          int ipoin = ent2poi(ielem,jj);
          int irank = idx_point[ipoin];
          if(irank < 0) continue; // Boundary point (or point not optimized)

          int ivar_dX = irank;
          int ivar_dD = noptim_points + irank;
          int ivar_t  = 2*noptim_points; 
          int icst_1 = irank;
          int icst_2 = noptim_points + irank;
          int icst_4 = MAE ? -1 : 2*noptim_points + ncoefglob + irank;

          // Sps a method computes the desired positions for the control points, stored in pos_ctrlp
          // pos_ctrlp is size 1:nvar, 1:gdim
          // METRIS_ASSERT(ioptim <= pos_ctrlp.get_n());

          if(msh.poi2tag(0, ipoin) < msh.tag[0]){
            msh.poi2tag(0, ipoin) = msh.tag[0];

            double opt_ctrlpt = pos_ctrlp(irank,icoor);
            //double ctrlpt = 0;
            double ctrlpt = coord(ipoin,icoor);

            dxmax[iwhich][icoor] = MAX(dxmax[iwhich][icoor],abs(opt_ctrlpt - ctrlpt));
            dxavg[iwhich][icoor] += abs(opt_ctrlpt - ctrlpt);

            //printf("debug ipoin %d icst_1 = %d icst_2 = %d \n",ipoin,icst_1, icst_2);

            // INEQUALITY CONSTRAINTS 
            // dD_i - dX_i >=  Pi - Pi^opt - |Pi - Pi*|
            solver.setConstraintMatrix(icst_1, ivar_dD, 1);
            solver.setConstraintMatrix(icst_1, ivar_dX, -1);
            solver.setConstraintLB(icst_1, ctrlpt - opt_ctrlpt - abs(ctrlpt - opt_ctrlpt));
            solver.setConstraintUB(icst_1, solver.INFINITY_p);

            //printf("row %d put at j = %d %d \n",icst_2, ivar_dD, ivar_dX);
            //printf("cstrt >= %14.7e\n",opt_ctrlpt - ctrlpt);
            // INEQUALITY CONSTRAINTS
            // dD_i + dX_i >= Pi^opt - Pi - |Pi - Pi*|
            solver.setConstraintMatrix(icst_2, ivar_dD, 1);
            solver.setConstraintMatrix(icst_2, ivar_dX, 1);
            solver.setConstraintLB(icst_2, opt_ctrlpt - ctrlpt - abs(ctrlpt - opt_ctrlpt));
            solver.setConstraintUB(icst_2, solver.INFINITY_p);

            
            if (!MAE){
              // INEQUALITY CONSTRAINTS
              // t - dD_i >= |Pi - Pi*|
              // set t = |Pi - Pi*| + dt
              // min t + dt <=> min dt
              // cstr becomes
              // dt - dDi >= 0
              // solver.setConstraintMatrix(row_idx_glob, noptim_points+col_idx, 1);
              solver.setConstraintMatrix(icst_4, ivar_t , 1);
              solver.setConstraintMatrix(icst_4, ivar_dD, -1);
              solver.setConstraintLB(icst_4, 0);
              solver.setConstraintUB(icst_4, solver.INFINITY_p);  
            }
          }
          

          for(int ii = 0; ii < ncoef; ii++){
            // INEQUALITY CONSTRAINTS
            // dN dXi  >= jtol - alpha_i
            //  ^^^
            int icst_3 = 2*noptim_points + ncoef*ielem + ii;
            solver.setConstraintMatrix(icst_3, ivar_dX, d_ccoef(ii, jj));
          }
          
        } // for nnode

        for(int ii = 0; ii < ncoef; ii++){
            // dN dXi  >= jtol - alpha_i
            //               ^^^
            int icst_3 = 2*noptim_points + ncoef*ielem + ii;
            solver.setConstraintLB(icst_3, jtol*vol*vfact - ccoef[ii]); 
            solver.setConstraintUB(icst_3, solver.INFINITY_p);
          }
        
      }
      double t1 = get_wall_time();
    
      //for(int ii = 0; ii < nrow; ii++){
      //  printf("%d: ",ii);
      //  for(int jj = 0; jj < ncol; jj++){
      //    printf(" %5f ", sparseget(solver.Acstr_ALG,ii,jj));
      //  }
      //  printf("\n");
      //}


      // Solver can always end at 0, we need to manually compute the error. 
      double obj_val = solver.optimize();
      obj_change[icoor] = abs(obj_prev[icoor] - obj_val);
      obj_prev[icoor] = obj_val;

      //printf("Solution vector:\n");
      //for(int ii = 0; ii < ncol; ii++) {
      //  printf("%d : %15.7e\n",ii,solver.x_ALG[ii]);
      //}
      //printf("debug coord ip 6 %f %f \n", msh.coord(6,0), msh.coord(6,1));

      double t2 = get_wall_time();


      double min_ccoef_before = getminccoef<gdim,ideg>(msh);
      solver.updateCoord(msh, idx_point, icoor, qupdt);
      double min_ccoef_after = getminccoef<gdim,ideg>(msh);


      min_ccoef = min_ccoef_after;

      printf("   - %d/%d f = %15.7e t asbl %f solve %f ccoef %f -> %f\n", 
             niter,icoor,obj_val,t1-t0,t2-t1,min_ccoef_before/vol0,min_ccoef_after/vol0);


      if(iverb >= 3){
        std::string fname = "bezGapsLP"+std::to_string(niter)+"."+std::to_string(icoor);
        writeMesh(fname,msh);
      }
    } // for icoor 


    bool istop = true;
    for(int icoor = 0; icoor < gdim; icoor++){
      dxavg[iwhich][icoor] /= noptim_points;
      if(iverb >= 2) printf(" - niter %d/%d dx[%d] avg = %15.7e max = %15.7e "
                            "davg = %15.7e dmax = %15.7e\n",
                            niter,miter,icoor,dxavg[iwhich][icoor],dxmax[iwhich][icoor],
                            abs(dxavg[iwhich][icoor] - dxavg[1-iwhich][icoor]),
                            abs(dxmax[iwhich][icoor] - dxmax[1-iwhich][icoor]));
      if(!MAE && dxmax[iwhich][icoor] < obj_change_tol) continue;
      if( MAE && dxavg[iwhich][icoor] < obj_change_tol) continue;
      if(!MAE && abs(dxmax[iwhich][icoor] - dxmax[1-iwhich][icoor]) < obj_change_tol
              && niter >= 1) continue;
      if( MAE && abs(dxavg[iwhich][icoor] - dxavg[1-iwhich][icoor]) < obj_change_tol
              && niter >= 1) continue;
      istop = false;

      //if(obj_change[icoor] < obj_change_tol) continue;
      //istop = false;
    }
    if(istop){
      if(iverb >= 1) printf(" - changes too small in objective, stop \n");
      return min_ccoef;
    }
    iwhich = 1 - iwhich;
    if(istop){
      if(iverb >= 1) printf(" - changes too small, stop \n");
      return min_ccoef;
    }

  } // for niter 

  //printf("debug coord ip 6 %f %f \n", msh.coord(6,0), msh.coord(6,1));
  //METRIS_THROW(AlgoExcept());
  return -1;
}

#define BOOST_PP_LOCAL_MACRO(ideg)\
template double bezGapsLP<2,2,ideg>(MeshBase &msh, const intAr1 &idx_point, \
                              const dblAr2 &pos_ctrlp,LPMethod method,  LPLib lib);\
template double bezGapsLP<3,3,ideg>(MeshBase &msh, const intAr1 &idx_point, \
                              const dblAr2 &pos_ctrlp,LPMethod method,  LPLib lib);
#define BOOST_PP_LOCAL_LIMITS     (2, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


}// namespace Metris
