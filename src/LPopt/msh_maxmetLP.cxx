//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See $METRIS_ROOT/License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_maxccoef.hxx"
#include "LPsolver.hxx"

#include "../codegen_lag2bez.hxx"
#include "../codegen_ccoef_d.hxx"
#include "../low_geo.hxx"
#include "../io_libmeshb.hxx"
#include "../aux_timer.hxx"
#include "../aux_utils.hxx"
#include "../linalg/det.hxx"

#include "../Mesh/Mesh.hxx"
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

template<class MFT, int gdim, int tdim, int ideg>
double maximizeMetCcoef(Mesh<MFT> &msh, OptDoF idofs, LPMethod method, 
                        LPLib lib, const bool cstrCcoef){
  METRIS_ASSERT_MSG(ideg==2, "maximizeMetCcoef not implemented for ideg = " << ideg);
  // The P2 and Metric specialization
  const bool MAE = false;
  const bool metOn = true;
  constexpr int jdeg = tdim * (ideg - 1);
  constexpr int ncoef = tdim == 2 ? facnpps[jdeg] : tetnpps[jdeg];
  constexpr int nnode = tdim == 2 ? facnpps[ideg] : tetnpps[ideg];
  constexpr int nnod1 = tdim + 1;
  constexpr int vol0 = ifact<tdim>();
  const intAr2 &ent2poi = tdim == 2 ? msh.fac2poi : msh.tet2poi;
        dblAr2 &coord   = msh.coord;

  msh.met.setSpace(MetSpace::Log);
  msh.met.setBasis(FEBasis::Bezier);
  constexpr int nnmet = (gdim*(gdim+1))/2;
  double metl[nnmet];

  constexpr AsDeg asdmet = AsDeg::P1;

  // Initialize Barycentric coordinates
  auto ordent = ORDELT(tdim);
  double bary[ncoef][tdim+1];
  for(int ii = 0 ;ii < ncoef; ii++){
    for(int jj = 0; jj < tdim + 1; jj++){
      bary[ii][jj] = ordent[jdeg][ii][jj] / ((double) jdeg);
    }
  }

  double jtol = msh.param->jtol;
  // jtol = -1;
  if(jtol < 0) jtol = 1;
  const double qupdt0 = 0.2;
  const double qupdt1 = 0.8;
  
  double ccoef[ncoef];
  double min_ccoef = getminccoef<gdim,ideg>(msh);
  double obj_value_prev = 0;

  printf("-- Enter maximizeMetCcoef with jtol = %e min ccoef = %e \n",jtol,min_ccoef);

  int nelems = msh.nentt(tdim); 
  int ncoefglob = nelems * ncoef;

  intAr1 idx_point(msh.npoin);
  idx_point.set_n(msh.npoin);
  int noptim_points = 0;

  /* Selecting the optim variables */
  if (idofs == OptDoF::Full){
      // As in the unit test, whole mesh
     for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      idx_point[ipoin] = -1;
      if(msh.poi2bpo[ipoin] >= 0) continue;
      idx_point[ipoin] = noptim_points;
      noptim_points++; 
    }
  }else if (idofs == OptDoF::HO){
    // Only edge control points 
    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      idx_point[ipoin] = -1;
    }
    for(int ientt = 0; ientt < nelems; ientt++){
      // Ordering is always P1 nodes first, then high order (HO). 
      for(int ii = nnod1; ii < nnode; ii++){
        int ipoin = ent2poi(ientt, ii);
        if((idx_point[ipoin] >= 0)||(msh.poi2bpo[ipoin] >= 0)) continue;
        idx_point[ipoin] = noptim_points;
        noptim_points++;
      }
    }
  }
  printf("\n NUMBER OF OPTIM = %d \n", noptim_points);

  // Virtual variable t, and 1 variable (coordinate) per point and 2 identities
  int ncol;
  int nrow;

  // DoFs are ordered as X, followed by the Y_i+, the Y_i-, followed by virtual 
  // variable t in the case of min max. 
  ncol = 1*noptim_points + 2*ncoefglob + !MAE;
  // Constraints are ordered as N_i - M_i = Y_i+ - Y_i-
  // then, if minmax (!MAE), the t >= Y_i+ + Y_i-
  // then, if constrained ccoef, the N_i >= jtol
  nrow = ncoefglob*(1 + !MAE + cstrCcoef);
  
  LPsolver solver(lib, method);
  solver.allocate(nrow,ncol);
  
  if(MAE){
    // The L1 norm does not depend on the X_i:
    for(int ii = 0; ii < noptim_points; ii++) solver.obj(ii) = 0;
    // but on the y+ and y-: (could be weighted in the future)
    for(int ii = noptim_points; ii < ncol; ii++) solver.obj(ii) = 1.0/ncoefglob; 
  }else{
    // If min-max, the cost function depends only on the virtual variable 
    for(int ii = 0; ii < ncol; ii++) solver.obj(ii) = 0;
    solver.obj(ncol - 1) = 1; // min t
  }
  
  // Constraint bounds. All constraints are written as Au <= 0.
  for(int ii = 0; ii < nrow; ii++) solver.setConstraintUB(ii, 0);
  solver.setConstraintLBinf();

  // These (dbl/int)Ar classes entail dynamic alloc so they should be initialized
  // the fewest amount of times (not in loops, ideally)
  dblAr2 d_ccoef(ncoef, nnode);
  d_ccoef.set_n(ncoef);

  double detm, detJ0;

  // For lag2bez 
  int lfld[ncoef];
  for(int ii = 0; ii < ncoef; ii++) lfld[ii] = ii;
  dblAr2 rfld0(ncoef, 1);
  dblAr2 rfld1(ncoef, 1);
  rfld0.set_n(ncoef);
  rfld1.set_n(ncoef);
  
  const int miter = 100;
  const double obj_change_tol = 1.0e-5;
  for(int niter = 0; niter < miter; niter++){

    double t = niter / (double) miter;

    const double qupdt = t * qupdt1 + (1.0 - t)*qupdt0;
    //const double qupdt = qupdt0;
    // Get max norm along x or y displacement 
    for(int icoor = 0; icoor < tdim ; icoor++){

      double t0 = get_wall_time();

      for(int ielem = 0; ielem < nelems; ielem++){
        if (isdeadent(ielem, ent2poi)) continue;

        // full case we need to scale the control coeffs and their derivatives
        double vol = getmeasentP1<gdim>(ent2poi[ielem], coord);

        int row = ielem;
        int idx_start_row = row * ncoef;

        // Gather the metric determinants at the Lagrange nodes into rfld0
        for(int jj = 0; jj < ncoef; jj++){
          msh.met.getMetBary(asdmet, DifVar::None, 
                             MetSpace::Exp, ent2poi[ielem], 
                             tdim, bary[jj], metl, NULL);
          detm = detsym<tdim, double>(metl);
          // For s.p.d. det(M^-1/2) = det(M)^-1/2 
          rfld0(jj,0) = pow(detm, -0.5);
        }
        // Bezier coeffs of determinants at nodes (from lagrange coeffs)
        // Right now we are computing lag2bez to get the bezier coeffs 
        // at the nodes (from the lagrange coeffs)
        lag2bez2<ideg,1>(lfld, rfld0, rfld1);
        

        if(tdim==2){
          d_ccoef_genbez2<ideg>(ent2poi, coord, ielem, ccoef, icoor, d_ccoef);
        }else if(tdim==3){
          d_ccoef_genbez3<ideg>(ent2poi, coord, ielem, icoor, ccoef, d_ccoef);
        }else METRIS_THROW_MSG(TODOExcept(), 
                               "derivatives not implemented for tdim = "<<tdim);
  
        for(int ii = 0; ii < ncoef; ii++){
            int row_idx = idx_start_row + ii;

            // If max-min, set the constraints t >= y+ + y- ? (why 0 for t ?)
            if (!MAE){
              solver.setConstraintMatrix(row_idx, ncol-1, 0);
              solver.setConstraintMatrix(row_idx + ncoefglob, ncol-1, -1);
            }
            
            for(int jj = 0; jj < nnode; jj++){
                int ipoin = ent2poi(ielem,jj);
                int irank = idx_point[ipoin];
                if(irank < 0) continue; // Boundary point (or point not optimized)
                int col_idx = irank;
                
                // Set X part of the constraint N_i - M_i = Y_i+ - Y_i-
                solver.setConstraintMatrix(row_idx, col_idx, d_ccoef(ii, jj));

                // If ccoef >= jtol enforced, set X part of N_i/vol >= jtol
                if (cstrCcoef){
                  if(MAE) solver.setConstraintMatrix(row_idx + ncoefglob, col_idx, 
                                                     d_ccoef(ii,jj)/vol);
                  else    solver.setConstraintMatrix(row_idx +2*ncoefglob, col_idx, 
                                                     d_ccoef(ii, jj)/vol);
                }
                
            }
            detJ0 = tdim == 2 ? sqrt(3)/2 : sqrt(2)/2;

            // rdlf0 stores the Lagrange Coeffs (determinants evaluated at the nodes)
            // rdlf1 stores the corresponding Bezier Coeffs (control coeffs)
            // Set LB/UB for constraints N_i - M_i = Y_i+ - Y_i-
            solver.setConstraintLB(row_idx, rfld1(ii, 0) * detJ0 - ccoef[ii]);
            solver.setConstraintUB(row_idx, rfld1(ii, 0) * detJ0 - ccoef[ii]);


            if (cstrCcoef){
              double row_index;
              if (MAE){  
                row_index = row_idx+ncoefglob;
              }else{
                row_index = row_idx + 2*ncoefglob;
              }
              // METRIS_ASSERT(rfld1(ii, 0) * detJ0 < 0);
              solver.setConstraintLB(row_index, -ccoef[ii]/vol + jtol);
              // if (rfld1(ii, 0) * detJ0 < 0) solver.setConstraintLB(row_index, -INFINITY);
              // solver.setConstraintLB(row_index, -INFINITY);
              solver.setConstraintUB(row_index, INFINITY);
            }
        }
      } // for ielem

      double t1 = get_wall_time();
      // I blocks
      for(int k = 0; k < ncoefglob; k++){
        if (MAE){
          solver.setConstraintMatrix(k, noptim_points+k, -1);
          solver.setConstraintMatrix(k, noptim_points+k+ncoefglob, 1);
        }else{
          solver.setConstraintMatrix(k, noptim_points+k, -1);
          solver.setConstraintMatrix(k, noptim_points+k+ncoefglob, 1);
          solver.setConstraintMatrix(k+ncoefglob, noptim_points+k, 1);
          solver.setConstraintMatrix(k+ncoefglob, noptim_points+k+ncoefglob, 1);
        }
      }

      if (metOn){
        minlpsetbci(solver.state, ncol-1, 0, INFINITY);
        for(int i = 0; i<noptim_points; i++) minlpsetbci(solver.state, i, -INFINITY, INFINITY); 
        for(int i = noptim_points; i<noptim_points+2*ncoefglob; i++) minlpsetbci(solver.state, i, 0, INFINITY);
      }else minlpsetbcall(solver.state, -INFINITY, INFINITY);

      double obj_val = solver.optimize();
      double obj_change = abs(obj_value_prev - obj_val);
      printf(" - f = %f prev = %f df = %15.7e\n", obj_val, obj_value_prev,
              obj_change);
      obj_value_prev = obj_val;

      double t2 = get_wall_time();

      printf(" - CPU time assembly %f solve %f \n",t1-t0,t2-t1);

      double min_ccoef_before = getminccoef<gdim,ideg>(msh);

      // Metric Reinterpolation
      for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
        int irnk1 = idx_point[ipoin];
        if(irnk1<0) continue;
        int iseed;
        if(gdim == 2){
          iseed = getpoifac(msh, ipoin);
        }else{
          iseed = getpoitet(msh, ipoin);
        }
        msh.interpMetBack(ipoin, msh.idim, iseed, -1, NULL);
      }

      solver.updateCoord(msh, idx_point, icoor, qupdt);

      double min_ccoef_after = getminccoef<gdim,ideg>(msh);
      min_ccoef = min_ccoef_after;

      printf(" %d - %d : ccoef %f -> %f\n",niter,icoor,min_ccoef_before/vol0,min_ccoef_after/vol0);
      // mxnrm = MAX(mxnrm,nrm);
      // mxdNm = MAX(mxdNm,min_ccoef_after - min_ccoef_before);
      
      if ((obj_change < obj_change_tol)){
        printf(" - gone over threshold\n");
        
        // reset update 
        // double t0s = get_wall_time();

        // double q0 = 0, q1 = 1;
        // double qp = qupdt; // q previous
        // double tol0 = 0.1;
        // int miters = 100;
        // bool iok = false;
        // double min_final;
        // for(int niters = 0; niters < miters; niters++){
        //   // Cancel previous 
        //   solver.updateCoord(msh, idx_point, icoor, -qp);


        //   qp = (q0 + q1)/2;
        //   // useful when qp -> 1
        //   // if(1-qp < tol0) qp = 1-qp;
        //   solver.updateCoord(msh, idx_point, icoor, qp);
          
        //   min_ccoef = getminccoef<gdim,ideg>(msh);

        //   printf(" - bisection iter %d / %d q = %e min = %e \n",niters,miters,qp,min_ccoef);

        //   if(min_ccoef >= jtol && min_ccoef < jtol*(1.0 + tol0)){
        //     iok = true;
        //     min_final = min_ccoef;
        //     break;
        //   }

        //   if(min_ccoef > jtol){
        //     q1 = qp;
        //   }else{
        //     q0 = qp;
        //   }

        // } // for niters
        // double t1s = get_wall_time();
        // printf("Bisection time %f \n",t1s-t0s);

        // if(!iok) METRIS_THROW(AlgoExcept());

        // // Got a good correction
        // double t1_tot = get_wall_time();
        // printf(" - Backtract exit min ccoef = %e > %e = jtol total time = %f \n",
        //                                           min_final,jtol,t1_tot-t0_tot);

        // return min_final;
        return min_ccoef;

      }
      msh.setBasis(FEBasis::Lagrange); // Vizir assumes Lagrange
      std::string fname;
      fname = "MetLP"+std::to_string(niter)+"."+std::to_string(icoor);
      writeMesh(fname,msh);
      msh.setBasis(FEBasis::Bezier);

    } // for icoor 


    // if(mxdNm < dNtol){
    //   printf("-> min ccoef changes too small = %e, exit\n",mxdNm);
    //   return min_ccoef;
    // }
  } // for niter 

  METRIS_THROW(AlgoExcept());
  return -1;
}


#define BOOST_PP_LOCAL_MACRO(ideg) \
template double maximizeMetCcoef<MetricFieldFE        ,2,2,ideg>\
  (Mesh<MetricFieldFE> &msh, OptDoF idofs, LPMethod method, LPLib lib, const bool cstrCcoef);\
template double maximizeMetCcoef<MetricFieldFE        ,3,3,ideg>\
  (Mesh<MetricFieldFE> &msh, OptDoF idofs, LPMethod method, LPLib lib, const bool cstrCcoef);\
template double maximizeMetCcoef<MetricFieldAnalytical,2,2,ideg>\
  (Mesh<MetricFieldAnalytical> &msh, OptDoF idofs, LPMethod method, LPLib lib, const bool cstrCcoef);\
template double maximizeMetCcoef<MetricFieldAnalytical,3,3,ideg>\
  (Mesh<MetricFieldAnalytical> &msh, OptDoF idofs, LPMethod method, LPLib lib, const bool cstrCcoef);
#define BOOST_PP_LOCAL_LIMITS (2, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()


}// namespace Metris
