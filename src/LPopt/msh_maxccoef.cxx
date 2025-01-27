//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_maxccoef.hxx"
#include "LPsolver.hxx"

#include "../low_ccoef.hxx"
#include "../low_geo.hxx"
#include "../io_libmeshb.hxx"
#include "../aux_utils.hxx"
#include "../aux_timer.hxx"

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
 
template<int gdim, int tdim, int ideg>
double maximizeCcoef(MeshBase &msh, OptDoF idofs, LPMethod method, LPLib lib){
  // METRIS_THROW_MSG(TODOExcept(), "maximizeCcoef not implemented for ideg = "<<ideg);
  msh.setBasis(FEBasis::Bezier); // Vizir assumes Lagrange
  constexpr int jdeg = tdim * (ideg - 1);
  constexpr int ncoef = tdim == 2 ? facnpps[jdeg] : tetnpps[jdeg];
  // Full node count 
  constexpr int nnode = tdim == 2 ? facnpps[ideg] : tetnpps[ideg];
  // Only P1 nodes 
  constexpr int nnod1 = tdim == 2 ? 3 : 4;
  constexpr int vol0 = ifact<tdim>();
  bool use_vol_scale = false;

  const intAr2 &ent2poi = tdim == 2 ? msh.fac2poi : msh.tet2poi;
        dblAr2 &coord   = msh.coord;

  if(msh.param->iverb >= 2){
    FEBasis ibas0 = msh.getBasis();
    msh.setBasis(FEBasis::Lagrange);
    writeMesh("maximizeCcoef.inp", msh);
    msh.setBasis(ibas0);
  }

  
  double jtol = msh.param->jtol;
  if(jtol < 0) jtol = 1;
  const double qupdt0 = 0.2;
  const double qupdt1 = 0.8;
  
  double t0_tot = get_wall_time();

  double ccoef[ncoef];
  double min_ccoef = getminccoef<gdim,ideg>(msh);

  printf("-- Enter maximizeCcoef with jtol = %e min ccoef = %e \n",
         jtol,min_ccoef);
  if(min_ccoef >= jtol){
    printf(" - Nothing to do: exit\n");
    return min_ccoef;
  }


  int nelems = tdim == 2 ? msh.nface : msh.nelem;
  int nrow = nelems * ncoef;
  
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
    use_vol_scale = true;
    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
      idx_point[ipoin] = -1;
    }
    for(int ielem = 0; ielem < nelems; ielem++){
      // Ordering is always P1 nodes first, then high order (HO). 
      for(int ii = nnod1; ii < nnode; ii++){
        int ipoin = ent2poi(ielem, ii);
        if((idx_point[ipoin] >= 0)||(msh.poi2bpo[ipoin] >= 0)) continue;
        idx_point[ipoin] = noptim_points;
        noptim_points++;
      }
    }
  }
  printf("\n NUMBER OF OPTIM = %d \n", noptim_points);

  // Virtual variable t, and 1 variable (coordinate) per point 
  int ncol = 1 * noptim_points + 1;


  LPsolver solver(lib, method);
  solver.allocate(nrow, ncol);

  
  for(int ii = 0; ii < ncol; ii++) solver.obj(ii) = 0;
  solver.obj(ncol - 1) = -1;

  for(int jj = 0; jj < nrow; jj++) solver.setConstraintLB(jj, 0.0);
  solver.setConstraintUBinf(); // infinity upper bound 

  // These (dbl/int)Ar classes entail dynamic alloc so they should be initialized
  // the fewest amount of times (not in loops, ideally)
  dblAr2 d_ccoef(ncoef, nnode);
  d_ccoef.set_n(ncoef);

  const int miter = 30;
  const double dNtol = 1.0e-3;
  for(int niter = 0; niter < miter; niter++){

    double t = niter / (double) miter;

    const double qupdt = t * qupdt1 + (1.0 - t)*qupdt0; // this changes the beha
    // const double qupdt = qupdt0;
    // Get max norm along x or y displacement 
    double mxdNm = -1;
    for(int icoor = 0; icoor < tdim; icoor++){

      double t0 = get_wall_time();

      for(int ielem = 0; ielem < nelems; ielem++){
        
        if (isdeadent(ielem, ent2poi)) continue;
        // full case we need to scale the control coeffs and their derivatives
        double vol = use_vol_scale ? getmeasentP1<gdim>(ent2poi[ielem], coord) : 1.0;

        int row = ielem;
        int idx_start_row = row * ncoef;

        getccoef_dcoord<tdim,ideg>(msh,ielem,icoor,ccoef,d_ccoef); 

        for(int ii = 0; ii < ncoef; ii++){
          int row_idx = idx_start_row + ii;
          
          solver.setConstraintMatrix(row_idx, ncol-1, -1);
          
          for(int jj = 0; jj < nnode; jj++){
            int ipoin = ent2poi(ielem,jj);
            int irank = idx_point[ipoin];
            if(irank < 0) continue; // Boundary point (or point not optimized)
            int col_idx = irank;

            solver.setConstraintMatrix(row_idx, col_idx, d_ccoef(ii, jj)/vol);

          }
          solver.setConstraintLB(row_idx,-ccoef[ii]/vol);
        }
        row++;
      }
      double t1 = get_wall_time();

      minlpsetbcall(solver.state, -INFINITY, INFINITY);

      solver.optimize();

      double t2 = get_wall_time();

      printf(" - CPU time assembly %f solve %f \n",t1-t0,t2-t1);

      double min_ccoef_before = getminccoef<gdim,ideg>(msh);
      solver.updateCoord(msh, idx_point, icoor, qupdt);
      double min_ccoef_after = getminccoef<gdim,ideg>(msh);

      min_ccoef = min_ccoef_after;

      printf(" %d - %d : ccoef %f -> %f\n",niter,icoor,min_ccoef_before/vol0,min_ccoef_after/vol0);
      mxdNm = MAX(mxdNm,abs(min_ccoef_after - min_ccoef_before));

      if(min_ccoef_after > jtol && min_ccoef_before < jtol){
        printf(" - gone over threshold\n");
        
        // reset update 
        double t0s = get_wall_time();
        double q0 = 0, q1 = 1;
        double qp = qupdt; // q previous
        double tol0 = 0.1;
        int miters = 100;
        bool iok = false;
        double min_final;
        for(int niters = 0; niters < miters; niters++){
          // Cancel previous 
          solver.updateCoord(msh, idx_point, icoor, -qp);

          qp = (q0 + q1)/2;

          solver.updateCoord(msh, idx_point, icoor, qp);

          min_ccoef = getminccoef<gdim,ideg>(msh);
          printf(" - bisection iter %d / %d q = %e min = %e \n",niters,miters,qp,min_ccoef);

          if(min_ccoef >= jtol && min_ccoef < jtol*(1.0 + tol0)){
            iok = true;
            min_final = min_ccoef;
            break;
          }

          if(min_ccoef > jtol){
            q1 = qp;
          }else{
            q0 = qp;
          }

        } // for niters
        double t1s = get_wall_time();
        printf("Bisection time %f \n",t1s-t0s);

        if(!iok) METRIS_THROW(AlgoExcept());

        // Got a good correction
        double t1_tot = get_wall_time();
        printf(" - Backtract exit min ccoef = %e > %e = jtol total time = %f \n",
                                                  min_final,jtol,t1_tot-t0_tot);
        return min_final;
      }
    } // for icoor 

    if(mxdNm < dNtol){
      printf("-> min ccoef changes too small = %e, exit\n",mxdNm);
      return min_ccoef;
    }
  } // for niter 

  METRIS_THROW(AlgoExcept());
  return -1;
}

#define BOOST_PP_LOCAL_MACRO(ideg)\
template double maximizeCcoef<2,2,ideg>(MeshBase &msh, OptDoF idofs, LPMethod method, LPLib lib); \
template double maximizeCcoef<3,3,ideg>(MeshBase &msh, OptDoF idofs, LPMethod method, LPLib lib);

#define BOOST_PP_LOCAL_LIMITS     (2, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

template<int gdim, int ideg>
double getminccoef(MeshBase &msh){
  constexpr int tdim = gdim;
  constexpr int jdeg = tdim * (ideg - 1);
  constexpr int ncoef = tdim == 2 ? facnpps[jdeg] : tetnpps[jdeg];
  int nelems = gdim == 2 ? msh.nface : msh.nelem; // TODO : tdim instead of gim
  const intAr2 &ent2poi = gdim == 2 ? msh.fac2poi : msh.tet2poi;
  double iret = 1.0e30;
  constexpr int vol0 = ifact<gdim>();
  double ccoef[ncoef];
  for(int ielem = 0; ielem < nelems; ielem++){
    if(isdeadent(ielem, ent2poi)) continue;
    double vol = getmeasentP1<gdim>(ent2poi[ielem], msh.coord);
    //if(vol < Defaults::vtol){
    //  printf(" NEGATIVE VOLUME ? vol = %23.15e vtol = %15.7e",vol,msh.param->vtol);

    //  bool iflat;
    //  vol = getmeasentP1<gdim,tdim>(ent2poi[ielem],msh.coord,msh.param->vtol, 
    //                NULL, &iflat);
    //  printf(" recomp vol = %23.15e iflat %d",vol,iflat);
    //}
    bool iflat;
    if(!(vol > 0.0)){
      printf("## FATAL vol = %15.7e \n", vol);
      vol = getmeasentP1<gdim,gdim>(msh, ent2poi[ielem], NULL, &iflat);
      printf("Recompute with tol %15.7e iflat %d \n",vol, iflat);
      printf("element is %d = ", ielem);
      intAr1(tdim + 1, ent2poi[ielem]).print();
      writeMesh("debug_vol", msh);
    }
    METRIS_ENFORCE_MSG(vol>0.0,"vol = "<<vol);
    getsclccoef<gdim,gdim,ideg>(msh,ielem,NULL,ccoef,&iflat);
    for(int ii = 0; ii < ncoef; ii++){
      iret = MIN(iret,ccoef[ii]/vol/vol0);
    }
  }
  return iret;
}

#define BOOST_PP_LOCAL_MACRO(ideg)\
template double getminccoef<2,ideg>(MeshBase &msh); \
template double getminccoef<3,ideg>(MeshBase &msh);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()



}// end namespace
