//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_interperr

#include "common_setup.hxx"

#include "../src/SolutionField/minInterpError.hxx"
#include "../src/SolutionField/interpError.hxx"
#include "../src/SolutionField/SolutionField.hxx"
#include "../src/utils/mprintf.hxx"
#include "../src/utils/CT_loop.hxx"
#include "../src/utils/aux_timer.hxx"
#include "../src/linalg/det.hxx"
#include "../src/utils/bernstein_prod.hxx"

#include <random>
#include <cmath>

namespace utf = boost::unit_test;

using namespace Metris;
typedef MetricFieldAnalytical MFT;



#if 0
// Simpler version that just samples a bunch of points and integrates that way
template<int idim, int ideg, int pdeg, int pnorm>
double interpErr_debug(const SolutionFieldAnalytical &sol, int ielem){
  static_assert(idim == 2); //change ordfac

  const MeshBase &msh = *(sol.msh);
  METRIS_ASSERT(idim == msh.get_tdim());
  METRIS_ASSERT(ideg == msh.curdeg);

  const intAr2& ent2poi = msh.ent2poi(idim);

  constexpr auto ordent = ORDELT(idim);

  // Compute integral exactly using Bernstein basis
  // The Jacobian determinant in the Bernstein basis is easy -> ccoef
  // The interpolation error is computed in Lagrange then lag2bez'd
  // The power, if needed, is obtained using square_bernstein
  // Then the product of that to the Jacobian is computed using prod_bernstein

  constexpr int ideg_jdet = idim*(ideg - 1); // Jacobian det degree
  constexpr int ideg_intp = pdeg;  // interpolate degree
  constexpr int ideg_max  = METRIS_MAX_DEG_ORDERING; // sampling nodes

  // Node counts. 
  constexpr int nnode_jdet = getnnode(idim,ideg_jdet);
  constexpr int nnode_intp = getnnode(idim,ideg_intp);
  constexpr int nnode_max  = getnnode(idim,ideg_max);


  // -- Compute the interpolate
  dblAr2 rfld_intp(nnode_intp,1);
  double bary[idim+1];
  for(int inode = 0; inode < nnode_intp; inode++){
    for(int ii = 0; ii < idim + 1; ii++) 
      bary[ii] = ordent[ideg_intp][inode][ii] / (double) ideg_intp;

    METRIS_ASSERT(idim == 2 && abs(bary[0] + bary[1] + bary[2] - 1) < 1.0e-12
               || idim == 3 && abs(bary[0] + bary[1] + bary[2] + bary[3] - 1) < 1.0e-12);

    rfld_intp(inode, 0) = sol.getSolBary(idim,ielem,bary);
  }


  // Evaluation/lag2bez indices (identity)
  intAr1 lfld(nnode_intp);
  for(int ii = 0; ii < nnode_intp; ii++) lfld[ii] = ii;
  constexpr auto eval_intp = idim == 1 ? eval1<idim,ideg_intp> 
                           : idim == 2 ? eval2<idim,ideg_intp> : eval3<idim,ideg_intp>;
  constexpr auto eval_elt  = idim == 1 ? eval1<idim,ideg> 
                           : idim == 2 ? eval2<idim,ideg> : eval3<idim,ideg>;
  double errLp = 0;
  printf("debug nnode_max %d \n",nnode_max);
  for(int inode = 0; inode < nnode_max; inode++){
    for(int ii = 0; ii < idim + 1; ii++) 
      bary[ii] = ordent[ideg_max][inode][ii] / (double) ideg_max;

    // f(F_K(\xi))
    double feval = sol.getSolBary(idim,ielem,bary);
    // Pi_K (f \circ F_K)
    double intrp; 
    eval_intp(rfld_intp, &lfld[0], FEBasis::Lagrange, DifVar::None, DifVar::None,
              bary, &intrp, NULL, NULL);

    double dum[idim];
    double jmat[idim*idim];
    eval_elt(msh.coord, ent2poi[ielem], msh.getBasis(), DifVar::Bary, DifVar::None,
             bary, dum, jmat, NULL);

    errLp += std::pow(abs(feval-intrp),pnorm)*detmat<idim>(jmat) / (double) nnode_max;
    //errLp += detmat<idim>(jmat) / (double) nnode_max;
  }
  errLp /= ifact<idim>();
  return errLp;
}
#endif


BOOST_AUTO_TEST_CASE(test_interperr) 
{


  std::vector<std::string> meshes = { METRIS_CASES_DIR "/unit/2D/misc/2tri2D.mesh"
                                     ,METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k"
                                     ,METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k -t 2"
                                    };

  constexpr int pnorm = 1;
  constexpr int pdeg  = 1;

  const double minsl = 0.5;

  constexpr int idim = 2;

  CT_FOR0_INC(1,2,pdeg){
    CT_FOR0_INC(1,2,pnorm){
      printf("-- Testing pnorm %d \n",pnorm);
      for(auto s : meshes)
      {
        printf("  - Reading mesh %s\n",s.c_str());
        cargHandler arg("-in " + s + " -anamet 1 -vdepth 0 -verb 0");
        MetrisRunner run(arg.c, arg.v);
        Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
        METRIS_ASSERT(idim == msh.idim);

        run.degElevate();
        msh.setBasis(FEBasis::Bezier);

        for(int isol = 1; isol <= 3; isol++){

          SolutionFieldAnalytical sol(msh);
          sol.setAnalyticalSolution(isol);

          CT_FOR0_INC(1,2,ideg){if(ideg == msh.curdeg){

            int nentt = msh.nentt(idim);
            const intAr2& ent2poi = msh.ent2poi(idim);

            double t0 = get_cpu_time();
            double errLp = 0;
            double errLp_dbg = 0;
            for(int ientt = 0; ientt < nentt; ientt++){
              if(isdeadent(ientt,ent2poi)) continue;
              double err = interpErr<idim,ideg,pdeg,pnorm>(sol, ientt);
              errLp += err;
              //printf("Debug ientt %d err %e \n",ientt, err);
            }
            double t1 = get_cpu_time();

            if(s != "../cases/2tri2D.mesh")
              printf("   - pnorm %d pdeg %d ideg %d isol %d got err %7.3e elt/s = %dk/s \n",
                      pnorm,pdeg,ideg,isol,errLp,(int)(nentt/(1000*(t1-t0))) );

            if(isol <= pdeg+1){
              BOOST_TEST(abs(errLp) <= 1.0e-15);
            }else if(s == "../cases/2tri2D.mesh" && pdeg == 1){
              // The mesh is linear so it is very simple to compute the error.
              // Function is x^2 + 2xy + y^2, domain is [0,1]^2
              // If pnorm = 1,
              // [1/3 + 1/2 + 1/3] = 7/6
              // pnorm = 2, 31/15
              double trueErrLp;
              if(pnorm == 1){
                // 1.0/24 in each element
                trueErrLp = 1.0/6;
              }else{
                trueErrLp = 1.0/30.0;            
              }
              BOOST_TEST(abs(trueErrLp - errLp) <= 1.0e-15 * trueErrLp);
              if(!(abs(trueErrLp - errLp) <= 1.0e-15 * trueErrLp)){
                printf("    - Debug trueErrLp %e got %e \n",trueErrLp, errLp);
                wait();
              }
            }

            if(isol > pdeg+1){
              //printf("-- Derivatives test isol %d \n",isol);
              // ---- Derivatives 
              constexpr int nhess = (idim*(idim+1))/2;
              const double dx0   = 1.0e-1;
              const double dx1   = 1.0e-10;
              const double qdx   = 10.0;
              const double tol   = 1.0e-9;
              int ndx = 0;
              for(double dx = dx0; dx > dx1; dx /= qdx) ndx++;
              dblAr1 d1err(ndx), d2err(ndx), logdx(ndx);

              double min_slop1 = 1.0e30;
              double max_slop1 = -1.0e30;
              double min_slop2 = 1.0e30;
              double max_slop2 = -1.0e30;
              for(int ientt = 0; ientt < nentt; ientt++){
                if(isdeadent(ientt,ent2poi)) continue;
                double epsent = getepsent<idim>(msh, idim, ientt);
                for(int inode = 0; inode < getnnode(idim, ideg); inode++){
                  // Those sized [idim+1][2] store at (+0, Ø), (+dx, -dx), ... 
                  double d1err_ana[idim+1][2][idim], d2err_ana[nhess];
                  double err[idim+1][2];
                  double d1err_fd[idim];
                  double d2err_fd[nhess];
                  err[0][1] = err[0][0] = interpErr<idim,ideg,pdeg,pnorm>(sol, ientt, inode, {d1err_ana[0][0], d2err_ana});
                  int ipoin = ent2poi(ientt,inode);
                  int idx = 0;


                  for(double dx_ = dx0; dx_ > dx1; dx_ /= qdx){
                    double dx = dx_*epsent;
                    for(int ii = 0; ii < idim; ii++){
                      double coor0 = msh.coord(ipoin,ii);
                      for(int isig = 0; isig <= 1; isig++){
                        // +- dx 
                        msh.coord(ipoin,ii) = coor0 + (2*isig-1)*dx;
                        err[ii+1][isig] = interpErr<idim,ideg,pdeg,pnorm>(sol, ientt, inode, 
                                       {d1err_ana[ii+1][isig]});
                        //if(isig == 1) 
                        //  printf("dbg dx %e d1err_ana[0+1][1][0] = %e \n",dx,d1err_ana[0+1][1][0]);
                        msh.coord(ipoin,ii) = coor0;
                      }
                      d1err_fd[ii] = (err[ii+1][1] - err[ii+1][0])/(2*dx);
                    }// for int ii 
                    for(int ii = 0; ii < idim; ii ++){
                      for(int jj = 0; jj < idim; jj++){
                        d2err_fd[sym2idx(ii,jj)] = (d1err_ana[ii+1][1][jj] - d1err_ana[ii+1][0][jj])/(2*dx);
                      }
                    }
                    d1err[idx] = 0.5*log(geterrl2<idim >(d1err_fd,d1err_ana[0][0]));
                    d2err[idx] = 0.5*log(geterrl2<nhess>(d2err_fd,d2err_ana   ));
                    logdx[idx] = log(dx);
                    //printf("Debug dx %15.7e d1ana %f %f d1fdf %f %f \n",
                    //      dx,d1err_ana[0][0][0],d1err_ana[0][0][1],d1err_fd[0],d1err_fd[1]);
                    //printf("Debug dx %15.7e d2ana %f %f %f d2fdf %f %f %f\n",
                    //      dx,d2err_ana[0],d2err_ana[1],d2err_ana[2],
                    //         d2err_fd[0],d2err_fd[1],d2err_fd[2]);
                    idx++;
                  }// for dx

                  double slopd1 = linearRegression(ndx,&logdx[0],&d1err[0]);
                  double slopd2 = linearRegression(ndx,&logdx[0],&d2err[0]);
                  double minerrd1 = 1.0e30;
                  double minerrd2 = 1.0e30;
                  for(int ii = 0; ii < ndx; ii++){
                    minerrd1 = MIN(minerrd1, exp(d1err[ii]));
                    minerrd2 = MIN(minerrd2, exp(d2err[ii]));
                  }
                  //BOOST_CHECK_MESSAGE(minsl < slopd1 || minerrd1 < tol,
                  //  " d1 slope "<<slopd1<<" under minimum "<<minsl<<" min = "<<minerrd1<<"\n");
                  //BOOST_CHECK_MESSAGE(minsl < slopd2 || minerrd2 < tol,
                  //  " d2 slope "<<slopd2<<" under minimum "<<minsl<<" min = "<<minerrd2<<"\n");
                  //if(!(minsl < slopd1 || minerrd1 < tol) || !(minsl < slopd2 || minerrd2 < tol)){
                  //  printf("## WAIT HERE pnorm %d pdeg %d ideg %d isol %d inode %d ielem %d\n",pnorm,pdeg,ideg,isol,inode,ientt);
                  //  wait();
                  //}
                  //if(minerrd1 >= tol) 
                    min_slop1 = MIN(min_slop1, slopd1);
                  //if(minerrd1 >= tol) 
                    max_slop1 = MAX(max_slop1, slopd1);
                  //if(minerrd2 >= tol) 
                    min_slop2 = MIN(min_slop2, slopd2);
                  //if(minerrd2 >= tol) 
                    max_slop2 = MAX(max_slop2, slopd2);
                }// for inode
              }// for ientt

              if(s != "../cases/2tri2D.mesh")
                printf("    -> %7.3f < d1 slope < %7.3f ; %7.3f < d2 slope < %7.3f\n",
                       min_slop1, max_slop1,min_slop2, max_slop2);
            }



            // Convergence study? 

          }}CT_FOR1(ideg);
        }
      }// for s meshes
    }CT_FOR1(pnorm);
  }CT_FOR1(pdeg);


}













