//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_eval_bezierfunc 

#include <boost/test/included/unit_test.hpp> 

#include "../src/utils/CT_loop.hxx"
#include "../src/low_eval.hxx"

#include "common_setup.hxx"

#include <random>



namespace Metris{


BOOST_AUTO_TEST_CASE(test_eval_bezierfunc) 
{

  const int nsamp = 100;

  std::uniform_real_distribution<double> unif(-10.0,10.0);
  std::default_random_engine rng(0);

  const double tol_zero = 1.0e-16;
  const double tol  = 1.0e-12;


  const double dx0   = 1.0e+3;
  const double dx1   = 1.0e-6;
  const double qdx   = 3.0;
  const double minsl = 0.5;
  int ndx = 0;
  for(double dx = dx0; dx > dx1; dx /= qdx) ndx++;
  

  CT_FOR0_INC(1,3,tdim){
    constexpr auto ordent = ORDELT(tdim);
    // Generate samples
    dblAr2 bary(nsamp,tdim+1);
    for(int isamp = 0; isamp < nsamp; isamp++){
      double sum = 0;
      do{
        for(int jj = 0; jj < tdim + 1; jj++){
          bary(isamp,jj) = unif(rng);
          sum += bary(isamp,jj);
        }
      }while(abs(sum) < 1.0e-16);

      for(int jj = 0; jj < tdim + 1; jj++) bary(isamp,jj) /= sum;
    }// for isamp

    double dfun[tdim];

    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ifun = 0; ifun < getnnode(tdim,ideg); ifun++){
          double feval = eval_bezierfunc<ideg, tdim>(ordent[ideg][ifun], bary[isamp], 1, dfun);
          double fbez = 1;
          for(int idim = 0; idim < tdim+1; idim++){
            fbez *= pow(bary[isamp][idim], ordent[ideg][ifun][idim]);
          }
          if constexpr(tdim == 1){
            fbez *= cbzedg.s[ideg][ifun];
          }else if(tdim == 2){
            fbez *= cbzfac.s[ideg][ifun];
          }else{
            fbez *= cbztet.s[ideg][ifun];
          }

          BOOST_REQUIRE(abs(fbez - feval) < tol * abs(feval)
             || abs(feval) < tol_zero && abs(fbez) < tol_zero);

          // Test derivative
          for(int ivar = 0; ivar < tdim; ivar++){
            dblAr1 errdiff(ndx), logdx(ndx);
            double minerr = 1.0e30;
            int idx = 0;
            for(double dx = dx0; dx > dx1; dx /= qdx, idx++){
              double bar1[tdim+1];
              for(int ii = 0; ii < tdim + 1; ii++) bar1[ii] = bary[isamp][ii];
              bar1[0]      -= dx;
              bar1[1+ivar] += dx;
              double feva1 = eval_bezierfunc<ideg, tdim>(ordent[ideg][ifun], bar1, 0, NULL);
              double diff_disc = (feva1 - feval) / dx;
              double err0  = abs(diff_disc - dfun[ivar]);
              minerr = MIN(minerr,err0);
              //printf("debug dx = %e err %e fun ref %e new %e dref %e ddf %e\n",dx,err0,feval,feva1,
              //  dfun[ivar],diff_disc);
              errdiff[idx] = err0;
              logdx[idx] = log(dx);
            }// for dx

            //printf("debug idim %d ideg %d ifun %d ivar %d minerr %e\n",tdim, ideg, ifun, ivar, minerr);
            // If not smallest less than tolerance, compute slope
            if(minerr >= tol){
              for(int idx = 0; idx < ndx; idx++) errdiff[idx] = log(errdiff[idx]);
              double slope = linearRegression(ndx,&logdx[0],&errdiff[0]);
              BOOST_CHECK_MESSAGE(slope >= minsl, 
                "Diff slope ivar "<<ivar<<" ifun "<<ifun<<" is small "<<slope<<" minerr = "<<minerr<<
                "\n ideg = "<<ideg<<" ifun = "<<ifun<<" ivar = "<<ivar);
              //BOOST_REQUIRE_MESSAGE(slope >= minsl, 
              //  "Diff slope ivar "<<ivar<<" ifun "<<ifun<<" is small "<<slope<<" minerr = "<<minerr<<
              //  "\n ideg = "<<ideg<<" ifun = "<<ifun<<" ivar = "<<ivar);
            }

          }// for ivar

        }// for ifun
      }// for isamp
    }CT_FOR1(ideg);
  }CT_FOR1(tdim);

}// end test case


} // end namespace