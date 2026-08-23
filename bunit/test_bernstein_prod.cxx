//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE bernstein_prod

#include "common_setup.hxx"

#include "utils/CT_loop.hxx"
#include "utils/bernstein_prod.hxx"

#include <random>
#include <cmath>

namespace utf = boost::unit_test;

using namespace Metris;
typedef MetricFieldAnalytical MFT;

BOOST_AUTO_TEST_CASE(bernstein_prod) 
{


  std::vector<std::string> meshes = {METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k"
                                     //,"../cases/curved_p2.meshb"
                                    };

  constexpr int pnorm = 2;
  constexpr int pdeg  = 1;

  const int pdeg_max = 3;

  printf("-- TEST: test square_bernstein\n");

  std::uniform_real_distribution<double> unif(-10.0,10.0);
  std::default_random_engine rng(0);
  // Maximum degree eval is instantiated for
  CT_FOR0_INC(1,3,idim){
    printf(" - Test dimension %d \n",idim);
    const int nsamb = 1;
    const int nsamp = 1; 
    // Generate bary wherever (not necessarily in elt)
    dblAr2 bary(nsamb,idim+1);
    for(int ii = 0; ii < nsamb; ii++){
      double sum = 0;
      do{
        for(int jj = 0; jj < idim+1; jj++){
          bary(ii,jj) = unif(rng);
          sum += bary(ii,jj);
        }
      }while(abs(sum) < 1.0e-16);

      for(int jj = 0; jj < idim+1; jj++){
        bary(ii,jj) /= sum;
      }
    }

    //CT_FOR0_INC(1,METRIS_MAX_DEG_JACOBIAN/2,ideg){
    CT_FOR0_INC(1,METRIS_MAX_DEG_JACOBIAN/2+1,ideg){
      printf("  - Test sq degree %d \n",ideg);
      int nnod1 = getnnode(idim,  ideg);
      int nnod2 = getnnode(idim,2*ideg);
      dblAr2 coef1(nnod1,1);
      dblAr2 coef2(nnod2,1);
      intAr1 lfld(nnod2);
      for(int ii = 0; ii < nnod2; ii++) lfld[ii] = ii;

      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnod1; ii++) coef1(ii,0) = unif(rng);
        square_bernstein<idim, idim, ideg>(coef1, coef2);

        dblAr2 coef2_2(nnod2,1); 
        prod_bernstein<idim, idim, ideg, ideg>(coef1, coef1, coef2_2);
        double err = 0;
        for(int ii = 0; ii < nnod2; ii++){
          err += abs(coef2(ii,0)-coef2_2(ii,0)) / METRIS_MAX(abs(coef2(ii,0)),1.0);
        }
        BOOST_REQUIRE(err < 1.0e-10);

        double f1, f2, fsq;
        if constexpr(ideg <= METRIS_MAX_DEG_JACOBIAN/2){
          constexpr auto eval_d = idim == 1 ? eval1<1,ideg> : 
                                  idim == 2 ? eval2<1,ideg> : eval3<1,ideg>;
          constexpr auto eval_2d = idim == 1 ? eval1<1,2*ideg> : 
                                   idim == 2 ? eval2<1,2*ideg> : eval3<1,2*ideg>;
          for(int isamb = 0; isamb < nsamb; isamb++){
            eval_d(coef1,&lfld[0], FEBasis::Bezier, DifVar::None, DifVar::None, bary[isamb], 
                   &f1, NULL, NULL);
            fsq = f1*f1; 
            eval_2d(coef2,&lfld[0], FEBasis::Bezier, DifVar::None, DifVar::None, bary[isamb], 
                   &f2, NULL, NULL);

            BOOST_TEST(abs(f2-fsq) < 1.0e-10*METRIS_MAX(f2,fsq));
          }
        }else{
          // Just check the squared values are the squares
          for(int ii = 0; ii < idim + 1; ii++){
            fsq = coef1(ii,0) * coef1(ii,0);
            f2 = coef2(ii,0);
            BOOST_TEST(abs(f2-fsq) < 1.0e-10*METRIS_MAX(f2,fsq));
            if(!(abs(f2-fsq) < 1.0e-10*METRIS_MAX(f2,fsq))){
              printf("Error ii %d f2 = %f fsq = %f \n",ii,f2,fsq);
              printf("coef1 = \n");
              coef1.print();
              printf("coef2 = \n");
              coef2.print();
              wait();
            }
          }
        }
      }
    }CT_FOR1(ideg);

    printf(" - Passed square test\n");


    //CT_FOR0_EXC(1,METRIS_MAX_DEG_JACOBIAN,ideg1){
    //  CT_FOR0_INC(1,METRIS_MAX_DEG_JACOBIAN - c_ideg1,ideg2){
    CT_FOR0_EXC(1,METRIS_MAX_DEG_JACOBIAN,ideg1){
      CT_FOR0_INC(1,METRIS_MAX_DEG_JACOBIAN - c_ideg1,ideg2){
        
        printf("   - Prod test ideg1 %d ideg2 %d \n",ideg1,ideg2);

        constexpr auto eval_d1 = idim == 1 ? eval1<1,ideg1> : 
                                 idim == 2 ? eval2<1,ideg1> : eval3<1,ideg1>;
        constexpr auto eval_d2 = idim == 1 ? eval1<1,ideg2> : 
                                 idim == 2 ? eval2<1,ideg2> : eval3<1,ideg2>;
        constexpr auto eval_ds = idim == 1 ? eval1<1,ideg1+ideg2> : 
                                 idim == 2 ? eval2<1,ideg1+ideg2> : eval3<1,ideg1+ideg2>;
        int nnod1 = getnnode(idim,ideg1);
        int nnod2 = getnnode(idim,ideg2);
        int nnods = getnnode(idim,ideg1+ideg2);
        dblAr2 coef1(nnod1,1);
        dblAr2 coef2(nnod2,1);
        dblAr2 coefs(nnods,1);
        intAr1 lfld(nnods);
        for(int ii = 0; ii < nnods; ii++) lfld[ii] = ii;

        for(int isamp = 0; isamp < nsamp; isamp++){
          for(int ii = 0; ii < nnod1; ii++) coef1(ii,0) = unif(rng);
          for(int ii = 0; ii < nnod2; ii++) coef2(ii,0) = unif(rng);
          prod_bernstein<idim, idim, ideg1, ideg2>(coef1, coef2, coefs);
          double f1, f2, fs;
          for(int isamb = 0; isamb < nsamb; isamb++){
            eval_d1(coef1,&lfld[0], FEBasis::Bezier, DifVar::None, DifVar::None, bary[isamb], 
                    &f1, NULL, NULL);
            eval_d2(coef2,&lfld[0], FEBasis::Bezier, DifVar::None, DifVar::None, bary[isamb], 
                    &f2, NULL, NULL);
            eval_ds(coefs,&lfld[0], FEBasis::Bezier, DifVar::None, DifVar::None, bary[isamb], 
                    &fs, NULL, NULL);

            BOOST_TEST(abs(fs-f1*f2) < 1.0e-10*METRIS_MAX(abs(fs),abs(f1*f2)));
            //printf("Debug ideg1 %d ideg2 %d \n",ideg1, ideg2);
            //printf("coef1: \n");coef1.print();
            //printf("coef2: \n");coef2.print();
            //printf("coefs: \n");coefs.print();
            //printf("Z21 comp %f should be %f idx %d \n",coefs(mul2nod(2,1),0),
            //  (2*coef1(0,0)*coef2(2,0) + coef1(1,0)*coef2(0,0))/3,mul2nod(2,1));

          }
        }
      }CT_FOR1(ideg2);
    }CT_FOR1(ideg1);
  }CT_FOR1(idim);


}