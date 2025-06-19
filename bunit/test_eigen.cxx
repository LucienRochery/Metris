//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 

#include <boost/test/included/unit_test.hpp> 
#include "common_setup.hxx"

#include <boost/timer/progress_display.hpp>
#include <random>

//#include "../src/utils/aux_utils.hxx"
#include "../src/low_geo.hxx"
#include "../src/linalg/det.hxx"
#include "../src/linalg/eigen.hxx"
#include "../src/linalg/explogmet.hxx"


#include <Eigen/Dense>


namespace Metris{


// -- Test metqua3_xi_d 
// Constant metric fields should yield reliable derivatives in all cases 
// In non constant metric fields, derivatives only defined for DoFs in back element
// interiors... 

template<int ndim>
void generate_metric(double aniso, double eigval[ndim], double eigvec[ndim][ndim],
  std::uniform_real_distribution<double>& unif, std::default_random_engine& rng);

BOOST_AUTO_TEST_CASE(test_eval3) 
{//METRIS_MAX_DEG

  const int nsamp = 1e5;
  const double tol = 1.0e-11;

  std::uniform_real_distribution<double> unif(1.0e-16,1.0);
  std::default_random_engine rng(0);

  double aniso_max = 2e7;
  double aniso_mul = 10;

  CT_FOR0_INC(2,3,ndim){
    constexpr int nnmet = (ndim*(ndim+1))/2;

    double eigval[ndim],eigvec[ndim][ndim];
    double eigva2[ndim],eigve2[ndim][ndim];
    double met[nnmet], met2[nnmet];

    printf("-- Tests for dim = %d \n",ndim);


    for(double anisorat = 2; anisorat <= aniso_max + 1; anisorat *= aniso_mul){

      double errminmet_aniso = 1.0e30;
      double errmaxmet_aniso = -1;
      double erravgmet_aniso = 0;

      double errminmetL_aniso = 1.0e30;
      double errmaxmetL_aniso = -1;
      double erravgmetL_aniso = 0;

      double errminmetE_aniso = 1.0e30;
      double errmaxmetE_aniso = -1;
      double erravgmetE_aniso = 0;

      double errmineig_aniso = 1.0e30;
      double errmaxeig_aniso = -1;
      double erravgeig_aniso = 0;

      double errmineigL_aniso = 1.0e30;
      double errmaxeigL_aniso = -1;
      double erravgeigL_aniso = 0;

      double errmineigE_aniso = 1.0e30;
      double errmaxeigE_aniso = -1;
      double erravgeigE_aniso = 0;

      double nerro_aniso = 0;

      dblAr2 eigval_samples(nsamp, ndim);
      dblAr2 eigvec_samples(nsamp, ndim*ndim);
      eigval_samples.set_n(nsamp);
      eigvec_samples.set_n(nsamp);

      dblAr2 met_samples(nsamp, nnmet);

      for(int isamp = 0; isamp < nsamp; isamp++){
        nerro_aniso++;

        generate_metric<ndim>(anisorat, eigval, eigvec, unif, rng);
        for(int ii = 0; ii < ndim; ii++){
          eigval_samples(isamp, ii) = eigval[ii];
          for(int jj = 0; jj < ndim; jj++){
            eigvec_samples(isamp, ndim*ii+jj) = eigvec[ii][jj];
          }
        }

        // Generate metric from eigvecs, eigvals
        eig2met<ndim,double>(eigval,eigvec[0],met);
        for(int ii = 0; ii < nnmet; ii++) met_samples(isamp, ii) = met[ii];


        // ----- Test geteigsym (SurrealS compatible dsyevq)

        // Re-decompose and compare
        geteigsym<ndim,double>(met,eigva2,eigve2[0]);

        // Test unit norm 
        for(int i = 0; i < ndim; i++){
          BOOST_TEST( abs(sqrt(getnrml2<ndim>(eigve2[0])) - 1) < tol);
        }
        // Test orthogonality 
        BOOST_TEST( abs(getprdl2<ndim>(eigve2[0],eigve2[1])) < tol);
        if constexpr (ndim > 2){
          BOOST_TEST( abs(getprdl2<ndim>(eigve2[0],eigve2[2])) < tol);
          BOOST_TEST( abs(getprdl2<ndim>(eigve2[1],eigve2[2])) < tol);
        }

        eig2met<ndim,double>(eigva2,eigve2[0],met2);

        double errmet = sqrt(geterrl2<nnmet>(met2,met) / getnrml2<nnmet>(met));
        errminmet_aniso = MIN(errminmet_aniso, errmet);
        errmaxmet_aniso = MAX(errmaxmet_aniso, errmet);
        erravgmet_aniso += errmet;

        double erreig = sqrt(geterrl2<ndim>(eigval,eigva2) / getnrml2<ndim>(eigval));
        errmineig_aniso = MIN(errmineig_aniso, erreig);
        errmaxeig_aniso = MAX(errmaxeig_aniso, erreig);
        erravgeig_aniso += erreig;



        // ----- Test LAPACK
        double rwork[10];
        geteigsym_LAPACK<ndim>(met,10,rwork,eigva2,eigve2[0]);

        // Test unit norm 
        for(int i = 0; i < ndim; i++){
          BOOST_TEST( abs(sqrt(getnrml2<ndim>(eigve2[0])) - 1) < tol);
        }
        // Test orthogonality 
        BOOST_TEST( abs(getprdl2<ndim>(eigve2[0],eigve2[1])) < tol);
        if constexpr (ndim > 2){
          BOOST_TEST( abs(getprdl2<ndim>(eigve2[0],eigve2[2])) < tol);
          BOOST_TEST( abs(getprdl2<ndim>(eigve2[1],eigve2[2])) < tol);
        }

        eig2met<ndim,double>(eigva2,eigve2[0],met2);

        double errmetL = sqrt(geterrl2<nnmet>(met2,met) / getnrml2<nnmet>(met));
        errminmetL_aniso = MIN(errminmetL_aniso, errmetL);
        errmaxmetL_aniso = MAX(errmaxmetL_aniso, errmetL);
        erravgmetL_aniso += errmetL;

        double erreigL = sqrt(geterrl2<ndim>(eigval,eigva2) / getnrml2<ndim>(eigval));
        errmineigL_aniso = MIN(errmineigL_aniso, erreigL);
        errmaxeigL_aniso = MAX(errmaxeigL_aniso, erreigL);
        erravgeigL_aniso += erreigL;




        // ---- Test Eigen

        typedef Eigen::Matrix<double,ndim,ndim> MatrixN;
        typedef Eigen::Vector<double,ndim> VectorN;
        MatrixN met_Eigen;

        for(int ii = 0; ii < ndim; ii++) 
          for(int jj = 0; jj < ndim; jj++) 
            met_Eigen(ii, jj) = met[sym2idx(ii,jj)];

        METRIS_ENFORCE(met_Eigen.isApprox(met_Eigen.transpose()))

        Eigen::SelfAdjointEigenSolver<MatrixN> solver(met_Eigen);

        VectorN eigenvalues  = solver.eigenvalues();
        MatrixN eigenvectors = solver.eigenvectors();

        for(int ii = 0; ii < ndim; ii++) eigva2[ii] = eigenvalues[ii];
        for(int ii = 0; ii < ndim; ii++) 
          for(int jj = 0; jj < ndim; jj++) 
            eigve2[jj][ii] = eigenvectors(ii,jj);

        eig2met<ndim,double>(eigva2,eigve2[0],met2);
        double errmetE = sqrt(geterrl2<nnmet>(met2,met) / getnrml2<nnmet>(met));
        errminmetE_aniso = MIN(errminmetE_aniso, errmetE);
        errmaxmetE_aniso = MAX(errmaxmetE_aniso, errmetE);
        erravgmetE_aniso += errmetE;

        double erreigE = sqrt(geterrl2<ndim>(eigval,eigva2) / getnrml2<ndim>(eigval));
        errmineigE_aniso = MIN(errmineigE_aniso, erreigE);
        errmaxeigE_aniso = MAX(errmaxeigE_aniso, erreigE);
        erravgeigE_aniso += erreigE;

      }// for isamp

      erravgmet_aniso /= nerro_aniso;
      erravgmetL_aniso /= nerro_aniso;
      erravgmetE_aniso /= nerro_aniso;
      erravgeig_aniso /= nerro_aniso;
      erravgeigL_aniso /= nerro_aniso;
      erravgeigE_aniso /= nerro_aniso;

      printf("   - metric error min = %e avg = %e max = %e\n",errminmet_aniso,erravgmet_aniso,errmaxmet_aniso);
      printf("   - eigenv error min = %e avg = %e max = %e\n",errmineig_aniso,erravgeig_aniso,errmaxeig_aniso);
      printf("   - metric error min = %e avg = %e max = %e (LAPACK)\n",errminmetL_aniso,erravgmetL_aniso,errmaxmetL_aniso);
      printf("   - eigenv error min = %e avg = %e max = %e (LAPACK)\n",errmineigL_aniso,erravgeigL_aniso,errmaxeigL_aniso);
      printf("   - metric error min = %e avg = %e max = %e (Eigen)\n",errminmetE_aniso,erravgmetE_aniso,errmaxmetE_aniso);
      printf("   - eigenv error min = %e avg = %e max = %e (Eigen)\n",errmineigE_aniso,erravgeigE_aniso,errmaxeigE_aniso);

      BOOST_TEST(errmaxmet_aniso <= tol);
      BOOST_TEST(errmaxeig_aniso <= tol);

      BOOST_TEST(errmaxmetL_aniso <= tol);
      BOOST_TEST(errmaxeigL_aniso <= tol);

      BOOST_TEST(errmaxmetE_aniso <= tol);
      BOOST_TEST(errmaxeigE_aniso <= tol);




      // ------------------------------------------------------------ Benchmark
      #ifdef NDEBUG
      double t0_DSYEVQ = get_wall_time();
      double dum_DSYEVQ = 0;
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met[ii] = met_samples(isamp, ii);
        double rwork[10];
        geteigsym<ndim,double>(met,eigva2,eigve2[0]);
        dum_DSYEVQ += eigva2[0];
      }
      double t1_DSYEVQ = get_wall_time();


      double t0_LAPACK = get_wall_time();
      double dum_LAPACK = 0;
      int rwork[10];
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met[ii] = met_samples(isamp, ii);
        double rwork[10];
        geteigsym_LAPACK<ndim>(met,10,rwork,eigva2,eigve2[0]);
        dum_LAPACK += eigva2[0];
      }
      double t1_LAPACK = get_wall_time();


      double t0_EIGEN = get_wall_time();
      double dum_EIGEN = 0;
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met[ii] = met_samples(isamp, ii);
        typedef Eigen::Matrix<double,ndim,ndim> MatrixN;
        typedef Eigen::Vector<double,ndim> VectorN;
        MatrixN met_Eigen;

        for(int ii = 0; ii < ndim; ii++) 
          for(int jj = 0; jj < ndim; jj++) 
            met_Eigen(ii, jj) = met[sym2idx(ii,jj)];

        Eigen::SelfAdjointEigenSolver<MatrixN> solver(met_Eigen);

        VectorN eigenvalues  = solver.eigenvalues();
        MatrixN eigenvectors = solver.eigenvectors();
        dum_EIGEN += eigenvalues[0];
      }
      double t1_EIGEN = get_wall_time();


      double t0_EIGEN2 = get_wall_time();
      double dum_EIGEN2 = 0;
      typedef Eigen::Matrix<double,ndim,ndim> MatrixN;
      typedef Eigen::Vector<double,ndim> VectorN;
      MatrixN met_Eigen2;
      VectorN eigenvalues;
      MatrixN eigenvectors;
      Eigen::SelfAdjointEigenSolver<MatrixN> solver;
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met[ii] = met_samples(isamp, ii);

        for(int ii = 0; ii < ndim; ii++) 
          for(int jj = ii; jj < ndim; jj++) 
            met_Eigen2(ii, jj) = met[sym2idx(ii,jj)];

        solver.compute(met_Eigen2.template selfadjointView<Eigen::Upper>());
        //Eigen::SelfAdjointEigenSolver<MatrixN> solver(met_Eigen2.template selfadjointView<Eigen::Upper>());

        eigenvalues  = solver.eigenvalues();
        eigenvectors = solver.eigenvectors();
        dum_EIGEN2 += eigenvalues[0];
      }
      double t1_EIGEN2 = get_wall_time();


      printf("  -- DONE benchmarks DSYEVQ time : %8.2e = %d op/s\n",
                  t1_DSYEVQ-t0_DSYEVQ,(int)(nerro_aniso/(t1_DSYEVQ-t0_DSYEVQ)));
      printf("                     LAPACK time : %8.2e = %d op/s, fac = %4.2fx\n",
                  t1_LAPACK-t0_LAPACK,(int)(nerro_aniso/(t1_LAPACK-t0_LAPACK)),
                  (t1_LAPACK-t0_LAPACK)/(t1_DSYEVQ-t0_DSYEVQ));
      printf("                      Eigen time : %8.2e = %d op/s, fac = %4.2fx\n",
                  t1_EIGEN-t0_EIGEN,(int)(nerro_aniso/(t1_EIGEN-t0_EIGEN)),
                  (t1_EIGEN-t0_EIGEN)/(t1_DSYEVQ-t0_DSYEVQ));
      printf("                     Eigen2 time : %8.2e = %d op/s, fac = %4.2fx\n",
                  t1_EIGEN2-t0_EIGEN2,(int)(nerro_aniso/(t1_EIGEN2-t0_EIGEN2)),
                  (t1_EIGEN2-t0_EIGEN2)/(t1_DSYEVQ-t0_DSYEVQ));


      #endif

      printf("  -- DONE with aniso ratio %5.1e\n",anisorat);
    }// for anisorat



    printf("-- All tests done for ndim = %d \n",ndim);
  }CT_FOR1(ndim);

}


template <int ndim>
void generate_metric(double anisorat, double eigval[ndim], double eigvec[ndim][ndim],
  std::uniform_real_distribution<double>& unif, std::default_random_engine& rng){

  for(int idim = 0; idim < ndim; idim++){
    eigval[idim] = 1.0/unif(rng);
  }

  sortupto8_inc(eigval, ndim);
  // eigval[0] is the smallest eigval, i.e. largest size. 
  // To simulate the kind of anisotropic metrics we get, we rather decrease
  // the smallest size.
  eigval[ndim-1] *= eigval[0] / eigval[ndim-1] * anisorat * anisorat;

  for(int ii = 0; ii <= ndim - 2; ii++){
    if(ii > 0 && eigval[ii] > eigval[ii+1]) eigval[ii] = sqrt(eigval[ii-1]*eigval[ii+1]);
    METRIS_ENFORCE(eigval[ii] <= eigval[ii+1]);
  }

  // Generate Frénet frame
  while(true){
    for(int jj = 0; jj < ndim; jj++){
     eigvec[0][jj] = unif(rng);
    }
    double nrm = sqrt(getnrml2<ndim>(eigvec[0]));

    if(nrm < 1.0e-12) continue;

    for(int jj = 0; jj < ndim; jj++){
     eigvec[0][jj] = eigvec[0][jj] / nrm;
    }

    bool ifnd = false;
    for(int jj = 0; jj < ndim; jj++){
      if( abs(eigvec[0][jj]) > 1.0e-6 ){
        for(int kk = 0; kk < ndim; kk++){
          eigvec[1][kk] = 0;
        }
        eigvec[1][(jj+1)%ndim] = -eigvec[0][jj];
        eigvec[1][jj]          = eigvec[0][(jj+1)%ndim];
        nrm = sqrt(getnrml2<ndim>(eigvec[1]));
        for(int jj = 0; jj < ndim; jj++){
         eigvec[1][jj] = eigvec[1][jj] / nrm;
        }

        ifnd = true;
        break;
      }
    }

    if(!ifnd) continue;

    METRIS_ENFORCE( abs(getprdl2<ndim>(eigvec[0],eigvec[1])) < 1.0e-12);

    if constexpr(ndim == 3){
      vecprod(eigvec[0],eigvec[1],eigvec[2]);
      double nrm = sqrt(getnrml2<ndim>(eigvec[2]));
      for(int jj = 0; jj < ndim; jj++){
       eigvec[2][jj] = eigvec[2][jj] / nrm;
      }

      METRIS_ENFORCE( abs(getnrml2<3>(eigvec[2]) - 1) < 1.0e-12);
      METRIS_ENFORCE( abs(getprdl2<3>(eigvec[0],eigvec[1])) < 1.0e-12);
      METRIS_ENFORCE( abs(getprdl2<3>(eigvec[2],eigvec[1])) < 1.0e-12);
      METRIS_ENFORCE( abs(getprdl2<3>(eigvec[2],eigvec[0])) < 1.0e-12);
    }

    break;
  }
}

template 
void generate_metric<2>(double aniso, double eigval[2], double eigvec[2][2],
  std::uniform_real_distribution<double>& unif, std::default_random_engine& rng);
template 
void generate_metric<3>(double aniso, double eigval[3], double eigvec[3][3],
  std::uniform_real_distribution<double>& unif, std::default_random_engine& rng);


}