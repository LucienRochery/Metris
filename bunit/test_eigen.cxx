//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_eigen

#include <boost/test/included/unit_test.hpp> 
#include "common_setup.hxx"

#include <boost/timer/progress_display.hpp>
#include <random>
#include <fstream>

//#include "../src/utils/aux_utils.hxx"
#include "../src/low_geo/misc.hxx"
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

BOOST_AUTO_TEST_CASE(test_eigen) 
{//METRIS_MAX_DEG

  const int nsamp = 1e5;
  const double tol = 1.0e-11;

  std::uniform_real_distribution<double> unif(1.0e-16,1.0);
  std::default_random_engine rng(0);

  double aniso_max = 2e7;
  double aniso_mul = 10;


  // i = idim - 2, j = 0 ? avg : max;
  scriptArrayString errseigD[2][2], errseigE[2][2], errseigL[2][2];
  scriptArrayString errseicD[2][2], errseicE[2][2], errseicL[2][2];

  scriptArrayString benchD[2], benchE[2], benchL[2];

  scriptArrayString anisos("aniso");
  for(double anisorat = 2; anisorat <= aniso_max + 1; anisorat *= aniso_mul){
    anisos += anisorat;
  }
  anisos.finish();

  for(int ii = 0; ii < 2; ii++){
    errseigD[ii][0].setName("erreig_avg_DSYEVQ"+std::to_string(ii+2));
    errseigD[ii][1].setName("erreig_max_DSYEVQ"+std::to_string(ii+2));

    #ifdef METRIS_USE_LAPACK
    errseigL[ii][0].setName("erreig_avg_LAPACK"+std::to_string(ii+2));
    errseigL[ii][1].setName("erreig_max_LAPACK"+std::to_string(ii+2));
    #endif

    errseigE[ii][0].setName("erreig_avg_EIGEN"+std::to_string(ii+2));
    errseigE[ii][1].setName("erreig_max_EIGEN"+std::to_string(ii+2));

    errseicD[ii][0].setName("erreic_avg_DSYEVQ"+std::to_string(ii+2));
    errseicD[ii][1].setName("erreic_max_DSYEVQ"+std::to_string(ii+2));

    #ifdef METRIS_USE_LAPACK
    errseicL[ii][0].setName("erreic_avg_LAPACK"+std::to_string(ii+2));
    errseicL[ii][1].setName("erreic_max_LAPACK"+std::to_string(ii+2));
    #endif

    errseicE[ii][0].setName("erreic_avg_EIGEN"+std::to_string(ii+2));
    errseicE[ii][1].setName("erreic_max_EIGEN"+std::to_string(ii+2));

    benchD[ii].setName("time_DSYEVQ"+std::to_string(ii+2));
    #ifdef METRIS_USE_LAPACK
    benchL[ii].setName("time_LAPACK"+std::to_string(ii+2));
    #endif
    benchE[ii].setName("time_EIGEN"+std::to_string(ii+2));
  }


  // Tests
  dblAr2 met_samples[4];
  met_samples[2].allocate(nsamp,(2*(2+1))/2);
  met_samples[2].set_n(nsamp);
  met_samples[3].allocate(nsamp,(3*(3+1))/2);
  met_samples[3].set_n(nsamp);
  CT_FOR0_INC(2,3,ndim){
    constexpr int nnmet = (ndim*(ndim+1))/2;

    double eigval[ndim],eigvec[ndim][ndim];
    double eigva2[ndim],eigve2[ndim][ndim];
    double met[nnmet], met2[nnmet];

    printf("-- Tests for dim = %d \n",ndim);

    dblAr1 errD_series, aniso_series;

    for(double anisorat = 2; anisorat <= aniso_max + 1; anisorat *= aniso_mul){

      MinMaxAvg errmetD, errmetL, errmetE;
      MinMaxAvg erreigD, erreigL, erreigE;
      MinMaxAvg erreicD, erreicL, erreicE;

      dblAr2 eigval_samples(nsamp, ndim);
      dblAr2 eigvec_samples(nsamp, ndim*ndim);
      eigval_samples.set_n(nsamp);
      eigvec_samples.set_n(nsamp);


      double err;

      for(int isamp = 0; isamp < nsamp; isamp++){

        generate_metric<ndim>(anisorat, eigval, eigvec, unif, rng);
        for(int ii = 0; ii < ndim; ii++){
          eigval_samples(isamp, ii) = eigval[ii];
          for(int jj = 0; jj < ndim; jj++){
            eigvec_samples(isamp, ndim*ii+jj) = eigvec[ii][jj];
          }
        }

        // Generate metric from eigvecs, eigvals
        eig2met<ndim,double>(eigval,eigvec[0],met);
        for(int ii = 0; ii < nnmet; ii++) met_samples[ndim](isamp, ii) = met[ii];


        // ----- Test geteigsym (SurrealS compatible dsyevq)

        // Re-decompose and compare
        geteigsym<ndim,double>(met,eigva2,eigve2[0]);

        sorteig<ndim, double>(eigva2,eigve2[0]);

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

        err = sqrt(geterrl2<nnmet>(met2,met) / getnrml2<nnmet>(met));
        errmetD += err;

        err = sqrt(geterrl2<ndim>(eigval,eigva2) / getnrml2<ndim>(eigval));
        erreigD += err;

        err = 0;
        for(int ii = 0; ii < ndim; ii++){
          err += abs(eigval[ii] - eigva2[ii])/abs(eigval[ii]);
        }
        erreicD += err;



        // ----- Test LAPACK
      #ifdef METRIS_USE_LAPACK
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

        err = sqrt(geterrl2<nnmet>(met2,met) / getnrml2<nnmet>(met));
        errmetL += err;

        err = sqrt(geterrl2<ndim>(eigval,eigva2) / getnrml2<ndim>(eigval));
        erreigL += err;

        err = 0;
        for(int ii = 0; ii < ndim; ii++){
          err += abs(eigval[ii] - eigva2[ii])/abs(eigval[ii]);
        }
        erreicL += err;
      #endif



        // ---- Test Eigen

        int ierro = geteigsym_Eigen<ndim, double>(met, eigva2, eigve2[0]);
        #if 0
        typedef Eigen::Matrix<double,ndim,ndim> MatrixN;
        typedef Eigen::Vector<double,ndim> VectorN;
        MatrixN met_Eigen;

        for(int ii = 0; ii < ndim; ii++) 
          for(int jj = 0; jj <= ii; jj++) 
            met_Eigen(ii, jj) = met[sym2idx(ii,jj)];

        //METRIS_ENFORCE(met_Eigen.isApprox(met_Eigen.transpose()))

        Eigen::SelfAdjointEigenSolver<MatrixN> solver(met_Eigen);

        VectorN eigenvalues  = solver.eigenvalues();
        MatrixN eigenvectors = solver.eigenvectors();

        for(int ii = 0; ii < ndim; ii++) eigva2[ii] = eigenvalues[ii];
        for(int ii = 0; ii < ndim; ii++) 
          for(int jj = 0; jj < ndim; jj++) 
            eigve2[jj][ii] = eigenvectors(ii,jj);
        #endif

        eig2met<ndim,double>(eigva2,eigve2[0],met2);
        err = sqrt(geterrl2<nnmet>(met2,met) / getnrml2<nnmet>(met));
        errmetE += err;

        err = sqrt(geterrl2<ndim>(eigval,eigva2) / getnrml2<ndim>(eigval));
        erreigE += err;

        err = 0;
        for(int ii = 0; ii < ndim; ii++){
          err += abs(eigval[ii] - eigva2[ii])/abs(eigval[ii]);
        }
        erreicE += err;

      }// for isamp

      errseigD[ndim-2][0] += erreigD.avg();
      errseigL[ndim-2][0] += erreigL.avg();
      errseigE[ndim-2][0] += erreigE.avg();

      errseicD[ndim-2][0] += erreicD.avg();
      errseicL[ndim-2][0] += erreicL.avg();
      errseicE[ndim-2][0] += erreicE.avg();

      errseigD[ndim-2][1] += erreigD.max();
      errseigL[ndim-2][1] += erreigL.max();
      errseigE[ndim-2][1] += erreigE.max();

      errseicD[ndim-2][1] += erreicD.max();
      errseicL[ndim-2][1] += erreicL.max();
      errseicE[ndim-2][1] += erreicE.max();

      printf("   - metric error min = %e avg = %e max = %e (DSYEV)\n",errmetD.min(),errmetD.avg(),errmetD.max());
      printf("   - eigvec error min = %e avg = %e max = %e (DSYEV)\n",erreigD.min(),erreigD.avg(),erreigD.max());
      printf("   - eigrel error min = %e avg = %e max = %e (DSYEV)\n",erreicD.min(),erreicD.avg(),erreicD.max());

      #ifdef METRIS_USE_LAPACK
      printf("   - metric error min = %e avg = %e max = %e (LAPACK)\n",errmetL.min(),errmetL.avg(),errmetL.max());
      printf("   - eigvec error min = %e avg = %e max = %e (LAPACK)\n",erreigL.min(),erreigL.avg(),erreigL.max());
      printf("   - eigrel error min = %e avg = %e max = %e (LAPACK)\n",erreicL.min(),erreicL.avg(),erreicL.max());
      #endif

      printf("   - metric error min = %e avg = %e max = %e (Eigen)\n",errmetE.min(),errmetE.avg(),errmetE.max());
      printf("   - eigvec error min = %e avg = %e max = %e (Eigen)\n",erreigE.min(),erreigE.avg(),erreigE.max());
      printf("   - eigrel error min = %e avg = %e max = %e (Eigen)\n",erreicE.min(),erreicE.avg(),erreicE.max());

      BOOST_TEST(errmetD <= tol);
      BOOST_TEST(erreigD <= tol);

      BOOST_TEST(errmetL <= tol);
      BOOST_TEST(erreigL <= tol);

      BOOST_TEST(errmetE <= tol);
      BOOST_TEST(erreigE <= tol);


      errD_series.stack(erreicD.avg());
      aniso_series.stack(anisorat);

      printf("  -- DONE with aniso ratio %5.1e\n",anisorat);
    }// for anisorat


    // Linear regression on dsyevq (what is used in Metris) 
    // check O(aniso^2) as empirically observed (note matrix conditioning is aniso^2)
    // and expected under 10^{-13} at aniso 10.
    int nseries = errD_series.get_n();
    for(int ii = 0; ii < nseries; ii++){
      errD_series[ii] = log(errD_series[ii]);
      aniso_series[ii] = log(aniso_series[ii]);
    }
    LinReg linreg(aniso_series,errD_series);
    BOOST_TEST(linreg.slope <= 2.5);

    double err10 = linreg.origin + linreg.slope*log(10);
    BOOST_TEST(err10 < log(1.0e-13));


    printf("-- All tests done for ndim = %d \n",ndim);
  }CT_FOR1(ndim);


  std::ofstream fil;
  std::string fname = "test_eigen.py";
  fil.open(fname.c_str(), std::ios::out);
  fil << "import matplotlib.pyplot as plt\n\n";
  fil << anisos.str() << "\n\n";
  int ifig = 0;
  for(int ndim = 0; ndim < 2; ndim++){
    for(int iar = 0; iar < 2; iar++){

      std::string errtype = iar == 0 ? "avg" : "max";

      errseigD[ndim][iar].finish();
      #ifdef METRIS_USE_LAPACK
      errseigL[ndim][iar].finish();
      #endif
      errseigE[ndim][iar].finish();
      fil << errseigD[ndim][iar].str() << "\n";
      fil << errseigL[ndim][iar].str() << "\n";
      fil << errseigE[ndim][iar].str() << "\n";
      fil << "plt.figure("<<std::to_string(++ifig)<<")\n";
      fil << "plt.loglog(aniso," << errseigD[ndim][iar].name()<<",label='DSYEVQ')\n";
      #ifdef METRIS_USE_LAPACK
      fil << "plt.loglog(aniso," << errseigL[ndim][iar].name()<<",label='LAPACK')\n";
      #endif
      fil << "plt.loglog(aniso," << errseigE[ndim][iar].name()<<",label='Eigen')\n";
      fil << "plt.title('Eigenvalue vector relative error (" << errtype << "), dim "
          << std::to_string(ndim+2)<< "')\n";
      fil << "plt.xlabel('anisotropic ratio')\n";
      fil << "plt.ylabel('relative error')\n";
      fil << "plt.legend()\n";
      fil << "plt.grid(True)\n";
      fil << "plt.savefig('errseig_eigen_"<<errtype<<std::to_string(ndim+2)<<".png')\n";
      fil << "\n\n";


      errseicD[ndim][iar].finish();
      #ifdef METRIS_USE_LAPACK
      errseicL[ndim][iar].finish();
      #endif
      errseicE[ndim][iar].finish();
      fil << errseicD[ndim][iar].str() << "\n";
      #ifdef METRIS_USE_LAPACK
      fil << errseicL[ndim][iar].str() << "\n";
      #endif
      fil << errseicE[ndim][iar].str() << "\n";
      fil << "plt.figure("<<std::to_string(++ifig)<<")\n";
      fil << "plt.loglog(aniso," << errseicD[ndim][iar].name()<<",label='DSYEVQ')\n";
      #ifdef METRIS_USE_LAPACK
      fil << "plt.loglog(aniso," << errseicL[ndim][iar].name()<<",label='LAPACK')\n";
      #endif
      fil << "plt.loglog(aniso," << errseicE[ndim][iar].name()<<",label='Eigen')\n";
      fil << "plt.title('Eigenvalue coordinate-wise relative error (" << errtype << "), dim "
          << std::to_string(ndim+2)<< "')\n";
      fil << "plt.xlabel('anisotropic ratio')\n";
      fil << "plt.ylabel('relative error')\n";
      fil << "plt.legend()\n";
      fil << "plt.grid(True, which = 'both', linestyle='--')\n";
      fil << "plt.savefig('errseic_eigen_"<<errtype<<std::to_string(ndim+2)<<".png')\n";
      fil << "\n\n";
    }
  }
  fil << "plt.show()\n";
  fil.close();



  // ------------------------------------------------------------ Benchmark
  CT_FOR0_INC(2,3,ndim){
    constexpr int nnmet = (ndim*(ndim+1))/2;

    double eigval[ndim],eigvec[ndim][ndim];
    double eigva2[ndim],eigve2[ndim][ndim];
    double met[nnmet], met2[nnmet];

    printf("-- Tests for dim = %d \n",ndim);

    dblAr1 errD_series, aniso_series;

    for(double anisorat = 2; anisorat <= aniso_max + 1; anisorat *= aniso_mul){


      #ifdef NDEBUG
      double t0_DSYEVQ = get_wall_time();
      double dum_DSYEVQ = 0;
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met[ii] = met_samples[ndim](isamp, ii);
        double rwork[10];
        geteigsym<ndim,double>(met,eigva2,eigve2[0]);
        dum_DSYEVQ += eigva2[0];
      }
      double t1_DSYEVQ = get_wall_time();


      double dum_LAPACK = 1;
    #ifdef METRIS_USE_LAPACK
      double t0_LAPACK = get_wall_time();
      int rwork[10];
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met[ii] = met_samples[ndim](isamp, ii);
        double rwork[10];
        geteigsym_LAPACK<ndim>(met,10,rwork,eigva2,eigve2[0]);
        dum_LAPACK += eigva2[0];
      }
      double t1_LAPACK = get_wall_time();
    #endif


      double t0_EIGEN = get_wall_time();
      double dum_EIGEN = 0;
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met[ii] = met_samples[ndim](isamp, ii);
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
        for(int ii = 0; ii < nnmet; ii++) met[ii] = met_samples[ndim](isamp, ii);

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
                  t1_DSYEVQ-t0_DSYEVQ,(int)(nsamp/(t1_DSYEVQ-t0_DSYEVQ)));
      #ifdef METRIS_USE_LAPACK
      printf("                     LAPACK time : %8.2e = %d op/s, fac = %4.2fx\n",
                  t1_LAPACK-t0_LAPACK,(int)(nsamp/(t1_LAPACK-t0_LAPACK)),
                  (t1_LAPACK-t0_LAPACK)/(t1_DSYEVQ-t0_DSYEVQ));
      #endif
      printf("                      Eigen time : %8.2e = %d op/s, fac = %4.2fx\n",
                  t1_EIGEN-t0_EIGEN,(int)(nsamp/(t1_EIGEN-t0_EIGEN)),
                  (t1_EIGEN-t0_EIGEN)/(t1_DSYEVQ-t0_DSYEVQ));
      printf("                     Eigen2 time : %8.2e = %d op/s, fac = %4.2fx\n",
                  t1_EIGEN2-t0_EIGEN2,(int)(nsamp/(t1_EIGEN2-t0_EIGEN2)),
                  (t1_EIGEN2-t0_EIGEN2)/(t1_DSYEVQ-t0_DSYEVQ));

      benchD[ndim-2] += nsamp/(t1_DSYEVQ - t0_DSYEVQ);
      #ifdef METRIS_USE_LAPACK
      benchL[ndim-2] += nsamp/(t1_LAPACK - t0_LAPACK);
      #endif
      benchE[ndim-2] += nsamp/(t1_EIGEN  - t0_EIGEN );

      #endif

      printf("  -- DONE with aniso ratio %5.1e\n",anisorat);
    }// for anisorat

    printf("-- All tests done for ndim = %d \n",ndim);
  }CT_FOR1(ndim);

  fname = "bench_eigen.py";
  fil.open(fname.c_str(), std::ios::out);
  fil << "import matplotlib.pyplot as plt\n";
  fil << "import numpy as np\n";
  fil << "\n";

  fil << anisos.str() << "\n\n";
  ifig = 0;
  for(int ndim = 0; ndim < 2; ndim++){

    benchD[ndim].finish();
    #ifdef METRIS_USE_LAPACK
    benchL[ndim].finish();
    #endif
    benchE[ndim].finish();
    fil << benchD[ndim].str() << "\n";
    #ifdef METRIS_USE_LAPACK
    fil << benchL[ndim].str() << "\n";
    #endif
    fil << benchE[ndim].str() << "\n";
    fil << "plt.figure("<<std::to_string(++ifig)<<")\n";
    #ifdef METRIS_USE_LAPACK
    fil << "plt.semilogx(aniso,np.array(" << benchD[ndim].name()<<")"
        << "/np.array(" << benchL[ndim].name() << ")" <<",label='DSYEVQ')\n";
    fil << "plt.semilogx(aniso,np.array(" << benchE[ndim].name()<<")"
        << "/np.array(" << benchL[ndim].name() << ")" <<",label='Eigen')\n";
    fil << "plt.ylabel('op/s divided by LAPACKs')\n";
    #else
    fil << "plt.semilogx(aniso,np.array(" << benchD[ndim].name()<<"))"
        << ",label='DSYEVQ')\n";
    fil << "plt.semilogx(aniso,np.array(" << benchE[ndim].name()<<"))"
        << ",label='Eigen')\n";
    fil << "plt.ylabel('op/s')\n";
    #endif
    fil << "plt.xlabel('anisotropic ratio')\n";
    fil << "plt.title('Benchmark dim "<< std::to_string(ndim+2)<< "')\n";
    fil << "plt.legend()\n";
    fil << "plt.grid(True)\n";
    fil << "plt.savefig('bench_eigen_"<<std::to_string(ndim+2)<<".png')\n";
    fil << "\n\n";
  }
  fil << "plt.show()\n";
  fil.close();


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