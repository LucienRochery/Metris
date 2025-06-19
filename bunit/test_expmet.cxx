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


#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>


namespace Metris{


// -- Test metqua3_xi_d 
// Constant metric fields should yield reliable derivatives in all cases 
// In non constant metric fields, derivatives only defined for DoFs in back element
// interiors... 

template<int ndim>
void generate_metric(double aniso, double eigval[ndim], double eigvec[ndim][ndim],
  MetSpace metspac,
  std::uniform_real_distribution<double>& unif, std::default_random_engine& rng);

BOOST_AUTO_TEST_CASE(test_eval3) 
{//METRIS_MAX_DEG

  const int nsamp = 1e5;
  const double tol = 1.0e-11;

  std::uniform_real_distribution<double> unif(1.0e-16,1.0);
  std::default_random_engine rng(0);

  double aniso_max = 2e6;
  double aniso_mul = 10;

  dblAr2 eigval_samples[2];
  dblAr2 eigvec_samples[2];
  dblAr2 met_samples[2];
  for(int idim = 2; idim <= 3; idim++){
    eigval_samples[idim-2].allocate(nsamp, idim);
    eigval_samples[idim-2].set_n(nsamp);
    
    eigvec_samples[idim-2].allocate(nsamp, idim*idim);
    eigvec_samples[idim-2].set_n(nsamp);

    met_samples[idim-2].allocate(nsamp,(idim*(idim+1))/2);
    met_samples[idim-2].set_n(nsamp);
  }


  CT_FOR0_INC(2,3,ndim){
    constexpr int nnmet = (ndim*(ndim+1))/2;
    typedef Eigen::Matrix<double,ndim,ndim> MatrixN;
    typedef Eigen::Vector<double,ndim> VectorN;

    double eigval[ndim],eigvec[ndim][ndim];
    double eigva2[ndim],eigve2[ndim][ndim];
    double met[nnmet], met2[nnmet], met3[nnmet];
    MatrixN met_Eigen, met2_Eigen;


    for(MetSpace metspac : {MetSpace::Exp, MetSpace::Log}){

      printf("\n-- Test dim %d",ndim);
      if(metspac == MetSpace::Exp) printf(" log\n");
      else                         printf(" exp\n");

      for(double anisorat = 2; anisorat < aniso_max; anisorat *= aniso_mul){

        double errmin1_aniso = 1.0e30;
        double errmax1_aniso = -1;
        double erravg1_aniso = 0;

        double errmin2_aniso = 1.0e30;
        double errmax2_aniso = -1;
        double erravg2_aniso = 0;

        double errmin3_aniso = 1.0e30;
        double errmax3_aniso = -1;
        double erravg3_aniso = 0;


        double nerro_aniso = 0;


        for(int isamp = 0; isamp < nsamp; isamp++){

          generate_metric<ndim>(anisorat, eigval, eigvec, metspac, unif, rng);
          for(int ii = 0; ii < ndim; ii++){
            eigval_samples[ndim-2](isamp, ii) = eigval[ii];
            for(int jj = 0; jj < ndim; jj++){
              eigvec_samples[ndim-2](isamp, ndim*ii+jj) = eigvec[ii][jj];
            }
          }

          // Generate metric from eigvecs, eigvals
          eig2met<ndim,double>(eigval,eigvec[0],met);
          for(int ii = 0; ii < nnmet; ii++) met_samples[ndim-2](isamp, ii) = met[ii];

          for(int ii = 0; ii < nnmet; ii++) met2[ii] = met[ii];
          if(metspac == MetSpace::Exp){
            getlogmet_inp<ndim,double>(met2);
            for(int ii = 0; ii < ndim; ii++) eigva2[ii] = log(eigval[ii]);
          }else{
            getexpmet_inp<ndim,double>(met2);
            for(int ii = 0; ii < ndim; ii++) eigva2[ii] = exp(eigval[ii]);
          }
          eig2met<ndim,double>(eigva2,eigvec[0],met3);

          double err1 = sqrt(geterrl2<nnmet>(met2,met3) / getnrml2<nnmet>(met3));
          errmin1_aniso = MIN(errmin1_aniso, err1);
          errmax1_aniso = MAX(errmax1_aniso, err1);
          erravg1_aniso += err1;
          nerro_aniso++;


          for(int ii = 0; ii < ndim; ii++) 
            for(int jj = 0; jj < ndim; jj++) 
              met_Eigen(ii, jj) = met[sym2idx(ii,jj)];

          if(metspac == MetSpace::Exp){
            met_Eigen.log().evalTo(met2_Eigen);
          }else{
            met_Eigen.exp().evalTo(met2_Eigen); 
          }

          for(int ii = 0; ii < ndim; ii++) 
            for(int jj = 0; jj < ndim; jj++) 
              met2[sym2idx(ii,jj)] = met2_Eigen(ii, jj);
          double err2 = sqrt(geterrl2<nnmet>(met2,met3) / getnrml2<nnmet>(met3));
          errmin2_aniso = MIN(errmin2_aniso, err2);
          errmax2_aniso = MAX(errmax2_aniso, err2);
          erravg2_aniso += err2;




          Eigen::SelfAdjointEigenSolver<MatrixN> solver(met_Eigen);

          VectorN eigenvalues  = solver.eigenvalues();
          MatrixN eigenvectors = solver.eigenvectors();

          if(metspac == MetSpace::Exp){
            eigenvalues = (eigenvalues.array().log()).matrix();
          }else{
            eigenvalues = (eigenvalues.array().exp()).matrix();
          }

          met2_Eigen = eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();
          for(int ii = 0; ii < ndim; ii++) 
            for(int jj = 0; jj < ndim; jj++) 
              met2[sym2idx(ii,jj)] = met2_Eigen(ii, jj);
          double err3 = sqrt(geterrl2<nnmet>(met2,met3) / getnrml2<nnmet>(met3));
          errmin3_aniso = MIN(errmin3_aniso, err3);
          errmax3_aniso = MAX(errmax3_aniso, err3);
          erravg3_aniso += err3;


        }// for isamp

        erravg1_aniso /= nerro_aniso;
        erravg2_aniso /= nerro_aniso;
        erravg3_aniso /= nerro_aniso;

        printf("  -- DONE with aniso ratio %5.1e\n",anisorat);
        printf("   - error min = %e avg = %e max = %e\n",errmin1_aniso,erravg1_aniso,errmax1_aniso);
        printf("   - error min = %e avg = %e max = %e (Eigen)\n",errmin2_aniso,erravg2_aniso,errmax2_aniso);
        printf("   - error min = %e avg = %e max = %e (Eigen2)\n",errmin3_aniso,erravg3_aniso,errmax3_aniso);

        BOOST_TEST(errmax1_aniso <= tol);
        BOOST_TEST(errmax2_aniso <= tol);
        BOOST_TEST(errmax3_aniso <= tol);

      }// for anisorat
    }// for MetSpace
 
  }CT_FOR1(ndim);


  // ---------------------------------------------------------------------
  // ----------------------------------------------------------- Benchmark
  #ifdef NDEBUG
  for(MetSpace metspac : {MetSpace::Exp, MetSpace::Log}){
    printf("\n\n");
    CT_FOR0_INC(2,3,ndim){
      constexpr int nnmet = (ndim*(ndim+1))/2;
      typedef Eigen::Matrix<double,ndim,ndim> MatrixN;
      typedef Eigen::Vector<double,ndim> VectorN;
  
      double eigval[ndim],eigvec[ndim][ndim];
      double eigva2[ndim],eigve2[ndim][ndim];
      double met[nnmet], met2[nnmet], met3[nnmet];
      MatrixN met_Eigen, met2_Eigen;


      for(double anisorat = 2; anisorat < aniso_max; anisorat *= aniso_mul){

        double t0_LAPACK = get_wall_time();
        double dum_LAPACK = 0;
        for(int isamp = 0; isamp < nsamp; isamp++){
          for(int ii = 0; ii < nnmet; ii++) met2[ii] = met_samples[ndim-2](isamp, ii);
          if(metspac == MetSpace::Exp){
            getlogmet_inp<ndim,double>(met2);
          }else{
            getexpmet_inp<ndim,double>(met2);
          }
          dum_LAPACK += met2[0];
        }
        double t1_LAPACK = get_wall_time();


        double t0_EIGEN = get_wall_time();
        double dum_EIGEN = 0;
        for(int isamp = 0; isamp < nsamp; isamp++){
          for(int ii = 0; ii < ndim; ii++) 
            for(int jj = 0; jj < ndim; jj++) 
              met_Eigen(ii, jj) = met_samples[ndim-2](isamp,sym2idx(ii,jj));

          if(metspac == MetSpace::Exp){
            met_Eigen.log().evalTo(met2_Eigen);
          }else{
            met_Eigen.exp().evalTo(met2_Eigen); 
          }

          for(int ii = 0; ii < ndim; ii++) 
            for(int jj = 0; jj < ndim; jj++) 
              met[sym2idx(ii,jj)] = met2_Eigen(ii, jj);
          dum_EIGEN += met[0];
        }
        double t1_EIGEN = get_wall_time();


        double t0_EIGEN2 = get_wall_time();
        double dum_EIGEN2 = 0;
        for(int isamp = 0; isamp < nsamp; isamp++){
          for(int ii = 0; ii < ndim; ii++) 
            for(int jj = 0; jj < ndim; jj++) 
              met_Eigen(ii, jj) = met_samples[ndim-2](isamp,sym2idx(ii,jj));

          Eigen::SelfAdjointEigenSolver<MatrixN> solver(met_Eigen);

          VectorN eigenvalues  = solver.eigenvalues();
          MatrixN eigenvectors = solver.eigenvectors();

          if(metspac == MetSpace::Exp){
            eigenvalues = (eigenvalues.array().log()).matrix();
          }else{
            eigenvalues = (eigenvalues.array().exp()).matrix();
          }

          met_Eigen = eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();

          for(int ii = 0; ii < ndim; ii++) 
            for(int jj = 0; jj < ndim; jj++) 
              met[sym2idx(ii,jj)] = met2_Eigen(ii, jj);
          dum_EIGEN += met[0];
        }
        double t1_EIGEN2 = get_wall_time(); 


        std::string funname = metspac == MetSpace::Exp ? "log" : "exp";
        printf("  -- DONE bench dim %1d aniso %5.1e %s LAPACK time : %8.2es = %dk op/s\n",
                    ndim,anisorat,funname.c_str(),
                    t1_LAPACK-t0_LAPACK,(int)(nsamp/(t1_LAPACK-t0_LAPACK)/1000));
        printf("                                         Eigen time : %8.2es = %dk op/s, fac = %4.2fx\n",
                    t1_EIGEN-t0_EIGEN,(int)(nsamp/(t1_EIGEN-t0_EIGEN)/1000),
                    (t1_EIGEN-t0_EIGEN)/(t1_LAPACK-t0_LAPACK));
        printf("                                        Eigen2 time : %8.2es = %dk op/s, fac = %4.2fx\n",
                    t1_EIGEN2-t0_EIGEN2,(int)(nsamp/(t1_EIGEN2-t0_EIGEN2)/1000),
                    (t1_EIGEN2-t0_EIGEN2)/(t1_LAPACK-t0_LAPACK));


      }// for anisorat
    }CT_FOR1(ndim);
  }// for MetSpace
  #endif

}


template <int ndim>
void generate_metric(double anisorat, double eigval[ndim], double eigvec[ndim][ndim],
  MetSpace metspac,
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

  if(metspac == MetSpace::Log){
    for(int ii = 0; ii < ndim; ii++) eigval[ii] = log(eigval[ii]);
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
  MetSpace metspac,
  std::uniform_real_distribution<double>& unif, std::default_random_engine& rng);
template 
void generate_metric<3>(double aniso, double eigval[3], double eigvec[3][3],
  MetSpace metspac,
  std::uniform_real_distribution<double>& unif, std::default_random_engine& rng);


}