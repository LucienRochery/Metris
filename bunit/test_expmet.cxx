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

BOOST_AUTO_TEST_CASE(test_expmet) 
{//METRIS_MAX_DEG

  const int nsamp = 1e5;
  const double tol = 1.0e-11;
  const double exptol = 1.0e-16;

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


  // i = idim - 2, j = 0 ? log : exp;
  scriptArrayString errsD[2][2], errsE1[2][2], errsE2[2][2];
  scriptArrayString benchD[2][2], benchE1[2][2], benchE2[2][2];

  scriptArrayString anisos("aniso");
  for(double anisorat = 2; anisorat <= aniso_max + 1; anisorat *= aniso_mul){
    anisos += anisorat;
  }
  anisos.finish();

  for(int ii = 0; ii < 2; ii++){
    errsD[ii][0].setName("err_max_log_DSYEVQ_"+std::to_string(ii+2));
    errsD[ii][1].setName("err_max_exp_DSYEVQ_"+std::to_string(ii+2));

    errsE1[ii][0].setName("err_max_log_EIGEN1_"+std::to_string(ii+2));
    errsE1[ii][1].setName("err_max_exp_EIGEN1_"+std::to_string(ii+2));

    errsE2[ii][0].setName("err_max_log_EIGEN2_"+std::to_string(ii+2));
    errsE2[ii][1].setName("err_max_exp_EIGEN2_"+std::to_string(ii+2));

    benchD[ii][0].setName("time_max_log_DSYEVQ_"+std::to_string(ii+2));
    benchD[ii][1].setName("time_max_exp_DSYEVQ_"+std::to_string(ii+2));

    benchE1[ii][0].setName("time_max_log_EIGEN1_"+std::to_string(ii+2));
    benchE1[ii][1].setName("time_max_exp_EIGEN1_"+std::to_string(ii+2));

    benchE2[ii][0].setName("time_max_log_EIGEN2_"+std::to_string(ii+2));
    benchE2[ii][1].setName("time_max_exp_EIGEN2_"+std::to_string(ii+2));
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

      for(double anisorat = 2; anisorat < aniso_max + 1; anisorat *= aniso_mul){

        MinMaxAvg errD, errE1, errE2;


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


          // --- Using DSYEVQ to compute eigvals, then log/exp eigval
          for(int ii = 0; ii < nnmet; ii++) met2[ii] = met[ii];
          if(metspac == MetSpace::Exp){
            getlogmet_inp<ndim,double>(met2);
            for(int ii = 0; ii < ndim; ii++) eigva2[ii] = log(eigval[ii]);
          }else{
            getexpmet_dsyevq<ndim,double>(met2);
            for(int ii = 0; ii < ndim; ii++) eigva2[ii] = exp(eigval[ii]);
          }
          eig2met<ndim,double>(eigva2,eigvec[0],met3);

          errD += sqrt(geterrl2<nnmet>(met2,met3) / getnrml2<nnmet>(met3));



          // --- Eigen using log() / exp()
          for(int jj = 0; jj < ndim; jj++) 
            for(int ii = 0; ii < ndim; ii++) 
              met_Eigen(ii, jj) = met[sym2idx(ii,jj)];

          if(metspac == MetSpace::Exp){
            met_Eigen.log().evalTo(met_Eigen);
          }else{
            met_Eigen.exp().evalTo(met_Eigen); 
          }

          for(int ii = 0; ii < ndim; ii++) 
            for(int jj = 0; jj < ndim; jj++) 
              met2[sym2idx(ii,jj)] = met_Eigen(ii, jj);
          errE1 += sqrt(geterrl2<nnmet>(met2,met3) / getnrml2<nnmet>(met3));


          // --- Eigen using eigvals then log/exp eigvals
          for(int jj = 0; jj < ndim; jj++) 
            for(int ii = 0; ii < ndim; ii++) 
              met_Eigen(ii, jj) = met[sym2idx(ii,jj)];
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
          errE2 += sqrt(geterrl2<nnmet>(met2,met3) / getnrml2<nnmet>(met3));


        }// for isamp

        int ifun = metspac == MetSpace::Exp ? 0 : 1;

        errsD[ndim-2][ifun]  += errD.max();
        errsE1[ndim-2][ifun] += errE1.max();
        errsE2[ndim-2][ifun] += errE2.max();


        printf("  -- DONE with aniso ratio %5.1e\n",anisorat);
        printf("   - error min = %e avg = %e max = %e (DSYEVQ)\n",errD.min(),errD.avg(),errD.max());
        printf("   - error min = %e avg = %e max = %e (Eigen direct)\n",errE1.min(),errE1.avg(),errE1.max());
        printf("   - error min = %e avg = %e max = %e (Eigen decomp)\n",errE2.min(),errE2.avg(),errE2.max());

        BOOST_TEST(errD <= tol);
        BOOST_TEST(errE1 <= tol);
        BOOST_TEST(errE2 <= tol);

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


      for(double anisorat = 2; anisorat <= aniso_max + 1; anisorat *= aniso_mul){

        // Using DSYEVQ
        double t0_DSYEVQ = get_wall_time();
        double dum_DSYEVQ = 0;
        for(int isamp = 0; isamp < nsamp; isamp++){
          for(int ii = 0; ii < nnmet; ii++) met2[ii] = met_samples[ndim-2](isamp, ii);
          if(metspac == MetSpace::Exp){
            getlogmet_inp<ndim,double>(met2);
          }else{
            getexpmet_inp<ndim,double>(met2);
          }
          dum_DSYEVQ += met2[0];
        }
        double t1_DSYEVQ = get_wall_time();


        double t0_EIGEN1 = get_wall_time();
        double dum_EIGEN1 = 0;
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
          dum_EIGEN1 += met[0];
        }
        double t1_EIGEN1 = get_wall_time();


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
          dum_EIGEN2 += met[0];
        }
        double t1_EIGEN2 = get_wall_time(); 


        std::string funname = metspac == MetSpace::Exp ? "log" : "exp";
        printf("  -- DONE bench dim %1d aniso %5.1e %s DSYEVQ time : %8.2es = %dk op/s\n",
                    ndim,anisorat,funname.c_str(),
                    t1_DSYEVQ-t0_DSYEVQ,(int)(nsamp/(t1_DSYEVQ-t0_DSYEVQ)/1000));
        printf("                                         Eigen time : %8.2es = %dk op/s, fac = %4.2fx\n",
                    t1_EIGEN1-t0_EIGEN1,(int)(nsamp/(t1_EIGEN1-t0_EIGEN1)/1000),
                    (t1_EIGEN1-t0_EIGEN1)/(t1_DSYEVQ-t0_DSYEVQ));
        printf("                                        Eigen2 time : %8.2es = %dk op/s, fac = %4.2fx\n",
                    t1_EIGEN2-t0_EIGEN2,(int)(nsamp/(t1_EIGEN2-t0_EIGEN2)/1000),
                    (t1_EIGEN2-t0_EIGEN2)/(t1_DSYEVQ-t0_DSYEVQ));

      int ifun = metspac == MetSpace::Exp ? 0 : 1;
      benchD[ndim-2][ifun]  += nsamp/(t1_DSYEVQ - t0_DSYEVQ);
      benchE1[ndim-2][ifun] += nsamp/(t1_EIGEN1 - t0_EIGEN1);
      benchE2[ndim-2][ifun] += nsamp/(t1_EIGEN2 - t0_EIGEN2);

      }// for anisorat
    }CT_FOR1(ndim);
  }// for MetSpace
  #endif

  std::ofstream fil;
  std::string fname = "test_explogmet.py";
  fil.open(fname.c_str(), std::ios::out);
  fil << "import matplotlib.pyplot as plt\n\n";
  fil << anisos.str() << "\n\n";
  int ifig = 0;
  for(int ndim = 0; ndim < 2; ndim++){
    for(int iar = 0; iar < 2; iar++){
      std::string funname = iar == 0 ? "log" : "exp";
      errsD[ndim][iar].finish();
      errsE1[ndim][iar].finish();
      errsE2[ndim][iar].finish();
      fil << errsD[ndim][iar].str() << "\n";
      fil << errsE1[ndim][iar].str() << "\n";
      fil << errsE2[ndim][iar].str() << "\n";
      fil << "plt.figure("<<std::to_string(++ifig)<<")\n";
      fil << "plt.loglog(aniso," << errsD[ndim][iar].name()<<",label='DSYEVQ')\n";
      fil << "plt.loglog(aniso," << errsE1[ndim][iar].name()<<",label='Eigen1')\n";
      fil << "plt.loglog(aniso," << errsE2[ndim][iar].name()<<",label='Eigen2')\n";
      fil << "plt.title('Output relative error "<<funname<<" (max), dim "
          << std::to_string(ndim+2)<< "')\n";
      fil << "plt.xlabel('anisotropic ratio')\n";
      fil << "plt.ylabel('relative error')\n";
      fil << "plt.legend()\n";
      fil << "plt.grid(True)\n";
      fil << "plt.savefig('errs_explogmet_"<<funname<<std::to_string(ndim+2)<<".png')\n";
      fil << "\n\n";
    }
  }
  fil << "plt.show()\n";
  fil.close();

  fname = "bench_explogmet.py";
  fil.open(fname.c_str(), std::ios::out);
  fil << "import matplotlib.pyplot as plt\n";
  fil << "import numpy as np\n";
  fil << "\n";

  fil << anisos.str() << "\n\n";
  ifig = 0;
  for(int ndim = 0; ndim < 2; ndim++){
    for(int iar = 0; iar < 2; iar++){
      std::string funname = iar == 0 ? "log" : "exp";
      benchD[ndim][iar].finish();
      benchE1[ndim][iar].finish();
      benchE2[ndim][iar].finish();
      fil << benchD[ndim][iar].str() << "\n";
      fil << benchE1[ndim][iar].str() << "\n";
      fil << benchE2[ndim][iar].str() << "\n";
      fil << "plt.figure("<<std::to_string(++ifig)<<")\n";
      fil << "plt.semilogx(aniso,np.array(" << benchE1[ndim][iar].name()<<")"
          << "/np.array(" << benchD[ndim][iar].name() << ")" <<",label='Eigen1')\n";
      fil << "plt.semilogx(aniso,np.array(" << benchE2[ndim][iar].name()<<")"
          << "/np.array(" << benchD[ndim][iar].name() << ")" <<",label='Eigen2')\n";
      fil << "plt.title('Benchmark "<<funname<<" dim "<< std::to_string(ndim+2)<< "')\n";
      fil << "plt.xlabel('anisotropic ratio')\n";
      fil << "plt.ylabel('op/s divided by DSYEVQs')\n";
      fil << "plt.legend()\n";
      fil << "plt.grid(True)\n";
      fil << "plt.savefig('bench_explogmet_"<<funname<<"_"<<std::to_string(ndim+2)<<".png')\n";
      fil << "\n\n";
    }
  }
  fil << "plt.show()\n";
  fil.close();

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