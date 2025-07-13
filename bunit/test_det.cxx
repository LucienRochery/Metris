//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 

#include <boost/test/included/unit_test.hpp> 
#include "common_setup.hxx"

#include <boost/timer/progress_display.hpp>
#include <random>
#include <fstream>

//#include "../src/utils/aux_utils.hxx"
#include "../src/low_geo.hxx"
#include "../src/linalg/det.hxx"
#include "../src/linalg/eigen.hxx"
#include "../src/linalg/explogmet.hxx"


#include <Eigen/Dense>

namespace Metris{

template<int ndim,typename ftype>
ftype detsym_Eigen_LDLT(const ftype *met);



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

  const double aniso_max = 2e7;
  const double aniso_mul = 10;


  // i = idim - 2, j = 0 ? avg : max;
  scriptArrayString errsN[2], errsNf4[2], errsELLT[2], errsEdet[2], errsELDLT[2], errsELDLTf4[2], errsL[2];

  scriptArrayString benchN[2], benchELLT[2], benchEdet[2], benchELDLTf4[2], benchELDLT[2], benchL[2], benchNf4[2];

  scriptArrayString anisos("aniso");
  for(double anisorat = 2; anisorat <= aniso_max + 1; anisorat *= aniso_mul){
    anisos += anisorat;
  }
  anisos.finish();

  for(int ii = 0; ii < 2; ii++){
    errsN[ii].setName("err_max_Naive"+std::to_string(ii+2));
    errsNf4[ii].setName("err_max_Naivef4"+std::to_string(ii+2));
    errsL[ii].setName("err_max_LAPACK"+std::to_string(ii+2));
    errsELLT[ii].setName("err_max_EIGENLLT"+std::to_string(ii+2));
    errsEdet[ii].setName("err_max_EIGENDET"+std::to_string(ii+2));
    errsELDLT[ii].setName("err_max_EigenLDLT"+std::to_string(ii+2));
    errsELDLTf4[ii].setName("err_max_EigenLDLTf4"+std::to_string(ii+2));

    benchN[ii].setName("time_Naive"+std::to_string(ii+2));
    benchNf4[ii].setName("time_Naivef4"+std::to_string(ii+2));
    benchL[ii].setName("time_LAPACK"+std::to_string(ii+2));
    benchELLT[ii].setName("time_EIGENLLT"+std::to_string(ii+2));
    benchEdet[ii].setName("time_EIGENDET"+std::to_string(ii+2));
    benchELDLT[ii].setName("time_EigenLDLT"+std::to_string(ii+2));
    benchELDLTf4[ii].setName("time_EigenLDLTf4"+std::to_string(ii+2));
  }

  dblAr2 met_samples[2];
  for(int idim = 2; idim <= 3; idim++){
    met_samples[idim-2].allocate(nsamp,(idim*(idim+1))/2);
    met_samples[idim-2].set_n(nsamp);
  }


#ifdef USE_MULTIPRECISION
  CT_FOR0_INC(2,3,ndim){
    constexpr int nnmet = (ndim*(ndim+1))/2;
    typedef Eigen::Matrix<double,ndim,ndim> MatrixN;
    typedef Eigen::Vector<double,ndim> VectorN;

    double eigval[ndim],eigvec[ndim][ndim];
    double eigva2[ndim],eigve2[ndim][ndim];
    double met[nnmet], met2[nnmet];

    printf("-- Tests for dim = %d \n",ndim);

    dblAr1 errLLT_series, aniso_series;

    for(double anisorat = 2; anisorat <= aniso_max + 1; anisorat *= aniso_mul){

      MinMaxAvg errN, err3, errL, errNf4, errD, errELLT, errE2, errELDLT, errELDLTf4;

      double nerro_aniso = 0;


      double err;

      for(int isamp = 0; isamp < nsamp; isamp++){
        nerro_aniso++;

        generate_metric<ndim>(anisorat, eigval, eigvec, unif, rng);

        // Generate metric from eigvecs, eigvals
        eig2met<ndim,double>(eigval,eigvec[0],met);
        for(int ii = 0; ii < nnmet; ii++) met_samples[ndim-2](isamp, ii) = met[ii];

        // Truth value
        float8 met8[nnmet];
        for(int ii = 0; ii < nnmet; ii++) met8[ii] = (float8) met[ii];
        float8 det0 = detsym<ndim, float8>(met8);

        // ----- Test detsym<float4> (naive computation)
        float4 met4[nnmet];
        for(int ii = 0; ii < nnmet; ii++) met4[ii] = (float4) met[ii];
        float4 det4 = detsym<ndim, float4>(met4);
        float8 err = abs( ((float8) det4) - det0 ) /det0;
        errNf4 += (double) err;

        // ----- Test detsym (naive computation)
        double det1 = detsym<ndim>(met);
        err = abs( ((float8) det1) - det0 ) /det0;
        errN += (double) err;

        #ifdef METRIS_USE_LAPACK
        // ----- Test detsym_LAPACK which uses LAPACK in 3D (same as previous in 2D)
        double detL = detsym_LAPACK<ndim>(met);
        err = abs( ((float8) detL) - det0 ) /det0;
        errL += (double) err;
        #endif

        // ----- Test detsym3 which does a kind of preconditioning 
        double det3 = detsym3<ndim>(met);
        err = abs( ((float8) det3) - det0 ) /det0;
        err3 += (double) err;

        // ----- Test eigen decomp which uses LAPACK in 3D (same as previous in 2D)
        geteigsym<ndim,double>(met,eigval,eigvec[0]);
        double detD = eigval[0];
        for(int ii = 1; ii < ndim; ii++) detD *= eigval[ii];
        err = abs( ((float8) detD) - det0 ) /det0;
        errD += (double) err;


        // ----- Test Eigen using a Cholesky decomposition
        double detE = detsym_Eigen_LLT<ndim>(met);
        err = abs( ((float8) detE) - det0 ) /det0;
        errELLT += (double) err;



        // ----- Test Eigen using determinant
        MatrixN met_Eigen;
        for(int ii = 0; ii < ndim; ii++) 
          for(int jj = 0; jj < ndim; jj++) 
            met_Eigen(ii, jj) = met[sym2idx(ii,jj)];
        detE = met_Eigen.determinant();
        err = abs( ((float8) detE) - det0 ) /det0;
        errE2 += (double) err;



        // ----- Test Eigen using LDLT
        double detE3 = detsym_Eigen_LDLT<ndim>(met);
        err = abs( ((float8) detE3) - det0 ) /det0;
        errELDLT += (double) err;

        // ----- Test Eigen using LDLT and float4
        float4 detE4 = detsym_Eigen_LDLT<ndim,float4>(met4);
        err = abs( ((float8) detE4) - det0 ) /det0;
        errELDLTf4 += (double) err;



      }// for isamp

      errsN[ndim-2] += errN.max();
      errsNf4[ndim-2] += errNf4.max();
      #ifdef METRIS_USE_LAPACK
      errsL[ndim-2] += errL.max();
      #endif
      errsELLT[ndim-2] += errELLT.max();
      errsEdet[ndim-2] += errE2.max();
      errsELDLT[ndim-2] += errELDLT.max();
      errsELDLTf4[ndim-2] += errELDLTf4.max();

      printf("  -- DONE with aniso ratio %5.1e dim %d\n",anisorat,ndim);
      printf("   - det error min = %e avg = %e max = %e (Naive)\n",errN.min(),errN.avg(),errN.max());
      printf("   - det error min = %e avg = %e max = %e (Naive f4)\n",errNf4.min(),errNf4.avg(),errNf4.max());
      #ifdef METRIS_USE_LAPACK
      printf("   - det error min = %e avg = %e max = %e (LAPACK)\n",errL.min(),errL.avg(),errL.max());
      #endif
      printf("   - det error min = %e avg = %e max = %e (Precond)\n",err3.min(),err3.avg(),err3.max());
      printf("   - det error min = %e avg = %e max = %e (DSYEVQ)\n",errD.min(),errD.avg(),errD.max());
      printf("   - det error min = %e avg = %e max = %e (Eigen LLT)\n",errELLT.min(),errELLT.avg(),errELLT.max());
      printf("   - det error min = %e avg = %e max = %e (Eigen det)\n",errE2.min(),errE2.avg(),errE2.max());
      printf("   - det error min = %e avg = %e max = %e (Eigen LDLT)\n",errELDLT.min(),errELDLT.avg(),errELDLT.max());
      printf("   - det error min = %e avg = %e max = %e (Eigen LDLT f4)\n",errELDLTf4.min(),errELDLTf4.avg(),errELDLTf4.max());

      //BOOST_TEST(errN <= tol);
      //BOOST_TEST(errL <= tol);
      //BOOST_TEST(errELLT <= tol);
      //BOOST_TEST(errELDLT <= tol);

      errLLT_series.stack(errELLT.avg());
      aniso_series.stack(anisorat);


    }// for anisorat

    // Linear regression on errLLT (what is used in Metris) 
    // check O(aniso^2) as empirically observed (note matrix conditioning is aniso^2)
    // and expected under 10^{-13} at aniso 10.
    int nseries = errLLT_series.get_n();
    for(int ii = 0; ii < nseries; ii++){
      errLLT_series[ii] = log(errLLT_series[ii]);
      aniso_series[ii] = log(aniso_series[ii]);
    }
    LinReg linreg(aniso_series,errLLT_series);
    BOOST_TEST(linreg.slope <= 2.5);

    double err10 = linreg.origin + linreg.slope*log(10);
    BOOST_TEST(err10 < log(1.0e-13));

  }CT_FOR1(ndim);


#else
  printf("## Unit test disabled as USE_MULTIPRECISION not defined\n");
#endif



  CT_FOR0_INC(2,3,ndim){
    constexpr int nnmet = (ndim*(ndim+1))/2;
    typedef Eigen::Matrix<double,ndim,ndim> MatrixN;
    typedef Eigen::Vector<double,ndim> VectorN;

    double eigval[ndim],eigvec[ndim][ndim];
    double eigva2[ndim],eigve2[ndim][ndim];
    double met[nnmet], met2[nnmet];
    #ifdef USE_MULTIPRECISION
    float4 met4[nnmet];
    #endif

    for(double anisorat = 2; anisorat <= aniso_max + 1; anisorat *= aniso_mul){
      // ------------------------------------------------------------ Benchmark
      double dum_tot = 0;

      #ifdef NDEBUG
      double t0_N = get_wall_time();
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met[ii] = met_samples[ndim-2](isamp, ii);
        double det = detsym<ndim>(met);
        dum_tot += det;
      }
      double t1_N = get_wall_time();

      #ifdef USE_MULTIPRECISION
      double t0_Nf4 = get_wall_time();
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met4[ii] = (float4) met_samples[ndim-2](isamp, ii);
        float4 det = detsym<ndim,float4>(met4);
        dum_tot += (double) det;
      }
      double t1_Nf4 = get_wall_time();
      #endif


      #ifdef METRIS_USE_LAPACK
      double t0_L = get_wall_time();
      int rwork[10];
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met[ii] = met_samples[ndim-2](isamp, ii);
        double det = detsym_LAPACK<ndim>(met);
        dum_tot += det;
      }
      double t1_L = get_wall_time();
      #endif


      double t0_ELLT = get_wall_time();
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met[ii] = met_samples[ndim-2](isamp, ii);
        double det = detsym_Eigen_LLT<ndim>(met);
        dum_tot += det;
      }
      double t1_ELLT = get_wall_time();

      double t0_Edet = get_wall_time();
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met[ii] = met_samples[ndim-2](isamp, ii);
        
        MatrixN met_Eigen;
        for(int ii = 0; ii < ndim; ii++) 
          for(int jj = 0; jj < ndim; jj++) 
            met_Eigen(ii, jj) = met[sym2idx(ii,jj)];
        double det = met_Eigen.determinant();
        dum_tot += det;
      }
      double t1_Edet = get_wall_time();

      double t0_ELDLT= get_wall_time();
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met[ii] = met_samples[ndim-2](isamp, ii);
        double det = detsym_Eigen_LDLT<ndim>(met);
        dum_tot+= det;
      }
      double t1_ELDLT= get_wall_time();

      #ifdef USE_MULTIPRECISION
      double t0_ELDLTf4 = get_wall_time();
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met4[ii] = (float4) met_samples[ndim-2](isamp, ii);
        float4 det = detsym_Eigen_LDLT<ndim,float4>(met4);
        dum_tot+= (double) det;
      }
      double t1_ELDLTf4 = get_wall_time();
      #endif


      printf("  -- DONE bench dim %1d aniso %5.1e dum %5.0e\n",ndim,anisorat,
             dum_tot);
      printf("  -- DONE benchmarks  Naive    time : %8.2e = %dk op/s\n",
                  t1_N-t0_N,(int)(nsamp/(t1_N-t0_N))/1000);
      #ifdef USE_MULTIPRECISION
      printf("                     Naivef4   time : %8.2e = %dk op/s, fac = %4.2fx\n",
                  t1_Nf4-t0_Nf4,(int)(nsamp/(t1_Nf4-t0_Nf4)/1000),
                  (t1_Nf4-t0_Nf4)/(t1_N-t0_N));
      #endif
      #ifdef METRIS_USE_LAPACK
      printf("                     LAPACK    time : %8.2e = %dk op/s, fac = %4.2fx\n",
                  t1_L-t0_L,(int)(nsamp/(t1_L-t0_L)/1000),
                  (t1_L-t0_L)/(t1_N-t0_N));
      #endif
      printf("                      ELLT     time : %8.2e = %dk op/s, fac = %4.2fx\n",
                  t1_ELLT-t0_ELLT,(int)(nsamp/(t1_ELLT-t0_ELLT)/1000),
                  (t1_ELLT-t0_ELLT)/(t1_N-t0_N));
      printf("                      Edet     time : %8.2e = %dk op/s, fac = %4.2fx\n",
                  t1_Edet-t0_Edet,(int)(nsamp/(t1_Edet-t0_Edet)/1000),
                  (t1_Edet-t0_Edet)/(t1_N-t0_N));
      printf("                      ELDLT   time : %8.2e = %dk op/s, fac = %4.2fx\n",
                  t1_ELDLT-t0_ELDLT,(int)(nsamp/(t1_ELDLT-t0_ELDLT)/1000),
                  (t1_ELDLT-t0_ELDLT)/(t1_N-t0_N));
      #ifdef USE_MULTIPRECISION
      printf("                      ELDLTf4 time : %8.2e = %dk op/s, fac = %4.2fx\n",
                  t1_ELDLTf4-t0_ELDLTf4,(int)(nsamp/(t1_ELDLTf4-t0_ELDLTf4)/1000),
                  (t1_ELDLTf4-t0_ELDLTf4)/(t1_N-t0_N));
      #endif

      benchN[ndim-2] += nsamp/(t1_N - t0_N);
      #ifdef USE_MULTIPRECISION
      benchNf4[ndim-2] += nsamp/(t1_Nf4 - t0_Nf4);
      #endif
      #ifdef METRIS_USE_LAPACK
      benchL[ndim-2] += nsamp/(t1_L - t0_L);
      #endif
      benchELLT[ndim-2] += nsamp/(t1_ELLT - t0_ELLT);
      benchEdet[ndim-2] += nsamp/(t1_Edet - t0_Edet);
      benchELDLT[ndim-2] += nsamp/(t1_ELDLT- t0_ELDLT);
      #ifdef USE_MULTIPRECISION
      benchELDLTf4[ndim-2] += nsamp/(t1_ELDLTf4- t0_ELDLTf4);
      #endif

      #endif

    }// for anisorat

  }CT_FOR1(ndim);


  std::ofstream fil;
  std::string fname = "test_det.py";
  fil.open(fname.c_str(), std::ios::out);
  fil << "import matplotlib.pyplot as plt\n\n";
  fil << anisos.str() << "\n\n";
  int ifig = 0;
  for(int ndim = 0; ndim < 2; ndim++){
    errsN[ndim].finish();
    #ifdef METRIS_USE_LAPACK
    errsL[ndim].finish();
    #endif
    errsELLT[ndim].finish();
    errsEdet[ndim].finish();
    errsELDLT[ndim].finish();
    fil << errsN[ndim].str() << "\n";
    #ifdef METRIS_USE_LAPACK
    fil << errsL[ndim].str() << "\n";
    #endif
    fil << errsELLT[ndim].str() << "\n";
    fil << errsEdet[ndim].str() << "\n";
    fil << errsELDLT[ndim].str() << "\n";
    fil << "plt.figure("<<std::to_string(++ifig)<<")\n";
    fil << "plt.loglog(aniso," << errsN[ndim].name()<<",'--o',label='Naive',)\n";
    #ifdef METRIS_USE_LAPACK
    fil << "plt.loglog(aniso," << errsL[ndim].name()<<",'--x',label='LAPACK')\n";
    #endif
    fil << "plt.loglog(aniso," << errsELLT[ndim].name()<<",'--+',label='Eigen (LLT)')\n";
    fil << "plt.loglog(aniso," << errsELDLT[ndim].name()<<",'--o',label='Eigen (LDLT)')\n";
    fil << "plt.loglog(aniso," << errsEdet[ndim].name()<<",'--^',label='Eigen (det)')\n";
    fil << "plt.title('Determinant relative error (max), dim "
        << std::to_string(ndim+2)<< "')\n";
    fil << "plt.xlabel('anisotropic ratio')\n";
    fil << "plt.ylabel('relative error')\n";
    fil << "plt.legend()\n";
    fil << "plt.grid(True)\n";
    fil << "plt.savefig('errs_det_"<<std::to_string(ndim+2)<<".png')\n";
    fil << "\n\n";
  }
  fil << "plt.show()\n";
  fil.close();

  fname = "bench_det.py";
  fil.open(fname.c_str(), std::ios::out);
  fil << "import matplotlib.pyplot as plt\n";
  fil << "import numpy as np\n";
  fil << "\n";

  fil << anisos.str() << "\n\n";
  ifig = 0;
  for(int ndim = 0; ndim < 2; ndim++){

    benchN[ndim].finish();
    #ifdef METRIS_USE_LAPACK
    benchL[ndim].finish();
    #endif
    benchELLT[ndim].finish();
    benchEdet[ndim].finish();
    benchELDLT[ndim].finish();
    fil << benchN[ndim].str() << "\n";
    #ifdef METRIS_USE_LAPACK
    fil << benchL[ndim].str() << "\n";
    #endif
    fil << benchELLT[ndim].str() << "\n";
    fil << benchEdet[ndim].str() << "\n";
    fil << benchELDLT[ndim].str() << "\n";
    fil << "plt.figure("<<std::to_string(++ifig)<<")\n";
    #ifdef METRIS_USE_LAPACK
    fil << "plt.loglog(aniso,np.array(" << benchL[ndim].name()<<")"
        << "/np.array(" << benchN[ndim].name() << ")" <<",label='LAPACK')\n";
    #endif
    fil << "plt.loglog(aniso,np.array(" << benchELLT[ndim].name()<<")"
        << "/np.array(" << benchN[ndim].name() << ")" <<",label='Eigen (LLT)')\n";
    fil << "plt.loglog(aniso,np.array(" << benchELDLT[ndim].name()<<")"
        << "/np.array(" << benchN[ndim].name() << ")" <<",label='Eigen (LDLT)')\n";
    fil << "plt.loglog(aniso,np.array(" << benchEdet[ndim].name()<<")"
        << "/np.array(" << benchN[ndim].name() << ")" <<",label='Eigen (det)')\n";
    //fil << "plt.loglog(aniso," << benchN[ndim].name()<<",label='DSYEVQ')\n";
    //fil << "plt.loglog(aniso," << benchL[ndim].name()<<",label='LAPACK')\n";
    //fil << "plt.loglog(aniso," << benchE[ndim].name()<<",label='Eigen')\n";
    fil << "plt.title('Benchmark dim "<< std::to_string(ndim+2)<< "')\n";
    fil << "plt.xlabel('anisotropic ratio')\n";
    fil << "plt.ylabel('op/s divided by Naives')\n";
    fil << "plt.legend()\n";
    fil << "plt.grid(True)\n";
    fil << "plt.savefig('bench_det_"<<std::to_string(ndim+2)<<".png')\n";
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




template<int ndim, typename ftype>
ftype detsym_Eigen_LDLT(const ftype *met){

  typedef Eigen::Matrix<ftype,ndim,ndim> MatrixN;
  typedef Eigen::Vector<ftype,ndim> VectorN;

  MatrixN met_Eigen;
  for(int ii = 0; ii < ndim; ii++) 
    for(int jj = 0; jj <= ii; jj++) 
      met_Eigen(ii, jj) = met[sym2idx(ii,jj)];

  Eigen::LDLT<MatrixN> ldlt(met_Eigen.template selfadjointView<Eigen::Lower>());
  if(ldlt.info() != Eigen::Success) METRIS_THROW(GeomExcept());

  VectorN Dv = ldlt.vectorD();

  ftype detE3 = Dv(0);
  for(int ii = 1; ii < ndim; ii++) detE3 *= Dv(ii);

  return detE3;
}

template  double detsym_Eigen_LDLT<2>(const double *met);
template  double detsym_Eigen_LDLT<3>(const double *met);
#ifdef USE_MULTIPRECISION
template  float4 detsym_Eigen_LDLT<2>(const float4 *met);
template  float4 detsym_Eigen_LDLT<3>(const float4 *met);
#endif


}