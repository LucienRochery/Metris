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
#include "../src/linalg/eigen.hxx"
#include "../src/linalg/det.hxx"
#include "../src/linalg/invmat.hxx"


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

BOOST_AUTO_TEST_CASE(test_invspd) 
{//METRIS_MAX_DEG

  const int nsamp = 1e5;
  const double tol = 1.0e-11;
  const double exptol = 1.0e-16;

  std::uniform_real_distribution<double> unif(1.0e-16,1.0);
  std::default_random_engine rng(0);

  double aniso_max = 2e6;
  double aniso_mul = 10;


  dblAr2 met_samples[2],invmet_samples[2];
  for(int idim = 2; idim <= 3; idim++){
    met_samples[idim-2].allocate(nsamp,(idim*(idim+1))/2);
    met_samples[idim-2].set_n(nsamp);

    invmet_samples[idim-2].allocate(nsamp,(idim*(idim+1))/2);
    invmet_samples[idim-2].set_n(nsamp);
  }


  // i = idim - 2, j = 0 ? log : exp;
  scriptArrayString spd_errsD[2], spd_errsE[2], spd_errsL[2];
  scriptArrayString spd_benchD[2], spd_benchE[2], spd_benchL[2];

  scriptArrayString inv_errsL[2], inv_errsPP[2], inv_errsFP[2], inv_errsN[2];
  scriptArrayString inv_benchL[2], inv_benchPP[2], inv_benchFP[2], inv_benchN[2];

  scriptArrayString anisos("aniso");
  for(double anisorat = 2; anisorat <= aniso_max + 1; anisorat *= aniso_mul){
    anisos += anisorat;
  }
  anisos.finish();

  for(int ii = 0; ii < 2; ii++){
    spd_errsD[ii].setName("err_max_spd_DSYEVQ_"+std::to_string(ii+2));
    spd_errsE[ii].setName("err_max_spd_EIGEN_"+std::to_string(ii+2));
    spd_errsL[ii].setName("err_max_spd_LAPACK_"+std::to_string(ii+2));

    spd_benchD[ii].setName("time_spd_DSYEVQ_"+std::to_string(ii+2));
    spd_benchE[ii].setName("time_spd_EIGEN_"+std::to_string(ii+2));
    spd_benchL[ii].setName("time_spd_LAPACK_"+std::to_string(ii+2));

    inv_errsN[ii].setName("err_max_inv_NAIVE_"+std::to_string(ii+2));
    inv_errsL[ii].setName("err_max_inv_LAPACK_"+std::to_string(ii+2));
    inv_errsFP[ii].setName("err_max_inv_EIGENFP_"+std::to_string(ii+2));
    inv_errsPP[ii].setName("err_max_inv_EIGENPP_"+std::to_string(ii+2));

    inv_benchN[ii].setName("time_inv_NAIVE_"+std::to_string(ii+2));
    inv_benchL[ii].setName("time_inv_LAPACK_"+std::to_string(ii+2));
    inv_benchFP[ii].setName("time_inv_EIGENFP_"+std::to_string(ii+2));
    inv_benchPP[ii].setName("time_inv_EIGENPP_"+std::to_string(ii+2));
  }


  #ifdef USE_MULTIPRECISION
    CT_FOR0_INC(2,3,ndim){
      constexpr int nnmet = (ndim*(ndim+1))/2;
      typedef Eigen::Matrix<double,ndim,ndim> MatrixN;
      typedef Eigen::Vector<double,ndim> VectorN;

      double eigval[ndim],eigvec[ndim][ndim];
      double eigva2[ndim],eigve2[ndim][ndim];
      double met[nnmet], met2[nnmet], invmet[nnmet];
      int ierro;
      MatrixN met_Eigen, met2_Eigen;


      for(double anisorat = 2; anisorat < aniso_max + 1; anisorat *= aniso_mul){

        MinMaxAvg spd_errD, spd_errL, spd_errE;

        for(int isamp = 0; isamp < nsamp; isamp++){

          generate_metric<ndim>(anisorat, eigval, eigvec, MetSpace::Exp, unif, rng);

          // Generate metric from eigvecs, eigvals
          eig2met<ndim,double>(eigval,eigvec[0],met);
          for(int ii = 0; ii < nnmet; ii++) met_samples[ndim-2](isamp, ii) = met[ii];

          // Get metric inverse
          double invmet[nnmet];
          float8 f8_invmet[nnmet];
          float8 f8_eigval[ndim], f8_eigvec[ndim*ndim];

          for(int ii = 0; ii < ndim; ii++) f8_eigval[ii] = 1.0 / (float8) eigval[ii];
          for(int ii = 0; ii < ndim; ii++) 
            for(int jj = 0; jj < ndim; jj++) 
              f8_eigvec[ii*ndim+jj] = (float8) eigvec[ii][jj];
          eig2met<ndim,float8>(f8_eigval, f8_eigvec, f8_invmet);
          for(int ii = 0; ii < nnmet; ii++){
            invmet[ii] = (double) f8_invmet[ii];
            invmet_samples[ndim-2](isamp, ii) = invmet[ii];
          }


          // --- Using LAPACK
          #ifdef USE_LAPACK
          for(int ii = 0; ii < nnmet; ii++) met2[ii] = met[ii];
          invspd_LAPACK(ndim, met2);
          spd_errL += sqrt(geterrl2<nnmet>(invmet,met2) / getnrml2<nnmet>(invmet));
          #endif

          // --- Using Eigen's LL^T
          ierro = invspd_Eigen<ndim,double>(met,met2);
          METRIS_ENFORCE_MSG(ierro == 0,"invspd_Eigen failed");
          spd_errE += sqrt(geterrl2<nnmet>(invmet,met2) / getnrml2<nnmet>(invmet));

          // --- Using eigendecomposition (DSYEVQ)
          double eigva2[ndim], eigve2[ndim*ndim];
          geteigsym<ndim>(met, eigva2, eigve2);
          for(int ii = 0; ii < ndim; ii++) eigva2[ii] = 1.0 / eigva2[ii];
          eig2met<ndim>(eigva2, eigve2, met2);
          spd_errD += sqrt(geterrl2<nnmet>(invmet,met2) / getnrml2<nnmet>(invmet));

        }// for isamp

        spd_errsD[ndim-2] += spd_errD.max();
        spd_errsE[ndim-2] += spd_errE.max();
        #ifdef USE_LAPACK
        spd_errsL[ndim-2] += spd_errL.max();
        #endif


        printf("  -- DONE with aniso ratio %5.1e dim %d\n",anisorat,ndim);
        #ifdef USE_LAPACK
        printf("   - error invspd min = %e avg = %e max = %e (LAPACK)\n",spd_errL.min(),spd_errL.avg(),spd_errL.max());
        #endif
        printf("   - error invspd min = %e avg = %e max = %e (Eigen LLT)\n",spd_errE.min(),spd_errE.avg(),spd_errE.max());
        printf("   - error invspd min = %e avg = %e max = %e (DSYEVQ)\n",spd_errD.min(),spd_errD.avg(),spd_errD.max());

      }// for anisorat

      printf("-- DONE spd tests\n\n");

      double mat[ndim*ndim], mat2[ndim*ndim], imat[ndim*ndim];
      for(double anisorat = 2; anisorat < aniso_max + 1; anisorat *= aniso_mul){

        MinMaxAvg inv_errFP, inv_errL, inv_errPP, inv_errN;

        for(int isamp = 0; isamp < nsamp; isamp++){

          for(int ii = 0; ii < ndim; ii++){
            for(int jj = 0; jj < ndim; jj++){
              mat[ndim*ii+jj] = met_samples[ndim-2](isamp, sym2idx(ii,jj));
              imat[ndim*ii+jj] = invmet_samples[ndim-2](isamp, sym2idx(ii,jj));
            }
          }


          // --- Using LAPACK
          #ifdef USE_LAPACK
          for(int ii = 0; ii < ndim*ndim; ii++) mat2[ii] = mat[ii];
          int ierroL = invmat_LAPACK(ndim, mat2);
          //METRIS_ENFORCE_MSG(ierro == 0, "invmat_LAPACK failed");
          //BOOST_CHECK(ierroL == 0);
          if(ierroL == 0)
            inv_errL += sqrt(geterrl2<nnmet>(imat,mat2) / getnrml2<nnmet>(imat));
          #endif

          // --- Using Naive
          if constexpr(ndim == 2){
            for(int ii = 0; ii < ndim*ndim; ii++) mat2[ii] = mat[ii];
            invmat_naive<ndim>(mat2);
            inv_errN += sqrt(geterrl2<nnmet>(imat,mat2) / getnrml2<nnmet>(imat));
          }

          // --- Using Eigen's PartialPivLU
          int ierroPP = invmat_EigenLUPP<ndim,double>(mat,mat2);
          //METRIS_ENFORCE_MSG(ierroPP == 0, "invmat_EigenLUPP failed");
          //BOOST_CHECK(ierroPP == 0);
          if(ierroPP == 0)
            inv_errPP += sqrt(geterrl2<nnmet>(imat,mat2) / getnrml2<nnmet>(imat));

          // --- Using Eigen's FullPivLU
          int ierroFP = invmat_EigenLUFP<ndim,double>(mat,mat2);
          //METRIS_ENFORCE_MSG(ierroFP == 0, "invmat_EigenLUFP failed");
          //BOOST_CHECK(ierroFP == 0);
          if(ierroFP == 0)
            inv_errFP += sqrt(geterrl2<nnmet>(imat,mat2) / getnrml2<nnmet>(imat));

        }// for isamp

        inv_errsFP[ndim-2] += inv_errFP.max();
        inv_errsPP[ndim-2] += inv_errPP.max();
        #ifdef USE_LAPACK
        inv_errsL[ndim-2] += inv_errL.max();
        #endif
        if(ndim == 2)
          inv_errsN[ndim-2] += inv_errN.max();


        printf("  -- DONE with aniso ratio %5.1e dim %d\n",anisorat,ndim);
        #ifdef USE_LAPACK
        printf("   - error invmat min = %e avg = %e max = %e (LAPACK)\n",inv_errL.min(),inv_errL.avg(),inv_errL.max());
        #endif
        printf("   - error invmat min = %e avg = %e max = %e (Eigen PP)\n",inv_errPP.min(),inv_errPP.avg(),inv_errPP.max());
        printf("   - error invmat min = %e avg = %e max = %e (Eigen FP)\n",inv_errFP.min(),inv_errFP.avg(),inv_errFP.max());
        if(ndim == 2) 
          printf("   - error invmat min = %e avg = %e max = %e (Naive)\n",inv_errN.min(),inv_errN.avg(),inv_errN.max());


      }// for anisorat
      printf("-- DONE invmat tests\n\n");

   
    }CT_FOR1(ndim);
  #else
    printf("## invmat tests skipped as USE_MULTIPRECISION disabled\n");
  #endif

  // ---------------------------------------------------------------------
  // ----------------------------------------------- Benchmark SPD
  #ifdef NDEBUG
  printf("\n\n");
  CT_FOR0_INC(2,3,ndim){
    constexpr int nnmet = (ndim*(ndim+1))/2;
    double met[nnmet], met2[nnmet], met3[nnmet];

    double dumtot = 0;
    for(double anisorat = 2; anisorat <= aniso_max + 1; anisorat *= aniso_mul){

      // --- Using LAPACK
      #ifdef USE_LAPACK
      double t0_L = get_wall_time();
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met2[ii] = met_samples[ndim-2](isamp, ii);
        for(int ii = 0; ii < nnmet; ii++) met2[ii] = met[ii];
        invspd_LAPACK(ndim, met2);
        dumtot += met2[0];
      }
      double t1_L = get_wall_time();
      #endif

      // --- Using Eigen's LL^T
      double t0_E = get_wall_time();
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met2[ii] = met_samples[ndim-2](isamp, ii);
        invspd_Eigen<ndim,double>(met,met2);
        dumtot += met2[0];
      }
      double t1_E = get_wall_time();

      // --- Using eigendecomposition (DSYEVQ)
      double t0_D = get_wall_time();
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < nnmet; ii++) met2[ii] = met_samples[ndim-2](isamp, ii);
        double eigva2[ndim], eigve2[ndim*ndim];
        geteigsym<ndim>(met, eigva2, eigve2);
        for(int ii = 0; ii < ndim; ii++) eigva2[ii] = 1.0 / eigva2[ii];
        eig2met<ndim>(eigva2, eigve2, met2);
        dumtot += met2[0];
      }
      double t1_D = get_wall_time();


      #ifdef USE_LAPACK
      printf("  -- DONE bench dim %1d aniso %5.1e LAPACK time : %8.2es = %dk op/s\n",
                  ndim,anisorat,
                  t1_L-t0_L,(int)(nsamp/(t1_L-t0_L)/1000));
      printf("                                     Eigen time : %8.2es = %dk op/s, fac = %4.2fx\n",
                  t1_E-t0_E,(int)(nsamp/(t1_E-t0_E)/1000),
                  (t1_E-t0_E)/(t1_L-t0_L));
      printf("                                    DSYEVQ time : %8.2es = %dk op/s, fac = %4.2fx\n",
                  t1_D-t0_D,(int)(nsamp/(t1_D-t0_D)/1000),
                  (t1_D-t0_D)/(t1_L-t0_L));
      spd_benchL[ndim-2] += nsamp/(t1_L - t0_L);
      #else
      printf("  -- DONE bench dim %1d aniso %5.1e Eigen time : %8.2es = %dk op/s\n",
                  ndim,anisorat,
                  t1_E-t0_E,(int)(nsamp/(t1_E-t0_E)/1000));
      printf("                                    DSYEVQ time : %8.2es = %dk op/s\n",
                  t1_D-t0_D,(int)(nsamp/(t1_D-t0_D)/1000));
      #endif
      spd_benchE[ndim-2] += nsamp/(t1_E - t0_E);
      spd_benchD[ndim-2] += nsamp/(t1_D - t0_D);

    }// for anisorat
    printf("dum = %e\n",dumtot);
  }CT_FOR1(ndim);


  // ---------------------------------------------------------------------
  // ----------------------------------------------- Benchmark invmat
  printf("\n\n");
  CT_FOR0_INC(2,3,ndim){
    constexpr int nnmet = (ndim*(ndim+1))/2;
    double mat[nnmet],mat2[nnmet];

    double dumtot = 0;
    for(double anisorat = 2; anisorat <= aniso_max + 1; anisorat *= aniso_mul){

      #ifdef USE_LAPACK
      // --- Using LAPACK
      double t0_L = get_wall_time();
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < ndim; ii++)
          for(int jj = 0; jj < ndim; jj++)
            mat[ndim*ii+jj] = met_samples[ndim-2](isamp, sym2idx(ii,jj));
        invmat_LAPACK(ndim, mat);
        dumtot += mat2[0];
      }
      double t1_L = get_wall_time();
      #endif

      // --- Using Eigen's LL^T
      double t0_PP = get_wall_time();
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < ndim; ii++)
          for(int jj = 0; jj < ndim; jj++)
            mat[ndim*ii+jj] = met_samples[ndim-2](isamp, sym2idx(ii,jj));
        invmat_EigenLUPP<ndim,double>(mat,mat2);
        dumtot += mat2[0];
      }
      double t1_PP = get_wall_time();

      // --- Using eigendecomposition (DSYEVQ)
      double t0_FP = get_wall_time();
      for(int isamp = 0; isamp < nsamp; isamp++){
        for(int ii = 0; ii < ndim; ii++)
          for(int jj = 0; jj < ndim; jj++)
            mat[ndim*ii+jj] = met_samples[ndim-2](isamp, sym2idx(ii,jj));
        invmat_EigenLUFP<ndim,double>(mat,mat2);
        dumtot += mat2[0];
      }
      double t1_FP = get_wall_time();


      #ifdef USE_LAPACK
      printf("  -- DONE invmat bench dim %1d aniso %5.1e LAPACK time : %8.2es = %dk op/s\n",
                  ndim,anisorat,
                  t1_L-t0_L,(int)(nsamp/(t1_L-t0_L)/1000));
      printf("                                               FP time : %8.2es = %dk op/s, fac = %4.2fx\n",
                  t1_FP-t0_FP,(int)(nsamp/(t1_FP-t0_FP)/1000),
                  (t1_FP-t0_FP)/(t1_L-t0_L));
      printf("                                               PP time : %8.2es = %dk op/s, fac = %4.2fx\n",
                  t1_PP-t0_PP,(int)(nsamp/(t1_PP-t0_PP)/1000),
                  (t1_PP-t0_PP)/(t1_L-t0_L));
      inv_benchL[ndim-2]  += nsamp/(t1_L - t0_L);
      #else
      printf("  -- DONE invmat bench dim %1d aniso %5.1e     FP time : %8.2es = %dk op/s\n",
                  ndim,anisorat,
                  t1_FP-t0_FP,(int)(nsamp/(t1_FP-t0_FP)/1000));
      printf("                                               PP time : %8.2es = %dk op/s\n",
                  t1_PP-t0_PP,(int)(nsamp/(t1_PP-t0_PP)/1000));
      #endif
      inv_benchFP[ndim-2] += nsamp/(t1_FP - t0_FP);
      inv_benchPP[ndim-2] += nsamp/(t1_PP - t0_PP);

    }// for anisorat
    printf("dum = %e\n",dumtot);
  }CT_FOR1(ndim);
  #endif

  std::ofstream fil;
  std::string fname = "test_explogmet.py";
  fil.open(fname.c_str(), std::ios::out);
  fil << "import matplotlib.pyplot as plt\n\n";
  fil << anisos.str() << "\n\n";
  int ifig = 0;
  for(int ndim = 0; ndim < 2; ndim++){
    spd_errsD[ndim].finish();
    spd_errsE[ndim].finish();
    #ifdef USE_LAPACK
    spd_errsL[ndim].finish();
    #endif
    fil << spd_errsD[ndim].str() << "\n";
    fil << spd_errsE[ndim].str() << "\n";
    #ifdef USE_LAPACK
    fil << spd_errsL[ndim].str() << "\n";
    #endif
    fil << "plt.figure("<<std::to_string(++ifig)<<")\n";
    fil << "plt.loglog(aniso," << spd_errsD[ndim].name()<<",label='DSYEVQ')\n";
    fil << "plt.loglog(aniso," << spd_errsE[ndim].name()<<",label='Eigen')\n";
    #ifdef USE_LAPACK
    fil << "plt.loglog(aniso," << spd_errsL[ndim].name()<<",label='LAPACK')\n";
    #endif
    fil << "plt.title('Output relative error (max), dim "
        << std::to_string(ndim+2)<< "')\n";
    fil << "plt.xlabel('anisotropic ratio')\n";
    fil << "plt.ylabel('relative error')\n";
    fil << "plt.legend()\n";
    fil << "plt.grid(True)\n";
    fil << "plt.savefig('errs_explogmet_"<<std::to_string(ndim+2)<<".png')\n";
    fil << "\n\n";
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
    spd_benchD[ndim].finish();
    spd_benchE[ndim].finish();
    #ifdef USE_LAPACK
    spd_benchL[ndim].finish();
    #endif
    fil << spd_benchD[ndim].str() << "\n";
    fil << spd_benchE[ndim].str() << "\n";
    #ifdef USE_LAPACK
    fil << spd_benchL[ndim].str() << "\n";
    #endif
    fil << "plt.figure("<<std::to_string(++ifig)<<")\n";
    fil << "plt.semilogx(aniso,np.array(" << spd_benchE[ndim].name()<<")"
        << "/np.array(" << spd_benchD[ndim].name() << ")" <<",label='Eigen')\n";
    #ifdef USE_LAPACK
    fil << "plt.semilogx(aniso,np.array(" << spd_benchL[ndim].name()<<")"
        << "/np.array(" << spd_benchD[ndim].name() << ")" <<",label='LAPACK')\n";
    #endif
    fil << "plt.title('Benchmark  dim "<< std::to_string(ndim+2)<< "')\n";
    fil << "plt.xlabel('anisotropic ratio')\n";
    fil << "plt.ylabel('op/s divided by DSYEVQs')\n";
    fil << "plt.legend()\n";
    fil << "plt.grid(True)\n";
    fil << "plt.savefig('bench_explogmet_"<<std::to_string(ndim+2)<<".png')\n";
    fil << "\n\n";
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