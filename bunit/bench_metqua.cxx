//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 

#include <boost/test/included/unit_test.hpp> 
#include "common_setup.hxx"

#include "../src/ho_constants.hxx"
#include "../src/utils/aux_misc.hxx"
#include "../src/quality/low_metqua.hxx"

#include <random>

using namespace Metris;


// P1 only and uses 1 eval with avg metric. This is current metqua implementation
// (05/2025)
template <typename MFT, int idim>
double metqua1(Mesh<MFT> &msh, int ientt, double difto);

// P1 only and uses tdim + 1 evals instead of avg metric
template <typename MFT, int idim>
double metqua2(Mesh<MFT> &msh, int ientt, double difto);

// Like metqua1 but using linearly averaged metric
template <typename MFT, int idim>
double metqua3(Mesh<MFT> &msh, int ientt, double difto);




BOOST_AUTO_TEST_CASE(bench_metqua) 
{

  // bool is whether straight
  std::vector<std::string> meshes = {
     METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k.meshb"
    ,METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k.meshb -t 2"
    ,METRIS_CASES_DIR "/unit/2D/square/iso.p1.100k.meshb"
    ,METRIS_CASES_DIR "/unit/2D/square/iso.p1.100k.meshb -t 2"
    };

  const int tarop = 1e6;

  // Accumulate into dummy to avoid optimizing out
  double dum = 0;

  for(auto meshname : meshes){

    CT_FOR0_INC(1,2,itest){
      
      std::string  metarg = itest == 1 ? " -anamet 1" : " ";
      cargHandler arg("-in " + meshname + metarg + " -opt-norm 2 -verb 0");
      MetrisRunner run(arg.c, arg.v);
      run.degElevate();
      using MFT = typename std::conditional<itest == 1, MetricFieldAnalytical, MetricFieldFE>::type;
      Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);

      std::string testname = itest == 1 ? "Analytical" : "FE";


      CT_FOR0_INC(1,3,gdim){if(gdim == msh.idim){
        CT_FOR0_INC(2,gdim,tdim){

          int nentt = msh.nentt(tdim);
          const int nloop = ceil(tarop / (double) nentt);

          // For slight mesh perturbations
          std::uniform_real_distribution<double> unif(0.0,1.0);
          std::default_random_engine rng(0);

          printf("\n-- Running mesh %s metric %s gdim %d tdim %d ideg %d nloop = %d\n",
                 meshname.c_str(), testname.c_str(), gdim, tdim, msh.curdeg, nloop);

          msh.met.setSpace(MetSpace::Exp);
          double tm0_exp = get_wall_time();
          for(int iloop = 0; iloop < nloop; iloop++){
            for(int ielem = 0; ielem < nentt; ielem++){
              dum += metqua<MFT,gdim,tdim>(msh, AsDeg::Pk, AsDeg::Pk, ielem, 1.0);
            }
          }// for iloop
          double tm1_exp = get_wall_time();

          int psm_exp = (int) ((((double)nentt) * nloop) / (tm1_exp - tm0_exp) / 1000);
          printf("-- Test end prod metqua = %dk/s\n",psm_exp);
          
          if(tdim == gdim && msh.curdeg == 1){

            msh.met.setSpace(MetSpace::Log);
            double t01 = get_wall_time();
            for(int iloop = 0; iloop < nloop; iloop++){
              for(int ielem = 0; ielem < nentt; ielem++){
                dum -= metqua1<MFT,gdim>(msh, ielem, 1);
              }
            }// for iloop
            double t11 = get_wall_time();


            msh.met.setSpace(MetSpace::Exp);
            double t02 = get_wall_time();
            for(int iloop = 0; iloop < nloop; iloop++){
              for(int ielem = 0; ielem < nentt; ielem++){
                dum += metqua2<MFT,gdim>(msh, ielem, 1);
              }
            }// for iloop
            double t12 = get_wall_time();


            // Accumulate into dummy to avoid optimizing out
            msh.met.setSpace(MetSpace::Exp);
            double t03 = get_wall_time();
            for(int iloop = 0; iloop < nloop; iloop++){
              for(int ielem = 0; ielem < nentt; ielem++){
                dum -= metqua3<MFT,gdim>(msh, ielem, 1);
              }
            }// for iloop
            double t13 = get_wall_time();

            int ps1 = (int) ((((double)nentt) * nloop) / (t11 - t01) / 1000);
            int ps2 = (int) ((((double)nentt) * nloop) / (t12 - t02) / 1000);
            int ps3 = (int) ((((double)nentt) * nloop) / (t13 - t03) / 1000);


              
            printf("          bench metqua1 = %dk/s\n",ps1);
            printf("          bench metqua2 = %dk/s ratio 2/1 = %e\n",ps2, (t12 - t02) / (t11 - t01));
            printf("          bench metqua3 = %dk/s ratio 3/1 = %e\n",ps3, (t13 - t03) / (t11 - t01));
          }
        
        }CT_FOR1(tdim);
      }}CT_FOR1(gdim);

    }CT_FOR1(itest);

  }// for meshname : meshes

  printf("\n\ndum = %e\n",dum);
}// end boost test case






// Same as metqua for P1 only and uses tdim + 1 evals instead of avg metric
template <typename MFT, int idim>
double metqua2(Mesh<MFT> &msh, int ientt, double difto){
  constexpr int tdim = idim;
  constexpr int gdim = idim;

  const intAr2 &ent2poi = msh.ent2poi(tdim);

  METRIS_ASSERT(!isdeadent(ientt, ent2poi));
  METRIS_ASSERT_MSG(msh.curdeg==1, "This test case is for Q1 and interior only.")
  METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);


  double bary[tdim+1];

  const int pnorm = msh.param->opt_pnorm;

  double qutet = 0; 
  // Performance impact should be zero
  constexpr auto quafun_xi = get_quafun_xi<MFT,gdim,tdim,QuaFun::Distortion,double>();

  
  // Value doesn't matter. (P1)
  for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0 / (tdim  + 1);
  METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);
  for(int iquad = 0; iquad < tdim + 1; iquad++){
    int ipoin = ent2poi(ientt, iquad);
    double qua0 = quafun_xi(msh,AsDeg::Pk,AsDeg::Pk,ent2poi[ientt],bary,msh.met[ipoin]);
    // You'd think this wouldn't be a bottleneck but it eats up 20% of optimization
    // time to run pow() here even if pnorm = 2 or 1.
    qua0 = abs(qua0 - difto);
    if(pnorm == 2){
      qua0 *= qua0;
    }else if(pnorm > 2){
      qua0 = pow(qua0, pnorm);
    }
    qutet += qua0 / (tdim + 1);
  }

  return qutet;
}
template double metqua2<MetricFieldAnalytical, 2>(Mesh<MetricFieldAnalytical> &msh, int ientt, double difto);
template double metqua2<MetricFieldAnalytical, 3>(Mesh<MetricFieldAnalytical> &msh, int ientt, double difto);
template double metqua2<MetricFieldFE, 2>(Mesh<MetricFieldFE> &msh, int ientt, double difto);
template double metqua2<MetricFieldFE, 3>(Mesh<MetricFieldFE> &msh, int ientt, double difto);


template <typename MFT, int idim>
double metqua1(Mesh<MFT> &msh, int ientt, double difto){
  constexpr int tdim = idim;
  constexpr int gdim = idim;

  const intAr2 &ent2poi = msh.ent2poi(tdim);

  METRIS_ASSERT(!isdeadent(ientt, ent2poi));
  METRIS_ASSERT_MSG(msh.curdeg==1, "This test case is for Q1 and interior only.")


  const int pnorm = msh.param->opt_pnorm;

  // Performance impact should be zero
  constexpr auto quafun_xi = get_quafun_xi<MFT,gdim,tdim,QuaFun::Distortion,double>();

  double bary[tdim+1];
  for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0 / (tdim  + 1);
  double qutet = quafun_xi(msh,AsDeg::Pk,AsDeg::Pk,ent2poi[ientt],bary,NULL);
  // You'd think this wouldn't be a bottleneck but it eats up 20% of optimization
  // time to run pow() here even if pnorm = 2 or 1.
  qutet = abs(qutet - difto);
  if(pnorm == 2){
    qutet *= qutet;
  }else if(pnorm > 2){
    qutet = pow(qutet, pnorm);
  }
  return qutet;
}
template double metqua1<MetricFieldAnalytical, 2>(Mesh<MetricFieldAnalytical> &msh, int ientt, double difto);
template double metqua1<MetricFieldAnalytical, 3>(Mesh<MetricFieldAnalytical> &msh, int ientt, double difto);
template double metqua1<MetricFieldFE, 2>(Mesh<MetricFieldFE> &msh, int ientt, double difto);
template double metqua1<MetricFieldFE, 3>(Mesh<MetricFieldFE> &msh, int ientt, double difto);



template <typename MFT, int idim>
double metqua3(Mesh<MFT> &msh, int ientt, double difto){
  constexpr int tdim = idim;
  constexpr int gdim = idim;

  const intAr2 &ent2poi = msh.ent2poi(tdim);

  METRIS_ASSERT(!isdeadent(ientt, ent2poi));
  METRIS_ASSERT(msh.met.getSpace() == MetSpace::Exp);
  METRIS_ASSERT_MSG(msh.curdeg==1, "This test case is for Q1 and interior only.")

  double bary[tdim+1];

  const int pnorm = msh.param->opt_pnorm;

  // Performance impact should be zero
  constexpr auto quafun_xi = get_quafun_xi<MFT,gdim,tdim,QuaFun::Distortion,double>();

  constexpr int nnmet = (idim*(idim+1))/2;
  double met[nnmet] = {};
  for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0 / (tdim  + 1);
  for(int ii = 0; ii < tdim + 1; ii++){
    int ipoin = ent2poi(ientt, ii);
    for(int jj = 0; jj < nnmet; jj++){
      met[jj] += msh.met(ipoin,jj) / (tdim + 1);
    }
  }
  double qutet = quafun_xi(msh,AsDeg::Pk,AsDeg::Pk,ent2poi[ientt],bary,met);
  // You'd think this wouldn't be a bottleneck but it eats up 20% of optimization
  // time to run pow() here even if pnorm = 2 or 1.
  qutet = abs(qutet - difto);
  if(pnorm == 2){
    qutet *= qutet;
  }else if(pnorm > 2){
    qutet = pow(qutet, pnorm);
  }
  return qutet;
}
template double metqua3<MetricFieldAnalytical, 2>(Mesh<MetricFieldAnalytical> &msh, int ientt, double difto);
template double metqua3<MetricFieldAnalytical, 3>(Mesh<MetricFieldAnalytical> &msh, int ientt, double difto);
template double metqua3<MetricFieldFE, 2>(Mesh<MetricFieldFE> &msh, int ientt, double difto);
template double metqua3<MetricFieldFE, 3>(Mesh<MetricFieldFE> &msh, int ientt, double difto);
