//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE My Test 

#include <boost/test/included/unit_test.hpp> 
#include <bunit/common_setup.hxx>

#include "../src/ho_constants.hxx"
#include "../src/utils/aux_misc.hxx"
#include "../src/quality/low_metqua.hxx"

#include <random>

namespace utf = boost::unit_test;

using namespace Metris;

typedef MetricFieldAnalytical MFT;


// P1 only and uses 1 eval with avg metric. This is current metqua implementation
// (05/2025)
template <int idim>
double metqua1(Mesh<MFT> &msh, int ientt, double difto);

// P1 only and uses tdim + 1 evals instead of avg metric
template <int idim>
double metqua2(Mesh<MFT> &msh, int ientt, double difto);

// Like metqua1 but using linearly averaged metric
template <int idim>
double metqua3(Mesh<MFT> &msh, int ientt, double difto);




BOOST_AUTO_TEST_CASE(bench_metqua) 
{

  // bool is whether straight
  std::vector<std::string> meshes = {
    METRIS_CASES_DIR "/unit/3D/cube/p1.2k.meshb"
    METRIS_CASES_DIR "/unit/2D/square/p1.100k.meshb"
    };

  const int tarop = 1e6;

  for(auto meshname : meshes){

    cargHandler arg("-in " + meshname + " -anamet 1 -opt-norm 2 -verb 0");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);



    int idim = msh.idim;
    METRIS_ASSERT(idim == msh.get_tdim());

    int nentt = msh.nentt(idim);
    const int nloop = ceil(tarop / (double) nentt);
    printf("-- Running mesh %s nloop = %d\n",meshname.c_str(), nloop);

    // For slight mesh perturbations
    std::uniform_real_distribution<double> unif(0.0,1.0);
    std::default_random_engine rng(0);

    // Accumulate into dummy to avoid optimizing out
    msh.met.setSpace(MetSpace::Log);
    double dum1 = 0;
    double t01 = get_wall_time();
    for(int iloop = 0; iloop < nloop; iloop++){
      for(int ielem = 0; ielem < nentt; ielem++){
        if(idim == 2){
          dum1 += metqua1<2>(msh, ielem, 1);
        }else{
          dum1 += metqua1<3>(msh, ielem, 1);
        }
      }
    }// for iloop
    double t11 = get_wall_time();


    // Accumulate into dummy to avoid optimizing out
    msh.met.setSpace(MetSpace::Exp);
    double dum2 = 0;
    double t02 = get_wall_time();
    for(int iloop = 0; iloop < nloop; iloop++){
      for(int ielem = 0; ielem < nentt; ielem++){
        if(idim == 2){
          dum2 += metqua2<2>(msh, ielem, 1);
        }else{
          dum2 += metqua2<3>(msh, ielem, 1);
        }
      }
    }// for iloop
    double t12 = get_wall_time();


    // Accumulate into dummy to avoid optimizing out
    msh.met.setSpace(MetSpace::Exp);
    double dum3 = 0;
    double t03 = get_wall_time();
    for(int iloop = 0; iloop < nloop; iloop++){
      for(int ielem = 0; ielem < nentt; ielem++){
        if(idim == 2){
          dum3 += metqua3<2>(msh, ielem, 1);
        }else{
          dum3 += metqua3<3>(msh, ielem, 1);
        }
      }
    }// for iloop
    double t13 = get_wall_time();

    int ps1 = (int) ((((double)nentt) * nloop) / (t11 - t01) / 1000);
    int ps2 = (int) ((((double)nentt) * nloop) / (t12 - t02) / 1000);
    int ps3 = (int) ((((double)nentt) * nloop) / (t13 - t03) / 1000);

    printf(" dum1 = %e dum2 = %e dum3 = %e\n",dum1,dum2,dum3);
    printf("-- Test end metqua1 time = %e metqua2 = %e metqua3 = %e\n"
           "ps 1: %dk/s %dk/s %dk/s\n"
           "ratio 2/1 = %e 3/1 = %e\n\n\n",
      t11 - t01, t12 - t02, t13 - t03, 
      ps1, ps2, ps3,
      (t12 - t02) / (t11 - t01), (t13 - t03) / (t11 - t01));

  }// for meshname : meshes

}// end boost test case






// Same as metqua for P1 only and uses tdim + 1 evals instead of avg metric
template <int idim>
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
template double metqua2<2>(Mesh<MFT> &msh, int ientt, double difto);
template double metqua2<3>(Mesh<MFT> &msh, int ientt, double difto);



template <int idim>
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
template double metqua1<2>(Mesh<MFT> &msh, int ientt, double difto);
template double metqua1<3>(Mesh<MFT> &msh, int ientt, double difto);



template <int idim>
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
template double metqua3<2>(Mesh<MFT> &msh, int ientt, double difto);
template double metqua3<3>(Mesh<MFT> &msh, int ientt, double difto);
