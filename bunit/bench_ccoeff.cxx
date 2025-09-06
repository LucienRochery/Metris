//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE bench_ccoef

#include "common_setup.hxx"


using namespace Metris;

typedef MetricFieldAnalytical MFT;

namespace utf = boost::unit_test;
// 1.0e-8 relative error -> 1.0e-6% utf::tolerance()
BOOST_AUTO_TEST_CASE(bench_ccoef, * utf::tolerance(double(1.0e-6)) )  
{ 
  #ifdef NDEBUG

  // bool is whether straight
  std::vector<std::pair<std::string,bool>> meshes = {
   {METRIS_CASES_DIR "/unit/2D/square/iso.p1.100k",true}
  ,{METRIS_CASES_DIR "/unit/2D/square/circmet.p2.500",false}
  ,{METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k",true}
  ,{METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k -t 2",true}
  ,{METRIS_CASES_DIR "/unit/3D/cube/curved.p2.2k",false}
  #if METRIS_MAX_DEG >= 3
  ,{METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k -t 3",true}
  ,{METRIS_CASES_DIR "/unit/3D/cube/curved.p2.2k -t 3",false}
  #if METRIS_MAX_DEG >= 4
  ,{METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k -t 4",true}
  ,{METRIS_CASES_DIR "/unit/3D/cube/curved.p2.2k -t 4",false}
  #endif
  #endif
  };


  double tol = 1.0e-12;
  for(auto testcase : meshes)
  {
    std::string s = testcase.first;
    bool istr8    = testcase.second;

    cargHandler arg("-in " + s + "  -anamet 1 -verb 0");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
    // For those meshes that have a higher target degree than 1 (does nothing to the others)
    run.degElevate();
    msh.cleanup();

    std::cout<<"\n\n-- Mesh "<<s<<" dim "<<msh.idim<<" deg "<<msh.curdeg<<"\n";

    int ps;

    CT_FOR0_INC(2,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
    CT_FOR0_INC(2,3,idim){if(idim == msh.idim){

      msh.setBasis(FEBasis::Bezier);

      constexpr auto ordent = ORDELT(idim);
      constexpr int jdeg = idim*(ideg-1);
      constexpr int nnodj = idim == 2 ? getnnod2(jdeg) : getnnod3(jdeg);
      const int vol0 = ifact<idim>();
      const intAr2 &ent2poi = msh.ent2poi(idim);
      const int nentt = msh.nentt(idim);
      double ccoef[nnodj], ccoef2[nnodj];

      constexpr auto ccoef_genbez = idim == 2 ? ccoef_genbez2<ideg> : ccoef_genbez3<ideg>;
      constexpr auto evalf =  idim == 2 ? eval2<idim,ideg> : eval3<idim,ideg>;
      constexpr auto evalj =  idim == 2 ? eval2<1   ,jdeg> : eval3<1   ,jdeg>;

      printf("-- Release target: benchmarks\n");
      std::cout<<"Mesh "<<s<<"\n";
      double t0,t1,dum[3]={0.0};
      double jmin,jmax;

      printf("--- Codegen \n");
      jmin = 1.0e30;
      jmax =-1.0e30;
      t0 = get_cpu_time();
      for(int ielem = 0; ielem < nentt; ielem++){
        ccoef_genbez(ent2poi,msh.coord,ielem,ccoef);
        dum[0] += ccoef[0];
        double vol = getmeasentP1<idim>(ent2poi[ielem],msh.coord)*vol0;
        for(int i = 0; i < nnodj; i++){
          jmin = MIN(ccoef[i] / vol, jmin);
          jmax = MAX(ccoef[i] / vol, jmax);
        }
      }
      t1 = get_cpu_time();
      ps = (int)(nentt/(t1-t0));
      ps /= 1000;
      printf(" %2.0f P%d Full elt coefs %dk/s \n",dum[0],ideg,ps);
      printf("   %15.8e < J_K < %15.8e \n",jmin,jmax);
      if(istr8){
        BOOST_TEST(jmin == 1.0);
        BOOST_TEST(jmax == 1.0);
      }
    
      printf("--- Manual compute at nodes \n");
      jmin = 1.0e30;
      jmax =-1.0e30;
      t0 = get_cpu_time();
      for(int ielem = 0; ielem < nentt; ielem++){
        ccoef_eval<idim,idim,ideg>(msh.getBasis(),ent2poi,msh.coord,ielem,NULL,ccoef);
        dum[0] += ccoef[0];
        double vol = getmeasentP1<idim>(ent2poi[ielem],msh.coord)*vol0;
        for(int i = 0; i < nnodj; i++){
          jmin = MIN(ccoef[i] / vol, jmin);
          jmax = MAX(ccoef[i] / vol, jmax);
        }
      }
      t1 = get_cpu_time();
      ps = (int)(nentt/(t1-t0));
      ps /= 1000;
      printf(" %2.0f P%d Full elt coefs %dk/s \n",dum[0],ideg,ps);
      printf("   %15.8e < J_K < %15.8e \n",jmin,jmax);
      

    }}CT_FOR1(idim);
    }}CT_FOR1(ideg);
  }
  

  #else

    std::cout<<"Skipping test_ccoeff in debug mode\n";
    BOOST_TEST(true);

  #endif
}