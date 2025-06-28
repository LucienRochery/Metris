//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#define BOOST_TEST_MODULE MyTest 


#include <egads.h>

#include <random>
#include <string>
#include <regex>

#include "common_setup.hxx"
#include "../src/utils/aux_pp_inc.hxx"
#include "../src/utils/mprintf.hxx"
#include "../src/Localization/low_localization.hxx"
#include "../src/low_ccoef.hxx"
#include "../src/low_normal.hxx"
#include "../src/linalg/invmat.hxx"

#include <filesystem>
#include <boost/hana.hpp> 
#include <nlopt.hpp>

namespace hana = boost::hana;
using namespace hana::literals;

namespace utf = boost::unit_test;

namespace Metris {

typedef MetricFieldAnalytical MFT;



BOOST_AUTO_TEST_CASE(test_inveval) 
{
  
  std::vector<std::string> meshes = {
    //"../cases/2D/square.circmet.50.curved.meshb",
    //"../cases/2D/square.circmet.5k.curved.meshb",
    "../cases/1200_p1.meshb -back ../cases/2400_p1.meshb",
    //"../cases/invevalP2_2",
  };

  const int nsamp = 100;
  dblAr2 bar1(nsamp,2),bar2(nsamp,3),bar3(nsamp,4);
  dblAr2 bar1_out(nsamp,2),bar2_out(nsamp,3),bar3_out(nsamp,4);
  std::uniform_real_distribution<double> unif(0.0,1.0);
  std::default_random_engine rng(0);
  for(int ii = 0; ii < nsamp; ii++){
    double sum1 = 0, sum2 = 0, sum3 = 0;
    do{
      for(int jj = 0; jj < 4; jj++){
        if(jj < 2) bar1(ii,jj) = unif(rng);
        if(jj < 3) bar2(ii,jj) = unif(rng);
                   bar3(ii,jj) = unif(rng);
        if(jj < 2) sum1 += bar1(ii,jj);
        if(jj < 3) sum2 += bar2(ii,jj);
                   sum3 += bar3(ii,jj);
      }
    }while(abs(sum1) < 1.0e-16 || abs(sum2) < 1.0e-16 || abs(sum3) < 1.0e-16);

    for(int jj = 0; jj < 4; jj++){
      if(jj < 2) bar1(ii,jj) /= sum1;
      if(jj < 3) bar2(ii,jj) /= sum2;
                 bar3(ii,jj) /= sum3;
    }
  }
  std::uniform_real_distribution<double> unif_out(-100.0,100.0);
  std::default_random_engine rng_out(0);
  for(int ii = 0; ii < nsamp; ii++){
    double sum1 = 0, sum2 = 0, sum3 = 0;
    do{
      for(int jj = 0; jj < 4; jj++){
        if(jj < 2) bar1_out(ii,jj) = unif_out(rng_out);
        if(jj < 3) bar2_out(ii,jj) = unif_out(rng_out);
                   bar3_out(ii,jj) = unif_out(rng_out);
        if(jj < 2) sum1 += bar1_out(ii,jj);
        if(jj < 3) sum2 += bar2_out(ii,jj);
                   sum3 += bar3_out(ii,jj);
      }
    }while(abs(sum1) < 1.0e-16 || abs(sum2) < 1.0e-16 || abs(sum3) < 1.0e-16);

    for(int jj = 0; jj < 4; jj++){
      if(jj < 2) bar1_out(ii,jj) /= sum1;
      if(jj < 3) bar2_out(ii,jj) /= sum2;
                 bar3_out(ii,jj) /= sum3;
    }
  }


  int nsucc0 = 0;
  int nerro0 = 0;

  std::filesystem::create_directory("tmp");
  const int iverb = 0;

  for(auto s : meshes)
  {
    std::cout << "Mesh " << s << std::endl;
    // We want to alter mesh somewhat so back and front don't coincide
    //cargHandler arg("-in " + s + " -sclmet 2 -verb 2 -vdepth 5 -adapt 2 -prefix tmp/");
    cargHandler arg("-in " + s + " -sclmet 2 -adapt 2 -prefix tmp/ -verb 0");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MetricFieldFE> &msh = *((Mesh<MetricFieldFE>*) run.msh_g);

    if(msh.idim == 2){
      msh.param->iverb = 0;
      run.adaptMesh();
      msh.param->iverb = iverb;
    }else{
      msh.param->iverb = 2;
      msh.param->ivdepth= 5;
    }

    INCVDEPTH(msh.param);

    msh.cleanup();

    int nsucc1 = 0;
    int nerro1 = 0;

    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
      CT_FOR0_INC(2,3,gdim){if(gdim == msh.idim){
        CT_FOR0_INC(1,gdim,tdim){

          printf("\n\n   - start tdim %d/%d\n",tdim,gdim);
          int nsucc2 = 0;
          int nerro2 = 0;

          const intAr2& ent2poi = msh.ent2poi(tdim);
          const intAr1& ent2ref = msh.ent2ref(tdim);
          int nentt = msh.nentt(tdim);
          auto evalf = tdim == 1 ? eval1<gdim,ideg> : 
                       tdim == 2 ? eval2<gdim,ideg> : 
                                   eval3<gdim,ideg>;

          dblAr2& bary = tdim == 1 ? bar1 :
                         tdim == 2 ? bar2 : bar3;
          int nlast = 0;
          for(int ientt = 0; ientt < nentt; ientt++){
            for(int isamp = 0; isamp < nsamp; isamp++){
              int ipnew = msh.newpoitopo(tdim, ientt);
              evalf(msh.coord,ent2poi[ientt],msh.getBasis(),DifVar::None,
                    DifVar::None,bary[isamp],msh.coord[ipnew],NULL,NULL);
              double algnd[gdim];
              if(tdim < gdim){
                msh.newbpotopo(ipnew, tdim, ientt);
                if(tdim == 1){
                  for(int ii = 0; ii < gdim; ii++){
                    algnd[ii] = msh.coord(msh.edg2poi(ientt,0),ii)
                              - msh.coord(msh.edg2poi(ientt,1),ii);
                  }
                }else{
                  getnorfacP1(ent2poi[ientt], msh.coord, algnd);
                }
              }

              int ierro = msh.interpMetBack(ipnew, tdim, ientt, ent2ref[ientt], algnd);

              if(ierro == 0) nsucc2++;
              else           nerro2++;


              msh.killpoint(ipnew);
              msh.set_npoin(msh.npoin-1);
            }// for isamp
            if(nlast++ == 100){
              msh.cleanup();
              nlast = 0;
            }

          }// for ientt

          nsucc1 += nsucc2;
          nerro1 += nerro2;
          printf("   - tdim %d/%d nsucc = %d nerro = %d = %f%%\n",
            tdim,gdim,nsucc2,nerro2,(100*nerro2/(double)(nerro2 + nsucc2)));
        }CT_FOR1(tdim);

      }}CT_FOR1(gdim);
    }}CT_FOR1(ideg);

    nsucc0 += nsucc1;
    nerro0 += nerro1;
    printf("-- Mesh %s nsucc = %d nerro = %d \n",s.c_str(),nsucc1, nerro1);

  }// for auto s : meshes
  printf("-- End all tests nsucc = %d nerro = %d \n",nsucc0, nerro0);
}

} //namespace Metris 









