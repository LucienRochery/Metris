//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#define BOOST_TEST_MODULE My Test 

#include <boost/test/included/unit_test.hpp> 
#include <random>

#include "gen_bary.hxx"

#include "../src/utils/CT_loop.hxx"
#include "../src/utils/mprintf.hxx"
#include "../src/quality/low_metqua.hxx"
#include "../src/quality/low_metqua_d.hxx"
#include "../src/metris_options.hxx"
#include "../src/MetrisRunner/MetrisRunner.hxx"
#include "../src/Mesh/Mesh.hxx"

#include "../src/quality/low_metqua.hxx"
#include "../src/quality/quafun_tradet.hxx"
#include "../src/quality/quafun.hxx"

#include "SANS/Surreal/SurrealS.h"

namespace Metris{

typedef MetricFieldAnalytical MFT;
typedef double ftype;
typedef std::pair<AsDeg,AsDeg> AsDegPair;

BOOST_AUTO_TEST_CASE(test_metqua_d) 
{

  std::vector<std::string> meshes = 
  {METRIS_CASES_DIR "/2D/square.p1.10.meshb",
   METRIS_CASES_DIR "/1200_p1.meshb",
   METRIS_CASES_DIR "/2D/square.circmet.5k.curved.meshb",
  };


  for(auto testcase : meshes){
    std::string mesh_name = testcase;
    std::cout<<"-- Test case mesh = "<<mesh_name<<std::endl;
    cargHandler arg("-in "+mesh_name+" -verb 0 -anamet 1 -sclmet 1");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
    msh.met.setSpace(MetSpace::Log);

    dblAr2 bary[4];
    const int nsamp = 20;
    genBary(nsamp, 2, bary[2]);
    genBary(nsamp, 3, bary[3]);

    GETVDEPTH(msh.param);

    CT_FOR0_INC(2,3,gdim){if(msh.idim == gdim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
      constexpr int nhess = (gdim*(gdim+1))/2;
      // Use Boost PP to loop between options in the future
      constexpr auto iquaf = QuaFun::Distortion;


      CT_FOR0_INC(2,gdim,tdim){if(tdim <= msh.get_tdim()){
        INCVDEPTH(msh.param);
        MPRINTF("-- Tests tdim = %d gdim = %d\n",tdim,gdim);
        constexpr auto quafun_xi   = get_quafun_xi<MFT,gdim,tdim,iquaf,ftype>();
        constexpr auto d_quafun_xi = get_d_quafun_xi<MFT,gdim,tdim,iquaf,ftype>();
        const intAr2 &ent2poi = msh.ent2poi(tdim);
        const int nnode = getnnode(tdim, ideg);

        {//INCVDEPTH
        INCVDEPTH(msh.param);
        MPRINTF("-- Test 1: quafun_xi == d_quafun_xi\n");
        for(AsDegPair asdeg : {AsDegPair({AsDeg::P1, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::Pk})}){
          std::string asdmsh_str = asdeg.first  == AsDeg::P1 ? "P1" : "Pk";
          std::string asdmet_str = asdeg.second == AsDeg::P1 ? "P1" : "Pk";
          double err_min = 1.0e30;
          double err_max = -1.0;
          double err_avg = 0;
          double ntest = 0;
          for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
            for(int isamp = 0; isamp < nsamp; isamp++){
              double qua = quafun_xi(msh, asdeg.first, asdeg.second, ent2poi[ientt], 
                                     bary[tdim][isamp], NULL);
              double dquael[gdim], hquael[nhess];
              for(int ivar = 0; ivar < nnode; ivar++){
                double qua_d = d_quafun_xi(msh, asdeg.first, asdeg.second, ent2poi[ientt], 
                                           bary[tdim][isamp], NULL,
                                           ivar, msh.getBasis(), DifVar::None, 
                                           dquael, hquael);
                double err = abs(qua - qua_d);
                err_min = MIN(err_min, err);
                err_max = MAX(err_max, err);
                err_avg += err;
                ntest += 1;
              }
            }// for isamp
          }// for ientt
          if(ntest == 0) continue;
          err_avg /= ntest;
          MPRINTF(" - DONE using asdmsh = %s, asdmet = %s:" 
                  " got err min = %e, avg = %e, max = %e\n",
                  asdmsh_str.c_str(),asdmet_str.c_str(),
                  err_min,err_avg,err_max);
          BOOST_TEST(err_max < 1.0e-12);
        }// for AsDegPair



        MPRINTF("-- Test 2.1: d_quafun_tradet_SurrealS derivatives (identity)\n");
        {// using
        constexpr int nvar = gdim;
        auto d_quafun = d_quafun_tradet_SurrealS<MFT,gdim,tdim,gdim,ftype>;
        auto   quafun =   quafun_tradet         <MFT,gdim,tdim,     ftype>;
        SANS::DLA::MatrixS<gdim,nvar,double> dpoint = SANS::DLA::Identity();
        //for(int ii = 0; ii < gdim; ii++){
        //  for(int jj = 0; jj < gdim; jj++){
        //    dpoint(ii,jj) = ii == jj;
        //  }
        //}
        for(AsDegPair asdeg : {AsDegPair({AsDeg::P1, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::Pk})}){
          double ntest = 0;
          for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
            for(int isamp = 0; isamp < nsamp; isamp++){
              SANS::SurrealS<nvar, ftype> tra, det;
              ftype tra0, det0;

              quafun(msh, asdeg.first, asdeg.second,
                     ent2poi[ientt], bary[tdim][isamp],
                     NULL, 
                     &tra0, &det0);
              for(int ivar = 0; ivar < nnode; ivar++){
                d_quafun(msh, asdeg.first, asdeg.second,
                         ent2poi[ientt], bary[tdim][isamp],
                         ivar, msh.getBasis(), DifVar::None,
                         NULL, 
                         tra, det, &dpoint);
                ntest += 1;
              }
            }// for isamp
          }// for ientt
        }// for AsDegPair
        }// using

        }// INCVDEPTH

      }}CT_FOR1(tdim);
    }}CT_FOR1(ideg);
    }}CT_FOR1(gdim);


  }// for testcase


}// BOOST_AUTO_TEST_CASE

}// namespace