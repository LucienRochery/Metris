//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE reg_2D_OAT15A_x05

#include "../../reg_common.hxx"


using namespace Metris;


BOOST_AUTO_TEST_CASE(reg_2D_OAT15A_x05)
{

  RegressionTestManager manager("2D/OAT15A", "reg_x0.5");

  std::vector<int> l_adp_opt_niter  = { -1, 5, 1, 0 };

  MetrisParameters param;
  param.setMetricScale(sqrt(2));
  param.adp_niter = 30;
  param.opt_niter = 10;
  param.iverb = 2;
  param.ivdepth = 1;

  for(int tardeg : {1,2}){
    for(int adp_opt_niter : l_adp_opt_niter){
      for(bool adp_geo_lines : {true, false}){

        std::string test_name = "Q1toQ" + std::to_string(tardeg) 
                              + "_cost" + std::to_string(adp_opt_niter);
        if(adp_geo_lines) test_name += "_geo";

        param.adp_line_adapt = adp_geo_lines;
        param.adp_opt_niter  = adp_opt_niter;
        param.usrTarDeg      = tardeg;

        try{
          manager.runTest(param, test_name, "OAT15A.meshb", "OAT15A.egads", "OAT15A.solb");
        }catch(const MetrisExcept& e){
          fmt::print(stderr,"## Test {} raised exception:\n{}\n", test_name, e.what());
        }
      }

    }
  }

}