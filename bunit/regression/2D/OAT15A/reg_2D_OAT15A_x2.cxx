//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE reg_2D_OAT15A_x2

#include "../../reg_common.hxx"


using namespace Metris;

typedef MetricFieldAnalytical MFT;

BOOST_AUTO_TEST_CASE(reg_2D_OAT15A_x2) 
{

  RegressionTestManager manager("2D/OAT15A", "reg_x2");

  std::vector<int> l_adp_opt_niter  = { -1, 5, 1, 0 };

  MetrisParameters param;
  param.setMetricScale(1/sqrt(2));
  param.adp_niter = 30;
  param.opt_niter = 10;
  param.iverb = 2;
  param.ivdepth = 1;

  for(int tardeg : {1,2}){
    for(int adp_opt_niter : l_adp_opt_niter){

      std::string test_name = "Q1toQ" + std::to_string(tardeg) + "_cost" + std::to_string(adp_opt_niter);

      fmt::print("Call with name {}\n"  , test_name);
      param.adp_opt_niter = adp_opt_niter;
      param.usrTarDeg     = tardeg;

      try{
        manager.runTest(param, test_name, "OAT15A.meshb", "OAT15A.egads", "OAT15A.solb");
      }catch(const MetrisExcept& e){
        fmt::print(stderr,"## Test {} raised exception:\n{}\n", test_name, e.what());
      }

    }
  }

}