//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE reg_2D_OAT15A_x2

#include "../../reg_common.hxx"


using namespace Metris;


BOOST_AUTO_TEST_CASE(reg_2D_OAT15A_x2) 
{

  RegressionTestManager manager("3D/cube", "reg_anamet");

  std::vector<int> l_adp_opt_niter  = {5, 1, 0 };
  std::vector<std::pair<int,double>> l_ianamet = { {1,0.5}, {2,0.05} };

  MetrisParameters param;
  param.adp_niter = 20;
  param.opt_niter = 10;
  param.iverb = 2;
  param.ivdepth = 1;
  

  for(int tardeg : {1,2}){
    for(int adp_opt_niter : l_adp_opt_niter){
      for(bool adp_geo_lines : {true, false}){

        for(std::pair<int,double> ianamet : l_ianamet){
          std::string test_name = "Q1toQ" + std::to_string(tardeg) 
                                + "_cost" + std::to_string(adp_opt_niter)
                                + "_met" + std::to_string(ianamet.first);
          if(adp_geo_lines) test_name += "_geo";

          param.setAnalyticalMetric(ianamet.first);
          param.setMetricScale(ianamet.second);
          param.adp_line_adapt = adp_geo_lines;
          param.adp_opt_niter  = adp_opt_niter;
          param.usrTarDeg      = tardeg;

          try{
            manager.runTest(param, test_name, "cube.meshb", "cube.egads", "");
          }catch(const MetrisExcept& e){
            fmt::print(stderr,"## Test {} raised exception:\n{}\n", test_name, e.what());
          }
        }
      }

    }
  }

}