//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE reg_3D_cube_anamet

#include "../../reg_common.hxx"


using namespace Metris;



BOOST_AUTO_TEST_CASE(reg_3D_cube_anamet_Q1toQ2) 
{

  RegressionTestManager manager("3D/cube", "reg_anamet_Q1toQ2", "3D_cube_anamet_Q1toQ2");

  std::vector<int> l_adp_opt_niter  = {5, 1, 0 };
  std::vector<std::pair<int,double>> l_ianamet = { {1,0.5}, {2,0.05} };

  MetrisParameters param;
  param.adp_niter = 20;
  param.opt_niter = 10;
  param.iverb = 2;
  param.ivdepth = 1;
  param.usrTarDeg = 2;

  for(int adp_opt_niter : l_adp_opt_niter){
    for(bool adp_geo_lines : {true, false}){

      for(std::pair<int,double> ianamet : l_ianamet){
        std::string test_name = "Q1toQ2_cost" + std::to_string(adp_opt_niter)
                              + "_met" + std::to_string(ianamet.first);
        if(adp_geo_lines) test_name += "_geo";

        param.setAnalyticalMetric(ianamet.first);
        param.setMetricScale(ianamet.second);
        param.adp_line_adapt = adp_geo_lines;
        param.adp_opt_niter  = adp_opt_niter;

        try{
          manager.runTest(param, test_name, "cube.meshb", "cube.egads", "");
        }catch(const MetrisExcept& e){
          fmt::print(stderr,"## Test {} raised exception:\n{}\n", test_name, e.what());
        }
      }
    }

  }

}


BOOST_AUTO_TEST_CASE(reg_3D_cube_anamet_Q1toQ1) 
{

  RegressionTestManager manager("3D/cube", "reg_anamet_Q1toQ1", "3D_cube_anamet_Q1toQ1");

  std::vector<int> l_adp_opt_niter  = {5, 1, 0 };
  std::vector<std::pair<int,double>> l_ianamet = { {1,0.5}, {2,0.05} };

  MetrisParameters param;
  param.adp_niter = 20;
  param.opt_niter = 10;
  param.iverb = 2;
  param.ivdepth = 1;
  

  for(int adp_opt_niter : l_adp_opt_niter){
    for(bool adp_geo_lines : {true, false}){

      for(std::pair<int,double> ianamet : l_ianamet){
        std::string test_name = "Q1toQ1_cost" + std::to_string(adp_opt_niter)
                              + "_met" + std::to_string(ianamet.first);
        if(adp_geo_lines) test_name += "_geo";

        param.setAnalyticalMetric(ianamet.first);
        param.setMetricScale(ianamet.second);
        param.adp_line_adapt = adp_geo_lines;
        param.adp_opt_niter  = adp_opt_niter;

        try{
          manager.runTest(param, test_name, "cube.meshb", "cube.egads", "");
        }catch(const MetrisExcept& e){
          fmt::print(stderr,"## Test {} raised exception:\n{}\n", test_name, e.what());
        }
      }
    }

  }

}



BOOST_AUTO_TEST_CASE(reg_3D_cube_anamet_Q2toQ2) 
{

  RegressionTestManager manager("3D/cube", "reg_anamet_Q2toQ2", "3D_cube_anamet_Q2toQ2");

  std::vector<int> l_adp_opt_niter  = {5, 1, 0 };
  std::vector<std::pair<int,double>> l_ianamet = { {1,0.5}, {2,0.05} };

  MetrisParameters param;
  param.adp_niter = 20;
  param.opt_niter = 10;
  param.iverb = 2;
  param.ivdepth = 1;

  MetrisParameters param_degelev;
  param_degelev.adp_niter = 0;
  param_degelev.opt_niter = 10;
  param_degelev.iverb = 2;
  param_degelev.ivdepth = 1;
  param_degelev.usrTarDeg = 2;
  param_degelev.setMeshIn(manager.getMeshIn("cube.meshb"));
  param_degelev.setCAD(manager.getCADIn("cube.egads"));
  param_degelev.setLogFile(manager.getOutDir() + "degelev.log");

  MetrisRunner run_degelev(NULL, NULL, param_degelev);
  run_degelev.degElevate();
  MetrisAPI HOdata(run_degelev);

  for(int adp_opt_niter : l_adp_opt_niter){
    for(bool adp_geo_lines : {true, false}){

      for(std::pair<int,double> ianamet : l_ianamet){
        std::string test_name = "Q1toQ2_cost" + std::to_string(adp_opt_niter)
                              + "_met" + std::to_string(ianamet.first);
        if(adp_geo_lines) test_name += "_geo";

        param.setAnalyticalMetric(ianamet.first);
        param.setMetricScale(ianamet.second);
        param.adp_line_adapt = adp_geo_lines;
        param.adp_opt_niter  = adp_opt_niter;

        try{
          MetrisAPI HOdata_current(HOdata);
          manager.runTest(param, test_name, &HOdata_current, NULL);
        }catch(const MetrisExcept& e){
          fmt::print(stderr,"## Test {} raised exception:\n{}\n", test_name, e.what());
        }
      }
    }

  }

}