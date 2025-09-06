//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE reg_2D_OAT15A_x05

#include "../../reg_common.hxx"


using namespace Metris;


BOOST_AUTO_TEST_CASE(reg_2D_OAT15A_x05_Q2toQ2) 
{

  RegressionTestManager manager("2D/OAT15A", "reg_x0.5_Q2toQ2", "2D_OAT15A_x0.5_Q2toQ2");

  std::vector<int> l_adp_opt_niter  = {0, 1, 5, -1};

  MetrisParameters param;
  param.setMetricScale(sqrt(2));
  param.adp_niter = 30;
  param.opt_niter = 10;
  //fmt::print("## DEBUG SET ITERATIONS FROM 30, 10 to 0 0\n");
  param.iverb = 2;
  param.ivdepth = 1;
  param.usrTarDeg = 2;
  
  MetrisParameters param_degelev;
  param_degelev.adp_niter = 0;
  param_degelev.opt_niter = 10;
  param_degelev.iverb = 2;
  param_degelev.ivdepth = 1;
  param_degelev.usrTarDeg = 2;
  param_degelev.setMeshIn(manager.getMeshIn("OAT15A.meshb"));
  param_degelev.setMetricFile(manager.getMetricIn("OAT15A.solb"));
  param_degelev.setCAD(manager.getCADIn("OAT15A.egads"));
  param_degelev.setLogFile(manager.getOutDir() + "degelev.log");

  MetrisRunner run_degelev(NULL, NULL, param_degelev);
  run_degelev.degElevate();
  MetrisAPI HOdata(run_degelev);


  for(int adp_opt_niter : l_adp_opt_niter){
    for(bool adp_geo_lines : {true, false}){

      std::string test_name = "Q2toQ2_cost" + std::to_string(adp_opt_niter);
      if(adp_geo_lines) test_name += "_geo";

      param.adp_line_adapt = adp_geo_lines;
      param.adp_opt_niter  = adp_opt_niter;
      param.usrTarDeg      = 2;

      try{
        // Hard copy because the runner moves the data.
        MetrisAPI HOdata_current(HOdata);
        manager.runTest(param, test_name, &HOdata_current, NULL);
      }catch(const MetrisExcept& e){
        fmt::print(stderr,"## Test {} raised exception:\n{}\n", test_name, e.what());
      }
    }

  }

}


BOOST_AUTO_TEST_CASE(reg_2D_OAT15A_x05_Q1toQ1) 
{

  RegressionTestManager manager("2D/OAT15A", "reg_x0.5_Q1toQ1", "2D_OAT15A_x0.5_Q1toQ1");

  std::vector<int> l_adp_opt_niter  = {0, 1, 5, -1};

  MetrisParameters param;
  param.setMetricScale(sqrt(2));
  param.adp_niter = 30;
  param.opt_niter = 10;
  param.iverb = 2;
  param.ivdepth = 1;

  for(int adp_opt_niter : l_adp_opt_niter){
    for(bool adp_geo_lines : {true, false}){

      std::string test_name = "Q1toQ1_cost" + std::to_string(adp_opt_niter);
      if(adp_geo_lines) test_name += "_geo";

      param.adp_line_adapt = adp_geo_lines;
      param.adp_opt_niter  = adp_opt_niter;

      try{
        manager.runTest(param, test_name, "OAT15A.meshb", "OAT15A.egads", "OAT15A.solb");
      }catch(const MetrisExcept& e){
        fmt::print(stderr,"## Test {} raised exception:\n{}\n", test_name, e.what());
      }
    }

  }

}


BOOST_AUTO_TEST_CASE(reg_2D_OAT15A_x05_Q1toQ2) 
{

  RegressionTestManager manager("2D/OAT15A", "reg_x0.5_Q1toQ2", "2D_OAT15A_x0.5_Q1toQ2");

  std::vector<int> l_adp_opt_niter  = {0, 1, 5, -1};

  MetrisParameters param;
  param.setMetricScale(sqrt(2));
  param.adp_niter = 30;
  param.opt_niter = 10;
  param.iverb = 2;
  param.ivdepth = 1;
  param.usrTarDeg = 2;

  for(int adp_opt_niter : l_adp_opt_niter){
    for(bool adp_geo_lines : {true, false}){

      std::string test_name = "Q1toQ2_cost" + std::to_string(adp_opt_niter);
      if(adp_geo_lines) test_name += "_geo";

      param.adp_line_adapt = adp_geo_lines;
      param.adp_opt_niter  = adp_opt_niter;

      try{
        manager.runTest(param, test_name, "OAT15A.meshb", "OAT15A.egads", "OAT15A.solb");
      }catch(const MetrisExcept& e){
        fmt::print(stderr,"## Test {} raised exception:\n{}\n", test_name, e.what());
      }
    }

  }

}

