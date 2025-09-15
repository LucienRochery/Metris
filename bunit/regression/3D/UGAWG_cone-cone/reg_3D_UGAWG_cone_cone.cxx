//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE reg_3D_UGAWG_cone_cone

#include "../../reg_common.hxx"


using namespace Metris;


BOOST_AUTO_TEST_CASE(reg_3D_UGAWG_cone_cone_uniform)
{

  RegressionTestManager manager("3D/UGAWG_cone-cone", "reg_uniform", "3D_cone_cone_uniform");

  int ianamet = 6;
  std::string metname = "uniform";

  MetrisParameters param;
  param.adp_niter = 20;
  param.opt_niter = 10;
  param.adp_opt_niter = 1;
  param.refineConventionsInp = true;

  for(bool adp_geo_lines : {true, false}){
    for(int adp_smoolen : {true, false}){

      std::string test_name = metname;
      if(adp_geo_lines) test_name += "_geo";
      if(adp_smoolen) test_name += "_smoolen";

      param.setAnalyticalMetric(ianamet);
      param.setMetricScale(1.0);
      param.adp_line_adapt = adp_geo_lines;
      param.adp_smoo_len = adp_smoolen;

      try{
        manager.runTest(param, test_name, "cone-cone.meshb", "cone-cone.egads", "");
      }catch(const MetrisExcept& e){
        fmt::print(stderr,"## Test {} raised exception:\n{}\n", test_name, e.what());
      }
    }
  }


}
