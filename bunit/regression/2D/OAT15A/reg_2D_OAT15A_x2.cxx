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

  std::vector<int> l_adp_opt_niter  = { -1, 5, 1, 0 };
  const double metScale = 1/sqrt(2);

  std::string time_stamp = time2str();

  std::string case_dir = METRIS_REGRESSION_DIR "/2D/OAT15A/";
  std::string out_dir  = case_dir + "/reg_x2/" + time_stamp + "/";
  std::string tmp_dir  = case_dir + "/reg_x2/tmp/";

  std::filesystem::create_directories(tmp_dir);
  std::filesystem::create_directories(out_dir);

  nlohmann::json json_all_tests;
  json_all_tests["metadata"]["time_stamp"] = time_stamp;
  json_all_tests["metadata"]["git_repo"] = METRIS_GIT_URL;
  json_all_tests["metadata"]["git_hash"] = METRIS_GIT_COMMIT_HASH;

  json_all_tests["runs"] = nlohmann::json::array();

  for(int adp_opt_niter : l_adp_opt_niter){

    for(int tardeg : {1,2}){

      std::string test_name = "Q1toQ" + std::to_string(tardeg) + "_cost" + std::to_string(adp_opt_niter);
      std::string outfile = out_dir + test_name;

      fmt::print("\n========================================\n");
      fmt::print("Running case {}\n",test_name);

      MetrisParameters param;
      param.setMeshIn(case_dir + "OAT15A.meshb");
      param.setMetricFile(case_dir + "OAT15A.solb");
      param.setCAD(case_dir + "OAT15A.egads");
      param.setMetricScale(metScale);

      param.adp_niter = 30;
      param.adp_opt_niter = adp_opt_niter;
      param.opt_niter = 10;
      param.iverb = 2;
      param.ivdepth = 1;
      param.outmPrefix = tmp_dir;
      param.usrTarDeg = tardeg;

      param.setLogFile(outfile + ".log");
      param.setMeshOut(outfile);


      //cargHandler arg("-in OAT15A -cad OAT15A -met OAT15A -verb 2 -vdepth 1 -prefix " + tmp_dir
      //  + " -adp-opt-niter " + std::to_string(adp_opt_niter) 
      //  + " -sclmet " + sclmet
      //  + " -adapt 30 -opt-niter 10 -out " + outfile + ".meshb"
      //  + " -log " + outfile + ".log"
      //  + tardeg
      //  );
      //MetrisRunner run(arg.c, arg.v);

      double t0 = get_wall_time();
      MetrisRunner run(NULL, param);

      run.runMetris();
      double t1 = get_wall_time();
      fmt::print("Total time: {:.2f} s\n", t1 - t0);
      
      MeshStat stat;
      run.statMesh(0, &stat);

      stat.print();

      nlohmann::json json_entry;
      json_entry["params"] = *run.param;
      json_entry["result"] = stat;
      json_entry["CPU_time"] = t1 - t0;
      json_entry["test_name"] = test_name;
      json_entry["logfile"] = outfile + ".log";

      json_all_tests["runs"].push_back(json_entry);
    }
    
  }

  std::ofstream json_file(out_dir + "runs.json");
  json_file << json_all_tests.dump(2);  // Pretty print with 2-space indent


}