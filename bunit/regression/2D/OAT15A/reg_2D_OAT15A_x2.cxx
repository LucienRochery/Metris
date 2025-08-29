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

  nlohmann::json json_current;
  json_current["metadata"]["time_stamp"] = time_stamp;
  json_current["metadata"]["git_repo"] = METRIS_GIT_URL;
  json_current["metadata"]["git_hash"] = METRIS_GIT_COMMIT_HASH;

  //json_current["runs"] = nlohmann::json::array();
  json_current["runs"] = nlohmann::json::object();

  std::vector<std::string> test_names;

  for(int adp_opt_niter : l_adp_opt_niter){

    for(int tardeg : {1,2}){

      std::string test_name = "Q1toQ" + std::to_string(tardeg) + "_cost" + std::to_string(adp_opt_niter);
      test_names.push_back(test_name);

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

      bool iexcept = false;
      std::string except_message;
      try{
        run.runMetris();
        //bool dummy = false;
        //METRIS_ENFORCE_MSG(dummy == true, "Dummy assertion\n");
      }catch(const MetrisExcept& e){
        iexcept = true;
        except_message = e.what();
        fmt::print(stderr,"## Test {} raised exception:\n{}\n", test_name, except_message);
      }
      double t1 = get_wall_time();
      fmt::print("Total time: {:.2f} s\n", t1 - t0);
      
      MeshStat stat;
      if(!iexcept){
        run.statMesh(0, &stat);
        stat.print();
      }

      nlohmann::json json_entry;
      json_entry["params"] = *run.param;
      if(!iexcept) json_entry["result"] = stat;
      else         json_entry["except"] = except_message;
      json_entry["CPU_time"] = t1 - t0;
      json_entry["logfile"]  = outfile + ".log";

      json_current["runs"][test_name] = json_entry;
    }
    
  }

  std::ofstream json_file(out_dir + "runs.json");
  json_file << json_current.dump(2);  // Pretty print with 2-space indent

  // Now look for the baseline json file and compare.
  std::string baseline_json_file = METRIS_REGRESSION_DIR "/2D/OAT15A/reg_2D_OAT15A_x2.baseline.json";
  if(std::filesystem::exists(baseline_json_file)){
    std::ifstream jfile(baseline_json_file);
    nlohmann::json json_baseline;
    jfile >> json_baseline;
    for(const std::string& test_name : test_names){
      fmt::print("========================================\n");
      fmt::print("Comparing case {}\n",test_name);

      BOOST_REQUIRE(json_baseline["runs"].contains(test_name));
      BOOST_REQUIRE(json_current["runs"].contains(test_name));

      auto json_baseline_run = json_baseline["runs"][test_name];
      auto json_current_run  = json_current["runs"][test_name];

      BOOST_REQUIRE(!json_baseline_run.contains("except"));

      BOOST_TEST(!json_current_run.contains("except"));
      if(json_current_run.contains("except")){
        std::string except_message = json_current_run["except"];
        fmt::print("Test {} raised exception:\n{}\nSkipping comparison.\n", 
                   test_name, except_message);
        continue;
      }

      MeshStat baseline_stat = json_baseline_run["result"];
      MeshStat current_stat = json_current_run["result"];

      bool iclose_pctunit, iclose_pctunit_bdry, iclose_avgqua, iclose_avgqua_bdry;
      bool icloseStat = isCloseMeshStat(baseline_stat, current_stat,
          &iclose_pctunit, &iclose_pctunit_bdry,
          &iclose_avgqua, &iclose_avgqua_bdry);

      if(!icloseStat){
        fmt::print(stderr, "Mesh stats diverge for test {}:\n",test_name);
        printUnmatchedMeshStat(baseline_stat, current_stat, 
                               iclose_pctunit, iclose_pctunit_bdry,
                               iclose_avgqua, iclose_avgqua_bdry);
        baseline_stat.print("baseline",stderr);
        current_stat.print("current",stderr);
      }
      BOOST_CHECK_MESSAGE(icloseStat, "MeshStat for test " << test_name << " differ from baseline");
      if(icloseStat)
        fmt::print("Current test statistics better or equal as baseline\n");

      double baseline_time = json_baseline_run["CPU_time"];
      double current_time  = json_current_run["CPU_time"];
      bool iclose_time = current_time <= 1.2 * baseline_time;
      BOOST_CHECK_MESSAGE(iclose_time, 
                         "CPU time for test " << test_name << " increased > 1.2x, from " 
                          << baseline_time << " s to " << current_time << " s");
      if(iclose_time)
        fmt::print("Current test time {:.2f}s ~= {:.2f}s of baseline\n", current_time, baseline_time);
    }

  }else{
    //fmt::print("No baseline file found at {}, cannot compare results.\n"
    //           "If the results are expected, please copy {} to that location.\n",
    //           baseline_json_file, out_dir + "runs.json");   
  }
}