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
  #ifdef METRIS_GIT_DIRTY
  bool update_baseline = false;
  #else
  bool update_baseline = true;
  #endif

  if(update_baseline) fmt::print("-- Clean git working tree, will update baseline.\n");
  else                fmt::print("## Uncommitted changes: will not update baseline.\n");
  
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

  HardwareID hwid;
  fmt::print("-- Running on {} \n", hwid.to_string());

  for(int tardeg : {1,2}){
    for(int adp_opt_niter : l_adp_opt_niter){


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

      double t0 = get_cpu_time();
      double t0b = get_wall_time();
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
      double t1 = get_cpu_time();
      double t1b = get_wall_time();
      fmt::print("Total time: wall = {:.2f}s, user = {:.2f}s\n", t1b - t0b, t1 - t0);
      
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
  std::string baseline_json_fname = METRIS_REGRESSION_DIR "/2D/OAT15A/reg_2D_OAT15A_x2.baseline.json";
  if(std::filesystem::exists(baseline_json_fname)){
    std::ifstream jfile(baseline_json_fname);
    nlohmann::json json_baseline;
    jfile >> json_baseline;
    for(const std::string& test_name : test_names){
      fmt::print("========================================\n");
      fmt::print("Comparing case {}\n",test_name);

      BOOST_REQUIRE(json_baseline["runs"].contains(test_name));
      BOOST_REQUIRE(json_current["runs"].contains(test_name));

      bool test_passes = true;

      nlohmann::json& json_baseline_run = json_baseline["runs"][test_name];
      nlohmann::json& json_current_run  = json_current["runs"][test_name];

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
        test_passes = false;
      }
      BOOST_CHECK_MESSAGE(icloseStat, "MeshStat for test " << test_name << " differ from baseline");
      if(icloseStat)
        fmt::print("Current test statistics better or equal as baseline\n");

      double baseline_time = get_baseline_CPU(json_baseline_run);
      double current_time  = json_current_run["CPU_time"];
      bool iclose_time = current_time <= baseline_time;
      BOOST_CHECK_MESSAGE(iclose_time, 
                         "CPU time for test " << test_name << " increased > 1.2x, from " 
                          << baseline_time << " s to " << current_time << " s");
      if(iclose_time)
        fmt::print("Current test time {:.2f}s <= {:.2f}s baseline bound\n", current_time, baseline_time);
      else{
        fmt::print(stderr, "CPU time {:.2f}s > {:.2f}s baseline bound\n", baseline_time, current_time);
        test_passes = false;
      }

      if(!test_passes || !update_baseline) continue;

      // Update the baseline CPU and stats
      json_baseline_run["CPU_times"].push_back(current_time);
      // Update stat and metadata if this case improved.
      if(current_stat.pctunit > baseline_stat.pctunit){
        json_baseline_run["result"] = current_stat;
        // metadata is stored per run on the baseline, as each case 
        // might be updated at a different time.
        json_baseline_run["metadata"] = json_current["metadata"]; 
      }


    }

  }else if(update_baseline){
    fmt::print("## No baseline file found at {}, creating it.\n",
               baseline_json_fname);   
    nlohmann::json json_baseline = json_current;
    nlohmann::json metadata = json_current["metadata"];
    // Remove CPU_time and log into CPU_times.
    for(const std::string& test_name : test_names){
      auto& json_run = json_baseline["runs"][test_name];
      double cpu_time = json_run["CPU_time"];
      json_run.erase("CPU_time");
      json_run["CPU_times"] = nlohmann::json::array();
      json_run["CPU_times"].push_back(cpu_time);
      // metadata stored per run on the baseline, as each case 
      // might be updated at a different time.
      json_run["metadata"] = metadata;
    }
    
    std::ofstream baseline_json_file(baseline_json_fname);
    baseline_json_file << json_baseline.dump(2);  // Pretty print with 2-space indent
  }

  if (boost::unit_test::results_collector.results(boost::unit_test::framework::current_test_case().p_id).passed() == false) {
    fmt::print("## Some tests have failed, skipping baseline update.\n");
    return;
  }
}