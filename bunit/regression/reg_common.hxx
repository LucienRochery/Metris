//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_REG_COMMON__
#define __METRIS_REG_COMMON__

#include "../common_setup.hxx"
#include "reg_hwid.hxx"

#include "nlohmann/json.hpp"

#include <filesystem>
#include <chrono>
#include <string>
#include <vector>

namespace Metris{



// Initialize a field json_baseline_allhwid[hwid] passed in as json_baseline
void initialize_baseline_hwid(const std::vector<std::string>& test_names,
                              std::string str_hwid,
                                    nlohmann::json& json_baseline_allhwid, 
                              const nlohmann::json& json_current){

  json_baseline_allhwid[str_hwid] = json_current;
  nlohmann::json metadata = json_current["metadata"];

  json_baseline_allhwid[str_hwid].erase("metadata");

  // Remove CPU_time and log into CPU_times.
  for(const std::string& test_name : test_names){
    nlohmann::json& json_run = json_baseline_allhwid[str_hwid]["runs"][test_name];
    double cpu_time = json_run["CPU_time"];
    json_run.erase("CPU_time");
    json_run["CPU_times"] = nlohmann::json::array();
    json_run["CPU_times"].push_back(cpu_time);
    // metadata stored per run on the baseline, as each case 
    // might be updated at a different time.
    json_run["metadata"] = metadata;
  }
  
}


// Returns an upper bound on CPU time
double get_baseline_CPU(nlohmann::json baseline_json){
  METRIS_ENFORCE(baseline_json.contains("CPU_times"));

  int ncpu = baseline_json["CPU_times"].size();
  METRIS_ENFORCE(ncpu > 0);

  dblAr1 cpu_times(ncpu);

  double cpu_avg = 0.0;
  int ii = 0;
  for(double cpu_time : baseline_json["CPU_times"]){
    cpu_times[ii++] = cpu_time;
    cpu_avg += cpu_time;
  }
  cpu_avg /= ncpu;

  // Start at 1.5x and gradually go down to 1x
  if(ncpu < 5) return cpu_avg * 1.5 - 0.1*cpu_avg;

  // Otherwise, use 2.5 sigma, which is a little under 99% confidence
  // interval (assuming normal distribution)
  double cpu_stdev = 0.0;
  for(double cpu_time : cpu_times){
    cpu_stdev += (cpu_time - cpu_avg) * (cpu_time - cpu_avg);
  }
  cpu_stdev = sqrt(cpu_stdev / ncpu);
  return cpu_avg + 2.5*cpu_stdev;
}


// Return if untracked changes excluding untracked files
// and the whole bunit/ directory. 
bool isGitDirty() {
  FILE* pipe = popen("git status --porcelain --untracked-files=no 2>/dev/null | grep -v 'bunit/'", "r");
  if (!pipe) return true;
  
  int c = fgetc(pipe);
  pclose(pipe);
  
  return (c != EOF);
}



bool isCloseMeshStat(const MeshStat& baseline, const MeshStat& current,
                     bolAr1 &iclose_pctunit, 
                     bolAr1 &iclose_avgquaP1,
                     bolAr1 &iclose_avgquaPk){
  int tdim = baseline.pctunit.get_n();
  // We can't throw here as this'll be called in the destructor 
  // of RegressionTestManager
  if( tdim != baseline.quaP1.get_n()
   || tdim != current.len.get_n()
   || tdim != current.quaP1.get_n()){
    fmt::print(stderr,"## ERROR in isCloseMeshStat, inconsistent tdim\n");
    baseline.print("baseline",stderr);
    current.print("current",stderr);
    return false;
   }

  bool iHO = baseline.quaPk.get_n() > 0 || baseline.quaPk.get_n() > 0;
  if(iHO && (baseline.quaPk.get_n() != tdim 
          || current.quaPk.get_n() != tdim)){
    fmt::print(stderr,"## ERROR in isCloseMeshStat, inconsistent tdim (HO)\n");
    baseline.print("baseline",stderr);
    current.print("current",stderr);
    return false;
  }

  iclose_pctunit.set_n(tdim);
  iclose_avgquaP1.set_n(tdim);
  if(iHO) iclose_avgquaPk.set_n(tdim);
  else    iclose_avgquaPk.set_n(0);

  bool iret = true;
  for(int ii = 0; ii < tdim; ii++){
    // Set pct unit close to 1%
    iret = iret && 
      (iclose_pctunit[ii]  = abs(baseline.pctunit[ii] - current.pctunit[ii]) < 1.0);
    // Quality (0,1) differs by 0.01
    iret = iret && 
      (iclose_avgquaP1[ii] = abs(baseline.quaP1[ii].avg() - current.quaP1[ii].avg()) < 1.0e-2);

    if(!iHO) continue;
    // Quality (0,1) differs by 0.01
    iret = iret && 
      (iclose_avgquaPk[ii] = abs(baseline.quaPk[ii].avg() - current.quaPk[ii].avg()) < 1.0e-2);
  }

  return iret;
}

void printUnmatchedMeshStat(const MeshStat& baseline, const MeshStat& current,
                            bolAr1 &iclose_pctunit, 
                            bolAr1 &iclose_avgquaP1,
                            bolAr1 &iclose_avgquaPk){
  int tdim = baseline.pctunit.get_n();
  // We can't throw here as this'll be called in the destructor 
  // of RegressionTestManager
  if( tdim != baseline.quaP1.get_n()
   || tdim != current.len.get_n()
   || tdim != current.quaP1.get_n()){
    fmt::print(stderr,"## ERROR in isCloseMeshStat, inconsistent tdim\n");
    baseline.print("baseline",stderr);
    current.print("current",stderr);
    return;
   }

  bool iHO = baseline.quaPk.get_n() > 0 || baseline.quaPk.get_n() > 0;
  if(iHO && (baseline.quaPk.get_n() != tdim 
          || current.quaPk.get_n() != tdim)){
    fmt::print(stderr,"## ERROR in isCloseMeshStat, inconsistent tdim (HO)\n");
    baseline.print("baseline",stderr);
    current.print("current",stderr);
    return;
  }

  for(int ii = 0; ii < tdim; ii++){
    if(!iclose_pctunit[ii]) fmt::print(stderr," - {}D % unit edges differ: baseline {} vs current {}\n",
        ii+1,baseline.pctunit[ii], current.pctunit[ii]);
    if(!iclose_avgquaP1[ii]) fmt::print(stderr," - {}D avg quality P1 differ: baseline {} vs current {}\n",
        ii+1,baseline.quaP1[ii].avg(), current.quaP1[ii].avg());
    if(!iclose_avgquaPk[ii]) fmt::print(stderr," - {}D avg quality Pk differ: baseline {} vs current {}\n",
        ii+1,baseline.quaPk[ii].avg(), current.quaPk[ii].avg());
  }
}






/* --------------------------------------------------------------- */
class RegressionTestManager{
public:
  RegressionTestManager() = delete;

  // case_subdir example "2D/OAT15A"
  // out_subdir example "reg_x2" (subdir of 2D/OAT15A)
  RegressionTestManager(std::string case_subdir, std::string out_subdir, std::string baseline_name_base){
  
    #ifndef METRIS_ROOT_DIR
      #error "METRIS_ROOT_DIR not defined, should be in CMakeLists.txt"
    #endif

    if(std::string(METRIS_ROOT_DIR).empty()){
      METRIS_THROW_MSG("METRIS_ROOT_DIR not defined, should be in CMakeLists.txt");
    }

    std::string metris_regression_dir = METRIS_ROOT_DIR "/bunit/regression/";

    t0_all = get_cpu_time();
    t0w_all = get_wall_time();

    update_baseline = !isGitDirty();
    update_baseline = true;
    fmt::print(stderr,"## WARNING FORCED update_baseline = true\n");

    if(update_baseline) fmt::print("-- Clean git working tree, will update baseline.\n");
    else                fmt::print("## Uncommitted changes: will not update baseline.\n");
    
    time_stamp = time2str();
    
    case_dir = METRIS_ROOT_DIR "/bunit/regression/" + case_subdir + "/";
    out_dir  = case_dir + out_subdir + "/" + time_stamp + "/";
    tmp_dir  = case_dir + out_subdir + "/tmp/";

    std::filesystem::create_directories(tmp_dir);
    std::filesystem::create_directories(out_dir);
    
    fmt::print("-- Writing results to {}\n", out_dir);

    json_current["metadata"]["time_stamp"] = time_stamp;
    json_current["metadata"]["git_repo"] = METRIS_GIT_URL;
    json_current["metadata"]["git_hash"] = METRIS_GIT_COMMIT_HASH;
    json_current["runs"] = nlohmann::json::object();

    fmt::print("-- Running on {} \n", hwid.to_string());


    baseline_json_fname = metris_regression_dir + case_subdir + "/"
                        + baseline_name_base + ".baseline.json";

    fmt::print("-- Baseline file: {}\n", baseline_json_fname);
  }

  std::string getMeshIn(const std::string &msh_base_name) const {
    return case_dir + msh_base_name;
  }
  std::string getMetricIn(const std::string &met_base_name) const {
    if(met_base_name.empty()) return "";
    return case_dir + met_base_name;
  }
  std::string getCADIn(const std::string &CAD_base_name) const {
    if(CAD_base_name.empty()) return "";
    return case_dir + CAD_base_name;
  }
  std::string getOutDir() const {
    return out_dir;
  }

  // Call this to run from files
  // CAD_base_name and met_base_name can be ""
  void runTest(MetrisParameters& param, 
               std::string test_name,
               std::string mesh_base_name,
               std::string CAD_base_name,
               std::string met_base_name){
    
    std::string outfile = out_dir + test_name;
    
    param.iverb = 2;
    param.ivdepth = 1;
    param.setMeshIn(getMeshIn(mesh_base_name));
    if(!CAD_base_name.empty()) param.setCAD(getCADIn(CAD_base_name));
    if(!met_base_name.empty()) param.setMetricFile(getMetricIn(met_base_name));
    param.outmPrefix = tmp_dir;
    param.setLogFile(outfile + ".log");
    param.setMeshOut(outfile);

    MetrisRunner run(NULL, param);

    runTest(test_name, run);
  }

  // Call this to run from existing MetrisAPI data
  // (e.g. if some preprocessing was done)
  // Example reg_2D_OAT15A_x2.cxx test reg_2D_OAT15A_x2_Q2toQ2
  // Note we could also have passed a preprocessing function handle,
  // but its type is likely to change with time, requiring more maintenance.
  // This approach is more flexible, if a little more cumbersome. 
  void runTest(MetrisParameters& param, 
               std::string test_name,
               MetrisAPI *data_front,
               MetrisAPI *data_back,
               std::string CAD_base_name,
               std::string met_base_name){

    param.iverb = 2;
    param.ivdepth = 1;
    std::string outfile = out_dir + test_name;
    param.outmPrefix = tmp_dir;
    param.setLogFile(outfile + ".log");
    param.setMeshOut(outfile);

    if(!CAD_base_name.empty()) param.setCAD(getCADIn(CAD_base_name));
    if(!met_base_name.empty()) param.setMetricFile(getMetricIn(met_base_name));

    MetrisRunner run(data_front, data_back, param);
    runTest(test_name, run);
  }

private:
  void runTest(std::string test_name,
               MetrisRunner& run){

    
    fmt::print("\n========================================\n");
    fmt::print("Running case {}\n",test_name);
    
    test_names.push_back(test_name);


    double t1_case, t0_case = get_cpu_time();
    double t1w_case,t0w_case = get_wall_time();
  
    std::string except_message;
    nlohmann::json json_entry;
    json_entry["params"] = *run.param;
    try{
      run.runMetris();
      t1_case = get_cpu_time();
      t1w_case = get_wall_time();
      fmt::print("Total time: wall = {:.2f}s, user = {:.2f}s\n", t1w_case - t0w_case, t1_case - t0_case);
      
      MeshStat stat;
      run.statMesh(0, &stat);
      stat.print();
      json_entry["result"] = stat;
      json_entry["CPU_time"] = t1_case - t0_case;
    }catch(const MetrisExcept& e){
      except_message = e.what();
      fmt::print(stderr,"## Test {} raised exception:\n{}\n", test_name, except_message);
      json_entry["except"] = except_message;
    }
 
    json_entry["logfile"]  = run.param->logFileName;
    json_current["runs"][test_name] = json_entry;

  }
public:

  ~RegressionTestManager(){

    t1_all = get_cpu_time();
    t1w_all = get_wall_time();

    fmt::print("\n========================================\n");
    fmt::print("-- All tests done.\n");
    fmt::print("-- Total time: wall = {:.2f}s, user = {:.2f}s\n", t1w_all - t0w_all, t1_all - t0_all);
    fmt::print("-- Results written to {}\n", out_dir);
    fmt::print("\n========================================\n");

    std::ofstream json_file(out_dir + "runs.json");
    json_file << json_current.dump(2);  // Pretty print with 2-space indent

    // Now look for the baseline json file and compare.

    nlohmann::json json_baseline_allhwid;
    if(std::filesystem::exists(baseline_json_fname)){
      fmt::print("-- Baseline file {} found.\n",baseline_json_fname);
      std::ifstream jfile(baseline_json_fname);
      jfile >> json_baseline_allhwid;
    }else if(update_baseline){
      fmt::print("## No baseline file {} found, creating.\n",baseline_json_fname);   
    }

    // True both if file didn't exist before (empty object)
    // or if it did but HWID didn't. 
    if(!json_baseline_allhwid.contains(hwid.to_string())){
      if(update_baseline){
        fmt::print("-- Baseline doesn't contain data for HWID {}, initializing.\n", hwid.to_string());

        json_baseline_allhwid[hwid.to_string()] = nlohmann::json::object();

        initialize_baseline_hwid(test_names, hwid.to_string(), json_baseline_allhwid, json_current);

        std::ofstream baseline_json_file(baseline_json_fname);
        baseline_json_file << json_baseline_allhwid.dump(2);  // Pretty print with 2-space indent
      }else{
        fmt::print("## Baseline doesn't exist or doesn't contain data for HWID, and should not be updated (uncommitted changes)\n");
        BOOST_FAIL("No baseline data to compare against and dirty Git tree.");
      }
      return;
    }
    

    bool updates_done = false;
    nlohmann::json& json_baseline = json_baseline_allhwid[hwid.to_string()];
    for(const std::string& test_name : test_names){
      fmt::print("========================================\n");
      fmt::print("Comparing case {}\n",test_name);

      if(!json_baseline["runs"].contains(test_name)){
        fmt::print(stderr, "## ERROR: Baseline missing test {}\n", test_name);
      }
      BOOST_REQUIRE(json_baseline["runs"].contains(test_name));
      if(!json_current["runs"].contains(test_name)){
        fmt::print(stderr, "## ERROR: Current missing test {}\n", test_name);
      }
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

      // An error that can innocently occur if the baseline name 
      // matches another test suite's, e.g. OAT15A_x2 vs OAT15A_x0.5
      // if one is built by copy-pasting another and modifying...
      MetrisParameters baseline_param = json_baseline_run["params"];
      MetrisParameters current_param  = json_current_run["params"];
      //fmt::print(stderr, "## WARNING DISABLED CHECK FOR PARAM\n");
      if(!(baseline_param == current_param)){
        fmt::print(stderr,"## ERROR: Parameters for test {} differ from baseline\n",test_name);
        baseline_param.printDifference(current_param,"baseline");
      }
      BOOST_REQUIRE(baseline_param == current_param);

      MeshStat baseline_stat = json_baseline_run["result"];
      MeshStat current_stat = json_current_run["result"];

      bolAr1 iclose_pctunit, iclose_avgquaP1, iclose_avgquaPk;
      bool icloseStat = isCloseMeshStat(baseline_stat, current_stat,
          iclose_pctunit, iclose_avgquaP1, iclose_avgquaPk);

      if(!icloseStat){
        fmt::print(stderr, "Mesh stats diverge for test {}:\n",test_name);
        printUnmatchedMeshStat(baseline_stat, current_stat, 
                               iclose_pctunit, iclose_avgquaP1, iclose_avgquaPk);
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

      if(!test_passes) fmt::print(stderr, "## Test {} failed.\n",test_name);
      else             fmt::print("-- Test {} passed.\n",test_name);

      if(!test_passes || !update_baseline) continue;

      updates_done = true;

      fmt::print("-- Updating baseline for test {}\n",test_name);

      // Update the baseline CPU and stats
      json_baseline_run["CPU_times"].push_back(current_time);
      // Update stat and metadata if this case improved.
      int tdim = current_stat.pctunit.get_n();
      if(tdim != baseline_stat.pctunit.get_n()){
        fmt::print(stderr,"## ERROR in updating baseline, inconsistent tdim\n");
        baseline_stat.print("baseline",stderr);
        current_stat.print("current",stderr);
        continue;
      }
      if(current_stat.pctunit[tdim-1] > baseline_stat.pctunit[tdim-1]){
        json_baseline_run["result"] = current_stat;
        // metadata is stored per run on the baseline, as each case 
        // might be updated at a different time.
        json_baseline_run["metadata"] = json_current["metadata"]; 
      }
    }

    if(updates_done){
      fmt::print("-- Updates done, overwriting baseline file {}\n",baseline_json_fname);
      std::fstream jfile(baseline_json_fname);
      jfile << json_baseline_allhwid.dump(2);
    }

  }
  
public:
// Provided
  std::string case_dir, out_dir, tmp_dir, baseline_json_fname;

// Internal
  HardwareID hwid;
  // Timings (process and wall) for full lifetime
  double t0_all, t0w_all, t1_all, t1w_all;
  bool update_baseline;
  std::string time_stamp;
  std::vector<std::string> test_names;
  nlohmann::json json_current;
};

} // namespace Metris

#endif



