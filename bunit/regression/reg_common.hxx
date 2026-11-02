//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_REG_COMMON__
#define __METRIS_REG_COMMON__

// Regression test manager.
//
// One BOOST_AUTO_TEST_CASE = one regression case = one ctest entry of the same
// name (see CMakeLists.txt). A case runs several variants ("runs"). Each case
// writes its results to <out>/<case>.json, and its meshes and logs to
// <out>/<case>/, where <out> is $METRIS_REG_OUT, else <build>/regression/out.
//
// Environment:
//  - METRIS_REG_OUT=<dir>        where results go (see above).
//  - METRIS_REG_BASELINE=<dir>   compare against <dir>/<hwid>/<case>.json.
//                                Unset: record only, the case fails only if a
//                                run throws. Set and no baseline for this
//                                hwid/case: the case fails.
//  - METRIS_UPDATE_BASELINE=1    write the baseline: runs that pass, are new or
//                                whose parameters changed replace theirs; runs
//                                that fail or throw keep the old one.
//                                Into $METRIS_REG_BASELINE, else into
//                                bunit/regression/baselines.
//    METRIS_UPDATE_BASELINE=all  same, but failed runs are replaced too.
//  - METRIS_REG_SKIP_TIME=1      no CPU time check (CPU time still recorded).
//  - METRIS_CI_HWID=<name>       overrides the hardware ID.
//
// Comparison (one-sided: improvements pass) happens in finish(), which every
// case must call last. The same rules are in reg_compare.py, used by tools
// that compare two result sets; keep the tolerances in sync.

#include "../common_setup.hxx"
#include "reg_hwid.hxx"

#include "nlohmann/json.hpp"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <exception>
#include <map>
#include <set>
#include <memory>

#ifdef __linux__
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif

#ifndef METRIS_REG_DEFAULT_OUT
  #error "METRIS_REG_DEFAULT_OUT not defined, should be in CMakeLists.txt"
#endif
#ifndef METRIS_BUILD_TYPE
  #error "METRIS_BUILD_TYPE not defined, should be in CMakeLists.txt"
#endif

namespace Metris{

// Result file format version. Bump on incompatible changes, and in reg_compare.py.
constexpr int REG_SCHEMA = 1;

// One-sided tolerances: current may be at most this much worse than baseline
// (see regStatChecked for which way is better).
constexpr double REG_TOL_PCTUNIT = 1.0;     // % unit edges, in points
constexpr double REG_TOL_QUALITY = 1.0e-2;  // average quality, in (0,1)

// Number of CPU time samples kept per run in a baseline.
constexpr int REG_MAX_CPU_SAMPLES = 20;

// Parameters set by the manager itself or irrelevant to results: not compared.
inline bool regParamIgnored(const std::string& key){
  return key == "iverb" || key == "ivdepth";
}


inline std::string regGetEnv(const char* name){
  const char* val = std::getenv(name);
  return val ? std::string(val) : std::string();
}

// Output of a shell command, trailing newline removed. ok = false on failure.
inline std::string regShell(const std::string& cmd, bool& ok){
  ok = false;
  FILE* pipe = popen(cmd.c_str(), "r");
  if(!pipe) return "";
  std::string out;
  char buf[256];
  while(fgets(buf, sizeof(buf), pipe)) out += buf;
  ok = (pclose(pipe) == 0);
  while(!out.empty() && (out.back() == '\n' || out.back() == '\r')) out.pop_back();
  return out;
}


// Instructions retired in user space by this process and the threads it
// starts after start(). Linux only; elsewhere, or when perf_event_open is not
// allowed (perf_event_paranoid > 2, most CI VMs), available() is false.
class RegInstructionCounter{
public:
  RegInstructionCounter(){
#ifdef __linux__
    perf_event_attr attr{};
    attr.type = PERF_TYPE_HARDWARE;
    attr.size = sizeof(attr);
    attr.config = PERF_COUNT_HW_INSTRUCTIONS;
    attr.disabled = 1;
    attr.exclude_kernel = 1;
    attr.exclude_hv = 1;
    attr.inherit = 1;
    fd = syscall(SYS_perf_event_open, &attr, 0, -1, -1, 0);
#endif
  }
  ~RegInstructionCounter(){
#ifdef __linux__
    if(fd >= 0) close(fd);
#endif
  }
  RegInstructionCounter(const RegInstructionCounter&) = delete;
  RegInstructionCounter& operator=(const RegInstructionCounter&) = delete;

  bool available() const { return fd >= 0; }

  void start(){
#ifdef __linux__
    if(fd < 0) return;
    ioctl(fd, PERF_EVENT_IOC_RESET, 0);
    ioctl(fd, PERF_EVENT_IOC_ENABLE, 0);
#endif
  }
  // Returns the count since start(), or -1 if not available.
  long long stop(){
#ifdef __linux__
    if(fd < 0) return -1;
    ioctl(fd, PERF_EVENT_IOC_DISABLE, 0);
    long long count = 0;
    if(read(fd, &count, sizeof(count)) != sizeof(count)) return -1;
    return count;
#else
    return -1;
#endif
  }
private:
  int fd = -1;
};


// Upper bound on CPU time from a baseline's samples.
inline double regCPUBound(const nlohmann::json& samples){
  int n = samples.size();
  double avg = 0.0;
  for(double t : samples) avg += t;
  avg /= n;
  // Fewer than 5 samples: 1.4x the average.
  if(n < 5) return 1.4*avg;
  // Otherwise 2.5 sigma, a little under 99% if normally distributed.
  double var = 0.0;
  for(double t : samples) var += (t - avg)*(t - avg);
  return avg + 2.5*std::sqrt(var/n);
}


// The checked statistics. sign = +1: higher is better (% unit edges);
// -1: lower is better (quality: the stats hold the conformity error |q - 1|,
// see getmetquamesh, so 0 is ideal).
inline bool regStatChecked(const std::string& key, double& tol, int& sign){
  if(key.rfind("pctunit", 0) == 0){ tol = REG_TOL_PCTUNIT; sign = 1; return true; }
  if(key.rfind("avgquaP1_", 0) == 0 || key.rfind("avgquaPk_", 0) == 0){
    tol = REG_TOL_QUALITY; sign = -1; return true;
  }
  return false;
}


/* --------------------------------------------------------------- */
class RegressionTestManager{
public:
  RegressionTestManager() = delete;
  RegressionTestManager(const RegressionTestManager&) = delete;
  RegressionTestManager& operator=(const RegressionTestManager&) = delete;

  // case_subdir: where the inputs are, relative to bunit/regression/, e.g. "2D/OAT15A".
  // The case name is the Boost test case's.
  explicit RegressionTestManager(std::string case_subdir){
    if(std::string(METRIS_ROOT_DIR).empty()){
      METRIS_THROW_MSG("METRIS_ROOT_DIR not defined, should be in CMakeLists.txt");
    }

    case_name = boost::unit_test::framework::current_test_case().p_name;

    t0_all  = get_cpu_time();
    t0w_all = get_wall_time();

    in_dir = std::string(METRIS_ROOT_DIR) + "/bunit/regression/" + case_subdir + "/";

    std::string out_root = regGetEnv("METRIS_REG_OUT");
    if(out_root.empty()) out_root = METRIS_REG_DEFAULT_OUT;
    out_root += "/";
    result_fname = out_root + case_name + ".json";
    out_dir = out_root + case_name + "/";
    tmp_dir = out_dir + "tmp/";
    // Results of an earlier run of this case must not mix with this one's.
    std::filesystem::remove_all(out_dir);
    std::filesystem::remove(result_fname);
    std::filesystem::create_directories(tmp_dir);

    baseline_root = regGetEnv("METRIS_REG_BASELINE");
    std::string upd = regGetEnv("METRIS_UPDATE_BASELINE");
    update_baseline = !upd.empty() && upd != "0";
    update_all      = (upd == "all");
    skip_time       = !regGetEnv("METRIS_REG_SKIP_TIME").empty()
                   &&  regGetEnv("METRIS_REG_SKIP_TIME") != "0";

    json_current["schema"]     = REG_SCHEMA;
    json_current["case"]       = case_name;
    json_current["hwid"]       = hwid.to_string();
    json_current["build_type"] = METRIS_BUILD_TYPE;
    setGitInfo();
    json_current["time_stamp"] = time2str();
    json_current["runs"]       = nlohmann::json::object();

    fmt::print("-- Case {} on {} ({}), commit {}{}\n", case_name, hwid.to_string(),
               METRIS_BUILD_TYPE, json_current["commit"].get<std::string>(),
               json_current["dirty"].get<bool>() ? " (dirty)" : "");
    fmt::print("-- Results: {}\n", result_fname);
    if(!counter.available())
      fmt::print("-- Instruction counts not available (perf_event_open).\n");
  }

  std::string getMeshIn(const std::string &msh_base_name) const {
    return in_dir + msh_base_name;
  }
  std::string getMetricIn(const std::string &met_base_name) const {
    if(met_base_name.empty()) return "";
    return in_dir + met_base_name;
  }
  std::string getCADIn(const std::string &CAD_base_name) const {
    if(CAD_base_name.empty()) return "";
    return in_dir + CAD_base_name;
  }
  std::string getOutDir() const {
    return out_dir;
  }

  // Run from files. CAD_base_name and met_base_name can be "".
  void runTest(MetrisParameters& param,
               std::string test_name,
               std::string mesh_base_name,
               std::string CAD_base_name,
               std::string met_base_name){
    param.setMeshIn(getMeshIn(mesh_base_name));
    if(!CAD_base_name.empty()) param.setCAD(getCADIn(CAD_base_name));
    if(!met_base_name.empty()) param.setMetricFile(getMetricIn(met_base_name));
    setOutputs(param, test_name);
    runGuarded(test_name, param, [&](){ return std::make_unique<MetrisRunner>(nullptr, param); });
  }

  // Run from existing MetrisAPI data (e.g. after some preprocessing).
  // Example: reg_2D_OAT15A_x2.cxx, reg_2D_OAT15A_x2_Q2toQ2.
  void runTest(MetrisParameters& param,
               std::string test_name,
               MetrisAPI *data_front,
               MetrisAPI *data_back,
               std::string CAD_base_name,
               std::string met_base_name){
    if(!CAD_base_name.empty()) param.setCAD(getCADIn(CAD_base_name));
    // Don't set metric file when data_front already carries a metric: the
    // constructor promotes data_front to data_back when data_back is NULL,
    // and MeshBack rejects having both data->imet and param.inpMet set.
    auto hasmet = [](MetrisAPI *d) -> bool {
      if(!d) return false;
      int idim,ideg,ncorn,npoin,nedge,nface,nelem; bool imet; FEBasis mb,mtb;
      d->getConstants(&idim,&ideg,&ncorn,&npoin,&imet,&nedge,&nface,&nelem,&mb,&mtb);
      return imet;
    };
    if(!met_base_name.empty() && !hasmet(data_front) && !hasmet(data_back))
      param.setMetricFile(getMetricIn(met_base_name));
    setOutputs(param, test_name);
    runGuarded(test_name, param, [&](){
      return std::make_unique<MetrisRunner>(data_front, data_back, param); });
  }

  // Writes the result file, then compares against the baseline if asked and
  // writes the baseline if asked. Call last in every case.
  void finish(){
    if(finished) return;
    finished = true;

    double t1_all = get_cpu_time(), t1w_all = get_wall_time();
    fmt::print("\n========================================\n");
    fmt::print("-- Case {} done: wall = {:.2f}s, cpu = {:.2f}s\n",
               case_name, t1w_all - t0w_all, t1_all - t0_all);

    writeResult();

    for(auto& [name, run] : json_current["runs"].items()){
      if(run["status"] != "ok"){
        BOOST_ERROR("Run " << name << " raised an exception: "
                    << run["error"].get<std::string>());
      }
    }

    nlohmann::json json_baseline;
    bool have_baseline = false;
    std::string why;
    if(!baseline_root.empty()){
      std::string fname = baselineFile(baseline_root);
      have_baseline = readJSON(fname, json_baseline);
      if(have_baseline && !validBaseline(json_baseline, why)){
        BOOST_ERROR("Baseline " << fname << " is malformed (" << why << ")");
        have_baseline = false;
      }else if(have_baseline){
        fmt::print("-- Comparing against baseline {}\n", fname);
        compare(json_baseline);
      }else{
        BOOST_ERROR("No baseline " << fname << " for hwid " << hwid.to_string()
                    << ": nothing to compare against");
      }
    }else{
      fmt::print("-- METRIS_REG_BASELINE not set: results recorded, not compared.\n");
      // Updating without comparing: still start from the existing baseline,
      // so CPU samples and the runs that threw this time are kept.
      if(update_baseline){
        std::string fname = baselineFile(baselineWriteRoot());
        have_baseline = readJSON(fname, json_baseline) && validBaseline(json_baseline, why);
        if(!have_baseline && std::filesystem::exists(fname))
          fmt::print(stderr, "## Existing baseline {} unreadable ({}): replaced\n", fname, why);
      }
    }

    if(update_baseline) writeBaseline(have_baseline ? json_baseline : nlohmann::json());
  }

  ~RegressionTestManager(){
    // Comparing here is not an option: a failed BOOST_REQUIRE throws, and
    // throwing from a destructor terminates. Only make sure the results are
    // on disk and flag the missing finish().
    if(finished) return;
    try{
      writeResult();
      if(!std::uncaught_exceptions())
        BOOST_ERROR("RegressionTestManager::finish() was not called for case " << case_name);
    }catch(...){}
  }

private:
  void setGitInfo(){
    bool ok;
    std::string root = METRIS_ROOT_DIR;
    std::string commit = regShell("git --no-optional-locks -C '" + root + "' rev-parse HEAD 2>/dev/null", ok);
    if(ok && !commit.empty()){
      std::string st = regShell("git --no-optional-locks -C '" + root + "' status --porcelain "
                                "--untracked-files=no 2>/dev/null", ok);
      json_current["commit"] = commit;
      // Unknown status counts as dirty: such a run is never logged as clean.
      json_current["dirty"]  = !ok || !st.empty();
      return;
    }
    // No git (e.g. CI container checkout without .git).
    std::string sha = regGetEnv("GITHUB_SHA");
    json_current["commit"] = !sha.empty() ? sha : std::string(METRIS_GIT_COMMIT_HASH);
    json_current["dirty"]  = sha.empty();
  }

  void setOutputs(MetrisParameters& param, const std::string& test_name){
    std::string outfile = out_dir + test_name;
    param.iverb = 2;
    param.ivdepth = 1;
    param.outmPrefix = tmp_dir;
    param.setLogFile(outfile + ".log");
    param.setMeshOut(outfile);
  }

  // Non-default parameters, file names reduced to their base name.
  static nlohmann::json paramsNonDefault(const MetrisParameters& param){
    nlohmann::json cur = param;
    MetrisParameters defparam;
    nlohmann::json def = defparam;
    nlohmann::json out = nlohmann::json::object();
    for(auto& [key, val] : cur.items()){
      if(regParamIgnored(key)) continue;
      nlohmann::json v = val;
      if(v.is_string()){
        std::string s = v.get<std::string>();
        v = s.substr(s.find_last_of("/\\") + 1);
      }
      if(def.contains(key) && def[key] == val) continue;
      out[key] = v;
    }
    return out;
  }

  template<typename MAKE_RUNNER>
  void runGuarded(const std::string& test_name, MetrisParameters& param, MAKE_RUNNER make_runner){
    fmt::print("\n========================================\n");
    fmt::print("Running {}/{}\n", case_name, test_name);

    if(json_current["runs"].contains(test_name)){
      BOOST_ERROR("Run name " << test_name << " used twice in case " << case_name);
    }

    nlohmann::json entry;
    entry["status"]            = "ok";
    entry["error"]             = nullptr;
    entry["cpu_s"]             = nullptr;
    entry["wall_s"]            = nullptr;
    entry["instructions"]      = nullptr;
    entry["params_nondefault"] = paramsNonDefault(param);
    entry["stats"]             = nullptr;

    double t0 = get_cpu_time(), t0w = get_wall_time();
    counter.start();
    try{
      auto run = make_runner();
      run->runMetris();
      long long ninstr = counter.stop();
      double t1 = get_cpu_time(), t1w = get_wall_time();
      fmt::print("Total time: wall = {:.2f}s, cpu = {:.2f}s\n", t1w - t0w, t1 - t0);
      entry["cpu_s"]  = t1 - t0;
      entry["wall_s"] = t1w - t0w;
      if(ninstr >= 0) entry["instructions"] = ninstr;
      MeshStat stat;
      run->statMesh(0, &stat);
      stat.print();
      entry["stats"] = stat;
    }catch(const std::exception& e){
      counter.stop();
      entry["status"] = "exception";
      entry["error"]  = std::string(e.what());
      fmt::print(stderr, "## Run {} raised exception:\n{}\n", test_name, e.what());
    }
    json_current["runs"][test_name] = entry;
  }

  // Serializes first, then writes a temporary file and renames it, so a
  // failure never leaves a truncated or empty file behind. Invalid UTF-8
  // (e.g. in an exception message) is replaced rather than thrown on.
  static void writeJSONFile(const std::string& fname, const nlohmann::json& jj){
    std::string text = jj.dump(2, ' ', false, nlohmann::json::error_handler_t::replace) + "\n";
    std::string tmp = fname + ".tmp";
    {
      std::ofstream ff(tmp);
      ff << text;
      if(!ff) throw std::runtime_error("Could not write " + tmp);
    }
    std::filesystem::rename(tmp, fname);
  }

  void writeResult() const {
    writeJSONFile(result_fname, json_current);
    fmt::print("-- Results written to {}\n", result_fname);
  }

  std::string baselineWriteRoot() const {
    return baseline_root.empty()
         ? std::string(METRIS_ROOT_DIR) + "/bunit/regression/baselines"
         : baseline_root;
  }

  // Structure checks before any field is read: nlohmann's const operator[]
  // on a missing key is undefined behaviour in Release.
  static bool validBaseline(const nlohmann::json& jj, std::string& why){
    why.clear();
    if(!jj.is_object())                                   why = "not an object";
    else if(!jj.contains("schema") || jj["schema"] != REG_SCHEMA)
      why = "schema is not " + std::to_string(REG_SCHEMA);
    else if(!jj.contains("runs") || !jj["runs"].is_object()) why = "no runs object";
    else for(auto& [name, run] : jj["runs"].items()){
      if(!run.is_object() || !run.contains("status") || !run["status"].is_string()){
        why = "run " + name + ": no status"; break;
      }
      if(run["status"] != "ok") continue;
      if(!run.contains("params_nondefault") || !run["params_nondefault"].is_object()
      || !run.contains("stats") || !run["stats"].is_object()){
        why = "run " + name + ": no params_nondefault or stats object"; break;
      }
      if(run.contains("cpu_s_samples")){
        bool okk = run["cpu_s_samples"].is_array();
        if(okk) for(auto& t : run["cpu_s_samples"]) okk = okk && t.is_number();
        if(!okk){ why = "run " + name + ": cpu_s_samples not an array of numbers"; break; }
      }
    }
    return why.empty();
  }

  std::string baselineFile(const std::string& root) const {
    return root + "/" + hwid.to_string() + "/" + case_name + ".json";
  }

  static bool readJSON(const std::string& fname, nlohmann::json& jj){
    if(!std::filesystem::exists(fname)) return false;
    try{
      std::ifstream ff(fname);
      ff >> jj;
      return true;
    }catch(const std::exception& e){
      fmt::print(stderr, "## Could not read {}: {}\n", fname, e.what());
      return false;
    }
  }

  // Compares runs against the baseline; records which runs passed in run_passed.
  void compare(const nlohmann::json& json_baseline){
    // Structure already checked by validBaseline.
    const nlohmann::json& base_runs = json_baseline["runs"];

    for(auto& [name, base] : base_runs.items()){
      if(!json_current["runs"].contains(name)){
        BOOST_ERROR("Run " << name << " is in the baseline but was not run");
      }
    }

    for(auto& [name, cur] : json_current["runs"].items()){
      fmt::print("========================================\n");
      fmt::print("Comparing {}\n", name);
      run_passed[name] = false;
      if(cur["status"] != "ok") continue; // already reported
      if(!base_runs.contains(name)){
        fmt::print("-- {} not in the baseline (new run): not compared.\n", name);
        run_passed[name] = true;
        continue;
      }
      const nlohmann::json& base = base_runs[name];
      bool pass = true;

      if(base["params_nondefault"] != cur["params_nondefault"]){
        fmt::print(stderr, "## {}: non-default parameters differ from baseline\n"
                           "   baseline {}\n   current  {}\n", name,
                   base["params_nondefault"].dump(), cur["params_nondefault"].dump());
        BOOST_ERROR("Run " << name << ": parameters differ from baseline");
        params_changed.insert(name);
        continue;
      }

      for(auto& [key, bval] : base["stats"].items()){
        double tol;
        int sign;
        if(!regStatChecked(key, tol, sign)) continue;
        if(!cur["stats"].contains(key)){
          BOOST_ERROR("Run " << name << ": stat " << key << " missing");
          pass = false;
          continue;
        }
        if(!bval.is_number()){
          // NaN in the baseline (written as null): nothing to compare against.
          fmt::print(" - {}: baseline {} is not a number, not compared\n", name, key);
          continue;
        }
        double b = bval, c = cur["stats"][key];
        if(!std::isfinite(c)){
          fmt::print(stderr, " - {}: {} is {}\n", name, key, c);
          BOOST_ERROR("Run " << name << ": " << key << " is not finite");
          pass = false;
          continue;
        }
        if(sign*(c - b) < -tol){
          fmt::print(stderr, " - {}: {} worse: baseline {:.4g}, current {:.4g} (tol {})\n",
                     name, key, b, c, tol);
          BOOST_ERROR("Run " << name << ": " << key << " " << c
                      << " worse than baseline " << b);
          pass = false;
        }
      }

      if(!skip_time && base.contains("cpu_s_samples") && !base["cpu_s_samples"].empty()){
        double bound = regCPUBound(base["cpu_s_samples"]);
        double tt = cur["cpu_s"];
        if(tt > bound){
          fmt::print(stderr, " - {}: CPU time {:.2f}s > {:.2f}s baseline bound\n", name, tt, bound);
          BOOST_ERROR("Run " << name << ": CPU time " << tt
                      << " s over the baseline bound " << bound << " s");
          pass = false;
        }else{
          fmt::print("-- {}: CPU time {:.2f}s <= {:.2f}s baseline bound\n", name, tt, bound);
        }
      }

      run_passed[name] = pass;
      fmt::print("-- {} {}\n", name, pass ? "passed" : "FAILED");
    }
  }

  void writeBaseline(nlohmann::json json_baseline){
    std::string fname = baselineFile(baselineWriteRoot());
    std::string why;
    if(!validBaseline(json_baseline, why)){
      json_baseline = nlohmann::json::object();
      json_baseline["runs"] = nlohmann::json::object();
    }
    nlohmann::json runs = json_baseline["runs"];
    json_baseline = json_current;
    json_baseline.erase("dirty");
    json_baseline.erase("time_stamp");
    json_baseline["runs"] = nlohmann::json::object();

    for(auto& [name, cur] : json_current["runs"].items()){
      bool has_old = runs.contains(name);
      bool passed  = !run_passed.count(name) || run_passed[name];
      // Changed parameters mean the case itself was changed: take the new run.
      bool replace = passed || update_all || params_changed.count(name);
      if(cur["status"] != "ok" || (has_old && !replace)){
        if(has_old) json_baseline["runs"][name] = runs[name];
        fmt::print("-- Baseline for {} {}\n", name, has_old ? "kept" : "not written");
        continue;
      }
      nlohmann::json samples = nlohmann::json::array();
      // Samples only carry over while the parameters are the same.
      if(has_old && runs[name].contains("cpu_s_samples")
         && runs[name]["params_nondefault"] == cur["params_nondefault"])
        samples = runs[name]["cpu_s_samples"];
      samples.push_back(cur["cpu_s"]);
      while((int)samples.size() > REG_MAX_CPU_SAMPLES) samples.erase(samples.begin());
      nlohmann::json entry = cur;
      entry["cpu_s_samples"] = samples;
      entry["commit"] = json_current["commit"];
      json_baseline["runs"][name] = entry;
      fmt::print("-- Baseline for {} updated\n", name);
    }
    // Runs that are in the old baseline but no longer run are dropped.

    std::filesystem::create_directories(std::filesystem::path(fname).parent_path());
    writeJSONFile(fname, json_baseline);
    fmt::print("-- Baseline written to {}\n", fname);
  }

public:
  std::string case_name, in_dir, out_dir, tmp_dir, result_fname, baseline_root;

private:
  HardwareID hwid;
  RegInstructionCounter counter;
  double t0_all, t0w_all;
  bool update_baseline, update_all, skip_time;
  bool finished = false;
  nlohmann::json json_current;
  std::map<std::string, bool> run_passed;
  std::set<std::string> params_changed;
};

} // namespace Metris

#endif
