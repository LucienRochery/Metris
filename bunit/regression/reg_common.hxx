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
#include <iomanip>
#include <string>
#include <sstream>

namespace Metris{

bool isGitDirty() {
  FILE* pipe = popen("git status --porcelain 2>/dev/null", "r");
  if (!pipe) return true;
  
  int c = fgetc(pipe);
  pclose(pipe);
  
  return (c != EOF);
}

// Outputs a string YYYY-MM-DD-HH-MM-SS
std::string time2str() {
  auto now = std::chrono::system_clock::now();
  auto time_t = std::chrono::system_clock::to_time_t(now);
  
  std::ostringstream oss;
  oss << std::put_time(std::localtime(&time_t), "%Y-%m-%d-%H-%M-%S");
  return oss.str();
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

bool isCloseMeshStat(const MeshStat& baseline, const MeshStat& current,
    bool* iclose_pctunit, bool* iclose_pctunit_bdry,
    bool* iclose_avgqua, bool* iclose_avgqua_bdry){
  *iclose_pctunit      = abs(baseline.pctunit      - current.pctunit     ) < 1.0;
  *iclose_pctunit_bdry = abs(baseline.pctunit_bdry - current.pctunit_bdry) < 1.0;
  *iclose_avgqua       = abs(baseline.avgqua       - current.avgqua      ) < 1.0e-2;
  *iclose_avgqua_bdry  = abs(baseline.avgqua_bdry  - current.avgqua_bdry ) < 1.0e-2;
  return *iclose_pctunit && *iclose_pctunit_bdry && *iclose_avgqua && *iclose_avgqua_bdry;
}

void printUnmatchedMeshStat(const MeshStat& baseline, const MeshStat& current,
    bool iclose_pctunit, bool iclose_pctunit_bdry,
    bool iclose_avgqua, bool iclose_avgqua_bdry){

  if(!iclose_pctunit) fmt::print(stderr," - pctunit differ: baseline {} vs current {}\n",
      baseline.pctunit, current.pctunit);
  if(!iclose_pctunit_bdry) fmt::print(stderr," - pctunit_bdry differ: baseline {} vs current {}\n",
      baseline.pctunit_bdry, current.pctunit_bdry);
  if(!iclose_avgqua) fmt::print(stderr," - avg quality differ: baseline {} vs current {}\n",
      baseline.avgqua, current.avgqua);
  if(!iclose_avgqua_bdry) fmt::print(stderr," - avg quality bdry differ: baseline {} vs current {}\n",
      baseline.avgqua_bdry, current.avgqua_bdry);
  }

} // namespace Metris

#endif



