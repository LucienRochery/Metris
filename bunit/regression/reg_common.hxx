//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_REG_COMMON__
#define __METRIS_REG_COMMON__

#include "../common_setup.hxx"

#include "nlohmann/json.hpp"

#include <filesystem>
#include <chrono>
#include <iomanip>
#include <string>
#include <sstream>

namespace Metris{


// Outputs a string YYYY-MM-DD-HH-MM-SS
std::string time2str() {
  auto now = std::chrono::system_clock::now();
  auto time_t = std::chrono::system_clock::to_time_t(now);
  
  std::ostringstream oss;
  oss << std::put_time(std::localtime(&time_t), "%Y-%m-%d-%H-%M-%S");
  return oss.str();
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



