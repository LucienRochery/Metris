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



} // namespace Metris

#endif



