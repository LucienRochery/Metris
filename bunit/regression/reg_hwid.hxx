//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_REG_HWID__
#define __METRIS_REG_HWID__

#include <string>
#include <thread>           // For std::thread::hardware_concurrency()
#include <cstdlib>          // For std::getenv

// Platform-specific includes
#ifdef __APPLE__
#include <sys/sysctl.h>
#include <sys/types.h>
#elif defined(__linux__)
#include <fstream>
#endif

namespace Metris {

#ifdef __APPLE__
std::string get_cpu_model(){
  char buffer[256];
  size_t size = sizeof(buffer);
  if(sysctlbyname("machdep.cpu.brand_string", &buffer, &size, NULL, 0) == 0){
    return std::string(buffer);
  }
  return "Unknown_Apple_CPU";
}


#elif defined(__linux__)
// Linux conveniently provides /proc/cpuinfo to cat.

std::string get_cpu_model(){
  std::ifstream cpuinfo("/proc/cpuinfo");
  std::string line;
  while(std::getline(cpuinfo, line)){
    // Line looks like:
    // model name      : Intel(R) Pentium(R) CPU G2030 @ 3.00GHz
    if(line.find("model name") == std::string::npos) continue;
    return line.substr(line.find(":") + 2);
  }
  return "Unknown_Linux_CPU";
}

#endif

// Hardware identification structure
struct HardwareID {
  std::string cpu_model;
  std::string cpu_architecture; 
  int num_cores;
  std::string os_name;

  HardwareID(){
    // CPU detection
    cpu_model = get_cpu_model();
    std::replace(cpu_model.begin(), cpu_model.end(), ' ', '_');

    #ifdef __APPLE__
    os_name = "macOS";
    #elif defined(__linux__)
    os_name = "Linux";
    #endif
    
    // Architecture detection
    #ifdef __x86_64__ 
    cpu_architecture = "x86_64";
    #elif defined(__aarch64__) || defined(__arm64__)
    cpu_architecture = "arm64";
    #elif defined(__i386__)
    cpu_architecture = "x86";
    #endif
    
    // Core count
    num_cores = std::thread::hardware_concurrency();
  }
  
  std::string to_string() const {
    const char* ci_hwid = std::getenv("METRIS_CI_HWID");
    if(ci_hwid) return std::string(ci_hwid);
    return fmt::format("{}_{}_{}c_{}",
                       cpu_model, cpu_architecture, num_cores, os_name);
  }
  
};

} // namespace Metris

#endif
