//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "aux_timer.hxx"
// Sourced from stack https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>   // For getrusage()
#include <sstream>
#include <iomanip>
#include <chrono>

namespace Metris{


double get_wall_time(){
  struct timeval time;
  if (gettimeofday(&time,NULL)){
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

// Process CPU time, similar to time command "user + system" value
double get_cpu_time(){
  struct rusage usage;
  if(getrusage(RUSAGE_SELF, &usage) == 0){
    return usage.ru_utime.tv_sec + usage.ru_utime.tv_usec * 1e-6 +
           usage.ru_stime.tv_sec + usage.ru_stime.tv_usec * 1e-6;
  }
  return 0.0;
}


// Outputs a string YYYY-MM-DD-HH-mm-ss
std::string time2str(){
  auto now = std::chrono::system_clock::now();
  auto time_t = std::chrono::system_clock::to_time_t(now);
  
  std::ostringstream oss;
  oss << std::put_time(std::localtime(&time_t), "%Y-%m-%d-%H-%M-%S");
  return oss.str();
}

} // End namespace
