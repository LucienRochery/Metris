//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "aux_timer.hxx"
// Sourced from stack https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>   // For getrusage()

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

} // End namespace
