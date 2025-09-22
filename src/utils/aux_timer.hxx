//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __SRC_AUX_TIMER__
#define __SRC_AUX_TIMER__

#include <string>

// Sourced from stack https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows

namespace Metris{
	
double get_wall_time();
double get_cpu_time();

// Outputs a string YYYY-MM-DD-HH-mm-ss
std::string time2str();

} // End namespace
#endif