//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../src/Metris.h"
#include "../src/libmeshb.hxx"
#include <fstream>
#include <string>
#include <algorithm>

#include <fstream>

using namespace Metris;


std::string ReplaceAll(std::string str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}


int main(int argc, char** argv){
  MetrisOptions opt(argc, argv);
  MetrisParameters param(opt);

  MetrisRunner run(argc,argv);

  Mesh<MetricFieldFE> &msh = (Mesh<MetricFieldFE> &) *(run.msh_g);
  GETVDEPTH(msh.param);

  if(!opt.count("rflag1") || !opt.count("rflag2") || !(opt.count("rflag3") && msh.idim == 3)){
    PRINTF("## Error: use rflag1, rflag2 and rflag3 to specify dx, dy, dz.\n");
    return 1;
  }

  double dx = msh.param->rflag1;
  double dy = msh.param->rflag2;
  double dz = msh.param->rflag3;

  PRINTF("Using dx = {}, dy = {}",dx,dy);
  if(msh.idim == 3) PRINTF(", dz = {}",dz);
  PRINTF("\n");

  for(int ipoin = 0; ipoin < msh.npoin; ipoin++){
    msh.coord(ipoin,0) += dx;
    msh.coord(ipoin,1) += dy;
    if(msh.idim == 3) msh.coord(ipoin,2) += dz;
  }

  std::string fname = param.meshFileName;
  fname = ReplaceAll(fname, ".meshb", "");
  fname = ReplaceAll(fname, ".mesh", "");
  fname += ".offset.meshb";

  writeMesh(fname, msh, true);

  return 0;
}