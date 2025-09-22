//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "Metris.h"
#include "libmeshb.hxx"
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

  #if 0
  msh.set_nelem(0);

  for(int iface = 0; iface < msh.nface; iface++){
    if(isdeadent(iface, msh.fac2poi)) continue;
    msh.fac2tet(iface,0) = -1;
    msh.fac2tet(iface,1) = -1;
    for(int ii = 0; ii < msh.nnode(2); ii++){
      int ipoin = msh.fac2poi(iface, ii);
      msh.poi2ent(ipoin, 0) = iface;
      msh.poi2ent(ipoin, 1) = 2;
    }
  }

  for(int iedge = 0; iedge < msh.nface; iedge++){
    if(isdeadent(iedge, msh.fac2poi)) continue;
    msh.fac2tet(iedge,0) = -1;
    msh.fac2tet(iedge,1) = -1;
    for(int ii = 0; ii < msh.nnode(2); ii++){
      int ipoin = msh.fac2poi(iedge, ii);
      msh.poi2ent(ipoin, 0) = iedge;
      msh.poi2ent(ipoin, 1) = 1;
    }
  }
  #endif

  std::string fname = param.meshFileName;
  fname = ReplaceAll(fname, ".meshb", "");
  fname = ReplaceAll(fname, ".mesh", "");
  fname += ".surf.meshb";

  writeMesh(fname, msh, true, 0, 0, msh.nelem);

  return 0;
}