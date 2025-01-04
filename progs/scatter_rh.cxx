//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../src/Metris.h"
#include "../libs/libmeshb.hxx"
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

  dblAr1 array_r(msh.nface);
  array_r.set_n(msh.nface);
  dblAr1 array_h(msh.nface);
  array_h.set_n(msh.nface);

  for(int iface = 0; iface < msh.nface; iface++){
    double coop[2] = {0};
    for(int ii = 0; ii < 2; ii++)
      coop[ii] += msh.coord(msh.fac2poi(iface,0),ii)
                + msh.coord(msh.fac2poi(iface,1),ii)
                + msh.coord(msh.fac2poi(iface,2),ii);

    double r = sqrt(getnrml2<2>(coop));

    double h = getmeasentP1<2>(msh.fac2poi[iface], msh.coord);
    h = sqrt(h);

    array_r[iface] = r;
    array_h[iface] = h;
  }

  std::string fname = param.meshFileName;
  fname = ReplaceAll(fname, ".meshb", "");
  fname = ReplaceAll(fname, ".mesh", "");


  std::ofstream myfile;
  myfile.open(fname + "_rh.csv");
  myfile << "r,h\n";
  for(int iface = 0; iface < msh.nface; iface++){
    myfile << array_r[iface] << "," << array_h[iface] <<"\n";
  }
  myfile.close();
  return 0;
}