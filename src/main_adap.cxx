//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "MetrisRunner/MetrisRunner.hxx"


namespace Metris{

int main_metris(int argc, char** argv){ 

  MetrisRunner run(argc,argv);

  run.runMetris();

  return 0;
}

} // End namespace

