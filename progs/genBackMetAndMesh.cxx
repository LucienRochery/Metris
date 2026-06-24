//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include <src/MetrisRunner/MetrisParameters.hxx>
#include <src/MetrisRunner/MetrisRunner.hxx>
#include <src/MetrisRunner/MeshStat.hxx>
#include <src/msh_metricCost.hxx>
#include <src/low_topo.hxx>
#include <src/Mesh/Mesh.hxx>
#include <src/metris_options.hxx>
#include <src/aux_exceptions.hxx>
#include <egads.h>
#include <src/aux_EGADSprinterr.hxx>
#include <fstream>
#include "libmeshb.hxx"
#include <string>
#include <algorithm>

#include "../src/io_libmeshb.hxx"

using namespace Metris;

int main(int argc, char** argv){

  MetrisOptions opt(argc, argv);
  MetrisParameters param(opt);

  MetrisRunner run(argc,argv);

  if (run.metricFE){
    // print mesh and metric
    Mesh<MetricFieldFE>& msh = static_cast<Mesh<MetricFieldFE>&>(*(run.msh_g));

    writeMesh("mshForBack.meshb",msh);
    msh.met.writeMetricFile("metInMeshForBack.solb");
  }
  else{
    // print mesh and metric
    Mesh<MetricFieldAnalytical>& msh = static_cast<Mesh<MetricFieldAnalytical>&>(*(run.msh_g));

    writeMesh("mshForBack.meshb",msh);
    msh.met.writeMetricFile("metInMeshForBack.solb");
  }

  return 0;
}