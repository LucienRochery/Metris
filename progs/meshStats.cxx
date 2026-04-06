//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include <src/MetrisRunner/MetrisParameters.hxx>
#include <src/MetrisRunner/MetrisRunner.hxx>
#include <src/MetrisRunner/MeshStat.hxx>
#include <src/Mesh/Mesh.hxx>
#include <src/metris_options.hxx>
#include <src/aux_exceptions.hxx>
#include <egads.h>
#include <src/aux_EGADSprinterr.hxx>
#include <fstream>
#include "libmeshb.hxx"
#include <string>
#include <algorithm>

using namespace Metris;

template <class MFT>
void outputStats(MetrisRunner& run) {

  Mesh<MFT> &msh = (Mesh<MFT> &) *(run.msh_g);

  MeshStat stat;
  run.statMesh(0,&stat);

  std::string filename = "meshStatsOutput.txt";

  bool fileExists = std::filesystem::exists(filename);

  std::ofstream foutput;
  foutput.open(filename, std::ios::out | std::ios::app);

  if (!foutput.good()) METRIS_THROW_MSG("meshStats -- Error opening file!");

  if (!fileExists) {
    foutput
            << std::setw(30) << "npoin"
            << std::setw(30) << "nedge"
            << std::setw(30) << "nface"
            << std::setw(30) << "ntet"
            << std::setw(30) << "Length 1D min"
            << std::setw(30) << "Length 1D avg"
            << std::setw(30) << "Length 1D max"
            << std::setw(30) << "Length 1D % unit"
            << std::setw(30) << "Length 2D min"
            << std::setw(30) << "Length 2D avg"
            << std::setw(30) << "Length 2D max"
            << std::setw(30) << "Length 2D % unit"
            << std::setw(30) << "Length 3D min"
            << std::setw(30) << "Length 3D avg"
            << std::setw(30) << "Length 3D max"
            << std::setw(30) << "Length 3D % unit"
            << std::setw(30) << "Quality 2D min"
            << std::setw(30) << "Quality 2D avg"
            << std::setw(30) << "Quality 2D max"
            << std::setw(30) << "Quality 3D min"
            << std::setw(30) << "Quality 3D avg"
            << std::setw(30) << "Quality 3D max"
            << std::endl;
  }

  msh.cleanup();

  foutput
          << std::setw(30) << msh.npoin
          << std::setw(30) << msh.nentt(1)
          << std::setw(30) << msh.nentt(2)
          << std::setw(30) << msh.nentt(3)
          << std::setw(30) << stat.len[0].min()
          << std::setw(30) << stat.len[0].avg()
          << std::setw(30) << stat.len[0].max()
          << std::setw(30) << stat.pctunit[0]
          << std::setw(30) << stat.len[1].min()
          << std::setw(30) << stat.len[1].avg()
          << std::setw(30) << stat.len[1].max()
          << std::setw(30) << stat.pctunit[1]
          << std::setw(30) << stat.len[2].min()
          << std::setw(30) << stat.len[2].avg()
          << std::setw(30) << stat.len[2].max()
          << std::setw(30) << stat.pctunit[2]
          << std::setw(30) << stat.quaP1[1].min()
          << std::setw(30) << stat.quaP1[1].avg()
          << std::setw(30) << stat.quaP1[1].max()
          << std::setw(30) << stat.quaP1[2].min()
          << std::setw(30) << stat.quaP1[2].avg()
          << std::setw(30) << stat.quaP1[2].max()
          << std::endl;

  foutput.close();
}

int main(int argc, char** argv){
  MetrisOptions opt(argc, argv);
  MetrisParameters param(opt);

  MetrisRunner run(argc,argv);

  if (run.metricFE) outputStats<MetricFieldFE        >(run);
  else              outputStats<MetricFieldAnalytical>(run);

  return 0;
}

