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

using namespace Metris;

template <class MFT>
void outputStats(MetrisRunner& run) {

  Mesh<MFT> &msh = (Mesh<MFT> &) *(run.msh_g);

  double actualCost = 0; // the actual cost is the integral of sqrt(det(M_K)) where M_K is the mesh-implied metric. that is equal to (n. of elements) * (volume of ideal element)
  double integratedCost = 0; // the integrated cost refers to sqrt(det(M)) where M is the target metric

  msh.cleanup();

  int tdim = msh.get_tdim();

  actualCost  = (double)msh.nentt(tdim);
  if (tdim == 2){
    actualCost *= sqrt(3.)/4.;
  }else{
    actualCost *= sqrt(2.)/12.;
  }
  if (tdim == 2) integratedCost = getMetricCost<MFT,2,2>(msh);
  else           integratedCost = getMetricCost<MFT,3,3>(msh);

  double avgBalance = 0;

  const int mball = 100;
  intAr1 lball(mball);

  for (int ipoin = 0; ipoin < msh.npoin; ipoin++){

    constexpr int ideg = 1;

    int ientt = getpoient(msh, ipoin, tdim);
    int iver = tdim == 2  ? msh.template getverfac<ideg>(ientt, ipoin)
                          : msh.template getvertet<ideg>(ientt, ipoin);
    METRIS_ENFORCE(iver >= 0);

    int ierro = 0;
    int iopen;
    if (msh.get_tdim() == 2){
      intAr1 dum;
      bool imani = false;
      ierro = ball2(msh,ipoin,ientt,lball,dum,&iopen,&imani,0);
      METRIS_ENFORCE(imani == true);
    }else{
      ierro = ball3(msh,ipoin,ientt,lball,&iopen,0);
    }
    METRIS_ENFORCE(ierro == 0);

    int nenttball = lball.get_n();
    avgBalance += nenttball;

  }
  avgBalance /= msh.npoin;

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
            << std::setw(30) << "Target Cost (integrated)"
            << std::setw(30) << "Actual Cost"
            << std::setw(30) << "Avg. Balance"
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
          << std::setw(30) << integratedCost
          << std::setw(30) << actualCost
          << std::setw(30) << avgBalance
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

