//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "MetrisRunner.hxx"

#include "../Mesh/Mesh.hxx"

#include "../utils/mprintf.hxx"
#include "../utils/aux_timer.hxx"

#include "../io_libmeshb.hxx"
#include "../msh_checktopo.hxx"
#include "../msh_metricCost.hxx"

#ifdef WRITEQUALFIELD
#include "../quality/quafun.hxx"
#include "../quality/low_metqua.hxx"
#include "../low_topo.hxx"
#endif

namespace Metris{

void MetrisRunner::runMetris(){

  std::string outputAdaptStats = param->outmPrefix;
  if (param->MOESS_adapt_it >= 0){

    outputAdaptStats += "outputAdaptStats_MOESS_a" + std::to_string(param->MOESS_adapt_it) + ".txt";

    std::string meshName = "mesh_MOESS_initial_a" + std::to_string(param->MOESS_adapt_it) + ".meshb";
    std::string metName  = "met_MOESS_initial_a"  + std::to_string(param->MOESS_adapt_it) + ".solb";

    writeMesh(meshName,*msh_g);
    Mesh<MetricFieldFE> &msh = (Mesh<MetricFieldFE> &)(*msh_g);
    msh.met.writeMetricFile(metName);
  }
  else{
    outputAdaptStats += "outputAdaptStats.txt";
  }

  foutputAdaptStats.open(outputAdaptStats, std::ios::app);
  METRIS_ASSERT_MSG(foutputAdaptStats.good(), "Error opening file: " + outputAdaptStats);

  foutputAdaptStats << std::setw(30) << "iter"
                    << std::setw(30) << "nTrySmoo"
                    << std::setw(30) << "nSuccSmoo"
                    << std::setw(30) << "nTryIns"
                    << std::setw(30) << "nSuccIns"
                    << std::setw(30) << "nTryCol"
                    << std::setw(30) << "nSuccCol"
                    << std::endl;

  GETVDEPTH(param);

  double t0, t1;
  t0 = get_cpu_time();

  try{

    if(DOPRINTS1()) statMesh();

    if(param->dbgfull) check_topo(*msh_g,0);

    //if(param->opt_unif){
    //  printf("## EXPERIMENTAL opt-unif ONLY \n");
    //  if(metricFE){
    //    Mesh<MetricFieldFE        > *msh = (Mesh<MetricFieldFE        > *) msh_g;
    //    rebalanceMesh<MetricFieldFE,2>(*msh);
    //  }else{
    //    Mesh<MetricFieldAnalytical> *msh = (Mesh<MetricFieldAnalytical> *) msh_g;
    //    rebalanceMesh<MetricFieldAnalytical,2>(*msh);
    //  }
    //  writeOutputs();
    //  return 0;
    //}

    #ifdef TESTQUALITYALGO
    int niter = 5;
    #else
    int niter = 1;
    #endif

    #ifdef WRITEQUALFIELD
    #ifdef STEPDISTANCE
    constexpr QuaFun iquaf = QuaFun::StepDistance;
    #else
    constexpr QuaFun iquaf = QuaFun::SizeShape;
    #endif

    Mesh<MetricFieldAnalytical>& msh = static_cast<Mesh<MetricFieldAnalytical>&>(*msh_g);
    double scl = msh.param->metScale;
    const int tdim = msh.get_tdim();
    METRIS_ENFORCE_MSG(tdim == 3, "WRITEQUALFIELD not supported in 2D");
    #endif

    int iter = 0;
    while (iter < niter){

      foutputAdaptStats << std::setw(30) << iter;

      #ifdef WRITEQUALFIELD

      msh.cleanup();

      std::string meshQualFieldNameStart = "meshWithQualField_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".meshb";
      std::string qualFieldNameStart     = "meshWithQualField_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".solb";

      const int nenttStart = msh.nentt(tdim);
      dblAr1 rfldStart(nenttStart);
      for (int ientt = 0; ientt < nenttStart; ientt++){

        if (isdeadent(ientt,msh.ent2poi(tdim))) continue;

        double quaent;
        if (tdim == 2) quaent = metqua<MetricFieldAnalytical,2,2,iquaf>(msh,AsDeg::P1,AsDeg::P1,ientt,1.);
        else           quaent = metqua<MetricFieldAnalytical,3,3,iquaf>(msh,AsDeg::P1,AsDeg::P1,ientt,1.);
        rfldStart[ientt] = quaent;

      }

      writeMesh(meshQualFieldNameStart,msh);
      writeField(qualFieldNameStart,   msh,SolTyp::P0Elt,rfldStart);

      intAr1 lball(100);
      constexpr auto quafun_xi = get_quafun_xi<MetricFieldAnalytical,3,3,iquaf,double>();
      const int npoinStart = msh.npoin;
      dblAr1 quaCGP1(npoinStart);
      intAr2& ent2poi = msh.ent2poi(tdim);
      for (int ipoin = 0; ipoin < npoinStart; ipoin++){

        int ientt = getpoient(msh, ipoin, tdim);

        int iopen;
        int ierro = ball3(msh,ipoin,ientt,lball,&iopen,0);

        int nenttBall = 0;
        double quaSum = 0;
        for (int ienttBall : lball){

          int localIndex = -1;
          for (int ivertex = 0; ivertex < 4; ivertex++){
            if (ent2poi(ienttBall,ivertex) == ipoin){
              localIndex = ivertex;
              break;
            }
          }

          double bary[4] = {0,0,0,0};
          bary[localIndex] = 1.;

          double qpt = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ienttBall], bary, msh.met[ipoin]);

          quaSum += qpt;
          nenttBall++;
        }

        quaCGP1[ipoin] = quaSum/(double)nenttBall;
      }

      // for (int ipoin = 0; ipoin < npoinStart; ipoin++) quaCGP1[ipoin] = 1./quaCGP1[ipoin];

      std::string meshPointwiseQualFieldNameStart = "meshWithPointWiseQual_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".meshb";
      std::string pointwiseQualFieldNameStart     = "meshWithPointWiseQual_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".solb";

      writeMesh(meshPointwiseQualFieldNameStart,msh);
      writeField(pointwiseQualFieldNameStart, msh, SolTyp::CG, quaCGP1, 1);

      #endif

      adaptMesh2();

      if(param->usrTarDeg > 1 || msh_g->curdeg > 1){
        int ielev = degElevate();

        if(param->dbgfull) check_topo(*msh_g,0);

        if(param->curveType > 0 && !ielev){ // Not really smoothing, rather metric based curving
          curveMesh();
        }

        if(param->dbgfull) check_topo(*msh_g,0);
      }

      // Run the requested optimizer once, after any degree elevation and
      // high-order curving.
      optimMesh();

      iter++;
    }

    // A successful high-order run must finish with a degree-aware validity
    // check even when the expensive debug checks are otherwise disabled.
    if(msh_g->curdeg > 1) check_topo(*msh_g,0);

    writeOutputs();

    #ifdef OUTPUTTIMEANDUNITINFO
    printUnit = true;
    #endif

    if(DOPRINTS1()) statMesh();

  }catch(const MetrisExcept &e){
    printf("## Exception thrown, print mesh\n");
    std::string time_stamp = time2str();
    std::string fname = "exception-" + time_stamp;
    writeMesh(fname,*msh_g);
    if(metricFE){
      Mesh<MetricFieldFE> &msh = (Mesh<MetricFieldFE> &)(*msh_g);
      msh.met.writeMetricFile(fname);
    }else{
      Mesh<MetricFieldAnalytical> &msh = (Mesh<MetricFieldAnalytical> &)(*msh_g);
      msh.met.writeMetricFile(fname);
    }
    throw(e);
  }

#ifdef OUTPUTTIMEANDUNITINFO
  foutputTimeUnit << std::setw(30) << t1-t0
                  << std::endl;
#endif


  double actualCost = 0;

  double integratedCost = 0;
  if(metricFE){
    Mesh<MetricFieldFE> &msh = (Mesh<MetricFieldFE> &)(*msh_g);
    actualCost  = (double)msh.nentt(msh.get_tdim());
    if (msh.get_tdim() == 2){
      actualCost *= sqrt(3.)/4.;
    }else{
      actualCost *= sqrt(2.)/12.;
    }
    if (msh.get_tdim() == 2) integratedCost = getMetricCost<MetricFieldFE,2,2>(msh);
    else                     integratedCost = getMetricCost<MetricFieldFE,3,3>(msh);
  }else{
    Mesh<MetricFieldAnalytical> &msh = (Mesh<MetricFieldAnalytical> &)(*msh_g);
    actualCost  = (double)msh.nentt(msh.get_tdim());
    if (msh.get_tdim() == 2){
      actualCost *= sqrt(3.)/4.;
    }else{
      actualCost *= sqrt(2.)/12.;
    }
    if (msh.get_tdim() == 2) integratedCost = getMetricCost<MetricFieldAnalytical,2,2>(msh);
    else                     integratedCost = getMetricCost<MetricFieldAnalytical,3,3>(msh);
  }

  std::cout << std::setprecision(16) << "Integrated target cost = " << integratedCost << std::endl;
  std::cout << std::setprecision(16) << "Actual cost = " << actualCost << std::endl;

  if (param->MOESS_adapt_it >= 0){
    std::string meshName = "mesh_MOESS_final_a" + std::to_string(param->MOESS_adapt_it) + ".meshb";
    std::string metName  = "met_MOESS_final_a"  + std::to_string(param->MOESS_adapt_it) + ".solb";

    writeMesh(meshName,*msh_g);
    Mesh<MetricFieldFE> &msh = (Mesh<MetricFieldFE> &)(*msh_g);
    msh.met.writeMetricFile(metName);
  }

  foutputAdaptStats.close();

  t1 = get_cpu_time();
  MPRINTF("\n-- END Metris total runtime {:.2e}s\n",t1-t0);

}

template <class MFT>
void MetrisRunner::runMetrisProgressive(){

  GETVDEPTH(param);

  std::string fileName = "scaleEvol_sclInpt_" + std::to_string(param->metScale) + ".txt";
  std::fstream foutputScaleEvol(fileName, std::fstream::out);
  METRIS_ENFORCE_MSG(foutputScaleEvol.good(), "Error opening file: " + fileName);

  foutputScaleEvol  << "# scaleInput = " << param->metScale << std::endl;

  double t0, t1;
  t0 = get_cpu_time();

  const int niterAdapt = 1;

  // helper to set target metric scale
  auto setMetScale = [&](double sclmet) -> void {

    METRIS_ENFORCE(sclmet > 0);
    // The 1D pass must be repeated for each progressive target metric.
    objectiveLineAdapted = false;
    param->setMetricScale(sclmet);
    Mesh<MFT> &msh = (Mesh<MFT> &)(*msh_g);
    msh.met.normalize(sclmet);
    bak.met.normalize(sclmet);

    for (int ipoin = 0; ipoin < msh.npoin; ipoin++){
      if constexpr(std::is_same<MFT, MetricFieldAnalytical>::value){
        msh.met.getMetPhys(DifVar::None,msh.met.getSpace(),msh.coord[ipoin],msh.met[ipoin],NULL);
      }else{

        int tdim = msh.get_tdim();
        int ientt = getpoient(msh,ipoin,tdim);
        intAr2& ent2poi = msh.ent2poi(tdim);

        // if (tdim == 2){
        //   int enttdim = msh.poi2ent(ipoin,1);
        //   if (enttdim == 2) ientt = msh.poi2ent(ipoin,0);
        //   else if (enttdim == 1){
        //     int iedge = msh.poi2ent(ipoin,0);
        //     ientt = msh.edg2fac[iedge];
        //   }
        //   else METRIS_THROW_MSG("poi2ent points to a tdim=0 entt");
        // }
        // else{
        //   int enttdim = msh.poi2ent(ipoin,1);
        //   if (enttdim == 3) ientt = msh.poi2ent(ipoin,0);
        //   else if (enttdim == 2){
        //     int ifac = msh.poi2ent(ipoin,0);
        //     ientt = msh.fac2tet(ifac,0);
        //   }
        //   else if (enttdim == 1){
        //     int iedge = msh.poi2ent(ipoin,0);
        //     int ifac = msh.edg2fac[iedge];
        //     ientt = msh.fac2tet(ifac,0);
        //   }
        //   else METRIS_THROW_MSG("poi2ent points to a tdim=0 entt");
        // }

        int ibary = -1;

        for (int iver = 0; iver < tdim+1; iver++){
          if (ent2poi(ientt,iver) == ipoin){
            ibary = iver;
            break;
          }
        }
        METRIS_ENFORCE(ibary >=0);

        double bary[tdim+1];
        for (int ii = 0; ii < tdim+1; ii++) bary[ii] = 0.;
        bary[ibary] = 1.;

        msh.met.getMetBary(AsDeg::P1,
                          DifVar::None,
                          msh.met.getSpace(),
                          ent2poi[ientt],
                          tdim,
                          bary,
                          msh.met[ipoin],
                          NULL);
      }
    }
  };

  // helper to do a mesh adaptation cycle
  auto runAdaptCycle = [&](double sclmet = -1.) -> void {

    if (sclmet > 0) setMetScale(sclmet);

    for (int iterAdapt = 0; iterAdapt < niterAdapt; iterAdapt++){

      statMesh();

      if(param->dbgfull) check_topo(*msh_g,0);

      adaptMesh2();
      optimMesh();
    }
  };

  const double scaleInpt0 = param->metScale;

  double currentMeshCost;
  double currentMeshScale;
  double targetCost;
  double targetScale;

  Mesh<MFT> &msh = (Mesh<MFT> &)(*msh_g);
  msh.cleanup();
  currentMeshCost  = (double)msh.nentt(msh.get_tdim());
  if (msh.get_tdim() == 2){
    currentMeshCost *= sqrt(3.)/4.;
  }else{
    currentMeshCost *= sqrt(2.)/12.;
  }
  currentMeshScale = pow(currentMeshCost,-1./msh.get_tdim());

  if (msh.get_tdim() == 2) targetCost = getMetricCost<MFT,2,2>(msh);
  else                     targetCost = getMetricCost<MFT,3,3>(msh);

  targetScale = pow(targetCost,-1./msh.get_tdim());

  // the target mesh is coarser than initial mesh, no need to approach progressively, just call runMetris
  if (targetScale >= currentMeshScale){

    runAdaptCycle();

    writeOutputs();

    statMesh();

    t1 = get_cpu_time();

    MPRINTF("\n-- END Metris total runtime {:.2e}s\n",t1-t0);
    return;
  }

  targetScale *= 0.85; // assume we are overestimating the target scale

  foutputScaleEvol << "# Initial mesh scale = "   << currentMeshScale << std::endl;
  foutputScaleEvol << "# Approx. target scale = " << targetScale << std::endl;

  foutputScaleEvol  << std::setw(6) << "# Iter"
                    << std::setw(30) << "currenMeshScale"
                    << std::setw(30) << "targetScale"
                    << std::setw(30) << "scaleRun"
                    << std::endl;

  // the target mesh is finer than initial mesh, approach progressively
  double scaleRun = currentMeshScale/targetScale * scaleInpt0; // let first generation has similar target scale as current scale
  double lastRatio;
  int iiter = 0;
  bool exitLoop = false;
  while (true){

    std::cout << "Running adaptation for scale = " << scaleRun * targetScale / scaleInpt0 << std::endl;
    std::cout << "scaleRun = " << scaleRun << std::endl;

    foutputScaleEvol << std::setw(6)  << iiter
                     << std::setw(30) << currentMeshScale
                     << std::setw(30) << targetScale
                     << std::setw(30) << scaleRun
                     << std::endl;

    writeMesh("meshSTART_iter" + std::to_string(iiter) + ".meshb",msh);

    runAdaptCycle(scaleRun);
    msh.cleanup();
    setMetScale(scaleInpt0);

    writeMesh("meshEND_iter" + std::to_string(iiter) + ".meshb",msh);

    currentMeshCost  = (double)msh.nentt(msh.get_tdim());
    if (msh.get_tdim() == 2){
      currentMeshCost *= sqrt(3.)/4.;
    }else{
      currentMeshCost *= sqrt(2.)/12.;
    }
    currentMeshScale = pow(currentMeshCost,-1./msh.get_tdim());

    if (msh.get_tdim() == 2) targetCost = getMetricCost<MFT,2,2>(msh);
    else                     targetCost = getMetricCost<MFT,3,3>(msh);

    targetScale = pow(targetCost,-1./msh.get_tdim());

    lastRatio = currentMeshScale/targetScale;

    if (exitLoop){
      foutputScaleEvol << "# Final mesh scale = "   << currentMeshScale << std::endl;
      foutputScaleEvol << "# Approx. target scale = " << targetScale << std::endl;
      break;
    }

    scaleRun = 0.85*currentMeshScale / targetScale * scaleInpt0;

    std::cout << "After adaptation " << iiter << ":" << std::endl;
    std::cout << "currentMeshScale = " << currentMeshScale << std::endl;
    std::cout << "targetScale = " << targetScale << std::endl;
    std::cout << "scaleInpt0 = " << scaleInpt0 << std::endl;
    std::cout << "scaleRun = " << scaleRun << std::endl;
    std::cout << "Last currentScale/targetScale ratio = " << lastRatio << std::endl;

    if (scaleRun < scaleInpt0){

      scaleRun = scaleInpt0;
      exitLoop = true;
    }

    iiter++;
  }

  foutputScaleEvol.close();

  writeOutputs();

  statMesh();

  t1 = get_cpu_time();

  MPRINTF("\n-- END Metris total runtime {:.2e}s\n",t1-t0);

  return;
}

template void MetrisRunner::runMetrisProgressive<MetricFieldAnalytical>();
template void MetrisRunner::runMetrisProgressive<MetricFieldFE>();

} // end namespace
