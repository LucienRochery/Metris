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

#ifdef WRITEQUALFIELD
#include "../quality/quafun.hxx"
#include "../quality/low_metqua.hxx"
#include "../low_topo.hxx"
#endif

namespace Metris{

void MetrisRunner::runMetris(){

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
    int niter = 10;
    #else
    int niter = 1;
    #endif

    #ifdef WRITEQUALFIELD
    Mesh<MetricFieldAnalytical>& msh = static_cast<Mesh<MetricFieldAnalytical>&>(*msh_g);
    double scl = msh.param->metScale;
    const int tdim = msh.get_tdim();
    METRIS_ENFORCE_MSG(tdim == 3, "WRITEQUALFIELD not supported in 2D");
    #endif

    int iter = 0;
    while (iter < niter){

      #ifdef WRITEQUALFIELD

      msh.cleanup();

      #ifdef ONEPOINTQUAL
      std::string meshQualFieldNameStart = "meshWithQualField1Point_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".meshb";
      std::string qualFieldNameStart     = "meshWithQualField1Point_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".solb";
      #elif defined(TDIM1POINTSQUAL)
      std::string meshQualFieldNameStart = "meshWithQualFieldTDIM1Points_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".meshb";
      std::string qualFieldNameStart     = "meshWithQualFieldTDIM1Points_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".solb";
      #elif defined(KEAST4QUAL)
      std::string meshQualFieldNameStart = "meshWithQualFieldKEAST4_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".meshb";
      std::string qualFieldNameStart     = "meshWithQualFieldKEAST4_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".solb";
      #endif

      const int nenttStart = msh.nentt(tdim);
      dblAr1 rfldStart(nenttStart);
      for (int ientt = 0; ientt < nenttStart; ientt++){

        if (isdeadent(ientt,msh.ent2poi(tdim))) continue;

        double quaent;
        if (tdim == 2) quaent = metqua<MetricFieldAnalytical,2,2,QuaFun::SizeShape>(msh,AsDeg::P1,AsDeg::P1,ientt,1.);
        else           quaent = metqua<MetricFieldAnalytical,3,3,QuaFun::SizeShape>(msh,AsDeg::P1,AsDeg::P1,ientt,1.);
        rfldStart[ientt] = quaent;

      }

      writeMesh(meshQualFieldNameStart,msh);
      writeField(qualFieldNameStart,   msh,SolTyp::P0Elt,rfldStart);

      intAr1 lball(100);
      constexpr auto quafun_xi = get_quafun_xi<MetricFieldAnalytical,3,3,QuaFun::SizeShape,double>();
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

      #if defined(ONEPOINTQUAL)
      std::string meshPointwiseQualFieldNameStart = "meshWithPointWiseQual1Point_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".meshb";
      std::string pointwiseQualFieldNameStart     = "meshWithPointWiseQual1Point_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".solb";
      #elif defined(TDIM1POINTSQUAL)
      std::string meshPointwiseQualFieldNameStart = "meshWithPointWiseQualTDIM1Points_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".meshb";
      std::string pointwiseQualFieldNameStart     = "meshWithPointWiseQualTDIM1Points_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".solb";
      #elif defined(KEAST4QUAL)
      std::string meshPointwiseQualFieldNameStart = "meshWithPointWiseQualKEAST4_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".meshb";
      std::string pointwiseQualFieldNameStart     = "meshWithPointWiseQualKEAST4_scl_" + std::to_string(scl) + "_StartIter" + std::to_string(iter) + ".solb";
      #endif

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

        if(param->smoo_type == 0){
          optimMesh();
        }

        if(param->dbgfull) check_topo(*msh_g,0);
      }

      optimMesh();

      iter++;
    }

    //if(param->anaSol && param->smoo_type == 1){
    //  optimMesh();
    //}
    writeOutputs();

    printUnit = true;
    if(DOPRINTS1()) statMesh();

    //#ifdef METRIS_USE_PETSC
    //  PetscCall(PetscFinalize());
    //#endif
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

  t1 = get_cpu_time();

#ifdef OUTPUTTIMEANDUNITINFO
  foutputTimeUnit << std::setw(30) << t1-t0
                  << std::endl;
#endif

  MPRINTF("\n-- END Metris total runtime {:.2e}s\n",t1-t0);

}



} // end namespace






