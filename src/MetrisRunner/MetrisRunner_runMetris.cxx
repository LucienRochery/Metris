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

    //if(param->anaSol && param->smoo_type == 1){
    //  optimMesh();
    //}
    writeOutputs();


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
  MPRINTF("\n-- END Metris total runtime {:.2e}s\n",t1-t0);

}



} // end namespace






