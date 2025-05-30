//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "common_includes.hxx"
#include "metris_options.hxx"
#include "utils/mprintf.hxx"

#include "MetrisRunner/MetrisRunner.hxx"
#include "msh_checktopo.hxx"
#include "io_libmeshb.hxx"

#ifdef USE_PETSC
  #include <petscsys.h>
#endif


namespace Metris{

// First API steps, mostly for use with Boost::test
int main_metris(int argc, char** argv){ 

  //#ifdef USE_PETSC
  //  static char help_PETSc[] = "PETSc Metris instance.\n\n";
  //  char help[] = "ok\0";
  //  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
  //#endif

  MetrisRunner run(argc,argv);
  MetrisParameters &param = run.param_;
  MetrisOptions opt = run.opt;


  try{

    if(param.iverb >= 1) run.statMesh();

    if(param.dbgfull) check_topo(*run.msh_g,0);

    //if(run.param.opt_unif){
    //  printf("## EXPERIMENTAL opt-unif ONLY \n");
    //  if(run.metricFE){
    //    Mesh<MetricFieldFE        > *msh = (Mesh<MetricFieldFE        > *) run.msh_g;
    //    rebalanceMesh<MetricFieldFE,2>(*msh);
    //  }else{
    //    Mesh<MetricFieldAnalytical> *msh = (Mesh<MetricFieldAnalytical> *) run.msh_g;
    //    rebalanceMesh<MetricFieldAnalytical,2>(*msh);
    //  }
    //  run.writeOutputs();
    //  return 0;
    //}


    run.adaptMesh();

    run.optimMesh();


    if(param.usrTarDeg > 1 || run.msh_g->curdeg > 1){
      int ielev = run.degElevate();

      if(param.dbgfull) check_topo(*run.msh_g,0);

      if(param.curveType > 0 && !ielev){ // Not really smoothing, rather metric based curving
        run.curveMesh();
      }

      if(param.smoo_type == 0){
        run.optimMesh();
      }
      
      if(param.dbgfull) check_topo(*run.msh_g,0);
    }

    if(param.anaSol && param.smoo_type == 1){
      run.optimMesh();
    }
    run.writeOutputs();


    if(param.iverb >= 1) run.statMesh();

    //#ifdef USE_PETSC
    //  PetscCall(PetscFinalize());
    //#endif
  }catch(const MetrisExcept &e){
    printf("## Exception thrown, print mesh\n");
    writeMesh("exception",*run.msh_g);
    throw(e);
  }

  return 0;
}

} // End namespace

