//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "common_includes.hxx"
#include "metris_options.hxx"

#include "MetrisRunner/MetrisRunner.hxx"
#include "msh_checktopo.hxx"

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
  MetrisParameters param = run.param;
  MetrisOptions opt = run.opt;


  int iverb = param.iverb;

  if(iverb >= 1) run.statMesh();

  if(param.dbgfull) check_topo(*run.msh_g);

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


  int ielev = run.degElevate();


  if(param.dbgfull) check_topo(*run.msh_g);

  
  if(param.curveType > 0 && !ielev){ // Not really smoothing, rather metric based curving
    run.curveMesh();
  }
  
  run.optimMesh();

  if(param.dbgfull) check_topo(*run.msh_g);

  run.writeOutputs();


  if(iverb >= 1) run.statMesh();

  //#ifdef USE_PETSC
  //  PetscCall(PetscFinalize());
  //#endif

	return 0;
}

} // End namespace

