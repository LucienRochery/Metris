//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "../src/main_adap.hxx"
#include "../src/metris_options.hxx"
#ifdef USE_PETSC
  #include <petscsys.h>
#endif

/*
(tet2poi|tet2tet)\((iele[m2]),([a-Z0-9]*)\)
(tet2poi|tet2tet)\((iele[m2]),(lno(fa|ed)3\[i(fa|ed)2?\]\[[0-9a-Z]\])\)
$1[$2][$3]
*/

using namespace Metris;


int main(int argc, char** argv){



  //char **argv2 = (char**) malloc(256*sizeof(char*));
  //int argc2;
  //gen_argv(&argc2,argv2,"-ksp_monitor -start_in_debugger --with-strict-petscerrorcode");
  cargHandler arg2("-ksp_monitor -start_in_debugger --with-strict-petscerrorcode");

  printf("Call: ");
  for(int ii = 0; ii < argc; ii++){
    printf(" %s ",argv[ii]);
  }
  printf("\n");
//  gen_argv(&argc2,argv2,"");

  #ifdef USE_PETSC
    PetscFunctionBeginUser;
    //PetscCall(PetscInitialize(&arg2.c,&arg2.v,(char *)NULL, "Default help message"));
    PetscCall(PetscInitialize(&argc,&argv,(char *)NULL, "Default help message"));
  
    PetscMPIInt MPI_Rank = 0;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &MPI_Rank));
  #endif

  int icod; 
  
#ifdef USE_PETSC
  if(MPI_Rank == 0){
#endif
    //Mesh msh, bak;
    try{
      icod = main_metris(argc, argv);
      //icod = main_metris(argc, argv, msh, bak);
    }catch(const MetrisExcept &e){
      printf("\n\n## MAIN_METRIS THROWS EXCEPTION:\n");
      std::cout<<"## Type: "<<e.what()<<std::endl;
  
    #ifndef NO_BOOST_EXCEPT
      if(std::string const * ms=boost::get_error_info<excMessage>(e) )
        std::cout<<"## Message: "<<*ms; 
      if(boost::stacktrace::stacktrace const * tr=boost::get_error_info<excStackTrace>(e) )
        std::cerr << "## Call stack: \n" << *tr;
    #endif
    }
#ifdef USE_PETSC
  }
#endif
  #ifdef USE_PETSC
    PetscCall(PetscFinalize());
  #endif
  //for(int ii = 0; ii < argc2; ii++) free(argv2[ii]);
  //free(argv2);

  return icod;
}


