//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __COMMON_SETUP_MSH__
#define __COMMON_SETUP_MSH__

#include <string>
#include <src/types.hxx>
#include <src/io_libmeshb.hxx>
#include <src/ho_constants.hxx>
#include <src/aux_topo.hxx>
#include <src/utils/aux_misc.hxx>
#include <src/low_topo.hxx>

#include <src/low_geo.hxx>
#include <src/linalg/matprods.hxx>
#include <src/msh_degelev.hxx>
#include <src/aux_exceptions.hxx>
#include <src/MetrisRunner/MetrisRunner.hxx>

#include <src/main_adap.hxx>
#include <fcntl.h>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp> 
#include <boost/timer/progress_display.hpp>

namespace Metris{


#if 0
template<class MetricFieldType>
struct 
MeshTestSetup{
	MeshTestSetup() = delete;
//	MeshTestSetup(std::string meshFileName, int tardeg = 1){
//
//  	int usrMaxDeg = METRIS_MAX_DEG, usrMinDeg = tardeg; // Empty values
//
//  // usrMaxDeg is the very maximum the user is allowing for storage. It is hard bounded by the constant METRIS_MAX_DEG
//  // usrMinDeg is the minimum degree the user wants. 
//  if(usrMaxDeg > METRIS_MAX_DEG){
//    std::cout<<"!! MAX DEG "<<usrMaxDeg<<" SET ABOVE HARD LIMIT OF "<<METRIS_MAX_DEG<<std::endl;
//    usrMaxDeg = METRIS_MAX_DEG;
//  }
//  if(usrMinDeg > METRIS_MAX_DEG) std::cout<<"!! TARGET MIN DEG "<<usrMinDeg<<" SET ABOVE HARD LIMIT OF "<<METRIS_MAX_DEG<<std::endl;
//  if(usrMinDeg > usrMaxDeg){
//    std::cout<<"!! TARGET MIN DEG "<<usrMinDeg<<" ABOVE MAX ADMISSIBLE "<<usrMaxDeg<<std::endl;
//    usrMinDeg = usrMaxDeg;
//  }
//
//  nt npoin  = 0, nedge = 0, nface = 0, nelem = 0;
//  nt maxDeg = 0; // This will be determined by the reader 
//
//  nt ilagMsh;
//  td::string CADName = "";
//  if( iniMesh(meshFileName , CADName, usrMinDeg , usrMaxDeg , msh) ) {
//    std::cout<<"## FAILED TO READ MESH"<<std::endl;
//    exit(1);
//  }
//
//
//
//  int ierro;
//     // This is a compile-time for loop
//  hana::while_(hana::less_equal.than(hana::int_c<METRIS_MAX_DEG>), 1_c, [&](auto ideg_c){
//  constexpr int ideg = ideg_c;
//  hana::while_(hana::less_equal.than(hana::int_c<METRIS_MAX_DEG>), ideg_c+1_c, [&](auto tdeg_c){
//    constexpr int tdeg = tdeg_c;
//    if(ideg == msh.curdeg && tdeg == usrMinDeg){
//      if(msh.ilag == 1){
//        deg_elevate<ideg,tdeg,1>(msh,  &ierro);
//      }else{
//        deg_elevate<ideg,tdeg,0>(msh,  &ierro);
//      }
//    }
//    return tdeg_c+1_c;});
//  return ideg_c+1_c;});
//
//
//  }

  // This one allows us to control the initializer more finely with the full range of options
  MeshTestSetup(int argc, char** argv) : run(argc,argv){
    iniFromArgs(argc,argv);
  }

  //MeshTestSetup(std::string meshFileName, int tardeg = 1){
  //  char **argv = (char**) malloc(256*sizeof(char*));
  //  int argc;
  //  
  //  gen_argv(&argc,argv,"-i "+meshFileName+" -tardeg "+std::to_string(tardeg)+" -nosort");

  //  iniFromArgs(argc,argv);

  //  for(int i = 0; i < argc; i++) free(argv[i]);
  //  free(argv);
  //}


  void iniFromArgs(int argc, char** argv){
    // In Release mode, get rid of all initialization prints.
    run.degElevate();

    //METRIS_ENFORCE( (run.runnerMetricFE         && std::is_same<MetricFieldType, MetricFieldFE        >::value
    //            || run.runnerMetricAnalytical && std::is_same<MetricFieldType, MetricFieldAnalytical>::value));



    msh = (Mesh<MetricFieldType>*) run.msh_g;


    //try{
    //  int ret = main_metris(argc,argv,msh);
    //}catch(const MetrisExcept &e){
    //  printf("## EXCEPTION CAUGHT BY iniFromArgs\n");

    //  std::cout<<"## Type: "<<e.what()<<std::endl;
  
    //  if(std::string const * ms=boost::get_error_info<excMessage>(e) )
    //    std::cout<<"## Message: "<<*ms; 
    //  if(boost::stacktrace::stacktrace const * tr=boost::get_error_info<excStackTrace>(e) )
    //    std::cerr << "## Call stack: \n" << *tr;
    //}
    //#ifdef NDEBUG
    //  fflush(stdout);
    //  dup2(stdout_fd, STDOUT_FILENO);
    //  close(stdout_fd);
    //#endif
  }

	~MeshTestSetup(){} // Mesh and MshConComps already handle the deallocation

	Mesh<MetricFieldType> *msh;
  MetrisRunner run;
};
#endif



} // end namespace

#endif