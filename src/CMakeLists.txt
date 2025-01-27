
# -- 
# libMetris is the main target: all files compiled
# It requires 3 codegen targets of which one requires the other two... and also 
# partial compilation of src/* files. 
# This proceeds in three steps: 
#   - stage 1 codegen needs nothing and generates files needed for 
#     stage 2 (third codegen target) codegen
#   - meshBASIC compiles the needed files (which require stage 1 codegen) for:
#   - stage 2 codegen compiles the last codegen target 
# -- 



# 1 --- CODEGEN TARGETS
# 1.1 - gen_lageval: generate Lagrange evaluation routines
if(REQ_CODEGEN)
  add_executable(gen_lageval ${PROJECT_SOURCE_DIR}/tools/gen_lageval.cxx ho_constants.hxx)
  setMetrisFlags(gen_lageval PUBLIC)
  target_link_libraries(gen_lageval PUBLIC  ${GINAC_LIBRARIES} Boost::headers SANS ${BOOST_TRACELIBS})
  target_include_directories(gen_lageval PUBLIC 
  	                         ${PROJECT_SOURCE_DIR}
                             ${PETSC_INCLUDE}
                             ${GINAC_INCLUDE_DIRS}
                             ${CLN_INCLUDE_DIRS}
                             ${Boost_INCLUDE_DIRS}
  	                         "/usr/local/include" ".")
  add_custom_command(OUTPUT  ${PROJECT_SOURCE_DIR}/src/codegen_lagrange.hxx
  	                 COMMAND make gen_lageval 
  	                 COMMAND ./gen_lageval
  	                 COMMAND cp codegen_* ${PROJECT_SOURCE_DIR}/src/)
endif()


# 1.2 - gen_ccoef: generate Jacobian determinant control coefficients
set(SRC_CCOEF_OUTPUT "")
foreach(ideg RANGE 2 ${METRIS_MAX_DEG})
	list(APPEND SRC_CCOEF_OUTPUT ${PROJECT_SOURCE_DIR}/src/codegen_ccoef3.0${ideg}.cxx)
	list(APPEND SRC_CCOEF_OUTPUT ${PROJECT_SOURCE_DIR}/src/codegen_ccoef2.0${ideg}.cxx)
  list(APPEND SRC_CCOEF_OUTPUT ${PROJECT_SOURCE_DIR}/src/codegen_ccoef3.0${ideg}_d.cxx)
  list(APPEND SRC_CCOEF_OUTPUT ${PROJECT_SOURCE_DIR}/src/codegen_ccoef2.0${ideg}_d.cxx)
  list(APPEND SRC_CCOEF_OUTPUT ${PROJECT_SOURCE_DIR}/src/codegen_ccoef2.0${ideg}_d_pt.cxx)
endforeach()



if(REQ_CODEGEN)
  add_executable(gen_ccoef ${PROJECT_SOURCE_DIR}/tools/gen_ccoef.cxx ${PROJECT_SOURCE_DIR}/tools/gen_ccoef_d.cxx ho_constants.hxx)
  setMetrisFlags(gen_ccoef PUBLIC)
  target_link_libraries(gen_ccoef PUBLIC SANS ${BOOST_TRACELIBS})
  target_include_directories(gen_ccoef PUBLIC 
  	                         ${PROJECT_SOURCE_DIR} 
                             ${Boost_INCLUDE_DIRS}
  	                         "/usr/local/include" ".")
  add_custom_command(OUTPUT  ${PROJECT_SOURCE_DIR}/src/codegen_ccoef.hxx ${PROJECT_SOURCE_DIR}/src/codegen_ccoef_d.hxx ${SRC_CCOEF_OUTPUT}
  	                 COMMAND make gen_ccoef
  	                 COMMAND ./gen_ccoef
  	                 COMMAND cp codegen_* ${PROJECT_SOURCE_DIR}/src/)

  # 1 - wrapping up: create dependency for later targets
  add_custom_target(codegen1 DEPENDS ${PROJECT_SOURCE_DIR}/src/codegen_ccoef.hxx
                                     ${PROJECT_SOURCE_DIR}/src/codegen_ccoef_d.hxx
                                     ${PROJECT_SOURCE_DIR}/src/codegen_lagrange.hxx)
endif()


# 2 --- meshBASIC (low_eval.cxx)
# Files needed by stage 2 codegen which don't need stage 2 codegen themselves...
add_library(meshBASIC low_eval.cxx ho_constants.hxx
                      Arrays/aux_msharrays1D.cxx 
                      Arrays/aux_msharrays2D.cxx )
target_include_directories(meshBASIC PRIVATE "." ".." )
target_include_directories(meshBASIC PUBLIC 
													 "libs"
					                 "/usr/local/include"
					                 ${LAPACK_INCLUDE_DIR}
	                         ${ESP_ROOT}/include
	                         ${PETSC_INCLUDE}
                           ${ABSL_INCLUDE_DIRS}
                           ${Boost_INCLUDE_DIRS}
	                         ${SRC_INCLUDE_DIRS})
target_link_libraries(meshBASIC PUBLIC absl::hash 
                                       absl::flat_hash_map
                                       ${BOOST_TRACELIBS}
                                       ) 
if(REQ_CODEGEN)
  add_dependencies(meshBASIC codegen1)


  # 3 --- CODEGEN TARGETS 2
  # gen_lag2bez: generate Lagrange-to-Bézier conversion routines
  add_executable(gen_lag2bez ${PROJECT_SOURCE_DIR}/tools/gen_lag2bezfuncs.cxx)
  setMetrisFlags(gen_lag2bez PUBLIC)
  target_link_libraries(gen_lag2bez PUBLIC meshBASIC 
                                           SANS
                                           ${BOOST_TRACELIBS}
  														        	   absl::hash absl::flat_hash_map)
  target_include_directories(gen_lag2bez PUBLIC 
                             ${PROJECT_BINARY_DIR}
                             ${ABSL_INCLUDE_DIRS}
  					                 ${PROJECT_SOURCE_DIR}
                             ${Boost_INCLUDE_DIRS}
  					                 "/usr/local/include")
  target_include_directories(gen_lag2bez PRIVATE ${SRC_INCLUDE_DIRS})
  target_include_directories(gen_lag2bez PRIVATE ".")
  add_custom_command(OUTPUT  ${PROJECT_SOURCE_DIR}/src/codegen_lag2bez.hxx ${PROJECT_SOURCE_DIR}/src/codegen_lag2bez.cxx
  	                 COMMAND make gen_lag2bez
  	                 COMMAND ./gen_lag2bez
  	                 COMMAND cp codegen_* ${PROJECT_SOURCE_DIR}/src/)
  add_custom_target(codegen2 DEPENDS codegen1
                                     ${PROJECT_SOURCE_DIR}/src/codegen_lag2bez.hxx)
endif()


if(METRIS_MAX_DEG)
  add_definitions(-DMETRIS_MAX_DEG=${METRIS_MAX_DEG})
else()
	message(FATAL_ERROR "Define METRIS_MAX_DEG")
endif()

#file (GLOB_RECURSE SRC_SOURCES CONFIGURE_DEPENDS "*.cxx")

set(SRC_SOURCES 
    ${PROJECT_SOURCE_DIR}/src/Arrays/aux_msharrays1D.cxx
    ${PROJECT_SOURCE_DIR}/src/Arrays/aux_msharrays2D.cxx
	  ${PROJECT_SOURCE_DIR}/src/aux_hashtab.cxx
	  ${PROJECT_SOURCE_DIR}/src/aux_timer.cxx
	  ${PROJECT_SOURCE_DIR}/src/aux_topo.cxx
    ${PROJECT_SOURCE_DIR}/src/aux_utils.cxx
    ${PROJECT_SOURCE_DIR}/src/mprintf.cxx
    ${PROJECT_SOURCE_DIR}/src/aux_histogram.cxx
	  ${PROJECT_SOURCE_DIR}/src/aux_EGADSprinterr.cxx
	  ${PROJECT_SOURCE_DIR}/src/io_libmeshb.cxx
	  ${PROJECT_SOURCE_DIR}/src/low_ccoef.cxx
	  ${PROJECT_SOURCE_DIR}/src/low_eval.cxx
    ${PROJECT_SOURCE_DIR}/src/low_geo.cxx
    ${PROJECT_SOURCE_DIR}/src/low_lenedg.cxx
	  ${PROJECT_SOURCE_DIR}/src/msh_lenedg.cxx
	  ##${PROJECT_SOURCE_DIR}/src/low_subdivide.cxx
	  ${PROJECT_SOURCE_DIR}/src/low_topo.cxx
	  ${PROJECT_SOURCE_DIR}/src/main_adap.cxx
	  ${PROJECT_SOURCE_DIR}/src/msh_anamet_2D.cxx
	  ${PROJECT_SOURCE_DIR}/src/msh_anamet_3D.cxx
	  ${PROJECT_SOURCE_DIR}/src/msh_inineigh.cxx
	  ${PROJECT_SOURCE_DIR}/src/msh_intrinsicmet.cxx
	  ${PROJECT_SOURCE_DIR}/src/msh_lag2bez.cxx
	  ${PROJECT_SOURCE_DIR}/src/msh_checktopo.cxx
    ${PROJECT_SOURCE_DIR}/src/msh_degelev.cxx
	  ${PROJECT_SOURCE_DIR}/src/opt_generic.cxx
    ${PROJECT_SOURCE_DIR}/src/Localization/low_localization.cxx
    ${PROJECT_SOURCE_DIR}/src/Localization/msh_localization.cxx
    ${PROJECT_SOURCE_DIR}/src/Boundary/msh_inisurf.cxx
    ${PROJECT_SOURCE_DIR}/src/Boundary/low_projsurf.cxx
    #${PROJECT_SOURCE_DIR}/src/Boundary/msh_optsurf.cxx
    #${PROJECT_SOURCE_DIR}/src/Boundary/low_geosurf.cxx
    ${PROJECT_SOURCE_DIR}/src/BezierOffsets/msh_bez2gaps.cxx
    ${PROJECT_SOURCE_DIR}/src/BezierOffsets/low_gaps.cxx
    #${PROJECT_SOURCE_DIR}/src/BezierOffsets/low_gaps_3D.cxx
    ${PROJECT_SOURCE_DIR}/src/BezierOffsets/msh_curve_offsets.cxx
		#${PROJECT_SOURCE_DIR}/src/quality/low_surfqua.cxx
    ${PROJECT_SOURCE_DIR}/src/quality/low_metqua.cxx
    ${PROJECT_SOURCE_DIR}/src/quality/quafun_distortion.cxx
    ${PROJECT_SOURCE_DIR}/src/quality/quafun_unit.cxx
    ${PROJECT_SOURCE_DIR}/src/quality/quafun_tradet.cxx
		${PROJECT_SOURCE_DIR}/src/quality/low_metqua_d.cxx
    ${PROJECT_SOURCE_DIR}/src/quality/msh_metqua.cxx
		##${PROJECT_SOURCE_DIR}/src/quality/low_metqua_2D.cxx
		#${PROJECT_SOURCE_DIR}/src/quality/low_metqua_3D.cxx
		#${PROJECT_SOURCE_DIR}/src/quality/opt_metqua_2D.cxx
		#${PROJECT_SOURCE_DIR}/src/quality/opt_metqua_3D.cxx
		#${PROJECT_SOURCE_DIR}/src/quality/opt_metqua_costfunc.cxx
		#${PROJECT_SOURCE_DIR}/src/quality/opt_metqua_tests.cxx
		${PROJECT_SOURCE_DIR}/src/linalg/dsyevq.cxx
		${PROJECT_SOURCE_DIR}/src/linalg/dsytrd.cxx
		${PROJECT_SOURCE_DIR}/src/linalg/eigen.cxx
		${PROJECT_SOURCE_DIR}/src/linalg/explogmet.cxx
		${PROJECT_SOURCE_DIR}/src/linalg/matprods.cxx
		${PROJECT_SOURCE_DIR}/src/linalg/invmat.cxx
		${PROJECT_SOURCE_DIR}/src/linalg/utils.cxx
    ${PROJECT_SOURCE_DIR}/src/cavity/check_cavity_topo.cxx
		${PROJECT_SOURCE_DIR}/src/cavity/correct_cavity_fast.cxx
		${PROJECT_SOURCE_DIR}/src/cavity/msh_cavity.cxx
		${PROJECT_SOURCE_DIR}/src/cavity/reconnect_lincav.cxx
    ${PROJECT_SOURCE_DIR}/src/cavity/reconnect_faccav.cxx
		${PROJECT_SOURCE_DIR}/src/cavity/update_cavity.cxx
	  ${PROJECT_SOURCE_DIR}/src/codegen_lag2bez.cxx
	  #${PROJECT_SOURCE_DIR}/src/Mesh/Mesh.cxx
    ${PROJECT_SOURCE_DIR}/src/Mesh/CADInfo.cxx
	  ${PROJECT_SOURCE_DIR}/src/Mesh/MeshBase.cxx
    ${PROJECT_SOURCE_DIR}/src/Mesh/MeshBase_init.cxx
	  ${PROJECT_SOURCE_DIR}/src/Mesh/MeshBase_topo.cxx
	  ${PROJECT_SOURCE_DIR}/src/Mesh/MeshBack.cxx
	  ${PROJECT_SOURCE_DIR}/src/Mesh/MeshMetric.cxx
    ${PROJECT_SOURCE_DIR}/src/Mesh/Mesh.cxx
    ${PROJECT_SOURCE_DIR}/src/Mesh/Mesh_cleanup.cxx
	  ${PROJECT_SOURCE_DIR}/src/MetricField/MetricFieldFE.cxx
    ${PROJECT_SOURCE_DIR}/src/MetricField/MetricFieldAnalytical.cxx
	  ${PROJECT_SOURCE_DIR}/src/MetricField/MetricField_getMetNodes.cxx
    ${PROJECT_SOURCE_DIR}/src/MetricField/msh_explogmet.cxx
	  ${PROJECT_SOURCE_DIR}/src/MetricField/msh_checkmet.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/MetrisParameters.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/MetrisRunner.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/MetrisRunner_iniMetris.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/MetrisRunner_curveMesh.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/msh_adaptMesh.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/msh_optimMesh.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/statMesh.cxx
    ${PROJECT_SOURCE_DIR}/src/adapt/msh_collapse.cxx
	  ${PROJECT_SOURCE_DIR}/src/adapt/low_collapse.cxx
    ${PROJECT_SOURCE_DIR}/src/adapt/msh_swap2D.cxx
    ${PROJECT_SOURCE_DIR}/src/adapt/low_swap2D.cxx
    ${PROJECT_SOURCE_DIR}/src/adapt/msh_insert2D.cxx
    ${PROJECT_SOURCE_DIR}/src/adapt/low_insert2D.cxx
    ${PROJECT_SOURCE_DIR}/src/adapt/low_increasecav.cxx
    ${PROJECT_SOURCE_DIR}/src/adapt/msh_lineadapt.cxx
    ${PROJECT_SOURCE_DIR}/src/adapt/msh_lineforce.cxx
    ${PROJECT_SOURCE_DIR}/src/adapt/low_delaunay.cxx
    ${PROJECT_SOURCE_DIR}/src/adapt/msh_reinsert_flat.cxx
    ${PROJECT_SOURCE_DIR}/src/smoothing/msh_smooball.cxx
    #${PROJECT_SOURCE_DIR}/src/smoothing/msh_smoofull.cxx
    #${PROJECT_SOURCE_DIR}/src/smoothing/msh_smoofull_TAO.cxx
    #${PROJECT_SOURCE_DIR}/src/smoothing/msh_smoofull_custom.cxx
    #${PROJECT_SOURCE_DIR}/src/smoothing/msh_smoofull_costfun.cxx
    #${PROJECT_SOURCE_DIR}/src/smoothing/low_smoodirect.cxx
    ${PROJECT_SOURCE_DIR}/src/smoothing/low_smooballdiff.cxx
    #${PROJECT_SOURCE_DIR}/src/smoothing/msh_uniformization.cxx
    ${PROJECT_SOURCE_DIR}/src/API/MetrisAPI.cxx
    #${PROJECT_SOURCE_DIR}/src/API/MetrisAPI.cxx
    #${PROJECT_SOURCE_DIR}/src/API/MetrisAPI_vertices.cxx
    #${PROJECT_SOURCE_DIR}/src/API/MetrisAPI_elements.cxx
    #${PROJECT_SOURCE_DIR}/src/API/MetrisAPI_metfield.cxx
    #${PROJECT_SOURCE_DIR}/src/low_ccoef_d.cxx
    ${PROJECT_SOURCE_DIR}/src/LPopt/msh_maxccoef.cxx
    ${PROJECT_SOURCE_DIR}/src/LPopt/msh_maxmetLP.cxx
    ${PROJECT_SOURCE_DIR}/src/LPopt/msh_bezgapLP.cxx
    ${PROJECT_SOURCE_DIR}/src/LPopt/LPsolver.cxx
		${SRC_CCOEF_OUTPUT}
)
#message("In SRC got CLP_INCLUDE_DIRS ${CLP_INCLUDE_DIRS}")
add_library(libMetris ${SRC_SOURCES})
setMetrisFlags(libMetris PUBLIC)
set_property(TARGET libMetris PROPERTY CXX_STANDARD 17)
#get_target_property(DEBUG_METRISFLAGS libMetris COMPILE_OPTIONS)
#MESSAGE("libMetris compile options: ${DEBUG_METRISFLAGS}")

#set(NLOPT_INCLUDE ${PROJECT_BINARY_DIR}/libs/nlopt)
target_include_directories(libMetris PRIVATE "." "..")
target_include_directories(libMetris PUBLIC 
                           "../SANS/"
													 "../libs"
                           "../libs/libkdtree/"
                           "../libs/alglib/src/"
					                 "/usr/local/include"
					                 ${LAPACK_INCLUDE_DIR}
	                         ${ESP_ROOT}/include
	                         ${PETSC_INCLUDE}
	                         ${SRC_INCLUDE_DIRS}
                           ${CLP_INCLUDE_DIRS}
                           ${ABSL_INCLUDE_DIRS}
                           ${NLOPT_INCLUDE}
                           ${Boost_INCLUDE_DIRS}
	                         #${LPLIB_INCLUDE_DIR}
	                         )
if(CLP_FOUND)
  target_include_directories(libMetris PUBLIC ${CLP_INCLUDE_DIRS})
endif()
message("BOOST_TRACELIBS = ${BOOST_TRACELIBS}")
target_link_libraries(libMetris PUBLIC libmeshb
                                        SANS 
                                        #nlopt
                                        Boost::program_options
                                        #${Boost_PROGRAM_OPTIONS_LIBRARY}
														        	  dl
                                        absl::hash absl::flat_hash_map #absl::node_hash_map
                                        lapack 
                                        libegads libegadslite
                                        #LP.3
                                        lplib3
                                        ${PETSC_LIBRARIES}
                                        ${MPI_LIBRARIES}
                                        alglib
                                        ${BOOST_TRACELIBS}
                                        nlopt
                                        )
if(CLP_FOUND)
  target_link_libraries(libMetris PUBLIC ${CLP_LIBRARIES})
endif()

if(REQ_CODEGEN)
  add_dependencies(libMetris codegen2)
endif()
   