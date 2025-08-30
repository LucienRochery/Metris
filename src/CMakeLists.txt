
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
  setMetrisCXXFlags(gen_lageval PUBLIC)
  target_link_libraries(gen_lageval PUBLIC  Boost::headers 
                                            libSANS_Metris 
                                            )
  target_include_directories(gen_lageval PUBLIC ${PROJECT_SOURCE_DIR})
  add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/src/codegen_lagrange.hxx
                     DEPENDS gen_lageval
                     COMMAND make gen_lageval 
                     COMMAND ./gen_lageval)
endif()


# 1.2 - gen_ccoef: generate Jacobian determinant control coefficients
set(CCOEF_OUT_SRC "")
foreach(ideg RANGE 2 ${METRIS_MAX_DEG})
  if(REQ_CODEGEN)
    list(APPEND CCOEF_OUT_SRC ${PROJECT_BINARY_DIR}/src/codegen_ccoef3.0${ideg}.cxx)
    list(APPEND CCOEF_OUT_SRC ${PROJECT_BINARY_DIR}/src/codegen_ccoef2.0${ideg}.cxx)
    list(APPEND CCOEF_OUT_SRC ${PROJECT_BINARY_DIR}/src/codegen_ccoef3.0${ideg}_d_coord.cxx)
    list(APPEND CCOEF_OUT_SRC ${PROJECT_BINARY_DIR}/src/codegen_ccoef2.0${ideg}_d_coord.cxx)
    list(APPEND CCOEF_OUT_SRC ${PROJECT_BINARY_DIR}/src/codegen_ccoef2.0${ideg}_d_pt.cxx)
  else()
    list(APPEND CCOEF_OUT_SRC ${PROJECT_SOURCE_DIR}/src/codegen_ccoef3.0${ideg}.cxx)
    list(APPEND CCOEF_OUT_SRC ${PROJECT_SOURCE_DIR}/src/codegen_ccoef2.0${ideg}.cxx)
    list(APPEND CCOEF_OUT_SRC ${PROJECT_SOURCE_DIR}/src/codegen_ccoef3.0${ideg}_d_coord.cxx)
    list(APPEND CCOEF_OUT_SRC ${PROJECT_SOURCE_DIR}/src/codegen_ccoef2.0${ideg}_d_coord.cxx)
    list(APPEND CCOEF_OUT_SRC ${PROJECT_SOURCE_DIR}/src/codegen_ccoef2.0${ideg}_d_pt.cxx)
  endif()
endforeach()

if(REQ_CODEGEN)
  add_executable(gen_ccoef ${PROJECT_SOURCE_DIR}/tools/gen_ccoef.cxx 
                           ${PROJECT_SOURCE_DIR}/tools/gen_ccoef_d_pt.cxx 
                           ${PROJECT_SOURCE_DIR}/tools/gen_ccoef_d_coord.cxx 
                           ho_constants.hxx)
  setMetrisCXXFlags(gen_ccoef PUBLIC)
  target_link_libraries(gen_ccoef PUBLIC libSANS_Metris)
  target_include_directories(gen_ccoef PUBLIC ${METRIS_DEPS_INCLUDE_DIRS}
                                              ${PROJECT_SOURCE_DIR} )

  add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/src/codegen_ccoef.hxx 
                            ${PROJECT_BINARY_DIR}/src/codegen_ccoef_d.hxx 
                            ${CCOEF_OUT_SRC}
                    DEPENDS gen_ccoef
                    COMMAND make gen_ccoef
                    COMMAND ./gen_ccoef)
                    #COMMAND cp codegen_* ${PROJECT_SOURCE_DIR}/src/)

  # 1 - wrapping up: create dependency for later targets
  add_custom_target(codegen1 DEPENDS ${PROJECT_BINARY_DIR}/src/codegen_ccoef.hxx 
                                     ${PROJECT_BINARY_DIR}/src/codegen_ccoef_d.hxx 
                                     ${CCOEF_OUT_SRC}
                                     ${PROJECT_BINARY_DIR}/src/codegen_lagrange.hxx)
endif()


# 2 --- meshBASIC (low_eval.cxx)
# Files needed by stage 2 codegen which don't need stage 2 codegen themselves...
add_library(meshBASIC STATIC low_eval.cxx ho_constants.hxx
                             Arrays/aux_msharrays1D.cxx 
                             Arrays/aux_msharrays2D.cxx 
                             utils/aux_timer.cxx )
setMetrisCXXFlags(meshBASIC PUBLIC)
set_property(TARGET meshBASIC PROPERTY CXX_STANDARD 17)
target_include_directories(meshBASIC PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
                           $<INSTALL_INTERFACE:include/src>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}> # For directories SANS/ libs/ etc.
                           $<INSTALL_INTERFACE:include>
                           ${METRIS_EXTERNAL_INCLUDE_DIRS})
target_link_libraries(meshBASIC PUBLIC ${ABSL_LIBRARIES}
 ${METRIS_DEPS_LIBRARIES}
                                       ) 
if(REQ_CODEGEN)
  add_dependencies(meshBASIC codegen1  gen_lageval)

  # 3 --- CODEGEN TARGETS 2
  # gen_lag2bez: generate Lagrange-to-Bézier conversion routines
  add_executable(gen_lag2bez ${PROJECT_SOURCE_DIR}/tools/gen_lag2bezfuncs.cxx)
  setMetrisCXXFlags(gen_lag2bez PUBLIC)
  target_link_libraries(gen_lag2bez PUBLIC meshBASIC 
                                           libSANS_Metris
                                           ${ABSL_LIBRARIES})
  target_include_directories(gen_lag2bez PUBLIC 
                             ${PROJECT_SOURCE_DIR})

  add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/src/codegen_lag2bez.hxx 
                            ${PROJECT_BINARY_DIR}/src/codegen_lag2bez.cxx
                     DEPENDS gen_lag2bez
                     COMMAND make gen_lag2bez
                     COMMAND ./gen_lag2bez)
                     #COMMAND cp codegen_* ${PROJECT_SOURCE_DIR}/src/)

  add_custom_target(codegen2 DEPENDS codegen1
                                     ${PROJECT_BINARY_DIR}/src/codegen_lag2bez.hxx
                                     ${PROJECT_BINARY_DIR}/src/codegen_lag2bez.cxx)
endif()

set(CODEGEN_OUT_SRC ${CCOEF_OUT_SRC})
if(REQ_CODEGEN)
  list(APPEND CODEGEN_OUT_SRC ${PROJECT_BINARY_DIR}/src/codegen_lag2bez.cxx)
  list(APPEND METRIS_DEPS_INCLUDE_DIRS ${PROJECT_BINARY_DIR}/src)
else()
  list(APPEND CODEGEN_OUT_SRC ${PROJECT_SOURCE_DIR}/src/codegen_lag2bez.cxx)
endif()

# These are large files that have been vetted
# We get compile timeouts after dozens of minutes compiling 
# with clang -fsanitize=address
foreach(ff ${CODEGEN_OUT_SRC})
  set_source_files_properties(${ff} PROPERTIES COMPILE_FLAGS "-fno-sanitize=address")
endforeach()

if(METRIS_MAX_DEG)
  add_definitions(-DMETRIS_MAX_DEG=${METRIS_MAX_DEG})
else()
  message(FATAL_ERROR "Define METRIS_MAX_DEG")
endif()

#file (GLOB_RECURSE SRC_SOURCES CONFIGURE_DEPENDS "*.cxx")


set(SRC_SOURCES 
    ${PROJECT_SOURCE_DIR}/src/aux_hashtab.cxx
    ${PROJECT_SOURCE_DIR}/src/aux_topo.cxx
    ${PROJECT_SOURCE_DIR}/src/aux_histogram.cxx
    ${PROJECT_SOURCE_DIR}/src/io_libmeshb.cxx
    ${PROJECT_SOURCE_DIR}/src/low_eval.cxx
    ${PROJECT_SOURCE_DIR}/src/msh_lenedg.cxx
    ${PROJECT_SOURCE_DIR}/src/low_topo.cxx
    ${PROJECT_SOURCE_DIR}/src/main_adap.cxx
    ${PROJECT_SOURCE_DIR}/src/msh_anamet_context.cxx
    ${PROJECT_SOURCE_DIR}/src/msh_anamet_2D.cxx
    ${PROJECT_SOURCE_DIR}/src/msh_anamet_3D.cxx
    ${PROJECT_SOURCE_DIR}/src/msh_inineigh.cxx
    ${PROJECT_SOURCE_DIR}/src/msh_lag2bez.cxx
    ${PROJECT_SOURCE_DIR}/src/msh_checktopo.cxx
    ${PROJECT_SOURCE_DIR}/src/msh_degelev.cxx

    ${PROJECT_SOURCE_DIR}/src/utils/aux_misc.cxx
    ${PROJECT_SOURCE_DIR}/src/utils/aux_timer.cxx
    ${PROJECT_SOURCE_DIR}/src/utils/mprintf.cxx
    ${PROJECT_SOURCE_DIR}/src/utils/binom.cxx
    ${PROJECT_SOURCE_DIR}/src/utils/bernstein_prod.cxx
    ${PROJECT_SOURCE_DIR}/src/utils/EGADSprinterr.cxx
    
    #${PROJECT_SOURCE_DIR}/src/utils/tuple_hashtable.cxx
    ${PROJECT_SOURCE_DIR}/src/Arrays/aux_msharrays1D.cxx
    ${PROJECT_SOURCE_DIR}/src/Arrays/aux_msharrays2D.cxx

    ${PROJECT_SOURCE_DIR}/src/low_geo/ccoef.cxx
    ${PROJECT_SOURCE_DIR}/src/low_geo/measure.cxx
    ${PROJECT_SOURCE_DIR}/src/low_geo/height.cxx
    ${PROJECT_SOURCE_DIR}/src/low_geo/normal.cxx
    ${PROJECT_SOURCE_DIR}/src/low_geo/misc.cxx
    ${PROJECT_SOURCE_DIR}/src/low_geo/lenedg.cxx

    ${PROJECT_SOURCE_DIR}/src/Optimization/Newton.cxx
    ${PROJECT_SOURCE_DIR}/src/Optimization/DIRECT.cxx
    ${PROJECT_SOURCE_DIR}/src/Optimization/nlopt_luksan_mssubs.cxx
    ${PROJECT_SOURCE_DIR}/src/Optimization/nlopt_luksan_pssubs.cxx

    ${PROJECT_SOURCE_DIR}/src/Localization/low_localization.cxx
    ${PROJECT_SOURCE_DIR}/src/Localization/msh_localization.cxx
    ${PROJECT_SOURCE_DIR}/src/Localization/msh_interpFrontBack.cxx

    ${PROJECT_SOURCE_DIR}/src/Boundary/msh_inisurf.cxx
    ${PROJECT_SOURCE_DIR}/src/Boundary/low_projsurf.cxx
    #${PROJECT_SOURCE_DIR}/src/Boundary/msh_optsurf.cxx
    #${PROJECT_SOURCE_DIR}/src/Boundary/low_geosurf.cxx

    ${PROJECT_SOURCE_DIR}/src/BezierOffsets/msh_bez2gaps.cxx
    ${PROJECT_SOURCE_DIR}/src/BezierOffsets/low_gaps.cxx
    ${PROJECT_SOURCE_DIR}/src/BezierOffsets/msh_curve_offsets.cxx

    #${PROJECT_SOURCE_DIR}/src/quality/low_surfqua.cxx
    ${PROJECT_SOURCE_DIR}/src/quality/low_metqua.cxx
    ${PROJECT_SOURCE_DIR}/src/quality/quafun_distortion.cxx
    ${PROJECT_SOURCE_DIR}/src/quality/quafun_unit.cxx
    ${PROJECT_SOURCE_DIR}/src/quality/quafun_tradet.cxx
    ${PROJECT_SOURCE_DIR}/src/quality/low_metqua_d.cxx
    ${PROJECT_SOURCE_DIR}/src/quality/msh_metqua.cxx

    ${PROJECT_SOURCE_DIR}/src/linalg/det.cxx
    ${PROJECT_SOURCE_DIR}/src/linalg/dsyevq.cxx
    ${PROJECT_SOURCE_DIR}/src/linalg/dsytrd.cxx
    ${PROJECT_SOURCE_DIR}/src/linalg/eigen.cxx
    ${PROJECT_SOURCE_DIR}/src/linalg/explogmet.cxx
    ${PROJECT_SOURCE_DIR}/src/linalg/matprods.cxx
    ${PROJECT_SOURCE_DIR}/src/linalg/invmat.cxx
    ${PROJECT_SOURCE_DIR}/src/linalg/utils.cxx

    ${PROJECT_SOURCE_DIR}/src/cavity/check_cavity_topo.cxx
    ${PROJECT_SOURCE_DIR}/src/cavity/correct_cavity.cxx
    ${PROJECT_SOURCE_DIR}/src/cavity/msh_cavity.cxx
    ${PROJECT_SOURCE_DIR}/src/cavity/reconnect_lincav.cxx
    ${PROJECT_SOURCE_DIR}/src/cavity/reconnect_faccav.cxx
    ${PROJECT_SOURCE_DIR}/src/cavity/reconnect_tetcav.cxx
    ${PROJECT_SOURCE_DIR}/src/cavity/update_bpois_newp.cxx
    ${PROJECT_SOURCE_DIR}/src/cavity/update_cavity.cxx

    ${PROJECT_SOURCE_DIR}/src/Mesh/CADInfo.cxx
    ${PROJECT_SOURCE_DIR}/src/Mesh/MeshBase.cxx
    ${PROJECT_SOURCE_DIR}/src/Mesh/MeshBase_init.cxx
    ${PROJECT_SOURCE_DIR}/src/Mesh/MeshBase_topo.cxx
    ${PROJECT_SOURCE_DIR}/src/Mesh/MeshBase_Work.cxx
    ${PROJECT_SOURCE_DIR}/src/Mesh/MeshBack.cxx
    ${PROJECT_SOURCE_DIR}/src/Mesh/MeshMetric.cxx
    ${PROJECT_SOURCE_DIR}/src/Mesh/Mesh_interpMetBack.cxx
    ${PROJECT_SOURCE_DIR}/src/Mesh/Mesh_cleanup.cxx
    ${PROJECT_SOURCE_DIR}/src/Mesh/Mesh.cxx
    ${PROJECT_SOURCE_DIR}/src/Mesh/WorkArray.cxx

    ${PROJECT_SOURCE_DIR}/src/MetricField/MetricFieldFE.cxx
    ${PROJECT_SOURCE_DIR}/src/MetricField/MetricFieldAnalytical.cxx
    ${PROJECT_SOURCE_DIR}/src/MetricField/MetricField_getMetNodes.cxx
    ${PROJECT_SOURCE_DIR}/src/MetricField/msh_explogmet.cxx
    ${PROJECT_SOURCE_DIR}/src/MetricField/msh_checkmet.cxx
    ${PROJECT_SOURCE_DIR}/src/MetricField/msh_intrinsicmet.cxx
    ${PROJECT_SOURCE_DIR}/src/MetricField/msh_snapmetsurf.cxx

    ${PROJECT_SOURCE_DIR}/src/SolutionField/interpError.cxx
    ${PROJECT_SOURCE_DIR}/src/SolutionField/minInterpError.cxx
    ${PROJECT_SOURCE_DIR}/src/SolutionField/SolutionField.cxx
    ${PROJECT_SOURCE_DIR}/src/SolutionField/anasol2D.cxx
    ${PROJECT_SOURCE_DIR}/src/SolutionField/anasol3D.cxx

    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/MetrisParameters.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/MeshStat.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/MetrisRunner.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/MetrisRunner_iniMetris.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/MetrisRunner_curveMesh.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/MetrisRunner_runMetris.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/MetrisRunner_statMesh.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/msh_adaptMesh.cxx
    ${PROJECT_SOURCE_DIR}/src/MetrisRunner/msh_optimMesh.cxx

    ${PROJECT_SOURCE_DIR}/src/Adaptation/msh_collapse.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/low_collapse.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/msh_swap.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/low_swapface.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/low_swaptetra.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/low_increasecav.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/msh_lineadapt.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/msh_lineforce.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/low_delaunay.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/msh_reinsert_flat.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/low_cavqual.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/Insertion/msh_insert.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/Insertion/low_insert.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/Insertion/aux_insert.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/Insertion/EdgeSeed.cxx
    ${PROJECT_SOURCE_DIR}/src/Adaptation/Insertion/low_Steiner.cxx

    ${PROJECT_SOURCE_DIR}/src/smoothing/msh_smooball.cxx
    ${PROJECT_SOURCE_DIR}/src/smoothing/msh_smoolen.cxx
    #${PROJECT_SOURCE_DIR}/src/smoothing/msh_smoofull.cxx
    #${PROJECT_SOURCE_DIR}/src/smoothing/msh_smoofull_TAO.cxx
    #${PROJECT_SOURCE_DIR}/src/smoothing/msh_smoofull_custom.cxx
    #${PROJECT_SOURCE_DIR}/src/smoothing/msh_smoofull_costfun.cxx
    #${PROJECT_SOURCE_DIR}/src/smoothing/low_smoodirect.cxx
    ${PROJECT_SOURCE_DIR}/src/smoothing/low_smooballdiff.cxx
    ${PROJECT_SOURCE_DIR}/src/smoothing/low_smoolen.cxx

    ${PROJECT_SOURCE_DIR}/src/API/MetrisAPI.cxx
    
    ${PROJECT_SOURCE_DIR}/src/LPopt/msh_maxccoef.cxx
    ${PROJECT_SOURCE_DIR}/src/LPopt/msh_maxmetLP.cxx
    ${PROJECT_SOURCE_DIR}/src/LPopt/msh_bezgapLP.cxx
    ${PROJECT_SOURCE_DIR}/src/LPopt/LPsolver.cxx
    ${CODEGEN_OUT_SRC}
)
if(EXISTS ${PROJECT_SOURCE_DIR}/src/Adaptation/msh_lineadapt2.cxx)
  set(SRC_SOURCES 
      ${SRC_SOURCES}
      ${PROJECT_SOURCE_DIR}/src/Adaptation/msh_lineadapt2.cxx
      ${PROJECT_SOURCE_DIR}/src/Adaptation/msh_genPointsCurve.cxx
      ${PROJECT_SOURCE_DIR}/src/Adaptation/msh_insPointsCurve.cxx
  )
endif()
add_library(libMetris STATIC ${SRC_SOURCES})
setMetrisCXXFlags(libMetris PUBLIC)
set_property(TARGET libMetris PROPERTY CXX_STANDARD 17)

message("-- METRIS_DEPS_INCLUDE_DIRS = ${METRIS_DEPS_INCLUDE_DIRS}")
target_include_directories(libMetris PUBLIC  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
                                             $<INSTALL_INTERFACE:include/src>
                                             $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}> # For directories SANS/ libs/ etc.
                                             $<INSTALL_INTERFACE:include>
                                             ${METRIS_EXTERNAL_INCLUDE_DIRS} )

target_link_libraries(libMetris PUBLIC ${METRIS_DEPS_LIBRARIES})

if(REQ_CODEGEN)
  add_dependencies(libMetris codegen1 codegen2)
endif()

SET_TARGET_PROPERTIES(libMetris PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)

list(APPEND CMAKE_BUILD_RPATH ${CMAKE_BINARY_DIR}/src/ ${CMAKE_BINARY_DIR}/SANS/)
list(APPEND CMAKE_INSTALL_RPATH ${CMAKE_BINARY_DIR}/src/ ${CMAKE_BINARY_DIR}/SANS/)

# C Prepro target for debugging (optional)
add_library(metrisPrepro STATIC ${SRC_SOURCES})
set_property(TARGET metrisPrepro PROPERTY CXX_STANDARD 17)
set_target_properties(metrisPrepro PROPERTIES EXCLUDE_FROM_ALL TRUE)
target_compile_options(metrisPrepro PUBLIC ${METRIS_CXX_FLAGS} -E)

target_include_directories(metrisPrepro PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
                                               $<INSTALL_INTERFACE:include/src>
                                               $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}> # For directories SANS/ libs/ etc.
                                               $<INSTALL_INTERFACE:include>
                                               ${METRIS_EXTERNAL_INCLUDE_DIRS} )
if(REQ_CODEGEN)
  add_dependencies(libMetris codegen1 codegen2)
endif() 
   
