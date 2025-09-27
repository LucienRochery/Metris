add_library(libSANS_Metris STATIC 
              tools/SANSException.cpp
              tools/demangle.cpp
              tools/BackTraceException.cpp
              LinearAlgebra/DenseLinAlg/StaticSize/MatMul/MatrixS_MatMul_Native_Real_Instantiate.cpp
              LinearAlgebra/DenseLinAlg/StaticSize/MatMul/MatrixS_MatMul_Native_Int_Instantiate.cpp
              LinearAlgebra/DenseLinAlg/StaticSize/MatMul/MatrixS_MatMul_Native_Surreal_Instantiate.cpp
              LinearAlgebra/DenseLinAlg/StaticSize/MatMul/MatrixS_MatMul_Native_MatrixS_Instantiate.cpp
              )

target_include_directories(libSANS_Metris PUBLIC  
                            ${Boost_INCLUDE_DIRS}
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
                            $<INSTALL_INTERFACE:include/libs>)

if(USE_GMP)
  target_link_libraries(libSANS_Metris PUBLIC ${GMP_LIBRARIES})
  target_include_directories(libSANS_Metris PUBLIC ${GMP_INCLUDE_DIRS})
endif()

install(TARGETS libSANS_Metris
        EXPORT libMetrisTargets 
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

export(TARGETS libSANS_Metris FILE libMetrisTargets.cmake APPEND)


#list(APPEND METRIS_DEPS_LIBRARIES libSANS_Metris)
target_link_libraries(metris_deps INTERFACE libSANS_Metris)
