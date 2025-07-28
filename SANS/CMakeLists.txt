add_library(libSANS_Metris tools/SANSException.cpp
                           tools/demangle.cpp
                           tools/BackTraceException.cpp
                           LinearAlgebra/DenseLinAlg/StaticSize/MatMul/MatrixS_MatMul_Native_Real_Instantiate.cpp
                           LinearAlgebra/DenseLinAlg/StaticSize/MatMul/MatrixS_MatMul_Native_Int_Instantiate.cpp
                           LinearAlgebra/DenseLinAlg/StaticSize/MatMul/MatrixS_MatMul_Native_Surreal_Instantiate.cpp
                           LinearAlgebra/DenseLinAlg/StaticSize/MatMul/MatrixS_MatMul_Native_MatrixS_Instantiate.cpp
                           )
target_include_directories(libSANS_Metris PUBLIC ${Boost_INCLUDE_DIRS}
                                                 $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                                                 $<INSTALL_INTERFACE:include>)

if(USE_GMP)
  target_link_libraries(libSANS_Metris PUBLIC ${GMP_LIBRARIES})
  target_include_directories(libSANS_Metris PUBLIC ${GMP_INCLUDE_DIRS})
endif()

list(APPEND METRIS_DEPS_LIBRARIES libSANS_Metris)
list(APPEND METRIS_INTERNAL_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/SANS)

set(METRIS_DEPS_LIBRARIES ${METRIS_DEPS_LIBRARIES} PARENT_SCOPE)
set(METRIS_INTERNAL_INCLUDE_DIRS ${METRIS_INTERNAL_INCLUDE_DIRS} PARENT_SCOPE)
