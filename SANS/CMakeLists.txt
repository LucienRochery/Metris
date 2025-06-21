add_library(SANS tools/SANSException.cpp
                 tools/demangle.cpp
                 tools/BackTraceException.cpp
                 LinearAlgebra/DenseLinAlg/StaticSize/MatMul/MatrixS_MatMul_Native_Real_Instantiate.cpp
                 LinearAlgebra/DenseLinAlg/StaticSize/MatMul/MatrixS_MatMul_Native_Int_Instantiate.cpp
                 LinearAlgebra/DenseLinAlg/StaticSize/MatMul/MatrixS_MatMul_Native_Surreal_Instantiate.cpp
                 LinearAlgebra/DenseLinAlg/StaticSize/MatMul/MatrixS_MatMul_Native_MatrixS_Instantiate.cpp
                 )
target_include_directories(SANS PUBLIC ${Boost_INCLUDE_DIRS} 
                                       ${PROJECT_SOURCE_DIR})

if(USE_GMP)
  target_link_libraries(SANS PUBLIC ${GMP_LIBRARIES})
  target_include_directories(SANS PUBLIC ${GMP_INCLUDE_DIRS})
endif()

list(APPEND METRIS_DEPS_LIBRARIES SANS)
list(APPEND METRIS_DEPS_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/SANS)
