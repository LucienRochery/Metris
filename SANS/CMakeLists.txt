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
list(APPEND METRIS_DEPS_LIBRARIES SANS)
list(APPEND METRIS_DEPS_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/SANS)
