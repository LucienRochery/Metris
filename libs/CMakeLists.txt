if(NOT TARGET libmeshb)
  add_library(libmeshb STATIC libmeshb7.c)
  setMetrisFlags(libmeshb PUBLIC)
endif()
list(APPEND METRIS_DEPS_LIBRARIES libmeshb)

add_library(lplib3 STATIC lplib3.c)
list(APPEND METRIS_DEPS_LIBRARIES lplib3)
list(APPEND METRIS_DEPS_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/libs/)



add_library(libkdtree INTERFACE)
target_include_directories(libkdtree INTERFACE "${PROJECT_SOURCE_DIR}/libs/libkdtree")
set(METRIS_LIBKDTREE_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/src/libs/libkdtree/)
message("-- METRIS_LIBKDTREE_INCLUDE_DIRS = ${METRIS_LIBKDTREE_INCLUDE_DIRS}")
list(APPEND METRIS_DEPS_INCLUDE_DIRS ${METRIS_LIBKDTREE_INCLUDE_DIRS})

file (GLOB_RECURSE ALGLIB_SRC CONFIGURE_DEPENDS "alglib-cpp/src/*.cpp")
add_library(alglib STATIC ${ALGLIB_SRC})
target_include_directories(alglib INTERFACE "${PROJECT_SOURCE_DIR}/libs/alglib-cpp/src")
list(APPEND METRIS_DEPS_LIBRARIES alglib)
set(ALGLIB_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/libs/alglib-cpp/src)
message("-- ALGLIB_INCLUDE_DIRS = ${ALGLIB_INCLUDE_DIRS}")
list(APPEND METRIS_DEPS_INCLUDE_DIRS ${ALGLIB_INCLUDE_DIRS})




#add_subdirectory(${PROJECT_SOURCE_DIR}/libs/nlopt)
## the nlopt_config.h file is found in here
#set(NLOPT_INCLUDE_DIRS ${PROJECT_BINARY_DIR}/libs/nlopt)
#message("-- NLOPT_INCLUDE_DIRS = ${NLOPT_INCLUDE_DIRS}")
#list(APPEND METRIS_DEPS_INCLUDE_DIRS ${NLOPT_INCLUDE_DIRS})
#list(APPEND METRIS_DEPS_LIBRARIES nlopt_Metris)


set(METRIS_DEPS_INCLUDE_DIRS ${METRIS_DEPS_INCLUDE_DIRS} PARENT_SCOPE)
set(METRIS_DEPS_LIBRARIES    ${METRIS_DEPS_LIBRARIES} PARENT_SCOPE)