if(NOT TARGET libmeshb)
  add_library(libmeshb STATIC libmeshb7.c)
endif()
add_library(lplib3  STATIC lplib3.c)

add_library(libkdtree INTERFACE)
target_include_directories(libkdtree INTERFACE "${PROJECT_SOURCE_DIR}/libs/libkdtree")
file (GLOB_RECURSE ALGLIB_SRC CONFIGURE_DEPENDS "alglib-cpp/src/*.cpp")
add_library(alglib STATIC ${ALGLIB_SRC})
target_include_directories(alglib INTERFACE "${PROJECT_SOURCE_DIR}/libs/alglib-cpp/src")

#add_library(nlopt INTERFACE)
#target_include_directories(nlopt INTERFACE "${PROJECT_SOURCE_DIR}/libs/nlopt")
#ADD_SUBDIRECTORY("${PROJECT_SOURCE_DIR}/libs/nlopt")

#add_library(nlopt SHARED IMPORTED) # or STATIC instead of SHARED
#set_target_properties(nlopt PROPERTIES
#  IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/lib/libbar.so"
#  INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/include/libbar"
#)

#find_library(NAMES nlopt PATHS ${PROJECT_SOURCE_DIR}/libs/nlopt)
#link_directories(${PROJECT_SOURCE_DIR}/libs/nlopt)
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/nlopt)
# the nlopt_config.h file is found in here
set(NLOPT_INCLUDE ${PROJECT_BINARY_DIR}/libs/nlopt PARENT_SCOPE)
