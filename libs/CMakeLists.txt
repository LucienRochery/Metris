list(APPEND METRIS_INTERNAL_INCLUDE_DIRS libs/)


# -- libMeshb and LPlib
if(NOT TARGET libmeshb)
  add_library(libmeshb STATIC libmeshb7/libmeshb7.c)
  setMetrisCFlags(libmeshb PUBLIC)
endif()
list(APPEND METRIS_DEPS_LIBRARIES libmeshb)

add_library(lplib3 STATIC lplib3/lplib3.c)
setMetrisCFlags(lplib3 PUBLIC)
target_link_libraries(lplib3 PUBLIC pthread)
list(APPEND METRIS_DEPS_LIBRARIES lplib3)



# -- libkdtree
add_library(libkdtree INTERFACE)
target_include_directories(libkdtree INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libs/libkdtree>
                                               $<INSTALL_INTERFACE:include/libs/libkdtree>)
list(APPEND METRIS_DEPS_LIBRARIES libkdtree)

# -- Alglib
file (GLOB_RECURSE ALGLIB_SRC CONFIGURE_DEPENDS "alglib-cpp/src/*.cpp")
add_library(alglib STATIC ${ALGLIB_SRC})
target_include_directories(alglib INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libs/alglib-cpp/src>
                                            $<INSTALL_INTERFACE:include/libs/alglib-cpp/src>)
list(APPEND METRIS_DEPS_LIBRARIES alglib)




set(METRIS_INTERNAL_INCLUDE_DIRS ${METRIS_INTERNAL_INCLUDE_DIRS} PARENT_SCOPE)
set(METRIS_DEPS_LIBRARIES        ${METRIS_DEPS_LIBRARIES} PARENT_SCOPE)