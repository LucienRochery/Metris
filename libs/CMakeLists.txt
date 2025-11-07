target_include_directories(metris_deps INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                            $<INSTALL_INTERFACE:include/libs>)

# -- libMeshb and LPlib
if(NOT TARGET libmeshb)
  add_library(libmeshb STATIC libmeshb7/libmeshb7.c)
  setMetrisCFlags(libmeshb PUBLIC)
endif()
if(METRIS_INSTALL)
  install(TARGETS libmeshb
          EXPORT libMetrisTargets
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(metris_deps INTERFACE libmeshb)

if(NOT TARGET lplib3)
  add_library(lplib3 STATIC lplib3/lplib3.c)
  target_link_libraries(lplib3 PUBLIC pthread)
  if(METRIS_INSTALL)
    install(TARGETS lplib3
            EXPORT libMetrisTargets
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()
endif()
target_link_libraries(metris_deps INTERFACE lplib3)




## -- libkdtree
#add_library(libkdtree INTERFACE)
#target_include_directories(libkdtree INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libs/libkdtree>
#                                               $<INSTALL_INTERFACE:include/libs/libkdtree>)
#list(APPEND METRIS_DEPS_LIBRARIES libkdtree)

# -- Alglib
file(GLOB_RECURSE ALGLIB_SRC CONFIGURE_DEPENDS "alglib-cpp/src/*.cpp")
add_library(alglib STATIC ${ALGLIB_SRC})
#target_include_directories(alglib INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libs/alglib-cpp/src>
#                                            $<INSTALL_INTERFACE:include/libs/alglib-cpp/src>)
target_link_libraries(metris_deps INTERFACE alglib)
if(METRIS_INSTALL)
  install(TARGETS alglib
          EXPORT libMetrisTargets
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

add_subdirectory(SANS)

