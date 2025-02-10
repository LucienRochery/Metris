
# arg1: target name 
# arg2: scope (PUBLIC INTERFACE etc)
function(setMetrisFlags arg1 arg2)
  target_compile_options(${arg1} ${arg2} $<$<CONFIG:MEMCHECK>:${METRIS_CXX_FLAGS_MEMCHECK}>)
  target_compile_options(${arg1} ${arg2} $<$<CONFIG:DEBUG>:${METRIS_CXX_FLAGS_DEBUG}>)
  target_compile_options(${arg1} ${arg2} $<$<CONFIG:RELEASE>:${METRIS_CXX_FLAGS_RELEASE}>)
  #target_compile_options(${arg1} ${arg2} $<$<CONFIG:DEBUG>:-fsanitize=address>)
  #target_link_options(${arg1} ${arg2} $<$<CONFIG:DEBUG>:-fsanitize=address>)
  target_link_options(${arg1} ${arg2} $<$<CONFIG:MEMCHECK>:${METRIS_CXX_FLAGS_MEMCHECK}>)
  target_link_options(${arg1} ${arg2} $<$<CONFIG:DEBUG>:${METRIS_CXX_FLAGS_DEBUG}>)
  target_link_options(${arg1} ${arg2} $<$<CONFIG:RELEASE>:${METRIS_CXX_FLAGS_RELEASE}>)
endfunction()


set(CBT_lower ${CMAKE_BUILD_TYPE})
string( TOLOWER "${CBT_lower}" CBT_lower )

if(CBT_lower MATCHES "memcheck")
  remove_definitions("-DNDEBUG")
endif()

message("Metris using build type = ${CMAKE_BUILD_TYPE}")

set(METRIS_WARNING_FLAGS -Wno-gnu-zero-variadic-macro-arguments  -Wno-logical-op-parentheses
    -Wno-gcc-compat -Wno-variadic-macros)  
set(METRIS_CXX_FLAGS ${METRIS_WARNING_FLAGS})
if(USE_TRACELIBS)
  set(METRIS_CXX_FLAGS ${METRIS_CXX_FLAGS} -DBOOST_STACKTRACE_USE_ADDR2LINE)
endif()

#Somehow, a straight comparison with EQUAL icc or EQUAL icc doesn't register here. Perhaps there's a space in there. 
if(SHORT_COMPILER_NAME STREQUAL icc OR SHORT_COMPILER_NAME STREQUAL icx OR CMAKE_CXX_COMPILER_ID STREQUAL IntelLLVM)
  message("Using Intel compiler ${CMAKE_C_COMPILER} ${SHORT_COMPILER_NAME}")
  # NDEBUG should be set in release mode by default, but somehow this doesn't work with icc
  #  -guide=4
  #set(METRIS_CXX_FLAGS_RELEASE  -DNDEBUG -fno-alias -funroll-all-loops -fno-fnalias -fast -fno-protect-parens -Ofast -flto -diag-disable=10441 -qopt-subscript-in-range)
  set(METRIS_CXX_FLAGS_RELEASE -DNDEBUG -O3 -fPIC )
  set(METRIS_CXX_FLAGS_DEBUG   -g -O3 -diag-disable=10441 -fPIC)
  set(METRIS_C_FLAGS_RELEASE  -O3)
  #set(METRIS_C_FLAGS_RELEASE   -fno-alias -fno-fnalias -fast -fno-protect-parens -Ofast -flto -diag-disable=10441 -qopt-subscript-in-range)
  set(METRIS_C_FLAGS_DEBUG ${METRIS_CXX_FLAGS_DEBUG})
  message(Using Intel compiler ${CMAKE_C_COMPILER} ${SHORT_COMPILER_NAME})
elseif(CMAKE_CXX_COMPILER_ID STREQUAL GNU)
  message("Using GNU compiler ${CMAKE_C_COMPILER} ${SHORT_COMPILER_NAME}")
  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fconstexpr-ops-limit=10000000 -fPIC)
  set(METRIS_CXX_FLAGS_RELEASE   -DNDEBUG -O3 )
  #set(METRIS_CXX_FLAGS_DEBUG   -Og -ggdb3 -Wall -Wextra -pedantic  -march=native -no-pie -fno-pie  -rdynamic) # -S -fverbose-asm
  set(METRIS_CXX_FLAGS_DEBUG    -Og -g -Wall -march=native ) #  -rdynamic # -S -fverbose-asm -ggdb3
  set(METRIS_CXX_FLAGS_MEMCHECK -Os -fsanitize=address -fno-omit-frame-pointer)

  set(METRIS_C_FLAGS_RELEASE  ${METRIS_CXX_FLAGS_RELEASE})
  set(METRIS_C_FLAGS_DEBUG ${METRIS_CXX_FLAGS_DEBUG})
  set(METRIS_C_FLAGS_MEMCHECK ${METRIS_CXX_FLAGS_MEMCHECK})
elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  message("Using Clang ${CMAKE_C_COMPILER} ${SHORT_COMPILER_NAME}")
  set(METRIS_CXX_FLAGS_RELEASE   -DNDEBUG  -O3 -fPIC)
  set(METRIS_CXX_FLAGS_DEBUG      -Og -g  -Wall -Wextra -pedantic  -march=native  -fno-pie  -fPIC) # -S -fverbose-asm -rdynamic -ggdb3
  #set(METRIS_CXX_FLAGS_DEBUG  -fsanitize=address  -fconstexpr-steps=10000000 -O0 -g3  -march=native -fno-pie ) # -S -fverbose-asm
  set(METRIS_CXX_FLAGS_MEMCHECK    -Wall -fsanitize=address -Og -g3 -fPIC) # -S -fverbose-asm
  #-fno-omit-frame-pointer
  set(METRIS_C_FLAGS_RELEASE  ${METRIS_CXX_FLAGS_RELEASE})
  set(METRIS_C_FLAGS_DEBUG ${METRIS_CXX_FLAGS_DEBUG})
  set(METRIS_C_FLAGS_MEMCHECK ${METRIS_CXX_FLAGS_MEMCHECK})
else()
  message(FATAL_ERROR "Unknown compiler ID = ${CMAKE_CXX_COMPILER_ID}, SHORT_COMPILER_NAME = ${SHORT_COMPILER_NAME}")
endif()

set(METRIS_CXX_FLAGS_RELEASE  ${METRIS_CXX_FLAGS_RELEASE} ${METRIS_CXX_FLAGS} )
set(METRIS_CXX_FLAGS_DEBUG    ${METRIS_CXX_FLAGS_DEBUG} ${METRIS_CXX_FLAGS} )
set(METRIS_CXX_FLAGS_MEMCHECK ${METRIS_CXX_FLAGS_MEMCHECK} ${METRIS_CXX_FLAGS} )


if(PREPRO STREQUAL "True")
  message(WARNING "preprocessor mode")
  set(METRIS_CXX_FLAGS_RELEASE  ${METRIS_CXX_FLAGS_RELEASE} -E )
  set(METRIS_CXX_FLAGS_DEBUG ${METRIS_CXX_FLAGS_DEBUG} -E )
  set(METRIS_CXX_FLAGS_DEBUG ${METRIS_CXX_FLAGS_MEMCHECK} -E )
endif()