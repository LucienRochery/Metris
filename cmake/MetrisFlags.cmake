
# arg1: target name 
# arg2: scope (PUBLIC INTERFACE etc)
function(setMetrisFlags arg1 arg2)
  target_compile_options(${arg1} ${arg2} $<$<CONFIG:DEBUG>:${METRIS_CXX_FLAGS_DEBUG}>)
  target_compile_options(${arg1} ${arg2} $<$<CONFIG:RELEASE>:${METRIS_CXX_FLAGS_RELEASE}>)
  #target_compile_options(${arg1} ${arg2} $<$<CONFIG:DEBUG>:-fsanitize=address>)
  #target_link_options(${arg1} ${arg2} $<$<CONFIG:DEBUG>:-fsanitize=address>)
  target_link_options(${arg1} ${arg2} $<$<CONFIG:DEBUG>:${METRIS_CXX_FLAGS_DEBUG}>)
  target_link_options(${arg1} ${arg2} $<$<CONFIG:RELEASE>:${METRIS_CXX_FLAGS_RELEASE}>)
endfunction()



#Somehow, a straight comparison with EQUAL icc or EQUAL icc doesn't register here. Perhaps there's a space in there. 
if(SHORT_COMPILER_NAME STREQUAL icc OR SHORT_COMPILER_NAME STREQUAL icx)
  message("Using Intel compiler ${CMAKE_C_COMPILER} ${SHORT_COMPILER_NAME}")
  # NDEBUG should be set in release mode by default, but somehow this doesn't work with icc
  #  -guide=4
  #set(METRIS_CXX_FLAGS_RELEASE  -DNDEBUG -fno-alias -funroll-all-loops -fno-fnalias -fast -fno-protect-parens -Ofast -flto -diag-disable=10441 -qopt-subscript-in-range)
  set(METRIS_CXX_FLAGS_RELEASE  -DNDEBUG -O3 )
  set(METRIS_CXX_FLAGS_DEBUG    -g -O3 -diag-disable=10441)
  set(METRIS_C_FLAGS_RELEASE   -O3)
  #set(METRIS_C_FLAGS_RELEASE   -fno-alias -fno-fnalias -fast -fno-protect-parens -Ofast -flto -diag-disable=10441 -qopt-subscript-in-range)
  set(METRIS_C_FLAGS_DEBUG ${METRIS_CXX_FLAGS_DEBUG})
  message(Using Intel compiler ${CMAKE_C_COMPILER} ${SHORT_COMPILER_NAME})
elseif(CMAKE_CXX_COMPILER_ID STREQUAL GNU)
  message("Using GNU compiler ${CMAKE_C_COMPILER} ${SHORT_COMPILER_NAME}")
  set(METRIS_CXX_FLAGS_RELEASE  -DNDEBUG -O3 -fPIE)
  #set(METRIS_CXX_FLAGS_DEBUG   -Og -ggdb3 -Wall -Wextra -pedantic -DBOOST_STACKTRACE_USE_ADDR2LINE -march=native -no-pie -fno-pie  -rdynamic) # -S -fverbose-asm
  set(METRIS_CXX_FLAGS_DEBUG   -fconstexpr-ops-limit=10000000 -Og -ggdb3 -Wall -DBOOST_STACKTRACE_USE_ADDR2LINE -march=native ) #  -rdynamic # -S -fverbose-asm
  set(METRIS_C_FLAGS_RELEASE  ${METRIS_CXX_FLAGS_RELEASE})
  set(METRIS_C_FLAGS_DEBUG ${METRIS_CXX_FLAGS_DEBUG})
elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  message("Using Clang ${CMAKE_C_COMPILER} ${SHORT_COMPILER_NAME}")
  set(METRIS_CXX_FLAGS_RELEASE -DNDEBUG  -O3 )
  #set(METRIS_CXX_FLAGS_DEBUG   -Og -ggdb3 -Wall -Wextra -pedantic -DBOOST_STACKTRACE_USE_ADDR2LINE -march=native -no-pie -fno-pie  -rdynamic) # -S -fverbose-asm
  #set(METRIS_CXX_FLAGS_DEBUG  -fsanitize=address  -fconstexpr-steps=10000000 -O0 -g3 -DBOOST_STACKTRACE_USE_ADDR2LINE -march=native -fno-pie ) # -S -fverbose-asm
  if(USE_SANITIZER)
    message(WARNING "Using sanitizer")
    set(METRIS_CXX_FLAGS_DEBUG -Wall -Wno-logical-op-parentheses -fsanitize=address -Og -g3) # -S -fverbose-asm
  else()
    set(METRIS_CXX_FLAGS_DEBUG -Wall -Wno-logical-op-parentheses -Og -g3) # -S -fverbose-asm
  endif()
  #-fno-omit-frame-pointer
  set(METRIS_C_FLAGS_RELEASE  ${METRIS_CXX_FLAGS_RELEASE})
  set(METRIS_C_FLAGS_DEBUG ${METRIS_CXX_FLAGS_DEBUG})
else()
  message(FATAL_ERROR "Unknown compiler ID = ${CMAKE_CXX_COMPILER_ID}, SHORT_COMPILER_NAME = ${SHORT_COMPILER_NAME}")
endif()


if(PREPRO STREQUAL "True")
  message(WARNING "preprocessor mode")
  set(METRIS_CXX_FLAGS_RELEASE  ${METRIS_CXX_FLAGS_RELEASE} -E )
  set(METRIS_CXX_FLAGS_DEBUG ${METRIS_CXX_FLAGS_DEBUG} -E )
endif()