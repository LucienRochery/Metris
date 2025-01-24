add_library(SANS tools/SANSException.cpp
				    tools/demangle.cpp
            tools/BackTraceException.cpp)
target_include_directories(libkdtree INTERFACE 
                             ${Boost_INCLUDE_DIRS}
                             ${Boost_INCLUDE_DIR})
