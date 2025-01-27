add_library(SANS tools/SANSException.cpp
				         tools/demangle.cpp
                 tools/BackTraceException.cpp)
target_include_directories(SANS INTERFACE ${Boost_INCLUDE_DIRS})
