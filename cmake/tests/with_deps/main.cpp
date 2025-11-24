#include <iostream>
#include "Metris.h"
#include <nlopt.h>
#include <Eigen/Dense>
// CMake will define BOOST_LIBRARIES and Boost::program_options if legacy compatibility works
#ifdef Boost_LIBRARIES
#pragma message("Boost_LIBRARIES is defined: " BOOST_LIBRARIES)
#endif
#ifdef Boost_INCLUDE_DIR
#pragma message("Boost_INCLUDE_DIR is defined: " Boost_INCLUDE_DIR)
#endif
int main(){
  std::cout << "Consumer with_deps: Metris + NLopt + Eigen included\n";
#ifdef Boost_LIBRARIES
  std::cout << "Boost_LIBRARIES: " << BOOST_LIBRARIES << std::endl;
#else
  std::cout << "Boost_LIBRARIES not defined" << std::endl;
#endif
#ifdef Boost_INCLUDE_DIR
  std::cout << "Boost_INCLUDE_DIR: " << Boost_INCLUDE_DIR << std::endl;
#else
  std::cout << "Boost_INCLUDE_DIR not defined" << std::endl;
#endif
  return 0;
}
