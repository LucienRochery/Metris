//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_ANASOL__
#define __METRIS_MSH_ANASOL__


#include <boost/preprocessor/iteration/local.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>

#include <initializer_list>

#define MAX_ANASOL3D_DEFINED 1
#define MAX_ANASOL2D_DEFINED 5

#define MAX_ANASOL_DEFINED(dim) (dim == 2 ? MAX_ANASOL2D_DEFINED : MAX_ANASOL3D_DEFINED)


namespace Metris{

// Number of derivatives is specified by the std::initializer_list length. 
// This holds d1f, d2f...
typedef double(*anasol_proto)(void*,const double*__restrict__,std::initializer_list<double*>);


// A wrapper to ensure a function pointer passed to the overloaded function in MetrisParameters 
// is never interpreted as an integer. 
// This has no use to anybody except MetrisParameters. 
class AnaSolFun{
public:
  explicit AnaSolFun(anasol_proto anasol_ptr){
    this->anasol_ptr = anasol_ptr;
  }
  operator int() const = delete;
  operator anasol_proto() const {
    return anasol_ptr;
  }
  operator anasol_proto() {
    return anasol_ptr;
  }
  anasol_proto anasol_ptr;
};

#define BOOST_PP_LOCAL_MACRO(n)\
double anasol3D_##n(void*,const double*__restrict__,std::initializer_list<double*>);
#define BOOST_PP_LOCAL_LIMITS (1, MAX_ANASOL3D_DEFINED)
#include BOOST_PP_LOCAL_ITERATE()

#define BOOST_PP_LOCAL_MACRO(n)\
double anasol2D_##n(void*,const double*__restrict__,std::initializer_list<double*>);
#define BOOST_PP_LOCAL_LIMITS (1, MAX_ANASOL2D_DEFINED)
#include BOOST_PP_LOCAL_ITERATE()


constexpr anasol_proto __ANASOL2D[MAX_ANASOL2D_DEFINED] ={
  anasol2D_1
  #if MAX_ANASOL2D_DEFINED > 1
  #define BOOST_PP_LOCAL_MACRO(n) ,anasol2D_##n
  #define BOOST_PP_LOCAL_LIMITS (2, MAX_ANASOL2D_DEFINED)
  #include BOOST_PP_LOCAL_ITERATE()
  #endif
};

constexpr anasol_proto __ANASOL3D[MAX_ANASOL3D_DEFINED] ={
  anasol3D_1
  #if MAX_ANASOL3D_DEFINED > 1
  #define BOOST_PP_LOCAL_MACRO(n) ,anasol3D_##n
  #define BOOST_PP_LOCAL_LIMITS (2, MAX_ANASOL3D_DEFINED)
  #include BOOST_PP_LOCAL_ITERATE()
  #endif
};


} // End namespace

#endif
