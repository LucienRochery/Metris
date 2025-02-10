//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_MSH_ANAMET__
#define __METRIS_MSH_ANAMET__


#include <boost/preprocessor/iteration/local.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>


#define MAX_ANAMET3D_DEFINED 4
#define MAX_ANAMET2D_DEFINED 6

#define MAX_ANAMET_DEFINED(dim) (dim == 2 ? MAX_ANAMET2D_DEFINED : MAX_ANAMET3D_DEFINED)
// Allowed prototype: 
//  void (*anamet)(void* ctx, double *crd, int idif1, double *met, double *dmet);
// Met stored 1 2 4  
//              3 5
//                6
// Dmet [i][j] = d_i M_j

namespace Metris{

typedef void(*anamet_proto)(void*,const double*__restrict__,double,int,double*,double*);


// A wrapper to ensure a function pointer passed to the overloaded function in MetrisParameters 
// is never interpreted as an integer. 
// This has no use to anybody except MetrisParameters. 
class AnaMetFun{
public:
  explicit AnaMetFun(anamet_proto anamet_ptr){
    this->anamet_ptr = anamet_ptr;
  }
  operator int() const = delete;
  operator anamet_proto() const {
    return anamet_ptr;
  }
  operator anamet_proto() {
    return anamet_ptr;
  }
  anamet_proto anamet_ptr;
};

#define BOOST_PP_LOCAL_MACRO(n)\
void anamet3D_##n(void *ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet);
#define BOOST_PP_LOCAL_LIMITS (1, MAX_ANAMET3D_DEFINED)
#include BOOST_PP_LOCAL_ITERATE()

#define BOOST_PP_LOCAL_MACRO(n)\
void anamet2D_##n(void *ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet);
#define BOOST_PP_LOCAL_LIMITS (1, MAX_ANAMET2D_DEFINED)
#include BOOST_PP_LOCAL_ITERATE()



constexpr anamet_proto __ANAMET2D[MAX_ANAMET2D_DEFINED] ={
  anamet2D_1
  #define BOOST_PP_LOCAL_MACRO(n) ,anamet2D_##n
  #define BOOST_PP_LOCAL_LIMITS (2, MAX_ANAMET2D_DEFINED)
  #include BOOST_PP_LOCAL_ITERATE()
};

constexpr anamet_proto __ANAMET3D[MAX_ANAMET3D_DEFINED] ={
  anamet3D_1
  #define BOOST_PP_LOCAL_MACRO(n) ,anamet3D_##n
  #define BOOST_PP_LOCAL_LIMITS (2, MAX_ANAMET3D_DEFINED)
  #include BOOST_PP_LOCAL_ITERATE()
};


} // End namespace

#endif
