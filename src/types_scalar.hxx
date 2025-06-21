//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#ifdef USE_MULTIPRECISION
  #ifdef USE_GMP
    // Boost MP promises greater performance with GMP
    #include <boost/multiprecision/gmp.hpp>
    #include <boost/multiprecision/cpp_bin_float.hpp>
    typedef  boost::multiprecision::number<boost::multiprecision::gmp_float<64>> float8;
    typedef  long double float4;
    //typedef  boost::multiprecision::number<boost::multiprecision::gmp_float<0>> float4;
    //typedef  boost::multiprecision::cpp_bin_float_double_extended float4;
    //typedef  boost::multiprecision::number<boost::multiprecision::gmp_float<12>> float4;
    //typedef  boost::multiprecision::mpf_float_100 float8;
    //typedef  boost::multiprecision::mpf_float_50 float4;
  #else
    #include <boost/multiprecision/cpp_bin_float.hpp>
    typedef  boost::multiprecision::cpp_bin_float_oct  float8;
    typedef  boost::multiprecision::cpp_bin_float_quad float4;
  #endif

  #define QUA_FTYPE_SEQ (double)(float4)(float8)
#else
  #define QUA_FTYPE_SEQ (double)
#endif