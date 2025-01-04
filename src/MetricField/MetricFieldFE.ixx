//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php




#if !BOOST_PP_IS_ITERATING
   #ifndef FILE_H_
      #define FILE_H_
      #include <boost/preprocessor/iteration/iterate.hpp>
      #define BOOST_PP_ITERATION_PARAMS_1 (3, (1, METRIS_MAX_DEG, "src/MetricField/MetricFieldFE.ixx"))
      #include BOOST_PP_ITERATE()
   #endif
#else
   #include <boost/preprocessor/repetition/repeat_from_to.hpp>
   #define INSTANTIATE1(z,n,t) \
template void MetricFieldFE::getMetBary0<2, 1+n, BOOST_PP_ITERATION()>\
                            (DifVar idiff,  MetSpace tarspac, \
                              const int*__restrict__ ent2pol, int tdimn, \
                             const double*__restrict__ bary,double*__restrict__ metl, double*__restrict__ dmet) ;\
template void MetricFieldFE::getMetBary0<3, 1+n, BOOST_PP_ITERATION()>\
                            ( DifVar idiff, MetSpace tarspac, \
                              const int*__restrict__ ent2pol, int tdimn,  \
                             const double*__restrict__ bary, double*__restrict__ metl, double*__restrict__ dmet) ;
   BOOST_PP_REPEAT_FROM_TO(0, BOOST_PP_ITERATION(), INSTANTIATE1, )
   // Mshdeg goes from 1 to METRIS_MAX_DEG - 1 (both included)
   // Tardeg goes from 2 to METRIS_MAX_DEG (both included)

   #undef INSTANTIATE1
#endif