
//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "quafun_sizeshape.hxx"
#include "quafun_tradet.hxx"

#include "../Mesh/Mesh.hxx"
#include "../metris_constants.hxx"
#include "../utils/aux_misc.hxx"

#include "../linalg/symidx.hxx"

#include "../utils/aux_pp_inc.hxx"

#include <limits>

namespace Metris{

// For some special barys (nodes), met is already known -> pass it in
template <class MFT, int gdim, int tdim, typename ftype>
ftype quafun_sizeshape(Mesh<MFT> &msh,
                        AsDeg asdmsh, AsDeg asdmet,
                        const int*__restrict__ ent2pol,
                        const double*__restrict__ bary,
                        const double*__restrict__ met_){

  static_assert(gdim == 2 || gdim == 3);
  static_assert(tdim <= gdim);

  METRIS_ASSERT(gdim == msh.idim);
  ftype tra, det;
  quafun_tradet<MFT,gdim,tdim,ftype>(msh,asdmsh,asdmet,ent2pol,bary,
                                     met_,&tra,&det,
                                     QualitySingularityPolicy::Reject);

  ftype size_shape_quality;
  if constexpr (tdim == 2){
    size_shape_quality = tra*tra*(1.+1./(det*det))/8.;
  }else{
    size_shape_quality = tra*tra*tra*(1.+1./(det*det))/54.;
  }

  ftype size_shape_error = size_shape_quality - 1.;
  constexpr double ideal_roundoff_tolerance
      = 32.0*std::numeric_limits<double>::epsilon();
  if(abs(size_shape_error) <= ideal_roundoff_tolerance){
    size_shape_error = 0.;
  }else if(size_shape_error < 0.){
    METRIS_THROW_MSG(
        "SizeShape quality below its ideal minimum: {:e}",
        size_shape_quality);
  }

  const double objective_p = msh.param->objective_p;
  METRIS_ASSERT(objective_p >= 1.0);
  if(objective_p == 1.0) return size_shape_error;
  return pow(size_shape_error,objective_p);
}


#define EXPAND_TEMPLATE(r,SEQ) \
                  INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ),\
                              BOOST_PP_SEQ_ELEM(1, SEQ))
#define MFT_SEQ (MetricFieldFE)(MetricFieldAnalytical)


#define INSTANTIATE(MFT_VAL,FTYPE)\
template FTYPE quafun_sizeshape< MFT_VAL , 2, 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary,\
                   const double*__restrict__ met); \
template FTYPE quafun_sizeshape< MFT_VAL , 3, 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary,\
                   const double*__restrict__ met); \
template FTYPE quafun_sizeshape< MFT_VAL , 3, 3,  FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary,\
                   const double*__restrict__ met);
BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TEMPLATE,\
                              (MFT_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE
#undef EXPAND_TEMPLATE

/*
NOTE: This does not yet use the metric field derivatives.
That is something we might not want, it'll be costlier as we can no longer
provide a metric, but instead must interpolate (which dominates CPU time
because of the matrix exponential).
Compute quality function and derivative with respect to node ivar
- gdim is geometric dimension: also topological dimension !
- mshdeg is mesh degree
- asdmet is AsDeg::P1 or AsDeg::Pk: MetricField handles its degree
- ftype is arithmetic precision (debug): double, float8...
- ivar is the DoF, specified Bézier or Lagrange by:
- dofbas is either FEBasis::Lagrange or FEBasis::Bezier -> whether the DoF is a
  Lagrange node or Bézier control point
- idifmet is either DifVar::None ("static" metric approximation) or DifVar::Phys
- quael is output quality and derivatives
*/
template <class MFT, int gdim, int tdim, typename ftype>
ftype d_quafun_sizeshape(Mesh<MFT> &msh,
                  AsDeg asdmsh, AsDeg asdmet,
                  const int* ent2pol,
                  const double*__restrict__ bary,
                  const double*__restrict__ met_,
                  int ivar,
                  FEBasis dofbas,
                  DifVar idifmet,
                  ftype*__restrict__ dquael,
                  ftype*__restrict__ hquael){

  constexpr int nhess = (gdim*(gdim+1))/2;
  ftype tra, det;
  ftype dtra[gdim], htra_[nhess];
  ftype ddet[gdim], hdet_[nhess];

  // This is not needed for the gradient as ivar pilots that
  ftype *htra = NULL, *hdet = NULL;
  if(hquael != NULL){
    htra = htra_;
    hdet = hdet_;
  }

  // get derivative of trace and determinant
  d_quafun_tradet<MFT,gdim,tdim,ftype>
      (msh,asdmsh,asdmet,ent2pol,bary,ivar,dofbas,idifmet,
       met_,
       &tra,dtra,htra,
       &det,ddet,hdet,
       QualitySingularityPolicy::Reject);

  const int dimension_power = iipow<tdim>(tdim);
  const ftype trace_power_dimension_minus_two
      = irpow<tdim-2,ftype>(tra);
  const ftype trace_power_dimension_minus_one
      = trace_power_dimension_minus_two*tra;
  const ftype trace_power_dimension
      = trace_power_dimension_minus_one*tra;

  const ftype size_shape_quality
      = trace_power_dimension*(1. + 1./(det*det))
       /(2.*dimension_power);
  ftype size_shape_error = size_shape_quality - 1.;
  constexpr double ideal_roundoff_tolerance
      = 32.0*std::numeric_limits<double>::epsilon();
  if(abs(size_shape_error) <= ideal_roundoff_tolerance){
    size_shape_error = 0.;
  }else if(size_shape_error < 0.){
    METRIS_THROW_MSG(
        "SizeShape quality below its ideal minimum: {:e}",
        size_shape_quality);
  }

  const double objective_p = msh.param->objective_p;
  METRIS_ASSERT(objective_p >= 1.0);
  const ftype psi = objective_p == 1.0
                  ? size_shape_error
                  : pow(size_shape_error,objective_p);
  if(ivar < 0) return psi;

  // Quality function base is 1/(2*d^d) * tra^d * (1 + 1/det^2)
  // derivative then:
  // 1/(2*d^d) * [ d * tra^(d-1) * dtra * (1 + 1/det^2) - 2 * tra^d/det^3 * ddet ]
  // hessian is to cumbersome to write...
  for(int ii = 0; ii < gdim; ii++){
    dquael[ii]
        = 1./(2.*dimension_power)
        * (tdim*trace_power_dimension_minus_one*dtra[ii]
                     *(1. + 1./(det*det))
           - 2.*trace_power_dimension/(det*det*det)*ddet[ii]);
  }
  if(hquael != NULL){ // asking for hessian as well
    for(int ii = 0; ii < gdim; ii++){
      for(int jj = ii; jj < gdim; jj++){
        hquael[sym2idx(ii,jj)]
            = 1./(2.*dimension_power)
            * (tdim
               *((tdim - 1)*trace_power_dimension_minus_two
                              *dtra[ii]*dtra[jj]
                              *(1. + 1./(det*det))
                 + trace_power_dimension_minus_one
                              *(htra[sym2idx(ii,jj)]
                                   *(1. + 1./(det*det))
                                - 2.*dtra[ii]/(det*det*det)*ddet[jj]))
               - 2.*(tdim*trace_power_dimension_minus_one*dtra[jj]
                              /(det*det*det)*ddet[ii]
                      + trace_power_dimension
                              *(-3./(det*det*det*det)*ddet[ii]*ddet[jj]
                                + hdet[sym2idx(ii,jj)]/(det*det*det))));
      }
    }
  }

  if(objective_p == 1.0) return psi;

  if(size_shape_error == 0.){
    for(int ii = 0; ii < gdim; ii++){
      dquael[ii] = 0.;
    }
    if(hquael != NULL){
      for(int ihessian = 0; ihessian < nhess; ihessian++){
        hquael[ihessian] = 0.;
      }
    }
    return psi;
  }

  ftype gradient_scale;
  ftype gradient_outer_product_scale;
  if(objective_p == 2.0){
    gradient_scale = 2.*size_shape_error;
    gradient_outer_product_scale = 2.;
  }else{
    gradient_scale
        = objective_p*pow(size_shape_error,objective_p - 1.0);
    gradient_outer_product_scale
        = objective_p*(objective_p - 1.0)
         *pow(size_shape_error,objective_p - 2.0);
  }

  if(hquael != NULL){
    for(int ii = 0; ii < gdim; ii++){
      for(int jj = ii; jj < gdim; jj++){
        const int ihessian = sym2idx(ii,jj);
        hquael[ihessian]
            = gradient_outer_product_scale*dquael[ii]*dquael[jj]
            + gradient_scale*hquael[ihessian];
      }
    }
  }

  for(int ii = 0; ii < gdim; ii++){
    dquael[ii] *= gradient_scale;
  }

  return psi;
}

// While cumbersome, this replaces a bunch of manual instantiations, about to
// be made worse the day we add tdimn as a template argument.
#define EXPAND_TEMPLATE(z,SEQ) \
                  INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ),\
                              BOOST_PP_SEQ_ELEM(1, SEQ))

#define INSTANTIATE(MFT_VAL,FTYPE)\
template FTYPE d_quafun_sizeshape< MFT_VAL , 2, 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int* ent2pol, \
                   const double*__restrict__ bary,\
                   const double*__restrict__ met_,\
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, \
                   FTYPE*__restrict__ hquael);\
template FTYPE d_quafun_sizeshape< MFT_VAL , 3, 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int* ent2pol, \
                   const double*__restrict__ bary,\
                   const double*__restrict__ met_,\
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, \
                   FTYPE*__restrict__ hquael);\
template FTYPE d_quafun_sizeshape< MFT_VAL , 3, 3, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int* ent2pol, \
                   const double*__restrict__ bary,\
                   const double*__restrict__ met_,\
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, \
                   FTYPE*__restrict__ hquael);
BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TEMPLATE,\
                              (MFT_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE
#undef EXPAND_TEMPLATE

#undef MFT_SEQ
#undef ASDEG_SEQ

#if 0
// While cumbersome, this replaces a bunch of manual instantiations, about to
// be made worse the day we add tdimn as a template argument.
#define EXPAND_TEMPLATE(z,gdim,SEQ) \
                  INSTANTIATE(z,gdim,BOOST_PP_SEQ_ELEM(0, SEQ),\
                                     BOOST_PP_SEQ_ELEM(1, SEQ),\
                                     BOOST_PP_SEQ_ELEM(2, SEQ),\
                                     BOOST_PP_SEQ_ELEM(3, SEQ))
#define REPEAT_GDIM(z,n,SEQ) BOOST_PP_REPEAT(2,EXPAND_TEMPLATE,(n)SEQ)
#define REPEAT_IDEG(r,SEQ)   BOOST_PP_REPEAT(METRIS_MAX_DEG,REPEAT_GDIM,SEQ)

#define INSTANTIATE(z,gdim,ideg,MFT_VAL,ASDMET,FTYPE)\
template FTYPE d_quafun_sizeshape< MFT_VAL , 2+gdim, 1+ideg, ASDMET, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   const int* ent2pol, \
                   const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, \
                   FTYPE*__restrict__ hquael);
BOOST_PP_SEQ_FOR_EACH_PRODUCT(REPEAT_IDEG,(MFT_SEQ)(ASDEG_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE
#undef EXPAND_TEMPLATE
#undef REPEAT_GDIM
#undef REPEAT_IDEG

#undef MFT_SEQ
#undef ASDEG_SEQ
#undef QUA_FTYPE_SEQ
#endif


}
