
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
  const int power = msh.param->opt_power;
  METRIS_ASSERT(power == 1 || power == -1);


  ftype tra, det;
  quafun_tradet<MFT,gdim,tdim,ftype>(msh,asdmsh,asdmet,ent2pol,bary,
                                     met_,&tra,&det);

  ftype quent;
  if constexpr (tdim == 2){
    if(power > 0){
      quent = tra*tra*(1.+1./(det*det))/8.;
    }else{
      quent = 8./(tra*tra*(1.+1./(det*det)));
    }
  }else{
    if(power > 0){
      quent = tra*tra*tra*(1.+1./(det*det))/54.;
    }else{
      quent = 54./(tra*tra*tra*(1.+1./(det*det)));
    }
  }

  return quent;
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

  const int power = msh.param->opt_power;
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
       &det,ddet,hdet);

  // This is used later on -> store it
  int dpowd = iipow<tdim>(tdim);              // n^n
  ftype trapowdm2 = irpow<tdim-2,ftype>(tra); // tra^(n-2)
  ftype trapowdm1 = trapowdm2*tra;            // tra^(n-1)
  ftype trapowd   = trapowdm1*tra;            // tra^n

  ftype quael1;
  quael1 = trapowd*(1. + 1./(det*det))/(2.*dpowd);

  ftype quael;
  const int ppower = abs(power);
  if(ivar < 0)
  {
    if (power < 0)
      quael1 = 1./quael1;

    quael = pow(quael1, ppower);

    return quael;
  }
  // From here, we compute derivatives.

  // Quality function base is 1/(2*d^d) * tra^d * (1 + 1/det^2)
  // derivative then:
  // 1/(2*d^d) * [ d * tra^(d-1) * dtra * (1 + 1/det^2) - 2 * tra^d/det^3 * ddet ]
  // hessian is to cumbersome to write...
  for(int ii = 0; ii < gdim; ii++){
    dquael[ii] = 1./(2.*dpowd) * (tdim * trapowdm1 * dtra[ii] * (1. + 1./(det*det))
                                - 2. * trapowd / (det*det*det) * ddet[ii]
                               );
  }
  if(hquael != NULL){ // asking for hessian as well
    for(int ii = 0; ii < gdim; ii++){
      for(int jj = ii; jj < gdim; jj++){
        hquael[sym2idx(ii,jj)] = 1./(2.*dpowd) * ( tdim * ( (tdim-1) * trapowdm2 * dtra[ii] * dtra[jj] * (1. + 1./(det*det))
                                                          + trapowdm1 * ( htra[sym2idx(ii,jj)] * (1. + 1./(det*det)) - 2. * dtra[ii] / (det*det*det) * ddet[jj])
                                                        )
                                                 -2 * ( tdim * trapowdm1 * dtra[jj] / (det*det*det) * ddet[ii] + trapowd * ( -3./(det*det*det*det) * ddet[ii]*ddet[jj]
                                                                                                                             + hdet[sym2idx(ii,jj)]/(det*det*det)
                                                                                                                           )
                                                      )
                                               );
      }
    }
  }
  if(power < 0){

    if(hquael != NULL){
      // modify Hessian first, as it uses the gradient
      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++)
          hquael[sym2idx(ii,jj)] = 2. / (quael1*quael1*quael1) * dquael[ii] * dquael[jj] - 1. / (quael1*quael1) * hquael[sym2idx(ii,jj)];
      }
    }

    // now modify the gradient
    for(int ii = 0; ii < gdim; ii++)
      dquael[ii] = -dquael[ii]/(quael1*quael1);

    // finaly modify value itself
    quael1 = 1./quael1;
  }

  quael = quael1;

  // if we have the final quality as the base Q raised to a positive power
  // modify value and derivatives accordingly
  if(ppower != 1){
    // Q^(abs(power)) , quael is Q^p, quael1 is just Q

    quael = pow(quael, ppower);

    if(hquael != NULL){
      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++)
          hquael[sym2idx(ii,jj)] = ppower * (ppower-1) * quael / (quael1*quael1) * dquael[ii] * dquael[jj] + ppower * quael/quael1 * hquael[sym2idx(ii,jj)];
      }
    }

    for(int ii = 0; ii < gdim; ii++){
      dquael[ii] = ppower * dquael[ii] * quael/quael1;
    }
  }

  return quael;
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