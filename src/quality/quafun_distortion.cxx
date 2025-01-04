//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "quafun_distortion.hxx"
#include "quafun_tradet.hxx"

#include "../Mesh/Mesh.hxx"
#include "../metris_constants.hxx"
#include "../aux_utils.hxx"

#include "../linalg/sym3idx.hxx"

#include "../aux_pp_inc.hxx"

namespace Metris{

// For some special barys (nodes), met is already known -> pass it in
template <class MFT, int gdim, int tdim, typename ftype>
ftype quafun_distortion(Mesh<MFT> &msh,
                        AsDeg asdmsh, AsDeg asdmet,
                        const int*__restrict__ ent2pol,  
                        const double*__restrict__ bary, int power, 
                        double*__restrict__ met_){

  static_assert(gdim == 2 || gdim == 3);
  static_assert(tdim <= gdim); 

  METRIS_ASSERT(gdim == msh.idim);

  METRIS_ASSERT(power != 0);

  ftype tra, det;
  #ifndef NDEBUG
  try{
  #endif
  quafun_tradet<MFT,gdim,tdim,ftype>(msh,asdmsh,asdmet,ent2pol,bary,
                                     power,met_,&tra,&det);
  #ifndef NDEBUG
  }catch(const MetrisExcept& e){
    printf("##quafun_distortion excpt ent2pol = \n");
    intAr1(facnpps[msh.curdeg],ent2pol).print();
    throw(e);
  }
  #endif


  ftype quent; 
  if constexpr (tdim == 2){
    if(power > 0){
      quent = pow(tra*tra/det/4,power);
    }else{
      quent = pow(4*det/(tra*tra),-power);
    }
  }else{
    if(power > 0){
      quent = pow(tra*tra*tra/det/27,power);
    }else{
      quent = pow(27*det/(tra*tra*tra),-power);
    }
  }

  return quent;
}


#define EXPAND_TEMPLATE(r,SEQ) \
                  INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ),\
                              BOOST_PP_SEQ_ELEM(1, SEQ))
#define MFT_SEQ (MetricFieldFE)(MetricFieldAnalytical) 


#define INSTANTIATE(MFT_VAL,FTYPE)\
template FTYPE quafun_distortion< MFT_VAL , 2, 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met); \
template FTYPE quafun_distortion< MFT_VAL , 3, 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met); \
template FTYPE quafun_distortion< MFT_VAL , 3, 3,  FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met); 
BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TEMPLATE,\
                              (MFT_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE
#undef EXPAND_TEMPLATE

#if 0
#define EXPAND_TEMPLATE(z,ideg,SEQ) \
                  INSTANTIATE(z,ideg,BOOST_PP_SEQ_ELEM(0, SEQ),\
                                     BOOST_PP_SEQ_ELEM(1, SEQ),\
                                     BOOST_PP_SEQ_ELEM(2, SEQ))
#define REPEAT_IDEG(r,SEQ)   BOOST_PP_REPEAT(METRIS_MAX_DEG,EXPAND_TEMPLATE,SEQ)
#define MFT_SEQ (MetricFieldFE)(MetricFieldAnalytical)
#define ASDEG_SEQ (AsDeg::P1)(AsDeg::Pk) 


#define INSTANTIATE(z,ideg,MFT_VAL,ASDMET,FTYPE)\
template FTYPE quafun_distortion< MFT_VAL , 2, 2, 1+ideg, ASDMET, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met); \
template FTYPE quafun_distortion< MFT_VAL , 3, 2, 1+ideg, ASDMET, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met); \
template FTYPE quafun_distortion< MFT_VAL , 3, 3, 1+ideg, ASDMET, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met); 
BOOST_PP_SEQ_FOR_EACH_PRODUCT(REPEAT_IDEG,(MFT_SEQ)(ASDEG_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE
#undef EXPAND_TEMPLATE
#undef REPEAT_IDEG
#endif





/*
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
template <class MFT, int gdim, typename ftype>
ftype d_quafun_distortion(Mesh<MFT> &msh, 
                  AsDeg asdmsh, AsDeg asdmet,
                  const int* ent2pol,
                  const double*__restrict__ bary, 
                  int power, 
                  int ivar, 
                  FEBasis dofbas, 
                  DifVar idifmet, 
                  ftype*__restrict__ dquael, 
                  ftype*__restrict__ hquael){


  constexpr int tdim  = gdim;
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

  d_quafun_tradet<MFT,gdim,ftype>
      (msh,asdmsh,asdmet,ent2pol,bary,power,ivar,dofbas,idifmet,
       &tra,dtra,htra,
       &det,ddet,hdet);


  // This is used later on -> store it 
  int dpowd = iipow<tdim>(tdim);
  ftype trapowdm2 = irpow<tdim-2,ftype>(tra);
  ftype trapowdm1 = trapowdm2*tra;
  ftype trapowd   = trapowdm1*tra;

  ftype quael;
  if(power > 0){
    quael = trapowd/(det*dpowd);
  }else{
    quael = (det*dpowd)/trapowd;
  }
  ftype quae1 = quael; // for derivatives
  quael = pow(quael, abs(power));



  // From here, we compute derivatives. 
  if(ivar < 0) return quael;
  // See docs/quality/qualityiff.pdf for details 


  // Quality function is ( d^d tra^d/det ) ^ power
  //int dpowdpowp = pow(dpowd,power); 
  if(power > 0){
    // We could use idpow<tdim-1,ftype>(tra) = trapowd/tra but is this as stable? 
    // idpow expands to tra*tra*... here either tra or tra*tra. This is not slower. 
    // We can factorize by idpow<tdim->,ftype>(tra)
    for(int ii = 0; ii < gdim; ii++){
      dquael[ii] = trapowdm1*( tdim*dtra[ii]*det
                             - tra*ddet[ii])
                  /(det*det*dpowd); 
    }
    if(hquael != NULL){
      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++){
          hquael[sym3idx(ii,jj)] = ( tdim*trapowdm2*det*det*( tra*htra[sym3idx(ii,jj)]
                                                            + (tdim - 1)*dtra[ii]*dtra[jj] )
                                   - trapowdm1*det*( tra*hdet[sym3idx(ii,jj)]
                                                   + tdim*dtra[ii]*ddet[jj]
                                                   + tdim*dtra[jj]*ddet[ii])
                                   + 2.0*trapowd*ddet[ii]*ddet[jj]
                                   )/(det*det*det);
          hquael[sym3idx(ii,jj)] /= dpowd;
        }
      }
    }
  }else{
    for(int ii = 0; ii < gdim; ii++){
      dquael[ii] = dpowd*(         tra   *ddet[ii]
                   - tdim*dtra[ii]*  det) 
                  /(trapowd*tra); 
    }
    if(hquael != NULL){
      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++){
          hquael[sym3idx(ii,jj)] = (-(tdim+1)*dtra[jj]*(ddet[ii]*tra - tdim*det*dtra[ii])
                                   + tra*( hdet[sym3idx(ii,jj)]*tra + ddet[ii]*dtra[jj] - tdim*ddet[jj]*dtra[ii]
                                         - tdim*det*htra[sym3idx(ii,jj)]) 
                                   )/(trapowd*tra*tra);
          hquael[sym3idx(ii,jj)] *= dpowd;
        }
      }
    } 
  }

  if(abs(power) != 1){
    // Q^(abs(power)) , quael is Q^p, quae1 is just Q
    for(int ii = 0; ii < gdim; ii++){
      dquael[ii] = abs(power)*dquael[ii]*quael/quae1;
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
template FTYPE d_quafun_distortion< MFT_VAL , 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int* ent2pol, \
                   const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, \
                   FTYPE*__restrict__ hquael);\
template FTYPE d_quafun_distortion< MFT_VAL , 3, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int* ent2pol, \
                   const double*__restrict__ bary, \
                   int power, \
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
template FTYPE d_quafun_distortion< MFT_VAL , 2+gdim, 1+ideg, ASDMET, FTYPE>\
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