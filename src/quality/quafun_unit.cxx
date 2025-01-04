//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "quafun_unit.hxx"
#include "quafun_tradet.hxx"

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../linalg/sym3idx.hxx"

#include "../metris_constants.hxx"
#include "../aux_pp_inc.hxx"

namespace Metris{

#define UNIT_QUAFUN_EXP

// For some special barys (nodes), met is already known -> pass it in
template <class MFT, int gdim, int tdim, typename ftype>
ftype quafun_unit(Mesh<MFT> &msh,
                AsDeg asdmsh, AsDeg asdmet,
                const int*__restrict__ ent2poi,  
                const double*__restrict__ bary, int power, 
                double*__restrict__ met_){

  static_assert(gdim == 2 || gdim == 3);
  static_assert(tdim <= gdim); 

  METRIS_ASSERT(gdim == msh.idim);


  ftype tra, det;
  quafun_tradet<MFT,gdim,tdim,ftype>
              (msh,asdmsh,asdmet,ent2poi,bary,power,met_,&tra,&det);



  double coef_det = msh.param->opt_coef_det; 
  double coef_tra = msh.param->opt_coef_tra; 

  #ifndef UNIT_QUAFUN_EXP
    int    powr_det = msh.param->opt_powr_det; 
    int    powr_tra = msh.param->opt_powr_tra; 
  #else
    int    powr_det = abs(msh.param->opt_powr_det); 
    int    powr_tra = abs(msh.param->opt_powr_tra); 
  #endif

  ftype quael;

  #ifndef UNIT_QUAFUN_EXP
    quael = coef_tra*pow(tra - tdim,powr_tra)
          + coef_det*pow(det -    1,powr_det);
  #else
    ftype x, y; 
    if(tra > tdim){
      x = tra - tdim;
    }else{
      x = 1.0 / tra - 1.0 / tdim;
    }
    if(det > 1){
      y = det - 1;
    }else{
      y = 1.0 / det - 1.0;
    }
    //printf("Debug x = %15.7e y = %15.7e tra = %15.7e det = %15.7e \n",
    //  (double)x,(double)y,(double)tra,(double)det);
    quael = coef_tra * exp(powr_tra * x)
          + coef_det * exp(powr_det * y);
  #endif

  return quael;
}


#define EXPAND_TEMPLATE(r,SEQ) \
                  INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ),\
                              BOOST_PP_SEQ_ELEM(1, SEQ))
#define MFT_SEQ (MetricFieldFE)(MetricFieldAnalytical) 


#define INSTANTIATE(MFT_VAL,FTYPE)\
template FTYPE quafun_unit< MFT_VAL , 2, 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh, AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2poi, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met); \
template FTYPE quafun_unit< MFT_VAL , 3, 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh, AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2poi, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met); \
template FTYPE quafun_unit< MFT_VAL , 3, 3, FTYPE>\
                  (Mesh< MFT_VAL > &msh, AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2poi, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met); 
BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TEMPLATE,(MFT_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE
#undef EXPAND_TEMPLATE






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
ftype d_quafun_unit(Mesh<MFT> &msh,
                    AsDeg asdmsh, AsDeg asdmet,
                    const int* ent2poi,
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
      (msh,asdmsh,asdmet,
       ent2poi,bary,power,ivar,dofbas,idifmet,
       &tra,dtra,htra,
       &det,ddet,hdet);


  double coef_det = msh.param->opt_coef_det; 
  double coef_tra = msh.param->opt_coef_tra; 

  #ifndef UNIT_QUAFUN_EXP
    int    powr_det = msh.param->opt_powr_det; 
    int    powr_tra = msh.param->opt_powr_tra; 
  #else
    int    powr_det = abs(msh.param->opt_powr_det); 
    int    powr_tra = abs(msh.param->opt_powr_tra); 
  #endif


  ftype quael;

  #ifndef UNIT_QUAFUN_EXP
    quael = coef_tra*pow(tra - tdim,powr_tra)
          + coef_det*pow(det -    1,powr_det);
  #else
    ftype x, y; 
    if(tra > tdim){
      x = tra - tdim;
    }else{
      x = 1.0 / tra - 1.0 / tdim;
    }
    if(det > 1){
      y = det - 1;
    }else{
      y = 1.0 / det - 1;
    }
    ftype exptra = exp(powr_tra * x);
    ftype expdet = exp(powr_det * y);
    quael = coef_tra * exptra
          + coef_det * expdet;
  #endif

  // From here, we compute derivatives. 
  if(ivar < 0) return quael;
  // See docs/quality/qualityiff.pdf for details 

  #ifndef UNIT_QUAFUN_EXP
    for(int ii = 0; ii < gdim; ii++){
      dquael[ii] = coef_tra*powr_tra*dtra[ii]*pow(tra-tdim,powr_tra-1) 
                 + coef_det*powr_det*ddet[ii]*pow(det-   1,powr_det-1);
    }
  #else
    ftype dx[gdim], dy[gdim];
    if(tra > tdim){
      for(int ii = 0; ii < gdim; ii++) dx[ii] = dtra[ii];
    }else{
      for(int ii = 0; ii < gdim; ii++) dx[ii] = -dtra[ii]/(tra*tra);
    }
    if(det > 1){
      for(int ii = 0; ii < gdim; ii++) dy[ii] = ddet[ii];
    }else{
      for(int ii = 0; ii < gdim; ii++) dy[ii] = -ddet[ii]/(det*det);
    }

    for(int ii = 0; ii < gdim; ii++){
      dquael[ii] = coef_tra*powr_tra*dx[ii]*exptra
                 + coef_det*powr_det*dy[ii]*expdet;
    }
  #endif

  if(hquael == NULL) return quael;

  #ifndef UNIT_QUAFUN_EXP
    ftype trapowdm2 = pow(tra-tdim,powr_tra-2);
    //ftype trapowdm1 = (tra-tdim)*trapowdm2;

    ftype detpowdm2 = pow(det-1   ,powr_det-2);
    //ftype detpowdm1 = (det-tdim)*detpowdm2;

    for(int ii = 0; ii < gdim; ii++){
      for(int jj = ii; jj < gdim; jj++){
        hquael[sym3idx(ii,jj)] = 
            powr_tra*coef_tra*( htra[sym3idx(ii,jj)]*(tra-tdim)
                     + (powr_tra-1)*dtra[jj]*dtra[ii] )*trapowdm2
          + powr_det*coef_det*( hdet[sym3idx(ii,jj)]*(det-1   )
                     + (powr_det-1)*ddet[jj]*ddet[ii] )*detpowdm2;
      }
    }
  #else
    ftype hx[nhess], hy[nhess];
    if(tra > tdim){
      for(int ii = 0; ii < gdim; ii++)
        for(int jj = ii; jj < gdim; jj++)
          hx[sym3idx(ii,jj)] = dtra[ii]*dtra[jj] + htra[sym3idx(ii,jj)];
    }else{
      for(int ii = 0; ii < gdim; ii++)
        for(int jj = ii; jj < gdim; jj++)
          hx[sym3idx(ii,jj)] = (  -htra[sym3idx(ii,jj)]
                                + 2*dtra[ii]*dtra[jj]/tra )/(tra*tra);
    }
    if(det > 1){
      for(int ii = 0; ii < gdim; ii++)
        for(int jj = ii; jj < gdim; jj++)
          hy[sym3idx(ii,jj)] = ddet[ii]*ddet[jj] + hdet[sym3idx(ii,jj)];
    }else{
      for(int ii = 0; ii < gdim; ii++)
        for(int jj = ii; jj < gdim; jj++)
          hy[sym3idx(ii,jj)] = (  -hdet[sym3idx(ii,jj)]
                                + 2*ddet[ii]*ddet[jj]/det)/(det*det);
    }
    for(int ii = 0; ii < gdim; ii++){
      for(int jj = ii; jj < gdim; jj++){
        hquael[sym3idx(ii,jj)] = 
          coef_tra*powr_tra*( hx[sym3idx(ii,jj)]
                            + dx[ii]*dx[jj]*powr_tra )*exptra 
        + coef_det*powr_det*( hy[sym3idx(ii,jj)]
                            + dy[ii]*dy[jj]*powr_det )*expdet ;
      }
    }
  #endif


  return quael;
}



// While cumbersome, this replaces a bunch of manual instantiations, about to 
// be made worse the day we add tdimn as a template argument. 
#define EXPAND_TEMPLATE(z,gdim,SEQ) \
                  INSTANTIATE(gdim,BOOST_PP_SEQ_ELEM(0, SEQ),\
                                   BOOST_PP_SEQ_ELEM(1, SEQ))
#define REPEAT_GDIM(r,SEQ) BOOST_PP_REPEAT(2,EXPAND_TEMPLATE,SEQ)

#define INSTANTIATE(gdim,MFT_VAL,FTYPE)\
template FTYPE d_quafun_unit< MFT_VAL , 2+gdim, FTYPE>\
                  (Mesh< MFT_VAL > &msh, AsDeg asdmsh, AsDeg asdmet,\
                   const int* ent2poi, \
                   const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, \
                   FTYPE*__restrict__ hquael); 
BOOST_PP_SEQ_FOR_EACH_PRODUCT(REPEAT_GDIM,(MFT_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE
#undef EXPAND_TEMPLATE
#undef REPEAT_GDIM

#undef MFT_SEQ 



}