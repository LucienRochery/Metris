//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "quafun_tradet.hxx"

#include "../Mesh/Mesh.hxx"
#include "../metris_constants.hxx"
#include "../aux_utils.hxx"

#include "../linalg/matprods.hxx"
#include "../linalg/det.hxx"

#include "../aux_pp_inc.hxx"

namespace Metris{

// For some special barys (nodes), met is already known -> pass it in
template <class MFT, int gdim, int tdim, typename ftype>
void quafun_tradet(Mesh<MFT> &msh,AsDeg asdmsh, AsDeg asdmet,
                   const int*__restrict__ ent2pol,  
                   const double*__restrict__ bary, int power, 
                   double*__restrict__ met_,
                   ftype*__restrict__ tra,
                   ftype*__restrict__ det){

  static_assert(gdim == 2 || gdim == 3);
  static_assert(tdim <= gdim); 

  METRIS_ASSERT(gdim == msh.idim);

  METRIS_ASSERT(power != 0);


  if(msh.met.getSpace() != MetSpace::Log) METRIS_THROW_MSG(WArgExcept(),
      "## SET MESH METRIC TO LOG BEFORE CALLING metqua2_xi");

  constexpr int nnmet = (gdim*(gdim+1))/2;

  double jmat[tdim*gdim],coopr[gdim];
  double met[nnmet];


  // Get Jacobian matrix at xi
  if(asdmsh == AsDeg::P1){
    if constexpr(tdim == 2){  
      eval2<gdim,1>(msh.coord,ent2pol,msh.getBasis(),
                       DifVar::Bary,DifVar::None,
                       bary,coopr,jmat,NULL);
    }else{
      eval3<gdim,1>(msh.coord,ent2pol,msh.getBasis(),
                       DifVar::Bary,DifVar::None,
                       bary,coopr,jmat,NULL);
    }
  }else{
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){

      if constexpr(tdim == 2){  
        eval2<gdim,ideg>(msh.coord,ent2pol,msh.getBasis(),
                         DifVar::Bary,DifVar::None,
                         bary,coopr,jmat,NULL);
      }else{
        eval3<gdim,ideg>(msh.coord,ent2pol,msh.getBasis(),
                         DifVar::Bary,DifVar::None,
                         bary,coopr,jmat,NULL);
      }

    }}CT_FOR1(ideg);
  }

  if(met_ == NULL){
    // Get metric interpolated at xi
    // The metric field class fetches geometric dimension from the mesh
    if(asdmet == AsDeg::P1){
      msh.met.getMetFullinfo(asdmet,DifVar::None,MetSpace::Exp,
                             ent2pol,tdim,bary,coopr,met,NULL);
    }else{
      msh.met.getMetFullinfo(asdmet,DifVar::None,MetSpace::Exp,
                             ent2pol,tdim,bary,coopr,met,NULL);
    }
  }else{
    for(int ii = 0; ii < nnmet; ii++) met[ii] = met_[ii];
  }

  // Compute J_0^{-T} J_K^T M J_K J_0^{-1}
  // Starting with J_K J_0^{-1}
  // Note that J_K is stored transposed w.r.t. above ! -> jmat[i,j] = d_i F_j 
  // whereas (J_K)_{ij} = d_j F_i .. 


  // Get J_0^{-T} J_K^T
  ftype invtJ_0tJ_K[tdim*gdim];

  matXmat<tdim,tdim,gdim>(Constants::invtJ_0[hana::type_c<ftype>][tdim],jmat,invtJ_0tJ_K);


  ftype J0tJtMJJ0_diag[tdim];
  matXsymXtmat_diag<gdim, tdim, double, ftype, ftype>(met, invtJ_0tJ_K, J0tJtMJJ0_diag);
  *tra = J0tJtMJJ0_diag[0] + J0tJtMJJ0_diag[1];
  if constexpr (tdim == 3) *tra += J0tJtMJJ0_diag[2];
  
  // This is an actual exception that should never theoretically happen. 
  if(*tra < 1.0e-16) METRIS_THROW(GeomExcept());


  if constexpr(tdim == gdim){
    ftype det1 = detmat<gdim>(invtJ_0tJ_K); // Error of 10^{-19} = 10^-14 relative compared to Matlab on wonky case
    //ftype tmp = detsym<gdim>(met); // Error of 10^5 ... // Matlab yields 7.346639223296765e+09 we get 7.346911345315383e+09 // Even in relative this is terrible
    ftype tmp = detsym2<gdim>(met); // Also wrong, same error. Note Matlab gets 10^-9 final quality relative error ! Our determinant is terribly bad
    *det = det1*det1*tmp; 
  }else{
    static_assert(tdim == 2);
    ftype J0tJtMJJ0[3];
    matXsymXtmat<2,3,double,ftype,ftype>(met,invtJ_0tJ_K,J0tJtMJJ0);
    *det = detsym2<2>(J0tJtMJJ0);
  }

  if(abs(*det) < 1.0e-16 && power > 0) 
     METRIS_THROW_MSG(GeomExcept(), "Singular J^TMJ det = "<<det);

   return;
}


#define EXPAND_TEMPLATE(r,SEQ) \
                  INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ),\
                              BOOST_PP_SEQ_ELEM(1, SEQ))
#define MFT_SEQ (MetricFieldFE)(MetricFieldAnalytical)


#define INSTANTIATE(MFT_VAL,FTYPE)\
template void quafun_tradet< MFT_VAL , 2, 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met,\
                   FTYPE*__restrict__ tra,\
                   FTYPE*__restrict__ det); \
template void quafun_tradet< MFT_VAL , 3, 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met,\
                   FTYPE*__restrict__ tra,\
                   FTYPE*__restrict__ det); \
template void quafun_tradet< MFT_VAL , 3, 3, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met,\
                   FTYPE*__restrict__ tra,\
                   FTYPE*__restrict__ det); 
BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TEMPLATE,\
                              (MFT_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE
#undef EXPAND_TEMPLATE
#undef REPEAT_IDEG

#if 0
#define EXPAND_TEMPLATE(z,ideg,SEQ) \
                  INSTANTIATE(z,ideg,BOOST_PP_SEQ_ELEM(0, SEQ),\
                                     BOOST_PP_SEQ_ELEM(1, SEQ),\
                                     BOOST_PP_SEQ_ELEM(2, SEQ))
#define REPEAT_IDEG(r,SEQ)   BOOST_PP_REPEAT(METRIS_MAX_DEG,EXPAND_TEMPLATE,SEQ)
#define MFT_SEQ (MetricFieldFE)(MetricFieldAnalytical)
#define ASDEG_SEQ (AsDeg::P1)(AsDeg::Pk) 


#define INSTANTIATE(z,ideg,MFT_VAL,ASDEG_VAL,FTYPE)\
template void quafun_tradet< MFT_VAL , 2, 2, 1+ideg, ASDEG_VAL, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met,\
                   FTYPE*__restrict__ tra,\
                   FTYPE*__restrict__ det); \
template void quafun_tradet< MFT_VAL , 3, 2, 1+ideg, ASDEG_VAL, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met,\
                   FTYPE*__restrict__ tra,\
                   FTYPE*__restrict__ det); \
template void quafun_tradet< MFT_VAL , 3, 3, 1+ideg, ASDEG_VAL, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   const int*__restrict__ ent2pol, \
                   const double*__restrict__ bary, int power,\
                   double*__restrict__ met,\
                   FTYPE*__restrict__ tra,\
                   FTYPE*__restrict__ det); 
BOOST_PP_SEQ_FOR_EACH_PRODUCT(REPEAT_IDEG,(MFT_SEQ)(ASDEG_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE
#undef EXPAND_TEMPLATE
#undef REPEAT_IDEG
#endif


/*
Compute quality function and derivative with respect to node ivar
- gdim is geometric dimension: also topological dimension !
- asdmsh is mesh as P1 or Pk
- asdmet is AsDeg::P1 or AsDeg::Pk: MetricField handles its degree
- ftype is arithmetic precision (debug): double, float8...
- ivar is the DoF, specified Bézier or Lagrange by: 
- dofbas is either FEBasis::Lagrange or FEBasis::Bezier -> whether the DoF is a 
  Lagrange node or Bézier control point 
- idifmet is either DifVar::None ("static" metric approximation) or DifVar::Phys
- quael is output quality and derivatives 
*/
template <class MFT, int gdim, typename ftype>
void d_quafun_tradet(Mesh<MFT> &msh, AsDeg asdmsh, AsDeg asdmet,
                     const int* ent2pol,
                     const double*__restrict__ bary, 
                     int power, 
                     int ivar, 
                     FEBasis dofbas, 
                     DifVar idifmet, 
                     ftype*__restrict__ tra_, 
                     ftype*__restrict__ dtra, 
                     ftype*__restrict__ htra, 
                     ftype*__restrict__ det_, 
                     ftype*__restrict__ ddet, 
                     ftype*__restrict__ hdet){


  static_assert(gdim == 2 || gdim == 3);
  METRIS_ASSERT(gdim == msh.idim);
  METRIS_ASSERT(power != 0);
  METRIS_ASSERT(msh.met.getSpace() == MetSpace::Log)
  // Differentiate or don't, but there is no barycentric derivative in this context 
  METRIS_ASSERT(idifmet == DifVar::None || idifmet == DifVar::Phys);
  if(idifmet != DifVar::None) METRIS_THROW_MSG(TODOExcept(), 
                           "Metric field derivative not implemented in quality")
  METRIS_ASSERT( !(ivar >= 0 && dofbas == FEBasis::Undefined) );
  if(dofbas == FEBasis::Bezier) METRIS_THROW_MSG(TODOExcept(), 
    "Ctrl pt dof not implemented -> do lag2bez derivatives of metric")

  constexpr int tdim  = gdim;
  constexpr int nnmet = (gdim*(gdim+1))/2;
  //constexpr int nhess = nnmet;

  ftype &tra = *tra_;
  ftype &det = *det_;

  // Get Jacobian matrix at xi  
  // Derivatives are not needed, we compute them ourselves, as they greatly simplify
  // see docs/quality/qualityiff.pdf
  double jmat[tdim*gdim],coopr[gdim];
  // Get Jacobian matrix at xi
  if(asdmsh == AsDeg::P1){
    if constexpr(tdim == 2){  
      eval2<gdim,1>(msh.coord,ent2pol,msh.getBasis(),
                       DifVar::Bary,DifVar::None,
                       bary,coopr,jmat,NULL);
    }else{
      eval3<gdim,1>(msh.coord,ent2pol,msh.getBasis(),
                       DifVar::Bary,DifVar::None,
                       bary,coopr,jmat,NULL);
    }
  }else{
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){

      if constexpr(tdim == 2){  
        eval2<gdim,ideg>(msh.coord,ent2pol,msh.getBasis(),
                         DifVar::Bary,DifVar::None,
                         bary,coopr,jmat,NULL);
      }else{
        eval3<gdim,ideg>(msh.coord,ent2pol,msh.getBasis(),
                         DifVar::Bary,DifVar::None,
                         bary,coopr,jmat,NULL);
      }

    }}CT_FOR1(ideg);
  }
  // Compute J_0^{-T} J_K^T M J_K J_0^{-1}
  // Note that J_K is stored transposed w.r.t. above ! -> jmat[i,j] = d_i F_j 
  // whereas (J_K)_{ij} = d_j F_i .. 

  // Get J_0^{-T} J_K^T 
  ftype invtJ0_tJ[tdim][gdim];
  matXmat<tdim,tdim,gdim>(Constants::invtJ_0[hana::type_c<ftype>][tdim],
                          jmat,invtJ0_tJ[0]);
  // Get metric interpolated at xi
  // The metric field class fetches geometric dimension from the mesh
  double met[nnmet],dmet[gdim][nnmet];
  msh.met.getMetFullinfo(asdmsh,idifmet,MetSpace::Exp,
                         ent2pol,tdim,bary,coopr,met,dmet[0]);

  // Compute the trace 
  tra = tra_matXsymXtmat<gdim, double, ftype, ftype>(met, invtJ0_tJ[0]);
  //tra = invtJ0_tJ_M_J_invJ0iag[0] + invtJ0_tJ_M_J_invJ0iag[1];
  //if constexpr (tdim == 3) tra += invtJ0_tJ_M_J_invJ0iag[2]; 
  // This is an actual exception that should never theoretically happen. 
  if(tra < 1.0e-16) 
    METRIS_THROW_MSG(GeomExcept(), "NEGATIVE J^TMJ trace "<<tra);



  ftype detM, det_invtJ0_tJ;
  if constexpr(tdim == gdim){
    // Error of 10^{-19} = 10^-14 relative compared to Matlab on wonky case
    det_invtJ0_tJ = detmat<gdim,ftype>(invtJ0_tJ[0]);
    // Error of 10^5 ... 
    // Matlab yields 7.346639223296765e+09 we get 7.346911345315383e+09 
    // Even in relative this is terrible>(invtJ0_tJ[0]); 
    //ftype tmp = detsym<gdim>(met); 
    // Also wrong, same error. Note Matlab gets 10^-9 final quality 
    // relative error ! Our determinant is terribly bad
    detM = detsym2<gdim>(met); 
    det  = det_invtJ0_tJ*det_invtJ0_tJ*detM; 
  }else{
    static_assert(tdim == 2);
    ftype J0tJtMJJ0[3];
    matXsymXtmat<2,3,double,ftype,ftype>(met,invtJ0_tJ[0],J0tJtMJJ0);
    det = detsym2<2>(J0tJtMJJ0);
  }

  if(abs(det) < 1.0e-16 && power > 0) 
     METRIS_THROW_MSG(GeomExcept(), "Singular J^TMJ det = "<<det);


  // This is used later on -> store it 
  //int dpowd = iipow<tdim>(tdim);
  //ftype trapowdm2 = irpow<tdim-2,ftype>(tra);//irpow.template operator()<tdim-2>(tra);
  //ftype trapowdm1 = trapowdm2*tra;
  //ftype trapowd   = trapowdm1*tra;


  if(ivar < 0) return;
  // See docs/quality/qualityiff.pdf for details 


  // Get the derivatives (d_k+1 - d_1) f
  double dfun[tdim];
  constexpr auto ordent = ORDELT(tdim);
  // multi-index of ivar: 
  int idx[tdim+1];
  for(int ii = 0; ii < tdim+1; ii++) idx[ii] = ordent[msh.curdeg][ivar][ii];

  // This is what we called \psi_k in the pdf document 
  if(asdmsh == AsDeg::P1){
    if(dofbas == FEBasis::Bezier){
      eval_bezierfunc<1,tdim>(idx,bary,1,dfun);
    }else{
      eval_lagrangefunc<1,tdim>(idx,bary,1,dfun);
    }
  }else{
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(ideg == msh.curdeg){
      if(dofbas == FEBasis::Bezier){
        eval_bezierfunc<ideg,tdim>(idx,bary,1,dfun);
      }else{
        eval_lagrangefunc<ideg,tdim>(idx,bary,1,dfun);
      }
    }}CT_FOR1(ideg);
  }

  // The derivative is much simplified, see pdf doc 
  ftype D_J_invJ0[tdim];

  // d_i(invtJ0_tJ)_{ij} = \sum_k \psi_k C_kj^T
  // In these notations, C_kj^T = Constants::invtJ_0[hana::type_c<ftype>][tdim][k][j]
  // BECAUSE THE JACOBIAN MATRICES ARE TRANSPOSED! 
  ftype sumDk2 = 0; // needed for second derivative of trace 
  for(int ii = 0; ii < tdim; ii++){
    D_J_invJ0[ii] = 0;
    for(int kk = 0; kk < tdim; kk++){
      D_J_invJ0[ii] 
        += dfun[kk]*Constants::invtJ_0[hana::type_c<ftype>][tdim][tdim*ii+kk];
    }
    sumDk2 += D_J_invJ0[ii]*D_J_invJ0[ii];
  }

  // Compute A^TM 
  ftype invtJ0_tJ_M[tdim][gdim];
  matXsym<gdim>(invtJ0_tJ[0],met,invtJ0_tJ_M[0]);



  // Compute trace derivatives
  // Given by d_i tr = 2 \sum_j (d_iA^T)_{ji} (A^TM)_{ji}
  //ftype dtra[gdim] = {0};//, htra[nhess] = {0};

  for(int ii = 0; ii < gdim; ii++){
    dtra[ii] = 0;
    for(int jj = 0; jj < gdim; jj++){
      dtra[ii] += 2*invtJ0_tJ_M[jj][ii]
                   *D_J_invJ0[jj];
    }
    //dtrapowd[ii] = tdim*dtra[ii]*trapowdm1;
  }
  if(htra != NULL){
    //ftype htra[nhess];
    for(int ii = 0; ii < gdim; ii++){
      for(int jj = ii; jj < gdim; jj++){
        //htra[sym2idx(ii,jj)] = 2.0*met[sym2idx(ii,jj)]*sumDk2;
        htra[sym2idx(ii,jj)] = 2.0*met[sym2idx(ii,jj)]*sumDk2;
        //htrapowd[sym2idx(ii,jj)] = trapowdm2*tdim*( htra[sym2idx(ii,jj)]*tra
        //                                          + (tdim-1)*dtra[ii]*dtra[jj]);
      }
    }
  }

  // Compute determinant derivatives 
  // Given by two terms:
  // d_i det(A^TMA) = det(A)det(M) ( T1 + T2 )
  // T1 = det(diA_(:1)   A_(:2)   A_(:3)) 
  //    + det(  A_(:1) diA_(:2)   A_(:3)) 
  //    + det(  A_(:1)   A_(:2) diA_(:3)) 
  // T2 = det(diA_(1:) + det(  A_(1:) + det(  A_(1:)                                     
  //            A_(2:)       diA_(2:)         A_(2:)                                 
  //            A_(3:))        A_(3:))      diA_(3:))
  //    = only one term, the one where di is in the i-th position. 
  // We'll need both invtJ0_tJ and its transpose:
  ftype J_invJ0[tdim][tdim];
  for(int ii = 0; ii < tdim; ii++){
    for(int jj = 0; jj < tdim; jj++){
      J_invJ0[jj][ii] = invtJ0_tJ[ii][jj];
    }
  }
  // Compute T2 first. 
  if constexpr (gdim == 2){
    ddet[0] = detvec2<ftype>(D_J_invJ0   ,  J_invJ0[1]);
    ddet[1] = detvec2<ftype>(  J_invJ0[0],D_J_invJ0   );
  }else{
    ddet[0] = detvec3<ftype>(D_J_invJ0   ,  J_invJ0[1],  J_invJ0[2]);
    ddet[1] = detvec3<ftype>(  J_invJ0[0],D_J_invJ0   ,  J_invJ0[2]);
    ddet[2] = detvec3<ftype>(  J_invJ0[0],  J_invJ0[1],D_J_invJ0   );
  }
  // Now T1. Normally the matrix is column-wise, and we eliminate (i,1), (i,2)...
  // Since subdetvec takes lines, we instead give it the lines of the transpose,
  // and so we eliminate (1,i), (2,i) and so on. 

  // At this stage, ddet stores grad of det(A). (if the full matrix is A^TMA). 
  // Now compute the second derivatives:
  if(hdet != NULL){
    for(int ii = 0; ii < gdim; ii++){
      for(int jj = ii; jj < gdim; jj++){
        hdet[sym2idx(ii,jj)] = 2.0*detM*ddet[ii]*ddet[jj];
      }
    }
  }
  for(int ii = 0; ii < gdim;ii++){
    ddet[ii] *= 2*detM*det_invtJ0_tJ;
  }

  return;
}



// While cumbersome, this replaces a bunch of manual instantiations, about to 
// be made worse the day we add tdimn as a template argument. 
#define EXPAND_TEMPLATE(r,SEQ) \
                  INSTANTIATE(BOOST_PP_SEQ_ELEM(0, SEQ),\
                              BOOST_PP_SEQ_ELEM(1, SEQ))

#define INSTANTIATE(MFT_VAL,FTYPE)\
template void d_quafun_tradet< MFT_VAL , 2, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int* ent2pol,\
                   const double*__restrict__ bary,\
                   int power,\
                   int ivar,\
                   FEBasis dofbas,\
                   DifVar idifmet,\
                   FTYPE*__restrict__ tra,\
                   FTYPE*__restrict__ dtra,\
                   FTYPE*__restrict__ htra,\
                   FTYPE*__restrict__ det,\
                   FTYPE*__restrict__ ddet,\
                   FTYPE*__restrict__ hdet);\
template void d_quafun_tradet< MFT_VAL , 3, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   AsDeg asdmsh, AsDeg asdmet,\
                   const int* ent2pol,\
                   const double*__restrict__ bary,\
                   int power,\
                   int ivar,\
                   FEBasis dofbas,\
                   DifVar idifmet,\
                   FTYPE*__restrict__ tra,\
                   FTYPE*__restrict__ dtra,\
                   FTYPE*__restrict__ htra,\
                   FTYPE*__restrict__ det,\
                   FTYPE*__restrict__ ddet,\
                   FTYPE*__restrict__ hdet); 
BOOST_PP_SEQ_FOR_EACH_PRODUCT(EXPAND_TEMPLATE,\
                              (MFT_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE
#undef EXPAND_TEMPLATE
#undef REPEAT_GDIM

#undef MFT_SEQ  


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

#define INSTANTIATE(z,gdim,ideg,MFT_VAL,ASDEG_VAL,FTYPE)\
template void d_quafun_tradet< MFT_VAL , 2+gdim, 1+ideg, ASDEG_VAL, FTYPE>\
                  (Mesh< MFT_VAL > &msh,\
                   const int* ent2pol,\
                   const double*__restrict__ bary,\
                   int power,\
                   int ivar,\
                   FEBasis dofbas,\
                   DifVar idifmet,\
                   FTYPE*__restrict__ tra,\
                   FTYPE*__restrict__ dtra,\
                   FTYPE*__restrict__ htra,\
                   FTYPE*__restrict__ det,\
                   FTYPE*__restrict__ ddet,\
                   FTYPE*__restrict__ hdet); 
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