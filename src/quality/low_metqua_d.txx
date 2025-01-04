//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __LOW_VOLQUA_D__
#define __LOW_VOLQUA_D__

#include "../aux_exceptions.hxx"

#include "../low_geo.hxx"
#include "../ho_quadrature.hxx"
//#include "../low_eval_d.hxx"
#include "../linalg/explogmet.hxx"

#include "../../SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"


namespace Metris{

//#include "../low_localization.hxx"

//#include "../codegen_ccoef.hxx"
//#include "../low_ccoef.hxx"

#if 0
template <class MFT, int gdim, int tdim, int mshdeg, 
          AsDeg metdeg, typename ftype>
void d_quafun_distortion(Mesh<MFT> &msh, int ielem, int power, 
                 const double*__restrict__ bary, int ivar, DifVar metderiv, ftype* qutet){

  static_assert(gdim == 2 || gdim == 3);
  static_assert(tdim <= gdim); 

  METRIS_ASSERT(gdim == msh.idim);

  METRIS_ASSERT(power != 0);

  if(msh.met.getSpace() != MetSpace::Log) METRIS_THROW_MSG(WArgExcept(),
      "## SET MESH METRIC TO LOG BEFORE CALLING metqua2_xi");

  if(msh.getBasis() != FEBasis::Bezier ) METRIS_THROW_MSG(WArgExcept(), 
      "## METQUA DIFF ONLY AVAILABLE IN BEZIER");

  // Differentiate or don't, but there is no barycentric derivative in this context 
  METRIS_ASSERT(metderiv == DifVar::None || metderiv == DifVar::Phys);

  constexpr int nnmet = (gdim*(gdim+1))/2;

  double jmat[tdim*gdim],met[nnmet],dmet[gdim][nnmet],coopr[gdim];

  int* ent2pol = tdim == 2 ? msh.fac2poi[ielem] : msh.tet2poi[ielem];

  // Get Jacobian matrix at xi
  if constexpr(tdim == 2){
    eval2<gdim,mshdeg>(msh.coord,ent2pol,msh.getBasis(),DifVar::Bary,DifVar::None,
                                                          bary,coopr,jmat,NULL);
  }else{
    eval3<gdim,mshdeg>(msh.coord,ent2pol,msh.getBasis(),DifVar::Bary,DifVar::None,
                                                          bary,coopr,jmat,NULL);
  }

  auto ordent = 
  double bezfun = eval_bezierfunc<ideg,tdim>()

  // Get metric interpolated at xi
  // The metric field class fetches geometric dimension from the mesh
  msh.met.template getMetFullinfo<metdeg,metdeg>(metderiv,MetSpace::Exp,
                                              ent2pol,tdim,bary,coopr,met,dmet[0]);

  // Compute J_0^{-T} J_K^T M J_K J_0^{-1}
  // Starting with J_K J_0^{-1}
  // Note that J_K is stored transposed w.r.t. above ! -> jmat[i,j] = d_i F_j 
  // whereas (J_K)_{ij} = d_j F_i .. 


  // Get J_0^{-T} J_K^T -> the Jacobian depends on the variable. 
  SANS::SurrealS<gdim, ftype> invtJ_0tJ_K[tdim*gdim];
  matXmat<tdim,tdim,gdim>(Constants::invtJ_0[hana::type_c<ftype>][tdim],
                                                              jmatS,invtJ_0tJ_K);


  ftype J0tJtMJJ0_diag[tdim];
  matXsymXtmat_diag<gdim, tdim, double, ftype, ftype>(met, invtJ_0tJ_K, J0tJtMJJ0_diag);
  ftype tra = J0tJtMJJ0_diag[0] + J0tJtMJJ0_diag[1];
  if(tdim == 3) tra += J0tJtMJJ0_diag[2];
  
  // This is an actual exception that should never theoretically happen. 
  if(tra < 1.0e-16) METRIS_THROW_MSG(GeomExcept(), "NEGATIVE J^TMJ trace "<<tra);


  ftype det ;
  if constexpr(tdim == gdim){
    ftype det1 = detmat<gdim>(invtJ_0tJ_K); // Error of 10^{-19} = 10^-14 relative compared to Matlab on wonky case
    //ftype tmp = detsym<gdim>(met); // Error of 10^5 ... // Matlab yields 7.346639223296765e+09 we get 7.346911345315383e+09 // Even in relative this is terrible
    ftype tmp = detsym2<gdim>(met); // Also wrong, same error. Note Matlab gets 10^-9 final quality relative error ! Our determinant is terribly bad
    det = det1*det1*tmp; 
  }else{
    static_assert(tdim == 2);
    ftype J0tJtMJJ0[3];
    matXsymXtmat<2,3,double,ftype,ftype>(met,invtJ_0tJ_K,J0tJtMJJ0);
    det = detsym2<2>(J0tJtMJJ0);
  }

  if(abs(det) < 1.0e-16 && power > 0) 
     METRIS_THROW_MSG(GeomExcept(), "Singular J^TMJ det = "<<det);

  if constexpr (tdim == 2){
    if(power > 0){
      *qutet = pow(tra*tra/det/4,power);
    }else{
      *qutet = pow(4*det/(tra*tra),-power);
    }
  }else{
    if(power > 0){
      *qutet = pow(tra*tra*tra/det/27,power);
    }else{
      *qutet = pow(27*det/(tra*tra*tra),-power);
    }
  }

}
#endif

template <class MFT, int gdim, int mshdeg, AsDeg metdeg, typename ftype>
d_quafun_distortion<MFT,gdim,mshdeg,metdeg,ivar,ftype>::d_quafun_distortion(Mesh<MFT> &msh, 
            int ielem, int power, 
            const double*__restrict__ bary, DifVar metderiv, 
            int ivar, 
            SANS::SurrealS<gdim,ftype>* qutet){
  static_assert(gdim == 2 || gdim == 3);
  METRIS_ASSERT(gdim == msh.idim);
  METRIS_ASSERT(power != 0);
  if(msh.met.getSpace() != MetSpace::Log) METRIS_THROW_MSG(WArgExcept(),
      "## SET MESH METRIC TO LOG BEFORE CALLING metqua2_xi");
  if(msh.getBasis() != FEBasis::Bezier ) METRIS_THROW_MSG(WArgExcept(), 
      "## METQUA DIFF ONLY AVAILABLE IN BEZIER");
  // Differentiate or don't, but there is no barycentric derivative in this context 
  METRIS_ASSERT(metderiv == DifVar::None || metderiv == DifVar::Phys);
  if(metderiv != DifVar::None) METRIS_THROW_MSG(TODOExcept(), 
                           "Metric field derivative not implemented in quality")

  constexpr int tdim  = gdim;
  constexpr int nnmet = (gdim*(gdim+1))/2;


  double jmat[tdim*gdim],met[nnmet],dmet[gdim][nnmet],coopr[gdim];

  int* ent2pol = tdim == 2 ? msh.fac2poi[ielem] : msh.tet2poi[ielem];

  // Get Jacobian matrix at xi
  if constexpr(tdim == 2){
    eval2<gdim,mshdeg>(msh.coord,ent2pol,msh.getBasis(),DifVar::Bary,DifVar::None,
                                                          bary,coopr,jmat,NULL);
  }else{
    eval3<gdim,mshdeg>(msh.coord,ent2pol,msh.getBasis(),DifVar::Bary,DifVar::None,
                                                          bary,coopr,jmat,NULL);
  }


  // Get metric interpolated at xi
  // The metric field class fetches geometric dimension from the mesh
  msh.met.getMetFullinfo(metdeg,metderiv,MetSpace::Exp,
                         ent2pol,tdim,bary,coopr,met,dmet[0]);

  // Compute J_0^{-T} J_K^T M J_K J_0^{-1}
  // Starting with J_K J_0^{-1}
  // Note that J_K is stored transposed w.r.t. above ! -> jmat[i,j] = d_i F_j 
  // whereas (J_K)_{ij} = d_j F_i .. 


  // Get J_0^{-T} J_K^T -> the Jacobian depends on the variable. 
  SANS::SurrealS<gdim, ftype> invtJ_0tJ_K[tdim*gdim];
  //ftype invtJ_0tJ_K[tdim*gdim];
  //matXmat<tdim,tdim,gdim>(Constants::invtJ_0[hana::type_c<ftype>][tdim],
  //                                                            jmatS,invtJ_0tJ_K);
  //HERE;

  METRIS_THROW_MSG(TODOExcept(),"ORDELT ARRAY REMOVED -> FIX d_quafun_distortion here");
  //constexpr auto ordent = ORDELT_ARRAY(tdim); 
  //constexpr std::array<int,tdim+1> idxP = ordent[mshdeg][ivar];
  //constexpr std::array<int,tdim+1> idx1 = ordent[mshdeg][ivar] - ordent[1][0];
  constexpr std::array<int,tdim+1> idxP = 0;
  constexpr std::array<int,tdim+1> idx1 = 0;
  //double dbezf[tdim];
  //constexpr int idx[tdim+1] = [&](){
  //  if constexpr(tdim == 1){
  //    return {ordedg.s[mshdeg][ivar][0],ordedg.s[mshdeg][ivar][1]};
  //  }else if(tdim ==2){
  //    return {ordfac.s[mshdeg][ivar][0],
  //            ordfac.s[mshdeg][ivar][1],
  //            ordfac.s[mshdeg][ivar][2]};
  //  }else{
  //    return {ordtet.s[mshdeg][ivar][0],
  //            ordtet.s[mshdeg][ivar][1],
  //            ordtet.s[mshdeg][ivar][2],
  //            ordtet.s[mshdeg][ivar][3]};
  //  }
  //}();
  //for(int ii = 0; ii < tdim+1; ii++){
  //  idx[ii] = ordent[mshdeg][ivar]; 
  //}
  double dbez0 = eval_bezierfunc<mshdeg,tdim>(ordent[mshdeg],bary,0,NULL);
  for(int ii = 0; ii < tdim; ii++){
    // Get multi-index idx - e_{ii+1}
    constexpr int idx2[tdim + 2] = ordent[mshdeg][ivar] - ordent[1][ii+1];
    double bezfun = eval_bezierfunc<mshdeg,tdim>(ordent[mshdeg],bary,0,NULL);

  }


  ftype J0tJtMJJ0_diag[tdim];
  matXsymXtmat_diag<gdim, tdim, double, ftype, ftype>(met, invtJ_0tJ_K, J0tJtMJJ0_diag);
  ftype tra = J0tJtMJJ0_diag[0] + J0tJtMJJ0_diag[1];
  if(tdim == 3) tra += J0tJtMJJ0_diag[2];
  
  // This is an actual exception that should never theoretically happen. 
  if(tra < 1.0e-16) METRIS_THROW_MSG(GeomExcept(), "NEGATIVE J^TMJ trace "<<tra);


  ftype det ;
  if constexpr(tdim == gdim){
    ftype det1 = detmat<gdim>(invtJ_0tJ_K); // Error of 10^{-19} = 10^-14 relative compared to Matlab on wonky case
    //ftype tmp = detsym<gdim>(met); // Error of 10^5 ... // Matlab yields 7.346639223296765e+09 we get 7.346911345315383e+09 // Even in relative this is terrible
    ftype tmp = detsym2<gdim>(met); // Also wrong, same error. Note Matlab gets 10^-9 final quality relative error ! Our determinant is terribly bad
    det = det1*det1*tmp; 
  }else{
    static_assert(tdim == 2);
    ftype J0tJtMJJ0[3];
    matXsymXtmat<2,3,double,ftype,ftype>(met,invtJ_0tJ_K,J0tJtMJJ0);
    det = detsym2<2>(J0tJtMJJ0);
  }

  if(abs(det) < 1.0e-16 && power > 0) 
     METRIS_THROW_MSG(GeomExcept(), "Singular J^TMJ det = "<<det);

  if constexpr (tdim == 2){
    if(power > 0){
      *qutet = pow(tra*tra/det/4,power);
    }else{
      *qutet = pow(4*det/(tra*tra),-power);
    }
  }else{
    if(power > 0){
      *qutet = pow(tra*tra*tra/det/27,power);
    }else{
      *qutet = pow(27*det/(tra*tra*tra),-power);
    }
  }

}



/*----------------
OLD IMPLEMENTATIONS
-----------------*/
///*
//  Differentiation of quality/low_volqua.hxx functions. 
//  New arguments:
//   - gdim (template param, optional): number of variables, typically 3 (physical) or 4 (barycentric). 
//   Default 3 and no dpoivar argument for physical coordinates. Could be 2 for surface-bound points. 
//   - dmetvar (mandatory argument): derivatives of ivar-th metric w.r.t. the gdims. 
//   - dpoivar (optional argument, mandatory if gdim != 3): derivatives of ivar-th control point/node 
//   w.r.t. the gdims. 
//*/
//template <int gdim, int tdim, int mshdeg, AsDeg metdeg, int ivar, int gdim, typename ftype>
//metqua3_xi_d<gdim,tdim,mshdeg,metdeg,ivar,gdim,ftype>
//::metqua3_xi_d(const Mesh &msh, int ielem, int power, 
//               const double*__restrict__ bary, 
//               const double*__restrict__ dmetvar, 
//               const double*__restrict__ dpoivar, 
//               SANS::SurrealS<gdim,ftype>* qutet){
//
//  const double exptol = 1.0e-12;
//
//  bool debug = true;
//
//  if(msh.ianamet <= 0){
//    if(dmetvar == NULL){
//      printf("## metqua3_xi_d: dpoivar == NULL admissible but not dmetvar == NULL\n");
//      exit(1);
//    }
//    if(msh.ilogmet != 1){
//      printf("## SET MESH METRIC TO LOG BEFORE CALLING metqua3_xi_d \n");
//      exit(1);
//    }
//  }
//  double jmat[9],djmat[3][9],lmet[6],met[6],dlmet[3][6],dmet[3][6],coopr[3],ddum[9];
//
//  // Get Jacobian matrix and derivatives at xi.
//  // eval3_d handles whether dpoivar == NULL or not itself. 
//  eval3_d<3,ideg,ilag,1,0,ivar,gdim>(msh.coord,msh.tet2poi[ielem],bary,coopr,jmat,NULL,
//                                     ddum,djmat[0],NULL,dpoivar);
//
//  SANS::SurrealS<gdim,ftype> jmat_S[9];
//  // Could we instead consider returning a Surreal directly from eval3?
//  for(int ii=0; ii<9; ii++){
//    jmat_S[ii].value()  =  (ftype) jmat[ii];
//    for(int idof = 0; idof < gdim; idof++)
//      jmat_S[ii].deriv(idof) = (ftype) djmat[idof][ii];
//  }
//
//  /*
//  Interpolate metric and derivatives at xi. 
//  Derivatives depend on dmetvar, the derivatives of the metric at ivar w.r.t. to the variables. 
//  */
//
//
//  if(msh.ianamet <= 0){
//    if(msh.ilagmet == 1){
//      eval3_d<6,ideg,1,0,0,ivar,gdim>(msh.met,msh.tet2poi[ielem],bary,lmet,NULL,NULL,dlmet[0],NULL,NULL,dmetvar);
//    }else{
//      eval3_d<6,ideg,0,0,0,ivar,gdim>(msh.met,msh.tet2poi[ielem],bary,lmet,NULL,NULL,dlmet[0],NULL,NULL,dmetvar);
//    }
//    getexpmet_cpy_d(lmet,dlmet[0],met,dmet[0]);
//  }else{
//    msh.anamet(NULL,coopr,1,met,dmet[0]);
//  }
//
//  SANS::SurrealS<gdim,ftype> met_S[6];
//  for(int ii=0; ii<6; ii++){
//    met_S[ii].value() = (ftype)met[ii];
//    for(int idof = 0; idof < gdim; idof++)
//      met_S[ii].deriv(idof) = (ftype)dmet[idof][ii];
//  }
//
//  // Note, our jmat is transposed:
//  // jmat[3*i+j] = d_i F_j
//
//  //printf("Debug ilag %d ideg %d jmat :\n",ilag,ideg);
//  //dblAr2(3,3,jmat).print();
//
//
//  // Compute J_0^{-T} J_K^T M J_K J_0^{-1}
//  // Starting with J_K J_0^{-1}
//  // Note that J_K is stored transposed w.r.t. above 
//  // The below J_0^{-1} is also transposed. 
//  ftype invtJ_0[3][3] = {
//      -0.57735026918962562,-0.57735026918962562 , 0                 ,
//      1                   ,-1                   , 0                 ,
//      -0.40824829046386302,-0.4082482904638630  , 1.2247448713915889};
//
//
//  // Get J_0^{-T} J_K^T
//  SANS::SurrealS<gdim,ftype> invtJ_0tJ_K[9];
//  matXmat<3                         ,
//                              ftype ,
//          SANS::SurrealS<gdim,ftype>,
//          SANS::SurrealS<gdim,ftype>>(invtJ_0[0],jmat_S,invtJ_0tJ_K);
//
//  // 2. Get whole product
//  //SANS::SurrealS<gdim,double> J0tJtMJJ0[6];
//  //matXsymXtmat<SANS::SurrealS<gdim,double>,
//  //     SANS::SurrealS<gdim,double>,
//  //     SANS::SurrealS<gdim,double>>(met_S, invtJ_0tJ_K, J0tJtMJJ0);
//  //SANS::SurrealS<gdim,double> tra = J0tJtMJJ0[0]
//  //                                + J0tJtMJJ0[2]
//  //                                + J0tJtMJJ0[5];
//  SANS::SurrealS<gdim,ftype> J0tJtMJJ0_diag[3];
//  matXsymXtmat_diag<3,SANS::SurrealS<gdim,ftype>,
//            SANS::SurrealS<gdim,ftype>,
//            SANS::SurrealS<gdim,ftype>>(met_S, invtJ_0tJ_K, J0tJtMJJ0_diag);
//
//  SANS::SurrealS<gdim,ftype> tra = J0tJtMJJ0_diag[0]
//                                 + J0tJtMJJ0_diag[1]
//                                 + J0tJtMJJ0_diag[2];
//
//  SANS::SurrealS<gdim,ftype> det1 
//   = detmat<3>(invtJ_0tJ_K); // <SANS::SurrealS<gdim,ftype>>
//  SANS::SurrealS<gdim,ftype> det 
//  = det1*det1*detsym<3>(met_S); // <SANS::SurrealS<gdim,ftype>>
//  //SANS::SurrealS<gdim,double> det 
//  //  = detsym<3,SANS::SurrealS<gdim,double>>(J0tJtMJJ0);
//
//  if(tra.value() < 1.0e-16) METRIS_THROW_MSG(GeomExcept(),
//    "NEGATIVE J^TMJ trace "<<tra.value());
//    
//  if(abs(det.value()) < 1.0e-16 && power > 0){
//    METRIS_THROW_MSG(GeomExcept(),
//    "Singular J^TMJ det = "<<det.value()<<" met det "<<detsym<3>(met)
//    <<" prod det "<<det1<<" met = "<<
//    met[0]<<" "<<met[1]<<" "<<met[2]<<" "<<
//    met[3]<<" "<<met[4]<<" "<<met[5]<<"\n"<<
//    "jmat "<<jmat[0]<<" "<<jmat[1]<<" "<<jmat[2]<<
//             jmat[3]<<" "<<jmat[4]<<" "<<jmat[5]<<
//             jmat[6]<<" "<<jmat[7]<<" "<<jmat[8]<<
//        "\n Use negative power to allow (almost-)singular elements.\n"<<power);
//  } 
//
//  if(power > 0){
//    *qutet = pow(tra*tra*tra/det/27,power);
//  }else{
//    *qutet = pow(27*det/(tra*tra*tra),-power);
//  }
//
//}
//
//
//
//
//template <int ideg, int ilag, int ivar, int gdim, typename ftype>
//metqua3_d<ideg,ilag,ivar,gdim,ftype>
//::metqua3_d(const Mesh & msh, int ielem, int power,
//            const double* __restrict__ dmetvar, 
//            const double* __restrict__  dpoivar, 
//            SANS::SurrealS<gdim,ftype> *qutet){
//
//  if constexpr(ideg > 1){
//    double bary[4];
//    *qutet = 0; // This takes care of derivatives
//  
//  
//    SANS::SurrealS<gdim,ftype> qua0;
//    int nquad = tetnpps[ideg - 1];
//    for(int iquad = 0; iquad < nquad; iquad++){
//      bary[0] = ordtet.s[ideg-1][iquad][0]/((double) (ideg - 1));
//      bary[1] = ordtet.s[ideg-1][iquad][1]/((double) (ideg - 1));
//      bary[2] = ordtet.s[ideg-1][iquad][2]/((double) (ideg - 1));
//      bary[3] = ordtet.s[ideg-1][iquad][3]/((double) (ideg - 1));
//      metqua3_xi_d<ideg,ilag,ivar,gdim,ftype>(
//         msh,ielem,power,bary,dmetvar,dpoivar,&qua0);
//      (*qutet) += qua0 / nquad;
//    }
//
//    //for(int iquad = 0; iquad < nquad_Keaton6; iquad++){
//    //  bary[0] = tuvquad_Keaton6[iquad][0];
//    //  bary[1] = tuvquad_Keaton6[iquad][1];
//    //  bary[2] = tuvquad_Keaton6[iquad][2];
//    //  bary[3] = 1 - bary[0] - bary[1] - bary[2];
//    //  metqua3_xi_d<ideg,ilag,ivar,gdim>(msh,ielem,bary,dmetvar,dpoivar,&qua0);
//    //  (*qutet) += qua0 * wtquad_Keaton6[iquad];
//    //}
//
//
//  }else{
//    double bary[4] = {0.25,0.25,0.25,0.25};
//    metqua3_xi_d<1,ilag,ivar,gdim,ftype>(msh,ielem,power,bary,dmetvar,dpoivar,qutet);
//  }
//}
//
//
//
//template <int ideg, int ilag, int gdim>
//metqua3_shell_d<ideg,ilag,gdim>
//::metqua3_shell_d(const Mesh &msh, int ipoin, int iele0, int power,
//                  int mshell, 
//                  int* __restrict__ nshell, int* __restrict__ lshell, 
//                  const double* __restrict__ dmetvar, 
//                  const double* __restrict__ dpoivar, 
//                  SANS::SurrealS<gdim,double>*  qushe){
//
//  int ierro;
//  if(*nshell <= 0){
//    int iopen;
//    int iver = getvertet<ideg>(iele0, msh.tet2poi, ipoin);
//    if(iver < 4 || iver >= 4 + 6*(edgnpps[ideg]-2))
//      METRIS_THROW_MSG(TopoExcept(),"VERTEX NOT ON EDGE")
//  
//    int ied = (iver - 4) / (edgnpps[ideg] - 2);
//  
//    int ipoi1 = msh.tet2poi(iele0,lnoed3[ied][0]);
//    int ipoi2 = msh.tet2poi(iele0,lnoed3[ied][1]);
//    shell3(msh,ipoi1,ipoi2,iele0,mshell,nshell,lshell,&iopen);
//    if(iopen >= 0) METRIS_THROW_MSG(TopoExcept(),"SHELL IS OPEN IN OPTIM")
//  }
//  
//  constexpr int nrfld = tetnpps[ideg];
//  auto nrfld_c = hana::int_c<nrfld>;
//
//  *qushe = 0.0;
//  for(int ishell = 0; ishell < *nshell; ishell++){
//    int ielem = lshell[ishell];
//    assert(ielem >= 0 && ielem < msh.nelem);
//
//    bool ifnd = false;
//    hana::while_(hana::less.than(nrfld_c), 4_c, [&](auto ivar_c){
//    constexpr int ivar = ivar_c;
//      if(msh.tet2poi(ielem,ivar) == ipoin){
//        ifnd = true;
//
//        SANS::SurrealS<3,double> qutet = 0;
//        metqua3_d<ideg,ilag,ivar,gdim>(msh,ielem,power,dmetvar,dpoivar,&qutet);
//        if(qutet < 0) METRIS_THROW_MSG(GeomExcept(),
//          "NEGATIVE ELEMENT QUALITY = "<<qutet<<" ielem = "<<ielem)
//    
//        (*qushe) += qutet;
//      }
//    return ivar_c+1_c;});
//    assert(ifnd == true);
//  }
//}






} // end namespace




#endif
