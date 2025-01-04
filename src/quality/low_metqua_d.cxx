//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "low_metqua_d.hxx"
#include "quafun.hxx"

#include "../linalg/sym3idx.hxx"
#include "../linalg/det.hxx"
#include "../linalg/matprods.hxx"
#include "../Mesh/Mesh.hxx"
#include "../aux_utils.hxx"

#include "../aux_pp_inc.hxx"

namespace Metris{


// START DEBUG 

//#define DEBUG_MACRO(r,SEQ)  toto(SEQ )
//BOOST_PP_SEQ_FOR_EACH_PRODUCT(DEBUG_MACRO,(MFT_SEQ)(ASDEG_SEQ)(FTYPE_SEQ))



template <class MFT, int gdim, QuaFun iquaf, typename ftype>
ftype d_metqua(Mesh<MFT> &msh, AsDeg asdmsh, AsDeg asdmet,
               int ielem, int power, 
               int ivar, FEBasis dofbas, DifVar idifmet, 
               ftype*__restrict__ dquael, ftype*__restrict__ hquael, 
               int pnorm, double difto){
  constexpr int tdim = gdim;
  int* ent2poi = tdim == 2 ? msh.fac2poi[ielem] : msh.tet2poi[ielem];
  return d_metqua0<MFT,gdim,iquaf,ftype>(msh,asdmsh,asdmet,
                                         ent2poi,power,
                                         ivar,dofbas,idifmet,
                                         dquael,hquael,
                                         pnorm,difto);
}

// Hessian is optional (pass in NULL)
template <class MFT, int gdim, QuaFun iquaf, typename ftype>
ftype d_metqua0(Mesh<MFT> &msh, AsDeg asdmsh, AsDeg asdmet,
                const int* ent2poi, int power,
                int ivar, FEBasis dofbas, DifVar idifmet, 
                ftype*__restrict__ dquael, ftype*__restrict__ hquael,
                int pnorm, double difto){
  static_assert(gdim==2 || gdim==3);
  METRIS_ASSERT(pnorm > 0);
  constexpr int tdim = gdim;
  constexpr int nhess = (gdim*(gdim+1))/2;

  double bary[tdim+1];

  ftype qutet = 0; 


  constexpr auto d_quafun_xi 
    = get_d_quafun_xi<MFT,gdim,tdim,iquaf,ftype>();

  const int ideg = msh.curdeg;
  if(asdmet == AsDeg::Pk && ideg > 1){

    const int idegj = SMOO_DEGJ(ideg);
    const int nnodj = tdim == 2 ? facnpps[idegj] : tetnpps[idegj];

    const auto ordelt = ORDELT(tdim);

    if(hquael != NULL && ivar >= 0)
      for(int ii = 0; ii < nhess; ii++) hquael[ii] = 0;
    if(ivar >= 0)
      for(int ii = 0; ii < gdim; ii++) dquael[ii] = 0;

    ftype qua0, dqua0[gdim], hqua0[nhess];
    for(int iquad = 0; iquad < nnodj; iquad++){

      for(int ii = 0; ii < tdim + 1; ii++)
        bary[ii] = ordelt[idegj][iquad][ii]/((double) (idegj));

      if(hquael == NULL){
        qua0 = d_quafun_xi(msh,asdmsh,asdmet,
                           ent2poi,bary,power,
                           ivar,dofbas,idifmet,
                           dqua0,NULL);
      }else{
        qua0 = d_quafun_xi(msh,asdmsh,asdmet,
                           ent2poi,bary,power,
                           ivar,dofbas,idifmet,
                           dqua0,hqua0);
      }
      ftype powm1 = pow(abs(qua0 - difto),pnorm-1);
      qutet += abs(qua0 - difto)*abs(powm1)/nnodj;

      if(ivar < 0) continue;

      // If p%2 = 1 and f = qua0 - difto < 0, 
      // notice d|f|^p = d(-f)^p = p(-f')(-f)^(p-1) = p(-1)^p f' f^(p-1)
      // but since p%2 = 1, (-1)^p = -1, so d|f|^p = -pf'f^(p-1) = -d(f)^p
      // note f^(p-1) = |f|^(p-1) as p-1 = 0 [2] 
      int sg = 1;
      if(qua0 - difto < 0) sg = -1;
      for(int ii = 0; ii < gdim; ii++){
        dquael[ii] += sg*pnorm*dqua0[ii]*powm1/nnodj;
      }

      if(hquael == NULL) continue;

      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++){
          hquael[sym3idx(ii,jj)] += sg*pnorm*hqua0[sym3idx(ii,jj)]*powm1;
        }
      }

      if(pnorm < 2) continue;

      ftype powm2 = pow(abs(qua0 - difto),pnorm-2);
      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++){
          // Only the gradient gets a sign difference
          hquael[sym3idx(ii,jj)] += 
            pnorm*(pnorm - 1)*dqua0[ii]*dqua0[jj]*powm2;
        }
      }

    }// for iquad 
  }else{
    for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0/(tdim  + 1);
    qutet = d_quafun_xi(msh,asdmsh,asdmet,
                        ent2poi,bary,power,
                        ivar,dofbas,idifmet,
                        dquael,hquael);
    ftype powm1 = pow(abs(qutet - difto),pnorm-1);
    qutet = powm1*abs(qutet - difto); 

    if(ivar < 0) return qutet;

    int sg = 1;
    if(qutet - difto < 0) sg = -1;
    for(int ii = 0; ii < gdim; ii++){
      dquael[ii] = sg*pnorm*dquael[ii]*powm1;
    }

    if(hquael == NULL) return qutet;

    for(int ii = 0; ii < gdim; ii++){
      for(int jj = ii; jj < gdim; jj++){
        hquael[sym3idx(ii,jj)] = sg*pnorm*hquael[sym3idx(ii,jj)]*powm1;
      }
    }
    if(pnorm >= 2){
      ftype powm2 = pow(abs(qutet - difto),pnorm-2);
      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++){
          // Only the gradient gets a sign difference
          hquael[sym3idx(ii,jj)] += pnorm*(pnorm - 1)*dquael[ii]*dquael[jj]*powm2;
        }
      }
    }
  }
  return qutet;
}

// While cumbersome, this replaces a bunch of manual instantiations, about to 
// be made worse the day we add tdimn as a template argument. 
#define EXPAND_TEMPLATE(z,gdim,SEQ) \
                  INSTANTIATE(gdim,BOOST_PP_SEQ_ELEM(0, SEQ),\
                                   BOOST_PP_SEQ_ELEM(1, SEQ),\
                                   BOOST_PP_SEQ_ELEM(2, SEQ))
#define REPEAT_GDIM(r,SEQ) BOOST_PP_REPEAT(2,EXPAND_TEMPLATE,SEQ)
#define MFT_SEQ (MetricFieldFE)(MetricFieldAnalytical)
#define QUAFUN_SEQ (QuaFun::Distortion)(QuaFun::Unit)
#define INSTANTIATE(gdim,MFT_VAL,QUAFUN,FTYPE)\
template FTYPE d_metqua< MFT_VAL , 2+gdim, QUAFUN,FTYPE>\
                  (Mesh< MFT_VAL > &msh, AsDeg asdmsh, AsDeg asdmet,\
                   int ielem, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, FTYPE*__restrict__ hquael, \
                   int pnorm,\
                   double difto);
BOOST_PP_SEQ_FOR_EACH_PRODUCT(REPEAT_GDIM,\
                              (MFT_SEQ)(QUAFUN_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE

#define INSTANTIATE(gdim,MFT_VAL,QUAFUN,FTYPE)\
template FTYPE d_metqua0< MFT_VAL , 2+gdim, QUAFUN, FTYPE>\
                  (Mesh< MFT_VAL > &msh, AsDeg asdmsh, AsDeg asdmet,\
                   const int* ent2poi, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, FTYPE*__restrict__ hquael, \
                   int pnorm, double difto);
BOOST_PP_SEQ_FOR_EACH_PRODUCT(REPEAT_GDIM,\
                              (MFT_SEQ)(QUAFUN_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE
#undef EXPAND_TEMPLATE
#undef REPEAT_GDIM






/*
------------------------- Full derivatives (all control points)
*/


/*
D_quafun_distortion : computes derivatives (up to second) wrt all control points in the 
element. 
- dquael is ftype* size N = gdim*entnpps[ideg] (entnpps -> gdim)
- hquael is ftype* size sym NxN. Ordering: 
  hquael(symidx(gdim*ip1 + i,gdim*ip2 + j)) = d_{ij}^{ip1,ip2} 
Refer to d_quafun_distortion for additional comments within the routine. 
*/
template <class MFT, int gdim, int ideg, AsDeg asdmet, typename ftype>
ftype D_quafun_distortion(Mesh<MFT> &msh, 
                  const int* ent2poi,
                  const double*__restrict__ bary, 
                  int power, 
                  FEBasis dofbas, 
                  DifVar idifmet, 
                  ftype*__restrict__ dquael, 
                  ftype*__restrict__ hquael){


  static_assert(gdim == 2 || gdim == 3);
  METRIS_ASSERT(gdim == msh.idim);
  METRIS_ASSERT(power != 0);
  METRIS_ASSERT(msh.met.getSpace() == MetSpace::Log)
  // Differentiate or don't, but there is no barycentric derivative in this context 
  METRIS_ASSERT(idifmet == DifVar::None || idifmet == DifVar::Phys);
  if(idifmet != DifVar::None) METRIS_THROW_MSG(TODOExcept(), 
                           "Metric field derivative not implemented in quality")
  //METRIS_ASSERT( !(idiff != DifVar::None && dofbas == FEBasis::Undefined) );
  if(dofbas == FEBasis::Bezier) METRIS_THROW_MSG(TODOExcept(), 
    "Ctrl pt dof not implemented -> do lag2bez derivatives of metric")
  //METRIS_ASSERT( !(idiff != DifVar::None && dquael == NULL) );

  constexpr int tdim  = gdim;
  constexpr int nnmet = (gdim*(gdim+1))/2;
  constexpr int nhess = nnmet;


  ftype quael;


  // Get Jacobian matrix at xi  
  // Derivatives are not needed, we compute them ourselves, as they greatly simplify
  // see docs/quality/qualityiff.pdf
  double jmat[tdim*gdim],coopr[gdim];
  if constexpr(tdim == 2){
    eval2<gdim,ideg>(msh.coord,ent2poi,msh.getBasis(),DifVar::Bary,DifVar::None,
                                                          bary,coopr,jmat,NULL);
  }else{
    eval3<gdim,ideg>(msh.coord,ent2poi,msh.getBasis(),DifVar::Bary,DifVar::None,
                                                          bary,coopr,jmat,NULL);
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
  msh.met.getMetFullinfo(asdmet,idifmet,MetSpace::Exp,
                         ent2poi,tdim,bary,coopr,met,dmet[0]);

  //printf("Debug metric = \n");
  //dblAr1(nnmet,met).print();

  // Compute the trace 
  ftype tra;
  //ftype invtJ0_tJ_M_J_invJ0iag[tdim];
  tra = tra_matXsymXtmat<gdim, double, ftype, ftype>(met, invtJ0_tJ[0]);
  //tra = invtJ0_tJ_M_J_invJ0iag[0] + invtJ0_tJ_M_J_invJ0iag[1];
  //if constexpr (tdim == 3) tra += invtJ0_tJ_M_J_invJ0iag[2]; 
  // This is an actual exception that should never theoretically happen. 
  if(tra < 1.0e-16) METRIS_THROW_MSG(GeomExcept(), "NEGATIVE J^TMJ trace "<<tra);



  ftype det ;
  ftype detM, det_invtJ0_tJ;
  if constexpr(tdim == gdim){
    det_invtJ0_tJ = detmat<gdim,ftype>(invtJ0_tJ[0]); // Error of 10^{-19} = 10^-14 relative compared to Matlab on wonky case
    //ftype tmp = detsym<gdim>(met); // Error of 10^5 ... // Matlab yields 7.346639223296765e+09 we get 7.346911345315383e+09 // Even in relative this is terrible
    detM = detsym2<gdim>(met); // Also wrong, same error. Note Matlab gets 10^-9 final quality relative error ! Our determinant is terribly bad
    det = det_invtJ0_tJ*det_invtJ0_tJ*detM; 
  }else{
    static_assert(tdim == 2);
    ftype J0tJtMJJ0[3];
    matXsymXtmat<2,3,double,ftype,ftype>(met,invtJ0_tJ[0],J0tJtMJJ0);
    det = detsym2<2>(J0tJtMJJ0);
  }

  if(abs(det) < 1.0e-16 && power > 0) 
     METRIS_THROW_MSG(GeomExcept(), "Singular J^TMJ det = "<<det);

  //constexpr auto irpow = [&]<int n>(ftype x) -> ftype{
  //  if constexpr(std::is_same<ftype,double>::value) return idpow<n>(x);
  //  if constexpr(std::is_same<ftype,float4>::value) return id4pow<n>(x);
  //  if constexpr(std::is_same<ftype,float8>::value) return id8pow<n>(x);
  //};

  // This is used later on -> store it 
  int dpowd = iipow<tdim>(tdim);
  ftype trapowdm2 = irpow<tdim-2,ftype>(tra);//irpow.template operator()<tdim-2>(tra);
  ftype trapowdm1 = trapowdm2*tra;
  ftype trapowd   = trapowdm1*tra;

  if(power > 0){
    quael = trapowd/(det*dpowd);
  }else{
    quael = (det*dpowd)/trapowd;
  }
  ftype quae1 = quael; // for derivatives
  quael = pow(quael, abs(power));



  //// From here, we compute derivatives. 
  //if(idiff == DifVar::None) return quael;
  //// See docs/quality/qualityiff.pdf for details 

  constexpr int nnode = tdim == 1 ? edgnpps[ideg] : 
                        tdim == 2 ? facnpps[ideg] : tetnpps[ideg] ;
  constexpr auto ordent = ORDELT(tdim);


  // Compute A^TM 
  ftype invtJ0_tJ_M[tdim][gdim];
  matXsym<gdim>(invtJ0_tJ[0],met,invtJ0_tJ_M[0]);
  ftype J_invJ0[tdim][tdim];
  for(int ii = 0; ii < tdim; ii++){
    for(int jj = 0; jj < tdim; jj++){
      J_invJ0[jj][ii] = invtJ0_tJ[ii][jj];
    }
  }

  // Compute D_J_invJ0 which depends on each node, 
  // then dtra and ddet
  ftype D_J_invJ0[nnode][tdim], ddet[nnode][gdim], ddetA[nnode][gdim], dtra[nnode][gdim];
  for(int inode = 0; inode < nnode; inode++){
    // Get the derivatives (d_k+1 - d_1) f
    double dfun[tdim];
    // multi-index of inode: 
    int idx[tdim+1];
    for(int ii = 0; ii < tdim+1; ii++) idx[ii] = ordent[ideg][inode][ii];

    // This is what we called \psi_k in the pdf document 
    if(dofbas == FEBasis::Bezier){
      eval_bezierfunc<ideg,tdim>(idx,bary,1,dfun);
    }else{
      eval_lagrangefunc<ideg,tdim>(idx,bary,1,dfun);
    }
    // The derivative is much simplified, see pdf doc 
    // d_i(invtJ0_tJ)_{ij} = \sum_k \psi_k C_kj^T
    // In these notations, C_kj^T = Constants::invtJ_0[hana::type_c<ftype>][tdim][k][j]
    // BECAUSE THE JACOBIAN MATRICES ARE TRANSPOSED! 
    for(int ii = 0; ii < tdim; ii++){
      D_J_invJ0[inode][ii] = 0;
      for(int kk = 0; kk < tdim; kk++){
        D_J_invJ0[inode][ii] += dfun[kk]*Constants::invtJ_0[hana::type_c<ftype>][tdim][tdim*ii+kk];
      }
    }

    for(int ii = 0; ii < gdim; ii++){
      dtra[inode][ii] = 0;
      for(int jj = 0; jj < gdim; jj++){
        dtra[inode][ii] += 2*invtJ0_tJ_M[jj][ii]
                            *D_J_invJ0[inode][jj];
      }
    }

    if constexpr (gdim == 2){
      ddetA[inode][0] = detvec2<ftype>(D_J_invJ0[inode],  J_invJ0[1]);
      ddetA[inode][1] = detvec2<ftype>(  J_invJ0[0]    ,D_J_invJ0[inode]);
    }else{
      ddetA[inode][0] = detvec3<ftype>(D_J_invJ0[inode],  J_invJ0[1]    ,  J_invJ0[2]    );
      ddetA[inode][1] = detvec3<ftype>(  J_invJ0[0]    ,D_J_invJ0[inode],  J_invJ0[2]    );
      ddetA[inode][2] = detvec3<ftype>(  J_invJ0[0]    ,  J_invJ0[1]    ,D_J_invJ0[inode]);
    }
    for(int ii = 0; ii < gdim;ii++){
      ddet[inode][ii] = 2*ddetA[inode][ii]*detM*det_invtJ0_tJ;
    }

    if(power > 0){
      for(int ii = 0; ii < gdim; ii++){
        dquael[inode*gdim + ii] = trapowdm1*( tdim*dtra[inode][ii]*det
                               - tra*ddet[inode][ii])
                    /(det*det*dpowd); 
        //dquael[inode*gdim + ii] = trapowdm1*( tdim*dtra[inode][ii]*det_invtJ0_tJ
        //                                    - 2.0*ddetA[inode][ii]*tra )/(det_invtJ0_tJ*det*dpowd); 
      }
    }else{
      for(int ii = 0; ii < gdim; ii++){
        dquael[inode*gdim+ii] = dpowd*(     2* tra   *ddetA[inode][ii]
                                      - tdim*dtra[inode][ii]*det_invtJ0_tJ)*detM*det_invtJ0_tJ/(trapowd*tra); 
        //printf("grad = %f Contrib  =",(double)dquael[inode*gdim+ii]);
        //std::cout<<ddetA[inode][ii]<<" "<<dtra[inode][ii]<<" "<<trapowd<<" "<<tra<<" \n";
      }
    }
    if(abs(power) != 1){
      // Q^(abs(power)) , quael is Q^p, quae1 is just Q
      for(int ii = 0; ii < gdim; ii++){
        dquael[inode*gdim+ii] = abs(power)*dquael[inode*gdim+ii]*quael/quae1;
      }
    }

  }


  // Sized for two points local Hessian 
  constexpr int nhess2 = ( (2*gdim)*(2*gdim + 1) )/ 2;
  // Compute Hessian
  //const int dpowdpowp = pow(dpowd,power); 
  for(int inode = 0; inode < nnode; inode++){
    ftype sumDk2 = 0; // needed for second derivative of trace 
    for(int ii = 0; ii < tdim; ii++){
      sumDk2 += D_J_invJ0[inode][ii]*D_J_invJ0[inode][ii];
    }

    ftype htra[nhess]; // These are temp 
    for(int ii = 0; ii < gdim; ii++){
      for(int jj = ii; jj < gdim; jj++){
        htra[sym3idx(ii,jj)] = 2.0*met[sym3idx(ii,jj)]*sumDk2;
      }
    }

    ftype hdet[nhess] = {0};//, hdet2[nhess2] = {0};
    for(int ii = 0; ii < gdim; ii++){
      for(int jj = ii; jj < gdim; jj++){
        hdet[sym3idx(ii,jj)] = 2.0*detM*ddetA[inode][ii]*ddetA[inode][jj];
      }
    }

    if(power > 0){
      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++){ 
          hquael[sym3idx(gdim*inode+ii,gdim*inode+jj)] = 
            ( tdim*trapowdm2*det*det*( tra*htra[sym3idx(ii,jj)]
                                     + (tdim - 1)*dtra[inode][ii]*dtra[inode][jj] )
            - trapowdm1*det*( tra*hdet[sym3idx(ii,jj)]
                            + tdim*dtra[inode][ii]*ddet[inode][jj]
                            + tdim*dtra[inode][jj]*ddet[inode][ii])
            + 2.0*trapowd*ddet[inode][ii]*ddet[inode][jj]
            )/(det*det*det);
          hquael[sym3idx(gdim*inode+ii,gdim*inode+jj)] /= dpowd;
        }
      }
    }else{
      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++){
          hquael[sym3idx(gdim*inode+ii,gdim*inode+jj)]  = 
            (-(tdim+1)*dtra[inode][jj]*(ddet[inode][ii]*tra - tdim*det*dtra[inode][ii])
            + tra*( hdet[sym3idx(ii,jj)]*tra + ddet[inode][ii]*dtra[inode][jj] 
                  - tdim*ddet[inode][jj]*dtra[inode][ii]
                  - tdim*det*htra[sym3idx(ii,jj)]) 
            )/(trapowd*tra*tra);
          hquael[sym3idx(gdim*inode+ii,gdim*inode+jj)]  *= dpowd;
        }
      }
    }// if power > 0
  }


  for(int inod1 = 0; inod1 < nnode; inod1++){
    for(int inod2 = inod1+1; inod2 < nnode; inod2++){
      ftype sumDkDk = 0;
      for(int ii = 0; ii < tdim; ii++){
        sumDkDk += D_J_invJ0[inod1][ii]*D_J_invJ0[inod2][ii];
      }

      ftype htra[nhess]; // These are temp 
      for(int ii = 0; ii < gdim; ii++){
        for(int jj = ii; jj < gdim; jj++){
          htra[sym3idx(ii,jj)] = 2.0*met[sym3idx(ii,jj)]*sumDkDk;
        }
      }

      // some extradiagonal terms are added due to d^2 detA != 0 now 
      ftype hdet[nhess2] = {0};//, hdet2[nhess2] = {0};
      const ftype detA = det_invtJ0_tJ; // Just for clarity 
      if constexpr (gdim == 2){
        // The crossed-point terms (also crossed coordinate)
        hdet[sym3idx(0*gdim+0,1*gdim+1)] = detvec2<ftype>(D_J_invJ0[inod1], D_J_invJ0[inod2])*detA;
        hdet[sym3idx(0*gdim+1,1*gdim+0)] = - hdet[sym3idx(0*gdim+0,1*gdim+1)];//detvec2<ftype>(D_J_invJ0[inod1], D_J_invJ0[inod2])*detA;
      }else{
        hdet[sym3idx(0*gdim+0,1*gdim+1)] = 
          detvec3<ftype>(D_J_invJ0[inod1], D_J_invJ0[inod2],  J_invJ0[2]      )*detA;
        hdet[sym3idx(0*gdim+0,1*gdim+2)] = 
          detvec3<ftype>(D_J_invJ0[inod1],   J_invJ0[1]    ,  D_J_invJ0[inod2])*detA;
        hdet[sym3idx(0*gdim+1,1*gdim+2)] = 
          detvec3<ftype>(  J_invJ0[0]    , D_J_invJ0[inod1],  D_J_invJ0[inod2])*detA;

        hdet[sym3idx(1*gdim+0,0*gdim+1)] = - hdet[sym3idx(0*gdim+0,1*gdim+1)];
        hdet[sym3idx(1*gdim+0,0*gdim+2)] = - hdet[sym3idx(0*gdim+0,1*gdim+2)];
        hdet[sym3idx(1*gdim+1,0*gdim+2)] = - hdet[sym3idx(0*gdim+1,1*gdim+2)];
      }

      for(int ii = 0; ii < gdim; ii++){
        for(int jj = 0; jj < gdim; jj++){
          hdet[sym3idx(0*gdim+ii,1*gdim+jj)] += ddetA[inod1][ii]*ddetA[inod2][jj];
          hdet[sym3idx(0*gdim+ii,1*gdim+jj)] *= 2.0*detM;
        }
      }
       
      if(power > 0){
        for(int ii = 0; ii < gdim; ii++){
          for(int jj = 0; jj < gdim; jj++){
            hquael[sym3idx(gdim*inod1+ii,gdim*inod2+jj)] = 
              ( tdim*trapowdm2*det*det*( tra*htra[sym3idx(ii,jj)]
                                       + (tdim - 1)*dtra[inod1][ii]*dtra[inod2][jj] )
              - trapowdm1*det*( tra*hdet[sym3idx(0*gdim+ii,1*gdim+jj)]
                              + tdim*dtra[inod1][ii]*ddet[inod2][jj]
                              + tdim*dtra[inod2][jj]*ddet[inod1][ii])
              + 2.0*trapowd*ddet[inod1][ii]*ddet[inod2][jj]
              )/(det*det*det);
            hquael[sym3idx(gdim*inod1+ii,gdim*inod2+jj)] /= dpowd;
          }
        }
      }else{
        for(int ii = 0; ii < gdim; ii++){
          for(int jj = 0; jj < gdim; jj++){
            hquael[sym3idx(gdim*inod1+ii,gdim*inod2+jj)] = 
              (-(tdim+1)*dtra[inod2][jj]*(ddet[inod1][ii]*tra - tdim*det*dtra[inod1][ii])
              + tra*( hdet[sym3idx(0*gdim+ii,1*gdim+jj)]*tra 
                    + ddet[inod1][ii]*dtra[inod2][jj] 
                    - tdim*ddet[inod2][jj]*dtra[inod1][ii]
                    - tdim*det*htra[sym3idx(ii,jj)]) 
              )/(trapowd*tra*tra);
            hquael[sym3idx(gdim*inod1+ii,gdim*inod2+jj)] *= dpowd;
          }
        }
      }
      //if(power > 0){
      //  for(int ii = 0; ii < gdim; ii++){
      //    for(int jj = 0; jj < gdim; jj++){
      //      hquael[sym3idx(gdim*inod1+ii,gdim*inod2+jj)] =    
      //        (- hquael[sym3idx(gdim*inod1+ii,gdim*inod2+jj)]*quael
      //        + 2.0*dquael[gdim*inod1+ii]*dquael[gdim*inod2+jj]) /(quael*quael*quael);
      //    }
      //  }
      //}

    }
  }



  return quael;
}


#define INSTANTIATE(gdim,ideg,MFT_VAL,ASDEG_VAL,FTYPE)\
template FTYPE D_quafun_distortion< MFT_VAL , 2+gdim, 1+ideg, ASDEG_VAL, FTYPE>\
                  (Mesh< MFT_VAL > &msh, \
                  const int* ent2poi,\
                  const double*__restrict__ bary, \
                  int power, /*DifVar idiff, \*/\
                  FEBasis dofbas, \
                  DifVar idifmet, \
                  FTYPE*__restrict__ dquael, \
                  FTYPE*__restrict__ hquael);

#define EXPAND_TEMPLATE(z,gdim,SEQ) \
                  INSTANTIATE(gdim,BOOST_PP_SEQ_ELEM(0, SEQ),\
                                   BOOST_PP_SEQ_ELEM(1, SEQ),\
                                   BOOST_PP_SEQ_ELEM(2, SEQ),\
                                   BOOST_PP_SEQ_ELEM(3, SEQ))
                  
#define REPEAT_GDIM(z,n,SEQ) BOOST_PP_REPEAT(2,EXPAND_TEMPLATE,(n)SEQ)
#define REPEAT_IDEG(r,SEQ)   BOOST_PP_REPEAT(METRIS_MAX_DEG,REPEAT_GDIM,SEQ)

#define ASDEG_SEQ (AsDeg::P1)(AsDeg::Pk)
BOOST_PP_SEQ_FOR_EACH_PRODUCT(REPEAT_IDEG,\
                             (MFT_SEQ)(ASDEG_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE




template <class MFT, int gdim, int ideg, AsDeg asdmet, typename ftype>
ftype D_metqua(Mesh<MFT> &msh, int ielem, int power, 
               FEBasis dofbas, DifVar idifmet, 
               ftype*__restrict__ dquael, ftype*__restrict__ hquael, 
               int pnorm, double difto){
  constexpr int tdim = gdim;
  int* ent2poi = tdim == 2 ? msh.fac2poi[ielem] : msh.tet2poi[ielem];
  return D_metqua<MFT,gdim,ideg,asdmet,ftype>(msh,ent2poi,power,
                                                dofbas,idifmet,
                                                dquael,hquael,
                                                pnorm,difto);
}


template <class MFT, int gdim, int ideg, AsDeg asdmet, typename ftype>
ftype D_metqua(Mesh<MFT> &msh, const int* ent2poi, int power,
              FEBasis dofbas, DifVar idifmet, 
              ftype*__restrict__ dquael, 
              ftype*__restrict__ hquael,
              int pnorm, double difto){
  static_assert(gdim==2 || gdim==3);
  METRIS_ASSERT(pnorm > 0);
  constexpr int tdim = gdim;

  double bary[tdim+1];

  ftype qutet = 0; 

  for(int ii = 0; ii < gdim; ii++) dquael[ii] = 0;

  if constexpr(asdmet == AsDeg::Pk && ideg > 1){

    constexpr int idegj = SMOO_DEGJ(ideg);
    constexpr int nnodj = tdim == 2 ? facnpps[idegj] : tetnpps[idegj];

    constexpr int nnode = tdim == 1 ? edgnpps[ideg] : 
                          tdim == 2 ? facnpps[ideg] : tetnpps[ideg] ;

    constexpr auto ordelt = ORDELT(tdim);

    if(dquael != NULL){
      for(int ii = 0; ii < gdim*nnode; ii++) dquael[ii] = 0;
    }
    if(hquael != NULL){
      for(int ii = 0; ii < gdim*nnode; ii++){
        for(int jj = ii; jj < gdim*nnode; jj++){
          hquael[sym3idx(ii,jj)] = 0;
        }
      }
    }
    constexpr int nhess = (nnode*gdim*(nnode*gdim + 1))/2;
    ftype dqua0[nnode*gdim], hqua0[nhess];
    for(int iquad = 0; iquad < nnodj; iquad++){
      for(int ii = 0; ii < tdim + 1; ii++){
        bary[ii] = ordelt[idegj][iquad][ii]/((double) (idegj));
      }
      ftype qua0 = D_quafun_distortion<MFT,gdim,ideg,asdmet,ftype>(msh,ent2poi,bary,power,
                                                       dofbas,idifmet,
                                                       dqua0,hqua0);
      ftype powm1 = pow(qua0 - difto,pnorm-1);
      qutet += (qua0 - difto)*powm1/nnodj;

      int sg = 1;
      if(qua0 - difto < 0 && pnorm % 2 == 1) sg = -1;


      if(hquael != NULL){

        for(int ii = 0; ii < gdim*nnode; ii++){
          for(int jj = ii; jj < gdim*nnode; jj++){
            hquael[sym3idx(ii,jj)] += 
             sg*pnorm*hqua0[sym3idx(ii,jj)]*powm1;
          }
        }

        if(pnorm >= 2){
          ftype powm2 = pow(qua0 - difto,pnorm-2);
          for(int ii = 0; ii < gdim*nnode; ii++){
            for(int jj = ii; jj < gdim*nnode; jj++){
              hquael[sym3idx(ii,jj)] +=
                sg*pnorm*(pnorm - 1) 
               *dqua0[ii]*dqua0[jj]*powm2;
            }
          }
        } // if(pnorm >= 2)
        
      } // if(hquael != NULL)



      for(int ii = 0; ii < gdim*nnode; ii++){
        dquael[ii] += sg*pnorm*dqua0[ii]*powm1/nnodj;
      }


    }
  }else{
    constexpr int nnode = tdim + 1;

    for(int ii = 0; ii < tdim + 1; ii++) bary[ii] = 1.0/(tdim  + 1);
    ftype qua0 = D_quafun_distortion<MFT,gdim,1,asdmet,ftype>(msh,ent2poi,bary,power,
                                                      dofbas,idifmet,
                                                      dquael,hquael);

    int sg = 1;
    if(qua0 - difto < 0 && pnorm % 2 == 1) sg = -1;


    ftype powm1 = pow(qua0 - difto,pnorm-1);
    qutet = sg*powm1*(qua0 - difto); 


    if(hquael != NULL){
      for(int ii = 0; ii < gdim*nnode; ii++){
        for(int jj = ii; jj < gdim*nnode; jj++){
          hquael[sym3idx(ii,jj)] = 
            sg*pnorm*hquael[sym3idx(ii,jj)]*powm1;
        }
      }
      if(pnorm >= 2){
        ftype powm2 = pow(qua0 - difto,pnorm-2);
        for(int ii = 0; ii < gdim*nnode; ii++){
          for(int jj = ii; jj < gdim*nnode; jj++){
            hquael[sym3idx(ii,jj)] +=
                  sg*pnorm*(pnorm - 1) 
                 *dquael[ii]*dquael[jj]*powm2;
          }
        }
      }// endif pnorm 
    }

    for(int ii = 0; ii < gdim*nnode; ii++){
      dquael[ii] = sg*pnorm*dquael[ii]*powm1;
    }

  }
  return qutet;
}

#define EXPAND_TEMPLATE(z,gdim,SEQ) \
                  INSTANTIATE(gdim,BOOST_PP_SEQ_ELEM(0, SEQ),\
                                   BOOST_PP_SEQ_ELEM(1, SEQ),\
                                   BOOST_PP_SEQ_ELEM(2, SEQ),\
                                   BOOST_PP_SEQ_ELEM(3, SEQ))
#define REPEAT_GDIM(z,n,SEQ) BOOST_PP_REPEAT(2,EXPAND_TEMPLATE,(n)SEQ)
#define REPEAT_IDEG(r,SEQ)   BOOST_PP_REPEAT(METRIS_MAX_DEG,REPEAT_GDIM,SEQ)
                  
#define INSTANTIATE(gdim,ideg,MFT_VAL,ASDEG_VAL,FTYPE)\
template FTYPE D_metqua< MFT_VAL , 2+gdim, 1+ideg, ASDEG_VAL, FTYPE>\
                  (Mesh< MFT_VAL > &msh, \
                   int ielem, \
                   int power, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, FTYPE*__restrict__ hquael, \
                   int pnorm,\
                   double difto); 
BOOST_PP_SEQ_FOR_EACH_PRODUCT(REPEAT_IDEG,(MFT_SEQ)(ASDEG_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE

#define INSTANTIATE(gdim,ideg,MFT_VAL,ASDEG_VAL,FTYPE)\
template FTYPE D_metqua< MFT_VAL , 2+gdim, 1+ideg, ASDEG_VAL, FTYPE>\
                  (Mesh< MFT_VAL > &msh, \
                   const int* ent2poi, \
                   int power, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   FTYPE*__restrict__ dquael, FTYPE*__restrict__ hquael, \
                   int pnorm, double difto); 
BOOST_PP_SEQ_FOR_EACH_PRODUCT(REPEAT_IDEG,(MFT_SEQ)(ASDEG_SEQ)(QUA_FTYPE_SEQ))
#undef INSTANTIATE







#undef EXPAND_TEMPLATE
#undef REPEAT_GDIM
#undef REPEAT_IDEG
#undef MFT_SEQ // note these two could go into headers 
#undef ASDEG_SEQ 



#if 0
#define BOOST_PP_LOCAL_MACRO(n)\
template double d_quafun_distortion<MetricFieldAnalytical, 2, n, AsDeg::P1, double>\
                  (Mesh<MetricFieldAnalytical> &msh, \
                   int ielem, const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   double*__restrict__ dquael);\
/* switch met  */ \
template double d_quafun_distortion<MetricFieldFE        , 2, n, AsDeg::P1, double>\
                  (Mesh<MetricFieldFE        > &msh, \
                   int ielem, const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   double*__restrict__ dquael);\
/* switch dim + met */ \
template double d_quafun_distortion<MetricFieldAnalytical, 3, n, AsDeg::P1, double>\
                  (Mesh<MetricFieldAnalytical> &msh, \
                   int ielem, const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   double*__restrict__ dquael);\
/* switch met  */ \
template double d_quafun_distortion<MetricFieldFE        , 3, n, AsDeg::P1, double>\
                  (Mesh<MetricFieldFE        > &msh, \
                   int ielem, const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   double*__restrict__ dquael);\
/* switch AsDeg + dim + met  */ \
template double d_quafun_distortion<MetricFieldAnalytical, 2, n, AsDeg::Pk, double>\
                  (Mesh<MetricFieldAnalytical> &msh, \
                   int ielem, const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   double*__restrict__ dquael);\
/* switch met  */ \
template double d_quafun_distortion<MetricFieldFE        , 2, n, AsDeg::Pk, double>\
                  (Mesh<MetricFieldFE        > &msh, \
                   int ielem, const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   double*__restrict__ dquael);\
/* switch dim + met */ \
template double d_quafun_distortion<MetricFieldAnalytical, 3, n, AsDeg::Pk, double>\
                  (Mesh<MetricFieldAnalytical> &msh, \
                   int ielem, const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   double*__restrict__ dquael);\
/* switch met  */ \
template double d_quafun_distortion<MetricFieldFE        , 3, n, AsDeg::Pk, double>\
                  (Mesh<MetricFieldFE        > &msh, \
                   int ielem, const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   double*__restrict__ dquael);
#define BOOST_PP_LOCAL_LIMITS     (1, __MAXEG_JACOBIAN__)
#include BOOST_PP_LOCAL_ITERATE()
#endif

#if 0
#define INSTANTIATE(z,gdim,ideg) \
template double d_quafun_distortion<MetricFieldAnalytical, 2+gdim, 1+ideg, AsDeg::P1, double>\
                  (Mesh<MetricFieldAnalytical> &msh, \
                   int ielem, const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   double*__restrict__ dquael);\
/* switch met  */ \
template double d_quafun_distortion<MetricFieldFE        , 2+gdim, 1+ideg, AsDeg::P1, double>\
                  (Mesh<MetricFieldFE        > &msh, \
                   int ielem, const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   double*__restrict__ dquael);\
/* switch AsDeg + dim + met  */ \
template double d_quafun_distortion<MetricFieldAnalytical, 2+gdim, 1+ideg, AsDeg::Pk, double>\
                  (Mesh<MetricFieldAnalytical> &msh, \
                   int ielem, const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   double*__restrict__ dquael);\
/* switch met  */ \
template double d_quafun_distortion<MetricFieldFE        , 2+gdim, 1+ideg, AsDeg::Pk, double>\
                  (Mesh<MetricFieldFE        > &msh, \
                   int ielem, const double*__restrict__ bary, \
                   int power, \
                   int ivar, \
                   FEBasis dofbas, \
                   DifVar idifmet, \
                   double*__restrict__ dquael);
#endif

////BOOST_PP_REPEAT(2,INSTANTIATE,~)
//#define REPEAT_GDIM(z,n,data) BOOST_PP_REPEAT(2,INSTANTIATE,n)
//BOOST_PP_REPEAT(METRIS_MAX_DEG,REPEAT_GDIM,~)


//#define LOOP_MFT(data) BOOST_PP_SEQ_FOR_EACH(REPEAT_IDEG,data,MFT_SEQ)
////LOOP_MFT(~)
//
//#define LOOP_MFT_r(r,data,tupelem) BOOST_PP_SEQ_FOR_EACH_R(r,REPEAT_IDEG,tupelem,MFT_SEQ)
//BOOST_PP_SEQ_FOR_EACH(LOOP_MFT_r,~,MFT_SEQ)
////BOOST_PP_SEQ_FOR_EACH(REPEAT_IDEG,~,MFT_SEQ)
//
//
//#define DEBUG_TAKES_PRODUCT(r,product) f(BOOST_PP_SEQ_ELEM(0,product), BOOST_PP_SEQ_ELEM(1,product))
//BOOST_PP_SEQ_FOR_EACH_PRODUCT(DEBUG_TAKES_PRODUCT,(MFT_SEQ)(ASDEG_SEQ))
//
//
//
//#define EXPAND_TEMPLATE2(z,gdim,SEQ) g(gdim,SEQ) h(gdim,BOOST_PP_SEQ_ELEM(0,SEQ),\
//                                                        BOOST_PP_SEQ_ELEM(1,SEQ),\
//                                                        BOOST_PP_SEQ_ELEM(2,SEQ))
//
////h(gdim,BOOST_PP_SEQ_ELEM(1,SEQ),\
//                                                 //       BOOST_PP_SEQ_ELEM(0,BOOST_PP_SEQ_ELEM(0,SEQ)),\
//                                                 //       BOOST_PP_SEQ_ELEM(1,BOOST_PP_SEQ_ELEM(0,SEQ))
//#define REPEAT_GDIM3(r,n,SEQ) BOOST_PP_REPEAT(2,EXPAND_TEMPLATE2, (n)SEQ)
//#define REPEAT_IDEG2(r,SEQ) BOOST_PP_REPEAT(METRIS_MAX_DEG,REPEAT_GDIM3,SEQ)
//
//
//BOOST_PP_SEQ_FOR_EACH_PRODUCT(REPEAT_IDEG2,(MFT_SEQ)(ASDEG_SEQ))



//#define BOOST_PP_LOCAL_MACRO(n) BOOST_PP_REPEAT(2,MACRO2,n)
//#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG_JACOBIAN)
//#include BOOST_PP_LOCAL_ITERATE()






} // End namespace
