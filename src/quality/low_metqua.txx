//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_LOW_METQUA_D__
#define __METRIS_LOW_METQUA_D__

#include "../aux_exceptions.hxx"

#include "../low_geo.hxx"
#include "../ho_quadrature.hxx"
#include "../low_eval_d.hxx"
#include "../linalg/explogmet.hxx"

#include "../../SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"




//#include "../low_localization.hxx"

//#include "../codegen_ccoef.hxx"
//#include "../low_ccoef.hxx"

/*
  Differentiation of quality/low_metqua.hxx functions. 
  New arguments:
   - nvar (template param, optional): number of variables, typically 3 (physical) or 4 (barycentric). 
   Default 3 and no dpoivar argument for physical coordinates. Could be 2 for surface-bound points. 
   - dmetvar (mandatory argument): derivatives of ivar-th metric w.r.t. the nvars. 
   - dpoivar (optional argument, mandatory if nvar != 3): derivatives of ivar-th control point/node 
   w.r.t. the nvars. 
   gdim is the physical dimension of the mesh
*/
template <class MFT, int ideg, int gdim, int ivar, int nvar, typename ftype>
d_quafun_distortion<MFT,ideg,idim,ivar,nvar,ftype>
::d_quafun_distortion(const Mesh<MFT> &msh, const int* ent2poi, int power, 
              const double* bary, 
              const double* __restrict__ dmetvar, 
              const double* __restrict__ dpoivar, 
              SANS::SurrealS<nvar,ftype>* qutet){

	METRIS_ENFORCE_MSG(msh.met.getSpace() == MetSpace::Log, "Set metric to log before d_quafun_distortion call");

  if constexpr (std::is_same<MFT,MetricFieldAnalytical>::value) 
                                                 METRIS_ASSERT(dmetvar != NULL);

  constexpr int nnmet = (gdim*(gdim+1))/2;

  double jmat[gdim*gdim],
         djmat[gdim][gdim*gdim],
         lmet[nnmet],
         met[nnmet],
         dlmet(gdim,nnmet),
         dmet[gdim][nnmet],
         coopr[gdim],
         ddum[gdim*gdim];

  // Get Jacobian matrix and derivatives at xi.
  // eval3_d handles whether dpoivar == NULL or not itself. 
  if constexpr(gdim == 2){
   eval2_d<2,ideg,ivar,nvar>(msh.coord,ent2poi,msh.getBasis(),DifVar::Bary,0,bary,coopr,jmat,NULL,
                             ddum,djmat[0],NULL,dpoivar);
  }else{
   eval3_d<3,ideg,ivar,nvar>(msh.coord,ent2poi,msh.getBasis(),DifVar::Bary,0,bary,coopr,jmat,NULL,
                             ddum,djmat[0],NULL,dpoivar);
  }

  SANS::SurrealS<nvar,ftype> jmat_S[gdim*gdim];
  // Could we instead consider returning a Surreal directly from eval3?
  for(int ii=0; ii<gdim*gdim ; ii++){
    jmat_S[ii].value()  =  (ftype) jmat[ii];
    for(int idof = 0; idof < nvar; idof++)
      jmat_S[ii].deriv(idof) = (ftype) djmat[idof][ii];
  }

  /*
  Interpolate metric and derivatives at xi. 
  Derivatives depend on dmetvar, the derivatives of the metric at ivar w.r.t. to the variables. 
  */

  msh.met.getMetFullInfo_d()

  if(msh.ianamet <= 0){
    if(msh.ilagmet == 1){
      eval3_d<6,ideg,1,0,0,ivar,nvar>(msh.met,ent2poi,bary,lmet,NULL,NULL,dlmet[0],NULL,NULL,dmetvar);
    }else{
      eval3_d<6,ideg,0,0,0,ivar,nvar>(msh.met,ent2poi,bary,lmet,NULL,NULL,dlmet[0],NULL,NULL,dmetvar);
    }
    getexpmet_cpy_d(lmet,dlmet[0],met,dmet[0]);
  }else{
    msh.anamet(NULL,coopr,1,met,dmet[0]);
  }

  SANS::SurrealS<nvar,ftype> met_S[nnmet];
  for(int ii=0; ii<nnmet ; ii++){
    met_S[ii].value() = (ftype)met[ii];
    for(int idof = 0; idof < nvar; idof++)
      met_S[ii].deriv(idof) = (ftype)dmet[idof][ii];
  }

  // Note, our jmat is transposed:
  // jmat[3*i+j] = d_i F_j

  //printf("Debug ilag %d ideg %d jmat :\n",ilag,ideg);
  //dblAr2(3,3,jmat).print();


  // Compute J_0^{-T} J_K^T M J_K J_0^{-1}
  // Starting with J_K J_0^{-1}
  // Note that J_K is stored transposed w.r.t. above 
  // The below J_0^{-1} is also transposed. 

  //ftype invtJ_0[gdim][gdim] = [&]() -> std::initializer_list<ftype> {
  //  if constexpr(gdim == 2){
  //    return { 1                  , 0                 ,
  //            -0.577350269189626 , 1.154700538379252};
  //  }else{
  //    return {
  //    -0.57735026918962562,-0.57735026918962562 , 0                 ,
  //    1                   ,-1                   , 0                 ,
  //    -0.40824829046386302,-0.4082482904638630  , 1.2247448713915889};
  //  } 
  //}();

  // Get J_0^{-T} J_K^T
  SANS::SurrealS<nvar,ftype> invtJ_0tJ_K[gdim*gdim];
  matXmat<gdim                      ,
                              ftype ,
          SANS::SurrealS<nvar,ftype>,
          SANS::SurrealS<nvar,ftype>>(Constants::invtJ_0[hana::type_c<ftype>][gdim],jmat_S,invtJ_0tJ_K);

  // 2. Get whole product
  //SANS::SurrealS<nvar,double> J0tJtMJJ0[6];
  //matXsymXtmat<SANS::SurrealS<nvar,double>,
  //     SANS::SurrealS<nvar,double>,
  //     SANS::SurrealS<nvar,double>>(met_S, invtJ_0tJ_K, J0tJtMJJ0);
  //SANS::SurrealS<nvar,double> tra = J0tJtMJJ0[0]
  //                                + J0tJtMJJ0[2]
  //                                + J0tJtMJJ0[5];
  SANS::SurrealS<nvar,ftype> J0tJtMJJ0_diag[gdim];
  matXsymXtmat_diag<gdim,SANS::SurrealS<nvar,ftype>,
             				   SANS::SurrealS<nvar,ftype>,
             				   SANS::SurrealS<nvar,ftype>>(met_S, invtJ_0tJ_K, J0tJtMJJ0_diag);

  SANS::SurrealS<nvar,ftype> tra = J0tJtMJJ0_diag[0]
                                 + J0tJtMJJ0_diag[1];
  if(gdim == 3) tra += J0tJtMJJ0_diag[2]; 


  SANS::SurrealS<nvar,ftype> det1 = detmat<gdim>(invtJ_0tJ_K); // <SANS::SurrealS<nvar,ftype>>
  SANS::SurrealS<nvar,ftype> det = det1*det1*detsym<gdim>(met_S); // <SANS::SurrealS<nvar,ftype>>
  //SANS::SurrealS<nvar,double> det 
  //  = detsym<3,SANS::SurrealS<nvar,double>>(J0tJtMJJ0);

  if(tra.value() < 1.0e-16) METRIS_THROW_MSG(GeomExcept(),
    "NEGATIVE J^TMJ trace "<<tra.value());
    
	if(abs(det.value()) < 1.0e-16 && power > 0) 
     METRIS_THROW_MSG(GeomExcept(), "Singular J^TMJ det = "<<det);

  if(gdim == 2){
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




template <int ideg, int ilag, int ivar, int nvar, typename ftype>
metqua_d<ideg,ilag,ivar,nvar,ftype>
::metqua_d(const Mesh & msh, int ielem, int power,
            const double* __restrict__ dmetvar, 
            const double* __restrict__  dpoivar, 
            SANS::SurrealS<nvar,ftype> *qutet){

  //int* ent2poi = tdim == 2 ? msh.fac2poi[ielem] : msh.tet2poi[ielem];
  int *ent2poi = msh.tet2poi[ielem];
  if constexpr(ideg > 1){
    double bary[4];
    *qutet = 0; // This takes care of derivatives
  
  
    SANS::SurrealS<nvar,ftype> qua0;
    int nquad = tetnpps[ideg - 1];
    for(int iquad = 0; iquad < nquad; iquad++){
      bary[0] = ordtet.s[ideg-1][iquad][0]/((double) (ideg - 1));
      bary[1] = ordtet.s[ideg-1][iquad][1]/((double) (ideg - 1));
      bary[2] = ordtet.s[ideg-1][iquad][2]/((double) (ideg - 1));
      bary[3] = ordtet.s[ideg-1][iquad][3]/((double) (ideg - 1));
      d_quafun_distortion<ideg,ilag,ivar,nvar,ftype>(
         msh,ielem,power,bary,dmetvar,dpoivar,&qua0);
      (*qutet) += qua0 / nquad;
    }

    //for(int iquad = 0; iquad < nquad_Keaton6; iquad++){
    //  bary[0] = tuvquad_Keaton6[iquad][0];
    //  bary[1] = tuvquad_Keaton6[iquad][1];
    //  bary[2] = tuvquad_Keaton6[iquad][2];
    //  bary[3] = 1 - bary[0] - bary[1] - bary[2];
    //  d_quafun_distortion<ideg,ilag,ivar,nvar>(msh,ielem,bary,dmetvar,dpoivar,&qua0);
    //  (*qutet) += qua0 * wtquad_Keaton6[iquad];
    //}


  }else{
    double bary[4] = {0.25,0.25,0.25,0.25};
    d_quafun_distortion<1,ilag,ivar,nvar,ftype>(msh,ielem,power,bary,dmetvar,dpoivar,qutet);
  }
}



template <int ideg, int ilag, int nvar>
metqua_shell_d<ideg,ilag,nvar>
::metqua_shell_d(const Mesh &msh, int ipoin, int iele0, int power,
                  int mshell, 
                  int* __restrict__ nshell, int* __restrict__ lshell, 
                  const double* __restrict__ dmetvar, 
                  const double* __restrict__ dpoivar, 
                  SANS::SurrealS<nvar,double>*  qushe){

  int ierro;
  if(*nshell <= 0){
    int iopen;
    int iver = msh.getvertet<ideg>(iele0, ipoin);
    if(iver < 4 || iver >= 4 + 6*(edgnpps[ideg]-2))
      METRIS_THROW_MSG(TopoExcept(),"VERTEX NOT ON EDGE")
  
    int ied = (iver - 4) / (edgnpps[ideg] - 2);
  
    int ipoi1 = msh.tet2poi(iele0,lnoed3[ied][0]);
    int ipoi2 = msh.tet2poi(iele0,lnoed3[ied][1]);
    shell3(msh,ipoi1,ipoi2,iele0,mshell,nshell,lshell,&iopen);
    if(iopen >= 0) METRIS_THROW_MSG(TopoExcept(),"SHELL IS OPEN IN OPTIM")
  }
  
  constexpr int nrfld = tetnpps[ideg];
  auto nrfld_c = hana::int_c<nrfld>;

  *qushe = 0.0;
  for(int ishell = 0; ishell < *nshell; ishell++){
    int ielem = lshell[ishell];
    assert(ielem >= 0 && ielem < msh.nelem);

    bool ifnd = false;
    hana::while_(hana::less.than(nrfld_c), 4_c, [&](auto ivar_c){
    constexpr int ivar = ivar_c;
      if(msh.tet2poi(ielem,ivar) == ipoin){
        ifnd = true;

        SANS::SurrealS<3,double> qutet = 0;
        metqua_d<ideg,ilag,ivar,nvar>(msh,ielem,power,dmetvar,dpoivar,&qutet);
        if(qutet < 0) METRIS_THROW_MSG(GeomExcept(),
          "NEGATIVE ELEMENT QUALITY = "<<qutet<<" ielem = "<<ielem)
    
        (*qushe) += qutet;
      }
    return ivar_c+1_c;});
    assert(ifnd == true);
  }
}






#endif