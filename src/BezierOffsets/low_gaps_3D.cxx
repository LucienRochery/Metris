//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#include "low_gaps.hxx"

#include "../aux_exceptions.hxx"
#include "../linalg/matprods.hxx"
#include "../linalg/invmat.hxx"


#include "../ho_quadrature.hxx"
#include "../low_localization.hxx"
#include "../linalg/eigen.hxx"
#include "../SANS/Surreal/SurrealS.h"

#include "../low_ccoef.hxx"


namespace Metris{
  
void d2unittensor(const Mesh<MetricFieldAnalytical> &msh, int ielem, double *tens3sym_){

  METRIS_THROW_MSG(TODOExcept(), "Read this over and consolidate with existing")
  
  double metP1[6],dmetP1[18]; 
  double eigval[3],eigvec[9],rwork[10];
  double srmet[6];

  // Get metric used to get orientation and scaling 
  // Use metP1 for this
  double bary[4] = {0.25,0.25,0.25,0.25};
  double coop[3];
  eval3<3,1>(msh.coord,msh.tet2poi[ielem],msh.getBasis(),DifVar::None,DifVar::None,bary,coop,NULL,NULL);
  msh.anamet(NULL,coop,1,metP1,dmetP1);

  // Compute square root 
  geteigsym<3>(metP1,eigval,eigvec);
  for(int i = 0; i < 3; i++) eigval[i] = sqrt(eigval[i]);
  eig2met<3>(eigval,eigvec,srmet); 

  // 1. Get P1 element scaling and R*J_0, also metric and derivatives at barycentre
  double scale, tJ0tR[9];
  scalrotJ0<3>(msh,ielem,srmet,&scale,tJ0tR);


  SANS::SurrealS<3,double>  metS[6];
  SANS::SurrealS<3,double>  eigvalS[3];
  SANS::SurrealS<3,double>  eigvecS[9];
  
  getmet_dbl2SurS<3,3>(metP1,dmetP1,metS);

  //for(int ii = 0; ii < 6; ii++){
  //  metS[ii] = metP1[ii];
  //}
  //for(int jj = 0; jj < 3; jj++){
  //  for(int ii = 0; ii < 6; ii++){
  //    metS[ii].deriv(jj) = dmetP1[6*jj+ii];
  //  }
  //}


  geteigsym<3>(metS,eigvalS,eigvecS);


  for(int ii = 0; ii < 3; ii++){
    eigvalS[ii] = 1.0 / sqrt(eigvalS[ii]);
  }

  eig2met<3>(eigvalS,eigvecS,metS);

  getmet_dbl2SurS<3,3>(metS,metP1,dmetP1);

  //for(int ii = 0; ii < 6; ii++){
  //  metP1[ii] = metS[ii].value();
  //}

  //for(int jj = 0; jj < 3; jj++){
  //  for(int ii = 0; ii < 6; ii++){
  //    dmetP1[6*jj + ii] = metS[ii].deriv(jj);
  //  }
  //}


  double tens3sym[18];

  // First compute T = dM x_1 M
  tens3sym1X1mat3sym<3>(dmetP1,metP1,tens3sym);

  // Next compute S = N x_2 T x_1 N^T
  // with N = J_0^T R^T 
  mat3X2tens3sym1X1tmat3(tJ0tR,tens3sym,tens3sym_);
  METRIS_THROW_MSG(TODOExcept(),"Changed to non symmetric verify symmetric");

}

/*
Use P2 element intrinsic metric at all stages
This is mainly for validation
*/
void d2unittensor2(const MeshBase &msh, int ielem, double *tens3sym_){
  if(msh.curdeg != 2) METRIS_THROW(WArgExcept());

  constexpr int ideg = 2;
  constexpr int gdim = 3;
  constexpr int tdim = 3; 

  METRIS_THROW_MSG(TODOExcept(), "Read this over and consolidate with existing")

  int iprt = 1;

  double tJ0tR[9];
  SANS::SurrealS<3,double> intmetS[6];
  double intmet[6];

  // --------------------------------------------------------
  // 1. Compute rotation and scaling
  double bary[4] = {0.25,0.25,0.25,0.25};
  // a) Get intrinsic metric and physical derivatives
  getintmetxi<gdim,gdim,ideg>(msh.coord,msh.tet2poi[ielem],msh.getBasis(),bary,intmetS);
  if(iprt > 0){
    printf("1. print intmetS and derivatives\n");
    MeshArray1D<SANS::SurrealS<3,double>>(6,intmetS).print();
  }
    // State is now intmetS = M
  // b) Get Jacobian to compute scale
	double eval[3], jmat[9];
  eval3<3,2>(msh.coord,msh.tet2poi[ielem],FEBasis::Bezier,DifVar::Bary,DifVar::None,bary,eval,jmat,NULL);
  // c) Get scale
  for(int ii = 0; ii < 6; ii++) intmet[ii] = intmetS[ii].value();
  double detJ = detmat<3>(jmat);
  double detM = detmat<3>(intmet);
  double detJ0 = 1 / sqrt(2);
  //double scale = pow(detJ * sqrt(detM) / detJ0, 1.0/3.0);
  double scale = 1;
  printf("### RESTORE SCALE \n");
    // State is not intmet = M
  // d) Get rotation; for this we need to square root the metric.
  double eigval[3], eigvec[9];
  geteigsym<3>(intmet,eigval,eigvec);
  for(int ii = 0; ii < 3; ii++) eigval[ii] = sqrt(eigval[ii]);
  eig2met<3>(eigval,eigvec,intmet);
  matXsym<3>(jmat,intmet,tJ0tR);
    // State is not intmet = M^{1/2}

  
  // Physical derivatives of M^{-1/2}
  double dpMm12[18];
  // --------------------------------------------------------
  // 2. Compute M^{-1/2} and get physical derivatives
  // a) Compute M^{-1/2}
  SANS::SurrealS<3,double> eigvalS[3], eigvecS[9];
  geteigsym<3>(intmetS,eigvalS,eigvecS);
  for(int ii = 0; ii < 3; ii++) eigvalS[ii] = 1.0 / sqrt(eigvalS[ii]);
  eig2met<3>(eigvalS,eigvecS,intmetS);
    // State is now intmetS = M^{-1/2}
  // b) Compute physical derivatives of M^{-1/2}, i.e. diff of M^{-1/2} \circ F_K^{-1}:
  // = dF_K^{-1} x_1 dM^{-1/2}
  //   jmat ^           ^ intmetS.deriv(:)
  // i) gather diff M^{-1/2}
  double dbMm12[18];
  for(int idif = 0; idif < 3; idif++){
    for(int ii = 0; ii < 6; ii++){
      dbMm12[6*idif + ii] = intmetS[ii].deriv(idif);
    }
  }
  // ii) get dF_K^{-1}
	double invjmat[9];
  for(int i = 0; i < 9; i++) invjmat[i] = jmat[i];
  invmat(3,invjmat);
  // iii) conclude: tensor product
  mat3X1tens3sym1(dbMm12,invjmat,dpMm12);

  

  // --------------------------------------------------------
  // 3. Assemble everything: 
  // Denoting N = J°0^T R^T (tJ0tR)
  //          T = d_phy M^{-1/2} x_1 M^{-1/2}
  // compute N X_2 MM X_1 N

  double tens3sym[18];

  // a) compute T
  // Recall intmet is M^{1/2} but intmetS is M^{-1/2}:
  for(int ii = 0; ii < 6; ii++) intmet[ii] = intmetS[ii].value();
  tens3sym1X1mat3sym<3>(dpMm12,intmet,tens3sym);

  printf("%f\n",tJ0tR[0]);
  printf("%f\n",tJ0tR[1]);
  printf("%f\n",tJ0tR[2]);
  printf("%f\n",tJ0tR[3]);
  printf("%f\n",tJ0tR[4]);
  printf("%f\n",tJ0tR[5]);
  printf("%f\n",tJ0tR[6]);
  printf("%f\n",tJ0tR[7]);
  printf("%f\n",tJ0tR[8]);

  printf("%f\n",tens3sym[0]);
  printf("%f\n",tens3sym[1]);
  printf("%f\n",tens3sym[2]);
  printf("%f\n",tens3sym[3]);
  printf("%f\n",tens3sym[4]);
  printf("%f\n",tens3sym[5]);
  printf("%f\n",tens3sym[6]);
  printf("%f\n",tens3sym[7]);
  printf("%f\n",tens3sym[8]);
  printf("%f\n",tens3sym[9]);
  printf("%f\n",tens3sym[10]);
  printf("%f\n",tens3sym[11]);
  printf("%f\n",tens3sym[12]);
  printf("%f\n",tens3sym[13]);
  printf("%f\n",tens3sym[14]);
  printf("%f\n",tens3sym[15]);
  printf("%f\n",tens3sym[16]);
  printf("%f\n",tens3sym[17]);

  // Next compute S = N x_2 T x_1 N^T
  double tens3[27];
  mat3X2tens3sym1X1tmat3(tens3sym,tJ0tR,tens3);

  printf("%f\n",tens3[0]);
  printf("%f\n",tens3[1]);
  printf("%f\n",tens3[2]);
  printf("%f\n",tens3[3]);
  printf("%f\n",tens3[4]);
  printf("%f\n",tens3[5]);
  printf("%f\n",tens3[6]);
  printf("%f\n",tens3[7]);
  printf("%f\n",tens3[8]);
  printf("%f\n",tens3[9]);
  printf("%f\n",tens3[10]);
  printf("%f\n",tens3[11]);
  printf("%f\n",tens3[12]);
  printf("%f\n",tens3[13]);
  printf("%f\n",tens3[14]);
  printf("%f\n",tens3[15]);
  printf("%f\n",tens3[16]);
  printf("%f\n",tens3[17]);
  printf("%f\n",tens3[18]);
  printf("%f\n",tens3[19]);
  printf("%f\n",tens3[20]);
  printf("%f\n",tens3[21]);
  printf("%f\n",tens3[22]);
  printf("%f\n",tens3[23]);
  printf("%f\n",tens3[24]);
  printf("%f\n",tens3[25]);
  printf("%f\n",tens3[26]);

  printf("Check 3 symmetry:");
  for(int kk = 0; kk < 3; kk++){
    for(int ii = 0; ii < 3; ii++){
      for(int jj = ii + 1; jj < 3; jj++){
        printf("%d %d / %d %d : %f %f \n",ii,jj,jj,ii,
                    tens3[9*ii+3*jj+kk],tens3[9*jj+3*ii+kk]);
      }
    }
  }
}




///*
//Compute element orientation and scale using P1 subjacent element.
//See scitech24 paper. 
//Also return metbar, the metric at the P1 barycentre. 
//And derivatives if dmetbar != NULL.
//
//Note: we could return RJ_0 directly, and thus never bother with J_0 at all. 
//*/
//void scalrot3(const Mesh &msh, int ielem , 
//              double* __restrict__ metbar,
//              double* __restrict__ dmetbar,
//              double* __restrict__ scale ,
//              double* __restrict__ rotmat){
//  if(msh.ilogmet != 1 && msh.ianamet <=0)
//    METRIS_THROW_MSG(WArgExcept(), "Set metric to log before calling scalrot3");
//
//  // Used in both
//  double jmat[3][3];
//  for(int i = 0; i < 3; i++){
//    jmat[0][i] = msh.coord[msh.tet2poi(ielem,1)][i]
//               - msh.coord[msh.tet2poi(ielem,0)][i];
//    jmat[1][i] = msh.coord[msh.tet2poi(ielem,2)][i]
//               - msh.coord[msh.tet2poi(ielem,0)][i];
//    jmat[2][i] = msh.coord[msh.tet2poi(ielem,3)][i]
//               - msh.coord[msh.tet2poi(ielem,0)][i];
//  }
//  double detJ1 = detmat<3>(jmat[0]);
//  double detJ0 = 1 / sqrt(2);
//
//
//  double bary[4] = {0.25,0.25,0.25,0.25};
//
//
//  /*Get orientation*/
//  getintmetxi<1,0>(msh.coord,msh.tet2poi[ielem],bary,metbar);
//  double eigval[3],eigvec[9],rwork[10];
//  geteigsym(metbar,10,rwork,eigval,eigvec);
//  for(int i = 0; i < 3; i++) eigval[i] = sqrt(eigval[i]);
//  eig2met(eigval,eigvec,metbar); 
//
//
//  const double invtJ_0[3][3] = {
//    -0.57735026918962562,-0.57735026918962562 , 0                 ,
//    1                   ,-1                   , 0                 ,
//    -0.40824829046386302,-0.4082482904638630  , 1.2247448713915889};
//  double jmat1T[3][3];
//  matXmat<3>(invtJ_0[0],jmat[0],jmat1T[0]);
//
//  sym3tmat(metbar,jmat1T[0],rotmat);
//
//
//
//  /*
//  Scale. If discrete metric, interpolate directly from front vertices as log-met, then expmet
//  Otherwise, get 
//  */
//  double detM;
//  if(msh.ianamet <= 0){
//    if(dmetbar == NULL){
//      eval3<6,1,1,0,0>(msh.met,msh.tet2poi[ielem],bary,metbar,NULL,NULL);
//      detM = exp(metbar[0] + metbar[2] + metbar[5]);
//      getexpmet_inp(metbar);
//    }else{
//      double buf1[6], buf2[18];
//      eval3<6,1,1,1,0>(msh.met,msh.tet2poi[ielem],bary,buf1,buf2,NULL);
//      detM = exp(buf1[0] + buf1[2] + buf1[5]);
//      getexpmet_cpy_d(buf1,buf2,metbar,dmetbar);
//    }
//  }else{
//    double coop[3];
//    eval3<3,1,1,0,0>(msh.coord,msh.tet2poi[ielem],bary,coop,NULL,NULL);
//    msh.anamet(NULL,coop,dmetbar!=NULL,metbar,dmetbar);
//    detM = detsym<3>(metbar);
//  }
//  // tr(log M) = log(det(M))
//  if(detM < 1.0e-16) METRIS_THROW_MSG(GeomExcept(),"SINGULAR METRIC " << detM
//  <<" log-metric "<<metbar[0]<<" "<<metbar[1]<<" "<<metbar[2]<<" "
//                  <<metbar[3]<<" "<<metbar[4]<<" "<<metbar[5]<<" ")
//  double detMm12 = 1.0 / sqrt(detM);
//
//  *scale = pow(detJ1 / (detMm12 * detJ0),1.0/3.0);
//
//}
//










} // End namespace

