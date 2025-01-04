//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "linalg/eigen.hxx"
#include "linalg/utils.hxx"
#include "../SANS/Surreal/SurrealS.h"

// Allowed prototype: 
//  void (*anamet)(void* ctx, double x, double y, double z, int idif1, double *met, double *dmet);
// Met stored 1 2 4  
//              3 5
//                6
// Dmet [i][j] = d_i M_j

namespace Metris{


void anamet3D_1(void *ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
  met[0] = 1.0/(scale*scale);
  met[1] = 0.0;
  met[2] = 1.0/(scale*scale);
  met[3] = 0.0;
  met[4] = 0.0;
  met[5] = 1.0/(scale*scale);

  if(idif1 > 0){
    for(int ii = 0; ii < 3; ii++){
      for(int jj = 0; jj < 6; jj++){
        dmet[6*ii + jj] = 0;
      }
    }
  }
}



/*
Sinusoidal BL cf "HIGH-ORDER METRIC INTERPOLATION FOR CURVED R-ADAPTION BY DISTORTION MINIMIZATION" 
Guillermo Aparicio-Estrems Abel Gargallo-Peiro Xevi Roca
*/
void anamet3D_2(void *ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
  double x0 = 0;
  double y0 = 0;
  double z0 = 0;

  SANS::SurrealS<3,double> x = crd[0] - x0;
  x.deriv(0) = 1;
  x.deriv(1) = 0;
  x.deriv(2) = 0;

  SANS::SurrealS<3,double> y = crd[1] - y0;
  y.deriv(0) = 0;
  y.deriv(1) = 1;
  y.deriv(2) = 0;

  SANS::SurrealS<3,double> z = crd[2] - z0;
  z.deriv(0) = 0;
  z.deriv(1) = 0;
  z.deriv(2) = 1;


  double hmin = 0.1;
  double gamma = 2;
  SANS::SurrealS<3,double> hh = hmin + gamma*abs(z);
  hh *= scale;

  SANS::SurrealS<3,double> eigval[3] = {1.0/(scale*scale), 1.0/(scale*scale), 1.0/(hh*hh)};

  SANS::SurrealS<3,double> eigvec[3][3];

  double pi = 3.1415;

//  SANS::SurrealS<3,double> phi[3] = {x, y, (10*z - cos(2*pi*x)*cos(2*pi*y))/(100.0 + 8.0*pi*pi)};
  SANS::SurrealS<3,double> d1phi3 = 2*pi*sin(2*pi*x)*cos(2*pi*y)/(100.0 + 8.0*pi*pi);
  SANS::SurrealS<3,double> d2phi3 = 2*pi*cos(2*pi*x)*sin(2*pi*y)/(100.0 + 8.0*pi*pi);
  SANS::SurrealS<3,double> d3phi3 = 10.0/(100.0 + 8.0*pi*pi);

  eigvec[0][0] = 1;
  eigvec[1][0] = 0; 
  eigvec[2][0] = 0; 

  eigvec[0][1] = 0;
  eigvec[1][1] = 1; 
  eigvec[2][1] = 0; 

  SANS::SurrealS<3,double> nrm = d1phi3*d1phi3
                               + d2phi3*d2phi3
                               + d3phi3*d3phi3;
  
  //std::cout<<"norm"<<nrm<<"\n";
  //for(int idif = 0; idif < 3; idif++){
  //  printf("d%d phi\n",idif);
  //  for(int i = 0; i< 3; i++){
  //    printf(" %15.7e ",phi[i].deriv(idif));
  //  }
  //  printf("\n");
  //}
  nrm = sqrt(nrm);
  eigvec[0][2] = d1phi3/nrm;
  eigvec[1][2] = d2phi3/nrm;
  eigvec[2][2] = d3phi3/nrm;

  //for(int idif = 0; idif < 3; idif ++){
  //printf("Debug eigvec d%d: \n",idif);
  //for(int i = 0; i < 3; i++){
  //  for(int j = 0; j < 3; j++){
  //    printf(" %15.7e ",eigvec[i][j].deriv(idif));
  //  }
  //  printf("\n");
  //}
  //}

  SANS::SurrealS<3,double> metS[6];
  eig2met<3,SANS::SurrealS<3,double>>(eigval,eigvec[0],metS);

  getmet_SurS2dbl<3>(metS,met,dmet);

  //for(int jj = 0; jj < 6; jj++){
  //  met[jj] = metS[jj].value();
  //}

  //if(idif1 > 0){
  //  for(int ii = 0; ii < 3; ii++){
  //    for(int jj = 0; jj < 6; jj++){
  //      dmet[6*ii + jj] = metS[jj].deriv(ii);
  //    }
  //  }
  //}



//  SANS::SurrealS<3,double> r = sqrt( x*x + y*y + z*z );
//
//  SANS::SurrealS<3,double> tmp = y / r;
//
//  SANS::SurrealS<3,double> theta = acos(tmp);
//
//  tmp = (crd[0] - x0) / sqrt( x*x + y*y );
//  SANS::SurrealS<3,double> phi = atan(tmp);
//  if(y.value() < 0) phi *= -1;
//
//  // For x = r sin theta cos phi
//  //     y = r sin theta sin phi
//  //     z = r cos theta
//  SANS::SurrealS<3,double> rot[3][3];
//
//  rot[0][0] = sin(theta)*cos(phi);
//  rot[0][1] = sin(theta)*sin(phi);
//  rot[0][2] = cos(theta);
//
//  rot[1][0] = cos(theta)*cos(phi);
//  rot[1][1] = cos(theta)*sin(phi);
//  rot[1][2] = -sin(phi);
//
//  rot[2][0] = -sin(phi);
//  rot[2][1] = cos(phi);
//  rot[2][2] = 0;
// 
//  SANS::SurrealS<3,double> eigval[3] = {1, 2, 3};
//
//  SANS::SurrealS<3,double> metS[6];
//  eig2met<SANS::SurrealS<3,double>>(eigval,rot[0],metS);
//
//  for(int jj = 0; jj < 6; jj++){
//    met[jj] = metS[jj].value();
//  }
//
//  if(idif1 > 0){
//    for(int ii = 0; ii < 3; ii++){
//      for(int jj = 0; jj < 6; jj++){
//        dmet[6*ii + jj] = metS[jj].deriv(ii);
//      }
//    }
//  }
//

}




/*
Cylindrical around axis z 
Centered around -1, -1, -1 to avoid singularity on common cube cases 
*/
void anamet3D_3(void *ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
  double x0 = -0.6;
  double y0 = -0.6;

  SANS::SurrealS<2,double> x = crd[0] - x0;
  x.deriv(0) = 1;
  x.deriv(1) = 0;

  SANS::SurrealS<2,double> y = crd[1] - y0;
  y.deriv(0) = 0;
  y.deriv(1) = 1;

  SANS::SurrealS<2,double> r = sqrt( x*x + y*y );

  SANS::SurrealS<2,double> tmp = y / r;

  SANS::SurrealS<2,double> theta = 2*3.141592653589793238462643383279502884*acos(tmp);


  // For x = r sin theta cos phi
  //     y = r sin theta sin phi
  //     z = r cos theta
  SANS::SurrealS<2,double> rot[3][3];

  rot[0][0] = cos(theta);
  rot[0][1] = -sin(theta);
  rot[0][2] = 0;

  rot[1][0] = sin(theta);
  rot[1][1] = cos(theta);
  rot[1][2] = 0;

  rot[2][0] = 0;
  rot[2][1] = 0;
  rot[2][2] = 1;
 
  SANS::SurrealS<2,double> eigval[3] = {1.0/(scale*scale), 2/(scale*scale), 10/(scale*scale)};

  SANS::SurrealS<2,double> metS[6];
  eig2met<3,SANS::SurrealS<2,double>>(eigval,rot[0],metS);



  for(int jj = 0; jj < 6; jj++){ 
    met[jj] = metS[jj].value();
  }


  if(idif1 > 0){
    for(int ii = 0; ii < 2; ii++){
      for(int jj = 0; jj < 6; jj++){
        dmet[6*ii + jj] = metS[jj].deriv(ii);
      }
    }

    for(int jj = 0; jj < 6; jj++){
      dmet[6*2 + jj] = 0;
    }
  }
}


/*
Cylindrical around axis z 
Centered around -1, -1, -1 to avoid singularity on common cube cases 
*/
void anamet3D_4(void *ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
  double x0 = 0;
  double y0 = 0;
  double z0 = 0;

  SANS::SurrealS<3,double> x = crd[0] - x0;
  x.deriv(0) = 1;
  x.deriv(1) = 0;
  x.deriv(2) = 0;

  SANS::SurrealS<3,double> y = crd[1] - y0;
  y.deriv(0) = 0;
  y.deriv(1) = 1;
  y.deriv(2) = 0;

  SANS::SurrealS<3,double> z = crd[2] - z0;
  z.deriv(0) = 0;
  z.deriv(1) = 0;
  z.deriv(2) = 1;

  SANS::SurrealS<3,double> r = 1 + x*x + y*y + z*z;
  r *= scale*scale;

  SANS::SurrealS<3,double> metS[6] = {1.0/r, 0, 1.0/r, 0, 0, 1.0/r};
  // SANS::SurrealS<3,double> metS[6] = {1/(scal*scal), 0, 1/(scal*scal), 0, 0, 1/(scal*scal)};


  getmet_SurS2dbl<3>(metS,met,dmet);

  //for(int jj = 0; jj < 6; jj++){
  //  met[jj] = metS[jj].value();
  //}

  //if(idif1 > 0){
  //  for(int ii = 0; ii < 3; ii++){
  //    for(int jj = 0; jj < 6; jj++){
  //      dmet[6*ii + jj] = metS[jj].deriv(ii);
  //    }
  //  }
  //}
}



// circle BL
void anamet3D_5(void *ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
  const double pi = 3.141592653589793238462643383279502884;
  SANS::SurrealS<2,double> X[2];
  X[0] = crd[0];
  X[0].deriv(0) = 1;
  X[0].deriv(1) = 0;

  X[1] = crd[1];
  X[1].deriv(0) = 0;
  X[1].deriv(1) = 1;

  ///const double r0 = 0.5;

  SANS::SurrealS<2,double> r = sqrt(X[0]*X[0] + X[1]*X[1]) + 0.01;// - r0 + 1.0e-2;

  double hy_min = 0.01;
  double hy_max = 0.1;
  double hx = scale*0.5;
  SANS::SurrealS<2,double> hy = scale * (r * hy_max + (1.0 - r) * hy_min);
  double hz = scale * 0.02;

  SANS::SurrealS<2,double> eigval[3] = {1.0/(hy*hy), 1.0/(hx*hx), 1.0/(hz*hz)};
  SANS::SurrealS<2,double> eigvec[9];

  SANS::SurrealS<2,double> y[2] = {X[0] / r, X[1] / r};
  SANS::SurrealS<2,double> theta;
  if(y[0].value() > 0){
    theta = atan(y[1]/y[0]);
  }else{
    theta = pi + atan(y[1]/y[0]);
  }

  if(ctx != NULL){
    printf("anamet2D_2 debug r = %.6f theta = %.6f\n", r.value(), theta.value());
  }

  // eig2met is in R^T D R format. Worst case we are using -theta. 
  eigvec[0] =  cos(theta);
  eigvec[1] =  sin(theta); 
  eigvec[2] =  0; 
  eigvec[3] = -sin(theta);
  eigvec[4] =  cos(theta);
  eigvec[5] =  0; 
  eigvec[6] = 0;
  eigvec[7] = 0;
  eigvec[8] = 1;

  SANS::SurrealS<2,double> metS[3];
  eig2met<3,SANS::SurrealS<2,double>>(eigval,eigvec,metS);

  getmet_SurS2dbl<3,2>(metS,met,dmet);

  //for(int ii = 0; ii < 3; ii ++) met[ii] = metS[ii].value();

 //if(idif1 > 0){
 //  for(int ii = 0; ii < 2; ii++){
 //    for(int jj = 0; jj < 3; jj++){
 //      dmet[3*ii + jj] = metS[jj].deriv(ii);
 //    }
 //  }
 //}
}



} // End namespace
