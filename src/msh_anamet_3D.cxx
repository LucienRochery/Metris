//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include "msh_anamet.hxx"

#include "linalg/eigen.hxx"
#include "linalg/utils.hxx"
#include "SANS/Surreal/SurrealS.h"


#include <cmath>

// Allowed prototype:
//  void (*anamet)(void* ctx, double x, double y, double z, int idif1, double *met, double *dmet);
// Met stored 1 2 4
//              3 5
//                6
// Dmet [i][j] = d_i M_j

namespace Metris{


void anamet3D_1([[maybe_unused]] const AnaMetCtx* ctx, [[maybe_unused]] const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
  // Not too coarse at the scale of 1
  double h0 = 0.05;

  met[0] = 1/(h0*h0*scale*scale);
  met[1] = 0.0;
  met[2] = 1/(h0*h0*scale*scale);
  met[3] = 0.0;
  met[4] = 0.0;
  met[5] = 1/(h0*h0*scale*scale);

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
void anamet3D_2([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
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
  SANS::SurrealS<3,double> hh = hmin + gamma*std::abs(z);
  hh *= scale;

  SANS::SurrealS<3,double> eigval[3] = {1.0/(scale*scale), 1.0/(scale*scale), 1.0/(hh*hh)};

  SANS::SurrealS<3,double> eigvec[3][3];

  double pi = 3.141592653589793;

//  SANS::SurrealS<3,double> phi[3] = {x, y, (10*z - cos(2*pi*x)*cos(2*pi*y))/(100.0 + 8.0*pi*pi)};
  SANS::SurrealS<3,double> d1phi3 = 2*pi*sin(2*pi*x)*cos(2*pi*y)/(100.0 + 8.0*pi*pi);
  SANS::SurrealS<3,double> d2phi3 = 2*pi*cos(2*pi*x)*sin(2*pi*y)/(100.0 + 8.0*pi*pi);
  SANS::SurrealS<3,double> d3phi3 = 10.0/(100.0 + 8.0*pi*pi);

  eigvec[0][0] = 1;
  eigvec[0][1] = 0;
  eigvec[0][2] = 0;

  eigvec[1][0] = 0;
  eigvec[1][1] = 1;
  eigvec[1][2] = 0;

  SANS::SurrealS<3,double> nrm = d1phi3*d1phi3
                               + d2phi3*d2phi3
                               + d3phi3*d3phi3;

  //std::cout<<"norm"<<nrm<<"\n";
  //for(int idif = 0; idif < 3; idif++){
  //  printf("d{} phi\n",idif);
  //  for(int i = 0; i< 3; i++){
  //    printf(" {:15.7e} ",phi[i].deriv(idif));
  //  }
  //  printf("\n");
  //}
  nrm = sqrt(nrm);
  eigvec[2][0] = d1phi3/nrm;
  eigvec[2][1] = d2phi3/nrm;
  eigvec[2][2] = d3phi3/nrm;

  //for(int idif = 0; idif < 3; idif ++){
  //printf("Debug eigvec d{}: \n",idif);
  //for(int i = 0; i < 3; i++){
  //  for(int j = 0; j < 3; j++){
  //    printf(" {:15.7e} ",eigvec[i][j].deriv(idif));
  //  }
  //  printf("\n");
  //}
  //}

  SANS::SurrealS<3,double> metS[6];
  eig2met<3,SANS::SurrealS<3,double>>(eigval,eigvec[0],metS);

  getmet_SurS2dbl<3>(metS,met,idif1 > 0 ? dmet : NULL);

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
void anamet3D_3([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
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


  getmet_SurS2dbl<3>(metS,met,idif1 > 0 ? dmet : NULL);

}


/*
Cylindrical around axis z
Centered around -1, -1, -1 to avoid singularity on common cube cases
*/
void anamet3D_4([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
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


  getmet_SurS2dbl<3>(metS,met,idif1 > 0 ? dmet : NULL);

}



/*
UGAWG (https://github.com/UGAWG/adapt-benchmarks/tree/master/cube)
cone-cone metric:
   + h_x^-2   0      0    +
M =|   0    h_y^-2   0    |
   +   0      0    h_z^-2 +

where

h_x = 0.1
h_y = 0.1
h_z = 0.1
*/
void anamet3D_5([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){

  const double h0 = 0.001;

  SANS::SurrealS<3,double> z = crd[2] - h0;
  z.deriv(0) = 0;
  z.deriv(1) = 0;
  z.deriv(2) = 1;

  double hx = 0.1*scale;
  double hy = 0.1*scale;
  double hz = 0.1*scale;

  SANS::SurrealS<3,double> metS[6] = {1.0/(hx*hx), 0, 1.0/(hy*hy), 0, 0, 1.0/(hz*hz)};

  getmet_SurS2dbl<3>(metS,met,idif1 > 0 ? dmet : NULL);

}


/*
UGAWG (https://github.com/UGAWG/adapt-benchmarks/tree/master/cube)
cube & cube-cylinder linear metric:
   + h_x^-2   0      0    +
M =|   0    h_y^-2   0    |
   +   0      0    h_z^-2 +

where

h_x = 0.1
h_y = 0.1
h_z = h0 + 2*(0.1-h0)*abs(z-0.5)
h0 = 0.001
*/
void anamet3D_6([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){

  const double h0 = 0.001;

  SANS::SurrealS<3,double> z = crd[2] - h0;
  z.deriv(0) = 0;
  z.deriv(1) = 0;
  z.deriv(2) = 1;

  double hx = 0.1*scale;
  double hy = 0.1*scale;
  SANS::SurrealS<3,double> hz = (h0 + 2*(0.1-h0)*std::abs(z-0.5))*scale;

  SANS::SurrealS<3,double> metS[6] = {1.0/(hx*hx), 0, 1.0/(hy*hy), 0, 0, 1.0/(hz*hz)};

  getmet_SurS2dbl<3>(metS,met,idif1 > 0 ? dmet : NULL);

}


/*
UGAWG (https://github.com/UGAWG/adapt-benchmarks/tree/master/cube)
cube-cylinder polar-1:
   + cos(t) -sin(t)   0    ++ h_r^-2   0      0    ++ cos(t)  sin(t)   0    +
M =| sin(t)  cos(t)   0    ||   0    h_t^-2   0    ||-sin(t)  cos(t)   0    |
   +   0       0      1    ++   0      0    h_z^-2 ++   0       0      1    +

where

r=sqrt(x^2+y^2)
t=atan2(y,x)
h_z = 0.1
h_t = 0.1
h_r = h0 + 2*(0.1-h0)*abs(r-0.5)
h0 = 0.001
# t is in the theta direction
# r is radial direction
*/
void anamet3D_7([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){

  double dx = ctx == NULL ? 0 : ctx->dx;
  double dy = ctx == NULL ? 0 : ctx->dy;
  double dz = ctx == NULL ? 0 : ctx->dz;

  SANS::SurrealS<3,double> x = std::abs(crd[0] + dx);
  x.deriv(0) = crd[0] + dx > 0 ? 1 : -1;
  x.deriv(1) = 0;
  x.deriv(2) = 0;

  SANS::SurrealS<3,double> y = std::abs(crd[1] + dy);
  y.deriv(0) = 0;
  y.deriv(1) = crd[1] + dy > 0 ? 1 : -1;
  y.deriv(2) = 0;

  SANS::SurrealS<3,double> z = std::abs(crd[2] + dz);
  z.deriv(0) = 0;
  z.deriv(1) = 0;
  z.deriv(2) = crd[2] + dz > 0 ? 1 : -1;

  SANS::SurrealS<3,double> r = sqrt(x*x + y*y);
  SANS::SurrealS<3,double> t = atan2(y,x);

  double hz = 0.1*scale;
  double ht = 0.1*scale;
  const double h0 = 0.001;
  SANS::SurrealS<3,double> hr = (h0 + 2*(0.1-h0)*std::abs(r-0.5))*scale;

  SANS::SurrealS<3,double> eigval[3], eigvec[3][3];
  eigval[0] = 1.0/(hr*hr);
  eigval[1] = 1.0/(ht*ht);
  eigval[2] = 1.0/(hz*hz);

  // eig2met computes R^T D R, so eigvec is the right side matrix
  eigvec[0][0] =  cos(t);
  eigvec[0][1] =  sin(t);
  eigvec[0][2] = 0;

  eigvec[1][0] = -sin(t);
  eigvec[1][1] =  cos(t);
  eigvec[1][2] = 0;

  eigvec[2][0] = 0;
  eigvec[2][1] = 0;
  eigvec[2][2] = 1;

  SANS::SurrealS<3,double> metS[6];
  eig2met<3,SANS::SurrealS<3,double>>(eigval, eigvec[0], metS);

  getmet_SurS2dbl<3>(metS,met,idif1 > 0 ? dmet : NULL);
}


/*
UGAWG (https://github.com/UGAWG/adapt-benchmarks/tree/master/cube)
cube-cylinder polar-2:
A modified polar-1 metric that is easier to satisfy with high-quality elements by refining along theta near the layer,
d = (0.6 - r)*10
h_t = (d < 0) ? (0.1) : (d*(1 / 40) + (1 - d)*0.1)
*/
void anamet3D_8([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){

  double dx = ctx == NULL ? 0 : ctx->dx;
  double dy = ctx == NULL ? 0 : ctx->dy;
  double dz = ctx == NULL ? 0 : ctx->dz;

  SANS::SurrealS<3,double> x = std::abs(crd[0] + dx);
  x.deriv(0) = crd[0] + dx > 0 ? 1 : -1;
  x.deriv(1) = 0;
  x.deriv(2) = 0;

  SANS::SurrealS<3,double> y = std::abs(crd[1] + dy);
  y.deriv(0) = 0;
  y.deriv(1) = crd[1] + dy > 0 ? 1 : -1;
  y.deriv(2) = 0;

  SANS::SurrealS<3,double> z = std::abs(crd[2] + dz);
  z.deriv(0) = 0;
  z.deriv(1) = 0;
  z.deriv(2) = crd[2] + dz > 0 ? 1 : -1;

  SANS::SurrealS<3,double> r = sqrt(x*x + y*y);
  SANS::SurrealS<3,double> t = atan2(y,x);

  double hz = 0.1*scale;
  SANS::SurrealS<3,double> d = (0.6 - r)*10;
  SANS::SurrealS<3,double> ht;
  if(d < 0) ht = 0.1;
  else      ht = d / 40 + (1 - d)*0.1;
  ht *= scale;
  const double h0 = 0.001;
  SANS::SurrealS<3,double> hr = (h0 + 2*(0.1-h0)*std::abs(r-0.5))*scale;

  SANS::SurrealS<3,double> eigval[3], eigvec[3][3];
  eigval[0] = 1.0/(hr*hr);
  eigval[1] = 1.0/(ht*ht);
  eigval[2] = 1.0/(hz*hz);

  // eig2met computes R^T D R, so eigvec is the right side matrix
  eigvec[0][0] =  cos(t);
  eigvec[0][1] =  sin(t);
  eigvec[0][2] = 0;

  eigvec[1][0] = -sin(t);
  eigvec[1][1] =  cos(t);
  eigvec[1][2] = 0;

  eigvec[2][0] = 0;
  eigvec[2][1] = 0;
  eigvec[2][2] = 1;

  SANS::SurrealS<3,double> metS[6];
  eig2met<3,SANS::SurrealS<3,double>>(eigval, eigvec[0], metS);

  getmet_SurS2dbl<3>(metS,met,idif1 > 0 ? dmet : NULL);
}


// Boundary-layer along x centered at 0
void anamet3D_9([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){

  SANS::SurrealS<3,double> X[3];
  const double x0 = 0.0;
  X[0] = std::abs(crd[0] - x0);
  X[0].deriv(0) = crd[0] >= x0 ? 1 : -1;
  X[0].deriv(1) = 0;
  X[0].deriv(2) = 0;

  X[1] = crd[1];
  X[1].deriv(0) = 0;
  X[1].deriv(1) = 1;
  X[1].deriv(2) = 0;

  X[2] = crd[2];
  X[2].deriv(0) = 0;
  X[2].deriv(1) = 0;
  X[2].deriv(2) = 1;


  double hx_min = 0.001;
  double hx_max = 0.1;
  SANS::SurrealS<3,double> hx = scale*(X[0]*(hx_max - hx_min) + hx_min);
  double hy = scale*0.1;
  double hz = scale*0.1;


  SANS::SurrealS<3,double> eigval[3] = {1.0/(hx*hx), 1.0/(hy*hy), 1.0/(hz*hz)};
  SANS::SurrealS<3,double> eigvec[9];

  eigvec[0] = 1.0;
  eigvec[1] = 0.0;
  eigvec[2] = 0.0;

  eigvec[3] = 0.0;
  eigvec[4] = 1.0;
  eigvec[5] = 0.0;

  eigvec[6] = 0.0;
  eigvec[7] = 0.0;
  eigvec[8] = 1.0;

  SANS::SurrealS<3,double> metS[6];
  eig2met<3,SANS::SurrealS<3,double>>(eigval,eigvec,metS);
  getmet_SurS2dbl<3>(metS,met,idif1 > 0 ? dmet : NULL);

}

// Boundary-layer along x centered at 0.5
// Merge with previous once (if) ctx is implemented.
void anamet3D_10([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){

  SANS::SurrealS<3,double> X[3];
  const double x0 = 0.5;
  X[0] = std::abs(crd[0] - x0);
  X[0].deriv(0) = crd[0] >= x0 ? 1 : -1;
  X[0].deriv(1) = 0;
  X[0].deriv(2) = 0;

  X[1] = crd[1];
  X[1].deriv(0) = 0;
  X[1].deriv(1) = 1;
  X[1].deriv(2) = 0;

  X[2] = crd[2];
  X[2].deriv(0) = 0;
  X[2].deriv(1) = 0;
  X[2].deriv(2) = 1;


  double hx_min = 0.001;
  double hx_max = 0.1;
  SANS::SurrealS<3,double> hx = scale*(X[0]*(hx_max - hx_min) + hx_min);
  double hy = scale*0.1;
  double hz = scale*0.1;


  SANS::SurrealS<3,double> eigval[3] = {1.0/(hx*hx), 1.0/(hy*hy), 1.0/(hz*hz)};
  SANS::SurrealS<3,double> eigvec[9];

  eigvec[0] = 1.0;
  eigvec[1] = 0.0;
  eigvec[2] = 0.0;

  eigvec[3] = 0.0;
  eigvec[4] = 1.0;
  eigvec[5] = 0.0;

  eigvec[6] = 0.0;
  eigvec[7] = 0.0;
  eigvec[8] = 1.0;

  SANS::SurrealS<3,double> metS[6];
  eig2met<3,SANS::SurrealS<3,double>>(eigval,eigvec,metS);
  getmet_SurS2dbl<3>(metS,met,idif1 > 0 ? dmet : NULL);

}

// Euclidean-anisotropic, i.e. constant metric but anisotropic
void anamet3D_11([[maybe_unused]] const AnaMetCtx* ctx, [[maybe_unused]] const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
  // Not too coarse at the scale of 1
  double hx = 0.05;
  double hy = hx;
  double anisoRatioZ = 16.;

  double hz = hx/anisoRatioZ;
  // hy = hz;

  met[0] = 1/(hx*hx*scale*scale);
  met[1] = 0.0;
  met[2] = 1/(hy*hy*scale*scale);
  met[3] = 0.0;
  met[4] = 0.0;
  met[5] = 1/(hz*hz*scale*scale);

  if(idif1 > 0){
    for(int ii = 0; ii < 3; ii++){
      for(int jj = 0; jj < 6; jj++){
        dmet[6*ii + jj] = 0;
      }
    }
  }
}

// BL from Yano's thesis, p. 88
// x: normal direction to BL
// y and z: tangent to BL

void anamet3D_12([[maybe_unused]] const AnaMetCtx* ctx, [[maybe_unused]] const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){

  SANS::SurrealS<3,double> X[3];

  X[0] = crd[0];
  X[0].deriv(0) = 1;
  X[0].deriv(1) = 0;
  X[0].deriv(2) = 0;

  X[1] = crd[1];
  X[1].deriv(0) = 0;
  X[1].deriv(1) = 1;
  X[1].deriv(2) = 0;

  X[2] = crd[2];
  X[2].deriv(0) = 0;
  X[2].deriv(1) = 0;
  X[2].deriv(2) = 1;

  // parameters
  double epsilon = 0.2;
  double p = 1.;
  double delta = epsilon * (p + 3./2.) * (1. - 1./(2*p*p + 7*p + 6));
  double k1 = 1./delta;

  SANS::SurrealS<3,double> h1 = scale * exp(k1*X[0]);

  double betay = 1.;
  double betaz = 1.;

  SANS::SurrealS<3,double> ARy = 1./(epsilon * pow(betay,1./(p+1.))) * exp(-X[0]/(epsilon*(p+1)));
  SANS::SurrealS<3,double> ARz = 1./(epsilon * pow(betaz,1./(p+1.))) * exp(-X[0]/(epsilon*(p+1)));

  SANS::SurrealS<3,double> h2 = ARy * h1;
  SANS::SurrealS<3,double> h3 = ARz * h1;

  SANS::SurrealS<3,double> metS[6];
  metS[0] = 1./(h1*h1);
  metS[1] = 0;
  metS[2] = 1./(h2*h2);
  metS[3] = 0.;
  metS[4] = 0.;
  metS[5] = 1./(h3*h3);

  getmet_SurS2dbl<3>(metS,met,idif1 > 0 ? dmet : NULL);

}

} // End namespace
