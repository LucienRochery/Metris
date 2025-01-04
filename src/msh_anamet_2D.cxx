//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_anamet.hxx"

#include "linalg/eigen.hxx"
#include "linalg/utils.hxx"

#include "../SANS/Surreal/SurrealS.h"
#include "aux_exceptions.hxx"
#include <cmath>


namespace Metris{

void anamet2D_1(void *ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
  met[0] = 1.0/(scale*scale);
  met[1] = 0.0;
  met[2] = 1.0/(scale*scale);

  if(idif1 > 0){
    for(int ii = 0; ii < 2; ii++){
      for(int jj = 0; jj < 3; jj++){
        dmet[3*ii + jj] = 0;
      }
    }
  }
}

// circle
void anamet2D_2(void *ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
  const double pi = 3.141592653589793238462643383279502884;
  SANS::SurrealS<2,double> X[2];
  X[0] = crd[0] - 0.5;
  //X[0] = crd[0] + 0.01;
  X[0].deriv(0) = 1;
  X[0].deriv(1) = 0;

  X[1] = crd[1] - 0.5;
  //X[1] = crd[1] + 0.01;
  X[1].deriv(0) = 0;
  X[1].deriv(1) = 1;


  double h1 = scale*0.1;
  double h2 = scale*0.5;

  SANS::SurrealS<2,double> eigval[2] = {1.0/(h1*h1), 1.0/(h2*h2)};
  //SANS::SurrealS<2,double> eigval[2] = {h1, h2};
  SANS::SurrealS<2,double> eigvec[4];

  SANS::SurrealS<2,double> r = sqrt(X[0]*X[0] + X[1]*X[1]) + 0.01;
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
  eigvec[2] = -sin(theta);
  eigvec[3] =  cos(theta);

  SANS::SurrealS<2,double> metS[3];
  eig2met<2,SANS::SurrealS<2,double>>(eigval,eigvec,metS);

  getmet_SurS2dbl<2>(metS,met,dmet);
}

// Boundary-layer mesh 
void anamet2D_3(void *ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){

  SANS::SurrealS<2,double> X[2];
  X[0] = crd[0];
  X[0].deriv(0) = 1;
  X[0].deriv(1) = 0;

  X[1] = crd[1];
  X[1].deriv(0) = 0;
  X[1].deriv(1) = 1;


  double hy_min = 0.001;
  double hy_max = 0.1;
  double hx = scale*0.5;
  SANS::SurrealS<2,double> hy = scale * (X[1] * hy_max + (1 - X[1] + 0.0) * hy_min);


  SANS::SurrealS<2,double> eigval[2] = {1.0/(hx*hx), 1.0/(hy*hy)};
  SANS::SurrealS<2,double> eigvec[4];

  eigvec[0] = 1.0;
  eigvec[1] = 0.0; 
  eigvec[2] = 0.0;
  eigvec[3] = 1.0;

  SANS::SurrealS<2,double> metS[3];
  eig2met<2,SANS::SurrealS<2,double>>(eigval,eigvec,metS);
  getmet_SurS2dbl<2>(metS,met,dmet);
}


// Boundary-layer mesh, slanted, wall = { x + y - 0.5 = 0 }
void anamet2D_4(void *ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){

  SANS::SurrealS<2,double> X;
  X = crd[0] + crd[1] - 0.5;
  X.deriv(0) = 1;
  X.deriv(1) = 1;

  //X[1] = crd[0] - crd[1];
  //X[1].deriv(0) = 1;
  //X[1].deriv(1) =-1;


  double hy_min = 0.001;
  double hy_max = 0.1;
  double hx = scale*0.5;
  SANS::SurrealS<2,double> hy = scale * (X * hy_max + (1 - X + 0.0) * hy_min);


  SANS::SurrealS<2,double> eigval[2] = {1.0/(hx*hx), 1.0/(hy*hy)};
  SANS::SurrealS<2,double> eigvec[4];

  eigvec[0] = sqrt(2);
  eigvec[1] =-sqrt(2); 
  eigvec[2] = sqrt(2);
  eigvec[3] = sqrt(2);

  SANS::SurrealS<2,double> metS[3];
  eig2met<2,SANS::SurrealS<2,double>>(eigval,eigvec,metS);
  getmet_SurS2dbl<2>(metS,met,dmet);
}




// circle BL
void anamet2D_5(void *ctx, const double*__restrict__ crd, double scale, int idif1, 
  double *met, double *dmet){
  const double pi = 3.141592653589793238462643383279502884;
  double x0 = 0.01;
  SANS::SurrealS<2,double> X[2];
  X[0] = crd[0] + x0;
  X[0].deriv(0) = 1;
  X[0].deriv(1) = 0;

  X[1] = crd[1] + x0;
  X[1].deriv(0) = 0;
  X[1].deriv(1) = 1;

  const double r0 = 0.5;

  SANS::SurrealS<2,double> r = sqrt(X[0]*X[0] + X[1]*X[1]) - r0;

  double hy_min = 0.001;
  double hy_max = 0.1;
  double hx = scale*0.5;
  SANS::SurrealS<2,double> hy = scale * (abs(r)*hy_max + (1 - abs(r))*hy_min);


  SANS::SurrealS<2,double> eigval[2] = {1.0/(hy*hy), 1.0/(hx*hx)};
  SANS::SurrealS<2,double> eigvec[4];

  SANS::SurrealS<2,double> y[2] = {X[0] / (abs(r) + 1.0e-6), X[1] / (abs(r) + 1.0e-6)};
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
  eigvec[2] = -sin(theta);
  eigvec[3] =  cos(theta);

  SANS::SurrealS<2,double> metS[3];
  eig2met<2,SANS::SurrealS<2,double>>(eigval,eigvec,metS);

  getmet_SurS2dbl<2>(metS,met,dmet);

  #ifndef NDEBUG
  if(dmet != NULL){
    //printf("debug r = %15.7e theta = %15.7e print met = %15.7e %15.7e %15.7e \n",
    //  r.value(), theta.value(), met[0],met[1],met[2]);
    for(int ii = 0; ii < 6; ii++){
      if(std::isnan(dmet[ii])){
        printf("## NAN METRIS IN ANAMET 5 ! coop = %f %f \n",crd[0],crd[1]);
        METRIS_THROW(GeomExcept());
      }
    }
  }
  #endif

  //for(int ii = 0; ii < 3; ii ++) met[ii] = metS[ii].value();

 //if(idif1 > 0){
 //  for(int ii = 0; ii < 2; ii++){
 //    for(int jj = 0; jj < 3; jj++){
 //      dmet[3*ii + jj] = metS[jj].deriv(ii);
 //    }
 //  }
 //}
}


// circle centered on 0 
void anamet2D_6(void *ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
  const double pi = 3.141592653589793238462643383279502884;
  SANS::SurrealS<2,double> X[2];
  X[0] = crd[0] + 0.01;
  X[0].deriv(0) = 1;
  X[0].deriv(1) = 0;

  X[1] = crd[1] + 0.01;
  X[1].deriv(0) = 0;
  X[1].deriv(1) = 1;


  double h1 = scale*0.1;
  double h2_ = scale*0.5;


  SANS::SurrealS<2,double> r = sqrt(X[0]*X[0] + X[1]*X[1]) + 0.01;
  SANS::SurrealS<2,double> y[2] = {X[0] / r, X[1] / r};

  SANS::SurrealS<2,double> h2 = h2_ * (1 + 10*r);

  SANS::SurrealS<2,double> eigval[2] = {1.0/(h1*h1), 1.0/(h2*h2)};
  //SANS::SurrealS<2,double> eigval[2] = {h1, h2};
  SANS::SurrealS<2,double> eigvec[4];

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
  eigvec[2] = -sin(theta);
  eigvec[3] =  cos(theta);

  SANS::SurrealS<2,double> metS[3];
  eig2met<2,SANS::SurrealS<2,double>>(eigval,eigvec,metS);

  getmet_SurS2dbl<2>(metS,met,dmet);
}


} // End namespace
