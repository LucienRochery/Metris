//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#include "msh_anamet.hxx"

#include "linalg/eigen.hxx"
#include "linalg/utils.hxx"

#include "../libs/SANS/Surreal/SurrealS.h"
#include "aux_exceptions.hxx"
#include <cmath>


#include "fmt/format.h"


namespace Metris{

void anamet2D_1([[maybe_unused]] const AnaMetCtx* ctx,
                [[maybe_unused]] const double*__restrict__ crd,
                double scale, int idif1, double *met, double *dmet){
  // Not too coarse at the scale of 1
  double h0 = 0.05;

  met[0] = 1/(h0*h0*scale*scale);
  met[1] = 0.0;
  met[2] = 1/(h0*h0*scale*scale);

  if(idif1 > 0){
    for(int ii = 0; ii < 2; ii++){
      for(int jj = 0; jj < 3; jj++){
        dmet[3*ii + jj] = 0;
      }
    }
  }
}

// circle
void anamet2D_2([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
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

  //if(ctx != NULL){
  //  fmt::print("anamet2D_2 debug r = {:.6f} theta = {:.6f}\n", r.value(), theta.value());
  //}

  // eig2met is in R^T D R format. Worst case we are using -theta.
  eigvec[0] =  cos(theta);
  eigvec[1] =  sin(theta);
  eigvec[2] = -sin(theta);
  eigvec[3] =  cos(theta);

  SANS::SurrealS<2,double> metS[3];
  eig2met<2,SANS::SurrealS<2,double>>(eigval,eigvec,metS);
  getmet_SurS2dbl<2>(metS,met,idif1 > 0 ? dmet : NULL);
}

// Boundary-layer along x centered at 0.5
void anamet2D_3([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){

  SANS::SurrealS<2,double> X[2];
  X[0] = std::abs(crd[0] - 0.5);
  X[0].deriv(0) = crd[0] >= 0.5 ? 1 : -1;
  X[0].deriv(1) = 0;

  X[1] = crd[1];
  X[1].deriv(0) = 0;
  X[1].deriv(1) = 1;


  double hx_min = 0.001;
  double hx_max = 0.1;
  SANS::SurrealS<2,double> hx = scale*(X[0]*(hx_max - hx_min) + hx_min);
  double hy = scale*0.1;


  SANS::SurrealS<2,double> eigval[2] = {1.0/(hx*hx), 1.0/(hy*hy)};
  SANS::SurrealS<2,double> eigvec[4];

  eigvec[0] = 1.0;
  eigvec[1] = 0.0;
  eigvec[2] = 0.0;
  eigvec[3] = 1.0;

  SANS::SurrealS<2,double> metS[3];
  eig2met<2,SANS::SurrealS<2,double>>(eigval,eigvec,metS);
  getmet_SurS2dbl<2>(metS,met,idif1 > 0 ? dmet : NULL);
}


// Boundary-layer mesh, slanted, wall = { x + y - 0.5 = 0 }
void anamet2D_4([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){

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
  SANS::SurrealS<2,double> hy = scale*(X*hy_max + (1 - X + 0.0)*hy_min);


  SANS::SurrealS<2,double> eigval[2] = {1.0/(hx*hx), 1.0/(hy*hy)};
  SANS::SurrealS<2,double> eigvec[4];

  eigvec[0] = sqrt(2);
  eigvec[1] =-sqrt(2);
  eigvec[2] = sqrt(2);
  eigvec[3] = sqrt(2);

  SANS::SurrealS<2,double> metS[3];
  eig2met<2,SANS::SurrealS<2,double>>(eigval,eigvec,metS);
  getmet_SurS2dbl<2>(metS,met,idif1 > 0 ? dmet : NULL);
}




// circle BL
void anamet2D_5([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1,
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
  SANS::SurrealS<2,double> hy = scale*(std::abs(r)*hy_max + (1 - std::abs(r))*hy_min);


  SANS::SurrealS<2,double> eigval[2] = {1.0/(hy*hy), 1.0/(hx*hx)};
  SANS::SurrealS<2,double> eigvec[4];

  SANS::SurrealS<2,double> y[2] = {X[0] / (std::abs(r) + 1.0e-6), X[1] / (std::abs(r) + 1.0e-6)};
  SANS::SurrealS<2,double> theta;
  if(y[0].value() > 0){
    theta = atan(y[1]/y[0]);
  }else{
    theta = pi + atan(y[1]/y[0]);
  }

  //if(ctx != NULL){
  //  fmt::print("anamet2D_2 debug r = {:.6f} theta = {:.6f}\n", r.value(), theta.value());
  //}

  // eig2met is in R^T D R format. Worst case we are using -theta.
  eigvec[0] =  cos(theta);
  eigvec[1] =  sin(theta);
  eigvec[2] = -sin(theta);
  eigvec[3] =  cos(theta);

  SANS::SurrealS<2,double> metS[3];
  eig2met<2,SANS::SurrealS<2,double>>(eigval,eigvec,metS);

  getmet_SurS2dbl<2>(metS,met,idif1 > 0 ? dmet : NULL);

  #ifndef NDEBUG
  if(idif1 > 0){
    //fmt::print("debug r = {:15.7e} theta = {:15.7e} print met = {:15.7e} {:15.7e} {:15.7e} \n",
    //  r.value(), theta.value(), met[0],met[1],met[2]);
    for(int ii = 0; ii < 6; ii++){
      METRIS_ASSERT_MSG(!std::isnan(met[ii]),
        "NAN METRIS IN ANAMET 5 ! coop = {} {} \n",crd[0],crd[1]);
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
void anamet2D_6([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){
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

  SANS::SurrealS<2,double> h2 = h2_*(1 + 10*r);

  SANS::SurrealS<2,double> eigval[2] = {1.0/(h1*h1), 1.0/(h2*h2)};
  //SANS::SurrealS<2,double> eigval[2] = {h1, h2};
  SANS::SurrealS<2,double> eigvec[4];

  SANS::SurrealS<2,double> theta;
  if(y[0].value() > 0){
    theta = atan(y[1]/y[0]);
  }else{
    theta = pi + atan(y[1]/y[0]);
  }

  //if(ctx != NULL){
  //  fmt::print("anamet2D_2 debug r = {:.6f} theta = {:.6f}\n", r.value(), theta.value());
  //}

  // eig2met is in R^T D R format. Worst case we are using -theta.
  eigvec[0] =  cos(theta);
  eigvec[1] =  sin(theta);
  eigvec[2] = -sin(theta);
  eigvec[3] =  cos(theta);

  SANS::SurrealS<2,double> metS[3];
  eig2met<2,SANS::SurrealS<2,double>>(eigval,eigvec,metS);

  getmet_SurS2dbl<2>(metS,met,idif1 > 0 ? dmet : NULL);
}

// optimal metric for corner singularity at 0
// M = C * r^{-2k} * I
// with k = 1 - (alpha + 1)/(p+2) ---- alpha = singularity strength, p = solution order
void anamet2D_7([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){

  const double alpha = 0.2;
  const int p = 1;

  const double k = 1. - (alpha + 1.)/(p + 2.);

  SANS::SurrealS<2,double> X[2];
  X[0] = crd[0];
  X[0].deriv(0) = 1;
  X[0].deriv(1) = 0;

  X[1] = crd[1];
  X[1].deriv(0) = 0;
  X[1].deriv(1) = 1;

  SANS::SurrealS<2,double> r2 = X[0]*X[0] + X[1]*X[1] + 1e-12; // to avoid blow up at the origin

  SANS::SurrealS<2,double> metricVal = 1./(scale * scale) * pow(r2,-k);

  SANS::SurrealS<2,double> metS[3];
  metS[0] = metricVal; // M11
  metS[1] = 0;         // M12
  metS[2] = metricVal; // M22

  getmet_SurS2dbl<2>(metS,met,idif1 > 0 ? dmet : NULL);
}

// optimal metric for boundary layer at x = 0

// x: normal to BL
// y: tangent to BL

// function: u(x,y) = exp(-x/epsilon) + beta/(p+1)! * y^(p+1)

// metric: diag(1/h1^2,1/h2^2)
// h1 = exp(k1*x)
// h2 = AR * h1
void anamet2D_8([[maybe_unused]] const AnaMetCtx* ctx, const double*__restrict__ crd, double scale, int idif1, double *met, double *dmet){

  const double epsilon = 0.1;
  const double beta = 1.;
  const int p = 1;

  const double delta = epsilon * (p + 3./2.) * (1. - 1./(4.*p*p + 12.*p + 9.));
  const double k1 = 1./delta;

  const double kR = -1./(epsilon*(p+1.));
  const double AR0 = 1./(pow(beta,1./(p+1.)) * epsilon);

  SANS::SurrealS<2,double> X[2];
  X[0] = crd[0];
  X[0].deriv(0) = 1;
  X[0].deriv(1) = 0;

  X[1] = crd[1];
  X[1].deriv(0) = 0;
  X[1].deriv(1) = 1;

  SANS::SurrealS<2,double> h1 = scale * exp(k1*X[0]);
  SANS::SurrealS<2,double> AR = AR0 * exp(kR*X[0]);
  SANS::SurrealS<2,double> h2 = AR * h1;

  SANS::SurrealS<2,double> metS[3];
  metS[0] = 1./(h1*h1); // M11
  metS[1] = 0;          // M12
  metS[2] = 1./(h2*h2); // M22

  getmet_SurS2dbl<2>(metS,met,idif1 > 0 ? dmet : NULL);
}

} // End namespace
