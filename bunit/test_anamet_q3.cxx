// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_anamet_q3
#include <boost/test/included/unit_test.hpp>

#include "aux_exceptions.hxx"
#include "msh_anamet.hxx"
#include "quality/quafun_connection.hxx"

#include <cmath>
#include <limits>

namespace Metris{
namespace{

double cubic_map(double s){
  return s*(1.5 - 0.25*s - 0.25*s*s);
}

} // namespace

BOOST_AUTO_TEST_CASE(q3_metric_registration_and_pullback){
  BOOST_REQUIRE(MAX_ANAMET2D_DEFINED >= 10);
  BOOST_CHECK(__ANAMET2D[9] == anamet2D_10);
  for(double scale : {1.,0.1,0.025}){
    for(int i = 0; i <= 64; i++){
      const double s = i/64.;
      const double d = 1.5 - 0.5*s - 0.75*s*s;
      const double point[2] = {cubic_map(s),0.37};
      double met[3], dmet[6], valueOnly[3];
      __ANAMET2D[9](nullptr,point,scale,1,met,dmet);
      anamet2D_10(nullptr,point,scale,0,valueOnly,nullptr);
      for(int m = 0; m < 3; m++) BOOST_CHECK_EQUAL(met[m],valueOnly[m]);
      BOOST_CHECK_GT(met[0],0.);
      BOOST_CHECK_GT(met[0]*met[2]-met[1]*met[1],0.);
      // DPsi^T M DPsi = C/scale^2 for the descending-diagonal grid.
      BOOST_CHECK_CLOSE_FRACTION(scale*scale*d*d*met[0],1.,1.e-12);
      BOOST_CHECK_CLOSE_FRACTION(scale*scale*d*met[1],0.5,1.e-12);
      BOOST_CHECK_CLOSE_FRACTION(scale*scale*met[2],1.,1.e-12);

      // Evaluate the exact cubic map's connection directly; no Q3 mesh
      // storage or degree-elevation capability is required for this check.
      const double jac[4] = {d,0.,0.,1.};
      const double hess[6] = {-0.5-1.5*s,0.,0.,0.,0.,0.};
      BOOST_CHECK_SMALL(connection_difference_norm2<2>(jac,hess,met,dmet),1.e-23);
    }
  }
}

BOOST_AUTO_TEST_CASE(q3_metric_physical_derivatives){
  for(double scale : {1.,0.1}){
    for(double s : {0.,0.07,0.31,0.68,0.93,1.}){
      double point[2] = {cubic_map(s),0.37};
      double met[3], dmet[6];
      anamet2D_10(nullptr,point,scale,1,met,dmet);
      for(int r = 0; r < 2; r++){
        const double initial = point[r], eps = 1.e-5;
        double plus[3], minus[3], plusHalf[3], minusHalf[3];
        point[r] = initial+eps;
        anamet2D_10(nullptr,point,scale,0,plus,nullptr);
        point[r] = initial-eps;
        anamet2D_10(nullptr,point,scale,0,minus,nullptr);
        point[r] = initial+eps/2;
        anamet2D_10(nullptr,point,scale,0,plusHalf,nullptr);
        point[r] = initial-eps/2;
        anamet2D_10(nullptr,point,scale,0,minusHalf,nullptr);
        point[r] = initial;
        for(int m = 0; m < 3; m++){
          const double coarse = (plus[m]-minus[m])/(2*eps);
          const double fine = (plusHalf[m]-minusHalf[m])/eps;
          const double numerical = (4*fine-coarse)/3;
          BOOST_CHECK_SMALL(dmet[3*r+m]-numerical,2.e-8*(1.+std::abs(numerical)));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(q3_metric_rejects_invalid_branch_and_scale){
  double point[2] = {0.5,0.5}, met[3];
  for(double scale : {0.,-1.,std::numeric_limits<double>::infinity()}){
    BOOST_CHECK_THROW(anamet2D_10(nullptr,point,scale,0,met,nullptr),MetrisExcept);
  }
  for(double x : {-3.,1.1,std::numeric_limits<double>::quiet_NaN()}){
    point[0] = x;
    BOOST_CHECK_THROW(anamet2D_10(nullptr,point,1.,0,met,nullptr),MetrisExcept);
  }
}

} // namespace Metris
