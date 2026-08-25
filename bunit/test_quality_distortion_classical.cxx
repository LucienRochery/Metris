// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_quality_distortion_classical

#include <boost/test/included/unit_test.hpp>

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "linalg/symidx.hxx"
#include "low_geo/measure.hxx"
#include "metris_options.hxx"
#include "quality/low_metqua.hxx"
#include "quality/low_metqua_d.hxx"
#include "quality/quafun_distortion.hxx"

#include <algorithm>
#include <cmath>

namespace Metris
{

namespace
{

void check_classical_distortion_derivatives(Mesh<MetricFieldFE> &msh,
                                            int power,
                                            double expected_value)
{
  constexpr int gdim = 2;
  constexpr int tdim = 2;
  constexpr int nhess = 3;

  msh.param->opt_power = power;

  const int *nodes = msh.fac2poi[0];
  const double bary[3] = {0.2,0.3,0.5};
  const double metric[3] = {1.7,0.25,0.9};
  const double fd_step = 1.e-6;

  const double value
      = quafun_distortion<MetricFieldFE,gdim,tdim,double>(
          msh,AsDeg::P1,AsDeg::P1,nodes,bary,metric);
  BOOST_CHECK_CLOSE_FRACTION(value,expected_value,2.e-15);

  for(int ivar = 0; ivar < tdim + 1; ivar++){
    double gradient[gdim];
    double hessian[nhess];
    const double differentiated_value
        = d_quafun_distortion<MetricFieldFE,gdim,tdim,double>(
            msh,AsDeg::P1,AsDeg::P1,nodes,bary,metric,
            ivar,msh.getBasis(),DifVar::None,gradient,hessian);
    BOOST_CHECK_CLOSE_FRACTION(differentiated_value,value,2.e-15);

    const int ipoin = nodes[ivar];
    for(int j = 0; j < gdim; j++){
      const double coordinate = msh.coord(ipoin,j);
      double gradient_plus[gdim];
      double gradient_minus[gdim];

      msh.coord(ipoin,j) = coordinate + fd_step;
      const double value_plus
          = quafun_distortion<MetricFieldFE,gdim,tdim,double>(
              msh,AsDeg::P1,AsDeg::P1,nodes,bary,metric);
      d_quafun_distortion<MetricFieldFE,gdim,tdim,double>(
          msh,AsDeg::P1,AsDeg::P1,nodes,bary,metric,
          ivar,msh.getBasis(),DifVar::None,gradient_plus,NULL);

      msh.coord(ipoin,j) = coordinate - fd_step;
      const double value_minus
          = quafun_distortion<MetricFieldFE,gdim,tdim,double>(
              msh,AsDeg::P1,AsDeg::P1,nodes,bary,metric);
      d_quafun_distortion<MetricFieldFE,gdim,tdim,double>(
          msh,AsDeg::P1,AsDeg::P1,nodes,bary,metric,
          ivar,msh.getBasis(),DifVar::None,gradient_minus,NULL);
      msh.coord(ipoin,j) = coordinate;

      const double gradient_fd
          = (value_plus - value_minus)/(2.0*fd_step);
      BOOST_CHECK_SMALL(
          gradient[j] - gradient_fd,
          2.e-5*std::max(1.0,std::abs(gradient_fd)));

      for(int i = 0; i <= j; i++){
        const double hessian_fd
            = (gradient_plus[i] - gradient_minus[i])/(2.0*fd_step);
        BOOST_CHECK_SMALL(
            hessian[sym2idx(i,j)] - hessian_fd,
            5.e-4*std::max(1.0,std::abs(hessian_fd)));
      }
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(test_classical_distortion_both_power_signs)
{
  MetrisParameters parameters;
  parameters.iverb = 0;

  Mesh<MetricFieldFE> msh;
  msh.idim = 2;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &parameters;
  msh.set_npoin(3);
  msh.set_nface(1);

  msh.coord(0,0) = 0.0;
  msh.coord(0,1) = 0.0;
  msh.coord(1,0) = 1.2;
  msh.coord(1,1) = 0.1;
  msh.coord(2,0) = 0.2;
  msh.coord(2,1) = 0.7;
  msh.fac2poi(0,0) = 0;
  msh.fac2poi(0,1) = 1;
  msh.fac2poi(0,2) = 2;

  constexpr double direct_value = 1.9620361990797364;
  constexpr double inverse_value = 0.50967459237960799;

  parameters.objective_p = 1.0;
  check_classical_distortion_derivatives(msh,1,direct_value);
  check_classical_distortion_derivatives(msh,-1,inverse_value);

  // The shared objective exponent must not affect Classical Distortion.
  parameters.objective_p = 3.25;
  check_classical_distortion_derivatives(msh,1,direct_value);
  check_classical_distortion_derivatives(msh,-1,inverse_value);

  BOOST_CHECK_CLOSE_FRACTION(direct_value*inverse_value,1.0,2.e-15);
}

BOOST_AUTO_TEST_CASE(test_classical_integrated_distortion_uses_log_metric)
{
  constexpr int gdim = 2;
  constexpr int tdim = 2;
  constexpr int nmetric = 3;

  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.opt_power = 1;
  parameters.opt_pnorm = 1;

  Mesh<MetricFieldFE> msh;
  msh.idim = gdim;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &parameters;
  msh.set_npoin(3);
  msh.set_nface(1);
  msh.met.forceBasisFlag(FEBasis::Lagrange);
  msh.met.forceSpaceFlag(MetSpace::Exp);

  msh.coord(0,0) = 0.0; msh.coord(0,1) = 0.0;
  msh.coord(1,0) = 1.2; msh.coord(1,1) = 0.1;
  msh.coord(2,0) = 0.2; msh.coord(2,1) = 0.7;
  for(int inode = 0; inode < 3; inode++) msh.fac2poi(0,inode) = inode;

  const double nodal_metrics[3][nmetric] = {
      {1.0,0.0,1.0},
      {64.0,0.0,1.0},
      {1.0,0.0,4.0}};
  for(int inode = 0; inode < 3; inode++){
    for(int imetric = 0; imetric < nmetric; imetric++){
      msh.met(inode,imetric) = nodal_metrics[inode][imetric];
    }
  }

  const double barycenter[3] = {1.0/3.0,1.0/3.0,1.0/3.0};
  double log_interpolated_metric[nmetric];
  msh.met.getMetBary(AsDeg::P1,DifVar::None,MetSpace::Exp,
                     msh.fac2poi[0],tdim,barycenter,
                     log_interpolated_metric,NULL);
  BOOST_CHECK_CLOSE_FRACTION(log_interpolated_metric[0],4.0,2.e-14);
  BOOST_CHECK_SMALL(log_interpolated_metric[1],2.e-14);
  BOOST_CHECK_CLOSE_FRACTION(
      log_interpolated_metric[2],std::pow(4.0,1.0/3.0),2.e-14);

  double expected
      = std::abs(
          quafun_distortion<MetricFieldFE,gdim,tdim,double>(
              msh,AsDeg::P1,AsDeg::P1,msh.fac2poi[0],barycenter,
              log_interpolated_metric)
          - 1.0);
#ifdef TESTQUALITYALGO
  double measure;
  isvalideltP1<gdim,tdim>(msh,0,NULL,&measure);
  expected *= measure;
#ifdef INTQUALINRIEMSPACE
  expected *= std::sqrt(
      log_interpolated_metric[0]*log_interpolated_metric[2]
      - log_interpolated_metric[1]*log_interpolated_metric[1]);
#endif
#endif

  const double value
      = metqua<MetricFieldFE,gdim,tdim,QuaFun::Distortion,double>(
          msh,AsDeg::P1,AsDeg::P1,0,1.0);
  const double differentiated_value
      = d_metqua<MetricFieldFE,gdim,tdim,QuaFun::Distortion,double>(
          msh,AsDeg::P1,AsDeg::P1,0,-1,msh.getBasis(),DifVar::None,
          NULL,NULL,1.0);
  BOOST_CHECK_CLOSE_FRACTION(value,expected,2.e-13);
  BOOST_CHECK_CLOSE_FRACTION(differentiated_value,expected,2.e-13);

  const double arithmetic_metric[nmetric] = {22.0,0.0,2.0};
  double arithmetic_expected
      = std::abs(
          quafun_distortion<MetricFieldFE,gdim,tdim,double>(
              msh,AsDeg::P1,AsDeg::P1,msh.fac2poi[0],barycenter,
              arithmetic_metric)
          - 1.0);
#ifdef TESTQUALITYALGO
  double arithmetic_measure;
  isvalideltP1<gdim,tdim>(msh,0,NULL,&arithmetic_measure);
  arithmetic_expected *= arithmetic_measure;
#ifdef INTQUALINRIEMSPACE
  arithmetic_expected *= std::sqrt(
      arithmetic_metric[0]*arithmetic_metric[2]
      - arithmetic_metric[1]*arithmetic_metric[1]);
#endif
#endif
  BOOST_CHECK_GT(
      std::abs(expected - arithmetic_expected),
      1.e-3*(1.0 + std::abs(expected)));
}

} // namespace Metris
