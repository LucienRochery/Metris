// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_metric_interpolation_log_euclidean

#include <boost/test/included/unit_test.hpp>

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "ho_constants.hxx"
#include "linalg/explogmet.hxx"

#include <algorithm>
#include <array>
#include <cmath>

namespace Metris
{

namespace
{

constexpr int gdim = 2;
constexpr int nnmet = 3;

void initialize_triangle(Mesh<MetricFieldFE> &msh,
                         MetrisParameters &parameters)
{
  msh.idim = gdim;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &parameters;
  msh.set_npoin(3);
  msh.set_nentt(2,1);
  msh.met.forceBasisFlag(FEBasis::Lagrange);
  msh.met.forceSpaceFlag(MetSpace::Exp);

  msh.coord(0,0) = 0.0; msh.coord(0,1) = 0.0;
  msh.coord(1,0) = 1.0; msh.coord(1,1) = 0.0;
  msh.coord(2,0) = 0.0; msh.coord(2,1) = 1.0;
  for(int inode = 0; inode < 3; inode++){
    msh.fac2poi(0,inode) = inode;
  }

  // SPD matrices with different eigenvectors. Their log-Euclidean
  // interpolation is detectably different from componentwise averaging.
  constexpr double metric[3][nnmet] = {
      {4.0, 0.0, 1.0},
      {2.5, 1.5, 2.5},
      {1.2,-0.4, 3.1}};
  for(int inode = 0; inode < 3; inode++){
    for(int imet = 0; imet < nnmet; imet++){
      msh.met(inode,imet) = metric[inode][imet];
    }
  }
}

std::array<double,nnmet>
expected_log_metric(const Mesh<MetricFieldFE> &msh,
                    const double *bary)
{
  std::array<double,nnmet> expected = {};
  for(int inode = 0; inode < 3; inode++){
    double nodal_log[nnmet];
    for(int imet = 0; imet < nnmet; imet++){
      nodal_log[imet] = msh.met(inode,imet);
    }
    if(msh.met.getSpace() == MetSpace::Exp){
      getlogmet_inp<gdim,double>(nodal_log);
    }
    for(int imet = 0; imet < nnmet; imet++){
      expected[imet] += bary[inode]*nodal_log[imet];
    }
  }
  return expected;
}

void check_close(const double *actual,
                 const double *expected,
                 double tolerance = 2.e-13)
{
  for(int imet = 0; imet < nnmet; imet++){
    BOOST_CHECK_SMALL(
        actual[imet]-expected[imet],
        tolerance*std::max(1.0,std::abs(expected[imet])));
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(fe_metric_interpolation_is_log_euclidean_in_any_storage)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  Mesh<MetricFieldFE> msh;
  initialize_triangle(msh,parameters);

  const double bary[3] = {0.2,0.3,0.5};
  const int *nodes = msh.fac2poi[0];
  const auto expected_log = expected_log_metric(msh,bary);
  auto expected_exp = expected_log;
  getexpmet_inp<gdim,double>(expected_exp.data());

  double metric_log[nnmet],metric_exp[nnmet];
  msh.met.getMetBary(AsDeg::P1,DifVar::None,MetSpace::Log,
                     nodes,2,bary,metric_log,NULL);
  msh.met.getMetBary(AsDeg::P1,DifVar::None,MetSpace::Exp,
                     nodes,2,bary,metric_exp,NULL);
  check_close(metric_log,expected_log.data());
  check_close(metric_exp,expected_exp.data());

  double arithmetic[nnmet] = {};
  for(int inode = 0; inode < 3; inode++){
    for(int imet = 0; imet < nnmet; imet++){
      arithmetic[imet] += bary[inode]*msh.met(inode,imet);
    }
  }
  double difference = 0.0;
  for(int imet = 0; imet < nnmet; imet++){
    difference = std::max(
        difference,std::abs(metric_exp[imet]-arithmetic[imet]));
  }
  BOOST_CHECK_GT(difference,1.e-2);

  msh.met.setSpace(MetSpace::Log);
  msh.met.getMetBary(AsDeg::P1,DifVar::None,MetSpace::Log,
                     nodes,2,bary,metric_log,NULL);
  msh.met.getMetBary(AsDeg::P1,DifVar::None,MetSpace::Exp,
                     nodes,2,bary,metric_exp,NULL);
  check_close(metric_log,expected_log.data());
  check_close(metric_exp,expected_exp.data());
}

BOOST_AUTO_TEST_CASE(fe_metric_physical_derivatives_follow_log_then_exp)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  Mesh<MetricFieldFE> msh;
  initialize_triangle(msh,parameters);

  const double bary[3] = {0.2,0.3,0.5};
  const int *nodes = msh.fac2poi[0];
  double metric[nnmet],derivative[gdim*nnmet];
  msh.met.getMetBary(AsDeg::P1,DifVar::Phys,MetSpace::Exp,
                     nodes,2,bary,metric,derivative);

  constexpr double step = 1.e-7;
  for(int idim = 0; idim < gdim; idim++){
    double bary_plus[3] = {bary[0],bary[1],bary[2]};
    double bary_minus[3] = {bary[0],bary[1],bary[2]};
    bary_plus[0] -= step;
    bary_minus[0] += step;
    bary_plus[idim + 1] += step;
    bary_minus[idim + 1] -= step;

    double metric_plus[nnmet],metric_minus[nnmet];
    msh.met.getMetBary(AsDeg::P1,DifVar::None,MetSpace::Exp,
                       nodes,2,bary_plus,metric_plus,NULL);
    msh.met.getMetBary(AsDeg::P1,DifVar::None,MetSpace::Exp,
                       nodes,2,bary_minus,metric_minus,NULL);
    for(int imet = 0; imet < nnmet; imet++){
      const double finite_difference
          = (metric_plus[imet]-metric_minus[imet])/(2.0*step);
      BOOST_CHECK_SMALL(
          derivative[idim*nnmet + imet]-finite_difference,
          2.e-7*std::max(1.0,std::abs(finite_difference)));
    }
  }
}

BOOST_AUTO_TEST_CASE(degree_elevation_metric_sampling_is_log_euclidean)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  Mesh<MetricFieldFE> msh;
  initialize_triangle(msh,parameters);

  constexpr int target_degree = 2;
  constexpr int target_nodes = getnnod2(target_degree);
  double sampled[target_nodes*nnmet];
  msh.met.getMetNodes<gdim,2,1,target_degree>(0,sampled);

  constexpr auto ordering = ORDELT(2);
  for(int inode = 0; inode < target_nodes; inode++){
    double bary[3];
    for(int ibary = 0; ibary < 3; ibary++){
      bary[ibary] = ordering[target_degree][inode][ibary]
                  / static_cast<double>(target_degree);
    }
    const auto expected_log = expected_log_metric(msh,bary);
    auto expected_exp = expected_log;
    getexpmet_inp<gdim,double>(expected_exp.data());
    check_close(&sampled[inode*nnmet],expected_exp.data());
  }

  msh.met.setSpace(MetSpace::Log);
  msh.met.getMetNodes<gdim,2,1,target_degree>(0,sampled);
  for(int inode = 0; inode < target_nodes; inode++){
    double bary[3];
    for(int ibary = 0; ibary < 3; ibary++){
      bary[ibary] = ordering[target_degree][inode][ibary]
                  / static_cast<double>(target_degree);
    }
    const auto expected_log = expected_log_metric(msh,bary);
    check_close(&sampled[inode*nnmet],expected_log.data());
  }
}

} // namespace Metris
