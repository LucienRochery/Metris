// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_quality_sizeshape_pointwise_derivatives

#include <boost/test/included/unit_test.hpp>

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "linalg/symidx.hxx"
#include "quality/quafun_sizeshape.hxx"

#include <cmath>

namespace Metris
{

namespace
{

void initialize_ideal_triangle(Mesh<MetricFieldFE> &msh,
                               MetrisParameters &parameters)
{
  msh.idim = 2;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &parameters;
  msh.set_npoin(3);
  msh.set_nface(1);

  msh.coord(0,0) = 0.0;
  msh.coord(0,1) = 0.0;
  msh.coord(1,0) = 1.0;
  msh.coord(1,1) = 0.0;
  msh.coord(2,0) = 0.5;
  msh.coord(2,1) = std::sqrt(3.0)/2.0;
  for(int inode = 0; inode < 3; inode++){
    msh.fac2poi(0,inode) = inode;
  }
}

void initialize_ideal_tetrahedron(Mesh<MetricFieldFE> &msh,
                                  MetrisParameters &parameters)
{
  msh.idim = 3;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &parameters;
  msh.set_npoin(4);
  msh.set_nelem(1);

  msh.coord(0,0) = 0.0;
  msh.coord(0,1) = 0.0;
  msh.coord(0,2) = 0.0;
  msh.coord(1,0) = -std::sqrt(3.0)/2.0;
  msh.coord(1,1) = 0.5;
  msh.coord(1,2) = 0.0;
  msh.coord(2,0) = -std::sqrt(3.0)/2.0;
  msh.coord(2,1) = -0.5;
  msh.coord(2,2) = 0.0;
  msh.coord(3,0) = -1.0/std::sqrt(3.0);
  msh.coord(3,1) = 0.0;
  msh.coord(3,2) = std::sqrt(2.0/3.0);
  for(int inode = 0; inode < 4; inode++){
    msh.tet2poi(0,inode) = inode;
  }
}

template<int gdim, int tdim>
double evaluate_sizeshape_value(Mesh<MetricFieldFE> &msh,
                                const int *nodes,
                                const double *bary,
                                const double *metric)
{
  return quafun_sizeshape<MetricFieldFE,gdim,tdim,double>(
      msh,AsDeg::P1,AsDeg::P1,nodes,bary,metric);
}

template<int gdim, int tdim>
double evaluate_sizeshape_derivatives(Mesh<MetricFieldFE> &msh,
                                      const int *nodes,
                                      const double *bary,
                                      const double *metric,
                                      int ivar,
                                      double *gradient,
                                      double *hessian)
{
  return d_quafun_sizeshape<MetricFieldFE,gdim,tdim,double>(
      msh,AsDeg::P1,AsDeg::P1,nodes,bary,metric,
      ivar,msh.getBasis(),DifVar::None,gradient,hessian);
}

template<int gdim, int tdim>
void check_sizeshape_finite_differences(Mesh<MetricFieldFE> &msh,
                                        const int *nodes,
                                        const double *bary,
                                        const double *metric,
                                        int ivar,
                                        double finite_difference_step,
                                        double gradient_tolerance,
                                        double hessian_tolerance)
{
  constexpr int nhessian = gdim*(gdim + 1)/2;
  double gradient[gdim];
  double hessian[nhessian];
  const double differentiated_value
      = evaluate_sizeshape_derivatives<gdim,tdim>(
          msh,nodes,bary,metric,ivar,gradient,hessian);
  const double value
      = evaluate_sizeshape_value<gdim,tdim>(msh,nodes,bary,metric);
  BOOST_CHECK_CLOSE_FRACTION(differentiated_value,value,2.0e-14);

  const int ipoin = nodes[ivar];
  for(int icomponent = 0; icomponent < gdim; icomponent++){
    const double coordinate = msh.coord(ipoin,icomponent);
    msh.coord(ipoin,icomponent) = coordinate + finite_difference_step;
    const double plus_value
        = evaluate_sizeshape_value<gdim,tdim>(msh,nodes,bary,metric);
    msh.coord(ipoin,icomponent) = coordinate - finite_difference_step;
    const double minus_value
        = evaluate_sizeshape_value<gdim,tdim>(msh,nodes,bary,metric);
    msh.coord(ipoin,icomponent) = coordinate;

    const double finite_difference_gradient
        = (plus_value - minus_value)/(2.0*finite_difference_step);
    const double gradient_error
        = std::abs(gradient[icomponent] - finite_difference_gradient);
    BOOST_CHECK_SMALL(
        gradient_error,
        gradient_tolerance*(1.0 + std::abs(finite_difference_gradient)));
  }

  for(int jcomponent = 0; jcomponent < gdim; jcomponent++){
    const double coordinate = msh.coord(ipoin,jcomponent);
    double plus_gradient[gdim];
    double minus_gradient[gdim];

    msh.coord(ipoin,jcomponent) = coordinate + finite_difference_step;
    evaluate_sizeshape_derivatives<gdim,tdim>(
        msh,nodes,bary,metric,ivar,plus_gradient,NULL);
    msh.coord(ipoin,jcomponent) = coordinate - finite_difference_step;
    evaluate_sizeshape_derivatives<gdim,tdim>(
        msh,nodes,bary,metric,ivar,minus_gradient,NULL);
    msh.coord(ipoin,jcomponent) = coordinate;

    for(int icomponent = 0; icomponent <= jcomponent; icomponent++){
      const double finite_difference_hessian
          = (plus_gradient[icomponent] - minus_gradient[icomponent])
           /(2.0*finite_difference_step);
      const double hessian_error
          = std::abs(hessian[sym2idx(icomponent,jcomponent)]
                     - finite_difference_hessian);
      BOOST_CHECK_SMALL(
          hessian_error,
          hessian_tolerance*(1.0 + std::abs(finite_difference_hessian)));
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(test_anisotropic_triangle_finite_differences)
{
  MetrisParameters parameters;
  parameters.iverb = 0;

  Mesh<MetricFieldFE> msh;
  initialize_ideal_triangle(msh,parameters);
  const double barycenter[3] = {1.0/3.0,1.0/3.0,1.0/3.0};
  const double metric[3] = {4.0,0.3,1.2};
  const double powers[3] = {1.0,2.0,1.5};

  for(int ipower = 0; ipower < 3; ipower++){
    parameters.objective_p = powers[ipower];
    check_sizeshape_finite_differences<2,2>(
        msh,msh.fac2poi[0],barycenter,metric,1,
        2.0e-6,2.0e-8,2.0e-7);
  }
}

BOOST_AUTO_TEST_CASE(test_anisotropic_tetrahedron_finite_differences)
{
  MetrisParameters parameters;
  parameters.iverb = 0;

  Mesh<MetricFieldFE> msh;
  initialize_ideal_tetrahedron(msh,parameters);
  const double barycenter[4] = {0.25,0.25,0.25,0.25};
  const double metric[6] = {3.0,0.2,1.4,-0.1,0.15,0.8};
  const double powers[3] = {1.0,2.0,1.5};

  for(int ipower = 0; ipower < 3; ipower++){
    parameters.objective_p = powers[ipower];
    check_sizeshape_finite_differences<3,3>(
        msh,msh.tet2poi[0],barycenter,metric,2,
        2.0e-6,5.0e-8,8.0e-7);
  }
}

BOOST_AUTO_TEST_CASE(test_exact_and_near_ideal_derivative_limits)
{
  MetrisParameters parameters;
  parameters.iverb = 0;

  Mesh<MetricFieldFE> msh;
  initialize_ideal_triangle(msh,parameters);
  const double barycenter[3] = {1.0/3.0,1.0/3.0,1.0/3.0};
  const double identity_metric[3] = {1.0,0.0,1.0};
  constexpr int gdim = 2;
  constexpr int nhessian = 3;

  parameters.objective_p = 1.0;
  check_sizeshape_finite_differences<2,2>(
      msh,msh.fac2poi[0],barycenter,identity_metric,1,
      2.0e-6,2.0e-9,2.0e-7);

  const double powers_above_one[2] = {2.0,1.5};
  for(int ipower = 0; ipower < 2; ipower++){
    parameters.objective_p = powers_above_one[ipower];
    double gradient[gdim];
    double hessian[nhessian];
    const double value = evaluate_sizeshape_derivatives<2,2>(
        msh,msh.fac2poi[0],barycenter,identity_metric,1,
        gradient,hessian);
    BOOST_CHECK_EQUAL(value,0.0);
    for(int icomponent = 0; icomponent < gdim; icomponent++){
      BOOST_CHECK_EQUAL(gradient[icomponent],0.0);
    }
    for(int ihessian = 0; ihessian < nhessian; ihessian++){
      BOOST_CHECK_EQUAL(hessian[ihessian],0.0);
    }
    check_sizeshape_finite_differences<2,2>(
        msh,msh.fac2poi[0],barycenter,identity_metric,1,
        5.0e-7,2.0e-8,1.0e-5);
  }

  msh.coord(1,0) += 1.0e-3;
  const double powers[3] = {1.0,2.0,1.5};
  for(int ipower = 0; ipower < 3; ipower++){
    parameters.objective_p = powers[ipower];
    check_sizeshape_finite_differences<2,2>(
        msh,msh.fac2poi[0],barycenter,identity_metric,1,
        1.0e-7,2.0e-8,5.0e-6);
  }
}

BOOST_AUTO_TEST_CASE(test_derivatives_are_independent_of_classical_opt_power)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_p = 1.5;

  Mesh<MetricFieldFE> msh;
  initialize_ideal_triangle(msh,parameters);
  const double barycenter[3] = {1.0/3.0,1.0/3.0,1.0/3.0};
  const double metric[3] = {4.0,0.3,1.2};
  constexpr int gdim = 2;
  constexpr int nhessian = 3;

  double direct_gradient[gdim];
  double direct_hessian[nhessian];
  parameters.opt_power = 1;
  const double direct_value = evaluate_sizeshape_derivatives<2,2>(
      msh,msh.fac2poi[0],barycenter,metric,1,
      direct_gradient,direct_hessian);

  double inverse_setting_gradient[gdim];
  double inverse_setting_hessian[nhessian];
  parameters.opt_power = -1;
  const double inverse_setting_value
      = evaluate_sizeshape_derivatives<2,2>(
          msh,msh.fac2poi[0],barycenter,metric,1,
          inverse_setting_gradient,inverse_setting_hessian);

  BOOST_CHECK_EQUAL(inverse_setting_value,direct_value);
  for(int icomponent = 0; icomponent < gdim; icomponent++){
    BOOST_CHECK_EQUAL(inverse_setting_gradient[icomponent],
                      direct_gradient[icomponent]);
  }
  for(int ihessian = 0; ihessian < nhessian; ihessian++){
    BOOST_CHECK_EQUAL(inverse_setting_hessian[ihessian],
                      direct_hessian[ihessian]);
  }
}

} // namespace Metris
