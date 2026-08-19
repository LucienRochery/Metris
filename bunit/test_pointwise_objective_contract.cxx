// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_pointwise_objective_contract

#include <boost/test/included/unit_test.hpp>

#include "Mesh/Mesh.hxx"
#include "metris_options.hxx"
#include "quality/pointwise_objective.hxx"

namespace Metris
{

namespace
{

template<QuaFun iquaf>
void check_current_objective_adapter(Mesh<MetricFieldFE> &msh)
{
  constexpr int gdim = 2;
  constexpr int tdim = 2;
  constexpr int nhessian = 3;

  const int *nodes = msh.fac2poi[0];
  const double bary[3] = {0.2,0.3,0.5};
  const double metric[3] = {1.7,0.25,0.9};
  const int ivar = 1;

  const double direct_value
      = get_quafun_xi<MetricFieldFE,gdim,tdim,iquaf,double>()(
          msh,AsDeg::P1,AsDeg::P1,nodes,bary,metric);
  const PointwiseObjectiveResult<gdim,double> value_result
      = evaluate_pointwise_objective_value<
          MetricFieldFE,gdim,tdim,iquaf,double>(
              msh,AsDeg::P1,AsDeg::P1,nodes,bary,metric);

  BOOST_CHECK_EQUAL(value_result.psi,direct_value);
  BOOST_CHECK(!value_result.has_gradient());
  BOOST_CHECK(!value_result.has_hessian());

  double direct_gradient[gdim];
  double direct_hessian[nhessian];
  const double direct_differentiated_value
      = get_d_quafun_xi<MetricFieldFE,gdim,tdim,iquaf,double>()(
          msh,AsDeg::P1,AsDeg::P1,nodes,bary,metric,
          ivar,msh.getBasis(),DifVar::None,
          direct_gradient,direct_hessian);

  const PointwiseObjectiveResult<gdim,double> gradient_result
      = evaluate_pointwise_objective_derivatives<
          MetricFieldFE,gdim,tdim,iquaf,double>(
              msh,AsDeg::P1,AsDeg::P1,nodes,bary,metric,
              ivar,msh.getBasis(),DifVar::None,
              PointwiseDerivativeOrder::Gradient);
  BOOST_CHECK_EQUAL(gradient_result.psi,direct_differentiated_value);
  BOOST_CHECK(gradient_result.has_gradient());
  BOOST_CHECK(!gradient_result.has_hessian());
  for(int i = 0; i < gdim; i++){
    BOOST_CHECK_EQUAL(gradient_result.gradient[i],direct_gradient[i]);
  }

  const PointwiseObjectiveResult<gdim,double> hessian_result
      = evaluate_pointwise_objective_derivatives<
          MetricFieldFE,gdim,tdim,iquaf,double>(
              msh,AsDeg::P1,AsDeg::P1,nodes,bary,metric,
              ivar,msh.getBasis(),DifVar::None,
              PointwiseDerivativeOrder::Hessian);
  BOOST_CHECK_EQUAL(hessian_result.psi,direct_differentiated_value);
  BOOST_CHECK(hessian_result.has_gradient());
  BOOST_CHECK(hessian_result.has_hessian());
  for(int i = 0; i < gdim; i++){
    BOOST_CHECK_EQUAL(hessian_result.gradient[i],direct_gradient[i]);
  }
  for(int i = 0; i < nhessian; i++){
    BOOST_CHECK_EQUAL(hessian_result.hessian[i],direct_hessian[i]);
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(test_result_state_and_current_objective_adapters)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.opt_power = 1;
  parameters.opt_pnorm = 1;
  parameters.objective_p = 1.5;
  parameters.step_distance_regularization = 1.e-8;

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

  check_current_objective_adapter<QuaFun::SizeShape>(msh);
  check_current_objective_adapter<QuaFun::StepDistance>(msh);
}

BOOST_AUTO_TEST_CASE(test_differentiated_adapter_rejects_value_order)
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
  msh.coord(1,0) = 1.0;
  msh.coord(1,1) = 0.0;
  msh.coord(2,0) = 0.0;
  msh.coord(2,1) = 1.0;
  msh.fac2poi(0,0) = 0;
  msh.fac2poi(0,1) = 1;
  msh.fac2poi(0,2) = 2;

  const double bary[3] = {1.0/3.0,1.0/3.0,1.0/3.0};
  const double metric[3] = {1.0,0.0,1.0};
  BOOST_CHECK_THROW(
      (evaluate_pointwise_objective_derivatives<
          MetricFieldFE,2,2,QuaFun::SizeShape,double>(
              msh,AsDeg::P1,AsDeg::P1,msh.fac2poi[0],bary,metric,
              0,msh.getBasis(),DifVar::None,
              PointwiseDerivativeOrder::Value)),
      MetrisExcept);
}

} // namespace Metris
