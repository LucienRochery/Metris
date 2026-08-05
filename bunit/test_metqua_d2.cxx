//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php


#define BOOST_TEST_MODULE test_metqua_d2

#include <boost/test/included/unit_test.hpp>
#include <algorithm>
#include <array>
#include <random>

#include "gen_bary.hxx"

#include "utils/CT_loop.hxx"
#include "utils/mprintf.hxx"
#include "quality/low_metqua.hxx"
#include "quality/low_metqua_d.hxx"
#include "metris_options.hxx"
#include "MetrisRunner/MetrisRunner.hxx"
#include "Mesh/Mesh.hxx"
#include "aux_badEntHandler.hxx"

#include "quality/low_metqua.hxx"
#include "quality/quafun_tradet.hxx"
#include "quality/quafun.hxx"
#include "quality/aux_volumeMeasure.hxx"

#include "libs/SANS/Surreal/SurrealS.h"

namespace Metris{

typedef MetricFieldAnalytical MFT;
typedef double ftype;
typedef std::pair<AsDeg,AsDeg> AsDegPair;

BOOST_AUTO_TEST_CASE(test_bad_ent_handler_global_totals_include_seen_entities)
{
  std::array<double,3> qualities = {3.0,2.0,1.0};
  const std::array<double,3> weights = {2.0,3.0,5.0};
  std::array<bool,3> alive = {true,true,true};

  BadEntHandler handler(1,100.0,0.00001);
  handler.setCallbacks(
      [&](int ientt){ return qualities[ientt]; },
      [&](int ientt){ return !alive[ientt]; });
  handler.setObjectiveWeightCallback(
      [&](int ientt){ return weights[ientt]; });
  handler.seedFromSortedIDs({0,1,2});

  intAr2 ent2ent(3,2);
  intAr2 ent2tag(2,3);
  for(int ientt = 0; ientt < 3; ientt++){
    for(int ii = 0; ii < 2; ii++) ent2ent(ientt,ii) = -1;
  }
  for(int ithrd = 0; ithrd < 2; ithrd++){
    for(int ientt = 0; ientt < 3; ientt++) ent2tag(ithrd,ientt) = 0;
  }

  BOOST_CHECK_EQUAL(handler.getQualityCount(),3);
  BOOST_CHECK_CLOSE_FRACTION(handler.getQualitySum(),6.0,1.e-15);
  BOOST_CHECK_CLOSE_FRACTION(handler.getObjectiveWeightSum(),10.0,1.e-15);
  BOOST_CHECK_CLOSE_FRACTION(handler.getWeightedQualitySum(),17.0,1.e-15);

  // Processing entity 1 moves the worse entity 0 out of K and into
  // seenEntts. Its aggregate contribution must nevertheless stay tracked.
  qualities[1] = 1.5;
  handler.affectedEnttsAlive[1] = qualities[1];
  handler.updateK(1,ent2ent,ent2tag,1,0);
  BOOST_CHECK_EQUAL(handler.getQualityCount(),3);
  BOOST_CHECK_CLOSE_FRACTION(handler.getQualitySum(),5.5,1.e-15);
  BOOST_CHECK_CLOSE_FRACTION(handler.getWeightedQualitySum(),15.5,1.e-15);

  // Updating an element that is already in seenEntts replaces its old global
  // contribution; it must not be counted as a newly created element.
  qualities[0] = 2.5;
  handler.affectedEnttsAlive[0] = qualities[0];
  handler.updateK(2,ent2ent,ent2tag,2,0);
  BOOST_CHECK_EQUAL(handler.getQualityCount(),3);
  BOOST_CHECK_CLOSE_FRACTION(handler.getQualitySum(),5.0,1.e-15);
  BOOST_CHECK_CLOSE_FRACTION(handler.getObjectiveWeightSum(),10.0,1.e-15);
  BOOST_CHECK_CLOSE_FRACTION(handler.getWeightedQualitySum(),14.5,1.e-15);

  // Deleting another seen element must remove its stored contribution.
  alive[1] = false;
  qualities[0] = 2.4;
  handler.deadEntts.push_back(1);
  handler.affectedEnttsAlive[0] = qualities[0];
  handler.updateK(0,ent2ent,ent2tag,3,0);
  BOOST_CHECK_EQUAL(handler.getQualityCount(),2);
  BOOST_CHECK_CLOSE_FRACTION(handler.getQualitySum(),3.4,1.e-15);
  BOOST_CHECK_CLOSE_FRACTION(handler.getObjectiveWeightSum(),7.0,1.e-15);
  BOOST_CHECK_CLOSE_FRACTION(handler.getWeightedQualitySum(),9.8,1.e-15);
}

BOOST_AUTO_TEST_CASE(test_stepdistance_arbitrary_p_derivatives)
{
  MetrisParameters param;
  param.iverb = 0;
  Mesh<MetricFieldFE> msh;
  msh.idim = 2;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &param;
  msh.set_npoin(3);
  msh.set_nface(1);
  msh.coord(0,0) = 0.0; msh.coord(0,1) = 0.0;
  msh.coord(1,0) = 1.2; msh.coord(1,1) = 0.1;
  msh.coord(2,0) = 0.2; msh.coord(2,1) = 0.7;
  msh.fac2poi(0,0) = 0;
  msh.fac2poi(0,1) = 1;
  msh.fac2poi(0,2) = 2;

  constexpr int gdim = 2;
  constexpr int tdim = 2;
  constexpr int nhess = 3;
  constexpr auto quafun =
      get_quafun_xi<MetricFieldFE,gdim,tdim,
                    QuaFun::StepDistance,double>();
  constexpr auto dquafun =
      get_d_quafun_xi<MetricFieldFE,gdim,tdim,
                      QuaFun::StepDistance,double>();

  const int iface = 0;
  const int* nodes = msh.fac2poi[iface];
  const double bary[3] = {1./3.,1./3.,1./3.};
  const double met[3] = {1.,0.,1.};
  const double fd_step = 2.e-6;
  msh.param->step_distance_regularization = 1.e-8;

  for(double power : {0.5,1.0,1.5,2.0,3.25}){
    msh.param->step_distance_p = power;
    for(int ivar = 0; ivar < 3; ivar++){
      double grad[2];
      double hess[nhess];
      const double value = quafun(
          msh,AsDeg::P1,AsDeg::P1,nodes,bary,met);
      const double differentiated_value = dquafun(
          msh,AsDeg::P1,AsDeg::P1,nodes,bary,met,
          ivar,msh.getBasis(),DifVar::None,grad,hess);
      BOOST_CHECK_SMALL(value-differentiated_value,1.e-12);

      const int ipoin = nodes[ivar];
      for(int j = 0; j < gdim; j++){
        const double coordinate = msh.coord(ipoin,j);
        double grad_plus[gdim];
        double grad_minus[gdim];

        msh.coord(ipoin,j) = coordinate+fd_step;
        const double value_plus = quafun(
            msh,AsDeg::P1,AsDeg::P1,nodes,bary,met);
        dquafun(msh,AsDeg::P1,AsDeg::P1,nodes,bary,met,
                ivar,msh.getBasis(),DifVar::None,grad_plus,NULL);

        msh.coord(ipoin,j) = coordinate-fd_step;
        const double value_minus = quafun(
            msh,AsDeg::P1,AsDeg::P1,nodes,bary,met);
        dquafun(msh,AsDeg::P1,AsDeg::P1,nodes,bary,met,
                ivar,msh.getBasis(),DifVar::None,grad_minus,NULL);
        msh.coord(ipoin,j) = coordinate;

        const double grad_fd = (value_plus-value_minus)/(2.*fd_step);
        BOOST_CHECK_SMALL(
            grad[j]-grad_fd,2.e-5*std::max(1.,std::abs(grad_fd)));

        for(int i = 0; i <= j; i++){
          const double hess_fd =
              (grad_plus[i]-grad_minus[i])/(2.*fd_step);
          BOOST_CHECK_SMALL(
              hess[sym2idx(i,j)]-hess_fd,
              5.e-4*std::max(1.,std::abs(hess_fd)));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_stepdistance_metric_volume_barrier_derivatives)
{
  constexpr int gdim = 2;
  constexpr int tdim = 2;
  constexpr int nhess = 3;
  const double met[3] = {2.0,0.2,1.3};
  const double gradN[2] = {-0.7,0.3};
  const double rho0 = 0.7;
  const double beta = 2.0;
  const double fd_step = 1.e-6;
  double Jreg_T[4] = {0.4,0.1,0.05,0.8};

  double rho;
  double barrier;
  double grad[gdim];
  double hess[nhess];
  VolumeMeasureHelpers::eval_metric_volume_barrier_fixed_metric_grad<
      gdim,tdim,double>(
          Jreg_T,met,gradN,rho0,beta,&rho,&barrier,grad);
  VolumeMeasureHelpers::
      eval_metric_volume_barrier_fixed_metric_hess_by_surreal<gdim,tdim>(
          Jreg_T,met,gradN,rho0,beta,hess);

  BOOST_TEST(rho < rho0);
  BOOST_TEST(barrier > 0.0);

  for(int j = 0; j < gdim; j++){
    double Jplus[4];
    double Jminus[4];
    std::copy(Jreg_T,Jreg_T+4,Jplus);
    std::copy(Jreg_T,Jreg_T+4,Jminus);
    for(int i = 0; i < tdim; i++){
      Jplus[i*gdim+j] += fd_step*gradN[i];
      Jminus[i*gdim+j] -= fd_step*gradN[i];
    }

    double rho_plus,rho_minus,barrier_plus,barrier_minus;
    double grad_plus[gdim],grad_minus[gdim];
    VolumeMeasureHelpers::eval_metric_volume_barrier_fixed_metric_grad<
        gdim,tdim,double>(
            Jplus,met,gradN,rho0,beta,
            &rho_plus,&barrier_plus,grad_plus);
    VolumeMeasureHelpers::eval_metric_volume_barrier_fixed_metric_grad<
        gdim,tdim,double>(
            Jminus,met,gradN,rho0,beta,
            &rho_minus,&barrier_minus,grad_minus);

    const double grad_fd = (barrier_plus-barrier_minus)/(2.*fd_step);
    BOOST_CHECK_SMALL(
        grad[j]-grad_fd,1.e-7*std::max(1.,std::abs(grad_fd)));
    for(int i = 0; i <= j; i++){
      const double hess_fd =
          (grad_plus[i]-grad_minus[i])/(2.*fd_step);
      BOOST_CHECK_SMALL(
          hess[sym2idx(i,j)]-hess_fd,
          2.e-6*std::max(1.,std::abs(hess_fd)));
    }
  }

  double inactive_barrier;
  double inactive_grad[gdim];
  double inactive_rho;
  VolumeMeasureHelpers::eval_metric_volume_barrier_fixed_metric_grad<
      gdim,tdim,double>(
          Jreg_T,met,gradN,0.1,beta,
          &inactive_rho,&inactive_barrier,inactive_grad);
  BOOST_CHECK_SMALL(inactive_barrier,1.e-15);
  BOOST_CHECK_SMALL(inactive_grad[0],1.e-15);
  BOOST_CHECK_SMALL(inactive_grad[1],1.e-15);
}

BOOST_AUTO_TEST_CASE(test_integrated_stepdistance_objective_derivatives)
{
  MetrisParameters param;
  param.iverb = 0;
  param.opt_pnorm = 1;
  param.step_distance_p = 1.0;
  param.step_distance_regularization = 1.e-8;
  param.step_distance_barrier_rho0 = 0.7;

  Mesh<MetricFieldFE> msh;
  msh.idim = 2;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &param;
  msh.set_npoin(3);
  msh.set_nface(1);
  msh.met.forceBasisFlag(FEBasis::Lagrange);
  msh.met.forceSpaceFlag(MetSpace::Exp);

  msh.coord(0,0) = 0.0; msh.coord(0,1) = 0.0;
  msh.coord(1,0) = 0.5; msh.coord(1,1) = 0.0;
  msh.coord(2,0) = 0.1; msh.coord(2,1) = 0.3;
  msh.fac2poi(0,0) = 0;
  msh.fac2poi(0,1) = 1;
  msh.fac2poi(0,2) = 2;
  for(int ipoin = 0; ipoin < 3; ipoin++){
    msh.met(ipoin,0) = 1.0;
    msh.met(ipoin,1) = 0.0;
    msh.met(ipoin,2) = 1.0;
  }

  constexpr int gdim = 2;
  constexpr int tdim = 2;
  constexpr int nhess = 3;
  const int ivar = 0;
  const int ipoin = msh.fac2poi(0,ivar);
  const double fd_step = 1.e-6;

  auto evaluate = [&](double beta,
                      double* grad,
                      double* hess){
    msh.param->step_distance_barrier_beta = beta;
    if(grad == NULL){
      return metqua<MetricFieldFE,gdim,tdim,
                    QuaFun::StepDistance,double>(
          msh,AsDeg::P1,AsDeg::P1,0,1.0);
    }
    return d_metqua<MetricFieldFE,gdim,tdim,
                    QuaFun::StepDistance,double>(
        msh,AsDeg::P1,AsDeg::P1,0,ivar,
        msh.getBasis(),DifVar::None,grad,hess,1.0);
  };

  double grad_active[gdim],grad_inactive[gdim];
  double hess_active[nhess],hess_inactive[nhess];
  const double value_active = evaluate(2.0,grad_active,hess_active);
  const double value_inactive = evaluate(0.0,grad_inactive,hess_inactive);
  BOOST_TEST(value_active > value_inactive);

  for(int j = 0; j < gdim; j++){
    const double coordinate = msh.coord(ipoin,j);
    msh.coord(ipoin,j) = coordinate+fd_step;
    const double barrier_plus = evaluate(2.0,NULL,NULL)-evaluate(0.0,NULL,NULL);
    double grad_plus_active[gdim],grad_plus_inactive[gdim];
    evaluate(2.0,grad_plus_active,NULL);
    evaluate(0.0,grad_plus_inactive,NULL);

    msh.coord(ipoin,j) = coordinate-fd_step;
    const double barrier_minus = evaluate(2.0,NULL,NULL)-evaluate(0.0,NULL,NULL);
    double grad_minus_active[gdim],grad_minus_inactive[gdim];
    evaluate(2.0,grad_minus_active,NULL);
    evaluate(0.0,grad_minus_inactive,NULL);
    msh.coord(ipoin,j) = coordinate;

    const double barrier_grad = grad_active[j]-grad_inactive[j];
    const double barrier_grad_fd =
        (barrier_plus-barrier_minus)/(2.*fd_step);
    BOOST_CHECK_SMALL(
        barrier_grad-barrier_grad_fd,
        1.e-7*std::max(1.,std::abs(barrier_grad_fd)));

    for(int i = 0; i <= j; i++){
      const double barrier_hess =
          hess_active[sym2idx(i,j)]-hess_inactive[sym2idx(i,j)];
      const double barrier_hess_fd =
          ((grad_plus_active[i]-grad_plus_inactive[i])
           -(grad_minus_active[i]-grad_minus_inactive[i]))/(2.*fd_step);
      BOOST_CHECK_SMALL(
          barrier_hess-barrier_hess_fd,
          2.e-6*std::max(1.,std::abs(barrier_hess_fd)));
    }
  }

  // Exercise the P1 target-volume normalized element value used by
  // CavityTargetAverage with nonuniform target densities. The metric samples
  // and their normalized weights are frozen in the derivative model.
  msh.met(0,0) = 1.0; msh.met(0,1) = 0.0; msh.met(0,2) = 1.0;
  msh.met(1,0) = 4.0; msh.met(1,1) = 0.0; msh.met(1,2) = 1.0;
  msh.met(2,0) = 1.0; msh.met(2,1) = 0.0; msh.met(2,2) = 9.0;

  msh.param->step_distance_cavity_target_average = true;
  double grad_target_average[gdim],hess_target_average[nhess];
  const double value_target_average =
      evaluate(2.0,grad_target_average,hess_target_average);
  double grad_target_average_no_barrier[gdim];
  double hess_target_average_no_barrier[nhess];
  const double value_target_average_no_barrier =
      evaluate(0.0,grad_target_average_no_barrier,
               hess_target_average_no_barrier);

  constexpr auto pointwise_step_distance =
      get_quafun_xi<MetricFieldFE,gdim,tdim,
                    QuaFun::StepDistance,double>();
  double expected_numerator = 0.0;
  double expected_denominator = 0.0;
  for(int iquad = 0; iquad < tdim + 2; iquad++){
    double bary[tdim + 1] = {};
    if(iquad < tdim + 1){
      bary[iquad] = 1.0;
    }else{
      for(int i = 0; i < tdim + 1; i++){
        bary[i] = 1.0/(tdim + 1);
      }
    }

    double met[3];
    if(iquad < tdim + 1){
      const int metric_point = msh.fac2poi(0,iquad);
      for(int i = 0; i < 3; i++) met[i] = msh.met(metric_point,i);
    }else{
      msh.met.getMetBary(AsDeg::P1,DifVar::None,MetSpace::Exp,
                         msh.fac2poi[0],tdim,bary,
                         met,NULL);
    }

    const double target_density =
        VolumeMeasureHelpers::eval_target_metric_volume_density<
            gdim,double>(met);
    const double phi = pointwise_step_distance(
        msh,AsDeg::P1,AsDeg::P1,msh.fac2poi[0],bary,met);
    expected_numerator += target_density*phi;
    expected_denominator += target_density;
  }
  BOOST_CHECK_CLOSE_FRACTION(
      value_target_average,
      expected_numerator/expected_denominator,
      2.e-14);

  BOOST_CHECK_SMALL(
      value_target_average-value_target_average_no_barrier,1.e-14);
  for(int i = 0; i < gdim; i++){
    BOOST_CHECK_SMALL(
        grad_target_average[i]-grad_target_average_no_barrier[i],1.e-14);
  }
  for(int i = 0; i < nhess; i++){
    BOOST_CHECK_SMALL(
        hess_target_average[i]-hess_target_average_no_barrier[i],1.e-14);
  }

  for(int j = 0; j < gdim; j++){
    const double coordinate = msh.coord(ipoin,j);
    msh.coord(ipoin,j) = coordinate+fd_step;
    const double value_plus = evaluate(2.0,NULL,NULL);
    double grad_plus[gdim];
    evaluate(2.0,grad_plus,NULL);

    msh.coord(ipoin,j) = coordinate-fd_step;
    const double value_minus = evaluate(2.0,NULL,NULL);
    double grad_minus[gdim];
    evaluate(2.0,grad_minus,NULL);
    msh.coord(ipoin,j) = coordinate;

    const double grad_fd = (value_plus-value_minus)/(2.*fd_step);
    BOOST_CHECK_SMALL(
        grad_target_average[j]-grad_fd,
        2.e-5*std::max(1.,std::abs(grad_fd)));

    for(int i = 0; i <= j; i++){
      const double hess_fd =
          (grad_plus[i]-grad_minus[i])/(2.*fd_step);
      BOOST_CHECK_SMALL(
          hess_target_average[sym2idx(i,j)]-hess_fd,
          5.e-4*std::max(1.,std::abs(hess_fd)));
    }
  }

  double moving_coordinates[4];
  for(int ip = 1; ip < 3; ip++){
    for(int j = 0; j < gdim; j++){
      moving_coordinates[(ip-1)*gdim+j] = msh.coord(ip,j);
      msh.coord(ip,j) *= 1.e-8;
    }
  }
  const double collapsed_value = evaluate(2.0,NULL,NULL);
  for(int ip = 1; ip < 3; ip++){
    for(int j = 0; j < gdim; j++){
      msh.coord(ip,j) = moving_coordinates[(ip-1)*gdim+j];
    }
  }
  BOOST_TEST(collapsed_value > value_target_average);

  const double global_sum = 5.0;
  const int global_count = 10;
  const int old_cavity_count = 2;
  const int new_cavity_count = 3;

  // Improving the cavity ratio alone does not imply improvement of the full
  // mesh ratio when the replacement changes the target-weight denominator.
  const double old_weighted_cavity_sum = 2.0;
  const double old_cavity_weight = 2.0;
  const double new_weighted_cavity_sum = 2.7;
  const double new_cavity_weight = 3.0;
  BOOST_TEST(new_weighted_cavity_sum/new_cavity_weight
             < old_weighted_cavity_sum/old_cavity_weight);
  const double new_weighted_global_objective =
      step_distance_replaced_region_objective(
          global_sum,old_weighted_cavity_sum,new_weighted_cavity_sum,
          10.0,old_cavity_weight,new_cavity_weight);
  BOOST_TEST(new_weighted_global_objective
             > global_sum/10.0);
  BOOST_CHECK_CLOSE_FRACTION(
      new_weighted_global_objective,5.7/11.0,1.e-15);

  StepDistanceObjectiveState cavity_average_state;
  cavity_average_state.numerator = global_sum;
  cavity_average_state.element_count = global_count;
  cavity_average_state.target_weight = 10.;
  cavity_average_state.best_objective = global_sum/10.;
  cavity_average_state.cavity_global_relative_tolerance = 1.e-6;
  cavity_average_state.cavity_global_gain_fraction = 0.05;

  // A substantial global worsening is rejected despite local improvement.
  BOOST_CHECK(!cavity_average_state.accepts_replacement(
      old_weighted_cavity_sum,old_cavity_count,old_cavity_weight,
      new_weighted_cavity_sum,new_cavity_count,new_cavity_weight));

  // A tiny global worsening supported by a strong mesh-scaled local gain is
  // accepted inside the best-so-far envelope.
  const double filtered_cavity_sum = 2.500001;
  BOOST_CHECK(cavity_average_state.accepts_replacement(
      old_weighted_cavity_sum,old_cavity_count,old_cavity_weight,
      filtered_cavity_sum,new_cavity_count,new_cavity_weight));

  cavity_average_state.best_objective = 0.499;
  BOOST_CHECK(!cavity_average_state.accepts_replacement(
      old_weighted_cavity_sum,old_cavity_count,old_cavity_weight,
      filtered_cavity_sum,new_cavity_count,new_cavity_weight));
  cavity_average_state.best_objective = global_sum/10.;

  BOOST_CHECK(!cavity_target_average_global_filter_accepts(
      1.0,0.999999,0.5,0.5000001,0.5,
      2.0,3.0,10.0,1.e-6,0.05));

  // Local worsening never reaches the global filter.
  BOOST_CHECK(!cavity_average_state.accepts_replacement(
      old_weighted_cavity_sum,old_cavity_count,old_cavity_weight,
      3.3,new_cavity_count,new_cavity_weight));

  const double best_before_replace = cavity_average_state.best_objective;
  cavity_average_state.replace(
      old_weighted_cavity_sum,old_cavity_count,old_cavity_weight,
      filtered_cavity_sum,new_cavity_count,new_cavity_weight);
  BOOST_CHECK_EQUAL(cavity_average_state.best_objective,best_before_replace);

  // Cavity-level target average. Regional aggregation uses
  // sum_K(D_K e_K)/sum_K(D_K).
  msh.set_npoin(4);
  msh.set_nface(2);
  msh.coord(3,0) = 0.7; msh.coord(3,1) = 0.4;
  msh.fac2poi(1,0) = 1;
  msh.fac2poi(1,1) = 3;
  msh.fac2poi(1,2) = 2;
  msh.met(3,0) = 16.0; msh.met(3,1) = 0.0; msh.met(3,2) = 0.25;

  msh.param->step_distance_cavity_target_average = true;

  auto evaluate_cavity_target_average = [&](double* grad,
                                             double* hess){
    double numerator = 0.;
    double denominator = 0.;
    for(int iface = 0; iface < 2; iface++){
      double element_gradient[gdim] = {};
      double element_hessian[nhess] = {};
      double element_value;
      if(iface == 0 && grad != NULL){
        element_value =
            d_metqua<MetricFieldFE,gdim,tdim,
                     QuaFun::StepDistance,double>(
                msh,AsDeg::P1,AsDeg::P1,iface,ivar,
                msh.getBasis(),DifVar::None,
                element_gradient,hess == NULL ? NULL : element_hessian,1.0);
      }else{
        element_value =
            metqua<MetricFieldFE,gdim,tdim,
                   QuaFun::StepDistance,double>(
                msh,AsDeg::P1,AsDeg::P1,iface,1.0);
      }

      const double element_weight =
          step_distance_element_target_weight<MetricFieldFE,gdim,tdim>(
              msh,AsDeg::P1,iface);
      numerator += step_distance_region_contribution(
          element_value,element_weight,true);
      denominator += element_weight;

      if(iface == 0 && grad != NULL){
        for(int i = 0; i < gdim; i++){
          grad[i] = element_weight*element_gradient[i];
        }
        if(hess != NULL){
          for(int i = 0; i < nhess; i++){
            hess[i] = element_weight*element_hessian[i];
          }
        }
      }
    }

    if(grad != NULL){
      for(int i = 0; i < gdim; i++) grad[i] /= denominator;
      if(hess != NULL){
        for(int i = 0; i < nhess; i++) hess[i] /= denominator;
      }
    }
    return step_distance_region_objective(
        numerator,denominator,true);
  };

  double cavity_gradient[gdim];
  double cavity_hessian[nhess];
  msh.param->step_distance_barrier_beta = 50.0;
  const double cavity_value =
      evaluate_cavity_target_average(cavity_gradient,cavity_hessian);
  msh.param->step_distance_barrier_beta = 0.0;
  const double cavity_value_no_barrier =
      evaluate_cavity_target_average(NULL,NULL);
  BOOST_CHECK_SMALL(cavity_value-cavity_value_no_barrier,1.e-14);

  const double element_value0 =
      metqua<MetricFieldFE,gdim,tdim,QuaFun::StepDistance,double>(
          msh,AsDeg::P1,AsDeg::P1,0,1.0);
  const double element_value1 =
      metqua<MetricFieldFE,gdim,tdim,QuaFun::StepDistance,double>(
          msh,AsDeg::P1,AsDeg::P1,1,1.0);
  const double element_weight0 =
      step_distance_element_target_weight<MetricFieldFE,gdim,tdim>(
          msh,AsDeg::P1,0);
  const double element_weight1 =
      step_distance_element_target_weight<MetricFieldFE,gdim,tdim>(
          msh,AsDeg::P1,1);
  const double expected_cavity_value =
      (element_weight0*element_value0 + element_weight1*element_value1)
      /(element_weight0 + element_weight1);
  BOOST_CHECK_CLOSE_FRACTION(
      cavity_value,expected_cavity_value,2.e-14);

  // Regression for cavity smoothing acceptance: a decreasing weighted
  // numerator is not sufficient when the target-volume denominator changes.
  // These are the values from the BL a0 cavity that originally triggered the
  // post-smoothing assertion.
  const double before_numerator = 2.13698184801736318e+02;
  const double before_denominator = 2.24266627232477333e+02;
  const double after_numerator = 2.13686305062136768e+02;
  const double after_denominator = 2.24237408380270040e+02;
  const double before_objective = step_distance_region_objective(
      before_numerator,before_denominator,true);
  const double after_objective = step_distance_region_objective(
      after_numerator,after_denominator,true);
  BOOST_CHECK(after_numerator < before_numerator);
  BOOST_CHECK(after_objective > before_objective);

  for(int j = 0; j < gdim; j++){
    const double coordinate = msh.coord(ipoin,j);
    msh.coord(ipoin,j) = coordinate+fd_step;
    const double value_plus = evaluate_cavity_target_average(NULL,NULL);
    double gradient_plus[gdim];
    evaluate_cavity_target_average(gradient_plus,NULL);

    msh.coord(ipoin,j) = coordinate-fd_step;
    const double value_minus = evaluate_cavity_target_average(NULL,NULL);
    double gradient_minus[gdim];
    evaluate_cavity_target_average(gradient_minus,NULL);
    msh.coord(ipoin,j) = coordinate;

    const double gradient_fd = (value_plus-value_minus)/(2.*fd_step);
    BOOST_CHECK_SMALL(
        cavity_gradient[j]-gradient_fd,
        2.e-5*std::max(1.,std::abs(gradient_fd)));
    for(int i = 0; i <= j; i++){
      const double hessian_fd =
          (gradient_plus[i]-gradient_minus[i])/(2.*fd_step);
      BOOST_CHECK_SMALL(
          cavity_hessian[sym2idx(i,j)]-hessian_fd,
          5.e-4*std::max(1.,std::abs(hessian_fd)));
    }
  }
}

BOOST_AUTO_TEST_CASE(test_metqua_d2)
{
  std::cout << "================================" << std::endl;
  std::cout << "Testing with quafun_distortion" << std::endl;
  std::cout << "================================" << std::endl;

  std::vector<std::string> meshes =
  {METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k"
  ,METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k"
  ,METRIS_CASES_DIR "/unit/2D/square/circmet.p2.5k"
  };


  for(auto testcase : meshes){
    std::string mesh_name = testcase;
    std::cout<<"-- Test case mesh = "<<mesh_name<<std::endl;
    cargHandler arg("-in "+mesh_name+" -verb 0 -anamet 1 -sclmet 1");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
    msh.met.setSpace(MetSpace::Log);

    dblAr2 bary[4];
    const int nsamp = 20;
    genBary(nsamp, 2, bary[2]);
    genBary(nsamp, 3, bary[3]);

    GETVDEPTH(msh.param);

    CT_FOR0_INC(2,3,gdim){if(msh.idim == gdim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
      constexpr int nhess = (gdim*(gdim+1))/2;
      // Use Boost PP to loop between options in the future
      constexpr auto iquaf = QuaFun::Distortion;


      CT_FOR0_INC(2,gdim,tdim){if(tdim <= msh.get_tdim()){
        INCVDEPTH(msh.param);
        MPRINTF("-- Tests tdim = {} gdim = {}\n",tdim,gdim);
        constexpr auto quafun_xi   = get_quafun_xi<MFT,gdim,tdim,iquaf,ftype>();
        constexpr auto d_quafun_xi = get_d_quafun_xi<MFT,gdim,tdim,iquaf,ftype>();
        const intAr2 &ent2poi = msh.ent2poi(tdim);
        const int nnode = getnnode(tdim, ideg);

        {//INCVDEPTH
        INCVDEPTH(msh.param);
        MPRINTF("-- Test 1: quafun_xi == d_quafun_xi\n");
        for(AsDegPair asdeg : {AsDegPair({AsDeg::P1, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::Pk})}){
          std::string asdmsh_str = asdeg.first  == AsDeg::P1 ? "P1" : "Pk";
          std::string asdmet_str = asdeg.second == AsDeg::P1 ? "P1" : "Pk";
          double err_min = 1.0e30;
          double err_max = -1.0;
          double err_avg = 0;
          double ntest = 0;
          for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
            for(int isamp = 0; isamp < nsamp; isamp++){
              double qua = quafun_xi(msh, asdeg.first, asdeg.second, ent2poi[ientt],
                                     bary[tdim][isamp], NULL);
              double dquael[gdim], hquael[nhess];
              for(int ivar = 0; ivar < nnode; ivar++){
                double qua_d = d_quafun_xi(msh, asdeg.first, asdeg.second, ent2poi[ientt],
                                           bary[tdim][isamp], NULL,
                                           ivar, msh.getBasis(), DifVar::None,
                                           dquael, hquael);
                double err = abs(qua - qua_d);
                err_min = MIN(err_min, err);
                err_max = MAX(err_max, err);
                err_avg += err;
                ntest += 1;
              }
            }// for isamp
          }// for ientt
          if(ntest == 0) continue;
          err_avg /= ntest;
          MPRINTF(" - DONE using asdmsh = {}, asdmet = {}:"
                  " got err min = {}, avg = {}, max = {}\n",
                  asdmsh_str.c_str(),asdmet_str.c_str(),
                  err_min,err_avg,err_max);
          BOOST_TEST(err_max < 1.0e-12);
        }// for AsDegPair



        MPRINTF("-- Test 2.1: d_quafun_tradet_SurrealS derivatives (identity)\n");
        {// using
        constexpr int nvar = gdim;
        auto d_quafun = d_quafun_tradet_SurrealS<MFT,gdim,tdim,gdim>;
        auto   quafun =   quafun_tradet         <MFT,gdim,tdim,ftype>;
        SANS::DLA::MatrixS<gdim,nvar,double> dpoint = SANS::DLA::Identity();
        //for(int ii = 0; ii < gdim; ii++){
        //  for(int jj = 0; jj < gdim; jj++){
        //    dpoint(ii,jj) = ii == jj;
        //  }
        //}
        for(AsDegPair asdeg : {AsDegPair({AsDeg::P1, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::Pk})}){
          double ntest = 0;
          for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
            for(int isamp = 0; isamp < nsamp; isamp++){
              SANS::SurrealS<nvar, ftype> tra, det;
              ftype tra0, det0;

              quafun(msh, asdeg.first, asdeg.second,
                     ent2poi[ientt], bary[tdim][isamp],
                     NULL,
                     &tra0, &det0);
              for(int ivar = 0; ivar < nnode; ivar++){
                d_quafun(msh, asdeg.first, asdeg.second,
                         ent2poi[ientt], bary[tdim][isamp],
                         ivar, msh.getBasis(), DifVar::None,
                         NULL,
                         tra, det, &dpoint);
                ntest += 1;
              }
            }// for isamp
          }// for ientt
        }// for AsDegPair
        }// using

        }// INCVDEPTH

      }}CT_FOR1(tdim);
    }}CT_FOR1(ideg);
    }}CT_FOR1(gdim);


  }// for testcase


}// BOOST_AUTO_TEST_CASE

BOOST_AUTO_TEST_CASE(test_quafun_sizeshape)
{

  std::cout << "================================" << std::endl;
  std::cout << "Testing with quafun_sizeshape" << std::endl;
  std::cout << "================================" << std::endl;

  std::vector<std::string> meshes =
  {METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k"
  // ,METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k"
  // ,METRIS_CASES_DIR "/unit/2D/square/circmet.p2.5k"
  };


  for(auto testcase : meshes){
    std::string mesh_name = testcase;
    std::cout<<"-- Test case mesh = "<<mesh_name<<std::endl;
    cargHandler arg("-in "+mesh_name+" -verb 0 -anamet 1 -sclmet 1");
    MetrisRunner run(arg.c, arg.v);
    Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
    msh.met.setSpace(MetSpace::Log);

    dblAr2 bary[4];
    const int nsamp = 20;
    genBary(nsamp, 2, bary[2]);
    genBary(nsamp, 3, bary[3]);

    GETVDEPTH(msh.param);

    CT_FOR0_INC(2,3,gdim){if(msh.idim == gdim){
    CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
      constexpr int nhess = (gdim*(gdim+1))/2;
      // Use Boost PP to loop between options in the future
      constexpr auto iquaf = QuaFun::SizeShape;


      CT_FOR0_INC(2,gdim,tdim){if(tdim <= msh.get_tdim()){
        INCVDEPTH(msh.param);
        MPRINTF("-- Tests tdim = {} gdim = {}\n",tdim,gdim);
        constexpr auto quafun_xi   = get_quafun_xi<MFT,gdim,tdim,iquaf,ftype>();
        constexpr auto d_quafun_xi = get_d_quafun_xi<MFT,gdim,tdim,iquaf,ftype>();
        const intAr2 &ent2poi = msh.ent2poi(tdim);
        const int nnode = getnnode(tdim, ideg);

        {//INCVDEPTH
        INCVDEPTH(msh.param);
        MPRINTF("-- Test 1: quafun_xi == d_quafun_xi\n");
        for(AsDegPair asdeg : {AsDegPair({AsDeg::P1, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::P1}),
                               AsDegPair({AsDeg::Pk, AsDeg::Pk})}){
          std::string asdmsh_str = asdeg.first  == AsDeg::P1 ? "P1" : "Pk";
          std::string asdmet_str = asdeg.second == AsDeg::P1 ? "P1" : "Pk";
          double err_min = 1.0e30;
          double err_max = -1.0;
          double err_avg = 0;
          double ntest = 0;
          for(int ientt = 0; ientt < msh.nentt(tdim); ientt++){
            for(int isamp = 0; isamp < nsamp; isamp++){
              double qua = quafun_xi(msh, asdeg.first, asdeg.second, ent2poi[ientt],
                                     bary[tdim][isamp], NULL);
              double dquael[gdim], hquael[nhess];
              for(int ivar = 0; ivar < nnode; ivar++){
                double qua_d = d_quafun_xi(msh, asdeg.first, asdeg.second, ent2poi[ientt],
                                           bary[tdim][isamp], NULL,
                                           ivar, msh.getBasis(), DifVar::None,
                                           dquael, hquael);
                double err = abs(qua - qua_d);
                err_min = MIN(err_min, err);
                err_max = MAX(err_max, err);
                err_avg += err;
                ntest += 1;
              }
            }// for isamp
          }// for ientt
          if(ntest == 0) continue;
          err_avg /= ntest;
          MPRINTF(" - DONE using asdmsh = {}, asdmet = {}:"
                  " got err min = {}, avg = {}, max = {}\n",
                  asdmsh_str.c_str(),asdmet_str.c_str(),
                  err_min,err_avg,err_max);
          BOOST_TEST(err_max < 1.0e-12);
        }// for AsDegPair

        MPRINTF("-- Test 2: Check quality value for knwon P1 cases (unit regular, scaled regular, anisotropic example)\n");
        if (gdim == 2 && tdim == 2 && ideg == 1 && msh.nentt(tdim) > 0) {
          const int ientt = 0;
          intAr2& ent2poi = msh.ent2poi(2);

          // save old coords to restore later
          dblAr2 oldCoord(3,2);
          for (int ipoi = 0; ipoi < 3; ipoi++) {
            oldCoord(ipoi,0) = msh.coord(ent2poi(ientt,ipoi),0);
            oldCoord(ipoi,1) = msh.coord(ent2poi(ientt,ipoi),1);
          }

          // ===================
          // Test unit regular
          // ===================

          // overwrite with unit equilateral triangle
          msh.coord(ent2poi(ientt,0),0) = 0.; msh.coord(ent2poi(ientt,0),1) = 0.;
          msh.coord(ent2poi(ientt,1),0) = 1.; msh.coord(ent2poi(ientt,1),1) = 0.;
          msh.coord(ent2poi(ientt,2),0) = 0.5; msh.coord(ent2poi(ientt,2),1) = sqrt(3.)/2.;

          // barycenter
          const double btri[3] = { 1./3., 1./3., 1./3. };

          // use identity metric for this test
          const double metI[3] = { 1., 0., 1.};

          double trueQuality = 1.;

          for (int powerSign : {+1, -1}) {
            msh.param->opt_power = powerSign;

            // value-only
            double qv = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);

            // returns value as well, and additionally compute derivatives
            double dquael_loc[gdim], hquael_loc[(gdim*(gdim+1))/2];
            double qv_d = d_quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri,
                                    metI, /*ivar=*/0, msh.getBasis(), DifVar::None,
                                    dquael_loc, hquael_loc);

            BOOST_CHECK_SMALL(qv - trueQuality, 1e-13);
            BOOST_CHECK_SMALL(qv_d - trueQuality, 1e-13);
          }

          // ===================
          // Test scaled regular
          // ===================

          const double scale = 8.;

          msh.coord(ent2poi(ientt,0),0) *= scale; msh.coord(ent2poi(ientt,0),1) *= scale;
          msh.coord(ent2poi(ientt,1),0) *= scale; msh.coord(ent2poi(ientt,1),1) *= scale;
          msh.coord(ent2poi(ientt,2),0) *= scale; msh.coord(ent2poi(ientt,2),1) *= scale;

          trueQuality = 1./2. * (pow(scale,4.) + pow(scale,-4.));

          for (int powerSign : {+1, -1}) {
            msh.param->opt_power = powerSign;

            // value-only
            double qv = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);

            // returns value as well, and additionally compute derivatives
            double dquael_loc[gdim], hquael_loc[(gdim*(gdim+1))/2];
            double qv_d = d_quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri,
                                    metI, /*ivar=*/0, msh.getBasis(), DifVar::None,
                                    dquael_loc, hquael_loc);

            if (powerSign == 1){
              BOOST_CHECK_SMALL(qv - trueQuality, 1e-12);
              BOOST_CHECK_SMALL(qv_d - trueQuality, 1e-12);
            }
            else{
              BOOST_CHECK_SMALL(qv - 1./trueQuality, 1e-12);
              BOOST_CHECK_SMALL(qv_d - 1./trueQuality, 1e-12);
            }
          }

          // ===================
          // Test anisotropic
          // ===================

          // vertices: (0,0) - (1,0) - (0.5,0.1)

          msh.coord(ent2poi(ientt,0),0) = 0.; msh.coord(ent2poi(ientt,0),1) = 0.;
          msh.coord(ent2poi(ientt,1),0) = 1.; msh.coord(ent2poi(ientt,1),1) = 0.;
          msh.coord(ent2poi(ientt,2),0) = 0.5; msh.coord(ent2poi(ientt,2),1) = 0.1;

          const dblAr1 e1({msh.coord(ent2poi(ientt,1),0) - msh.coord(ent2poi(ientt,0),0), msh.coord(ent2poi(ientt,1),1) - msh.coord(ent2poi(ientt,0),1)});
          const dblAr1 e2({msh.coord(ent2poi(ientt,2),0) - msh.coord(ent2poi(ientt,0),0), msh.coord(ent2poi(ientt,2),1) - msh.coord(ent2poi(ientt,0),1)});

          double tra = 4./3. * (e1[0]*e1[0] + e1[1]*e1[1] +  e2[0]*e2[0] + e2[1]*e2[1] - e1[0]*e2[0] - e1[1]*e2[1]);
          double det = 4./3. * (e1[0]*e2[1] - e1[1]*e2[0])*(e1[0]*e2[1] - e1[1]*e2[0]);

          trueQuality = 1./(2*2*2) * tra * tra * (1. + 1./(det*det));

          for (int powerSign : {+1, -1}) {
            msh.param->opt_power = powerSign;

            // value-only
            double qv = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);

            // returns value as well, and additionally compute derivatives
            double dquael_loc[gdim], hquael_loc[(gdim*(gdim+1))/2];
            double qv_d = d_quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri,
                                    metI, /*ivar=*/0, msh.getBasis(), DifVar::None,
                                    dquael_loc, hquael_loc);

            if (powerSign == 1){
              BOOST_CHECK_SMALL(qv - trueQuality, 1e-12);
              BOOST_CHECK_SMALL(qv_d - trueQuality, 1e-12);
            }
            else{
              BOOST_CHECK_SMALL(qv - 1./trueQuality, 1e-12);
              BOOST_CHECK_SMALL(qv_d - 1./trueQuality, 1e-12);
            }
          }

          // Restore original coords
          for (int ipoin = 0; ipoin < 3; ipoin++) {
            msh.coord(ent2poi(ientt,ipoin),0) = oldCoord(ipoin,0);
            msh.coord(ent2poi(ientt,ipoin),1) = oldCoord(ipoin,1);
          }


        }
        MPRINTF("-- Test 3: Check quality derivatives using 1st order FD\n");
        if (gdim == 2 && tdim == 2 && ideg == 1 && msh.nentt(tdim) > 0) {

          // use identity metric for this test
          const double metI[3] = {1., 0., 1.};
          const double btri[3]  = {1./3., 1./3., 1./3.};

          const int nnode = getnnode(tdim, ideg);

          // check for several elements (not entire mesh)
          const int ncheck_ent = std::min(10, msh.nentt(tdim));

          // step sizes
          const int nh = 7;
          const double hvec[7] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7};

          // to save max error
          double err_max_by_h[7] = {0,0,0,0,0,0,0};

          intAr2& ent2poi = msh.ent2poi(2);

          for (int ientt = 0; ientt < ncheck_ent; ientt++) {

            // Compute analytic value and gradient
            double qua0;
            dblAr1 grad_analytic(6); // 3 nodes with 2 components each

            // we will need the coords backup to restore after each perturbation
            dblAr2 coord_backup(3,2);
            for (int ipoin = 0; ipoin < 3; ipoin++) {
              for (int ii = 0; ii < gdim; ii++)
                coord_backup(ipoin,ii) = msh.coord(ent2poi(ientt,ipoin), ii);
            }

            // get base value once
            qua0 = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);

            // fill analytic gradient per node
            for (int ivar = 0; ivar < nnode; ivar++) {
              double dqa[gdim], hqa[(gdim*(gdim+1))/2];
              double qua_d = d_quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI,
                                        ivar, msh.getBasis(), DifVar::None,
                                        dqa, hqa);
              // consistency : qua_d == qua0
              BOOST_CHECK_SMALL(qua_d - qua0, 1e-12);
              for (int ii = 0; ii < gdim; ii++)
                grad_analytic[ivar*gdim + ii] = dqa[ii];
            }

            // now finite differences for each coordinate
            for (int ivar = 0; ivar < nnode; ivar++) {
              for (int ii = 0; ii < gdim; ii++) {

                double prev_err = std::numeric_limits<double>::quiet_NaN();

                for (int ih = 0; ih < nh; ih++) {
                  const double h = hvec[ih];

                  // perturb x(ivar,k) by +h
                  msh.coord(ent2poi(ientt,ivar), ii) = coord_backup(ivar,ii) + h;

                  // recompute quality value
                  double quah = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);

                  // forward finite difference
                  const double fd = (quah - qua0) / h;

                  // error against analytic
                  const double an = grad_analytic[ivar*gdim + ii];

                  const double err = std::abs(fd - an);
                  err_max_by_h[ih] = std::max(err_max_by_h[ih], err);

                  // crude rate check: first order
                  if (ih > 0 && std::isfinite(prev_err) && prev_err > 0) {
                    const double h_prev = hvec[ih-1];
                    const double order = std::log(prev_err/err) / std::log(h_prev/h);
                  }

                  // restore the perturbed coordinate
                  msh.coord(ent2poi(ientt,ivar), ii) = coord_backup(ivar,ii);
                  prev_err = err;
                } // hvec
              } // ii
            } // ivar
          } // ientt

          // print max errors per h
          MPRINTF(" - Gradient 1st order FD max errors by h:\n");
          for (int ih = 0; ih < nh; ih++) MPRINTF("  - h = {:.1e}, max err = {:.1e}\n",hvec[ih],err_max_by_h[ih]);

        }

        MPRINTF("-- Test 4: Check quality derivatives using 2nd order FD\n");
        if (gdim == 2 && tdim == 2 && ideg == 1 && msh.nentt(tdim) > 0) {

          // use identity metric for this test
          const double metI[3] = {1., 0., 1.};
          const double btri[3]  = {1./3., 1./3., 1./3.};

          const int nnode = getnnode(tdim, ideg);

          // check for several elements (not entire mesh)
          const int ncheck_ent = std::min(10, msh.nentt(tdim));

          // step sizes
          const int nh = 7;
          const double hvec[7] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7};

          // to save max error
          double err_max_by_h[7] = {0,0,0,0,0,0,0};

          intAr2& ent2poi = msh.ent2poi(2);

          for (int ientt = 0; ientt < ncheck_ent; ientt++) {

            // Compute analytic value and gradient
            double qua0;
            dblAr1 grad_analytic(6); // 3 nodes with 2 components each

            // we will need the coords backup to restore after each perturbation
            dblAr2 coord_backup(3,2);
            for (int ipoin = 0; ipoin < 3; ipoin++) {
              for (int ii = 0; ii < gdim; ii++)
                coord_backup(ipoin,ii) = msh.coord(ent2poi(ientt,ipoin), ii);
            }

            // get base value once
            qua0 = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);

            // fill analytic gradient per node
            for (int ivar = 0; ivar < nnode; ivar++) {
              double dqa[gdim], hqa[(gdim*(gdim+1))/2];
              double qua_d = d_quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI,
                                        ivar, msh.getBasis(), DifVar::None,
                                        dqa, hqa);
              // consistency : qua_d == qua0
              BOOST_CHECK_SMALL(qua_d - qua0, 1e-12);
              for (int ii = 0; ii < gdim; ii++)
                grad_analytic[ivar*gdim + ii] = dqa[ii];
            }

            // now finite differences for each coordinate
            for (int ivar = 0; ivar < nnode; ivar++) {
              for (int ii = 0; ii < gdim; ii++) {

                double prev_err = std::numeric_limits<double>::quiet_NaN();

                for (int ih = 0; ih < nh; ih++) {
                  const double h = hvec[ih];

                  // perturb x(ivar,k) by +h
                  msh.coord(ent2poi(ientt,ivar), ii) = coord_backup(ivar,ii) + h;
                  // compute quality value
                  double quahplus = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);
                  // perturb x(ivar,k) by -h
                  msh.coord(ent2poi(ientt,ivar), ii) = coord_backup(ivar,ii) - h;
                  // compute quality value
                  double quahminus = quafun_xi(msh, AsDeg::P1, AsDeg::P1, ent2poi[ientt], btri, metI);

                  // forward finite difference
                  const double fd = (quahplus - quahminus) / (2.*h);

                  // error against analytic
                  const double an = grad_analytic[ivar*gdim + ii];

                  const double err = std::abs(fd - an);
                  err_max_by_h[ih] = std::max(err_max_by_h[ih], err);

                  // crude rate check: first order
                  if (ih > 0 && std::isfinite(prev_err) && prev_err > 0) {
                    const double h_prev = hvec[ih-1];
                    const double order = std::log(prev_err/err) / std::log(h_prev/h);
                  }

                  // restore the perturbed coordinate
                  msh.coord(ent2poi(ientt,ivar), ii) = coord_backup(ivar,ii);
                  prev_err = err;
                } // hvec
              } // ii
            } // ivar
          } // ientt

          // print max errors per h
          MPRINTF(" - Gradient 2nd order FD max errors by h:\n");
          for (int ih = 0; ih < nh; ih++) MPRINTF("  - h = {:.1e}, max err = {:.1e}\n",hvec[ih],err_max_by_h[ih]);

        }
        }// INCVDEPTH

      }}CT_FOR1(tdim);
    }}CT_FOR1(ideg);
    }}CT_FOR1(gdim);

  }// for testcase

}// BOOST_AUTO_TEST_CASE

}// namespace
