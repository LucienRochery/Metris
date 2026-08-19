// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_quality_sizeshape_p1_integrated_derivatives

#include <boost/test/included/unit_test.hpp>

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "linalg/det.hxx"
#include "linalg/symidx.hxx"
#include "low_eval.hxx"
#include "low_geo/measure.hxx"
#include "quality/low_metqua.hxx"
#include "quality/low_metqua_d.hxx"
#include "quality/quafun_sizeshape.hxx"
#include "quality/simplex_quadrature.hxx"

#include <array>
#include <cmath>
#include <type_traits>

namespace Metris
{

namespace
{

void test_analytical_metric_2d(
    const AnaMetCtx *,
    const double *coordinate,
    double scale,
    int derivative_order,
    double *metric,
    double *metric_derivative)
{
  metric[0] = scale*(1.4 + 0.20*coordinate[0]);
  metric[1] = scale*(0.08 + 0.03*coordinate[1]);
  metric[2] = scale*(0.9 + 0.05*coordinate[0] + 0.15*coordinate[1]);

  if(derivative_order == 0) return;
  for(int iderivative = 0; iderivative < 6; iderivative++){
    metric_derivative[iderivative] = 0.0;
  }
  metric_derivative[0*3 + 0] = scale*0.20;
  metric_derivative[0*3 + 2] = scale*0.05;
  metric_derivative[1*3 + 1] = scale*0.03;
  metric_derivative[1*3 + 2] = scale*0.15;
}

void test_analytical_metric_3d(
    const AnaMetCtx *,
    const double *coordinate,
    double scale,
    int derivative_order,
    double *metric,
    double *metric_derivative)
{
  metric[0] = scale*(1.5 + 0.12*coordinate[0]);
  metric[1] = scale*(0.05 + 0.02*coordinate[1]);
  metric[2] = scale*(1.1 + 0.10*coordinate[1]);
  metric[3] = scale*(-0.03 + 0.01*coordinate[2]);
  metric[4] = scale*(0.04 + 0.015*coordinate[0]);
  metric[5] = scale*(0.85 + 0.08*coordinate[2]);

  if(derivative_order == 0) return;
  for(int iderivative = 0; iderivative < 18; iderivative++){
    metric_derivative[iderivative] = 0.0;
  }
  metric_derivative[0*6 + 0] = scale*0.12;
  metric_derivative[0*6 + 4] = scale*0.015;
  metric_derivative[1*6 + 1] = scale*0.02;
  metric_derivative[1*6 + 2] = scale*0.10;
  metric_derivative[2*6 + 3] = scale*0.01;
  metric_derivative[2*6 + 5] = scale*0.08;
}

template<class MFT, int gdim>
void initialize_element(Mesh<MFT> &msh,
                        MetrisParameters &parameters)
{
  static_assert(gdim == 2 || gdim == 3);
  msh.idim = gdim;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &parameters;
  msh.set_npoin(gdim + 1);
  msh.set_nentt(gdim,1);
  msh.met.forceBasisFlag(FEBasis::Lagrange);
  msh.met.forceSpaceFlag(MetSpace::Exp);

  if constexpr(gdim == 2){
    msh.coord(0,0) = 0.05;
    msh.coord(0,1) = 0.10;
    msh.coord(1,0) = 1.15;
    msh.coord(1,1) = 0.18;
    msh.coord(2,0) = 0.25;
    msh.coord(2,1) = 0.92;
  }else{
    msh.coord(0,0) = 0.05;
    msh.coord(0,1) = 0.10;
    msh.coord(0,2) = 0.08;
    msh.coord(1,0) = 1.05;
    msh.coord(1,1) = 0.18;
    msh.coord(1,2) = 0.04;
    msh.coord(2,0) = 0.22;
    msh.coord(2,1) = 0.98;
    msh.coord(2,2) = 0.16;
    msh.coord(3,0) = 0.14;
    msh.coord(3,1) = 0.28;
    msh.coord(3,2) = 0.91;
  }

  intAr2 &element_to_point = msh.ent2poi(gdim);
  for(int inode = 0; inode < gdim + 1; inode++){
    element_to_point(0,inode) = inode;
  }

  if constexpr(std::is_same<MFT,MetricFieldAnalytical>::value){
    if constexpr(gdim == 2){
      parameters.setAnalyticalMetric(
          AnaMetFun(test_analytical_metric_2d));
    }else{
      parameters.setAnalyticalMetric(
          AnaMetFun(test_analytical_metric_3d));
    }
    msh.met.setAnalyticalMetric(parameters);
    for(int ipoin = 0; ipoin < gdim + 1; ipoin++){
      msh.met.getMetPhys(DifVar::None,MetSpace::Exp,
                         msh.coord[ipoin],msh.met[ipoin],NULL);
    }
  }else if constexpr(gdim == 2){
    const double nodal_metrics[3][3] = {
      {1.55, 0.10,0.82},
      {2.05, 0.14,1.12},
      {0.95,-0.04,1.75}
    };
    for(int ipoin = 0; ipoin < 3; ipoin++){
      for(int imetric = 0; imetric < 3; imetric++){
        msh.met(ipoin,imetric) = nodal_metrics[ipoin][imetric];
      }
    }
  }else{
    const double nodal_metrics[4][6] = {
      {1.55, 0.05,1.02, 0.02,-0.03,0.82},
      {2.05, 0.10,1.22,-0.04, 0.02,0.92},
      {1.15,-0.05,1.82, 0.03, 0.06,1.32},
      {0.95, 0.02,1.12,-0.01, 0.04,2.05}
    };
    for(int ipoin = 0; ipoin < 4; ipoin++){
      for(int imetric = 0; imetric < 6; imetric++){
        msh.met(ipoin,imetric) = nodal_metrics[ipoin][imetric];
      }
    }
  }
}

void initialize_embedded_surface_triangle(
    Mesh<MetricFieldFE> &msh,
    MetrisParameters &parameters)
{
  msh.idim = 3;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &parameters;
  msh.set_npoin(3);
  msh.set_nface(1);
  msh.met.forceBasisFlag(FEBasis::Lagrange);
  msh.met.forceSpaceFlag(MetSpace::Exp);

  msh.coord(0,0) = 0.05;
  msh.coord(0,1) = 0.10;
  msh.coord(0,2) = 0.02;
  msh.coord(1,0) = 1.10;
  msh.coord(1,1) = 0.16;
  msh.coord(1,2) = 0.12;
  msh.coord(2,0) = 0.22;
  msh.coord(2,1) = 0.94;
  msh.coord(2,2) = 0.18;
  for(int inode = 0; inode < 3; inode++){
    msh.fac2poi(0,inode) = inode;
  }

  const double nodal_metrics[3][6] = {
    {1.45, 0.05,1.02, 0.02,-0.03,0.88},
    {1.95, 0.09,1.25,-0.04, 0.02,0.96},
    {1.10,-0.04,1.72, 0.03, 0.05,1.28}
  };
  for(int ipoin = 0; ipoin < 3; ipoin++){
    for(int imetric = 0; imetric < 6; imetric++){
      msh.met(ipoin,imetric) = nodal_metrics[ipoin][imetric];
    }
  }
}

// Advertise CAD presence without supplying usable CAD entities. Objective
// integration must not inspect this token; any legacy nordev evaluation would.
class ScopedCadPresence
{
public:
  explicit ScopedCadPresence(CADInfo &cad) : cad_(cad)
  {
    METRIS_ENFORCE(cad_.EGADS_model == NULL);
    cad_.EGADS_model = reinterpret_cast<ego>(&presence_token_);
  }

  ~ScopedCadPresence()
  {
    cad_.EGADS_model = NULL;
  }

private:
  CADInfo &cad_;
  char presence_token_ = 0;
};

template<QuaFun iquaf>
void check_objective_ignores_surface_quality_parameters(
    Mesh<MetricFieldFE> &msh)
{
  constexpr int gdim = 3;
  constexpr int tdim = 2;
  constexpr int nhessian = gdim*(gdim + 1)/2;

  msh.param->qua_surf_wt_quality = 0.75;
  msh.param->qua_surf_wt_normal = 1.25;

  const double baseline_value
      = metqua<MetricFieldFE,gdim,tdim,iquaf,double>(
            msh,AsDeg::P1,AsDeg::P1,0,1.0);
  double baseline_gradient[gdim];
  double baseline_hessian[nhessian];
  const double baseline_derivative_value
      = d_metqua<MetricFieldFE,gdim,tdim,iquaf,double>(
            msh,AsDeg::P1,AsDeg::P1,0,0,
            FEBasis::Lagrange,DifVar::None,
            baseline_gradient,baseline_hessian,1.0);
  const double baseline_value_only_derivative_entry
      = d_metqua<MetricFieldFE,gdim,tdim,iquaf,double>(
            msh,AsDeg::P1,AsDeg::P1,0,-1,
            FEBasis::Lagrange,DifVar::None,NULL,NULL,1.0);

  msh.param->qua_surf_wt_quality = 17.0;
  msh.param->qua_surf_wt_normal = 1.e40;

  const double changed_value
      = metqua<MetricFieldFE,gdim,tdim,iquaf,double>(
            msh,AsDeg::P1,AsDeg::P1,0,1.0);
  double changed_gradient[gdim];
  double changed_hessian[nhessian];
  const double changed_derivative_value
      = d_metqua<MetricFieldFE,gdim,tdim,iquaf,double>(
            msh,AsDeg::P1,AsDeg::P1,0,0,
            FEBasis::Lagrange,DifVar::None,
            changed_gradient,changed_hessian,1.0);
  const double changed_value_only_derivative_entry
      = d_metqua<MetricFieldFE,gdim,tdim,iquaf,double>(
            msh,AsDeg::P1,AsDeg::P1,0,-1,
            FEBasis::Lagrange,DifVar::None,NULL,NULL,1.0);

  BOOST_CHECK_EQUAL(changed_value,baseline_value);
  BOOST_CHECK_EQUAL(changed_derivative_value,baseline_derivative_value);
  BOOST_CHECK_EQUAL(
      changed_value_only_derivative_entry,
      baseline_value_only_derivative_entry);
  BOOST_CHECK_EQUAL(baseline_derivative_value,baseline_value);
  BOOST_CHECK_EQUAL(baseline_value_only_derivative_entry,baseline_value);
}

template<int gdim>
struct FrozenQuadratureSamples
{
  static constexpr int nquad = gdim + 2;
  static constexpr int nnmet = gdim*(gdim + 1)/2;

  std::array<std::array<double,gdim + 1>,nquad> barycentric_points{};
  std::array<std::array<double,nnmet>,nquad> metrics{};
  std::array<double,nquad> weights{};
};

template<class MFT, int gdim>
FrozenQuadratureSamples<gdim> capture_frozen_samples(Mesh<MFT> &msh)
{
  FrozenQuadratureSamples<gdim> samples;
  const SimplexQuadratureView<gdim> quadrature
      = get_vertex_barycenter_quadrature<gdim>();
  const int *nodes = msh.ent2poi(gdim)[0];

  double measure = 1.0;
  #ifdef TESTQUALITYALGO
  BOOST_REQUIRE((isvalideltP1<gdim,gdim>(msh,0,NULL,&measure)));
  #endif

  for(int iquad = 0; iquad < quadrature.size(); iquad++){
    const SimplexQuadraturePointView<gdim> quadrature_point
        = quadrature[iquad];
    for(int ibary = 0; ibary < gdim + 1; ibary++){
      samples.barycentric_points[iquad][ibary]
          = quadrature_point.bary[ibary];
    }

    double *metric = samples.metrics[iquad].data();
    if(iquad < gdim + 1){
      const int ipoin = nodes[iquad];
      for(int imetric = 0;
          imetric < FrozenQuadratureSamples<gdim>::nnmet;
          imetric++){
        metric[imetric] = msh.met(ipoin,imetric);
      }
    }else if constexpr(
        std::is_same<MFT,MetricFieldAnalytical>::value){
      double physical_point[gdim];
      if constexpr(gdim == 2){
        eval2<gdim,1>(msh.coord,nodes,msh.getBasis(),
                      DifVar::None,DifVar::None,
                      quadrature_point.bary,physical_point,NULL,NULL);
      }else{
        eval3<gdim,1>(msh.coord,nodes,msh.getBasis(),
                      DifVar::None,DifVar::None,
                      quadrature_point.bary,physical_point,NULL,NULL);
      }
      msh.met.getMetPhys(DifVar::None,MetSpace::Exp,
                         physical_point,metric,NULL);
    }else{
      msh.met.getMetBary(AsDeg::P1,DifVar::None,MetSpace::Exp,
                         nodes,gdim,quadrature_point.bary,metric,NULL);
    }

    double weight = quadrature_point.weight*measure;
    #ifdef TESTQUALITYALGO
    #ifdef INTQUALINRIEMSPACE
    weight *= std::sqrt(detsym<gdim>(metric));
    #endif
    #endif
    samples.weights[iquad] = weight;
  }
  return samples;
}

template<class MFT, int gdim>
double evaluate_frozen_value(
    Mesh<MFT> &msh,
    const FrozenQuadratureSamples<gdim> &samples)
{
  const int *nodes = msh.ent2poi(gdim)[0];
  double value = 0.0;
  for(int iquad = 0;
      iquad < FrozenQuadratureSamples<gdim>::nquad;
      iquad++){
    value += samples.weights[iquad]
           * quafun_sizeshape<MFT,gdim,gdim,double>(
                 msh,AsDeg::P1,AsDeg::P1,nodes,
                 samples.barycentric_points[iquad].data(),
                 samples.metrics[iquad].data());
  }
  return value;
}

template<class MFT, int gdim>
double evaluate_frozen_derivatives(
    Mesh<MFT> &msh,
    const FrozenQuadratureSamples<gdim> &samples,
    int ivar,
    double *gradient,
    double *hessian)
{
  constexpr int nhessian = gdim*(gdim + 1)/2;
  const int *nodes = msh.ent2poi(gdim)[0];
  for(int icomponent = 0; icomponent < gdim; icomponent++){
    gradient[icomponent] = 0.0;
  }
  if(hessian != NULL){
    for(int ihessian = 0; ihessian < nhessian; ihessian++){
      hessian[ihessian] = 0.0;
    }
  }

  double value = 0.0;
  for(int iquad = 0;
      iquad < FrozenQuadratureSamples<gdim>::nquad;
      iquad++){
    double point_gradient[gdim];
    double point_hessian[nhessian];
    const double point_value
        = d_quafun_sizeshape<MFT,gdim,gdim,double>(
              msh,AsDeg::P1,AsDeg::P1,nodes,
              samples.barycentric_points[iquad].data(),
              samples.metrics[iquad].data(),
              ivar,msh.getBasis(),DifVar::None,
              point_gradient,hessian == NULL ? NULL : point_hessian);
    const double weight = samples.weights[iquad];
    value += weight*point_value;
    for(int icomponent = 0; icomponent < gdim; icomponent++){
      gradient[icomponent] += weight*point_gradient[icomponent];
    }
    if(hessian == NULL) continue;
    for(int ihessian = 0; ihessian < nhessian; ihessian++){
      hessian[ihessian] += weight*point_hessian[ihessian];
    }
  }
  return value;
}

template<class MFT, int gdim>
void check_integrated_derivatives(Mesh<MFT> &msh)
{
  constexpr int nhessian = gdim*(gdim + 1)/2;
  const FrozenQuadratureSamples<gdim> samples
      = capture_frozen_samples<MFT,gdim>(msh);
  const double integrated_value
      = metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
            msh,AsDeg::P1,AsDeg::P1,0,1.0);
  const double frozen_value
      = evaluate_frozen_value<MFT,gdim>(msh,samples);
  BOOST_CHECK_CLOSE_FRACTION(integrated_value,frozen_value,3.0e-14);

  const double finite_difference_step = 2.0e-6;
  for(int ivar = 0; ivar < gdim + 1; ivar++){
    double integrated_gradient[gdim];
    double integrated_hessian[nhessian];
    const double differentiated_value
        = d_metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
              msh,AsDeg::P1,AsDeg::P1,0,ivar,
              msh.getBasis(),DifVar::Phys,
              integrated_gradient,integrated_hessian,-5.0);
    BOOST_CHECK_CLOSE_FRACTION(
        differentiated_value,integrated_value,3.0e-14);

    double gradient_only[gdim];
    const double gradient_only_value
        = d_metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
              msh,AsDeg::P1,AsDeg::P1,0,ivar,
              msh.getBasis(),DifVar::Phys,
              gradient_only,NULL,4.0);
    BOOST_CHECK_CLOSE_FRACTION(
        gradient_only_value,integrated_value,3.0e-14);
    for(int icomponent = 0; icomponent < gdim; icomponent++){
      BOOST_CHECK_EQUAL(
          gradient_only[icomponent],integrated_gradient[icomponent]);
    }

    double reconstructed_gradient[gdim];
    double reconstructed_hessian[nhessian];
    const double reconstructed_value
        = evaluate_frozen_derivatives<MFT,gdim>(
              msh,samples,ivar,
              reconstructed_gradient,reconstructed_hessian);
    BOOST_CHECK_CLOSE_FRACTION(
        reconstructed_value,integrated_value,3.0e-14);
    for(int icomponent = 0; icomponent < gdim; icomponent++){
      BOOST_CHECK_CLOSE_FRACTION(
          integrated_gradient[icomponent],
          reconstructed_gradient[icomponent],3.0e-14);
    }
    for(int ihessian = 0; ihessian < nhessian; ihessian++){
      BOOST_CHECK_CLOSE_FRACTION(
          integrated_hessian[ihessian],
          reconstructed_hessian[ihessian],3.0e-14);
    }

    const int ipoin = msh.ent2poi(gdim)(0,ivar);
    for(int icomponent = 0; icomponent < gdim; icomponent++){
      const double coordinate = msh.coord(ipoin,icomponent);
      msh.coord(ipoin,icomponent)
          = coordinate + finite_difference_step;
      const double plus_value
          = evaluate_frozen_value<MFT,gdim>(msh,samples);
      msh.coord(ipoin,icomponent)
          = coordinate - finite_difference_step;
      const double minus_value
          = evaluate_frozen_value<MFT,gdim>(msh,samples);
      msh.coord(ipoin,icomponent) = coordinate;

      const double finite_difference_gradient
          = (plus_value - minus_value)/(2.0*finite_difference_step);
      const double gradient_error
          = std::abs(integrated_gradient[icomponent]
                     - finite_difference_gradient);
      BOOST_CHECK_SMALL(
          gradient_error,
          8.0e-8*(1.0 + std::abs(finite_difference_gradient)));
    }

    for(int jcomponent = 0; jcomponent < gdim; jcomponent++){
      const double coordinate = msh.coord(ipoin,jcomponent);
      double plus_gradient[gdim];
      double minus_gradient[gdim];

      msh.coord(ipoin,jcomponent)
          = coordinate + finite_difference_step;
      evaluate_frozen_derivatives<MFT,gdim>(
          msh,samples,ivar,plus_gradient,NULL);
      msh.coord(ipoin,jcomponent)
          = coordinate - finite_difference_step;
      evaluate_frozen_derivatives<MFT,gdim>(
          msh,samples,ivar,minus_gradient,NULL);
      msh.coord(ipoin,jcomponent) = coordinate;

      for(int icomponent = 0;
          icomponent <= jcomponent;
          icomponent++){
        const double finite_difference_hessian
            = (plus_gradient[icomponent] - minus_gradient[icomponent])
             /(2.0*finite_difference_step);
        const double hessian_error
            = std::abs(
                integrated_hessian[sym2idx(icomponent,jcomponent)]
                - finite_difference_hessian);
        BOOST_CHECK_SMALL(
            hessian_error,
            1.0e-6*(1.0 + std::abs(finite_difference_hessian)));
      }
    }
  }

  const double value_only_derivative_entry
      = d_metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
            msh,AsDeg::P1,AsDeg::P1,0,-1,
            msh.getBasis(),DifVar::Phys,NULL,NULL,9.0);
  BOOST_CHECK_EQUAL(value_only_derivative_entry,integrated_value);

  #ifndef TESTQUALITYALGO
  double baseline_pnorm_gradient[gdim];
  double baseline_pnorm_hessian[nhessian];
  const double baseline_pnorm_value
      = d_metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
            msh,AsDeg::P1,AsDeg::P1,0,0,
            msh.getBasis(),DifVar::None,
            baseline_pnorm_gradient,baseline_pnorm_hessian,1.0);
  msh.param->opt_pnorm = 3;
  double changed_pnorm_gradient[gdim];
  double changed_pnorm_hessian[nhessian];
  const double changed_pnorm_value
      = d_metqua<MFT,gdim,gdim,QuaFun::SizeShape,double>(
            msh,AsDeg::P1,AsDeg::P1,0,0,
            msh.getBasis(),DifVar::None,
            changed_pnorm_gradient,changed_pnorm_hessian,1.0);
  BOOST_CHECK_CLOSE_FRACTION(
      changed_pnorm_value,baseline_pnorm_value,3.0e-14);
  for(int icomponent = 0; icomponent < gdim; icomponent++){
    BOOST_CHECK_EQUAL(
        changed_pnorm_gradient[icomponent],
        baseline_pnorm_gradient[icomponent]);
  }
  for(int ihessian = 0; ihessian < nhessian; ihessian++){
    BOOST_CHECK_EQUAL(
        changed_pnorm_hessian[ihessian],
        baseline_pnorm_hessian[ihessian]);
  }
  msh.param->opt_pnorm = 1;
  #endif
}

template<class MFT, int gdim>
void run_integrated_derivative_case()
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_p = 1.5;
  parameters.opt_pnorm = 1;
  parameters.opt_power = -1;

  Mesh<MFT> msh;
  initialize_element<MFT,gdim>(msh,parameters);
  check_integrated_derivatives<MFT,gdim>(msh);
}

} // namespace

BOOST_AUTO_TEST_CASE(test_fe_triangle_integrated_derivatives)
{
  run_integrated_derivative_case<MetricFieldFE,2>();
}

BOOST_AUTO_TEST_CASE(test_fe_tetrahedron_integrated_derivatives)
{
  run_integrated_derivative_case<MetricFieldFE,3>();
}

BOOST_AUTO_TEST_CASE(test_analytical_triangle_integrated_derivatives)
{
  run_integrated_derivative_case<MetricFieldAnalytical,2>();
}

BOOST_AUTO_TEST_CASE(test_analytical_tetrahedron_integrated_derivatives)
{
  run_integrated_derivative_case<MetricFieldAnalytical,3>();
}

BOOST_AUTO_TEST_CASE(test_objectives_exclude_cad_normal_deviation)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  parameters.objective_p = 1.5;
  parameters.step_distance_barrier_rho0 = 2.0;
  parameters.step_distance_barrier_beta = 0.2;

  Mesh<MetricFieldFE> msh;
  initialize_embedded_surface_triangle(msh,parameters);
  ScopedCadPresence cad_presence(msh.CAD);

  BOOST_TEST_CONTEXT("SizeShape")
  {
    check_objective_ignores_surface_quality_parameters<QuaFun::SizeShape>(msh);
  }

  parameters.step_distance_shape_volume = false;
  parameters.step_distance_cavity_target_average = false;
  BOOST_TEST_CONTEXT("StepDistance")
  {
    check_objective_ignores_surface_quality_parameters<
        QuaFun::StepDistance>(msh);
  }

  parameters.step_distance_cavity_target_average = true;
  BOOST_TEST_CONTEXT("StepDistance cavity target average")
  {
    check_objective_ignores_surface_quality_parameters<
        QuaFun::StepDistance>(msh);
  }

  parameters.step_distance_cavity_target_average = false;
  parameters.step_distance_shape_volume = true;
  BOOST_TEST_CONTEXT("StepDistance ShapeVolume")
  {
    check_objective_ignores_surface_quality_parameters<
        QuaFun::StepDistance>(msh);
  }
}

} // namespace Metris
