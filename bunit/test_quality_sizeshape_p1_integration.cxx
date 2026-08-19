// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_quality_sizeshape_p1_integration

#include <boost/test/included/unit_test.hpp>

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "linalg/det.hxx"
#include "low_geo/measure.hxx"
#include "quality/low_metqua.hxx"
#include "quality/msh_metqua.hxx"
#include "quality/quafun_tradet.hxx"
#include "quality/simplex_quadrature.hxx"

#include <algorithm>
#include <cmath>

namespace Metris
{

namespace
{

void initialize_triangle(Mesh<MetricFieldFE> &msh,
                         MetrisParameters &parameters)
{
  msh.idim = 2;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &parameters;
  msh.set_npoin(3);
  msh.set_nface(1);
  msh.met.forceBasisFlag(FEBasis::Lagrange);
  msh.met.forceSpaceFlag(MetSpace::Exp);

  msh.coord(0,0) = 0.0;
  msh.coord(0,1) = 0.0;
  msh.coord(1,0) = 1.1;
  msh.coord(1,1) = 0.1;
  msh.coord(2,0) = 0.2;
  msh.coord(2,1) = 0.8;
  for(int inode = 0; inode < 3; inode++){
    msh.fac2poi(0,inode) = inode;
  }

  const double nodal_metrics[3][3] = {
    {1.5, 0.10,0.8},
    {2.0, 0.15,1.1},
    {0.9,-0.05,1.8}
  };
  for(int ipoin = 0; ipoin < 3; ipoin++){
    for(int imetric = 0; imetric < 3; imetric++){
      msh.met(ipoin,imetric) = nodal_metrics[ipoin][imetric];
    }
  }
}

void initialize_tetrahedron(Mesh<MetricFieldFE> &msh,
                            MetrisParameters &parameters)
{
  msh.idim = 3;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.param = &parameters;
  msh.set_npoin(4);
  msh.set_nelem(1);
  msh.met.forceBasisFlag(FEBasis::Lagrange);
  msh.met.forceSpaceFlag(MetSpace::Exp);

  msh.coord(0,0) = 0.0;
  msh.coord(0,1) = 0.0;
  msh.coord(0,2) = 0.0;
  msh.coord(1,0) = 1.0;
  msh.coord(1,1) = 0.1;
  msh.coord(1,2) = 0.0;
  msh.coord(2,0) = 0.2;
  msh.coord(2,1) = 0.9;
  msh.coord(2,2) = 0.1;
  msh.coord(3,0) = 0.1;
  msh.coord(3,1) = 0.2;
  msh.coord(3,2) = 0.8;
  for(int inode = 0; inode < 4; inode++){
    msh.tet2poi(0,inode) = inode;
  }

  const double nodal_metrics[4][6] = {
    {1.5, 0.05,1.0, 0.02,-0.03,0.8},
    {2.0, 0.10,1.2,-0.04, 0.02,0.9},
    {1.1,-0.05,1.8, 0.03, 0.06,1.3},
    {0.9, 0.02,1.1,-0.01, 0.04,2.1}
  };
  for(int ipoin = 0; ipoin < 4; ipoin++){
    for(int imetric = 0; imetric < 6; imetric++){
      msh.met(ipoin,imetric) = nodal_metrics[ipoin][imetric];
    }
  }
}

template<int gdim, int tdim>
double reconstruct_vertex_barycenter_sizeshape(
    Mesh<MetricFieldFE> &msh,
    double objective_p)
{
  constexpr int nnmet = gdim*(gdim + 1)/2;
  const int *nodes = msh.ent2poi(tdim)[0];
  const SimplexQuadratureView<tdim> quadrature
      = get_vertex_barycenter_quadrature<tdim>();

  double measure = 1.0;
  #ifdef TESTQUALITYALGO
  BOOST_REQUIRE((isvalideltP1<gdim,tdim>(msh,0,NULL,&measure)));
  #endif

  double expected = 0.0;
  for(int iquad = 0; iquad < quadrature.size(); iquad++){
    const SimplexQuadraturePointView<tdim> quadrature_point
        = quadrature[iquad];
    double metric[nnmet];
    if(iquad < tdim + 1){
      const int ipoin = nodes[iquad];
      for(int imetric = 0; imetric < nnmet; imetric++){
        metric[imetric] = msh.met(ipoin,imetric);
      }
    }else{
      msh.met.getMetBary(AsDeg::P1,DifVar::None,MetSpace::Exp,
                         nodes,tdim,quadrature_point.bary,
                         metric,NULL);
    }

    double trace;
    double determinant;
    quafun_tradet<MetricFieldFE,gdim,tdim,double>(
        msh,AsDeg::P1,AsDeg::P1,nodes,quadrature_point.bary,metric,
        &trace,&determinant,QualitySingularityPolicy::Reject);

    double size_shape_quality;
    if constexpr(tdim == 2){
      size_shape_quality
          = trace*trace*(1.0 + 1.0/(determinant*determinant))/8.0;
    }else{
      size_shape_quality
          = trace*trace*trace
           *(1.0 + 1.0/(determinant*determinant))/54.0;
    }
    const double size_shape_error
        = std::max(0.0,size_shape_quality - 1.0);
    const double psi = objective_p == 1.0
                     ? size_shape_error
                     : std::pow(size_shape_error,objective_p);

    double weight = quadrature_point.weight*measure;
    #ifdef TESTQUALITYALGO
    #ifdef INTQUALINRIEMSPACE
    weight *= std::sqrt(detsym<gdim>(metric));
    #endif
    #endif
    expected += weight*psi;
  }
  return expected;
}

template<int gdim, int tdim>
void check_p1_sizeshape_integration(Mesh<MetricFieldFE> &msh)
{
  msh.param->objective_p = 1.0;
  msh.param->opt_power = 1;
  msh.param->opt_pnorm = 1;

  const double historical_direct_value
      = reconstruct_vertex_barycenter_sizeshape<gdim,tdim>(msh,1.0);
  const double integrated_value
      = metqua<MetricFieldFE,gdim,tdim,QuaFun::SizeShape,double>(
          msh,AsDeg::P1,AsDeg::P1,0,1.0);
  BOOST_CHECK_CLOSE_FRACTION(
      integrated_value,historical_direct_value,2.0e-14);

  bool invalid_mesh = false;
  double mesh_minimum;
  double mesh_maximum;
  double mesh_average;
  const double mesh_objective_pnorm_one
      = getmetquamesh<MetricFieldFE,QuaFun::SizeShape>(
          msh,tdim,AsDeg::P1,AsDeg::P1,
          &invalid_mesh,&mesh_minimum,&mesh_maximum,&mesh_average,NULL);
  msh.param->opt_pnorm = 3;
  const double mesh_objective_pnorm_three
      = getmetquamesh<MetricFieldFE,QuaFun::SizeShape>(
          msh,tdim,AsDeg::P1,AsDeg::P1,
          &invalid_mesh,&mesh_minimum,&mesh_maximum,&mesh_average,NULL);
  BOOST_CHECK_CLOSE_FRACTION(
      mesh_objective_pnorm_one,integrated_value,2.0e-14);
  BOOST_CHECK_CLOSE_FRACTION(
      mesh_objective_pnorm_three,integrated_value,2.0e-14);
  msh.param->opt_pnorm = 1;

  const double different_difto_value
      = metqua<MetricFieldFE,gdim,tdim,QuaFun::SizeShape,double>(
          msh,AsDeg::P1,AsDeg::P1,0,-7.0);
  BOOST_CHECK_EQUAL(different_difto_value,integrated_value);

  msh.param->opt_power = -1;
  const double inverse_setting_value
      = metqua<MetricFieldFE,gdim,tdim,QuaFun::SizeShape,double>(
          msh,AsDeg::P1,AsDeg::P1,0,3.0);
  BOOST_CHECK_EQUAL(inverse_setting_value,integrated_value);

  #ifndef TESTQUALITYALGO
  msh.param->opt_pnorm = 3;
  const double different_pnorm_value
      = metqua<MetricFieldFE,gdim,tdim,QuaFun::SizeShape,double>(
          msh,AsDeg::P1,AsDeg::P1,0,1.0);
  BOOST_CHECK_EQUAL(different_pnorm_value,integrated_value);
  #endif

  msh.param->opt_pnorm = 1;
  msh.param->objective_p = 2.0;
  const double expected_squared_objective
      = reconstruct_vertex_barycenter_sizeshape<gdim,tdim>(msh,2.0);
  const double squared_objective
      = metqua<MetricFieldFE,gdim,tdim,QuaFun::SizeShape,double>(
          msh,AsDeg::P1,AsDeg::P1,0,1.0);
  BOOST_CHECK_CLOSE_FRACTION(
      squared_objective,expected_squared_objective,2.0e-14);
}

} // namespace

BOOST_AUTO_TEST_CASE(test_triangle_p1_sizeshape_integration)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  Mesh<MetricFieldFE> msh;
  initialize_triangle(msh,parameters);
  check_p1_sizeshape_integration<2,2>(msh);
}

BOOST_AUTO_TEST_CASE(test_tetrahedron_p1_sizeshape_integration)
{
  MetrisParameters parameters;
  parameters.iverb = 0;
  Mesh<MetricFieldFE> msh;
  initialize_tetrahedron(msh,parameters);
  check_p1_sizeshape_integration<3,3>(msh);
}

} // namespace Metris
