// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2026, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_metric_discontinuous_storage
#include <boost/test/included/unit_test.hpp>

#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "Adaptation/Insertion/msh_insert.hxx"
#include "Adaptation/msh_collapse.hxx"
#include "Adaptation/msh_swap.hxx"
#include "msh_degelev.hxx"
#include "smoothing/msh_smooball.hxx"
#include "smoothing/msh_smoolen.hxx"

#include <utility>

namespace Metris{
namespace{

struct TwoTriangles{
  MetrisParameters parameters;
  Mesh<MetricFieldFE> mesh;

  TwoTriangles(){
    parameters.iverb = 0;
    mesh.param = &parameters;
    mesh.idim = 2;
    mesh.curdeg = mesh.strdeg = 2;
    mesh.forceBasisFlag(FEBasis::Lagrange);
    mesh.set_npoin(9);
    mesh.set_nentt(2,2);
    mesh.met.forceBasisFlag(FEBasis::Lagrange);
    const int nodes[2][6] = {{0,1,2,4,5,6},{1,3,2,7,4,8}};
    const double coordinates[9][2] = {
      {0.,0.},{1.,0.},{0.,1.},{1.,1.},
      {0.5,0.5},{0.,0.5},{0.5,0.},{0.5,1.},{1.,0.5}
    };
    for(int face = 0; face < 2; face++)
      for(int node = 0; node < 6; node++) mesh.fac2poi(face,node) = nodes[face][node];
    for(int point = 0; point < 9; point++){
      for(int k = 0; k < 2; k++) mesh.coord(point,k) = coordinates[point][k];
      mesh.met(point,0) = 2. + point;
      mesh.met(point,1) = 0.125;
      mesh.met(point,2) = 1. + point;
    }
  }
};

} // namespace

BOOST_FIXTURE_TEST_CASE(continuous_layout_preserves_shared_geometry_indexing,TwoTriangles){
  BOOST_CHECK(mesh.met.getLayout() == MetricFieldLayout::Continuous);
  BOOST_CHECK_EQUAL(mesh.met.getDegree(),2);
  BOOST_CHECK_EQUAL(mesh.met.getCoefficientCount(),9);
  for(int face = 0; face < 2; face++) for(int node = 0; node < 6; node++){
    BOOST_CHECK_EQUAL(mesh.met.getElementCoefficientIndex(face,node),mesh.fac2poi(face,node));
    BOOST_CHECK(mesh.met.elementCoefficient(face,node) == mesh.met[mesh.fac2poi(face,node)]);
  }
  BOOST_CHECK(mesh.met.elementCoefficient(0,1) == mesh.met.elementCoefficient(1,0));
  BOOST_CHECK(mesh.met.elementCoefficient(0,3) == mesh.met.elementCoefficient(1,4));
  MetricFieldFE copy(mesh.met);
  copy(1,0) = 42.;
  BOOST_CHECK_EQUAL(mesh.met(1,0),3.);
  // Continuous metric degree follows the existing geometry-degree convention.
  mesh.curdeg = 1;
  BOOST_CHECK_EQUAL(mesh.met.getDegree(),1);
}

BOOST_FIXTURE_TEST_CASE(discontinuous_coefficients_have_independent_vertex_ownership,TwoTriangles){
  mesh.met.initializeDiscontinuous(1);
  BOOST_CHECK(mesh.met.getLayout() == MetricFieldLayout::Discontinuous);
  BOOST_CHECK(mesh.met.getBasis() == FEBasis::Lagrange);
  BOOST_CHECK(mesh.met.getSpace() == MetSpace::Log);
  BOOST_CHECK_EQUAL(mesh.met.getDegree(),1);
  BOOST_CHECK_EQUAL(mesh.curdeg,2);
  BOOST_CHECK_EQUAL(mesh.met.getCoefficientCount(),6);
  BOOST_CHECK_EQUAL(mesh.npoin,9);
  BOOST_CHECK_EQUAL(mesh.fac2poi(0,1),mesh.fac2poi(1,0));
  BOOST_CHECK_EQUAL(mesh.fac2poi(0,3),mesh.fac2poi(1,4));
  for(int face = 0; face < 2; face++) for(int node = 0; node < 3; node++){
    for(int component = 0; component < 3; component++)
      BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(face,node)[component],0.);
  }
  BOOST_CHECK_NE(mesh.met.getElementCoefficientIndex(0,1),mesh.met.getElementCoefficientIndex(1,0));
  mesh.met.elementCoefficient(0,1)[0] = 0.7;
  mesh.met.elementCoefficient(1,0)[0] = -0.3;
  const MetricFieldFE &met = mesh.met;
  BOOST_CHECK_EQUAL(met.elementCoefficient(0,1)[0],0.7);
  BOOST_CHECK_EQUAL(met.elementCoefficient(1,0)[0],-0.3);
  // Geometry point storage can grow without resizing or aliasing DG coefficients.
  mesh.set_npoin(12);
  BOOST_CHECK_EQUAL(met.getCoefficientCount(),6);
  BOOST_CHECK_EQUAL(met.elementCoefficient(0,1)[0],0.7);
  BOOST_CHECK_EQUAL(met.elementCoefficient(1,0)[0],-0.3);
  BOOST_CHECK_THROW(met.elementCoefficient(0,3),MetrisExcept);
}

BOOST_FIXTURE_TEST_CASE(discontinuous_assignment_preserves_layout_and_owns_coefficients,TwoTriangles){
  TwoTriangles destination;
  mesh.met.initializeDiscontinuous();
  mesh.met.elementCoefficient(0,1)[0] = 0.7;
  mesh.met.elementCoefficient(1,0)[0] = -0.3;
  destination.mesh.met = mesh.met;
  BOOST_CHECK(destination.mesh.met.getLayout() == MetricFieldLayout::Discontinuous);
  BOOST_CHECK_EQUAL(destination.mesh.met.getDegree(),1);
  BOOST_CHECK_EQUAL(destination.mesh.met.getCoefficientCount(),6);
  BOOST_CHECK_EQUAL(destination.mesh.met.elementCoefficient(0,1)[0],0.7);
  BOOST_CHECK_EQUAL(destination.mesh.met.elementCoefficient(1,0)[0],-0.3);
  destination.mesh.met.elementCoefficient(0,1)[0] = 1.2;
  BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(0,1)[0],0.7);
  const MetricFieldFE &self = destination.mesh.met;
  destination.mesh.met = self;
  BOOST_CHECK_EQUAL(destination.mesh.met.elementCoefficient(0,1)[0],1.2);

  TwoTriangles reordered;
  std::swap(reordered.mesh.fac2poi(0,0),reordered.mesh.fac2poi(0,1));
  BOOST_CHECK_THROW(reordered.mesh.met = mesh.met,MetrisExcept);
  BOOST_CHECK(reordered.mesh.met.getLayout() == MetricFieldLayout::Continuous);

  TwoTriangles continuous;
  destination.mesh.met = continuous.mesh.met;
  BOOST_CHECK(destination.mesh.met.getLayout() == MetricFieldLayout::Continuous);
  BOOST_CHECK_EQUAL(destination.mesh.met.getDegree(),2);
  BOOST_CHECK_EQUAL(destination.mesh.met.getCoefficientCount(),9);
  BOOST_CHECK_EQUAL(destination.mesh.met(1,0),continuous.mesh.met(1,0));
}

BOOST_FIXTURE_TEST_CASE(dead_elements_have_no_discontinuous_coefficients,TwoTriangles){
  mesh.set_nentt(2,3);
  mesh.fac2poi(2,0) = -1;
  mesh.met.initializeDiscontinuous();
  BOOST_CHECK_EQUAL(mesh.met.getCoefficientCount(),6);
  BOOST_CHECK_THROW(mesh.met.elementCoefficient(2,0),MetrisExcept);
  BOOST_CHECK_THROW(mesh.met.elementCoefficient(-1,0),MetrisExcept);
  BOOST_CHECK_THROW(mesh.met.elementCoefficient(0,-1),MetrisExcept);
  BOOST_CHECK_THROW(mesh.met.elementCoefficient(3,0),MetrisExcept);
}

BOOST_FIXTURE_TEST_CASE(discontinuous_storage_rejects_ambiguous_legacy_operations,TwoTriangles){
  BOOST_CHECK_THROW(mesh.met.initializeDiscontinuous(2),MetrisExcept);
  BOOST_CHECK(mesh.met.getLayout() == MetricFieldLayout::Continuous);
  BOOST_CHECK_EQUAL(mesh.met(1,0),3.);
  mesh.met.initializeDiscontinuous();
  BOOST_CHECK_THROW(mesh.met.initializeDiscontinuous(),MetrisExcept);
  BOOST_CHECK_THROW(mesh.met[1],MetrisExcept);
  BOOST_CHECK_THROW(mesh.met(1,0),MetrisExcept);
  BOOST_CHECK_THROW(mesh.met.setSpace(MetSpace::Undefined),MetrisExcept);
  BOOST_CHECK_THROW(mesh.met.setBasis(FEBasis::Undefined),MetrisExcept);
  BOOST_CHECK_THROW(mesh.met.forceSpaceFlag(MetSpace::Exp),MetrisExcept);
  BOOST_CHECK_THROW(mesh.met.forceBasisFlag(FEBasis::Bezier),MetrisExcept);
  BOOST_CHECK_THROW(mesh.met.normalize(0.),MetrisExcept);
  const double bary[3] = {0.2,0.3,0.5};
  double met[3];
  BOOST_CHECK_THROW(mesh.met.getMetBary(AsDeg::Pk,DifVar::None,MetSpace::Log,
                                      mesh.fac2poi[0],2,bary,met,nullptr),MetrisExcept);
}

BOOST_AUTO_TEST_CASE(copying_an_empty_metric_preserves_its_state){
  Mesh<MetricFieldFE> mesh;
  const MetricFieldFE copy(mesh.met);
  BOOST_CHECK(copy.getLayout() == MetricFieldLayout::Continuous);
  BOOST_CHECK_EQUAL(copy.getCoefficientCount(),0);
  BOOST_CHECK(copy.getBasis() == FEBasis::Undefined);
}

BOOST_FIXTURE_TEST_CASE(dg_rejects_unsupported_support_before_modifying_coefficients,TwoTriangles){
  for(int degree : {1,3}){
    mesh.curdeg = degree;
    BOOST_CHECK_THROW(mesh.met.initializeDiscontinuous(),MetrisExcept);
    BOOST_CHECK(mesh.met.getLayout() == MetricFieldLayout::Continuous);
  }
  mesh.curdeg = 2;
  mesh.idim = 3;
  BOOST_CHECK_THROW(mesh.met.initializeDiscontinuous(),MetrisExcept);
  mesh.idim = 2;
  mesh.forceBasisFlag(FEBasis::Undefined);
  BOOST_CHECK_THROW(mesh.met.initializeDiscontinuous(),MetrisExcept);
  mesh.forceBasisFlag(FEBasis::Lagrange);
  mesh.met.initializeDiscontinuous();
  mesh.met.elementCoefficient(0,0)[0] = 0.7;
  for(int degree : {1,3}){
    mesh.curdeg = degree;
    BOOST_CHECK_THROW(mesh.met.setSpace(MetSpace::Exp),MetrisExcept);
    BOOST_CHECK_THROW(mesh.met.normalize(2.),MetrisExcept);
    BOOST_CHECK_THROW(mesh.met.setBasis(FEBasis::Bezier),MetrisExcept);
    BOOST_CHECK_THROW(mesh.met.elementCoefficient(0,0),MetrisExcept);
  }
  mesh.curdeg = 2;
  mesh.idim = 3;
  BOOST_CHECK_THROW(mesh.met.validateDiscontinuous(),MetrisExcept);
  mesh.idim = 2;
  mesh.set_nentt(2,0);
  BOOST_CHECK_THROW(mesh.met.validateDiscontinuous(),MetrisExcept);
  mesh.set_nentt(2,2);
  BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(0,0)[0],0.7);
  BOOST_CHECK(mesh.met.getSpace() == MetSpace::Log);
  BOOST_CHECK(mesh.met.getBasis() == FEBasis::Lagrange);
}

BOOST_FIXTURE_TEST_CASE(dg_smoothing_rejects_wrong_objectives_and_geometry_basis,TwoTriangles){
  // Continuous fields retain their existing smoothing choices.
  BOOST_CHECK_NO_THROW(validateSmoothingMetric(mesh,mesh.met,QuaFun::StepDistance));
  mesh.met.initializeDiscontinuous();
  for(MetSpace space : {MetSpace::Log,MetSpace::Exp}){
    mesh.met.setSpace(space);
    for(FEBasis basis : {FEBasis::Lagrange,FEBasis::Bezier}){
      mesh.met.setBasis(basis);
      BOOST_CHECK_NO_THROW(validateSmoothingMetric(mesh,mesh.met,QuaFun::Connection));
    }
  }
  const auto wrongObjective = [](const MetrisExcept &error){
    return std::string(error.what()).find("DG smoothing requires the connection objective") != std::string::npos;
  };
  for(QuaFun objective : {QuaFun::StepDistance,QuaFun::SizeShape,QuaFun::Distortion}){
    BOOST_CHECK_EXCEPTION(smoothInterior_Ball(mesh,objective,0,1),MetrisExcept,wrongObjective);
    BOOST_CHECK_EXCEPTION((smoothInterior_Ball0<MetricFieldFE,2,2>(mesh,objective,0,1)),MetrisExcept,wrongObjective);
  }
  mesh.forceBasisFlag(FEBasis::Bezier);
  BOOST_CHECK_NO_THROW(mesh.met.validateDiscontinuous());
  BOOST_CHECK_THROW(smoothInterior_Ball(mesh,QuaFun::Connection,0,1),MetrisExcept);
  mesh.forceBasisFlag(FEBasis::Lagrange);
  BOOST_CHECK_EQUAL(mesh.coord(4,0),0.5);
  BOOST_CHECK_EQUAL(mesh.coord(4,1),0.5);
}

BOOST_FIXTURE_TEST_CASE(dg_topology_operations_reject_before_mutating_mesh,TwoTriangles){
  mesh.met.initializeDiscontinuous();
  mesh.met.elementCoefficient(0,1)[0] = 0.7;
  int count = 17;
  const int tag = mesh.tag[0];
  const int stride = mesh.fac2poi.get_stride();
  const auto requiresContinuous = [](const MetrisExcept &error){
    return std::string(error.what()).find("requires a continuous metric layout") != std::string::npos;
  };
  BOOST_CHECK_EXCEPTION((insertLongEdges<MetricFieldFE,2,2>(mesh,2,&count,0,1)),MetrisExcept,requiresContinuous);
  BOOST_CHECK_EXCEPTION((collapseShortEdges<MetricFieldFE,2,2>(mesh,2,1.,&count,0,1,2,3)),MetrisExcept,requiresContinuous);
  BOOST_CHECK_EXCEPTION((swapMesh<MetricFieldFE,2,2>(mesh,swapOptions(parameters),&count,0,1,2)),MetrisExcept,requiresContinuous);
  BOOST_CHECK_EXCEPTION(mesh.cleanup(),MetrisExcept,requiresContinuous);
  BOOST_CHECK_EXCEPTION(mesh.newpoitopo(PointType::Vertex,2,0),MetrisExcept,requiresContinuous);
  BOOST_CHECK_EXCEPTION((deg_elevate<MetricFieldFE,1,2>(mesh)),MetrisExcept,requiresContinuous);
  BOOST_CHECK_EXCEPTION(smoothMeshLength(mesh,2,0,1),MetrisExcept,requiresContinuous);
  BOOST_CHECK_EQUAL(count,17);
  BOOST_CHECK_EQUAL(mesh.tag[0],tag);
  BOOST_CHECK_EQUAL(mesh.fac2poi.get_stride(),stride);
  BOOST_CHECK_EQUAL(mesh.nface,2);
  BOOST_CHECK_EQUAL(mesh.npoin,9);
  BOOST_CHECK_EQUAL(mesh.met.elementCoefficient(0,1)[0],0.7);
  BOOST_CHECK(mesh.met.getSpace() == MetSpace::Log);
  const int expected[2][6] = {{0,1,2,4,5,6},{1,3,2,7,4,8}};
  for(int ientt = 0; ientt < 2; ientt++) for(int node = 0; node < 6; node++)
    BOOST_CHECK_EQUAL(mesh.fac2poi(ientt,node),expected[ientt][node]);
}

} // namespace Metris
