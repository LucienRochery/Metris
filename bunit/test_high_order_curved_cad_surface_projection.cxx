// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_high_order_curved_cad_surface_projection

#include "common_setup.hxx"

#include "API/MetrisAPI.hxx"
#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "MetrisRunner/MetrisRunner.hxx"
#include "low_topo.hxx"

#include <cmath>
#include <filesystem>
#include <memory>

using namespace Metris;

namespace
{

struct CurvedSurfaceCase
{
  ego context = nullptr;
  ego model = nullptr;
  std::unique_ptr<MetrisRunner> runner;

  ~CurvedSurfaceCase()
  {
    runner.reset();
    if(context != nullptr) EG_close(context);
  }
};

void require_egads_success(int status, const char *operation)
{
  METRIS_ENFORCE_MSG(
      status == EGADS_SUCCESS,"{} failed with EGADS status {}",
      operation,status);
}

double squared_distance(const double *first, const double *second)
{
  double distance = 0.0;
  for(int component = 0; component < 3; component++){
    const double difference = first[component] - second[component];
    distance += difference*difference;
  }
  return distance;
}

int mesh_point_for_cad_node(ego cadNode,
                            const double coordinates[4][3])
{
  double nodeCoordinates[4];
  int objectClass = 0;
  int topologyType = 0;
  int childCount = 0;
  ego geometry = nullptr;
  ego *children = nullptr;
  int *senses = nullptr;
  require_egads_success(
      EG_getTopology(
          cadNode,&geometry,&objectClass,&topologyType,nodeCoordinates,
          &childCount,&children,&senses),
      "EG_getTopology(NODE)");
  BOOST_REQUIRE_EQUAL(objectClass,NODE);

  int closestPoint = -1;
  double closestDistance = 1.0;
  for(int point = 0; point < 4; point++){
    const double distance
        = squared_distance(nodeCoordinates,coordinates[point]);
    if(distance >= closestDistance) continue;
    closestDistance = distance;
    closestPoint = point;
  }
  BOOST_REQUIRE_GE(closestPoint,0);
  BOOST_REQUIRE_LT(closestDistance,1.e-24);
  return closestPoint;
}

std::unique_ptr<CurvedSurfaceCase> make_curved_surface_case()
{
  auto testCase = std::make_unique<CurvedSurfaceCase>();
  require_egads_success(EG_open(&testCase->context),"EG_open");
  EG_setOutLevel(testCase->context,0);

  double sphereData[10] = {
      0.0,0.0,0.0,
      1.0,0.0,0.0,
      0.0,1.0,0.0,
      1.0};
  ego sphere;
  require_egads_success(
      EG_makeGeometry(
          testCase->context,SURFACE,SPHERICAL,nullptr,nullptr,
          sphereData,&sphere),
      "EG_makeGeometry(SPHERICAL)");
  const double faceRange[4] = {0.2,1.4,-0.4,0.8};
  ego face;
  require_egads_success(
      EG_makeFace(sphere,SFORWARD,faceRange,&face),"EG_makeFace");
  ego body;
  require_egads_success(
      EG_makeTopology(
          testCase->context,nullptr,BODY,FACEBODY,nullptr,
          1,&face,nullptr,&body),
      "EG_makeTopology(BODY)");
  require_egads_success(
      EG_makeTopology(
          testCase->context,nullptr,MODEL,0,nullptr,
          1,&body,nullptr,&testCase->model),
      "EG_makeTopology(MODEL)");

  const double vertexParameters[4][2] = {
      {faceRange[0],faceRange[2]},
      {faceRange[1],faceRange[2]},
      {faceRange[1],faceRange[3]},
      {faceRange[0],faceRange[3]}};
  double coordinates[4][3];
  for(int point = 0; point < 4; point++){
    double evaluation[18];
    require_egads_success(
        EG_evaluate(face,vertexParameters[point],evaluation),
        "EG_evaluate(FACE vertex)");
    for(int component = 0; component < 3; component++){
      coordinates[point][component] = evaluation[component];
    }
  }

  MetrisAPI input(
      3,1,4,8,6,4,true,4,2,0,
      FEBasis::Lagrange,FEBasis::Lagrange,MetSpace::Exp);
  const double identityMetric[6] = {1.0,0.0,1.0,0.0,0.0,1.0};
  for(int point = 0; point < 4; point++){
    input.setCoord(point,coordinates[point]);
    input.setMetric(point,identityMetric);
  }

  int cadNodeCount = 0;
  ego *cadNodes = nullptr;
  require_egads_success(
      EG_getBodyTopos(body,nullptr,NODE,&cadNodeCount,&cadNodes),
      "EG_getBodyTopos(NODE)");
  BOOST_REQUIRE_EQUAL(cadNodeCount,4);
  for(int cadNode = 0; cadNode < cadNodeCount; cadNode++){
    const int point
        = mesh_point_for_cad_node(cadNodes[cadNode],coordinates);
    const int reference = EG_indexBodyTopo(body,cadNodes[cadNode]) - 1;
    BOOST_REQUIRE_GE(reference,0);
    input.setCorner(reference,point);
  }
  EG_free(cadNodes);

  constexpr int boundaryPointPairs[4][2]
      = {{0,1},{1,2},{2,3},{3,0}};
  int cadEdgeCount = 0;
  ego *cadEdges = nullptr;
  require_egads_success(
      EG_getBodyTopos(body,nullptr,EDGE,&cadEdgeCount,&cadEdges),
      "EG_getBodyTopos(EDGE)");
  BOOST_REQUIRE_EQUAL(cadEdgeCount,4);

  int geometricEdgePoints[8][2];
  double geometricEdgeParameters[8][2];
  for(int meshEdge = 0; meshEdge < 4; meshEdge++){
    const int firstPoint = boundaryPointPairs[meshEdge][0];
    const int secondPoint = boundaryPointPairs[meshEdge][1];
    int matchingCadEdge = -1;
    for(int cadEdge = 0; cadEdge < cadEdgeCount; cadEdge++){
      double firstParameter[2] = {0.0,0.0};
      double secondParameter[2] = {0.0,0.0};
      double firstEvaluation[18];
      double secondEvaluation[18];
      const int firstStatus = EG_invEvaluate(
          cadEdges[cadEdge],coordinates[firstPoint],
          firstParameter,firstEvaluation);
      const int secondStatus = EG_invEvaluate(
          cadEdges[cadEdge],coordinates[secondPoint],
          secondParameter,secondEvaluation);
      if(firstStatus != EGADS_SUCCESS
         || secondStatus != EGADS_SUCCESS) continue;
      if(squared_distance(firstEvaluation,coordinates[firstPoint]) > 1.e-24
         || squared_distance(secondEvaluation,coordinates[secondPoint])
                > 1.e-24) continue;
      matchingCadEdge = cadEdge;
      geometricEdgePoints[2*meshEdge][0] = firstPoint;
      geometricEdgePoints[2*meshEdge][1] = meshEdge;
      geometricEdgeParameters[2*meshEdge][0] = firstParameter[0];
      geometricEdgeParameters[2*meshEdge][1] = 0.0;
      geometricEdgePoints[2*meshEdge + 1][0] = secondPoint;
      geometricEdgePoints[2*meshEdge + 1][1] = meshEdge;
      geometricEdgeParameters[2*meshEdge + 1][0]
          = secondParameter[0];
      geometricEdgeParameters[2*meshEdge + 1][1] = 0.0;
      break;
    }
    BOOST_REQUIRE_GE(matchingCadEdge,0);
    const int reference
        = EG_indexBodyTopo(body,cadEdges[matchingCadEdge]) - 1;
    BOOST_REQUIRE_GE(reference,0);
    input.setElement(
        1,meshEdge,boundaryPointPairs[meshEdge],reference);
  }
  EG_free(cadEdges);
  input.setVerticesOnGeometricEdges(
      0,8,&geometricEdgePoints[0][0],
      &geometricEdgeParameters[0][0]);

  const int facePoints[2][3] = {{0,1,2},{0,2,3}};
  const int faceReference = EG_indexBodyTopo(body,face) - 1;
  BOOST_REQUIRE_GE(faceReference,0);
  input.setElement(2,0,facePoints[0],faceReference);
  input.setElement(2,1,facePoints[1],faceReference);

  const int geometricFacePoints[6][2]
      = {{0,0},{1,0},{2,0},{0,1},{2,1},{3,1}};
  double geometricFaceParameters[6][3];
  for(int record = 0; record < 6; record++){
    const int point = geometricFacePoints[record][0];
    geometricFaceParameters[record][0] = vertexParameters[point][0];
    geometricFaceParameters[record][1] = vertexParameters[point][1];
    geometricFaceParameters[record][2] = 0.0;
  }
  input.setVerticesOnGeometricTriangles(
      0,6,&geometricFacePoints[0][0],
      &geometricFaceParameters[0][0]);
  input.setCADModel(testCase->context,testCase->model);

  MetrisParameters parameters;
  parameters.usrTarDeg = 2;
  parameters.adp_niter = 0;
  parameters.opt_niter = 0;
  parameters.iverb = 0;
  parameters.outmPrefix
      = (std::filesystem::temp_directory_path()
         / "metris_high_order_curved_cad_surface_projection").string()
        + "/";
  testCase->runner = std::make_unique<MetrisRunner>(
      &input,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(testCase->runner->degElevate(),1);
  return testCase;
}

std::string curved_surface_output_directory(const char *stage)
{
  const std::filesystem::path directory
      = std::filesystem::temp_directory_path()
      / "metris_high_order_curved_cad_surface_projection"
      / stage;
  std::filesystem::create_directories(directory);
  return directory.string() + "/";
}

void check_curved_surface_node_invariants(Mesh<MetricFieldFE> &mesh)
{
  mesh.setBasis(FEBasis::Lagrange);
  BOOST_REQUIRE_EQUAL(mesh.curdeg,2);
  BOOST_REQUIRE(mesh.getBasis() == FEBasis::Lagrange);

  const int diagonalLocalEdge = getedgfac(mesh,0,0,2);
  BOOST_REQUIRE_GE(diagonalLocalEdge,0);
  const int surfacePoint = mesh.fac2poi(0,3 + diagonalLocalEdge);
  BOOST_REQUIRE_GE(surfacePoint,0);
  const int primaryRecord = mesh.poi2bpo[surfacePoint];
  BOOST_REQUIRE_GE(primaryRecord,0);
  BOOST_REQUIRE_EQUAL(mesh.bpo2ibi(primaryRecord,1),2);
  const int surfaceRecord = mesh.poi2ebp(surfacePoint,2,0,-1);
  BOOST_REQUIRE_GE(surfaceRecord,0);

  const double expectedParameters[2] = {0.8,0.2};
  BOOST_CHECK_SMALL(
      mesh.bpo2rbi(surfaceRecord,0) - expectedParameters[0],1.e-14);
  BOOST_CHECK_SMALL(
      mesh.bpo2rbi(surfaceRecord,1) - expectedParameters[1],1.e-14);

  const int reference = mesh.fac2ref[0];
  const ego cadFace = mesh.CAD.cad2fac[reference];
  double evaluation[18];
  require_egads_success(
      EG_evaluate(cadFace,expectedParameters,evaluation),
      "EG_evaluate(elevated surface point)");
  for(int component = 0; component < 3; component++){
    BOOST_CHECK_SMALL(
        mesh.coord(surfacePoint,component) - evaluation[component],1.e-13);
  }

  double chordMidpoint[3];
  for(int component = 0; component < 3; component++){
    chordMidpoint[component]
        = 0.5*(mesh.coord(0,component) + mesh.coord(2,component));
  }
  BOOST_CHECK_GT(
      std::sqrt(squared_distance(mesh.coord[surfacePoint],chordMidpoint)),
      0.1);
  BOOST_CHECK_GE(mesh.poi2bak[surfacePoint],0);

  const int curvePoint = mesh.edg2poi(0,2);
  const int curveRecord = mesh.poi2ebp(curvePoint,1,0,-1);
  BOOST_REQUIRE_GE(curveRecord,0);
  BOOST_REQUIRE_EQUAL(mesh.bpo2ibi(mesh.poi2bpo[curvePoint],1),1);
  const ego cadEdge = mesh.CAD.cad2edg[mesh.edg2ref[0]];
  const double curveParameter = mesh.bpo2rbi(curveRecord,0);
  require_egads_success(
      EG_evaluate(cadEdge,&curveParameter,evaluation),
      "EG_evaluate(elevated curve point)");
  for(int component = 0; component < 3; component++){
    BOOST_CHECK_SMALL(
        mesh.coord(curvePoint,component) - evaluation[component],1.e-13);
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(p2_internal_surface_edge_node_uses_elevated_uv)
{
  auto testCase = make_curved_surface_case();
  auto &mesh = static_cast<Mesh<MetricFieldFE>&>(
      *testCase->runner->msh_g);
  check_curved_surface_node_invariants(mesh);
}

BOOST_AUTO_TEST_CASE(p2_curved_cad_surface_native_input_regression)
{
  auto testCase = make_curved_surface_case();
  auto &elevatedMesh = static_cast<Mesh<MetricFieldFE>&>(
      *testCase->runner->msh_g);
  check_curved_surface_node_invariants(elevatedMesh);

  MetrisAPI nativeP2Data(*testCase->runner);
  MetrisParameters parameters;
  parameters.usrTarDeg = 2;
  parameters.adp_niter = 0;
  parameters.opt_niter = 0;
  parameters.iverb = 0;
  parameters.outmPrefix
      = curved_surface_output_directory("native_p2");
  testCase->runner
      = std::make_unique<MetrisRunner>(&nativeP2Data,nullptr,parameters);

  auto &nativeMesh = static_cast<Mesh<MetricFieldFE>&>(
      *testCase->runner->msh_g);
  // A complete runMetris() lifecycle for tdim=2, gdim=3 remains part of the
  // deferred embedded-surface scope; native construction and localization are
  // the contracts qualified here.
  BOOST_CHECK_EQUAL(testCase->runner->degElevate(),0);
  check_curved_surface_node_invariants(nativeMesh);
}
