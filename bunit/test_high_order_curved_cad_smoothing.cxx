// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#define BOOST_TEST_MODULE test_high_order_curved_cad_smoothing

#include "common_setup.hxx"

#include "API/MetrisAPI.hxx"
#include "Mesh/Mesh.hxx"
#include "MetrisRunner/MetrisParameters.hxx"
#include "MetrisRunner/MetrisRunner.hxx"
#include "low_geo/validity.hxx"
#include "smoothing/low_smooballdiff.hxx"
#include "smoothing/msh_smooball.hxx"

#include <cmath>
#include <filesystem>
#include <memory>

using namespace Metris;

namespace
{

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double circle_start = pi/6.0;
constexpr double circle_end = 5.0*pi/6.0;

struct CurvedPlanarCadCase
{
  ego context = nullptr;
  ego model = nullptr;
  std::unique_ptr<MetrisRunner> runner;

  ~CurvedPlanarCadCase()
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

std::unique_ptr<CurvedPlanarCadCase> make_curved_planar_cad_case()
{
  auto testCase = std::make_unique<CurvedPlanarCadCase>();
  require_egads_success(EG_open(&testCase->context),"EG_open");
  EG_setOutLevel(testCase->context,0);

  const double endpoint_x = std::cos(circle_start);
  const double endpoint_y = std::sin(circle_start);
  double coordinates[3][3] = {
      { endpoint_x,endpoint_y,0.0},
      {-endpoint_x,endpoint_y,0.0},
      {0.0,-1.0,0.0}};

  ego nodes[3];
  for(int inode = 0; inode < 3; inode++){
    require_egads_success(
        EG_makeTopology(
            testCase->context,nullptr,NODE,0,
            coordinates[inode],0,nullptr,nullptr,
            &nodes[inode]),
        "EG_makeTopology(NODE)");
  }

  ego edges[3];
  double circleData[10] = {
      0.0,0.0,0.0,
      1.0,0.0,0.0,
      0.0,1.0,0.0,
      1.0};
  ego circle;
  require_egads_success(
      EG_makeGeometry(
          testCase->context,CURVE,CIRCLE,nullptr,nullptr,
          circleData,&circle),
      "EG_makeGeometry(CIRCLE)");
  double curvedRange[2] = {circle_start,circle_end};
  ego curvedNodes[2] = {nodes[0],nodes[1]};
  require_egads_success(
      EG_makeTopology(
          testCase->context,circle,EDGE,TWONODE,curvedRange,
          2,curvedNodes,nullptr,&edges[0]),
      "EG_makeTopology(curved EDGE)");

  const int lineEndpoints[2][2] = {{1,2},{2,0}};
  for(int iline = 0; iline < 2; iline++){
    const int first = lineEndpoints[iline][0];
    const int second = lineEndpoints[iline][1];
    double lineData[6] = {
        coordinates[first][0],
        coordinates[first][1],
        coordinates[first][2],
        coordinates[second][0] - coordinates[first][0],
        coordinates[second][1] - coordinates[first][1],
        coordinates[second][2] - coordinates[first][2]};
    ego line;
    require_egads_success(
        EG_makeGeometry(
            testCase->context,CURVE,LINE,nullptr,nullptr,lineData,&line),
        "EG_makeGeometry(LINE)");
    double lineRange[2] = {0.0,1.0};
    ego lineNodes[2] = {nodes[first],nodes[second]};
    require_egads_success(
        EG_makeTopology(
            testCase->context,line,EDGE,TWONODE,lineRange,
            2,lineNodes,nullptr,&edges[iline + 1]),
        "EG_makeTopology(line EDGE)");
  }

  int edgeSenses[3] = {SFORWARD,SFORWARD,SFORWARD};
  ego loop;
  require_egads_success(
      EG_makeTopology(
          testCase->context,nullptr,LOOP,CLOSED,nullptr,
          3,edges,edgeSenses,&loop),
      "EG_makeTopology(LOOP)");
  ego face;
  require_egads_success(
      EG_makeFace(loop,SFORWARD,nullptr,&face),"EG_makeFace");
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

  MetrisAPI input(
      2,1,3,6,0,3,false,3,1,0,
      FEBasis::Lagrange,FEBasis::Undefined,MetSpace::Exp);
  for(int ipoint = 0; ipoint < 3; ipoint++){
    input.setCoord(ipoint,coordinates[ipoint]);
    input.setCorner(ipoint,ipoint);
  }

  for(int iedge = 0; iedge < 3; iedge++){
    const int curvedEdgePoints[2] = {0,1};
    const int *points
        = iedge == 0 ? curvedEdgePoints : lineEndpoints[iedge - 1];
    const int reference = EG_indexBodyTopo(body,edges[iedge]) - 1;
    BOOST_REQUIRE_GE(reference,0);
    input.setElement(1,iedge,points,reference);
  }
  const int facePoints[3] = {0,1,2};
  input.setElement(2,0,facePoints,0);

  const int geometricEdgePoints[6][2] = {
      {0,0},{1,0},
      {1,1},{2,1},
      {2,2},{0,2}};
  const double geometricEdgeParameters[6][2] = {
      {circle_start,0.0},{circle_end,0.0},
      {0.0,0.0},{1.0,0.0},
      {0.0,0.0},{1.0,0.0}};
  input.setVerticesOnGeometricEdges(
      0,6,&geometricEdgePoints[0][0],
      &geometricEdgeParameters[0][0]);
  input.setCADModel(testCase->context,testCase->model);

  MetrisParameters parameters;
  parameters.setAnalyticalMetric(1);
  parameters.usrTarDeg = 2;
  parameters.adp_niter = 0;
  parameters.opt_niter = 0;
  parameters.iverb = 0;
  parameters.outmPrefix
      = (std::filesystem::temp_directory_path()
         / "metris_high_order_curved_cad_smoothing").string() + "/";
  testCase->runner = std::make_unique<MetrisRunner>(
      &input,nullptr,parameters);
  BOOST_REQUIRE_EQUAL(testCase->runner->degElevate(),1);
  return testCase;
}

void check_curved_boundary_smoothing(QuaFun objective)
{
  auto testCase = make_curved_planar_cad_case();
  auto &mesh = static_cast<Mesh<MetricFieldAnalytical>&>(
      *testCase->runner->msh_g);
  mesh.met.setSpace(MetSpace::Exp);
  mesh.setBasis(FEBasis::Lagrange);

  constexpr int curvedMeshEdge = 0;
  const int controlPoint = mesh.edg2poi(curvedMeshEdge,2);
  BOOST_REQUIRE_GE(controlPoint,0);
  const int boundaryRecord
      = mesh.poi2ebp(controlPoint,1,curvedMeshEdge,-1);
  BOOST_REQUIRE_GE(boundaryRecord,0);
  BOOST_CHECK_EQUAL(mesh.bpo2ibi(boundaryRecord,1),1);

  const int reference = mesh.edg2ref[curvedMeshEdge];
  const ego cadEdge = mesh.CAD.cad2edg[reference];
  const double elevatedParameter = mesh.bpo2rbi(boundaryRecord,0);
  BOOST_CHECK_SMALL(
      elevatedParameter - 0.5*(circle_start + circle_end),1.e-14);
  double evaluation[18];
  require_egads_success(
      EG_evaluate(cadEdge,&elevatedParameter,evaluation),
      "EG_evaluate(elevated control point)");
  BOOST_CHECK_SMALL(mesh.coord(controlPoint,0) - evaluation[0],1.e-13);
  BOOST_CHECK_SMALL(mesh.coord(controlPoint,1) - evaluation[1],1.e-13);

  const int firstEndpoint = mesh.edg2poi(curvedMeshEdge,0);
  const int secondEndpoint = mesh.edg2poi(curvedMeshEdge,1);
  const double chordMidpoint[2] = {
      0.5*(mesh.coord(firstEndpoint,0) + mesh.coord(secondEndpoint,0)),
      0.5*(mesh.coord(firstEndpoint,1) + mesh.coord(secondEndpoint,1))};
  BOOST_CHECK_GT(
      std::hypot(
          mesh.coord(controlPoint,0) - chordMidpoint[0],
          mesh.coord(controlPoint,1) - chordMidpoint[1]),
      0.25);

  int seedFace = -1;
  int localEdge = -1;
  for(int candidateEdge = 0; candidateEdge < 3; candidateEdge++){
    if(mesh.fac2poi(0,3 + candidateEdge) != controlPoint) continue;
    seedFace = 0;
    localEdge = candidateEdge;
    break;
  }
  BOOST_REQUIRE_GE(localEdge,0);
  intAr1 region(2);
  buildEdgeControlPointSmoothingRegion(
      mesh,2,seedFace,localEdge,region);
  BOOST_REQUIRE_EQUAL(region.get_n(),1);

  const double perturbedParameter = elevatedParameter - 0.25;
  require_egads_success(
      EG_evaluate(cadEdge,&perturbedParameter,evaluation),
      "EG_evaluate(perturbed control point)");
  mesh.bpo2rbi(boundaryRecord,0) = perturbedParameter;
  mesh.coord(controlPoint,0) = evaluation[0];
  mesh.coord(controlPoint,1) = evaluation[1];
  BOOST_REQUIRE((
      classify_element_validity<2,2>(mesh,seedFace)
          .accepted_conservatively()));

  double objectiveBefore = 0.0;
  double maximumBefore = 0.0;
  double objectiveAfter = 0.0;
  double maximumAfter = 0.0;
  const int status
      = smooballdiff_boundary<MetricFieldAnalytical,2,2>(
          mesh,controlPoint,1,region,
          &objectiveBefore,&maximumBefore,
          &objectiveAfter,&maximumAfter,objective);
  BOOST_REQUIRE_EQUAL(status,0);
  BOOST_CHECK_LE(objectiveAfter,objectiveBefore);

  const double optimizedParameter = mesh.bpo2rbi(boundaryRecord,0);
  BOOST_CHECK_GT(
      std::abs(optimizedParameter - perturbedParameter),1.e-10);
  BOOST_CHECK_GT(optimizedParameter,circle_start);
  BOOST_CHECK_LT(optimizedParameter,circle_end);
  require_egads_success(
      EG_evaluate(cadEdge,&optimizedParameter,evaluation),
      "EG_evaluate(optimized control point)");
  BOOST_CHECK_SMALL(mesh.coord(controlPoint,0) - evaluation[0],1.e-12);
  BOOST_CHECK_SMALL(mesh.coord(controlPoint,1) - evaluation[1],1.e-12);
  BOOST_CHECK((
      classify_element_validity<2,2>(mesh,seedFace)
          .accepted_conservatively()));
}

} // namespace

BOOST_AUTO_TEST_CASE(p2_curved_cad_edge_sizeshape_smoothing)
{
  check_curved_boundary_smoothing(QuaFun::SizeShape);
}

BOOST_AUTO_TEST_CASE(p2_curved_cad_edge_stepdistance_smoothing)
{
  check_curved_boundary_smoothing(QuaFun::StepDistance);
}
