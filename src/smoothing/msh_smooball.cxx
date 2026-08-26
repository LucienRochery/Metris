//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

/*
Routine for "direct" smoothing as P1. From each (facet, metric) pair, generate remaining vertex to be unit. Then average over ball.
Simplest possible approach.
*/

#include "../Mesh/Mesh.hxx"
#include "../MetrisRunner/MetrisParameters.hxx"

#include "../smoothing/msh_smooball.hxx"
#include "../smoothing/low_smooballdiff.hxx"
#include "../smoothing/smoothing_progress.hxx"

#include "../aux_topo.hxx"
#include "../utils/aux_timer.hxx"
#include "../low_topo.hxx"
#include "../utils/mprintf.hxx"
#include "../quality/low_metqua.hxx"
#include "../low_geo/measure.hxx"
#include "../low_geo/validity.hxx"
#include "../io_libmeshb.hxx"

#include "lplib3/lplib3.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#undef USE_LPLIB_SMOOTHINTERIOR


namespace Metris{

namespace {

struct BoundaryParameterSnapshot
{
  std::vector<int> records;
  std::vector<std::array<double,2>> parameters;
};

BoundaryParameterSnapshot captureBoundaryParameters(
    const MeshBase &msh,
    int point)
{
  BoundaryParameterSnapshot snapshot;
  int record = msh.poi2bpo[point];
  int iteration = 0;
  while(record >= 0 && msh.bpo2ibi(record,0) == point){
    METRIS_ENFORCE_MSG(
        iteration++ < METRIS_MAX_WHILE,
        "Boundary-record cycle while capturing point {}",point);
    snapshot.records.push_back(record);
    snapshot.parameters.push_back({
        msh.bpo2rbi(record,0),msh.bpo2rbi(record,1)});
    record = msh.bpo2ibi(record,3);
  }
  return snapshot;
}

void restoreBoundaryParameters(
    MeshBase &msh,
    const BoundaryParameterSnapshot &snapshot)
{
  METRIS_ASSERT(snapshot.records.size() == snapshot.parameters.size());
  for(std::size_t index = 0; index < snapshot.records.size(); index++){
    const int record = snapshot.records[index];
    msh.bpo2rbi(record,0) = snapshot.parameters[index][0];
    msh.bpo2rbi(record,1) = snapshot.parameters[index][1];
  }
}

// A curve-owned point can also carry one face record for each incident CAD
// face. Keep those secondary (u,v) parameters synchronized with the primary
// curve parameter before the move is counted as accepted.
void synchronizeCurvePointFaceParameters(MeshBase &msh, int point)
{
  const int primaryRecord = msh.poi2bpo[point];
  METRIS_ENFORCE_MSG(primaryRecord >= 0,
                     "Missing boundary record for point {}",point);
  METRIS_ENFORCE_MSG(msh.bpo2ibi(primaryRecord,1) == 1,
                     "Point {} is not primarily curve-owned",point);

  const int edge = msh.bpo2ibi(primaryRecord,2);
  METRIS_ENFORCE_MSG(edge >= 0 && edge < msh.nedge,
                     "Invalid owner edge {} for point {}",edge,point);
  const int edgeReference = msh.edg2ref[edge];
  METRIS_ENFORCE_MSG(
      edgeReference >= 0 && edgeReference < msh.CAD.ncaded,
      "Invalid CAD edge reference {} for mesh edge {}",
      edgeReference,edge);
  const ego cadEdge = msh.CAD.cad2edg[edgeReference];
  const double parameter = msh.bpo2rbi(primaryRecord,0);

  int record = primaryRecord;
  int iteration = 0;
  while(record >= 0 && msh.bpo2ibi(record,0) == point){
    METRIS_ENFORCE_MSG(
        iteration++ < METRIS_MAX_WHILE,
        "Boundary-record cycle while synchronizing point {}",point);
    if(msh.bpo2ibi(record,1) == 2){
      const int face = msh.bpo2ibi(record,2);
      METRIS_ENFORCE_MSG(face >= 0 && face < msh.nface,
                         "Invalid incident face {} for point {}",face,point);
      const int faceReference = msh.fac2ref[face];
      METRIS_ENFORCE_MSG(
          faceReference >= 0 && faceReference < msh.CAD.ncadfa,
          "Invalid CAD face reference {} for mesh face {}",
          faceReference,face);
      const ego cadFace = msh.CAD.cad2fac[faceReference];
      double parameters[2];
      const int status
          = EG_getEdgeUV(cadFace,cadEdge,0,parameter,parameters);
      METRIS_ENFORCE_MSG(
          status == EGADS_SUCCESS,"EG_getEdgeUV failed: {}",status);
      msh.bpo2rbi(record,0) = parameters[0];
      msh.bpo2rbi(record,1) = parameters[1];
    }
    record = msh.bpo2ibi(record,3);
  }
}

template<class MFT, int idim, int ideg>
void enforceAcceptedBoundarySmoothingInvariants(
    const Mesh<MFT> &msh,
    int point,
    const intAr1 &region)
{
  constexpr double coordinateTolerance = 1.e-10;
  const double squaredTolerance
      = coordinateTolerance*coordinateTolerance;

  int record = msh.poi2bpo[point];
  METRIS_ENFORCE_MSG(record >= 0,
                     "Accepted boundary point {} has no CAD record",point);
  int iteration = 0;
  while(record >= 0 && msh.bpo2ibi(record,0) == point){
    METRIS_ENFORCE_MSG(
        iteration++ < METRIS_MAX_WHILE,
        "Boundary-record cycle while verifying point {}",point);
    const int cadDimension = msh.bpo2ibi(record,1);
    if(cadDimension == 1 || cadDimension == 2){
      const int entity = msh.bpo2ibi(record,2);
      ego cadEntity = nullptr;
      if(cadDimension == 1){
        METRIS_ENFORCE_MSG(entity >= 0 && entity < msh.nedge,
                           "Invalid owner edge {} for point {}",entity,point);
        const int reference = msh.edg2ref[entity];
        METRIS_ENFORCE_MSG(
            reference >= 0 && reference < msh.CAD.ncaded,
            "Invalid CAD edge reference {} for mesh edge {}",
            reference,entity);
        cadEntity = msh.CAD.cad2edg[reference];
        METRIS_ENFORCE_MSG(
            std::isfinite(msh.bpo2rbi(record,0)),
            "Non-finite CAD curve parameter for point {}",point);
      }else{
        METRIS_ENFORCE_MSG(entity >= 0 && entity < msh.nface,
                           "Invalid owner face {} for point {}",entity,point);
        const int reference = msh.fac2ref[entity];
        METRIS_ENFORCE_MSG(
            reference >= 0 && reference < msh.CAD.ncadfa,
            "Invalid CAD face reference {} for mesh face {}",
            reference,entity);
        cadEntity = msh.CAD.cad2fac[reference];
        METRIS_ENFORCE_MSG(
            std::isfinite(msh.bpo2rbi(record,0))
            && std::isfinite(msh.bpo2rbi(record,1)),
            "Non-finite CAD surface parameters for point {}",point);
      }

      double evaluation[18];
      const int status
          = EG_evaluate(cadEntity,&msh.bpo2rbi(record,0),evaluation);
      METRIS_ENFORCE_MSG(
          status == EGADS_SUCCESS,
          "EG_evaluate failed for accepted point {} record {}: {}",
          point,record,status);
      double squaredError = 0.;
      for(int component = 0; component < idim; component++){
        const double difference
            = msh.coord(point,component) - evaluation[component];
        squaredError += difference*difference;
      }
      METRIS_ENFORCE_MSG(
          squaredError <= squaredTolerance,
          "Accepted boundary point {} is inconsistent with CAD record {}: "
          "distance {}",
          point,record,std::sqrt(squaredError));
    }
    record = msh.bpo2ibi(record,3);
  }

  for(const int element : region){
    if constexpr(ideg == 1){
      METRIS_ENFORCE_MSG(
          (isvalideltP1<idim,idim>(msh,element)),
          "Accepted boundary move invalidated P1 element {}",element);
    }else{
      const ElementValidityResult validity
          = classify_element_validity<idim,ideg>(msh,element);
      METRIS_ENFORCE_MSG(
          validity.accepted_conservatively(),
          "Accepted boundary move left element {} uncertified or invalid",
          element);
    }
  }
}

using Point2 = std::array<double,2>;

struct QuadraticBezierSegment
{
  std::array<Point2,3> control;
  double parameter0 = 0.;
  double parameter1 = 1.;
};

struct Box2
{
  Point2 minimum;
  Point2 maximum;
};

struct SegmentProximity
{
  double squaredDistance = std::numeric_limits<double>::infinity();
  double firstParameter = 0.;
  double secondParameter = 0.;
};

double dot2(const Point2 &first, const Point2 &second)
{
  return first[0]*second[0] + first[1]*second[1];
}

double cross2(const Point2 &first, const Point2 &second)
{
  return first[0]*second[1] - first[1]*second[0];
}

Point2 difference2(const Point2 &first, const Point2 &second)
{
  return {first[0] - second[0],first[1] - second[1]};
}

Point2 midpoint2(const Point2 &first, const Point2 &second)
{
  return {0.5*(first[0] + second[0]),
          0.5*(first[1] + second[1])};
}

Box2 quadraticBezierBox(const QuadraticBezierSegment &segment)
{
  Box2 box{segment.control[0],segment.control[0]};
  for(int control = 1; control < 3; control++){
    for(int component = 0; component < 2; component++){
      box.minimum[component]
          = std::min(box.minimum[component],segment.control[control][component]);
      box.maximum[component]
          = std::max(box.maximum[component],segment.control[control][component]);
    }
  }
  return box;
}

bool boxesAreSeparated(const Box2 &first,
                       const Box2 &second,
                       double tolerance)
{
  for(int component = 0; component < 2; component++){
    if(first.maximum[component] < second.minimum[component] - tolerance
       || second.maximum[component]
              < first.minimum[component] - tolerance){
      return true;
    }
  }
  return false;
}

double pointSegmentSquaredDistance(const Point2 &point,
                                   const Point2 &first,
                                   const Point2 &second,
                                   double &parameter)
{
  const Point2 direction = difference2(second,first);
  const double squaredLength = dot2(direction,direction);
  if(squaredLength <= std::numeric_limits<double>::min()){
    parameter = 0.;
    const Point2 offset = difference2(point,first);
    return dot2(offset,offset);
  }
  parameter = dot2(difference2(point,first),direction)/squaredLength;
  parameter = std::max(0.,std::min(1.,parameter));
  const Point2 projection = {
      first[0] + parameter*direction[0],
      first[1] + parameter*direction[1]};
  const Point2 offset = difference2(point,projection);
  return dot2(offset,offset);
}

SegmentProximity segmentProximity(const Point2 &first0,
                                  const Point2 &first1,
                                  const Point2 &second0,
                                  const Point2 &second1)
{
  SegmentProximity proximity;
  const Point2 firstDirection = difference2(first1,first0);
  const Point2 secondDirection = difference2(second1,second0);
  const Point2 origins = difference2(second0,first0);
  const double denominator = cross2(firstDirection,secondDirection);
  const double denominatorScale
      = std::sqrt(dot2(firstDirection,firstDirection)
                  *dot2(secondDirection,secondDirection));
  const double parallelTolerance
      = 32.*std::numeric_limits<double>::epsilon()
       *std::max(1.,denominatorScale);
  if(std::abs(denominator) > parallelTolerance){
    const double firstParameter
        = cross2(origins,secondDirection)/denominator;
    const double secondParameter
        = cross2(origins,firstDirection)/denominator;
    if(firstParameter >= 0. && firstParameter <= 1.
       && secondParameter >= 0. && secondParameter <= 1.){
      proximity.squaredDistance = 0.;
      proximity.firstParameter = firstParameter;
      proximity.secondParameter = secondParameter;
      return proximity;
    }
  }

  double parameter = 0.;
  double distance = pointSegmentSquaredDistance(
      first0,second0,second1,parameter);
  proximity = {distance,0.,parameter};
  distance = pointSegmentSquaredDistance(
      first1,second0,second1,parameter);
  if(distance < proximity.squaredDistance){
    proximity = {distance,1.,parameter};
  }
  distance = pointSegmentSquaredDistance(
      second0,first0,first1,parameter);
  if(distance < proximity.squaredDistance){
    proximity = {distance,parameter,0.};
  }
  distance = pointSegmentSquaredDistance(
      second1,first0,first1,parameter);
  if(distance < proximity.squaredDistance){
    proximity = {distance,parameter,1.};
  }
  return proximity;
}

double quadraticBezierFlatness(const QuadraticBezierSegment &segment)
{
  double parameter = 0.;
  return std::sqrt(pointSegmentSquaredDistance(
      segment.control[1],segment.control[0],segment.control[2],parameter));
}

std::array<QuadraticBezierSegment,2> subdivideQuadraticBezier(
    const QuadraticBezierSegment &segment)
{
  const Point2 firstMidpoint
      = midpoint2(segment.control[0],segment.control[1]);
  const Point2 secondMidpoint
      = midpoint2(segment.control[1],segment.control[2]);
  const Point2 curveMidpoint = midpoint2(firstMidpoint,secondMidpoint);
  const double parameterMidpoint
      = 0.5*(segment.parameter0 + segment.parameter1);
  return {{
      {{{segment.control[0],firstMidpoint,curveMidpoint}},
       segment.parameter0,parameterMidpoint},
      {{{curveMidpoint,secondMidpoint,segment.control[2]}},
       parameterMidpoint,segment.parameter1}}};
}

bool isAllowedSharedEndpointContact(
    const std::array<int,2> &firstEndpoints,
    const std::array<int,2> &secondEndpoints,
    double firstParameter,
    double secondParameter)
{
  constexpr double parameterTolerance = 1.e-8;
  for(int firstEndpoint = 0; firstEndpoint < 2; firstEndpoint++){
    for(int secondEndpoint = 0; secondEndpoint < 2; secondEndpoint++){
      if(firstEndpoints[firstEndpoint] != secondEndpoints[secondEndpoint]){
        continue;
      }
      if(std::abs(firstParameter - firstEndpoint) <= parameterTolerance
         && std::abs(secondParameter - secondEndpoint)
                <= parameterTolerance){
        return true;
      }
    }
  }
  return false;
}

bool quadraticBezierSegmentsIntersect(
    const QuadraticBezierSegment &first,
    const QuadraticBezierSegment &second,
    const std::array<int,2> &firstEndpoints,
    const std::array<int,2> &secondEndpoints,
    double tolerance,
    int depth)
{
  if(boxesAreSeparated(
         quadraticBezierBox(first),quadraticBezierBox(second),tolerance)){
    return false;
  }

  const double firstFlatness = quadraticBezierFlatness(first);
  const double secondFlatness = quadraticBezierFlatness(second);
  constexpr int maximumDepth = 30;
  if(depth >= maximumDepth
     || (firstFlatness <= tolerance && secondFlatness <= tolerance)){
    const SegmentProximity proximity = segmentProximity(
        first.control[0],first.control[2],
        second.control[0],second.control[2]);
    if(proximity.squaredDistance > 16.*tolerance*tolerance) return false;
    const double firstParameter
        = first.parameter0
        + proximity.firstParameter*(first.parameter1 - first.parameter0);
    const double secondParameter
        = second.parameter0
        + proximity.secondParameter*(second.parameter1 - second.parameter0);
    return !isAllowedSharedEndpointContact(
        firstEndpoints,secondEndpoints,firstParameter,secondParameter);
  }

  if(firstFlatness >= secondFlatness){
    const auto halves = subdivideQuadraticBezier(first);
    return quadraticBezierSegmentsIntersect(
               halves[0],second,firstEndpoints,secondEndpoints,
               tolerance,depth + 1)
        || quadraticBezierSegmentsIntersect(
               halves[1],second,firstEndpoints,secondEndpoints,
               tolerance,depth + 1);
  }
  const auto halves = subdivideQuadraticBezier(second);
  return quadraticBezierSegmentsIntersect(
             first,halves[0],firstEndpoints,secondEndpoints,
             tolerance,depth + 1)
      || quadraticBezierSegmentsIntersect(
             first,halves[1],firstEndpoints,secondEndpoints,
             tolerance,depth + 1);
}

QuadraticBezierSegment getPlanarP2BoundaryBezier(
    const MeshBase &msh,
    int edge)
{
  const int firstPoint = msh.edg2poi(edge,0);
  const int secondPoint = msh.edg2poi(edge,1);
  const int midpoint = msh.edg2poi(edge,2);
  QuadraticBezierSegment segment;
  for(int component = 0; component < 2; component++){
    segment.control[0][component] = msh.coord(firstPoint,component);
    segment.control[2][component] = msh.coord(secondPoint,component);
    segment.control[1][component]
        = 2.*msh.coord(midpoint,component)
        - 0.5*(msh.coord(firstPoint,component)
               + msh.coord(secondPoint,component));
  }
  return segment;
}

bool edgeContainsP2Point(const MeshBase &msh, int edge, int point)
{
  for(int localNode = 0; localNode < 3; localNode++){
    if(msh.edg2poi(edge,localNode) == point) return true;
  }
  return false;
}

} // namespace

bool planarP2BoundaryIsIntersectionFreeAroundPoint(
    const MeshBase &msh,
    int movedPoint)
{
  METRIS_ENFORCE_MSG(msh.idim == 2 && msh.get_tdim() == 2,
                     "Planar boundary intersection check requires a 2D mesh");
  METRIS_ENFORCE_MSG(msh.curdeg == 2,
                     "Planar P2 boundary check received degree {}",msh.curdeg);
  METRIS_ENFORCE_MSG(
      msh.getBasis() == FEBasis::Lagrange,
      "Planar P2 boundary check requires Lagrange geometry");
  METRIS_ENFORCE_MSG(movedPoint >= 0 && movedPoint < msh.npoin,
                     "Invalid moved boundary point {}",movedPoint);

  for(int firstEdge = 0; firstEdge < msh.nedge; firstEdge++){
    if(isdeadent(firstEdge,msh.edg2poi)
       || !edgeContainsP2Point(msh,firstEdge,movedPoint)) continue;
    const QuadraticBezierSegment first
        = getPlanarP2BoundaryBezier(msh,firstEdge);
    const std::array<int,2> firstEndpoints
        = {msh.edg2poi(firstEdge,0),msh.edg2poi(firstEdge,1)};

    for(int secondEdge = 0; secondEdge < msh.nedge; secondEdge++){
      if(secondEdge == firstEdge || isdeadent(secondEdge,msh.edg2poi)){
        continue;
      }
      if(edgeContainsP2Point(msh,secondEdge,movedPoint)
         && secondEdge < firstEdge) continue;

      const QuadraticBezierSegment second
          = getPlanarP2BoundaryBezier(msh,secondEdge);
      const std::array<int,2> secondEndpoints
          = {msh.edg2poi(secondEdge,0),msh.edg2poi(secondEdge,1)};
      Point2 coordinateMinimum = first.control[0];
      Point2 coordinateMaximum = first.control[0];
      double coordinateMagnitude = 1.;
      for(const Point2 &control : first.control){
        for(int component = 0; component < 2; component++){
          coordinateMinimum[component] = std::min(
              coordinateMinimum[component],control[component]);
          coordinateMaximum[component] = std::max(
              coordinateMaximum[component],control[component]);
        }
        coordinateMagnitude = std::max(
            coordinateMagnitude,
            std::max(std::abs(control[0]),std::abs(control[1])));
      }
      for(const Point2 &control : second.control){
        for(int component = 0; component < 2; component++){
          coordinateMinimum[component] = std::min(
              coordinateMinimum[component],control[component]);
          coordinateMaximum[component] = std::max(
              coordinateMaximum[component],control[component]);
        }
        coordinateMagnitude = std::max(
            coordinateMagnitude,
            std::max(std::abs(control[0]),std::abs(control[1])));
      }
      const double geometricScale = std::max(
          coordinateMaximum[0] - coordinateMinimum[0],
          coordinateMaximum[1] - coordinateMinimum[1]);
      const double tolerance
          = 1.e-10*std::max(1.,geometricScale)
          + 64.*std::numeric_limits<double>::epsilon()
               *coordinateMagnitude;
      if(quadraticBezierSegmentsIntersect(
             first,second,firstEndpoints,secondEndpoints,tolerance,0)){
        return false;
      }
    }
  }
  return true;
}

void buildEdgeControlPointSmoothingRegion(
    const MeshBase &msh,
    int tdim,
    int seed_entity,
    int local_edge,
    intAr1 &region){
  METRIS_ENFORCE_MSG(tdim == 2 || tdim == 3,
                     "Unsupported smoothing-region dimension {}",tdim);
  METRIS_ENFORCE_MSG(seed_entity >= 0 && seed_entity < msh.nentt(tdim),
                     "Invalid smoothing-region seed {}",seed_entity);

  region.allocate(tdim == 2 ? 2 : 10);
  region.set_n(0);
  region.stack(seed_entity);

  if(tdim == 2){
    METRIS_ENFORCE_MSG(local_edge >= 0 && local_edge < 3,
                       "Invalid triangle local edge {}",local_edge);
    const int neighbor = msh.fac2fac(seed_entity,local_edge);
    if(neighbor >= 0) region.stack(neighbor);
    return;
  }

  METRIS_ENFORCE_MSG(local_edge >= 0 && local_edge < 6,
                     "Invalid tetrahedron local edge {}",local_edge);
  const int ip1 = msh.tet2poi(seed_entity,lnoed3[local_edge][0]);
  const int ip2 = msh.tet2poi(seed_entity,lnoed3[local_edge][1]);
  intAr1 boundary_faces;
  int open_shell;
  shell3(msh,ip1,ip2,seed_entity,region,boundary_faces,&open_shell);
}

void reactivateSmoothingRegionGeometry(
    MeshBase &msh,
    int tdim,
    int geometry_degree,
    const intAr1 &region,
    int ithread)
{
  METRIS_ENFORCE_MSG(tdim == 2 || tdim == 3,
                     "Unsupported smoothing dimension {}",tdim);
  METRIS_ENFORCE_MSG(geometry_degree >= 1,
                     "Invalid smoothing geometry degree {}",geometry_degree);
  METRIS_ENFORCE_MSG(ithread >= 0 && ithread < METRIS_MAXTAGS,
                     "Invalid smoothing tag thread {}",ithread);

  const int nodes_per_element = getnnode(tdim,geometry_degree);
  const int reactivated_tag = msh.tag[ithread] - 1;
  const intAr2 &element_points = msh.ent2poi(tdim);
  for(const int element : region){
    if(element < 0 || element >= msh.nentt(tdim)
       || isdeadent(element,element_points)) continue;
    for(int local_node = 0; local_node < nodes_per_element; local_node++){
      const int point = element_points(element,local_node);
      METRIS_ASSERT(point >= 0 && point < msh.npoin);
      if(msh.isdeadpoint(point)) continue;
      msh.poi2tag(ithread,point) = reactivated_tag;
    }
  }
}

template<class MFT>
double smoothing_region_objective(const Mesh<MFT>& msh,
                                  QuaFun iquaf,
                                  double numerator,
                                  double targetWeight){
  if(iquaf == QuaFun::StepDistance
     && msh.param->step_distance_cavity_target_average){
    METRIS_ENFORCE_MSG(targetWeight > 0.,
                       "Nonpositive StepDistance smoothing target weight");
    return numerator/targetWeight;
  }
  return numerator;
}

template<class MFT>
double smoothInterior_Ball(Mesh<MFT> &msh, QuaFun iquaf, int ithrd1, int ithrd2){

  int tdimn = msh.get_tdim();

  METRIS_ASSERT_MSG(tdimn > 1, "TODO: edge smooth interior ball");

  // Geo and topo dimn must match otherwise surface specific
  METRIS_ASSERT(tdimn == msh.idim);
  double noper;
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    if(tdimn == 2){
      noper = smoothInterior_Ball0<MFT,2,ideg>(msh,iquaf,ithrd1,ithrd2);
    }else{
      noper = smoothInterior_Ball0<MFT,3,ideg>(msh,iquaf,ithrd1,ithrd2);
    }
  }}CT_FOR1(ideg);

  return noper;
}

template double smoothInterior_Ball<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh,
                                             QuaFun iquaf, int ithrd1, int ithrd2);
template double smoothInterior_Ball<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh,
                                             QuaFun iquaf, int ithrd1, int ithrd2);


// =============================================================================================== //
// =============================================================================================== //

// idim: gdim = tdim
template<class MFT, int idim, int ideg>
double smoothInterior_Ball0(Mesh<MFT> &msh, QuaFun iquaf,
                            int ithrd1, int ithrd2){
  GETVDEPTH(msh.param);

  constexpr int tdim = idim;
  //constexpr int gdim = idim;
  if(ideg > idim + 1){
    PRINTF("## SMOOTHING DISABLED FOR DEGREE {} and dim {} \n",ideg,idim);
    return -1.0;
  }


  int nentt = msh.nentt(tdim);
  const intAr2 &ent2poi = msh.ent2poi(tdim);


  #ifdef USE_LPLIB_SMOOTHINTERIOR
    // LPlib init
    int nproc = msh.param->nproc;
    int nthread = GetNumberOfCores();
    if(nthread <= 0){
      CPRINTF1("## WARNING: LPlib function GetNumberOfCores() returned "
               "negative threads. Set to default {}.\n",METRIS_MAXTAGS);
      nthread = METRIS_MAXTAGS;
    }else{
      CPRINTF2("-- LPlib found ncore = {} \n",nthread);
      if(nthread > METRIS_MAXTAGS){
        CPRINTF1("## WARNING: must verify nthread <= METRIS_MAXTAGS = {}."
                " Increase in metris_constants.hxx.\n", METRIS_MAXTAGS);
        nthread = METRIS_MAXTAGS;
      }
    }
    if(nproc > 0) nthread = MIN(nthread, nproc);
    int64_t LibIdx = InitParallel(nthread);
    int LP_elt = NewType(LibIdx, nentt);
    int LP_poi = NewType(LibIdx, msh.npoin);
    float LP_stat[2];
    BeginDependency(LibIdx, LP_elt, LP_poi);
    for(int ientt = 0; ientt < nentt; ientt++){
      for(int ii = 0; ii < tdim + 1; ii++)
        AddDependency(LibIdx, ientt+1, ent2poi(ientt,ii)+1);
    }
    EndDependency(LibIdx, LP_stat);
    // END LPlib init

    int itag_shared = ithrd1;
    // ithrd2 can be used freely as it is for elements, whose collisions are
    // avoided by LPlib
  #endif

  //msh.met.setSpace(MetSpace::Log);



  // Eventually move all constants to MetrisParameters
  // L2 conformity error from 0 to 1
  const double difto = 1.0;
  const int miter = msh.param->opt_smoo_niter;
  //const double maxwt = 20.0;
  //const double qrthr = 2.0;
  const double tolavg = msh.param->opt_smoo_tol;
  const double tolmax = msh.param->opt_smoo_tol;

  dblAr1 work;
  if(msh.param->iflag2 != 0){
    MPRINTF("\n\n##WARNING EXPERIMENTAL SMOOTHING FUNCTION 2\n");
  }

  // 1 -> no maximum quality increase allowed
  //const double maxinc_worst = 1.00;

  constexpr int nnmet = (idim*(idim+1))/2;

  METRIS_ASSERT(msh.param->opt_power == 1 || msh.param->opt_power == -1);
  // Otherwise not only edge nodes
  METRIS_ENFORCE(ideg <= tdim + 1);


  #ifndef USE_LPLIB_SMOOTHINTERIOR
  const int mball = 100;
  intAr1 lball(mball);
  #endif

  msh.tag[ithrd1]++;

  auto quafun = get_quafun<MFT,idim,idim>(iquaf);

  const bool useGlobalStepDistance =
      iquaf == QuaFun::StepDistance
      && msh.param->step_distance_cavity_target_average;
  StepDistanceObjectiveState globalStepDistance;
  if(useGlobalStepDistance){
    globalStepDistance = step_distance_global_objective_state<MFT,idim,idim>(
        msh,AsDeg::Pk,AsDeg::Pk);
  }

  auto accumulateStepDistanceRegion = [&](const intAr1& region,
                                           double& numerator,
                                           int& count,
                                           double& targetWeight){
    numerator = 0.;
    count = 0;
    targetWeight = 0.;
    for(const int ient2 : region){
      if(isdeadent(ient2,ent2poi)) continue;
      const double quality = quafun(
          msh,AsDeg::Pk,AsDeg::Pk,ient2,difto);
      if(iquaf == QuaFun::StepDistance
         && msh.param->step_distance_cavity_target_average){
        const double weight =
            step_distance_element_target_weight<MFT,idim,idim>(
                msh,AsDeg::Pk,ient2);
        numerator += weight*quality;
        targetWeight += weight;
      }else{
        numerator += quality;
      }
      count++;
    }
  };

  double noper = 0;
  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh.param);

    #if 0
    for(int ii = 0; ii < msh.npoin; ii++){
      rpoqe[ii] = 0.0;
      ipone[ii] = 0;
    }
    #endif

    double qmin = 1.0e30,qmax = -1.0e30, qnrm = 0.0;
    int imax = -1;
    int navg = 0;
    for(int ientt = 0; ientt < nentt; ientt++){

      if(isdeadent(ientt,ent2poi)) continue;

      navg++;

      double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,ientt,difto);

      qnrm += quael;
      qmin = MIN(qmin,quael);
      if(qmax < quael){
        imax = ientt;
        qmax = quael;
      }
    }

    qnrm /= navg;
    double t0 = get_cpu_time();
    CPRINTF1(" - smoo iter {:3} init {:10.6e} < q < {:10.6e} (at {}), avg = {:10.6e} "
                   "(p = {})\n",niter,qmin,qmax,imax,qnrm,msh.param->opt_pnorm);
    //if(iverb >= 2 && qmax >= 1e10){
    //  printf("## HIGH QMAX mshdeg = {} \n",msh.curdeg);
    //  std::string fname = "qmax"+std::to_string(imax);
    //  writeMesh(fname,msh);
    //  //wait();
    //}

    int nsucc = 0;
    int nmov  = 0;

    #if 0
    #ifdef USE_LPLIB_SMOOTHINTERIOR
    void (*inerloop_LPlib)(int,int,int,Mesh<MFT>*,int,int,double*,
          intAr1*,intAr1*,double,double)
    = [] (int ipoi0, int ipoi1, int ithread, Mesh<MFT> *msh,
          int itag_shared, int itag2, double *qmax,
          intAr1 *nsuccthr, intAr1 *nmovthr,
          double tolavg, double tolmax){

      const int mball = 100;
      constexpr int nnmet = (idim*(idim+1))/2;
      const int iverb = msh->param->iverb;
      intAr2 &ent2poi = tdim == 2 ? msh->fac2poi : msh->tet2poi;

      intAr1 lball(mball);
      dblAr1 qball(mball);

      for(int ipoin = ipoi0 - 1; ipoin < ipoi1; ipoin++){

        if(msh->poi2tag(itag_shared,ipoin) >= msh->tag[itag_shared]) continue;

        int ib = msh->poi2bpo[ipoin];
        if(ib >= 0) continue;

        int ientt = getpoient(*msh, ipoin, tdim);

        //double qpoin = rpoqe[ipoin] / ipone[ipoin];
        //if(qpoin > qnrm / qrthr){
          if(iverb >= 3){
            //printf("   - smoo pt {} seed elt {} quapt = {:10.6e}"
            //  " qthrs = {:10.6e} qnrm = {:10.6e}\n",
            //  ipoin,ientt,qpoin,qrthr * qnrm,qnrm);
            printf("   - smoo pt {} seed elt {} \n", ipoin,ientt);
          }

          int iopen;
          bool imani = false;
          int ierro = 0,itmp = 0;
          if constexpr (idim == 2){
            intAr1 dum;
            ierro = ball2(*msh,ipoin,ientt,lball,dum,&iopen,&imani,itag2);
          }else{
            ierro = ball3(*msh,ipoin,ientt,lball,&iopen,itag2);
          }
          METRIS_ASSERT(ierro == 0);
          METRIS_ASSERT(iopen == 0);
          METRIS_ASSERT(imani == true);

          double coor0[idim];
          double met0[nnmet];
          for(int ii = 0; ii < idim; ii++) coor0[ii] = msh->coord(ipoin,ii);
          for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh->met(ipoin,ii);
          double qnrm0, qmax0, qnrm1, qmax1;
          try{
            //ierro = smooballdirect<MFT,idim,ideg>(msh,ipoin,lball,qball,
            //                       &qnrm0,&qmax0,&qnrm1,&qmax1,
            //                       qpower,qpnorm,difto,maxwt,inorm,iverb,itag2);
            ierro = smooballdiff<MFT,idim,ideg>(*msh,ipoin,lball,qball,
                                   &qnrm0,&qmax0,&qnrm1,&qmax1);//maxinc_worst,
            if(qmax1 > *qmax){
              if(iverb >= 2) printf("  - reject move, worst above last worst "
                " {:15.7e} > {:15.7e}\n", qmax1, *qmax);
              for(int ii = 0; ii < idim; ii++)
                msh->coord(ipoin,ii) = coor0[ii];
              for(int ii = 0; ii < nnmet;ii++)
                msh->met(ipoin,ii)   =  met0[ii];
              ierro = 1;
            }
          }catch(const MetrisExcept &e){
            printf("## FAILED  smooballdirect : METRIC INVALID\n");
            writeMesh("smooth_error.meshb",*msh);
            throw(e);
          }
          if(ierro == 0){
            (*nsuccthr)[ithread]++;
            if(iverb >= 3) printf("   - success smoothing {} q avg"
                                     " {:10.6e} -> {:10.6e} max {:10.6e} -> {:10.6e}\n",
                                     ipoin,qnrm0,qnrm1,qmax0,qmax1);

            const bool imov = smoothing_neighborhood_should_be_reactivated(
                qnrm0,qnrm1,qmax0,qmax1,tolavg,tolmax);
            if(imov){
              (*nmovthr)[ithread]++;
              reactivateSmoothingRegionGeometry(
                  *msh,idim,ideg,lball,itag_shared);
            }else{
              msh->poi2tag(itag_shared,ipoin) = msh->tag[itag_shared]; // deactivate
            }
          }else{
            msh->poi2tag(itag_shared,ipoin) = msh->tag[itag_shared]; // deactivate
          }

        //} // if qpoin
      } // for ientt
      // Control sizes here if provided (hmin hmax)
    };

    float acc = LaunchParallelMultiArg(LibIdx, LP_elt, LP_poi,
                                       (void*)inerloop_LPlib, 8,
                                       &msh, itag_shared, ithrd2, &qmax,
                                       &nsuccthr, &nmovthr,
                                       tolavg, tolmax);
    CPRINTF1("Smoothing accel = {} \n",acc);
    for(int ii = 0; ii < nthread; ii++){
      nsucc += nsuccthr[ii];
      nmov  += nmovthr[ii];
    }

    #endif
    #endif


    #ifndef USE_LPLIB_SMOOTHINTERIOR

    for(int ipoin = 0; ipoin < msh.npoin; ipoin++){

      // if (ipoin != 15361) continue;

      // skip if dead point or tagged
      if(msh.isdeadpoint(ipoin)) continue;
      if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
      INCVDEPTH(msh.param);

      // check if boundary point
      int ibpoin = msh.poi2bpo[ipoin];
      bool pointOnEdge = false;
      bool pointOnFace = false;
      if (ibpoin >= 0){
        METRIS_ASSERT(msh.bpo2ibi(ibpoin,0) == ipoin);

        // A topological boundary record without a CAD model provides no
        // parameter space in which the point can be moved. Keep that boundary
        // geometry fixed; in particular, reactivating a neighboring P2 control
        // point must not turn it into an unconstrained physical-space variable.
        if(!msh.CAD()) continue;

        #if !defined(SMOOTHEDGES) && !defined(SMOOTFACES)
        continue;
        #endif

        if (msh.bpo2ibi(ibpoin,1) == 0) continue;           // point is a CAD corner, don't move

        #ifdef SMOOTHEDGES
        if (msh.bpo2ibi(ibpoin,1) == 1) pointOnEdge = true; // point is in CAD edge
        #else
        if (msh.bpo2ibi(ibpoin,1) == 1) continue;
        #endif

        #ifdef SMOOTHFACES
        if (msh.bpo2ibi(ibpoin,1) == 2) pointOnFace = true; // point is in CAD face
        #else
        if (msh.bpo2ibi(ibpoin,1) == 2) continue;
        #endif
      }
      const bool interiorPoint = !pointOnEdge && !pointOnFace;

      int ientt = getpoient(msh, ipoin, tdim);
      int iver = tdim == 2 ? msh.template getverfac<ideg>(ientt, ipoin)
                           : msh.template getvertet<ideg>(ientt, ipoin);

      METRIS_ASSERT(iver >= 0);


      CPRINTF1(" - smoo pt {} seed elt {} \n", ipoin,ientt);
      int ierro = 0;

      if(iver < tdim+1){
        int iopen;
        // Vertex case
        if constexpr (idim == 2){
          intAr1 dum;
          bool imani = false;
          ierro = ball2(msh,ipoin,ientt,lball,dum,&iopen,&imani,ithrd2);
          METRIS_ASSERT(imani == true);
        }else{
          ierro = ball3(msh,ipoin,ientt,lball,&iopen,ithrd2);
        }
        // METRIS_ASSERT(iopen == 0);
      }else{
        // HO node
        int nppe = getnnod1(ideg) - 2;
        int ied = (iver - (tdim + 1)) / nppe;
        buildEdgeControlPointSmoothingRegion(
            msh,tdim,ientt,ied,lball);

      }
      METRIS_ASSERT(ierro == 0);

      double oldObjectiveNumerator = 0.;
      double oldObjectiveTargetWeight = 0.;
      int oldObjectiveCount = 0;
      accumulateStepDistanceRegion(
          lball,oldObjectiveNumerator,oldObjectiveCount,
          oldObjectiveTargetWeight);

      double coor0[idim];
      double met0[nnmet];
      for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipoin,ii);
      for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh.met(ipoin,ii);
      const BoundaryParameterSnapshot boundaryParameters0
          = captureBoundaryParameters(msh,ipoin);

      // std::cout << "printing qualities of elements in ball BEFORE smoobaldiff" << std::endl;

      // for (int ieleInBall : lball){

      //   double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,ieleInBall,difto);
      //   std::cout << "ele = " << ieleInBall << " q = " << quael << std::endl;
      // }

      double qnrm0, qmax0, qnrm1, qmax1;
      double newObjectiveNumerator = 0.;
      double newObjectiveTargetWeight = 0.;
      int newObjectiveCount = 0;
      bool objectiveReplacementAccepted = true;
      try{
        //ierro = smooballdirect<MFT,idim,ideg>(msh,ipoin,lball,qball,
        //                       &qnrm0,&qmax0,&qnrm1,&qmax1,
        //                       qpower,qpnorm,difto,maxwt,inorm,iverb,ithrd2);
        if(msh.param->iflag2 == 0){
          if (pointOnEdge)   ierro = smooballdiff_boundary<MFT,idim,ideg>(msh,ipoin,1,lball,&qnrm0,&qmax0,&qnrm1,&qmax1,iquaf);
          if (pointOnFace)   ierro = smooballdiff_boundary<MFT,idim,ideg>(msh,ipoin,2,lball,&qnrm0,&qmax0,&qnrm1,&qmax1,iquaf);
          if (interiorPoint) ierro = smooballdiff<MFT,idim,ideg>(msh,ipoin,lball,&qnrm0,&qmax0,&qnrm1,&qmax1,iquaf);
        }else{
          ierro = smooballdiff_luksan<MFT,idim,ideg>(msh,ipoin,lball,
                                     &qnrm0,&qmax0,&qnrm1,&qmax1,work,iquaf);
        }
        if(ierro == 0){
          accumulateStepDistanceRegion(
              lball,newObjectiveNumerator,newObjectiveCount,
              newObjectiveTargetWeight);
          if(useGlobalStepDistance){
            objectiveReplacementAccepted =
                globalStepDistance.accepts_replacement(
                    oldObjectiveNumerator,oldObjectiveCount,
                    oldObjectiveTargetWeight,
                    newObjectiveNumerator,newObjectiveCount,
                    newObjectiveTargetWeight);
          }else{
            const double objectiveOld = smoothing_region_objective(
                msh,iquaf,oldObjectiveNumerator,
                oldObjectiveTargetWeight);
            const double objectiveNew = smoothing_region_objective(
                msh,iquaf,newObjectiveNumerator,
                newObjectiveTargetWeight);
            objectiveReplacementAccepted =
                objective_strictly_improves(objectiveNew,objectiveOld);
          }
        }
        if(!objectiveReplacementAccepted){
          CPRINTF1(" - reject move: configured StepDistance replacement objective did not improve\n");
          for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii];
          for(int ii = 0; ii < nnmet;ii++) msh.met(ipoin,ii)   =  met0[ii];
          restoreBoundaryParameters(msh,boundaryParameters0);
          ierro = 1;
        // Objective-driven paths have already applied their configured
        // regional or global acceptance contract. Do not layer the legacy
        // worst-element veto on top of that decision.
        }else if(!isObjectiveDrivenSmoothing(iquaf) && qmax1 > qmax){
          CPRINTF1(" - reject move, worst above last worst {:15.7e} > {:15.7e}\n",
                   qmax1, qmax);
          for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii];
          for(int ii = 0; ii < nnmet;ii++) msh.met(ipoin,ii)   =  met0[ii];
          restoreBoundaryParameters(msh,boundaryParameters0);
          ierro = 1;
        }else if(ierro == 0 && pointOnEdge){
          synchronizeCurvePointFaceParameters(msh,ipoin);
        }
        if(ierro == 0 && (pointOnEdge || pointOnFace)){
          bool intersectionFree = true;
          if constexpr(idim == 2 && ideg == 2){
            intersectionFree
                = planarP2BoundaryIsIntersectionFreeAroundPoint(msh,ipoin);
          }
          if(!intersectionFree){
            CPRINTF1(
                " - reject move: curved planar boundary intersection\n");
            for(int ii = 0; ii < idim; ii++){
              msh.coord(ipoin,ii) = coor0[ii];
            }
            for(int ii = 0; ii < nnmet; ii++){
              msh.met(ipoin,ii) = met0[ii];
            }
            restoreBoundaryParameters(msh,boundaryParameters0);
            ierro = 4;
          }else{
            enforceAcceptedBoundarySmoothingInvariants<MFT,idim,ideg>(
                msh,ipoin,lball);
          }
        }
      }catch(const MetrisExcept &e){
        PRINTF("## FAILED  smoothing\n");
        MPRINTF("ipoin = {}", ipoin);
        writeMesh("smooth_error.meshb",msh);
        msh.met.writeMetricFile("smooth_error_metric.solb");
        for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii];
        for(int ii = 0; ii < nnmet;ii++) msh.met(ipoin,ii)   =  met0[ii];
        restoreBoundaryParameters(msh,boundaryParameters0);
        writeMesh("smooth_error_0.meshb",msh);
        msh.met.writeMetricFile("smooth_error_metric_0.solb");
        ierro = 1;
        throw(e);
      }
      if(ierro == 0){
        if(useGlobalStepDistance){
          globalStepDistance.replace(
              oldObjectiveNumerator,oldObjectiveCount,
              oldObjectiveTargetWeight,
              newObjectiveNumerator,newObjectiveCount,
              newObjectiveTargetWeight);
        }
        nsucc++;
        CPRINTF1(" - success smoothing {} q avg {:10.6e} -> {:10.6e} "
                 "max {:10.6e} -> {:10.6e}\n",ipoin,qnrm0,qnrm1,qmax0,qmax1);


        const bool imov = smoothing_neighborhood_should_be_reactivated(
            qnrm0,qnrm1,qmax0,qmax1,tolavg,tolmax);
        if(imov){
          nmov ++;
          reactivateSmoothingRegionGeometry(
              msh,tdim,ideg,lball,ithrd1);
        }else{
          msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1]; // deactivate
        }
      }else{
        msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1]; // deactivate
      }

    } // for ipoin // for ientt
    #endif

    double t1 = get_cpu_time();
    CPRINTF1(" - Iteration end time = {:.2e}s nsuccess = {} nmov = {} \n",
                          t1-t0,nsucc,nmov);
    noper += nmov;

    if(nmov == 0) break;
  } // end for niter

  return noper / (double) nentt;
}


#define BOOST_PP_LOCAL_MACRO(n)\
template double smoothInterior_Ball0<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh,\
                                        QuaFun iquaf, int ithrd1, int ithrd2);\
template double smoothInterior_Ball0<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh,\
                                        QuaFun iquaf, int ithrd1, int ithrd2);\
template double smoothInterior_Ball0<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh,\
                                        QuaFun iquaf, int ithrd1, int ithrd2);\
template double smoothInterior_Ball0<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh,\
                                        QuaFun iquaf, int ithrd1, int ithrd2);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

// =============================================================================================== //
// =============================================================================================== //

template<class MFT>
double smoothElement_Ball(Mesh<MFT> &msh, const int ientt, BadEntHandler& handler, QuaFun iquaf, int ithrd1, int ithrd2){

  int tdimn = msh.get_tdim();

  METRIS_ASSERT_MSG(tdimn > 1, "TODO: edge smooth interior ball");

  // Geo and topo dimn must match otherwise surface specific
  METRIS_ASSERT(tdimn == msh.idim);
  double noper;
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    if(tdimn == 2){
      noper = smoothElement_Ball0<MFT,2,ideg>(msh,ientt,handler,iquaf,ithrd1,ithrd2);
    }else{
      noper = smoothElement_Ball0<MFT,3,ideg>(msh,ientt,handler,iquaf,ithrd1,ithrd2);
    }
  }}CT_FOR1(ideg);

  return noper;
}

template double smoothElement_Ball<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, const int ientt, BadEntHandler& handler,
                                             QuaFun iquaf, int ithrd1, int ithrd2);
template double smoothElement_Ball<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, const int ientt, BadEntHandler& handler,
                                             QuaFun iquaf, int ithrd1, int ithrd2);

// =============================================================================================== //
// =============================================================================================== //


// idim: gdim = tdim
template<class MFT, int idim, int ideg>
double smoothElement_Ball0(Mesh<MFT> &msh, const int ientt, BadEntHandler& handler, QuaFun iquaf,
                           int ithrd1, int ithrd2){
  // TODO: for now just vertex smoothing. Add HO nodes in the future
  GETVDEPTH(msh.param);

  constexpr int tdim = idim;
  //constexpr int gdim = idim;
  if(ideg > idim + 1){
    PRINTF("## SMOOTHING DISABLED FOR DEGREE {} and dim {} \n",ideg,idim);
    return -1.0;
  }

  int nentt = msh.nentt(tdim);
  const intAr2 &ent2poi = msh.ent2poi(tdim);
  const intAr2 &ent2ent = msh.ent2ent(tdim);

  // Eventually move all constants to MetrisParameters
  // L2 conformity error from 0 to 1
  const double difto = 1.0;
  const int miter = msh.param->opt_smoo_niter;
  //const double maxwt = 20.0;
  //const double qrthr = 2.0;
  const double tolavg = msh.param->opt_smoo_tol;
  const double tolmax = msh.param->opt_smoo_tol;

  dblAr1 work;
  if(msh.param->iflag2 != 0){
    MPRINTF("\n\n##WARNING EXPERIMENTAL SMOOTHING FUNCTION 2\n");
  }

  constexpr int nnmet = (idim*(idim+1))/2;

  // METRIS_ENFORCE(msh.param->opt_power < 0); // Otherwise rework the mins / maxs
  // Otherwise not only edge nodes
  METRIS_ENFORCE(ideg <= tdim + 1);

  const int mball = 100;
  intAr1 lball(mball);

  msh.tag[ithrd1]++;

  auto quafun = get_quafun<MFT,idim,idim>(iquaf);

  const bool useGlobalStepDistance =
      iquaf == QuaFun::StepDistance
      && msh.param->step_distance_cavity_target_average;
  StepDistanceObjectiveState globalStepDistance;
  if(useGlobalStepDistance){
    globalStepDistance.element_count = handler.getQualityCount();
    globalStepDistance.numerator = handler.getQualitySum();
    globalStepDistance.target_weight = handler.getQualityCount();
  }

  double noper = 0;
  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh.param);

    #if 0
    for(int ii = 0; ii < msh.npoin; ii++){
      rpoqe[ii] = 0.0;
      ipone[ii] = 0;
    }
    #endif

    int nsucc = 0;
    int nmov  = 0;

    // P2 edge control points are production smoothing variables too. Higher
    // degrees retain the historical vertex-only targeted pass until their
    // face- and volume-interior variable policies are qualified.
    constexpr int candidateNodeCount
        = ideg == 2 ? getnnode(tdim,ideg) : tdim + 1;
    for(int ii = 0; ii < candidateNodeCount; ii++){

      const int ipoin = ent2poi(ientt, ii);

      // skip if dead point or tagged
      if(msh.isdeadpoint(ipoin)) continue;
      if(msh.poi2tag(ithrd1,ipoin) >= msh.tag[ithrd1]) continue;
      INCVDEPTH(msh.param);

     // check if boundary point
      int ibpoin = msh.poi2bpo[ipoin];
      bool pointOnEdge = false;
      bool pointOnFace = false;
      if (ibpoin >= 0){
        METRIS_ASSERT(msh.bpo2ibi(ibpoin,0) == ipoin);

        // Boundary points can only be smoothed through their CAD parameters.
        // Without CAD, reactivation means "reconsider affected geometry", not
        // permission to release a boundary control point into physical space.
        if(!msh.CAD()) continue;

        #if !defined(SMOOTHEDGES) && !defined(SMOOTFACES)
        continue;
        #endif

        if (msh.bpo2ibi(ibpoin,1) == 0) continue;           // point is a CAD corner, don't move

        #ifdef SMOOTHEDGES
        if (msh.bpo2ibi(ibpoin,1) == 1) pointOnEdge = true; // point is in CAD edge
        #else
        if (msh.bpo2ibi(ibpoin,1) == 1) continue;
        #endif

        #ifdef SMOOTHFACES
        if (msh.bpo2ibi(ibpoin,1) == 2) pointOnFace = true; // point is in CAD face
        #else
        if (msh.bpo2ibi(ibpoin,1) == 2) continue;
        #endif
      }
      const bool interiorPoint = !pointOnEdge && !pointOnFace;

      CPRINTF1(" - smoo pt {} seed elt {} \n", ipoin,ientt);
      int ierro = 0;

      if(ii < tdim + 1){
        int iopen;
        if constexpr(idim == 2){
          intAr1 dum;
          bool imani = false;
          ierro = ball2(
              msh,ipoin,ientt,lball,dum,&iopen,&imani,ithrd2);
          METRIS_ASSERT(imani == true);
        }else{
          ierro = ball3(msh,ipoin,ientt,lball,&iopen,ithrd2);
        }
      }else{
        if constexpr(ideg == 2){
          constexpr int edgeNodesPerEdge = getnnod1(ideg) - 2;
          const int localEdge = (ii - (tdim + 1))/edgeNodesPerEdge;
          buildEdgeControlPointSmoothingRegion(
              msh,tdim,ientt,localEdge,lball);
        }else{
          METRIS_THROW_MSG(
              "Only P2 targeted high-order smoothing is qualified");
        }
      }
      #ifndef NDEBUG
      // if (!pointOnEdge)  METRIS_ASSERT(iopen == 0);
      // else               METRIS_ASSERT(iopen == 1);
      METRIS_ASSERT(ierro == 0);
      #endif

      // we now have the ball for ipoin, get initial qualities

      double qmin = 1.0e30,qmax = -1.0e30, qnrm = 0.0;
      double qsum = 0.;
      double targetWeightSum = 0.;
      int imax = -1;
      int navg = 0;
      for(auto ient2 : lball){

        if(isdeadent(ient2,ent2poi)) continue;

        navg++;

        double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,ient2,difto);

        qnrm += quael;
        double regionWeight = 1.;
        if(iquaf == QuaFun::StepDistance
           && msh.param->step_distance_cavity_target_average){
          regionWeight =
              step_distance_element_target_weight<MFT,idim,idim>(
                  msh,AsDeg::Pk,ient2);
          targetWeightSum += regionWeight;
        }
        qsum += regionWeight*quael;
        qmin = MIN(qmin,quael);
        if(qmax < quael){
          imax = ient2;
          qmax = quael;
        }
      }

      qnrm /= navg;
      double t0 = get_cpu_time();
      CPRINTF1(" - smoo iter {:3}, ipoin {} (ientt {}, point {}), ball init {:10.6e} < q < {:10.6e} (at {}), avg = {:10.6e} "
                     "(p = {})\n",niter,ipoin,ientt,ii,qmin,qmax,imax,qnrm,msh.param->opt_pnorm);

      const bool protectEdgeObjective = idim == 2
                                     && pointOnEdge
                                     && msh.CAD()
                                     && msh.param->adp_line_adapt;
      auto incidentEdgeLengthObjective = [&]() -> double{
        if constexpr(idim != 2){
          return 0.0;
        }else{
          if(!protectEdgeObjective) return 0.0;
          double objective = 0.0;
          for(int iedge = 0; iedge < msh.nedge; iedge++){
            if(isdeadent(iedge,msh.edg2poi)) continue;
            if(msh.edg2poi(iedge,0) != ipoin
               && msh.edg2poi(iedge,1) != ipoin) continue;
            if(iquaf == QuaFun::StepDistance){
              objective += metqua1_length<MFT,2,QuaFun::StepDistance>(
                  msh,msh.edg2poi[iedge]);
            }else{
              objective += metqua1_length<MFT,2,QuaFun::SizeShape>(
                  msh,msh.edg2poi[iedge]);
            }
          }
          return objective;
        }
      };
      const double edgeObjectiveOld = incidentEdgeLengthObjective();

      double coor0[idim];
      double met0[nnmet];
      for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipoin,ii);
      for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh.met(ipoin,ii);
      const BoundaryParameterSnapshot boundaryParameters0
          = captureBoundaryParameters(msh,ipoin);
      double qnrm0, qmax0, qnrm1, qmax1;
      double qsumNew = 0.;
      double targetWeightSumNew = 0.;
      try{
        if(msh.param->iflag2 == 0){
          if (pointOnEdge)   ierro = smooballdiff_boundary<MFT,idim,ideg>(msh,ipoin,1,lball,&qnrm0,&qmax0,&qnrm1,&qmax1,iquaf);
          if (pointOnFace)   ierro = smooballdiff_boundary<MFT,idim,ideg>(msh,ipoin,2,lball,&qnrm0,&qmax0,&qnrm1,&qmax1,iquaf);
          if (interiorPoint) ierro = smooballdiff<MFT,idim,ideg>(msh,ipoin,lball,&qnrm0,&qmax0,&qnrm1,&qmax1,iquaf);
        }else{
          ierro = smooballdiff_luksan<MFT,idim,ideg>(msh,ipoin,lball,
                                     &qnrm0,&qmax0,&qnrm1,&qmax1,work,iquaf);
        }
        double qmaxNew = 0.;
        for(auto ient2 : lball){

          if(isdeadent(ient2,ent2poi)) continue;

          double quael = quafun(msh,AsDeg::Pk,AsDeg::Pk,ient2,difto);
          double regionWeight = 1.;
          if(iquaf == QuaFun::StepDistance
             && msh.param->step_distance_cavity_target_average){
            regionWeight =
                step_distance_element_target_weight<MFT,idim,idim>(
                    msh,AsDeg::Pk,ient2);
            targetWeightSumNew += regionWeight;
          }
          qsumNew += regionWeight*quael;
          if (quael > qmaxNew) qmaxNew = quael;

        }
        const double objectiveSum = smoothing_region_objective(
            msh,iquaf,qsum,targetWeightSum);
        const double objectiveSumNew = smoothing_region_objective(
            msh,iquaf,qsumNew,targetWeightSumNew);
        bool improveQuaSum =
            handler.checkSuccess(objectiveSumNew,objectiveSum);
        bool improveGlobalObjective = true;
        if(useGlobalStepDistance && improveQuaSum){
          const double globalObjective = globalStepDistance.value();
          const double globalObjectiveNew = globalStepDistance.replaced_value(
              qsum,navg,targetWeightSum,
              qsumNew,navg,targetWeightSumNew);
          improveGlobalObjective =
              objective_strictly_improves(
                  globalObjectiveNew,globalObjective);
        }
        bool improveQuaMax = true;
        #ifdef IMPROVEMAXQUAL
          if(!isObjectiveDrivenSmoothing(iquaf)){
            improveQuaMax = qmaxNew < qmax;
          }
        #endif
        const double edgeObjectiveNew = incidentEdgeLengthObjective();
        const double edgeTolerance =
            64.0*std::numeric_limits<double>::epsilon()
            * MAX(1.0,std::abs(edgeObjectiveOld));
        const bool improveEdgeObjective = !protectEdgeObjective
            || edgeObjectiveNew <= edgeObjectiveOld + edgeTolerance;
        if(!improveQuaSum || !improveGlobalObjective || !improveQuaMax
           || !improveEdgeObjective){
          CPRINTF1(" - reject move, quality error increased: {:15.7e} > {:15.7e}\n",
                   objectiveSumNew, objectiveSum);
          if(!improveEdgeObjective){
            CPRINTF2(" - reject move, 1D length objective increased: "
                     "{:15.7e} > {:15.7e}\n",
                     edgeObjectiveNew,edgeObjectiveOld);
          }
          for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii];
          for(int ii = 0; ii < nnmet;ii++) msh.met(ipoin,ii)   =  met0[ii];
          restoreBoundaryParameters(msh,boundaryParameters0);

          ierro = 1;
        }else if(ierro == 0 && pointOnEdge){
          synchronizeCurvePointFaceParameters(msh,ipoin);
        }
        if(ierro == 0 && (pointOnEdge || pointOnFace)){
          bool intersectionFree = true;
          if constexpr(idim == 2 && ideg == 2){
            intersectionFree
                = planarP2BoundaryIsIntersectionFreeAroundPoint(msh,ipoin);
          }
          if(!intersectionFree){
            CPRINTF1(
                " - reject move: curved planar boundary intersection\n");
            for(int ii = 0; ii < idim; ii++){
              msh.coord(ipoin,ii) = coor0[ii];
            }
            for(int ii = 0; ii < nnmet; ii++){
              msh.met(ipoin,ii) = met0[ii];
            }
            restoreBoundaryParameters(msh,boundaryParameters0);
            ierro = 4;
          }else{
            enforceAcceptedBoundarySmoothingInvariants<MFT,idim,ideg>(
                msh,ipoin,lball);
          }
        }
      }catch(const MetrisExcept &e){
        PRINTF("## FAILED  smooballdirect\n");
        writeMesh("smooth_error.meshb",msh);
        for(int ii = 0; ii < idim; ii++) msh.coord(ipoin,ii) = coor0[ii];
        for(int ii = 0; ii < nnmet;ii++) msh.met(ipoin,ii)   = met0[ii];
        restoreBoundaryParameters(msh,boundaryParameters0);
        throw(e);
      }
      if(ierro == 0){
        if(useGlobalStepDistance){
          globalStepDistance.replace(
              qsum,navg,targetWeightSum,
              qsumNew,navg,targetWeightSumNew);
        }
        nsucc++;
        CPRINTF1(" - success smoothing {} q avg {:10.6e} -> {:10.6e} "
                 "max {:10.6e} -> {:10.6e}\n",ipoin,qnrm0,qnrm1,qmax0,qmax1);

        nmov++;
        reactivateSmoothingRegionGeometry(
            msh,tdim,ideg,lball,ithrd1);
        for(int ient2 : lball){
          // inform affected entities: id and new quality
          handler.affectedEnttsAlive[ient2] = quafun(msh,AsDeg::Pk,AsDeg::Pk,ient2,difto);
        }
      }else{
        msh.poi2tag(ithrd1,ipoin) = msh.tag[ithrd1]; // deactivate
      }

      double t1 = get_cpu_time();
      CPRINTF1(" - Iteration end time = {:.2e}s nsuccess = {} nmov = {} \n",
                            t1-t0,nsucc,nmov);
      noper += nmov;
      if(nmov == 0) break;
    } // ii (vertex in element)
  } // iter

  return noper / (double) nentt;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template double smoothElement_Ball0<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh,\
                                        const int ientt, BadEntHandler& handler, QuaFun iquaf, int ithrd1, int ithrd2);\
template double smoothElement_Ball0<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh,\
                                        const int ientt, BadEntHandler& handler, QuaFun iquaf, int ithrd1, int ithrd2);\
template double smoothElement_Ball0<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh,\
                                        const int ientt, BadEntHandler& handler, QuaFun iquaf, int ithrd1, int ithrd2);\
template double smoothElement_Ball0<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh,\
                                        const int ientt, BadEntHandler& handler, QuaFun iquaf, int ithrd1, int ithrd2);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

// =============================================================================================== //
// =============================================================================================== //

template<class MFT>
double smoothCavity(Mesh<MFT> &msh, MshCavity& cav, BadEntHandler& handler, QuaFun iquaf,
                    const double quaCav0, const double quaMaxCav0,
                    const double targetWeightCav0,
                    double& quaCav1, double& quaMaxCav1,
                    double& targetWeightCav1,
                    int ithrd1, int ithrd2){

  int tdimn = msh.get_tdim();

  METRIS_ASSERT_MSG(tdimn > 1, "TODO: edge smooth interior ball");

  // Geo and topo dimn must match otherwise surface specific
  METRIS_ASSERT(tdimn == msh.idim);
  double noper;
  CT_FOR0_INC(1,METRIS_MAX_DEG,ideg){if(msh.curdeg == ideg){
    if(tdimn == 2){
      noper = smoothCavity0<MFT,2,ideg>(msh,cav,handler,iquaf,
          quaCav0,quaMaxCav0,targetWeightCav0,
          quaCav1,quaMaxCav1,targetWeightCav1,ithrd1,ithrd2);
    }else{
      noper = smoothCavity0<MFT,3,ideg>(msh,cav,handler,iquaf,
          quaCav0,quaMaxCav0,targetWeightCav0,
          quaCav1,quaMaxCav1,targetWeightCav1,ithrd1,ithrd2);
    }
  }}CT_FOR1(ideg);

  return noper;
}

template double smoothCavity<MetricFieldAnalytical>(Mesh<MetricFieldAnalytical> &msh, MshCavity& cav, BadEntHandler& handler, QuaFun iquaf,
                                                    const double quaCav0, const double quaMaxCav0,
                                                    const double targetWeightCav0,
                                                    double& quaCav1, double& quaMaxCav1,
                                                    double& targetWeightCav1,
                                                    int ithrd1, int ithrd2);
template double smoothCavity<MetricFieldFE        >(Mesh<MetricFieldFE        > &msh, MshCavity& cav, BadEntHandler& handler, QuaFun iquaf,
                                                    const double quaCav0, const double quaMaxCav0,
                                                    const double targetWeightCav0,
                                                    double& quaCav1, double& quaMaxCav1,
                                                    double& targetWeightCav1,
                                                    int ithrd1, int ithrd2);

// =============================================================================================== //
// =============================================================================================== //


// idim: gdim = tdim
template<class MFT, int idim, int ideg>
double smoothCavity0(Mesh<MFT> &msh, MshCavity& cav, BadEntHandler& handler, QuaFun iquaf, const double quaCav0, const double quaMaxCav0, const double targetWeightCav0, double& quaCav1, double& quaMaxCav1, double& targetWeightCav1,
                           int ithrd1, int ithrd2){
  // TODO: for now just vertex smoothing. Add HO nodes in the future
  GETVDEPTH(msh.param);

  constexpr int tdim = idim;
  //constexpr int gdim = idim;
  if(ideg > idim + 1){
    PRINTF("## SMOOTHING DISABLED FOR DEGREE {} and dim {} \n",ideg,idim);
    return -1.0;
  }

  int nentt = msh.nentt(tdim);
  const intAr2 &ent2poi = msh.ent2poi(tdim);
  const intAr2 &ent2ent = msh.ent2ent(tdim);

  // Eventually move all constants to MetrisParameters
  // L2 conformity error from 0 to 1
  const double difto = 1.0;
  // const int miter = msh.param->opt_smoo_niter;
  const int miter = 1;
  //const double maxwt = 20.0;
  //const double qrthr = 2.0;
  const double tolavg = msh.param->opt_smoo_tol;
  const double tolmax = msh.param->opt_smoo_tol;

  constexpr int nnmet = (idim*(idim+1))/2;

  METRIS_ENFORCE(ideg <= tdim + 1);

  msh.tag[ithrd1]++;

  double noper = 0;
  for(int niter = 0; niter < miter; niter++){
    INCVDEPTH(msh.param);

    int nsucc = 0;
    int nmov  = 0;

    const int ipins = cav.ipins;
    quaCav1 = quaCav0;
    quaMaxCav1 = quaMaxCav0;
    targetWeightCav1 = targetWeightCav0;

    CPRINTF1(" - smoo cavity for insertion pt {} \n", ipins);
    int ierro = 0;

    double t0 = get_cpu_time();

    double coor0[idim];
    double met0[nnmet];
    for(int ii = 0; ii < idim; ii++) coor0[ii] = msh.coord(ipins,ii);
    for(int ii = 0; ii < nnmet; ii++) met0[ii] = msh.met(ipins,ii);
    const int ibpoin = msh.poi2bpo[ipins];
    const double tparam0 = ibpoin >= 0 ? msh.bpo2rbi(ibpoin,0) : 0.;
    double qnrm0, qmax0, qnrm1, qmax1;
    try{

      if(ibpoin >= 0){
        int cadDim = msh.bpo2ibi(ibpoin,1);
        ierro = smoocavdiff_boundary<MFT,idim,ideg>(
            msh,cav,cadDim,quaCav1,quaMaxCav1,targetWeightCav1,
            iquaf,ithrd1);
      }else{
        ierro = smoocavdiff<MFT,idim,ideg>(
            msh,cav,quaCav1,quaMaxCav1,targetWeightCav1,
            iquaf,ithrd1);
      }

      // std::cout << "In smoothCavity0" << std::endl;
      // std::cout << "quaCav1 = " << quaCav1 << std::endl;
      // std::cout << "quaMaxCav1 = " << quaMaxCav1 << std::endl;
      // std::cout << "quaCav0 = " << quaCav0 << std::endl;
      // std::cout << "quaMaxCav0 = " << quaMaxCav0 << std::endl;

      if (ierro == 0){
        const double objectiveCav0 = smoothing_region_objective(
            msh,iquaf,quaCav0,targetWeightCav0);
        const double objectiveCav1 = smoothing_region_objective(
            msh,iquaf,quaCav1,targetWeightCav1);
        bool improveCavSum = handler.checkSuccess(objectiveCav1,objectiveCav0);
        bool improveCavMax = true;
        bool pointMoved = true;
        if constexpr(ideg == 2){
          double displacementSquared = 0.;
          double coordinateScale = 1.;
          for(int component = 0; component < idim; component++){
            const double displacement
                = msh.coord(ipins,component) - coor0[component];
            displacementSquared += displacement*displacement;
            coordinateScale = MAX(
                coordinateScale,
                MAX(std::abs(msh.coord(ipins,component)),
                    std::abs(coor0[component])));
          }
          const double movementTolerance
              = 128.*std::numeric_limits<double>::epsilon()*coordinateScale;
          pointMoved
              = displacementSquared > movementTolerance*movementTolerance;
        }
        #ifdef IMPROVEMAXQUAL
        if(!isObjectiveDrivenSmoothing(iquaf)){
          improveCavMax = quaMaxCav1 <= quaMaxCav0;
        }
        #endif
        if(!improveCavSum || !improveCavMax || !pointMoved){
          CPRINTF1(" - reject move, quality error increased ({:15.7e} > {:15.7e}) or max error increased ({:15.7e} > {:15.7e})\n",
                    quaCav1, quaCav0, quaMaxCav1, quaMaxCav0);
          for(int ii = 0; ii < idim; ii++) msh.coord(ipins,ii) = coor0[ii];
          for(int ii = 0; ii < nnmet;ii++) msh.met(ipins,ii)   =  met0[ii];
          if(ibpoin >= 0) msh.bpo2rbi(ibpoin,0) = tparam0;

          quaCav1 = quaCav0;
          quaMaxCav1 = quaMaxCav0;
          targetWeightCav1 = targetWeightCav0;
          ierro = 1;
        }else if constexpr(ideg == 2){
          // A successful interior relocation is no longer a subdivision point
          // of the original edge. Future growth and final reconstruction must
          // therefore use the same released affine-spoke family optimized
          // above. CAD-constrained points retain endpoint provenance because
          // their child geometry remains CAD-owned.
          if(ibpoin < 0) cav.clear_split_edge_geometry();
        }
      }
      else{
        for(int ii = 0; ii < idim; ii++) msh.coord(ipins,ii) = coor0[ii];
        for(int ii = 0; ii < nnmet;ii++) msh.met(ipins,ii)   =  met0[ii];
        if(ibpoin >= 0) msh.bpo2rbi(ibpoin,0) = tparam0;
        targetWeightCav1 = targetWeightCav0;
      }
    }catch(const MetrisExcept &e){
      PRINTF("## FAILED  smoocavdiff\n");
      writeMesh("smooth_error.meshb",msh);
      throw(e);
    }
    if(ierro == 0){
      nsucc++;
      CPRINTF1(" - success smoothing {} q avg {:10.6e} -> {:10.6e} "
                "max {:10.6e} -> {:10.6e}\n",ipins,quaCav0,quaCav1,quaMaxCav0,quaMaxCav1);

    }

    double t1 = get_cpu_time();
    CPRINTF1(" - Iteration end time = {:.2e}s nsuccess = {} nmov = {} \n",
                          t1-t0,nsucc,nmov);
    noper += nsucc;
    if(nsucc == 0) break;
  } // iter

  return noper / (double) nentt;
}

#define BOOST_PP_LOCAL_MACRO(n)\
template double smoothCavity0<MetricFieldAnalytical,2,n>(Mesh<MetricFieldAnalytical> &msh,\
                                        MshCavity& cav, BadEntHandler& handler, QuaFun iquaf, const double quaCav0, const double quaMaxCav0, const double targetWeightCav0, double& quaCav1, double& quaMaxCav1, double& targetWeightCav1, \
                           int ithrd1, int ithrd2);\
template double smoothCavity0<MetricFieldFE        ,2,n>(Mesh<MetricFieldFE        > &msh,\
                                        MshCavity& cav, BadEntHandler& handler, QuaFun iquaf, const double quaCav0, const double quaMaxCav0, const double targetWeightCav0, double& quaCav1, double& quaMaxCav1, double& targetWeightCav1, \
                           int ithrd1, int ithrd2);\
template double smoothCavity0<MetricFieldAnalytical,3,n>(Mesh<MetricFieldAnalytical> &msh,\
                                        MshCavity& cav, BadEntHandler& handler, QuaFun iquaf, const double quaCav0, const double quaMaxCav0, const double targetWeightCav0, double& quaCav1, double& quaMaxCav1, double& targetWeightCav1, \
                           int ithrd1, int ithrd2);\
template double smoothCavity0<MetricFieldFE        ,3,n>(Mesh<MetricFieldFE        > &msh,\
                                        MshCavity& cav, BadEntHandler& handler, QuaFun iquaf, const double quaCav0, const double quaMaxCav0, const double targetWeightCav0, double& quaCav1, double& quaMaxCav1, double& targetWeightCav1, \
                           int ithrd1, int ithrd2);
#define BOOST_PP_LOCAL_LIMITS     (1, METRIS_MAX_DEG)
#include BOOST_PP_LOCAL_ITERATE()

} // end namespace
