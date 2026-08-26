//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

/*
Ball-based as-P1 smoothing.
Low-level drivers:
 - low_smoodirect -> facets propose points, average
 - low_smooballnewton -> Newton
*/

#ifndef __METRIS_MSH_SMOOBAL__
#define __METRIS_MSH_SMOOBAL__

#include "../Mesh/MeshFwd.hxx"
#include "../quality/quafun.hxx"
#include "../aux_badEntHandler.hxx"
#include "../cavity/msh_cavity.hxx"

namespace Metris{

// Production uses its configured objective at every supported geometry
// degree. The geometric arguments remain explicit so future policies can make
// a dimension- or degree-specific choice without changing callers.
inline constexpr QuaFun productionSmoothingObjective(
    [[maybe_unused]] int geometricDimension,
    [[maybe_unused]] int geometryDegree) noexcept
{
#ifdef TESTQUALITYALGO
  #ifdef STEPDISTANCE
  return QuaFun::StepDistance;
  #else
  return QuaFun::SizeShape;
  #endif
#else
  return QuaFun::Distortion;
#endif
}

// Objective-driven smoothing owns its acceptance policy through the complete
// regional or global objective. The legacy worst-element veto is a separate
// non-objective policy and must not be layered on top of these objectives.
inline constexpr bool isObjectiveDrivenSmoothing(QuaFun objective) noexcept
{
  return objective == QuaFun::SizeShape
      || objective == QuaFun::StepDistance;
}

// Build the element region affected by moving a high-order edge control point.
// In 2D this is the seed triangle plus its across-edge neighbor when one
// exists; a boundary edge therefore produces a one-triangle region. In 3D it
// is the complete tetrahedral shell around the edge.
void buildEdgeControlPointSmoothingRegion(
    const MeshBase &msh,
    int tdim,
    int seed_entity,
    int local_edge,
    intAr1 &region);

// Make every geometric degree of freedom belonging to an affected element
// eligible for another smoothing visit. P1 reactivates only vertices; P2 also
// reactivates the shared edge-interior control points.
void reactivateSmoothingRegionGeometry(
    MeshBase &msh,
    int tdim,
    int geometry_degree,
    const intAr1 &region,
    int ithread);

// Check the quadratic polynomial boundary edges affected by moving one point
// in a full-dimensional planar P2 mesh. Intersections at a shared topological
// endpoint are allowed; every other contact with another boundary edge is
// rejected. The mesh must use Lagrange geometry so its edge-interior node is
// the physical midpoint sample of the quadratic map.
bool planarP2BoundaryIsIntersectionFreeAroundPoint(
    const MeshBase &msh,
    int moved_point);

// Returns number of operations as double; this is because it may exceed element
// count, something that cannot be anticipated easily by selecting integer
// types.
template<class MetricFieldType>
double smoothInterior_Ball(Mesh<MetricFieldType> &msh,
             QuaFun iquaf, int ithrd1, int ithrd);


// idim: gdim = tdim
template<class MetricFieldType, int idim, int ideg>
double smoothInterior_Ball0(Mesh<MetricFieldType> &msh,
             QuaFun iquaf, int ithrd1, int ithrd);


template<class MetricFieldType>
double smoothElement_Ball(Mesh<MetricFieldType> &msh, const int ientt, BadEntHandler& handler,
                          QuaFun iquaf, int ithrd1, int ithrd);

template<class MetricFieldType, int idim, int ideg>
double smoothElement_Ball0(Mesh<MetricFieldType> &msh, const int ientt, BadEntHandler& handler,
                           QuaFun iquaf, int ithrd1, int ithrd);

// to do cavity smoothing: optimize the insertion point position
template<class MetricFieldType>
double smoothCavity(Mesh<MetricFieldType> &msh, MshCavity& cav, BadEntHandler& handler, QuaFun iquaf,
                    const double quaCav0, const double quaMaxCav0,
                    const double targetWeightCav0,
                    double& quaCav1, double& quaMaxCav1,
                    double& targetWeightCav1,
                    int ithrd1, int ithrd2);

template<class MetricFieldType, int idim, int ideg>
double smoothCavity0(Mesh<MetricFieldType> &msh, MshCavity& cav, BadEntHandler& handler, QuaFun iquaf, const double quaCav0, const double quaMaxCav0, const double targetWeightCav0, double& quaCav1, double& quaMaxCav1, double& targetWeightCav1,
                           int ithrd1, int ithrd2);


} // end namespace
#endif
