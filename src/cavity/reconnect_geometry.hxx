//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1

#ifndef __METRIS_CAVITY_RECONNECT_GEOMETRY_HXX__
#define __METRIS_CAVITY_RECONNECT_GEOMETRY_HXX__

#include "msh_cavity.hxx"

#include "../Mesh/Mesh.hxx"
#include "../low_eval.hxx"

namespace Metris{

enum class QuadraticConeSpokePolicy
{
  PreserveSplitEdgeGeometry,
  ReleasedAffine
};

template<int tdim>
constexpr int quadratic_simplex_edge_count()
{
  static_assert(tdim == 2 || tdim == 3);
  return tdim == 2 ? 3 : 6;
}

template<int tdim>
int quadratic_simplex_edge_vertex(int edge, int endpoint)
{
  static_assert(tdim == 2 || tdim == 3);
  if constexpr(tdim == 2) return lnoed2[edge][endpoint];
  return lnoed3[edge][endpoint];
}

template<int tdim>
int quadratic_simplex_edge_node(int edge)
{
  int index[tdim + 1] = {};
  index[quadratic_simplex_edge_vertex<tdim>(edge,0)] = 1;
  index[quadratic_simplex_edge_vertex<tdim>(edge,1)] = 1;
  return mul2nod(tdim,index);
}

// Initialize the single interior coefficient of a quadratic child edge when
// an existing non-CAD polynomial edge is split. The result is the exact
// restriction of the parent curve in the mesh's current basis. Returning
// false means this is a genuinely new edge, for which the caller retains its
// ordinary affine initialization.
template<class MFT, int gdim>
bool initialize_quadratic_split_child_coefficient(
    const MeshMetric<MFT>& msh,
    const MshCavity& cav,
    int endpoint0,
    int endpoint1,
    double* coefficient){
  if(!cav.preserve_split_edge_geometry
  || cav.split_edge_points.get_n() != getnnode(1,2)) return false;

  int parent_endpoint = -1;
  if(endpoint0 == cav.ipins) parent_endpoint = endpoint1;
  if(endpoint1 == cav.ipins) parent_endpoint = endpoint0;
  if(parent_endpoint < 0) return false;

  // Metris stores the two vertices first and the P2 edge-interior coefficient
  // last, rather than in increasing one-dimensional parameter order.
  const int parent0 = cav.split_edge_points[0];
  const int parent2 = cav.split_edge_points[1];
  const int parent1 = cav.split_edge_points[2];
  if(parent_endpoint != parent0 && parent_endpoint != parent2) return false;

  const double bary0 = cav.split_edge_barycentric[0];
  const double bary1 = cav.split_edge_barycentric[1];
  if(msh.getBasis() == FEBasis::Lagrange){
    const double parent_bary[2] = {
        0.5*(bary0 + (parent_endpoint == parent0 ? 1. : 0.)),
        0.5*(bary1 + (parent_endpoint == parent2 ? 1. : 0.))};
    eval1<gdim,2>(msh.coord,&cav.split_edge_points[0],
                  FEBasis::Lagrange,DifVar::None,DifVar::None,
                  parent_bary,coefficient,nullptr,nullptr);
    return true;
  }

  METRIS_ENFORCE(msh.getBasis() == FEBasis::Bezier);
  // De Casteljau at the split parameter. bary0 is the weight of parent0.
  // The adjacent first-level control point is the middle Bezier coefficient
  // of the corresponding restricted child curve.
  for(int component = 0; component < gdim; component++){
    if(parent_endpoint == parent0){
      coefficient[component]
          = bary0*msh.coord(parent0,component)
          + bary1*msh.coord(parent1,component);
    }else{
      coefficient[component]
          = bary0*msh.coord(parent1,component)
          + bary1*msh.coord(parent2,component);
    }
  }
  return true;
}

// Evaluate the coefficient at the midpoint of a CAD parameter interval shared
// by two endpoints. A parent-edge child may use a CAD curve; a genuinely new
// boundary spoke is owned only by a common CAD surface.
template<class MFT>
bool evaluate_common_cad_midpoint(const MeshMetric<MFT>& msh,
                                  int point0,
                                  int point1,
                                  int minimum_owner_dimension,
                                  double* coefficient)
{
  if(!msh.CAD() || msh.getBasis() != FEBasis::Lagrange) return false;

  for(int owner_dimension = minimum_owner_dimension;
      owner_dimension <= 2 && owner_dimension < msh.idim;
      owner_dimension++){
    for(int boundary_point = msh.poi2bpo[point0]; boundary_point >= 0;
        boundary_point = msh.bpo2ibi(boundary_point,3)){
      if(msh.bpo2ibi(boundary_point,1) != owner_dimension) continue;
      const int owner_entity = msh.bpo2ibi(boundary_point,2);
      if(owner_entity < 0 || owner_entity >= msh.nentt(owner_dimension)){
        continue;
      }
      const int owner_reference = msh.ent2ref(owner_dimension)[owner_entity];
      const int other_boundary_point = msh.poi2ebp(
          point1,owner_dimension,-1,owner_reference);
      if(other_boundary_point < 0) continue;

      double parameter[2] = {0.,0.};
      for(int component = 0; component < owner_dimension; component++){
        parameter[component]
            = 0.5*(msh.bpo2rbi(boundary_point,component)
                  + msh.bpo2rbi(other_boundary_point,component));
      }
      const ego owner = owner_dimension == 1
          ? msh.CAD.cad2edg[owner_reference]
          : msh.CAD.cad2fac[owner_reference];
      double result[18];
      const int error = EG_evaluate(owner,parameter,result);
      METRIS_ENFORCE_MSG(error == 0,
          "EG_evaluate failed with error {} for local P2 cavity probe",
          error);
      for(int component = 0; component < msh.idim; component++){
        coefficient[component] = result[component];
      }
      return true;
    }
  }
  return false;
}

// Complete one quadratic cone element from its already installed vertices.
// The edge opposite cav.ipins is inherited from source_element. Spokes are
// CAD-owned where applicable; otherwise they either retain exact split-edge
// restriction or are released to the affine P2 representation selected by the
// caller.
template<class MFT, int gdim, int tdim>
void complete_quadratic_cone_element(
    MeshMetric<MFT>& msh,
    const MshCavity& cav,
    int source_element,
    int candidate_element,
    QuadraticConeSpokePolicy spoke_policy
        = QuadraticConeSpokePolicy::PreserveSplitEdgeGeometry)
{
  static_assert(gdim == tdim);
  METRIS_ASSERT(msh.curdeg == 2);
  intAr2& element_to_point = msh.ent2poi(tdim);

  for(int candidate_edge = 0;
      candidate_edge < quadratic_simplex_edge_count<tdim>();
      candidate_edge++){
    const int local0
        = quadratic_simplex_edge_vertex<tdim>(candidate_edge,0);
    const int local1
        = quadratic_simplex_edge_vertex<tdim>(candidate_edge,1);
    const int point0 = element_to_point(candidate_element,local0);
    const int point1 = element_to_point(candidate_element,local1);
    const int coefficient_point = element_to_point(
        candidate_element,
        quadratic_simplex_edge_node<tdim>(candidate_edge));

    if(point0 != cav.ipins && point1 != cav.ipins){
      bool copied_existing_edge = false;
      for(int source_edge = 0;
          source_edge < quadratic_simplex_edge_count<tdim>(); source_edge++){
        const int source0 = element_to_point(
            source_element,
            quadratic_simplex_edge_vertex<tdim>(source_edge,0));
        const int source1 = element_to_point(
            source_element,
            quadratic_simplex_edge_vertex<tdim>(source_edge,1));
        if(!((source0 == point0 && source1 == point1)
          || (source0 == point1 && source1 == point0))) continue;
        const int source_coefficient = element_to_point(
            source_element,quadratic_simplex_edge_node<tdim>(source_edge));
        for(int component = 0; component < gdim; component++){
          msh.coord(coefficient_point,component)
              = msh.coord(source_coefficient,component);
        }
        copied_existing_edge = true;
        break;
      }
      METRIS_ENFORCE_MSG(copied_existing_edge,
          "Local P2 cavity probe could not find inherited edge {} {} "
          "in source element {}",point0,point1,source_element);
      continue;
    }

    const int other_point = point0 == cav.ipins ? point1 : point0;
    const bool is_parent_endpoint
        = cav.split_edge_points.get_n() >= 2
       && (other_point == cav.split_edge_points[0]
        || other_point == cav.split_edge_points[1]);

    // CAD owns physical boundary geometry under either spoke policy. A curve
    // is authoritative only for a true child of the split edge; a common CAD
    // surface may own any newly created boundary spoke.
    const int minimum_cad_dimension = is_parent_endpoint ? 1 : 2;
    bool initialized = evaluate_common_cad_midpoint(
        msh,cav.ipins,other_point,minimum_cad_dimension,
        msh.coord[coefficient_point]);
    if(!initialized
       && spoke_policy
          == QuadraticConeSpokePolicy::PreserveSplitEdgeGeometry){
      initialized = initialize_quadratic_split_child_coefficient<MFT,gdim>(
          msh,cav,point0,point1,msh.coord[coefficient_point]);
    }
    if(!initialized){
      for(int component = 0; component < gdim; component++){
        msh.coord(coefficient_point,component)
            = 0.5*(msh.coord(point0,component)
                  + msh.coord(point1,component));
      }
    }
  }
}

} // namespace Metris

#endif
