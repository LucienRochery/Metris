//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1

#ifndef __METRIS_CAVITY_RECONNECT_GEOMETRY_HXX__
#define __METRIS_CAVITY_RECONNECT_GEOMETRY_HXX__

#include "msh_cavity.hxx"

#include "../Mesh/Mesh.hxx"
#include "../low_eval.hxx"

namespace Metris{

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

} // namespace Metris

#endif
