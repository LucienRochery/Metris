//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_resolve_edge_refs

#include "common_setup.hxx"

#include <array>
#include <map>
#include <set>
#include <vector>

using namespace Metris;

typedef std::array<double,2> Point2;
typedef std::array<double,3> Point3;
typedef std::array<int,2>    Edge2;
typedef std::array<int,3>    Tri3;


// Surface mesh of the unit cube with a distinct face-ref per box face. Twelve
// triangles (two per face), eight vertices. Triangle vertex orderings are
// chosen so each face's outward normal points away from the cube interior;
// the choice doesn't matter for the connectivity work resolveEdgeRefs does.
static std::vector<Point3> cube_coord(){
  return {
    Point3{0.0,0.0,0.0}, // 0
    Point3{1.0,0.0,0.0}, // 1
    Point3{1.0,1.0,0.0}, // 2
    Point3{0.0,1.0,0.0}, // 3
    Point3{0.0,0.0,1.0}, // 4
    Point3{1.0,0.0,1.0}, // 5
    Point3{1.0,1.0,1.0}, // 6
    Point3{0.0,1.0,1.0}  // 7
  };
}

static std::vector<Tri3> cube_faces(){
  return {
    // z=0 (ref 1)
    Tri3{0,2,1}, Tri3{0,3,2},
    // z=1 (ref 2)
    Tri3{4,5,6}, Tri3{4,6,7},
    // y=0 (ref 3)
    Tri3{0,1,5}, Tri3{0,5,4},
    // y=1 (ref 4)
    Tri3{3,7,6}, Tri3{3,6,2},
    // x=0 (ref 5)
    Tri3{0,4,7}, Tri3{0,7,3},
    // x=1 (ref 6)
    Tri3{1,2,6}, Tri3{1,6,5}
  };
}

static std::vector<int> cube_face_refs(){
  return {1,1, 2,2, 3,3, 4,4, 5,5, 6,6};
}


// Build a MeshBack from the supplied connectivity. Mirrors test_periodic_faces.
static MeshBack *new_surface_mesh(const std::vector<Point3> &coord,
                                  const std::vector<Edge2> &edges,
                                  const std::vector<int>   &edgref,
                                  const std::vector<Tri3>  &faces,
                                  const std::vector<int>   &facref,
                                  const std::vector<int>   &corners){
  BOOST_REQUIRE(edges.size() == edgref.size());
  BOOST_REQUIRE(faces.size() == facref.size());

  MetrisAPI *data = new MetrisAPI(3, 1,
                                  (int) corners.size(), 0, 0,
                                  (int) coord.size(), true,
                                  (int) edges.size(), (int) faces.size(), 0,
                                  FEBasis::Lagrange, FEBasis::Lagrange, MetSpace::Exp);

  for(int ipoin = 0; ipoin < (int) coord.size(); ipoin++){
    data->setCoord(ipoin, &coord[ipoin][0]);
    double met[6] = {1.0,0.0,1.0,0.0,0.0,1.0};
    data->setMetric(ipoin, met);
  }
  for(int iedge = 0; iedge < (int) edges.size(); iedge++){
    data->setElement(1, iedge, &edges[iedge][0], edgref[iedge]);
  }
  for(int iface = 0; iface < (int) faces.size(); iface++){
    data->setElement(2, iface, &faces[iface][0], facref[iface]);
  }
  for(int icorn = 0; icorn < (int) corners.size(); icorn++){
    data->setCorner(icorn, corners[icorn]);
  }

  MetrisParameters *param = new MetrisParameters;
  param->iverb = 0;
  param->ivdepth = 0;

  MeshBack *msh = new MeshBack;
  msh->initialize(data, *param);
  return msh;
}


static bool is_corner(const MeshBase &msh, int ipoin){
  int ibpoi = msh.poi2bpo[ipoin];
  return ibpoi >= 0 && msh.bpo2ibi(ibpoi,1) == 0;
}

// Walk the edge graph from a seed edge and return the set of edge indices in
// the connected component. A component is one curve: it stops where there is
// no neighbour (edg2edg < 0) and at corners, which bound curves. Matches what
// resolveEdgeRefs walks internally; duplicating the walk in the test lets us
// check the post-condition without trusting the SUT to report it.
static std::set<int> connected_component(const MeshBase &msh, int ie0){
  std::set<int> seen;
  std::vector<int> stack;
  stack.push_back(ie0);
  seen.insert(ie0);
  while(!stack.empty()){
    int ie = stack.back();
    stack.pop_back();
    for(int ive = 0; ive < 2; ive++){
      if(is_corner(msh, msh.edg2poi(ie, ive))) continue;
      int ien = msh.edg2edg(ie, 1-ive);
      if(ien < 0) continue;
      if(seen.find(ien) != seen.end()) continue;
      seen.insert(ien);
      stack.push_back(ien);
    }
  }
  return seen;
}

// Verify the post-condition resolveEdgeRefs is supposed to enforce:
//  1. every active edge has a non-negative ref;
//  2. within each curve (connected component up to corners), all edges share
//     the same ref.
static void check_post_condition(const MeshBase &msh){
  std::vector<int> visited(msh.nedge, 0);
  for(int ie = 0; ie < msh.nedge; ie++){
    if(isdeadent(ie, msh.edg2poi)) continue;
    BOOST_TEST_REQUIRE(msh.edg2ref[ie] >= 0);
    // The convention walkPaintEdgeComponent relies on to know which point it
    // is about to walk across: slot ive holds the neighbour through the
    // vertex of index 1 - ive.
    for(int ive = 0; ive < 2; ive++){
      int ien = msh.edg2edg(ie, ive);
      if(ien < 0) continue;
      int ipoin = msh.edg2poi(ie, 1 - ive);
      BOOST_TEST_REQUIRE((msh.edg2poi(ien,0) == ipoin || msh.edg2poi(ien,1) == ipoin));
    }
  }
  for(int ie0 = 0; ie0 < msh.nedge; ie0++){
    if(isdeadent(ie0, msh.edg2poi)) continue;
    if(visited[ie0]) continue;
    std::set<int> comp = connected_component(msh, ie0);
    int ref0 = msh.edg2ref[ie0];
    for(int ie : comp){
      visited[ie] = 1;
      BOOST_TEST_REQUIRE(msh.edg2ref[ie] == ref0);
    }
  }
}


// Build a 2D MeshBack: triangles are the elements, edges are the boundary.
static MeshBack *new_2d_mesh(const std::vector<Point2> &coord,
                             const std::vector<Edge2>  &edges,
                             const std::vector<int>    &edgref,
                             const std::vector<Tri3>   &faces,
                             const std::vector<int>    &facref,
                             const std::vector<int>    &corners){
  BOOST_REQUIRE(edges.size() == edgref.size());
  BOOST_REQUIRE(faces.size() == facref.size());

  MetrisAPI *data = new MetrisAPI(2, 1,
                                  (int) corners.size(), 0, 0,
                                  (int) coord.size(), true,
                                  (int) edges.size(), (int) faces.size(), 0,
                                  FEBasis::Lagrange, FEBasis::Lagrange, MetSpace::Exp);

  for(int ipoin = 0; ipoin < (int) coord.size(); ipoin++){
    data->setCoord(ipoin, &coord[ipoin][0]);
    double met[3] = {1.0,0.0,1.0};
    data->setMetric(ipoin, met);
  }
  for(int iedge = 0; iedge < (int) edges.size(); iedge++){
    data->setElement(1, iedge, &edges[iedge][0], edgref[iedge]);
  }
  for(int iface = 0; iface < (int) faces.size(); iface++){
    data->setElement(2, iface, &faces[iface][0], facref[iface]);
  }
  for(int icorn = 0; icorn < (int) corners.size(); icorn++){
    data->setCorner(icorn, corners[icorn]);
  }

  MetrisParameters *param = new MetrisParameters;
  param->iverb = 0;
  param->ivdepth = 0;

  MeshBack *msh = new MeshBack;
  msh->initialize(data, *param);
  return msh;
}


// The unit square: four sides, four refs, meeting at four corners of degree 2.
// This is the case the cube misses entirely — there every ridge endpoint is a
// triple junction, so edg2edg is already negative and the walk cannot cross it.
// At a degree-2 corner edg2edg links straight through, so the walk has to stop
// on the corner itself or it wanders into the next side and its different ref.
// Run with and without the corners declared in the input: undeclared, they get
// reconstructed from the ref change.
static void check_square(bool declare_corners){
  std::vector<Point2> coord = {
    Point2{0.0,0.0}, Point2{1.0,0.0}, Point2{1.0,1.0}, Point2{0.0,1.0}
  };
  std::vector<Edge2> edges  = {Edge2{0,1}, Edge2{1,2}, Edge2{2,3}, Edge2{3,0}};
  std::vector<int>   edgref = {0, 1, 2, 3};
  std::vector<Tri3>  faces  = {Tri3{0,1,2}, Tri3{0,2,3}};
  std::vector<int>   facref = {0, 0};
  std::vector<int>   corners;
  if(declare_corners) corners = {0,1,2,3};

  MeshBack *msh = new_2d_mesh(coord, edges, edgref, faces, facref, corners);

  BOOST_REQUIRE_EQUAL(msh->nedge, 4);
  check_post_condition(*msh);

  // Every side keeps the ref it came in with: no propagation across a corner.
  for(int ie = 0; ie < msh->nedge; ie++){
    BOOST_TEST_REQUIRE(msh->edg2ref[ie] == edgref[ie]);
  }
  for(int ipoin = 0; ipoin < 4; ipoin++){
    BOOST_TEST_REQUIRE(is_corner(*msh, ipoin));
  }
}

BOOST_AUTO_TEST_CASE(square_four_refs_no_corners)
{
  check_square(false);
}

BOOST_AUTO_TEST_CASE(square_four_refs_with_corners)
{
  check_square(true);
}


// Case 1: cube surface, no input edges, no input corners. resolveEdgeRefs has
// to fabricate refs for every reconstructed edge from scratch. With six
// distinct face refs the cube has twelve ridge edges, each its own connected
// component (every ridge endpoint is a triple junction → corner), so every
// edge ends up with a distinct fresh ref.
BOOST_AUTO_TEST_CASE(no_input_edges_no_corners)
{
  std::vector<Point3> coord  = cube_coord();
  std::vector<Edge2>  edges  = {};
  std::vector<int>    edgref = {};
  std::vector<Tri3>   faces  = cube_faces();
  std::vector<int>    facref = cube_face_refs();
  std::vector<int>    corners = {};

  MeshBack *msh = new_surface_mesh(coord, edges, edgref, faces, facref, corners);

  // 12 box edges should be reconstructed (one per pair of adjacent face refs).
  BOOST_REQUIRE_EQUAL(msh->nedge, 12);

  check_post_condition(*msh);

  // Each ridge connects two cube corners (degree-3 junctions), so each edge
  // is its own component → 12 distinct refs.
  std::set<int> refs_used;
  for(int ie = 0; ie < msh->nedge; ie++){
    if(isdeadent(ie, msh->edg2poi)) continue;
    refs_used.insert(msh->edg2ref[ie]);
  }
  BOOST_TEST_REQUIRE(refs_used.size() == 12u);
}


// Case 2: a few seed edges with coherent refs, NO explicit corners. The seeds
// each carry a distinct ref; resolveEdgeRefs should keep those on the seeded
// component (one edge per component, between corner endpoints) and allocate
// fresh refs to the remaining nine reconstructed ridges.
BOOST_AUTO_TEST_CASE(some_input_edges_no_corners)
{
  std::vector<Point3> coord = cube_coord();
  std::vector<Edge2> edges = {
    Edge2{0,1},   // bottom front edge   ref 100
    Edge2{2,3},   // bottom back edge    ref 101
    Edge2{4,5}    // top    front edge   ref 102
  };
  std::vector<int> edgref = {100, 101, 102};
  std::vector<Tri3> faces = cube_faces();
  std::vector<int> facref = cube_face_refs();
  std::vector<int> corners = {};

  MeshBack *msh = new_surface_mesh(coord, edges, edgref, faces, facref, corners);

  BOOST_REQUIRE_EQUAL(msh->nedge, 12);
  check_post_condition(*msh);

  // The three pre-existing refs should be preserved verbatim somewhere in
  // edg2ref. Each ridge between cube corners is its own connected component,
  // so the three seeds remain isolated and each propagates onto exactly one
  // edge.
  std::map<int,int> ref_count;
  for(int ie = 0; ie < msh->nedge; ie++){
    if(isdeadent(ie, msh->edg2poi)) continue;
    ref_count[msh->edg2ref[ie]]++;
  }
  BOOST_TEST_REQUIRE(ref_count.count(100) == 1u);
  BOOST_TEST_REQUIRE(ref_count.count(101) == 1u);
  BOOST_TEST_REQUIRE(ref_count.count(102) == 1u);
  BOOST_TEST_REQUIRE(ref_count[100] == 1);
  BOOST_TEST_REQUIRE(ref_count[101] == 1);
  BOOST_TEST_REQUIRE(ref_count[102] == 1);

  // 12 components total → 12 distinct refs.
  BOOST_TEST_REQUIRE(ref_count.size() == 12u);

  // Freshly allocated refs must avoid colliding with the pre-existing 100-102.
  for(auto &kv : ref_count){
    if(kv.first == 100 || kv.first == 101 || kv.first == 102) continue;
    BOOST_TEST_REQUIRE(kv.first >= 103);
  }
}


// Case 3: a few seed edges with coherent refs, AND every cube corner declared
// explicitly as a Metris corner. The reconstructed-corner detection should
// then be redundant; resolveEdgeRefs's component walk should still see each
// ridge as its own component and produce the same outcome as Case 2.
BOOST_AUTO_TEST_CASE(some_input_edges_with_corners)
{
  std::vector<Point3> coord = cube_coord();
  std::vector<Edge2> edges = {
    Edge2{0,1},
    Edge2{2,3},
    Edge2{4,5}
  };
  std::vector<int> edgref = {100, 101, 102};
  std::vector<Tri3> faces = cube_faces();
  std::vector<int> facref = cube_face_refs();
  std::vector<int> corners = {0,1,2,3,4,5,6,7};

  MeshBack *msh = new_surface_mesh(coord, edges, edgref, faces, facref, corners);

  BOOST_REQUIRE_EQUAL(msh->nedge, 12);
  check_post_condition(*msh);

  std::map<int,int> ref_count;
  for(int ie = 0; ie < msh->nedge; ie++){
    if(isdeadent(ie, msh->edg2poi)) continue;
    ref_count[msh->edg2ref[ie]]++;
  }
  BOOST_TEST_REQUIRE(ref_count[100] == 1);
  BOOST_TEST_REQUIRE(ref_count[101] == 1);
  BOOST_TEST_REQUIRE(ref_count[102] == 1);
  BOOST_TEST_REQUIRE(ref_count.size() == 12u);
}
