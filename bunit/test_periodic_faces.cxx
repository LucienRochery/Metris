//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_periodic_faces

#include "common_setup.hxx"

#include <array>
#include <vector>

using namespace Metris;

typedef std::array<double,3> Point3;
typedef std::array<int,2> Edge2;
typedef std::array<int,3> Tri3;


MeshBack *new_surface_mesh(const std::vector<Point3> &coord,
                           const std::vector<Edge2> &edges,
                           const std::vector<int> &edgref,
                           const std::vector<Tri3> &faces,
                           const std::vector<int> &facref){
  BOOST_REQUIRE(edges.size() == edgref.size());
  BOOST_REQUIRE(faces.size() == facref.size());

  // MeshBase keeps a raw param pointer, and initialization moves MetrisAPI
  // storage into the mesh arrays. Keep this minimal fixture alive for the
  // whole test process rather than tearing it down piecemeal.
  MetrisAPI *data = new MetrisAPI(3,1,0,0,0,
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

  MetrisParameters *param = new MetrisParameters;
  param->iverb = 0;
  param->ivdepth = 0;

  MeshBack *msh = new MeshBack;
  msh->initialize(data, *param);
  return msh;
}


std::vector<Point3> cube_coord(){
  return {
    Point3{0.0,0.0,0.0},
    Point3{1.0,0.0,0.0},
    Point3{1.0,1.0,0.0},
    Point3{0.0,1.0,0.0},
    Point3{0.0,0.0,1.0},
    Point3{1.0,0.0,1.0},
    Point3{1.0,1.0,1.0},
    Point3{0.0,1.0,1.0}
  };
}


BOOST_AUTO_TEST_CASE(test_periodic_faces_open_surface)
{
  std::vector<Point3> coord = cube_coord();
  std::vector<Edge2> edges;
  std::vector<int> edgref;
  std::vector<Tri3> faces = {
    Tri3{0,2,1}, Tri3{0,3,2},
    Tri3{0,1,5}, Tri3{0,5,4},
    Tri3{3,6,2}, Tri3{3,7,6},
    Tri3{0,4,7}, Tri3{0,7,3},
    Tri3{1,2,6}, Tri3{1,6,5}
  };
  std::vector<int> facref = {
    0,0,
    1,1,
    2,2,
    3,3,
    4,4
  };

  MeshBack *msh = new_surface_mesh(coord, edges, edgref, faces, facref);

  BOOST_REQUIRE(msh->isperiodic_face.get_n() == 5);
  BOOST_REQUIRE(msh->nperiodic_face == 0);
  for(int iref = 0; iref < msh->isperiodic_face.get_n(); iref++){
    BOOST_REQUIRE(!msh->isperiodic_face[iref]);
  }

}


BOOST_AUTO_TEST_CASE(test_periodic_faces_same_ref_cube)
{
  std::vector<Point3> coord = cube_coord();
  std::vector<Edge2> edges = {
    Edge2{0,1}, Edge2{1,2}, Edge2{2,3}, Edge2{3,0},
    Edge2{4,5}, Edge2{5,6}, Edge2{6,7}, Edge2{7,4},
    Edge2{0,4}, Edge2{1,5}, Edge2{2,6}, Edge2{3,7}
  };
  std::vector<int> edgref(edges.size(),0);
  std::vector<Tri3> faces = {
    Tri3{0,2,1}, Tri3{0,3,2},
    Tri3{4,5,6}, Tri3{4,6,7},
    Tri3{0,1,5}, Tri3{0,5,4},
    Tri3{3,7,6}, Tri3{3,6,2},
    Tri3{0,4,7}, Tri3{0,7,3},
    Tri3{1,2,6}, Tri3{1,6,5}
  };
  std::vector<int> facref(faces.size(),0);

  MeshBack *msh = new_surface_mesh(coord, edges, edgref, faces, facref);

  BOOST_REQUIRE(msh->isperiodic_face.get_n() == 1);
  BOOST_REQUIRE(msh->isperiodic_face[0]);
  BOOST_REQUIRE(msh->nperiodic_face == 1);

}


BOOST_AUTO_TEST_CASE(test_periodic_faces_nonmanifold)
{
  std::vector<Point3> coord = {
    Point3{0.0,0.0,0.0},
    Point3{1.0,0.0,0.0},
    Point3{0.0,1.0,0.0},
    Point3{0.0,0.0,1.0},
    Point3{0.0,-1.0,0.0}
  };
  std::vector<Edge2> edges;
  std::vector<int> edgref;
  std::vector<Tri3> faces = {
    Tri3{0,1,2},
    Tri3{1,0,3},
    Tri3{0,1,4}
  };
  std::vector<int> facref(faces.size(),0);

  MeshBack *msh = new_surface_mesh(coord, edges, edgref, faces, facref);

  BOOST_REQUIRE(msh->isperiodic_face.get_n() == 1);
  BOOST_REQUIRE(msh->isperiodic_face[0]);
  BOOST_REQUIRE(msh->nperiodic_face == 1);

  int iedge = getedgglo(*msh,0,1);
  BOOST_REQUIRE(iedge >= 0);
  intAr1 lshell(4);
  shell2_nm(*msh, iedge, lshell);
  BOOST_REQUIRE(lshell.get_n() == 3);

}
