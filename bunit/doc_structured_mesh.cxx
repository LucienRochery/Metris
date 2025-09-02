//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#define BOOST_TEST_MODULE test_doc_structured_mesh

#include "common_setup.hxx"
#include "doc_structured_mesh.hxx"

using namespace Metris;
typedef MetricFieldAnalytical MFT;

BOOST_AUTO_TEST_CASE(test_doc_structured_mesh)
{

  const int ithread = 0;

  // mesh params
  double xmin = 0; double xmax = 1;
  double ymin = 0; double ymax = 1;

  int nx = 10; int ny = 10;

  // create connectivity arrays
  intAr2 fac2poi;
  intAr1 fac2ref;
  intAr2 edg2poi;
  intAr1 edg2ref;
  dblAr2 coord;
  intAr2 fac2fac;
  intAr2 edg2edg;
  intAr1 geonodes;
  doc_structured_mesh(xmin,xmax,ymin,ymax,nx,ny,fac2poi,fac2ref,edg2poi,edg2ref,coord,fac2fac,edg2edg,geonodes);

  // set up for MetrisAPI
  int idim = 2;
  int ideg = 1;
  int ncorn = 4;
  int ngpoe = 0;
  int ngpof = 0;
  int npoin = coord.get_n(); // points
  bool imet = false;
  int nedge = edg2edg.get_n(); // bnd edges
  int nface = fac2fac.get_n(); // faces (triangles)
  int nelem = 0; // tetra
  FEBasis meshbasis = FEBasis::Lagrange;
  FEBasis metbasis = FEBasis::Undefined;
  MetSpace metspace = MetSpace::Undefined;

  // init API and pass outputs of doc_structured_mesh
  MetrisAPI myAPI(idim,ideg,ncorn,ngpoe,ngpof,npoin,imet,nedge,nface,nelem,meshbasis,metbasis,metspace);
  myAPI.setCoord(std::move(coord)); // pass coord array
  myAPI.setElement(2, std::move(fac2poi), std::move(fac2ref)); // pass face connectivity
  myAPI.setElement(1, std::move(edg2poi), std::move(edg2ref));  // pass edge connectivity
  for (int ii = 0; ii < geonodes.get_n(); ii++) myAPI.setCorner(ii,geonodes[ii]); // pass corners

  MetrisParameters param;

  // MetrisRunner init and run
  MetrisRunner run(&myAPI,&myAPI,param);
  MeshBase &myMsh = *(run.msh_g);

  // visualize
  writeMesh("doc_structured_mesh", myMsh);
}