//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

// Tests for locMesh (mesh-level walking point location).
// We use convex domains (unit square / unit cube) so that seed-independence
// holds: the barycentric walk always reaches the unique containing element
// regardless of where it starts.

#define BOOST_TEST_MODULE test_locMesh

#include "common_setup.hxx"

using namespace Metris;

typedef MetricFieldAnalytical MFT;

// ===========================================================================
// Test 1: 2D localization on the unit square (iso.p1.10k, convex, P1)
//
//  a) Out-of-bounding-box rejection
//  b) Centroid of every triangle located from a fixed seed (elem 0)
//     → LOC_ERR_NOERR, distance < tol, bary ≈ [1/3,1/3,1/3]
//  c) Every mesh vertex located from a fixed seed
//     → LOC_ERR_NOERR, distance < tol
//  d) Seed independence: centroid of 200 elements located from seed 0 and
//     from the last live element → same found element both times
// ===========================================================================
BOOST_AUTO_TEST_CASE(test_locMesh_2D)
{
  cargHandler arg("-in " METRIS_CASES_DIR "/unit/2D/square/iso.p1.10k"
                  " -anamet 1 -verb 0");
  MetrisRunner run(arg.c, arg.v);
  Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
  msh.param->iverb = 0;

  const double tol = 1.0e-10;
  double coopr[2], bary[3];

  // --- a) Out-of-BB rejection -----------------------------------------------
  {
    double coop_out[2] = {-2.0, 0.5};
    int ielem = 0;
    int ierr = locMesh<2,2,1>(msh, &ielem, coop_out, 2,
                              nullptr, -1, nullptr, coopr, bary, tol);
    BOOST_CHECK_EQUAL(ierr, LOC_ERR_OUTBB);
  }

  // --- b) Locate every triangle centroid from seed 0 ------------------------
  {
    int n_loc_fail  = 0;
    int n_dist_fail = 0;
    int n_bary_fail = 0;

    for(int i = 0; i < msh.nface; i++){
      if(isdeadent(i, msh.fac2poi)) continue;

      double cx = 0, cy = 0;
      for(int v = 0; v < 3; v++){
        cx += msh.coord(msh.fac2poi(i, v), 0);
        cy += msh.coord(msh.fac2poi(i, v), 1);
      }
      cx /= 3.0; cy /= 3.0;
      double coop[2] = {cx, cy};

      int ielem = 0;
      int ierr = locMesh<2,2,1>(msh, &ielem, coop, 2,
                                nullptr, -1, nullptr, coopr, bary, tol);
      if(ierr != LOC_ERR_NOERR){ n_loc_fail++;  continue; }

      double dist = std::sqrt(geterrl2<2>(coop, coopr));
      if(dist > 1.0e-6){ n_dist_fail++; continue; }

      // At its own centroid, bary must be [1/3, 1/3, 1/3]
      for(int v = 0; v < 3; v++){
        if(std::abs(bary[v] - 1.0/3.0) > 1.0e-6){ n_bary_fail++; break; }
      }
    }
    BOOST_CHECK_EQUAL(n_loc_fail,  0);
    BOOST_CHECK_EQUAL(n_dist_fail, 0);
    BOOST_CHECK_EQUAL(n_bary_fail, 0);
  }

  // --- c) Locate every vertex from seed 0 -----------------------------------
  // (works because the unit-square domain is convex)
  {
    int n_vtx_fail = 0;
    for(int ip = 0; ip < msh.npoin; ip++){
      if(msh.isdeadpoint(ip)) continue;
      double coop[2] = {msh.coord(ip, 0), msh.coord(ip, 1)};
      int ielem = 0;
      int ierr = locMesh<2,2,1>(msh, &ielem, coop, 2,
                                nullptr, -1, nullptr, coopr, bary, tol);
      if(ierr != LOC_ERR_NOERR){ n_vtx_fail++; continue; }
      double dist = std::sqrt(geterrl2<2>(coop, coopr));
      if(dist > 1.0e-6) n_vtx_fail++;
    }
    BOOST_CHECK_EQUAL(n_vtx_fail, 0);
  }

  // --- d) Seed independence on convex domain --------------------------------
  // Find the last live element to use as an alternative seed.
  int seed_last = -1;
  for(int i = msh.nface - 1; i >= 0; i--){
    if(!isdeadent(i, msh.fac2poi)){ seed_last = i; break; }
  }
  BOOST_REQUIRE(seed_last > 0); // must differ from seed 0

  {
    int n_indep_fail = 0;
    int ntested = 0;
    for(int i = 0; i < msh.nface && ntested < 200; i++){
      if(isdeadent(i, msh.fac2poi)) continue;
      ntested++;

      double cx = 0, cy = 0;
      for(int v = 0; v < 3; v++){
        cx += msh.coord(msh.fac2poi(i, v), 0);
        cy += msh.coord(msh.fac2poi(i, v), 1);
      }
      cx /= 3.0; cy /= 3.0;
      double coop[2] = {cx, cy};

      double bary0[3], bary1[3], coopr0[2], coopr1[2];
      int ielem0 = 0, ielem1 = seed_last;
      int ierr0 = locMesh<2,2,1>(msh, &ielem0, coop, 2,
                                 nullptr, -1, nullptr, coopr0, bary0, tol);
      int ierr1 = locMesh<2,2,1>(msh, &ielem1, coop, 2,
                                 nullptr, -1, nullptr, coopr1, bary1, tol);
      if(ierr0 != LOC_ERR_NOERR || ierr1 != LOC_ERR_NOERR){
        n_indep_fail++;
        continue;
      }
      if(ielem0 != ielem1) n_indep_fail++;
    }
    BOOST_CHECK_EQUAL(n_indep_fail, 0);
  }
}


// ===========================================================================
// Test 2: 3D localization on the unit cube (iso.p1.2k, convex, P1)
//
//  a) Out-of-bounding-box rejection
//  b) Centroid of every tet located from seed 0
//     → LOC_ERR_NOERR, distance < tol, bary ≈ [1/4,1/4,1/4,1/4]
//  c) Every mesh vertex located from seed 0
//     → LOC_ERR_NOERR, distance < tol
// ===========================================================================
BOOST_AUTO_TEST_CASE(test_locMesh_3D)
{
  cargHandler arg("-in " METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k"
                  " -anamet 1 -verb 0");
  MetrisRunner run(arg.c, arg.v);
  Mesh<MFT> &msh = *((Mesh<MFT>*) run.msh_g);
  msh.param->iverb = 0;

  const double tol = 1.0e-10;
  double coopr[3], bary[4];

  // --- a) Out-of-BB rejection -----------------------------------------------
  {
    double coop_out[3] = {-2.0, 0.5, 0.5};
    int ielem = 0;
    int ierr = locMesh<3,3,1>(msh, &ielem, coop_out, 3,
                              nullptr, -1, nullptr, coopr, bary, tol);
    BOOST_CHECK_EQUAL(ierr, LOC_ERR_OUTBB);
  }

  // --- b) Locate every tet centroid from seed 0 -----------------------------
  {
    int n_loc_fail  = 0;
    int n_dist_fail = 0;
    int n_bary_fail = 0;

    for(int i = 0; i < msh.nelem; i++){
      if(isdeadent(i, msh.tet2poi)) continue;

      double cx = 0, cy = 0, cz = 0;
      for(int v = 0; v < 4; v++){
        cx += msh.coord(msh.tet2poi(i, v), 0);
        cy += msh.coord(msh.tet2poi(i, v), 1);
        cz += msh.coord(msh.tet2poi(i, v), 2);
      }
      cx /= 4.0; cy /= 4.0; cz /= 4.0;
      double coop[3] = {cx, cy, cz};

      int ielem = 0;
      int ierr = locMesh<3,3,1>(msh, &ielem, coop, 3,
                                nullptr, -1, nullptr, coopr, bary, tol);
      if(ierr != LOC_ERR_NOERR){ n_loc_fail++;  continue; }

      double dist = std::sqrt(geterrl2<3>(coop, coopr));
      if(dist > 1.0e-6){ n_dist_fail++; continue; }

      // At its centroid, bary must be [1/4, 1/4, 1/4, 1/4]
      for(int v = 0; v < 4; v++){
        if(std::abs(bary[v] - 0.25) > 1.0e-6){ n_bary_fail++; break; }
      }
    }
    BOOST_CHECK_EQUAL(n_loc_fail,  0);
    BOOST_CHECK_EQUAL(n_dist_fail, 0);
    BOOST_CHECK_EQUAL(n_bary_fail, 0);
  }

  // --- c) Locate every vertex from seed 0 -----------------------------------
  {
    int n_vtx_fail = 0;
    for(int ip = 0; ip < msh.npoin; ip++){
      if(msh.isdeadpoint(ip)) continue;
      double coop[3] = {msh.coord(ip,0), msh.coord(ip,1), msh.coord(ip,2)};
      int ielem = 0;
      int ierr = locMesh<3,3,1>(msh, &ielem, coop, 3,
                                nullptr, -1, nullptr, coopr, bary, tol);
      if(ierr != LOC_ERR_NOERR){ n_vtx_fail++; continue; }
      double dist = std::sqrt(geterrl2<3>(coop, coopr));
      if(dist > 1.0e-6) n_vtx_fail++;
    }
    BOOST_CHECK_EQUAL(n_vtx_fail, 0);
  }

  // --- d) Seed independence on convex domain --------------------------------
  int seed_last_3d = -1;
  for(int i = msh.nelem - 1; i >= 0; i--){
    if(!isdeadent(i, msh.tet2poi)){ seed_last_3d = i; break; }
  }
  BOOST_REQUIRE(seed_last_3d > 0);

  {
    int n_indep_fail = 0;
    int ntested = 0;
    for(int i = 0; i < msh.nelem && ntested < 200; i++){
      if(isdeadent(i, msh.tet2poi)) continue;
      ntested++;

      double cx = 0, cy = 0, cz = 0;
      for(int v = 0; v < 4; v++){
        cx += msh.coord(msh.tet2poi(i, v), 0);
        cy += msh.coord(msh.tet2poi(i, v), 1);
        cz += msh.coord(msh.tet2poi(i, v), 2);
      }
      cx /= 4.0; cy /= 4.0; cz /= 4.0;
      double coop[3] = {cx, cy, cz};

      double bary0[4], bary1[4], coopr0[3], coopr1[3];
      int ielem0 = 0, ielem1 = seed_last_3d;
      int ierr0 = locMesh<3,3,1>(msh, &ielem0, coop, 3,
                                 nullptr, -1, nullptr, coopr0, bary0, tol);
      int ierr1 = locMesh<3,3,1>(msh, &ielem1, coop, 3,
                                 nullptr, -1, nullptr, coopr1, bary1, tol);
      if(ierr0 != LOC_ERR_NOERR || ierr1 != LOC_ERR_NOERR){
        n_indep_fail++;
        continue;
      }
      if(ielem0 != ielem1) n_indep_fail++;
    }
    BOOST_CHECK_EQUAL(n_indep_fail, 0);
  }
}
