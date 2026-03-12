//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2025, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

// Tests for:
//  - Low-level geometric predicates (isintetP1, isinfacP1, getmeasentP1)
//  - Neighbour consistency after iniMeshNeighbours on a structured 2D mesh
//  - Cavity point insertion in 2D (interior face, structured mesh)
//  - Cavity point insertion in 3D (interior tet, file-loaded mesh)
//  - Intrinsic metric computation on the structured 2D mesh

#define BOOST_TEST_MODULE test_cavity

#include "common_setup.hxx"

#include "Metris.h"

using namespace Metris;

typedef MetricFieldFE MFT;

// ---------------------------------------------------------------------------
// Helper: build a structured N×N triangle mesh of the unit square [0,1]².
// Each N×N cell is split into 2 triangles (lower-left and upper-right).
// No CAD, no boundary points beyond corners.  Requires idim set first.
// ---------------------------------------------------------------------------
static void build_structured_2D(Mesh<MFT> &msh, MetrisParameters &param, int N)
{
  msh.idim   = 2;
  msh.curdeg = 1;
  msh.strdeg = 1;
  msh.forceBasisFlag(FEBasis::Lagrange);
  msh.ndomn  = 0;
  msh.param  = &param;

  const int npoin = (N+1)*(N+1);
  const int nface = 2*N*N;

  // idim must be set before set_npoin so MeshMetric allocates met.rfld with
  // correct nnmet = idim*(idim+1)/2.
  msh.set_npoin(npoin);
  msh.set_nbpoi(0);        // allocate bpo2ibi even when empty
  msh.set_nface(nface);
  msh.set_nedge(0);        // boundary edges created by iniMeshNeighbours

  // Vertices: (i,j) → ipoin = i*(N+1)+j, coords (i/N, j/N)
  for(int i = 0; i <= N; i++){
    for(int j = 0; j <= N; j++){
      int ip = i*(N+1) + j;
      msh.coord(ip, 0) = (double)i / N;
      msh.coord(ip, 1) = (double)j / N;
      msh.poi2bpo[ip]  = -1;
    }
  }
  // poi2bak must be -1 (no back mesh)
  msh.poi2bak.fill(-1);

  // Triangles: each square (i,j) → two triangles
  //   A (lower-left):  v00, v10, v01
  //   B (upper-right): v11, v01, v10
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      int sq  = i*N + j;
      int v00 = i    *(N+1) + j;
      int v10 = (i+1)*(N+1) + j;
      int v01 = i    *(N+1) + j+1;
      int v11 = (i+1)*(N+1) + j+1;

      msh.fac2poi(2*sq,   0) = v00;
      msh.fac2poi(2*sq,   1) = v10;
      msh.fac2poi(2*sq,   2) = v01;
      msh.fac2ref[2*sq]      = 1;

      msh.fac2poi(2*sq+1, 0) = v11;
      msh.fac2poi(2*sq+1, 1) = v01;
      msh.fac2poi(2*sq+1, 2) = v10;
      msh.fac2ref[2*sq+1]    = 1;
    }
  }

  // Zero tag arrays, edg2fac, poi2ent_ etc.
  msh.zeroArrays();

  // Build fac2fac, boundary edges and corners
  iniMeshNeighbours<1>(msh);

  // Set poi2ent_ so isdeadpoint() returns false.
  // Faces first (dim 2), then boundary edges override with dim 1.
  for(int i = 0; i < msh.nface; i++){
    if(isdeadent(i, msh.fac2poi)) continue;
    for(int v = 0; v < 3; v++)
      msh.set_poi2ent(Vertex{msh.fac2poi(i,v)}, 2, i);
  }
  for(int ie = 0; ie < msh.nedge; ie++){
    if(isdeadent(ie, msh.edg2poi)) continue;
    for(int v = 0; v < 2; v++)
      msh.set_poi2ent(Vertex{msh.edg2poi(ie,v)}, 1, ie);
  }
}


// ===========================================================================
// Test 1: Low-level geometric predicates — no Mesh object required
// ===========================================================================
BOOST_AUTO_TEST_CASE(test_geo_predicates)
{
  // --- isintetP1: reference tet {(0,0,0),(1,0,0),(0,1,0),(0,0,1)}
  {
    double p0[3] = {0,0,0}, p1[3] = {1,0,0}, p2[3] = {0,1,0}, p3[3] = {0,0,1};

    double pin [3] = {0.25, 0.25, 0.25}; // centroid — inside
    double pout[3] = {1.0,  1.0,  1.0};  // outside
    double pfac[3] = {0.5,  0.5,  0.0};  // on face p0-p1-p2 — boundary (inside with default tol)
    double pcor[3] = {1.0,  0.0,  0.0};  // vertex p1 — exactly on tet

    BOOST_CHECK( isintetP1(p0,p1,p2,p3, pin));
    BOOST_CHECK(!isintetP1(p0,p1,p2,p3, pout));
    BOOST_CHECK( isintetP1(p0,p1,p2,p3, pfac));
    BOOST_CHECK( isintetP1(p0,p1,p2,p3, pcor));
  }

  // --- isinfacP1: right triangle in 3D plane z=0
  {
    double q0[3] = {0,0,0}, q1[3] = {1,0,0}, q2[3] = {0,1,0};

    double qin [3] = {0.25, 0.25, 0.0}; // inside
    double qout[3] = {1.0,  1.0,  0.0}; // outside (sum of bary coords > 1)
    double qedg[3] = {0.5,  0.0,  0.0}; // midpoint of q0-q1 — on edge

    BOOST_CHECK( isinfacP1(q0,q1,q2, qin));
    BOOST_CHECK(!isinfacP1(q0,q1,q2, qout));
    BOOST_CHECK( isinfacP1(q0,q1,q2, qedg));
  }

  // --- getmeasentP1<2>: area of right triangle (0,0)-(1,0)-(0,1) = 0.5
  {
    dblAr2 coord(3, 2);
    coord(0,0) = 0; coord(0,1) = 0;
    coord(1,0) = 1; coord(1,1) = 0;
    coord(2,0) = 0; coord(2,1) = 1;
    int ep[3] = {0, 1, 2};
    double area = getmeasentP1<2>(ep, coord);
    BOOST_CHECK_CLOSE(area, 0.5, 1e-10);
  }

  // --- getmeasentP1<2>: negative for flipped winding
  {
    dblAr2 coord(3, 2);
    coord(0,0) = 0; coord(0,1) = 0;
    coord(1,0) = 0; coord(1,1) = 1; // swapped vs previous
    coord(2,0) = 1; coord(2,1) = 0;
    int ep[3] = {0, 1, 2};
    BOOST_CHECK(getmeasentP1<2>(ep, coord) < 0.0);
  }

  // --- getmeasentP1<3>: volume of reference tet = 1/6
  {
    dblAr2 coord(4, 3);
    coord(0,0) = 0; coord(0,1) = 0; coord(0,2) = 0;
    coord(1,0) = 1; coord(1,1) = 0; coord(1,2) = 0;
    coord(2,0) = 0; coord(2,1) = 1; coord(2,2) = 0;
    coord(3,0) = 0; coord(3,1) = 0; coord(3,2) = 1;
    int ep[4] = {0, 1, 2, 3};
    double vol = getmeasentP1<3>(ep, coord);
    BOOST_CHECK_CLOSE(vol, 1.0/6.0, 1e-10);
  }

  // --- getmeasentP1<3>: negative for flipped tet
  {
    dblAr2 coord(4, 3);
    coord(0,0) = 0; coord(0,1) = 0; coord(0,2) = 0;
    coord(1,0) = 0; coord(1,1) = 1; coord(1,2) = 0; // swapped p1/p2
    coord(2,0) = 1; coord(2,1) = 0; coord(2,2) = 0;
    coord(3,0) = 0; coord(3,1) = 0; coord(3,2) = 1;
    int ep[4] = {0, 1, 2, 3};
    BOOST_CHECK(getmeasentP1<3>(ep, coord) < 0.0);
  }
}


// ===========================================================================
// Test 2: fac2fac symmetry and shared-vertex consistency after
//         iniMeshNeighbours on a 5×5 structured 2D mesh
// ===========================================================================
BOOST_AUTO_TEST_CASE(test_neighbours_2D)
{
  MetrisParameters param;
  param.iverb = 0;

  Mesh<MFT> msh;
  const int N = 5; // 5×5 grid → 50 triangles, 36 vertices
  build_structured_2D(msh, param, N);

  const int nface = msh.nface;

  // --- fac2fac symmetry: if fac2fac(i,e)==j then ∃f: fac2fac(j,f)==i
  int symm_fail = 0;
  for(int i = 0; i < nface; i++){
    for(int e = 0; e < 3; e++){
      int j = msh.fac2fac(i, e);
      if(j < 0) continue;
      bool found = false;
      for(int f = 0; f < 3; f++)
        if(msh.fac2fac(j, f) == i){ found = true; break; }
      if(!found) symm_fail++;
    }
  }
  BOOST_CHECK_EQUAL(symm_fail, 0);

  // --- shared vertices: if fac2fac(i,e)==j, the two vertices of i not at
  //     local index e must both appear in j's vertex list
  int vtx_fail = 0;
  for(int i = 0; i < nface; i++){
    for(int e = 0; e < 3; e++){
      int j = msh.fac2fac(i, e);
      if(j < 0) continue;
      int iv0 = msh.fac2poi(i, (e+1)%3);
      int iv1 = msh.fac2poi(i, (e+2)%3);
      bool has0 = false, has1 = false;
      for(int v = 0; v < 3; v++){
        if(msh.fac2poi(j, v) == iv0) has0 = true;
        if(msh.fac2poi(j, v) == iv1) has1 = true;
      }
      if(!has0 || !has1) vtx_fail++;
    }
  }
  BOOST_CHECK_EQUAL(vtx_fail, 0);

  // --- boundary edge count: N×N grid has 4×N boundary edges, each
  //     corresponding to one (face, local-edge) pair with fac2fac == -1
  int nbdry = 0;
  for(int i = 0; i < nface; i++)
    for(int e = 0; e < 3; e++)
      if(msh.fac2fac(i, e) == -1) nbdry++;
  BOOST_CHECK_EQUAL(nbdry, 4*N);

  // --- all boundary-edge (face,edge) pairs have a corresponding global edge
  int missing_edge = 0;
  for(int i = 0; i < nface; i++){
    for(int e = 0; e < 3; e++){
      if(msh.fac2fac(i, e) != -1) continue;
      int ip1 = msh.fac2poi(i, lnoed2[e][0]);
      int ip2 = msh.fac2poi(i, lnoed2[e][1]);
      if(getedgglo(msh, ip1, ip2) < 0) missing_edge++;
    }
  }
  BOOST_CHECK_EQUAL(missing_edge, 0);
}


// ===========================================================================
// Test 3: cavity point insertion in 2D — insert a vertex at the centroid
//         of a fully interior face of the structured mesh.
//         Expected: one face removed, three new faces created, all valid.
// ===========================================================================
BOOST_AUTO_TEST_CASE(test_cavity_insert_2D)
{
  MetrisParameters param;
  param.iverb = 0;

  Mesh<MFT> msh;
  const int N = 5;
  build_structured_2D(msh, param, N);

  // Find a fully interior face (all three neighbours present)
  int iface_int = -1;
  for(int i = 0; i < msh.nface; i++){
    if(msh.fac2fac(i,0) >= 0 && msh.fac2fac(i,1) >= 0 && msh.fac2fac(i,2) >= 0){
      iface_int = i;
      break;
    }
  }
  BOOST_REQUIRE_MESSAGE(iface_int >= 0, "No interior face found");

  // Centroid of the chosen face
  double cx = 0, cy = 0;
  for(int v = 0; v < 3; v++){
    cx += msh.coord(msh.fac2poi(iface_int, v), 0);
    cy += msh.coord(msh.fac2poi(iface_int, v), 1);
  }
  cx /= 3.0; cy /= 3.0;

  // Add new interior point
  const int ipnew = msh.npoin;
  msh.set_npoin(ipnew + 1);
  msh.coord(ipnew, 0) = cx;
  msh.coord(ipnew, 1) = cy;
  msh.poi2bpo[ipnew]  = -1;
  msh.poi2bak[ipnew]  = -1;
  msh.set_poi2ent(Vertex{ipnew}, 0, -1); // unconnected

  const int nfac_before = msh.nface;

  MshCavity  cav(0, 10, 0);
  cav.ipins = ipnew;
  cav.inewp = 0; // brand-new point
  cav.lcfac.stack(iface_int);

  CavOprOpt opts;
  opts.allow_topological_correction = false;
  opts.skip_topo_checks             = false;
  opts.qmax_nec = -1; opts.qmax_suf = -1; opts.qmax_iff = -1;

  CavWrkArrs work;
  CavOprInfo info;

  int ierr = cavity_operator<MFT, 1>(msh, cav, opts, work, info, 0);

  BOOST_CHECK_EQUAL(ierr, CAV_NOERR);
  BOOST_CHECK(info.done);

  // 1 face marked dead, 3 new faces appended → nface += 3
  BOOST_CHECK_EQUAL(msh.nface, nfac_before + 3);

  // All three new faces must contain ipnew and have positive area
  int n_new = 0;
  for(int i = nfac_before; i < msh.nface; i++){
    if(isdeadent(i, msh.fac2poi)) continue;
    bool has_new = false;
    for(int v = 0; v < 3; v++)
      if(msh.fac2poi(i, v) == ipnew){ has_new = true; break; }
    if(!has_new) continue;
    n_new++;
    double area = getmeasentP1<2>(msh.fac2poi[i], msh.coord);
    BOOST_CHECK_GT(area, 0.0);
  }
  BOOST_CHECK_EQUAL(n_new, 3);
}


// ===========================================================================
// Test 4: cavity point insertion in 3D — load a simple tetrahedral mesh,
//         insert a vertex at the centroid of the first fully interior tet.
//         Expected: one tet removed, four new tets created, all valid.
// ===========================================================================
BOOST_AUTO_TEST_CASE(test_cavity_insert_3D)
{
  typedef MetricFieldAnalytical MFT3;

  cargHandler arg("-in " METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k"
                  " -anamet 1 -verb 0");
  MetrisRunner run(arg.c, arg.v);
  Mesh<MFT3> &msh = *static_cast<Mesh<MFT3>*>(run.msh_g);

  // Find a fully interior tet (all four face-neighbours >= 0)
  int itet_int = -1;
  for(int i = 0; i < msh.nelem; i++){
    if(isdeadent(i, msh.tet2poi)) continue;
    if(msh.tet2tet(i,0) >= 0 && msh.tet2tet(i,1) >= 0 &&
       msh.tet2tet(i,2) >= 0 && msh.tet2tet(i,3) >= 0){
      itet_int = i;
      break;
    }
  }
  BOOST_REQUIRE_MESSAGE(itet_int >= 0, "No interior tet found");

  // Centroid
  double cx=0, cy=0, cz=0;
  for(int v = 0; v < 4; v++){
    int ip = msh.tet2poi(itet_int, v);
    cx += msh.coord(ip, 0);
    cy += msh.coord(ip, 1);
    cz += msh.coord(ip, 2);
  }
  cx /= 4.0; cy /= 4.0; cz /= 4.0;

  // Add new interior point
  const int ipnew = msh.npoin;
  msh.set_npoin(ipnew + 1);
  msh.coord(ipnew, 0) = cx;
  msh.coord(ipnew, 1) = cy;
  msh.coord(ipnew, 2) = cz;
  msh.poi2bpo[ipnew]  = -1;
  msh.poi2bak[ipnew]  = -1;
  msh.set_poi2ent(Vertex{ipnew}, 0, -1);

  // Metric at the new point: interpolate from the tet (needed even with
  // quality disabled to avoid uninitialised reads inside reconnect_tetcav)
  {
    const int nnmet = 6; // 3D
    for(int j = 0; j < nnmet; j++)
      msh.met(ipnew, j) = 0.25*(msh.met(msh.tet2poi(itet_int,0), j)
                               + msh.met(msh.tet2poi(itet_int,1), j)
                               + msh.met(msh.tet2poi(itet_int,2), j)
                               + msh.met(msh.tet2poi(itet_int,3), j));
  }

  const int ntet_before = msh.nelem;

  MshCavity  cav(10, 0, 0);
  cav.ipins = ipnew;
  cav.inewp = 0;
  cav.lctet.stack(itet_int);

  CavOprOpt opts;
  opts.allow_topological_correction = false;
  opts.skip_topo_checks             = false;
  opts.qmax_nec = -1; opts.qmax_suf = -1; opts.qmax_iff = -1;

  CavWrkArrs work;
  CavOprInfo info;

  int ierr = cavity_operator<MFT3, 1>(msh, cav, opts, work, info, 0);

  BOOST_CHECK_EQUAL(ierr, CAV_NOERR);
  BOOST_CHECK(info.done);

  // 1 tet marked dead, 4 new tets appended → nelem += 4
  BOOST_CHECK_EQUAL(msh.nelem, ntet_before + 4);

  // All four new tets must contain ipnew and have positive volume
  int n_new = 0;
  for(int i = ntet_before; i < msh.nelem; i++){
    if(isdeadent(i, msh.tet2poi)) continue;
    bool has_new = false;
    for(int v = 0; v < 4; v++)
      if(msh.tet2poi(i, v) == ipnew){ has_new = true; break; }
    if(!has_new) continue;
    n_new++;
    double vol = getmeasentP1<3>(msh.tet2poi[i], msh.coord);
    BOOST_CHECK_GT(vol, 0.0);
  }
  BOOST_CHECK_EQUAL(n_new, 4);
}


// ===========================================================================
// Test 5: intrinsic metric on the structured 2D mesh
//   - getMetMesh completes without error
//   - At every non-dead vertex the metric trace (m11+m22 in log space) is
//     positive, meaning the geometric mean eigenvalue of M is > 1 (as
//     expected for a mesh with h ≈ 1/N < 1, so M eigenvalues ≈ N²).
//   - The distribution of traces is reasonably uniform for a uniform mesh
//     (max/min ratio < 5).
// ===========================================================================
BOOST_AUTO_TEST_CASE(test_intrinsic_metric_2D)
{
  MetrisParameters param;
  param.iverb = 0;
  param.hmin  = 1.0e-30;
  param.hmax  = 1.0e30;

  Mesh<MFT> msh;
  const int N = 5;
  build_structured_2D(msh, param, N);

  getMetMesh<MFT, 1>(param, msh);

  // In log space: met(ip,0)=m11_log, met(ip,1)=m12_log, met(ip,2)=m22_log
  // trace_log = m11_log + m22_log  (= log(det M) for an isotropic metric)
  double tr_min =  1e30;
  double tr_max = -1e30;
  int n_finite_fail = 0;

  for(int ip = 0; ip < msh.npoin; ip++){
    if(msh.isdeadpoint(ip)) continue;
    double m11 = msh.met(ip, 0);
    double m22 = msh.met(ip, 2);
    if(!std::isfinite(m11) || !std::isfinite(m22)){ n_finite_fail++; continue; }
    double tr = m11 + m22;
    tr_min = std::min(tr_min, tr);
    tr_max = std::max(tr_max, tr);
  }

  BOOST_CHECK_EQUAL(n_finite_fail, 0);
  // For N=5 (h=0.2), eigenvalues ≈ 25, log ≈ 3.22 each → trace ≈ 6.44
  BOOST_CHECK_GT(tr_min, 0.0);
  // Uniform mesh: trace spread should be moderate (allow factor 5 for corners)
  BOOST_CHECK_LT(tr_max / tr_min, 5.0);
}


// ===========================================================================
// Test 6: intrinsic metric on a 3D tetrahedral mesh.
//   - getMetMesh completes without error (3D code path)
//   - At every live vertex the 3D metric trace m11+m22+m33 (in log space,
//     = log(det M)) is positive: for h < 1, eigenvalues > 1 → log > 0.
//   - Trace spread is moderate for a roughly uniform mesh.
// ===========================================================================
BOOST_AUTO_TEST_CASE(test_intrinsic_metric_3D)
{
  typedef MetricFieldFE MFT3;

  cargHandler arg("-in " METRIS_CASES_DIR "/unit/3D/cube/iso.p1.2k"
                  " -verb 0");
  MetrisRunner run(arg.c, arg.v);
  Mesh<MFT3> &msh = *static_cast<Mesh<MFT3>*>(run.msh_g);

  MetrisParameters &param = *msh.param;
  param.hmin = 1.0e-30;
  param.hmax = 1.0e30;

  getMetMesh<MFT3, 1>(param, msh);

  // For 3D, sym2idx(i,i): m11→0, m22→2, m33→5
  // trace_log = met(ip,0) + met(ip,2) + met(ip,5) = log(det M)
  double tr_min =  1e30;
  double tr_max = -1e30;
  int n_finite_fail = 0;

  for(int ip = 0; ip < msh.npoin; ip++){
    if(msh.isdeadpoint(ip)) continue;
    double m11 = msh.met(ip, 0);
    double m22 = msh.met(ip, 2);
    double m33 = msh.met(ip, 5);
    if(!std::isfinite(m11) || !std::isfinite(m22) || !std::isfinite(m33)){
      n_finite_fail++;
      continue;
    }
    double tr = m11 + m22 + m33;
    tr_min = std::min(tr_min, tr);
    tr_max = std::max(tr_max, tr);
  }

  BOOST_CHECK_EQUAL(n_finite_fail, 0);
  // h ≈ 1/7 → eigenvalues ≈ 49, log(det M) = 3*log(49) ≈ 11.5
  BOOST_CHECK_GT(tr_min, 0.0);
  // Allow factor 8 spread for corners of a tetrahedral mesh
  BOOST_CHECK_LT(tr_max / tr_min, 8.0);
}


// ===========================================================================
// Test 7: getmeasentP1 degenerate cases — zero area and zero volume.
//   Collinear triangle → area exactly 0.
//   Coplanar tet       → volume exactly 0.
// ===========================================================================
BOOST_AUTO_TEST_CASE(test_measure_degenerate)
{
  // --- 2D: three collinear points → area = 0
  {
    dblAr2 coord(3, 2);
    coord(0,0) = 0.0; coord(0,1) = 0.0;
    coord(1,0) = 0.5; coord(1,1) = 0.5; // midpoint of 0→2
    coord(2,0) = 1.0; coord(2,1) = 1.0;
    int ep[3] = {0, 1, 2};
    BOOST_CHECK_SMALL(getmeasentP1<2>(ep, coord), 1.0e-15);
  }

  // --- 2D: two identical points → area = 0
  {
    dblAr2 coord(3, 2);
    coord(0,0) = 0.0; coord(0,1) = 0.0;
    coord(1,0) = 1.0; coord(1,1) = 0.0;
    coord(2,0) = 0.0; coord(2,1) = 0.0; // duplicate of vertex 0
    int ep[3] = {0, 1, 2};
    BOOST_CHECK_SMALL(getmeasentP1<2>(ep, coord), 1.0e-15);
  }

  // --- 3D: four coplanar points → volume = 0
  // All four in the z=0 plane
  {
    dblAr2 coord(4, 3);
    coord(0,0) = 0.0; coord(0,1) = 0.0; coord(0,2) = 0.0;
    coord(1,0) = 1.0; coord(1,1) = 0.0; coord(1,2) = 0.0;
    coord(2,0) = 0.0; coord(2,1) = 1.0; coord(2,2) = 0.0;
    coord(3,0) = 0.5; coord(3,1) = 0.5; coord(3,2) = 0.0; // in plane
    int ep[4] = {0, 1, 2, 3};
    BOOST_CHECK_SMALL(getmeasentP1<3>(ep, coord), 1.0e-15);
  }

  // --- 3D: two identical vertices → volume = 0
  {
    dblAr2 coord(4, 3);
    coord(0,0) = 0.0; coord(0,1) = 0.0; coord(0,2) = 0.0;
    coord(1,0) = 1.0; coord(1,1) = 0.0; coord(1,2) = 0.0;
    coord(2,0) = 0.0; coord(2,1) = 1.0; coord(2,2) = 0.0;
    coord(3,0) = 0.0; coord(3,1) = 0.0; coord(3,2) = 0.0; // duplicate of vertex 0
    int ep[4] = {0, 1, 2, 3};
    BOOST_CHECK_SMALL(getmeasentP1<3>(ep, coord), 1.0e-15);
  }
}
