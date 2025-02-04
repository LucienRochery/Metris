//Metris: high-order metric-based non-manifold tetrahedral remesher
//Copyright (C) 2023-2024, Massachusetts Institute of Technology
//Licensed under The GNU Lesser General Public License, version 2.1
//See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php



#ifndef __METRIS_CONSTANTS__
#define __METRIS_CONSTANTS__


#include <array>
#include <boost/hana/map.hpp>
#include <boost/hana/pair.hpp>
#include "../SANS/LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"
#include "metris_defaults.hxx"
namespace hana = boost::hana;

// Maximum while loop counter for non-manifold structures and other linked lists.
#ifndef METRIS_MAX_WHILE
#define METRIS_MAX_WHILE 1000
#endif

// How many tag arrays per entity. 
// The idea is that a function that uses tags in its body should use layer i+1 with i being the highest layer of its callees. 
// Thus a function calling no other function that uses tags should use layer 0. 
// The cavity operator needs layer 1.
// We need additional tags for multi-threading. 
#ifndef METRIS_MAXTAGS
#define METRIS_MAXTAGS 8
#endif


#ifndef METRIS_MAX_DEG
#define METRIS_MAX_DEG 3
#endif



#ifndef METRIS_MAX_DEG_JACOBIAN
#define METRIS_MAX_DEG_JACOBIAN (3*(METRIS_MAX_DEG - 1))
#endif


#ifdef USE_MULTIPRECISION
  #include <boost/multiprecision/cpp_bin_float.hpp>
  typedef  boost::multiprecision::cpp_bin_float_oct  float8;
  typedef  boost::multiprecision::cpp_bin_float_quad float4;

  #define QUA_FTYPE_SEQ (double)(float4)(float8)
#else
  #define QUA_FTYPE_SEQ (double)
#endif


namespace Metris{

// Evaluation related 
enum class FEBasis {Undefined=-1, Lagrange, Bezier};

enum class DifVar {None=-1, Bary, Phys};

enum class DifOrd {None=0 , Grad=1, Hess=2};

enum class MetSpace {Undefined=-1, Log, Exp};

enum class AsDeg {P1, Pk};

// API Related
enum class APIDefaultOrderings{ExodusII,Metris};


constexpr int lnoed2[3][2] = {{1,2},
                              {2,0},
                              {0,1}};
constexpr int lnoed3[6][2] = {{0,1},
                              {1,2},
                              {2,0},
                              {0,3},
                              {1,3},
                              {2,3}};
constexpr int lnofa3[4][3] = {{1,3,2},
                              {2,3,0},
                              {3,1,0},
                              {0,1,2}};
constexpr int ledfa3[4][3] = {{5,1,4},
                              {3,2,5},
                              {0,3,4},
                              {1,2,0}};



namespace Constants{

  // We'll fit the number of edges to an isotropic mesh of a square case (slight overestimation wr.t. sphere, good) with:
  // 5.4k verts
  // 10.6k triangles
  // 250 edges

  // This writes E = E_0 Nvert^(1/2) i.e. E_0 = 250 * 5400^(-1/2), T_0 = Constants::memfitCoeff12 below:
  // this is slightly under 3.5. 
  const double memfitCoeff12 = 3.5;
  // We fit the number of triangles using a case with 19k verts, 104k tet, 6400 triangles
  // Use another example? 3.271M verts, 206k triangles
  // This writes T = T_0 Nvert^(2/3) i.e. T_0 = 6400 * 19000^(-2/3), T_0 = fitCoeff23 below:
  // These are both quite close, and to 9. The largest is 9.5. Let's just set it at 10. 
  //double fitCoeff23 = 6400.0/(pow(19000.0,2.0/3.0)); 
  const double memfitCoeff23 = 10.0;
  const double memfitCoeff13 = 10500 / pow(3271000.0,1.0/3.0);


  // vK0[i] = Volume of K0 in dimension i
  const double vK0[4] = {-1, 1, 0.433012701892219298, 0.1178511301977579};


  // Tolerance for barycentric coordinates. 
  // Barycentric coordinates considered valid if in range: 
  // -baryTol < bary[ii] < 1 + baryTol
  // And sum between 1 - (idim+1) baryTol and 1 + (idim+1) baryTol
  const double baryTol = 1.0e-9;

  // Norm below which we refuse to divide by the square root. 
  const double vecNrmTol = 1.0e-16;

  // Norm of distance gradient before stop
  const double projedgTol = 1.0e-12;

  // Threshold below which dot product between two vectors is indicative
  // of them not being aligned. Introduced for CAD normals so some 
  // slack is given. 
  const double dtprdMisAlign = 0.4;

  const double detTol = 1.0e-32;

  // Threshold for t or (u,v) coordinates to be considered different at a scale of 1
  const double CADparamTol = 1.0e-10;


  // Expected number of entities per vertex (not ctrl pt)
  // See below for "proofs"
  // Note edges, triangles are not geometric, rather interior
  // Geometric entities are not proportional to number of vertices ! 
  // -> ^2/3 triangles, ^1/3 edges in 3D, 
  // -> ^1/2 edges in 2D !
  // See Mesh::allocate(). 
  // 3D
  const double tetpver3 = 6.0;
  const double facpver3 = 2.0*tetpver3;
  const double edgpver3 = 6.5;
  // 2D 
  const double facpver2 = 2.0;
  const double edgpver2 = 3.0;

  // Vertices per point 3D
  const std::array<double,METRIS_MAX_DEG_JACOBIAN+1> 
  verppoi3{[](){
    std::array<double,METRIS_MAX_DEG_JACOBIAN+1> ret{};
    for(int ideg = 0; ideg <= METRIS_MAX_DEG_JACOBIAN; ideg++){
      const int c1=(ideg-1)*(ideg>1), // pure edge 
                c2=(ideg-1)*(ideg-2)*(ideg>2), // pure triangle 
                c3=(ideg-1)*(ideg-2)*(ideg-3)*(ideg>3); // pure tetra 
      ret[ideg] = 1.0 / (1.0 + edgpver3*c1 + facpver3*c2 + tetpver3*c3); 
    }
    return ret;
  }()};
  // Vertices per point 2D
  const std::array<double,METRIS_MAX_DEG_JACOBIAN+1> 
  verppoi2{[](){
    std::array<double,METRIS_MAX_DEG_JACOBIAN+1> ret{};
    for(int ideg = 0; ideg <= METRIS_MAX_DEG_JACOBIAN; ideg++){
      const int c1=(ideg-1)*(ideg>1), // pure edge 
                c2=(ideg-1)*(ideg-2)*(ideg>2); // pure triangle 
      ret[ideg] = 1.0 / (1.0 + edgpver2*c1 + facpver2*c2); 
    }
    return ret;
  }()};


  static double r8invtJ_0[4][3*3] = {
    {0},
    {1},
    {1,0,-0.577350269189626,1.154700538379252},
    { -0.57735026918962562,-0.57735026918962562 , 0                 ,
      1                   ,-1                   , 0                 ,
      -0.40824829046386302,-0.4082482904638630  , 1.2247448713915889}
  };

  #ifdef USE_MULTIPRECISION
    static float4 r16invtJ_0[4][3*3] = {
      {0},
      {1},
      {1,0,-0.577350269189626,1.154700538379252},
      { -0.57735026918962562,-0.57735026918962562 , 0                 ,
        1                   ,-1                   , 0                 ,
        -0.40824829046386302,-0.4082482904638630  , 1.2247448713915889}
    };
    static float8 r32invtJ_0[4][3*3] = {
      {0},
      {1},
      {1,0,-0.577350269189626,1.154700538379252},
      { -0.57735026918962562,-0.57735026918962562 , 0                 ,
        1                   ,-1                   , 0                 ,
        -0.40824829046386302,-0.4082482904638630  , 1.2247448713915889}
    };
  #endif

  // This container is very simply accessed using the type in operator[]. 
  // Example invtJ_0[hana::type_c<double>] is r8invtJ_0
  static auto invtJ_0 = 
  hana::make_map(
    hana::make_pair(hana::type_c<double>, r8invtJ_0)
    #ifdef USE_MULTIPRECISION
    ,
    hana::make_pair(hana::type_c<float4>, r16invtJ_0),
    hana::make_pair(hana::type_c<float8>, r32invtJ_0)
    #endif
  );

  static SANS::DLA::MatrixS<4,9,double> r8MSinvtJ_0(r8invtJ_0[0],36);
  #ifdef USE_MULTIPRECISION
    static SANS::DLA::MatrixS<4,9,float4> r16MSinvtJ_0(r16invtJ_0[0],36);
    static SANS::DLA::MatrixS<4,9,float8> r32MSinvtJ_0(r32invtJ_0[0],36);
  #endif

  static auto invtJ_0_MS = 
  hana::make_map(
    hana::make_pair(hana::type_c<double>, r8MSinvtJ_0)
    #ifdef USE_MULTIPRECISION
    ,
    hana::make_pair(hana::type_c<float4>, r16MSinvtJ_0),
    hana::make_pair(hana::type_c<float8>, r32MSinvtJ_0)
    #endif
  );


}

};



/*
  Computation of average element counts per point.
*/

///  --- 3D 
// Consider a uniform mesh of R^3 made up of cubes in a lattice. 
// Each point is associated the single cube whose lll (0,0,0 in the case of [0,1]^3) corner it is. 
// We consider two possible cube splits and will average out quantities in the end. 

// 1.Diamond lattice, 5 tetra per cube, no interior edge, new edge per subcube face (diagonal)
//               5 tetra per interior point
//   Each vertex has 4*3 = 12 incident cube faces, half of which are split with this point as a diagonal edge extremity
//   These edges are accounted for twice, once per extremity, yielding 12/2/2 = 3 edges per point. 
//   We can also think of these as one per cube face, that is 3 unique per cube. Cubes are 1 cube = 1 point. 
//   Moreover, we have the cube edges, 6 counted twice, contributing 3 to the point. 
//   Alternatively, we have 12 edges per cube, each shared between 4 cubes, so 3 edges unique per cube. 
//   This leads to:
//               6 edges per interior point 

// 2. "Marching tetrahedra", 6 tetra per cube, 1 subcube interior edge + new edge per subcube face (diag)
//               6 tetra per interior point
//    We have, again, 3 unique edges per cube, plus 3 face diagonals, plus one interior diagonal. 
//               7 edges per interior point

// Average: 5.5 tetra per point
//          6.5 edges per point

// 0.Faces follow from a different reasoning: introduce half-faces HF, we have HF = 2F and HF = 4T thus 
//               2 faces per tetra = 11 faces per point (on average). 

// Replacing, 
// We have Nvert * 1        pure vertex control points ( = vertices ...)
//         6.5 * Nvert * (d-1)    pure edge control points 
//          11 * Nvert * (d-1)(d-2)/2 pure face control points
//         5.5 * Nvert * (d-1)(d-2)(d-3)/6 pure tet ctr points
// This means for N points in a degree d mesh, we need space for 
//  N points
//  5.5 Nvert tetrahedra 
// conservatively set to 6 


///  --- 2D 
// Half edges: 2 per edge. 3 per triangle. 
// Thus 2Nedg = 3Ntri
// Regular grid, 1 vertex -> 2 triangles. 
// Thus per point, 2 triangles, 3 edges. 


#endif