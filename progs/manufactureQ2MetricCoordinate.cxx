// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2026, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#include <libmeshb7/libmeshb7.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

struct Vertex
{
  double x;
  double y;
  double xi;
  double eta;
};

struct Triangle
{
  std::array<int,3> vertex;
  std::array<int,3> midpoint;
};

struct BoundaryEdge
{
  std::array<int,2> vertex;
  int midpoint;
  int reference;
};

using EdgeKey = std::pair<int,int>;

void checkCall(int result, const std::string& operation)
{
  if (result == 0) throw std::runtime_error(operation);
}

EdgeKey edgeKey(int vertex0, int vertex1)
{
  return std::minmax(vertex0,vertex1);
}

Vertex physicalPoint(double xi, double eta, double a)
{
  // Psi maps the global metric-coordinate square to the physical square.
  // Its inverse first coordinate is the manufactured scalar function u.
  return {(1.0 + a)*xi - a*xi*xi,eta,xi,eta};
}

std::array<int,2> selectGridDimensions(int targetVertexCount)
{
  // The compatible metric uses one square computational grid.  Equal numbers
  // of intervals in both directions and one fixed diagonal make the triangles
  // equilateral after applying J0 transformation.
  const int pointCount = std::max(
      2,static_cast<int>(std::llround(std::sqrt(targetVertexCount))));
  return {{pointCount,pointCount}};
}

int addMidpoint(int vertex0, int vertex1, double a,
                std::vector<Vertex>& vertices,
                std::map<EdgeKey,int>& edgeMidpoint)
{
  const EdgeKey key = edgeKey(vertex0,vertex1);
  const auto existing = edgeMidpoint.find(key);
  if (existing != edgeMidpoint.end()) return existing->second;

  const Vertex& endpoint0 = vertices[vertex0];
  const Vertex& endpoint1 = vertices[vertex1];
  const double xi = 0.5*(endpoint0.xi + endpoint1.xi);
  const double eta = 0.5*(endpoint0.eta + endpoint1.eta);
  const int midpoint = static_cast<int>(vertices.size());
  vertices.push_back(physicalPoint(xi,eta,a));
  edgeMidpoint.emplace(key,midpoint);
  return midpoint;
}

void addTriangle(std::array<int,3> triangleVertices,
                 std::vector<Triangle>& triangles)
{
  triangles.push_back({triangleVertices,{{-1,-1,-1}}});
}

void manufacture(int targetVertexCount,
                 const std::filesystem::path& outputPath,
                 double a,
                 int geometryDegree,
                 const std::string& diagonal)
{
  if (targetVertexCount < 4)
    throw std::runtime_error("Target vertex count must be at least four");
  if (!(std::abs(a) < 1.0))
    throw std::runtime_error("Parameter a must satisfy |a| < 1");
  if (geometryDegree != 1 && geometryDegree != 2)
    throw std::runtime_error("Geometry degree must be one or two");
  if (diagonal != "descending" && diagonal != "ascending")
    throw std::runtime_error("Diagonal must be descending or ascending");
  if (std::filesystem::exists(outputPath))
    throw std::runtime_error("Refusing to overwrite: " + outputPath.string());

  const std::array<int,2> grid = selectGridDimensions(targetVertexCount);
  const int nxPoint = grid[0];
  const int nyPoint = grid[1];
  const int nxCell = nxPoint - 1;
  const int nyCell = nyPoint - 1;
  const int topologicalVertexCount = nxPoint*nyPoint;

  const auto vertexIndex = [nxPoint](int ix, int iy)
  {
    return iy*nxPoint + ix;
  };

  std::vector<Vertex> vertices;
  vertices.reserve(topologicalVertexCount +
                   nxCell*nyPoint + nyCell*nxPoint + nxCell*nyCell);
  for (int iy = 0; iy < nyPoint; ++iy)
  {
    const double eta = static_cast<double>(iy)/nyCell;
    for (int ix = 0; ix < nxPoint; ++ix)
    {
      const double xi = static_cast<double>(ix)/nxCell;
      vertices.push_back(physicalPoint(xi,eta,a));
    }
  }

  std::vector<Triangle> triangles;
  triangles.reserve(2*nxCell*nyCell);
  for (int iy = 0; iy < nyCell; ++iy)
  {
    for (int ix = 0; ix < nxCell; ++ix)
    {
      const int lowerLeft = vertexIndex(ix,iy);
      const int lowerRight = vertexIndex(ix + 1,iy);
      const int upperLeft = vertexIndex(ix,iy + 1);
      const int upperRight = vertexIndex(ix + 1,iy + 1);

      // By default use the lower-right-to-upper-left diagonal.  Under the
      // computational metric C = J0^T J0, both triangles are exact copies of
      // the equilateral ideal triangle.  Alternating diagonals would
      // destroy that property in every other cell.
      if (diagonal == "descending")
      {
        addTriangle({{lowerLeft,lowerRight,upperLeft}},triangles);
        addTriangle({{lowerRight,upperRight,upperLeft}},triangles);
      }
      else
      {
        // Opposite diagonal: bottom-left to top-right, both triangles CCW.
        addTriangle({{lowerLeft,lowerRight,upperRight}},triangles);
        addTriangle({{lowerLeft,upperRight,upperLeft}},triangles);
      }
    }
  }

  std::vector<BoundaryEdge> boundaryEdges;
  boundaryEdges.reserve(2*nxCell + 2*nyCell);
  // Match examples/2D/square/square.egads:
  //   1 = left, 2 = bottom, 3 = right, 4 = top.
  // Keeping these references consistent lets us recover the CAD edge
  // parameters from the physical boundary coordinates when --cad is given.
  for (int ix = 0; ix < nxCell; ++ix)
    boundaryEdges.push_back(
        {{{vertexIndex(ix,0),vertexIndex(ix + 1,0)}},-1,2});
  for (int iy = 0; iy < nyCell; ++iy)
    boundaryEdges.push_back(
        {{{vertexIndex(nxCell,iy),vertexIndex(nxCell,iy + 1)}},-1,3});
  for (int ix = nxCell; ix > 0; --ix)
    boundaryEdges.push_back(
        {{{vertexIndex(ix,nyCell),vertexIndex(ix - 1,nyCell)}},-1,4});
  for (int iy = nyCell; iy > 0; --iy)
    boundaryEdges.push_back(
        {{{vertexIndex(0,iy),vertexIndex(0,iy - 1)}},-1,1});

  if (geometryDegree == 2)
  {
    std::map<EdgeKey,int> edgeMidpoint;
    for (Triangle& triangle : triangles)
    {
      // Ordering: (0,1,1), (1,0,1), (1,1,0).
      triangle.midpoint[0] = addMidpoint(
          triangle.vertex[1],triangle.vertex[2],a,vertices,edgeMidpoint);
      triangle.midpoint[1] = addMidpoint(
          triangle.vertex[0],triangle.vertex[2],a,vertices,edgeMidpoint);
      triangle.midpoint[2] = addMidpoint(
          triangle.vertex[0],triangle.vertex[1],a,vertices,edgeMidpoint);
    }
    for (BoundaryEdge& edge : boundaryEdges)
    {
      const auto midpoint =
          edgeMidpoint.find(edgeKey(edge.vertex[0],edge.vertex[1]));
      if (midpoint == edgeMidpoint.end())
        throw std::runtime_error("A boundary edge is absent from all triangles");
      edge.midpoint = midpoint->second;
    }
  }

  if (!outputPath.parent_path().empty())
    std::filesystem::create_directories(outputPath.parent_path());
  constexpr int meshVersion = 2;
  constexpr int dimension = 2;
  const int64_t output = GmfOpenMesh(
      outputPath.string().c_str(),GmfWrite,meshVersion,dimension);
  if (output == 0)
    throw std::runtime_error("Could not create output mesh: " + outputPath.string());

  checkCall(GmfSetKwd(output,GmfVertices,
                      static_cast<int64_t>(vertices.size())),
            "Could not create vertices");
  for (const Vertex& vertex : vertices)
    checkCall(GmfSetLin(output,GmfVertices,vertex.x,vertex.y,0),
              "Could not write a vertex");

  if (geometryDegree == 2)
  {
    checkCall(GmfSetKwd(output,GmfTrianglesP2Ordering,6),
              "Could not create triangle ordering");
    checkCall(GmfSetLin(output,GmfTrianglesP2Ordering,2,0,0),
              "Could not write triangle ordering");
    checkCall(GmfSetLin(output,GmfTrianglesP2Ordering,0,2,0),
              "Could not write triangle ordering");
    checkCall(GmfSetLin(output,GmfTrianglesP2Ordering,0,0,2),
              "Could not write triangle ordering");
    checkCall(GmfSetLin(output,GmfTrianglesP2Ordering,0,1,1),
              "Could not write triangle ordering");
    checkCall(GmfSetLin(output,GmfTrianglesP2Ordering,1,0,1),
              "Could not write triangle ordering");
    checkCall(GmfSetLin(output,GmfTrianglesP2Ordering,1,1,0),
              "Could not write triangle ordering");

    checkCall(GmfSetKwd(output,GmfTrianglesP2,
                        static_cast<int64_t>(triangles.size())),
              "Could not create quadratic triangles");
    for (const Triangle& triangle : triangles)
      checkCall(GmfSetLin(output,GmfTrianglesP2,
                          triangle.vertex[0] + 1,triangle.vertex[1] + 1,
                          triangle.vertex[2] + 1,triangle.midpoint[0] + 1,
                          triangle.midpoint[1] + 1,triangle.midpoint[2] + 1,1),
                "Could not write a quadratic triangle");

    checkCall(GmfSetKwd(output,GmfEdgesP2Ordering,3),
              "Could not create edge ordering");
    checkCall(GmfSetLin(output,GmfEdgesP2Ordering,2),
              "Could not write edge ordering");
    checkCall(GmfSetLin(output,GmfEdgesP2Ordering,0),
              "Could not write edge ordering");
    checkCall(GmfSetLin(output,GmfEdgesP2Ordering,1),
              "Could not write edge ordering");

    checkCall(GmfSetKwd(output,GmfEdgesP2,
                        static_cast<int64_t>(boundaryEdges.size())),
              "Could not create quadratic boundary edges");
    for (const BoundaryEdge& edge : boundaryEdges)
      checkCall(GmfSetLin(output,GmfEdgesP2,
                          edge.vertex[0] + 1,edge.vertex[1] + 1,
                          edge.midpoint + 1,edge.reference),
                "Could not write a quadratic boundary edge");
  }
  else
  {
    checkCall(GmfSetKwd(output,GmfTriangles,
                        static_cast<int64_t>(triangles.size())),
              "Could not create linear triangles");
    for (const Triangle& triangle : triangles)
      checkCall(GmfSetLin(output,GmfTriangles,
                          triangle.vertex[0] + 1,triangle.vertex[1] + 1,
                          triangle.vertex[2] + 1,1),
                "Could not write a linear triangle");

    checkCall(GmfSetKwd(output,GmfEdges,
                        static_cast<int64_t>(boundaryEdges.size())),
              "Could not create linear boundary edges");
    for (const BoundaryEdge& edge : boundaryEdges)
      checkCall(GmfSetLin(output,GmfEdges,
                          edge.vertex[0] + 1,edge.vertex[1] + 1,
                          edge.reference),
                "Could not write a linear boundary edge");
  }

  checkCall(GmfSetKwd(output,GmfCorners,4),"Could not create corners");
  checkCall(GmfSetLin(output,GmfCorners,vertexIndex(0,0) + 1),
            "Could not write a corner");
  checkCall(GmfSetLin(output,GmfCorners,vertexIndex(nxCell,0) + 1),
            "Could not write a corner");
  checkCall(GmfSetLin(output,GmfCorners,vertexIndex(nxCell,nyCell) + 1),
            "Could not write a corner");
  checkCall(GmfSetLin(output,GmfCorners,vertexIndex(0,nyCell) + 1),
            "Could not write a corner");
  checkCall(GmfCloseMesh(output),"Could not close output mesh");

  std::cout << "Manufactured " << outputPath << '\n'
            << "  requested topological vertices: " << targetVertexCount << '\n'
            << "  grid points:                    " << nxPoint << " x "
            << nyPoint << '\n'
            << "  actual topological vertices:    " << topologicalVertexCount << '\n'
            << "  geometry points:                " << vertices.size() << '\n'
            << "  triangles:                      " << triangles.size() << '\n'
            << "  geometry degree:                " << geometryDegree << '\n'
            << "  computational diagonal:         " << diagonal << '\n'
            << "  parameter a:                    " << a << '\n';
}

}

int main(int argc, char** argv)
{
  if (argc < 3 || argc > 6)
  {
    std::cerr << "Usage: " << argv[0]
              << " TARGET_TOPOLOGICAL_VERTICES OUTPUT.meshb [a] "
              << "[geometry_degree] [descending|ascending]\n";
    return 2;
  }

  try
  {
    const int targetVertexCount = std::stoi(argv[1]);
    const double a = argc >= 4 ? std::stod(argv[3]) : 0.5;
    const int geometryDegree = argc >= 5 ? std::stoi(argv[4]) : 2;
    const std::string diagonal = argc >= 6 ? argv[5] : "descending";
    manufacture(targetVertexCount,argv[2],a,geometryDegree,diagonal);
  }
  catch (const std::exception& error)
  {
    std::cerr << "manufactureQ2MetricCoordinate: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
