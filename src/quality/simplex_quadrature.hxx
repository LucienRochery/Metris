// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_SIMPLEX_QUADRATURE__
#define __METRIS_SIMPLEX_QUADRATURE__

#include "../aux_exceptions.hxx"

namespace Metris
{

// Non-owning access to one reference-simplex quadrature point. Barycentric
// storage is owned by the quadrature rule that created the view.
template <int tdim>
struct SimplexQuadraturePointView
{
  static_assert(tdim == 2 || tdim == 3);

  const double *bary;
  double weight;
};

// Lightweight, non-owning view of a triangle or tetrahedron quadrature rule.
// Barycentric coordinates are stored point-major with tdim + 1 entries per
// point. Rule construction and storage are deliberately separate from this
// interface so vertex-barycenter and tabulated positive rules can share the
// same consumer.
template <int tdim>
class SimplexQuadratureView
{
public:
  static_assert(tdim == 2 || tdim == 3);

  using PointView = SimplexQuadraturePointView<tdim>;

  constexpr SimplexQuadratureView(int nquad,
                                  const double *bary,
                                  const double *weights) noexcept
      : nquad_(nquad), bary_(bary), weights_(weights)
  {
  }

  constexpr int size() const noexcept
  {
    return nquad_;
  }

  constexpr PointView operator[](int iquad) const noexcept
  {
    return {bary_ + iquad * (tdim + 1), weights_[iquad]};
  }

private:
  int nquad_;
  const double *bary_;
  const double *weights_;
};

namespace detail
{

template <int tdim>
struct VertexBarycenterQuadratureStorage;

template <>
struct VertexBarycenterQuadratureStorage<2>
{
  static constexpr int nquad = 4;

  inline static constexpr double bary[nquad * 3] = {
      1.0,       0.0,       0.0,
      0.0,       1.0,       0.0,
      0.0,       0.0,       1.0,
      1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};

  inline static constexpr double weights[nquad] = {
      1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0};
};

template <>
struct VertexBarycenterQuadratureStorage<3>
{
  static constexpr int nquad = 5;

  inline static constexpr double bary[nquad * 4] = {
      1.0,       0.0,       0.0,       0.0,
      0.0,       1.0,       0.0,       0.0,
      0.0,       0.0,       1.0,       0.0,
      0.0,       0.0,       0.0,       1.0,
      1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0};

  inline static constexpr double weights[nquad] = {
      1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0};
};

template <int tdim, int order>
struct PositiveSimplexQuadratureStorage;

// Positive degree-2 and degree-3 rules imported from SANS
//   src/Quadrature/QuadratureArea_Triangle.cpp
//   src/Quadrature/QuadratureVolume_Tetrahedron.cpp
// SANS stores triangle points as (x, y) and tetrahedron points as (x, y, z).
// The tables below preserve its point ordering and convert each point to the
// Metris barycentric convention (1 - sum(x_i), x_1, ...). SANS normalizes the
// weights of these reference-simplex rules to sum to one.
template <>
struct PositiveSimplexQuadratureStorage<2, 2>
{
  static constexpr int nquad = 3;

  inline static constexpr double bary[nquad * 3] = {
      1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0,
      1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0,
      2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0};

  inline static constexpr double weights[nquad] = {
      1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
};

template <>
struct PositiveSimplexQuadratureStorage<2, 3>
{
  static constexpr int nquad = 6;
  inline static constexpr double a
      = 0.65902762237409221517838077125540;
  inline static constexpr double b
      = 0.23193336855303057249678456117469;
  inline static constexpr double c
      = 0.10903900907287721232483466756991;

  inline static constexpr double bary[nquad * 3] = {
      c, a, b,
      a, b, c,
      b, c, a,
      b, a, c,
      c, b, a,
      a, c, b};

  inline static constexpr double weights[nquad] = {
      1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0,
      1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};
};

template <>
struct PositiveSimplexQuadratureStorage<3, 2>
{
  static constexpr int nquad = 4;
  inline static constexpr double a
      = 0.13819660112501051517954131656344;
  inline static constexpr double b
      = 0.58541019662496845446137605030969;

  inline static constexpr double bary[nquad * 4] = {
      a, b, a, a,
      a, a, b, a,
      a, a, a, b,
      b, a, a, a};

  inline static constexpr double weights[nquad] = {
      1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0};
};

template <>
struct PositiveSimplexQuadratureStorage<3, 3>
{
  static constexpr int nquad = 8;
  inline static constexpr double a = 3.2805469671142667e-01;
  inline static constexpr double b = 1.5835909865719922e-02;
  inline static constexpr double c = 1.0695227393293068e-01;
  inline static constexpr double d = 6.7914317820120795e-01;

  inline static constexpr double bary[nquad * 4] = {
      a, a, a, b,
      a, a, b, a,
      a, b, a, a,
      b, a, a, a,
      c, c, c, d,
      c, c, d, c,
      c, d, c, c,
      d, c, c, c};

  inline static constexpr double weights[nquad] = {
      2.3087994418643690e-02 * 6.0,
      2.3087994418643690e-02 * 6.0,
      2.3087994418643690e-02 * 6.0,
      2.3087994418643690e-02 * 6.0,
      1.8578672248022978e-02 * 6.0,
      1.8578672248022978e-02 * 6.0,
      1.8578672248022978e-02 * 6.0,
      1.8578672248022978e-02 * 6.0};
};

} // namespace detail

// Historical quality integration rule: all reference-simplex vertices followed
// by the barycenter, with equal normalized weights. Static rule storage keeps
// every returned non-owning view valid for the lifetime of the program.
template <int tdim>
constexpr SimplexQuadratureView<tdim>
get_vertex_barycenter_quadrature() noexcept
{
  return {detail::VertexBarycenterQuadratureStorage<tdim>::nquad,
          detail::VertexBarycenterQuadratureStorage<tdim>::bary,
          detail::VertexBarycenterQuadratureStorage<tdim>::weights};
}

// Compile-time access to the positive rules.
template <int tdim, int order>
constexpr SimplexQuadratureView<tdim>
get_positive_simplex_quadrature() noexcept
{
  static_assert(order == 2 || order == 3);
  return {detail::PositiveSimplexQuadratureStorage<tdim, order>::nquad,
          detail::PositiveSimplexQuadratureStorage<tdim, order>::bary,
          detail::PositiveSimplexQuadratureStorage<tdim, order>::weights};
}

// Runtime selection entry point shared by every objective-driven consumer.
// Unsupported requests must never fall back silently.
template <int tdim>
SimplexQuadratureView<tdim> get_objective_quadrature(int order)
{
  METRIS_ENFORCE_MSG(
      order == 0 || order == 2 || order == 3,
      "Unsupported objective quadrature order {}: available orders are "
      "0, 2, and 3",
      order);
  if (order == 0)
  {
    return get_vertex_barycenter_quadrature<tdim>();
  }
  if (order == 2)
  {
    return get_positive_simplex_quadrature<tdim, 2>();
  }
  return get_positive_simplex_quadrature<tdim, 3>();
}

} // namespace Metris

#endif // __METRIS_SIMPLEX_QUADRATURE__
