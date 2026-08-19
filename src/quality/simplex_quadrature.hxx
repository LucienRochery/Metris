// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef __METRIS_SIMPLEX_QUADRATURE__
#define __METRIS_SIMPLEX_QUADRATURE__

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
// interface so legacy and tabulated positive rules can share the same consumer.
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

} // namespace Metris

#endif // __METRIS_SIMPLEX_QUADRATURE__
