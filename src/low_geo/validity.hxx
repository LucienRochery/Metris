// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under the GNU Lesser General Public License, version 2.1

#ifndef __METRIS_LOW_GEO_VALIDITY_HXX__
#define __METRIS_LOW_GEO_VALIDITY_HXX__

#include <limits>

namespace Metris {

class MeshBase;

// A failed Bernstein lower bound is not, by itself, proof that an element is
// invalid.  Keep the proof outcome separate from the conservative acceptance
// policy used by optimization and adaptation.
enum class ElementValidityStatus {
  Certified,
  Invalid,
  Uncertified
};

constexpr const char *element_validity_status_name(
    ElementValidityStatus status) noexcept
{
  switch(status){
  case ElementValidityStatus::Certified: return "Certified";
  case ElementValidityStatus::Invalid: return "Invalid";
  case ElementValidityStatus::Uncertified: return "Uncertified";
  }
  return "Unknown";
}

struct ElementValidityResult {
  ElementValidityStatus status = ElementValidityStatus::Uncertified;

  // Minimum normalized Bernstein coefficient. For a full-dimensional element
  // this is a rigorous lower bound for the normalized Jacobian determinant
  // over the complete reference simplex.
  double normalized_lower_bound =
      std::numeric_limits<double>::quiet_NaN();

  // Zero-based local position of the minimum coefficient in the array returned
  // by getccoef: 0..5 for a P2 triangle and 0..19 for a P2 tetrahedron. This is
  // not a mesh point index. A value of -1 means that no bound was computed.
  int lower_bound_coefficient_index = -1;

  // A direct normalized Jacobian evaluation below the configured threshold is
  // a witness of invalidity.  These remain unset when no witness was evaluated
  // or when a failed coefficient bound is merely inconclusive.
  double normalized_witness = std::numeric_limits<double>::quiet_NaN();

  // Zero-based local position in the deterministic barycentric sample set used
  // by the validity classifier. This is not a mesh point index. A value of -1
  // means that no direct failing witness was found or evaluated.
  int witness_sample_index = -1;

  constexpr bool is_certified() const noexcept {
    return status == ElementValidityStatus::Certified;
  }

  constexpr bool is_invalid() const noexcept {
    return status == ElementValidityStatus::Invalid;
  }

  constexpr bool is_uncertified() const noexcept {
    return status == ElementValidityStatus::Uncertified;
  }

  // Initial Phase 3 policy: only a positive certificate is accepted.  Both an
  // invalid witness and an inconclusive bound are rejected conservatively.
  constexpr bool accepted_conservatively() const noexcept {
    return is_certified();
  }
};

// Classify one full-dimensional element against the configured normalized
// Jacobian threshold. The P1 corner orientation is a prerequisite for a
// certificate. A failed prerequisite or coefficient bound is Uncertified
// unless a direct high-order Jacobian sample witnesses invalidity.
template<int gdim, int ideg>
ElementValidityResult classify_element_validity(const MeshBase &msh,
                                                 int element);

} // namespace Metris

#endif
