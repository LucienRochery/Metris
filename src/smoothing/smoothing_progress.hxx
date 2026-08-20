// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef METRIS_SMOOTHING_PROGRESS_HXX
#define METRIS_SMOOTHING_PROGRESS_HXX

#include <cmath>

namespace Metris
{

// Quality-based objectives are nonnegative and are minimized.  Expressing the
// stopping decision as a relative reduction keeps it unchanged under any
// positive rescaling of the objective.
inline bool smoothing_reduction_is_substantive(double value_before,
                                                double value_after,
                                                double relative_tolerance)
{
  if (!std::isfinite(value_before) || !std::isfinite(value_after)
      || !std::isfinite(relative_tolerance))
  {
    return false;
  }
  if (value_before <= 0.0 || value_after < 0.0
      || relative_tolerance < 0.0)
  {
    return false;
  }

  return value_after < (1.0 - relative_tolerance) * value_before;
}

inline bool smoothing_neighborhood_should_be_reactivated(
    double average_before,
    double average_after,
    double maximum_before,
    double maximum_after,
    double average_relative_tolerance,
    double maximum_relative_tolerance)
{
  return smoothing_reduction_is_substantive(
             average_before,average_after,average_relative_tolerance)
      || smoothing_reduction_is_substantive(
             maximum_before,maximum_after,maximum_relative_tolerance);
}

} // namespace Metris

#endif
