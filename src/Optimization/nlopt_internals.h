// Metris: high-order metric-based non-manifold tetrahedral remesher
// Copyright (C) 2023-2025, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See /License.txt or http://www.opensource.org/licenses/lgpl-2.1.php

// This file is adapted from the NLopt library, specifically the luksan_pnet algorithm.
// It is verbatim copy-paste of code contained in headers which is not in the public NLopt API.
// We use this for the modified luksan_pnet which avoids dynamic allocations.
#ifndef NLOPT_INTERNALS_HXX
#define NLOPT_INTERNALS_HXX

// Copied from libs/.nlopt/src/nlopt_stopping.h
/* stopping criteria */
typedef struct
{
  unsigned n;
  double minf_max;
  double ftol_rel;
  double ftol_abs;
  double xtol_rel;
  const double *xtol_abs;
  const double *x_weights;
  int *nevals_p, maxeval;
  double maxtime, start;
  int *force_stop;
  char **stop_msg; /* pointer to msg string to update */
} nlopt_stopping;


#endif // NLOPT_INTERNALS_HXX