# Request-20 connection-smoothing experiment

Experiment directory: `/Users/renat/MIT/HOMeshing/ExactQ2MetricField/ConnectionSmoothingFromAdaptedQ1/requested_000020_20260905`.

| Mesh | Geometry nodes | CG P1 solution DOFs | Triangles | SANS L2 projection error |
|---|---:|---:|---:|---:|
| Fully adapted Q1, StepDistance SV | 16 | 16 | 18 | 4.2080154620e-03 |
| Elevated Q2, all-node connection smoothing | 49 | 16 | 18 | 4.5208586107e-10 |

The error decreases by approximately 9.31 million times. These are numerical
errors, not a claim of exact zero.

## Setup

- Same manufactured field as the earlier experiments: analytical metric 9,
  Q2Coordinate scalar solution, parameter a=0.5, unit-square CAD.
- The existing manifest maps request 20 to metric scale 1/3 and a manufactured
  grid count of 16 topological vertices. The adapted mesh also has 16 vertices.
- Started directly from the regular square with 9 vertices and 8 triangles.
- Q1 adaptation: 20-iteration limit, StepDistance SV, objective p=1,
  regularization 1e-8, no barrier, quadrature order 5; all other controls are
  recorded in the experiment's run.sh. Adaptation reached its stopping criterion.
- Elevated without changing connectivity; 33 edge nodes were added.
- Simultaneous L-BFGS smoothing of interior vertex/edge-node coordinates and
  boundary vertex/edge-node CAD parameters (70 scalar optimization variables).
  Corners were fixed. The existing CAD-parameter trust intervals preserve
  boundary ordering. Invalid elements and boundary intersections are rejected.
- Objective: sum over elements of the quadrature approximation to
  integral_K ||Gamma_K - Gamma(M)||_M^2 dV_M. Uses the same order-5 quadrature
  and Riemannian volume weighting as StepDistance SV.
- Gradient approximation matches StepDistance SV: freeze M, its spatial first
  derivatives, and the complete quadrature volume factor, including det J.
  Differentiate only the mesh-implied connection. Re-evaluate the actual target
  and integration weights for trial values; retain the best valid state.
- Projection: existing SANS solver, CG P1 solution space on both geometries,
  UMFPACK, unchanged manufactured solution and projection settings. The driver
  now accepts an optional expected mesh count; existing four-level behavior is
  retained by default. This run passes 1.

## Smoothing and validity

- Initial integrated connection energy: 2.05213689824186.
- Final integrated connection energy: 1.3751040476211822e-14.
- 30 optimizer evaluations; NLopt result 4 (coordinate-tolerance convergence).
- Final mesh passed the conservative element-validity and boundary checks.
- Minimum Jacobian determinant lower bound: 0.057807189779950127.
- Minimum normalized Jacobian lower bound: 0.74125684616695198.
- The frozen gradient differs from a finite difference of the fully varying
  integral, as expected; initial directional relative mismatch was 1.95e-4.
  Separate frozen-surrogate gradient tests pass in 2D and 3D.

## Artifacts

The experiment directory contains `adapted_Q1.meshb`,
`elevated_unsmoothed_Q2.meshb`, `connection_smoothed_Q2.meshb`,
`adaptation.log`, `adaptation_stdout.log`, `smoothing.log`,
`optimization_history.csv`, `projection_Q1.dat`, `projection_Q2.dat`,
both SANS logs, staged SANS input links, and `run.sh`.
Two initial setup-failure logs were retained (output prefix and the former
four-level projection requirement); neither failure altered the final results.

Reproduction entry point: `progs/runConnectionRequest20.sh` in Metris.
Use a fresh EXPERIMENT directory and TAG to avoid overwriting existing results.
