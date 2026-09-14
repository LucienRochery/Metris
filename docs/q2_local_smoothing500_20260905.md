# Q2 node-by-node smoothing of the inspected request-500 Q1 mesh

Experiment root:
`/Users/renat/MIT/HOMeshing/ExactQ2MetricField/Q1Inspection/requested_000500_20260905/Q2LocalSmoothing`.

## Final meshes

- Shared starting elevation: `elevated_Q2.meshb`.
- StepDistance SV: `step-distance-sv/smoothed_Q2.meshb`.
- Connection: `connection_production_stop/smoothed_Q2.meshb`.

Both final meshes descend independently from the same elevation of
`../adapted_Q1.meshb`. Both have 485 topological vertices, 1395 edge nodes,
1880 geometry points, and 911 triangles. Exact P2 connectivity is unchanged.
No insertion, collapse, swapping, global simultaneous optimization, or SANS
projection was performed.

## Algorithm and stopping

Both use `smoothInterior_Ball` and its production single-node interior/CAD
Newton solves, geometry validity checks, regional acceptance, and neighbor
reactivation. Vertices and edge nodes are eligible; corners remain fixed.
The active set persists across the sweeps of each production call.

The user selected **production stopping**, not numerical-stagnation restarts:

- StepDistance SV stopped after 170 sweeps: final sweep had 7 accepted moves,
  but 0 moves large enough to trigger reactivation (`nmov = 0`).
- Connection initially ran 2000 sweeps. It was then resumed from that saved
  production mesh for one uninterrupted call, stopping after another 3905
  sweeps: final sweep had 0 accepted moves and `nmov = 0`.
- Reactivation uses the configured 0.5% relative improvement threshold.
  Regional acceptance retains production's 1e-7 relative improvement test;
  local Newton termination retains its 1% relative objective tolerance.
  These are algorithmic stopping conditions, not proof of mathematical minima.
- Experimental extra-restart jobs in `until_stalled` were terminated when the
  user chose production stopping. Their checkpoints were **not** used for the
  final meshes above. `connection_preliminary_200` is also preliminary, not final.
- The connection result's `stalled:false` log field refers only to the unused
  outer numerical-stagnation test across a complete batch; its production loop
  demonstrably ended with `nsuccess = 0 nmov = 0`.

## Objectives and approximations

Metric 9; scale 0.046814989529053899; order-5 quadrature; objective p=1;
StepDistance SV regularization 1e-8, barrier zero.

Connection integrates the squared target norm of the connection difference
with exactly the distance quadrature and Riemannian volume measure. Gradient
and Hessian freeze M, its physical first derivatives, and the full quadrature
weight including det J. Only the mesh connection is differentiated. Trial
values re-evaluate target and weights. The exact frozen geometry Hessian was
added for production's local Newton solve; no finite differences are used in
the runtime evaluator. CAD chain rules use the connection Hessian's sym2idx
packing; the legacy StepDistance CAD path is unchanged.

| Objective (not comparable across rows) | Elevated | Final |
|---|---:|---:|
| StepDistance SV integral | 133.04862208119678 | 106.71701374138539 |
| Connection integral | 2.052800941376423 | 0.007914148029221604 |

## Verification

| Check | StepDistance SV | Connection |
|---|---:|---:|
| Vertices moved (>1e-13) | 477 | 481 |
| Edge nodes moved (>1e-13) | 1386 | 1395 |
| Maximum displacement from elevation | 0.01472254298 | 0.0008876562054 |
| Minimum Jacobian determinant bound | 9.13156733e-5 | 2.80488528e-4 |
| Minimum normalized Jacobian bound | 0.140387719 | 0.921694427 |

All corner coordinates are unchanged, all boundary nodes remain on the square,
and both meshes have positive conservative Jacobian bounds. See `validation.json`.

Tests passed: connection value/gradient/Hessian and quadrature (12 cases),
node locality/reactivation (5), high-order StepDistance (8), curved-CAD
smoothing including connection (6). `git diff --check` passed.

Reproduction sources in Metris: `progs/runQ1InspectionLocalSmoothing.py`,
`progs/finishQ2ConnectionProduction.py`, `progs/auditQ2LocalSmoothing.py`.
Exact commands and logs are preserved in the respective experiment folders.
