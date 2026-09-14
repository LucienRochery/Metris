# Production adaptation and connection postprocessing: requests 500–4000

Results: `/Users/renat/MIT/HOMeshing/ExactQ2MetricField/ConnectionPostprocessingSequence/start500_20260905`.

## Projection results

SANS CG P1 solution space on all three geometry variants; DOF means actual
topological vertex count, not Q2 geometry-node count. Projection settings are
unchanged from the existing Q2Coordinate experiment (quadrature order 4).

| Request | Q1 DOFs | Q2 DOFs (both variants) | Q1 L2 error | Original Q2 L2 error | Connection-smoothed Q2 L2 error |
|---:|---:|---:|---:|---:|---:|
| 500 | 485 | 481 | 1.1601432785e-4 | 1.0438218032e-3 | 1.0436640100e-3 |
| 1000 | 1029 | 993 | 4.4290595430e-5 | 6.8743881362e-4 | 6.8744452914e-4 |
| 2000 | 2013 | 1953 | 2.2220376745e-5 | 5.0893237554e-4 | 5.0937420301e-4 |
| 4000 | 4046 | 3920 | 1.0854855921e-5 | 3.5461225130e-4 | 3.5461810570e-4 |

Q1 has lower projection error at every level. Q2 postprocessing changes the
projection error by less than 0.09%; the two Q2 curves nearly overlap.
These runs do not establish the behavior of a fully minimized connection
objective: every smoothing run stopped on coordinate tolerance with a large
remaining objective and gradient.

## Mesh generation and smoothing

- Independent continued Q1 and Q2 production sequences, requests 500, 1000,
  2000, 4000; initial square mesh has 9 vertices and 8 triangles.
- Analytical metric 9, usual production scale `1/(sqrt(request)-1)`;
  StepDistance SV, objective p=1, regularization 1e-8, barrier zero, order-5
  quadrature. Full production adaptation/optimization settings and executable
  commands are in `run.py` and each `*.completed.json` receipt.
- Only the original adapted mesh seeds the next resolution. Q2 connection
  smoothing is a separate postprocessing branch without topology changes.
- Simultaneous smoothing of all free geometry nodes: interior vertices and
  edge nodes, with boundary vertex/edge nodes sliding on CAD and corners fixed.
  Existing validity checks and boundary-ordering trust intervals are retained.
- Connection objective uses the same quadrature and Riemannian measure as the
  distance objective. Gradient freezes the target metric, its physical first
  derivatives, and the entire quadrature volume weight (including det J).
  Only the mesh-implied connection is differentiated; trial values re-evaluate
  the target and volume weight.
- First pass (`Q2Connection`) used the existing 1e-8 coordinate tolerance.
  It stopped after small steps and barely changed projection errors.
- Reported results use `Q2ConnectionRefined`: coordinate tolerance 1e-12,
  fixed objective/gradient scaling by max(1, initial objective), and maximum
  5000 evaluations. This changes numerical conditioning/stopping only, not the
  objective, frozen-gradient approximation, or admissible node motions.
  Both passes and all original meshes are preserved.

| Request | Initial connection integral | Final connection integral | Reduction | Optimizer evaluations |
|---:|---:|---:|---:|---:|
| 500 | 5.1051567036e6 | 4.6393130835e5 | 90.91% | 748 |
| 1000 | 3.2454077815e6 | 2.0934776299e6 | 35.49% | 71 |
| 2000 | 5.6030359016e6 | 2.2889223753e5 | 95.91% | 1194 |
| 4000 | 2.4520128632e6 | 1.1971331589e6 | 51.18% | 350 |

All four have NLopt result 4 (coordinate-tolerance stop); saved outputs retain
the best valid state. Large residuals and many rejected invalid trials remain.
The connection integral reduction is therefore not evidence of convergence to
a connection-matched mesh, nor a general test against connection conformity.

## Artifacts and reproduction

- `convergence.png`, `convergence.pdf`: log-log error versus actual DOF^-1/2.
- `convergence.csv`: exact results, counts, metric scales, mesh paths.
- `Q1`, `Q2`, `Q2ConnectionRefined`: reported meshes, logs and command receipts.
- `Q2Connection`: preserved first-pass smoothing results.
- `projection_*.dat`, `projection_*.log`, `SANSInput`: SANS data and staged meshes.

Entry point: `/Users/renat/Metris/progs/runConnectionSequence.py`.
Reported pipeline options: `--smoothing-variant Q2ConnectionRefined --tight-smoothing`.
The saved `run.py` is a snapshot of that entry point. Plotting uses the existing
`ExactQ2MetricField/plot_l2_projection_errors.py` styling and matplotlib.
To reproduce in a new location, set a fresh `--root` and a unique `TAG` in the
script; the latter avoids colliding with preserved SANS output directories.
