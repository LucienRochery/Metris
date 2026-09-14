# Q1 adaptation and local Q2 smoothing: targets 100–3200

Experiment root:
`/Users/renat/MIT/HOMeshing/ExactQ2MetricField/Q1ElevatedQ2Sequences/start100_n06_20260905`.

## Results

SANS CG P1 L2 projection of the unchanged manufactured Q2Coordinate scalar.
All three variants have the same actual solution DOFs and connectivity at each
resolution. Q2 geometry nodes are not counted as solution DOFs.

| Request | Actual DOFs | Triangles | Q1 error | Q2 + SV error | Q2 + connection error |
|---:|---:|---:|---:|---:|---:|
| 100 | 100 | 166 | 5.3630253570e-4 | 2.1291527207e-3 | 2.2644968172e-16 |
| 200 | 201 | 354 | 2.5155590438e-4 | 1.1253983308e-3 | 2.1799083204e-16 |
| 400 | 409 | 747 | 1.1453446694e-4 | 8.8243966650e-4 | 2.1460284827e-16 |
| 800 | 801 | 1505 | 5.7493425731e-5 | 5.5325202746e-4 | 2.1677316910e-16 |
| 1600 | 1616 | 3093 | 2.8073940577e-5 | 3.8651234932e-4 | 2.1630762082e-16 |
| 3200 | 3226 | 6251 | 1.3755528570e-5 | 2.6545002331e-4 | 2.1715873645e-16 |

Connection-smoothed Q2 reaches roundoff-level projection error at every level.
Its flat convergence curve reflects that numerical floor, not an observed
algebraic convergence rate. SV smoothing has higher projection error than Q1
at every level in this experiment. These are results of the current production
node-by-node optimizer, not claims about all possible objective minima.

## Mesh locations

Under the experiment root, each mesh is named `mesh.meshb`:

- `Q1/requested_NNNNNN/mesh.meshb`: continued production Q1 adaptation.
- `Q2_Elevated/requested_NNNNNN/mesh.meshb`: shared unsmoothed elevation.
- `Q2_SV/requested_NNNNNN/mesh.meshb`: SV-smoothed Q2.
- `Q2_Connection/requested_NNNNNN/mesh.meshb`: connection-smoothed Q2.

`NNNNNN` is `000100`, `000200`, `000400`, `000800`, `001600`, or `003200`.

## Setup and verification

- Six target levels, 100, 200, 400, 800, 1600, 3200. Q1 starts from the regular
  2x2 square (9 vertices, 8 triangles); each next level starts from the preceding
  adapted Q1 mesh. No smoothed Q2 mesh feeds the Q1 sequence.
- Analytic metric 9, scale `1/(sqrt(request)-1)`, production StepDistance SV,
  objective p=1, order-5 quadrature, regularization 1e-8, barrier zero.
  Production adaptation/optimization flags are recorded in each command receipt.
- At each level, both smoothing branches start from the same saved elevation.
  Their initial snapshots are hash-identical. Vertices and edge nodes are
  optimized one at a time; CAD corners remain fixed, boundary nodes slide.
- Each smoothing branch uses a single uninterrupted production call with
  persistent node tags, neighbor reactivation, and production stopping
  (`nmov = 0`). The 100000-sweep safety cap was never reached. No extra full-node
  restart loop, topology modifications, or simultaneous global solve is used.
- Connection derivatives freeze M, its physical first derivatives, and the
  entire quadrature weight; trial values re-evaluate target and weights.
- All twelve final Q2 meshes have exactly the elevated connectivity and node
  counts, preserve the square boundary, and have positive conservative Jacobian
  bounds. Minimum normalized bounds range from 0.1996 to 0.4538 for SV and from
  0.8307 to 0.9611 for connection.
- Projection uses the existing SANS CG P1 implementation, UMFPACK, and unchanged
  order-4 projection quadrature. All 18 solves succeeded; mesh hashes were
  unchanged by projection.

| Request | SV sweeps | Connection sweeps | Final connection integral |
|---:|---:|---:|---:|
| 100 | 85 | 321 | 1.8975329252e-27 |
| 200 | 31 | 1535 | 1.5438270400e-26 |
| 400 | 43 | 303 | 3.6985891704e-26 |
| 800 | 89 | 337 | 2.0436010559e-25 |
| 1600 | 67 | 993 | 9.9989122537e-25 |
| 3200 | 108 | 1280 | 3.0327885576e-24 |

## Artifacts and reproduction

- `convergence.png` and `convergence.pdf`: error versus actual DOF^-1/2.
- `convergence.csv` and `results.json`: exact error values and absolute mesh paths.
- `validation.json`: per-mesh counts, validity bounds, objectives, and sweep counts.
- `manifest.json`: scaling, initial mesh, executable hashes, and method settings.
- `projection_*.dat` and `projection_*.log`: SANS outputs.
- Every mesh directory contains its generation log and exact command receipt.
- `run.py` and `runConnectionSequence.py`: saved reproduction scripts.

Source entry point: `/Users/renat/Metris/progs/runQ1ElevatedQ2Sequence.py`.
Use fresh ROOT and TAG values in that script for a new experiment. Existing
completed runs can be inspected/resumed by stage without overwriting logs.
