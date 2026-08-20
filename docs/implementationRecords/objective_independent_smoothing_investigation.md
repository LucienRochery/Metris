# Objective-independent smoothing investigation

## Scope

This record investigates why smoothing behaves differently for the
objective-driven SizeShape and StepDistance variants, and defines controlled
tests for optimizer tuning.  It does not change the Classical objective or its
adaptation algorithm.

## Paper-run evidence

The authoritative raw paper data are in:

`IMR-2027-paper/runAndResults/moess_runs/BL`

The following totals sum every row of every `outputAdaptStats_MOESS_a*.txt`
file for each BL target size.  They count the element-ball fallback attempted
after all eligible insertion and collapse operations on an element fail.

| Target | SizeShape attempts / successes | StepDistanceBarrier attempts / successes | StepDistanceSV attempts / successes |
|---:|---:|---:|---:|
| 500 | 224,529 / 79,713 | 0 / 0 | 0 / 0 |
| 1,000 | 86,136 / 34,071 | 0 / 0 | 0 / 0 |
| 2,000 | 174,246 / 73,102 | 0 / 0 | 0 / 0 |
| 4,000 | 290,651 / 125,279 | 0 / 0 | 114 / 80 |
| 8,000 | 799,718 / 332,177 | 361 / 197 | 114 / 76 |
| 16,000 | 2,634,634 / 1,025,476 | 39,849 / 16,259 | 1,686 / 937 |
| 32,000 | 2,524,245 / 1,074,395 | 34,147 / 16,112 | 6,601 / 3,847 |

The success rates of attempted smooths are not poor.  The dominant discrepancy
is whether smoothing is attempted at all.

## Code audit findings

1. The objective-driven adaptation loop processes elements from largest
   objective value (worst) to smallest (best).
2. After all out-of-range edge operations fail, fallback smoothing is attempted
   only when the raw elemental objective is less than `0.01`.
3. The same `0.01` is used for every objective even though their numerical
   scales differ.  This explains the paper-run discrepancy above.
4. Dividing by the minimum objective is not a robust replacement: every
   objective has ideal value zero, so a good mesh can have a zero or arbitrarily
   small minimum.
5. Mesh-wide quality smoothing accepts only objective-decreasing point moves,
   but its neighborhood-reactivation test historically required an absolute
   reduction greater than `opt_smoo_tol` (`0.005`).  This is another
   objective-scale-dependent decision.
6. The low-level Newton stopping tolerance is relative.  Most cavity acceptance
   gates are also relative or strict comparisons, so they are not the primary
   source of the observed discrepancy.
7. The manuscript's objective-driven algorithm description goes directly from
   a failed prioritized-element traversal to mesh-wide optimization.  It does
   not describe the element-level fallback, although the archived runs show
   that the fallback was active.

## Implemented experiment controls

- `opt-smoo-tol` is interpreted as the relative reduction required to
  reactivate neighboring vertices in quality-based mesh-wide smoothing.
- `--adp-quality-smoothing` enables the element-level fallback without changing
  the mesh-wide optimization pass.  The fallback is disabled by default.
- Unit tests verify that the reactivation decision is unchanged when an
  objective is multiplied by factors from `1e-12` through `1e12`.

## Controlled tuning experiments

Run the same saved input mesh and metric for each objective.  Record wall time,
all operation attempts and successes, initial/final objective sum and maximum,
validity, final element count, and the metric-edge-length distribution.

1. **Fallback A/B test**: fallback enabled versus disabled, with all mesh-wide
   optimizer settings fixed.
2. **Outer-pass sweep**: `opt-niter` in `{1, 3, 5, 10}` with the fallback
   disabled.
3. **Inner-sweep sweep**: `opt-smoo-niter` in `{1, 3, 5, 10}`.
4. **Relative revisit tolerance sweep**: `opt-smoo-tol` in
   `{1e-2, 5e-3, 1e-3, 1e-4}`.
5. Repeat the smallest useful matrix on 500, 4,000, and 32,000 targets so that
   a tuning choice is not inferred from a single resolution.

The existing `test_stepdistance_smoothing_scan` diagnostic can be pointed at a
saved case with `METRIS_SMOOTHING_CASE_DIR`, `METRIS_SMOOTHING_CAD`, and
`METRIS_SMOOTHING_ITERATION`.  `METRIS_SMOOTHING_MAX_ATTEMPTS` limits isolated
point scans for quick checks.  Its definition of a substantive move now uses
the same relative criterion as production smoothing.

A one-point diagnostic smoke test on the archived 500-target iteration-16 mesh
passed in isolated, locally committed, and globally committed modes.  The full
StepDistance production smoothing pass made 443 substantive moves and reduced
the configured whole-mesh objective from `525.11817` to `477.42895`, confirming
that disabling the adaptation fallback does not suppress mesh-wide smoothing.

Initial A/B tests on archived BL iteration-16 inputs found:

| Target | Fallback | Wall time | Attempts / successes | Final triangles | SizeShape average | StepDistance average |
|---:|:---:|---:|---:|---:|---:|---:|
| 500 | on | 1.80 s | 5,141 / 2,113 | 1,020 | 7.20349013e-2 | 1.97683197e-1 |
| 500 | off | 0.47 s | 0 / 0 | 1,020 | 7.20058235e-2 | 1.97631062e-1 |
| 4,000 | on | 2.58 s | 2,383 / 921 | 7,429 | 7.17229223e-2 | 2.08831050e-1 |
| 4,000 | off | 2.06 s | 0 / 0 | 7,429 | 7.17238561e-2 | 2.08818366e-1 |

Both comparisons produced the same insertion/collapse successes and final
element count.  The retained objective differences after mesh-wide smoothing
are negligible and neither mode is uniformly better.  This justifies making
the fallback opt-in while larger tests continue.

The next tuning work should focus on `opt-niter`, `opt-smoo-niter`, and the
relative revisit tolerance.  Objective-specific parameter defaults should be
introduced only if the same relative
criteria still produce materially different convergence behavior.

### Preliminary optimizer sweep

On the same archived 500-target, iteration-16 input with fallback smoothing
disabled, the SizeShape outer-pass sweep gave:

| `opt-niter` | Metris runtime | SizeShape average | StepDistance average |
|---:|---:|---:|---:|
| 1 | 0.271 s | 7.21281984e-2 | 1.97754986e-1 |
| 3 | 0.398 s | 7.20276971e-2 | 1.97658212e-1 |
| 5 | 0.457 s | 7.20058235e-2 | 1.97631062e-1 |
| 10 | 0.452 s | 7.20058235e-2 | 1.97631062e-1 |

Thus five outer passes are sufficient for this case; ten is identical because
the optimizer terminates early.

With `opt-niter=5`, changing the maximum number of smoothing sweeps performed
before each swap pass showed that one sweep can be preferable to ten:

| Objective used by run | `opt-smoo-niter` | Runtime | Final triangles | Objective average | Objective maximum |
|:---|---:|---:|---:|---:|---:|
| SizeShape | 1 | 0.408 s | 1,020 | 6.94972913e-2 | 6.08788097e-1 |
| SizeShape | 10 | 0.376 s | 1,020 | 7.20058235e-2 | 6.08788097e-1 |
| StepDistance | 1 | 9.95 s (debug) | 1,014 | 1.53199236e-1 | 5.12964418e-1 |
| StepDistance | 10 | 14.4 s (debug) | 1,008 | 1.53400959e-1 | 5.27652673e-1 |

The StepDistance timings come from the same debug build and are comparable to
each other, but not to the SizeShape release timings.  In both objectives, one
smoothing sweep before returning to swapping produced a lower final average;
for StepDistance it also produced a lower maximum.  This is evidence that
smoothing/swap interleaving matters.  It is not yet sufficient to change the
default: repeat at 4,000 and 32,000 targets before doing so.

## Deferred alternative fallback gate

If the fallback proves useful, replace the raw objective threshold with a
rank/quantile-based eligibility rule.  A quantile is invariant under positive
rescaling and monotone reparameterization of an objective and is more robust
than normalization by a near-zero minimum.  This should be implemented only
after the fallback A/B test establishes that the operation is worth retaining.
