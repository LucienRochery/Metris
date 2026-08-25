# High-Order Mesh Generation Implementation Plan

## Purpose and provenance

This document records the high-order implementation plan developed immediately
after the initial audit of Metris's P2 support. It is intended to be the starting
reference for a fresh implementation task.

The plan predates the completed objective-quadrature refactor. Consequently,
file locations and individual restrictions must be verified against the current
source before changing code. The architectural intent remains applicable: first
restore and qualify the existing Classic P2 path, then extend the common
objective-driven integration machinery to P2.

## Audit baseline

The audit established that P2 support is substantial rather than merely
experimental. The historical Classic workflow was designed to:

1. adapt a P1 mesh;
2. elevate it to P2;
3. treat P2 edge control points as smoothing variables;
4. evaluate the size-invariant distortion using the P2 mapping; and
5. reject trial configurations whose scaled Jacobian Bernstein coefficients
   fall below `jtol`.

At the time of the audit, the principal blockers were:

- 2D boundary P2 control-point smoothing assumed two adjacent triangles;
- the objective-driven SizeShape and StepDistance paths were restricted to P1;
- final metric-cost reporting was P1-only and could make a successful P2 run
  end with an exception;
- an already-P2 3D adaptation could reach an unimplemented P2 triangular-face
  projection;
- several local smoothing and topological acceptance paths were P1-only or
  evaluated only the underlying P1 quality;
- some P2 regressions caught and printed exceptions without failing; and
- the top-level high-order workflow invoked optimization twice.

These are historical audit findings, not assertions that every issue remains in
the current source. Phase 1 begins by checking them again.

## Governing principles

1. P2 is the initial supported high-order target. P3 and higher degrees remain
   outside the implementation scope until P2 is complete and qualified.
2. Quadrature accuracy and element validity are logically independent.
   Quadrature samples the objective; it must never be used as proof that a
   high-order element is valid.
3. All pointwise objectives should consume the same degree-aware geometric and
   quadrature data. Objective choice should not create a separate high-order
   integration path.
4. High-order validity should be enforced by one common certificate in trial
   moves, accepted moves, cavity operations, reconnection, and final checks.
5. Each phase should be implemented as small, independently tested changes.

## Phase 1: Make the baseline trustworthy

1. Re-audit the current source against every blocker listed above and record
   which items are fixed, still present, or changed by later refactors.
2. Make every high-order regression fail on an unexpected `MetrisExcept`.
3. Generalize or bypass the P1-only metric-cost postprocessor so a successful
   P2 run does not end with an unrelated exception.
4. Remove any duplicate high-order optimization invocation.
5. Establish four small deterministic baseline mesh workflows:

   - 2D P1 to P2 elevation;
   - 2D native P2 input;
   - 3D P1 to P2 elevation; and
   - 3D native P2 input.

6. Run these baseline meshes in debug and release builds before modifying
   smoothing or objective behavior.

## Phase 2: Recover Classic P2 smoothing

1. Use only the Classic size-invariant distortion objective.
2. Confirm that interior P2 edge control points are selected as smoothing
   variables.
3. Fix the 2D boundary-edge control-point neighborhood. A boundary edge has one
   adjacent triangle and must not be forced through a two-triangle patch.
4. Exercise the existing CAD-edge boundary optimizer for a boundary P2 control
   point.
5. Add focused tests for:

   - one interior 2D P2 edge control point;
   - one boundary 2D P2 edge control point;
   - a complete small 2D P2 mesh; and
   - one interior 3D P2 edge control point and its tetrahedral edge shell.

6. Establish that Classic P2 smoothing succeeds before involving objective
   variants, CAD surfaces, or full adaptation.

## Phase 3: Consolidate high-order validity

### Full-dimensional elements

For a degree-`p` simplex in topological dimension `d`,

\[
\deg(\det DF_K) \leq d(p-1).
\]

Writing the determinant in the Bernstein basis,

\[
\det DF_K(\xi) = \sum_\alpha c_\alpha B_\alpha(\xi),
\]

the nonnegative partition-of-unity property gives

\[
\det DF_K(\xi) \geq \min_\alpha c_\alpha.
\]

Strictly positive Bernstein coefficients therefore give a sufficient
certificate of a positive Jacobian determinant everywhere. This condition is
not necessary: a nonpositive coefficient means that the element is not
certified by the current bound, not necessarily that it is inverted.

### Incremental validity work

1. Give the validity result explicit semantics:

   - `certified`;
   - `invalid`; and
   - `uncertified`.

2. Unit-test coefficient generation against direct Jacobian evaluation for P2
   triangles and tetrahedra. Retain the existing mesh-wide coefficient test,
   but add deterministic affine and curved single-element cases covering
   Lagrange and Bezier storage and off-lattice evaluations.
3. Build one full-dimensional validity routine that owns the P1 orientation
   prerequisite, normalization, coefficient classification, invalidity
   witnesses, and diagnostics. Add focused cases for all three validity
   outcomes here, where the production classifier rather than test-only logic
   can be exercised. Keep this API degree-generic while qualifying P2 first.
4. Integrate the routine first into Classic P2 smoothing trials and accepted
   moves in 2D, followed by the existing 3D tetrahedral edge-shell path.
5. Use the same validity routine in completed-element checks for:

   - cavity construction and correction;
   - reconnection; and
   - final topology and validity checks.

   Provisional elements that have only their P1 skeleton retain their P1 check;
   certification is required once all high-order nodes have been assigned.
6. Keep the final high-order validity check unconditional rather than dependent
   on `dbgfull`. This was completed in Phase 1 and will be migrated to the new
   result type rather than reintroduced conditionally.
7. Treat embedded P2 surface triangles separately. Do not claim a polynomial
   certificate from samples of the surface-Jacobian magnitude; use the oriented
   polynomial described below when surface certification is implemented.
8. Add Bernstein subdivision later to refine inconclusive bounds and reduce
   conservative false rejections.

**Scope decision (2026-08-25).** With the full-dimensional P2 triangle and
tetrahedron path qualified, items 7 and 8 are deliberately deferred. The
immediate priority is the planar 2D high-order path, retaining 3D coverage when
the extension is straightforward. Embedded surface certification and
Bernstein subdivision will be resumed later rather than treated as gates for
Phase 4.

### Embedded surface elements

The magnitude of a surface Jacobian contains a square root and is not a
polynomial, so samples of that magnitude cannot be converted into rigorous
Bernstein bounds.

For a curved surface triangle, instead certify the polynomial

\[
g(\xi) = n_0 \cdot
\left(
\partial_{\xi_1}F_K \times \partial_{\xi_2}F_K
\right),
\]

where `n_0` is a fixed, consistently oriented reference normal. Positive
Bernstein coefficients of `g` certify a nonzero, consistently oriented surface
Jacobian. A future subdivision implementation may choose a separate reference
normal on each subpatch.

## Phase 4: Generalize SizeShape to P2

1. Remove the P1 restriction for SizeShape value evaluation.
2. Evaluate the P2 geometry and metric at the shared quadrature points. Keep
   the metric-degree choice explicit: the current compatibility contract uses
   `AsDeg::P1` for FE metrics even on a P2 mesh, whereas analytical metrics are
   evaluated directly at the physical P2 sample location. Moving FE metric
   interpolation to `AsDeg::Pk` is separate work and must not happen
   implicitly as part of geometry generalization.
3. Compare element values with an independent dense reference sampler.
4. Generalize differentiated evaluation for an active P2 edge control point.
5. Check gradients against both automatic differentiation of the frozen-sample
   value and centered finite differences of that same functional.
6. Check Hessians against automatic differentiation of the validated gradient
   and centered finite differences of the gradient.
7. Cover Lagrange and Bezier geometry representations where both are supported.
8. Enable SizeShape P2 smoothing only after value, gradient, Hessian, and
   validity tests pass.

The shared objective-quadrature traversal should remain the only integration
path. Degree awareness belongs in sample preparation and basis gradients, while
the SizeShape policy continues to provide its own pointwise value and
derivatives.

## Phase 5: Generalize StepDistance to P2

1. Replace any remaining hard-coded degree-one geometry evaluation.
2. Remove assumptions that the active variable must be a P1 vertex.
3. Implement and validate the basic StepDistance variant first.
4. Add the collapse-barrier variant.
5. Add ShapeVolume.
6. Add CavityTargetAverage as its deliberately distinct aggregation convention.
7. Verify each variant independently in 2D before enabling it in 3D.
8. For each variant, test values, gradients, Hessians, validity rejection, and
   smoothing behavior for P2 edge control points.

Objective selection should occur only in the lower-level pointwise policy. The
integrator should consume the common objective sample and `psi` abstraction in
the same way for SizeShape and every StepDistance variant.

## Phase 6: CAD and boundary support

1. Complete 2D curved CAD-boundary control-point smoothing.
2. Implement P2 triangular-face projection for 3D localization.
3. Test the polynomial orientation certificate for P2 surface triangles.
4. Add small CAD-curve and CAD-surface regression meshes.
5. Verify that CAD parameters, projected coordinates, and high-order validity
   remain consistent after every accepted move.
6. Add checks for curved-boundary self-intersection where local Jacobian
   positivity alone is insufficient.

## Phase 7: Generalize adaptation operations

Proceed one operation at a time, retaining a P1 compatibility test beside each
new P2 test.

1. Make insertion acceptance use both the P2 objective and P2 validity.
2. Do the same for collapse.
3. Generalize face and edge swaps.
4. Implement P2 length smoothing or disable it explicitly; it must not silently
   return while appearing to be supported.
5. Generalize cavity-targeted smoothing to high-order variables.
6. Ensure all neighboring high-order control points and affected elements are
   reactivated after a successful move.
7. Confirm that newly created cavity elements share P2 edge control points
   correctly and pass the common validity certificate.

## Phase 8: End-to-end qualification

The final qualification matrix should cover:

- 2D and 3D;
- analytical and FE metrics;
- P1 to P2 elevation and native P2 input;
- interior, CAD-curve, and CAD-surface variables;
- Classic, SizeShape, and every StepDistance variant;
- historical and real objective quadrature;
- optimization-only and full adaptation workflows; and
- debug and release builds.

For every configuration, record:

- successful completion and exit status;
- minimum high-order validity coefficient;
- objective decrease and accepted/rejected smoothing moves;
- topology-operation attempts and successes;
- output mesh degree and entity counts; and
- final independent validity verification.

## Initial restoration milestone

The shortest first milestone is a trustworthy Classic P2 generator:

1. repair 2D boundary control-point smoothing;
2. remove the P1-only final metric-cost failure;
3. eliminate duplicate high-order optimization;
4. make the P2 regressions hard-failing; and
5. demonstrate the four baseline workflows with unconditional high-order
   validity checks.

Only after that milestone should the objective-driven SizeShape and
StepDistance policies be enabled for P2 smoothing.

## Phase 1 current-source audit and execution record

This section is the living, numbered implementation record for the high-order
work. The immediate target is P2 in full-dimensional 2D. The same step covers
full-dimensional 3D whenever the implementation and focused test extend
naturally. P3 is a design constraint, not a current acceptance requirement;
explicit P2 or 2D code remains acceptable when a premature generalization
would obstruct the first working path.

1. **Done (2026-08-20): re-audited the historical blockers.** The current
   classifications are:

   - **Still present:** a 2D P2 edge control point constructs its smoothing
     region by unconditionally requiring the triangle across the associated
     edge. `smoothInterior_Ball0` asserts that this neighbor is nonnegative,
     so a boundary P2 edge control point cannot use its correct one-triangle
     region. This is the first Phase 2 repair.
   - **Still present by design:** high-order SizeShape uses the historical
     nodal compatibility integration, while StepDistance enters the common
     objective traversal only for effective degree one. These restrictions
     remain assigned to Phases 4 and 5.
   - **Still present:** `projptfac<ideg>` throws for `ideg > 1`. A native P2
     3D workflow that requires localization on a triangular face can therefore
     still reach an unimplemented projection. The Phase 1 3D baseline mesh
     excludes CAD and does not claim to cover this path.
   - **Still present:** several reconnection, swapping, and cavity-acceptance
     paths evaluate `AsDeg::P1`, call `isvalideltP1`, or both. Cavity correction
     does contain degree-aware Bernstein-coefficient checks, but the acceptance
     system is not yet consolidated. This remains Phase 3 and Phase 7 work.
   - **Still present:** high-order length smoothing reports that adjacent
     high-order nodes are not updated. It must eventually be implemented or
     explicitly disabled as described in Phase 7.
   - **Changed since the original audit:** the regression manager records a
     caught `MetrisExcept` and later compares that state with the baseline.
     It was strengthened in step 2 below so an exception fails immediately,
     including while a new baseline is being initialized.
   - **Confirmed present, then fixed in steps 3 and 4:** duplicate top-level
     optimization and the P1-only final metric-cost postprocessor.
   - **Confirmed present, then fixed in step 5:** the final high-order validity
     check depended on `dbgfull`.

2. **Done (2026-08-20): made regression exceptions hard failures.**
   `RegressionTestManager` still records the exception text in `runs.json`, but
   now also emits an immediate Boost test error. An unexpected exception can no
   longer disappear when the first baseline for a machine is initialized.

3. **Done (2026-08-20): removed duplicate high-order optimization.** The main
   workflow now calls `optimMesh()` once after adaptation, degree elevation,
   and optional curving. The previous call inside the high-order branch was
   removed. The separate optimization used internally by curvature mode 3 is
   retained because it is part of that curving algorithm rather than the
   duplicated top-level pass.

4. **Done (2026-08-20): generalized final metric-cost integration to the mesh
   degree.** The postprocessor now dispatches on `msh.curdeg` and reuses the
   shared degree-aware quadrature-sample preparation to evaluate physical
   geometry and the target metric. It selects the same positive quadrature rule
   as objective integration through `objective_quadrature_order`, including the
   automatic dimension-dependent defaults.

   **Important metric-degree boundary:** degree-aware metric-cost integration
   currently applies to the **geometry**, not to the polynomial degree of a
   discrete FE metric. On a P2 or higher-order element, the physical position
   and Jacobian use the full `msh.curdeg` geometry, but an FE metric is still
   requested with `AsDeg::P1`. Consequently, only corner metric values enter
   its interpolation; metric values at high-order control points are not used.
   An analytical metric behaves differently: at nonvertex quadrature points it
   is evaluated directly at the physical point produced by the high-order
   geometry. This asymmetry is the intentional Phase 1 compatibility contract,
   not a claim that P2 FE metric interpolation is supported. Switching FE
   metrics to `AsDeg::Pk` is deferred, must be coordinated with the objective
   kernels, and requires dedicated value and derivative tests.

   The metric is evaluated in exponential SPD space and the caller's original
   metric storage space is restored, which is necessary because degree
   elevation can leave the front metric in logarithmic storage. Straight P1
   and elevated P2 baseline meshes produce the same cost.

5. **Done (2026-08-20): made the final full-dimensional high-order check
   unconditional.** Every successful `runMetris()` with `curdeg > 1` now calls
   `check_topo` before output, independently of `dbgfull`. For the present P2
   2D and 3D volume scope, this checks the scaled Jacobian Bernstein
   coefficients against `jtol` in addition to the underlying P1 orientation
   and topology. Phase 3 will replace this boolean legacy interface with the
   explicit `certified`/`invalid`/`uncertified` validity contract.
   A focused negative test displaces one P2 edge control point while leaving
   the P1 triangle unchanged and verifies that the top-level run throws at the
   final high-order check with `dbgfull` disabled.

6. **Done (2026-08-20): established four small deterministic baseline
   workflows.** `test_high_order_phase1_baselines` covers:

   - a two-triangle 2D P1 mesh elevated to P2;
   - the resulting two-triangle mesh re-entering as native P2 data;
   - a one-tetrahedron 3D P1 mesh elevated to P2; and
   - the resulting one-tetrahedron mesh re-entering as native P2 data.

   Each workflow runs the top-level driver with adaptation and optimization
   disabled, checks the degree and exact entity counts, exercises final metric
   cost, and reaches the unconditional high-order validity check. The baseline
   meshes pass in both release and debug builds. They intentionally establish
   the trustworthy Phase 1 floor; they do not yet claim that smoothing, CAD,
   or adaptation operations work at P2.

7. **Done in the Phase 2 record below:** begin with the one-triangle boundary
   neighborhood repair and focused control-point smoothing tests. The 2D
   boundary case remains the primary acceptance case; the existing 3D
   tetrahedral edge-shell construction is covered by the same test series.

## Phase 2 Classic P2 smoothing execution record

1. **Done (2026-08-20): made the edge-control-point smoothing region explicit.**
   `buildEdgeControlPointSmoothingRegion` is now the common production and test
   path. For a 2D interior edge it returns the seed and across-edge triangles;
   for a 2D boundary edge it returns only the seed triangle. The 3D branch
   retains the complete tetrahedral `shell3` traversal.

2. **Done (2026-08-20): removed a stale debug-only P1 assertion from the
   generic quality traversal.** With `TESTQUALITYALGO` enabled, both value and
   derivative traversal asserted `ideg_eff == 1` even when Classic
   `QuaFun::Distortion` was explicitly requested. This disabled the intended
   Classic P2 path only in assertion-enabled builds. The generic assertions
   were removed; the separate SizeShape and StepDistance restrictions remain
   in their objective-specific branches and are still assigned to Phases 4
   and 5.

   Focused regression reruns also found that the objective-adaptation
   dispatcher did not honor `adp_niter == 0`, unlike the Classic dispatcher.
   The common `adaptMesh2` entry point now returns immediately in that case, so
   Phase 1 and Phase 2 tests can reliably disable adaptation while isolating
   degree elevation and smoothing.

3. **Done (2026-08-20): added focused 2D edge-control-point tests.** On the
   elevated two-triangle mesh, the shared interior P2 control point is
   recognized by both triangles, receives the two-triangle region, and moves
   under one Classic distortion smoothing pass without increasing the regional
   objective. A boundary P2 control point receives exactly one triangle.

4. **Done (2026-08-20): exercised the CAD-edge boundary optimizer.** A small
   square CAD mesh is elevated to P2, all points except one boundary edge
   control point are frozen, and the production Classic smoother completes on
   the one-triangle region. The resulting CAD parameter and physical coordinate
   are independently checked with `EG_evaluate`.

5. **Done (2026-08-20): covered complete 2D and focused 3D behavior.** One
   Classic pass completes on the full small square P2 mesh. On an elevated cube
   mesh, a genuine interior P2 edge control point is found, its multi-tetrahedron
   shell is built through the common region helper, and the production smoother
   visits it without increasing the shell objective.

6. **Qualification.** The five focused cases in
   `test_high_order_classic_smoothing` pass in debug and release builds. They
   establish the Classic P2 smoothing floor in 2D and retain immediate 3D shell
   coverage. They do not enable the SizeShape or StepDistance P2 paths, and do
   not broaden the CAD work beyond an existing 2D CAD edge. Phase 3 validity
   consolidation is next.

## Phase 3 high-order validity execution record

1. **Done (2026-08-25): defined the explicit validity-result contract.**
   `ElementValidityStatus` distinguishes `Certified`, `Invalid`, and
   `Uncertified`. `ElementValidityResult` carries the normalized Bernstein lower
   bound and coefficient index, plus an optional direct normalized-Jacobian
   witness and sample index. Its initial conservative policy accepts only
   `Certified`; both a witnessed invalid element and an inconclusive bound are
   rejected. The contract is intentionally independent of coefficient
   generation and production call sites, which are addressed by the following
   Phase 3 steps. `test_high_order_validity_contract` fixes the default state,
   status predicates, diagnostic sentinels, and conservative acceptance policy.

2. **Done (2026-08-25): added deterministic P2 Jacobian-coefficient tests.**
   `test_high_order_p2_jacobian_coefficients` constructs one triangle and one
   tetrahedron directly in memory, each in both affine and genuinely curved
   form. It converts the P2 geometry from Lagrange to Bezier storage, compares
   the generated coefficient formulas with `ccoef_eval` from both bases, and
   evaluates the resulting degree-2 triangle and degree-3 tetrahedron
   determinant polynomials at four and five off-lattice barycentric points,
   respectively. At every point, the three coefficient routes reproduce the
   determinant of the directly evaluated geometry Jacobian, and the Lagrange
   and Bezier geometry values and Jacobians agree. The affine cases additionally
   require a constant coefficient field; the curved cases require a positive,
   nonconstant field.

3. **Done (2026-08-25): added the full-dimensional validity classifier.**
   `classify_element_validity<gdim,ideg>` owns the P1 corner-orientation
   prerequisite, normalization by the absolute P1 corner determinant,
   coefficient generation and lower-bound diagnostics, and deterministic
   direct-Jacobian witnesses. It certifies only when the P1 prerequisite passes
   and the finite normalized coefficient minimum is at least `jtol`. Otherwise
   it evaluates the Jacobian at the degree-`gdim*(ideg-1)` simplex lattice and
   reports `Invalid` when a finite normalized value is below `jtol`; without
   such a witness it reports `Uncertified`. A zero or nonfinite P1 normalization
   is also `Uncertified`, with unavailable diagnostics left at their sentinels.
   `test_high_order_validity_classifier` qualifies this routine for P2
   triangles and tetrahedra in both Lagrange and Bezier storage. It covers an
   affine certificate, a curved inversion with a valid P1 skeleton, a positive
   determinant whose whole-element Bernstein bound is inconclusive, explicit
   P1-orientation gating, and degenerate P1 normalization. Production call-site
   migration remains the next Phase 3 step.

4. **Done (2026-08-25): integrated the classifier into Classic smoothing.**
   The interior-coordinate, CAD-edge, and CAD-face branches of
   `smooballdiff` now share one degree-aware region predicate. P1 meshes retain
   `isvalideltP1`; for P2 and higher, every element in the smoothing region
   must return `Certified` from `classify_element_validity`. Consequently,
   both `Invalid` and `Uncertified` Newton trials are rejected under the
   conservative policy.

   The optimizer's stored final candidate is classified again after its
   coordinates and interpolated metric have been installed. This check is
   unconditional: a non-certified candidate returns failure and restores the
   original coordinates, metric, and, for boundary moves, CAD parameters. It
   replaces the previous final P1-only assertion and debug-only legacy
   coefficient check. The same implementation covers P2 triangles in 2D and
   the existing complete tetrahedral edge-shell region in 3D. This step is
   intentionally confined to the Classic (`iflag2 == 0`) smoother; the Luksan
   alternative and cavity smoothing remain later migrations.

   `test_high_order_classic_smoothing` now contains six focused cases. In
   addition to the Phase 2 neighborhood and CAD checks, every successful 2D
   smoothing case independently classifies its affected triangles, the 3D
   case starts and ends with a certified multi-tetrahedron shell and verifies
   actual control-point movement, and a direct low-level 2D case verifies that
   a deliberately non-certified final candidate returns failure without being
   accepted. The six cases pass in debug and release builds, together with the
   validity classifier and contract tests.

5. **Done (2026-08-25): migrated completed cavity elements and final topology.**
   The cavity lifecycle now makes the provisional/completed distinction
   explicit. `reconnect_faccav` and `reconnect_tetcav` create only a P1
   skeleton; their high-order entries are deliberately unset, so their
   immediate orientation checks correctly remain `isvalideltP1`. After
   `correct_cavity` assigns, projects, and metric-interpolates the high-order
   nodes, every new full-dimensional P2 triangle or tetrahedron must be
   `Certified` by `classify_element_validity`. Both `Invalid` and
   `Uncertified` completed elements are added to `lbad` and cause the cavity
   operation to be rejected or enlarged through its existing correction
   policy.

   `check_topo` now dispatches the same classifier at the mesh degree for every
   live full-dimensional element, including P1. Its failure message reports
   the named status, normalized lower bound and coefficient index, and any
   direct witness and sample index. This replaces the separate P1 check plus
   boolean `getsclccoef` check, so the unconditional final high-order audit now
   uses the common three-state result type.

   Embedded P2 surface triangles in 3D remain deliberately separate:
   `correct_cavity` retains their historical compatibility check and does not
   call the full-dimensional classifier. Their oriented polynomial certificate
   is explicitly deferred under the scope decision above.

   Focused tests exercise the completed-element policy directly for certified,
   invalid, and positive-but-uncertified P2 triangles and tetrahedra. Two
   end-to-end cavity insertions additionally split one P2 triangle into three
   certified triangles and one P2 tetrahedron into four certified tetrahedra.
   The final topology test constructs the known positive-but-uncertified P2
   map and verifies that `check_topo` rejects it specifically as
   `Uncertified`. While adding the end-to-end tests, a reversed diagnostic for
   `cav.inewp` was corrected: `1` denotes a newly inserted point and `0` an
   existing point.

## Phase 4 SizeShape execution record

1. **Done (2026-08-25): removed the P1 restriction from value dispatch.**
   `metqua` now sends `SizeShape` value evaluation through the shared objective
   quadrature traversal at the effective geometry degree: `AsDeg::P1` retains
   the corner map, while `AsDeg::Pk` selects the complete mesh degree. The
   runtime dispatch is degree-generic and is qualified first for P2 triangles
   and tetrahedra. `StepDistance` remains explicitly restricted to P1 for its
   separate Phase 5 work, and differentiated `SizeShape` evaluation is
   unchanged.

   `test_high_order_sizeshape_value` constructs curved Lagrange P2 elements
   directly in 2D and 3D for both FE and analytical metric field types. It
   verifies that the public `metqua` entry point exactly matches the selected
   degree-2 shared-integrator specialization, that the P1 compatibility call
   still matches the degree-1 specialization, and that the two geometry
   choices produce observably different values. This is deliberately a
   dispatch contract only. Independent verification of P2 geometry and metric
   sampling, followed by dense-reference integration, remains steps 2 and 3.

2. **Done (2026-08-25): qualified P2 geometry and metric sampling at the
   shared quadrature points.** `test_high_order_sizeshape_value` now visits
   every point of the order-five objective quadrature rule for curved P2
   triangles and tetrahedra, with both FE and analytical metric fields. It
   independently reconstructs the physical point and canonical Jacobian from
   the individual P2 Lagrange basis functions, then checks the common sample's
   geometry and integration measure against those values.

   The metric-degree contract is tested separately from the geometry degree.
   For an FE field, the expected metric is independently formed by
   log-Euclidean P1 interpolation of the corner metrics. Changing every P2
   edge-node metric leaves the sample unchanged, while changing a corner
   metric changes it. For an analytical field, the expected metric is the
   prescribed spatially varying function evaluated at the reconstructed P2
   physical point; changing all stored nodal metrics leaves the sample
   unchanged. Thus this step explicitly preserves `AsDeg::P1` metric
   interpolation on a P2 geometry instead of implicitly promoting the metric
   to `AsDeg::Pk`.

   This is a pointwise sample-contract test, not an independent integral.
   Dense-reference integration remains Phase 4 step 3.

3. **Done (2026-08-25): compared P2 SizeShape values with an independent
   dense reference integral.** The focused high-order value test now builds a
   separate Duffy tensor-product Gauss reference with 64 points on the
   triangle and 512 points on the tetrahedron. Its oracle explicitly evaluates
   the P2 Lagrange basis and gradients, reconstructs the P1 FE or analytical
   metric according to the established contract, forms the regular-simplex
   Jacobian-metric matrix algebra, evaluates the complete SizeShape objective,
   and accumulates the physical measure. It does not call the production
   element evaluator, shared sample preparer, pointwise objective, or objective
   quadrature traversal.

   The test covers curved P2 triangles and tetrahedra with both FE and
   analytical metric fields, using the non-integer objective exponent
   `objective_p = 1.5`. Production orders two through five are reported against
   the dense reference. Order five must improve on order two and agree within
   `2e-5` in the test's absolute-or-relative scale. The observed order-five
   errors are approximately `8e-6` in 2D and `1e-7` in 3D.

4. **Done (2026-08-25): generalized differentiated SizeShape dispatch for an
   active P2 edge control point.** `d_metqua` now mirrors the value path's
   effective-geometry-degree dispatch and sends P2 SizeShape calls through the
   shared differentiated objective traversal. The call remains degree-generic
   at runtime, while StepDistance remains explicitly restricted to its P1 path
   for Phase 5.

   The focused high-order SizeShape test activates the first Lagrange P2 edge
   control point of a curved triangle or tetrahedron for both FE and analytical
   metric fields. It verifies exact agreement between the public `d_metqua`
   result and the selected degree-two differentiated-integrator specialization,
   agreement of the returned value with `metqua` to `3e-14`, and a finite,
   nonzero edge gradient. The established derivative contract keeps the metric
   and physical integration measure frozen at each quadrature sample. This
   step qualifies dispatch and edge activation only; centered finite-difference
   verification of that frozen-sample gradient remains step 5, and independent
   Hessian verification remains step 6.

5. **Done (2026-08-25): verified the active P2 edge gradient with automatic
   differentiation and centered finite differences of the frozen-sample
   functional.** At the unperturbed
   curved element, the test captures every order-five quadrature point, its
   P1 FE or analytical metric tensor, and its complete weight
   `quadrature_weight * theta`. It then perturbs one coordinate of the first P2
   Lagrange edge control point at a time and reevaluates only the pointwise P2
   SizeShape objective with those samples held fixed. This excludes metric and
   physical-measure motion exactly as required by the production derivative
   contract. A separate first-order `SurrealS` evaluation seeds the active
   control point's contribution to the regular Jacobian and automatically
   differentiates the independently assembled SizeShape value. Its result is
   checked directly against the production analytic gradient.

   Centered differences use a step of `2e-6` and are checked in every physical
   coordinate for curved triangles and tetrahedra with both FE and analytical
   metrics at `objective_p = 1.5`. The error limit is
   `1e-7 * (1 + abs(finite_difference_gradient))`. The largest observed errors
   are approximately `8.1e-10` in 2D and `4.2e-11` in 3D. The automatic
   gradient comparison uses the tighter scale `5e-14 * (1 + abs(gradient))`.
   The same AD-of-frozen-value check was backfilled into the legacy P1
   objective matrix for both SizeShape and all currently covered StepDistance
   variants. StepDistance remains P1-only here; its P2 extension is Phase 5.
   No production change was required for this qualification step. Independent
   Hessian verification remains Phase 4 step 6.

6. **Done (2026-08-25): verified the active P2 edge Hessian by automatic
   differentiation and centered finite differences of the frozen analytic
   gradient.** A test-only, scalar-generic implementation reconstructs the
   SizeShape gradient from the regular Jacobian, the frozen metric, and the P2
   regular-reference basis gradient. Seeding the active edge control point in
   `SurrealS` and differentiating that gradient produces an independent
   Hessian without invoking the production Hessian formula. Its primal
   gradient is also checked against the already qualified production gradient.

   The production integrated Hessian is separately reconstructed from the
   pointwise differentiated SizeShape routine at the captured order-five
   samples, confirming the common traversal and weighting. Centered differences
   then perturb each physical coordinate of the active edge point by `2e-6`
   and reevaluate the frozen analytic gradient. The metric, theta, quadrature
   coordinates, and weights remain fixed throughout both oracles.

   Curved triangles and tetrahedra are covered with FE and analytical metric
   fields at `objective_p = 1.5`. Analytic-versus-AD Hessian entries use the
   scale `5e-14 * (1 + abs(H_AD))`; the largest observed absolute discrepancy
   is approximately `1.6e-14`. Analytic-versus-finite-difference entries use
   `1e-7 * (1 + abs(H_FD))`; the largest observed absolute discrepancies are
   approximately `3.1e-8` in 2D and `5.9e-10` in 3D. No production change was
   required. Phase 4 step 7, Lagrange and Bezier representation coverage, is
   next; StepDistance remains P1-only until Phase 5.

7. **Done (2026-08-25): qualified equivalent P2 Lagrange and Bezier geometry
   representations.** The focused SizeShape test converts each curved P2
   triangle or tetrahedron from its Lagrange nodal coordinates to Bezier
   control coefficients through the production basis conversion. At every
   order-five sample it verifies invariant regular Jacobians, frozen P1 FE or
   analytical metrics, and complete objective weights. The integrated value
   and differentiated entry-point value are invariant as well.

   For the first quadratic edge degree of freedom, the conversion satisfies
   `B_edge = 2*L_edge - (L_vertex_0 + L_vertex_1)/2`. The test therefore checks
   the exact derivative transformation laws `g_L = 2*g_B` and
   `H_L = 4*H_B`. Across the 2D/3D FE/analytical matrix, the largest absolute
   sample discrepancy is approximately `3.3e-16`, the largest gradient-law
   discrepancy is `2.2e-15`, and the largest Hessian-law discrepancy is
   `2.9e-14`.

   The Phase 4 gradient and Hessian tests now also execute their independent
   AD and centered-finite-difference oracles with a Bezier active edge control
   point, rather than relying only on cross-representation scaling. All checks
   pass without production changes. Phase 4 step 8, enabling SizeShape P2
   smoothing after the completed value, derivative, and validity
   qualifications, is next. StepDistance remains P1-only until Phase 5.

8. **Done (2026-08-25): enabled SizeShape for production planar P2
   smoothing.** The production smoothing-objective policy now preserves the
   configured default objective for P1 meshes but selects SizeShape for degree
   two in geometric dimension two. With the current StepDistance build this is
   an explicit planar P2 fallback: StepDistance continues to own the P1 path
   and remains unelevated until its Phase 5 value, derivative, validity, and
   smoothing qualification is complete.

   The surrounding legacy optimizer diagnostics and topology-changing
   operations retain their existing objective contracts; this step changes
   only the ball-smoothing objective. The end-to-end high-order smoothing suite
   was migrated from classical Distortion to the production SizeShape
   selection. It checks an active interior P2 edge control point and a complete
   small mesh in 2D, the public production optimizer on an already elevated
   planar P2 mesh, CAD-constrained boundary smoothing, an interior edge shell
   in 3D, and conservative rejection of a noncertified final candidate. All
   accepted regions remain Bernstein-certified after smoothing.

   The isolated 3D SizeShape edge-shell test remains enabled and passes, but a
   full generated 3D optimizer case exposed a later interaction with legacy
   swap/cavity operations that can reach a singular configuration after the
   new smoothing pass. Consequently the 3D production selector deliberately
   retains its previous objective. This follows the planar-first scope: Phase 4
   is complete for production 2D, while 3D production integration is recorded
   as follow-up work rather than being enabled on the strength of an isolated
   test alone. Phase 5 is the separate StepDistance P2 generalization.

## Phase 5 StepDistance execution record

1. **Done (2026-08-25): replaced the degree-one geometry assumptions in the
   StepDistance value path.** `metqua` now dispatches both objective-driven
   value integrals at the effective geometry degree. `AsDeg::P1` continues to
   select the corner map, while `AsDeg::Pk` selects the complete mesh degree;
   metric interpolation remains the independent `asdmet` choice and therefore
   retains the established P1 FE-metric compatibility contract on P2
   geometry. The StepDistance pointwise value routine now makes the same
   degree-aware choice when it reconstructs the physical point and Jacobian,
   removing its internal `eval2<...,1>` and `eval3<...,1>` assumptions.

   `test_high_order_stepdistance_value` constructs curved Lagrange P2
   triangles and tetrahedra for both FE and analytical metric fields. For the
   basic StepDistance variant, it checks exact agreement between the public
   `metqua` entry point and the selected degree-two shared-integrator
   specialization, retains the equivalent P1 compatibility check, and
   requires the P1 and P2 geometry choices to produce observably different
   integrals. A separate frozen-metric pointwise probe perturbs only a P2 edge
   node and verifies that the P2 StepDistance value changes while the P1 value
   is exactly invariant. This specifically guards the pointwise geometry
   dispatch rather than testing only the outer integrator.

   Differentiated StepDistance evaluation remains explicitly P1-only. Removing
   its P1-vertex assumption and qualifying an active P2 edge control point is
   Phase 5 step 2. The collapse barrier, ShapeVolume, and
   CavityTargetAverage variants are not qualified by this step and retain
   their separate roadmap entries.

2. **Done (2026-08-25): removed the P1-vertex restriction from differentiated
   StepDistance dispatch.** `d_metqua` now follows the same effective geometry
   degree as the value path for both objective-driven policies. The pointwise
   StepDistance derivative reconstructs the geometry at that degree and
   accepts any local control-point index in the corresponding simplex. Its
   active basis gradient is evaluated as either a Lagrange or Bezier basis
   function and transformed into the regular reference frame.

   The regular-reference basis-gradient helper previously owned only by the
   shared quadrature traversal now lives in `objective_basis.hxx`. Both the
   traversal's integration-level terms and the pointwise StepDistance
   derivative use this same implementation, preventing the P2 objective and
   collapse-barrier contributions from silently differentiating different
   basis functions.

   The focused StepDistance test activates the first Lagrange P2 edge control
   point of a curved triangle or tetrahedron for FE and analytical metrics. It
   checks exact agreement between the public `d_metqua` result and the selected
   degree-two differentiated-integrator specialization, agreement of the
   returned value with `metqua` to `3e-14`, and a finite, nonzero edge
   gradient. This step qualifies differentiated dispatch and edge activation
   only. Independent automatic-differentiation and centered-finite-difference
   checks belong to the basic-variant validation in Phase 5 step 3; the later
   StepDistance variants retain their separate roadmap entries.

3. **Done (2026-08-25): implemented and independently validated the basic P2
   StepDistance variant.** The focused test now uses a spatially varying
   positive-definite metric and curved P2 triangles and tetrahedra. For both
   FE and analytical metric fields, an independent eight-point-per-coordinate
   Duffy quadrature reconstructs the explicit P2 Lagrange geometry, the
   regular Jacobian, the spectral StepDistance value, and the geometric
   integration weight without calling the production sample preparation or
   objective integration. The FE oracle explicitly retains log-Euclidean P1
   corner-metric interpolation. Configured objective quadrature orders two
   through five are compared with this dense reference; order five improves
   on order two in every case and remains within the `2e-5` scaled tolerance.

   The differentiated checks capture order-five samples and freeze their
   metric, quadrature point, and complete quadrature-times-theta weight. For
   the first P2 edge control point, in both Lagrange and Bezier bases, the
   production analytic gradient is compared with automatic differentiation
   of an independently assembled scalar StepDistance value and with centered
   differences of the frozen value. The analytic-versus-AD tolerance is
   `5e-14 * (1 + abs(gradient))`; the largest observed absolute discrepancy is
   approximately `1.33e-15`. The largest centered-difference gradient error is
   approximately `5.86e-11`.

   A separate scalar-generic oracle reconstructs the analytic spectral
   gradient through `A^{-1} log(A)` and differentiates it automatically. This
   checks the production Hessian independently of its AD-of-production-
   gradient implementation. Centered differences of the frozen analytic
   gradient provide the second Hessian oracle. The largest production-versus-
   AD Hessian discrepancy is approximately `1.78e-15`, under the same `5e-14`
   scaled tolerance, and the largest centered-difference Hessian error is
   approximately `3.67e-10`, under the `1e-7` scaled tolerance.

   This qualification covers the basic formula only, with the collapse
   barrier disabled and `ShapeVolume` and `CavityTargetAverage` false. Their
   separate Phase 5 steps remain pending, as do basic-variant validity
   rejection and production smoothing enablement.

4. **Done (2026-08-25): qualified the P2 metric-volume collapse barrier.**
   The barrier already entered through the shared objective-quadrature
   policies after the degree-aware value and derivative dispatch work. This
   step verifies that integration contract directly: the active-barrier
   result is subtracted from the identical `beta = 0` basic StepDistance
   result, isolating the additive barrier in both the value-only and
   differentiated public entry points. The independent reconstruction weights
   the barrier by the quadrature weight alone, without the geometric `theta`
   factor, and freezes the metric exactly as required by the optimizer.

   A test-owned explicit quadratic polynomial construction evaluates the
   active edge basis gradient for both Lagrange and Bezier representations and
   transforms it to the regular reference frame. It agrees with the shared
   production degree-two helper at the `5e-14` scaled tolerance for every
   order-five quadrature point in curved triangles and tetrahedra. The same
   independent gradient drives an AD-of-barrier-value oracle and a separately
   coded analytic barrier-gradient oracle. Differentiating the latter provides
   the independent AD-of-gradient Hessian comparison; centered differences of
   the frozen barrier value and gradient provide the numerical comparisons.

   The matrix covers 2D and 3D, FE and analytical metrics, and Lagrange and
   Bezier geometry. The largest isolated production-versus-AD gradient error
   is approximately `4.44e-16`, and the largest production-versus-independent-
   AD-of-gradient Hessian error is approximately `1.78e-15`. The largest
   centered-difference errors are approximately `1.82e-10` for the gradient
   and `3.98e-9` for the Hessian, within their `1e-7` scaled tolerances.

   ShapeVolume remains disabled and is Phase 5 step 5.
