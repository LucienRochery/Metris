# High-Order Mesh Generation Implementation Record

## Scope

The primary target is P2 in full-dimensional planar meshes. Full-dimensional
3D support is included when the extension follows the same implementation
without introducing a separate policy. P3 remains a design constraint, not a
current acceptance requirement.

The work started from the historical Classic P2 workflow, then generalized the
common objective and adaptation machinery to complete P2 geometry. This file is
the concise implementation record. Mathematical derivations and detailed test
formulations live in the complementary local high-order notes.

## Phase overview

| Phase | Status | Summary |
| --- | --- | --- |
| 1. Trustworthy baseline | Complete | Re-audited historical blockers, made regressions hard-failing, corrected the top-level lifecycle and metric-cost integration, and established deterministic 2D/3D elevation and native-P2 baselines. |
| 2. Classic P2 smoothing | Complete | Restored edge-control-point smoothing in planar interiors and boundaries, retained the tetrahedral edge-shell path, and qualified CAD-curve motion. |
| 3. High-order validity | Complete for full-dimensional P2 | Introduced one three-state Bernstein-based validity contract and applied it to smoothing, completed cavities, and final topology. Embedded surfaces and Bernstein subdivision are deferred. |
| 4. SizeShape P2 | Complete | Generalized value, gradient, and Hessian evaluation to P2; qualified quadrature, metric sampling, Lagrange/Bezier equivalence, and smoothing. |
| 5. StepDistance P2 | Complete | Generalized all StepDistance variants to P2 in 2D and 3D, qualified their derivatives and validity behavior, and enabled the configured production objective. |
| 6. CAD and boundaries | Complete for the planar-P2 scope | Made elevated nodes exact CAD images, generalized surface localization, audited accepted boundary moves, and added a planar quadratic boundary-intersection guard. |
| 7. Adaptation operations | Complete | Generalized insertion, growth, collapse, swaps, cavity smoothing, dependency reactivation, and completed-cavity conformity to P2. |
| 8. End-to-end qualification | Production follow-up | Focused tests are complete. Initial production runs are acceptable; broader production qualification and investigation of unusual StepDistance behavior are postponed. |

## Cross-cutting contracts

These contracts apply throughout the implementation and should not be changed
implicitly by later work.

1. **Geometry degree.** `AsDeg::P1` selects the corner map. `AsDeg::Pk`
   selects the complete mesh geometry at `msh.curdeg`.
2. **FE metric degree.** P2 geometry currently retains `AsDeg::P1` FE-metric
   interpolation. Only corner metric coefficients participate. An analytical
   metric is evaluated directly at each physical high-order sample. P2 FE
   metric interpolation is separate future work.
3. **Objective quadrature.** SizeShape, StepDistance, and final metric-cost
   integration use `objective_quadrature_order`. Order `0`, the
   vertex-barycenter rule, is the default in every dimension.
4. **Derivative convention.** Objective gradients and Hessians differentiate
   the geometry while holding the sampled metric and physical integration
   measure `theta` fixed. Tests reproduce this frozen-sample functional.
5. **Validity.** Quadrature samples do not prove element validity. A completed
   full-dimensional P2 element is accepted only when the common Bernstein
   classifier returns `Certified`.
6. **Cavity geometry.** Existing edges retain their stored polynomial trace;
   exact children of a split edge restrict the parent trace; genuinely new
   non-CAD edges are affine P2; CAD-owned nodes are exact CAD images.
7. **CAD ownership.** Every geometric node classified on a CAD entity is
   represented and optimized in that entity's parameter space. Lower-
   dimensional ownership takes precedence.
8. **Transactional acceptance.** Rejected smoothing and cavity candidates
   restore coordinates, metrics, CAD parameters, geometry provenance, and
   objective state.

## Detailed implementation record

### Phase 1: Trustworthy baseline

1. **Done — Source audit (2026-08-20).** The audit confirmed the boundary
   edge-shell bug, P1-only objective paths, P1-only final metric cost,
   incomplete CAD-face projection, mixed P1/P2 cavity acceptance, swallowed
   regression exceptions, duplicate optimization, and a debug-only final
   high-order check. Each issue was assigned to a later phase.

2. **Done — Hard-failing regressions (2026-08-20).** An unexpected
   `MetrisExcept` now produces an immediate Boost failure while its diagnostic
   remains recorded in `runs.json`.

3. **Done — Single top-level optimization pass (2026-08-20).** `runMetris()`
   now calls `optimMesh()` once, after adaptation, degree elevation, and
   optional curving. The optimization internal to curvature mode 3 remains
   part of that curving algorithm.

4. **Done — Degree-aware final metric cost (2026-08-20).** Final cost uses the
   full geometry degree and the objective quadrature rule. It evaluates the
   metric in exponential SPD space and restores the caller's storage space.
   The P1 FE-metric interpolation contract remains explicit.

5. **Done — Unconditional final validity check (2026-08-20).** Every
   successful full-dimensional high-order run checks topology and high-order
   validity before output, independently of `dbgfull`.

6. **Done — Four deterministic baselines (2026-08-20).**
   `test_high_order_phase1_baselines` covers P1-to-P2 elevation and native-P2
   input for a two-triangle planar mesh and a one-tetrahedron volume mesh. The
   runs verify degree, entity counts, final metric cost, and final validity in
   Release and Debug builds, with adaptation and optimization disabled.

### Phase 2: Classic P2 smoothing

1. **Done — Edge-control-point regions (2026-08-20).**
   `buildEdgeControlPointSmoothingRegion` returns two triangles for an interior
   planar edge, one triangle for a boundary edge, and the complete tetrahedral
   shell for a volume edge.

2. **Done — Classic P2 traversal (2026-08-20).** A stale debug-only P1
   assertion was removed from the generic Distortion traversal. The
   objective-adaptation dispatcher was also made to honor `adp_niter == 0`.

3. **Done — Planar control-point smoothing (2026-08-20).** Focused tests
   confirm that an interior edge node is shared by both triangles, receives
   the correct region, moves under Classic smoothing, and does not increase
   the regional objective. A boundary edge node receives exactly one triangle.

4. **Done — CAD-curve smoothing (2026-08-20).** An elevated boundary edge node
   is optimized in curve-parameter space. Its final physical coordinate is
   checked independently with `EG_evaluate`.

5. **Done — Complete planar and focused volume coverage (2026-08-20).** A full
   small planar P2 mesh completes a Classic pass. A genuine interior edge node
   in an elevated tetrahedral mesh receives its complete shell and is visited
   without increasing the shell objective.

### Phase 3: High-order validity

For a degree-`p` simplex of topological dimension `d`,
`deg(det DF) <= d(p-1)`. If

```text
det DF(xi) = sum_alpha c_alpha B_alpha(xi),
```

then the Bernstein partition of unity gives
`det DF(xi) >= min_alpha c_alpha`. Positive coefficients are therefore a
sufficient, but not necessary, certificate.

1. **Done — Three-state result (2026-08-25).** `ElementValidityResult`
   distinguishes `Certified`, `Invalid`, and `Uncertified`. It stores the
   normalized coefficient lower bound and index and, when available, a direct
   failing Jacobian witness and sample index. Conservative acceptance admits
   only `Certified`.

2. **Done — Coefficient qualification (2026-08-25).** Deterministic affine and
   curved P2 triangles and tetrahedra verify coefficient generation against
   direct off-lattice Jacobian evaluation in both Lagrange and Bezier storage.

3. **Done — Common full-dimensional classifier (2026-08-25).**
   `classify_element_validity<gdim,ideg>` owns the P1 orientation prerequisite,
   normalization, Bernstein lower bound, deterministic witness sampling, and
   diagnostics. Tests cover all three outcomes and degenerate normalization in
   2D and 3D.

4. **Done — Smoothing integration (2026-08-25).** Every completed high-order
   smoothing trial and accepted candidate must be certified. Rejection restores
   coordinates, metric data, and CAD parameters. The same predicate covers
   planar regions and tetrahedral edge shells.

5. **Done — Cavity and topology integration (2026-08-25).** Provisional cavity
   elements retain their P1 orientation prerequisite. After `correct_cavity`
   assigns all high-order nodes, every new full-dimensional element must pass
   the common classifier. `check_topo` uses the same result and reports its
   status, coefficient bound, and witness diagnostics.

6. **Done — Final-check migration (2026-08-25).** The unconditional final
   high-order check established in Phase 1 now uses the common three-state
   classifier.

7. **Deferred — Embedded surface certificate.** An embedded curved triangle
   requires a consistently oriented polynomial such as
   `n0 dot (dF/dxi1 cross dF/dxi2)`. Full-dimensional boundary faces already
   inherit local regularity from an incident certified tetrahedron. Resume this
   work with embedded 2D-in-3D meshes.

8. **Deferred — Bernstein subdivision.** Subdivision may refine an
   inconclusive whole-element lower bound and reduce conservative false
   rejection. It is not required for the current planar-P2 scope.

### Phase 4: SizeShape P2

1. **Done — Degree-aware value dispatch (2026-08-25).** `metqua` sends
   SizeShape through the shared objective traversal at the effective geometry
   degree while preserving explicit P1 compatibility.

2. **Done — Shared sample contract (2026-08-25).** Curved P2 triangles and
   tetrahedra verify physical points, regular Jacobians, integration weights,
   P1 log-Euclidean FE metrics, and analytical metrics at every quadrature
   sample.

3. **Done — Independent value reference (2026-08-25).** A separate dense
   Duffy quadrature reconstructs P2 geometry and the metric without calling the
   production evaluator. Orders two through five converge to this reference in
   2D and 3D for FE and analytical metrics.

4. **Done — Edge-control derivative dispatch (2026-08-25).** `d_metqua`
   accepts an active P2 edge node and dispatches differentiated SizeShape at
   the effective geometry degree.

5. **Done — Gradient qualification (2026-08-25).** The analytic gradient is
   compared with AD of an independently assembled frozen-sample value and with
   centered differences of that value. The AD comparison uses
   `5e-14 * (1 + abs(g_AD))`; the largest observed error is at machine
   precision. Equivalent AD checks were added to the legacy P1 SizeShape and
   StepDistance tests.

6. **Done — Hessian qualification (2026-08-25).** The production Hessian is
   compared with AD of an independently assembled analytic gradient and with
   centered differences of that gradient under the same frozen-sample
   convention.

7. **Done — Lagrange/Bezier equivalence (2026-08-25).** Values and samples are
   invariant under basis conversion. For a quadratic edge degree of freedom,
   the tests verify `g_L = 2 g_B` and `H_L = 4 H_B`, including independent AD
   and finite-difference checks in both bases.

8. **Done — Production smoothing integration (2026-08-25).** Planar P2
   SizeShape smoothing was enabled after value, derivative, and validity
   qualification. Tests cover an interior edge node, a complete small mesh,
   CAD-constrained motion, a tetrahedral edge shell, and conservative rejection
   of a noncertified candidate. The temporary production selection was later
   superseded by the completed StepDistance work in Phase 5.

### Phase 5: StepDistance P2

1. **Done — Degree-aware value path (2026-08-25).** Hard-coded degree-one
   geometry evaluation was removed from the outer and pointwise StepDistance
   value paths. P1 and Pk geometry remain explicit choices.

2. **Done — Degree-aware differentiated path (2026-08-25).** `d_metqua` and
   the pointwise derivative accept any local P2 control point. The shared
   regular-reference basis-gradient helper is used by both integration and
   pointwise policies.

3. **Done — Basic variant (2026-08-25).** Independent dense integration,
   AD-of-value gradients, AD-of-gradient Hessians, and centered differences
   qualify curved P2 triangles and tetrahedra for FE and analytical metrics in
   Lagrange and Bezier bases. The gradient comparison is at machine precision.

4. **Done — Collapse barrier (2026-08-25).** The isolated additive barrier is
   verified with its reference-measure weighting, frozen metric, independent
   active-basis construction, AD, and finite differences. The matrix covers
   both dimensions, metric types, and geometry bases.

5. **Done — ShapeVolume (2026-08-25).** Dense value references and independent
   derivative oracles qualify the centered log-shape coordinates and injective
   volume coordinate. Singular trials return the finite rejection sentinel
   with zero derivatives. The separate collapse barrier is disabled by this
   policy.

6. **Done — CavityTargetAverage (2026-08-25).** The element objective is the
   unit-weight reference-simplex average and the regional objective is the
   arithmetic mean by live-element count. Dense references, independent
   derivatives, count-changing replacement, and mesh-wide state updates verify
   this distinct aggregation contract.

7. **Done — Separate dimensional gates (2026-08-25).** The focused test is
   divided into matching 2D and 3D suites. Each covers value and derivative
   dispatch, Basic, CollapseBarrier, ShapeVolume, and CavityTargetAverage for
   FE and analytical metrics and both geometry bases.

8. **Done — Smoothing, validity, and production selection (2026-08-25).** All
   four variants move an interior P2 edge node without increasing the
   configured objective, preserve conservative validity, and restore rejected
   noncertified candidates in both 2D and 3D.

   Some valid StepDistance regions have an indefinite Hessian. High-order
   StepDistance smoothing therefore opts into a bounded steepest-descent
   fallback when Cholesky fails or the Newton direction is not a descent
   direction. No analytic StepDistance Hessian is claimed: the production
   Hessian is obtained by automatic differentiation of the analytic gradient.
   The fallback does not imply that SizeShape or P1 Hessians are always
   positive definite.

   Objective-driven smoothing now relies only on its configured regional or
   global objective plus conservative validity. The legacy worst-element guard
   remains exclusive to Distortion and Unit paths. In a StepDistance build,
   `productionSmoothingObjective` selects StepDistance at P1 and P2 in both 2D
   and 3D.

### Phase 6: CAD and boundary support

1. **Done — Planar curved-curve creation and smoothing (2026-08-25).** Degree
   elevation first converts geometry to Lagrange form, then sets every
   curve-owned edge-interior coordinate to the CAD image of its elevated
   parameter: `x = C(t)`. Smoothing uses the CAD tangent and curvature and
   constrains `t` to the topological edge interval. SizeShape and StepDistance
   may use the qualified bounded descent fallback in parameter space.

2. **Done — CAD-surface placement and localization (2026-08-25).** P1 vertex
   parameters are elevated as scalar Lagrange fields. A face-owned P2 edge node
   receives the interpolated `(u,v)`, is placed at `S(u,v)`, and is localized
   in the background metric. Curve ownership takes precedence over surface
   ownership.

3. **Deferred — Embedded surface orientation (2026-08-25).** A separate
   oriented surface certificate is unnecessary for a boundary face of a
   certified tetrahedron but remains necessary for embedded 2D-in-3D meshes.

4. **Done — Reusable CAD regressions (2026-08-25).** Self-contained EGADS
   regression cases construct a planar triangle bounded by one circular arc
   and two lines, and a two-triangle mesh on a spherical patch. They cover
   elevation and native-P2 re-ingestion without external binary data.

5. **Done — Accepted boundary-move audit (2026-08-25).** Before a move is
   counted, every CAD parameter must be finite, each CAD evaluation must match
   the physical coordinate, and every affected full-dimensional element must
   pass the degree-aware validity policy. Curve-owned points in 3D synchronize
   secondary surface records with `EG_getEdgeUV`. Rejection or exception
   restores the complete boundary-record chain.

6. **Done — Planar quadratic boundary intersections (2026-08-25).** Every
   changed boundary edge is compared globally with the other live boundary
   edges. Shared topological endpoints are allowed; crossing, overlap, or any
   other contact is rejected transactionally. The test converts the P2
   Lagrange midpoint to a quadratic Bezier control point and uses recursive de
   Casteljau subdivision with convex-hull bounds. This detects global contacts
   that local positive-Jacobian certificates cannot detect.

### Phase 7: Adaptation operations

1. **Done — Insertion and cavity growth (2026-08-26).** P1 topology and
   orientation remain prerequisites, but initial validity growth, local
   objective growth, and final acceptance use completed P2 reconnections.
   `correct_cavity` assigns all high-order geometry and CAD data before the
   final validity and objective callback; rejection rolls back before topology
   commit.

   Growth compares only the affected local patch. Configuration A keeps the
   outside neighbor and reconnects its incident cavity elements. Configuration
   B absorbs that neighbor and reconnects the enlarged local patch. The current
   configuration is already P2-valid by construction; an uncertified enlarged
   configuration is rejected, otherwise the exact local objective replacement
   decides growth. Unaffected cavity elements are neither rebuilt nor
   integrated.

   Existing edges retain their coefficients. Non-CAD children of a split edge
   are exact parent restrictions, using Lagrange evaluation or de Casteljau as
   appropriate. CAD-owned children use interpolated parameters and exact CAD
   evaluation. Genuinely new edges are affine P2. SizeShape uses additive local
   integrals; CavityTargetAverage applies the same local replacement to its
   maintained global numerator and element count.

2. **Done — Collapse (2026-08-26).** Collapse reuses the completed cavity
   machinery over the shell of the removed edge plus the endpoint balls. No
   split child survives, so split-edge provenance is empty. Existing cavity-
   boundary edges retain their geometry; new non-CAD replacement spokes are
   affine P2; CAD-owned coefficients remain exact CAD images. Final completed
   validity and objective acceptance precede commit.

3. **Done — Swaps (2026-08-26).** Planar edge flips and tetrahedral face/edge
   swaps use P1 topology proposals followed by completed-P2 reconstruction,
   conservative validity, and the configured completed-P2 objective. Surviving
   boundary edges retain their coefficients and genuinely new edges are affine
   P2. Rejected candidates restore topology and mesh-wide objective state.

4. **Done — Length smoothing policy (2026-08-26).** Legacy metric-length
   smoothing remains Classic-P1-only. Global P2 length smoothing exits before
   point traversal, and P2 cavity construction skips the legacy length move.
   Objective-driven insertion no longer has a length-based fallback or its
   obsolete arguments and statistics. Objective-driven ball and cavity
   smoothing remain enabled.

5. **Done — Cavity-targeted smoothing (2026-08-26).** Cavity smoothing moves
   the provisional apex and rebuilds every incident non-CAD P2 spoke as
   affine—the selected released-affine-spoke policy. Existing cavity-boundary
   edges remain fixed. CAD ownership overrides the affine rule and keeps nodes
   exact on their CAD entities.

   The derivative includes the dependent spoke motion. For an affine midpoint
   `C_0j(x) = (x + X_j)/2`,
   `dF/dx = partial F/partial x + 1/2 sum_j partial F/partial C_0j`. CAD-curve
   trials apply the corresponding parameter-space chain rule. The Hessian is a
   centered, validity-aware finite difference of this exact chain gradient.
   Acceptance is transactional: success releases non-CAD split provenance;
   rejection restores it exactly.

6. **Done — Dependency reactivation (2026-08-26).** After an accepted move,
   every geometric node of every affected element is reactivated. P1 therefore
   retains corner-only behavior, while P2 reactivates vertices and shared edge
   nodes. Targeted smoothing scans the complete P2 node set and refreshes every
   affected `BadEntHandler` value with the configured full-P2 objective.
   Boundary nodes remain constrained to CAD parameter space; boundary records
   without a CAD model do not become free physical variables.

7. **Done — Completed-cavity conformity (2026-08-26).** Reconnection tables
   already shared one global P2 control point per topological edge in planar
   and tetrahedral cavities. An independent oracle now builds both
   `edge -> control point` and `control point -> edge` maps over the live mesh.
   It requires identical edges to share one live identifier, distinct edges not
   to alias, and every newly completed element to be `Certified`.

   The audit runs before and after insertion commit and after accepted planar
   collapse and swap. It covers new-new and inherited-new interfaces in 2D and
   3D while retaining adjacent P1 compatibility tests.

### Phase 8: End-to-end qualification

Focused implementation qualification is complete. Initial production tests
with the circular analytical metric are acceptable. Unusual StepDistance
behavior observed in those runs is postponed as separate follow-up work rather
than silently folded into the completed Phases 1--7.

The remaining production qualification is intentionally small rather than a
Cartesian product of every already-qualified kernel option:

1. **Follow-up — Planar native-P2 adaptation.** Exercise the complete driver with
   real insertion, collapse, swap, and smoothing activity, then reload the
   output and independently verify degree, topology, edge-control ownership,
   CAD consistency, boundary intersections, and element validity.
2. **Follow-up — P1-to-P2 lifecycle.** Keep this distinct from native-P2
   adaptation: the standard driver adapts the P1 mesh before elevation and
   performs high-order optimization afterward.
3. **Follow-up — Metric and quadrature coverage.** Use the circular analytical
   metric and an FE field sampled from it. Retain the P1 FE interpolation
   contract. Exercise the default objective quadrature and one higher order.
4. **Follow-up — Full-dimensional 3D smoke test.** Qualify a small native-P2
   tetrahedral adaptation without making it a blocker for planar readiness.
5. **Follow-up — Build and output audit.** Run the selected workflows in Debug
   and Release and record exit status, operation counts, objective behavior,
   entity counts, minimum validity bound, and independent final verification.

## Qualification map

| Phase | Principal focused executables |
| --- | --- |
| 1 | `test_high_order_phase1_baselines` |
| 2 | `test_high_order_classic_smoothing` |
| 3 | `test_high_order_validity_contract`, `test_high_order_p2_jacobian_coefficients`, `test_high_order_validity_classifier`, `test_high_order_phase1_baselines` |
| 4 | `test_high_order_sizeshape_value`, `test_high_order_classic_smoothing` |
| 5 | `test_high_order_stepdistance_value`, `test_high_order_stepdistance_smoothing`, `test_high_order_classic_smoothing` |
| 6 | `test_high_order_curved_cad_smoothing`, `test_high_order_curved_cad_surface_projection`, `test_high_order_curved_boundary_intersection` |
| 7 | `test_high_order_phase1_baselines`, `test_high_order_cavity_smoothing`, `test_high_order_smoothing_reactivation`, `test_high_order_curved_cad_smoothing` |

The focused high-order suites were run directly in both Release and Debug
configurations as their phases were completed. Phase 8 adds complete-driver
production runs rather than duplicating these kernel and operation matrices.

## Explicitly deferred work

- P3 and higher geometry.
- P2 interpolation of FE metric coefficients (`AsDeg::Pk` metrics).
- Embedded 2D meshes in 3D and their oriented surface certificate.
- Bernstein subdivision for inconclusive validity bounds.
- General global self-intersection checks for curved surfaces in 3D.
- Higher-degree curved-boundary intersection checks.
