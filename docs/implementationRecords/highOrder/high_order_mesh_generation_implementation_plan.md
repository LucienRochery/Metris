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

1. **Done (2026-08-25): complete 2D curved CAD-boundary control-point
   smoothing.** See the Phase 6 execution record below.
2. **Done (2026-08-25): implement P2 triangular-face projection for 3D
   localization.** See the Phase 6 execution record below.
3. **Deferred (2026-08-25): test the polynomial orientation certificate for
   embedded P2 surface triangles.** Full-dimensional 3D meshes obtain local
   boundary-face regularity from the incident tetrahedron certificate. Resume
   this item when embedded 2D meshes in 3D are brought into scope.
4. **Done (2026-08-25): add small CAD-curve and CAD-surface regression
   meshes.** See the Phase 6 execution record below.
5. **Done (2026-08-25): verify that CAD parameters, projected coordinates, and
   high-order validity remain consistent after every accepted move.** See the
   Phase 6 execution record below.
6. **Done (2026-08-25): reject planar P2 curved-boundary intersections that
   local Jacobian positivity cannot detect.** See the Phase 6 execution record
   below. Phase 6 is complete for the agreed planar-P2 priority and the 3D
   extensions that followed naturally; embedded surfaces, global 3D surface
   intersection, and higher-degree checks remain explicit future work.

## Phase 7: Generalize adaptation operations

Proceed one operation at a time, retaining a P1 compatibility test beside each
new P2 test.

1. **Done (2026-08-26): make insertion acceptance use both the completed P2
   objective and P2 validity.** See the Phase 7 execution record below.
2. **Done (2026-08-26): apply the completed-P2 validity and objective contract
   to collapse.** See the Phase 7 execution record below.
3. **Done (2026-08-26): generalize face and edge swaps to completed-P2
   validity, geometry, and objective acceptance.** See the execution record.
4. **Done (2026-08-26): retain length smoothing as an explicitly Classic-P1-only
   heuristic and remove it completely from objective-driven insertion.** See the
   execution record below.
5. **Done (2026-08-26): generalize cavity-targeted smoothing to completed P2
   geometry using the released-affine-spoke policy.** See the execution record
   below.
6. **Done (2026-08-26): reactivate every geometric node and refresh every
   affected element after a successful smoothing move.** See the execution
   record below.
7. **Done (2026-08-26): confirm that newly created cavity elements share one
   global P2 control point per topological edge and pass the common validity
   certificate.** See the execution record below. Phase 7 is complete.

### Phase 7 execution record

1. **Done (2026-08-26): completed-geometry P2 insertion acceptance.** The
   provisional P1 skeleton remains useful for topology construction and as an
   orientation prerequisite, but every full-dimensional validity-growth and
   objective-growth decision is now made on the corresponding completed P2
   local reconnection. A final cavity callback also runs after
   `correct_cavity` has created every
   high-order coefficient, projected CAD-owned nodes, interpolated their
   metrics, and applied the common conservative P2 validity classifier. Only
   after that callback accepts does `update_cavity` commit the topology; a
   rejection returns `CAV_ERR_OBJECTIVE`, rolls all cavity-appended entities
   and control points back, and lets the insertion caller discard the proposed
   vertex through its ordinary failure cleanup.

   Concisely, relative to the **legacy Classic P2 insertion** (not merely to
   the objective-driven path), the completed implementation changes the
   following contract:

   | Concern | Legacy Classic P2 insertion | Completed Phase 7.1 insertion |
   | --- | --- | --- |
   | Early cavity decisions | P1 reconstruction, orientation, length, and quality checks drove the operation. | P1 topology and orientation remain prerequisites. Initial validity growth and later objective growth both construct the completed P2 local reconnection before deciding. |
   | Existing cavity-boundary edges | Their stored high-order coefficients were already reused. | This behavior is retained explicitly, so the old cavity boundary is unchanged. |
   | Children of a split curved edge | Newly generated non-CAD child coefficients followed the generic affine/chord initialization; CAD correction could subsequently place CAD-owned nodes on the CAD model. | Non-CAD children are exact polynomial restrictions of the parent; CAD-owned children are evaluated from interpolated CAD parameters. |
   | Genuinely new edges | Constructed as straight edges. | Still straight, but explicitly represented consistently as affine P2 edges. |
   | Final acceptance | `correct_cavity` completed the high-order mesh and applied its high-order checks, but Classic acceptance had no completed-P2 objective comparison. | The common conservative P2 validity result is followed by a completed-P2 SizeShape or StepDistance comparison before topology commit. |
   | Failure and quality state | No distinct post-completion objective rejection existed; quality bookkeeping could retain the earlier P1 geometry contract. | Callback rejection rolls the cavity back with `CAV_ERR_OBJECTIVE`, and objective/handler state consistently uses P2 geometry with frozen P1 metric interpolation. |

   Thus the new path does not discard the useful legacy P2 reconstruction: it
   preserves inherited coefficients and CAD correction, while making split
   geometry exact and adding a transactional completed-geometry validity and
   objective gate.

   Initial validity growth now establishes the invariant needed by the later
   objective-growth loop. Starting from the edge shell, each boundary cone is
   completed at P2 and passed through the common conservative classifier; a
   merely `Uncertified` cone is rejected just like an `Invalid` cone, so its
   outside neighbor is absorbed through the existing growth mechanism. The
   P1 cone check is retained only for P1 meshes and for the explicitly deferred
   embedded-surface compatibility path. Consequently, once objective growth
   starts, the current reconnected full-dimensional cavity is P2-valid by
   construction.

   Each subsequent neighbor decision is strictly local. Let the candidate
   outside neighbor be `N`, and let `I(N)` be only the current cavity elements
   incident to it. Configuration A keeps `N` unchanged and reconnects
   `I(N)`; configuration B absorbs `N` and reconnects the enlarged local patch
   `I(N) union {N}`. Both reconnections use completed P2 elements, but no
   unaffected element elsewhere in the cavity is rebuilt or integrated.
   Configuration A is asserted valid by the preceding invariant. If B is not
   conservatively certified, the neighbor is left outside; otherwise the two
   local objectives are compared. For additive SizeShape this is the direct
   local comparison. For mesh-average StepDistance the same local numerator
   and element-count replacement is applied to the maintained global state,
   avoiding a full-cavity recomputation while remaining algebraically exact.

   The scratch completion used by both growth stages follows the same edge
   geometry contract as final reconnection: inherited edges copy their stored
   coefficient; a child of the insertion edge is the exact restriction of the
   parent polynomial; a CAD-owned new coefficient is evaluated on its common
   curve or surface from interpolated CAD parameters; and a genuinely new
   non-CAD edge is affine but still represented as P2. Geometry is evaluated
   with `AsDeg::Pk`, the frozen FE metric contract remains `AsDeg::P1`, and
   the configured objective quadrature is used.

   The final objective uses `AsDeg::Pk` for geometry and retains the established
   frozen P1 FE-metric interpolation contract (`AsDeg::P1`). SizeShape compares
   the additive old and new cavity integrals through the configured
   `BadEntHandler` improvement threshold. StepDistance uses the same local
   completed-P2 values and, for CavityTargetAverage, replaces the old regional
   numerator and element count in an exact mesh-wide P2 objective state before
   deciding. The handler is now seeded and updated with that same P2-geometry,
   P1-metric contract, so later operations do not inherit stale P1 qualities.

   Reconnection now also satisfies the geometric prerequisite needed for a
   meaningful comparison. Existing cavity-boundary edges reuse their original
   high-order coefficients. When a non-CAD P2 edge is split, the insertion
   records the parent edge and split barycentric coordinate; both child-edge
   coefficients are then the exact polynomial restrictions of the parent. In
   Lagrange form they are evaluations at the child midpoint images, while in
   Bezier form they are obtained by the corresponding de Casteljau step.
   Genuinely new edges remain straight but are represented as P2. CAD-owned
   children deliberately override discrete restriction: their interpolated
   curve/surface parameters are evaluated by `correct_cavity`, so every new
   geometric node lies exactly on its owning CAD entity.

   The insertion-specific sharing checks above are a prerequisite of this
   step, not an early completion of Phase 7 item 7. Item 7 remains the later
   cross-operation audit after collapse and swap reconnection have also been
   generalized. At the completion of insertion item 1, `CAVSMOOTHING` was
   therefore still bypassed for P2. Phase 7 item 5 below removes that temporary
   bypass using the released-affine-spoke completed-P2 policy.

   The self-contained regression verifies that the final callback observes
   completed and certified P2 candidates for both a planar triangle and a
   tetrahedron, and that callback rejection restores the pre-cavity counts and
   topology. A curved two-triangle planar case verifies both exact child
   restrictions, reuse of all four inherited outer-edge coefficients, affine
   P2 construction of the two genuinely new spokes, and certification of all
   four children. A paired P1 cavity insertion verifies the unchanged optional-
   hook path. A direct Bezier-basis regression independently verifies that both
   restricted child coefficients are the expected first-level de Casteljau
   points. Two additional strip regressions qualify the growth stages
   themselves. One independently reconstructs both local configurations,
   proves that the reported objectives contain exactly the unchanged neighbor
   plus the single incident cavity cone, and verifies that the completed-P2
   values differ from their P1 surrogates. The other constructs a boundary
   cone whose P1 map is valid while its inherited curved P2 edge makes the
   completed element non-certified, and verifies that initial validity growth
   absorbs precisely that neighbor and leaves every final boundary cone
   certified. Release and Debug builds pass all eleven insertion-oriented
   Phase 1/7 cases; the seven-case legacy cavity suite also passes in Release.
   The older general cavity executable retains unrelated Debug-only
   mesh-construction assertions and is therefore claimed here only in Release.

2. **Done (2026-08-26): completed-geometry P2 collapse.** Collapse already
   enters the insertion cavity machinery as a special replacement: the seed
   starts with the shell of the selected edge, then appends the complete balls
   of both endpoints. The P2 validity-growth and local objective-growth work
   from item 1 therefore applies without a separate algorithm. The remaining
   insertion-only restriction was the final objective callback; it is now a
   completed-P2 cavity-replacement callback and is installed for collapse as
   well as insertion. Thus SizeShape and StepDistance both evaluate the final
   post-`correct_cavity` elements with `AsDeg::Pk` geometry, frozen
   `AsDeg::P1` metric interpolation, and configured objective quadrature before
   topology commit.

   The geometry provenance is intentionally simpler than insertion. The
   collapsed seed edge and both original endpoints lie inside the removed
   cavity, so no child trace survives. `aux_bisecPointLen` therefore clears and
   leaves empty the split-edge provenance for collapse. Existing edges on the
   boundary of the union of endpoint balls retain their stored P2 coefficient;
   every new non-CAD edge from the replacement point to that boundary is an
   affine P2 spoke. CAD-owned new coefficients continue to be evaluated on the
   owning CAD entity. Assertions at the collapse setup make it impossible to
   accidentally reuse insertion's exact-child restriction path.

   A self-contained planar regression uses a very short, genuinely curved P2
   interior edge. It verifies that the endpoint balls are replaced
   successfully, both original vertices disappear, split provenance remains
   empty, every new replacement-point spoke has its coefficient at the affine
   midpoint, all new triangles are conservatively certified, and the global
   completed-P2 SizeShape sum decreases. A paired P1 case performs the same
   collapse without high-order coefficients and protects the historical path.
   With these additions, all thirteen high-order Phase 1/7 baseline cases pass
   in Debug and Release, the seven-case legacy cavity suite passes in Release,
   and the broader quality-insertion diagnostic target compiles successfully.

3. **Done (2026-08-26): objective-driven completed-P2 swaps.** Swap topology
   is still proposed from the P1 skeleton, and P1 orientation remains a
   prerequisite. For a high-order mesh, however, the temporary P1 objective
   is no longer allowed to accept or reject a candidate. The common cavity
   operator first reconstructs the complete high-order replacement, reuses
   every surviving boundary coefficient, creates genuinely new non-CAD edges
   as affine P2 edges, applies CAD ownership where present, and runs the common
   conservative validity classifier. A swap-specific final callback then
   evaluates the configured SizeShape or StepDistance objective with
   `AsDeg::Pk` geometry, frozen `AsDeg::P1` metric interpolation, and configured
   quadrature. Only an improving completed replacement is committed.

   This contract covers planar edge flips, tetrahedral face-to-edge swaps, and
   tetrahedral edge-to-face swaps through one completed-swap acceptance helper.
   Because a swap removes its old diagonal, face, or edge rather than splitting
   it, there is no parent-child restriction provenance. Existing edges on the
   cavity boundary retain their exact stored P2 coefficients; a new planar
   diagonal or other genuinely new non-CAD edge is affine P2. The mesh-wide
   CavityTargetAverage state is initialized with P2 geometry and is updated
   only after a successful topology commit. Rejected completed candidates roll
   back without modifying that state.

   The planar regression flips an edge in a short-edge patch and checks the
   configured completed-P2 objective decrease, agreement between the updated
   and independently recomputed mesh-wide StepDistance state when applicable,
   certification of both new triangles, exact reuse of every inherited edge
   coefficient, and affine placement of the new diagonal coefficient. A paired
   P1 flip protects the historical path. A full-dimensional 3D regression
   reconstructs a valid completed-P2 three-to-two tetrahedral edge swap whose
   configured objective is worse; it verifies objective rejection, complete
   rollback, and preservation of all three original certified tetrahedra. The
   high-order baseline suite now contains sixteen cases.

4. **Done (2026-08-26): explicit Classic-P1-only length-smoothing policy.** Two
   distinct smoothing families are kept separate. `smoothMeshLength` and
   `movePointCavLen` are legacy metric-length heuristics: they move a
   topological vertex toward a damped, length-weighted average of its neighbors
   and use P1 edge lengths and P1 validity. They remain available without
   behavioral changes in Classic P1 adaptation. They are not a substitute for
   SizeShape or StepDistance optimization.

   A P2 implementation would have to choose and update the adjacent edge-,
   face-, and volume-interior geometry variables consistently. Moving only the
   topological vertex changes every incident polynomial trace while leaving its
   other coefficients fixed, so the legacy update is not a defined P2
   smoothing policy. Consequently, an explicitly requested global P2 length
   pass now exits at `smoothMeshLength` with a clear diagnostic before scanning
   points, and Classic P2 cavity construction skips `movePointCavLen` while
   retaining its existing insertion point and reports that policy once. The
   defensive low-level routines assert P1 rather than silently pretending to
   support high order.

   Objective-driven adaptation now has exactly one insertion contract. The
   `lengthBased` and `worsenPctg` insertion arguments, their alternate
   `setCavityInsertion3` plus `checkCavityQuality` branch, the obsolete final P1
   quality helper, and the always-zero `nTryInsLen`/`nSuccInsLen` statistics
   columns have been removed. A failed objective insertion proceeds to the
   remaining objective-driven operations; it is never retried with Classic
   length acceptance. This also closes the dormant route through
   `insertLongEdges` that could previously select the fallback in a
   `TESTQUALITYALGO` build.

   This decision does **not** disable ordinary SizeShape/StepDistance ball
   smoothing or `CAVSMOOTHING`. The latter is objective-based cavity-targeted
   smoothing and is completed separately in Phase 7 item 5 below. A paired
   regression confirms that the Classic P1 length stage remains callable while
   the same global stage on an elevated P2 mesh returns zero without changing
   any coordinate.

5. **Done (2026-08-26): released-affine P2 cavity-targeted smoothing.** This
   operation has a different geometric family from ordinary ball smoothing.
   Ordinary smoothing activates one stored mesh control point and freezes the
   others. Cavity smoothing instead moves the provisional reconnection apex.
   Under the selected option-b policy, every trial rebuilds each non-CAD P2
   spoke incident to that apex as an affine high-order edge. Thus, if an
   insertion originally split a curved non-CAD edge into exact polynomial
   children, that split provenance is used for the unsmoothed candidate but is
   released as soon as an actual smoothing move is accepted. Collapse has no
   split-child provenance and enters the same released family directly.

   Existing edges on the boundary of the cavity are not incident to the moving
   apex and retain their stored P2 coefficients exactly. CAD ownership has
   higher authority than the affine-spoke rule: a CAD-curve apex is evaluated
   at its trial curve parameter, and a CAD-owned child coefficient is evaluated
   on that curve at the interpolated child parameter. Consequently every
   accepted CAD node remains exactly on its CAD entity. The implementation
   covers planar interior and CAD-curve cavities and extends directly to
   full-dimensional interior tetrahedral cavities. The pre-existing unsupported
   3D CAD-boundary cavity modes remain outside this step; ordinary CAD-surface
   ball smoothing is unaffected.

   The geometry dependence is included in the objective derivatives. For an
   apex coordinate `x` and fixed opposite vertex `X_j`, an affine P2 spoke has
   midpoint coefficient `C_0j(x) = (x + X_j)/2`. Therefore the exact frozen-
   metric, frozen-measure objective gradient is
   `dF/dx = partial F/partial x + (1/2) sum_j partial F/partial C_0j`.
   Equivalently, the completed P2 map varies as
   `delta X_h(lambda) = lambda_0 delta x`, which is the affine-reproduction
   identity for the quadratic basis. On a CAD curve `x = C(t)`, the chain rule
   uses `C'(t)` for the apex, one half of `C'(t)` for affine non-CAD spokes, and
   `C'((t+t_j)/2)/2` for CAD child coefficients. The optimizer forms a centered,
   validity-aware finite difference of this exact chain gradient for its
   Hessian and symmetrizes the result; no analytic StepDistance Hessian is
   claimed. SizeShape and StepDistance both use configured objective
   quadrature, complete P2 geometry, the established frozen P1 metric and
   quadrature-measure derivative contract, and conservative P2 validity at
   every trial. StepDistance target-average normalization is applied once,
   after assembling the cavity numerator derivatives.

   Acceptance is transactional. The released family is optimized from the
   current apex, but the resulting completed-P2 objective is compared with the
   preserved pre-smoothing cavity. A P2 success requires an actual coordinate
   or CAD-parameter displacement, conservative certification, and the existing
   objective-improvement rule. Rejection restores the point, metric, CAD
   parameter, and split provenance exactly. An accepted interior move clears
   non-CAD split provenance; an accepted CAD move retains its endpoint
   provenance so subsequent CAD child construction keeps the correct parameter
   interval. Near-miss cavity growth similarly snapshots and restores this
   state when a smoothed enlarged local candidate is later rejected.

   `test_high_order_cavity_smoothing` qualifies the released completion and
   chain rule in planar triangles and tetrahedra, checks frozen-sample analytic
   gradients against centered differences, exercises both SizeShape and
   StepDistance evaluation, and verifies accepted release versus exact rollback
   on rejection. The curved-CAD regression moves a planar boundary insertion
   in curve-parameter space and independently confirms the apex and both child
   coefficients against CAD evaluations. The existing sixteen-case insertion,
   growth, collapse, and swap suite was updated to reconstruct its local P2
   objective at callback time, because cavity smoothing can now legitimately
   move the apex before objective growth.

6. **Done (2026-08-26): degree-aware smoothing dependency reactivation.** Let
   `N_p(K)` denote the complete degree-`p` geometric-node set of element `K`.
   After a substantive accepted move on smoothing region `R`, the active set is
   now the union `A_p(R) = union_{K in R} N_p(K)`. Thus P1 retains its exact
   vertex-only contract, while a P2 triangle or tetrahedron reactivates its
   vertices and every shared edge-interior control point. Shared global point
   identifiers provide the required deduplication across incident elements.
   The same helper is used by the serial mesh-wide pass, its dormant parallel
   implementation, and targeted bad-element smoothing; this also repairs the
   dormant parallel loop's historical assignment to the moved point instead of
   the neighboring point.

   Targeted smoothing now scans the complete P2 node set of its seed element.
   A vertex receives its ordinary ball, and an edge-interior control point
   receives the shell of its owning edge. Higher-degree face- and
   volume-interior candidate policies are deliberately not inferred here. On a
   targeted success, every live element in the affected ball or shell is
   inserted into `BadEntHandler::affectedEnttsAlive` with a freshly recomputed
   configured objective value using `AsDeg::Pk` geometry. The handler can
   therefore replace stale qualities and priority-queue entries for the whole
   dependency region, rather than only the seed element.

   Reactivation does not relax boundary ownership. A boundary point with CAD
   ownership remains movable only in its CAD parameter space. A topological
   boundary record in a mesh without a CAD model supplies no such parameter
   space and is consequently kept fixed when neighboring geometry reactivates
   it; it is never treated as an unconstrained physical-space variable.
   Rejected or numerically null moves continue to deactivate only the attempted
   point and do not propagate reactivation.

   `test_high_order_smoothing_reactivation` contains four regressions: exact
   P1 corner-only tags, complete six-node P2 triangle tags with an unaffected
   point, complete ten-node P2 tetrahedron tags, and an accepted targeted P2
   edge-control move. The end-to-end case verifies movement, reactivation of
   every node in both incident triangles, exact preservation of a reactivated
   no-CAD boundary control point, conservative validity, and exact agreement
   (to `5e-14`) between every handler entry and an independent full-P2
   SizeShape evaluation. The classic, StepDistance, cavity-smoothing, and
   sixteen-case topology baseline suites remain green.

7. **Done (2026-08-26): identity-level P2 cavity conformity qualification.**
   The reconnection implementation already had the required ownership
   mechanisms, so this step did not introduce a second construction path. In
   planar cavities, the first occurrence of a new apex edge creates its P2
   control point and records the edge in the cavity edge table; every later
   triangle with the same endpoints copies that global point identifier.
   In tetrahedral cavities, the external-edge hash preserves old boundary-edge
   identifiers, the apex-edge hash shares newly created spoke identifiers even
   between tetrahedra reached through different boundary faces, and the face
   hash copies the complete shared face. For P2 there are no face-interior or
   volume-interior geometry nodes, so shared vertices plus shared edge-control
   identifiers are the complete conformity condition.

   The qualification oracle is deliberately independent of those construction
   tables. It scans every live full-dimensional element, reduces each local
   edge to the sorted endpoint key `{min(i,j),max(i,j)}`, and builds both maps
   `edge -> control point` and `control point -> edge`. It requires every
   occurrence of one edge to reference the same live global control point and
   requires distinct edges not to alias one control point. At the same time,
   every newly completed element is passed through
   `classify_element_validity<d,2>` and must be `Certified`. Conformity and
   validity are therefore checked together but remain logically independent:
   positive Jacobians do not prove a continuous trace, and shared traces do not
   prove positive Jacobians.

   The centroid-insertion regression applies this audit inside the completed-
   element acceptance callback, before topology commit, and repeats it on the
   committed mesh. The planar case creates three triangles and the direct
   full-dimensional 3D case creates four tetrahedra. The accepted planar
   collapse and face-swap regressions repeat the same audit over their final
   live meshes, covering the other completed-cavity clients without duplicating
   reconnection logic. The adjacent P1 insertion, collapse, and swap tests
   remain unchanged and green. This strengthens the existing sixteen-case
   topology baseline suite without adding another test case.

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
   - **Resolved in Phase 7.4:** legacy metric-length smoothing is explicitly
     Classic-P1-only. P2 exits before point traversal, and objective-driven
     insertion no longer contains a length-based fallback or compatibility
     statistics.
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

5. **Done (2026-08-25): implemented and independently validated P2
   ShapeVolume.** The dense Duffy oracle now evaluates the complete
   ShapeVolume embedding independently: centered log-eigenvalues represent
   shape, while `det(A) - 1/det(A)` supplies the injective volume coordinate.
   Curved P2 triangles and tetrahedra are checked for FE and analytical metric
   fields at objective quadrature orders two through five. Order five improves
   on order two in every case and differs from the eight-point-per-coordinate
   reference by at most approximately `4.40e-6` in 2D and `3.70e-9` in 3D.

   Frozen order-five samples then validate the first P2 edge control point in
   both Lagrange and Bezier bases. The production analytic gradient agrees
   with AD of an independently assembled ShapeVolume value to a maximum
   absolute discrepancy of approximately `1.55e-15`, under the `5e-14`
   scaled tolerance. A separately coded scalar-generic ShapeVolume gradient
   is differentiated to obtain the independent AD-of-gradient Hessian oracle;
   its maximum discrepancy from the production AD Hessian is approximately
   `8.88e-16`. Centered differences of the frozen value and analytic gradient
   have maximum errors of approximately `3.56e-11` and `1.25e-9`,
   respectively, within their `1e-7` scaled tolerances.

   The ShapeVolume matrix deliberately sets the ordinary collapse-barrier
   coefficient to `1e6`. Agreement with the oracle, which contains no
   additive barrier, verifies that this separate barrier is disabled by the
   ShapeVolume policy. Numerically singular P2 maps are also exercised in 2D
   and 3D for both metric-field and geometry-basis choices. Both public value
   and derivative entry points return the finite ShapeVolume rejection
   sentinel, with zero gradient and Hessian, before invalid geometric weights
   can enter the integral.

   This sentinel test qualifies ShapeVolume's numerical SPD-boundary rejection
   policy. Bernstein validity rejection and production smoothing behavior
   remain part of the later per-variant integration work. Phase 5 step 6 is
   the deliberately different CavityTargetAverage aggregation convention.

6. **Done (2026-08-25): qualified P2 CavityTargetAverage as a distinct
   reference-space and regional aggregation policy.** Each elemental value is
   the unit-weight reference-simplex average of the basic StepDistance
   integrand: objective quadrature remains configured by
   `objective_quadrature_order`, but sample `theta` is exactly one. The
   geometric or metric volume measure is therefore absent. The policy also
   suppresses the separate collapse barrier and remains mutually exclusive
   with ShapeVolume.

   An independent dense Duffy oracle checks curved P2 triangles and
   tetrahedra for FE and analytical metrics. The order-five errors are at most
   approximately `5.70e-6` in 2D and `1.86e-8` in 3D. With both Lagrange and
   Bezier active edge bases, the production analytic gradient agrees with AD
   of the independently assembled frozen value to `2.55e-15`; the production
   AD Hessian agrees with AD of the test-owned analytic gradient to
   `3.55e-15`. The maximum centered-difference errors are approximately
   `1.69e-10` for the gradient and `7.92e-10` for the Hessian. The test sets
   the nominal collapse-barrier coefficient to `1e6`, so agreement with the
   barrier-free oracle directly checks its suppression.

   A separate disconnected two-element P2 mesh checks the second level of the
   contract. `step_distance_global_objective_state`, regional contribution,
   `getmetquamesh`, and count-changing replacement all reproduce the
   arithmetic mean by live-element count; arbitrary legacy target weights are
   ignored, and changing `opt_pnorm` from one to four leaves the result
   unchanged. For fixed topology, the regional gradient and Hessian are the
   elemental derivatives divided by the element count. Centered differences
   on the two-element FE case agree to approximately `1.30e-10` and
   `9.51e-10`, respectively.

   This step qualifies the CavityTargetAverage value, frozen elemental
   derivatives, and aggregation algebra. Per-variant Bernstein rejection and
   production smoothing behavior remain in Phase 5 step 8; consistent with
   the roadmap, they are not silently skipped here.

7. **Done (2026-08-25): made the planar qualification an explicit gate before
   the matching 3D extension.** The focused executable is now organized into
   `p2_stepdistance_2d` and `p2_stepdistance_3d` Boost suites. Each suite owns
   the same 14 independently named cases: degree-aware value and derivative
   dispatch; basic value, gradient, and Hessian; collapse-barrier
   qualification; ShapeVolume value, gradient, Hessian, and numerical
   singularity rejection; and CavityTargetAverage value, gradient, Hessian,
   and arithmetic aggregation. Every case retains both FE and analytical
   metric coverage where applicable, and every differentiated variant retains
   both Lagrange and Bezier active edge bases.

   The 2D suite was selected and run alone first, with the complete 3D suite
   explicitly disabled; all 14 planar cases passed. The 3D suite was then
   selected and run alone with the planar suite disabled; all 14 direct
   extensions passed. CTest now registers
   `test_high_order_stepdistance_value_2d` and
   `test_high_order_stepdistance_value_3d` separately so CI can report the
   dimensional gates independently. The original combined CTest entry remains
   intact for compatibility and full-matrix regression coverage.

   This step changes qualification structure only. It does not select a
   StepDistance variant for production P2 smoothing, and it does not claim
   Bernstein rejection coverage for basic, barrier, or CavityTargetAverage.
   Those validity and smoothing integrations remain the explicit work of
   Phase 5 step 8.

8. **Done (2026-08-25): qualified every P2 StepDistance variant for
   conservative smoothing and enabled the production objective.** The new
   `test_high_order_stepdistance_smoothing` executable perturbs an interior P2
   edge control point, first in the two-triangle planar mesh and then in a
   complete tetrahedral edge shell. For Basic, CollapseBarrier, ShapeVolume,
   and CavityTargetAverage it requires a nonzero production ball-smoothing
   move, a nonincreasing configured mesh objective, and a Bernstein-certified
   result for every affected element. It then constructs a noncertified
   position for the same control point and verifies that the low-level
   smoother rejects it, returns a nonzero status, and restores the candidate
   coordinates exactly. The matching `p2_stepdistance_2d` and
   `p2_stepdistance_3d` suites are registered as independent CTest entries.

   This completes the per-variant matrix together with the value, analytic
   gradient, and AD/finite-difference Hessian evidence from steps 3--7. During
   the smoothing qualification, the exact Basic, CollapseBarrier, and
   CavityTargetAverage regional Hessians were observed not to be positive
   definite at otherwise valid P2 configurations. The objective derivatives
   are not modified. Instead, the Newton driver now has an opt-in bounded
   steepest-descent fallback when Cholesky fails or the computed direction is
   not a descent direction. Only high-order StepDistance ball smoothing
   enables it. Its fallback step and spherical trust radius are both one
   tenth of the shortest P1 corner-edge length in the affected region. The
   existing Newton behavior is unchanged when the option is disabled. This
   qualification scope is not a convexity claim: neither StepDistance nor
   SizeShape has a general positive-definite-Hessian guarantee with respect to
   a moving control point, at P1 or P2. P2 StepDistance is where an indefinite
   regional Hessian was exposed by the production smoothing regressions;
   extending the fallback to other objectives or degrees requires its own
   optimizer qualification rather than an assumption that their Hessians are
   always positive definite.

   CavityTargetAverage revealed a second integration conflict: its global
   arithmetic mean accepted a move while the legacy local worst-element guard
   rejected it. The subsequent acceptance-policy audit showed that the same
   guard was also being layered on top of SizeShape and the regional
   StepDistance variants. All objective-driven smoothers now rely exclusively
   on their configured regional or global objective plus conservative
   validity; the legacy worst-element veto is retained only for the
   non-objective-driven Distortion and Unit paths. This separation is
   centralized in
   `isObjectiveDrivenSmoothing` and also covers the dormant
   `IMPROVEMAXQUAL` element and cavity checks if they are re-enabled.
   Conservative Bernstein rejection remains mandatory for all four
   StepDistance variants.

   The temporary production exception that selected SizeShape for planar P2
   has consequently been removed. In a StepDistance build,
   `productionSmoothingObjective` now returns StepDistance at P1 and P2 in
   both 2D and 3D. The Phase 4 SizeShape regressions request SizeShape
   explicitly and remain active, while the public planar P2 optimizer test now
   exercises the production StepDistance selection. Phase 5 is complete;
   CAD-boundary specialization remains the explicit subject of Phase 6.

   Release and Debug builds both pass the complete 28-case StepDistance value,
   gradient, and Hessian executable, the new eight-case StepDistance
   smoothing/validity executable, and the nine-case high-order smoothing
   regression containing the public production optimizer. The dimension-
   specific CTest registrations are present. A full fixture-driven CTest
   attempt did not reach them because the pre-test `deploy_cases` fixture
   failed independently while generating the 100k-square mesh, at the
   analytical-metric assertion `eigval[i] > 0`; the focused executables were
   therefore also run directly in both configurations.

## Phase 6 CAD and boundary execution record

1. **Done (2026-08-25): completed planar P2 curved-CAD-edge creation and
   smoothing.** Degree elevation now converts the elevated polynomial geometry
   to Lagrange form and then makes every topological edge-interior coordinate
   the exact CAD image of its already elevated edge parameter, before metric
   localization. This ordering is essential: Metris holds the incoming P1
   geometry in Bezier form, so evaluating a CAD point before the
   Bezier-to-Lagrange conversion would incorrectly treat that physical point
   as a Bezier coefficient. The resulting invariant is explicit:
   `coordinate = CAD_curve(parameter)`. The same creation path applies to
   topological CAD-curve nodes in 3D without a separate implementation.

   High-order curve-node smoothing continues to use the existing parameter-
   space chain rule through the CAD tangent and curvature. It now additionally
   requires Lagrange geometry and constrains the Newton trust region to remain
   strictly between the two endpoint parameters of the classified mesh edge.
   This prevents a high-order node from crossing a topological endpoint while
   remaining on the untrimmed underlying curve. SizeShape and StepDistance
   also opt into the already qualified bounded steepest-descent fallback when
   their parameter-space Hessian is not positive definite; classic smoothing
   receives the geometric interval constraint without changing its Hessian
   policy.

   The new regression constructs a genuinely curved planar EGADS boundary in
   memory: a circular arc and two straight segments bound one triangle. It
   verifies that P1-to-P2 elevation places the arc node at the CAD image of the
   interpolated parameter rather than at the chord midpoint. After perturbing
   the parameter, both SizeShape and StepDistance must produce a nonzero
   parameter-space move, not increase their objective, keep the parameter in
   the topological edge interval, restore the exact CAD coordinate, and leave
   the P2 triangle Bernstein-certified. No external binary CAD or mesh test
   data is needed.

   Release and Debug builds pass both new curved-CAD cases, the nine-case
   high-order classic/production smoothing executable, the eight-case P2
   StepDistance smoothing/validity executable, and the five-case Phase 1
   elevation/validity baseline. Phase 6 step 2 remains the distinct extension
   for P2 nodes classified on triangular CAD faces in 3D; surface orientation,
   broader CAD regression data, per-move invariant auditing, and curved-edge
   self-intersection remain steps 3--6 and are not claimed here.

2. **Done (2026-08-25): implemented exact P2 CAD-surface placement before
   3D background localization.** During degree elevation, the existing face
   path treats the P1 vertex parameters as two scalar Lagrange fields and
   evaluates them at every target face node. Thus, for a P1-to-P2 interior
   triangulation edge with endpoint parameters `(u0,v0)` and `(u1,v1)`, the
   new face-classified node receives their parameter-space midpoint. After the
   elevated polynomial geometry is converted from Bezier to Lagrange form,
   the runner now sets its physical coordinate to the exact CAD image
   `S(u,v)` before FE metric localization.

   The placement pass follows the lowest-dimensional ownership rule. A node
   whose primary boundary classification is a CAD curve is deliberately
   skipped by the surface pass and remains at `C(t)`, even though it also has
   an associated face record. Only nodes whose primary classification is a
   CAD face are placed with `S(u,v)`. The loop naturally includes the edge
   nodes needed by P2 and is also ready to place true face-interior Lagrange
   nodes once P3 is enabled.

   The new self-contained FE regression constructs a spherical EGADS patch,
   triangulates its rectangular parameter domain with two P1 triangles, and
   elevates them to P2. The control point on the shared diagonal is not on a
   topological CAD curve and is therefore the required face-classified case.
   The test verifies its exact interpolated `(u,v)`, proves that its physical
   coordinate is `S(u,v)` rather than the spatial chord midpoint, and checks
   that background localization assigns it a valid `poi2bak` seed. A boundary
   P2 node is checked in the same model to prove that the subsequent surface
   pass does not overwrite the stricter curve constraint.

   Release and Debug builds pass the new spherical-surface regression, both
   curved-CAD-edge objective cases, the nine-case high-order classic and
   production smoothing executable, and the five-case Phase 1 baseline.

3. **Deferred (2026-08-25): embedded surface-triangle orientation
   certification.** For a full-dimensional 3D mesh, a strictly positive
   incident-tetrahedron Jacobian implies that its boundary-face restriction is
   locally regular: the two restricted tangent vectors cannot become dependent
   while the three-dimensional Jacobian remains nonsingular. A separate
   surface certificate would therefore duplicate the current local validity
   guard for volume meshes. It remains necessary for embedded 2D meshes in 3D,
   where no incident tetrahedron exists, and will be implemented when that mesh
   class enters scope. This does not replace the later global curved-boundary
   self-intersection work in Phase 6 step 6.

4. **Done (2026-08-25): established small reusable CAD-curve and CAD-surface
   regression meshes.** The self-contained planar case consists of one
   triangle bounded by a circular arc and two straight CAD segments. The
   self-contained 3D surface case consists of two triangles covering a
   spherical CAD patch, with a shared interior triangulation edge. Both are
   constructed in memory so the regression has no external binary CAD or mesh
   dependency.

   The existing focused checks now share explicit invariant routines, and each
   regression covers both P1-to-P2 elevation and native P2 re-ingestion through
   `MetrisAPI`. After re-ingestion, the tests convert internal Bezier storage to
   Lagrange form before checking the geometric-node contract. The planar curve
   case also runs the complete zero-adaptation `runMetris()` lifecycle and
   checks the invariants afterward; it verifies the elevated parameter, exact
   `C(t)` coordinate, non-chord geometry, and full-dimensional P2 validity. The
   surface case verifies that native P2 input requests no further elevation,
   then checks the shared-edge `(u,v)`, exact `S(u,v)` coordinate, background
   localization, and preservation of the lower-dimensional curve ownership.

   A complete `runMetris()` lifecycle is intentionally not claimed for the
   embedded surface case. Its final metric-cost report currently assumes that
   a topological 2D mesh also has geometric dimension 2; supporting the
   `tdim=2`, `gdim=3` dispatch belongs with the deferred embedded-mesh work, not
   the planar-first Phase 6 scope.

   Exercising the full topology audit also qualified the regression inputs
   themselves. Straight EGADS lines use unit directions with physical-length
   parameter ranges, endpoint parameters are obtained by inverse evaluation,
   and the surface mesh supplies a parameter record for every
   vertex--incident-face pair. Release and Debug builds pass all three curve
   cases and both surface cases.

5. **Done (2026-08-25): made the accepted boundary-move invariant explicit in
   production smoothing.** Both production ball-smoothing acceptance paths now
   perform the same final audit after their objective policy accepts a move and
   before the move is counted as successful. Every curve or surface parameter
   in the moved point's complete boundary-record chain must be finite, and
   evaluating its CAD entity must reproduce the stored physical coordinate to
   the established `1e-10` geometric tolerance. Every affected
   full-dimensional element must also pass the degree-aware validity policy:
   the P1 orientation check for linear meshes and the conservative Bernstein
   certificate for P2 and higher meshes.

   For a curve-owned point in 3D, the primary `t` remains authoritative. Before
   the final audit, each secondary incident-face record is synchronized with
   `EG_getEdgeUV`, so its `(u,v)` evaluates to the same physical point. The
   synchronization and audit are shared by complete interior-ball smoothing
   and element-targeted smoothing rather than duplicated in the two acceptance
   loops. The implementation is dimension- and degree-generic, covering the
   planar P2 priority and the same full-dimensional 3D boundary path without a
   separate policy.

   Boundary parameters are snapshotted as a complete linked record chain before
   each trial. Ordinary rejection restores coordinates, metric data, and all
   parameters; an exception during synchronization or final verification does
   the same before propagating the diagnostic. Thus a partially updated
   secondary CAD record cannot survive a failed acceptance path.

   The curved planar regression now includes the production
   `smoothInterior_Ball` path, with every other point tagged out so that only
   the perturbed P2 circular-arc node can move. It requires a nonzero accepted
   parameter update, checks every resulting CAD record against the physical
   coordinate, and re-certifies the P2 triangle. Release and Debug builds pass
   all four curved-CAD cases, the nine-case classic/production smoothing suite,
   and the eight-case P2 StepDistance smoothing and validity suite. Embedded
   `tdim=2`, `gdim=3` production smoothing remains under the explicit deferred
   scope recorded in items 3 and 4.

6. **Done (2026-08-25): added a global planar P2 boundary-intersection guard to
   both production smoothing acceptance paths.** After a candidate boundary
   move has passed its objective and local element-validity checks, every
   topological boundary edge containing the moved point is compared with every
   other live boundary edge. Only contact at a common topological endpoint is
   legal. Any crossing, overlap, or other contact rejects the candidate as an
   ordinary smoothing failure and restores its coordinate, metric, and complete
   CAD-parameter snapshot before the move is counted.

   The check operates on the actual quadratic polynomial edge, not on its
   endpoint chord. For the current Lagrange representation, the physical
   midpoint sample `L1` is converted to the equivalent quadratic Bezier control
   point `B1 = 2 L1 - (B0 + B2)/2`. Recursive de Casteljau subdivision then
   supplies convex-hull bounding boxes for separation and reduces possible
   contacts to tolerance-aware straight-segment proximity tests once both
   subcurves are sufficiently flat. Because only edges incident to the moved
   point can have changed, the local set is compared globally without repeating
   an unrelated all-pairs boundary audit after every move.

   The focused regression elevates two disjoint planar triangles to P2, verifies
   their initially straight boundary is intersection-free, and then bows one
   edge into the other component. Both P2 triangles retain positive Bernstein
   Jacobian certificates, proving that the existing local validity policy alone
   cannot detect the global contact; the new polynomial-edge check rejects it.
   Restoring the straight edge is accepted again, which also exercises legal
   shared endpoints between adjacent edges.

   Release and Debug builds pass the new one-case intersection regression, all
   four curved-CAD smoothing cases, both curved-CAD surface-projection cases,
   the nine-case classic/production P2 smoothing suite, and the eight-case P2
   StepDistance smoothing and validity suite. This completes Phase 6 for the
   agreed priority: full-dimensional planar P2, together with the
   full-dimensional 3D CAD support that extended naturally in earlier steps.
   Global self-intersection of embedded surfaces, general curved surfaces in
   3D, and P3-or-higher boundary curves remain separate future extensions.
