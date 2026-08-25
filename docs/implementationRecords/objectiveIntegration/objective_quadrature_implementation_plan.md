# Metris objective and quadrature implementation plan

This record captures the finalized incremental implementation plan and the
state reached at the end of the quadrature work.

## Scope and conventions

- Follow the generic formulation in the IMR paper: a pointwise objective
  supplies `psi(A)`, a nonnegative function with unique ideal value zero that
  already includes the objective exponent `p`.
- Let the objective integration layer see `psi`, not a raw quality followed by
  an external shift, absolute value, or power.
- Require each objective implementation to define and compute its complete
  `psi`, gradient, and Hessian explicitly, including its own use of `p`,
  regularization, and ideal-state limits. Do not introduce a shared scalar
  transform, generic power-chain wrapper, or common internal machinery for
  converting a raw `f` into `psi`, even when formulas look similar.
- Share only the pointwise result contract, quadrature-sample preparation, and
  integration traversal. The integrator must consume final `psi` derivatives
  and must never apply `p` or reinterpret an objective's pointwise result.
- Use one quadrature selection and one quadrature traversal for SizeShape and
  all StepDistance variants.
- Keep quadrature points and weights on the canonical reference simplex.
- Normalize weights to sum to one, preserving the current reference-element
  average convention.
- Evaluate geometry, metric, and the complete objective integrand separately
  at every quadrature point. Never form an elemental metric by averaging its
  tensor components.
- Use the paper's sample-level structure `w_r[theta_r psi(A_r) + B_r]`.
  Geometry and integration-metric density belong to `theta`; an explicit
  collapse barrier remains additive and deliberately outside `theta`.
- Select the production `theta` convention independently of the pointwise
  objective: without `INTQUALINRIEMSPACE`, use physical measure; with it, use
  physical measure times `sqrt(det(M))`. Apply this to SizeShape, original
  StepDistance, and ShapeVolume. Keep CavityTargetAverage as the explicit
  normalized-reference-average exception because it is not the paper's
  integral formulation.
- Require every supported production quadrature rule to have strictly positive
  weights.
- Freeze the integration factor `theta` and the target metric at each
  quadrature point when differentiating, as in the paper. Retain derivatives of
  explicit collapse barriers, which remain outside `theta`.
- Keep the objective exponent as a runtime `double` that is fixed for a run.
  Do not make it a template parameter: Metris uses C++17 and StepDistance
  supports noninteger exponents.
- The new objective-based SizeShape is always based on
  `f_SS(A) = q_SS(A) - 1`. It does not support the inverse `opt_power == -1`
  convention.
- Preserve the historical Classical Distortion path, including both existing
  `opt_power` signs and its historical sampling, unless a later explicitly
  opt-in quadrature option is selected.

## Incremental implementation

1. **Done (2026-08-19).** Added standalone triangle and tetrahedron tests for
   the historical vertex-barycenter rule. They verify point layout, equal
   positive weights, sum-to-one, simplex membership, and exactness for
   constants and affine functions.

2. **Done (2026-08-19).** Added a lightweight, non-owning
   `SimplexQuadratureView<tdim>` and point view. It exposes a common `size()` and
   indexed barycentric-point/weight interface without owning or choosing rule
   storage. The step-1 tests now consume this interface while retaining their
   test-local vertex-barycenter data.

3. **Done (2026-08-19).** Added a production vertex-barycenter-rule factory
   backed by static-lifetime triangle and tetrahedron tables. The tests now
   exercise this production helper directly; their duplicate test-local rule
   construction was removed. No objective evaluation path changed. Renamed the
   rule from “legacy” to the descriptive “vertex-barycenter” terminology.

4. **Done (2026-08-19).** Refactored StepDistance value evaluation to obtain
   its point count, barycentric coordinates, and normalized weights from the
   shared vertex-barycenter helper. The integrand and metric-evaluation logic
   remain unchanged. Focused StepDistance tests pass before and after the
   refactor, including comparisons with the unchanged differentiated path.

5. **Done (2026-08-19).** Refactored StepDistance differentiated evaluation to
   obtain its point count, barycentric coordinates, and normalized weights
   from the shared vertex-barycenter helper. All seven focused StepDistance
   tests pass before and after the refactor, including value consistency,
   rejection behavior, 2D and 3D gradient/Hessian finite differences, and the
   metric-volume barrier checks.

6. **Done (2026-08-19).** Routed P1 SizeShape value evaluation through the
   shared vertex-barycenter helper. Each point now evaluates its own metric and
   complete quality-error integrand; metric components are no longer averaged
   before evaluation. Vertex metrics are sampled directly, while barycenter
   metrics use physical analytical evaluation or FE barycentric interpolation.
   The existing physical/Riemannian integration density is preserved. Both
   objective-enabled and release builds compile, and the quadrature,
   StepDistance regression, and objective-comparison tests pass. Existing
   high-order nodal sampling remains unchanged.

7. **Done (2026-08-19).** Established and documented the compatibility
   boundary in `/tmp/metris_objective_compatibility_audit.md`. The audit maps
   all SizeShape, StepDistance, and Distortion dispatch paths and classifies
   the ownership of `opt_power`, `opt_pnorm`, and `step_distance_p`. It also
   identifies the shared `quafun_tradet` singularity-policy coupling that must
   be separated later. Added a self-contained Classical Distortion regression
   on a non-ideal anisotropic P1 triangle. It pins direct and inverse values,
   their reciprocal relation, agreement of value and differentiated entry
   points, and every local-node gradient and Hessian against centered finite
   differences for both `opt_power` signs. The focused test passes in both the
   release/Classical and objective-enabled build configurations.

8. **Done (2026-08-19).** Added the shared runtime `double` parameter
   `objective_p` with the paper-aligned default `1.0` and validation
   `objective_p >= 1`. Added `--objective-p` as the canonical command-line
   option and retained `--step-distance-p` as a deprecated compatibility
   alias; identical dual specification is accepted and conflicting values are
   rejected. Routed all StepDistance value, gradient, Hessian, one-dimensional,
   and statistics paths to `objective_p`. Added focused default, canonical
   option, compatibility alias, conflict, and validation tests. The seven
   existing StepDistance value/derivative tests pass at their fixed legal
   exponents. Strengthened the Classical Distortion regression to prove that
   changing `objective_p` does not affect either `opt_power` sign. Both
   objective-enabled and release/Classical builds and tests pass.

9. **Done (2026-08-19).** Added an owning pointwise-objective result contract
   carrying final `psi`, gradient, packed Hessian, and explicit derivative
   availability. Added temporary value and differentiated adapters, restricted
   to SizeShape and StepDistance, which reproduce the current objective entry
   points without exposing raw `f`, applying `p`, adding a barrier, or using
   any shared scalar transform. The production integrators remain untouched.
   Focused tests verify exact agreement with the current pointwise values,
   gradients, and Hessians for both objectives, result-state behavior, and
   rejection of an invalid differentiated-value request. The tests pass in
   both objective-enabled and release/Classical build configurations, and the
   existing objective-parameter, Classical Distortion, and focused
   StepDistance regressions remain unchanged.

10. **Done (2026-08-19).** Changed the SizeShape pointwise value to the paper
    definition `f_SS = q_SS - 1` and `psi_SS,p = f_SS^p`, computed explicitly
    inside `quafun_sizeshape`. The value no longer reads `opt_power` or offers
    inverse SizeShape behavior. It clamps only ideal-scale roundoff, using 32
    double-precision epsilons, and rejects substantially negative `f_SS` as an
    invalid state. Replaced the trace/determinant helpers' implicit
    `opt_power` singularity decision with an explicit policy for value,
    differentiated, and SurrealS entry points. SizeShape always requests
    rejection; Distortion and Unit explicitly preserve their historical
    direct/inverse policy. Added self-contained 2D and 3D pointwise tests for
    ideal zero, nonnegativity, anisotropic states, `p = 1`, `p = 2`, and
    noninteger `p = 1.5`, independence from `opt_power`, invalid indefinite
    metrics, singular elements, and preserved Classical inverse behavior.
    The focused SizeShape, pointwise-contract, and Classical regressions pass
    in objective-enabled and release builds. This deliberately stages only
    the value definition: the historical differentiated SizeShape algebra
    remains until step 11, and historical integration post-processing remains
    until step 12.

11. **Done (2026-08-19).** Replaced the historical differentiated SizeShape
    algebra with an explicit implementation of the complete
    `psi_SS,p = (q_SS - 1)^p` value, gradient, and Hessian inside
    `d_quafun_sizeshape`. The implementation first computes the derivatives of
    `q_SS`, then applies the objective's own explicit `p` formulas without a
    shared scalar-transform helper. It has dedicated `p = 1` and `p = 2`
    branches and an explicit ideal-state limit: `p = 1` retains the gradient
    and Hessian of `q_SS`, while every `p > 1` returns zero gradient and Hessian
    at the ideal state. Removed all remaining `opt_power` and inverse behavior
    from differentiated SizeShape. Added self-contained triangle and
    tetrahedron gradient/Hessian finite-difference tests for `p = 1`, `p = 2`,
    and `p = 1.5`, including exact-ideal and near-ideal states and explicit
    independence from `opt_power`. Updated the older mesh-backed SizeShape
    expectations to zero-ideal semantics. The new value, derivative,
    pointwise-contract, Classical Distortion, and focused StepDistance tests
    pass in the objective-enabled build; the derivative, contract, and
    Classical tests also pass in the release build.

12. **Done (2026-08-19).** Simplified P1 SizeShape value integration to
    accumulate the complete pointwise `psi` directly. The P1 path no longer
    shifts by `difto`, takes an absolute value, applies `opt_pnorm`, or changes
    behavior with `opt_power`; the separate CAD normal-deviation term remains
    unchanged. Added self-contained triangle and tetrahedron regressions that
    independently reconstruct the vertex-barycenter integral from
    `quafun_tradet`. They verify agreement with the step-6 historical direct
    value for `p = 1`, invariance to `difto`, `opt_power`, and (in the release
    configuration) `opt_pnorm`, and correct direct integration for `p = 2`.
    The new tests and the focused SizeShape value/derivative,
    pointwise-contract, Classical Distortion, and StepDistance regressions pass
    in the objective-enabled build; the new tests, derivative, contract, and
    Classical regressions also pass in the release build.

13. **Done (2026-08-19).** Routed P1 SizeShape value, gradient, and Hessian
    evaluation in `d_metqua` through the same vertex-barycenter points,
    per-sample analytical or FE metrics, and weights as the value path. The
    differentiated integrator now consumes the shared pointwise-objective
    result contract and directly accumulates the complete `psi` derivatives;
    it applies no `difto`, absolute-value, `opt_pnorm`, or `opt_power`
    transformation. The target metric and integration weight are explicitly
    frozen for differentiation, including physical/Riemannian density when
    enabled. Added self-contained triangle and tetrahedron tests for both FE
    and spatially varying analytical metrics. They verify value agreement with
    `metqua`, independent quadrature reconstruction, every local vertex's
    gradient and Hessian against frozen-sample finite differences, caller
    `idifmet` isolation, value-only differentiated entry, and release-build
    independence from `opt_pnorm`. The new tests and all focused SizeShape,
    pointwise-contract, Classical Distortion, and StepDistance regressions pass
    in the objective-enabled build; the new tests and the focused SizeShape,
    contract, and Classical regressions pass in the release build. The CAD
    normal-deviation penalty remains a separate legacy term; its normalization
    and plumbing are intentionally left to step 14.

14. **Done (2026-08-19).** Normalized StepDistance to the shared pointwise
    `psi` contract in both value and differentiated P1 integration. The
    specialized path is now selected by `QuaFun::StepDistance` in every build,
    rather than being hidden behind the `STEPDISTANCE` compile definition, and
    it consumes `PointwiseObjectiveResult` for value, gradient, and Hessian.
    Renamed the remaining local `phi` terminology to `psi`. Removed the common
    StepDistance scalar-power helper and wrote the complete regularized
    `(d^2 + eps^2)^(p/2) - eps^p` value and derivative factors explicitly in
    the original and ShapeVolume implementations. Elemental StepDistance
    integration no longer reads or asserts on `opt_pnorm`, ignores `difto`,
    and excludes the CAD normal-deviation penalty; ShapeVolume rejection and
    the additive metric-volume barrier remain distinct integration policies.
    Whole-mesh SizeShape and StepDistance aggregation now returns the complete
    additive objective (or the explicitly selected StepDistance mean) without
    applying a second `opt_pnorm` root; Classical Distortion and Unit retain
    their historical aggregation. Moved the debug-only Classical
    `opt_pnorm == 1` assertion behind the returning P1 objective paths. Added
    regressions for elemental StepDistance value/gradient/Hessian independence
    from `difto` and `opt_pnorm`, whole-mesh StepDistance and SizeShape
    independence from `opt_pnorm`, and preservation of the explicit
    CavityTargetAverage mean. The seven focused StepDistance tests and the
    quadrature, SizeShape value/derivative/integration, pointwise-contract, and
    Classical Distortion suites pass in both objective-enabled and release
    builds.

15. **Done (2026-08-19).** Added a shared, owning
    `ObjectiveQuadratureSample<gdim,tdim,mshdeg>` and preparation function.
    Each sample now carries the barycentric coordinates and normalized rule
    weight, physical coordinates, canonical- and regular-reference Jacobian
    transposes, frozen analytical or FE target metric, and one explicit
    integration factor `theta`. The compile-time `mshdeg` is used directly by
    `eval2`/`eval3`, so geometry evaluation is degree-aware even though the
    first production consumers remain P1. Exact rule vertices preserve the
    historical direct nodal metric sampling; nonvertex analytical metrics are
    evaluated at the physical sample point and nonvertex FE metrics by
    barycentric interpolation. Made the currently different integration
    conventions explicit through reference-average, physical, physical-metric,
    and regular-reference-metric `theta` modes, preserving all existing
    SizeShape and StepDistance values while isolating the convention choice for
    later unification. Added `theta_is_valid` so preparation reports a
    singular geometry or metric without preempting ShapeVolume's finite
    rejection sentinel; consumers enforce the status only after pointwise
    rejection. Routed all P1 SizeShape and StepDistance value and
    differentiated paths through this shared preparation layer, removing their
    duplicated geometry, metric, Jacobian, and weight setup. Added focused
    triangle and tetrahedron FE tests, analytical-metric vertex/nonvertex
    tests, all four `theta` modes, invalid-`theta` reporting, and a curved P2
    test that verifies degree-2 physical coordinates and Jacobians against
    direct quadratic evaluation and distinguishes them from P1 preparation.
    The new tests and the focused StepDistance, SizeShape integration and
    derivative, pointwise-contract, Classical Distortion, and quadrature suites
    pass in both objective-enabled and release builds.

16. **Done (2026-08-19).** Added a common, degree-templated value traversal in
    `objective_quadrature_value.hxx`. It alone iterates the selected rule,
    prepares geometry and frozen-metric samples, validates `theta`, and
    accumulates the weighted result. Explicit SizeShape and StepDistance
    policies select their current `theta` convention and return the complete
    pointwise `psi` together with any integration-level rejection or additive
    barrier contribution. ShapeVolume retains rejection before invalid-theta
    enforcement, CavityTargetAverage remains a reference average, and the
    metric-volume barrier remains outside `theta`. The common accumulator also
    preserves the historical arithmetic ordering for ordinary samples and for
    samples with an additive barrier, avoiding value/derivative roundoff drift.
    Replaced both duplicated P1 value loops in `low_metqua.cxx` with calls to
    this traversal; the Classical paths remain unchanged. The focused 32-case
    StepDistance, SizeShape integration/derivative/pointwise, sample,
    quadrature, contract, and Classical suites pass in both objective-enabled
    and release builds.

17. **Done (2026-08-19).** Added the degree-templated common differentiated
    traversal in `objective_quadrature_derivatives.hxx`. It owns rule
    iteration, sample preparation, derivative-order selection, output
    initialization, `theta` value/gradient/Hessian handling, and the complete
    product-rule accumulation of `theta*psi + B`. Explicit SizeShape and
    StepDistance derivative policies provide their complete pointwise
    value/gradient/Hessian results and objective-specific rejection and
    metric-volume barrier data while reusing the value policies' `theta`
    conventions and validity checks. Replaced both duplicated P1 loops in
    `low_metqua_d.cxx`; ShapeVolume rejection ordering, CavityTargetAverage,
    frozen metric and frozen-theta behavior, and optional geometric-theta
    differentiation remain unchanged. Added a degree-aware regular-reference
    basis-gradient helper for Lagrange nodes and Bezier control points. Focused
    perturbation checks verify that gradient against the prepared regular
    Jacobian for P1 triangle and tetrahedron vertices and a P2 triangle control
    point in both bases, establishing the interface needed by future
    high-order traversal activation. The focused 32-case StepDistance,
    SizeShape integration/derivative/pointwise, sample, quadrature, contract,
    and Classical suites pass in both objective-enabled and release builds.

18. **Done (2026-08-19).** Removed the historical CAD-normal deviation term
    from every active objective-driven path. P1 SizeShape now returns through
    the common value and differentiated traversals before `nordev`, just as
    ordinary StepDistance, CavityTargetAverage, and ShapeVolume already do.
    Consequently `qua_surf_wt_quality`, `qua_surf_wt_normal`, and CAD-normal
    evaluation no longer alter objective values returned through either
    `metqua` or `d_metqua`. The Classical and current high-order compatibility
    paths retain their historical `nordev` behavior. Added an embedded-surface
    regression that advertises CAD presence without supplying usable CAD
    topology, exercises value-only, gradient, and Hessian entry points for
    SizeShape and all three StepDistance modes, and verifies invariant scalar
    results under widely different surface weights. The focused 33-case
    objective, quadrature, contract, and Classical suites pass in both the
    objective-enabled and release builds.

    **Deferred issue noted:** for an embedded surface (`gdim = 3`, `tdim = 2`),
    repeated SizeShape gradient/Hessian evaluations can produce different
    components, and optimized builds can produce non-finite components. The
    scalar SizeShape value remains stable. This is separate from `nordev` and
    will not be investigated during the present unification work.

19. **Done (2026-08-19).** Unified objective dispatch in `low_metqua.cxx` and
    `low_metqua_d.cxx`. A single compile-time branch now recognizes SizeShape
    and StepDistance, selects the common vertex-barycenter rule, and passes
    `iquaf` directly to the shared value or differentiated traversal. Both
    objectives therefore have the same outer consumption path; their
    pointwise `psi`, derivatives, `theta` convention, additive barrier, and
    rejection behavior remain selected below through their explicit policies.
    Moved historical compatibility scratch arrays, function pointers, `pnorm`,
    and `nordev` setup below the objective return. Preserved the current support
    boundary exactly: P1 uses the objective traversal, high-order SizeShape
    continues through its compatibility path, and StepDistance continues to
    assert its P1 requirement. The focused 33-case objective, quadrature,
    contract, and Classical suites pass in both objective-enabled and release
    builds.

20. **Done (2026-08-19).** Added dedicated SizeShape metric-sampling contracts
    for FE and analytical metrics on both triangles and tetrahedra. Each case
    first uses its varying anisotropic metric to verify that averaging the
    pointwise integrand is measurably different from averaging the sampled
    metric tensor and evaluating `psi` once. It then installs an explicitly
    constant anisotropic metric and verifies that the complete quadrature value,
    every control-point gradient, and every Hessian reduce to one pointwise
    evaluation multiplied by the sum of the integration weights. Finally, it
    installs the identity metric on an ideal regular simplex and verifies zero
    integrated value for `p = 1`, `1.5`, and `2`, together with zero gradient
    and Hessian for `p = 1.5`. These complement the existing frozen-sample and
    finite-difference checks for varying FE and analytical metrics. The focused
    suite now contains 37 cases and passes in both objective-enabled and release
    builds.

21. **Done (2026-08-19).** Added the shared integer parameter and command-line
    option `--objective-quadrature-order`, initially with default value `0` and
    later changed in step 27 to use `-1` as the automatic-default sentinel.
    Added a
    single runtime `get_objective_quadrature` selector and routed both the value
    and differentiated objective paths through it. Order `0` returns the
    historical vertex-barycenter rule exactly. Parameter validation and the
    selector both reject negative values and not-yet-available positive orders,
    preventing silent fallback while the positive tables are staged. The
    parameter participates automatically in copying, equality, and JSON
    serialization through the central parameter-field list. Added command-line,
    default, invalid-order, and order-zero selection tests. The focused suite
    now contains 40 cases and passes in both objective-enabled and release
    builds.

22. **Done (2026-08-19).** Imported compact degree-2 and degree-3 positive
    triangle and tetrahedron tables from the identified SANS source revision,
    with source-file and commit attribution in the production header. Converted
    SANS Cartesian reference coordinates to Metris barycentric coordinates
    while preserving point order and the normalized weight convention. Added a
    compile-time accessor for the four staged rules and an availability smoke
    test. The runtime objective selector deliberately remains restricted to
    order 0. The quadrature target builds and all four cases pass in both the
    objective-enabled and release configurations.

23. **Done (2026-08-19).** Added a complete contract test for each imported
    triangle and tetrahedron rule at degrees 2 and 3. Every rule is checked for
    strictly positive weights, normalized weight sum, nonnegative barycentric
    coordinates bounded by one, and barycentric coordinate sum equal to one.
    Exactness is checked against the analytic normalized-simplex moment for
    every Cartesian monomial of total degree up to the advertised degree,
    including all mixed terms. All eight quadrature cases pass in both the
    objective-enabled and release configurations.

24. **Done (2026-08-19).** Enabled objective quadrature orders 2 and 3 in the
    shared runtime selector and parameter validation, and updated the
    command-line documentation. The selector returns the corresponding
    positive rule exactly and continues to reject negative orders and
    unsupported orders 1 and 4 without fallback. Because the common value and
    differentiated traversals already consume this selector, the change
    activates the rules
    for SizeShape and every StepDistance variant without objective-specific
    integration edits. Added command-line, parameter-validation, and exact
    selector-equivalence coverage. Seven self-contained focused suites totaling
    29 cases pass directly in both the objective-enabled and release builds;
    order 0 remains the default and the Classical regression is unchanged.

25. **Done (2026-08-19).** Added a systematic 48-configuration integration
    matrix covering SizeShape, original StepDistance, StepDistance ShapeVolume,
    and StepDistance CavityTargetAverage; quadrature orders 0, 2, and 3; varying
    FE and analytical target metrics; and triangles and tetrahedra. Each case
    independently selects the explicit rule table, captures the baseline
    quadrature samples, and freezes both target metric and `theta` according to
    the paper's derivative convention. It checks the parameter-selected
    elemental value and every local node's gradient and packed Hessian against
    the frozen reconstruction, then checks gradients and Hessians by centered
    finite differences of the frozen functional, including active additive
    barrier derivatives for original StepDistance. The expanded suite contains
    13 cases and passes 19,266 assertions in the objective-enabled build and
    19,286 assertions in the release build.

26. **Done (2026-08-19).** Imported and attributed the positive SANS degree-5
    rules (7 triangle points and 14 tetrahedron points), enabled order 5 in the
    shared selector, parameter validation, and command-line interface, and
    extended positivity, normalization, barycentric-validity, and complete
    degree-5 monomial exactness tests. Then added the actual positive SANS
    degree-4 triangle rule (6 points) and exposed tetrahedron order 4 by reusing
    the degree-5 rule, matching SANS where the separate order-4 table is
    disabled because it is identical to the order-5 rule. Enabled order 4 in
    the same selector and interfaces and added complete degree-4 contract and
    exactness coverage, including an explicit tetrahedron order-4/order-5
    identity check. Extended the objective integration matrix to 80
    configurations by adding orders 4 and 5. Added diagnostic comparisons
    using a positive 8-point-per-axis Duffy-Gauss reference, damped-Newton
    moving-vertex experiments, and release timing without timing-based
    assertions. Degree 4 sharply reduced the triangle errors while using the
    same six points as degree 3; tetrahedron orders 4 and 5 are necessarily
    identical in value and cost. Every optimization run decreased its selected
    objective and reached a closely grouped vertex location. The complete
    parameter, quadrature, and objective suites pass in both objective-enabled
    and release builds. Detailed measurements are recorded in
    `/tmp/metris_objective_quadrature_order_comparison.md`.

27. **Done (2026-08-19).** Selected the balanced dimension-dependent defaults:
    degree 4 for triangles and degree 3 for tetrahedra. Changed the stored
    parameter default to `-1`, with the explicit meaning “automatic,” and made
    the shared selector resolve it from `tdim`; this lets explicit order `0`
    retain its historical vertex-barycenter meaning. Updated parameter
    validation and command-line help to document `-1`, the two automatic
    choices, and all explicit orders. Added tests proving the default parameter
    state, explicit command-line orders `-1` and `0`, automatic
    triangle/order-4 and tetrahedron/order-3 selection, and rejection of values
    below `-1`. The full parameter, quadrature-contract, and 80-configuration
    objective/derivative suites build and pass in the release configuration.
    Keep order 5 as the explicit high-accuracy 3D choice and order 0 for
    historical reproduction.

28. **Done (2026-08-19).** Removed the compile-time `ONEPOINTQUAL`,
    `TDIM1POINTSQUAL`, and `KEAST4QUAL` switches and every value, derivative,
    diagnostic-name, and dead-code branch associated with them. The former
    one-point Classical P1 calculation is now ordinary, unconditional
    Classical compatibility code, preserving its barycenter evaluation with
    the arithmetic mean of the vertex metrics; the objective-driven path
    remains entirely on the runtime positive-rule selector. Removed the old
    vertex-plus-barycenter experimental implementation and the modified
    11-point Keast implementation rather than leaving disabled code behind.
    Quality-field diagnostic filenames are now rule-neutral. Updated older
    SizeShape regressions that independently reconstruct the historical rule
    to request order `0` explicitly. A repository-wide search finds no
    remaining references to the three switches. Both the release and
    `TESTQUALITYALGO` libraries rebuild; in both configurations, the pinned
    Classical Distortion direct/inverse value and derivative test and the
    focused objective contract, SizeShape integration, and 80-configuration
    value/gradient/Hessian suites all pass.

**Final integration-policy correction (Done, 2026-08-19).** Centralized the
production objective `theta` selection above the pointwise value and derivative
policies. SizeShape, original StepDistance, and StepDistance ShapeVolume now
all use physical measure when `INTQUALINRIEMSPACE` is absent and physical
measure multiplied by `sqrt(det(M))` when it is present. CavityTargetAverage
deliberately remains the sole `ReferenceAverage` exception, and the
incompatible ShapeVolume-plus-CavityTargetAverage combination is still
rejected. Both common traversals call the same selector, so value, gradient,
and Hessian paths cannot choose different conventions. Added a focused
mode-selection regression and updated independent SizeShape and StepDistance
reconstructions for the common convention. The objective-sample, SizeShape
integration, 80-configuration value/gradient/Hessian matrix, and seven focused
StepDistance derivative cases pass in both release and `TESTQUALITYALGO`
configurations.

29. **Future extension (not part of the present implementation).** Offer the
    same real-quadrature option to Classical size-invariant Distortion only as
    an explicit opt-in. Its default must remain the historical sampling, and
    `opt_power == +1` and `opt_power == -1` must retain their existing
    meanings. Do not implement this until it is requested separately.

## Candidate rules from SANS

- Degree 2: 3 triangle points, 4 tetrahedron points.
- Degree 3: 6 triangle points, 8 tetrahedron points.
- Degree 4: 6 triangle points, 14 tetrahedron points; the tetrahedron rule is
  the degree-5 rule exposed as order 4 because SANS identifies them as the same
  rule.
- Degree 5: 7 triangle points, 14 tetrahedron points.

The selected automatic defaults are degree 4 in 2D and degree 3 in 3D. The
historical vertex-barycenter rule remains available explicitly as order 0.

## Rule for staging changes

Each numbered item should be independently reviewable and tested. Before an
objective contract changes, pin the current behavior that must be preserved.
Do not combine the SizeShape definition change, derivative change, and shared
integration-loop refactor in one commit. Do not alter Classical behavior while
refactoring the objective-driven path. The common traversal must be
degree-aware so the subsequent high-order implementation can reuse it without
another architectural split. Similar scalar formulas may be duplicated across
objective implementations deliberately: do not replace them with a shared
power-chain or scalar-transform abstraction.
