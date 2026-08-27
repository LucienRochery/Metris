# Objective and Quadrature Implementation Record

## Scope

This record covers the objective-integration refactor that preceded the
high-order implementation. It defines the pointwise objective contract,
quadrature sample preparation, common value and derivative traversals, runtime
positive quadrature rules, and the production integration policy.

Classical Distortion remains a compatibility path. SizeShape and every
StepDistance variant use the common objective-driven architecture.

## Phase overview

| Phase | Status | Summary |
| --- | --- | --- |
| 1. Compatibility baseline | Complete | Qualified the historical vertex-barycenter rule, introduced a non-owning rule view and production factory, migrated existing P1 consumers, and pinned Classical behavior. |
| 2. Pointwise objective contract | Complete | Introduced runtime `objective_p`, made each objective own its complete `psi` and derivatives, normalized SizeShape and StepDistance semantics, and removed integration-level scalar transforms. |
| 3. Shared integration architecture | Complete | Centralized degree-aware sample preparation and the value and differentiated traversals, unified objective dispatch, and qualified metric-sampling behavior. |
| 4. Runtime positive quadrature | Complete | Imported and verified positive simplex rules through degree five, enabled runtime selection, qualified the objective matrix, and retained vertex-barycenter order `0` as the universal default. |
| 5. Policy cleanup | Complete | Removed obsolete compile-time rule switches and centralized the production integration-density convention. |
| 6. Classical real quadrature | Deferred | A real-quadrature option for Classical Distortion remains an explicit future opt-in; its historical default is unchanged. |

## Final production contract

1. **Pointwise objective.** Each objective returns its complete nonnegative
   `psi(A)`, gradient, and packed Hessian, including its own exponent,
   regularization, ideal-state limits, and rejection policy. The integrator
   never shifts, powers, or otherwise reinterprets this result.
2. **Shared architecture.** SizeShape and every StepDistance variant use one
   quadrature selector, one sample-preparation layer, and common value and
   differentiated traversals. Objective-specific formulas remain in explicit
   lower-level policies.
3. **Rule convention.** Points are barycentric coordinates on the canonical
   reference simplex. Weights are strictly positive and normalized to sum to
   one, preserving the reference-average convention.
4. **Pointwise sampling.** Geometry, metric, and objective are evaluated at
   every quadrature point. An elemental metric is never formed by averaging
   tensor components.
5. **Integrated form.** Each sample contributes
   `w_r [theta_r psi(A_r) + B_r]`. An explicit collapse barrier `B_r` remains
   additive and outside `theta`.
6. **Integration density.** SizeShape, Basic StepDistance, and ShapeVolume use
   physical measure without `INTQUALINRIEMSPACE` and physical measure times
   `sqrt(det(M))` with it. CavityTargetAverage is the sole
   normalized-reference-average exception.
7. **Derivative convention.** The target metric and `theta` are frozen at each
   quadrature sample. Explicit additive-barrier derivatives remain active.
8. **Runtime exponent.** `objective_p` is a run-time `double`, must satisfy
   `objective_p >= 1`, and supports noninteger values. It is not a template
   parameter.
9. **SizeShape definition.** Objective-driven SizeShape uses
   `f_SS(A) = q_SS(A) - 1` and `psi_SS,p = f_SS(A)^p`. It does not support the
   inverse `opt_power == -1` convention.
10. **Classical compatibility.** Classical Distortion preserves both
    `opt_power` signs, historical aggregation, and historical sampling unless
    a future option explicitly requests otherwise.

## Production quadrature rules

| Requested order | Triangle rule | Tetrahedron rule |
| --- | --- | --- |
| `0` | Historical vertex-barycenter rule | Historical vertex-barycenter rule |
| `1` | One-point barycenter rule | One-point barycenter rule |
| `2` | Degree 2, 3 points | Degree 2, 4 points |
| `3` | Degree 3, 6 points | Degree 3, 8 points |
| `4` | Degree 4, 6 points | Degree 5, 14 points, exposed as order 4 |
| `5` | Degree 5, 7 points | Degree 5, 14 points |

The degree-2 through degree-5 rules were imported from SANS. Exact source-file
and revision attribution is retained in the production header. Tetrahedron
orders 4 and 5 intentionally use the same rule because the separate SANS
degree-4 table is disabled as identical.

## Detailed implementation record

### Phase 1: Compatibility baseline

1. **Done — Historical rule tests (2026-08-19).** Standalone triangle and
   tetrahedron tests pin the vertex-barycenter point layout, equal positive
   weights, normalized weight sum, simplex membership, and exact integration
   of constants and affine functions.

2. **Done — Non-owning rule view (2026-08-19).**
   `SimplexQuadratureView<tdim>` and its point view expose a common indexed
   barycentric-point and weight interface without owning or selecting storage.

3. **Done — Production vertex-barycenter factory (2026-08-19).** Static-lifetime
   triangle and tetrahedron tables replaced duplicate test-local construction.
   The descriptive vertex-barycenter name replaced the ambiguous legacy name.

4. **Done — StepDistance value migration (2026-08-19).** P1 StepDistance value
   evaluation obtains its rule points and weights from the shared factory. Its
   integrand and metric-evaluation semantics remain unchanged.

5. **Done — StepDistance derivative migration (2026-08-19).** The differentiated
   path uses the same shared rule. Existing value consistency, rejection,
   2D/3D gradient and Hessian finite differences, and metric-volume barrier
   tests remained unchanged.

6. **Done — SizeShape value migration (2026-08-19).** P1 SizeShape samples the
   metric and complete integrand independently at every vertex-barycenter
   point. Vertex metrics are sampled directly; nonvertex analytical metrics
   use physical evaluation and FE metrics use log-Euclidean barycentric
   interpolation. Tensor-component averaging was removed.

7. **Done — Compatibility audit (2026-08-19).** The dispatch and ownership of
   `opt_power`, `opt_pnorm`, and the StepDistance exponent were mapped before
   changing semantics. A self-contained anisotropic P1 regression pins
   Classical Distortion direct and inverse values, their reciprocal relation,
   and all local gradients and Hessians in Release and objective-enabled
   builds.

### Phase 2: Pointwise objective contract

8. **Done — Shared runtime exponent (2026-08-19).** Added `objective_p` with
   default `1.0`, validation `objective_p >= 1`, and the canonical
   `--objective-p` option. `--step-distance-p` remains a deprecated alias;
   conflicting dual specifications are rejected. Classical Distortion is
   independent of this parameter.

9. **Done — Pointwise result type (2026-08-19).**
   `PointwiseObjectiveResult` owns final `psi`, gradient, packed Hessian, and
   derivative availability. Initial SizeShape and StepDistance adapters
   reproduced the existing pointwise entry points without exposing raw values
   or applying a shared scalar transform.

10. **Done — SizeShape value definition (2026-08-19).**
    `quafun_sizeshape` explicitly computes
    `psi_SS,p = (q_SS - 1)^p`. It clamps only ideal-scale roundoff, rejects a
    substantially negative `q_SS - 1`, and no longer reads `opt_power`.
    Singularity handling became an explicit per-objective policy, preserving
    Classical direct and inverse behavior.

11. **Done — SizeShape derivatives (2026-08-19).**
    `d_quafun_sizeshape` explicitly applies the objective's `p` formulas to the
    derivatives of `q_SS`. At the ideal state, `p = 1` retains the derivatives
    of `q_SS`; every `p > 1` returns zero gradient and Hessian. Triangle and
    tetrahedron finite differences cover `p = 1`, `1.5`, and `2`.

12. **Done — Direct SizeShape integration (2026-08-19).** P1 integration now
    accumulates complete pointwise `psi` values. It no longer applies `difto`,
    absolute value, `opt_pnorm`, or `opt_power`. Independent triangle and
    tetrahedron reconstructions pin historical direct behavior at `p = 1` and
    the new explicit behavior at `p = 2`.

13. **Done — SizeShape differentiated integration (2026-08-19).** Value,
    gradient, and Hessian use the same vertex-barycenter samples, per-sample
    FE or analytical metrics, and normalized weights as the value path. Tests
    verify every local derivative against frozen-sample finite differences and
    isolate caller metric-derivative settings.

14. **Done — StepDistance normalization (2026-08-19).** Every build now selects
    the specialized StepDistance path by `QuaFun::StepDistance`. Value and
    derivatives consume `PointwiseObjectiveResult`; regularized Basic and
    ShapeVolume formulas own their complete exponent algebra. Element and
    mesh aggregation no longer apply `difto` or a second `opt_pnorm` transform.
    ShapeVolume rejection, the additive collapse barrier, and
    CavityTargetAverage's arithmetic mean remain explicit distinct policies.

### Phase 3: Shared integration architecture

15. **Done — Degree-aware sample preparation (2026-08-19).**
    `ObjectiveQuadratureSample<gdim,tdim,mshdeg>` stores barycentric coordinates,
    normalized rule weight, physical coordinates, canonical and regular
    Jacobians, frozen target metric, `theta`, and `theta_is_valid`. Exact rule
    vertices preserve direct nodal sampling; nonvertex metrics follow their FE
    or analytical policy. Reference, physical, physical-metric, and regular-
    reference-metric density modes were made explicit.

16. **Done — Common value traversal (2026-08-19).** The degree-templated value
    traversal owns rule iteration, sample preparation, `theta` validation, and
    weighted accumulation. SizeShape and StepDistance policies supply complete
    `psi`, rejection, density choice, and additive barrier data. ShapeVolume
    rejection precedes invalid-density enforcement, and CavityTargetAverage
    remains a reference average. The accumulator preserves the established
    arithmetic order for ordinary and additive-barrier samples.

17. **Done — Common differentiated traversal (2026-08-19).** The derivative
    traversal owns initialization, derivative-order selection, sample
    preparation, and product-rule accumulation of `theta*psi + B`. It shares
    value-policy density and validity rules and provides degree-aware regular-
    reference basis gradients for Lagrange and Bezier geometry.

18. **Done — Objective/CAD-normal separation (2026-08-19).** The historical CAD
    normal-deviation penalty was removed from every active objective-driven
    value and derivative path. Classical and then-current high-order
    compatibility paths retained their historical behavior.

    **Deferred observation:** repeated SizeShape derivatives on embedded
    `gdim = 3`, `tdim = 2` surfaces can vary and may become nonfinite in
    optimized builds although the scalar value remains stable. This is
    independent of the removed normal penalty and belongs with embedded-surface
    work.

19. **Done — Unified outer dispatch (2026-08-19).** `low_metqua.cxx` and
    `low_metqua_d.cxx` use one objective-driven branch for SizeShape and
    StepDistance and pass `iquaf` to the shared traversal. Objective-specific
    formulas, density, barriers, and rejection remain below that common layer.
    At this stage production activation remained P1; the later high-order work
    reused the degree-aware traversal without another architectural split.

20. **Done — Metric-sampling contracts (2026-08-19).** Varying FE and
    analytical metrics prove that integrating pointwise objectives differs
    from averaging metric tensors first. Constant anisotropic metrics reduce
    the complete value and derivatives to one pointwise result times the weight
    sum. Identity metrics on ideal regular simplices give the expected zero
    SizeShape value and, for `p > 1`, zero derivatives.

### Phase 4: Runtime positive quadrature

21. **Done — Runtime order parameter (2026-08-19).** Added
    `objective_quadrature_order`, its command-line option, validation, and the
    common `get_objective_quadrature` selector. Order `0` selects the historical
    rule exactly.

22. **Done — Degree-2 and degree-3 rule import (2026-08-19).** Positive SANS
    triangle and tetrahedron tables were converted from Cartesian reference
    coordinates to Metris barycentric coordinates while preserving order and
    normalized weights. Source attribution remains in the production header.

23. **Done — Rule contract tests (2026-08-19).** Every imported rule is checked
    for positive weights, normalized sum, valid barycentric coordinates, and
    exact normalized-simplex moments for every monomial through its advertised
    degree, including mixed terms.

24. **Done — Orders 2 and 3 enabled (2026-08-19).** Runtime selection,
    parameter validation, and command-line documentation expose both rules.
    Unsupported orders are rejected without fallback. The common traversals
    activate the rules for every objective without objective-specific edits.

25. **Done — Initial objective matrix (2026-08-19).** A 48-configuration matrix
    covers SizeShape, Basic StepDistance, ShapeVolume, and
    CavityTargetAverage; orders `0`, `2`, and `3`; FE and analytical metrics;
    and triangles and tetrahedra. Values and every local gradient and Hessian
    are checked against explicit frozen-sample reconstruction and centered
    differences, including collapse-barrier derivatives.

26. **Done — Orders 4 and 5 (2026-08-19).** Imported the positive degree-5
    triangle and tetrahedron rules and the degree-4 triangle rule. Tetrahedron
    order 4 deliberately reuses the degree-5 rule. Contract tests include full
    monomial exactness and explicit tetrahedron order-4/order-5 identity. The
    objective matrix expanded to 80 configurations. Positive Duffy-Gauss
    references and moving-vertex experiments confirmed improved accuracy and
    objective decrease; timing remained diagnostic only.

27. **Done — Final default and one-point extension.** Order `0`, the historical
    vertex-barycenter rule, is the default for triangles and tetrahedra. The
    temporary dimension-dependent `-1` selection was removed on 2026-08-27;
    valid runtime orders are now exactly `0` through `5`. Order `5` remains the
    explicit high-accuracy 3D choice.

    On 2026-08-25, runtime order `1` was added as a normalized one-point
    barycenter rule in triangles and tetrahedra. FE metric evaluation at that
    point uses the common log-Euclidean interpolation policy. Parameter
    validation, selector-equivalence tests, and the objective
    value/gradient/Hessian matrix cover the additional order, expanding that
    matrix from 80 to 96 configurations. The default remains order `0`.

### Phase 5: Policy cleanup

28. **Done — Rule-switch removal and density-policy cleanup (2026-08-19).**
    Deleted
    `ONEPOINTQUAL`, `TDIM1POINTSQUAL`, and `KEAST4QUAL` and all associated
    branches. Classical P1 keeps its unconditional historical barycenter
    calculation and arithmetic vertex-metric mean. Objective-driven paths use
    only the runtime rule selector. Dead experimental branches, including the
    modified Keast implementation, were removed rather than left disabled;
    the production vertex-barycenter rule remains order `0`.

    One selector above the pointwise value and derivative policies now chooses
    `theta` for both common traversals. SizeShape, Basic StepDistance, and
    ShapeVolume share the physical/physical-metric convention;
    CavityTargetAverage is the sole `ReferenceAverage` exception. ShapeVolume
    combined with CavityTargetAverage remains invalid. Focused mode-selection
    and independent reconstruction tests verify that value and derivative
    paths cannot diverge.

### Phase 6: Deferred Classical extension

29. **Deferred — Optional real quadrature for Classical Distortion.** If
    requested, expose the same runtime quadrature rules only through an
    explicit opt-in. Historical sampling must remain the default, and
    `opt_power == +1` and `opt_power == -1` must preserve their current
    meanings.

## Qualification map

| Area | Principal focused executables |
| --- | --- |
| Rule views, tables, selection, and sample preparation | `test_quality_quadrature_vertex_barycenter`, `test_objective_quadrature_sample`, `test_objective_parameters` |
| Log-Euclidean FE metric sampling | `test_metric_interpolation_log_euclidean` |
| Pointwise SizeShape | `test_quality_sizeshape_pointwise_value`, `test_quality_sizeshape_pointwise_derivatives` |
| Integrated SizeShape and objective matrix | `test_quality_sizeshape_p1_integration`, `test_quality_sizeshape_p1_integrated_derivatives` |
| Pointwise result contract | `test_pointwise_objective_contract` |
| StepDistance derivatives | `test_stepDistance_SurrealS` |
| Classical compatibility | `test_quality_distortion_classical` |

The focused suites passed in both Release and objective-enabled
`TESTQUALITYALGO` configurations. Timing measurements were diagnostic and were
never used as correctness assertions.

## Explicitly deferred work

- Optional real quadrature for Classical Distortion.
- Embedded-surface SizeShape derivative instability for `gdim = 3`,
  `tdim = 2`.
