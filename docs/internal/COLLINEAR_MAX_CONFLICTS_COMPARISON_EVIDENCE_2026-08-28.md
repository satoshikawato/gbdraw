# Collinear `max_conflicts` consumer and comparison evidence

Status: comparison complete; no fresh default selected

## Result

Explicit `max_conflicts=0` and `max_conflicts=1` now reach the Collinear
cluster-merge decision through the typed Linear request path. A fixed
five-anchor fixture produces different block membership at the one-conflict
boundary. Both values retain every input anchor.

This report does not recommend either value. The fresh defaults, Product
Contract, and `PD-OI-007` status are unchanged. `PD-OI-007` remains
`EVIDENCE_REQUIRED`.

## Evidence identity

| Item | Value |
| --- | --- |
| Repository | `satoshikawato/gbdraw` |
| Base branch | `origin/dev` |
| Base SHA | `db19778648c30fbff921951760dd9c1f9c70f3a6` |
| Base merge | PR #421, prerequisite report for the missing consumer |
| Product concern | `diagram-generation.collinear-max-conflicts-default` |
| Product record | `PD-OI-007`, scenario revision `1`, status `EVIDENCE_REQUIRED` |
| gbdraw | `0.14.0b0` |
| Python | `3.13.3` |
| Biopython | `1.85` |
| pandas | `2.3.0` |
| Evidence date | `2026-08-28` |

## Consumer proof

The production path is:

```text
LinearDiagramRequest
  -> LinearDiagramOptions.collinearity_params
  -> build_request_diagram()
  -> _build_linear_diagram()
  -> assemble_linear_diagram_from_records()
  -> build_orthogroup_collinearity_blocks()
  -> cluster_lossless_collinearity_anchors()
  -> _lossless_clusters_can_merge()
  -> params.max_conflicts
```

`test_typed_request_max_conflicts_reaches_real_collinearity_consumer` follows
this path for both explicit values. Only the external protein-search process is
replaced with fixed reciprocal hit tables. The test checks that the same typed
parameter object survives request resolution and that the algorithm returns
the value-specific block counts below.

The consumer counts retained anchors located strictly inside both order
intervals between two compatible cluster boundaries. Merge proceeds when that
count is less than or equal to `max_conflicts`. Retained singleton anchors are
not absorbed or discarded by the merge.

## Product impact preflight

The consumer correction is classified `IMPLEMENT_EXISTING_AUTHORITY`.
`PD-OI-007` requires explicit values `0` and `1` to execute while the fresh
default is pending. The base implementation violated that requirement because
the algorithm ignored both explicit values after validation.

The impact class is `PRODUCT_CHANGE`, authority resolution is
`STATIC_AUTHORITY`, and the corrected head is `CONFORMING`. On inputs with one
interior conflict, an explicit `1` can now merge two compatible clusters where
an explicit `0` keeps them separate. Existing CLI and Web fresh settings use
the numeric value `1`, so their output can change on affected inputs even
though no surface default is edited. Anchor retention is unchanged.

`EVIDENCE_REQUIRED` still governs selection of the fresh cross-surface default.
This pull request does not perform that selection. `PRODUCT_DECISION_REQUIRED`
does not apply to the consumer correction because the Product Contract already
requires both explicit alternatives to execute. `NOT_ALLOWED` does not apply
because the correction preserves both alternatives and changes no authority,
default, or compatibility promise.

## Deterministic fixtures

All fixtures use two record indexes, fixed anchor order, `min_anchors=1`,
`max_diagonal_drift=0`, and `merge_orientation="strand"`. The only comparison
variable is `max_conflicts`. The direct fixtures and typed-request fixture use
the same algorithm consumer.

| Fixture | Interior conflicts | `max_conflicts=0` | `max_conflicts=1` | Retained anchors |
| --- | ---: | --- | --- | ---: |
| `zero-interior-conflicts` | `0` | `[q0,q1,q3,q4]`, `[q2]` | `[q0,q1,q3,q4]`, `[q2]` | `5/5` for both |
| `one-interior-conflict` | `1` | `[q0,q1]`, `[q2]`, `[q3,q4]` | `[q0,q1,q3,q4]`, `[q2]` | `5/5` for both |
| `two-interior-conflicts` | `2` | `[q0,q1]`, `[q2]`, `[q3]`, `[q4,q5]` | `[q0,q1]`, `[q2]`, `[q3]`, `[q4,q5]` | `6/6` for both |

The one-conflict fixture is also run through `LinearDiagramRequest`. Its two
records each contain five CDS features. The third subject CDS uses the reverse
strand, and fixed reciprocal search rows pair features by order.

The observed difference is limited to cluster membership. At `0`, one interior
singleton prevents the two compatible two-anchor clusters from merging. At
`1`, those clusters merge and the singleton remains a separate block. The
zero- and two-conflict controls place the tested values below and above the
same threshold.

## Repeatability and checks

The focused comparison command was run twice without source changes. Each run
reported `4 passed, 102 deselected`; the normalized block memberships were
identical. The complete Collinear test module reported `106 passed`. No fixture
run raised an error. These small in-memory fixtures were used for correctness,
not performance comparison.

The repository fast suite reported `2953 passed, 17 skipped, 11 deselected`.
Its passing tests included all 16 tracked SVG reference-output comparisons,
the public API/CLI contract snapshot, and the current Session codec checks.
The Web change-policy checker reported `Gate: PASS` and `Review: CLEAR` with
no architecture differential.

```bash
python -m pytest tests/test_collinearity.py -q -k 'max_conflicts'
python -m pytest tests/test_collinearity.py -q
ruff check gbdraw/ tests/test_collinearity.py
python -m pytest tests/ -v -m 'not slow'
python -m pytest tests/test_output_comparison.py::TestOutputComparison -q
node tools/check-web-change-budget.mjs
```

## Scope and architecture review

The change adds one conflict-counting helper under the existing Collinear
algorithm owner and reads the existing typed field at the existing merge
decision. It does not add a semantic owner, production path, compatibility
path, public type, or persisted representation. The reviewed `OE`, `PE`, and
`CB` sets are unchanged.

No GUI, Session, cache, Worker, renderer, workflow, checker, numeric fresh
default, or Product Contract file is changed. The earlier prerequisite report
remains in the repository as the record of why the comparison could not run
before this consumer was restored.

## Interpretation limit

These fixtures establish parameter reachability and the merge-threshold effect
for zero, one, and two interior conflicts. They do not select a fresh default
or claim that either value is preferable for other datasets.
