# Collinear `max_conflicts` prerequisite evidence

Status: stopped at the PR-03 consumer-proof gate

## Result

The requested comparison cannot proceed on base
`8bc18e12f235ba6798c2ad7021d66bb73688c781`. Explicit `max_conflicts=0` and
`max_conflicts=1` pass typed validation, survive the current request codec,
and produce distinct cache/provenance identities. The Collinear clustering and
merge code does not read the field after validation. Comparing fixture output
would therefore measure the same algorithm twice.

This finding triggers the stop condition in
`PR-03_COMPARE_MAX_CONFLICTS_EVIDENCE_ONLY_INSTRUCTION_PROMPT.md`. No default,
Product authority, runtime selection, expected output, or production code is
changed by this report. `PD-OI-007` remains `EVIDENCE_REQUIRED`.

## Evidence identity

| Item | Value |
| --- | --- |
| Repository | `satoshikawato/gbdraw` |
| Base branch | `origin/dev` |
| Base SHA | `8bc18e12f235ba6798c2ad7021d66bb73688c781` |
| Base merge | PR #420, `Create approved Option Integrity Product Contract` |
| Product concern | `diagram-generation.collinear-max-conflicts-default` |
| Product record | `PD-OI-007`, scenario revision `1`, status `EVIDENCE_REQUIRED` |
| Change class and role | `STANDARD`, `EVIDENCE_ONLY` |
| gbdraw | `0.14.0b0` |
| Python | `3.13.3` |
| Biopython | `1.85` |
| pandas | `2.3.0` |
| Evidence date | `2026-08-28` |

The package versions describe the local evidence environment. The consumer
finding comes from the tracked Python source at the stated base SHA.

## Current paths and defaults

The following inventory records the state before this report was added.

| Surface or stage | Current behavior | Current fresh or omission value |
| --- | --- | --- |
| Typed parameter | `LosslessCollinearityParameters.validate()` accepts every integer greater than or equal to zero. | `max_conflicts=0` |
| Beginner/typed Python execution | `LinearDiagramOptions.collinearity_params` may be `None`; `_resolve_lossless_params()` then constructs `LosslessCollinearityParameters()`. | effective omission value `0` |
| Linear CLI | `_get_args()` sets the hidden `--collinear_max_conflicts_in_merge_gap` argument, then `linear_main()` constructs `LosslessCollinearityParameters`. | `1` |
| Web active config | `createDefaultLosat()` and current-session normalization own `collinearMaxConflictsInMergeGap`. | `1` |
| Web canonical request | `generatedProteinSettings()` writes `collinearityParams.parameters.maxConflicts`. Restore projects it back to `collinearMaxConflictsInMergeGap`. | fallback `1` |
| Web worker/helper | The worker forwards `collinearMaxConflictsInMergeGap`; the Python helper constructs `LosslessCollinearityParameters(max_conflicts=...)`. | fallback `1` |
| Canonical request codec | `_encode_pipeline_value()` and `_decode_collinearity_params()` preserve the typed field. | no independent default; `None` remains omission |
| Derived cache/provenance | `_collinearity_parameter_identity()` records `maxConflicts` and `maxConflictsInMergeGap`. | resolved typed value, normally `0` after Python omission |
| Collinear algorithm | `call_collinearity_blocks()` delegates to `cluster_lossless_collinearity_anchors()`. Validation reads the field; initial clustering and `_lossless_clusters_can_merge()` do not. | no execution effect |

The mismatched surface defaults are observations only. This report does not
select or normalize them.

## Existing fixture inventory

The repository already contains suitable candidates for a later comparison:

- `tests/test_collinearity.py::_anchor` supports small order, orientation, gap,
  and competing-anchor inputs. Existing cases cover no-conflict runs,
  plus/minus orientation, repeated or competing edges, and block-to-diagram row
  conversion.
- `docs/internal/SCENARIO_EVIDENCE.md` scenario `H-CLI-08` uses the five tracked
  `BGC0000708` through `BGC0000713` records for a reproducible Collinear CLI
  figure.
- Gallery entries `hepatoplasmataceae_collinear` and
  `vibrio-harveyi-group-collinear` provide representative multi-record
  Collinear datasets and checked-in reproduction commands.

None was promoted to a PR-03 comparison fixture. Until the field reaches the
consumer, block counts or diagram rows from these datasets cannot answer the
evidence question.

## Consumer proof

An AST inventory of `gbdraw/analysis/collinearity.py` found one
`max_conflicts` attribute read, in
`LosslessCollinearityParameters.validate()` at line 164 on the base SHA. A
post-validation access probe then called
`cluster_lossless_collinearity_anchors()` directly and recorded zero reads.

The dynamic probe used a fixed five-anchor order-sensitive input:

```text
(q0,s0), (q1,s1), (q2,s4), (q3,s3), (q4,s4)
```

All other parameters were fixed at `min_anchors=1`, `max_unit_gap=1`,
`max_diagonal_drift=1`, and `merge_orientation="order"`.

| Fixture identity | Explicit value | Actual consumer proof | Normalized result | Repeatability | Performance/error observation | Interpretation limit |
| --- | ---: | --- | --- | --- | --- | --- |
| `pr03-consumer-probe` | `0` | Typed validation passed; codec round trip passed; provenance recorded `0`; post-validation algorithm reads: `0` | `block_0001 plus [q0:s0,q1:s1]`; `block_0002 minus [q2:s4,q3:s3]`; `block_0003 plus [q4:s4]`; no unblocked anchors | `10/10` identical in each of two unchanged command runs | No error; runtime and memory were not measured because the consumer gate failed | This is a reachability probe, not representative evidence for choosing a default |
| `pr03-consumer-probe` | `1` | Typed validation passed; codec round trip passed; provenance recorded `1`; post-validation algorithm reads: `0` | Same normalized result as explicit `0` | `10/10` identical in each of two unchanged command runs | No error; runtime and memory were not measured because the consumer gate failed | Equality cannot show scientific equivalence when the compared field is unused |

The codec and provenance distinguish the two requests while the algorithm does
not. A cache miss or differing metadata would report requested intent, but the
derived block result would still come from identical clustering semantics.

## Exact commands

Base and authority checks:

```bash
git fetch origin
git rev-parse origin/dev
git log --oneline --decorate -1 origin/dev
git grep -n "PD-OI-007" origin/dev -- docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md
```

Current and historical attribute-use inventories:

```bash
rg -n "max_conflicts" gbdraw tests -g '*.py' -g '*.js' -g '*.json'
git grep -n 'max_conflicts' c727b04c^ -- gbdraw/analysis/collinearity.py tests/test_collinearity.py
```

The following probe was run twice without modification. Each run completed
successfully and produced the same normalized JSON observations.

```bash
python - <<'PY'
from __future__ import annotations

import ast
import json
from pathlib import Path

from gbdraw.analysis.collinearity import (
    CollinearityAnchor,
    LosslessCollinearityParameters,
    call_collinearity_blocks,
    cluster_lossless_collinearity_anchors,
)
from gbdraw.api.request_render import _collinearity_parameter_identity
from gbdraw.session_request_codec import (
    _decode_collinearity_params,
    _encode_pipeline_value,
)


def anchor(query_order: int, subject_order: int) -> CollinearityAnchor:
    return CollinearityAnchor(
        query_protein_id=f"q{query_order}",
        subject_protein_id=f"s{subject_order}",
        query_record_index=0,
        subject_record_index=1,
        query_order=query_order,
        subject_order=subject_order,
        query_start=query_order * 10 + 1,
        query_end=query_order * 10 + 9,
        subject_start=subject_order * 10 + 1,
        subject_end=subject_order * 10 + 9,
        identity=90.0,
        evalue=1e-30,
        bitscore=100.0,
        alignment_length=100,
        query_feature_svg_id=f"fq{query_order}",
        subject_feature_svg_id=f"fs{subject_order}",
        source="pr03-consumer-probe",
        query_unit_id=f"qu{query_order}",
        subject_unit_id=f"su{subject_order}",
        query_unit_kind="cds",
        subject_unit_kind="cds",
        query_locus_id=None,
        subject_locus_id=None,
        query_display_name=f"q{query_order}",
        subject_display_name=f"s{subject_order}",
    )


def normalized(result):
    return {
        "blocks": [
            {
                "id": block.block_id,
                "orientation": block.orientation,
                "members": [
                    [item.query_protein_id, item.subject_protein_id]
                    for item in block.anchors
                ],
            }
            for block in result.blocks
        ],
        "unblocked": [
            [item.query_protein_id, item.subject_protein_id]
            for item in result.unblocked_anchors
        ],
    }


tree = ast.parse(Path("gbdraw/analysis/collinearity.py").read_text(encoding="utf-8"))
parents = {}
for node in ast.walk(tree):
    for child in ast.iter_child_nodes(node):
        parents[child] = node
reads = []
for node in ast.walk(tree):
    if isinstance(node, ast.Attribute) and node.attr == "max_conflicts":
        current = node
        scopes = []
        while current in parents:
            current = parents[current]
            if isinstance(current, (ast.FunctionDef, ast.ClassDef)):
                scopes.append(current.name)
        reads.append({"line": node.lineno, "scope": ".".join(reversed(scopes))})

anchors = tuple(
    anchor(query, subject)
    for query, subject in ((0, 0), (1, 1), (2, 4), (3, 3), (4, 4))
)
observations = {}
for value in (0, 1):
    params = LosslessCollinearityParameters(
        min_anchors=1,
        max_unit_gap=1,
        max_diagonal_drift=1,
        max_conflicts=value,
        merge_orientation="order",
    )
    params.validate()
    encoded = _encode_pipeline_value("collinearity_params", params)
    decoded, legacy = _decode_collinearity_params(
        encoded,
        path="probe.collinearityParams",
        allow_standard=False,
    )
    assert decoded == params and legacy is None
    runs = [normalized(call_collinearity_blocks(anchors, params=params)) for _ in range(10)]
    assert all(run == runs[0] for run in runs)
    observations[str(value)] = {
        "codec": encoded,
        "provenance": _collinearity_parameter_identity(params),
        "normalized": runs[0],
        "repeatability": "10/10 identical",
    }


class AccessProbe:
    min_anchors = 1
    max_unit_gap = 1
    max_diagonal_drift = 1
    merge_orientation = "order"

    def __init__(self):
        self.max_conflicts_reads = 0

    @property
    def max_conflicts(self):
        self.max_conflicts_reads += 1
        return 1

    def validate(self):
        return None


probe = AccessProbe()
cluster_lossless_collinearity_anchors(anchors, params=probe)
assert probe.max_conflicts_reads == 0
assert observations["0"]["normalized"] == observations["1"]["normalized"]
print(json.dumps({
    "algorithmAttributeReads": reads,
    "postValidationConsumerReads": probe.max_conflicts_reads,
    "observations": observations,
}, sort_keys=True, indent=2))
PY
```

## Observation and interpretation boundary

Observed:

- Explicit `0` and `1` are valid typed values.
- Both values survive the current lossless request codec.
- Cache identity and provenance record the requested value.
- The current clustering and merge implementation never reads the value after
  validation.
- The fixed probe is deterministic, but its equality result is caused by the
  missing consumer path.

Not established:

- block merging, conflict consumption, retained anchors, diagram rows, runtime,
  or memory differences between functioning `0` and `1` semantics;
- whether either value violates a named scientific or domain rule;
- which representative fixtures expose meaningful differences; or
- a preferred fresh default.

The repository history contains an older
`_conflicts_between_blocks()`/`_blocks_can_merge()` path before commit
`c727b04c`, but history is evidence rather than current Product or scientific
authority. Restoring that code without reviewing its fit with the current
lossless algorithm would select execution semantics without current proof.

## Smallest correctness prerequisite

A separate implementation PR must make `max_conflicts` affect the current
lossless merge decision before the PR-03 comparison can resume. That PR needs
to:

1. define a conflict at the current block/anchor boundary and count it in one
   canonical Collinear consumer;
2. prove at the consumer that explicit `0` and `1` take distinct merge paths on
   a focused resolvable-conflict fixture;
3. preserve non-negative validation, request/Session round trip, and truthful
   provenance;
4. leave every fresh default and `PD-OI-007` unchanged; and
5. add a regression test that fails if either explicit value becomes inert.

After that prerequisite is merged into `dev`, restart the evidence comparison
from the new base. The restarted work can select the small fixture matrix,
measure representative and Gallery-quality inputs, and report block and
diagram consequences without choosing the Product default.

## Verification

The evidence probe above completed twice with identical observations. Focused
contract tests and repository checks also passed:

```text
pytest tests/test_collinearity.py tests/test_session_request_codec.py tests/test_api_request_render.py -v
261 passed in 5.44s

ruff check gbdraw/
All checks passed!

node tests/web/architecture-contracts.test.mjs
133 passed, 0 failed

node tools/check-web-change-budget.mjs --base origin/dev --head HEAD
Gate: PASS; Review: CLEAR; production additions and deletions: 0

git diff --cached --check
No whitespace errors
```

Reference-output tests were not run. This report stops before choosing a
functioning comparison fixture and neither generates nor changes a diagram or
tracked reference output. Running the reference-output suite here would not
add evidence about the missing consumer path.

## Architecture and rollback

This report is not architecture-bearing. It adds no semantic owner, canonical
path, compatibility path, production dependency, runtime branch, or persisted
format. `OE`, `PE`, and `CB` are unchanged. Rollback is removal of this report;
runtime behavior and Product authority would remain unchanged.
