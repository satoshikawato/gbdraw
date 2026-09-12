# Final 05A4 disposition

Starting `origin/dev`: `532e423cd8404df8748fbc6d309fbeed735408a1`.
Decision owner: `satoshikawato`. Decision date: 2026-09-12.

The explicit Product Decision Owner receipt is `satoshikawato`'s 2026-09-12
message following their review of [PR #514](https://github.com/satoshikawato/gbdraw/pull/514)
at `972967a6d603d8174d6dbee73bc30a4c06db2665`. The earlier SESSION 05A12-D request
supplied wording for preparation; this post-review message supplies the receipt.
This record preserves the approved Choice B semantics and the legend contract
established in SESSION 05A11-P. The retained Q01 observation
is a known performance limitation, not a release-blocking correctness defect.
All known 05A4 findings now have a final disposition; confirmed implementation
bugs remain closed. This is a disposition record, not a new acceptance run.

## Generated legend contract: 05A4-04

**CLOSED — EXPECTED BEHAVIOR / NOT A BUG.** Generated legends partition features
by specific color rules actually used in the current diagram:

| Case | Current membership | Generated caption(s) |
|---|---|---|
| L-C1 | All CDS use the default color | `CDS` |
| L-C2 | At least one specific CDS rule is used and default CDS remain | Specific caption(s) + `other proteins` |
| L-C3 | A specific CDS rule exists but matches no current feature | `CDS` |
| L-C4 | Remove/undo the used specific rule and regenerate | Default `CDS` semantics return |
| Non-CDS control | A specific gene rule is used and default genes remain | Specific caption(s) + `other genes` |

This is intended semantic partitioning, not stale legend state. The existing
[feature-rule reference](../REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md)
owns the public explanation. The small
`test_legend_partition_tracks_actual_rule_membership` cases in
[`tests/test_feature_visibility.py`](../../tests/test_feature_visibility.py)
exercise `preprocess_color_tables`, `precompute_used_color_rules`, and
`prepare_legend_table` on biological features. Each case recomputes the legend
after removing and restoring rules; L-C2 covers the L-C4 return to `CDS`.
Renderer behavior and legend ownership are unchanged.

## Q01 Product Decision receipt

**05A4-Q01: DEFERRED / NON-BLOCKING FOR v0.14.**

```text
PRODUCT_DECISION
Concern: performance.circular.external_labels
Scenario revision: 2
Choice: B defer_dense_performance_keep_controls
Rationale:
  The retained bacterial all-feature-label workloads are diagnostic stress
  cases and do not represent a useful mandatory figure-production target.
  The product should not introduce an arbitrary label-count boundary,
  automatic timeout, or numerical SLO from those measurements.
  v0.14 may ship without optimizing this stress case.
  Future optimization should target useful, readable figures and must
  improve performance without degrading label placement quality.
Must preserve:
  - existing outer-label controls;
  - explicit user label selections;
  - complete Result admission;
  - no silent label dropping;
  - Generate progress and Cancel;
  - prior Result preservation on failure/cancel;
  - canonical request/session semantics;
  - export semantics;
  - local browser execution;
  - source and feature identity.
May retire:
  - treating the retained bacterial all-feature-label adversarial workload
    as a v0.14 release-blocking acceptance requirement.
Accepted residual risk:
  A user may explicitly request a very dense external-label layout that
  requires a very long computation until completion or cancellation.
  No completion-time guarantee is introduced.
Owner: satoshikawato
Decision date: 2026-09-12
```

There is no numerical performance SLO or label-count boundary for Q01, and no
automatic timeout is introduced. The 100-label-per-record / 1,100-label
diagnostic ladder sets neither a supported maximum nor a Product limit.
The 1,800-second observation window establishes neither acceptable latency nor
a timeout. This receipt does not declare large Circular diagrams unsupported,
dense placement fast enough, or any public support guarantee. Users retain the
ability to request dense outer labels.

Bacterial whole-genome/all-feature-label jobs are not routine tests or release
acceptance requirements. Normal correctness, admission, cancellation and
recovery requirements remain in force. No performance or scientific-output
gate is waived, and no existing CI test or policy is changed.

## Decision scope and evidence

The [Product Impact Ratchet](PRODUCT_IMPACT_RATCHET.md) supplies the human receipt
format and procedural intake. On the starting SHA, the Product Impact map has
no Circular label-duration concern, `tools/web-product-decisions.json` has no
active decisions, and the preauthorized Option Integrity Product Contract is
absent. Web runtime and Session contracts retain their existing authority.
The post-review receipt above resolves the release-disposition question left
`PRODUCT_DECISION_REQUIRED` in the SESSION 05A11-P scenario-revision-2 packet.

This is an evidence-specific release disposition with no runtime transition or
new public support contract. It is not eligible for the mapped, exact-head
current-decision JSON route and does not add a `BD-###` or global registry entry.
The PR machine representation therefore keeps `decisions: []`; this note is not
a parallel evaluator or authority for future runtime changes. Human review is
required for Product-decision work under the Product Impact Ratchet and
[Web change policy](WEB_CHANGE_POLICY.md), independently of a passing Gate.

The local `gbdraw_v014_session05a11_p_evidence_2026-09-12` archive remains evidence
input outside this PR. Its report, scenario-revision-2 decision packet and
reviewed `changes.patch` supplied the prior diagnosis and two legend additions.
Historical pending classifications are superseded by this post-review receipt;
the raw archive is preserved without rerunning its measurements.

## Future optimization criteria

Evaluate any future Circular external-label optimization on useful, readable
figures. Performance improvement alone is insufficient: review wall time, text
overlap, leader-line crossing, clipping, label-to-feature association,
requested-label preservation, and visual stability/readability. Placement
quality must be preserved or improved. The 24,945-label stress case is not a
required acceptance fixture. Algorithm selection remains separate future work.

## Final finding dispositions and remaining release work

| Finding | Final disposition |
|---|---|
| 05A4-01 | CLOSED |
| 05A4-02 | CLOSED |
| 05A4-03 | CLOSED |
| 05A4-04 | CLOSED — EXPECTED BEHAVIOR / NOT A BUG |
| 05A4-05 | CLOSED |
| 05A4-06 | CLOSED |
| 05A4-07 | CLOSED |
| 05A4-08 | CLOSED |
| 05A4-09 | CLOSED |
| 05A4-10 | CLOSED |
| 05A4-11 | CLOSED |
| 05A4-Q01 | DEFERRED / NON-BLOCKING FOR v0.14 |

The confirmed-bug closures are inherited, not broadly retested here. Production
runtime changes: none. Session schema 41, canonical request schema 7, feature
catalog schema 3, Web file binding schema 2, CI policy, and PR smoke budget 10
remain unchanged. No generated artifact or reference output changes.

Full adversarial re-rerun: **NOT STARTED**. S11 reacceptance: **NOT READY**.
S12: **BLOCKED**. `TECHNICAL_BASELINE_SHA` replacement: **NOT ASSIGNED**.
`PUBLICATION_CANDIDATE_SHA`: **UNASSIGNED**. The next session is the full
48-journey adversarial re-rerun; this session does not start it.
