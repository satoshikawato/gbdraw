# PR-CLOSE-09: Option Integrity recovery completion report

Date: 2026-08-30

Closure base: `dev` merge `b76971110734a53144fd4a63bc4d4839bc11d99e`
(PR #442)

## Conclusion

The recovery baseline and every active Option Integrity implementation outcome
are present on the closure base with sufficient merged focused, required-CI,
browser, production-size, and staging evidence. The only authorized deferral is
the narrow direct non-adjacent display-link capability in `PD-OI-014`; no
implementation `GAP` or `CONFLICT` remains.

The current Gallery publication result is not wholly green. Its exact closure
disposition is:

- Gallery common 9: green;
- Gallery performance projection: green;
- Gallery Vibrio: known non-blocking CI reliability timeout, not green.

The maintainer decision excludes that classified
Vibrio timeout from this recovery's runtime and Product acceptance blockers. It
is not a Product acceptance failure, and this report does not describe Gallery
publication as green.

## Evidence inventory

| Product outcome | Current implementation/evidence | Status | Gap |
|---|---|---|---|
| Reconciled Gallery baseline (`REC-01B`) | PR #429 refreshed the official Vibrio assets and synchronized the independently established `579` count; its focused publication checks, two regenerations, visual review, maximum lane, required PR checks, and post-merge Gallery run `33249580558` passed. | `CONFORMING` | None. |
| Raw option integrity (`PD-OI-001`, `PD-OI-002`, `PD-OI-006`, raw part of `PD-OI-005/016`) | PR #430 made Candidate limit, mode, and Collinear scope agree across validation, LOSATP jobs, raw identity, current Session data, canonical request, provenance, and admission. Focused Python/Node evidence, its maximum Vibrio lane, required PR checks, and post-merge Gallery run `33255434740` passed. | `CONFORMING` | None. |
| Derived option integrity (`PD-OI-003` through `PD-OI-007`, derived part of `PD-OI-014/016`) | PR #431 connected every derived field to the Python/Web consumer, derived identity, current Session, and provenance while preserving raw identity. Focused Python/Node/browser/output evidence passed; PR #433's later maximum lane confirmed derived-only reuse on the integrated path. | `CONFORMING` | None. |
| Web comparison controls (`PD-OI-001` through `PD-OI-007`, `PD-OI-015`, `PD-OI-017`) | PR #433 restored separate analysis controls plus reachable ribbon/curve and height controls. Seven focused browser cases, accessibility and Session evidence, scientific-key invariance, a production-size Vibrio lane, and required PR checks passed. | `CONFORMING` | None. |
| Transactional Session behavior (`PD-OI-008`, `PD-OI-013`, `PD-OI-016`, `PD-OI-017`) | PR #434 added explicit import classification/resolution and full-result admission; PR #435 aligned the two staging expectations with that current contract. Focused state/history/browser failure cases, the designated production-size lane, and both PRs' required checks passed. | `CONFORMING` | None. |
| GUI-unmanaged configuration preservation (`PD-OI-009`, `PD-OI-015`) | PR #436 added the typed bounded override path and disclosure/reset journey; PRs #439 and #440 completed current-writer inventory and current Circular Session projection. Focused Python/Node/browser checks and required PR checks passed. | `CONFORMING` | None. |
| Circular presentation behavior (`PD-OI-010`, `PD-OI-012`, `PD-OI-013`, `PD-OI-015`) | PR #437 delivered deterministic single-record selection, crop, reverse display, label, and subtitle through the existing request/render path. Focused typed, Node, browser, visual, and required-CI evidence passed. | `CONFORMING` | None. |
| Circular canonical record identity preservation | PR #441 preserved saved canonical record keys through first Generate and replacement-input separation; focused request/state/browser and tobacco Gallery parity checks passed. | `CONFORMING` | None. |
| Lazy Circular record discovery preservation | PR #442 kept saved identities pending until the input was actually discovered, retained rollback state, and passed its focused sparse-depth, direct-edit, and tobacco parity checks plus required PR checks. | `CONFORMING` | None. |
| Independent Linear typography (`PD-OI-011`, `PD-OI-015`) | PR #438 retained independent public scale/ruler values with fresh linked UI state, explicit unlink/relink, current Session and History preservation, renderer evidence, scientific-cache invariance, realistic visual review, and required PR checks. | `CONFORMING` | None. |
| Current Session compatibility corrective | PRs #439–#442 corrected current-writer test inventory, inactive-mode full-config projection, Circular canonical record identity, and lazy discovery without changing the Session schema, Gallery assets, expected SVGs, Product authority, or CI. The full `dev` staging run `33300329391` passed on PR #442's merge. | `CONFORMING` | None. |
| Direct arbitrary non-adjacent Similarity/Collinear display links (`PD-OI-014`) | Product authority preserves Similarity all-vs-all evidence, both Collinear evidence scopes, uploaded comparison edges, and adjacent/all supported journeys; only direct arbitrary non-adjacent display links remain outside this program. | `AUTHORIZED_DEFERRED` | No implementation gap in the accepted program. |

No inventory row is `GAP` or `CONFLICT`.

## Integrated acceptance matrix disposition

- `REC-001`–`BASE-005`: satisfied by PRs #427–#430 and their recorded
  recovery, publication, visual, maximum-fixture, and Gallery evidence.
- `RAW-001`–`RAW-007`: satisfied by PR #430 and the unchanged relevant paths
  covered by later required CI and staging.
- `DER-001`–`DER-009`: satisfied by PR #431, with integrated raw-reuse and
  production-size continuation in PR #433.
- `UI-001`–`UI-008`: satisfied by PR #433.
- `SES-001`–`SES-012`: satisfied by PRs #434 and #435, with current Session
  continuation verified by PRs #439–#442.
- `CFG-001`–`CFG-005`: satisfied by PRs #436, #439, and #440.
- `CIR-001`–`CIR-006`: satisfied by PR #437, with canonical identity and lazy
  discovery continuations in PRs #441 and #442.
- `TYP-001`–`TYP-007`: satisfied by PR #438.
- `ARC-001`–`ARC-008` and `CI-001`–`CI-007`: satisfied by the changed-scope
  owner/path reports, compatibility dispositions, focused policy evidence,
  required PR checks, and current full `dev` staging recorded below.
- `CLOSE-001`–`CLOSE-005`: this inventory, the documentation-only closure diff,
  and the current evidence split satisfy the evidence conditions. Final
  admission is complete only when the closure PR's required checks pass; the
  merged PR record is the authority for the final checked head and merge result.

## Merged recovery sequence

| Stage | PR | Merge commit | Result |
|---|---:|---|---|
| Recovery boundary | #427 | `2fd1d5a087bcefb809de83a49d56eb5b27c4096d` | Reverted the two invalid implementation PRs without rewriting history. |
| Gallery cardinality round trip | #428 | `fd286c5e664e31f03e986f96c06c14d70580a963` | Preserved materialized record cardinality. |
| `REC-01B` | #429 | `64ba172f4c6ba9ffc4749983b10f721f4ee2cdcd` | Reconciled the current Vibrio baseline. |
| `OPT-RAW-02` | #430 | `aa2c7db176ed24ffdebe619a637b9c5ca26b5d7c` | Completed raw option conformance. |
| `OPT-DERIVED-03` | #431 | `2f86af97e81abe3f381fc8d9e728f039aeb1f3e2` | Completed derived option conformance. |
| `UI-04` | #433 | `b5192dd9cb9658a2d50d5a6149e8c76ccdc05c33` | Restored comparison controls and appearance. |
| `SESSION-05` | #434 | `9c6846eaca70df1bcafd151ee5eb2b2f20d3f10d` | Preserved imported intent and failure isolation. |
| Session staging correction | #435 | `9cb3c1dbd54c4cc752df356edc6c39f0a28584aa` | Aligned test expectations with the current contract; no runtime change. |
| `CFG-06` | #436 | `63a9f91dbc5b7e86941623a539db0993e1f223f2` | Preserved GUI-unmanaged configuration. |
| `CIRCULAR-07` | #437 | `721f01d8d451a94ff6fc2b9f77db48c266c398ea` | Completed Circular record presentation. |
| `TYPO-08` | #438 | `35cfaadb8005eb4cab41092f96501db25c047af0` | Preserved independent Linear typography. |
| Current-writer staging correction | #439 | `fe43ff9dd104f72c2bd963af732f8aa8138eea56` | Aligned the active-config test inventory; no runtime change. |
| Current Circular Session config correction | #440 | `1233dd3f6787d5dd3c4fb348ff4fa28cef714af9` | Preserved valid active/shared config leaves. |
| Current Circular Session identity correction | #441 | `b1a1e4ab77b34ab3f9c6f8e88b546d24f2186ae0` | Preserved canonical record identity on Generate. |
| Current Circular lazy-discovery correction | #442 | `b76971110734a53144fd4a63bc4d4839bc11d99e` | Preserved lazy discovery while retaining saved identity. |

PR #432 was an unrelated agent-harness documentation change and is not relied
on for closure evidence.

## Evidence relied on

### Focused and production-size evidence

- PR #429: publication/session/manifest checks, focused Vibrio Gallery browser
  check, two independent `16 / 579 / 2` regenerations, asset review, and a
  maximum lane without page close, crash, OOM, or worker error.
- PR #430: typed/raw-key/helper/current-option/Session/admission checks and one
  final maximum lane. It completed two Generates, reused all 47 raw jobs, and
  preserved `16 / 579 / 2`.
- PR #431: 411 focused Python checks plus Node, browser, output, architecture,
  and recipe evidence for derived consumers and identities. PR #433 supplied
  the later integrated maximum proof after the scoped expected-ID correction.
- PR #433: seven comparison-control browser cases and a 16.1-minute maximum
  lane whose derived-only mutation reused all 47 raw jobs.
- PR #434: focused import-resolution, state, History, Session, and four browser
  journeys plus the designated maximum transaction lane.
- PRs #436–#438: focused typed, Node, browser, architecture, Session, cache,
  and realistic visual evidence for configuration, Circular presentation, and
  Linear typography respectively.
- PRs #439–#442: bounded current Session regeneration, tobacco first-Generate
  parity, direct-edit identity, lazy-discovery, and sparse-depth round-trip
  evidence. No Gallery asset, expected SVG, Product authority, schema, or
  workflow was changed.

Every relevant merged PR from #429 through #442 passed both protected required
checks, `Web base policy (trusted base)` and `PR / gate`. PR #432 is unrelated
and omitted from this claim.

### Current `dev` staging

The full staging workflow
[`33300329391`](https://github.com/satoshikawato/gbdraw/actions/runs/33300329391)
passed on `b76971110734a53144fd4a63bc4d4839bc11d99e`. All four Playwright
functional shards and `Dev staging / gate` were successful, along with the
selected core, slow, browser, Gallery, recipe, lint, performance, cache, and
architecture jobs. Skipped PR-only jobs are not substituted or treated as
missing evidence.

### Current Gallery publication

The current Gallery workflow is
[`33300329387`](https://github.com/satoshikawato/gbdraw/actions/runs/33300329387)
on the same merge SHA.

- `Gallery browser (common 9)`: green.
- `Gallery publication performance (projection)`: green.
- `Gallery browser (Vibrio)`: known non-blocking CI reliability timeout, not
  green. Each of the three completed attempts remained at
  `losat-cache-preparation-start` during the second derived-only mutation until
  the 20-minute test timeout, then ended with `page-close`. None reported an
  assertion failure for Circular Session identity, tobacco Gallery parity, or
  SVG semantic parity.

The same official Vibrio job passed immediately beforehand in Gallery run
[`33299443823`](https://github.com/satoshikawato/gbdraw/actions/runs/33299443823)
for PR #441's merge. That prior pass does not make the current Gallery workflow
green; it is relevant only to the classified CI-reliability disposition.

No additional Vibrio rerun, manual Vibrio execution, or timeout change belongs
to this closure.

## Architecture closure summary

| Capability | Semantic owner | Canonical path | Compatibility branch | Remaining debt |
|---|---|---|---|---|
| Raw option integrity | `app/current-option-values.js` validates active Web values and `app/run-analysis.js` consumes raw-job semantics. Duplication was reduced by removing hidden Candidate substitution and invalid fallbacks. | Existing state -> `services/session-request.js` -> `app/run-analysis.js` -> existing LOSAT Worker/helper; path duplication was reduced. | Safe raw-cache miss/recompute when invocation identity differs; no reader or migration added. | None. |
| Derived option integrity | `api/request_render.py` owns typed request/provenance, `analysis/collinearity.py` owns collinearity consumption, and `app/run-analysis.js` owns Web derived execution/admission. Ignored and conflated decisions were removed. | Existing canonical request -> raw evidence -> derived consumer -> admitted Result; one path, reduced duplication. | Existing current Session fields and safe derived miss/recompute; no compatibility reader added. | None. |
| Comparison controls | `app/comparison-ui.js` owns availability and presentation projection; typed defaults and scientific semantics remain outside the UI. Duplication is unchanged. | Existing control state -> `services/session-request.js`; one projection path, unchanged. | Existing bounded missing-value fallback for older current Sessions; no new branch. | None. |
| Transactional Session behavior | Private `services/imported-comparison-intent.js` owns import classification/resolution; `services/config.js` owns committed Session state; `app/run-analysis.js` owns final admission. Partial projection and direct-assignment paths were reduced. | Import classification -> explicit draft resolution -> candidate request -> Worker -> complete Result admission -> committed replacement; one path. | One optional field in the current version-40 writer; no dual writer or branch-only migration. | None. |
| GUI-unmanaged configuration | Typed config/options remain the validity owners; `web_support/config_overrides.py` is the bounded Web adapter. The duplicate unsafe-key predicate was reduced. | Session import -> typed validation -> flat overlay -> GUI-precedence request projection -> existing typed Generate -> current writer. | Existing released readers remain bounded; current writer stays session 40/request 6. | None. |
| Circular presentation and identity | Existing discovery/form state, `services/session-request.js`, typed `RecordInput`/`RegionSpec`, and the Circular definition renderer retain their respective owners. Accepted-but-unused title and duplicate reverse paths were reduced. | Discovery -> stable identity/transform projection -> canonical request -> existing Worker renderer -> admitted Result; one path. | Current Session identity is carried only while discovery targets the same input; no schema field, reader, or compatibility framework added. | None. |
| Linear typography | `app/linear-typography.js` owns only the local linked-control relation; the two renderer fields remain independent semantic owners. The former request alias was removed. | Two existing state values -> `services/session-request.js` -> their separate renderer consumers; one path per public value. | Current Session stores the two values and local link state without a new version or reader. | None. |
| Current Session corrective | `web_support/config_overrides.py`, `services/config.js`, `services/session-request.js`, and `app/run-analysis.js` retain config, import, request, and discovery ownership. Duplication is unchanged. | Current Session import -> typed projection -> pending identity -> real discovery -> canonical Generate; one bounded path. | Existing current-session path only; no new schema, migration, fallback framework, or writer. | None. |

## Compatibility status and residual limits

- The current writer remains session version 40 and canonical request schema 6.
  Existing accepted session versions 27–33 and 39–40 and request schemas 1, 2,
  5, and 6 are unchanged.
- No dependency, Session/cache schema, compatibility reader, dual writer,
  workflow, checker, or branch-protection behavior changes in closure.
- The accepted `PD-OI-001` resource risk remains: an unbounded Candidate request
  can be expensive, while finite values remain exact.
- The accepted `PD-OI-008`, `PD-OI-009`, `PD-OI-010`, `PD-OI-014`, and
  `PD-OI-015` limits remain: some imports require explicit resolution, some
  valid settings are disclosed but not editable, arbitrary Circular
  multi-record subset editing is unavailable, and direct arbitrary
  non-adjacent Similarity/Collinear display links are deferred.
- The current Vibrio CI timeout remains a known non-blocking reliability limit.
  It is neither a Product acceptance failure nor a green Gallery result.

## Documentation closure

The closure corrects only the stale public statements found in the affected
documentation owners:

- `docs/FAQ.md` now states the fresh/Reset Collinear scope is **Adjacent pairs**;
- the two GUI protein-comparison Tutorials and the comparison reference use the
  restored always-reachable appearance controls directly instead of a Pairwise
  mode workaround;
- the GUI Collinear Tutorial no longer instructs users to change an already
  adjacent fresh scope to **Adjacent pairs**;
- the 0.14.0b0 release notes identify current canonical request schema 6 and its
  input-cardinality field.

No public or internal documentation promise was found that requires an
implementation change.

## Promotion-readiness conclusion

The Option Integrity recovery is ready for a separate promotion review. This
conclusion relies on the maintainer's explicit non-blocking classification for
the current Vibrio CI reliability timeout and must travel with the exact
Gallery split above; it is not a statement that Gallery publication is green.
