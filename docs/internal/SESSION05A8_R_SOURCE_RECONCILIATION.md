# Source replacement reconciliation (05A4-07 and 05A4-10)

Starting `dev`: `cc16bd33dbb7f0fe79c4c68535d3e9eced30634d`.
The remote was fetched and matched this SHA. Work uses a clean isolated worktree;
the original dirty workspace is preserved.

## Diagnosis before production changes

Both retained reproducers fail on this exact source. The original M20 harness
and journey module were copied with only their worktree root changed; its M13
function matches the retained M13 sequence. Original fixtures, snapshots,
downloaded SVGs/sessions, traces, and runtime logs are retained separately in
`/tmp/gbdraw-session05a8-evidence`.

- M13: load the mitochondrial Gallery session, Generate, hide `f24a47546`
  (`s-rRNA`), upload the retained `tobacco.gb`, Generate, Save/fresh Load,
  upload retained `mito.gb`, Generate. The override stays `{f24a47546: off}`;
  manual visibility rules stay empty. The target is absent from tobacco's
  catalog, remains dormant there, survives serialization, and hides the
  mitochondrial feature on return (36 rendered features instead of 37).
- M20: fresh Linear app, retained `lambda.gb`, labels All, Generate; replace
  with retained `tobacco.gb`, Generate; replace with `lambda.gb`, Generate.
  All 73 lambda feature geometries return, but the legend retains
  `repeat_region`, `tRNA`, and `rRNA`. There are no custom legend entries,
  color overrides, or added captions in this reproducer. Original order remains
  `[CDS]` while the editor inventory includes tobacco's categories.
- Both journeys maintain Result/mounted/export agreement and complete Save.
  Neither reports console errors, page errors, or attempted external requests.

These are **independent roots**, so fixes require sequential PRs. Visibility's
candidate mutation planner skips missing rendered targets without pruning its
durable override map; selector-cache preservation retains their rules. Legend
extraction initializes `originalLegendOrder` only once; subsequent generated
categories are then misclassified as manual additions by `candidate-render.js`.
They do not share an existing reconciliation operation whose one correction
would fix both defects.

## Identity and lifetime

Source bytes and ordinary/composite File views belong to
`session-resource-backing.js`; the current `matchesSessionResourceDescriptor`
operation includes 05A4-08's ordered composite/combined-resource equivalence.
Ordinary native Files are bound by the successful run's object identity in
`record-display-options.js`. Resource IDs and file names are not biological
feature identity. Save/Load and CLI imports reconstruct those same supported
backings; schema 1/2 binding readers and the current writer remain unchanged.

Python `features/source.py` captures the original biological catalog before
crop and visibility filtering, using `features/ids.py` identities (record ID,
feature type, original ordered coordinates, strand, and duplicate ordinal).
`services/feature-catalog.js` admits that catalog and its rendered projection.
The existing feature hash utilities distinguish rendered record suffixes from
the underlying identity. Mode navigation only projects the retained catalog;
it is not a successful source replacement.

The Generate caller uses `recordDisplayControls.isCurrentFeature` once per
previous source record to distinguish retained bindings from replacement.
Only a changed binding activates visibility reconciliation. Matching against
the candidate's original catalog then preserves valid targets even when a new
File has a different name. Hidden rendered instances use the visibility
owner's existing cached hash selector when their old rendered row is absent.
No identity bytes, descriptor comparison, or persistent index is duplicated.

| State | Lifetime and disposition |
| --- | --- |
| `featureVisibilityOverrides` | Feature-bound individual intent; discard only absent biological targets after a successful Generate. |
| `featureVisibilityManualRules` | Session-wide matcher rules, with explicit record/type/qualifier scope; preserve. An exact-product rule is a matcher, not one feature's identity. |
| `featureVisibilitySelectorCache` | Generated-artifact-only selector projection; existing refresh follows reconciled overrides. |
| `labelTextFeatureOverrides` | Feature-bound text intent; preserve its existing owner. |
| `labelTextBulkOverrides` | Source-bound source-text matching; preserve its existing context reconciliation. |
| `labelTextFeatureOverrideSources` | Feature-bound source-text context; preserve its existing owner. |
| `labelVisibilityOverrides` | Feature-bound label intent; preserve its existing owner. |
| `featureColorOverrides` | Feature-bound biological-key intent; existing candidate normalization already drops unmatched targets. |
| `featurePlacementOverrides` | Record/feature-bound placement; preserve the recent source-backing capability contract. |
| Catalog and feature selection | Generated-artifact-only catalog projection and mode-bound transient selection. |
| Legend entries/original inventory | Generated-artifact-only category inventory, with separately owned user edits; diagnosis only in the first PR. |
| Legend color/stroke/rename/delete edits | User intent keyed by category; do not treat every legend field as disposable. |
| Layout preferences/legend typography | Session-wide preferences with existing mode-specific slots; preserve. |
| Palette and color rules | Session-wide preferences/matchers; preserve. |
| Generated legend metadata and saved Result | Generated-artifact-only; admitted, mounted, exported, and replayed through existing boundaries. |

Classification for the individual visibility defect is
`IMPLEMENT_EXISTING_AUTHORITY`: the supported “This feature” scope and the
current request's authoritative biological catalog select the target. The
user's explicit mode/source distinction agrees with the retained 05A2-01
contract. Manual matcher lifetimes are not changed. No schema migration,
source fingerprint, second catalog, or reset-all operation is justified.

Generation already captures mutable editor intent in its rollback handle.
Reconciliation must publish only with a successfully admitted replacement and
must remain covered by that rollback and the existing History operation.
The canonical render request continues to describe how its immutable Result
was generated; active editor state determines the next Generate.

## Visibility implementation and verification

The visibility owner gains one reconciliation operation; Generate calls it
after candidate admission, activation, readiness, and cancellation checks.
It removes absent individual targets, refreshes the existing selector cache,
then captures the generated artifact identity. Save and the History checkpoint
therefore receive reconciled active intent. Existing rollback restores both
the map and cache if commitment fails.

The production path remains source binding (`record-display-options.js` and
`session-resource-backing.js`) → admitted catalog (`feature-catalog.js`) →
visibility owner (`feature-visibility.js`), coordinated by `run-analysis.js`
and wired in `app-setup.js`. There was no previous visibility reconciliation
path to retire. Owner/path excess and compatibility paths do not increase;
no architecture exception or persisted-format change is introduced. Rollback
is a normal revert with no data migration.

Local evidence includes the retained original M13 before/after replay, all
V1–V8 cases, 32 focused browser cases covering the closed neighboring defects,
412 focused Python tests, 137 architecture contracts, and 402 fast Web
contracts. The additional duplicate-feature control exposed an early candidate
that pruned intent on unchanged-source regeneration; the source-binding check
and existing selector cache correct it. The final browser set includes repeated
Generate and a renamed File carrying the same duplicate target.

Runtime checks block external requests in fresh browser contexts and use the
locally prepared browser wheel. The original reproducer also compares mounted,
Result, and exported SVG state and records console/page errors. An initial
concurrent browser run lost its shared test server; the corrected run uses a
separately owned local server. A sandbox Node run could not exercise subprocess
contracts normally; executable evidence comes from the escalated rerun.

The trusted-base change gate passes; the added public operation requires human
review under `WEB_CHANGE_POLICY.md`. Hosted and post-merge results are recorded
in the session handoff rather than inferred from local checks. CI tiers,
ten-case PR smoke budget, workflows, timeouts, Gallery artifacts, and unrelated
05A4 findings remain unchanged.

Visibility PR #510 was reviewed by the user at
`d83f46fd86d9b42fdefa288afc5afce2ac44ce42`, passed all PR gates, and merged as
`b67f9e505acd3a3a16320d6caa930bd30c1a2040`. Its final local head passed 403 fast
Web contracts and all 33 selected browser cases. The independent legend branch
starts at that fetched exact `dev`; retained M20 still fails there before
legend production changes.

## Legend implementation and verification

`legend/entry-actions.js` remains the owner of editor legend extraction and its
generated original-category inventory. Its existing extraction scan now
reconciles `originalLegendOrder` with the admitted diagram instead of only
initializing it once. Surviving categories retain their default relative order;
explicitly deleted rows retain their restore intent. New generated categories
are admitted into that same inventory. No second category cache or feature
catalog is introduced.

Manual additions use the existing `data-legend-owner="direct-editor"` marker
at their creation, including when a later renderer produces the same caption.
Extraction excludes those rows from the generated inventory, so they keep their
manual lifetime across Generate and replacement. Renderer categories are keyed
by caption, not biological feature ID: lambda and tobacco have disjoint feature
identities yet legitimately share `CDS`.

Source replacement can also remove a customized category. The existing source
binding predicate, already wired for visibility, supplies one `sourceReplaced`
fact to the candidate mutation planner. Old generated-category styles and
renames may be absent in that replacement. The unchanged-source malformed
binding control remains strict. Explicit deletion is idempotent while its
category is absent and keeps the existing user-facing restore contract.
`svg-result-ingestion.js` applies these plan capabilities using its existing
lazy legend index; it adds no parsing or scan to an empty mutation plan.

This is `IMPLEMENT_EXISTING_AUTHORITY`: current generated-category membership,
the existing manual Add/Restore contract, and the session's source-replacement
invariant determine these outcomes. Palette colors, valid category styles,
legend position, font size, layout preferences, visibility rules, labels, and
placement intent keep their existing owners. There are no persisted fields or
reader/writer changes.

The successful candidate's existing mounted binder updates the inventory, and
the Generate rollback handle already captures all affected legend state. A
rejected replacement therefore preserves the prior diagram and editor state;
History restores the corresponding source, Result, and legend inventory.
There is no watcher-order workaround or additional delay.

The old initialize-only inventory update is replaced in place. Ownership and
production paths remain Generate → existing source binding and candidate plan
→ SVG admission → existing legend extraction; owner/path excess and persisted
compatibility paths do not increase. Normal revert is sufficient for rollback.

The retained M20 passes after the inventory correction. Added browser coverage
exercises L1–L8, common category/different feature identity, manual rows,
customized absent categories, rejected replacement, and Save/fresh Load.
Deterministic tests cover inventory replacement, default-order retention,
manual ownership, deleted-row intent, absent-category admission, and strict
unchanged-source rejection. Local verification passes both original journeys,
the two L1–L8/manual-row browser cases, seven visibility/composition browser
controls, 405 fast Web contracts, and 137 architecture contracts. The policy
result is Gate PASS, Review CLEAR. Final exact-head and hosted results are
recorded in the session handoff.
