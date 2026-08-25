# Implementation plan: explicit LOSAT opt-out

Date: 2026-07-31

Issue:

- [#316: Allow not performing LOSAT searches](https://github.com/satoshikawato/gbdraw/issues/316)

Status: planned

> Superseded decision (2026-08-10): [Work package B](WORK_PACKAGE_B_LOSAT_COMPARISON_UI_IMPLEMENTATION_PLAN_2026-08-10.md)
> replaces this plan's new-session and Reset default of `adjacent + losat`.
> Fresh Linear state and **Reset Settings** now use `mode='none'`. Reset still
> retains uploaded BLAST TSV files and custom raw-result names as inactive
> drafts. The supported legacy-session migration rules and the normalizer's
> invalid-input fallback are unchanged.
>
> User amendment (2026-08-11): the Work package B UI projects the existing
> program through three `LOSAT Mode` buttons (LOSATN / LOSATP / TLOSATX) and,
> when LOSATP is active, projects `losat.blastp.mode` through the `LOSATP mode`
> menu (Similarity groups / Collinear blocks / Pairwise matches).
> This does not add persisted state or change this plan's selected-pair
> compatibility rules.
> Similarity groups always uses all-vs-all search evidence and exposes no
> evidence-scope control. Collinear blocks uses `collinearSearchScope`; its
> fresh and Reset default is All records, while a saved explicit Adjacent
> pairs value is restored unchanged.

## Current baseline and recorded decisions

This plan is based on the current `0.14.0b1` branch in
[#318](https://github.com/satoshikawato/gbdraw/pull/318) at `7aab08f8`. The
current implementation already has:

- Web and Python session writers at version 40;
- canonical render request schema 5;
- current-writer file bytes stored once in `resources` and UI file slots bound
  through `webFiles`; and
- comparison authority split among `blastSource`,
  `linearRecordLayoutEnabled`, `linearComparisons`,
  `linearSeqs[].blast`, and `linearSeqs[].losat_filename`.

Version 40 is branch-only at this baseline. Version 39 is evidenced in the
first-parent history of `main`; there is no `0.14.0b1` release tag.

The following decisions are part of this plan:

1. Integrate issue #316 before the current version-40 branch is merged. Keep
   the final writer at version 40, rewrite all branch-owned version-40
   artifacts to the final shape, and do not add a reader or migration for the
   superseded branch-only version-40 shape.
2. Reset Settings restores the active comparison mode to
   `adjacent + losat`, but preserves user-uploaded BLAST TSV files and custom
   LOSAT result filenames as inactive drafts. This keeps the existing promise
   that uploaded files are not removed by Reset Settings.

## Scope and product contract

Add `No comparison` as a first-class Web mode for Linear diagrams. The
comparison plan is independent of record placement and supports three modes:

- `No comparison`: resolve no comparison edge.
- `All adjacent`: derive the current positional adjacent topology.
- `Selected pairs`: use only explicitly stored edges.

The detailed behavior is:

- `none` performs no LOSAT job, consumes no uploaded comparison, and emits no
  active comparison descriptor or resource in the canonical `renderRequest`.
  A saved session may still retain an inactive uploaded TSV in its general
  `resources` inventory.
- `adjacent + losat` runs LOSAT for every positional adjacent edge.
- `adjacent + upload` uses only adjacent gaps that have an uploaded BLAST TSV.
  An empty upload gap remains a valid way to skip that comparison, preserving
  the behavior described in issue #316.
- `selected` uses exactly the stored drafts with `included=true` and allows
  LOSAT, upload, and omitted pairs in one plan. An included upload edge without
  an active file is incomplete input and produces a direct validation error.
- Selected-pair LOSAT supports LOSATN, TLOSATX, and LOSATP when LOSATP is in
  Pairwise mode.
- LOSATP Similarity groups and Collinear blocks retain their current all-record
  evidence expansion for `adjacent + losat`. A selected plan containing a
  LOSAT edge is rejected in either non-Pairwise mode instead of being silently
  expanded. A selected plan containing only uploaded edges remains valid
  because it does not request a LOSAT expansion.
- A valid one-record Linear plan resolves no comparison edge and must not
  revive a generated protein descriptor. A stored self-edge remains invalid.

For multi-row layouts, `adjacent` means the existing positional, or zipped,
pairing between adjacent rows. The existing `All adjacent-row pairs` batch
action remains a cross-product operation that materializes a `selected` plan;
it does not redefine the default adjacent topology.

Changing modes or global sources is non-destructive. Inactive uploaded files
and custom raw-result filenames remain in the editable draft, but the resolver
must never consume a dormant payload. New sessions continue to default to
`adjacent + losat`.

This is a Web planning, persistence, and orchestration change. The CLI and
Python comparison contracts do not change. The canonical renderer already
accepts an empty comparison list.

## Canonical Web state

Replace all current-writer comparison authority with one state object:

```js
{
  mode: 'none' | 'adjacent' | 'selected',
  defaultSource: 'losat' | 'upload',
  edges: [
    {
      id,
      queryUid,
      subjectUid,
      included,
      fileActive,
      losatFilenameActive,
      source: 'losat' | 'upload',
      file,
      losatFilename
    }
  ]
}
```

`file` is the in-memory File value. Current-writer config serializes its edge
identity and activation metadata, while the bytes and binding follow the
`resources`/`webFiles` ownership described below.

Every entry is an endpoint-keyed draft. `included` is its authoritative
selected-mode membership. `fileActive` and `losatFilenameActive` independently
control whether its retained file and custom filename may participate in an
active plan. Separate flags prevent reusing a filename from silently
reactivating a retained TSV, or vice versa.

The pure resolver applies these rules:

- `none` returns no active edges without deleting drafts.
- `adjacent + losat` derives every positional adjacent edge and applies
  `defaultSource='losat'`. It uses a matching custom filename only when that
  draft has `losatFilenameActive=true`; dormant uploads and names do not change
  its output.
- `adjacent + upload` derives the same topology, then activates only gaps with
  both `fileActive=true` and an uploaded file.
- `selected` returns exactly the drafts with `included=true`, in stored order,
  after validation. A draft with `included=false` cannot reappear merely
  because it retains a file.
- Omitting an edge sets `included=false`. A payload-free omitted draft may be
  pruned; a payload-bearing draft remains available for explicit reuse. Do not
  store `source: 'none'`.
- A persistent `id` owns UI editing and file binding. A derived directional
  `edgeKey` owns endpoint/result association and survives record reorder. A
  dense `ordinal` is assigned only at the job and render boundaries.
- Raw LOSAT cache identity remains content-, program-, and argument-based.
  Neither `id`, `edgeKey`, nor `ordinal` becomes part of the reusable raw-cache
  identity.
- Stable record UIDs identify endpoints. Array indexes are derived only at the
  execution and request boundaries.
- Duplicate, self, missing-UID, same-row, and invalid non-adjacent edges are
  rejected under the existing topology rules.

Record layout controls placement only. It no longer decides whether explicit
comparisons are authoritative.

The current writer must stop emitting and reading the following as active
state:

- `blastSource`;
- `config.linearRecordLayout.comparisons`;
- current-writer `linearSeqs[].blast`;
- current-writer `linearSeqs[].losat_filename`; and
- duplicate endpoint or source metadata in a file binding.

Supported older sessions may supply those fields to the direct migration
boundary.

## Ownership boundaries

| Concern | Owner |
| --- | --- |
| Editable comparison mode, endpoints, sources, and custom filenames | `config.linearComparisonPlan` |
| Uploaded BLAST TSV bytes | Deduplicated entries in top-level `resources`; multiple edge bindings may reference one resource |
| Edge-to-file association | Stable edge-ID bindings in `webFiles.bindings.linearComparisons` |
| Last successfully rendered comparison descriptors | Canonical `renderRequest` and its resources |
| Reusable raw LOSAT result | Existing content-addressed LOSAT cache |
| Record placement | `config.linearRecordLayout` without comparison membership |

Edits to the plan mode, membership, endpoints, source, file/name activation, or
LOSAT program invalidate
`files.linearCanonicalComparisons` and edge-associated `losatCacheInfo`.
Reusable raw LOSAT cache entries remain available when their content and
settings still match. Use explicit plan mutation methods as the invalidation
owner; do not add a deep watcher that duplicates mutation behavior.
Here, `files.linearCanonicalComparisons` names in-memory state; it does not
restore the retired top-level session `files` object.

Canonical request schema 5 remains materialized:

- `mode='none'` becomes `comparisons: []`;
- uploaded and resolved LOSATN/TLOSATX edges become `nucleotideBlast`
  descriptors with explicit endpoints; and
- selected LOSATP Pairwise edges use the existing explicit generated-pair
  contract.

No Web comparison-plan field is added to the Python render request.

## Implementation work

### 1. Build one pure comparison-plan resolver

Edit:

- `gbdraw/web/js/app/linear-comparisons.js`
- `gbdraw/web/js/state.js`

Changes:

- Add constants, normalization, and constructors for
  `none | adjacent | selected`.
- Reconcile drafts and selected edges by stable record UID and edge ID.
- Normalize `included`, `fileActive`, and `losatFilenameActive` independently
  so membership and one retained payload never activate another implicitly.
- Resolve adjacent topology independently of the record-layout enabled flag.
- Preserve the current positional adjacent-row algorithm and keep the
  cross-product algorithm as a selected-plan batch operation.
- Return normalized edges with endpoint UIDs, endpoint indexes, `edgeKey`,
  source, payload metadata, and dense `ordinal`.
- Distinguish the valid empty gaps of `adjacent + upload` from incomplete
  stored upload edges in `selected`.
- Add pure validation for duplicate, self, missing, same-row, non-adjacent, and
  missing-upload cases.
- Add a pure Pairwise LOSAT job-spec builder shared by LOSATN, TLOSATX, and
  LOSATP Pairwise.
- Derive `hasLinearComparisonIntent`, active LOSAT intent, and active upload
  intent from the same resolved result.
- Default new state to `adjacent + losat`.
- Move legacy per-record uploaded files and custom raw-result names into edge
  drafts once during supported-session migration.
- Remove runtime aliases for `blastSource`, `linearSeqs[].blast`, and
  `linearSeqs[].losat_filename` after callers and tests use the plan API.

### 2. Make the plan explicit in the UI

Edit:

- `gbdraw/web/index.html`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/js/app/losat-settings.js`

Changes:

- Add `No comparison` beside the current global `Run LOSAT` and
  `Upload BLAST TSV` choices.
- Treat those global choices as mode/source actions:
  - `No comparison` selects `mode='none'`;
  - `Run LOSAT` selects `mode='adjacent'` and
    `defaultSource='losat'`; and
  - `Upload BLAST TSV` selects `mode='adjacent'` and
    `defaultSource='upload'`.
- Preserve dormant edge payloads when a global action changes the active mode
  or source.
- Show retained-but-inactive files and custom filenames distinctly and provide
  separate reuse actions that set `fileActive=true` or
  `losatFilenameActive=true`. A mode/source switch alone must not reactivate a
  payload disabled by Reset Settings.
- Show per-gap `No comparison`, `Run LOSAT`, and `Upload BLAST TSV` actions.
  Editing one derived gap materializes the active adjacent edges as a selected
  plan, applies that edit, and leaves omitted gaps absent.
- Move the explicit comparison editor outside the
  `linearRecordLayoutEnabled` conditional.
- Retain custom endpoints and the two existing batch operations. Make their
  zipped versus cross-product semantics explicit.
- Rename `Run LOSAT protein pair` because selected LOSAT is no longer
  LOSATP-only.
- Associate raw-result download controls and custom filenames with resolved
  LOSAT edges rather than sequential record cards.
- Build default raw-result names from each edge's endpoint records and look up
  `losatCacheInfo` by `edgeKey`, not by `pairIndex`.
- Show LOSAT settings only when the resolved plan contains LOSAT work. Show
  comparison-only controls only when the plan resolves at least one
  comparison.
- Update Pairwise Match controls, comparison legends, `orthogroup_top`, and
  clicked orthogroup action availability to use resolved-plan intent instead
  of `blastSource`.
- Reject a selected LOSAT edge with Similarity groups or Collinear blocks using
  a direct user-facing error. Do not reject an upload-only selected plan.
- Compute job estimates and worker limits in `losat-settings.js` from resolved
  LOSAT work. Keep the current Orthogroup and Collinear expansion calculations
  for `adjacent + losat`.

### 3. Plan comparison work before initializing LOSAT

Edit:

- `gbdraw/web/js/app/run-analysis.js`

Inspect but do not change unless verification exposes a missing generic field
pass-through:

- `gbdraw/web/js/services/losat.js`

Changes:

- Resolve and validate one immutable comparison-plan snapshot before LOSAT
  runtime preparation or comparison file access.
- Derive `useLosat`, upload work, annotation comparison intent, and comparison
  rendering from that snapshot.
- For `none` or an otherwise empty resolved plan, skip only the comparison
  subpipeline before LOSAT warm-up, comparison FASTA extraction, raw-cache
  lookup, uploaded TSV reads, or job construction. Continue the common
  validation, annotation, depth, and diagram-rendering path.
- Build Pairwise job specs only for resolved LOSAT edges.
- Process resolved upload edges independently of record-layout state.
- Carry `edgeKey`, endpoint indexes, and dense `ordinal` on each job/result.
  Use:
  - `edgeKey` for plan, UI, raw-download, and result association;
  - endpoint indexes for record access; and
  - `ordinal` for temporary filenames and render order.
- Keep raw-cache lookup content-addressed. Do not replace its existing
  `cacheKey` with `edgeKey`.
- Replace array-position assumptions in `losatCacheInfo` with edge-keyed
  lookup.
- Preserve cancellation, stale-run guards, serial/threaded fallback, cache
  validation, and successful-render transaction behavior.
- Leave the LOSAT service unchanged if its existing generic job/result spread
  already carries `edgeKey`.

### 4. Make canonical request materialization obey the same snapshot

Edit:

- `gbdraw/web/js/services/session-request.js`
- the `buildCanonicalRenderRequest` call sites in
  `gbdraw/web/js/app/run-analysis.js`,
  `gbdraw/web/js/services/config.js`, and
  `gbdraw/web/js/services/gallery-session-migration.js`

Changes:

- Require a resolved comparison-plan snapshot when materializing a fresh
  Linear canonical request. Do not infer permission from layout or global
  source state.
- Put a hard `none => []` guard before reading uploaded files, persisted
  canonical comparison artifacts, or generated protein metadata.
- Iterate resolved edges once in dense ordinal order. Join uploaded files and
  executed LOSAT outputs to those edges by `edgeKey`; do not append uploads and
  generated results in separate orders.
- Consume a file only when its resolved edge has `source='upload'`. A stale
  file retained on a LOSAT draft must not become a nucleotide descriptor.
- Emit generated LOSATP descriptors only for resolved selected Pairwise edges.
- Remove the fallback that interprets an empty explicit list plus global LOSAT
  state as permission to regenerate a protein comparison.
- Keep LOSATN and TLOSATX outputs as `nucleotideBlast` resources with exact
  endpoints.
- Use the same snapshot for pre-run annotation intent and post-run comparison
  data. `none`, empty selected, and empty `adjacent + upload` plans must not
  preload comparison-dependent annotations.
- Keep the saved canonical `renderRequest` as the last committed render. Later
  editable plan changes must not rewrite it until a new render succeeds.
- When Gallery promotion must materialize a new request, resolve the migrated
  plan once and pass that exact snapshot.

The existing canonical
`generatedProteinComparison.mode='none'` means that a resolved artifact is
being reused. It belongs to a different namespace and must not be interpreted
as the new Web plan's global opt-out.

### 5. Persist one final version-40 plan

Edit:

- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/history-snapshot.js`
- `gbdraw/web/js/services/reset.js`
- `gbdraw/web/js/services/session-resources.js`
- `gbdraw/web/js/services/session-authority.js`
- `gbdraw/web/js/services/gallery-session-migration.js`
- `gbdraw/session_io.py`
- `docs/SESSION_COMPATIBILITY.md`

Changes:

- Keep Web `SESSION_VERSION` and Python `CURRENT_SESSION_VERSION` at 40. Keep
  their supported sets identical: versions 27–33, 39, and the current version
  40.
- Extend the branch-owned version-40 writer in place. Do not add version 41,
  advertise the superseded version-40 shape, or add a compatibility path for
  it.
- Write `config.linearComparisonPlan` separately from
  `config.linearRecordLayout`. The layout object retains placement fields only.
- Remove current-writer `blastSource`, nested comparisons, per-record BLAST
  files, and per-record LOSAT filenames from both config and UI-state output.
- Store uploaded comparison bytes once in top-level `resources`.
- Store only stable edge-ID-to-resource bindings in
  `webFiles.bindings.linearComparisons`. Do not write top-level legacy `files`
  or duplicate endpoint/source metadata in the binding.
- Restore plan metadata first, then attach restored File objects by edge ID.
  Reconcile endpoints only after records have stable UIDs.
- Make current-writer Web authority validation require the independent plan
  when a Web comparison draft is present and reject retired comparison
  authorities. Do not reject a CLI-only session merely because it has no Web
  draft. All branch-owned old-shape version-40 artifacts must be regenerated in
  the same implementation change.
- Sequence regeneration before enabling the final strict version-40 authority
  check, or use a temporary local converter to feed old branch state to the
  final writer. Remove that converter before commit; it must not become a
  shipped version-40 compatibility path.
- Store the plan and file bindings in history so undo/redo restores exact
  membership, sources, filenames, and uploads without duplicating bytes.
- Keep canonical artifact adoption subordinate to the editable plan and the
  committed render boundary.
- Reset active mode/source to `adjacent + losat`. Remove purely generated
  selected-edge drafts, but retain entries needed to bind user-uploaded TSVs or
  custom LOSAT filenames with `included=false`, `fileActive=false`, and
  `losatFilenameActive=false`. Each retained value requires its own explicit
  reuse action before the resolver may consume it. Do not clear results or
  reusable raw-cache entries.
- Update the Python promotion path to write the same final plan and resource
  ownership as the Web writer. Do not clone retired nested comparison fields
  into a version-40 envelope.

### 6. Migrate supported pre-40 sessions directly

Version 39 is the immediate supported predecessor from first-parent `main`.
Versions 32–33 also contain nested explicit comparison state, while earlier
supported Web sessions rely more heavily on global source and per-record
fields. Migrate each supported pre-40 Web comparison draft directly to the
final version-40 plan. CLI-only replay inputs retain their existing
compatibility path and do not gain a synthetic Web draft.

| Supported persisted state | Final version-40 plan |
| --- | --- |
| Layout disabled or absent, source LOSAT | `adjacent + losat`; import custom-name drafts with `losatFilenameActive=true` |
| Layout disabled or absent, source upload | `adjacent + upload`; import only populated gaps with `fileActive=true` |
| Layout enabled with explicit entries | `selected`; set `included=true`; activate each retained file/name independently |
| Layout enabled with an empty explicit list | `none` |
| Missing comparison fields in a context that previously used defaults | `adjacent + losat` |

Additional migration rules:

- Resolve legacy array positions to stable UIDs before constructing edge IDs.
- Import populated `linearSeqs[].blast` and
  `linearSeqs[].losat_filename` values only at this migration boundary so
  user-owned dormant values are not lost. Set each activation flag only when
  the legacy mode treated that value as active. An authoritative empty
  explicit list still maps to `none`, with any retained values inactive.
- Preserve explicit endpoint direction and per-edge source; do not reconstruct
  either from file-list order.
- Map each legacy `linearSeqs[].losat_filename` only to its original positional
  adjacent gap draft. Do not guess a selected edge from list ordinal or attach
  the name to a non-adjacent endpoint.
- Move legacy file bytes into `resources` and bind them by edge ID.
- Do not infer the new plan from
  `generatedProteinComparison.mode='none'`.
- Migrate directly from the source version; do not chain through an old
  branch-only version-40 shape.
- Keep a representative, genuine version-39 positive fixture sourced from
  first-parent `main`. Retain the existing genuine older fixture and add
  focused version-32/33 cases for nested selected comparisons.
- Replace tests that obtain “version 39” by relabeling current output when the
  comparison-state shape matters.

Update `docs/SESSION_COMPATIBILITY.md` to describe the final independent plan,
the `resources`/`webFiles` ownership split, and the unchanged accepted version
set. Do not document the superseded branch-only version-40 shape.

### 7. Refresh all branch-owned sessions and affected Gallery assets

Session refresh is not limited to Linear examples. The current invariant keeps
all bundled sessions on the current writer, so transactionally regenerate:

- all 11 files under `gbdraw/web/gallery/sessions/`; and
- both `tests/test_inputs/*.gbdraw-session.json` files.

Use `tools/refresh_gallery_sessions.py`; do not edit generated JSON or
compressed bytes manually. Run its validation before replacing outputs. When
the normal refresh changes rendered or catalog assets, review:

- `gbdraw/web/gallery/sources/`;
- `gbdraw/web/gallery/examples/`;
- `gbdraw/web/gallery/thumbnails/`; and
- `gbdraw/web/gallery/examples.json`.

Capture the genuine version-39 fixture before replacing current generated
sessions. After every old-shape branch artifact has been rewritten, enable the
strict final-writer validation and rerun the refresh checks with no temporary
old-v40 converter present.

The one-record `lambda_basic_linear` session must finish with zero comparison
descriptors, demonstrating that the retired global fallback cannot create an
edge without two records.

Update tutorial setup and captures:

- replace the direct `app.blastSource` assignment in
  `gbdraw/web/gallery/tutorials/majanivirus_orthogroup.json`;
- refresh the comparison-panel captures in the BGC,
  Hepatoplasmataceae Orthogroup, and Majanivirus tutorials;
- convert or recapture the stale non-declarative comparison-panel image in the
  Hepatoplasmataceae Collinear tutorial; and
- inspect the Vibrio evidence-scope crop, but recapture it only if the changed
  controls are visible.

Do not recapture unrelated Circular Pairwise panels solely because the Linear
selector changed. Follow
`.agents/skills/web-gallery-screenshot-maintenance/SKILL.md` for capture
maintenance and verify declared crop geometry and load checks.

### 8. Update user guidance without changing CLI contracts

Update:

- `docs/FAQ.md`
- `docs/TUTORIALS/2_Comparative_Genomics.md`
- `docs/TUTORIALS/4_Protein_Comparisons.md`
- `docs/TUTORIALS/7_Linear_Layout.md`

Explain:

- global `No comparison`;
- the difference between adjacent LOSAT, sparse adjacent uploads, and selected
  mixed edges;
- selected LOSAT support for LOSATN, TLOSATX, and LOSATP Pairwise; and
- the unchanged all-record scope of LOSATP Similarity groups and Collinear
  blocks.

Keep these as clearly labeled Web notes. Do not change the documented CLI
adjacent-BLAST or Python `pairs` contracts. `docs/RECIPES.md`, README, CLI
Reference, and Python API documentation are out of scope unless implementation
review finds a directly incorrect Web statement.

## Test plan

### Pure JavaScript and orchestration

- `tests/web/linear-comparisons.test.mjs`
  - Normalize all three modes and preserve inactive drafts without treating
    file presence as selected membership.
  - Reactivate a retained file and custom filename independently.
  - Resolve `none` and a one-record adjacent plan to zero edges.
  - Resolve positional adjacent rows separately from cross-product batch
    materialization.
  - Preserve an empty middle gap in `adjacent + upload`.
  - Reconcile reorder and record removal by UID.
  - Build exact LOSATN, TLOSATX, and LOSATP Pairwise job specs.
  - Validate duplicate, self, missing-UID, same-row, non-adjacent, and missing
    selected-upload cases independently.
  - Reject selected LOSAT with Similarity groups or Collinear blocks while
    accepting upload-only selected plans.
- `tests/web/losat-settings.test.mjs`
  - Estimate zero work for `none`.
  - Count only LOSAT edges in a mixed selected plan.
  - Preserve Orthogroup and Collinear expansion counts for
    `adjacent + losat`.
- `tests/web/run-analysis-simple-path.test.mjs`
  - Prove that an empty resolved plan does not initialize Pyodide/LOSAT, read
    comparison files, extract comparison FASTA, inspect the raw cache, or
    create workers/jobs.
  - Prove that normal record, annotation, depth, and render work continues.
- `tests/web/session-request.test.mjs`
  - Make `none` ignore legacy per-record files, explicit edge files, persisted
    canonical comparisons, and generated protein metadata.
  - Join mixed LOSAT/upload outputs by `edgeKey` and emit dense ordinal order.
  - Emit exact LOSATN/TLOSATX endpoints and selected LOSATP Pairwise pairs.
  - Ignore a stale file attached to a LOSAT edge.
  - Emit no generated descriptor for one record.
  - Preserve the last committed request until a new render succeeds.
- `tests/web/annotations.test.mjs`
  - Resolve comparison intent as false for `none`, empty selected, and empty
    adjacent-upload plans.
- `tests/web/history.test.mjs`
  - Restore mode, endpoints, sources, filenames, and uploaded files.
- `tests/web/mode-profiles.test.mjs`
  - Verify Reset Settings clears `included`, `fileActive`, and
    `losatFilenameActive` while retaining the File and custom filename values.
- `tests/web/session-draft-authority.test.mjs`
  - Keep the editable plan independent from the last committed
    `renderRequest`.
- `tests/web/session-resources.test.mjs`
- `tests/web/session-active-files.test.mjs`
  - Store each uploaded comparison once and bind it only by stable edge ID.
- `tests/web/session-authority.test.mjs`
  - Require the final version-40 plan for Web drafts and reject retired
    current-writer Web authorities without rejecting CLI-only sessions.
- `tests/web/session-losat-cache-validation.test.mjs`
  - Keep raw-cache validation content-based after edge identity changes.
- `tests/web/gallery-session-migration.test.mjs`
  - Cover direct migration from supported old global, per-record, and nested
    comparison shapes.
  - Retain serialized dormant files and filenames with their corresponding
    activation flags false.
- `tests/web/gallery-session-regeneration.playwright.spec.js`
  - Verify promoted sessions and regenerated render resources agree.

### Browser and packaging

- `tests/web/linear-multi-record.playwright.spec.js`
  - Global `No comparison` renders without starting LOSAT.
  - Three records with an empty middle upload gap use only the supplied edge.
  - One selected plan mixes LOSAT, upload, and an omitted pair exactly.
  - An uploaded selected edge works with row arrangement disabled.
  - Selected LOSAT works for LOSATN, TLOSATX, and LOSATP Pairwise.
  - Raw-result naming and download lookup survive record reorder.
  - Comparison-only controls and legends are absent for `none`.
  - Save/load preserves `none`, selected plans, and inactive upload drafts.
  - A retained TSV and custom filename can be reused independently without
    another upload.
- `tests/test_web_packaging.py`
  - Replace save/load checks that directly set `app.blastSource` with the plan
    API and verify the final session shape.
- `tests/web/depth-track-session.playwright.spec.js`
  - Replace the ad hoc `app.blastSource='none'` setup with
    `linearComparisonPlan.mode='none'`.
- `tests/web/losat-cache-migration.playwright.spec.js`
- `tests/run_losat_cache_browser_acceptance.py`
  - Keep their session-version constants and plan setup aligned, and verify
    cache reuse without restoring retired runtime aliases.
- `tests/web/gallery-tutorial.playwright.spec.js`
  - Restore the declared plan and validate affected captures.

### Python and generated-artifact boundaries

- `tests/test_session_io.py`
  - Keep Web/Python current and supported version sets equal.
  - Promote genuine version-39 state directly to the final version-40 plan.
  - Reject a superseded old-shape current-writer Web envelope.
- `tests/test_session_compat.py`
  - Preserve supported replay behavior without adding a branch-only migration.
- `tests/test_refresh_gallery_sessions.py`
  - Require all 11 Gallery and two tracked test-input sessions to use the final
    current-writer shape.
  - Verify committed session Results match generated source SVGs.
- `tests/test_gallery_session_semantics.py`
  - Preserve BGC, Hepatoplasmataceae, Majanivirus, and Vibrio comparison
    semantics.
  - Require `lambda_basic_linear` to have zero comparison descriptors.
- Existing selected-pair coverage in
  `tests/test_linear_multi_record_comparisons.py` remains green.
- No new Python renderer comparison contract is required.

### Focused verification

```bash
node tests/web/linear-comparisons.test.mjs
node tests/web/losat-settings.test.mjs
node tests/web/run-analysis-simple-path.test.mjs
node tests/web/session-request.test.mjs
node tests/web/annotations.test.mjs
node tests/web/history.test.mjs
node tests/web/mode-profiles.test.mjs
node tests/web/session-draft-authority.test.mjs
node tests/web/session-resources.test.mjs
node tests/web/session-active-files.test.mjs
node tests/web/session-authority.test.mjs
node tests/web/session-losat-cache-validation.test.mjs
node tests/web/gallery-session-migration.test.mjs
pytest tests/test_web_run_analysis_simple_path.py -v
pytest tests/test_session_io.py -v
pytest tests/test_session_compat.py -v
pytest tests/test_refresh_gallery_sessions.py -v
pytest tests/test_gallery_session_semantics.py -v
pytest tests/test_linear_multi_record_comparisons.py -v
node --check gbdraw/web/js/app/linear-comparisons.js
node --check gbdraw/web/js/app/losat-settings.js
node --check gbdraw/web/js/app/run-analysis.js
node --check gbdraw/web/js/services/session-request.js
node --check gbdraw/web/js/services/config.js
```

Run the focused Playwright specs after checking both the Node and Python
Playwright installations described in `AGENTS.md`. If Chromium hits the known
sandbox restriction, rerun the same local check with the required escalation.
After focused checks pass, run the fast Web/Python suite and review generated
session and SVG diffs separately from production-code and test diffs.

## Acceptance criteria

- Global `No comparison` succeeds with any number of Linear records without
  initializing LOSAT or reading a comparison input.
- After a successful `none` render, its canonical request has
  `comparisons: []` and no active comparison resource, ribbon, legend entry,
  or raw-result control. A session may retain an inactive TSV resource without
  exposing it to rendering.
- `adjacent + upload` skips every gap without a file and does not report it as
  invalid.
- A selected upload edge without its required file reports a direct validation
  error.
- One selected plan can execute LOSAT and upload edges while omitting another
  pair, with exact endpoint and render ordering.
- LOSATN, TLOSATX, and LOSATP Pairwise execute exactly the selected LOSAT
  edges.
- A selected plan with a LOSAT edge is rejected for Similarity groups and
  Collinear blocks without silent all-record expansion. Upload-only selected
  plans remain valid.
- Similarity groups and Collinear behavior is otherwise unchanged for
  `adjacent + losat`, and both are disabled by `none`.
- Record reorder preserves endpoints and raw-result association by stable
  identity. Removing a record removes active and dormant drafts that reference
  it.
- Mode, membership, endpoint, source, file/name activation, and program edits
  invalidate stale canonical comparison artifacts while leaving reusable
  content-addressed raw cache entries available.
- Reset Settings restores `adjacent + losat` and retains uploaded TSVs and
  custom filenames only as drafts with membership and both activation flags
  disabled. Reusing one retained value does not reactivate the other.
- Saving and loading `none` never revives legacy per-record files, nested
  comparisons, or the implicit generated-protein fallback.
- Supported pre-40 sessions retain their previous comparison intent and
  serialized dormant comparison inputs after one direct migration.
- The final current writer remains version 40 in Web and Python. All 13
  branch-owned bundled sessions use the final shape, and no compatibility path
  for the superseded branch-only version-40 shape remains.
- `lambda_basic_linear` contains no comparison descriptor.
- New sessions still default to all positional adjacent LOSAT comparisons.
- CLI and Python comparison behavior and canonical request schema 5 remain
  unchanged.
- Existing cancellation, stale-result, cache-reuse, and successful-render
  transaction tests remain green.

## Delivery and coordination

Integrate issue #316 into the current version-40 branch before that branch is
merged. Do not first publish the superseded branch-only version-40 session
shape. The implementation and the regenerated version-40 artifacts must land
together so the current-writer authority is internally consistent.

If the #311/#315 scale-visibility work is still active, land the smaller
overlapping change first and rebase before final verification. Otherwise,
re-audit the target branch immediately before implementation. Do not assume
that conflicts are limited to the original four files: current shared
ownership includes state, UI, run orchestration, config/session resources,
Gallery promotion, generated sessions, and their tests.

Review production-code, tests, generated sessions, rendered assets, and public
documentation as separate diff groups. This plan does not authorize changing
tracked reference diagrams unless an intentional renderer geometry change is
identified and reviewed.
