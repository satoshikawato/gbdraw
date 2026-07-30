# Implementation plan: explicit LOSAT opt-out

Date: 2026-07-30

Issue:

- [#316: Allow not performing LOSAT searches](https://github.com/satoshikawato/gbdraw/issues/316)

Status: planned

## Scope and product contract

Add `No comparison` as a first-class Web mode for Linear diagrams. This mode
must stop the comparison pipeline before any LOSAT program-specific planning.
It is valid regardless of whether the retained LOSAT settings name LOSATN,
TLOSATX, or LOSATP, including LOSATP Pairwise, Similarity groups, and Collinear
blocks.

Retain the issue's partial-skip use case. Pairwise workflows can select which
record pairs use LOSAT or an uploaded BLAST TSV. An omitted pair performs no
comparison.

The target behavior is:

- `No comparison`: run no LOSAT job and add no comparison resource.
- `All adjacent`: preserve the current behavior for LOSAT or uploaded BLAST
  inputs.
- `Selected pairs`: use only the listed edges. Missing edges mean no
  comparison.
- Selected-pair LOSAT is supported by LOSATN, TLOSATX, and LOSATP when LOSATP
  is in Pairwise mode.
- Similarity groups and Collinear blocks keep their current evidence scope.
  They can be disabled globally with `No comparison`, but they do not accept a
  partial selected-edge plan.

This is a Web planning and orchestration change. The CLI and Python APIs already
support diagrams without comparisons. The canonical renderer also accepts an
empty comparison list.

## Canonical Web state

Replace the current split authority formed by `blastSource`,
`linearRecordLayoutEnabled`, and `linearComparisons` with one current-writer
state:

```js
{
  mode: 'none' | 'adjacent' | 'selected',
  defaultSource: 'losat' | 'upload',
  edges: [
    {
      id,
      queryUid,
      subjectUid,
      source: 'losat' | 'upload',
      file
    }
  ]
}
```

Rules:

- `none` resolves to no edges.
- `adjacent` derives the current adjacent record or adjacent-row edges and uses
  `defaultSource`.
- `selected` resolves only the stored edges.
- Do not store `source: 'none'`. Removing an edge is the canonical
  representation of skipping that pair.
- Stable record UIDs identify endpoints. Array indexes are resolved only at the
  job and request boundaries.
- Record layout controls placement only. It no longer decides whether a
  comparison list is authoritative.
- Changing plan mode, endpoints, program, or source invalidates stale
  `files.linearCanonicalComparisons`. Reusable raw LOSAT cache entries may
  remain keyed by their content identity.

The schema-5 canonical render request remains materialized:

- `mode='none'` becomes `comparisons: []`.
- Uploaded or resolved nucleotide edges become `nucleotideBlast` entries with
  explicit endpoints.
- Selected LOSATP Pairwise edges use the existing explicit generated-pair
  contract.
- No Web comparison-plan field is added to the Python render request.

## Implementation work

### 1. Build one pure comparison-plan resolver

Edit:

- `gbdraw/web/js/app/linear-comparisons.js`
- `gbdraw/web/js/state.js`

Changes:

- Add constants and normalization for `none | adjacent | selected`.
- Add reconciliation by stable record UID.
- Resolve adjacent topology independently of the record-layout enabled flag.
- Return normalized edges with endpoint UIDs, endpoint indexes, a stable edge
  key, source, and a dense render ordinal.
- Add pure validation for duplicate, self, missing, same-row, and non-adjacent
  edges.
- Add a pure Pairwise LOSAT job-spec builder used by all three programs.
- Default new state to `adjacent + losat`, preserving current behavior.
- Keep legacy `blastSource` only in session migration input. Do not retain it as
  a second current-writer authority.

### 2. Make `No comparison` and selected pairs visible in the UI

Edit:

- `gbdraw/web/index.html`
- `gbdraw/web/js/app/app-setup.js`

Changes:

- Add `No comparison` beside the existing `Run LOSAT` and
  `Upload BLAST TSV` choices.
- Treat those choices as bulk actions:
  - `No comparison` selects `mode='none'`.
  - `Run LOSAT` selects `mode='adjacent'` and `defaultSource='losat'`.
  - `Upload BLAST TSV` selects `mode='adjacent'` and
    `defaultSource='upload'`.
- Show per-gap Pairwise controls with `No comparison`, `Run LOSAT`, and
  `Upload BLAST TSV`.
- When one automatically derived gap is changed, materialize the remaining
  adjacent edges and switch to `selected`. Choosing `No comparison` then
  removes only that edge.
- Move the explicit comparison editor outside the
  `linearRecordLayoutEnabled` conditional.
- Keep custom row endpoints and batch-add operations.
- Rename `Run LOSAT protein pair` because selected LOSAT also applies to LOSATN
  and TLOSATX.
- Hide LOSAT settings, raw-result downloads, and comparison-only controls when
  the plan is `none`.
- Reject `selected` with LOSATP Similarity groups or Collinear blocks using a
  direct user-facing error. Do not expand the selection silently.

### 3. Plan jobs only from resolved edges

Edit:

- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/services/losat.js` only if its public job identity lacks the
  stable edge key

Changes:

- Resolve and validate the comparison plan once before LOSAT initialization.
- Set `useLosat` from the resolved edge set, not from record-layout state or a
  global source fallback.
- Return early for `mode='none'` before warming a LOSAT runtime, extracting
  comparison FASTA, reading the comparison cache, or creating jobs.
- In Pairwise workflows, create job specs only for resolved LOSAT edges.
- Process resolved upload edges independently of record-layout state.
- Split the current overloaded `pairIndex` responsibilities:
  - stable edge key for cache and UI identity;
  - endpoint indexes for records;
  - dense ordinal for temporary filenames and render ordering.
- Preserve cancellation, stale-result guards, threaded/serial fallback, and
  cache content validation.
- Leave Orthogroup and Collinear job expansion unchanged when the plan is
  `adjacent + losat`.

### 4. Make the canonical request obey the same plan

Edit:

- `gbdraw/web/js/services/session-request.js`

Changes:

- Pass the resolved plan to `buildComparisons`.
- For `mode='none'`, return no comparison descriptors and ignore stale
  persisted comparison artifacts.
- Remove the fallback that interprets an empty explicit list plus global LOSAT
  state as permission to regenerate protein comparisons.
- Emit generated LOSATP pairs only for resolved selected Pairwise edges.
- Keep LOSATN and TLOSATX output as resolved `nucleotideBlast` resources.
- Ensure `hasLinearComparisonIntent` and annotation loading use the resolved
  plan.

This step fixes the current mismatch where orchestration can decide that no
LOSAT work is required while request projection can still add a generated
protein-comparison descriptor.

### 5. Persist the plan independently of record layout

Edit:

- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/history-snapshot.js`
- `gbdraw/web/js/services/reset.js`
- `gbdraw/web/js/services/gallery-session-migration.js`
- `docs/SESSION_COMPATIBILITY.md`

Changes:

- Write `config.linearComparisonPlan` separately from
  `config.linearRecordLayout`.
- Store edge metadata in the plan and uploaded file payloads in
  `files.linearComparisons`, keyed by stable edge ID.
- Restore and reconcile the complete plan for session load, history undo/redo,
  reset, and canonical artifact adoption.
- Reset to `adjacent + losat`, matching the current default.
- Advance the Web session writer from version 39 to version 40. Version 39 is
  release-backed by tag `0.14.0b1` and requires an explicit reader.
- Add version 39 to the supported migration set and keep one representative
  positive fixture.
- Write only the new independent plan in version 40. Do not continue writing
  the retired nested comparison list.

Version-39 migration:

| Persisted state | Version-40 plan |
| --- | --- |
| Layout disabled and `blastSource='losat'` | `adjacent + losat` |
| Layout enabled with nonempty comparison entries | `selected` with those edges |
| Layout enabled with an empty comparison list | `none` |
| Upload mode with one or more uploaded adjacent files | `selected` with those upload edges |
| Upload mode with no uploaded comparison file | `none` |

Canonical request schema 5 does not change. Update
`docs/SESSION_COMPATIBILITY.md` to list session versions 39 and 40 accurately.

### 6. Refresh generated Gallery sessions and user guidance

Gallery sessions are generator-owned outputs. Promote or regenerate the bundled
Linear sessions with the version-40 writer instead of editing their JSON bytes
by hand.

Check:

- `gbdraw/web/gallery/sessions/`
- `gbdraw/web/gallery/tutorials/`
- `gbdraw/web/gallery/examples.json`

Update:

- `docs/FAQ.md`
- `docs/TUTORIALS/2_Comparative_Genomics.md`
- `docs/TUTORIALS/4_Protein_Comparisons.md`
- `docs/TUTORIALS/7_Linear_Layout.md`

Explain global `No comparison`, selected Pairwise edges, and the unchanged
Orthogroup/Collinear evidence behavior. Update tutorial setup scripts that
write `app.blastSource` to use the comparison-plan API; do not retain a runtime
alias for the retired state. Several Gallery tutorials capture the Pairwise
Comparisons panel; follow
`docs/skills/web-gallery-screenshot-maintenance/SKILL.md` and regenerate only
affected captures.

## Test plan

Pure JavaScript:

- `tests/web/linear-comparisons.test.mjs`
  - Normalize all three modes.
  - Resolve an empty `none` plan.
  - Reconcile, deduplicate, reorder, and remove selected edges by UID.
  - Build exact LOSATN, TLOSATX, and LOSATP Pairwise job specs.
  - Reject selected Similarity groups and Collinear blocks.
- `tests/web/session-request.test.mjs`
  - `none` emits no descriptor or comparison resource.
  - Selected LOSATN/TLOSATX emit only executed nucleotide endpoints.
  - Selected LOSATP Pairwise emits only requested pairs.
  - Stale canonical artifacts do not return after switching to `none`.
  - Canonical projection preserves current plan semantics.
- `tests/web/annotations.test.mjs`
  - Comparison intent is false for `none` and empty selected plans.
- `tests/web/history.test.mjs`
  - Undo/redo restores mode, endpoints, sources, and uploaded files.
- `tests/web/gallery-session-migration.test.mjs`
  - Each version-39 mapping produces the expected current plan.

Browser:

- `tests/web/linear-multi-record.playwright.spec.js`
  - Global `No comparison` generates a diagram without starting LOSAT.
  - Omitting the middle gap from three records runs only the retained edge.
  - Selection works without `Arrange in rows`.
  - Save/load preserves `none` and selected plans.
- `tests/web/losat-cache-migration.playwright.spec.js`
  - Version-39 sessions migrate without changing their previous comparison
    intent.
- `tests/web/gallery-tutorial.playwright.spec.js`
  - Affected tutorials still restore their declared LOSAT state.

Static and Python boundary checks:

- `tests/test_gallery_session_semantics.py`
  - Regenerated Gallery sessions use the current version and retain their
    expected comparison descriptors.
- Existing selected-pair coverage in
  `tests/test_linear_multi_record_comparisons.py` remains green.
- No new Python comparison contract is required.

Focused verification:

```bash
node tests/web/linear-comparisons.test.mjs
node tests/web/session-request.test.mjs
node tests/web/annotations.test.mjs
node tests/web/history.test.mjs
node tests/web/gallery-session-migration.test.mjs
pytest tests/test_gallery_session_semantics.py -v
pytest tests/test_linear_multi_record_comparisons.py -v
node --check gbdraw/web/js/app/linear-comparisons.js
node --check gbdraw/web/js/app/run-analysis.js
node --check gbdraw/web/js/services/session-request.js
node --check gbdraw/web/js/services/config.js
```

Run the focused Playwright specs after confirming the available Node or Python
Playwright installation.

## Acceptance criteria

- Global `No comparison` succeeds with two or more Linear records without
  initializing LOSAT.
- The resulting canonical request has `comparisons: []` and no comparison
  resource.
- No pairwise ribbon or comparison legend entry is rendered.
- Omitting one selected pair causes no FASTA extraction, cache lookup, LOSAT
  job, conversion, or resource for that pair.
- LOSATN, TLOSATX, and LOSATP Pairwise execute exactly the selected edges.
- Similarity groups and Collinear behavior is unchanged unless the global plan
  is `none`.
- Record reorder preserves endpoints by UID. Removing a record removes its
  edges. Invalid non-adjacent topology is reported instead of silently changed.
- Version-39 sessions retain their prior intent after migration.
- Saving and loading `none` never revives the old implicit all-adjacent
  fallback.
- Default new sessions still use all adjacent LOSAT comparisons.

## Delivery and coordination

Implement this as a separate PR from #311/#315. This plan changes more Web and
session ownership code, so merge the scale-visibility PR first and rebase this
branch before final verification. Resolve shared-file conflicts manually in
`gbdraw/web/index.html`, `gbdraw/web/js/state.js`,
`gbdraw/web/js/services/session-request.js`, and their tests.
