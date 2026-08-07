# Implementation plan: place Linear comparison controls between record rows

Date: 2026-08-07

Status: planned

## Baseline and history

This plan is based on the `docs_renovation` working tree at `4a9f96e6`.
The working tree already contains unrelated and overlapping uncommitted work,
including changes to `gbdraw/web/js/app/linear-comparisons.js` and
`tests/web/linear-multi-record.playwright.spec.js`. Reinspect those diffs
before implementation and preserve changes that remain correct for the final
design.

The explicit Linear comparison plan introduced by `5b127905` is the current
authority. Before that change, each record card except the last contained the
controls for its comparison to the next record, including the raw LOSAT
filename and download action. That placement made adjacent pairs easy to
locate, but the underlying values were tied to a sequential record and pair
index. The current implementation correctly owns comparisons by stable record
UIDs and directional `edgeKey` values. It moved the controls into separate
consolidated **Adjacent gaps**, **Selected pairs and retained drafts**, and
**Raw LOSAT results** sections.

This work changes the presentation, not the comparison authority. Keep the
edge-based model from the explicit-plan implementation and place each active
or editable comparison at the visual boundary between its endpoint display
rows.

## Product decision

The Linear input flow will read in the same order as the diagram topology:

```text
Display row 1
  Record #1
  Record #2, when the row contains more than one record

Comparison boundary: display row 1 to display row 2
  #1 -> #3  [No comparison | Run LOSAT | Upload BLAST TSV]
  #2 -> #4  [No comparison | Run LOSAT | Upload BLAST TSV]
  Raw LOSAT TSV controls or the uploaded BLAST TSV live in each pair card

Display row 2
  Record #3
  Record #4
```

With the default one-record-per-row layout, the boundary appears directly
between record cards. With a custom multi-record layout, one boundary appears
between display-row groups and contains every comparison edge across that row
boundary. A selected cross-product therefore produces several pair cards in
one boundary rather than a detached list elsewhere in the sidebar.

The top-level **No comparison**, **Run LOSAT**, and **Upload BLAST TSV**
choices remain. Label them as an action for all adjacent gaps so they are
understood as a bulk default. LOSAT program, execution, thread, and threshold
settings remain global because they are shared by the active LOSAT work. Do
not duplicate those settings inside every pair card.

Each pair card owns all pair-specific controls:

- endpoint numbers and readable record names;
- `No comparison`, `Run LOSAT`, and `Upload BLAST TSV`;
- the uploaded BLAST TSV and its retained/inactive state;
- the optional raw LOSAT filename and its retained/inactive state;
- raw-result availability; and
- **Save Raw LOSAT TSV**, disabled until that edge has a downloadable cache
  result.

Remove the consolidated **Adjacent gaps** and **Raw LOSAT results** sections.
Do not restore `linearSeqs[].blast`, `linearSeqs[].losat_filename`, pair-index
lookup, or any other per-record comparison authority.

## Detailed UX contract

### Global comparison action

- Keep the global radio group above the record timeline.
- Add a visible group label such as **Apply to all adjacent gaps** and retain
  the existing radio labels exactly.
- The global actions continue to call `setLinearComparisonGlobalAction()` and
  preserve dormant files and custom filenames.
- When the plan is `selected`, show the global group as a bulk action rather
  than inventing a fourth selected radio value. The currently selected edges
  remain visible in their row boundaries.
- Keep the comparison-resolution error next to the global action so an
  invalid selected plan is reported before generation.

### Record and display-row order

- Derive display rows from the same effective layout used by the comparison
  resolver. When custom layout is disabled, every record receives its current
  positional row.
- Sort display-row groups by row number and keep records within a row in
  `linearSeqs` order.
- Preserve the stable sequence number shown to users. Reordering the DOM for
  display-row grouping must not renumber a record or mutate `linearSeqs`.
- Keep each record card keyed by `seq.uid` so a presentation reorder does not
  recreate file inputs or lose in-memory `File` values.
- Move the **Record Layout** editor before the timeline. A row assignment
  change should immediately move the affected record and comparison boundary
  without requiring generation.

### Comparison boundary and pair cards

- Render a boundary only between neighboring populated display rows.
- Show the row identities in the boundary heading and show exact endpoint
  numbers plus readable record labels in each pair legend.
- Include the positional zipped candidate edges so a user can change one gap
  from the bulk default.
- Also include active selected edges, selected cross-product edges, and
  endpoint-valid retained drafts. Merge entries by `edgeKey`; never render a
  resolved edge and its draft as two pair cards.
- Preserve endpoint direction. A selected lower-to-upper edge belongs to the
  boundary shared by the rows but keeps its directional arrow and `edgeKey`.
- Keep edits keyed by `edgeKey` or persistent edge ID. Array position and DOM
  order must not select a file, filename, cache result, or download.
- Editing one derived adjacent gap continues to materialize a selected plan
  through the existing action. The visual move must not change this behavior.

### Pair-specific payload and result controls

- For `source='upload'`, show the BLAST TSV uploader inside that pair card.
- For `source='losat'`, show one compact **Raw LOSAT TSV** subsection in the
  same card. It contains the filename, automatic default name, retained-name
  action when applicable, cache readiness, and download button.
- For `source='none'`, hide active upload and raw-result controls. If the edge
  retains a dormant file or filename, show a compact retained-input indicator
  with an explicit reuse action; do not reactivate it merely because the card
  moved.
- Continue to look up result availability through
  `linearLosatCacheInfoByEdgeKey` and download through
  `downloadLosatPair(edgeKey, filename)`. Do not fall back to ordinal or record
  index.
- Keep automatic filename generation based on the edge endpoints. A record
  reorder may change displayed numbers but must not attach the result to a
  different edge.

### Selected pairs and retained drafts

- Replace the always-open consolidated editor with a collapsed **Advanced
  pair setup** disclosure.
- Keep **Add**, **Zipped adjacent pairs**, **All adjacent-row pairs
  (cross-product)**, and **Omit all** available there.
- Endpoint-valid selected edges and retained drafts are edited in their row
  boundary. The advanced disclosure must not repeat their source, upload, raw
  filename, or download controls.
- Use the advanced disclosure for drafts that cannot be placed in a boundary,
  including missing, duplicate, same-row, and non-adjacent endpoints. Keep
  endpoint selectors and the existing validation message available so these
  drafts can be repaired or omitted.
- If **Add** creates a valid default pair, place and focus that pair in its
  boundary. Do not leave a second copy in the advanced disclosure.

### Accessibility and narrow layouts

- Use a `fieldset` and `legend` for every pair instead of relying on an
  unlabelled cluster of radios.
- Give filename and download controls edge-specific accessible names, for
  example `Raw LOSAT filename for #1 to #2` and
  `Save Raw LOSAT TSV for #1 to #2`. Keep the visible labels concise.
- Add stable DOM hooks such as `data-linear-display-row`,
  `data-linear-comparison-boundary`, and `data-edge-key`. Documentation and
  browser tests must scope repeated controls to an edge hook instead of using
  `.first`, `.nth()`, or a global count to identify the intended pair.
- At the normal narrow sidebar width, record labels may truncate but the
  endpoint numbers, source choice, filename, and action buttons must remain
  readable. Radios may wrap. The sidebar must not gain horizontal overflow.
- Check keyboard order from a record card into its following comparison
  boundary and then into the next display row.

## Surface ownership and exclusions

| Surface | Planned treatment |
| --- | --- |
| Comparison authority | Keep `linearComparisonPlan`, stable record UIDs, persistent edge IDs, and `edgeKey` unchanged. |
| Presentation projection | Add one pure display-row/boundary projection in `gbdraw/web/js/app/linear-comparisons.js`; expose one computed result from `app-setup.js`. |
| Web markup | Reorder the Linear record/layout markup and render pair controls inside row boundaries in `gbdraw/web/index.html`. |
| LOSAT execution | No behavior change. Continue consuming the frozen `linearComparisonResolution` snapshot. |
| Canonical render request | No schema or materialization change. The same UI state must produce the same comparison descriptors and resources. |
| Session/history/reset | No new field and no migration. Version 40 remains current; saved bytes may differ only through unrelated baseline work. |
| Raw LOSAT cache | No cache-key or payload change. UI availability remains an `edgeKey` lookup over the existing cache info. |
| Python API and CLI | Out of scope. Their comparison contracts and help remain unchanged. |
| Renderer and SVG | Out of scope. Do not change geometry, semantic IDs, or tracked reference SVGs. |
| Documentation and Gallery | Update label-dependent prose, executable selectors, affected control screenshots, and one representative inter-record gap capture. |

`gbdraw/web/js/services/session-request.js`,
`gbdraw/web/js/services/config.js`, `gbdraw/web/js/app/run-analysis.js`, and
Python request/session code are inspection-only for this work. Edit them only
if a focused regression proves that the presentation refactor exposed an
existing boundary error. Do not add a session version or compatibility path
for a DOM-only change.

## Implementation work

### 1. Freeze the current edge behavior

Inspect and run the focused comparison tests before editing:

- `tests/web/linear-comparisons.test.mjs`;
- the first four comparison scenarios in
  `tests/web/linear-multi-record.playwright.spec.js`;
- `tests/web/history.test.mjs`;
- `tests/web/session-request.test.mjs`; and
- `tests/web/session-losat-cache-validation.test.mjs`.

Record the current request descriptors, selected-plan transitions, raw
download names, and save/load results. These are the behavior oracle. A failing
baseline caused by overlapping working-tree changes must be resolved or
reported before the UI move is credited with a result.

### 2. Build one pure presentation projection

Edit:

- `gbdraw/web/js/app/linear-comparisons.js`
- `gbdraw/web/js/app/app-setup.js`
- `tests/web/linear-comparisons.test.mjs`

Add a pure helper, named according to the final module vocabulary, that
returns display rows with an optional boundary after each row. Its inputs are
the sequence list, effective record layout, normalized plan, and resolved
comparison snapshot. Its output contains presentation references only:

```js
[
  {
    row: 1,
    records: [{ uid, index, sequence }],
    boundaryAfter: {
      upperRow: 1,
      lowerRow: 2,
      pairs: [
        {
          edgeKey,
          edgeId,
          queryUid,
          subjectUid,
          queryIndex,
          subjectIndex,
          source,
          active,
          candidateKind,
          draft,
          resolved
        }
      ]
    }
  }
]
```

The exact property names may follow existing conventions, but the ownership
rules are fixed:

- reuse the module's row normalization and `adjacentRowPairs()` logic;
- do not reproduce comparison topology in `index.html`;
- produce positional zipped candidates first;
- add resolved and retained selected edges that are not zipped candidates;
- deduplicate by directional `edgeKey`;
- preserve source edge order where the selected plan makes order meaningful;
- classify invalid or duplicate drafts as unplaced rather than dropping them;
  and
- return no mutable or persisted presentation state.

Replace `linearAdjacentComparisonRows` and the view-only use of
`linearResolvedLosatEdges` with this projection. Keep
`linearLosatCacheInfoByEdgeKey` as the cache lookup owner.

Unit-test:

- two default records produce record, boundary, record;
- five default records produce four ordered boundaries;
- two records in one display row and two in the next produce two zipped pair
  cards in one boundary;
- a selected cross-product produces every distinct edge in that same
  boundary;
- a reverse-direction selected edge retains its direction;
- record reorder updates display indexes while preserving `edgeKey` and
  attached draft;
- an inactive retained filename remains on its edge without becoming active;
- duplicate, missing, same-row, and non-adjacent drafts remain available as
  unplaced drafts; and
- one record and `mode='none'` produce no boundary controls.

### 3. Render the record-row timeline

Edit:

- `gbdraw/web/index.html`
- `gbdraw/web/js/app/app-setup.js`
- `tests/test_web_ux_profile.py`

Changes:

- Add a visible label/help text to the global comparison action.
- Move **Record Layout** above the record timeline.
- Replace the single `v-for="(seq, idx) in linearSeqs"` with nested
  display-row and record loops while retaining each record card's current
  content and `seq.uid` key.
- Render the boundary after its row and move the current per-gap source,
  upload, retained payload, raw filename, and raw-download controls into each
  pair card.
- Remove the consolidated **Adjacent gaps** and **Raw LOSAT results** markup.
- Convert the selected-pair area into **Advanced pair setup** and ensure an
  edge has one active editing surface.
- Keep the Add/Remove sequence actions outside the row loop.
- Add the stable accessibility labels and DOM hooks described above.
- Use existing local utility classes or narrowly scoped styles in
  `index.html`; do not add a build step or external dependency.

Do not duplicate the large record-card template for default and custom row
layouts. Both paths must render from the same display-row projection.

### 4. Preserve pair actions, focus, and lifecycle behavior

Edit only as needed:

- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/js/app/linear-comparisons.js`
- `tests/web/linear-multi-record.playwright.spec.js`

Keep the existing mutation owners:

- `setLinearComparisonGlobalAction()`;
- `setLinearComparisonGapAction()`;
- `setLinearComparisonFile()` and its reuse/deactivate actions;
- the resolved and selected raw-filename actions;
- `downloadLosatPair()`; and
- the selected-pair endpoint and batch actions.

If the template currently needs separate resolved-edge and draft actions,
provide a small edge-card adapter that selects the existing action by
`edgeKey`/edge ID. Do not introduce a second mutation path. Focus a newly
created valid pair after Vue renders it, without storing focus state in the
session or history.

Browser coverage must prove:

- the DOM order is record #1, #1-to-#2 boundary, record #2 for the default
  layout;
- a per-gap source edit still materializes the expected selected plan;
- an uploaded TSV remains attached to the correct edge after record reorder;
- the filename and raw download remain attached to the same `edgeKey` after
  reorder and session reload;
- changing display rows moves the pair card to the correct row boundary;
- zipped and cross-product multi-record boundaries show the exact expected
  edge set once;
- `No comparison` hides active payload/result controls without deleting
  retained inputs;
- there is no consolidated **Adjacent gaps** or **Raw LOSAT results** heading;
- the narrow sidebar has no horizontal overflow; and
- keyboard focus follows the visual record/boundary order.

### 5. Update user guidance and executable documentation

Update label-dependent prose:

- `docs/FAQ.md`
- `docs/HOW_TO/GUI/use-uploaded-blast-results.md`
- `docs/HOW_TO/GUI/use-tlosatx.md`
- `docs/TUTORIALS/GUI/compare-genomes-losatn.md`
- `docs/TUTORIALS/GUI/compare-proteins-losatp.md`

Replace instructions that direct the reader to **Adjacent gaps** or **Raw
LOSAT results** with the exact endpoint boundary, for example, "In the
comparison between sequence 1 and sequence 2." Keep exact visible control
labels such as **Raw LOSAT filename** and **Save Raw LOSAT TSV**.

Update documentation capture flows to locate a specific edge container before
its repeated controls:

- `docs/capture/flows/tutorials/gui_losatn.py`
- `docs/capture/flows/bgc_losatp.py`
- `docs/capture/flows/how_to/nucleotide_comparisons.py`
- any shared LOSATP capture helper whose settings screenshot includes the
  changed record timeline
- `tests/test_documentation_capture_contracts.py`
- `tests/test_gui_protein_comparison_capture_contracts.py`
- `tests/test_gui_nucleotide_comparison_capture_contracts.py`

At minimum, inspect and recapture the changed settings/plan images for:

- `T-GUI-03`: `03-losatn-settings.png`;
- `H-GUI-04`: `comparison-plan.png`;
- `H-GUI-05`: `tlosatx-settings.png`;
- `T-GUI-04`: `03-losatp-settings.png`;
- `T-GUI-08`: `03-collinear-settings.png`;
- `H-GUI-07`: `group-settings.png`; and
- `H-GUI-08`: `collinear-settings.png`.

Compare each old and new image at the same displayed size. Keep an image only
when the changed controls are outside its visible viewport and the operation
remains truthful.

### 6. Refresh affected Gallery operation media

Update:

- `docs/internal/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md`
- `gbdraw/web/gallery/tutorials/BGC0000708-BGC0000713.json`
- `gbdraw/web/gallery/tutorials/hepatoplasmataceae_orthogroup.json`
- `gbdraw/web/gallery/tutorials/hepatoplasmataceae_collinear.json`
- `gbdraw/web/gallery/tutorials/majanivirus_orthogroup.json`
- affected files under each tutorial's
  `gbdraw/web/gallery/media/<example-id>/`
- capture-contract tests when selectors or media counts change

Add a new register subsection before replacing media. Recapture the existing
comparison-panel images if the new bulk-action label or crop extent is
visible:

- `BGC0000708-BGC0000713/manual-03-01-open-pairwise.webp`;
- `hepatoplasmataceae_orthogroup/manual-03-01-browser-losat.webp`;
- `hepatoplasmataceae_collinear/manual-03-01-open-pairwise.webp`; and
- `majanivirus_orthogroup/manual-03-01-open-pairwise.webp`.

Add one focused, data-dependent BGC operation image for the first comparison
boundary. It must use the exact BGC session, assert the endpoint record UIDs or
file identities, assert `source='losat'`, and keep the first pair's raw
filename/download controls inside the crop. The caption should tell the reader
that pair-specific evidence and export live between the two record rows.

Inspect the Vibrio **Adjacent pairs** dropdown image, but keep it when its
focused crop does not contain the moved controls. Do not recapture Circular
comparison panels, final previews, thumbnails, source SVGs, or Gallery
sessions solely for this DOM change. The render request and saved session
shape are unchanged.

## Verification

### Focused source and unit checks

```bash
node --check gbdraw/web/js/app/linear-comparisons.js
node --check gbdraw/web/js/app/app-setup.js
node tests/web/linear-comparisons.test.mjs
node tests/web/history.test.mjs
node tests/web/session-request.test.mjs
node tests/web/session-losat-cache-validation.test.mjs
python -m pytest tests/test_web_ux_profile.py -v
python -m pytest tests/test_documentation_capture_contracts.py -v
python -m pytest tests/test_gui_protein_comparison_capture_contracts.py -v
python -m pytest tests/test_gui_nucleotide_comparison_capture_contracts.py -v
```

### Browser checks

Check both Playwright installations first, as required by `AGENTS.md`:

```bash
command -v playwright && playwright --version
python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
node -e "console.log(require.resolve('@playwright/test'))"
```

Run the focused Linear browser spec through Node Playwright when available.
If it is not available, reproduce the same cases with Python Playwright. If
Chromium hits the known sandbox restriction, rerun the same command with the
required sandbox escalation.

```bash
npx --no-install playwright test tests/web/linear-multi-record.playwright.spec.js --project=chromium
```

Manually inspect at least:

- two records, default layout, at the normal desktop sidebar width;
- five records with four LOSAT boundaries;
- two records per display row with zipped selected pairs;
- the same layout with cross-product selected pairs;
- a mixture of LOSAT, uploaded TSV, omitted edge, and retained dormant input;
- record reorder after a raw result is cached; and
- a narrow/mobile viewport with the sidebar scrolled through every boundary.

### Gallery checks

```bash
python tools/capture_gallery_tutorial_screenshots.py --example BGC0000708-BGC0000713 --check
python tools/capture_gallery_tutorial_screenshots.py --example hepatoplasmataceae_orthogroup --check
python tools/capture_gallery_tutorial_screenshots.py --example hepatoplasmataceae_collinear --check
python tools/capture_gallery_tutorial_screenshots.py --example majanivirus_orthogroup --check
python -m pytest tests/test_gallery_capture_contracts.py -v
```

Review every replacement beside its previous image at the same size and check
desktop and mobile Gallery rendering. Confirm that the operated pair, source,
filename, and button are readable in the rendered tutorial, not only in the
full-size WebP.

### Documentation capture checks

After regenerating the affected committed images, verify each scenario:

```bash
python docs/capture/run_all.py --scenario T-GUI-03 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-04 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-05 --tier extended --check
python docs/capture/run_all.py --scenario T-GUI-04 --tier extended --check
python docs/capture/run_all.py --scenario T-GUI-08 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-07 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-08 --tier extended --check
```

These scenarios run real or qualified comparison workflows and may be slow.
Run focused source/selector checks before recapturing them.

### Broader regression gate

```bash
python -m pytest tests/ -v -m "not slow"
```

Review production, test, documentation, Gallery JSON, documentation PNG, and
Gallery WebP diffs as separate groups. `tests/reference_outputs/` remains
read-only. A DOM-only placement change must not require reference SVG
regeneration.

## Acceptance criteria

- In the default Linear layout, every positional comparison control appears
  between its two record cards.
- In a multi-record layout, comparisons appear once in the boundary between
  the two endpoint display rows, including selected cross-product edges.
- A pair's source, upload, raw filename, cache readiness, and raw download are
  edited in one pair card.
- The consolidated **Adjacent gaps** and **Raw LOSAT results** sections are
  gone.
- Shared LOSAT settings remain global and are not repeated per edge.
- Global actions retain their current mode/source semantics and preserve
  dormant payloads.
- Selected-pair editing, zipped batch creation, cross-product creation, omit,
  and repair of invalid drafts remain available without creating two active
  editing surfaces for one edge.
- Reorder and row reassignment change presentation only. Stable UIDs,
  `edgeKey`, files, custom filenames, raw cache results, and downloads remain
  associated with the same biological endpoints.
- The same pre- and post-change state produces the same canonical comparison
  descriptors, resources, LOSAT work, rendered SVG, and exports.
- Session version 40, current writer fields, history/reset semantics, and
  supported migrations remain unchanged.
- CLI, Python API, Circular comparison UI, and renderer geometry remain
  unchanged.
- Repeated controls have edge-specific accessible names and deterministic
  capture hooks.
- The sidebar has no horizontal overflow at the supported narrow width, and
  keyboard order matches visual order.
- Public instructions name the new endpoint-boundary location, executable
  capture scripts select the intended edge rather than a global ordinal, and
  affected screenshots show the real current UI.

## Delivery notes

Implement this as one Web presentation change with its tests, user guidance,
and affected screenshots. Do not refresh saved Gallery sessions or tracked
reference diagrams unless verification identifies a separate semantic change
and that expansion is reviewed explicitly.

Before handoff, re-run `git status --short` and audit overlapping dirty files.
Provide an English commit title and summary, but do not commit, push, deploy,
or publish without explicit authorization.
