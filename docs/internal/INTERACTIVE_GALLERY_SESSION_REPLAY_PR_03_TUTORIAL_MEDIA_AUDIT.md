# PR 3 instruction: Gallery tutorial and media audit after session replay repair

- Parent plan: [Interactive Gallery session replay remediation plan](./INTERACTIVE_GALLERY_SESSION_REPLAY_REMEDIATION_PLAN_2026-08-22.md)
- Suggested branch: `docs/gallery-replay-tutorial-media`
- Base: `origin/dev` after PR 2
- Dependencies: refreshed sessions/assets and passing publication/first-Generate gates from PR 2
- Scope class: Gallery tutorial JSON, executable capture metadata, operation register, referenced WebP media, browser presentation tests
- Generated session/source/example/thumbnail artifacts: read-only in this PR

## 1. Objective

Make every affected Gallery tutorial and screenshot truthful against the repaired checked-in session. Remove capture actions and prose that depend on broken active state, strengthen exact-session assertions, and recapture only the operations whose visible state or result changed.

This PR completes the user-facing remediation. A passing session replay test is not enough if tutorial media still depicts or silently repairs the old state.

## 2. Required reading and starting checks

Before editing, reread:

- `.agents/skills/web-gallery-screenshot-maintenance/SKILL.md`;
- `CLAUDE.md`;
- `gbdraw/web/CLAUDE.md`;
- `docs/internal/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md`;
- `docs/internal/WEB_GALLERY_DIAGRAM_LAYOUT_RECAPTURE_PLAN.md`;
- the five target tutorial JSON files and every referenced media file.

Start the branch and run strict validation:

```bash
git fetch origin
git switch --no-track -c docs/gallery-replay-tutorial-media origin/dev
git status --short

python tools/capture_gallery_tutorial_screenshots.py --all --check --strict
python -m pytest tests/test_gallery_capture_contracts.py -q
npx playwright test tests/web/gallery-tutorial.playwright.spec.js --project=chromium
```

Create a temporary contact sheet for each affected tutorial. Enumerate media from tutorial JSON first; do not use directory contents as the reference inventory.

## 3. Primary target tutorials

Audit these as one semantic set:

- `gbdraw/web/gallery/tutorials/HmmtDNA_ATskew.json`
- `gbdraw/web/gallery/tutorials/majanivirus_orthogroup.json`
- `gbdraw/web/gallery/tutorials/vibrio-harveyi-group-collinear.json`
- `gbdraw/web/gallery/tutorials/WSSV_genome_comparison.json`
- `gbdraw/web/gallery/tutorials/tobacco-chloroplast.json`

Run strict checks for all11 tutorials after the focused audit because PR 2 normalized every Gallery active state.

## 4. Operation register update

Update `docs/internal/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md` before replacing images.

Add a dated section for the session replay remediation. For every affected operation record:

- tutorial and media path;
- old assumption;
- repaired required app state;
- decision: `keep`, `recrop`, `replace`, `add`, or `remove`;
- exact capture assertions;
- final dimensions, DSF, quality, and review status.

Remove or supersede statements that say:

- Vibrio's editable session remains CLI-only with `linearComparisonPlan.mode = none`;
- capture actions may repair the session into adjacent LOSAT state;
- Gallery refresh synchronizes only legend/track/palette through individual helpers.

Replace them with the new publication contract: exact sessions restore their own reproducible active state before any capture action.

Do not rewrite historical completed sections to imply the old behavior never existed. Add a superseding section and adjust directly contradictory current-contract prose.

## 5. Capture metadata requirements

### 5.1 General rule

Every data-dependent operation must:

- use its own example session;
- declare `dataDependent: true`;
- declare exact `capture.session`;
- assert the state values named by the operation;
- list visible controls/text inside the crop;
- fail before capture if the restored session is wrong;
- avoid cross-example reuse for every example-specific file, organism, accession, result, or setting.

A genuinely data-independent control crop may be shared only when it contains no example-specific state and declares `genericMedia: true`. Do not use that escape hatch to weaken an exact-session assertion.

Capture scripts may open panels, scroll, add temporary selectors, or execute the documented UI action. They may not silently repair missing data/settings and then present the repaired state as if the session restored it.

For ordered or movable controls such as track slots, record rows, and comparison pairs, replay the documented setup from reset/default state and compare the resulting order with the repaired saved session. Document every required move explicitly; do not rely on a pre-arranged capture state.

### 5.2 Vibrio comparison operations

The current tutorial contains repeated calls to:

```javascript
app.setLinearComparisonGlobalAction('losat')
app.setLinearComparisonLosatMode('blastp')
app.setLinearComparisonLosatpMode('collinear')
```

Audit each occurrence.

- If the call exists only to compensate for saved `mode: none`, remove it.
- If the operation intentionally demonstrates pressing/selecting the control, first assert that the exact restored session already has the declared adjacent LOSAT/collinear recipe; then perform only the documented idempotent UI action.
- Do not let an after-action `assertAppState` substitute for a restored-state precondition.
- Prefer capture metadata that opens the real panel and shows its already selected value.

Update or replace `test_vibrio_capture_activates_losat_without_changing_the_cli_only_session` in `tests/test_gallery_capture_contracts.py`. The new contract must assert:

- checked-in Vibrio session itself restores adjacent LOSAT intent;
- captures do not turn `none` into `adjacent` as a hidden repair;
- declared collinear mode, search scope, color mode, filters, record rows, and gap are visible/verified.

### 5.3 AT skew

For input/track/result operations assert:

- exact session;
- enabled slot order `features,gc_content,gc_skew,a_skew_2,ticks`;
- `a_skew_2` dinucleotide `AT`;
- positive/negative colors and caption `AT skew`;
- `label_in_tick_out`;
- palette, label controls, and left legend state used by the figure.

The final preview media must show the complete custom stack and readable legend. Do not use a generated preview as the screenshot for a track-input operation.

### 5.4 Majanivirus

Assert before capture:

- bitscore 100 and identity 20;
- committed arrow geometry;
- active similarity-group/LOSATP intent;
- nine record identities/order;
- three custom caption/rule rows including `WSSV-like proteins`;
- generated preview retains committed matches and legend entries.

Comparison settings screenshots must show the real controls. Final preview appears only at the final/result step unless a later operation is specifically inspecting a popup or editor.

### 5.5 WSSV

Assert:

- exact reference and 20 ordered FASTA/resource identities;
- active conservation source and labels/colors;
- cache-backed Generate succeeds without an upload repair;
- 20 rings in the final preview;
- no operation asks the user to supply an unavailable external prepared FASTA.

If tutorial prose currently implies the JSON is only a preview and cannot regenerate, correct it. Preserve provenance notes for bundled prepared inputs.

The expensive Generate truth remains owned by PR 1/2's dedicated WSSV browser spec. Capture metadata should assert the restored inputs, cache identity, and displayed result; it must not duplicate the 20-ring renderer test merely to validate a crop.

### 5.6 Tobacco chloroplast

Assert:

- session import succeeds before capture;
- no retired `cli_circular_track_*` field;
- `plastome_regions` and four annotations;
- current three-slot stack;
- current color/priority resources;
- upper-left legend and final preview semantics.

Do not retain a capture workaround that removes invalid fields or copies track slots at capture time.

## 6. Writing and structured-content audit

Audit operation text, captions, alt text, and repeated-field data together.

- Use user-facing labels such as `Run LOSAT`, `LOSATP`, `Collinear blocks`, and `Adjacent pairs`.
- Remove implementation phrases such as `active draft`, `canonical request`, or `SessionResourceFileView` from tutorials.
- Separate pre-generation setup, generated-result inspection, and optional post-generation editing.
- Use tables for repeated file mappings, ring settings, record rows, track slots, or color rules.
- Do not repeat the same values in a table, operation body, and caption.
- Captions describe the action/state shown; alt text describes the visible control/result.
- Do not add generic Requirements or Files-tab steps unless needed for the actual workflow.
- Do not describe automatic Generate output as a manual editor task.

Run a sibling scan when changing a repeated-field representation so the five target tutorials remain structurally consistent.

## 7. Decide whether each image changes

Classify referenced images before capture:

- `keep`: exact restored state and caption remain true;
- `recrop`: state is correct but crop/readability is weak;
- `replace`: visible selected value/result is stale;
- `add`: a single overloaded image cannot show the required action clearly;
- `remove`: duplicate operation/media adds no information.

Do not recapture a final preview solely because session JSON bytes changed. Compare PR 2's source/result semantics first.

Likely high-risk images:

- Vibrio comparison command/mode/adjacent-pair/advanced settings operations;
- Majanivirus comparison filters and final orthogroup preview;
- AT skew slot settings and final preview;
- WSSV conservation input/track/final preview;
- chloroplast track and final preview operations.

The actual audit result, not this likelihood list, determines the final media diff.

After tutorial references are correct, scan `gbdraw/web/gallery/media/` for unreferenced WebP files. In this PR, delete only files made unreferenced by these edits or stale files inside the five audited example directories; report unrelated candidates for later ownership review instead of broadening the diff. Before deletion, confirm that no tutorial, focused test, documentation page, or register entry still owns the file and update intentional owners in the same change. Report any referenced manual-only screenshot that still lacks deterministic capture metadata instead of implying that `--all` regenerated it.

## 8. Capture procedure

For every `recrop`, `replace`, or `add`:

1. Copy the old WebP to `/tmp`.
2. Add/fix declarative capture metadata and its test first.
3. Run `--check --strict` and prove the exact session state.
4. Capture with the operation filter, quality 94, DSF 3 unless the existing semantic crop uses a documented different density.
5. Keep temporary PNG/intermediate output in `/tmp`; commit WebP only.
6. Compare old/new side-by-side at the same Gallery rendered size and at natural size.
7. Reject a wider or taller crop that adds unrelated UI without improving operation identification.
8. Verify visible selected values, labels, and controls are inside the final crop.
9. Update caption/alt text in the same change if image truth changed.
10. Record final evidence in the operation register.

Example focused commands:

```bash
python tools/capture_gallery_tutorial_screenshots.py \
  --example vibrio-harveyi-group-collinear \
  --check --strict

python tools/capture_gallery_tutorial_screenshots.py \
  --example vibrio-harveyi-group-collinear \
  --operation "adjacent" \
  --quality 94 \
  --device-scale-factor 3
```

Repeat for each explicitly approved operation. Do not run an unreviewed broad rewrite of every media file.

## 9. Tests

Update `tests/test_gallery_capture_contracts.py` to assert:

- affected data-dependent operations use their own exact session;
- restored state is sufficient before correction actions;
- Vibrio no longer relies on `none -> adjacent` repair;
- required state paths and visible controls/text exist;
- result captures use exact-session semantic anchors;
- removed media references are absent;
- any new tables or capture fields satisfy current schema.

Update browser tests where exact image count/src or tutorial rendering changes.

For final preview captures, rely on PR 2's first-Generate parity as the session correctness owner. Capture tests verify composition/readability and must not duplicate the shared SVG comparator or full renderer semantics.

## 10. Verification commands

Focused checks per target:

```bash
python -m json.tool gbdraw/web/gallery/tutorials/HmmtDNA_ATskew.json >/tmp/HmmtDNA_ATskew.json.check
python -m json.tool gbdraw/web/gallery/tutorials/majanivirus_orthogroup.json >/tmp/majanivirus.json.check
python -m json.tool gbdraw/web/gallery/tutorials/vibrio-harveyi-group-collinear.json >/tmp/vibrio.json.check
python -m json.tool gbdraw/web/gallery/tutorials/WSSV_genome_comparison.json >/tmp/WSSV.json.check
python -m json.tool gbdraw/web/gallery/tutorials/tobacco-chloroplast.json >/tmp/chloroplast.json.check

python tools/capture_gallery_tutorial_screenshots.py --example HmmtDNA_ATskew --check --strict
python tools/capture_gallery_tutorial_screenshots.py --example majanivirus_orthogroup --check --strict
python tools/capture_gallery_tutorial_screenshots.py --example vibrio-harveyi-group-collinear --check --strict
python tools/capture_gallery_tutorial_screenshots.py --example WSSV_genome_comparison --check --strict
python tools/capture_gallery_tutorial_screenshots.py --example tobacco-chloroplast --check --strict
```

All-example and browser-presentation gates:

```bash
python tools/capture_gallery_tutorial_screenshots.py --all --check --strict
python -m pytest tests/test_gallery_capture_contracts.py -q
node --check gbdraw/web/gallery/gallery.js
npx playwright test tests/web/gallery-tutorial.playwright.spec.js --project=chromium
```

Also run PR 2's DOM-free all11 publication gate once as a read-only session guard. Do not unconditionally rerun the WSSV, Vibrio, and remaining full-render browser suites for a tutorial/media-only diff: their PR 2 result remains the renderer correctness evidence. Rerun a focused full-render owner only when this PR changes a production/session input unexpectedly or the capture audit exposes a reproducible replay defect; such a defect must be fixed through the implementation owner before capture resumes.

Check desktop and mobile Gallery layouts manually. Verify thumbnail/result composition, caption specificity, atomic tokens in tables, and lightbox natural-size behavior.

## 11. Files that must not change

- `gbdraw/web/gallery/sessions/`
- `gbdraw/web/gallery/sources/`
- `gbdraw/web/gallery/examples/`
- `gbdraw/web/gallery/thumbnails/`
- `gbdraw/web/gallery/examples.json`
- `tests/reference_outputs/`
- `examples/gbdraw_social_preview.png`

If a session or generated-artifact defect is discovered, stop the media capture, fix it through the PR 2 owner path in a follow-up implementation PR, rerun the publication/replay gates, then resume this PR. Do not compensate in screenshot scripts.

## 12. Acceptance criteria

- [ ] Operation register records all five target audits and final dispositions.
- [ ] No current prose says Vibrio session restores `none` when it now restores adjacent LOSAT.
- [ ] Capture actions do not silently repair broken session state.
- [ ] Every target data-dependent capture uses its own exact session and strong state assertions.
- [ ] Required visible controls/text are inside each crop.
- [ ] Only audited stale/unreadable images are changed.
- [ ] Every changed image has an old/new same-size review.
- [ ] Captions and alt text match the changed image.
- [ ] All11 strict capture checks pass.
- [ ] Desktop and mobile Gallery browser tests pass.
- [ ] PR 2's DOM-free all11 publication guard remains passing; any full-render suite selected by change scope also passes.
- [ ] No session/source/example/thumbnail/social-preview/reference artifact changes.

## 13. Diff review and handoff

Review separately:

```bash
git diff -- docs/internal/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md
git diff -- gbdraw/web/gallery/tutorials/
git diff -- tests/test_gallery_capture_contracts.py tests/web/gallery-tutorial.playwright.spec.js
git diff --stat -- gbdraw/web/gallery/media/
git status --short
```

For each changed WebP, include old/new dimensions, file size, operation, exact session, DSF/quality, and review disposition.

Suggested commit title:

```text
Align Gallery tutorials with reproducible session state
```

Suggested summary:

```text
- make tutorial capture metadata assert the repaired session state directly
- remove hidden Vibrio comparison-state repair from capture workflows
- recapture only affected Gallery operations and verify all tutorials strictly
```
