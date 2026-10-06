# Web typed boundaries, phase 2: the remaining modules

Status: proposed, 2026-10-06. The Owner approved phase 2 on 2026-10-06. It
continues [the phase 1 plan](WEB_TYPED_BOUNDARIES_IMPLEMENTATION_PLAN_2026-10-06.md),
whose decisions D1-D10, guard (Section 4), and common slice rules (Section 5)
still apply.
Baseline: `dev` at `8c5e3b9f` (merge of #855). The phase-2 error counts are
the same at `dc521b07` (merge of #858, phase 1 B3).
Scope: the 83 modules that phase 1 left in `UNCHECKED_MODULES`, one header
line of `tools/generate_mode_profiles.py`, and runtime fixes of the design
findings in Section 2.5. Workers stay out of scope (D8).

## 1. Objective

Check every module under `gbdraw/web/js/` except `workers/` with `tsc`. Phase 2
is complete when `UNCHECKED_MODULES` in `tests/web/typed-boundaries.test.mjs`
is `new Set([])`. The empty literal stays, so R14 keeps its registered
allowlist; the Gate reads the empty set as a contraction (checked in a scratch
commit at `8c5e3b9f`: Gate PASS, `contraction: R14 …#UNCHECKED_MODULES`).

The compiler options stay those of phase 1: `strict: false` plus
`strictFunctionTypes`, `strictBindCallApply`, `noImplicitThis`, and
`noFallthroughCasesInSwitch`. `strictNullChecks` and `noImplicitAny` are left
to a phase 3 (Section 10). As in phase 1, a slice changes shipped code only by
comments and JSDoc-cast parentheses. A design problem that typing exposes is
fixed in the slice only when it takes a few lines (D10). Otherwise it gets its
own runtime pull request (Section 5.2).

## 2. Evidence

### 2.1 Start state

| Set | Count | Modules |
| --- | ---: | --- |
| `UNCHECKED_MODULES` at `8c5e3b9f` | 100 | |
| `UNCHECKED_MODULES` at `dc521b07` | 98 | B3 (#858) checked `services/config.js` and `services/gallery-session-migration.js` |
| Phase 1 still open | 15 | A2: `app/feature-editor.js`, `app/feature-editor/*` (9), `app/svg-styles.js`, `app/rule-matching.js`. A4: `app/app-setup.js`, `app/run-analysis.js`, `app/watchers.js` |
| Phase 2 | 83 | Section 2.6 |

Phase 2 starts at 83 entries, after A4 merges (Section 6).

### 2.2 Error counts

The method is the one in phase 1 Section 2.3. A scratch copy of
`gbdraw/web/js` gets `// @ts-check` on the 83 modules. `tsc` then runs with
the options of `tests/web/types/tsconfig.json` and TypeScript 7.0.2 from the
shared `node_modules`. One program run takes 0.23 s.

| Configuration | 83 phase-2 modules | Whole tree, `checkJs: true` (166 modules) |
| --- | ---: | ---: |
| Phase 1 options | 225 in 30 modules (53 have none) | 364 (225, plus 139 in the 17 phase-1 modules open at `8c5e3b9f`) |
| + `strictNullChecks` | 685 | 1,794 (+1,430) |
| + `noImplicitAny` | 3,718 | 9,800 (+9,436; 6,154 are TS7006) |
| `strict: true` | 3,834 | 10,287 |
| + `noUnusedLocals` | +6 | +59 (+69 with `noUnusedParameters`) |

The 225 errors by code: TS2339 158, TS2353 32, TS2345 12, TS2554 9, TS2322 7,
TS2741 3, and one each of TS2739, TS2349, TS2769, and TS2810. They fall into
four kinds:

1. Inference artifacts of parameters with defaults (phase 1, kind 2). These
   are most of the errors.
   - A destructured binding without a default is dropped from the inferred
     options type. Examples: `services/losat.js` (23),
     `services/error-normalization.js` (26, from one inferred
     `{ code, stage }` shape), and `app/depth-track-state.js` (14).
   - A default function narrows the port to its own arity or return type.
     `drawnPlacement = () => null` makes `drawnPlacement(target)` a TS2554
     (`app/annotations/table-codec.js:173`). `isCurrent = () => true` makes a
     `boolean` result not assignable to `true` (`app/preview-runtime.js`).
     The same holds for `mutationAvailability = () => null`
     (`services/history.js:204`).
   - A declared receiver type removes them.
2. Object shapes built from `{}` or a literal and extended later. Examples:
   `createLinearSeq(overrides = {})` in `state.js` (15), the match descriptor
   in `app/pairwise-match-popup.js` (10), and `app/run-info.js`.
3. String unions narrowed to one literal by a default (`app/comparison-ui.js`,
   `app/linear-comparisons.js`).
4. `Error` values carrying extra properties (`error.code`, `error.canceled`).

None of the 225 is a port-drift defect: every TS2554 site was checked against
its default and its callers. The unread arguments that do exist are reported
by `noUnusedLocals` and `noUnusedParameters` instead (Section 2.5).

Grouped by the type the message names, the 225 errors form 77 clusters, each
removed by one declared type.

### 2.3 Size estimates

The estimate for a slice adds these lines:

- 1 per module (the pragma);
- 4 per exported factory that takes parameters;
- 1 per factory parameter, counting each destructured name;
- 6 per error cluster;
- 1 per error.

Phase 1 slices landed at 0.5 to 3.4 times their own estimates (net: A1 +185,
A3 +286, B1 +407, B2 +104, B3 +126, C1 +34). The slices in Section 2.6
therefore keep the estimate at or under about 120 lines. A 2.5-times overrun
still stays near the cap of about 300 net lines.

### 2.4 Special modules

| Module | Finding | Consequence |
| --- | --- | --- |
| `mode-profiles.generated.js` | `tools/generate_mode_profiles.py:18-30` writes it. `tests/test_web_mode_profiles.py` runs the generator with `--check`. The guard needs `// @ts-check` on line 1, so the output cannot stay unchanged. The module has 0 errors. | P2 adds the pragma to the generator's header string and regenerates; the file changes by that line only (D12). A `tools/` path makes the PR `full` in ci-impact (`FULL_BY_DEFAULT`). |
| `services/standalone-interactivity-assets.js` | The two template literals are embedded in exported SVG and in the Gallery. Python reads them by name (`gbdraw/render/interactive_svg.py:1238-1250`). The browser wheel ships the file (`gbdraw/_build_support.py:28,101,220`). 0 errors. | Pragma only, outside the literals (D15). |
| `app/python-helpers.js` | One template literal of Python source, imported only by the render worker. Two tests split the file on its first backtick (`tests/test_protein_colinearity.py:177`, `tests/test_collinearity.py:2321`); one splits on the export line (`tests/test_record_metadata.py:124`). 0 errors. | Pragma only. No comment containing a backtick above the literal (D15). |
| `config.js` | Rewritten by regular expressions anchored on its export lines (`gbdraw/_build_support.py:20-21`, `gbdraw/_web_assets.py:377-379`). 0 errors. | Pragma only. |
| `state.js` | Vue comes from `window.Vue`, which `tests/web/types/web-globals.d.ts` declares `any`, so every `ref` and `computed` value is `any`. | D13. |

CI classes (`classifyPath` in `tools/ci-impact-policy.mjs`): every phase-2
module is `web-runtime`, except three `session-persistence` modules
(`app/history-inputs.js`, `app/history-shortcuts.js`,
`app/session-feature-metadata.js`) and the five LOSAT modules
(`losat-integration`). The `web-runtime` tier runs no pytest job, so the
Python readers of the modules above run locally (Section 5.1).

### 2.5 Design findings

The Owner asked on 2026-10-06 that coupling and separation-of-concerns
problems found while typing be fixed. The measurement already shows these.
Slice agents add more from F-08 on (Section 5.3).

| ID | Where | Problem | Plan |
| --- | --- | --- | --- |
| F-01 | `state.js` and 14 modules under `services/` | 79 runtime import edges point up the R13 layers ("A module depends only on lower layers") into 39 modules under `app/`. `services/config.js` imports 23 of them, `services/session-request.js` 21, `services/session-active-config-contract.js` 7, and `state.js` 7 (`color-utils`, `feature-selector`, `feature-visibility`, `layout-preferences`, `linear-comparisons`, `match-sequences`, `plot-title-position`). No detector and no Gate rule checks runtime import direction. R14 assertion 5 forbids only the type form. | Not fixed in phase 2. A lower layer that needs a type defined under `app/` declares the narrow shape it reads (D14). Fixing it moves pure helpers out of `app/`: new paths, rows of the ownership table in `gbdraw/web/CLAUDE.md`, and a new detector with a baseline. That is a separate plan after phase 2 (OD-1). |
| F-02 | `app/annotations/record-catalog.js:162-167`, `app/app-setup.js:1125-1135,2457,3086`, `app/run-analysis.js:2057,2172-2177` | `buildAnnotationRecordCatalog` destructures `inputType` and `loadComparison` and reads neither. `getAnnotationRecordCatalog(loadComparisonOverride, …)` in the composition root still computes both. Two ports carry the dead value from `app/run-analysis.js`: `validateAnnotationTargets({ loadComparison })` and the first parameter of `prepareLinearRecordCatalog(loadComparison, …)`. | Runtime PR F-02 before P5: remove the two options, the dead parameters of both ports and of `getAnnotationRecordCatalog`, and the values computed only for them. No behavior change. |
| F-03 | `app/pairwise-match-popup.js:358,746-784` | Six module-private declarations are never referenced: `featureOrthogroupId`, `getRenderedFeatureForMember`, `getFeatureForMember`, `getGroupMemberForFeatureSvgId`, `getOrthogroupForMatch`, and `buildFallbackOrthogroup`. That is about 45 lines, more than a D10 fix. | Runtime PR F-03 before P7, together with F-04 and F-06. |
| F-04 | `app/orthogroups.js:184` | `memberFastaText(member, sequenceKind, orthogroupId)` never reads `orthogroupId`. Both callers (`:628`, `:636`) pass it. | In PR F-03. |
| F-05 | `services/history.js` `createHistoryManager`, `services/history-snapshot.js` `createHistorySnapshotService` | Each takes `fileStore` whole. These are the two `whole-object-port` subjects of the owner-graph baseline that Phase E kept by design. | P3 types `fileStore` as the narrow set of functions each receiver calls. If a receiver calls one function, P3 reports a port split as a runtime PR. Otherwise it records why the object stays. |
| F-06 | Checked phase-1 modules | `noUnusedLocals` and `noUnusedParameters` also report dead code there: 3 unused declarations and 1 unused inner function in `app/circular-track-slots.js` (`:66,497,547,2339`); 2 unused functions in `services/standalone-interactivity.js` (`:239,261`); the unread `filesData` parameter and the unused `textToBase64` import in `services/session-request.js` (`:2218,147`); the unread `request` parameter of `applyPlan` in `app/similarity-alignment.js:1026`; and the unused `newY` in `app/legend/entry-actions.js:327`. | In PR F-03. `newY` waits for OV-47, which edits `app/legend/entry-actions.js`. |
| F-07 | `state.js:466`, `services/history-snapshot.js:982,1111`, `services/config.js:3495,3556`, `services/reset.js:104`, `app/feature-editor/label-actions.js:445-511`, `app/app-setup.js:391,5116`, `index.html:6879` | `labelOverrideBuildWarning` is only ever cleared, restored, or captured, and never set to a non-empty value, so its template warning never shows (observed by Phase E on 2026-10-06). `services/config.js:3495` is `captureSessionImportTransientState`, which is not persisted. | Runtime PR F-07 before P2 and P3: remove the state and its readers. No behavior change (OD-4). |

Line numbers are those of `dev` at `8c5e3b9f`.

### 2.6 Slices

The order of the rows is the merge order of Section 5.4, not the numbering.
Errors are shown as phase 1 options / `strict: true`.

| Slice | Modules | Lines | Factories (params) | Errors | Clusters | Estimate | CI class | Waits for |
| --- | ---: | ---: | --- | ---: | ---: | ---: | --- | --- |
| P4 Errors, downloads, utilities | 10 | 1,430 | 2 (1) | 27 / 163 | 5 | ~70 | web-runtime | A4 |
| P8 Comparisons | 3 | 1,512 | 2 (9) | 24 / 222 | 11 | ~110 | web-runtime | A4 |
| P2 Shared state and profiles | 5 | 1,619 | 8 (4) | 15 / 92 | 1 | ~60 | full (generator) | F-07 |
| P3 History | 4 | 2,849 | 4 (38) | 6 / 477 | 3 | ~80 | session-persistence | F-07 |
| P5 Annotations | 7 | 1,046 | 4 (6) | 3 / 186 | 3 | ~50 | web-runtime | F-02 |
| P7 Match popups | 5 | 3,745 | 2 (2) | 35 / 614 | 12 | ~120 | web-runtime | F-03 |
| P6 Feature selection and visibility | 8 | 3,342 | 3 (5) | 23 / 501 | 12 | ~120 | session-persistence | #857 |
| P9 LOSAT | 5 | 1,852 | 2 (2) | 32 / 236 | 8 | ~100 | losat-integration | A4 |
| P1 Generate path and Result mount | 7 | 2,505 | 5 (7) | 22 / 437 | 8 | ~100 | web-runtime | OV-46, OV-47 |
| P10 Linear layout and records | 10 | 1,591 | 2 (8) | 3 / 246 | 3 | ~50 | web-runtime | A4 |
| P11 Track slots and option state | 11 | 2,902 | 5 (2) | 19 / 317 | 7 | ~80 | web-runtime | A4 |
| P12 Shell, run info, assets, Vue entry | 8 | 11,346 | 4 (11) | 16 / 343 | 4 | ~80 | web-runtime | all other slices |
| Total | 83 | 35,739 | 43 (95) | 225 / 3,834 | 77 | ~1,020 | | |

Module lists:

- P1: `app/candidate-render.js`, `app/feature-dom.js`, `app/preview-runtime.js`,
  `app/results.js`, `services/svg-result-ingestion.js`,
  `services/svg-result-normalization.js`, `services/svg-serialization.js`.
- P2: `state.js`, `config.js`, `mode-profiles.js`,
  `mode-profiles.generated.js`, `web-ux-profile.js`.
- P3: `services/history.js`, `services/history-snapshot.js`,
  `app/history-inputs.js`, `app/history-shortcuts.js`.
- P4: `services/error-normalization.js`, `services/export.js`,
  `services/pdf-fonts.js`, `services/text-download.js`, `utils/clipboard.js`,
  `utils/feature-rendering.js`, `utils/optional-positive-number.js`,
  `utils/png.js`, `utils/tsv-cell.js`, `utils/zip.js`.
- P5: `app/annotations.js`, `app/annotations/*` (6: `record-catalog`,
  `record-selector`, `state`, `table-codec`, `target-actions`, `validation`).
- P6: `app/feature-selection.js`, `app/feature-selector.js`,
  `app/feature-utils.js`, `app/feature-visibility.js`,
  `app/feature-metadata-extraction.js`, `app/session-feature-metadata.js`,
  `app/specific-color-rules.js`, `app/color-utils.js`.
- P7: `app/pairwise-match-popup.js`, `app/match-sequences.js`,
  `app/orthogroups.js`, `app/feature-sequence-fasta.js`,
  `app/record-source-coordinates.js`.
- P8: `app/comparison-ui.js`, `app/linear-comparisons.js`,
  `app/conservation-series.js`.
- P9: `app/losat-cache.js`, `app/losat-normalization.js`,
  `app/losat-settings.js`, `services/losat.js`, `services/losat-thread-plan.js`.
- P10: `app/linear-label-visibility.js`, `app/linear-record-layout.js`,
  `app/linear-record-selector.js`, `app/linear-sources.js`,
  `app/linear-typography.js`, `app/record-discovery.js`, `app/record-groups.js`,
  `app/record-options.js`, `app/genbank-header.js`, `app/file-imports.js`.
- P11: `app/depth-track-state.js`, `app/depth-tracks.js`,
  `app/track-slot-colors.js`, `app/track-slot-display.js`,
  `app/track-slot-validation.js`, `app/auto-value-display.js`,
  `app/current-option-values.js`, `app/definition-line-style-state.js`,
  `app/layout-preferences.js`, `app/plot-title-position.js`, `app/palettes.js`.
- P12: `app/ui.js`, `app/right-drawer.js`, `app/run-info.js`,
  `app/python-helpers.js`, `services/pyodide-assets.js`,
  `services/standalone-interactivity-assets.js`, `components.js`, `app.js`.

Compared with the starting proposal, LOSAT is split from the comparisons (P9),
because it is a separate CI class and the two together would estimate about
210 lines. The match popups are split from feature selection (P7), and
`app/feature-dom.js` and `app/preview-runtime.js` join the Generate path
(P1), whose ingestion and mount they serve.

## 3. Decisions

D1-D10 of phase 1 apply unchanged. All of D11-D16 are Owner-delegated: each
recommended option is adopted under the standing instruction of 2026-09-29.
Section 9 lists the ones the Owner may want to override.

- **D11 Order.** Slices start after A4 merges, which ends phase 1. Within
  that, lower layers go first: `services/error-normalization.js` is imported
  by 32 modules, and its declared types remove artifacts in later slices.
  Slices whose `UNCHECKED_MODULES` lines are adjacent merge one after another
  (Section 5.4). A runtime fix that a slice depends on merges before the
  slice.
- **D12 The generator emits the pragma.** `tools/generate_mode_profiles.py`
  writes `// @ts-check` above its "Generated by" line. The alternative, a guard
  exemption for generated files, is an authority change that loosens R14.
- **D13 No whole-state typedef.** `state.js` declares the shapes it builds
  (`LinearSeq`, the default rule and draft objects) and the parameters of its
  exported factories. It does not type the reactive state object: its values
  are `any` because Vue is `any`, so such a typedef would not be enforced.
- **D14 Narrow local shapes in lower layers.** When `state.js`, a module
  under `services/`, or one under `utils/` needs a type defined under `app/`,
  it declares the narrow shape it reads in its own module. It does not import
  the type (R14 assertion 5) and does not move the module in a slice (F-01).
- **D15 Literal-carrying modules get only the pragma.** Nothing changes
  inside a template literal. In `app/python-helpers.js`, no comment containing
  a backtick goes above the literal. The slice that checks these modules runs
  the Python readers locally (Section 5.1).
- **D16 Provider-side casts.** A slice may add JSDoc to a module that is
  already checked when the new receiver types expose an artifact at a provider
  (most often `app/app-setup.js` after A4). The change stays comment-only, and
  the body lists those files. A real defect at a provider is a D10 fix when it
  takes a few lines; otherwise it is a runtime PR that merges first.

## 4. Guard

No change. The guard, `tests/web/types/tsconfig.json`, and R14 stay as T0
landed them. The last slice (P12) leaves `const UNCHECKED_MODULES = new
Set([]);`.

## 5. PR sequence

### 5.1 Slice rules

Phase 1 Section 5 and `tb-common.md` apply: one STANDARD pull request per
slice, a Sonnet agent per slice, and "This is not architecture-bearing" in the
body. The additions are these:

- Before opening a slice, merge every open slice branch into a scratch
  worktree based on `origin/dev` and run
  `node --test tests/web/typed-boundaries.test.mjs` there. JSDoc in one slice
  can change inference in another slice's modules.
- P2 also runs `python tools/generate_mode_profiles.py --check` and
  `pytest tests/test_web_mode_profiles.py tests/test_web_packaging.py`.
- P12 also compares the Python view of the literals before and after:
  `_load_standalone_assets()` from `gbdraw/render/interactive_svg.py`, and the
  text between the first and last backtick of `app/python-helpers.js`, must be
  byte-identical (SHA-256). It then runs
  `pytest tests/test_record_metadata.py tests/test_web_packaging.py tests/test_interactive_svg_cli_format.py -m "not slow"`
  and the helper tests of `tests/test_collinearity.py` and
  `tests/test_protein_colinearity.py`.
- The comment-only check of phase 1 Section 5 must report only the D10 fixes
  that the body names. P2's generator line is outside `gbdraw/web/js`, and the
  body names it.

### 5.2 Runtime PRs for the findings

| PR | Content | Files | Waits for | Before |
| --- | --- | --- | --- | --- |
| F-02 | Remove the unread `inputType` and `loadComparison` options of `buildAnnotationRecordCatalog`, and the dead `loadComparison` parameters of `getAnnotationRecordCatalog`, `validateAnnotationTargets`, and `prepareLinearRecordCatalog` | `app/annotations/record-catalog.js`, `app/app-setup.js`, `app/run-analysis.js` | A4, OV-47 (both edit `app/app-setup.js`) | P5 |
| F-03 | Remove the declarations and parameters of F-03, F-04, and F-06 that are never read | `app/pairwise-match-popup.js`, `app/orthogroups.js`, `app/circular-track-slots.js`, `services/standalone-interactivity.js`, `services/session-request.js`, `app/similarity-alignment.js`, `app/legend/entry-actions.js` | OV-47 for `app/legend/entry-actions.js` (or drop that line) | P7 |
| F-07 | Remove `labelOverrideBuildWarning` | the files of F-07 in Section 2.5 | A4 | P2, P3 |

Each is STANDARD and follows Phase E practice:

- R13 applies.
- The owner-graph baseline only shrinks; compare
  `node tools/report-web-owner-graph.mjs` before and after.
- Behavior does not change. If it would, the PR is a bug fix instead (OV-48
  onward, with Before and After screenshots).

`services/standalone-interactivity.js` and `services/session-request.js` put
F-03 in the `session-persistence` CI class.

### 5.3 Findings during slices

Slice agents record each design problem they find as F-08 onward in
`/home/kawato/gbdraw-baselines/typed-boundaries-phase2-20261006/HANDOFF.md`,
with location, problem, and fix or the reason it is not fixed:

- A fix of a few lines in a module the slice checks goes into the slice as a
  named D10 fix.
- A larger fix becomes a runtime PR like those in Section 5.2. The slice
  waits for it only when its types depend on the fix.
- A fix that adds a baseline or allowlist entry, or relaxes rule text, waits
  for the Owner.
- A real defect becomes OV-48 onward. Minor ones are fixed on the spot.

### 5.4 Order

`UNCHECKED_MODULES` is sorted, so two slices whose entries are adjacent in
the list conflict. These waves contain no such pair. Slices in one wave can
be open at the same time; a wave can start before the previous one has fully
merged when its own conditions hold.

| Wave | Slices and runtime PRs | Conditions |
| --- | --- | --- |
| 0 | phase 1 A2, A4 (other session) | — |
| 1 | P4, P8; F-02, F-03, F-07 | A4 merged; F-02 and F-03 after OV-47 |
| 2 | P2, P3, P5, P7 | F-07 before P2 and P3; F-02 before P5; F-03 before P7 |
| 3 | P6, P9 | #857 merged (P6 types `app/feature-visibility.js`) |
| 4 | P1, P10 | OV-46 and OV-47 merged (P1 types `app/candidate-render.js` and `services/svg-result-ingestion.js`) |
| 5 | P11 | — |
| 6 | P12 | all other slices merged; it empties the list |

After a merge that shrinks `UNCHECKED_MODULES`, an open PR whose "Web base
policy (trusted base)" check has not passed yet gets `gh pr update-branch`
and auto-merge again (phase 1 lesson).

## 6. Interaction with in-flight work

| In-flight | Overlap with phase 2 | What waits |
| --- | --- | --- |
| A2 (feature editor, 12 modules) and A4 (`app/app-setup.js`, `app/run-analysis.js`, `app/watchers.js`), phase 1 | A4 checks the provider side of most phase-2 factories (`createHistoryManager`, `createPreviewRuntime`, `createFeatureSelection`, and others) | Every slice PR and F-02/F-07 wait for A4 to merge. Planning and unopened investigation may run before. |
| #857 (Legend rerender, OV-42/43/44) | `app/feature-visibility.js` (P6) | P6 |
| OV-46 (`fix/web-legend-only-row-any-result-ov46`) | `app/candidate-render.js`, `services/svg-result-ingestion.js` (P1) | P1 |
| OV-47 (uncommitted in `.worktrees/ov47`) | `app/candidate-render.js` (P1), `app/app-setup.js` (F-02, F-07), `app/legend/entry-actions.js` (F-03) | P1, F-02, the F-03 line in `app/legend/entry-actions.js` |
| Owner-graph v2 governance PR (baseline to v2; follows #856, waits for Owner approval) | no Web module | Runtime PRs of Section 5.2 re-read the owner-graph baseline after it merges |

Check overlaps again before opening each PR: `gh pr list` and
`gh pr diff --name-only <n>`.

## 7. Acceptance

- `UNCHECKED_MODULES` is `new Set([])`, and
  `node --test tests/web/typed-boundaries.test.mjs` passes on `dev`.
- Every exported factory with parameters declares them (guard assertion 3);
  every phase-2 port is a function type declared by its receiver.
- Every slice PR passes the comment-only check, apart from the D10 fixes and
  the P2 generator line that its body names.
- The literals of `services/standalone-interactivity-assets.js` and
  `app/python-helpers.js` are byte-identical to `8c5e3b9f`.
- F-02, F-03 (with F-04 and F-06), and F-07 are merged. The owner-graph
  baseline has not grown.
- The final report lists the merged PRs, the count of `UNCHECKED_MODULES`
  from 83 (or the count at the start) to 0, the D10 fixes, OV-48 onward, the
  F-xx findings, the Owner-delegated choices, and a recommendation on phase 3.

## 8. Rollback

As in phase 1 Section 8. A slice is fixed forward, because reverting it would
expand a registered literal together with runtime paths. A runtime PR of
Section 5.2 reverts as an ordinary STANDARD PR.

## 9. Owner decisions

None blocks. The Owner may want to override these delegated choices:

1. **OD-1 F-01, the layer inversion.**
   - (A, adopted) Record it now and write a separate plan after phase 2: a
     checker-only PR adds an import-direction detector, an authority PR
     registers its baseline, and runtime PRs move pure helpers out of `app/`.
     Each moved module is checked in its move PR (D3).
   - (B) Move modules during phase 2. Slices would stop being comment-only,
     and the ownership table in `gbdraw/web/CLAUDE.md` is an authority file.
   - (C) Accept the inversion and narrow R13's text. That relaxes a rule,
     which only the Owner can do.
2. **OD-2 D12, the generator header.** The alternative is a guard exemption
   for generated modules (authority change).
3. **OD-3 D13, no whole-state typedef.** The alternative types the reactive
   state object at about one line per field (several hundred lines), which
   `tsc` would not enforce while Vue is `any`.
4. **OD-4 F-07, removing `labelOverrideBuildWarning`.** The alternative keeps
   it and types it. That keeps a template branch that cannot show.

## 10. Non-goals and phase 3 input

- No change to the compiler options, the guard, or R14 text.
- No `.ts` source, build step, or `.d.ts` under `gbdraw/web/`.
- Workers and the contract between the Vue template and the setup return stay
  unchecked.
- No type for Python-owned option fields (R7).
- Phase 3 input, measured at `8c5e3b9f` on the whole tree:
  - `strictNullChecks` +1,430: TS2339 687, TS2322 370, TS18047 304; the top
    modules are `app/run-analysis.js` 279, `app/app-setup.js` 236, and
    `app/preview-runtime.js` 136.
  - `noImplicitAny` +9,436: TS7006 6,154, TS7031 903, TS7005 731.
  - `noUnusedLocals` +59 (+69 with `noUnusedParameters`). After F-03 removes
    the dead declarations, this would cost little. The report at the end of
    phase 2 gives a recommendation and a cost estimate.
