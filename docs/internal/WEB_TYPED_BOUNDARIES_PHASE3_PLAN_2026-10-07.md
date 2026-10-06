# Web typed boundaries, phase 3: unused declarations and strictNullChecks

Status: proposed, 2026-10-07. It continues
[the phase 2 plan](WEB_TYPED_BOUNDARIES_PHASE2_PLAN_2026-10-06.md), whose
Section 10 gave the input, and
[the phase 1 plan](WEB_TYPED_BOUNDARIES_IMPLEMENTATION_PLAN_2026-10-06.md).
Their decisions D1-D16 still apply unless a decision below replaces one.
Baseline: `dev` at `b58176bd` (merge of #893). Phase 2 ended at `207ba99e`
(#890): every module under `gbdraw/web/js/` except `workers/` starts with
`// @ts-check`, and `UNCHECKED_MODULES` is `new Set([])`.
Scope: the compiler options of the R14 guard (`tests/web/types/tsconfig.json`),
the runtime and JSDoc changes that those options need, and a temporary guard
for `strictNullChecks`. Workers stay out of scope (D8).

## 1. Objective

The guard runs `tsc` with `strict: false` plus `strictFunctionTypes`,
`strictBindCallApply`, `noImplicitThis`, and `noFallthroughCasesInSwitch`.
It does not check unused declarations or `null`. Phase 3 has two goals:

1. `noUnusedLocals` and `noUnusedParameters` are on in the guard config, with
   0 diagnostics.
2. `strictNullChecks` is on in the guard config, with 0 diagnostics. If the
   planned pull requests do not reach 0, a per-file count ratchet holds the
   remainder and phase 4 continues from it.

`noImplicitAny` is not a goal. Phase 3 measures it again at the end and
recommends whether phase 4 should take it on.

## 2. Evidence

### 2.1 Start state

Measured by running `tsc` (TypeScript 7.0.2, about 0.5 s per run) with the
guard config plus one stricter option at a time, on a copy of `gbdraw/web/js`
at `b58176bd`. The numbers are the same at `207ba99e`.

| Configuration | Diagnostics | Files | Notes |
| --- | ---: | ---: | --- |
| Guard config | 0 | 0 | guard green |
| + `noUnusedLocals` | 40 | 6 | all TS6133: `app/app-setup.js` 21, `app/run-analysis.js` 14, `app/feature-editor/color-actions.js` 2, `app/feature-editor/svg-actions.js` 1, `app/feature-editor/visibility-actions.js` 1, `app/watchers.js` 1 |
| + `noUnusedParameters` | 4 | 3 | `app/app-setup.js` 2, `app/feature-editor/rule-actions.js` 1, `app/match-sequences.js` 1 |
| + `strictNullChecks` | 1,114 | 88 | app 729, services 364, root 14, utils 7; 49 files have 5 or fewer |
| + `noImplicitAny` | 8,811 | 158 | mostly TS7006 (untyped parameters) |
| `strict: true` | 9,219 | 158 | |

### 2.2 Unused declarations

Of the 44, 33 are names destructured from `state` and never read
(`app/app-setup.js` 19, `app/run-analysis.js` 12, one each in
`app/feature-editor/color-actions.js` and `app/watchers.js`). Destructuring a
reactive object has no side effect outside an effect, so these bindings can go.
The others are three import specifiers (`normalizeCollinearSearchScope` in
`app/run-analysis.js`, `filterFeatureFillTargets` in
`app/feature-editor/svg-actions.js`, `resolveDisplayProteinId` in
`app/feature-editor/visibility-actions.js`), one local function
(`getLiveLegendColor` in `app/feature-editor/color-actions.js`), and four
parameters: `description` of `copyRunInfoCommand` and `options` of
`runAnalysis` (both in `app/app-setup.js`), `currentFeat` of
`findExistingColorForCaption` (`app/feature-editor/rule-actions.js`), and the
index of a `forEach` callback (`app/match-sequences.js`).

A name that the setup function returns to the Vue template counts as read,
because `tsc` sees the return object. The template cannot read a local that
the setup function does not return, so a TS6133 local is not read by the
template either.

### 2.3 strictNullChecks diagnostics

By code: TS2339 305, TS2322 284, TS18047 235, TS2345 150, TS2349 27,
TS18048 27, TS2538 19, TS2531 17, TS2532 14, and 36 in 12 other codes.

Most come from four inference patterns in JavaScript files:

1. `let x = null;` gets the type `null`. A later `x = value` is TS2322 (for
   example "Type 'number' is not assignable to type 'null'", 23 times). A
   property read after a truthiness check is TS2339 on type `never` (299 of
   the 305 TS2339).
2. `const list = [];` gets `never[]` (TS2322 and TS2345 with `never[]`).
3. A parameter whose default is `null` gets the type `null`. A call of it is
   TS2349 "This expression is not callable".
4. Values that can really be absent: a `querySelector` result, a `Map.get`
   result, an optional field. These give TS18047, TS18048, TS2531, and TS2532.

Patterns 1-3 are fixed by declaring the true type at the declaration, with no
executable change. Some TS18047 errors cluster on one declaration:
`losatTiming` alone has 64, all in `app/run-analysis.js`. Pattern 4 includes
both guarded cases and real null paths (Section 3, D18).

Largest files: `app/run-analysis.js` 195, `app/feature-editor/svg-actions.js`
102, `services/config.js` 85, `services/history-snapshot.js` 76,
`app/app-setup.js` 61, `services/session-request.js` 51,
`services/diagram-generation.js` 40, `app/similarity-alignment.js` 35.

### 2.4 Paths that move during phase 3

The layering plan of the owner-coupling session (Phase E, R13 import
direction) moves about 40 modules from `app/` into `services/` and `utils/`
in four runtime pull requests, A to D, starting on 2026-10-07:

| PR | Moves and splits with strictNullChecks diagnostics on `b58176bd` |
| --- | --- |
| A | `services/error-normalization.js` 4 -> `utils/`; `app/conservation-series.js` 1, `app/genbank-header.js` 1, `app/definition-line-style-state.js` 3, `app/linear-record-layout.js` 2, `app/specific-color-rules.js` 10, `app/depth-track-state.js` 9, `app/track-slot-validation.js` 10 -> `services/`; `app/circular-track-slots/measure-editor.js` 3 -> `services/circular-track-measure.js`. It also edits `app/legend-layout/composition-actions.js` 18, `app/legend/stroke-actions.js` 1, and `services/config.js` 85 (a port, no move) |
| B | `app/linear-comparisons.js` 19, `app/match-sequences.js` 11, `app/feature-visibility.js` 2 -> `services/`; slices of `app/right-drawer.js` 1, `app/run-info.js` 10, `app/linear-typography.js`, and `app/rule-matching.js` |
| C | slices of `app/circular-track-slots.js` 16 and `app/linear-track-slots.js` 4 -> two new `services/` modules |
| D | `app/record-discovery.js` 3, `app/feature-metadata-extraction.js` 9 -> `services/`; `app/session-feature-metadata.js` 10 -> `services/session-feature-recovery.js`; slice of `app/record-display-options.js` 1; `app/legend/utils.js` 1 -> `services/legend-svg.js`; `app/feature-search/preview-svg.js` 1 merged into `services/svg-serialization.js` |

The Gate classifies a registered map of counts (`kind: count-map`) key by key.
A key that the trusted base lacks is an expansion
(`tools/check-web-change-budget.mjs`, `classifyAllowlistDelta`), and an
expansion beside production runtime paths fails the Gate with
`design-rule.co-change` (WEB_CHANGE_POLICY.md "Design-rule co-change"). With a
per-file baseline in place, a move of a file with a nonzero count would
therefore fail: its old key no longer matches a file, and its new key is an
expansion. A split that carries diagnostics into a new module fails the same
way. This decides the order (D21).

The same session also has OV-65 (#896) and OV-80 to OV-82 in flight. They edit
`services/config.js`, `state.js`, `app/app-setup.js`, and the Legend modules.

### 2.5 Slices

One slice is one commit and one agent. Slices of the same kind share a pull
request (D20).

| PR | Slice | Files | Diagnostics | Content |
| --- | --- | ---: | ---: | --- |
| U | U | 7 | 44 (unused) | Section 2.2 |
| S1 | N1 | 19 | 52 | `utils/` (2), `mode-profiles.js`, `state.js`, and the `services/` files with 5 or fewer |
| S1 | N2 | 7 | 81 | `services/history.js` 16, `services/svg-result-ingestion.js` 16, `services/losat.js` 13, `services/bounded-json-transport.js` 11, `services/session-import-client.js` 10, `services/legacy-similarity-alignment.js` 8, `services/feature-edit-migration.js` 7 |
| S2 | N3 | 1 | 85 | `services/config.js` |
| S2 | N4 | 1 | 76 | `services/history-snapshot.js` |
| S2 | N5 | 2 | 91 | `services/session-request.js` 51, `services/diagram-generation.js` 40 |
| S3 | N6 | 15 | 117 | Feature editor and selection: `app/feature-editor.js`, `app/feature-editor/*` except `svg-actions.js`, `app/feature-selection.js` 25, `app/feature-visibility.js`, `app/feature-metadata-extraction.js`, `app/feature-search/*`, `app/svg-styles.js`, `app/specific-color-rules.js`, `app/session-feature-metadata.js` |
| S3 | N7 | 8 | 60 | Legend and Legend layout: `app/legend-layout.js`, `app/legend-layout/*`, `app/legend/*` |
| S3 | N8 | 26 | 114 | Tracks, records, comparisons, and annotations: `app/circular-track-slots.js` 16, `app/linear-comparisons.js` 19, `app/match-sequences.js` 11, `app/track-slot-validation.js` 10, `app/depth-track-state.js` 9, and 21 files with 4 or fewer |
| S3 | N9 | 6 | 80 | Generate path and shell: `app/similarity-alignment.js` 35, `app/preview-runtime.js` 13, `app/candidate-render.js` 10, `app/run-info.js` 10, `app/ui.js` 10, `app/watchers.js` 2 |
| S4 | N10 | 1 | 195 | `app/run-analysis.js` |
| S4 | N11 | 1 | 102 | `app/feature-editor/svg-actions.js` |
| S4 | N12 | 1 | 61 | `app/app-setup.js` |
| | Total | 88 | 1,114 | |

Declarations shared across slices change inference in other slices. A return
type declared in a `services/` module can expose a missing null check in an
`app/` caller, so lower layers go first (D20), and each slice is measured
against the whole tree.

## 3. Decisions

D1-D16 apply. All of D17-D23 are Owner-delegated: each recommended option is
adopted under the standing instruction of 2026-09-29. Section 9 lists the ones
the Owner may want to override.

- **D17 Baseline shape.** The ratchet is `STRICT_NULL_BASELINE` in a new guard
  test, `tests/web/strict-null-ratchet.test.mjs`. It is an object literal of
  `{ '<path>': <count> }`, sorted by path. Paths are relative to
  `gbdraw/web/js/`, like `UNCHECKED_MODULES`. A count is the number of `tsc`
  diagnostics in that file under `tests/web/types/tsconfig.strict-null.json`,
  which extends the guard config and adds only `strictNullChecks: true`.
  The test compares both ways, like `TYPE_DEBT_BASELINE`
  (`tests/test_type_check_ratchet.py`) and `LAYER_IMPORT_BASELINE`:
  - a count above its entry fails and lists the diagnostics;
  - a count below its entry fails until the same pull request lowers or
    removes the entry;
  - an entry for a missing file fails;
  - a file without an entry, including a new file, has no allowance.

  The baseline is registered in R14 in `tools/web-design-rule-guards.json`
  with `kind: count-map`. Lowering an entry or removing one is then a
  contraction that a runtime pull request may carry. Adding or raising an
  entry is an expansion that needs an authority-only pull request.
  The test also requires the strict config to be exactly the two lines above,
  so that the ratchet cannot be loosened without a guard change.
- **D18 Classification of fixes.** Every change in an S pull request is one of
  three kinds, and the body says which:
  - (a) Types only: comments and JSDoc-cast parentheses. It declares a true
    type (`/** @type {T | null} */` on `let x = null`, `/** @type {T[]} */` on
    `[]`, `@param {T | null} [name]` on a `null` default), or types a value
    that the code already guarantees: after an existing check, an element of
    the module's own template, or a typedef field that every writer sets.
  - (b) D10: a change of a few lines that states existing null handling in
    code without changing behavior, for example an early return that already
    happened implicitly. The body lists each one.
  - (c) A real null path: a reachable path where `null` makes the code throw
    or misbehave. It gets an OV number (OV-70 to OV-79, then OV-90 onward;
    Phase E uses OV-80 to OV-89). A fix of a few lines lands in the slice with
    a test. A larger fix is a separate runtime pull request, and the slice
    keeps that diagnostic.
- **D19 Cast rule.** A cast that only removes `null`
  (`/** @type {T} */ (y)`) needs a one-line comment above it that says why `y`
  is not null there. Without that line the site is (b) or (c). Casts to `any`,
  `*`, or `Object`, declaring a value `any`, and weakening another module's
  declared type to hide a diagnostic are not used. R14 already forbids
  `@ts-ignore`, `@ts-expect-error`, and `@ts-nocheck`.
- **D20 Slices.** Slices follow the layers, lowest first: `utils/`, root, and
  the small `services/` files (S1), then the large `services/` files (S2), then
  `app/` files with fewer diagnostics (S3), then the three largest `app/` files
  (S4). Each pull request carries one kind of slice, one commit per slice,
  four pull requests in all. This follows the CI load rules (Section 5.1),
  as phase 2's combined #890 did. A slice ends with every one of its files at 0
  and no other file's count higher. When its types expose a diagnostic in
  another file, it fixes that file with types only or narrows its own
  declaration (D16).
- **D21 Order with the other sessions.** `STRICT_NULL_BASELINE` (N0) merges
  after layering PR D, so no layering move or split meets a per-file key
  (Section 2.4). Until then, the S pull requests lower the diagnostics without
  a ratchet. Section 5.2 gives the sequence.
  - Pull requests that touch the same files open one after another. Before
    opening one, this session sends its file list to the other sessions with
    `SendMessage`.
  - Each S branch is rebased onto `dev` after the previous merge, and its
    counts are measured again.
  - If every file reaches 0 before N0 could merge, N0 is dropped, and the
    final pull request turns `strictNullChecks` on directly.
- **D22 Compiler options in their own pull requests.** The tsconfig change
  for the unused checks follows the runtime pull request U as its own
  GOVERNANCE pull request. It only tightens the guard. The final change to
  `strictNullChecks` is a separate GOVERNANCE pull request too (Section 4.3).
- **D23 Runtime fixes stay behavior-neutral.** An S pull request changes
  behavior only for a (c) fix. That fix has a test and an OV number, and the
  body lists it. Locally, an S pull request runs only the Playwright specs
  that exercise its (b) and (c) files. CI runs the full functional suite in
  shards.

## 4. Guard

### 4.1 Unused declarations (pull request G-U)

`tests/web/types/tsconfig.json` gains `"noUnusedLocals": true` and
`"noUnusedParameters": true`. The existing assertion "tsc reports no
diagnostic" then covers them. Nothing else changes.

### 4.2 strictNullChecks ratchet (pull request N0, authority-only)

- `tests/web/types/tsconfig.strict-null.json`:
  `{ "extends": "./tsconfig.json", "compilerOptions": { "strictNullChecks": true } }`.
- `tests/web/strict-null-ratchet.test.mjs` asserts:
  1. the strict config is exactly that;
  2. one `tsc` run reports no diagnostic outside `gbdraw/web/js/`, and every
     line of its output parses;
  3. `STRICT_NULL_BASELINE` is sorted, and every key is a module under
     `gbdraw/web/js/` outside `workers/` with a positive count;
  4. per file, the count equals its entry, or 0 when the file has no entry
     (D17).

  It runs the pinned `typescript`, as `tests/web/typed-boundaries.test.mjs`
  does.
- `tools/web-design-rule-guards.json`, R14: the new guard path, and the
  allowlist `{ path: tests/web/strict-null-ratchet.test.mjs, symbol:
  STRICT_NULL_BASELINE, kind: count-map }`.
- R14 in `gbdraw/web/CLAUDE.md`: one paragraph that says
  `strictNullChecks` is enabled per module, that the baseline may only shrink,
  that a module without an entry has no diagnostics, and that an addition is
  an authority-only change. The guard sentence names the new test. Nothing
  else in the rule changes.
- The baseline is generated from `dev` immediately before the pull request
  opens. If `dev` changes counts while the pull request waits for approval,
  the counts are regenerated, and the body says which entries changed.

The pull request adds a registered baseline, which is an authority change. It
is not auto-merged; it waits for the Owner's approval. After it merges, every
other session is told that a pull request that raises a count now fails, and
how to fix it: declare the type, or handle the null; an entry is never raised.

### 4.3 Final (pull request T, authority-only)

When every file is at 0:

- `tests/web/types/tsconfig.json` gains `"strictNullChecks": true`;
- the ratchet test, its baseline, and `tsconfig.strict-null.json` are
  deleted;
- R14 drops the paragraph of Section 4.2, and the registry drops the guard
  and its allowlist.

The existing typed-boundaries guard then enforces 0 for every module. The
pull request removes a registered guard. The guard it removes checks a weaker
condition than the one that replaces it, so the pull request contracts and
can auto-merge. The Gate's own classification decides: if the Gate asks for
Review, the pull request waits for the Owner.

If N0 was dropped (D21), pull request T only changes `tsconfig.json`.

## 5. PR sequence

### 5.1 Common rules

- Each runtime pull request: STANDARD, "This is not architecture-bearing",
  R13 applies, the owner-graph baseline only shrinks.
- Before opening a pull request:
  - `node --test tests/web/typed-boundaries.test.mjs tests/web/owner-graph-baseline.test.mjs`,
    plus the ratchet test once it exists;
  - `node tools/check-web-change-budget.mjs --base origin/dev`: Gate PASS;
  - the trusted-base check: `node tools/check-web-change-budget.mjs --base
    origin/dev --head <sha>` in a `dev` worktree. It diffs the `dev` tree
    against the head tree (two-dot), so a baseline that `dev` lowered after
    the branch point reads as an expansion. When that happens, rebase onto
    `dev`;
  - the fast Web suite once, under the shared heavy lock;
  - the Python tests that name a changed Web file, the `tools/` and `gbdraw/`
    Python files that name one (checked by reading how they parse it), and
    once per pull request the whole non-slow Python suite (phase 2:
    `tools/benchmark_protein_comparison.py` failed in CI);
  - the comment-only check: in a slice, every file whose code changed is a (b)
    or (c) entry in the body;
  - `node tools/check-pr-language.mjs`.
- CI load rules (Owner, 2026-10-07): at most 2-3 pull requests of this session
  in CI at once. Pull requests that touch the same files are serialized:
  stack the next on the previous, and open it after the previous merges.
  `gh pr update-branch` only on a conflict, to take a needed fix, or after a
  failed trusted-base check.
- Stacks are linear (rebased or cherry-picked commits, no merge commits).
  Merge commits gave GitHub two merge bases in phase 2 and a false conflict.
- A body edit re-runs the trusted-base check against the `dev` of that moment,
  so the body is final before the last push.
- STANDARD and contraction pull requests are auto-merged. A pull request that
  adds a baseline or allowlist entry, loosens rule text, or removes a guard
  registration waits for the Owner, except T as Section 4.3 describes.

### 5.2 Order

| Step | Pull request | Class | Waits for | Auto-merge |
| --- | --- | --- | --- | --- |
| 1 | This plan | STANDARD, docs only | — | yes |
| 2 | U: remove the 44 unused declarations | STANDARD, runtime | — | yes |
| 3 | G-U: unused options on (Section 4.1) | GOVERNANCE | U merged | yes, unless the Gate asks for Review |
| 4 | S1: N1, N2 | STANDARD, runtime | before layering A, which rewrites the import lines of `services/error-normalization.js` importers | yes |
| 5 | S2: N3, N4, N5 | STANDARD, runtime | S1; serialized with Phase E pull requests that edit `services/config.js` (layering A and B, OV-80 to OV-82) | yes |
| 6 | S3: N6, N7, N8, N9 | STANDARD, runtime | layering A merged (S3 then edits the moved paths) | yes |
| 7 | S4: N10, N11, N12 | STANDARD, runtime | S3; serialized with Phase E pull requests that edit `app/app-setup.js` | yes |
| 8 | N0: ratchet (Section 4.2) | GOVERNANCE, authority | layering D merged and at least one file still above 0; dropped otherwise | no, Owner |
| 9 | T: final (Section 4.3) | GOVERNANCE, authority | every file at 0 | Section 4.3 |

If S4 does not reach 0, phase 3 stops after N0 and reports the remaining files
and counts for phase 4.

## 6. Interaction with in-flight work

| In flight | Overlap | Handling |
| --- | --- | --- |
| Phase E layering A-D (Section 2.4) | moves, splits, or edits 24 files with 150 diagnostics; rewrites import lines in most `app/` and `services/` modules | D21: N0 after D. S1 before A. S3 and S4 after A, on the moved paths. File lists exchanged before each pull request |
| Phase E OV-65 (#896), OV-80 to OV-82 | `services/config.js`, `state.js`, `app/app-setup.js`, Legend modules | serialized with S1 (`state.js`), S2 (`services/config.js`), S3 (N7 Legend), S4 (`app/app-setup.js`) |
| Phase E removals of unused imports | U removes 3 import specifiers that layering A and B rewrite | Phase E drops them when it integrates A and B (agreed 2026-10-07) |

Check overlaps again before opening each pull request: `gh pr list` and
`gh pr diff --name-only <n>`.

## 7. Acceptance

- `tests/web/types/tsconfig.json` has `noUnusedLocals` and
  `noUnusedParameters`, and the typed-boundaries guard passes on `dev`.
- Either `tests/web/types/tsconfig.json` has `strictNullChecks: true` and the
  guard passes on `dev`, or N0 is merged and the final report names the
  remaining files and counts.
- Every S pull request lists its (b) and (c) changes. Every (c) has an OV
  number and either a test or a separate runtime pull request.
- No cast in the S pull requests removes `null` without a reason line (D19).
- The owner-graph baseline and `LAYER_IMPORT_BASELINE` have not grown.
- The final report lists the merged pull requests, the counts (unused 44 to 0;
  strictNullChecks 1,114 through each pull request to the end), the number of
  (a) changes, the (b) and (c) lists, the Owner-delegated choices, lessons for
  the CI practice, and a phase 4 recommendation for `noImplicitAny` with a new
  measurement and cost.

## 8. Rollback

- U and the S pull requests revert as ordinary STANDARD pull requests while
  no ratchet exists. After N0, a revert would raise counts, so they are fixed
  forward.
- G-U, N0, and T revert only with the Owner's approval, because each revert
  relaxes a check.
- A (c) fix reverts with its test.

## 9. Owner decisions

None blocks. The Owner may want to override these delegated choices:

1. **OD-5 D21, ratchet after the layering moves.**
   - (A, adopted) N0 after layering D. The S pull requests before it run
     without a ratchet, so another pull request could add diagnostics
     unnoticed until the next measurement.
   - (B) N0 now, and every layering pull request first brings its moved files
     to 0. That adds typing work to behavior-neutral moves.
   - (C) Teach the Gate to read a moved key as unchanged: a checker change in
     a checker-only pull request, then an authority pull request.
2. **OD-6 D20, four combined pull requests.** The alternative is one pull
   request per slice, 12 in all, each with a full CI run.
3. **OD-7 D21, dropping N0 when every file reaches 0 first.** The alternative
   lands N0 anyway for the interval before T.
4. **OD-8 Section 4.3, T auto-merges when the Gate allows it.** T removes a
   guard registration, which needs the Owner's approval as a rule. Here the
   removed guard is superseded by a stricter check in the guard that stays.
   The alternative is to wait for the Owner's approval.

## 10. Non-goals and phase 4 input

- `noImplicitAny` and the rest of `strict` stay off.
- No `.ts` source, build step, or `.d.ts` under `gbdraw/web/`.
- Workers and the contract between the Vue template and the setup return stay
  unchecked.
- No type for Python-owned option fields (R7).
- Phase 4 input at `b58176bd`: `noImplicitAny` adds 8,811 diagnostics in 158
  files. The largest are `services/session-request.js` 582,
  `services/config.js` 466, `app/app-setup.js` 406, `app/run-analysis.js` 404,
  and `app/circular-track-slots.js` 397. The final report measures it again.
