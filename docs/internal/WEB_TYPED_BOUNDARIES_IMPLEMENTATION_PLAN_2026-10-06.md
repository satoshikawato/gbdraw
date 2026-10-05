# Web typed boundaries implementation plan

Status: proposed, 2026-10-06. The Owner approved the approach on 2026-10-06:
JSDoc types and `// @ts-check`, checked by `tsc --noEmit`, no build step, no
TypeScript migration, introduced as a ratchet in which the set of checked
modules only grows.
Baseline: `dev` at `bcbf1d7f` (merge of #828). Slice (a) is also measured at
`e83b0b4d`, the tip of the unmerged Phase E stack (E3 `4fdda97b`, E5
`844cbc15`, E4 `797bf4f4`, E9a, E9b), because E9 rewrites the factory
parameters that slice (a) types.
Scope: `gbdraw/web/js/**` except `workers/`, `gbdraw/web/CLAUDE.md`,
`tools/web-design-rule-guards.json`, `tests/web/typed-boundaries.test.mjs`,
`tests/web/types/`, `package.json`, `package-lock.json`.

## 1. Objective

Have a compiler, not review, check the boundaries the owner-coupling work
introduced (ports between owners, R13) and the persisted contracts (Session 45,
canonical request schema 9, feature catalog schema 5). After this plan, a
missing port, an argument the receiver no longer reads, a port with the wrong
signature, or a writer that emits a key its contract lacks fails a fast node
test. Types are JSDoc comments in the `.js` module that owns them. Shipped files
change only by comments and JSDoc-cast parentheses.

The owner-coupling plan's non-goal "No TypeScript migration" still holds: no
`.ts` source and no build step are added. Its detectors see the shape of the
owner graph (who receives which owner or port). They do not see signatures
(the names, arity, and argument types at a port), and this plan adds that check.

## 2. Evidence

### 2.1 Tooling, Gate, and CI

| Question | Finding | Evidence |
| --- | --- | --- |
| Is TypeScript available? | No. The shared `node_modules` holds `@playwright`, `@types`, `playwright`, `playwright-core`, `undici-types`. `gbdraw/web/js` has no `@param`, `@typedef`, or `@ts-check`. | `ls /home/kawato/gbdraw-work/node_modules`; `grep -rn '@param' gbdraw/web/js` → 0 |
| Registry | `latest` is 7.0.2, the native compiler; `npm ci` installs one `@typescript/typescript-<platform>` optional package. 6.0.3 is the last JavaScript release. | `npm view typescript dist-tags` (2026-10-06) |
| Lockfile churn | `npm install --package-lock-only --save-dev --save-exact` in a scratch copy: 7.0.2 adds 377 lockfile lines (20 platform entries), 6.0.3 adds 16. `package.json` gains one `devDependencies` line. | scratch copy, removed |
| Web change Gate | A `devDependencies` addition is listed under "Dependency changes" and does not fail the Gate. Only additions to `dependencies`, `optionalDependencies`, or `peerDependencies`, and bare production imports, fail it ("new production dependencies are not allowed"). No other npm or lockfile rule exists. | `tools/check-web-change-budget.mjs:1613-1661,1820-1822`; `WEB_CHANGE_POLICY.md:146` |
| CI impact class | `package.json` and `package-lock.json` are `packaging`, which runs the full PR tier. | `tools/ci-impact-policy.mjs:26,42,174` |
| Do the jobs that run `tests/web/*.test.mjs` install devDependencies? | Yes. `Browser` (dev push) and `Web contracts` (PR) both run `npm ci` without `--omit` or `NODE_ENV`, then run every `tests/web/*.test.mjs` through `find`. `Web change budget` installs nothing and runs only `architecture-contracts.test.mjs`. `deploy_web.yml` runs one Gallery test. A new `tests/web/typed-boundaries.test.mjs` therefore runs in both jobs, and no workflow change is needed. | `.github/workflows/test.yml:382-385,420-426,480-483,518-525,97-145`; `grep -rn 'omit\|NODE_ENV' .github/workflows` → none |
| Where may typedefs live? | A `.d.ts` or `tsconfig.json` under `gbdraw/web/` is copied to Cloudflare Pages (`shutil.copytree(WEB_ROOT, …)`). It is excluded from the wheel and sdist, which take only `*.js` under `web/js`. ci-impact classifies a non-`.js` path under `gbdraw/web/js/` as `full`. | `tools/prepare_cloudflare_pages.py:150`; `gbdraw/_build_support.py:31-43`; `MANIFEST.in:9`; `tools/ci-impact-policy.mjs:192-193,210` |
| Do JSDoc type imports create edges? | No. The import scanner masks comments, so `@import` and `import('…')` in JSDoc add no import edge, cause no cycle, and make no request (CSP is unaffected). The owner-graph detectors also mask comments. | `tools/web-change-source.mjs:88-106`; `tools/web-owner-graph-detectors.mjs:14` |
| Registered-allowlist kinds | `count`, `count-map`, `set`, `writer-map`. A `set` contracts when entries are removed, and that contraction may ship with runtime changes. Adding an entry is an expansion and needs an authority-only PR. | `tools/check-web-change-budget.mjs:479,638-673`; `WEB_CHANGE_POLICY.md:337-365` |

### 2.2 Timing (local, 32 cores, Node 26)

| Run | TS 7.0.2 | TS 6.0.3 |
| --- | ---: | ---: |
| Whole tree (166 modules), `checkJs: true` | 0.33 s, 260 MB | 2.2 s, 430 MB |
| Prototype config, 7 modules with `// @ts-check` | 0.12 s | 0.83 s |
| Same, 29 slice-(a) modules | 0.27-0.29 s | 1.6 s |
| `node --test` guard that spawns `tsc` (7 modules) | 0.18 s wall | 0.89 s wall |

`tsc` follows imports, so the program always parses all 166 modules. The cost
barely depends on how many modules are checked. The fast-suite step limit is 5
minutes.

### 2.3 Error counts

With `allowJs` and `checkJs: false`, only modules whose first line is
`// @ts-check` report errors. Imported unchecked modules are parsed for type
inference only: in the prototype, 7 checked modules reported 3 errors, all in
checked modules. A pragma placed after code is ignored, and one placed after a
leading comment is honored. TypeScript 6 and later default to `strict: true`,
so the config must set `strict` explicitly.

The whole tree with `checkJs: true` (every module as if checked), workers
excluded:

| Configuration | `dev` `bcbf1d7f` | E9 tip `e83b0b4d` |
| --- | ---: | ---: |
| `strict: false` | 575 (in 65 of 166 modules; 101 have none) | 575 |
| `strict: true` | 10,737 (6,126 are TS7006 implicit-any parameters) | 10,753 |
| `strict: false` + `strictNullChecks` | — | 2,150 |
| `strict: false` + `noImplicitAny` | — | 10,266 |
| `strict: false` + `noUnusedLocals` | — | 634 |
| `strict: false` + `strictFunctionTypes`, `strictBindCallApply`, `noImplicitThis`, `noFallthroughCasesInSwitch` | — | 575 (+0) |

TS 6.0.3 reports nearly the same sites: 572 at `dev`; at the E9 tip, 564 sites
are common, 11 are reported only by TS 7, and 8 only by TS 6. `lib: esnext` removes the two `Uint8Array.setFromBase64`
errors in `services/byte-utils.js`, which the module feature-detects.

The `strict: false` errors fall into four kinds:

1. Port drift, which are real defects. At `dev`, `app/feature-editor.js:42-45` passes
   `nextTick` to `createFeatureColorActions`, which never reads it, and
   `app/legend.js:22` passes `{ state }` to `createLegendLayoutActions()`, which
   takes no argument.
2. Inference artifacts. In a parameter written `({ a, b = 1 } = {})`, a binding
   without a default is dropped from the inferred type, so a caller that passes
   it fails (`app/legend-layout.js:106` at `dev` → `captureDecorationContinuity({ canonical, … })`).
   A JSDoc type on the receiver removes the artifact.
3. DOM narrowing (`Element` versus `HTMLElement` for `.focus()` and
   `.dataset`), fixed by a JSDoc cast.
4. Object shapes inferred from a literal and extended later, which account for
   most of the 125 errors in `services/session-request.js`.

A demonstration with a receiver typedef of ports under `strict: false` plus
`strictFunctionTypes` caught four cases: a missing port (TS2345), an extra
argument (TS2353), a port with an extra required parameter (TS2322), and a port
used with the wrong argument type (TS2365).

### 2.4 Slices

Counts are at `e83b0b4d`. Errors are shown as `strict: false` / `strict: true`.
"Params" is the number of parameter names of exported `create*` and `setup*`
factories, and each name is one JSDoc `@property` line.

| PR | Slice | Modules | Factories (params) | Errors | Lines | Estimated net additions |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| A1 | (a) Legend owners | 14 | 10 (43) | 33 / 654 | 5,458 | ~120 |
| A2 | (a) Feature-editor owners | 12 | 11 (90) | 17 / 991 | 7,087 | ~150 |
| A3 | (a) Track slots, record display, alignment, feature search | 12 | 11 (63) | 34 / 1,489 | 8,532 | ~130 |
| A4 | (a) Composition root and Generate | 3 | 4 (58) | 94 / 1,265 | 11,299 | ~150 |
| B1 | (b) Feature catalog, Feature placement rows, session readers | 8 | 6 (2) | 0 / 552 | 4,944 | ~120 |
| B2 | (b) Canonical request schema 9 | 2 | 0 | 130 / 798 | 5,670 | ~200 |
| B3 | (b) Session 45 | 2 | 0 | 30 / 621 | 5,613 | ~90 |
| C1 | (c) Pure services (no DOM) | 30 | 8 (6) | 12 / 568 | 4,588 | ~60 |
| | Total | 83 of 166 | | 350 / 6,938 | | |

Many A4 and B errors are inference artifacts of kind 2, which disappear once
the receivers are typed. The estimates assume that. Module lists:

- A1: `app/legend.js`, `app/legend-layout.js`, `app/legend/*` (6),
  `app/legend-layout/*` (6).
- A2: `app/feature-editor.js`, `app/feature-editor/*` (9), `app/svg-styles.js`,
  `app/rule-matching.js`.
- A3: `app/track-slot-edits.js`, `app/circular-track-slots.js`,
  `app/circular-track-slots/*` (2), `app/linear-track-slots.js`,
  `app/record-display-options.js`, `app/record-display/*` (2),
  `app/similarity-alignment.js`, `app/feature-search/*` (3).
- A4: `app/app-setup.js`, `app/run-analysis.js`, `app/watchers.js`.
- B1: `services/feature-catalog.js`, `services/standalone-interactivity.js`,
  `services/feature-placement.js`, `services/feature-edit-migration.js`,
  `services/session-file.js`, `services/session-authority.js`,
  `services/session-active-config-contract.js`,
  `services/gallery-session-publication.js`.
- B2: `services/session-request.js`, `services/session-resources.js`.
- B3: `services/config.js`, `services/gallery-session-migration.js`.
- C1: the services without DOM use, outside B: `bounded-json-transport`, `byte-utils`,
  `canonical-comparisons`, `canonical-resource-references`,
  `comparison-warnings`, `current-worker-result-source`, `depth-file-codec`,
  `diagram-generation`, `diagram-resource-staging`, `diagram-worker-protocol`,
  `feature-identity`, `feature-override-identity`, `file-content-cache`,
  `history-files`, `imported-comparison-intent`, `json-clone`,
  `legacy-similarity-alignment`, `losat-runtime`,
  `main-session-comparison-frame`, `orthogroup-feature-metadata`, `reset`,
  `resource-payload-owner`, `result-normalization`, `runtime-capabilities`,
  `runtime-test-hooks`, `safe-object-keys`, `session-feature-metadata`,
  `session-import-client`, `session-resource-backing`, `svg-sanitization`.

Python already declares the persisted key sets for two contracts:
`_TOP_LEVEL_FIELDS_V5` (`gbdraw/session_request_codec.py:194-197`) and
`CURRENT_SESSION_TOP_LEVEL_FIELDS` (`gbdraw/session_io.py:62-84`).

## 3. Decisions

All are Owner-delegated: each recommended option is adopted under the standing
instruction of 2026-09-29. Section 9 lists the ones the Owner may want to
override before T0.

- **D1 New rule R14 "Typed boundaries"**, not an R13 extension. Its scope
  includes persisted contracts and services, beyond owner layering, and it has
  its own guard and allowlist. R13 keeps its owner-graph semantics.
- **D2 Ratchet as a shrink-only `UNCHECKED_MODULES` set.** The literal is
  registered with kind `set`. Removing entries is the contraction that the
  existing checker already allows alongside runtime changes. A grow-only
  "checked" literal would need a new allowlist kind, which means a checker-only
  PR followed by an authority PR, for the same effect.
- **D3 A new module is checked from its first commit.** A module that moves
  gets a new path and is treated as new: type it in the move PR, or register
  the new path first in an authority-only PR.
- **D4 TypeScript 7.0.2, exact pin** in `devDependencies`. It is about 7x
  faster than 6.0.3, its diagnostics nearly match on this code, and it is npm
  `latest`. An upgrade is its own PR, because a new compiler can add
  diagnostics.
- **D5 `strict: false`** plus the four strict-family flags that cost 0 today:
  `strictFunctionTypes` (which matters for function-typed ports),
  `strictBindCallApply`, `noImplicitThis`, and `noFallthroughCasesInSwitch`.
  `strictNullChecks`, `noImplicitAny`, and `noUnusedLocals` are deferred
  (Section 10). Boundary typing comes from declared JSDoc types, which are
  enforced without `noImplicitAny`. Guard assertion 3 makes declaring them
  mandatory for factories.
- **D6 Typedef location.** A port typedef lives in the receiving owner's
  module. A persisted-contract typedef lives in its writer, next to the version
  constant: `GbdrawSession` beside `SESSION_VERSION` in `services/config.js`,
  `CanonicalRenderRequest` beside `CANONICAL_REQUEST_SCHEMA` in
  `services/session-request.js`, and `FeatureCatalog` beside
  `FEATURE_CATALOG_SCHEMA` in `services/feature-catalog.js`. Other modules use
  `/** @import { T } from './x.js' */`. There is no `types.js` hub, which would
  be a cross-owner module that every owner imports, and no `.d.ts` under
  `gbdraw/web/`. The one exception is
  `tests/web/types/web-globals.d.ts`, outside the shipped tree. It declares the
  globals that vendored scripts install (`window.Vue`, `window.jspdf`,
  `DOMPurify`) and the test hooks (`__GBDRAW_HISTORY__`, `__GBDRAW_TEST_HOOKS__`,
  `__GBDRAW_LAST_LOSAT_TELEMETRY__`), all typed `any`.
- **D7 Persisted typedefs describe the current writer format only.** No
  `SessionV44`-style types exist (root `CLAUDE.md`, "Persisted-format
  compatibility"). A reader takes `Record<string, any>`, and a migrator returns
  the current typedef. Python owns the option field set (R7), so
  `diagramOptions` is `Record<string, any>`, except for the rows JavaScript
  builds itself (`featureOverrides` and `featurePlacements` from
  `services/feature-placement.js`).
- **D8 Workers stay out of scope.** They need `lib: webworker`, which
  conflicts with `dom` in one program.
- **D9 A missing `typescript` fails the guard** with "run `npm ci`". There is
  no skip, because a guard that skips is not a guard. `typescript` must be
  installed once in the parent clone's `node_modules` (Section 5, T0).
- **D10 Slice PRs are comment-only.** The only code change allowed is the
  parentheses of a JSDoc cast. A defect `tsc` exposes is logged as OV-xx. It is
  fixed in the slice PR only when it takes a few lines in a module the PR
  already checks (for example, removing a dead port argument), and the body
  names it.

## 4. Rule and guard (landed by T0)

### 4.1 R14 text for `gbdraw/web/CLAUDE.md`

```markdown
### R14: Typed boundaries

Modules under `js/` except `workers/` are checked by `tsc --noEmit` with JSDoc
types; there is no build step and no TypeScript source. A checked module starts
with the line `// @ts-check`. The modules not yet checked are listed in
`UNCHECKED_MODULES` in `tests/web/typed-boundaries.test.mjs`, which may only
shrink; a new module is checked from its first commit.

- A checked module declares the type of every parameter of its exported
  `create*` and `setup*` factories. A port (R13) is a function type declared by
  the owner that receives it; a composition root imports it with `@import`. A
  type naming another owner's factory result (`ReturnType<typeof createX>`) is
  a whole-object port.
- A persisted contract (Session 45, canonical request schema 9, feature catalog
  schema 5) has one typedef, in its writer, for the current format only. A
  reader takes unvalidated data and returns that typedef.
- Types are JSDoc in the module that owns them, never a `.d.ts` or a types
  module under `gbdraw/web/`. A type import follows the layers of R13.
- `@ts-ignore`, `@ts-expect-error`, and `@ts-nocheck` are not used: fix the
  type, cast with JSDoc, or fix the code.

Guard: `tests/web/typed-boundaries.test.mjs` (the checked set, the compiler
run, declared factory parameters, no suppression, type-import direction).
```

### 4.2 Guard assertions (`tests/web/typed-boundaries.test.mjs`)

1. Module set: every `.js` under `gbdraw/web/js/` except `workers/` either has
   `// @ts-check` as its first line or is in `UNCHECKED_MODULES`, and never
   both. Every entry names an existing module. Failure messages say what to
   change: "remove from `UNCHECKED_MODULES`" or "add `// @ts-check`".
2. Compiler: one run of
   `tsc -p tests/web/types/tsconfig.json --listFiles --pretty false`, spawned
   with `process.execPath` on the `bin/tsc` resolved from
   `typescript/package.json`. It must exit 0, print no diagnostic line, and
   list every checked module in the program. The last condition catches an
   `include` typo that would make the check vacuous.
3. Declared boundaries: in a checked module, every exported `create*` or
   `setup*` function has a JSDoc block with one `@param {T}` per declared
   parameter, where `T` is not `any`, `*`, `object`, or `Object`. Outside the
   composition roots (`WEB_OWNER_GRAPH_DEFAULTS.compositionRoots`, reused from
   `tools/web-owner-graph-detectors.mjs`), no JSDoc names
   `ReturnType<typeof create…>`.
4. No `@ts-ignore`, `@ts-expect-error`, or `@ts-nocheck` under `gbdraw/web/js/`.
5. Type imports (`@import … from '…'` and `import('…')` inside comments)
   resolve to an existing module under `gbdraw/web/js/`. `state.js`,
   `services/**`, and `utils/**` do not import types from `app/**`, and only a
   composition root imports types from a composition root.

### 4.3 `tests/web/types/tsconfig.json`

```json
{
  "compilerOptions": {
    "allowJs": true, "checkJs": false, "noEmit": true,
    "target": "es2022", "module": "esnext", "moduleResolution": "bundler",
    "lib": ["esnext", "dom", "dom.iterable"], "types": [], "skipLibCheck": true,
    "strict": false, "strictFunctionTypes": true, "strictBindCallApply": true,
    "noImplicitThis": true, "noFallthroughCasesInSwitch": true
  },
  "include": ["../../../gbdraw/web/js/**/*.js", "web-globals.d.ts"],
  "exclude": ["../../../gbdraw/web/js/workers/**"]
}
```

## 5. PR sequence

Common rules for every slice PR (A1-A4, B1-B3, C1):

- Change class STANDARD. "This is not architecture-bearing": the PR adds
  comments and cast parentheses, and moves no owner, path, or persisted format.
- Add `// @ts-check` as line 1 of each module. Type every exported factory
  parameter (ports as function types) and fix the remaining errors with JSDoc
  types or `/** @type {X} */ (expr)` casts. Remove the modules from
  `UNCHECKED_MODULES`; this is a contraction, so the Gate passes.
- Verification:
  ```bash
  node --test tests/web/typed-boundaries.test.mjs
  node tools/check-web-change-budget.mjs --base origin/dev   # Gate PASS; R14 contraction listed; owner-graph report unchanged
  node --input-type=module -e "
  import { execFileSync as x } from 'node:child_process';
  import { maskJavaScript as m } from './tools/web-change-source.mjs';
  const [b, h] = ['origin/dev', 'HEAD'];
  const code = (r, p) => { try { return m(x('git', ['show', r + ':' + p], { encoding: 'utf8', maxBuffer: 1 << 26 }), { strings: false }).replace(/[\s()]/g, ''); } catch { return null; } };
  const paths = x('git', ['diff', '--name-only', b, h, '--', 'gbdraw/web/js'], { encoding: 'utf8' }).split('\n').filter(Boolean);
  const changed = paths.filter((p) => code(b, p) !== code(h, p));
  console.log(changed.length ? 'code changed: ' + changed.join(' ') : 'comment-only: ' + paths.length);
  "   # expected: comment-only, or exactly the D10 fixes named in the body
  flock /home/kawato/gbdraw-baselines/owner-coupling-phase-e-20261005/heavy.lock sh -c "find tests/web -maxdepth 1 -type f -name '*.test.mjs' ! -name 'architecture-contracts.test.mjs' ! -name 'gallery-session-publication.test.mjs' -print0 | xargs -0 node --test --test-concurrency=8 > <log> 2>&1"
  ```
  No local Playwright run is needed, because behavior does not change. CI
  still runs the Web runtime jobs.
- Size review: REQUIRED is expected on net additions (JSDoc lines) and, for C1,
  on file count. The size checker does not fail the Gate.
- After merge: run `gh pr update-branch` on every open PR and re-arm
  auto-merge. The merge contracted a registered literal, so an open PR that
  is re-evaluated against the new `dev` sees its older `UNCHECKED_MODULES` as
  an expansion (Phase E HANDOFF, lesson of 2026-10-06).

### T0: register R14 and the guard (GOVERNANCE, authority-only)

Files:

- `gbdraw/web/CLAUDE.md`: heading "Design rules R1-R14", the R14 section
  (4.1), and in "Local build and verification",
  `npm ci` and `node --test tests/web/typed-boundaries.test.mjs`.
- `tools/web-design-rule-guards.json`: add
  `{ "id": "R14", "heading": "R14: Typed boundaries", "guards": ["tests/web/typed-boundaries.test.mjs"], "allowlists": [{ "path": "tests/web/typed-boundaries.test.mjs", "symbol": "UNCHECKED_MODULES", "kind": "set" }] }`.
- `tests/web/typed-boundaries.test.mjs` (new): the guard (4.2), with
  `UNCHECKED_MODULES` holding every non-worker module of the base (166 at
  `bcbf1d7f`). The policy requires a new registered literal to land in the same
  authority-only PR as its registration (`WEB_CHANGE_POLICY.md:362-365`).
- `tests/web/types/tsconfig.json` and `tests/web/types/web-globals.d.ts`
  (new).
- `package.json` (`"typescript": "7.0.2"`) and `package-lock.json`.

No production path changes, and no checker or workflow change is needed.

Verification:

```bash
npm ci
node --test tests/web/typed-boundaries.test.mjs tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs --base origin/dev   # Gate PASS; Review REQUIRED (governance); "Dependency changes" lists devDependencies typescript
# Negative checks in a scratch commit, then drop it:
#   `// @ts-check` on app/legend.js -> assertion 1 names it; after removing it from the list, assertion 2 reports the dead `{ state }` argument
#   `// @ts-check` on services/json-clone.js without removing it from the list -> assertion 1 says "remove from UNCHECKED_MODULES"
```

After merge: install `typescript` in the shared `node_modules`. Either run
`npm ci` in `/home/kawato/gbdraw-work` from a checkout that contains T0, or
run `npm install --no-save typescript@7.0.2` there. Then run
`gh pr update-branch` on every open PR, because T0 changes
`gbdraw/web/CLAUDE.md`.

Proposed title: `Add R14 typed boundaries with a shrink-only unchecked-module list`

### A1-A4: ports, slice (a) (STANDARD)

| PR | Modules (2.4) | Specific content |
| --- | --- | --- |
| A1 | Legend owners | Typedefs for the options of `createLegendManager`, `createLegendLayout`, and the eight sub-owner factories. Ports include `commitLegendRowRules`, `beginHistoryTransaction`, `commitHistoryTransaction`, `commitActiveResultEdit`, and `readActiveResultIdentity`. Remove the dead `{ state }` argument of `createLegendLayoutActions` (D10). |
| A2 | Feature-editor owners | `createFeatureEditor` options, including `projectPaletteAndRules` and `projectFeatureEdits`; the `editorPorts` bag as a typedef of `applyFeatureVisibilityToLabels` and `syncLabelEditor`; the rule, color, label, SVG, visibility, placement (`changeTrackLayout` consumer), and table factories; `createSvgStyles`; `createRulePreparation`. Remove the dead `nextTick` argument of `createFeatureColorActions` (D10). |
| A3 | Track slots, record display, alignment, search | `changeTrackLayout` (R10) at both slot editors and `track-slot-edits.js`; `runRecordRotation` (`app/record-display/feature-record-rotation.js`); `runRecordAlignment` (`app/similarity-alignment.js`); `createRecordDisplayControls`; `createPreviewFeatureSearch`. |
| A4 | Composition root and Generate | `app/app-setup.js` becomes checked, so the provider side of every port is verified, including `legendRowRulePorts` and `paletteRulePorts`; `createRunAnalysis` (31 parameters) and `setupWatchers`. |

Dependency: A1-A4 start after E9 merges (Section 6). A4 comes last because its
errors shrink as A1-A3 type the receivers.

### B1-B3: persisted contracts, slice (b) (STANDARD, session-persistence CI class)

| PR | Content | Parity |
| --- | --- | --- |
| B1 | `FeatureCatalog` (schema 5) in `services/feature-catalog.js`; `migrateLegacyFeatureCatalog`, `admitFeatureCatalog`, and `validateFeatureCatalog` return it, and the catalog writer in `services/standalone-interactivity.js` uses it. Feature placement and override row types in `services/feature-placement.js`. The active-config contract type in `services/session-active-config-contract.js`. The reader input types of `services/session-file.js`, `services/session-authority.js`, and `services/gallery-session-publication.js`. | — |
| B2 | `CanonicalRenderRequest` (schema 9) and the `{ renderRequest, resources, webFiles }` envelope in `services/session-request.js`. `tsc` checks the writer literal (`session-request.js:2656-2675`) against them. Readers and the schema promotion return the current type. | New pytest beside the codec tests: the typedef's property names equal `_TOP_LEVEL_FIELDS_V5`. |
| B3 | `GbdrawSession` (version 45) in `services/config.js`, referring to `CanonicalRenderRequest`, the active-config type, and `FeatureCatalog` instead of repeating render fields. The save literal (`config.js:4112-4157`) is checked. `migrateSessionDataToCurrent` and `services/gallery-session-migration.js` return `GbdrawSession`. | The same pytest: the property names equal `CURRENT_SESSION_TOP_LEVEL_FIELDS`. |

Dependency: after the OV-38/OV-40 and OV-39 fixes merge. B1 also needs E9,
because E5 and E9 change `services/feature-placement.js`.

### C1: pure services, slice (c) (STANDARD)

C1 covers the 30 modules in 2.4: 8 factories, 12 errors at `es2023`, and 10 at
`esnext`. Most modules need only the pragma. No in-flight PR touches them, so
C1 may land right after T0.

Order: priority is a, then b, then c. C1 may go first while E9 is still
landing, because it is independent.

## 6. Interaction with in-flight work

| In-flight | Overlap | What waits |
| --- | --- | --- |
| E3, E5, E4, E9a/E9b (stack at `e83b0b4d`) | The stack changes 24 Web modules, nearly all of them slice-(a) factories; E5 and E9 change `services/feature-placement.js` | A1-A4 and B1 start after E9 merges. Typing earlier would be rewritten in the same lines. T0 should merge after the stack, because an authority merge forces `update-branch` on every stacked PR. |
| #827 and `governance/r10-q3-change-track-layout-port` (R10 and R3 text) | both edit `gbdraw/web/CLAUDE.md` | T0 merges after both and is rebased onto them. Only one authority change is in flight at a time. |
| OV-38/OV-40 (session legacy readers), OV-39 (`app/session-feature-metadata.js`, `services/error-normalization.js`) | Session readers | B1-B3 wait for OV-38/OV-40. OV-39 touches no slice module. |
| OV-36 (after E5), Gallery refresh, #829 | none | — |
| Any PR that adds a module after T0 | D3 | The new module starts with `// @ts-check` and passes the guard. No in-flight branch adds a module today (`git diff --diff-filter=A` on the stack and on OV-39). |

## 7. Acceptance

- T0: the guard runs in `Web contracts` and `Browser` in under 2 s, and the two
  negative checks fail with the documented messages.
- After A1-A4: the 41 slice-(a) modules are checked; every port-receiving
  factory declares its ports; `app/app-setup.js` is checked; the dead
  arguments are gone.
- After B1-B3: the three persisted contracts have one current-format typedef
  each, at the writer; readers and migrators return it; the two pytest parity
  checks pass.
- After C1: `UNCHECKED_MODULES` has 83 entries (166 - 83), and the rest are
  input for later work.
- Every slice PR passes the comment-only check, apart from the D10 fixes its
  body names.

## 8. Rollback

- A slice PR is fixed forward. Reverting one would add entries back to a
  registered literal together with runtime paths, and the Gate fails that
  (`design-rule.co-change`). That is the ratchet working as intended: no
  runtime PR un-checks a module.
- To retire the mechanism, revert T0 in an authority-only PR. The checker
  treats the removal of the literal as a contraction. The pragmas and JSDoc
  that remain are inert comments, and removing them is an optional later
  runtime PR.
- If TS 7 misbehaves, pin 6.0.3 in one PR (`package.json`,
  `package-lock.json`). It reports nearly the same sites with these options.

## 9. Owner decisions

None blocks. The Owner may want to override these delegated choices before T0:

1. **D4 compiler.** 7.0.2 (native, 0.1-0.3 s, 377 lockfile lines) or 6.0.3
   (JavaScript, 0.8-1.6 s, 16 lockfile lines, stable compiler API).
2. **D6 one ambient `.d.ts`** in `tests/web/types/`. It is the only file in
   TypeScript syntax, it is outside the shipped tree, and it declares only
   vendored globals and test hooks. The alternative is JSDoc casts at about 34
   use sites.
3. **D3 new modules checked from their first commit.** The alternative is to
   let new modules stay unchecked, which needs a grow-only literal kind (a
   checker-only PR and then an authority PR) to keep the ratchet.

## 10. Non-goals

- No `.ts` source, transpile, emitted file, or bundler. Shipped JavaScript
  changes only by comments and cast parentheses.
- No `strict: true`, `strictNullChecks` (+1,575 errors), `noImplicitAny`, or
  `noUnusedLocals` (+59) in this plan. `tsc` has no per-file strictness, so a
  later plan would add a second program whose diagnostics are filtered to a
  shrinking set.
- Workers (5 modules) and the contract between the Vue template and the setup
  return are not checked.
- No type for Python-owned option fields (R7).
- No change to workflows or to selective CI. A comment-only Web change still
  runs the Web runtime jobs.
- No editor integration such as a root `jsconfig.json`. VS Code checks
  `// @ts-check` files with its default inferred project.
