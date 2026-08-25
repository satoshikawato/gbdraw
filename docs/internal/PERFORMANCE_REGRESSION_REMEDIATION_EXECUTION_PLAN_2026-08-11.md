# Performance regression remediation execution ledger

- Date: 2026-08-11
- Authoritative baseline: `3bca0e8d98d7743e5df4f30dfb33a0fa34ddbf2f`
- Scope executed here: Phase A P0-3/P0-4/P0-5 guardrails, followed by P0-1
  History ownership and P0-2 SVG-ingestion remediation
- Source audit: `GBDRAW_PERFORMANCE_REGRESSION_AUDIT_2026-08-11.html`
- Task briefs: `CODEX_TASK_PHASE_A_PERFORMANCE_GUARDRAILS.md` and
  `CODEX_TASK_P0_WEB_PERFORMANCE_REGRESSION_REMEDIATION.md`

This ledger reconciles the audit with the current HEAD and worktree. A status
describes the underlying production or verification contract, not whether the
audit once documented an incident. P0-1 and P0-2 were executed after the
validated Phase A work; later P1/P2 packages remain proposals.

## Authoritative baseline

- `git status --short --untracked-files=all`: only the user-provided,
  untracked `docs/internal/CODEX_TASK_PHASE_A_PERFORMANCE_GUARDRAILS.md` existed.
- `git diff --stat`: no tracked changes.
- `git rev-parse HEAD`: `3bca0e8d98d7743e5df4f30dfb33a0fa34ddbf2f`,
  which is also the audit baseline.
- Node Playwright discovery found 83 tests in 11 specifications. The existing
  Browser CI job installed Playwright but invoked no Playwright runner.
- The existing 25,000-feature performance specification passed locally before
  implementation: 5 tests passed with one worker and zero retries. Its
  mutation-count assertions remained unchanged.
- The old offline GUI mega-smoke failed locally before implementation after
  189.99 seconds at the final Linear cache wait in
  `tools/verify_gui_offline.py`, not in its palette assertions or export checks.
  The helper wrote the retired `app.blastSource` field, while the current
  `linearComparisonPlan.mode` stayed `none`. It therefore completed a
  no-comparison render and waited for a LOSAT cache entry that could not exist.
  This verifies the audit's late-stage diagnosis but does not prove a product
  LOSAT runtime hang.

## Pre-existing CI failures

These failures belong to the unmodified baseline and are recorded separately
from post-change verification.

- Latest main Tests run: GitHub Actions run `31466230978`, HEAD
  `3bca0e8d98d7743e5df4f30dfb33a0fa34ddbf2f`, conclusion `failure`.
- Failed checks: Slow Python 3.10 job `93699585890`, Slow Python 3.11 job
  `93699585856`, and Slow Python 3.12 job `93699585894`. Each failed step 6,
  `Run slow tests`.
- All other 14 jobs passed, including Browser Python 3.11 job `93699585912`.
  That Browser job's passing status did not represent Playwright coverage,
  because its step inventory ran only Node unit tests and Python browser tests.
- Public check annotations expose only exit code 1; unauthenticated GitHub log
  download is unavailable. The exact offline-smoke traceback is independently
  reproduced above. A separate pre-existing warning reports the Node 20 action
  runtime deprecation for `actions/checkout@v4` and `actions/setup-python@v5`.

## Audit incident reconciliation

"CI" below describes the workflow after this Phase A patch. Where that differs
from the authoritative baseline, the cell says so explicitly.

| Item | Status | Current owner | Existing regression test | CI | Missing structural invariant | Target CI lane | Proposed PR-sized work package | Dependencies | Expected generated-artifact impact |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| I-01 | Fixed for the audited Circular coordinate path; not a cross-mode claim | `gbdraw/svg/circular_tracks.py`; `gbdraw/analysis/skew.py` | GC-content/regression/output comparison tests; quick real-assembler smoke | Core 3.10/3.11/3.12; benchmark itself is not CI | Deterministic N/2N work or timing bound for Circular GC/skew | Quick structural Core; standard main/nightly benchmark | Add synthetic scale cases and base/head JSON | P2-2 benchmark lane | Benchmark JSON only; no reference SVG change expected |
| I-02 | Fixed for expensive preprocessing, weakly guarded | `gbdraw/labels/filtering.py`; feature factory and Circular/Linear label candidate owners | Label-override semantics, precalculated-label reuse, and no-label short-circuit tests | Core matrix | Filtering prepared once per diagram and label selection/measurement bounded once per feature and purpose | Quick structural Core | Introduce one prepared label-filter/candidate owner and call-count tests | Coordinate with I-03 prepared-record boundary | None expected; output comparison must remain unchanged |
| I-03 | Fixed | `gbdraw/diagrams/linear/precalc.py`, Linear assembler, and `render/groups/linear/seq_record.py` | `test_linear_feature_layers_are_built_once_per_record`; reference outputs | Core matrix | Fallback branches still permit rebuilding when prepared inputs are omitted | Core architecture/structural | Require one typed prepared-linear-record aggregate and retire rebuild fallbacks | Coordinate with I-02 | None expected |
| I-04 | Fixed with strong fixed-fixture guards | `gbdraw/layout/spatial.py`; Circular/Linear label and feature-track integrations | Spatial index tests, indexed/full-scan parity, collision-call caps, origin-span tests | Core matrix | Generated N/2N/4N candidate-work growth bound | Quick structural Core; standard main/nightly benchmark | Add deterministic sparse/dense scaling fixtures and persist benchmark report | P2-2 | Benchmark JSON only |
| I-05 | Primary listener/lookup symptom fixed; index lifecycle partially fixed | `web/js/app/feature-editor/svg-actions.js`, `state.js`, and `app/watchers.js` | Feature selection unit test and real feature-click Playwright coverage | Unit test in Browser; Playwright coverage now runs in the functional lane (dormant at baseline) | One listener set and one DOM-index build per root mount; same-root metadata updates must not rebuild | Node structural plus Playwright functional/performance | Stop forced same-root rebuilds and add listener/query-count spies | P0-3; aligns with P1-2 | None |
| I-06 | Partially fixed | `web/js/services/diagram-generation.js`, diagram Worker, `app/run-analysis.js`, and LOSAT service | Cancellation, stale-result, rollback, and browser render-cancel tests | Browser Node lane and dedicated LOSAT browser acceptance | Warm Worker reuse after cancellation and a latest-only broker/queue | Node lifecycle plus canonical functional/LOSAT acceptance | Add an explicit single-Worker job state machine and cooperative cancellation | P0-5 owner guard; automatic rerender still depends on P0-1/P0-2 | None |
| I-07 | Partially fixed | `_build_support.py`, `MANIFEST.in`, `tools/prepare_cloudflare_pages.py`, Cloudflare Worker, and `deploy_web.yml` | Package exclusion, Cloudflare manifest/bundle, immutable-ref, and route tests | Core/Browser, although some packaging coverage is duplicated on main | One canonical deploy owner; expanded-size/max-file budget; remote fetch/range/cache/error contract | Canonical 3.11 packaging plus deploy smoke/monitor | Make the Cloudflare builder the sole deploy path and add remote-path tests/budgets | Canonical deployment-owner decision and immutable retention | Cloudflare bundle manifest/deploy artifact only |
| I-08 | Fixed by P0-1 | `web/js/services/history.js`, `history-snapshot.js`, and `app/history-inputs.js` | Intent/command/checkpoint, failure-atomicity, file-retention, no-op, composition, and legend tests | Browser Node plus Playwright functional/performance | No remaining P0 invariant gap; retain compact-intent and explicit-checkpoint ownership guards | Node structural plus Playwright performance | Implemented here | Stable artifact/session authority and performance instrumentation | None |
| I-09 | Production fix present; CI guard fixed by Phase A | `web/js/services/standalone-interactivity-assets.js`; `gbdraw/render/interactive_svg.py` | 25,000-feature Playwright test with apply/clear mutation cap, exact two-mutation navigation, and elapsed ceilings | Dedicated retry-free Playwright performance lane; unreachable at baseline | No remaining hot-path invariant gap; retain required-lane inventory | Playwright performance | P0-3 lane activation and ownership inventory | Single-worker/retry-free config | None; fixture SVG is temporary |
| I-10 | Partially fixed | `session-file.js`, `session-resources.js`, `result-normalization.js`, `session-authority.js`, and Gallery refresh tooling | Compression/dedupe/result normalization and Vibrio size/component gates | Browser Node and Gallery/Core coverage; large Gallery checks still repeat on supported main versions | Reader-cap tests and generalized expanded/component/copy/transfer-byte budgets for every session | Node plus canonical Gallery 3.11; large trend nightly/release | P1-3 schema inspector and size/copy report | Canonical session writer; P0-1 for history-copy accounting | Guard only: none; a schema change must canonically regenerate affected Gallery sessions |
| I-11 | Fixed with strong guards | `web/js/app/losat-cache.js`, `app/run-analysis.js`, `analysis/protein_colinearity.py`, Gallery refresh validator | Compact-ID grammar, cache validation, Vibrio exact-pair/size gates, migration Playwright spec | Node/Gallery plus dedicated LOSAT browser acceptance | Row-count multiplied by string-width budget for each repeated metadata column | Node/session hard gate, Gallery budget, LOSAT acceptance | Extend P1-3 with per-column and copy budgets | Stable runtime-handle/manifest schema and I-12 | None for guard; schema work regenerates affected Gallery sessions canonically |
| I-12 | Python path fixed; Web path still partially fixed | `gbdraw/session_io.py`; `web/js/app/losat-cache.js` and `services/config.js` | Python materialization call-count=1; Web validation correctness tests | Python Core and Browser Node lanes | One prepared Web manifest context reused by raw and derived entries, with JS call count=1 | Node hard gate, retaining Python Core test | Prepare manifest validation once in `losat-cache.js` and consume it from config | I-11 manifest schema | None |
| I-13 | Fixed by P0-2 | `services/svg-result-ingestion.js`, `app/candidate-render.js`, `state.js`, and `services/svg-sanitization.js` | Sanitization security/profile, runtime trust, and Worker/reflow/session admission tests | Browser Node plus functional/performance | No remaining P0 invariant gap; persisted data cannot forge the module-private admission marker | Node security/structural plus Playwright functional/performance | Implemented here | P0-5 ingress inventory; stable artifact identity/history sequence | None; security profile unchanged |
| I-14 | Matching-session regression fixed; broader P1-1 remains | `web/js/app/pyodide.js`, session/legend restore in `app-setup.js`, `index.html`, and export service | Pyodide startup and no-op guard tests | Browser Node plus Playwright performance | Lazy export and production Vue/Tailwind asset work remain outside P0 | Node structural, Playwright performance, packaging/offline | Retain the no-op guard; execute remaining P1-1 packages separately | P0-5 allowlist; P0-2/session stability | None for this correction |
| I-15 | P0 duplicate serialization/remount fixed; combined color traversal remains P1 | `web/js/app/svg-styles.js`, `watchers.js`, feature color actions, `feature-dom.js`, and `preview-runtime.js` | One dirty flush, zero edit remount, and one feature-index/handler build per root | Browser Node plus Playwright performance | A single compiled palette/specific traversal remains a later optimization | Node counter test plus Playwright performance | P1-2 may combine palette/specific resolution after profiling | P0-2 commit boundary is complete | None |
| I-16 | Partially fixed | Python composition planner/render metadata and Circular/Linear assemblers; Web composition actions | Layout/render/assembler contracts, runtime parity, and quick real-assembler benchmark | Core and Browser Node; Playwright composition now functional-lane reachable | No production SVG serialize/reparse for bounds; no second current-schema policy/default owner | Core/Node architecture; standard main/nightly benchmark | Add owner/source contracts and promote benchmark reporting | P2-2; P1-5 lane ownership | Benchmark JSON only; no tracked SVG change expected |
| I-17 | Partially fixed; expensive browser duplication removed by Phase A | `test.yml`, `tests/test_web_packaging.py`, and `tools/verify_gui_offline.py` | CI selectors plus four independently named offline GUI contracts | Browser behavior once on canonical 3.11; 3.10/3.12 retain non-slow version-sensitive bridge coverage | Repository-wide test cost/contract/marker/runner inventory and feedback-time budget | Canonical 3.11 browser plus version-sensitive Core/bridge matrix | P1-5 owner inventory and duration budget | P0-4 provides the split | None |
| I-18 | Fixed by Phase A | `test.yml`, Playwright configs/scripts, and 12 Playwright specs | 86-test discovery; performance specs include structural mutation guards | 77 functional tests, 7 performance tests, and 2 dedicated LOSAT migration tests each have an always-reporting lane | Branch-protection required-check configuration is not observable without repository authentication | Separate functional and retry-free performance jobs | P0-3 implemented here; later add an inventory meta-contract if desired | None | Failure-only trace artifacts; no tracked output |
| I-19 | Fixed by Phase A | Four contracts in `verify_gui_offline.py`/`test_web_packaging.py`; canonical Browser job | Offline init/network, palette, exports, and Linear cache contracts | Slow browser contracts once on Python 3.11; excluded from three-version slow matrix | No mega-test remains; future additions must keep one subsystem per contract and semantic settlement diagnostics | Canonical Browser 3.11 | P0-4 implemented here | Prepared browser wheel and current comparison-plan API | None |
| I-20 | Still present | Repository review/process surface (`.github` has no performance template or validator) | None | No | Required N/complexity/thread/copy/cancellation/before-after/guard fields for performance-sensitive paths | Lightweight PR metadata check | P2-1 template, sensitive-path manifest, pure validator, and tests | P1-5 owner inventory; P1-4/P2-2 can supply evidence | None |

## Recommendation reconciliation

| Recommendation | Status | Current owner | Existing regression test | CI | Missing structural invariant | Target CI lane | Proposed PR-sized work package | Dependencies | Expected generated-artifact impact |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P0-1 | Fixed | One History stack with compact intent/command entries and explicit artifact checkpoints | History/input/legend/composition tests | Browser Node plus Playwright functional/performance | None for P0; full checkpoints remain only at explicit artifact boundaries | Node plus Playwright performance | Implemented here | Stable session/artifact identity preserved | None |
| P0-2 | Fixed | `services/svg-result-ingestion.js` owns admission; PreviewRuntime owns mounted edits/flush | Ingestion/security/preview tests | Browser Node plus functional/performance | None for P0; one minimum sanitize/parse/serialize/mount remains | Node security/structural plus functional/performance | Implemented here | P0-5 owner inventory retained | None |
| P0-3 | Fixed by Phase A | Functional/performance configs, package scripts, and workflow jobs | 77 functional plus 7 performance tests; 25k structural assertions unchanged | Separate always-reporting jobs | Authenticated branch-protection setting remains externally unverifiable | Playwright functional and performance | Implemented here | Dedicated LOSAT runner retains its two tests | Failure-only trace artifacts |
| P0-4 | Fixed by Phase A | Four offline GUI contracts and workflow selectors | Four named fresh-context browser tests | Canonical Browser 3.11 only; slow matrix excludes browser | Keep future contracts independently settled and diagnosed | Canonical Browser 3.11 | Implemented here | Prepared browser wheel | None |
| P0-5 | Fixed by Phase A | Existing diagram-generation service/Worker plus architecture contract | Static import/owner/canonical-dispatch tests and existing runtime behavior tests | Browser Node glob | Allowlist updates must accompany any intentional new owner | Fast Node architecture lane | Implemented here | None; wrapper remains shared but callable only from Worker | None |
| P1-1 | Partially fixed; not implemented here | Main Pyodide manager, restore/startup, export and page assets | Palette init=0 unit test | Browser Node | No-op restore init=0, lazy export, production asset and cold-start contracts | Node, Playwright performance, packaging | Three small packages: restore, export, production assets | P0-5 and packaging/CSP inventory | Only production-asset package changes vendored inventory |
| P1-2 | Partially fixed; not implemented here | Standalone/preview search, SVG editor index, color-rule/style owners | 25k search guard, preview dirty-flush and feature action tests | Node plus new performance lane | All bulk work proportional to changed items; no same-root index rebuild | Node counters plus performance | Repair index lifecycle and unify bulk color pass | P0-2 for single commit | None |
| P1-3 | Partially fixed; not implemented here | Session writer/reader/resources, LOSAT cache, Gallery validator | Session dedupe, compact IDs, Vibrio hard limits | Node/Gallery/LOSAT lanes | General expanded/component/column/copy/transfer-byte report | Canonical Gallery plus nightly/release trend | Schema/budget inspector and size-diff artifact | P0-1 copy accounting; stable schema | CI JSON; schema work may regenerate Gallery sessions |
| P1-4 | Still present; not implemented here | No single production owner yet | Existing elapsed and structural tests do not provide observer/base-head evidence | No dedicated observer lane | Long-task/LoAF marks, stage counters, same-runner median relative and absolute gates | Playwright performance | Test-injected observer and JSON/trace reporter | P0-2 and P1-2 instrumentation points | CI JSON/traces only |
| P1-5 | Partially fixed by P0-4; not implemented in full | CI workflow and test ownership surface | Workflow selectors | Browser duplication is removed; no duration budget | Contract/cost/marker/runner inventory and PR feedback budget | Lightweight CI meta-check | Add owner inventory and duration regression check | P0-4 split | None |
| P2-1 | Still present; not implemented here | Repository process surface | None | No | Empty performance-evidence fields must fail for sensitive paths | Lightweight required PR metadata job | Template, path manifest, validator, unit tests | P1-5; optional P1-4/P2-2 artifacts | None |
| P2-2 | Partially present; not implemented here | `tools/benchmark_diagram_layout.py` and benchmark tests | Quick and standard profiles exist | Quick contracts in Core; standard base/head is not CI | Stable 3-warmup/21-sample main/nightly comparison and stored median/MAD | Main/nightly benchmark; quick PR smoke | Promote existing runner and persist JSON | CI owner/budget decisions | Benchmark JSON only |

## Phase A evidence

### P0-3: required Playwright execution

- `playwright.functional.config.js` selected 76 tests during Phase A. The final
  combined tree selects 77 after adding the P0-2 frozen-v39 composition
  compatibility regression. It excludes the seven current performance tests
  and the two LOSAT cache migration tests, whose existing dedicated runner
  remains authoritative.
- `playwright.perf.config.js` selects the five Phase A tests in
  `interactive-svg-search-performance.playwright.spec.js` plus the two P0
  ownership/security tests in `webapp-performance.playwright.spec.js`, and fixes
  `retries=0`, `workers=1`, `fullyParallel=false`, and
  `trace=retain-on-failure`.
- Phase A validation ran 76 functional tests and five performance tests. Final
  combined-tree validation is recorded below. The 25,000-feature mutation
  limits and exact two-mutation navigation assertion were not edited.
- Both jobs have unconditional workflow reachability for the configured push,
  pull-request, and manual triggers; no changed-file filter can omit their
  checks. Failure-only trace directories are uploaded for seven days.
- Exact branch-protection required-check configuration could not be inspected:
  the local `gh` client is unauthenticated and the public endpoint returns 401.

### P0-4: narrow offline GUI contracts

- Pytest collection selects exactly four `slow and browser` tests:
  offline initialization/network isolation, palette preview, three-format
  export, and Linear LOSAT generation/cache readiness.
- Each test dispatches only one contract in a fresh browser context. The
  backwards-compatible `smoke-test` CLI aggregates the same four fresh runs.
- The Linear contract uses `setLinearComparisonGlobalAction('losat')`, asserts
  adjacent/active LOSAT intent before and after the run, and waits for the
  current run sentinel, current Result, or explicit error plus cache growth.
  A timeout reports run, comparison, cache, and runtime state. The retired
  `blastSource` assignment and 180-second blind wait were removed.
- Each focused node passed. The final exact canonical CI selector passed all
  four in 41.62 seconds; the slowest was Linear/cache at 14.17 seconds.
- `tests/test_web_packaging.py -m "browser and not slow"`: 3 passed, confirming
  version-sensitive non-slow packaging/browser coverage remains separate.
- The Browser Python 3.11 job now owns the four slow browser contracts. The
  Python 3.10/3.11/3.12 slow matrix selects `slow and not browser` and no longer
  installs Chromium.

### P0-5: architecture contracts

- `tests/web/architecture-contracts.test.mjs`: 5 subtests passed. The test
  traverses literal first-party imports from `web/js/app.js`, asserts the main
  graph reaches no `workers/` module, and allowlists Worker constructors,
  diagram client importers, `loadPyodide`, canonical request builder callers,
  the sole render dispatch payload, and the embedded Python wrapper lookup.
- The existing shared `app/python-helpers.js` still installs the wrapper in the
  main Pyodide runtime, but the new contract proves that only the diagram Worker
  retrieves it. Splitting that payload belongs to P1-1, not this guard-only PR.
- Phase A validation passed 64/64 Node test files, including the new
  architecture file and existing runtime/canonical behavior coverage. The
  final combined tree passes 65/65 below.

### Cross-checks and artifacts

- The exact workflow YAML parses successfully and contains 10 jobs, including
  distinct `playwright-functional` and `playwright-performance` jobs.
- `python -m py_compile tools/verify_gui_offline.py tests/test_web_packaging.py`:
  passed.
- `python tools/benchmark_diagram_layout.py run --quick`: passed; medians were
  74.35 ms multi-Circular, 7.06 ms multi-Linear, and 94.99 ms single-Circular.
- `ruff check gbdraw/`: passed. A broader direct check of the touched helper
  reports its two baseline violations (`E402` after deliberate `sys.path`
  setup and the existing HTTP-handler `E731`); running Ruff against the file at
  `HEAD` reproduces both, so they are not introduced or repaired in this scope.
- `git diff --check`: passed.
- No tracked generated artifact changed. The ignored browser wheel was prepared
  with `python tools/prepare_browser_wheel.py` solely for offline verification.

## P0-1/P0-2 execution evidence

### Worktree boundary and confirmed root causes

- P0 work started from baseline commit
  `3bca0e8d98d7743e5df4f30dfb33a0fa34ddbf2f` with the validated Phase A
  workflow, split Playwright configs, offline-contract split, packaging tests,
  architecture test, and this ledger already modified or untracked. Those
  changes were preserved and reviewed separately from production changes.
- The obsolete ordinary History path called `buildHistorySnapshot()` at both
  input begin and commit. That cloned Results, feature catalogs, editor state,
  sequence-derived state, and run artifacts; signed the complete object; then
  limit enforcement signed and sized retained snapshots again. Commit
  `e3f6d9ec` introduced this full-snapshot topology.
- SVG Results crossed two current admission paths: candidate commit sanitized,
  parsed, transformed, and serialized, then the `state.svgContent` computed
  value sanitized the same full string again before `v-html`. Live editor
  persistence rewrote `results[].content`, which made that computed value an
  edit event bus and caused root replacement, watcher work, and index/listener
  rebuilding. The duplicate topology traces to `7aab08f8`.

Verified previous call graphs were:

```text
ordinary setting edit
  -> history.begin/commit
  -> buildHistorySnapshot twice
  -> serialize Results / clone generated domains
  -> JSON signature
  -> retained-stack re-sign/re-size during limit enforcement

Worker Result
  -> optional pairwise-legend parse/serialize
  -> candidate sanitize/parse/override/serialize
  -> state computed sanitize
  -> v-html mount -> watcher/index/listener setup

session/imported Result
  -> applyResultsData(raw string)
  -> state computed sanitize
  -> v-html mount
  -> live layout/legend persistence -> Result rewrite -> sanitize/remount

feature color/stroke/visibility edit
  -> mutate live SVG -> serialize -> Result rewrite
  -> computed sanitize -> root replacement -> watcher/index/listener setup
```

### P0-1 ownership now

- `history.js` retains one Undo/Redo stack with three entry kinds: compact
  intent patches for ordinary edits, existing semantic commands for targeted
  editor actions, and artifact checkpoints for explicit multi-domain
  boundaries. There is no `buildHistorySnapshot`/`applyHistorySnapshot`
  compatibility fallback.
- `history-snapshot.js` builds configuration, file descriptors, inactive-mode
  settings, editor overrides, compact composition deltas, and other user intent
  without reading Results, feature catalogs, sequence-source artifacts, LOSAT
  caches/manifests, orthogroup result sets, or run state.
- Generate and confirmed multi-domain Reset use `runUndoableCheckpoint`.
  Session load establishes an explicit baseline. Generation cancellation and
  stale rollback retain stable runtime Result identity outside History.
- Entry signatures and byte sizes are computed once when an entry is created;
  retained-entry bytes are updated incrementally. Command size callbacks run
  once. Undo/Redo moves an entry only after its apply/revert succeeds.
- Manual configuration, color, visibility, stroke, label, legend, and
  composition edits remain in the same History stack. Intent restoration
  reconciles the existing mounted root and marks it dirty; it does not restore
  or remount a copied SVG.

Structural evidence for ten ordinary edits after the former large-session restore:

```text
artifact checkpoint builds added = 0
SVG bytes copied/serialized for History = 0
Result serialization for History = 0
previous entries re-signed = 0
previous entries re-sized = 0
no-op entries added = 0
artifact runtime identity changes across config Undo/Redo = 0
```

The browser structural probe and History unit test also cover failure atomicity,
30-action and 200 MiB limits, one-time command sizing, file release, explicit
Generate checkpoint Undo/Redo, and omission of every generated domain from an
ordinary intent.

### P0-2 ownership now

- `services/svg-result-ingestion.js` is the single Result admission boundary.
  Worker, reflow, supported session, and imported SVG strings use the existing
  DOMPurify profile, one detached parse, commit-time editor/legend transforms,
  one canonical serialization, and a module-private runtime marker. JSON and
  session data cannot forge that marker.
- `state.svgContent` displays only admitted Results and does no sanitization.
  Pairwise-legend suppression and candidate editor overrides are part of the
  detached admission transaction; the old standalone pre-parse was removed.
- PreviewRuntime owns the mounted Result identity. Live editor changes mark the
  root dirty; persistence may update Result text once without replacing the
  current root. Result switch, session save, and explicit boundaries flush the
  root. `feature-dom.js` owns the one WeakMap feature index shared by palette,
  handlers, and PreviewRuntime.
- A matching current-session legend is compared from the mounted sanitized DOM
  before any main-thread Python helper is requested. The no-op path initializes
  Pyodide zero times.

Observed large-session admission and ordinary edit counters:

```text
sanitize / parse / serialize / Result commit = 1 / 1 / 1 / 1
SVG-root mount / feature-index build / handler setup = 1 / 1 / 1
displaying committed Result: additional sanitizes = 0
ordinary feature visibility edit: parse 0, remount 0, index 0, handler setup 0
no-op action: parse 0, serialize 0, Result update 0, remount 0
matching-session main-thread Pyodide initialization = 0
```

The malicious browser fixture retained the existing removal of scripts,
`foreignObject`, event attributes, and `javascript:` URLs while preserving safe
interactive feature metadata.

### Browser timing evidence

- Audited before values for the same former large-session fixture/runtime were approximately
  7,000 ms restore and 2,890 ms maximum long task; History begin was roughly
  99--103 ms and commit 143--190 ms.
- Five already-collected controlled post-change large-session restore samples measured
  1,813, 1,823, 1,834, 1,849, and 1,863 ms (median 1,834 ms), a 73.8% improvement
  from the audited approximately 7,000 ms baseline. The final post-fix probe
  measured a 491 ms maximum long task, an 83.0% improvement from approximately
  2,890 ms. Its ten ordinary edits measured begin p95 2.4 ms, commit p95
  5.9 ms, early median 6.4 ms, and later median 6.3 ms. All absolute and
  non-growth targets passed. Per user direction, no additional three-repeat
  performance run was added; the final seven-test performance suite was run
  once after the implementation stabilized.

### Residual main-thread cost classification

- **A, removed:** eager main-thread Pyodide initialization for an unchanged
  matching session legend (`7aab08f8`, `app/app-setup.js`). Genuine explicitly
  requested Python-backed helper actions remain **B** (`63970c94`,
  `app/pyodide.js`).
- **A, removed:** duplicate committed-Result sanitize/parse/display work from
  `7aab08f8`. The surviving security sanitize, detached transform parse, one
  serialization, and browser mount are **C**, the minimum under the current
  secure DOM/SVG model.
- **C:** browser style/layout and one `v-html` mount for the admitted large SVG;
  this rendering model has existed since `1ff0c843` (`web/index.html`). Extra
  application post-mount reconciliation remains **B** (`63970c94`, `2554176a`).
- **B:** remaining initial feature/comparison metadata scans per new root
  (`2a0cf199`, `2554176a`, `3e88d5a4`). Same-root duplicate indexing is removed.
- **B:** Vue and Tailwind startup existed in the initial Web app (`1ff0c843`);
  `c5753bde` vendored the same assets rather than introducing the cost.

### P0 validation and artifact impact

- `node --test tests/web/*.test.mjs`: 65/65 test files passed in 6.46 s.
- `npm run test:web:functional-smoke`: 77/77 passed in 54.1 s. This includes
  History/session/editor/export behavior and the real frozen-v39 admission and
  first-authorized-layout-edit regression.
- `npm run test:web:perf-smoke`: 7/7 passed in 25.5 s with one worker and zero
  retries. The large-session probe reported one sanitize, parse, serialization, Result
  update, mount, feature-index build, and feature-handler setup; main-thread
  Pyodide initialization was zero. The malicious SVG security case also passed.
- `npx playwright test tests/web/losat-cache-migration.playwright.spec.js
  --project=chromium --workers=1`: 2/2 passed in 1.1 minutes.
- Focused Pairwise SVG admission and failed-generation rollback Playwright:
  2/2 passed in 15.7 s.
- `python tools/verify_gui_offline.py smoke-test --contract all`: passed all
  four fresh-context contracts. The exact CI selector
  `python -m pytest tests/test_web_packaging.py -m "slow and browser"
  --durations=30` also passed 4/4 in 44.38 s.
- `python -m pytest tests/test_web_packaging.py -m "browser and not slow" -q`:
  3 passed, 38 deselected in 7.23 s.
- Workflow YAML parsed with 10 jobs. Both functional and performance jobs
  prepare the generated browser wheel before Worker-ready browser tests.
- `ruff check gbdraw/` and `python -m py_compile` passed.
- `python tools/benchmark_diagram_layout.py run --quick`: passed (85.90 ms
  multi-Circular, 7.75 ms multi-Linear, 98.40 ms single-Circular).
- No SVG reference, Gallery session, screenshot, public showcase, wheel, or
  other tracked generated artifact changed. The prepared browser wheel remains
  ignored and untracked.

### CI follow-up

- The first pushed combined-tree run exposed one stale Python architecture
  assertion. `tests/test_web_mode_profiles.py` still required pairwise-legend
  DOM ownership in `run-analysis.js`, although P0-2 deliberately moved the
  shared selector and removal into `candidate-render.js` inside the single
  detached SVG-ingestion transaction. The test now checks the suppression
  intent at `run-analysis.js` and the selector/removal at its canonical owner;
  no production path was restored. The exact test file passes 2/2.
- The first CI performance run passed all 25,000-feature structural mutation
  gates but reported 1,254 ms against the 1,000 ms apply ceiling. The old
  Node-side `Date.now()` interval included Playwright locator resolution,
  browser IPC, assertion polling, and a second browser round trip. It did not
  isolate application work. The test now starts `performance.now()` in the
  browser immediately before the native click and ends after the same two-frame
  application settlement. The 1,000 ms apply and 250 ms navigation ceilings,
  match assertions, mutation caps, and exact two-mutation navigation gate are
  unchanged. A CI-mode focused run measured 153.9 ms apply and 6.7 ms
  navigation.
- `CI=1 TMPDIR=/tmp npm run test:web:perf-smoke`: 7/7 passed in 25.5 s with
  one worker and zero retries. Its large-session probe reported 1,980 ms restore, 523 ms
  maximum long task, begin p95 2.4 ms, commit p95 6.3 ms, one of each ingestion
  ownership operation, and zero main-thread Pyodide initialization.
