# Codex Task: gbdraw P0 Web Performance Regression Remediation

Repository: `satoshikawato/gbdraw`

## Objective

Remove the measured root causes of the current user-visible Web GUI performance regression.

This is a production-quality root-cause remediation. It is **not** a workaround, timeout adjustment, feature reduction, cache added over an incorrect ownership model, or temporary compatibility path.

Implement only:

- **P0-1 — History ownership remediation:** ordinary edits must store compact user intent instead of rebuilding, cloning, signing, and sizing generated artifacts.
- **P0-2 — SVG ingestion ownership remediation:** each untrusted SVG Result must cross one sanitization/parse/commit boundary and must not be fully sanitized or remounted again merely to display or persist an ordinary edit.

Preserve the current Phase A guardrail work in the worktree. Do not redo P0-3/P0-4/P0-5 except for a narrowly required correction.

The primary goal is to make the current application usable again and prepare validated commits suitable for fast-forwarding to `main`.

## Read first

Read in this order:

1. `AGENTS.md`
2. `CLAUDE.md`
3. `gbdraw/web/CLAUDE.md`, if present
4. `docs/internal/GBDRAW_PERFORMANCE_REGRESSION_AUDIT_2026-08-11.html`
5. `docs/internal/PERFORMANCE_REGRESSION_REMEDIATION_EXECUTION_PLAN_2026-08-11.md`
6. The current internal plans for:
   - Web app performance remediation
   - Web preview edit performance
   - session compatibility
   - SVG sanitization and Result normalization
   - Worker/Pyodide ownership
   - canonical request ownership

Current HEAD and the current worktree are authoritative. Verify the current call paths before changing them.

## Baseline and worktree safety

Before modifying files, run:

```bash
git status --short --untracked-files=all
git rev-parse HEAD
git diff --stat
```

Record which changes belong to the completed/in-progress Phase A work and which are unrelated.

Do not:

- reset, discard, overwrite, or broadly reformat existing changes;
- introduce a second rendering, request, History, state, or compatibility path;
- manually edit generated SVGs, sessions, Gallery assets, screenshots, or reference outputs;
- weaken thresholds, reduce functionality, reduce History scope, add retries, or increase timeouts as the primary remedy;
- commit unrelated files;
- push or open a PR during implementation.

## Preserve these contracts

- no-build, same-origin, offline-capable Web application;
- browser-local biological data processing;
- diagram generation and generation-time feature extraction remain module-Worker owned;
- current canonical typed request/session-request ownership;
- current supported session readers and writers;
- current SVG sanitization security profile;
- public Python and CLI behavior;
- one logical diagram per Result;
- one resource payload copy per current session;
- current Gallery/session compatibility;
- all current feature, label, legend, layout, palette, History, Undo/Redo, export, and session behavior.

## Establish a reproducible baseline

Use fixed representative fixtures and the same browser/runtime for before/after measurements. Include at least:

1. loading a representative large Gallery session;
2. one ordinary number/select/text setting edit after loading it;
3. ten repeated ordinary setting edits;
4. Generate-to-preview Result commit;
5. one feature color or stroke edit;
6. one no-op edit;
7. a malicious SVG sanitization fixture.

Instrument, preferably through test-only hooks:

- History snapshot builds;
- artifact checkpoint builds;
- History signature computations;
- History byte-estimate computations;
- copied or serialized SVG bytes used by History;
- full-string sanitizer calls by Result identity;
- full SVG parses;
- full SVG serializations;
- Result-content updates;
- SVG-root mounts;
- feature-handler/index builds;
- main-thread Pyodide initialization count;
- long tasks and elapsed time where supported.

Do not add permanent telemetry or production console spam.

Record the verified current call graphs for:

```text
ordinary setting edit → History begin/commit
Worker Result → preview
session Result → preview
imported SVG Result → preview
feature color/stroke edit → persistence
```

## P0-1 — History ownership remediation

Inspect at minimum:

- `gbdraw/web/js/services/history.js`
- `gbdraw/web/js/services/history-snapshot.js`
- `gbdraw/web/js/services/history-inputs.js`
- `gbdraw/web/js/services/history-files.js`
- existing users of command History
- configuration/session serialization helpers
- Generate, Reset, session-load, Result-replacement, and Result-switch boundaries

### Required architecture

1. Keep one logical Undo/Redo stack.
2. Reuse and complete the existing command-entry mechanism. Do not create another History framework.
3. Ordinary form/config/draft/UI edits must store compact before/after values or commands.
4. An ordinary edit must not clone, serialize, sign, or size:
   - SVG Results;
   - feature catalogs;
   - sequence sources;
   - LOSAT caches;
   - protein identity manifests;
   - embedded resources;
   - other generated artifacts.
5. Full artifact checkpoints are allowed only at explicit artifact-changing boundaries, such as:
   - Generate;
   - session load;
   - multi-domain Reset;
   - Result replacement;
   - Result switch only where an artifact flush is genuinely required.
6. Stable generated artifacts must be referenced by stable runtime identity rather than copied into every ordinary History entry.
7. Compute an entry's signature and byte estimate once when the entry is created and retain those values on the entry.
8. Maintain incremental total-byte accounting on push, pop, redo clearing, and eviction.
9. Do not repeatedly `JSON.stringify()` all retained snapshots while enforcing History limits.
10. Preserve:
    - existing Undo/Redo order;
    - 30-action behavior;
    - 200 MiB protection;
    - file-store retain/release semantics;
    - failure atomicity;
    - inactive-mode state;
    - editor overrides;
    - session-save correctness.
11. Keep async generation cancellation, stale-result rejection, and rollback ownership separate from History storage.
12. Remove or make unreachable the superseded ordinary-edit full-artifact path. Do not leave it as a silent fallback.

### Hard structural gates

An ordinary setting edit must produce:

```text
artifact snapshot builds = 0
SVG bytes cloned/serialized for History = 0
full Result serialization for History = 0
```

Ten ordinary edits must satisfy:

```text
previous entries re-signed = 0
previous entries re-sized = 0
work does not grow with History depth
```

A converted command edit must satisfy:

```text
buildHistorySnapshot() calls = 0
applyHistorySnapshot() calls = 0
```

Also verify:

- a no-op edit creates no History entry;
- Undo failure leaves stack position unchanged;
- Generate → config edit → Undo/Redo preserves artifact identity and configuration;
- session load → editor command → Undo/Redo → session save remains internally consistent;
- file resources are released only when no retained entry references them.

## P0-2 — SVG ingestion ownership remediation

Inspect at minimum:

- `gbdraw/web/js/state.js`
- `gbdraw/web/js/app/candidate-render.js`
- `gbdraw/web/js/services/svg-sanitization.js`
- `gbdraw/web/js/services/svg-serialization.js`
- Result normalization and Result commit modules
- Worker-Result ingestion
- session Result restoration
- imported SVG handling
- preview mounting and watchers
- color, visibility, stroke, legend, and layout edit paths

### Required architecture

1. Define one explicit SVG/Result ingestion boundary for:
   - diagram Worker Results;
   - supported session Results;
   - imported SVG Results.
2. Every untrusted SVG must pass through the existing shared DOMPurify profile before it can enter the live DOM.
3. Do not weaken, bypass, duplicate, or fork the sanitizer policy.
4. A committed in-memory Result may carry runtime-only trusted/sanitized state, but persisted or imported data must never be able to spoof that state.
5. Do not full-string sanitize the same committed Result again merely to display it.
6. Apply commit-time editor overrides in one non-live parse/commit transaction.
7. For one Result identity, ingress-to-display may perform at most:
   - one full-string sanitize;
   - one full parse;
   - one final Result commit;
   - one SVG-root mount.
8. One user action may perform at most:
   - one complete SVG serialization;
   - one Result-content update.
9. A no-op action performs zero parse, serialization, Result update, and remount work.
10. Build feature handlers and DOM indexes once per new SVG root, not whenever an equivalent Result string is rewritten.
11. Small fill, visibility, and stroke edits must not remount the SVG root merely to persist the edit.
12. Keep live-DOM mutation and persisted Result text ownership explicit. Do not use `results[].content` as a broad edit event bus.
13. Remove or make unreachable the duplicate committed-Result sanitization/remount path. Do not retain it as a fallback for current Results.
14. Preserve:
    - XSS and DOM-clobbering defenses;
    - supported session compatibility;
    - offline operation;
    - export behavior;
    - interactive feature/match metadata;
    - composition metadata;
    - editor overrides;
    - one Result per logical diagram.

### Hard structural gates

For Worker, session, and imported Results separately:

```text
full-string sanitize count per Result identity <= 1
full SVG parse count per Result identity <= 1
SVG-root mount count per Result identity = 1
```

Displaying an already committed Result:

```text
additional full-string sanitize count = 0
```

One color, specific-rule, visibility, or stroke action:

```text
full SVG serialization count <= 1
Result-content update count <= 1
SVG-root remount count = 0
```

No-op action:

```text
parse = 0
serialize = 0
Result update = 0
remount = 0
```

Feature handlers/indexes:

```text
build count = 1 per SVG root
```

All malicious SVG security fixtures must remain rejected or sanitized with the same effective security contract.

## Classify remaining main-thread costs

After P0-1 and P0-2 structural gates pass, profile the same large-session and ordinary-edit flows.

For each remaining cost, classify it with evidence as:

- **A — regression/reintroduced regression:** previously absent or previously removed, then reintroduced;
- **B — longstanding architectural cost:** present before the recent regression;
- **C — unavoidable minimum cost:** required by the current DOM/SVG/browser model;
- **D — unresolved:** insufficient evidence.

Classify at least:

1. main-thread Pyodide initialization;
2. DOM style/layout work for large SVGs;
3. remaining full-feature scans;
4. Tailwind/Vue runtime startup cost;
5. the minimum necessary large-SVG parse/mount cost.

Use current code plus `git log`, `git blame`, `git log -S`, or `git log -G` where useful. Record supporting commits/files in the execution ledger.

This classification is evidence work, not permission for a broad rewrite.

### Narrowly permitted additional correction

If an unchanged/no-op session restore still initializes main-thread Pyodide, and the initialization is unnecessary for the restored palette/legend state, remove that unnecessary initialization in this task.

Required invariant:

```text
unchanged session restore with already-matching saved palette/legend
main-thread Pyodide initializations = 0
```

Do not move diagram rendering or generation-time feature extraction back to the main thread.

Do not expand this task into Tailwind replacement, Vue production bundling, DOM virtualization, Canvas/WebGL rendering, or a general front-end rewrite. Record those as separate work only if profiling shows they remain important after P0-1/P0-2.

## Performance evidence and acceptance targets

Structural counters are hard gates. Also report before/after browser measurements from the same fixture and environment.

Targets:

```text
ordinary large-session setting edit:
  History begin p95 < 50 ms
  History commit p95 < 50 ms

ten repeated ordinary edits:
  later median <= early median * 1.20

large-session restore:
  median improves by >= 20%
  maximum main-thread long task improves by >= 40%

unchanged matching session restore:
  main-thread Pyodide initialization count = 0
```

Do not relax a target to declare success.

If a target is missed:

1. profile the remaining largest task;
2. determine whether it belongs to P0-1, P0-2, or the narrowly permitted no-op Pyodide path;
3. fix it only if it belongs to this scope;
4. otherwise report it explicitly as the next blocker and classify it A/B/C/D.

Completion requires both:

- removal of the measured obsolete hot paths; and
- passing structural ownership gates.

Adding tests without removing the hot paths is not completion.

## Validation

Run focused tests first, then the actual Phase A CI commands. Adapt names to the current repository rather than creating duplicate configurations.

At minimum consider:

```bash
node --test tests/web/*.test.mjs
npx playwright test
npx playwright test --config <performance-config>
python -m pytest <focused History/session/Web packaging/browser tests>
python tools/benchmark_diagram_layout.py run --quick
ruff check gbdraw/
git diff --check
```

Run the split offline checks relevant to initialization, palette, export, and Linear LOSAT. Do not hide failures with larger timeouts.

Verify at the end:

```bash
git status --short --untracked-files=all
git diff --stat
git diff --check
```

No broad generated-artifact refresh is expected. If the canonical wire/output format genuinely changes, use the canonical generator and list every generated file and command.

## Update the execution ledger

Update only:

`docs/internal/PERFORMANCE_REGRESSION_REMEDIATION_EXECUTION_PLAN_2026-08-11.md`

Record:

- P0-1 and P0-2 status;
- confirmed root causes;
- removed paths;
- structural counter evidence;
- before/after measurements;
- residual A/B/C/D classification;
- tests and CI lanes;
- unmet targets and remaining risk.

Do not create another overlapping implementation plan.

## Commit and main-branch handoff

After all required gates pass, create two reviewable commits where practical:

1. `fix(web): make ordinary history edits artifact-independent`
2. `fix(web): unify sanitized SVG result ingestion`

A narrowly permitted no-op main-thread Pyodide correction may be included in the second commit or a third focused commit.

Do not push automatically.

At completion report:

1. verified previous call graphs;
2. new ownership boundaries;
3. root causes removed;
4. files changed;
5. structural counters before/after;
6. browser measurements before/after;
7. residual main-thread costs classified A/B/C/D;
8. tests and exact outcomes;
9. pre-existing failures;
10. generated artifacts changed, if any;
11. exact commit hashes;
12. clean/dirty worktree status;
13. the exact commands to fast-forward the validated branch into `main` and push it.

A successful result must be ready for an informed, immediate main-branch push without relying on unverified claims.
