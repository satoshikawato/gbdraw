# Codex Task: gbdraw Performance Regression Guardrails — Phase A

Repository: `satoshikawato/gbdraw`

## Objective

Reconcile the performance-regression audit against the current HEAD, then implement only these executable guardrails:

- **P0-3:** Run existing Playwright functional and performance specs in required CI.
- **P0-4:** Split the offline GUI mega-smoke into narrow, diagnosable browser contracts and stop repeating the same expensive browser flow across Python 3.10/3.11/3.12.
- **P0-5:** Enforce Worker, main-thread Pyodide, and canonical-request ownership with lightweight architecture tests.

Do **not** implement P0-1 History redesign or P0-2 SVG-ingestion redesign in this task, except for minimal test instrumentation.

## Read first

Read in this order:

1. `AGENTS.md`
2. `CLAUDE.md`
3. `gbdraw/web/CLAUDE.md`, if present
4. `docs/internal/GBDRAW_PERFORMANCE_REGRESSION_AUDIT_2026-08-11.html`
5. Current internal plans concerning:
   - Web app performance remediation
   - Web preview edit performance
   - CI test reliability/runtime
   - Worker/Pyodide ownership
   - Session compatibility
   - Canonical request ownership

The audit baseline was `3bca0e8d98d7`, but current HEAD and worktree are authoritative. Verify every finding before acting.

## Baseline and safety

Before modifying files, run:

```bash
git status --short --untracked-files=all
git rev-parse HEAD
git diff --stat
```

Inspect the latest main-branch GitHub Actions run with `gh` if available. Record pre-existing failures separately.

Do not reset, discard, overwrite, broadly reformat, commit, push, or open a PR.

Preserve:

- no-build, same-origin, offline-capable Web app
- browser-local biological processing
- current session compatibility
- canonical typed request ownership
- SVG sanitization profile
- Worker-based diagram generation
- public Python/CLI behavior
- generated-artifact ownership
- one logical diagram per Result
- one resource payload copy per current session
- Gallery/session compatibility

Do not weaken thresholds, increase timeouts as the primary fix, add retries to hide failures, introduce alternate render/request paths, or manually edit generated artifacts.

## Required execution ledger

Create:

`docs/internal/PERFORMANCE_REGRESSION_REMEDIATION_EXECUTION_PLAN_2026-08-11.md`

Reconcile every audit item `I-01` through `I-20` and every P0/P1/P2 recommendation. For each item record:

- status: fixed / partially fixed / still present / obsolete
- current owner module
- existing regression test
- whether that test actually runs in CI
- missing structural invariant
- target CI lane
- proposed PR-sized work package
- dependencies
- expected generated-artifact impact

Do not copy the audit mechanically. Verify code, tests, workflows, and current CI.

## P0-3: Playwright must run in CI

Requirements:

- Keep Node unit tests and Playwright tests as separate commands/jobs.
- Add an explicit Playwright functional lane.
- Add a dedicated Playwright performance config with:
  - `retries: 0`
  - `workers: 1`
  - `fullyParallel: false`
  - trace retained on failure
- Ensure `tests/web/interactive-svg-search-performance.playwright.spec.js` is executed in CI.
- Preserve its mutation-count/structural assertions; do not rely only on elapsed time.
- Do not move it into an unreachable, optional, or perpetually skipped job.
- Required checks must complete even when changed-file filtering skips expensive work.

## P0-4: Split the offline GUI mega-smoke

Inspect the current equivalent of:

`tests/test_web_packaging.py::test_offline_gui_smoke_test_covers_palette_preview_behavior`

Verify whether its timeout still occurs later in Linear LOSAT readiness rather than palette preview.

Replace the mega-smoke with independently named contracts covering at least:

1. offline initialization and no external network
2. palette preview behavior
3. SVG/PNG/PDF export
4. Linear LOSAT generation and cache readiness

Requirements:

- Test names must identify the actual subsystem.
- Use semantic ready/error conditions, not arbitrary long sleeps.
- Run the expensive browser behavior flow once on canonical Python 3.11.
- Python 3.10 and 3.12 should retain only genuinely version-sensitive Python/core/bridge coverage.
- Do not solve this by raising the 180-second timeout, adding retries, or merely marking the old mega-test slow.
- Remove dead/redundant test paths after replacement.

## P0-5: Executable architecture contracts

Add lightweight static Node or Python tests that fail when:

- a second diagram-generation Worker client is introduced;
- Worker-only render modules are imported into the main-thread app path;
- `loadPyodide` or equivalent main-thread initialization appears outside an explicit allowlist;
- Web rendering bypasses the canonical typed request/session-request boundary;
- production code adds an alternate direct rendering path.

Tests must inspect semantic ownership/import boundaries, not fragile line numbers. Reuse existing test tooling; do not add a heavyweight architecture framework.

Diagram rendering and generation-time feature extraction must remain owned by the existing diagram-generation module Worker boundary.

## Validation

Run the narrowest relevant commands and the exact new CI commands. Include, as applicable:

```bash
node --test tests/web/*.test.mjs
npx playwright test
npx playwright test --config <performance-config>
python -m pytest <focused browser/packaging tests>
python tools/benchmark_diagram_layout.py run --quick
ruff check gbdraw/
```

Adapt to the repository's actual configuration; do not create duplicate configs unnecessarily.

No broad generated-artifact refresh is expected. If a generated artifact must change, use its canonical generator and list every changed artifact.

## Completion report

Report:

1. reconciled audit findings
2. files changed
3. commands/tests run and outcomes
4. pre-existing failures
5. structural invariants now enforced
6. CI jobs now executing previously dormant tests
7. remaining P0-1 and P0-2 work
8. generated artifacts changed and regeneration command
9. risks and uncertainties

Keep the patch minimal and PR-sized. Success means high-value tests run in the correct required lane, failures identify the correct subsystem, and ownership regressions fail immediately.
