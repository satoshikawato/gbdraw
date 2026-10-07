# Web periodic audit and promotion checklist

The 2026-09-30 Web GUI audit
([remediation proposal, section 8.3](web-gui-audit-20260930/01_REMEDIATION_PROPOSAL.md#83-ワークフロー);
[W8, section 3](web-gui-audit-20260930/remediation/W8_verification_workflow.md))
decided that periodic audits use no new gate or workflow. Each audit is a
checklist item in the `dev` to `main` promotion pull request. Steps 5 and 6 are
Phase F of the
[Web owner-coupling prevention plan](WEB_OWNER_COUPLING_PREVENTION_IMPLEMENTATION_PLAN_2026-10-05.md#phase-f--continuous-structural-audit)
in this checklist form, which the Owner chose on 2026-10-06 over a scheduled
workflow. This page holds the
procedure and the checklist. The admission rules for a `PROMOTION` stay in
[`WEB_CHANGE_POLICY.md`](WEB_CHANGE_POLICY.md); this page does not change them.

## When to audit

- Before each `dev` to `main` promotion.
- Or earlier, when about 30 runtime PRs have merged into `dev` since the last
  audit.

## Input

1. The capabilities changed since the last promotion. Run this from the
   repository root after `git fetch origin main dev`:

   ```bash
   git diff --name-status --no-renames origin/main...origin/dev | node --input-type=module -e "
   import { readFileSync } from 'node:fs';
   import { classifyChanges } from './tools/ci-impact-policy.mjs';
   const changes = readFileSync(0, 'utf8').trim().split('\n').map((line) => {
     const [status, ...paths] = line.split('\t');
     return { status, paths };
   });
   console.log(classifyChanges(changes).capabilities.join('\n'));"
   ```

2. The previous audit folder (`docs/internal/web-gui-audit-<date>/`) and the
   follow-ups it left open.

## Procedure

1. Dispatch `Tests` on `dev` with `tier=release` (the existing release stage;
   see [`SELECTIVE_CI.md`](SELECTIVE_CI.md)); it includes the
   `Vibrio full generation` job that otherwise only runs on `main`.
2. Run every recipe check, including the LOSATP recipes that the `recipe`
   pytest marker skips:

   ```bash
   PYTHONPATH=$PWD python docs/recipes/run_cli_scenarios.py --all --check
   PYTHONPATH=$PWD python docs/recipes/run_python_scenarios.py --all --check
   ```

   A stale artifact is regenerated with the same runner without `--check`. A
   recipe that fails before the comparison step is a product or documentation
   defect: record it as an audit row instead of regenerating.
3. Run the sweeps in [`tools/audit/`](../../tools/audit/README.md): parity
   sweep and replay, XSS sweep, viewport sweep, Gallery round trip, and session
   load timing.
4. Audit the changed areas by hand, with a fixed time limit per area. Pick the
   areas from the capabilities in the input: input, Session/History,
   Generate/options, comparisons, Feature editing, legend/Preview/export, and
   tracks/layout. Reproduce each suspected defect on `dev` before you record it.

5. Compare the owner graph of `main` and `dev`. Run from the repository root
   after `git fetch origin main dev`; it prints every frozen metric on both
   sides and exits 1 when one is higher on `dev`:

   ```bash
   D=$(mktemp -d)
   node tools/report-web-owner-graph.mjs --at origin/main --json > "$D/main.json"
   node tools/report-web-owner-graph.mjs --at origin/dev --json > "$D/dev.json"
   node -e "
   const { readFileSync } = require('node:fs');
   const [main, dev] = ['main', 'dev'].map((name) => JSON.parse(readFileSync(process.argv[1] + '/' + name + '.json', 'utf8')).summary);
   const flat = ({ projectionShapes, ...counts }) => ({
     ...counts,
     ...Object.fromEntries(Object.entries(projectionShapes).map(([domain, count]) => ['shapes:' + domain, count]))
   });
   const [before, after] = [flat(main), flat(dev)];
   const rows = Object.keys({ ...before, ...after }).map((metric) => [metric, before[metric] ?? 0, after[metric] ?? 0]);
   for (const [metric, onMain, onDev] of rows) console.log(metric.padEnd(34), String(onMain).padStart(4), String(onDev).padStart(4), onDev > onMain ? 'HIGHER ON DEV' : '');
   process.exitCode = rows.some(([, onMain, onDev]) => onDev > onMain) ? 1 : 0;
   " "$D"
   ```

   A metric that is higher on `dev` blocks the promotion, unless the R13
   baseline (`tests/web/owner-graph-baseline.test.mjs` on `dev`) records the
   higher count and an authority PR admitted it. To find which merge changed a
   metric, list the first-parent merges between the two:

   ```bash
   node tools/report-web-owner-graph.mjs --range origin/main..origin/dev --first-parent
   ```

   A `*` after a count marks a value that differs from the previous row, so
   the first row that carries it names the merge. The per-PR `Gate` and the push CI of `dev` already enforce
   the R13 baseline; this step adds the trend across the merges of one
   promotion.
6. Run the seeded live-vs-Generate random walk
   ([`live-generate-random-walk.promotion.spec.js`](../../tests/web/live-generate-random-walk.promotion.spec.js),
   run by `tests/web/playwright/promotion.config.js`; no PR or dev workflow collects it).
   Use a fixed seed and the default budget (20 steps on each of a Circular
   Result, a Linear Result, and a two-Result Circular batch), and rebuild the
   browser wheel first when Python under `gbdraw/` changed
   (`python tools/prepare_browser_wheel.py`):

   ```bash
   GBDRAW_RANDOM_WALK_SEED=20261006 npx playwright test --config=tests/web/playwright/promotion.config.js --workers=1
   ```

   After each step (feature fill with each scope choice, label text and
   visibility, Feature Visibility and specific color rules, legend rename,
   color, and sort, Undo, Redo, Result switch, Auto Reflow), the walk requires
   the displayed Result to equal the Result Generate draws from the same draft
   (R3, PD-OI-066), minus `tests/web/contracts/live-generate-parity-allowed.json`.
   The log starts with the seed; a failure names the seed, the fixture, the
   step index, the steps so far, and the differences. The same seed and
   `GBDRAW_RANDOM_WALK_STEPS` replay the same walk; `-g "<fixture name>"`
   replays one fixture. Each mismatch is a finding (OV-xx): log it, and either
   fix it or mark the matching case `test.fail` in
   `tests/web/live-generate-parity.playwright.spec.js` per R3 (naming the
   finding) before the promotion. Pick a new seed for each promotion and
   record it in the promotion pull request.

## Carrying evidence forward

A light change can move `dev` after promotion evidence was collected on an
earlier commit E. Evidence from E counts for the promotion head H when
`node tools/ci-impact.mjs classify --base <E> --head <H>` reports the matching
verdict true ([`SELECTIVE_CI.md`](SELECTIVE_CI.md#carrying-evidence-to-a-later-commit)).
Every verdict requires E to be an ancestor of H.

| Evidence | Verdict |
| --- | --- |
| `Tests` dispatched on `dev` with `tier=release` | `releaseEvidenceCarries` |
| Gallery refresh, Gallery tutorial media, and docs GUI capture checks, and the Gallery artifact manifest | `generatedArtifactChecksCarry` |
| Recipe runs (step 2), `TestOutputComparison`, the `tools/audit/` sweeps (step 3), and the `main`-written Session fixture tests | `localTestEvidenceCarries` |

The promotion body names E, links E's evidence, and includes the `classify`
output. Rerun on H the evidence whose verdict is false. The hand audit (step 4)
names the `dev` SHA it audited and needs no carry rule.

This rule carries checklist evidence only. `Promotion / gate` still requires
exact-SHA `Dev staging / gate` and `Gallery readiness / gate` evidence on H.
Each push to `dev` earns that evidence; after a light move, the runs inherit
their direct parent's evidence and run only the jobs the change needs
([light-change inheritance](SELECTIVE_CI.md#light-change-inheritance)).

## Output

Write `docs/internal/web-gui-audit-<date>/README.md` in the same format as
[the 2026-09-30 audit](web-gui-audit-20260930/README.md): ID, severity (P1, P2,
P3), symptom, and cause, with the commit audited. For every confirmed defect
class, record one of two decisions:

- add a guard (name the test and the PR stage that runs it); or
- add no guard, and say why.

The audit is finished when every P1 is fixed, or the Owner has given an
explicit waiver for it.

## Promotion PR checklist

Copy these items into the `PROMOTION` section of the `dev` to `main` pull
request and complete them.

```markdown
- Periodic audit: `docs/internal/web-gui-audit-<date>/` (audited `dev` SHA: <sha>)
  - [ ] Changed capabilities since `main` <sha>: <list>
  - [ ] `Tests` dispatched on `dev` with `tier=release`: <run URL>
  - [ ] `run_cli_scenarios.py --all --check` and `run_python_scenarios.py --all --check` pass, or each failure is an audit row: <rows>
  - [ ] `tools/audit/` sweeps run; findings recorded: <rows or "none">
  - [ ] Changed areas audited by hand, with the time spent per area
  - [ ] Each confirmed defect class has "guard added" (test, PR) or "no guard" (reason)
  - [ ] Every P1 is fixed, or has an Owner waiver: <links>
  - [ ] Owner graph, `origin/main` against `origin/dev` (Procedure step 5): no frozen metric is higher on `dev`, or each one is recorded in the R13 baseline by an authority PR: <PR>
  - [ ] Random walk (Procedure step 6): seed <seed>, 20 steps per fixture, no mismatch, or each mismatch is an OV-xx that is fixed or marked `test.fail` in the parity spec: <rows>
  - [ ] Evidence carried from an ancestor: none, or E <sha>, the `classify --base <E> --head <H>` output, and the verdicts used
```

### First promotion after the 2026-09-30 audit

That promotion carries user-visible changes that the release notes must state.
Check each one on the promotion head, or on an ancestor under
[Carrying evidence forward](#carrying-evidence-forward), and link the evidence:

- [ ] PV-08 and N-01 (D-24): Web-default Circular output changes. Long species
  lines wrap at word boundaries, and the wrap width uses an approximate glyph
  size. `definition_font_size` 18 is the default, so 18 counts as not explicit
  and may wrap; any other value is explicit and never wraps. Check the Gallery
  and the tutorial images.
- [ ] CO-05 (D-23): outfmt 6 tables with more than 12 columns (for example
  `-outfmt "6 std qlen slen"`) are read by their first 12 columns, which must
  have the right types. Columns 13 and later are dropped with an INFO log;
  `main` misread these tables silently. A CLI run given a missing or unreadable
  comparison file now fails instead of skipping it.
- [ ] Q-FRAME (D-18): comparison tables and CLI `-b` are read in the search
  frame (after record selection and crop, original strand). CLI users who
  combined reverse complement with `-b` see their coordinates mean something
  different. Main-saved Sessions are converted once on load.
- [ ] D-17: a Web PDF is 75% of its former physical size (CSS px converted to
  pt), which now matches CLI PDF and the PNG DPI.
