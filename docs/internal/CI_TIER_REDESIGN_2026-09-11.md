# CI tier redesign measurements — 2026-09-11

Starting `origin/dev`: `8d6c914ecf5bb12437e1c8370d35a28ea4b0af6a`.
Scope: CI workflows, routing policy, tests/tags, and CI documentation only.
No runtime production source, fixture, generated public figure, or release
publication/source checker changed. Remaining 05A4 defects are unchanged.

## Historical cost

Source: GitHub Actions run/job/step timestamps, collected on 2026-09-11 for
PRs #502–#505 (including earlier #503 attempts) and the four latest dev runs.
The [per-job CSV](CI_TIER_TIMINGS_2026-09-11.csv) includes every non-skipped job
in these Tests runs and the first redesign rollout, with direct evidence URLs. Partial/cancelled runs are identified;
they are not counted as successful verification. Setup includes checkout and runtime
setup; dependency time is the grouped install step; tests includes checks/planning;
teardown includes uploads and post steps. Unattributed time captures runner overhead
and timestamp rounding. Combined install steps do not expose exact independent
pip/npm/Chromium or embedded wheel durations, so those are not falsely separated.
Runner minutes are the sum of job wall durations, not rounded billed minutes.

| PR | Successful Tests run | Web PR smoke | Playwright phase | Tests critical path | Tests runner minutes |
| --- | --- | ---: | ---: | ---: | ---: |
| #502 | [34543831540](https://github.com/satoshikawato/gbdraw/actions/runs/34543831540) | 14:47 | 12:10 | 17:17 | 31.75 |
| #503 | [34548658285](https://github.com/satoshikawato/gbdraw/actions/runs/34548658285) | 13:21 | 10:59 | 15:46 | 30.63 |
| #504 | [34550113922](https://github.com/satoshikawato/gbdraw/actions/runs/34550113922) | 14:06 | 11:37 | 16:23 | 31.00 |
| #505 | [34551490388](https://github.com/satoshikawato/gbdraw/actions/runs/34551490388) | 14:22 | 11:49 | 16:31 | 30.23 |

Each successful PR ran eight Tests jobs, plus the separate trusted-base check
and three CodeQL language jobs (12 actual jobs total). Aggregate names and those
four external jobs are preserved. Critical path above runs from the first Tests
job start to PR gate completion, including inter-job scheduling; external checks
are not added serially.

Median phase seconds across sampled successful jobs (individual values and all
matrix members remain in CSV):

| Job | Setup | Dependencies | Wheel | Tests/checks | Teardown | Wall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Web PR smoke | 34.5 | 40 | 2.5 | 771.5 | 1 | 854 |
| Core PR | 30 | 14 | 0 | 203 | 0 | 253 |
| Recipes standard | 30 | 18 | 0 | 89 | 0 | 138 |
| Gallery | 29 | 14 | 0 | 210 | 0 | 254 |
| Web change budget | 75 | 0 | 0 | 131 | 0 | 208 |
| Browser | 26.5 | 42 | 2 | 154 | 1 | 229 |
| Playwright functional shard 1/4 | 35 | 47.5 | 3 | 755.5 | 0 | 845 |
| Playwright functional shard 2/4 | 30 | 41 | 2 | 239 | 0.5 | 312 |
| Playwright functional shard 3/4 | 34 | 42 | 2 | 458.5 | 0.5 | 538 |
| Playwright functional shard 4/4 | 25 | 45.5 | 2 | 549.5 | 0.5 | 629 |
| Recipe acceptance 3.10 / 3.12 | 27.5 / 28.5 | 15.5 / 20 | 0 | 83.5 / 96.5 | 0 | 133.5 / 147 |
| Gallery acceptance 3.10 / 3.12 | 30.5 / 27.5 | 18.5 / 20 | 0 | 214 / 200 | 1 / 0 | 266 / 244 |
| Browser acceptance 3.10 / 3.12 | 33.5 / 27.5 | 44.5 / 45 | included in install | 49 / 54 | 1 / 0.5 | 126 / 130.5 |
| Slow 3.10 / 3.11 / 3.12 | 32 / 28.5 / 29 | 18 / 18 / 19 | 2.5 / 3 / 3 | 74.5 / 71.5 / 77.5 | 1 / 0 / 0 | 129 / 120.5 / 131.5 |

### Why passing smoke hit 15 minutes

The original job serialized dependency installation, wheel preparation, roughly
17–19 seconds of JS contracts, 48–52 seconds of Python browser tests, and about
11–12 minutes of Playwright. Smoke had grown to 27 expanded cases at the starting
SHA (28 in #502, 30 in early #503). The previous regex counted only top-level
`tag: '@pr-smoke'` declarations, missing title tags, loops, and nested contracts.
Several cases perform repeated Generate, History, Save/Load, or large Vibrio
resource operations; they are integration regressions, not short boot probes.

[Job 103097229281](https://github.com/satoshikawato/gbdraw/actions/runs/34545278484/job/103097229281)
lasted 15:12 and logged **30 passed (12.6m)** immediately before cancellation.
[Job 103102959643](https://github.com/satoshikawato/gbdraw/actions/runs/34547226783/job/103102959643)
lasted 15:08. These runs exhausted the serialized job budget; assertion success
could not rescue a cancelled required job. A separate 6:27 rerun was interrupted
and is not counted as a timeout. The primary remedy removes serial work and
reduces smoke selection. Its job timeout decreases from 15 to 10 minutes.

### Unique failures and duplicated setup

The successful heads of #502–#505 had no failing assertion in any required PR job.
Earlier #503 PR gates failed because smoke was cancelled. Sampled dev Browser
jobs [103096237619](https://github.com/satoshikawato/gbdraw/actions/runs/34545113530/job/103096237619)
and [103110884127](https://github.com/satoshikawato/gbdraw/actions/runs/34549970407/job/103110884127)
uniquely failed `test_offline_gui_linear_losat_generation_populates_cache` while
core, recipes, Gallery, version-acceptance, slow, and all four functional shards
passed. That Browser contract remains mandatory on dev. No unique failure was
observed in those other jobs in this sample; that is not proof of redundancy.

Pip/npm caches already existed. Web dependency installation usually cost 38–44
seconds, and wheel preparation 2–3 seconds. Each functional shard repeats about
40–48 seconds of setup dependencies, but sharing a 2–3-second wheel via an extra
artifact job would add coordination/startup cost. No shared mutable state, wheel
artifact pipeline, or Chromium cache is introduced. Splitting JS/Python contracts
adds one roughly 75-second setup while removing their 65–70-second work from the
browser critical path. Removing 17 browser journeys saves much more.

Trusted-base checkouts previously cost about 25 seconds each in planning and the
PR gate. Those two helpers require only `tools/`; their checkouts now use sparse
checkout. This preserves exact base revisions and helper imports without copying
candidate/runtime state. Actual savings are reported separately from projections.

## Before / after

| Measure | Before | After design / local evidence |
| --- | --- | --- |
| Web-only PR Tests jobs | Planner, architecture, core, recipes, Gallery, lint, serialized Web smoke, gate | Planner, architecture, Web contracts, Playwright smoke, gate |
| Number of Tests jobs | 8 | 5 selective; 9 for this rollout/full-fallback PR |
| Including trusted-base + three CodeQL jobs | 12 | 9 selective; 13 for rollout/full fallback |
| Required browser cases | 27 | 10 |
| Full functional browser inventory | 167 | 167, same cases after stripping smoke tags |
| PR Tests wall | 15:46–17:17 | Projected 5–7 minutes, pending hosted measurement |
| PR Tests runner minutes | 30.23–31.75 | Projected 10–13 selective; 20–24 full fallback |
| Dev Tests jobs / runner minutes at starting SHA | 25 / 95.65 | 16 jobs; roughly 72 minutes plus retained package-build checks, projected |
| Dev Tests wall at starting SHA | 16:49 | Same four full browser shards; approximately 15–20 minutes, projected |
| Dev functional coverage | Core, recipes, Gallery, browser, LOSAT/cache, functional/performance Playwright, offline GUI | Retained; package builds stay on 3.11; multi-version core retained |
| Supported-version and exhaustive slow acceptance | Every dev push | Explicit release profile; same matrices and commands mandatory under Release gate |
| Release source and package publication checks | Existing exact source/admission/build/install path | Unchanged |

The projection uses historical durations of **only the ten retained cases**:
138.6 / 130.7 / 136.4 / 138.0 seconds in #502–#505, plus setup, planning, gate and
scheduling. It does not use the faster local machine as a hosted-run prediction.
Local Chromium on Linux/WSL, Playwright 1.61.1, zero retries: 10/10 passed with
one worker in 87.09 seconds; two workers passed 30/30 over three repetitions in
153.24 seconds total, with no skips or flakes. CI keeps one worker because the
parallel experiment did not run repeatedly on GitHub-hosted hardware.

Actual changed-path replays from the historical PRs (with successful exact-base
staging evidence available):

| PR | Affected capabilities | New Tests jobs including planner/gate |
| --- | --- | ---: |
| #502 | Web runtime + session + Gallery tests | 9; the independent capabilities jointly require every PR suite |
| #503 | Session persistence and history tests | 8; Gallery deferred to dev |
| #504 | Documentation + Gallery + packaging verifier | 9; packaging requires full coverage |
| #505 | Web runtime + session/history | 8; Gallery deferred to dev |

The five-job Web example is an editor/runtime change plus its normal Web tests.
Session or packaging contributions must not disappear merely to reduce job count.
All four historical replays benefit from the shorter browser critical path.

### First hosted rollout measurement

[Tests run 34563716800](https://github.com/satoshikawato/gbdraw/actions/runs/34563716800)
on `8512e5650edebabf2b006ba312e9bdfc1fd70a97` passed all nine jobs in **6:20**,
using **23.05 runner minutes**. The base's original full route was enforced:
core, recipes, Gallery, lint, architecture, and browser smoke all ran, alongside
the separated contracts. The independent trusted-base and CodeQL checks passed.

Web smoke took **3:38**, with **2:20** in Playwright; contracts took **2:24**.
Trusted helper checkout took **1 second** in both planner and PR gate. Core
(4:13), rather than browser smoke, set the full-route test critical path.
The full-route wall time fell about 62% against #505's 16:31. This measurement
precedes the final code-plus-test routing refinement and evidence-only report
update; final-head verification is recorded on the PR. Selective five-job timing
and post-merge dev measurements must not be claimed as observed from this run.

## Tests moved to later tiers

The following 17 expanded browser cases leave PR smoke and remain mandatory in
the full functional suite on dev and release. Their bodies and timeouts are
unchanged:

- `annotation-download.playwright.spec.js` — annotation TSV download/re-import and Python round-trip offline (1440px).
- `annotation-download.playwright.spec.js` — annotation TSV download/re-import and Python round-trip offline (390px).
- `composite-session-resources.playwright.spec.js` — Session export Vibrio composite resources survive minimal Save, fresh Load and Generate.
- `contracts/active-result-edit-transaction.playwright.spec.js` — visible tRNA fill scope remains coherent through History, Session, Generate, and SVG export.
- `draft-placement-capability.playwright.spec.js` — circular draft placement capability survives history and dirty session restore.
- `draft-placement-capability.playwright.spec.js` — linear draft placement capability survives history and dirty session restore.
- `history-inputs.playwright.spec.js` — Label Mode keyboard edit has one Undo step and survives generation and fresh Load.
- `history-inputs.playwright.spec.js` — Label Mode pointer edit has one Undo step and survives generation and fresh Load.
- `joint-display-placement.playwright.spec.js` — joint rotation placement and saved drafts in circular.
- `joint-display-placement.playwright.spec.js` — joint rotation placement and saved drafts in linear.
- `linear-typography.playwright.spec.js` — independent Linear typography follows linked, imported, and History journeys.
- `mode-record-identity.playwright.spec.js` — committed Circular placement retains source identity through mode history and Save/Load.
- `mode-record-identity.playwright.spec.js` — ungenerated Circular placement retains source identity through mode history and Save/Load.
- `mode-transition-editor-state.playwright.spec.js` — circular label intent survives mode Undo/Redo and Save/Load/Generate.
- `mode-transition-editor-state.playwright.spec.js` — linear label intent survives mode Undo/Redo and Save/Load/Generate.
- `preview-navigation.playwright.spec.js` — background pan preserves screen displacement after text selection and zoom.
- `preview-navigation.playwright.spec.js` — preview pan leaves feature and match gestures available.

For a Web-subsystem PR (including its Web tests), the entire core, recipe, Gallery, and Python lint jobs
also move to dev; their suites still run in PRs that affect those owners, and in
full fallback. Python renderer PRs defer recipes/Gallery to dev. Session changes
retain core and recipe contracts in PR; Gallery changes retain Gallery tests.
The complete Python/browser contract suite stays in the parallel Web contracts
job whenever browser smoke is required.

Dev→release movement is explicit: `acceptance-supported-main` retains recipe,
Gallery and non-slow browser tests on Python 3.10/3.12; `slow-main` retains every
non-browser slow test on 3.10/3.11/3.12. The three package build tests also run on
dev Python 3.11. Offline slow browser tests remain on dev and release Python 3.11.

## Review evidence

CI path classification and job requirements remain owned only by
`tools/ci-impact-policy.mjs`; the adapter and workflows consume its result.
The coarse strongest-class decision is replaced by the capability union in the
same owner. No additional policy registry, runtime path, or persisted migration
reader is added. Changed-scope owner excess, path excess and compatibility burden
do not increase. Schema 1 remains executable only in the unchanged trusted PR
base during rollout; the candidate helper accepts schema 2 only.

Architecture/Product detectors, rules, privileged allowlists, accepted Product
decisions, trusted-base workflow, release source gate, publication workflow and
branch protection remain unchanged. Product outcome: no diagram/session/editor
behavior change. This CI-governance change requires normal maintainer review;
passing candidate tests is not that approval.

Local validation: 58 CI policy/adapter/gate/inventory contracts; 137 architecture
contracts; expanded functional inventory equality; smoke tags-only comparison;
10-case smoke and three two-worker repetitions. Hosted results and the final
integration state are recorded in the PR/handoff. S11 is not reaccepted here.

The retained JS contract command passed 390 tests. Python browser checks passed
23/24 initially; the replay failure used the globally installed older CLI
(session 40 versus candidate 41). Installing this checkout in an isolated venv
made that test pass without source edits. All three retained package-build
integration tests passed (29.68 seconds). Generated artifacts remain uncommitted.
