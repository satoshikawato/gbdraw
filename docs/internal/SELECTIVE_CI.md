# CI validation tiers and impact routing

`tools/ci-impact-policy.mjs` owns path classification and required job registries.
`tools/ci-impact.mjs` reads Git changes, verifies inherited evidence, and validates
aggregate results. Workflows execute the selected jobs; they do not duplicate the
classification rules.

## Validation tiers

| Tier | Trigger | Required coverage |
| --- | --- | --- |
| PR | Every PR into `dev` | Changed subsystem, cross-layer smoke, architecture/Product policy; all functional Playwright in eight shards for Web runtime, session, Gallery, LOSAT, shared-test, and unknown changes; Python 3.11 primary |
| Integrated dev | Every push to `dev`; `Tests` dispatch with `tier=dev` | Core on Python 3.10/3.11/3.12; recipes, Gallery, browser/package integration, offline GUI, LOSAT cache, all functional Playwright in eight shards, performance smoke on 3.11 |
| Release / S11 | Explicit `Tests` dispatch on `dev` with `tier=release` | Every dev functional job, additional recipe/Gallery/browser acceptance on 3.10/3.12, exhaustive non-browser slow tests on all three versions, Vibrio full generation (`vibrio-generate-release`), exact-candidate Gallery readiness |

PR feedback targets 5–8 minutes without functional Playwright and 15–20 minutes
with it. Integrated functional validation targets 15–20 minutes where practical.
S11 remains intentionally expensive. The separate
`Gallery publication` workflow retains common-nine browser parity and its
three-trial performance projection, with the complete-refresh budget available
through its existing explicit dispatch. Deployment stress checks also remain.

`Release / gate` is CI evidence for S11, not certification of the entire release
process. S11 still requires its candidate freeze, package install/offline and
cross-platform evidence, complete-refresh/performance budgets, and source checks.
The tag publication workflow and `tools/check_release_source.py` are unchanged;
passing a tier never authorizes tagging, publication, or deployment.

### Python-version coverage

The inexpensive core matrix remains on 3.10/3.11/3.12 because language support,
TOML loading (`tomli` on 3.10), dependency resolution, and package/Python APIs can
differ by host interpreter. Browser rendering itself uses the packaged Pyodide
interpreter, not the runner's Python version. Host-version repetition still tests
packaging and browser-test adapters, so exhaustive recipe/Gallery/browser version
coverage remains mandatory in release acceptance. Python 3.11 already runs each
of those surfaces in the shared functional jobs.

All non-browser slow tests move to S11, except the three package-build integration
checks in `test_web_packaging.py`, which also run in the dev Browser job. Offline
GUI browser contracts remain on dev, including real Linear LOSAT generation;
recent integrated runs found failures uniquely in that path.

`npm run test:web:vibrio-generate` (one real Vibrio generation test, up to 20
minutes) runs in the release-only `vibrio-generate-release` job and in the
main-push `deploy_web.yml` verification. Before this, only main saw it, so a
stale expectation (T13, request schema 7) surfaced after promotion.

## Changed paths → capabilities → required jobs

Plans contain the ordered union of affected `capabilities`. `impact` is the
highest-ranked display label; it does not discard other capabilities. The job
registry unions their contributions, then emits unique jobs in registry order.

| Capability | Representative paths | Selective PR jobs |
| --- | --- | --- |
| `metadata` | Existing development-tool allowlist, `.gitignore`, citation/licenses | No test job |
| `documentation` | Root Markdown and `docs/`, excluding the five exact policy documents below | Recipes |
| `policy-documentation` | Architecture, Web, Product Impact, Option Integrity, and selective CI policy documents | Web change budget |
| `python-core` | Known Python core/API/config/I/O owners and data | Architecture, core, lint, Web contracts, browser smoke |
| `renderer` | Render/SVG/diagram/canvas/features/labels/layout owners | Architecture, core, lint, Web contracts, browser smoke |
| `web-runtime` | `gbdraw/web/index.html`, first-party JS | Architecture, Web contracts, browser smoke, functional Playwright |
| `session-persistence` | Python session codecs; JS session/config/history owners | Architecture, core, recipes, lint, Web contracts, browser smoke, functional Playwright |
| `gallery` | Gallery inputs and their generator owners | Architecture, Gallery, Web contracts, browser smoke, functional Playwright |
| `losat-integration` | Comparison owners, LOSAT JS/workers/Wasm | Architecture, core, lint, Web contracts, browser smoke, functional Playwright |
| `tests-only` | Shared fixtures/harness and unclassified tests | Full PR tier, functional Playwright |
| `packaging` | Dependencies, manifests, vendored runtime, packaging tools | Full PR tier |
| `ci-only` | Workflows, planner/tests, and Playwright configuration | Full PR tier |
| `full` | Unknown or invalid paths | Full PR tier, functional Playwright |

The `policy-documentation` paths are exactly `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`,
`docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`,
`docs/internal/PRODUCT_IMPACT_RATCHET.md`, `docs/internal/SELECTIVE_CI.md`,
and `docs/internal/WEB_CHANGE_POLICY.md`. Other `docs/` paths remain ordinary
documentation, including working templates and acceptance reports.

The full PR tier job IDs are `web-change-budget`, `core-pr`, `recipes-standard`,
`gallery`, `lint`, `web-contracts-pr`, and `web-pr-smoke`. `playwright-functional`
is the eighth PR job. The `web-runtime`, `session-persistence`, `gallery`,
`losat-integration`, `tests-only`, and `full` capabilities add it under both the
selective and the full decision, so an evidence fallback or the
`architecture-change` label keeps it whenever one of those capabilities is
present. Other full-tier routes do not run it before merge: `ci-only`,
`packaging`, and fallbacks of `metadata`, `python-core`, or `renderer` plans.
Integrated `dev` staging runs it on every push that does not inherit parent
evidence. On pull requests it is the same job, shard matrix, and command
as in `dev` staging.

`web-contracts-pr` groups the existing fast JS contracts and non-slow Python
browser suite, which together took approximately 65–70 seconds historically.
`web-pr-smoke` runs `npm run test:web:pr-smoke`, the `@pr-smoke` Playwright cases
counted under [Smoke inventory](#smoke-inventory-and-regression-retention). When a
PR plan requires `web-pr-smoke`, the `gallery` job also runs the common-nine Gallery
first-Generate parity command, `npm run test:web:gallery-publication`;
`tests/ci/ci-impact-cli.test.mjs` keeps that command out of `web-pr-smoke`.
The `web-contracts-pr` and `web-pr-smoke` jobs execute in parallel instead of
sharing one 15-minute serialized budget. Both jobs keep independent working
directories, dependency installs, wheels, and browser state. The contracts job
also runs when the trusted plan requires `web-pr-smoke`, so the pre-redesign base
planner still requires all three original suites during rollout.

Web unit/browser tests and their helpers route to the Web, session, Gallery, or
LOSAT capability they exercise. Adding a normal Web regression therefore does not
force a full PR route. Shared fixtures and unclassified tests remain conservative
because they can affect several suites. Unknown production roots also remain full;
recognizing known subsystem owners does not make arbitrary future paths safe.
Renames/copies classify both endpoints, deletions classify their old path, and
malformed/empty diffs or missing Git objects fail closed.

### Evidence required for selection

Documentation-only PRs into `dev`, including exact policy documents and accompanying
allowlisted metadata, run only the changed documentation or policy checks. Their
`DOCUMENTATION_ONLY_PR` plan contains no inherited evidence and does not query base
staging. Missing, failed, unfinished, or inaccessible base staging therefore cannot
block the PR or expand it to the full tier. The required jobs and trusted-base
checks still have to pass.

Other narrower PR routes require successful exact-SHA `Dev staging / gate`
evidence at the PR base. No older successful SHA can substitute; unavailable
evidence selects the complete tier. Mixing runtime, CI, packaging, or unknown
paths into documentation changes retains those capabilities' required coverage.
`architecture-change`, control-plane changes, dependencies, unknown paths,
unclassified/shared test inputs, and explicit dispatches require the full tier.

Every runtime/subsystem change runs complete integrated dev and Gallery coverage.
Metadata, documentation, and policy-documentation can inherit direct-parent
staging evidence. Metadata runs no test jobs; ordinary documentation runs recipes;
policy documentation runs Web change budget in `Tests`. When the direct parent has
no successful `Dev staging / gate` (missing, failed, cancelled, or unfinished), such
a `dev` push runs the complete dev tier with basis `INHERITED_EVIDENCE_UNAVAILABLE`
instead of failing its plan. Documentation-only Gallery
publication skips browser and performance when direct-parent Gallery readiness is
successful. A missing or unfinished direct-parent result still fails
documentation-only Gallery planning rather than running the full tier.

A later `dev` push does not cancel a running `Tests` run. One run executes, the
newest later push waits, and GitHub cancels an older waiting run when a newer
push replaces it. Every started staging run therefore finishes with exact-SHA
evidence. A SHA whose waiting run was replaced has no evidence of its own; the
newer run tests a tree that contains it, and a documentation-only push on top of
it runs the complete dev tier. Pull request runs and manual dispatches still
cancel the run they supersede.

## Smoke inventory and regression retention

Smoke retains app/Worker boot and Circular generation, Linear multi-record
rendering, mounted Feature/Label/Legend edits, export startup, Save/Load/Generate,
Circular and Linear mode transitions, invalid-annotation preservation, and invalid
composite-session rejection. Full functional discovery includes every smoke case.

`tests/ci/playwright-inventory.test.mjs` uses Playwright's expanded `--list` output.
It counts nested and parameterized cases and rejects fewer than 8 or more than 19;
counting tag strings in source missed the earlier growth to 27 cases. No assertions,
case bodies, test timeouts, or full-functional exclusions changed in this redesign.
The exact relocated cases and before/after counts are in the
[measurement report](CI_TIER_REDESIGN_2026-09-11.md). Issue #601 later raised the
upper bound from 13 to 19 for six existing cases
([inventory note](issue-601-pr-smoke-inventory-20260927.md)). The smoke config
collects 19 cases, so a new `@pr-smoke` case requires raising that bound.

The smoke config retains one worker and zero retries. Three local two-worker
repetitions passed, but these are not repeated GitHub runner evidence. One worker
already meets the projected target, so CI does not adopt unverified parallelism.
Full functional CI runs eight shards with two retries and two workers per shard,
on pull requests and in `dev` staging alike. Playwright's `--shard` splits by
case count in spec file order; it cannot split by duration. With equal counts,
four hosted shards took 13 to 34 minutes on `52be34a7` and `a5deed4a`, because
the heaviest composite and Session cases sit in a few files. Each matrix entry
therefore runs its own list of spec files, passed as Playwright CLI file
filters, from `tests/ci/functional-shards.json`. The file records the measured
minutes of each spec file and the eight lists.
`tools/balance-functional-shards.mjs` fills the lists longest file first into
the shard with the least work; within a shard, the two workers split the cases.

The lists replaced Playwright's undocumented `PWTEST_SHARD_WEIGHTS` case-count
weights. In the green runs 37010713759, 37025041812, 37026343792, 37027373572
and 37028524962 (512 cases in 80 files), the slowest weighted shard ran its
cases in 14.9 to 16.5 minutes (mean 15.7) and the fastest in 8.7 to 12.1
minutes. Over 20 green `dev` and pull request runs, the slowest shard job took
15.7 to 17.8 minutes, including about 1.5 minutes of setup. Replaying each
run's case durations on two workers, the file lists give a slowest shard of
13.2 to 14.3 minutes (13.7 to 15.1 when the replayed run is left out of the
measurement that builds the lists); equal case counts give 19.7 to 21.0
minutes. The largest file, `linear-multi-record.playwright.spec.js` (about 17.5
minutes on one worker), fits in one shard because its two workers share its
cases.

Rebalance when shard times drift apart, or when the inventory check reports a
new spec file. Download the shard artifacts of a green run with
`gh run download <run-id> -p 'playwright-functional-shard-*' -D <dir>`, then run
`node tools/balance-functional-shards.mjs <dir>/*/functional-report.json`.
Without arguments, the tool keeps the recorded minutes and estimates a new file
from its case count. The inventory check fails when a spec file is unassigned,
assigned twice, or no longer exists. It also lists every shard with its file
filters to prove that the matrix partitions the full suite without omissions or
duplicate execution.

The job limit is 45 minutes; case timeouts are unchanged. Hosted
shard 1 still exceeded 30 minutes with two workers on `9dd6359c`: its composite
placement took 7.2 minutes, a composite Session round trip took 4.4 minutes,
and a Definition completion wait and retries consumed the remaining budget.
The Definition test now triggers an actual edit; the job budget also allows
for the measured composite workload, setup and existing retries. The line
reporter names the running case, and the JSON report records durations,
retries and final results alongside traces on every completed job. A job killed
at its limit may not reach upload; its last running case remains in the log.

## Trust and aggregate checks

PR planning and `PR / gate` execute the PR base's helper under `.ci-trusted-base`.
Candidate helpers are test inputs only. The documentation-only PR exemption
does not supply integrated-dev or promotion evidence. A missing or invalid base
helper fails;
there is no candidate fallback. `Web base policy (trusted base)` remains a separate
`pull_request_target` check that treats candidate files as Git data.

The redesign PR therefore receives the old base's complete PR route. New selective
routing can take effect only after reviewed integration into protected `dev`.
Additional contract failures also block the old gate's unknown-job validation.
No branch protection or required external status name changes.

Pull request functional Playwright is introduced the same way. The workflow
change that lets `playwright-functional` run on pull requests and adds it to
`PR / gate` merges first. The base planner does not select the job yet, so the
job is skipped and the gate accepts it as a skipped unknown job. The planner
change that selects it merges next and applies to pull requests based on that
commit.

`PR / gate`, `Dev staging / gate`, `Gallery readiness / gate`, and `Promotion / gate`
remain stable. `Release / gate` is separate and requires every release-profile
job, including both exhaustive matrices, plus exact-candidate Gallery readiness.
A release dispatch cannot emit ordinary `Dev staging / gate` success.

Plans use schema 2 and include the capability union. The active helper validates
its own exact schema, profile, workflow SHA, required job list, and inherited
evidence. It rejects missing/skipped/failed/cancelled required jobs and unexpected
failures in optional/additional jobs. Matrix aggregates must succeed. Workflows
start on all relevant events and filter at job level; they never use workflow
`paths` filters that could omit an aggregate.

## Verification and rollback

Install locked Node dependencies, then run `node --test tests/ci/*.test.mjs` and
`node --test tests/web/architecture-contracts.test.mjs`. For browser execution,
prepare the browser wheel and run `npm run test:web:pr-smoke` plus
`npm run test:web:gallery-publication`. Full functional coverage remains
`npm run test:web:functional-full`.

The [timing CSV](CI_TIER_TIMINGS_2026-09-11.csv) records each historical job's setup,
dependency, wheel, test, teardown, and wall-clock times with direct job URLs.
Plans report capabilities, required jobs, evidence links, and fallback reasons.

Rollback this coherent CI change together: policy, adapter, workflow, smoke tags,
and policy/inventory contracts. Preserve trusted-base evaluation and stable gate
names. Restoring exhaustive matrices to dev changes cost, not release criteria.
