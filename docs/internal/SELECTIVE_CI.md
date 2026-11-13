# CI validation tiers and impact routing

`tools/ci-impact-policy.mjs` owns path classification and required job registries.
`tools/ci-impact.mjs` reads Git changes, verifies inherited evidence, and validates
aggregate results. Workflows execute the selected jobs; they do not duplicate the
classification rules.

## Validation tiers

| Tier | Trigger | Required coverage |
| --- | --- | --- |
| PR | Every PR into `dev` | Changed subsystem, cross-layer smoke, architecture/Product policy; [examples freshness](#examples-freshness) for changes that can alter rendering; all functional Playwright in eight shards, plus the performance smoke, for Web runtime, session, Gallery, LOSAT, shared-test, and unknown changes; Core on Python 3.10, and on 3.10–3.14 for [version-dependent changes](#python-version-coverage); other jobs on 3.11 |
| Integrated dev | Every push to `dev`; `Tests` dispatch with `tier=dev` | Core on Python 3.10–3.14, with line coverage on 3.14; recipes, Gallery, browser/package integration, offline GUI, LOSAT cache, all functional Playwright in eight shards, performance smoke on 3.11 |
| Release / S11 | Explicit `Tests` dispatch on `dev` with `tier=release` | Every dev functional job, additional recipe/Gallery/browser acceptance on 3.10/3.12, exhaustive non-browser slow tests on all three versions, Vibrio full generation (`vibrio-generate-release`), the render-stage universe (`render-stage-universe-release`), exact-candidate Gallery readiness, complete Web UI translations (`node tools/i18n/sync.mjs status --require-complete`); outside `Release / gate`, the promotion audit jobs (every recipe check, report-only `tools/audit/` sweeps, the user journeys, the random walk, and their summary) |

PR feedback targets 20 minutes without functional Playwright and 35 minutes
with it: the median wall time, from the first job's start to the last job's
end, of the latest green `Tests` run on the head commit of each of the last 20
pull requests merged into `dev`. The [weekly CI inventory](#weekly-ci-inventory)
reports both medians next to these targets and warns, without failing, when
one is above its target. Integrated functional validation targets 15–20
minutes where practical.
S11 remains intentionally expensive. The separate
`Gallery publication` workflow retains common-nine browser parity and its
three-trial performance projection, with the complete-refresh budget available
through its existing explicit dispatch. Deployment stress checks also remain.

`Release / gate` is CI evidence for S11, not certification of the entire release
process. S11 still requires its candidate freeze, package install/offline and
cross-platform evidence, complete-refresh/performance budgets, and source checks.
The tag publication workflow admits its source with `tools/check_release_source.py`
and publishes to PyPI only when a release-tier run's evidence carries to the
tagged commit
([Release evidence check](WEB_PERIODIC_AUDIT.md#release-evidence-check));
passing a tier never authorizes tagging, publication, or deployment.

### Line coverage

`Core (Python 3.14)` runs its selection under pytest-cov, writes the total line
coverage of `gbdraw/` (without `gbdraw/web/`) to the run summary, and uploads
`coverage.xml` and `htmlcov/` as the `coverage-<run id>-<attempt>` artifact for
30 days. `patch = ["subprocess"]` counts the CLI that tests start as
subprocesses. Its test step fails when the total falls below `fail_under` in
`[tool.coverage.report]` of `pyproject.toml`. The floor is set from dev-tier
results and only raised; lowering it needs the Owner. A pull request measures
it only when its plan runs `core` ([Python-version coverage](#python-version-coverage)).

### Python-version coverage

The inexpensive core matrix remains on 3.10 to 3.14 because language support,
TOML loading (`tomli` on 3.10), dependency resolution, and package/Python APIs can
differ by host interpreter. Browser rendering itself uses the packaged Pyodide
interpreter, not the runner's Python version. Host-version repetition still tests
packaging and browser-test adapters, so exhaustive recipe/Gallery/browser version
coverage remains mandatory in release acceptance. Python 3.11 already runs each
of those surfaces in the shared functional jobs.

`Core PR` runs on Python 3.10, the oldest supported version, so a pull request
catches a newer-only interpreter option or API (such as `python -P`) before it
reaches `dev`; mypy already checks for 3.10 but cannot see interpreter options.
Python 3.11 still runs on pull requests in the recipe, Gallery, Web contract,
and Playwright jobs, and the `dev` Core matrix covers 3.10 to 3.14.

A pull request whose plan has the `python-versions` capability also runs that
`dev` job, `core`: Python 3.10 to 3.14, with line coverage on 3.14. A changed
path adds the capability when `PYTHON_VERSION_PATHS` in
`tools/ci-impact-policy.mjs` lists it: `pyproject.toml`, and each package module
that imports a standard-library module whose output or availability differs
between these versions (`gzip`, `zlib`, `bz2`, `lzma`, `tarfile`, `zipfile`,
`tomllib` or `tomli`, `locale`) or that reads `sys.version_info`.
`tests/ci/ci-impact-policy.test.mjs` fails when a package module imports one of
those and is not listed. The label `all-python-versions` adds the capability to
any pull request; use it for a change that the list cannot see, such as one to
the Core job or its Python setup. The path keeps its own capability, so `Core
PR` and the other jobs still run.

Before this rule, a failure on one version reached `dev` first: a gzip header
test added by #1625 expected the header byte that Python 3.10 and 3.13+ write,
and `dev` Core failed on 3.11 and 3.12 for 4.9 hours on 2026-11-05. Replayed
over the last 20 merged pull requests that changed `gbdraw/` or
`pyproject.toml` (#1617 to #1650), the rule selects 7 of them, including
#1625, in 7 runs. At the median `dev` Core times of 12 green runs (3.10 17.5,
3.11 14.8, 3.12 18.2, 3.13 13.7, 3.14 12.6 minutes), each such run uses about
77 more runner-minutes, 537 over the 20. The Core jobs end before the slowest
functional shard, so the run takes longer only when the extra jobs wait for a
runner.

All non-browser slow tests move to S11, except the three package-build integration
checks in `test_web_packaging.py`, which also run in the dev Browser job. Offline
GUI browser contracts remain on dev, including real Linear LOSAT generation;
recent integrated runs found failures uniquely in that path.

`npm run test:web:vibrio-generate` (one real Vibrio generation test, up to 20
minutes) runs in the release-only `vibrio-generate-release` job and in the
main-push `deploy_web.yml` verification. Before this, only main saw it, so a
stale expectation (T13, request schema 7) surfaced after promotion.

The render-stage universe (`--render-stage-universe -m render_stage_universe`:
every option state of the six cases in `tests/utils/render_stage_equivalence.py`,
warm render against cold, with the declaration and expected-build checks) runs in
the release-only `render-stage-universe-release` job on Python 3.11, three
workers with `--dist loadscope` (the chains of a module share one module-scoped
fixture), up to 120 minutes. `slow-main` collects these tests but skips them
without the flag; running them there would cost about an hour on each of three
Python versions on every `dev` push. Before this job, no run used the flag, and
three stages that the base of a case could not cache went unseen (OV-3112). Core
renders each case base once (`test_each_case_base_has_no_undeclared_read_or_uncacheable_stage`)
so a base that cannot cache a stage fails on the pull request.

## Changed paths → capabilities → required jobs

Plans contain the ordered union of affected `capabilities`. `impact` is the
highest-ranked display label; it does not discard other capabilities. The job
registry unions their contributions, then emits unique jobs in registry order.

| Capability | Representative paths | Selective PR jobs |
| --- | --- | --- |
| `metadata` | Existing development-tool allowlist, `.gitignore`, citation/licenses | No test job |
| `documentation` | Root Markdown, changelog fragments under `changelog.d/`, and `docs/`, excluding the five exact policy documents below | Recipes |
| `policy-documentation` | Architecture, Web, Product Impact, Option Integrity, and selective CI policy documents | Web change budget |
| `python-core` | Known Python core/API/config/I/O owners and data | Architecture, core, examples, lint, Web contracts, browser smoke |
| `renderer` | Render/SVG/diagram/canvas/features/labels/layout owners | Architecture, core, examples, lint, Web contracts, browser smoke |
| `web-runtime` | `gbdraw/web/index.html`, first-party JS | Architecture, Web contracts, browser smoke, functional Playwright |
| `session-persistence` | Python session codecs; JS session/config/history owners | Architecture, core, recipes, examples, lint, Web contracts, browser smoke, functional Playwright |
| `gallery` | Gallery inputs and their generator owners | Architecture, Gallery, examples, Web contracts, browser smoke, functional Playwright |
| `losat-integration` | Comparison owners, LOSAT JS/workers/Wasm | Architecture, core, examples, lint, Web contracts, browser smoke, functional Playwright |
| `tests-only` | Shared fixtures/harness and unclassified tests | Full PR tier, functional Playwright |
| `python-versions` | `pyproject.toml` and the version-dependent package modules, beside their own capability; any path with the `all-python-versions` label | `core` ([Python-version coverage](#python-version-coverage)) |
| `packaging` | Dependencies, manifests, vendored runtime, packaging tools | Full PR tier |
| `ci-only` | Workflows, planner/tests, and Playwright configuration | Full PR tier |
| `full` | Unknown or invalid paths | Full PR tier, functional Playwright |

The `policy-documentation` paths are exactly `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`,
`docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`,
`docs/internal/PRODUCT_IMPACT_RATCHET.md`, `docs/internal/SELECTIVE_CI.md`,
and `docs/internal/WEB_CHANGE_POLICY.md`. Other `docs/` paths remain ordinary
documentation, including working templates and acceptance reports.

The full PR tier job IDs are `web-change-budget`, `core-pr`, `recipes-standard`,
`gallery`, `examples`, `lint`, `web-contracts-pr`, and `web-pr-smoke`.
`playwright-functional` is the ninth PR job. `core` joins a pull request only through the
`python-versions` capability. The `web-runtime`, `session-persistence`, `gallery`,
`losat-integration`, `tests-only`, and `full` capabilities add it under both the
selective and the full decision, so an evidence fallback or the
`architecture-change` label keeps it whenever one of those capabilities is
present. Other full-tier routes do not run it before merge: `ci-only`,
`packaging`, and fallbacks of `metadata`, `python-core`, or `renderer` plans.
Integrated `dev` staging runs it on every push that inherits neither parent
evidence nor the evidence of the pull request it merges. On pull requests it is the same job, shard matrix, and command
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

### Examples freshness

The `examples` job (`Examples (Python 3.11)`) runs on pull requests into `dev`
only, for every capability that can change what `tools/reproduce_examples.py`
renders: `python-core`, `renderer`, `session-persistence` (Session variant
figures), `gallery` (figures drawn from Gallery Sessions), `losat-integration`
(LOSATP comparison figures), and the full PR tier (`packaging`, `ci-only`,
`tests-only`, `full`; `examples/` and `tools/` paths classify as `full`). It
runs `python tools/reproduce_examples.py --check --jobs "$(nproc)"`, which
renders every figure of `tools/reproduce_examples_manifest.py` into a
temporary root, compares each output with the tracked file under `examples/`,
lists every figure that differs, and prints the command that refreshes them (`python tools/reproduce_examples.py --output-root . --figure
<id> ...`). An SVG is compared byte for byte. A PNG's pixels depend on the
machine's FreeType and fonts, so each PNG the tool writes carries, in its
`gbdraw-example-source` text chunk, the sha256 of what it was rasterized from:
the SVG bytes and sizing options gbdraw hands to CairoSVG's `svg2png`, recorded
at that call (inside the gbdraw command for a CLI PNG) for each panel, the
recipe, and, for a composite, the code that builds it. The check compares that digest, and a tracked PNG without one is
stale, so a refresh on any machine passes. A failed run uploads what it
rendered as the `examples-<run>` artifact. Every input a figure reads, BLAST
results included, is in the repository or the test-data Release, so the check
runs no BLAST+ and fails when a figure cannot render. `dev` pushes do not run the job: the pull request that changes a
drawing refreshes its own figures.

The 16 `@comparison-contract` Playwright cases run once per CI run. Their pytest
wrapper, `tests/test_linear_comparison_browser_contracts.py`, runs in
`web-contracts-pr` only when the trusted plan does not require
`playwright-functional`; otherwise the job deselects it and the functional shards
run the same cases. The `dev` `browser` job always deselects the wrapper, because
every `dev` push has functional shard evidence, run or inherited.

Web unit/browser tests and their helpers route to the Web, session, Gallery, or
LOSAT capability they exercise. Adding a normal Web regression therefore does not
force a full PR route. Shared fixtures and unclassified tests remain conservative
because they can affect several suites. Unknown production roots also remain full;
recognizing known subsystem owners does not make arbitrary future paths safe.
Renames/copies classify both endpoints, deletions classify their old path, and
malformed/empty diffs or missing Git objects fail closed. The one exception is a
`dev` or Gallery push whose tree is proven identical to its evidenced parent
([Light-change inheritance](#light-change-inheritance)).

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
A `dev` push that leaves the tree unchanged, or changes only metadata,
documentation, policy documentation, and leaf tests, inherits its direct parent's
evidence and runs only the jobs it needs
([Light-change inheritance](#light-change-inheritance)). A `dev` push that merges
a pull request whose green run tested the same tree inherits the jobs that run
ran with the dev-tier command
([Merged pull request inheritance](#merged-pull-request-inheritance)).

A later `dev` push does not cancel a running `Tests` or `Gallery publication`
run. For each workflow, one run executes, the newest later push waits, and
GitHub cancels an older waiting run when a newer push replaces it. Every started
push run therefore finishes with exact-SHA evidence. A SHA whose waiting run was
replaced has no evidence of its own; the newer run tests a tree that contains
it, and a light push on top of it inherits nothing. Pull request runs and manual
dispatches still cancel the run they supersede.

Fail-fast on pull requests: when `Core PR`, `Web contracts`,
or `Recipes standard` fails in a pull request run, a small
dependent job, `Fail fast (<job>)`, runs `gh run cancel` on that run, so the
`Playwright functional` shards and other jobs still running stop. Its condition
is `if: failure() && github.event_name == 'pull_request'`. The failed job keeps its
`failure` conclusion; on a green run the cancel jobs are skipped. Only the
cancel jobs get `actions: write`; the test jobs, which run third-party
packages, keep the read-only workflow token. A `dev` push or a manual dispatch
is never cancelled and still leaves its exact-SHA evidence. A shard cancelled
this way loses no evidence: the commit has already failed, and `PR / gate`,
which does not depend on the cancel jobs, still runs and fails. The three jobs
are the ones whose failures came first in the measured pull request runs; the
weekly CI inventory shows when another job should join them. `Web change budget`
is excluded: the topology contract in `tests/web/architecture-contracts.test.mjs`
keeps every job independent of it.
`tests/ci/ci-impact-cli.test.mjs` pins the cancel jobs and their permissions.

### Light-change inheritance

A `dev` push inherits evidence only from its direct parent
(`github.event.before`), and only when that parent has its own successful
exact-SHA aggregate for the same workflow: `Dev staging / gate` for `Tests`, and
`Gallery readiness / gate` for `Gallery publication`. With that evidence, three
kinds of push run only these jobs and inherit every other job:

| Push | Basis | `Tests` jobs | `Gallery publication` jobs |
| --- | --- | --- | --- |
| Identical tree: `before^{tree}` equals `HEAD^{tree}` | `IDENTICAL_TREE_WITH_DIRECT_PARENT_EVIDENCE` | None | None |
| Only metadata, documentation, and policy documentation | `LIGHT_CHANGE_WITH_DIRECT_PARENT_EVIDENCE` | `recipes-standard` for documentation, `web-change-budget` for policy documentation, none for metadata | None |
| Only those paths and leaf tests, with at least one leaf test | `LEAF_TEST_CHANGE_WITH_DIRECT_PARENT_EVIDENCE` | The jobs of the row above for the non-test paths, plus the [jobs that run the changed leaf tests](#jobs-for-changed-leaf-tests) | `browser` when it runs a changed leaf test; otherwise none |

A leaf-test plan lists the changed leaf tests and the jobs they select. Any
other changed path, such as runtime, a Gallery input, packaging, the control
plane, shared test code, or an unknown path, keeps the complete tier.

The fail-closed rules stay:

- Without successful direct-parent evidence (missing, failed, cancelled, or
  unfinished), `Tests` runs the complete dev tier with basis
  `INHERITED_EVIDENCE_UNAVAILABLE` instead of failing its plan. `Gallery
  publication` fails the plan of a push that changes only documentation and
  policy documentation, and runs its complete profile for every other light
  push, including identical trees and leaf-test changes.
- The planner proves an identical tree by comparing the two tree object IDs.
  An empty diff without that proof, a malformed diff, and a missing Git object
  select the complete tier.
- Manual dispatches and the release profile always run their complete tier
  (`MANUAL_FULL_RUN`). Pull requests keep the routes above; the
  `architecture-change` label still forces their complete tier.

#### Leaf tests

A leaf test is a test module or spec that no other file uses. A changed path is
a leaf test when all of these hold:

1. It matches `tests/test_*.py`, `tests/web/**/*.test.mjs`,
   `tests/web/**/*.playwright.spec.js`, or `tests/web/contracts/**/*.spec.js`.
2. It is not a control-plane test that is already `ci-only`
   (`tests/web/architecture*.test.mjs`, `tests/web/product-impact*.test.mjs`,
   `tests/web/promotion-readiness*.test.mjs`).
3. No other tracked file names it, by file name or, for Python, by module name
   (`tests.test_<name>`). The check reads the head tree. It ignores metadata,
   documentation, and policy-documentation paths and the runner lists
   `tests/ci/**`, `.github/workflows/**`, `playwright.config.js`,
   `tests/web/playwright/*.config.js`, and `package.json`.

Condition 3 makes a test file that another file imports, runs, or reads a
helper, whatever its name. On `dc255c4a`, 361 of the 404 files that meet
conditions 1 and 2 are leaf tests. Examples of the other 43:
`tests/test_protein_s078.py` imports `tests/test_collinearity_units.py`;
`tests/test_api_session.py` and `tools/measure_gallery_publication_performance.py`
run `tests/web/session-request.test.mjs`; `tools/web-product-impact-map.json`
names mapped behavior contracts such as `tests/web/feature-catalog.test.mjs`; and
`tools/check-web-change-budget.mjs` names `tests/web/web-promotion-context.test.mjs`.

Shared test code is never a leaf test: `tests/conftest.py`, `tests/utils/**`,
`tests/fixtures/**`, `tests/test_inputs/**`, `tests/reference_outputs/**`,
`tests/web/helpers/**`, `tests/web/fixtures/**`, `*.cjs` and other support
modules, runners such as `tests/run_*.py`, and Playwright configurations. These
keep their current routes, which select the complete tier on `dev`. So does a
file that meets the three conditions but no row of the job table below, such as
an unlisted `tests/web/contracts/*.serial.spec.js`, because no job runs it.

An added, deleted, or renamed leaf test (both endpoints of a rename) is a leaf
change only when no runner list changes with it. pytest collects
`tests/test_*.py`, the `browser` job finds `tests/web/*.test.mjs`, and the
performance configuration matches `*performance.playwright.spec.js`, so these
files can be added or removed as leaf changes. A changed functional spec is a
leaf only when it exists in the head and `tests/ci/functional-shards.json`
measures it, so adding, deleting, or renaming a functional spec keeps the
complete tier until a refresh measures it. A Playwright configuration names the
Gallery parity and Vibrio specs; changing it is a `ci-only` change. A deleted
leaf test selects the job that ran it, so a runner list that still names it
fails in that job.

#### Jobs for changed leaf tests

| Changed leaf test | Runs in | `Tests` job on `dev` |
| --- | --- | --- |
| `tests/web/**/*.test.mjs` | `node --test` steps of `browser` | `browser` |
| A spec measured in `tests/ci/functional-shards.json` | Functional Playwright shards | `playwright-functional`, restricted to the changed specs |
| `tests/web/*performance.playwright.spec.js` | `npm run test:web:perf-smoke` | `playwright-performance` |
| `tests/web/contracts/gallery-publication-parity.serial.spec.js` | Gallery publication `browser`, pull request `gallery`, `deploy_web.yml` | None; Gallery publication runs `browser` |
| `tests/web/contracts/vibrio-full-generation.serial.spec.js` | Release-only `vibrio-generate-release`, `deploy_web.yml` | None; the plan records it as release-only |
| `tests/test_*.py` | pytest jobs, selected by marker | `core`, plus `recipes-standard`, `gallery`, and `browser` for each of the `recipe`, `gallery`, and `browser` markers named in the file; all four when the markers cannot be read from the file |

Each of the eight shards runs `node tools/ci-impact.mjs shard-files <shard>`
before its setup. It reads the run's planned lists
([Functional shards](#functional-shards)) and, for a restricted functional run,
intersects its list with the changed specs and passes only those files to
Playwright. A shard with no changed spec skips its setup and succeeds, so the
shard matrix and its aggregate stay unchanged. Every other plan runs the whole
shard list. A changed spec in no planned list fails the shard instead of
running nowhere. On a pull request the shard runs the trusted base's
`tools/ci-impact.mjs`, the helper whose policy made the plan, as `PR / gate`
does. The plan then matches the policy that checks it, even when the pull
request adds a job or changes the plan schema.

`Gallery publication` runs `browser` (`npm run test:web:gallery-publication`)
only when a changed leaf test is a spec it runs, which today is only the Gallery
parity spec. Its `performance` job runs
`tools/measure_gallery_publication_performance.py`, not a test file, so a
leaf-test change always inherits it.

Tests marked `slow` run only in the release tier, except the packaging checks
that `browser` runs from `tests/test_web_packaging.py`, so a leaf-test push does
not run them on `dev`.

### Merged pull request inheritance

With merge commits and no up-to-date requirement on `dev`, a pull request run
often tests the exact tree that its merge commit later pushes to `dev`: in 26
of 59 merges between 2026-10-11 and 2026-10-14, the last green run started
after `dev` reached the merge commit's first parent. The `dev` run then inherits
that run's evidence for the jobs both tiers run with the same command:

| `Tests` job on `dev` | Why the pull request job is the same |
| --- | --- |
| `recipes-standard`, `lint` | Same steps on both events |
| `gallery` | Same pytest step; the pull request job may add the Gallery parity steps |
| `playwright-functional` | The same tree gives the same eight planned lists; a pull request plan has no leaf tests, so each shard runs its whole list |

The other dev jobs always run: `core` (Python 3.10 to 3.14; pull requests
run `core-pr` on 3.10, and `core` only for version-dependent changes),
`browser`, `playwright-performance` (a pull request that runs functional
Playwright also runs it, but `dev` keeps its own run),
`losat-cache-browser-acceptance`, and `web-change-budget` (it checks the
`dev` push range in its integrated mode, not the pull request range). A job
is inherited only when the pull request capabilities of the change select it,
so the plan's required jobs follow from its capabilities
(`pullRequestEvidenceJobsFor` in `tools/ci-impact-policy.mjs`). The basis is
`IDENTICAL_TREE_WITH_PR_EVIDENCE`; its inherited evidence names the pull
request run, its `PR / gate` job, and the pull request head.

The `CI impact plan` job of every pull request run uploads an artifact named
`ci-impact-pr-plan.<tree>` followed by `.<job ID>` for each planned job, where
`<tree>` is the tree of the merge commit the run tested. The `dev` planner reads it through
the Actions REST API with `actions: read`, and inherits only when all of these
hold:

1. The `dev` profile classifies the push as a full change with no control-plane
   (`ci-only`) path, and the change selects at least one inheritable job.
   Lighter pushes keep the direct-parent routes above.
2. The push is one merge commit whose first parent is `github.event.before`.
   Its second parent is the pull request head.
3. The newest `Tests` `pull_request` run for that head comes from the same
   repository, completed successfully, and its current attempt has a
   successful `PR / gate`
   (`verifyWorkflowEvidence` in `tools/check-promotion-readiness.mjs`).
4. The run has exactly one record for the merge commit's tree, and that record
   lists every inherited job.

Otherwise the plan stays `FULL_CHANGE` and its summary names the missing
evidence (`NOT_A_PULL_REQUEST_MERGE`, `TESTED_TREE_NOT_RECORDED`,
`PULL_REQUEST_JOBS_NOT_PLANNED`, or the verifier's code). This covers a fork,
a later head push that was not tested, a run that failed or is still running,
a re-run that is in progress or planned differently, a tree that changed because
`dev` moved, and unreadable evidence. `Gallery publication` keeps its own
routes: a pull request run does not execute its jobs.

### Carrying evidence to a later commit

`node tools/ci-impact.mjs classify --base <E> --head <H>` reads
`git diff --name-status -z --find-renames <E> <H>` from the local repository. It
queries no network service and no workflow evidence. It prints each changed path
with its capability, reason, and leaf-test jobs, then the capability union and
three verdicts. A verdict is true only when E is an ancestor of H
(`git merge-base --is-ancestor`) and every changed path is in the verdict's set
below. An identical tree meets every verdict; an invalid diff meets none. Test
paths count only as leaf tests; their subject capability, such as `web-runtime`
for a Web test, does not matter.

| Verdict | Changed paths allowed | Evidence it carries from E to H |
| --- | --- | --- |
| `releaseEvidenceCarries` | Metadata; documentation and policy documentation outside `docs/recipes/` and `docs/capture/`; leaf tests whose tests run only in jobs the dev tier also runs on the same Python versions: `.test.mjs` files, functional and performance specs, the Gallery parity spec (Gallery publication runs it on H), and `tests/test_*.py` files that name none of the `slow`, `recipe`, `gallery`, and `browser` markers | The `Tests` release-tier run |
| `generatedArtifactChecksCarry` | Metadata; documentation and policy documentation outside `docs/images/`, `docs/capture/`, and `docs/recipes/`; `ci-only` paths; leaf tests | The Gallery refresh, Gallery tutorial media, and docs GUI capture checks, and the Gallery artifact manifest |
| `localTestEvidenceCarries` | Metadata; documentation and policy documentation outside `docs/recipes/`; leaf tests, except `tests/test_*.py` files with a `slow` marker or unreadable markers | Recipe runs (`run_cli_scenarios.py` and `run_python_scenarios.py` with `--all --check`), `TestOutputComparison`, the `tools/audit/` sweeps, and the `main`-written Session fixture tests |

The release tier runs the Vibrio full-generation spec (`vibrio-generate-release`),
the render-stage universe (`render-stage-universe-release`; its test files spell
their `slow` mark, so a change to them carries no release evidence),
tests marked `slow` (`slow-main`), tests marked `recipe`, `gallery`, or
`browser` on Python 3.10 and 3.12 (`acceptance-supported-main`), and the
promotion audit jobs (`promotion-audit-*`: the recipe runners under
`docs/recipes/`, the `tools/audit/` sweeps and user journeys, and the random
walk). The dev tier repeats none of these, so changing them does not carry
release evidence. When
the markers of a Python test cannot be read from the file, the verdict is false.
A changed leaf test that is itself one of the local checks reruns on H in the
dev tier, unless it is marked `slow`.
Shared test code feeds the generators, the sweeps, and the local checks, so a
change to it meets no verdict. `docs/recipes/` and `docs/capture/` hold
executable scenario runners and expected outputs that recipe and capture tests
read, so the verdicts whose evidence reads them exclude them. The periodic audit
uses these verdicts for promotion evidence
([carry-forward rule](WEB_PERIODIC_AUDIT.md#carrying-evidence-forward)), and the
release workflow uses `releaseEvidenceCarries` before PyPI publication
([Release evidence check](WEB_PERIODIC_AUDIT.md#release-evidence-check)).

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

### Functional shards

Full functional CI runs eight shards with two retries and two workers per shard,
on pull requests and in `dev` staging alike. Playwright's `--shard` splits by
case count in spec file order; it cannot split by duration. With equal counts,
four hosted shards took 13 to 34 minutes on `52be34a7` and `a5deed4a`, because
the heaviest composite and Session cases sit in a few files. Each matrix entry
therefore runs its own list of spec files, passed as Playwright CLI file
filters. `tests/ci/functional-shards.json` records the shard count and the
measured seconds per case of each spec file; it holds no lists. The `CI impact
plan` job plans the lists once per run, after `npm ci`, with
`node tools/balance-functional-shards.mjs plan`: it lists the head's functional
cases, weighs each spec file by its current case count times its seconds per
case (a file without a measurement by the mean of the measured cases), and
fills the lists longest file first into the shard with the least work. Within a
shard, the two workers split the cases. The lists are the job's
`functional-shards` output, so all eight shard jobs of a run, and a rerun of
one of them, read the same lists. A new spec file or a spec file that gains
cases therefore moves to a lighter shard on its own; no pull request edits the
measurements for it.

The lists are head data, as the spec files are: the head's tool plans them on
pull requests too, and the trusted base's `shard-files` only applies the plan's
leaf-test filter to them (see [Trust and aggregate checks](#trust-and-aggregate-checks)).

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

Since the diagram Worker boots from a Pyodide snapshot, the lists were rebuilt from
the green `dev` runs 38944092373 and 38946503427 (all eight shards executed). The
recorded minutes fell by about 55 percent in total, the largest file
`live-generate-parity-legend.playwright.spec.js` from 17.4 to 7.1 minutes. The
previous lists would have given a slowest shard of 10.6 minutes against a mean of
9.1 (ratio 1.18); the new lists give 9.1 minutes in every shard.

Refresh the measurements when the weekly inventory asks for it (below), that is,
when the seconds per case have moved or spec files have none. Download the shard
artifacts of two green `dev` runs with
`gh run download <run-id> -p 'playwright-functional-shard-*' -D <dir>`, then run
`node tools/balance-functional-shards.mjs refresh <dir>/*/*/functional-report.json`,
which rewrites the measurements from the median of each case and prints the plan
they give. `tests/ci/playwright-inventory.test.mjs` plans the lists the same way
and lists every shard with its file filters to prove that the matrix partitions
the full suite without omissions or duplicate execution.

The job limit is 45 minutes; case timeouts are unchanged. Hosted
shard 1 still exceeded 30 minutes with two workers on `9dd6359c`: its composite
placement took 7.2 minutes, a composite Session round trip took 4.4 minutes,
and a Definition completion wait and retries consumed the remaining budget.
The Definition test now triggers an actual edit; the job budget also allows
for the measured composite workload, setup and existing retries. The line
reporter names the running case, and the JSON report records durations,
retries and final results alongside traces on every completed job. A job killed
at its limit may not reach upload; its last running case remains in the log.

## Pyodide snapshot

A Pyodide memory snapshot lets the diagram Worker skip its package installs in
functional tests. Each `playwright-functional` shard builds one after `Prepare
browser wheel` into `$RUNNER_TEMP/pyodide-snapshot` and passes that directory
to the tests; a build failure fails the job. Other jobs build none.
`node tools/build_pyodide_snapshot.mjs <dir>` boots Pyodide
in Node, runs the Worker's package recipe (`installRuntimePackages`: Pyodide
loads micropip and tzdata from its distribution, then one micropip call
installs the dependency wheels and the gbdraw wheel),
runs the Worker helpers, checks that every installed distribution's
requirements are installed, and writes `snapshot.bin` (the Wasm memory),
`state.json` (the side-module linker state and the validity key) and `fs.bin`
(the site-packages files), about 170 MB, in 5 to 10 seconds. It needs a Node
with WebAssembly JSPI (`WebAssembly.Suspending`), as Chromium has, such as the
version in `.node-version`: without JSPI Pyodide's function table starts one
entry shorter, Chromium cannot restore the snapshot, and the builder refuses
to run. With
`GBDRAW_TEST_PYODIDE_SNAPSHOT_DIR` set, `playwright.config.js` starts its
test server `tools/serve_web_tests.py` with the snapshot, so it also serves
that directory at `/__pyodide-snapshot__/`
(`tests/web/helpers/pyodide-snapshot-test-server.cjs`) and `config.js` with
`PYODIDE_SNAPSHOT_URL` pointing at it. The Worker then restores the snapshot
(about 0.5 seconds) instead of installing its packages (about 5 seconds
locally), runs `_start_web_worker_state()` so each Worker has its own caches
and cache epoch, and calls numpy and pandas once. Every restored Worker
shares the string-hash secret of the build, so functional CI is less likely
to expose a bug that depends on set or dict order of strings. On draft PR #1115, building it in each functional shard
cut functional test time from 297.8 to 139.8 test-minutes. The repository's
`config.js` keeps `PYODIDE_SNAPSHOT_URL = null`: users, smoke, browser and
release acceptance jobs boot cold.

The Worker boots cold whenever a restore fails: a missing file or one that
takes longer than 60 seconds to download, another recipe format, wheel or
helper source (the key in `state.json` covers the format version and the bytes
of every wheel and of the helper source), another Pyodide build (its
`BUILD_ID`, which Pyodide checks), a table slot or JS reference the builder
could not name, or an error while relinking, starting the Worker state or
calling numpy and pandas. The page
then discards that Worker, starts one without the snapshot, and boots every
later Worker cold. The page reports each boot as the `workerBoot` structural
metric: `snapshot`, `cold`, or `cold after snapshot error: <reason>`.
`tests/web/pyodide-snapshot-boot.playwright.spec.js` requires `snapshot` and a
working fallback, so a snapshot that stops restoring fails CI instead of
silently costing time. A case that asserts the cold boot (the wheel request)
calls `forceColdPyodideBoot(page)`, which aborts the snapshot requests.

To run specs locally with a snapshot, prepare the wheel, build the snapshot
outside the repository, and set the variable with a free port: Playwright
reuses a server already listening on the port locally, and a plain server
serves no snapshot.

```bash
python tools/prepare_browser_wheel.py
node tools/build_pyodide_snapshot.mjs /tmp/pyodide-snapshot
GBDRAW_TEST_PYODIDE_SNAPSHOT_DIR=/tmp/pyodide-snapshot GBDRAW_WEB_TEST_PORT=<free port> \
  npx playwright test --config=tests/web/playwright/functional.config.js <spec>
```

Rebuild the snapshot after every wheel or helper change; a stale one only makes
the Worker boot cold. Without the variable, the server and the Worker are
unchanged.

The restore depends on Pyodide and Emscripten internals
(`gbdraw/web/js/workers/pyodide-snapshot.js`): `makeMemorySnapshot` saves only
the Wasm memory, so the restore re-instantiates the side modules at the memory
and table bases of their dlopen handles (`struct dso` offsets) and rebuilds the
table slots, the GOT and the linker tables through the internals Pyodide
exports on its module (`LDSO`, `GOT`, `GOTHandler`, `getDylinkMetadata`,
`wasmTable`, `newDSO`, `mergeLibSymbols`, `setWasmTableEntry`,
`updateTableMap`). A Pyodide or vendored-package upgrade must rebuild the
snapshot and pass `pyodide-snapshot-boot.playwright.spec.js` before it merges,
and a builder failure must fail the job that builds it.

## Weekly CI inventory

The `CI inventory` workflow (`.github/workflows/ci-inventory.yml`) runs
`tools/ci-inventory.mjs` every Monday at 00:17 UTC and on manual dispatch. It
checks out `dev`, needs no session, and never runs on a pull request. GitHub
runs scheduled workflows from the default branch, so the schedule starts once
the workflow reaches `main`. Its report covers the last 7 days, except where an
item names its own window:

- runner-minutes by workflow and event, and by job name (matrix jobs together);
- `Tests` wall time (run start to its last job) and queue time (job created to
  job started), p50 and p90, for pull request and `dev` push runs;
- pull request feedback time over the last 20 pull requests merged into `dev`,
  whenever they merged: the median (nearest rank) for those whose run had a
  functional Playwright shard and for the others, next to the targets in
  [Validation tiers](#validation-tiers);
- the functional shards of the latest green `dev` push `Tests` run that kept
  every shard report: test-minutes per shard, the slowest shard against the
  mean, the case count, the 15 spec files with the most test-minutes, and the
  plan that measurements refreshed from that run would give;
- pull request CI, from the pull requests merged into `dev` in the window and
  the pull request `Tests` runs created in it (below);
- changes since the previous week in functional cases, functional
  test-minutes, runner-minutes, jobs in `test.yml`, and the four pull request
  CI numbers.

Pull request CI has four numbers:

- **Merged pull requests that passed on one run**: the share of the pull
  requests merged into `dev` in the window, among those with a `Tests` run,
  that had exactly one, and that run passed at its first attempt. A pull
  request's runs are the pull request `Tests` runs on its head branch (of its
  head repository) created between its opening and its merge, give or take
  five minutes, including runs before the window.
- **`Tests` runs per merged pull request**: the mean number of those runs.
- **Pull request runs replaced by a newer run or cancelled by hand**: the
  share of the pull request `Tests` runs created in the window that were
  cancelled without a failed job. Most were replaced by a newer run of the
  same pull request: `test.yml` also runs on `labeled` and `unlabeled`, and
  its concurrency group cancels the running run for a push and for a label
  change alike.
- **Pull request runs that failed**: the share of those runs that failed,
  timed out, or did not start, together with the runs that a `Fail fast (...)`
  job cancelled after a failed job.

GitHub reports both kinds of cancelled run as `cancelled`. When the report has
read a cancelled run's job list, the run failed if a job of its attempt other
than `PR / gate` failed, and was replaced or cancelled by hand otherwise. A
cancelled run without a job list counts as replaced when the next run on its
branch was created at most 120 seconds after it ended (or before it ended),
as failed when the gap is longer, and as replaced when no later run exists.
The other conclusions (such as `skipped`) count in the total but in neither
share.

The report goes to the job summary and to one comment a week on the open issue
labelled `ci-inventory`; the workflow creates the label and the issue when they
are missing, and a second run in the same week replaces that week's comment.
Each comment ends with a data block that the next week's run reads back; only
the workflow's own comments (`github-actions[bot]`) count, and a block with a
malformed field is ignored. The workflow's token is the first in this
repository with `issues: write`. Its 1,000 API requests an hour belong to the
whole repository, so every workflow shares them, `CI impact plan`, `PR / gate`
and the `Fail fast (...)` jobs among them. The tool therefore reads job lists
only for `Tests` and `Gallery publication` runs, the cancelled pull request `Tests` runs first (so the four pull request numbers use job lists), then the other `Tests` runs, and stops
reading them while 500 requests of the hour are left (`RATE_RESERVE`); after
its other reads (about 150), a scheduled run reads the job lists of at most
about 350 runs, fewer when other workflows use the hour, against about 750
before the reserve rose from 100; a busy week has about 1,200 `Tests` runs. Other runs, and the runs left when the job
lists stop, count by their own duration in rows marked "(run duration)".
Three figures come from the job lists read: the
runner-minutes by job (and the job-level split of runner-minutes by
workflow), the `Tests` wall and queue time, and the decision for a cancelled
pull request run (above). The runner-minutes by job are a sample of the budget left after the cancelled pull request runs; the decision covers all of them (about 190 a week fit in the 350). When the API refuses the run
listing (HTTP 403 or 429), the report covers the runs read so far and says it
is partial. A filtered run listing returns at most 1,000 runs, so the tool
reads the window a day at a time and a busier day in halves. Every read,
the artifact downloads included, is tried up to three times when GitHub
answers with a server error (HTTP 500 and up) or a body cut short; a write is
never repeated.

CI plans the shards from the recorded seconds per case and the current case
counts, so new spec files and new cases never unbalance them. A report asks for
fresh measurements when the slowest shard has more than 1.2 times the mean
test-minutes (the seconds per case have moved), or when a functional spec file
has no measurement in `tests/ci/functional-shards.json` (it is planned at the
mean case). A session that reads a flagged report:

1. Refreshes the measurements with the reports of two green `dev` runs
   ([Functional shards](#functional-shards)).
2. Reviews the spec files and jobs that grew most: merges cases that repeat a
   journey another case already proves, moves checks that need no browser to
   node or pytest, and removes cases and jobs whose behaviour is covered
   elsewhere. Removing an acceptance case follows "Smoke inventory and
   regression retention".

A report warns when either pull request feedback median is above its target.
The warning never fails the workflow. A session that reads it names the jobs
that set the wall time of the slow runs and asks the Owner whether to plan
work that shortens them.

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
The promotion audit jobs are not release-profile jobs and no gate reads them;
`Promotion audit / summary` tabulates their findings for the promotion
checklist ([`WEB_PERIODIC_AUDIT.md`](WEB_PERIODIC_AUDIT.md)).
A release dispatch cannot emit ordinary `Dev staging / gate` success.

Plans use schema 3 and include the capability union. A leaf-test plan also lists
each changed leaf test with its kind and jobs (`leafTests`; `null` in every other
plan), and an identical-tree plan has the single capability `none` and no changed
path. A merged pull request plan names a pull request run as its evidence. The
active helper validates
its own exact schema, profile, workflow SHA, required job list, and inherited
evidence. It rejects missing/skipped/failed/cancelled required jobs and unexpected
failures in optional/additional jobs. Matrix aggregates must succeed. Workflows
start on all relevant events and filter at job level; they never use workflow
`paths` filters that could omit an aggregate.

### Dev state on pull requests

dev is red when the newest `Tests` run of a push to `dev` that has a result
failed (`failure`, `timed_out`, or `startup_failure`). A run's result is that
of its newest attempt that finished with one. An attempt still running, or one
that ended `cancelled`, `skipped`, `neutral`, `stale`, or `action_required`,
does not count, and the attempt before it decides; a run with no such attempt
is skipped. A failed run therefore stays red until an attempt of it passes,
also while a rerun runs or after a rerun is cancelled. `PR / gate` reads this through the Actions REST API with
`actions: read` (`readDevTestsState` in `tools/check-promotion-readiness.mjs`,
called by the `gate` command of the trusted `tools/ci-impact.mjs`). While dev is
red it fails with `DEV_RED`, naming the run and its head SHA, unless the pull
request carries the label `dev-red-fix`; `AGENTS.md` ("When dev is red") limits
that label to the pull request that fixes the failure. The gate fails closed
with `DEV_STATE_UNAVAILABLE` and the reader's code when the run list or an
attempt cannot be read or does not match the run (`API_REQUEST_FAILED`,
`RUN_IDENTITY_MISMATCH`, and the other verifier codes), when no run among the
newest 100 has a result (`NO_DEV_STATE`), and when any listed run or read
attempt has an unknown conclusion (`UNKNOWN_RUN_CONCLUSION`). The step summary shows the dev
state and the label either way. The `dev`, release, and Gallery gates do not
read it.

The gate reads dev once, when it runs. A pull request held this way passes
after dev turns green by rerunning its `PR / gate` job (about a minute); adding
or removing a label starts a new run. The gate does not hold a pull request
whose gate passed before dev turned red, nor one whose `pull_request.base.sha`
predates this check because its head or base was not updated since, including
reruns of older runs: such a run uses that base's helper. A new run started by
`labeled`, `unlabeled`, or `reopened` may keep the old base SHA too. Updating
the pull request from `dev` (`gh pr update-branch`) gives it a run with the
current helper; the rule text still asks for a check before merging.

The check took effect in two steps, because a checker changes separately from
the workflow that feeds it ([`WEB_CHANGE_POLICY.md`](WEB_CHANGE_POLICY.md)).
The workflow change passed the repository, labels, and token to `PR / gate`
first; the base helper ignored them. The pull request that added the check ran
that base helper too, so it was not held by the red `dev` it was written for.

### Repository growth on pull requests

The step summary of `PR / gate` shows how many MB of new Git objects merging
the pull request adds to `dev`. It is a report and never fails a job. The CI
impact plan job measures it on the merge commit: the size of a thin pack of
the objects reachable from the merge commit and not from `dev`, with objects
already in `dev` as delta bases. This is close to what the merge adds to the
repository's storage. Its full-history checkout has no old file contents, so
the step first fetches, in one request, the blobs the merge adds and `dev`'s
versions of the paths it modifies; Git fetches any other delta base it needs by
itself. The step gives up after 3 minutes, because every job waits for the CI
impact plan job. A Gallery regeneration that rewrites every Session shows as
tens of MB; most pull requests add well under 1 MB.

### Job times against their limits

`PR / gate` and `Dev staging / gate` report the jobs of their own run attempt
that used 80% or more of their `timeout-minutes`. They report before they
validate the results, so a job that ran out of time is reported too.
`tools/ci-impact-job-times.mjs` reads each finished job's start and end through
the Actions REST API (`actions: read`) and its limit from the run's own
`test.yml`; a job without `timeout-minutes` has GitHub's 360 minutes.
`PR / gate` checks out that one file of the tested merge commit and reads it as
data. Each such job gets a warning annotation and a line in the step summary;
when there is none, the summary names the job with the highest share. The
report never fails a gate, and a failed read is a warning. A job near its limit
is the time to shorten or split it, or to raise its limit with measurements,
before it fails. Over the 12 green `dev` runs before the report, Core 3.12 took
up to 19.4 of its 25 minutes (77%); run 42146106065 had Web change budget at
7.7 of its 10 minutes (77%).

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
