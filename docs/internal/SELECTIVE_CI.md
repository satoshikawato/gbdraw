# CI validation tiers and impact routing

`tools/ci-impact-policy.mjs` owns path classification and required job registries.
`tools/ci-impact.mjs` reads Git changes, verifies inherited evidence, and validates
aggregate results. Workflows execute the selected jobs; they do not duplicate the
classification rules.

## Validation tiers

| Tier | Trigger | Required coverage |
| --- | --- | --- |
| PR | Every PR into `dev` | Changed subsystem, cross-layer smoke, architecture/Product policy; Python 3.11 primary |
| Integrated dev | Every push to `dev`; `Tests` dispatch with `tier=dev` | Core on Python 3.10/3.11/3.12; recipes, Gallery, browser/package integration, offline GUI, LOSAT cache, all functional Playwright in four shards, performance smoke on 3.11 |
| Release / S11 | Explicit `Tests` dispatch on `dev` with `tier=release` | Every dev functional job, additional recipe/Gallery/browser acceptance on 3.10/3.12, exhaustive non-browser slow tests on all three versions, exact-candidate Gallery readiness |

PR feedback targets 5–8 minutes. Integrated functional validation targets 15–20
minutes where practical. S11 remains intentionally expensive. The separate
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

## Changed paths → capabilities → required jobs

Plans contain the ordered union of affected `capabilities`. `impact` is the
highest-ranked display label; it does not discard other capabilities. The job
registry unions their contributions, then emits unique jobs in registry order.

| Capability | Representative paths | Selective PR jobs |
| --- | --- | --- |
| `metadata` | Existing development-tool allowlist, `.gitignore`, citation/licenses | No test job |
| `documentation` | Root Markdown and `docs/`, excluding policy authority | Recipes |
| `python-core` | Known Python core/API/config/I/O owners and data | Architecture, core, lint, Web contracts, browser smoke |
| `renderer` | Render/SVG/diagram/canvas/features/labels/layout owners | Architecture, core, lint, Web contracts, browser smoke |
| `web-runtime` | `gbdraw/web/index.html`, first-party JS | Architecture, Web contracts, browser smoke |
| `session-persistence` | Python session codecs; JS session/config/history owners | Architecture, core, recipes, lint, Web contracts, browser smoke |
| `gallery` | Gallery inputs and their generator owners | Architecture, Gallery, Web contracts, browser smoke |
| `losat-integration` | Comparison owners, LOSAT JS/workers/Wasm | Architecture, core, lint, Web contracts, browser smoke |
| `tests-only` | Shared fixtures/harness and unclassified tests | Full PR tier |
| `packaging` | Dependencies, manifests, vendored runtime, packaging tools | Full PR tier |
| `ci-only` | Workflows, planner/tests, Playwright configuration, architecture/Product authority | Full PR tier |
| `full` | Unknown or invalid paths | Full PR tier |

The full PR job IDs are `web-change-budget`, `core-pr`, `recipes-standard`,
`gallery`, `lint`, `web-contracts-pr`, and `web-pr-smoke`.

`web-contracts-pr` groups the existing fast JS contracts and non-slow Python
browser suite, which together took approximately 65–70 seconds historically.
`web-pr-smoke` runs only the ten selected Playwright cases. They execute in
parallel jobs instead of sharing one 15-minute serialized budget. Both jobs keep
independent working directories, dependency installs, wheels, and browser state.
The contracts job also runs when the trusted plan requires `web-pr-smoke`, so the
pre-redesign base planner still requires all three original suites during rollout.

Web unit/browser tests and their helpers route to the Web, session, Gallery, or
LOSAT capability they exercise. Adding a normal Web regression therefore does not
force a full PR route. Shared fixtures and unclassified tests remain conservative
because they can affect several suites. Unknown production roots also remain full;
recognizing known subsystem owners does not make arbitrary future paths safe.
Renames/copies classify both endpoints, deletions classify their old path, and
malformed/empty diffs or missing Git objects fail closed.

### Evidence required for selection

A narrower PR route requires successful exact-SHA `Dev staging / gate` evidence
at the PR base. No older successful SHA can substitute. API failures, missing or
unfinished runs, cancellation, and failed/malformed evidence select the complete
PR tier. `architecture-change`, control-plane changes, dependencies, unknown
paths, unclassified/shared test inputs, and explicit dispatches also require the full relevant tier.

Every runtime/subsystem change runs complete integrated dev and Gallery coverage.
Only metadata/documentation changes can inherit direct-parent staging evidence:
metadata runs no test jobs, documentation runs recipes in `Tests`, and neither
reruns Gallery publication when exact parent Gallery readiness is successful.
Consecutive pushes whose direct parent is unfinished/cancelled run the full tier.

## Smoke inventory and regression retention

Smoke retains app/Worker boot and Circular generation, Linear multi-record
rendering, mounted Feature/Label/Legend edits, export startup, Save/Load/Generate,
Circular and Linear mode transitions, invalid-annotation preservation, and invalid
composite-session rejection. Full functional discovery includes every smoke case.

`tests/ci/playwright-inventory.test.mjs` uses Playwright's expanded `--list` output.
It counts nested and parameterized cases and rejects fewer than 8 or more than 12;
counting tag strings in source missed the earlier growth to 27 cases. No assertions,
case bodies, test timeouts, or full-functional exclusions changed in this redesign.
The exact relocated cases and before/after counts are in the
[measurement report](CI_TIER_REDESIGN_2026-09-11.md).

The smoke config retains one worker and zero retries. Three local two-worker
repetitions passed, but these are not repeated GitHub runner evidence. One worker
already meets the projected target, so CI does not adopt unverified parallelism.
Full functional CI retains four shards and its existing retry policy.

## Trust and aggregate checks

PR planning and `PR / gate` execute the PR base's helper under `.ci-trusted-base`.
Candidate helpers are test inputs only. A missing or invalid base helper fails;
there is no candidate fallback. `Web base policy (trusted base)` remains a separate
`pull_request_target` check that treats candidate files as Git data.

The redesign PR therefore receives the old base's complete PR route. New selective
routing can take effect only after reviewed integration into protected `dev`.
Additional contract failures also block the old gate's unknown-job validation.
No branch protection or required external status name changes.

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
prepare the browser wheel and run `npm run test:web:pr-smoke`. Full functional
coverage remains `npm run test:web:functional-full`.

The [timing CSV](CI_TIER_TIMINGS_2026-09-11.csv) records each historical job's setup,
dependency, wheel, test, teardown, and wall-clock times with direct job URLs.
Plans report capabilities, required jobs, evidence links, and fallback reasons.

Rollback this coherent CI change together: policy, adapter, workflow, smoke tags,
and policy/inventory contracts. Preserve trusted-base evaluation and stable gate
names. Restoring exhaustive matrices to dev changes cost, not release criteria.
