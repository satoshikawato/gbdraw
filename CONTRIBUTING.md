# Contributing to gbdraw

Thank you for your interest in gbdraw. gbdraw is a maintainer-led open-source
project. Use, inspection, bug reports, reproducible examples, research and
publication use cases, scientific and domain review, testing in additional
environments, documentation corrections, and focused proposals are welcome.

The maintainer retains final authority over product behavior, scientific and
technical terminology, architecture, scope, compatibility policy, and
releases. Bug reports, reproducible examples, and questions about existing
behavior may be opened directly. Small documentation corrections may also go
directly to a pull request. Anything larger must begin with an issue and wait
for explicit agreement on the intended behavior and whether a pull request is
the appropriate next step.

External pull requests are optional proposals, not the project's default
development path or an entitlement to merge. An accepted idea may be revised,
split, independently implemented, or adopted without merging the submitted
patch. Contributing does not imply roadmap priority, repository access, a
merge, or maintainer status.

## Most useful contributions

The project particularly values:

- reproducible bug reports;
- exact commands or Python examples;
- minimal, shareable inputs that expose an edge case;
- concrete research or publication workflows;
- review of scientific terminology and figure semantics;
- verification on platforms, browsers, and datasets not covered by the
  maintainer;
- focused documentation corrections; and
- small, agreed changes with regression tests.

A feature request is evidence for a product decision. It is not a commitment
to implement the feature or accept a particular design.

## Governance and decision rights

User reports and external review inform the project, but gbdraw is not governed
by contributor voting or consensus. The maintainer decides:

- whether a problem belongs within gbdraw's scope;
- authoritative user-facing behavior;
- scientific and technical terminology;
- public API, CLI, session, and file-format contracts;
- architecture and dependency changes;
- compatibility and deprecation policy;
- release contents and timing; and
- whether a change is maintainable enough to merge.

A technically functioning change may still be declined if it would introduce:

- ambiguous behavior;
- a parallel implementation path;
- duplicated semantic ownership;
- an unjustified compatibility burden;
- excessive maintenance cost;
- an unsupported scientific interpretation; or
- scope drift.

## Before you start

- Search [existing issues](https://github.com/satoshikawato/gbdraw/issues)
  and [pull requests](https://github.com/satoshikawato/gbdraw/pulls) first.
- Bug reports, reproducible examples, and questions about existing behavior
  may be opened without prior approval.
- Small typo, broken-link, and narrowly scoped documentation corrections may
  go directly to a pull request.
- Start anything larger with an issue. Do not begin implementation until the
  maintainer explicitly confirms both the intended behavior and that a pull
  request is an appropriate next step. Agreement that a problem belongs in
  gbdraw does not approve a particular design.
- Obtain prior agreement for changes to rendering geometry, public APIs, CLI
  options, persisted formats, dependencies, architecture, or feature scope.
- Unsolicited large features, broad refactors, generated rewrites, and
  architecture changes may be closed without line-by-line review.
- Read [`CODE_OF_CONDUCT.md`](./CODE_OF_CONDUCT.md). It applies to issues,
  pull requests, and other project interactions.

These boundaries keep reports and evidence easy to share while avoiding
substantial work on an approach the project cannot maintain.

## Reporting a bug

Open a [GitHub issue](https://github.com/satoshikawato/gbdraw/issues/new) with:

- the `gbdraw` version (`gbdraw -h` prints it) and installation method
  (Bioconda, PyPI, or source);
- the affected interface (Web, CLI, or Python) and diagram mode (Circular or
  Linear), when relevant;
- the exact command or Python snippet you ran;
- what you expected and what happened instead;
- the operating system or browser for environment-specific failures;
- whether the problem is reproducible with the current stable release; and
- a minimal input file or small excerpt, when possible.

Do not upload confidential, restricted, or otherwise non-shareable genomic
data. A synthetic or redacted input that preserves the failure is welcome.

## Proposing a change

Begin with the user problem rather than a preferred implementation. Describe:

- the affected interface and diagram mode;
- the current behavior;
- the expected behavior;
- a minimal example;
- why existing options are insufficient;
- possible effects on existing figures, saved sessions, requests, commands,
  or Python code;
- relevant scientific terminology or interpretation; and
- plausible behavioral alternatives when more than one answer is reasonable.

Substantial implementation must wait for maintainer confirmation. Agreement
that a problem exists is not approval of a patch. The maintainer may first ask
for a smaller evidence-producing change or an explicit behavior decision.

Before implementation, maintainers and developers apply the
[Product Impact Ratchet](./docs/internal/PRODUCT_IMPACT_RATCHET.md) to both
registered architecture-subject changes and unmapped material user-effect
ambiguity. Unresolved evidence or Product judgment uses the non-normative
[Product Decision Pack template](./docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md);
the completed working template does not itself authorize an outcome.

## Development setup

```bash
git clone https://github.com/satoshikawato/gbdraw.git
cd gbdraw
python -m pip install -e ".[dev]"
```

Add the optional export dependency if your change touches PNG/PDF/EPS/PS
output:

```bash
python -m pip install -e ".[dev,export]"
```

Non-SVG export also needs system Cairo/Pango libraries; see
[Installation](./docs/INSTALL.md#4-source-installation-for-development) for
platform-specific package names.

## Running the checks CI runs

```bash
# Lint (blocks CI on failure)
ruff check gbdraw/

# Fast test suite (excludes slow/regression-heavy tests)
pytest tests/ -v -m "not slow"

# Full suite, including slow tests (runs in integrated staging and release CI)
pytest tests/ -v
```

Run a narrower slice while iterating:

```bash
pytest tests/test_regression.py -v
pytest tests/ -v -k "test_circular_basic"
pytest tests/ -v -m "circular"      # or -m "linear"
```

For Web changes, run the focused JavaScript tests and the relevant browser
scenario. Browser tests that depend on the packaged application require the
generated browser wheel:

```bash
node --test tests/web/<focused-test>.test.mjs
python tools/prepare_browser_wheel.py
npx playwright test tests/web/<focused-test>.playwright.spec.js
```

The JavaScript Playwright specifications require Node's `@playwright/test`.
When only Python Playwright is available, use an equivalent focused browser
check rather than omitting browser verification.

## Changing rendered SVG output

Many tests compare generated SVGs against tracked references in
`tests/reference_outputs/`. If your change intentionally alters diagram
geometry:

```bash
# 1. Confirm the differences are the ones you intended
pytest tests/test_output_comparison.py::TestOutputComparison -v

# 2. Regenerate references only after reviewing that output
pytest tests/test_output_comparison.py::TestGenerateReferences \
  --update-reference-outputs -v

# 3. Review the diff and re-verify
git diff -- tests/reference_outputs/
pytest tests/test_output_comparison.py::TestOutputComparison -v
```

Do not update reference outputs merely to make a failing test pass. First
confirm that the new geometry is intentional and correct; the reference files
record the project's intended rendering behavior.

## Changing documentation

Use [`docs/DOCS.md`](./docs/DOCS.md) as the stable index for
[Tutorials](./docs/TUTORIALS/README.md),
[Technical documentation](./docs/REFERENCE/README.md),
[FAQ](./docs/FAQ.md), and [Gallery](./docs/GALLERY.md).

Before adding a page, identify which existing document owns the reader's
question. A new interface or test scenario does not automatically require a
new top-level page. Technical documentation owns exact contracts; FAQ pages
answer choices and troubleshooting by linking to those contracts. Internal
plans, audits, and executable evidence registries belong under
`docs/internal/`.

Tutorial commands, public code blocks, and internal evidence scenarios are
checked against the current CLI, Web, and Python interfaces:

```bash
pytest tests/ -v -k documentation
```

## Style

- Python: use type hints throughout (`from __future__ import annotations`),
  `snake_case` for functions and modules, `PascalCase` for classes,
  `UPPER_SNAKE_CASE` for constants, and the `logging` module for diagnostics
  instead of bare `print` calls.
- Extend an existing configurator, drawer, group, or typed boundary instead of
  adding a parallel path. The CLI, Web, and Python interfaces should converge
  on shared internals.
- Do not add speculative options, fallbacks, or abstractions for hypothetical
  future use. Keep changes within the agreed problem.
- Architecture-bearing changes must follow the
  [architecture fitness-function ratchet](./docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)
  and remove superseded owners and paths in the same change.

## Branch roles

- Approved work starts from the latest merged `dev`.
- Work-branch pull requests target `dev` and use the fast implementation
  admission checks.
- `dev` is the integration branch. Heavy supported-version, slow, browser, and
  Gallery checks run after merge against the exact integrated `dev` SHA.
- `main` is the release branch. A maintainer promotion pull request from `dev`
  contains no implementation changes and uses trusted exact-SHA readiness and
  tree evidence.
- A push to `main` performs release verification, build, and deployment against
  the exact release SHA.
- Direct `main` hotfixes are unsupported. Urgent fixes still go through `dev`
  staging and promotion.

See the normative
[Web change, review, and delivery policy](./docs/internal/WEB_CHANGE_POLICY.md)
for lifecycle ownership and trust boundaries.

For a fork, keep your fork as `origin`, add this repository as `upstream`, and
create a non-tracking work branch from `upstream/dev`:

```bash
git remote add upstream https://github.com/satoshikawato/gbdraw.git
git fetch upstream
git switch --no-track -c <branch-name> upstream/dev
```

If this repository is already your `origin`, use:

```bash
git fetch origin
git switch --no-track -c <branch-name> origin/dev
```

## Pull request change classes and review results

Select exactly one change class in the pull request template:

- `STANDARD` for implementation, tests, documentation, refactoring, or cleanup
  that does not need a governance, architecture-exception, or promotion packet;
- `GOVERNANCE` for authority, policy, CI evidence-producer, template, or
  repository-governance changes;
- `ARCHITECTURE_EXCEPTION` only when architecture debt increases, a new
  persisted compatibility path is introduced, multiple owners or paths remain,
  an accepted deterministic violation is added, or a time-bounded hard-invariant
  waiver is requested; and
- `PROMOTION` only for a no-implementation `dev` to `main` release-admission
  pull request.

An ordinary architecture-bearing `STANDARD` change supplies concise owner/path
evidence and removes superseded locations; it does not need an empty OE/PE/CB
packet or a dedicated approval permalink. Full sets, arithmetic, alternatives,
removal conditions, and an explicit maintainer decision are reserved for
`ARCHITECTURE_EXCEPTION`.

Apply the `architecture-change` label to a Web implementation that changes an
architecture owner, canonical path, compatibility path, privileged boundary, or
lifecycle responsibility. The label selects the architecture size-review
profile and routes attention. It does not waive deterministic gates.

Policy reports distinguish `Gate: PASS | FAIL` from `Review: CLEAR | REQUIRED`.
Gate controls the CI exit status. `Review: REQUIRED` means a human must examine
the identified risk; by itself it exits zero and is not a CI failure.

Record pull request baselines in proportion to the change. Ordinary
non-architecture pull requests use the normal concise template. Product-impact
work additionally records the applicable Product authority and affected
journey effects. Architecture-bearing work additionally records the relevant
semantic owners, canonical paths, compatibility paths, and changed-scope
`OE`/`PE`/`CB` conclusion. Mark non-applicable fields `N/A`; they do not require
an expanded packet.

## Submitting a pull request

Submit a non-trivial pull request only after its scope has been agreed in the
corresponding issue. Then:

1. Create a branch from the latest merged `dev` as described above.
2. Keep the patch limited to the agreed problem. Avoid unrelated cleanup,
   speculative abstractions, and drive-by refactors.
3. Add regression tests for the agreed behavior.
4. Run the relevant local checks, including at minimum:

   ```bash
   ruff check gbdraw/
   pytest tests/ -v -m "not slow"
   ```

5. Target `dev` and link the agreed issue.
6. Select exactly one change class in the
   [pull request template](./.github/pull_request_template.md), then explain the
   problem, agreed behavior, implementation boundary, validation, compatibility
   effects, scientific-output effects, and any intentional rendering changes.
7. For an ordinary architecture-bearing change, provide the concise owner/path,
   behavior-verification, deterministic-check, and rollback evidence required by
   the [architecture fitness-function ratchet](./docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md).
   Complete the full OE/PE/CB packet only for `ARCHITECTURE_EXCEPTION`.
8. If SVG output changes, review the actual output and the reference diff, and
   identify every intentional rendering change.

Passing CI is necessary but not sufficient for merge. The maintainer also
evaluates product behavior, scientific meaning, architecture, compatibility,
scope, and long-term maintenance cost. A pull request may be narrowed, split,
substantially revised, or replaced by a separate implementation. A valid
proposal may be adopted without merging the submitted patch.

By contributing, you agree that your contribution is licensed under the
project's [MIT License](./LICENSE.txt).
