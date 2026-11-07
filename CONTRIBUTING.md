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

Where another document owns a rule, this guide links to it instead of
repeating it. The [norms index](./docs/internal/NORMS_INDEX.md) lists the
design rules and Product Contract records with the tests and CI jobs that check
them.

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
- compatibility and deprecation policy (the current one is in
  [Compatibility and deprecation](./docs/REFERENCE/compatibility-and-deprecation.md));
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

Report a security vulnerability privately as [`SECURITY.md`](./SECURITY.md)
describes, not in a public issue.

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

### Large test inputs

Input files of 2 MiB or more under `tests/test_inputs/` and `examples/` are
not in git. They are assets of a test-data GitHub Release, listed with their
size and sha256 in `tests/test_data.json`. pytest places them at the start of
a run, from your cache (`~/.cache/gbdraw/test-data`, or the folder that
`GBDRAW_TEST_DATA_CACHE` names), from git history, or from the Release. Tools
that read them outside pytest need one run of:

```bash
python tools/test_data.py fetch
```

When an input cannot be placed, a pytest run in CI (where the `CI`
environment variable is set) stops before any test runs and names the missing
files. A local run continues but skips each test module whose source names a
missing file; a test that reaches a missing file in another way fails.
`python tools/test_data.py fetch --require` exits with status 1 while an input
is missing.

To add or change such an input, put it at its path and run
`python tools/test_data.py add <path>` (it uploads to the Release, so it
needs write access), then commit `tests/test_data.json` and `.gitignore`.
If you do not have write access, say in the issue where the file can be
downloaded; the maintainer uploads it. Generated figures and Sessions stay in
git.

The test-data Release follows these rules:

- `tests/test_data.json` is the authority. A copy from the cache, git history,
  or a download is used only when its size and sha256 match the manifest.
- The bytes of an asset never change. An asset's name starts with the first 12
  hexadecimal digits of its sha256, so `add` uploads a changed input as a new
  asset, and older commits still fetch the asset their manifest names.
- Never delete a test-data Release, or one of its assets, that
  `tests/test_data.json` names in a commit of `dev` or `main`.
- Only the maintainer creates a test-data Release. Its tag is
  `test-data-YYYY-MM`, the year and month it was created (for example
  `test-data-2026-11`), and `RELEASE_TAG` in `tools/test_data.py` names the
  Release that `add` uploads to. The tag points at a commit that holds no
  workflow files, because publishing a Release starts the release workflow of
  the tagged commit. The Release is a prerelease and is not marked as the
  latest release, so the newest gbdraw release stays the latest.

## Running the checks CI runs

```bash
# Lint (blocks CI on failure)
ruff check gbdraw/

# Fast test suite (excludes slow/regression-heavy tests)
pytest tests/ -v -m "not slow"

# Full suite, including slow tests (CI runs it in release acceptance only)
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

Which Playwright installation runs which check, and what to use when Node's
`@playwright/test` is missing: [Browser / Playwright Checks](./AGENTS.md#browser--playwright-checks).

## Changing rendered SVG output

Many tests compare generated SVGs against tracked references in
`tests/reference_outputs/`. If your change intentionally alters diagram
geometry:

```bash
# 1. Confirm the differences are the ones you intended. Actual output is
#    retained under pytest's temporary directory, not reference_outputs/.
pytest tests/test_output_comparison.py::TestOutputComparison -v

# 2. Regenerate references only after reviewing that output
pytest tests/test_output_comparison.py::TestGenerateReferences \
  --update-reference-outputs -v

# 3. Review the diff and re-verify
git diff -- tests/reference_outputs/
pytest tests/test_output_comparison.py::TestOutputComparison -v
```

Without `--update-reference-outputs`, reference-generation tests are skipped and
normal test runs do not write to `tests/reference_outputs/`.

Do not update reference outputs merely to make a failing test pass. First
confirm that the new geometry is intentional and correct; the reference files
record the project's intended rendering behavior.

## Changing documentation

Public documentation has four teaching routes:
[Tutorials](./docs/TUTORIALS/README.md),
[Technical documentation](./docs/REFERENCE/README.md),
[FAQ](./docs/FAQ.md), and [Gallery](./docs/GALLERY.md).
[`docs/DOCS.md`](./docs/DOCS.md) is their stable index and the navigation
authority.

- Use the fewest public pages that answer distinct reader questions. Coverage
  does not require one page per capability, workflow, or interface; a new
  interface or test scenario does not automatically require a new top-level
  page.
- Before adding a public page, identify which existing document owns the
  reader's question, and record the reader question, the existing owner, and a
  `keep`, `merge`, `delete`, or `new` disposition. Add a page only when editing
  an existing owner cannot answer the question clearly.
- A Tutorial needs a deliberate learning progression to a finished result.
  Separate GUI, CLI, or Python evidence does not by itself justify separate
  public pages.
- Technical documentation owns exact behavior and contracts. FAQ owns concise
  decisions and troubleshooting, and links to the technical owner. Gallery is
  for discovering finished outcomes.
- Executable scenarios may support one public page or no public page. Do not
  turn an evidence inventory into a public-page inventory. Internal plans,
  audits, and executable evidence registries belong under `docs/internal/`.
- Keep user- and maintainer-facing prose concise, concrete, and technically
  accurate. Optimize internal implementation plans for execution fidelity.
- Preserve exact technical terms, UI labels, CLI options, identifiers, and
  scientifically necessary qualifications.

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
- Extend an existing boundary instead of adding a parallel path, so that the
  CLI, Web, and Python interfaces converge on shared internals:
  [Key Architecture](./CLAUDE.md#key-architecture) in `CLAUDE.md`.
- Do not add speculative options, fallbacks, or abstractions for hypothetical
  future use. Keep changes within the agreed problem.
- Architecture-bearing changes must follow the
  [architecture fitness-function ratchet](./docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)
  and remove superseded owners and paths in the same change.

## Branch roles

Approved work starts from the latest merged `dev`, and a maintainer opens the
promotion pull request from `dev` to `main`. The roles of work branches,
`dev`, and `main`, what each admits and runs, and the rule that direct `main`
hotfixes are unsupported are in
[Delivery lifecycle](./docs/internal/WEB_CHANGE_POLICY.md#delivery-lifecycle)
of the normative Web change, review, and delivery policy, which also defines
lifecycle ownership and trust boundaries.

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

An ordinary architecture-bearing `STANDARD` change removes every superseded
owner or path in the same change
([Merge acceptance](./docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md#merge-acceptance)).
What it records, and what only `ARCHITECTURE_EXCEPTION` needs:
[Required author evidence](./docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md#required-author-evidence).

Apply the `architecture-change` label to a Web implementation that changes an
architecture owner, canonical path, compatibility path, privileged boundary, or
lifecycle responsibility. The label selects the architecture size-review
profile and routes attention. It does not waive deterministic gates.

Policy reports give `Gate: PASS | FAIL` and `Review: CLEAR | REQUIRED`
separately:
[Gate and Review are independent](./docs/internal/WEB_CHANGE_POLICY.md#gate-and-review-are-independent).

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
7. For an architecture-bearing change, provide the evidence that
   [Required author evidence](./docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md#required-author-evidence)
   lists for its change class; an ordinary change removes every superseded
   owner or path in the same change.
8. If SVG output changes, review the actual output and the reference diff, and
   identify every intentional rendering change.

Passing CI is necessary but not sufficient for merge. The maintainer also
evaluates product behavior, scientific meaning, architecture, compatibility,
scope, and long-term maintenance cost. A pull request may be narrowed, split,
substantially revised, or replaced by a separate implementation. A valid
proposal may be adopted without merging the submitted patch.

By contributing, you agree that your contribution is licensed under the
project's [MIT License](./LICENSE.txt).
