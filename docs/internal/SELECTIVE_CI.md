# Selective CI impact planning

Status: pull-request, `dev` staging, and Gallery readiness routing are active.

## Purpose

gbdraw classifies each change before deciding which existing CI jobs are needed. The
classifier is deliberately coarse and conservative. It can avoid unrelated work only
when the changed paths are clearly lightweight and the directly preceding trusted
revision already has successful exact-SHA aggregate evidence.

The `Tests` workflow routes pull-request jobs from a plan produced by the trusted base
revision. Pushes to `dev` use the planner from the protected branch and inherit only
successful evidence from the direct parent commit. The `Gallery publication` workflow
uses the same protected-branch planner and direct-parent rule for its browser and
performance jobs.

## Impact classes

The three impact classes, from least to most consequential, are:

1. `metadata`: repository and development-support files that do not affect the
   application, package, Web bundle, Gallery output, or public operating instructions.
2. `documentation`: root Markdown files and files below `docs/`. These may participate
   in documented-contract tests, so they are not treated as no-ops.
3. `full`: runtime, tests, dependencies, build and packaging inputs, Web and Gallery
   files, workflow and planner control files, and every unknown path.

`tools/ci-impact-policy.mjs` is the only owner of the path allowlist and profile job
registries. The initial metadata allowlist is:

- `.agents/**`, `.claude/**`, `.codex/**`, and `.cursor/**`
- `.github/pull_request_template.md`
- `.dockerignore`, `.gitattributes`, and `.gitignore`
- `CITATION.cff`, `LICENSE.txt`, and `LICENSE_LIBERATION_FONTS.txt`

Root `*.md` files and `docs/**` are documentation. Everything else is full by default.
For multiple paths, the strongest class wins:

```text
metadata < documentation < full
```

Rename and copy records classify both the old and new path. Deletes classify the old
path. Unknown Git statuses, invalid paths, malformed or empty diffs, and missing Git
objects all fail closed to `full`.

## Plans and full fallback

The planner emits one versioned `ImpactPlan` JSON object. It records the profile,
impact, decision, stable basis code, change and workflow SHAs, required job IDs,
changed-path count, and inherited evidence when present. Job IDs, not display names or
human-readable reason text, form the routing contract.

The policy defines three profiles:

- `pr`: the `Tests` workflow for a pull request into `dev`
- `dev`: the `Tests` workflow for a push to `dev`
- `gallery`: the `Gallery publication` workflow for a push to `dev`

## Active pull-request routing

When the pull-request base has successful exact-SHA `Dev staging / gate` evidence, the
`pr` profile selects these project-owned jobs:

| Pull-request impact | Jobs that run |
| --- | --- |
| `metadata` | `CI impact plan`, `PR / gate` |
| `documentation` | `CI impact plan`, `Recipes standard (Python 3.11)`, `PR / gate` |
| `full` | `CI impact plan`, `Web change budget`, `Core PR (Python 3.11)`, `Recipes standard (Python 3.11)`, `Gallery (Python 3.11)`, `Lint`, `Web PR smoke`, `PR / gate` |

Documentation changes run the recipe suite because it owns documented-contract,
tutorial, CLI reference, and session compatibility checks. A full plan keeps the six
existing pull-request test commands unchanged.

`Web base policy (trusted base)` remains a separate required check. It evaluates Web
change policy from `pull_request_target`; it does not select tests or replace
`PR / gate`.

A full-impact change always produces a full plan without calling the GitHub Actions
API. Manual dispatch and a pull request carrying the `architecture-change` label also
force a full plan.

A metadata or documentation candidate can be selective only after inherited evidence
is verified. Missing, queued, in-progress, failed, cancelled, or malformed evidence;
API errors; and insufficient token access produce `decision=full` with
`basis=INHERITED_EVIDENCE_UNAVAILABLE`. These operational failures do not crash the
planner. Invalid inputs, malformed schemas, unsafe output writes, and programming errors
remain control-plane failures and do not get hidden by a full fallback.

## Active dev staging routing

For a push to `dev`, the planner compares `github.event.before` with the current SHA.
When the direct parent has successful exact-SHA `Dev staging / gate` evidence, the
`dev` profile selects these project-owned jobs:

| `dev` impact | Jobs that run |
| --- | --- |
| `metadata` | `CI impact plan`, `Dev staging / gate` |
| `documentation` | `CI impact plan`, `Recipes standard (Python 3.11)`, `Dev staging / gate` |
| `full` | `CI impact plan`, all 11 staging job IDs, `Dev staging / gate` |

The 11 staging job IDs and their existing commands, matrices, timeouts, and artifact
steps are unchanged. `workflow_dispatch` always produces a full plan.

The direct-parent requirement handles rapid consecutive pushes. If a later push
cancels the parent's run, or reaches the planner before the parent completes, the
current run cannot inherit that evidence and runs the full staging profile. This keeps
`cancel-in-progress: true` without admitting untested runtime changes.

`Dev staging / gate` validates the plan and every result from the current workflow run
with the protected branch's `tools/ci-impact.mjs`. It always creates aggregate evidence
for the current SHA. The promotion verifier continues to require that exact current-SHA
job; it does not accept the parent's job directly.

## Active Gallery readiness routing

For a push to `dev`, the Gallery planner compares `github.event.before` with the current
SHA. When the direct parent has successful exact-SHA `Gallery readiness / gate`
evidence, the `gallery` profile selects these jobs:

| Gallery impact | Jobs that run |
| --- | --- |
| `metadata` | `CI impact plan`, `Gallery readiness / gate` |
| `documentation` | `CI impact plan`, `Gallery readiness / gate` |
| `full` | `CI impact plan`, both `Gallery browser` matrix entries, `Gallery publication performance (projection)`, `Gallery readiness / gate` |

Documentation tests remain in the `Tests` workflow's recipe suite. Public Markdown and
files under `docs/` do not change the packaged Gallery sessions, Web bundle, browser
generation path, or performance projection, so the Gallery profile does not rerun
either Gallery test job for those changes.

The `browser` and `performance` job IDs retain their existing commands, matrices,
timeouts, and performance baseline. Every `workflow_dispatch` run is full. When
`complete_refresh=true`, the full performance job still runs its three-trial projection
and complete-refresh gates.

`Gallery readiness / gate` validates the Gallery plan and current-run job results with
the protected branch's shared gate validator. A selective run therefore still creates
successful aggregate evidence for the current SHA. Promotion continues to check that
current-SHA job rather than the inherited parent job.

## Direct evidence only

Selective planning inherits evidence from exactly one revision:

| Profile | Evidence revision | Workflow | Aggregate job |
| --- | --- | --- | --- |
| `pr` | pull-request base SHA | `.github/workflows/test.yml` | `Dev staging / gate` |
| `dev` | `github.event.before` | `.github/workflows/test.yml` | `Dev staging / gate` |
| `gallery` | `github.event.before` | `.github/workflows/gallery-publication.yml` | `Gallery readiness / gate` |

The planner reuses `verifyWorkflowEvidence()` from
`tools/check-promotion-readiness.mjs`. It does not copy GitHub API pagination or
workflow/job identity logic.

Only the direct base or parent is eligible because it proves that every unchanged
surface in the current revision was covered immediately before the lightweight change.
Searching for an older green run could cross an intervening unverified runtime change.
If the direct run was cancelled by concurrency, is still running, or did not succeed,
the current revision runs the full profile.

## Stable aggregate checks

The following display names are external admission contracts and stay fixed:

- `PR / gate`
- `Dev staging / gate`
- `Gallery readiness / gate`
- `Promotion / gate`

Each workflow continues to start on every relevant event and creates its aggregate job
for the current SHA. Pull-request, `dev` staging, and Gallery selection happen at job
level. Keeping the aggregate names and current-SHA jobs stable avoids branch-protection
changes and allows the existing promotion verifier to continue checking exact staging
evidence.

Workflow-level `paths` and `paths-ignore` filters are not used. Such filters can prevent
a required check from being created, leaving a pull request pending, and would remove
the current-SHA aggregate evidence required for promotion.

## Pull-request trust boundary

A pull request is untrusted input. Candidate code is tested, but it does not decide its
own admission route. The workflow checks out the pull-request base SHA under
`.ci-trusted-base` and runs that revision's `tools/ci-impact.mjs` for both planning and
aggregate validation. The planner reads Git history from the candidate repository root.
Candidate versions of `tools/ci-impact*.mjs` remain unit-test inputs only.

The bootstrap pull request ran the candidate planner only because its base had no
helper. That exception is no longer used. If the trusted base helper is absent,
unusable, or emits an invalid plan, `ci-impact` or `PR / gate` fails instead of guessing
a lightweight route. Changes below `.github/workflows/**`, `tools/ci-impact*.mjs`, and
the CI planner tests classify as full.

The `pull_request_target` workflow remains separate. It checks candidate Git data with
trusted base code and never executes the candidate checkout.

## Dev trust boundary

Code merged to protected `dev` is the staging control plane. Push and manual-dispatch
plans run the current checkout's helper. A change to a workflow, planner, policy, or
planner test classifies as `full`, so the commit that changes routing must run the full
staging and Gallery profiles before their aggregates can succeed.

This trust boundary differs from pull requests, where the candidate helper is tested
but the base revision decides the route. The Gallery workflow creates its own `gallery`
plan and gate; it does not consume the `Tests` workflow's `dev` plan.

## Aggregate validation

The shared gate validator fails closed unless:

- the plan schema, profile, required job list, and workflow SHA are valid;
- `ci-impact` succeeded;
- every required known job exists and succeeded;
- every unrequired known job exists and either succeeded or was skipped; and
- any additional job either succeeded or was skipped.

An unexpected failure or cancellation remains blocking even for a job the plan did not
require. Malformed plan or `needs` JSON is also blocking.

`PR / gate` runs the validator from `.ci-trusted-base`. `Dev staging / gate` and
`Gallery readiness / gate` run the same validator from the current protected-branch
checkout with their expected profile and the current workflow SHA.

## Rollout

Selective CI uses four separately reviewed stages:

1. Implemented: deterministic policy, adapter, validator, tests, documentation, and
   shadow summary.
2. Active: selective pull-request routing with the trusted base planner and validator.
3. Active: evidence-aware selection for exact `dev` staging runs.
4. Active: evidence-aware selection for exact Gallery readiness runs.

Each stage keeps the aggregate check names stable and can be reverted independently.
Reverting a routing stage restores the previous full-suite conditions without a
repository-settings change.

To roll back only Gallery selection, restore the `browser` and `performance` jobs to
their unconditional event behavior and restore the fixed two-job readiness check.
`Gallery readiness / gate` keeps the same display name, so repository settings do not
need to change. The planner and the other two profiles can remain active.

Selective CI does not alter `.github/workflows/deploy_web.yml`, deployment triggering,
or GitHub's CodeQL setup. Those systems remain outside the impact planner.

## Observability

The `CI impact plan` summary includes:

- profile, impact, decision, and stable basis;
- change base/head and workflow SHAs;
- changed-path count and a bounded, escaped path sample;
- planned required job IDs;
- inherited run and aggregate links when evidence succeeds; and
- the fallback code and bounded reason when evidence is unavailable.

The `pr`, `dev`, and `gallery` summaries state that routing is active and name the
applicable trust boundary. Tokens are never included in the plan or summary and are
redacted from CLI diagnostics.

## Troubleshooting

### A lightweight change produces a full plan

Check `basis` in the `CI impact plan` summary.

- `FULL_CHANGE`: at least one path is outside the two light allowlists.
- `UNKNOWN_OR_INVALID_CHANGE`: the diff was empty, malformed, unavailable, or contained
  an unsupported status/path.
- `MANUAL_FULL_RUN`: manual dispatch is intentionally full.
- `ARCHITECTURE_CHANGE`: the label intentionally forces full verification.
- `INHERITED_EVIDENCE_UNAVAILABLE`: inspect the listed evidence code and the direct
  base/parent workflow run. Do not substitute an older successful run.

### The planner job fails

Run the network-free focused tests first:

```bash
node --test tests/ci/*.test.mjs
```

Then verify that every required environment value is a correctly typed string and that
base, head, and workflow SHAs are complete object IDs. A failure to write
`GITHUB_OUTPUT` or `GITHUB_STEP_SUMMARY` is a control-plane error, not a reason to emit a
selective plan.

### The aggregate validator fails

Compare the plan profile and workflow SHA with the gate environment. Confirm that every
known job is present in `needs`, required jobs succeeded, and unrequired jobs did not
fail or get cancelled unexpectedly. Do not weaken the validator to accept a missing or
skipped required job.

### Consecutive dev pushes run a full profile

Inspect `basis` and the evidence failure in `CI impact plan`. If the direct parent's
`Tests` or `Gallery publication` workflow is queued, running, cancelled, failed, or
absent, the current push must use `INHERITED_EVIDENCE_UNAVAILABLE` and run the full
profile for that workflow. Wait for a successful direct parent before using a
lightweight commit when selective staging is important. Do not replace the direct
parent with an older green SHA or disable concurrency.

### Pull-request jobs are all skipped

Open `CI impact plan` first. If the planner failed or its `plan` output is missing,
`PR / gate` must fail even though the six test jobs are skipped. If the plan is valid,
compare its `requiredJobs` array with each job's result in the aggregate summary.

### The trusted helper cannot run

Confirm that the pull-request base SHA contains `tools/ci-impact.mjs`,
`tools/ci-impact-policy.mjs`, and `tools/check-promotion-readiness.mjs`. Do not run the
candidate helper as a fallback. Restore the trusted base prerequisite or keep the
control job failing.

## Policy-extension review checklist

Before making a path or profile less conservative, verify all of the following:

- The structural reason that the path cannot change runtime, packaging, Web, Gallery,
  public instructions, or another independently tested surface is documented.
- Historical changes and failures support the proposed classification.
- The unchanged surface is covered by the exact direct-parent/base evidence contract.
- A current-revision test owns every changed surface that will no longer use the full
  profile.
- Rename, copy, delete, mixed-change, unknown-path, and API-failure cases still fall back
  safely.
- The policy remains in `tools/ci-impact-policy.mjs`; no duplicate YAML or JSON owner is
  introduced.
- Existing full-profile commands, aggregate names, and promotion evidence contracts are
  unchanged.
- Focused policy, CLI, gate, workflow architecture, and applicable full regression tests
  pass without network access for unit coverage.
