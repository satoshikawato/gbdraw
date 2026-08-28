# Selective CI impact planning

Status: staged rollout; the planner currently runs in shadow mode.

## Purpose

gbdraw classifies each change before deciding which existing CI jobs are needed. The
classifier is deliberately coarse and conservative. It can avoid unrelated work only
when the changed paths are clearly lightweight and the directly preceding trusted
revision already has successful exact-SHA aggregate evidence.

The first rollout stage records the plan in the GitHub Actions job summary but does not
route test jobs. All existing pull-request and `dev` staging test conditions still run.

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
changed-path count, and inherited evidence when present. Job IDs—not display names or
human-readable reason text—form the routing contract.

Profiles cover the three eventual consumers:

- `pr`: the `Tests` workflow for a pull request into `dev`
- `dev`: the `Tests` workflow for a push to `dev`
- `gallery`: the `Gallery publication` workflow for a push to `dev`

A full-impact change always produces a full plan without calling the GitHub Actions
API. Manual dispatch and a pull request carrying the `architecture-change` label also
force a full plan.

A metadata or documentation candidate can be selective only after inherited evidence
is verified. Missing, queued, in-progress, failed, cancelled, or malformed evidence;
API errors; and insufficient token access produce `decision=full` with
`basis=INHERITED_EVIDENCE_UNAVAILABLE`. These operational failures do not crash the
planner. Invalid inputs, malformed schemas, unsafe output writes, and programming errors
remain control-plane failures and do not get hidden by a full fallback.

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
for the current SHA. Selection happens at job level in later rollout stages. Keeping the
aggregate names and current-SHA jobs stable avoids branch-protection changes and allows
the existing promotion verifier to continue checking exact staging evidence.

Workflow-level `paths` and `paths-ignore` filters are not used. Such filters can prevent
a required check from being created, leaving a pull request pending, and would remove
the current-SHA aggregate evidence required for promotion.

## Pull-request trust boundary

A pull request is untrusted input. Candidate code is tested, but it must not decide its
own admission route. Once selective pull-request routing is enabled, the workflow will
execute the planner and gate validator from the pull-request base SHA in a separate
trusted checkout. Candidate versions of the planner remain test inputs only.

The initial shadow-mode pull request is the bootstrap exception: no planner exists on
its base, so the candidate planner produces an observational summary while all existing
jobs still run. Changes below `.github/workflows/**`, `tools/ci-impact*.mjs`, and the CI
planner tests classify as full.

The `pull_request_target` workflow remains separate. It checks candidate Git data with
trusted base code and never executes the candidate checkout.

## Aggregate validation

The shared gate validator fails closed unless:

- the plan schema, profile, required job list, and workflow SHA are valid;
- `ci-impact` succeeded;
- every required known job exists and succeeded;
- every unrequired known job exists and either succeeded or was skipped; and
- any additional job either succeeded or was skipped.

An unexpected failure or cancellation remains blocking even for a job the plan did not
require. Malformed plan or `needs` JSON is also blocking.

The validator is implemented and tested during shadow mode. Existing shell aggregate
checks remain in place until the rollout stage that connects selective routing.

## Rollout

Selective CI is introduced in four separately reviewed pull requests:

1. Add the deterministic policy, adapter, validator, tests, documentation, and shadow
   summary. Keep all existing test routing.
2. Enable selective routing for pull requests using the trusted base planner.
3. Enable evidence-aware routing for exact `dev` staging runs.
4. Enable evidence-aware routing for Gallery readiness runs.

Each stage keeps the aggregate check names stable and can be reverted independently.
Reverting a routing stage restores the previous full-suite conditions without a
repository-settings change.

## Observability

The `CI impact plan` summary includes:

- profile, impact, decision, and stable basis;
- change base/head and workflow SHAs;
- changed-path count and a bounded, escaped path sample;
- planned required job IDs;
- inherited run and aggregate links when evidence succeeds; and
- the fallback code and bounded reason when evidence is unavailable.

In shadow mode the summary explicitly states that routing is observational. Tokens are
never included in the plan or summary and are redacted from CLI diagnostics.

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
