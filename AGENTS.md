# Agent guidance for gbdraw

Shared repository rules for Codex, Claude Code, and other coding agents.
`gbdraw` is a Python 3.10+ genome-diagram tool with a static browser UI.

## Work and completion

- Inspect the working tree before editing. Preserve unrelated changes. An
  in-scope dirty implementation may be replaced when needed for correctness;
  inspect its diff first. A requested Gallery refresh may overwrite dirty
  generator-owned sessions, sources, examples, thumbnails, and `examples.json`.
- Prefer the fewest change points, execution paths, and duplicated owners that
  satisfy the request. Remove superseded paths in the same change.
- Complete the requested scope through implementation, relevant verification,
  and correction of failures caused by the change. A review-only or planning
  request remains review-only or planning; a scoped phase is not the whole plan.
- Continue authorized local edits, builds, and disposable-fixture checks without
  approval at each step. Reuse authorization within its stated scope. Ask before
  an unauthorized publication, deployment, tag, push, external message, or
  destructive action. Before retrying a remote mutation, check whether it succeeded.
- Do not weaken contracts, expected results, or required gates to pass a check.
  Diagnose the failed boundary; report a concrete blocker if it cannot be resolved.

## Branches

- For new work that may create a branch or commit, fetch `origin` and start from
  the latest `origin/dev`: `git switch --no-track -c <branch-name> origin/dev`.
  Continue an explicitly requested existing work branch for a resumed task.
- Agent commits belong on work branches derived from `origin/dev`, without an
  upstream of `main` or `dev`. Verify branch and upstream before committing or
  pushing; publish only to the same-named remote work branch.
- Never commit or push directly to `main` or `dev` unless the user explicitly
  authorizes that exact direct target in the current request.

## Read for the task

Read only the relevant sections and references; a wording fix does not require
an architecture or runtime audit.

| Task | Guidance |
| --- | --- |
| Python/API architecture, persisted compatibility, public-doc ownership | [CLAUDE.md](CLAUDE.md) |
| Web runtime, UI, sessions, packaging, browser checks | [gbdraw/web/CLAUDE.md](gbdraw/web/CLAUDE.md) |
| Architecture-bearing owner/path changes | [Architecture ratchet](docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md) |
| Web runtime or normative Web behavior contracts | [Product Impact Ratchet](docs/internal/PRODUCT_IMPACT_RATCHET.md) |
| Reproducible procedural documentation | [.agents/skills/love-me-love-my-docs/SKILL.md](.agents/skills/love-me-love-my-docs/SKILL.md) |
| Gallery tutorial text, JSON, or screenshots | [.agents/skills/web-gallery-screenshot-maintenance/SKILL.md](.agents/skills/web-gallery-screenshot-maintenance/SKILL.md) |
| Creating, editing, or reviewing PR wording | [.agents/skills/write-clear-pull-request/SKILL.md](.agents/skills/write-clear-pull-request/SKILL.md) |

Architecture changes normally need concise owner/path evidence; full before/after
OE, PE, and CB sets apply only to the ratchet's defined exceptions.

For Web behavior changes, preserve every jointly required user outcome and use
base-branch authority. Run the policy's developer preflight when a material user
effect may change without a registered concern. Unresolved materially different
outcomes require the Product Decision route and response template in that policy;
continue unaffected work. Candidate authority never authorizes its own runtime.

Use `execute-plan-with-evidence` only when explicitly invoked or when the user
explicitly asks to execute the whole plan with evidence. A named plan or an
ordinary instruction prompt alone does not activate it.

## Verification

Choose checks for the changed behavior and the repository's required gates.
Prefer shared contract tests, representative workflows, and continuing benchmarks
over a new test for every symptom. Reuse valid evidence when the relevant code,
inputs, environment, and required acceptance conditions have not changed.

Review production, test, documentation, and generated diffs separately in one
final pass. Repeat a check or review only for a new change, failure, or unresolved
concern. Instructions-only edits normally need skill/frontmatter/link checks,
not diagram generation or a full application suite; required CI still applies.

```bash
pytest tests/ -v -m "not slow"  # broader Python gate when warranted
pytest tests/ -v               # full suite when required
ruff check gbdraw/
python -m build
```

Allow at least 30 minutes before treating a test command as timed out; monitor
long runs incrementally. Keep test-owned timeout assertions unchanged unless
requested. Browser setup and sandbox recovery are in the Web guide's
"Local build and verification" section. Put temporary browser specs in a dedicated
task directory, not directly under `/tmp` where discovery can scan unrelated files.

## Generated outputs and figures

- Normal tests must not write `tests/reference_outputs/`. For an intentional,
  reviewed geometry change, run `TestGenerateReferences` with
  `--update-reference-outputs`, inspect the SVG diff, then rerun
  `pytest tests/test_output_comparison.py::TestOutputComparison -v`.
- Do not hand-edit `dist/`, `gbdraw.egg-info/`, or the gitignored browser wheel.
  Use `python tools/prepare_browser_wheel.py` when needed; add
  `--refresh-cache-bust` only for a deployable bundle. Never commit the wheel.
- `examples/gbdraw_social_preview.png` is owner-maintained; agents must not edit
  or regenerate it.
- Public feature figures must be finished examples. Start from a realistic
  Gallery recipe when available and retain useful labels, metadata, legends,
  tracks, and comparison context. Minimal smoke diagrams belong in tests.
  Reproduce and visually inspect changed public figures at a readable scale.

## Handoff

Report the result, relevant checks, and remaining limitations. After implementation,
provide one proposed English commit title and a short English summary for the
session's changes. This handoff does not itself require a commit, push, or PR.
