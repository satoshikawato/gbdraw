# AGENTS.md

Guidance for automated agents working in this repository.

Where another document owns a rule, this file links to it instead of repeating
it. The [norms index](docs/internal/NORMS_INDEX.md) lists the design rules and
Product Contract records with the tests and CI jobs that check them.

## Read First
- コード1行増やすごとに技術負債が増えると心得よ。負債は行数ではなく、増えた変更点・分岐・重複で数える。
- See `CLAUDE.md` for project-wide guidance.
- If working on the web UI, also read `gbdraw/web/CLAUDE.md`.

## Agent Branching

- Before starting work that may create a branch or commit, fetch `origin` and
  create a fresh work branch from the latest `origin/dev`:
  `git switch --no-track -c <branch-name> origin/dev`.
- All agent-authored commits must be made on a work branch derived from
  `origin/dev`. Do not base agent work on `main` or `master`.
- Do not configure a work branch to track `main` or `dev`. Before committing or
  pushing, verify the current branch and its upstream; publish only to the
  same-named remote work branch.
- Never commit directly to, or push directly to, `main` or `dev` unless the user
  explicitly authorizes that exact direct target in the current request.
- Open every pull request against `dev`. Open a pull request that depends on
  another only after that one has merged, rebased onto `origin/dev`. Never merge
  a pull request into another work branch; confirm `base` is `dev` before
  enabling auto-merge or merging.
- Before opening or merging a pull request, confirm that every file its
  Changes section names is in its diff (`gh pr diff <n> --name-only`).
- Do not enable auto-merge on, or merge, a pull request that expands authority
  (Web design rule text that permits more, a guard registration removed, or an
  entry added to a registered baseline or allowlist). Wait for the Owner's
  approval.

## When dev is red

- dev is red when the newest `Tests` run of a push to `dev` that finished
  with a result failed (`failure`, `timed_out`, or `startup_failure`). Runs
  that ended `cancelled`, `skipped`, `neutral`, `stale`, or `action_required`
  do not count. A run being rerun counts by its previous attempt, so a failed
  run stays red while it reruns. Read the state before you open a pull
  request, arm auto-merge on one, or merge one (command in the
  `gbdraw-ship-pr` skill, "Merge").
- While dev is red, merge only a pull request that fixes the failure: it
  carries the label `dev-red-fix` and names the failing job and run in its
  body. Other pull requests may run CI but stay unmerged. `PR / gate` enforces
  the label ([Dev state on pull requests](docs/internal/SELECTIVE_CI.md#dev-state-on-pull-requests)).
- The session that finds dev red reruns a job it believes flaky once. If the
  job fails again, it fixes the failure or tells the session that merged its
  cause (or the Owner), and records this in its task log.
- A dev-tier job that fails on a pull request's merge commit is that pull
  request's author's to fix first.

## Workspace Cleanup

- GitHub deletes a pull request's head branch when the pull request merges.
  Your local worktree and branch stay; remove them yourself.
- When your pull request merges and you no longer need its worktree, run
  `git status` in it, then `git worktree remove <path>` and
  `git branch -D <branch>`. Keep the worktree while it holds uncommitted work
  you still need.
- Remove only the worktrees and branches you created. Merged or idle does not
  mean unused: another agent may still work in its worktree. To clean up
  across agents, ask each owner first.
- Keep logs, evidence, scratch scripts, virtual environments, and backups in
  one folder per task outside the repository, never inside `.worktrees/`.
- When a task ends, delete its remaining worktrees, local branches, virtual
  environments, caches, and backups. Keep its report and handoff notes. PR body
  drafts and session logs need no archive. A backup script keeps at most three
  snapshots.

## Project Summary

- `gbdraw` is a Python 3.10+ bioinformatics tool for publication-quality genome diagrams.
- Outputs: SVG, PNG, PDF, EPS, PS.
- Main package lives in `gbdraw/`; tests in `tests/`; web UI in `gbdraw/web/index.html` with JS modules under `gbdraw/web/js/` (app entry points in `gbdraw/web/js/app/` with submodules grouped by feature).

## Common Commands

- Test, lint, and reference-output commands:
  [Running the checks CI runs](CONTRIBUTING.md#running-the-checks-ci-runs) and
  [Changing rendered SVG output](CONTRIBUTING.md#changing-rendered-svg-output)
  in `CONTRIBUTING.md`.
- Browser wheel, hosted bundle, and build commands:
  [Local build and verification](gbdraw/web/CLAUDE.md#local-build-and-verification)
  in `gbdraw/web/CLAUDE.md`.
- Allow at least 30 minutes for test commands before treating them as timed
  out, and monitor long runs incrementally. Keep shorter test-owned timeout
  assertions unchanged unless the task explicitly requires changing them.
- When monitoring CI or other remote status, poll no more often than once every
  five minutes unless the user explicitly requests a different interval. Do not
  use 10-second polling loops.

## Browser / Playwright Checks

- Do not conclude that browser testing is unavailable just because `node_modules/`, `package.json`, or `@playwright/test` is missing at the repo root. This workspace may have Playwright installed through Python/conda instead.
- Check both paths when browser verification matters:
  - `command -v playwright && playwright --version`
  - `python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"`
- The JavaScript specs under `tests/web/*.playwright.spec.js` require Node's `@playwright/test`. Verify it with `node -e "console.log(require.resolve('@playwright/test'))"` before trying to run those specs.
- If Node's `@playwright/test` is unavailable, use Python Playwright for targeted browser checks instead of skipping browser verification.
- In Codex/agent sandboxes, Chromium may fail with `sandbox_host_linux.cc ... Operation not permitted`. When that happens, rerun the same local browser check with the required sandbox escalation rather than reporting that Playwright is unavailable.

## Expectations When Editing

- Audit production-code and test diffs separately before handoff. Optimize for
  fewer change points, branches, and duplicated behavior rather than raw line counts.
- Follow the [architecture fitness-function ratchet](docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)
  for every architecture-bearing change. An ordinary change removes every
  superseded owner or path in the same change
  ([Merge acceptance](docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md#merge-acceptance));
  the evidence each change class records is in
  [Required author evidence](docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md#required-author-evidence).
- Web UI layout (a single-page app with no build step, ES modules under
  `gbdraw/web/js/`, `create*` entry points in `app/*.js` with larger helpers in
  feature subfolders): [Non-negotiable properties](gbdraw/web/CLAUDE.md#non-negotiable-properties)
  and [Module ownership](gbdraw/web/CLAUDE.md#module-ownership) in `gbdraw/web/CLAUDE.md`.
- Do not add a CDN or other runtime network dependency; vendor the dependency
  and keep the CSP aligned with same-origin needs:
  [Assets, CSP, and privacy](gbdraw/web/CLAUDE.md#assets-csp-and-privacy).
- `tests/reference_outputs/` is read-only in normal test runs; regenerate it
  only for an intentional, reviewed change:
  [Changing rendered SVG output](CONTRIBUTING.md#changing-rendered-svg-output).
- Do not manually edit generated artifacts under `dist/` or `gbdraw.egg-info/`.
- The browser wheel `gbdraw/web/gbdraw-<version>-py3-none-any.whl` is generated
  and gitignored; do not commit it:
  [Non-negotiable properties](gbdraw/web/CLAUDE.md#non-negotiable-properties).

## Owners and Boundaries in Every Fix

Every fix or other change to code under `gbdraw/` applies these rules to the
code and data it touches. They do not widen a task's scope or authorize
cleanup beyond it.

- Before fixing, name the channel: which owner writes the value, which code
  reads it, and by which key. Fix the channel, not only the reader that showed
  the bug.
- One owner per value and per transformation. A value shown on two surfaces
  (the drawing and an editor, the live preview and Generate, Python and the
  Web) comes from one owner function. Never add another derivation. Remove the
  second one in the same change when it is within the task's scope; otherwise
  list it as remaining work. Fix a live versus Generate difference in the
  shared compile path (R3 in `gbdraw/web/CLAUDE.md`) and add a parity case.
- Find a row, slot, record, or feature by its identity (slot id, row id,
  feature id), never by its caption, display text, or kind when an identity
  exists.
- Type the data that crosses a boundary you touch (between modules, between
  Python and the Web, in a Session): a named JSDoc typedef, dataclass, or
  `TypedDict` that says which fields are required. Use no `any` or `Any`, and no
  optional chaining that hides a missing required field. A type baseline never
  rises without the Owner's approval (R14 in `gbdraw/web/CLAUDE.md`,
  [Type checking](docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md#type-checking)
  in the architecture ratchet).
- Do not hand-copy a Python-owned value or rule into JavaScript. Carry it
  through an existing generator with a `--check` test (for example
  `tools/generate_web_contracts.py`); a twin that cannot be avoided follows R9.
- When the conforming fix needs a structure `dev` does not have yet, take the
  largest step that does not need it and list the rest as remaining work.
- The step handoff and the pull request body include a short "Boundary
  analysis": the channel, the rules the fix applies, what it removed, any twin
  that remains, and the remaining work. Describe each rule in words; a design
  rule number such as R3 may follow the words, but campaign-internal IDs may
  not.

## Product Behavior Semantics Review

For changes affecting Web runtime or normative Web behavior contracts:

- Follow [`docs/internal/PRODUCT_IMPACT_RATCHET.md`](docs/internal/PRODUCT_IMPACT_RATCHET.md).
- Identify every reachable user-visible difference or continuation and separate
  the product outcome from the owner or path used to implement it.
- Search explicit supported behavior, compatibility commitments, authoritative
  domain rules, and accepted base-branch decisions. Current code and tests are
  evidence, not automatic product authority. Cite a `BD-###` only when it exists
  in the base branch.
- Before implementation, use the same policy's developer preflight when a
  proposed change may alter a material user effect but no registered concern is
  triggered. Use
  [`docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md`](docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md)
  when evidence or Product judgment is required.
- Compare option realization with AND-of-OR requirements. Matching option IDs do
  not prove that independent behavior contributions remain satisfied.
- Keep the common path automatic for changes with no registered Product Impact
  delta or unresolved outcome.
- Stop the affected convergence when materially different user-visible outcomes
  remain unresolved. Do not select a product option autonomously.
- Present the Product Decision Owner with stable outcome choices, decision
  routes, and the documented `PRODUCT_DECISION` response template.
- After an explicit human choice, serialize only that outcome and show the
  generated machine representation for review. Do not infer missing rationale,
  retirement intent, or accepted residual risk, and do not broaden the choice.
- Candidate machine authority (Product Impact map, durable decisions,
  architecture rules, checker) never authorizes the same candidate runtime.
  A static Product Contract change may merge with its implementation only
  through the reviewed co-change route in
  [`docs/internal/WEB_CHANGE_POLICY.md`](docs/internal/WEB_CHANGE_POLICY.md#static-product-contract-co-change).
  The explicit human receipt authorizes it, not the candidate text.

## Showcase Figure Quality

- Treat `examples/gbdraw_social_preview.png` as an owner-maintained asset. Automated agents must not regenerate, replace, or edit it.
- Treat every public figure that demonstrates a feature as a finished example, even when the underlying change is small.
- Start from a realistic Gallery-quality recipe or session when one exists, and change only the setting being demonstrated where practical.
- Keep the labels, legend, color rules, record metadata, quantitative tracks, and comparison context that make the source figure useful. Remove an element only when the example has a clear reason to omit it.
- Keep minimal smoke diagrams in tests. Do not use them as tutorial images, Gallery entries, release examples, or other public showcases.
- Render and visually inspect the final artifact at a readable scale. Confirm that its documented command, session, or reproduction recipe generates the displayed figure.

## Pull request communication

- Before creating, changing, or reviewing PR wording,
  read and apply `.claude/skills/write-clear-pull-request/SKILL.md`.
- Prepare the title and body file. Run
  `node tools/check-pr-language.mjs --title "<title>" --body-file <path>` once
  before `gh pr create`, using the same wording.
- Rerun only after a material rewrite before merge; never for ordinary shell commands.
- Do not use `gh pr create --fill`, `--fill-first`, or `--fill-verbose`.

## Completion Handoff

- After completing an implementation, treat the session as one commit and provide a proposed commit title and summary in English.
