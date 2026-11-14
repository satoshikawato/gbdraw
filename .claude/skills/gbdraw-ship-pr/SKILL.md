---
name: gbdraw-ship-pr
description: Open, watch, and merge a gbdraw pull request into dev - how many PRs and CI slots, body size and screenshots, the language check, auto-merge, when to update other branches, and how to read CI failures. Use when a branch is ready for a PR, while its CI runs, and when it merges.
---

# Ship a gbdraw pull request

This skill covers the procedure only; the content of the PR is decided
elsewhere. `AGENTS.md` and
[`docs/internal/WEB_CHANGE_POLICY.md`](../../../docs/internal/WEB_CHANGE_POLICY.md)
("Merging pull requests into `dev`") are the authority. Where they differ from
this skill, they win.

## Before opening

1. **Open few PRs.** Each runtime PR runs the full functional Playwright tier
   (about 40 minutes) on a job limit that every session shares. Fold small
   fixes into an open PR that touches the same files, and fold
   behaviour-neutral refactors into the PR they prepare. Keep a PR separate
   only when it is already in CI, when other work builds on it, or when its
   class must stay alone (authority, Product Contract, governance).
2. **Use at most two or three CI slots per session.** Serialize PRs that touch
   the same files: open the next one only after the previous one merges, then
   rebase it on `origin/dev`. Stack linearly; never stack through merge
   commits. Changelog entries never serialize PRs: each PR writes its own
   fragment under `changelog.d/` (`changelog.d/README.md`) and leaves
   `CHANGELOG.md` alone.
3. **Check version literals.** If the change bumps a schema or version
   constant, grep `tests/`, `docs/capture/`, `docs/recipes/`, and tool fixtures
   for the old value, including suites that run only on `main`, nightly, or by
   dispatch. Prefer asserting the exported constant.
4. **Predict the CI tier.**
   `node tools/ci-impact.mjs classify --base origin/dev --head HEAD`.
5. **Show visible GUI changes.** Take Before (base) and After (branch)
   screenshots with the same fixture, viewport, and state, cropped to the
   change; add a phone width when a popup or drawer changes. Push the PNGs to
   the orphan branch `pr-screenshots` as `pr/<PR#>/<name>-before.png` and
   `-after.png`, and embed
   `https://raw.githubusercontent.com/satoshikawato/gbdraw/pr-screenshots/pr/<PR#>/<file>`.
   Never commit screenshots that the Owner pasted.

## Write and open

- Follow `.claude/skills/write-clear-pull-request/SKILL.md` and choose the
  change class in the PR template.
- Keep the body under about 60,000 characters. GitHub refuses a body over
  65,536 characters, and "Web base policy (trusted base)" accepts any body
  GitHub accepts (`DEFAULT_MAX_BODY_BYTES`). Move long tables into a PR
  comment.
- In a bundled PR, list the bug IDs and the former PRs it contains.
- Run `node tools/check-pr-language.mjs --title "<title>" --body-file <file>`
  once, then `gh pr create --base dev`.
- A Product Contract change is named after its pull request
  ([Recording a revision](../../../docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md#recording-a-revision)).
  Create the PR as a draft, write `Revision #<number>` into its entries and
  `PD-OI-<number>` into a record it adds, run
  `node --test tests/ci/contract-revisions.test.mjs`, and push at once, so the
  push replaces the PR's first run within minutes; mark the PR ready afterwards.
- Where the machine defines a CI admission command, run every command that
  starts a pull request `Tests` run (a push to a branch with an open PR,
  `gh pr create`, `gh pr update-branch`, `gh run rerun`) through it, as one
  background command, and end your turn until it exits. Its pre-check runs CI's
  fast checks on the head merged with dev and refuses the push when one fails
  on the branch or the head conflicts with dev: fix the failure or rebase onto
  `origin/dev`, then run it again. A failure that dev also has does not stop the
  push. Skip the pre-check only with `CI_WAIT_SKIP_PRECHECK=<reason>` (a
  `dev-red-fix`, a check that is broken on dev). A push that stops a running
  `Tests` run states its reason (`CI_WAIT_WHY=conflict|failed-check|review|other`).

## Merge

- Before arming auto-merge, read dev's state (`AGENTS.md`, "When dev is
  red"). This prints the newest dev `Tests` push run with a result, or a
  rerun run (attempt above 1) without one:
  `gh run list -w Tests -b dev -e push -L 100 --json status,conclusion,attempt,url -q '[.[] | select((.status == "completed" and (.conclusion == "success" or .conclusion == "failure" or .conclusion == "timed_out" or .conclusion == "startup_failure")) or .attempt > 1)][0]'`.
  For a rerun run, read its earlier attempts newest first
  (`gh run view <id> --attempt <n> --json conclusion`), passing over those
  that ended cancelled or skipped, until one has a result; if none has, go on
  to the next run. The trusted `PR / gate` check is the authority
  (`docs/internal/SELECTIVE_CI.md`, "Dev state on pull requests"). While dev
  is red, arm auto-merge only on the fix.
- When the PR does not expand authority (see `AGENTS.md`), arm auto-merge:
  `gh pr merge <n> --auto --merge --match-head-commit <full sha>`. The Owner
  allows this for agent PRs into `dev`.
- An authority expansion waits for the Owner's approval. An
  `ARCHITECTURE_EXCEPTION` needs one ordinary Owner approval (a review or a
  comment); it stays valid across merges and fixes that leave the exception
  rows unchanged. Never post an approval yourself.
- Run `gh pr update-branch` on other open PRs only under the rules in
  `WEB_CHANGE_POLICY.md`, on a conflict, or to pick up a fix they need, never
  just to stay current. When such a merge lands while other PRs have
  auto-merge armed, disarm theirs before updating.

## Watch CI

- Wait with the machine's PR wait command (see the local instructions), as one
  background command per session. It exits at the first watched PR that needs
  its owner and prints one line per watched PR, `WAITING` for those still
  waiting; start it again for those. Never write a watcher script per task. On
  a machine without that command, `python3 .claude/skills/gbdraw-ship-pr/scripts/watch_prs.py --prs <n> ...`
  is the fallback; it polls every five minutes, as `AGENTS.md` requires.
- Act on each line's state.
- Re-pushes, rebases, and CI fixes whose log names the cause go to
  `gbdraw-mechanic`; resume the implementer only for a fix that bears on the
  design.
- Before you act on an event line (rebase, update-branch, merge, or report the
  PR done), re-read the PR with
  `gh pr view <n> --json state,mergeStateStatus,statusCheckRollup`. The Owner
  may have merged it already, and a green line covers only the checks that
  existed at that poll, not those of a workflow that starts later.
- A run whose jobs all passed but which ended `cancelled` is recovered with
  `gh run rerun <id> --failed`, run through the admission command.
- When a local check predicts a CI failure, settle it now with CI's own
  comparison, then fix it or state that nothing is needed. "Web change budget"
  and "Web base policy (trusted base)" measure the PR's own change, from the
  merge base of `pull_request.base.sha` and the head; reproduce them with
  `--base <PR base.sha> --head <head sha>`. A baseline lowered or a guard
  changed on `dev` after the branch forked does not fail them, so do not
  rebase for it.
- Read only failing logs, into a file:
  `gh run view <id> --log-failed > "$TASK_DIR/logs/<id>.log"`. Quote the
  failing lines, not the whole log.

## After the merge

Remove the worktree and the branch (skill `gbdraw-worktree`) and record the
merge in the task's `STATE.md` (skill `gbdraw-campaign`).
