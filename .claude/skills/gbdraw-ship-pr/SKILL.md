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
   commits.
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
- Keep the body under about 60 KB. "Web base policy (trusted base)" fails above
  64 KiB (`DEFAULT_MAX_BODY_BYTES`), counted after screenshot URLs expand.
  Move long tables into a PR comment.
- In a bundled PR, list the bug IDs and the former PRs it contains.
- Run `node tools/check-pr-language.mjs --title "<title>" --body-file <file>`
  once, then `gh pr create --base dev`.

## Merge

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

- Run one watcher per session, not one per agent:
  `python3 .claude/skills/gbdraw-ship-pr/scripts/watch_prs.py --prs <n> ...`
  under Monitor or in the background, and act on its event lines (merged,
  closed, conflict, failed, cancelled, green). It polls every five minutes, as
  `AGENTS.md` requires.
- A run whose jobs all passed but which ended `cancelled` is recovered with
  `gh run rerun <id> --failed`.
- When a local check predicts a CI failure, settle it now with CI's own
  comparison, then fix it or state that nothing is needed. "Web change budget"
  diffs `pull_request.base.sha` against the head (two-dot), so reproduce it
  with `--base <PR base.sha>`, not `origin/dev`. "Web base policy (trusted
  base)" compares the `dev` tree with the PR head, so a branch that lacks a
  later `dev` baseline reduction must merge or rebase `dev`.
- Read only failing logs, into a file:
  `gh run view <id> --log-failed > "$TASK_DIR/logs/<id>.log"`. Quote the
  failing lines, not the whole log.

## After the merge

Remove the worktree and the branch (skill `gbdraw-worktree`) and record the
merge in the task's `STATE.md` (skill `gbdraw-campaign`).
