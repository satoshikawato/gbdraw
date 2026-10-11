---
name: gbdraw-worktree
description: Set up, use, and remove a gbdraw git worktree for one task - create it from origin/dev, keep temporary files and logs outside the repository, run tests against the worktree's own package without reinstalling, push finished branches, and clean up after the PR merges. Use when starting work that needs its own branch or worktree, and when finishing it.
---

# gbdraw worktree lifecycle

This skill covers mechanics only. It does not decide what the task includes.

The concrete paths for this machine come from the local instructions
(`CLAUDE.local.md` or `~/.claude/CLAUDE.md`). The steps below name them:

- `$WORK_ROOT`: the repository clone that holds `.worktrees/`.
- `$TASK_DIR`: the task's folder outside the repository, for logs, evidence,
  scripts, screenshots, and `STATE.md`.

## Create

```bash
git -C "$WORK_ROOT" fetch origin
git -C "$WORK_ROOT" worktree add --detach ".worktrees/<name>" origin/dev
cd "$WORK_ROOT/.worktrees/<name>"
git switch --no-track -c <prefix>/<topic> origin/dev   # fix, feat, refactor, docs, test
export TMPDIR="$TASK_DIR/tmp"; mkdir -p "$TMPDIR"
```

- Keep logs, probes, screenshots, and virtual environments out of `.worktrees/`.
- Do not run `git stash`. Every worktree of the clone shares one stash list.

## Run code from the worktree

- Do not run `pip install` from a worktree, editable or not. The Python
  environment is shared by every session, and an editable install that points
  at a worktree breaks once the worktree is removed.
- `python -m pytest` run from the worktree root imports the worktree's
  `gbdraw`. For other Python entry points, set `PYTHONPATH="$PWD"`.
- Tests that call the `gbdraw` command need the worktree's CLI first on `PATH`:
  a shim script named `gbdraw` that runs
  `PYTHONPATH="$GBDRAW_SRC" python -c 'from gbdraw.cli import main; main()' "$@"`,
  with `GBDRAW_SRC` set to the worktree root.
- Use the port range assigned to the task for local Web servers and
  Playwright, so parallel tasks do not collide.
- Run long suites through `.claude/scripts/run_quiet.py` so the full log stays
  in a file.

## Push and finish

- Push each verified commit at once: `git push -u origin <branch>`; after a
  rebase, `--force-with-lease`. A push to a work branch starts no CI and is the
  off-machine backup. "Do not open the PR yet" never means "do not push".
- After the PR merges (GitHub deletes the remote branch), run `git status` in
  the worktree, then `git worktree remove <path>` and `git branch -D <branch>`.
  Remove only worktrees you created; ask the owner of any other one first.
- When the task ends, delete its remaining worktrees, local branches, virtual
  environments, caches, and backups. Keep the report and the handoff files.
