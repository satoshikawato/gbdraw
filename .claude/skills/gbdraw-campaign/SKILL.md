---
name: gbdraw-campaign
description: Keep the state of a multi-session gbdraw task in files - start, resume, hand off, report progress, and finish a campaign with STATE.md, LOG.md, DECISIONS.md, and progress.json in the task folder. Use when starting or resuming work that spans sessions or PRs, before compacting or ending such a session, and when asked how far the work has progressed.
---

# Campaign state files

A campaign is work that spans sessions, PRs, or agents. Its files live in one
task folder outside the repository (`$TASK_DIR`; the local instructions say
where). The files carry the state, not the conversation. This skill does not
decide the campaign's scope.

| File | Content | Rule |
| --- | --- | --- |
| `STATE.md` | Current state only: goal and done criteria, open PRs (number, branch, head SHA, CI), worktrees, the next steps, open questions, assigned ID blocks | Overwrite; at most 150 lines. Read it first after a resume or compaction. |
| `LOG.md` | What happened, in time order | Append only. Search it (`grep`, `tail`); never read it whole. |
| `DECISIONS.md` | Each Owner decision and each Owner-delegated recommendation: id, date, question, choice, source | Append only. Cite the ids in PR bodies and the final report. |
| `progress.json` | Work items with weight and status | Update when an item changes. |

Templates are in `templates/`.

## Start or resume

1. If the folder is new, copy `templates/` into it and fill `STATE.md` and
   `progress.json`.
2. Bind this session to the folder, so that the SessionStart hook reloads
   `STATE.md` after a compaction or a resume:
   `python3 .claude/skills/gbdraw-campaign/scripts/bind_session.py "$TASK_DIR"`.
3. Read `STATE.md` and `DECISIONS.md`. Search `LOG.md` only for what you need.

## While working

- After each merge, abandoned approach, or decision, update `STATE.md` and
  append to `LOG.md` or `DECISIONS.md` in the same turn.
- Record each newly found bug at once under the next ID of the campaign's
  block, with steps, expected result, actual result, and cause.
- Answer "how far along is it?" by running
  `python3 .claude/skills/gbdraw-campaign/scripts/progress.py "$TASK_DIR"`
  (`--gh` also reads PR states from GitHub), not from memory.

## Split the session

End the session and continue in a new one when any of these holds:

- a PR merged and the next item touches other files;
- the context is past about 300k tokens;
- a second compaction would be needed;
- the session's goal is done;
- the session sat idle for more than an hour, so the prompt cache is gone.

Before ending, update `STATE.md` and write `NEXT-SESSION-PROMPT.md` from the
template: under 10 KB, pointing to the files above instead of copying them.
Give the Owner its path.

## Finish

Write the final report from `DECISIONS.md`, `progress.json`, and `LOG.md`.
Mark the campaign complete in `STATE.md`, and remove the session bindings with
`bind_session.py --clear "$TASK_DIR"`.
