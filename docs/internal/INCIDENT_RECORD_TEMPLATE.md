# Incident <YYYY-MM-DD>: <what failed, in a few words>

Status: working template

Copy this file for one incident or one root-cause record: the hosted site
broke after a deploy, a release reached users with a defect, a check passed
on something it should have stopped, or work had to be recovered (as in the
`*_RECOVERY_*` documents). A periodic-audit defect class that reached `main`
gets one too. Save it as `docs/internal/INCIDENT_<YYYY-MM-DD>_<short-name>.md`,
or inside the audit or task folder that found it.

In the copy, delete this status line and these two paragraphs, and replace
every prompt in angle brackets. Write in the past tense, link commits, runs,
and pull requests, and write "unknown" instead of a guess.

- Status: <open | mitigated | closed>
- Detected: <UTC time>, by <check or run, report, or person>
- Mitigated: <UTC time>, by <rollback, fix, or other action>
- Affected: <surfaces (hosted site, PyPI package, CLI, Python API, Web app), versions or commits, and who could see it>

## What happened

<Two to four sentences that a reader who was not there can follow.>

## Timeline (UTC)

| Time | Event |
| --- | --- |
| | |

## Impact

<What users could not do or got wrong, for how long, and whether saved files,
Sessions, or published figures were affected.>

## Cause

- Trigger: <the change, commit, release, or outside event that started it>
- Root cause: <the owner that wrote the wrong value or allowed the wrong state,
  and why; name the channel: who writes the value, who reads it, and by which key>
- Why no check stopped it: <the missing, skipped, or too weak guard>

## Recovery

<What was done (rollback, fix forward, data repair), in order, and how recovery
was confirmed, with run links.>

## Guards

For each cause, one of:

- guard added: <test or check>, run by <CI job or step>, in <pull request>; or
- no guard: <reason>.

## Follow-up

| Item | Who | Link |
| --- | --- | --- |
| | | |
