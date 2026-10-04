# Web periodic audit and promotion checklist

The 2026-09-30 Web GUI audit
([remediation proposal, section 8.3](web-gui-audit-20260930/01_REMEDIATION_PROPOSAL.md#83-ワークフロー);
[W8, section 3](web-gui-audit-20260930/remediation/W8_verification_workflow.md))
decided that periodic audits use no new gate or workflow. Each audit is a
checklist item in the `dev` to `main` promotion pull request. This page holds the
procedure and the checklist. The admission rules for a `PROMOTION` stay in
[`WEB_CHANGE_POLICY.md`](WEB_CHANGE_POLICY.md); this page does not change them.

## When to audit

- Before each `dev` to `main` promotion.
- Or earlier, when about 30 runtime PRs have merged into `dev` since the last
  audit.

## Input

1. The capabilities changed since the last promotion. Run this from the
   repository root after `git fetch origin main dev`:

   ```bash
   git diff --name-status --no-renames origin/main...origin/dev | node --input-type=module -e "
   import { readFileSync } from 'node:fs';
   import { classifyChanges } from './tools/ci-impact-policy.mjs';
   const changes = readFileSync(0, 'utf8').trim().split('\n').map((line) => {
     const [status, ...paths] = line.split('\t');
     return { status, paths };
   });
   console.log(classifyChanges(changes).capabilities.join('\n'));"
   ```

2. The previous audit folder (`docs/internal/web-gui-audit-<date>/`) and the
   follow-ups it left open.

## Procedure

1. Dispatch `Tests` on `dev` with `tier=release` (the existing release stage;
   see [`SELECTIVE_CI.md`](SELECTIVE_CI.md)); it includes the
   `Vibrio full generation` job that otherwise only runs on `main`.
2. Run every recipe check, including the LOSATP recipes that the `recipe`
   pytest marker skips:

   ```bash
   PYTHONPATH=$PWD python docs/recipes/run_cli_scenarios.py --all --check
   PYTHONPATH=$PWD python docs/recipes/run_python_scenarios.py --all --check
   ```

   A stale artifact is regenerated with the same runner without `--check`. A
   recipe that fails before the comparison step is a product or documentation
   defect: record it as an audit row instead of regenerating.
3. Run the sweeps in [`tools/audit/`](../../tools/audit/README.md): parity
   sweep and replay, XSS sweep, viewport sweep, Gallery round trip, and session
   load timing.
4. Audit the changed areas by hand, with a fixed time limit per area. Pick the
   areas from the capabilities in the input: input, Session/History,
   Generate/options, comparisons, Feature editing, legend/Preview/export, and
   tracks/layout. Reproduce each suspected defect on `dev` before you record it.

## Carrying evidence forward

A light change can move `dev` after promotion evidence was collected on an
earlier commit E. Evidence from E counts for the promotion head H when
`node tools/ci-impact.mjs classify --base <E> --head <H>` reports the matching
verdict true ([`SELECTIVE_CI.md`](SELECTIVE_CI.md#carrying-evidence-to-a-later-commit)).
Every verdict requires E to be an ancestor of H.

| Evidence | Verdict |
| --- | --- |
| `Tests` dispatched on `dev` with `tier=release` | `releaseEvidenceCarries` |
| Gallery refresh, Gallery tutorial media, and docs GUI capture checks, and the Gallery artifact manifest | `generatedArtifactChecksCarry` |
| Recipe runs (step 2), `TestOutputComparison`, the `tools/audit/` sweeps (step 3), and the `main`-written Session fixture tests | `localTestEvidenceCarries` |

The promotion body names E, links E's evidence, and includes the `classify`
output. Rerun on H the evidence whose verdict is false. The hand audit (step 4)
names the `dev` SHA it audited and needs no carry rule. Until the planner change
that adds `classify` merges, compare `git diff --name-status <E> <H>` with the
verdict sets in `SELECTIVE_CI.md` by hand and include that comparison instead.

This rule carries checklist evidence only. `Promotion / gate` still requires
exact-SHA `Dev staging / gate` and `Gallery readiness / gate` evidence on H.
Each push to `dev` earns that evidence; after a light move, the runs inherit
their direct parent's evidence and run only the jobs the change needs
([light-change inheritance](SELECTIVE_CI.md#light-change-inheritance)).

## Output

Write `docs/internal/web-gui-audit-<date>/README.md` in the same format as
[the 2026-09-30 audit](web-gui-audit-20260930/README.md): ID, severity (P1, P2,
P3), symptom, and cause, with the commit audited. For every confirmed defect
class, record one of two decisions:

- add a guard (name the test and the PR stage that runs it); or
- add no guard, and say why.

The audit is finished when every P1 is fixed, or the Owner has given an
explicit waiver for it.

## Promotion PR checklist

Copy these items into the `PROMOTION` section of the `dev` to `main` pull
request and complete them.

```markdown
- Periodic audit: `docs/internal/web-gui-audit-<date>/` (audited `dev` SHA: <sha>)
  - [ ] Changed capabilities since `main` <sha>: <list>
  - [ ] `Tests` dispatched on `dev` with `tier=release`: <run URL>
  - [ ] `run_cli_scenarios.py --all --check` and `run_python_scenarios.py --all --check` pass, or each failure is an audit row: <rows>
  - [ ] `tools/audit/` sweeps run; findings recorded: <rows or "none">
  - [ ] Changed areas audited by hand, with the time spent per area
  - [ ] Each confirmed defect class has "guard added" (test, PR) or "no guard" (reason)
  - [ ] Every P1 is fixed, or has an Owner waiver: <links>
  - [ ] Evidence carried from an ancestor: none, or E <sha>, the `classify --base <E> --head <H>` output, and the verdicts used
```

### First promotion after the 2026-09-30 audit

That promotion carries user-visible changes that the release notes must state.
Check each one on the promotion head, or on an ancestor under
[Carrying evidence forward](#carrying-evidence-forward), and link the evidence:

- [ ] PV-08 and N-01 (D-24): Web-default Circular output changes. Long species
  lines wrap at word boundaries, and the wrap width uses an approximate glyph
  size. `definition_font_size` 18 is the default, so 18 counts as not explicit
  and may wrap; any other value is explicit and never wraps. Check the Gallery
  and the tutorial images.
- [ ] CO-05 (D-23): outfmt 6 tables with more than 12 columns (for example
  `-outfmt "6 std qlen slen"`) are read by their first 12 columns, which must
  have the right types. Columns 13 and later are dropped with an INFO log;
  `main` misread these tables silently. A CLI run given a missing or unreadable
  comparison file now fails instead of skipping it.
- [ ] Q-FRAME (D-18): comparison tables and CLI `-b` are read in the search
  frame (after record selection and crop, original strand). CLI users who
  combined reverse complement with `-b` see their coordinates mean something
  different. Main-saved Sessions are converted once on load.
- [ ] D-17: a Web PDF is 75% of its former physical size (CSS px converted to
  pt), which now matches CLI PDF and the PNG DPI.
