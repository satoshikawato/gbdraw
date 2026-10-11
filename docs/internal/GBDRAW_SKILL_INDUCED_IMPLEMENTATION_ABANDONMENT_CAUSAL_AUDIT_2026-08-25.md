# gbdraw Skill-induced implementation abandonment causal audit

- Audit date: 2026-08-25
- Repository: `satoshikawato/gbdraw`
- Audited repository baseline: `dev` at `a0be554d296d1278f1753e1bc4baa01b3c3e5820`
- GitHub disposition cutoff: 2026-08-25 04:22 UTC, when PR #368 was closed
- Skill-install intervention used for the equal-window comparison: 2026-08-04 04:03 UTC
- Research mode: read-only for product code, branches, worktrees, stashes, reflogs, GitHub state, and session logs

## 1. Executive verdict

The evidence does not support a repository-wide increase in abandonment caused primarily by Codex Skills. In the 21 days before and after `execute-plan-with-evidence` was installed, all delivery PRs increased from 25 to 49. Merged PRs increased from 24 to 45. The raw closed-unmerged rate moved from 4.0% to 8.2%, but the samples are small, the work is not independent, and Fisher's exact test does not distinguish the windows (`p = 0.657`). A feature-branch-only view gives 4.0% versus 7.7% and `p = 1.0`.

The post-installation failures are concentrated in two exceptionally large implementation attempts rather than spread across ordinary work:

- PR #322 changed 162 files and added/deleted 36,601 lines. It expanded a local Feature-editor wiring regression into new identity, session, request, History, multi-Result, label, stroke, Gallery, and compatibility contracts. The user stopped that branch and replaced it with the focused, merged PR #323.
- PR #337 changed 75 files and added/deleted 15,999 lines. It was closed after an independent review found a second general transaction/state layer, duplicated Legend layout ownership, and rollback and async regressions despite green CI.

The literal H1 claim is contradicted: the current `execute-plan-with-evidence` text does not treat repository drift or a different failure signature as a hard stop, and archived runs repeatedly adapted to both. It did, however, make a broad named plan more salient in PR #322 and contributed to an initial "out of scope" misclassification in the unfinished PR #368 session. That is a routing risk, not evidence of a systematic hard-stop mechanism.

The strongest Skill-specific mechanism is the former trigger of `change-gbdraw-rendering-surface`. Its description matched almost any setting, Web control, label, legend, session, Gallery example, or SVG change. In PR #322's planning session, the agent explicitly used that Skill to extend the matrix through canonical requests, sessions, renderer parity, Gallery, and reference outputs. A second Gallery attempt again loaded the all-surface workflow, grew across at least 20 paths, and was rolled back at the user's direction after ownership expanded. The Skill was therefore a material amplifier, though not the primary cause in either event.

One minimal correction is applied in this audit:

- `change-gbdraw-rendering-surface` now routes only genuine shared or persisted rendering-contract changes.
- A bug confined to one surface, an internal refactor, a test/doc-only change, or a Gallery refresh with unchanged semantics no longer triggers it.
- Its lifecycle is conditional on actual owners; a surface matrix cannot create new owners or authorize scope expansion.
- A lightweight Eval specification covers local routing, shared-contract routing, plan drift, unexpected tests, extra files, partial work, unresolved UX, and size heuristics.

No product code is changed. `execute-plan-with-evidence`, `maintain-python-api`, the architecture ratchet, and active hook configuration are left unchanged because their modification thresholds were not met. The numeric Web size gate is already advisory on `dev`; `main` still uses the old hard gate until the normal `dev` promotion reaches it.

## 2. Scope and definitions

### 2.1 Evidence inspected

The audit used:

- all 278 GitHub PRs, including PR bodies, comments, reviews, checks, base/head SHAs, and branch dispositions;
- all local and remote branches, current worktrees, stashes, reflogs, `git fsck`, and commit ancestry against `main` and `dev`;
- current and historical `.agents/skills`, global Codex Skills, `AGENTS.md`, `CLAUDE.md`, Web guidance, `.claude` settings, Git hooks, GitHub workflows, PR templates, and Web policy documents;
- archived Codex sessions for the relevant implementation, planning, correction, interruption, and rollback events;
- current branch-specific behavior of the Web size checker on `dev` and `main`.

No branch, worktree, stash, reflog, hook, or generated artifact was deleted or rewritten.

### 2.2 Units

An **implementation attempt** is a PR intended to deliver code, tests, documentation, governance, or an integration change to `main` or `dev`. The report shows both all delivery PRs and a narrower feature-branch subset. The latter excludes `main`/`dev` promotion heads.

An **abandonment event** is an implementation attempt that did not reach `main` or `dev`, or a substantial implementation added and later removed, reverted, or stashed without integration. It remains an event even when the disposition was justified.

A plan-only branch is not counted as a product implementation attempt unless product code was started from it. A closed PR is not classified as abandonment when the same head or objective reached the target through a replacement PR.

Classifications are:

- A: justified abandonment, duplicate, or complete replacement;
- B: Skill or guardrail was the primary cause;
- C: Skill amplified the attempt or its abandonment risk;
- D: implementation quality or design failure;
- E: external cause;
- F: unresolved from available evidence.

### 2.3 Causal standard

A Skill is called causal only when the following align:

1. the Skill or guardrail existed before the event;
2. a concrete clause, trigger, hook, or check supplied a plausible mechanism;
3. the event log, PR record, or diff shows that mechanism operating;
4. counterevidence and successor work do not better explain the disposition.

Temporal correlation alone is not enough. A causal finding is weaker when the user explicitly requested the broader scope, implementation quality independently justified rejection, or a later interruption has no recorded reason.

### 2.4 Comparison window and limitations

The equal windows are:

- before: 2026-07-14 04:03 UTC through 2026-08-04 04:03 UTC;
- after: 2026-08-04 04:03 UTC through 2026-08-25 04:03 UTC.

The intervention date is the filesystem creation time of the global Skill, corroborated by its installation session. The Skill is not versioned, so its creation transcript and current byte-identical WSL/Windows copies are the available provenance.

The windows contain overlapping changes in branch strategy, project complexity, CI, the `dev` staging flow, numeric size policy, and the architecture ratchet. PRs also depend on one another, so the Fisher tests are descriptive checks, not randomized causal estimates.

A trend chart is intentionally omitted. Two short endpoints with several intervention dates would suggest a continuous trend that the data cannot establish; exact tables preserve the evidence without that implication.

## 3. Skill and guardrail inventory

| Candidate | Location and execution point | Current or historical mechanism | Audit finding |
|---|---|---|---|
| `execute-plan-with-evidence` | `/home/kawato/.codex/skills/execute-plan-with-evidence/SKILL.md`; byte-identical Windows copy; selected by Codex | Global `AGENTS.md` requires it for any named plan. Metadata also matches checklists, roadmaps, remediation, and migration. The body says the plan is a proposed route, compares drift, diagnoses actual failures, and limits stopping to changed outcome/authority/destructive scope or genuine blockers. | Auto-routing is broad. The alleged plan-mismatch hard stop is absent and contradicted by runs. Minor plan-salience risk; unchanged. |
| `change-gbdraw-rendering-surface` | `.agents/skills/change-gbdraw-rendering-surface/`; implicit Codex Skill | Before this audit, metadata matched individual settings, labels, legends, Web controls, sessions, Gallery examples, and SVGs, while the body requested an applicability matrix. | Material scope amplifier in two events. Trigger and conditional lifecycle narrowed in this audit. |
| `maintain-python-api` | `.agents/skills/maintain-python-api/SKILL.md` | Dirty in-scope code may be replaced when a different implementation is needed; unrelated changes stay outside the diff. Similar policy exists in `CLAUDE.md`. | Plausible loss mechanism, but no verified unauthorized loss event. Changing the Skill alone would leave the authoritative duplicate. Unchanged pending evidence. |
| `web-gallery-screenshot-maintenance` | `.agents/skills/web-gallery-screenshot-maintenance/SKILL.md` | Generator-owned outputs may be regenerated rather than merged. | Recorded Gallery removals/regenerations were user-authorized or generator-owned. No causal loss found. |
| `refactor-gbdraw-web-safely` | Historical branch-only Skill in PR #337 | Required base oracles, exhaustive ownership matrices, adversarial review, split thresholds, a Stop hook, pre-push hook, CI, and a 17-part completion report. | Protective but over-broad. Existed for about two hours inside an unmerged branch and was deleted before the branch's final commit. No active removal needed. |
| Web numeric size policy | `tools/check-web-change-budget.mjs`, `docs/internal/WEB_CHANGE_POLICY.md`, trusted-base workflow | Initially hard-failed ordinary `8/800/100` and architecture `12/1500/400` profiles. `dev` now reports size separately and exits nonzero only for integrity/architecture failures. | Hazard existed, but timing and outcomes do not support abandonment causation. `main` still has the hard form pending normal promotion. |
| Architecture fitness ratchet | `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`, checker/rules/workflows | Complete before/after owner/path sets, no self-authorization, exact-head human approval, and owner/compatibility convergence. Explicitly says file and line count are not architecture goals. | Adds review latency and an approval dependency. It is not a size gate. One unresolved event (#368) is adjacent to the approval step, but closure cause is unknown. |
| Claude/Codex hooks and settings | `.claude/settings*.json`, `.git/config`, `.githooks`, `.git/hooks`, global configs | `core.hooksPath=.githooks`, but `.githooks` is empty. `.git/hooks` contains samples only. No active Stop hook. Permission rules allowing a command are not automatic instructions. | No active hook mechanism explains abandonment or revert events. |

The repository-local Skill inventory at the audit baseline contains `avoid-ai-writing`, `change-gbdraw-rendering-surface`, `love-me-love-my-docs`, `maintain-python-api`, and `web-gallery-screenshot-maintenance`. There is no current `.claude/skills` tree and no repository `.codex/skills` content.

## 4. Skill introduction timeline

| Date and time (JST) | Commit or source | Change | Causal relevance |
|---|---|---|---|
| 2026-01-11 and earlier | `1c4f97f`, `635079a`, and earlier reverts | Explicit revert commits predate the audited Skills. | Revert behavior existed before the intervention. |
| 2026-07-05 01:19 | `b44bee88` | Initial Gallery screenshot Skill. | Early specialized workflow. |
| 2026-07-14 18:23 | `4d57953e` | Initial Python API Skill under `docs/skills`. | API maintenance workflow begins. |
| 2026-07-30 07:01 | `ec8059ce` | Adds the dirty in-scope replacement policy to the Python API Skill, Gallery Skill, and `CLAUDE.md`. | H4 mechanism predates `execute-plan`. |
| 2026-07-30 17:32 | `1cdb0d60` | Adds more cross-surface and reference-resolution guidance. | Expands contract auditing. |
| 2026-08-04 13:03 | installation session and file birth time | Installs `execute-plan-with-evidence` globally in WSL and Windows; global `AGENTS.md` auto-routes named plans. | Main before/after intervention. |
| 2026-08-04 13:06 | `3d7c82da` | Adds `change-gbdraw-rendering-surface`; moves maintained Skills to `.agents/skills`. | H2 trigger begins. |
| 2026-08-13 | PR #322 / `30cfb387` | Large Feature-style plan implementation closed; narrow PR #323 succeeds. | First strong scope-amplification event. |
| 2026-08-16 10:38 | `5accc0fc` | Adds `refactor-gbdraw-web-safely`, Stop/pre-push/CI guards, templates, and a cross-trigger. | Transient maximum guardrail load in PR #337. |
| 2026-08-16 12:38 | `9bd55935` | Deletes that Skill, hooks, workflow, template, and bundle plans while continuing product rework. | Guardrail lifetime: 2h00m32s. |
| 2026-08-16 18:19 | `fdfa378d` | First hard Web numeric budget. | Occurs after PRs #322 and #337; cannot cause them. |
| 2026-08-16 22:54 | `915c2226` | Finite hard ordinary and architecture profiles. | H3 intervention. |
| 2026-08-17 to 2026-08-18 | `bf2c0517`, `9378f20e` | Makes `dev` a first-class target and adds trusted-base handling. | Changes PR denominator and routing. |
| 2026-08-18 22:35 onward | `79638639` through `1909ada8` | Defines and progressively enforces the architecture ratchet. | Adds owner/approval gates after #322/#337. |
| 2026-08-25 10:25 | `038a1049` | Converts numeric Web size limits to review-only on `dev`. | Direct H3 remediation. |
| 2026-08-25 10:44 | `a0be554d` / PR #367 | Merges advisory size behavior to `dev`. | Audit baseline. |

The current `change-gbdraw-rendering-surface` file was otherwise byte-equivalent to its initial `3d7c82da` version before this audit. There was no gradual trigger drift to explain the August events.

## 5. Implementation-attempt and abandonment-event ledger

### 5.1 Complete closed-unmerged PR screen

GitHub contains 278 PRs: 260 merged, 18 closed without merge, and no open PRs. The following table screens every closed-unmerged PR rather than treating only recent examples as the denominator.

| Event ID | Date | Branch / PR / commit | Intended objective | Work completed | Final disposition and reason | Relevant Skill or guardrail | Evidence | Class / confidence |
|---|---|---|---|---|---|---|---|---|
| GH-014 | 2025-07-17 | PR #14, `main` / `947c5671` | Reverse synchronization to Streamlit | Head was already on `main` | Closed duplicate/reverse-sync PR | None of the audited Skills | Head already reachable from `main` | A / high |
| GH-015 | 2025-07-17 | PR #15, `alert-autofix-1` | Code-scanning command fix | 1 file, `+7/-2` | Autofix proposal closed; manual command fix reached PR #19 | Pre-Skill | Successor implementation | A / medium |
| GH-016 | 2025-07-17 | PR #16, `alert-autofix-2` | Code-scanning command fix | 1 file, `+7/-17` | Autofix proposal closed; manual command fix reached PR #19 | Pre-Skill | Successor implementation | A / medium |
| GH-017 | 2025-07-17 | PR #17, `alert-autofix-3` | Code-scanning path fix | 1 file, `+10/-1` | Autofix proposal closed; manual fix reached PR #19 | Pre-Skill | Successor implementation | A / medium |
| GH-018 | 2025-07-18 | PR #18, `streamlit` / `b6c301f6` | Streamlit update | 5 files, `+368/-10` | Exact head included in merged PR #19 | Pre-Skill | Commit-to-PR ancestry | A / high |
| GH-036 | 2025-09-04 | PR #36, `streamlit` / `73d7a47a` | Streamlit-to-main implementation | 4 files, `+84/-89` | Closed; branch later deleted; no reason or exact successor proof | Pre-Skill | GitHub disposition only | F / medium |
| GH-039 | 2025-09-08 | PR #39, `label_whitelist` / `2134fe5d` | Label whitelist | 19 files, `+449551/-231` | Same head merged 64 seconds later through PR #40 | Pre-Skill | Exact head match | A / high |
| GH-045 | 2025-10-08 | PR #45, `GCskew_linear` | GC-skew Linear change | 12 files, `+346/-180` | Core objective already merged in PR #44; later extras discarded | Pre-Skill | Timing and overlap | A / medium |
| GH-078 | 2025-10-26 | PR #78, `dev` / `2dfb0b0f` | Development integration | 12 files, `+1400/-1279` | Same head already merged in PR #77 | Pre-Skill | Exact head match | A / high |
| GH-166 | 2026-03-12 | PR #166, `custom_tracks` | Custom tracks | 36 files, `+5345/-613` | Closed after 21 days; a later independent implementation merged as PR #222, but no explicit link | Pre-Skill | Possible successor only | F / medium |
| GH-178 | 2026-04-06 | PR #178, `gui_improvement` | Label displacement and wheel bundle | 21 files, `+752/-152` | Overlapping PR #179 opened and merged before closure | Pre-Skill | Successor timing and overlap | A / high |
| GH-215 | 2026-05-10 | PR #215, `mobile` | gbdraw icon | 6 files, `+98/-2` | Icon reimplemented in merged PR #218 | Pre-Skill | Successor implementation | A / medium |
| GH-235 | 2026-06-07 | PR #235, `cli_table_arguments` | CLI table arguments | 33 files, `+5833/-211` | Closed after 23 days; later TSV PR series likely superseded it, without an explicit closure link | Pre-Skill | Possible successor only | F / medium |
| GH-309 | 2026-07-20 | PR #309, `2026-07-20` / `520ba228` | Automatic region highlights | 3 commits; branch later advanced | Commit became an ancestor of merged PR #312 on the same head branch | Pre-Skill equal-window control | Ancestry and same branch | A / high |
| GH-322 | 2026-08-13 | PR #322, `bugfix-20260813` / `30cfb387` | Restore Feature-style scope editing and atomic legend behavior | 162 files, `+31721/-4880`; new schema/session/identity/command systems; several focused gates passed, broad gates remained red | User explicitly archived the branch and prohibited merge/cherry-pick; focused PR #323 opened before closure and later merged | Broad rendering-surface trigger, named plan, automatic execute-plan routing | Planning, implementation, and stop sessions; failed Browser/Recipes/Playwright checks; 11-path overlap with #323 | C, with D/A / high |
| GH-337 | 2026-08-15 | PR #337, `GBDRAW_WEB_MATURITY...` / `476fb10a` | Centralize live SVG edit commit and replay | 75 files, `+13712/-2287`; green reported CI; four commits | Closed with explicit comment citing a second transaction/state layer, duplicated Legend ownership, rollback/async/maintainability risk | Named plan; transient refactor Skill/hooks/CI; not numeric size policy | Closure comment and independent review with six concrete defects | D / high |
| GH-343 | 2026-08-17 | PR #343, `dev` / `6f084dca` | Promote an already merged integration | Zero effective diff | Closed after 15 seconds as redundant; head reaches both targets | Delivery workflow | Zero diff and ancestry | A / high |
| GH-368 | 2026-08-25 | PR #368, `fix/linear-arrange-in-rows` / `2ebf26c8` | Restore Linear shared-row authoring and row-label semantics | 8 files, `+1086/-43`; all completed checks passed | Owner closed without a reason after the user corrected one scope assumption; branch and clean worktree retained | Execute-plan auto-route; architecture approval field was blank; numeric size check passed | GitHub checks and local session timestamps | F / high disposition, medium causal |

Screening totals are A = 12, D = 1, C = 1, F = 4, B = 0, E = 0. A GitHub-only classification would put GH-322 in D because the Skill link comes from local session evidence rather than the PR record.

### 5.2 Revert, local WIP, stash, and plan-only events

| Event ID | Date | Branch / PR / commit | Intended objective | Work completed | Final disposition and reason | Relevant Skill or guardrail | Evidence | Class / confidence |
|---|---|---|---|---|---|---|---|---|
| LOC-001 | 2026-01-14 | `bugfix_webapp`, `8ee6c44d` | Historical Web bug fix | Commit subject/body `abadoned` | Local branch retained; no successor proof in this audit | Predates audited Skills | Local branch and commit | F / low |
| LOC-002 | 2026-04-29 | `cancel_mid_run`, `2098512f` | Historical cancellation work | Commit `under construction`; remote gone | Local WIP retained | Predates audited Skills | Local branch and stash history | F / low |
| REV-324 | 2026-08-13 | PR #324, `bc18a3cf` then revert `376da0f1` | Adopt current sessions with lazy resources | 13-file optimization added and reverted 30 minutes later | User instructions explicitly required preserving a follow-up branch and non-destructively removing this separate optimization from PR #324; equivalent capability landed through PR #329 | Execute-plan used, but instruction came from the user | Attachment, session, revert commit, later PR | A / high |
| LOC-003 | 2026-08-17 | `gallery-parity-implementation`, uncommitted WIP | Gallery first-Generate parity across examples | At least 20 paths edited; plan-required surface work plus a self-added XML/CSS parser and dependency | User challenged ownership growth and ordered return to the no-new-owner state; initial plan commit reset and revised plan recommitted | Execute-plan, rendering-surface Skill, broad named plan | Session log, patch log, reflog | C, with D / high |
| LOC-004 | 2026-08-24 to 25 | stash `3ca6fb60` on `gallery-replay-tutorial-media` | Gallery PR3 tutorial/media publication | 52 files, `+1023/-7977`; captures, tutorial JSON, tests, and user-authorized WSSV removal | Turn interrupted; work later stashed, not deleted; no PR | Three documentation/Gallery Skills, but interruption and user scope change dominate | Session, explicit deletion authorization, stash | F / high event, low Skill causality |
| LOC-005 | 2026-08-25 | `/tmp/gbdraw-pr1`, `fix/linear-arrange-in-rows` | Correct PR #368 after user clarified one-file/multi-record row behavior | Clean at `2ebf26c8`; promised correction did not begin | Both root and worker turns interrupted; PR later closed without comment | Execute-plan scope salience; architecture approval pending | Root/worker session and worktree | F / high event, medium Skill causality |
| PLAN-001 | 2026-08-17 to 22 | `gallery-parity-implementation`, `revise-gallery-session-replay-plan` | Architecture and Gallery planning | 5 docs-only commits / 2 files / `+3592`; 1 docs-only commit / 8 files / `+1839` | Plans were referenced or superseded; implementation PRs #364/#365 merged | Execute-plan appears in plan text | Branch diffs and successor PRs | A / high; excluded from product-attempt rate |

Other current worktrees do not add abandonment events. `/tmp/gbdraw-browser-guard-fix` is a prunable worktree for merged PR #366, and the detached Gallery baseline is a disposable read-only baseline. `git fsck --unreachable --no-reflogs` found no post-2026-08-04 dangling implementation commits.

## 6. Before/after comparison

### 6.1 Equal 21-day cohorts

| Metric | Before | After | Interpretation |
|---|---:|---:|---|
| All delivery PR attempts | 25 | 49 | Activity increased 96%. |
| All merged | 24 | 45 | Completed volume increased 88%. |
| All raw closed-unmerged | 1 | 4 | The before closure was superseded; the after set contains one redundant sync, two design/scope events, and one unresolved event. |
| All raw completion rate | 96.0% | 91.8% | Difference is 4.2 percentage points; Fisher `p = 0.657`. |
| Feature-branch attempts | 25 | 39 | Activity increased 56% after excluding promotion heads. |
| Feature-branch merged | 24 | 36 | Completed feature volume increased 50%. |
| Feature-branch raw closed-unmerged | 1 | 3 | PRs #322, #337, and #368. |
| Feature-branch raw completion rate | 96.0% | 92.3% | Fisher `p = 1.0`; no detectable cohort difference. |
| Explicit implementation revert inside a cohort PR | 0 / 25 (0.0%) | 1 / 49 (2.0%) | The one revert was an authorized scope split inside merged PR #324. |
| Median feature PR merge latency | 0.27 h | 0.53 h | Both medians are short; branching/review policy changed between windows. |
| Median changed files per feature PR | 13 | 5 | Typical post-window PRs were smaller, not larger. |

The post-window median hides two extreme abandoned attempts. PR #322 has 162 changed files and PR #337 has 75, while the post-window median is 5. This concentration supports a plan/design scope mechanism more strongly than a general decline in delivery capability.

The pre-window churn mean is not useful: PR #298 reports more than 1.35 million additions because its reused branch/diff topology makes the API statistic incomparable. Medians and adjudicated event histories are therefore preferred.

### 6.2 Hard-size-policy period

From PR #340 through PR #366, while the hard numeric policy was active, 27 PRs closed: 26 merged and only PR #343 did not. PR #343 had zero effective diff and was a redundant integration. The adjudicated completion rate is therefore 100% for substantive attempts, with no explicit implementation revert.

PR #365 landed at 1,499 of a 1,500-line architecture-profile threshold. That is evidence of threshold pressure, but no evidence shows necessary behavior being deleted, unnatural compression, or abandonment. The policy was nevertheless correctly demoted to an advisory on `dev` because such pressure can distort future changes.

### 6.3 Architecture-ratchet period

From PR #353 through PR #368, 15 of 16 PRs merged. The sole unresolved disposition is PR #368. Its architecture decision permalink remained blank, but all automated checks, including the advisory size check and trusted-base check, passed. Because the owner closed it without a comment after an interrupted correction session, the audit cannot distinguish approval friction, behavior disagreement, timing, or another reason.

### 6.4 Plan revisions and PR count per feature

The repository does not provide a reliable single denominator for plan revisions. Several plans remained untracked during execution, some were recommitted after a reset, and some were referenced by fixed SHA without merging. The audit therefore reports concrete plan histories rather than a synthetic rate:

- PR #322's plan grew from 967 to 1,268 lines and its implementation reached 162 files before a focused replacement PR was requested.
- The Gallery parity plan generated five docs-only commits, a reset/recommit, and a later eight-document plan before PRs #364 and #365 delivered narrower units.
- The architecture ratchet intentionally used several governance PRs; those are delivery design, not repeated failed product attempts.

## 7. Candidate-by-candidate causal assessment

| Hypothesis or candidate | Verdict | Supporting evidence | Counterevidence / alternative | Confidence |
|---|---|---|---|---|
| H1: `execute-plan-with-evidence` turns plan/repository mismatch into a hard stop | Not supported as stated; broad auto-routing is a minor contributor | Global `AGENTS.md` auto-selects it for any named plan. In PR #368, the agent initially treated a valid one-file/multi-record behavior as outside the planned journey. | Current Skill text calls the plan a proposed route, says to compare drift, diagnose actual failures, and not stop for long work, a dirty tree, or sandbox approval. PR #322 adapted a stale fixture assumption; PR #367 treated dirty-plan drift as non-blocking; many named-plan PRs merged. | High for no hard-stop text; medium for minor plan salience |
| H2: the all-surface Skill broadens local work until it is incomplete | Material contributor, not primary | PR #322 planning explicitly extended a local wiring regression through request/session/renderer/Gallery/reference surfaces under this Skill. Gallery WIP again loaded the workflow and crossed at least 20 paths before owner-growth rollback. | The user requested structural prevention and detailed plans. The Skill body prohibited duplicate renderers, and the Gallery parser duplication was an implementation error. | Medium-high |
| H3: numeric PR limits cause deletion, compression, or abandonment | Not supported for observed abandonment; protective but previously over-broad | Hard limits existed and PR #365 sat one line below a threshold. `main` still has the hard version. | The policy postdates #322/#337; the hard-policy substantive completion rate was 100%; PR #368 ran after `dev` became advisory and its size check passed. | High on timing, medium on future risk |
| H4: dirty in-scope replacement rules destroyed partial work | Insufficient evidence | `maintain-python-api` and `CLAUDE.md` contain a concrete replacement mechanism. | Every material revert/restore/stash located was user-authorized, generator-owned, or retained. No post-Skill dangling commit was found. No session proves unauthorized loss. | Medium |
| H5: stop conditions embedded in plans are the primary cause | Literal hard-stop claim not supported; over-broad plan scope is a material contributor | The Feature-style and Gallery plans prescribed large owner/surface matrices. The architecture plan contains many hard-stop clauses. | No abandonment event was triggered by a located architecture-plan stop clause; the architecture sequence completed. PR #322 was stopped by the user because of scope and remaining regressions, not by a plan stop. | High for plan-scope contribution; medium for literal negative |
| H6: work volume and difficulty changed after Skill installation | Material confound, but not a complete explanation | All PR attempts rose 96%. The two definite failures were exceptionally large and architecturally difficult. CI, `dev`, and ratchet workflows changed simultaneously. | Typical post-window PR size was smaller, and the exceptional breadth was itself partly plan-induced. | High |
| Historical `refactor-gbdraw-web-safely` | Protective but over-broad; process/churn contributor | Added a 463-line workflow plus Stop/pre-push/CI/template requirements inside PR #337; all were deleted two hours later. | Independent correctness findings justified rejecting PR #337. The Skill caught defects and never reached `main` or `dev`. | Medium |
| Architecture ratchet | Protective with approval overhead | Exact-head human decision and renewal after later commits can add latency; PR #368 had no decision permalink. | It rejects raw line/file minimization, and 15 of 16 PRs in its short window merged. #368 has no recorded closure reason. | Low-medium for #368 causation |

## 8. Primary mechanisms

The supported failure pattern is scope amplification, not an automatic stop:

```text
local or bounded defect
    -> broad Skill routing or broad plan inventory
    -> plan gains schema, compatibility, persistence, Gallery, and all-mode work
    -> execute-plan faithfully treats the named plan as the active work package
    -> implementation adds owners and interacting state machines
    -> real-browser or adversarial review finds regressions/duplication
    -> user requests a narrow rewrite or rejects the branch
```

### 8.1 Surface-matrix completion became scope

The former rendering-surface metadata made an individual Web control, label, legend, Gallery example, or SVG sufficient for selection. Its body said to omit a surface only after confirming non-applicability. In PR #322, the planning report explicitly cited that check when adding canonical request, session, renderer parity, Gallery, and reference-output work.

The audit does not claim that every added surface was unnecessary. It finds that the routing rule lacked a first gate separating "the shared contract is changing" from "one surface is incorrectly wired to an already-correct contract." Without that gate, an applicability audit could become an implementation inventory.

### 8.2 A named plan gained authority through routing

The global rule automatically selects `execute-plan-with-evidence` for any named plan. The Skill itself is adaptive, but automatic selection increases the salience of the plan's scope and phase boundaries. PR #368 shows the practical risk: the agent first defended an existing validation as outside PR1, then acknowledged that it had mistaken current code and plan wording for product intent after the user clarified the behavior.

### 8.3 Temporary duplicated guardrails increased process and code churn

PR #337 briefly carried the same refactor guard through a Skill, Claude Stop hook, pre-push hook, CI workflow, PR template, static tests, and prose. Those artifacts were removed two hours later. The final branch still failed on substantive architecture and runtime behavior, so duplicated enforcement was not the primary cause, but it added thousands of changed lines and another review surface.

### 8.4 Numeric size limits created pressure without a located failure

The hard checker combined size violations with integrity violations and exited nonzero. That mechanism could encourage splitting, deletion, or threshold gaming. The observed period does not show an abandonment caused by it, and `dev` now reports size separately. The residual risk is branch-specific: trusted-base checks for `main` still execute the hard checker from the `main` base until promotion.

### 8.5 No active automatic Stop or destructive hook exists

The historical Stop hook blocked a first stop when Web refactor guards failed and avoided an infinite loop on a repeated stop. It existed only in PR #337. The audit baseline has an empty `.githooks`, no `.claude/hooks` files, sample-only `.git/hooks`, and no Codex/Claude configuration key that automatically stops, reverts, restores, or cleans work.

## 9. Counterevidence and uncertainty

Several facts prevent a stronger "Skills caused abandonment" conclusion:

- 45 of 49 all delivery PRs merged in the post-installation window.
- There are no open PRs, so no stale-open backlog was reclassified as failure.
- Of 42 retained non-default GitHub branches, 37 already reach `main` or `dev`. The five that do not are PRs #322/#337/#368 and two docs-only planning branches.
- `execute-plan-with-evidence` adapted when PR #322's Gallery fixture differed from the plan, when PR #367 began with unrelated untracked plans, when CI failures had a different cause, and during several merged architecture and Gallery phases.
- PR #324's revert was explicitly requested, preserved on a follow-up branch, and later reimplemented successfully.
- PR #337 has a sufficient non-Skill explanation: independent review found concrete correctness and ownership defects.
- The Gallery WIP contained a self-added XML/CSS engine that directly violated the repository's convergence rules; the Skill did not require that implementation.
- No evidence links the H4 dirty-replacement clauses to an unauthorized loss.
- PR #368's closure has no comment or review. A blank architecture decision field is adjacent evidence, not a recorded reason.
- The 52-file PR3 stash followed a turn interruption and an explicit user request to remove WSSV publication material. It is suspended work, not an unauthorized revert.
- Two historical local WIP markers and several explicit reverts predate the August Skills.

The session archive is unusually rich but not a complete causal instrument. An agent can follow a Skill without naming every clause that influenced a decision, and a user can close a PR outside the recorded conversation. Confidence therefore refers to the evidence available, not to the impossibility of another explanation.

## 10. Recommended remediation

### Applied now

1. Narrow `change-gbdraw-rendering-surface` to intentional shared or persisted contract changes.
2. Explicitly exclude single-surface bugs, internal refactors, test/doc-only work, and content refreshes with unchanged semantics.
3. State that multiple changed files do not make a bug cross-surface.
4. Make each lifecycle step conditional on an actual owner and prohibit adding schema, migrations, APIs, Gallery assets, or references merely to complete a matrix.
5. Protect a partially correct path from removal solely to match a planned file set or size target.
6. Add lightweight routing and continuation Evals without introducing a test framework.

### Keep unchanged now

- Keep `execute-plan-with-evidence` unchanged. Its broad automatic route should be observed, but its stop text is already adaptive and the audit did not find multiple hard-stop events.
- Keep the current architecture ratchet. It addresses owner/path debt and explicitly rejects raw size as the architectural objective.
- Keep `maintain-python-api`, Gallery replacement guidance, and the duplicated `CLAUDE.md` rule unchanged until an actual unauthorized-loss event is captured. A speculative edit to one copy would create inconsistent policy.
- Keep hook state unchanged; there is no active Stop/pre-push hook to remove.

### Complete through normal delivery

- Promote the advisory numeric size behavior from `dev` to `main` through the repository's normal staging/release flow. Do not bypass the branch policy for this audit.

### Add minimal observability before the next causal audit

For future agent-authored attempts, record in the PR body or handoff:

- Skills actually selected and why;
- whether the requested behavior or only an implementation plan changed;
- each plan-drift decision and whether work continued;
- any surface added beyond the initial defect and the contract that requires it;
- interruption, user-requested stop, supersession, or close reason;
- whether retained work exists as a branch, worktree, or stash.

Re-audit after 30 additional feature PRs or after two new plausible B/C events, whichever occurs first. A second event directly showing `execute-plan-with-evidence` treating a feasible drift as terminal would justify narrowing its global auto-trigger.

## 11. Exact files and clauses to modify

### 11.1 Applied Skill changes

| File and clause | Before | After | Why |
|---|---|---|---|
| `.agents/skills/change-gbdraw-rendering-surface/SKILL.md`, frontmatter description | Any setting, track, label, legend, geometry, API, CLI flag, Web control, session, Gallery example, or SVG could trigger the Skill. | Trigger requires an intentional shared or persisted rendering-contract change; local bugs/refactors/content-only work are excluded. | Prevent recurrence of the PR #322 local-bug-to-all-surface route. |
| Same file, new applicability section | No precondition distinguished shared contract change from one-surface miswiring. | Requires shared/persisted semantics, says file count does not define applicability, and forbids creating surface owners to fill the matrix. | Keep adjacent-boundary inspection without turning it into scope. |
| Same file, lifecycle | Typed model, API/CLI, Web, persistence, docs/Gallery, and cleanup read as one default sequence. | Each step is conditional on whether that owner exposes or persists the changed contract. | Preserve parity where real while avoiding artificial owners. |
| Same file, removal rule | Retired paths were removed in the lifecycle. | Superseded paths are removed only as part of the complete implementation; partially correct work is not removed to match a planned file set or size target. | Preserve valid work while retaining convergence. |
| `.agents/skills/change-gbdraw-rendering-surface/agents/openai.yaml` | UI text promised work across every public surface and APIs/CLI/Web/sessions/tests. | UI text says shared contracts and only the owning/persisted surfaces that apply. | Metadata is the first routing boundary and must match the Skill. |
| `.agents/skills/change-gbdraw-rendering-surface/EVALS.md` | No Skill-local Eval. | Adds two routing controls and six continuation controls. | Covers the failure mechanisms without a large framework. |

The key routing diff is:

```diff
-description: Implement or review a gbdraw rendering option across every surface that owns or persists it. Use when adding, changing, renaming, or removing diagram settings, track behavior, labels, legends, geometry, comparison data, output topology, Python APIs, CLI flags, Web controls, typed render requests, saved sessions, Gallery examples, or reference SVGs.
+description: Implement or review a gbdraw rendering contract that intentionally changes multiple owning or persistence surfaces, such as shared option semantics, Python/CLI/Web request parity, session serialization, or renderer/output topology. Do not invoke for a bug confined to one surface, an internal refactor, a test or documentation-only change, or a Gallery content refresh when the shared contract is unchanged.

-Before editing, trace the concept through the surfaces it actually uses:
+After confirming that the Skill applies, trace the concept through the surfaces it actually uses:

-Omit a surface only after confirming that the concept does not apply there.
+Mark non-applicable surfaces explicitly; do not create an owner merely to fill the matrix.

-default_prompt: "Use $change-gbdraw-rendering-surface to implement this rendering change across gbdraw APIs, CLI, Web, sessions, and tests."
+default_prompt: "Use $change-gbdraw-rendering-surface to update this shared rendering contract across the owning and persisted surfaces that actually apply."
```

No product file, test product contract, generated artifact, user branch, worktree, or stash is changed.

### 11.2 Clauses considered but not changed

| Candidate | Exact clause | Proposed change if future evidence crosses the threshold | Why not now |
|---|---|---|---|
| Global execute-plan routing | `/home/kawato/.codex/AGENTS.md`: use the Skill for a named implementation plan; Skill description matches plans/checklists/roadmaps/remediation/migration | Require explicit user invocation or explicit designation of a reviewed plan as execution-authoritative. State that a named plan alone does not trigger it. Update both WSL and Windows copies together. | Broad routing is real, but multiple terminal plan-mismatch events were not found. Current body already instructs adaptation. |
| Dirty in-scope replacement | `CLAUDE.md` dirty-tree section; `maintain-python-api` dirty-tree section; Gallery generator guidance | Protect user/other-agent partial work by default; integrate it or ask only when an outcome conflict cannot be resolved. Update all authoritative duplicates together. | No verified unauthorized-loss event. |
| Numeric size enforcement on `main` | `main:tools/check-web-change-budget.mjs` combines budget and integrity violations | Promote the already merged `dev` advisory behavior. | The fix exists and is validated on `dev`; this audit does not authorize bypassing the normal promotion path. |
| Architecture plan hard-stop list | `docs/internal/GBDRAW_ARCHITECTURE_FITNESS_FUNCTION_RATCHET_IMPLEMENTATION_PLAN_2026-08-17.md`, stop-conditions section | Reword repository drift and file-set changes as reconciliation events unless behavior/authority/safety becomes unresolved. | No located abandonment was triggered by those clauses; the plan sequence completed. |

## 12. Skill Eval scenarios

The new `.agents/skills/change-gbdraw-rendering-surface/EVALS.md` defines observable decisions rather than matching prose.

| Scenario | Expected result |
|---|---|
| Local Feature picker calls a raw setter while all shared semantics are already correct | Do not invoke the Skill; repair the Web owner and focused tests only. |
| New persisted label-layout option spans typed request, API, CLI, Web, session, and two render modes | Invoke the Skill; map and update the real owners. |
| Plan names function A, current owner is B | Record drift and continue at B. |
| Expected test X passes but Y demonstrates the defect | Replace the diagnostic hypothesis and continue toward the same behavior. |
| A fourth file is necessary beyond the planned three | Explain and change it; do not delete behavior to preserve file count. |
| User or another agent has a partial implementation | Inspect and integrate it; do not reset/restore/checkout/clean without explicit authorization. |
| Two reasonable UX choices remain unresolved | Ask for the exact decision; preserve non-conflicting progress. |
| Complete fix exceeds a size guideline | Consider a cohesive split; retain the complete change if splitting creates an incomplete owner interval. |

Validation performed:

- Skill frontmatter and structure: `quick_validate.py` passed.
- `agents/openai.yaml`: parsed successfully and remains consistent with the Skill.
- `git diff --check` for the Skill directory: passed.
- Independent forward routing: an evaluator made both judgments before reading `EVALS.md`. It excluded the local Feature-picker wiring defect and invoked the Skill for the shared persisted label-layout option; its owner map also rejected automatic schema migration, Gallery regeneration, and duplicate mode-specific owners. These judgments matched the Eval expectations.
- Runtime metadata caveat: this already-running Codex session retains the pre-edit Skill description in its loaded catalog. The on-disk trigger is correct, but automatic routing must be reloaded in a later session before treating the new metadata as active.

## 13. Residual risks

1. `main` still hard-enforces numeric Web size profiles. A `main`-target PR can receive different behavior from the same change targeting `dev`.
2. Global `execute-plan-with-evidence` selection remains broad. Although its body is adaptive, future agents can still treat a detailed plan as more authoritative than the user outcome.
3. The dirty in-scope replacement rule remains in several authoritative locations. No event justified changing it, but the textual mechanism is real.
4. Narrower rendering-surface routing can create a false negative if a seemingly local bug actually changes persisted semantics. The revised Skill handles this by allowing escalation after concrete shared-contract evidence is found.
5. PR #368 remains causally unresolved. Its branch and clean worktree preserve the implementation, but the correction requested at 04:00 UTC is absent.
6. The large PR3 stash is retained but has no PR. Future work must distinguish intentional WSSV removal from incomplete remaining tutorial publication.
7. The two docs-only planning branches can be mistaken for failed implementation unless product-code ancestry is checked.
8. The comparison window is short and contains several governance interventions. Rate estimates should not be extrapolated.
9. Closed-PR reason is optional. Four historical events remain F because no reliable reason or exact successor was recorded.
10. The current Codex turn loaded the old rendering-surface description before the edit. A fresh metadata load is needed for the narrowed automatic trigger to take effect operationally.

## 14. Final implementation-readiness decision

```text
Primary cause:
Over-broad implementation-plan scope combined with implementation complexity and duplicated ownership in the two definite post-installation failures. No Skill hard stop is established as the primary cause.

Material contributors:
- The former broad trigger and quasi-mandatory surface matrix in change-gbdraw-rendering-surface.
- Plan scope growth in the Feature-style and Gallery parity work.
- Implementation-specific design errors in PRs #322 and #337.

Minor contributors:
- Automatic execute-plan routing, which increases plan salience even though the Skill text is adaptive.
- Transient duplicated refactor guardrails in PR #337.
- Possible architecture-approval friction around PR #368; the cause is unrecorded.

Not supported:
- A statistically or operationally established repository-wide increase in Skill-caused abandonment.
- execute-plan-with-evidence treating ordinary repository drift or a different failure signature as a hard stop.
- Numeric size limits causing the observed PR #322, #337, or #368 dispositions.
- Verified unauthorized loss from dirty in-scope replacement rules.
- Any current active Stop, pre-push, or destructive hook causing abandonment.

Still uncertain:
- Why PR #368 was closed after all checks passed.
- Final intended disposition of the retained PR3 stash.
- Four older closed-unmerged PRs without exact successor or close rationale.

Recommended action:
- Keep unchanged: execute-plan-with-evidence, maintain-python-api, current architecture ratchet, and inactive hook state.
- Narrow auto-trigger: applied to change-gbdraw-rendering-surface through discriminating metadata and an applicability gate.
- Rewrite stop conditions: not applied; no evidenced hard-stop clause caused an event.
- Convert hard limit to review heuristic: already applied on dev; promote normally to main.
- Add Eval: applied in the rendering-surface Skill directory.
- Deprecate: historical refactor-gbdraw-web-safely is already absent from integrated branches.
- Remove: nothing further.

Implementation readiness:
Safe to apply now for the repository-local Skill correction. Re-audit is required before changing global execute-plan routing or dirty-work replacement policy. Main's size-policy correction should arrive through the normal dev-to-main promotion.
```

The applied correction raises the chance of completing valuable work by changing the first routing decision: a local defect stays local unless evidence shows that the shared contract itself is wrong. It does not weaken safety stops, persistence compatibility, owner convergence, regression evidence, or explicit user authority.
