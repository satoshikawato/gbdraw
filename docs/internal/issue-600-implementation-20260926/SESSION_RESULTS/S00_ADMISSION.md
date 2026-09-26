# S00 follow-up — Separate trusted-base admission delivery

Status: local delivery candidates prepared; publication and dev integration pending.
The four approved outcomes are unchanged. Runtime is unchanged; S01 has not started.

## Inputs and actual remote state

- Independent clone: `/tmp/gbdraw-issue600-admission-4oPuA1`.
- Original S00 branch: `fix/issue-600-annotations-styles-20260926`.
- Original S00 SHA: `4eb3179757780b71bb8e73357f3e5f6182dcb006`.
- Actual fetched `origin/dev`: `4e84ac364beb7bea85ae0933120246944db18f39` (`B`).
- Both remote tips still matched those SHAs at candidate preparation. The original
  local/remote branch and [historical S00 report](S00.md) are preserved.
- Base inventory: contract revision 17, 39 unique active `PD-OI` records.
  All four issue #600 concerns are absent. The mapped concerns remain canonical
  request and saved-session continuity; durable `BD` decisions are empty.
- Proposed authority: exactly the original S00 contract bytes, revision 18,
  adding `PD-OI-040`–`043`. No existing record is deleted or renumbered.

## Diagnosed admission boundary

The checker was extracted unchanged from the actual base, not the candidate:

```bash
mkdir -p /tmp/gbdraw-issue600-admission-4oPuA1-evidence/base-checker
git archive 4e84ac364beb7bea85ae0933120246944db18f39 tools --output=/tmp/gbdraw-issue600-admission-4oPuA1-evidence/base-tools.tar
tar -xf /tmp/gbdraw-issue600-admission-4oPuA1-evidence/base-tools.tar -C /tmp/gbdraw-issue600-admission-4oPuA1-evidence/base-checker
node /tmp/gbdraw-issue600-admission-4oPuA1-evidence/base-checker/tools/check-web-change-budget.mjs --base 4e84ac364beb7bea85ae0933120246944db18f39 --head 4eb3179757780b71bb8e73357f3e5f6182dcb006
```

Observed exit **1**, Gate **FAIL**, Review **REQUIRED**. The sole blocking
violation is `Product Contract authority changes must be isolated from other
changed paths`. The complete output is `original-s00.log` in the artifact directory.

`check-web-change-budget.mjs` derives `diffRefs = [base, head]` and runs
`git diff --no-renames --name-status base head --`. It compares complete endpoint
trees, not one commit or a merge-base range. When the exact Product Contract path
changes, **every other changed path** is an invalid companion. The original
range includes CLI, plans, evidence and report, plus pre-existing differences in
`SELECTIVE_CI`, CI tools/tests and browser tests from its older branch base.
Splitting commits within that same PR range does not satisfy isolation.

The trusted workflow checks out the event's actual base SHA, fetches head as Git
data and executes base checker code. Its path protection does not parse static
`PD-OI` receipts or audit their Markdown links. Its `VALID` authority lines concern
mapped Product Impact JSON validation and must not be described as proof of the
four static outcomes. Receipt, active-record and link checks below supply that
separate evidence. No checker, policy, workflow or acceptance criterion is changed.

## Minimal delivery candidates and ownership

| Candidate | Parent / comparison | Changed scope |
| --- | --- | --- |
| `docs/issue-600-admission-evidence-20260926` (`D`) | Actual `B` → `D` | CLI plus the existing issue #600 plan, approvals, evidence, prompts and historical S00 report; this admission report. Product Contract stays byte-identical to `B`. |
| `docs/issue-600-contract-authority-20260926` (`A`) | Locally stacked on `D`; `D` → `A` is a sequence rehearsal | Only `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`; exact S00 payload. All authority source/acceptance documents already exist in its local parent `D`. |

Both branches are created with `--no-track` from the fetched `origin/dev`.
The authority branch then fast-forwards locally to `D` before its single contract
commit. This creates no merge commit and changes no original S00 ref. Neither
candidate is pushed in this session.

The existing [MASTER_PLAN](../MASTER_PLAN.md) remains the only implementation
plan owner; it now states docs/evidence → static authority → runtime. The original
[S00 prompt](../sessions/S00_AUTHORITY.md) is marked historical rather than left
as an active instruction to republish a mixed candidate. The original S00 report,
approvals, baseline evidence and S01–S05 prompts remain byte-identical.

CLI keeps the approved scopes and explicitly states static authority integration
and runtime implementation are pending. Its new links point to the approval
sections already delivered in `D`, avoiding links to absent `PD-OI-040/043` base
anchors. CLI is explanatory documentation, not competing Product authority.

Static Product authority remains in its existing sole contract owner. Runtime
owners, execution paths, compatibility readers and schemas are unchanged. There
is no new registry, evaluator, test store or compatibility path. Architecture
exception conditions are not triggered. Rollback of local preparation does not
supersede any approved outcome; preserve the original S00 branch and choose no
alternative Product behavior.

## Verification and evidence limits

The disposable validator lives outside the repository at
`/tmp/gbdraw-issue600-admission-4oPuA1-evidence/validate.py`. It checks:

- All nine JSON fields against both approved JSON and human `PRODUCT_DECISION`,
  including exact JSON bytes; five complete outcome clauses per record and all
  TSV/SEL/CLR/PX acceptance references.
- 43 unique active record IDs and concern/revision pairs, preserving all 39 base
  records, cross-surface clauses and the complete existing acceptance catalog.
  No mapped concern or `BD` duplication; the owner is in the base allowlist.
- Approval, historical evidence/report and untouched prompts against source S00
  bytes; master-plan acceptance definitions remain unchanged.
- Relative links and anchors in the plan package, new CLI links, and authority
  source/acceptance links against the documents delivered by `D`.
- Separate production, test, checker/workflow, generated and documentation diffs;
  candidate scope and whitespace. Final exact-commit checks are in the external
  handoff, avoiding a report added to the isolated contract candidate.

The receipt/preservation checks passed before candidate creation. Baseline runtime
observations in [EVIDENCE.md](../EVIDENCE.md) and historical S00 evidence are reused
only for their unchanged source and criteria. No TSV/SEL/CLR/PX/INT runtime
acceptance is claimed as PASS. Runtime suites, browser, wheel/build and reference
regeneration are outside this session; unchanged production/tests/generated files
are verified through endpoint diffs.

Local validation snapshots are disposable Git-index projections, not integrations
or trusted base commits. Their observed results are recorded below; they precede the final
candidate commits and contain the earlier report draft. The unchanged checker always comes from actual `B`.

| Evaluated range | Exit | Gate / Review | Meaning |
| --- | ---: | --- | --- |
| Actual `B` → docs snapshot `341cae9e43ebbd2372cd99032e8860a8d870c89f` | 0 | PASS / CLEAR | Docs/evidence alone is admissible to the fetched base. |
| Actual `B` → stacked authority snapshot `b157783397cf5cdb01c56d8c63430a920b928c5d` | 1 | FAIL / REQUIRED | Docs are not integrated yet; the full candidate still fails isolation. |
| Docs snapshot `341cae9e43ebbd2372cd99032e8860a8d870c89f` → authority snapshot `b157783397cf5cdb01c56d8c63430a920b928c5d` | 0 | PASS / REQUIRED | Contract-only local sequence rehearsal. The docs snapshot is **not trusted base** and this is **not actual admission**. |

Each command is `node /tmp/gbdraw-issue600-admission-4oPuA1-evidence/base-checker/tools/check-web-change-budget.mjs --base <table base> --head <table head>`.
Logs are `snapshot-docs-actual-base.log`, `snapshot-authority-actual-base.log`
and `snapshot-authority-sequence-rehearsal.log`. The exact final commits are checked
again by the same extracted base code; final results belong to the external handoff.
The authority snapshot's parent contains every referenced approval/acceptance
file; those files are absent from actual `B`. Checker path-only PASS by itself
therefore cannot establish this source-document integration boundary.

| Artifact | SHA-256 |
| --- | --- |
| Base `check-web-change-budget.mjs` | `c7888ffb9c6d136284b567712f6e11a474e6b08f54d4daaee7d78e1a402d7f1c` |

Extraction covers the entire base `tools/` tree so transitive checker modules
also come from the actual trusted base. Source/payload and log digests are retained
with the final handoff outside the isolated authority candidate.


Final actual candidate SHAs, endpoint diffs, checker commands, log digests and
results are recorded in `/tmp/gbdraw-issue600-admission-4oPuA1-evidence/FINAL_HANDOFF.md`.
That artifact is outside both candidates; no report is added to `A`. A candidate
SHA is not embedded as a self-reference in its own committed report.

## Publication and integration order

1. Review and explicitly authorize publication of the concrete `D` branch/SHA,
   then its PR/integration into dev. This session performs only local work.
2. After docs/evidence is actually integrated, fetch `origin/dev`. Verify its
   source documents match `D`, inspect inventory again and preserve any newer
   independent authority. Existing identical issue #600 records must not be
   re-registered; ID/revision collisions require adjusting only inventory against
   the real new base while keeping every approved field unchanged.
3. Re-evaluate the isolated contract candidate against that **actual dev SHA**
   using checker code freshly extracted from it. Its full endpoint diff must be
   exactly the Product Contract. If docs were squash/merge-integrated, reconstruct
   the contract commit on a fresh branch from that actual base rather than calling
   the unintegrated local `D` trusted or relying on the rehearsal PASS. Preserve
   every existing local/remote result; do not reset or force-push it.
4. Present the resulting concrete authority branch/SHA, one-file diff and real-base
   validation for separate publication and integration authorization. A new
   candidate SHA does not inherit permission tied to an earlier candidate.
5. After authority integration, fetch and confirm all four complete records and
   nine receipt fields in `origin/dev` exactly match
   [APPROVED_PRODUCT_DECISIONS](../APPROVED_PRODUCT_DECISIONS.md). Only that
   confirmation opens S01. No S01 work occurs here.

There is no direct push to main/dev, no PR creation, merge, deployment or tag.
The previous push authorization applies only to the original S00 candidate.
