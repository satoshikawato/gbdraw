# Issue 598 — Accepted Product Decisions

Status: ACCEPTED by the Product Decision Owner; admission as runtime base authority requires merge into `dev`.

- Product Decision Owner / author: `satoshikawato`
- Decision date: 2026-09-26
- Issue: https://github.com/satoshikawato/gbdraw/issues/598
- Prepared against dev commit: `d457b7189b137185a8dec800819a312c30b969fa`
- Authority branch: `product/issue-598-decisions-20260926`
- Implementation branch: `fix/issue-598-alignment-direction-reset-20260926`
- Approval: the owner explicitly accepted the complete direction and Reset choices and authorized Product Decision updates.
- Normative owner: [Option Integrity Product Contract](../OPTION_INTEGRITY_PRODUCT_CONTRACT.md).
  The linked concern-specific Packs preserve the exact approval receipts; this index and those Packs
  are not a second decision store or evaluator.

## Approved behavior

Display-direction selection is exclusive: Keep, all selected known-strand alignment features right-facing,
all left-facing, or Custom per-record Keep/right/left. The exact reference participates in direction choices;
its feature center retains its pre-Align canvas x. Source strands are unchanged. Each display reversal acts
on an entire record, including its features, labels and comparison ribbons.

Reset offers positions only or positions plus directions actually changed by the latest Align. Combined
Reset restores absolute before-Align directions only for those records, including a changed reference.
A target's later manual direction edit is overwritten only by the explicit combined scope and is disclosed
before Apply. Other settings and pending edits survive. Save/Load preserves trustworthy restoration receipts;
missing old direction evidence disables only direction restoration with an explicit reason.

BUG-17 (a new linear crop/rotation alternative) is out of scope. Existing circular rotation, linear region
editing, deterministic anchor resolution and narrow-screen palette coverage retain their accepted contracts.

## 関心別の承認済みDecision Packs

各Packは一つのconcernとその完全なowner receiptを保存する。責任・検証・実装routeを分け、総合計画書は依存関係だけを統合する。

| Decision | 独立したPack | 責任 |
| --- | --- | --- |
| PD-OI-027 | [表示方向と幾何変換](decisions/01_DISPLAY_DIRECTION_AND_GEOMETRY.md) | transform semantics |
| PD-OI-029 | [ResetとHistory](decisions/02_RESET_AND_HISTORY.md) | restoration/history |
| PD-OI-031 | [操作面と自動Align](decisions/03_OPERATION_SURFACE.md) | discoverability/default flow |
| PD-OI-034 | [draft・検証・再試行](decisions/04_LOCAL_REVIEW_AND_RETRY.md) | local intent/validation/recovery |



## Authority and implementation separation

The authority branch changes only the active Product Contract. This documentary approval record is stored
on the implementation branch with the plan. The authority change does not change
runtime, checker/workflow code, fixtures, mapped evidence or output baselines. Implementation begins only
when the four exact accepted records exist on the latest `origin/dev` base. Candidate records on a work
branch cannot authorize that same candidate runtime. No `BD-###` record or new authority store is created.
