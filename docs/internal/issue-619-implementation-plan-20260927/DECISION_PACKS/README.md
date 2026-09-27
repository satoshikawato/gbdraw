# Product Decision Packs

Issue #619の表示不具合を修正し、Circular Width/Radiusを数値＋単位controlsにするための独立した製品仕様を扱う。

- [01 — 入力表現](01_INPUT_REPRESENTATION.md): 通常の単位選択と旧percent/suffix入力。
- [02 — 単位変更](02_UNIT_CHANGE.md): 数値維持か、geometry維持か。
- [03 — Auto単位の寿命](03_AUTO_UNIT_LIFECYCLE.md): transientか、History/Sessionへ保存するか。

推奨は各A。各responseはrationale、Must preserve、May retire、riskまで記入済みで、Owner署名だけが未記入。日付が異なる場合は更新する。方向としての数値＋単位分離を、三つの未署名の詳細仕様の承認へ拡張してはならない。

このディレクトリはレビュー用の候補文書であり、新しいauthority registryではない。署名後はS00が既存のstatic Product Contractへ正確にserializeし、authority-only統合を確認する。別Packは別concernとして決裁する。

## S00 evidence handoff (2026-09-27)

3 Pack は引き続き未署名。最新 authority search の基点は
`origin/dev` = `d313b70b9f97c2c1d70f9ae885edbead80b62021`。
対象の署名済み receipt / authority はなく、推奨 A は候補のまま。
[S00結果](../SESSION_RESULTS/S00.md) に既存 authority と実測境界を記録する。

既存 typed Web draft の `{value:"1.",unit:"px"}` /
`{value:"1e-3",unit:"factor"}` は normalization、validation、request projection、
実 Save/Load、native reader、History で受理・保持された。
不正 numeric text と不正 unit は draft/History で保持し、Web Save/Load と
request projection で拒否する。typed value の空文字は Auto ではなく不正。
Auto は既存の literal null / 空欄を使う。

native reader による config 保持は Web admission の証明ではない。
内部の nonfinite number は History で null に変わり得るため、不正編集を
その値へ変換してはならない。Boolean の既存 coercion も有効 domain の
承認ではない。比較結果と再現入力は S00 の evidence にあり、
新 schema や invalid draft の null 化を必要とする失敗解消は行っていない。

Pack02-B の geometry 換算、Pack03-B の永続 namespace は未検証。
本 evidence はどの Product option も選択せず、署名の代わりにならない。

settings-only の実 Save は typed draft と無編集 fresh-app control の双方で
既存 writer/validator 不整合により拒否された。通常のrendered Sessionでの
表現受理を全Session形態の受理へ広げない。runtimeは本S00で変更せず、
この境界の解決と実Save/Load/native証拠もS01開始前の前提として残す。
