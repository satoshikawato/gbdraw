# Product Decision Packs

Issue #619の表示不具合を修正し、Circular Width/Radiusを数値＋単位controlsにするための独立した製品仕様を扱う。

- [01 — 入力表現](01_INPUT_REPRESENTATION.md): 通常の単位選択と旧percent/suffix入力。
- [02 — 単位変更](02_UNIT_CHANGE.md): 数値維持か、geometry維持か。
- [03 — Auto単位の寿命](03_AUTO_UNIT_LIFECYCLE.md): transientか、History/Sessionへ保存するか。

推奨は各A。各responseはrationale、Must preserve、May retire、riskまで記入済みで、Owner署名だけが未記入。日付が異なる場合は更新する。方向としての数値＋単位分離を、三つの未署名の詳細仕様の承認へ拡張してはならない。

このディレクトリはレビュー用の候補文書であり、新しいauthority registryではない。署名後はS00が既存のstatic Product Contractへ正確にserializeし、authority-only統合を確認する。別Packは別concernとして決裁する。
