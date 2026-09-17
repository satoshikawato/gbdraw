# INSTRUCTION PROMPT — S01: 基準・再現・契約確認

このセッションでは、Collinear / Similarity groups の性能改修の比較基準と再現資産を作成してください。後続の本番algorithm変更は実装せず、この範囲を測定・記録・引き継ぎまで完了してください。

実行結果: [S01 handoff](results/S01.md)、[authority / consumer inventory](results/S01_INVENTORY.md)。
採用した単一entryは `tools/benchmark_protein_comparison.py`（`run` / `compare` / `browser`）。
後続は handoff の固定設定と入力hashを用いる。

## 必読と前提

- [総合計画書](MASTER_PLAN.md) 全体。特に §2、§3、§7〜§9。
- repository `AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`。
- `docs/internal/report.md` が存在すれば読む。過去の計測値・一時script・未追跡状態を最新baseの証拠に置き換えない。
- 最新baseの `OPTION_INTEGRITY_PRODUCT_CONTRACT.md`、comparison semantics、session compatibility、architecture/Product ratchet。

## 担当範囲

- 既存benchmark runnerの調査と今回の最小のrunner / synthetic fixture生成 / 比較処理。
- 入力・baseline・authority・consumer・外部依存のinventory。
- この計画パックの該当箇所と `results/S01.md`。
- 本番コードはread-only。恒久的な計測hookが不可欠と判明した場合は境界と理由を先に記録し、最小限に限定する。

## 実行内容

1. working treeとlatest origin/devを確認する。計画時SHAと異なる変更を分類し、既に解決した問題を再実装対象から外す。
2. `47f5cebd` 相当のバッチidentity処理と `7db0539a` 相当のraw retryを確認する。key生成のfull manifest validationが何回か、direction/searchContext/順序が維持されるかを測る。
3. optional Collinear inference、limit retention、source batchingのauthorityと実装状態を照合する。optional inference/limit記憶は文書保存時の `65f231af` に統合済みである。欠けるcheckoutでは依存revisionを確認し、性能改修内に並行実装を作らない。
4. Galleryの保存rawとGenBankから report のcaseを再構築する。record/runtime ID一致を確認する。生物学的レコード数、source数、方向table数、実source job数を区別する。
5. `MASTER_PLAN.md` §8の合成caseを再生成可能にする。大規模な旧全経路列挙は実行せず、小規模の実測と式による値を分ける。
6. 既存 `tools/benchmark_diagram_layout.py` 等の方式を確認し、適切な既存runnerを拡張するか、最小の `tools/benchmark_protein_comparison.py` を作る。採用した単一entryと引数を記録する。
7. 別source-root / 別processで同一fixtureを評価できるようにする。JSON結果にsource revision、設定、入力hash、各sample、operation count、結果比較を残す。
8. wall-clockとprofile/memoryを分離する。計測回数、warmup、noise評価、regression判定基準を変更前に決める。native Pythonと実browserの転送・Worker/Pyodide時間を区別する。
9. `_build_ortholog_paths()` から Python API / typed resource / SVG / catalog / popup / Session までを検索し、S02用のconsumer inventoryを開始する。
10. focused baseline testsを実行する。既存failureは再現command・原因boundaryを記録し、skipやglobal stubで合格に見せない。

## 検証

- runnerを実際に実行し、同じseed/hash/settingsから同じsemantic結果になる。
- source-root取り違えや別Python importの混入を検出できる。
- HSP、sparse、cache、path、mergeのcounterと測定対象が明示される。
- baselineの「現在の観測」と「維持すべきauthority」が区別される。
- batch1回validationとキー一致は既存改修として記録される。
- raw search、post-search、metadata、serialization、render、転送の時間を混同しない。

## 完了条件と引き継ぎ

- reportに依存した一時scriptを再発見できなくても、後続が実行できる再現資産がリポジトリにある。
- exact command、fixture生成、出力先、比較方法、limitsを `results/S01.md` に保存する。
- S02にconsumer inventory、S03にHSP oracle、S04にsparse oracle、S05にcache/lifecycle baselineを渡す。
- 未解決のauthority conflictや外部依存を明記する。本番改善を実装済みと報告しない。
- 総合計画書 §9 の形式で、英語commit title/summaryを含めて報告する。
