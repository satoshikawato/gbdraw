# INSTRUCTION PROMPT — S03: HSP 集計を一回走査へ変更

Collinear / Similarity groups の共通HSP集計を、結果を維持して改善してください。実装、差分比較、性能・メモリ測定、handoffまでをこのセッションで完了してください。

## 必読と前提

- [総合計画書](MASTER_PLAN.md)、特に §3、§5.2、§7〜§9。
- `results/S01.md` と採用されたrunner・fixture・基準値。
- 最新の `protein_colinearity.py`、対象tests、repositoryガイド。

S01が必要です。S02のProduct選択を待つ必要はありません。

## 担当範囲

- `gbdraw/analysis/protein_colinearity.py` のHSP集計と直接必要なprivate helper。
- `tests/test_protein_colinearity.py` とS01の関連benchmark case。
- `results/S03.md`。
- 所属algorithm、cache、検索設定、path/schema変更は担当外です。

## 実装手順

1. callerからnumeric coercion、member limit、HSP集計、coverage filter、score modelまでの順序を確認する。
2. `hits.groupby(..., sort=False)` ごとの小DataFrame / `itertuples()`生成を、table全体の一回のrow iterationとペア別accumulatorへ置き換える。
3. 既存の代表HSP順位、coverage interval/clamp/unionのhelperを再利用する。代表行、区間、count、length合計を保持し、不要な全rowの二重materializationを避ける。
4. ペア順、ペア内の元行順によるtie-break、出力column/dtype、空tableを維持する。NaNを文字列IDへ変換して別の有効pairを作らない。
5. missing/unknown ID、長さ0、不正座標、NaN/非有限値の現在のaccepted/rejected境界を維持する。既存authorityと矛盾する場合はS01の分類に従い、偶然の挙動を暗黙に新契約にしない。
6. 旧本番ループを削除する。比較用oracleはtests/benchmarkに限定する。

## 必要な検証

- overlap/disjoint/隣接区間、reverse座標、範囲外clamp、複数HSPのunion。
- 同点score/evalue/identity/length/coordinate、完全同点の元行順。
- 空入力、欠損・未知ID、短い/不正protein、numeric coercionと例外。
- finite/unbounded member limitが選んだpairについて全HSPを維持する。
- 集計DataFrameのcolumn/dtype/値/順序と、最終group / CollinearityResultがbaselineと一致する。
- Gallery両modeとsyntheticで、行iterationのtable単位化、native時間、peak memoryを確認する。
- `pytest tests/test_protein_colinearity.py tests/test_collinearity.py tests/test_collinearity_units.py -v` を実行する。追加のgateは変更範囲に応じる。

## 合格条件

- coverageや代表HSPを変えず、ペアごとのpandas処理を除去できた。
- 時間改善とmemory tradeoffを実測し、小規模caseの回帰も確認した。
- 必要なinterval sortを残しており、全処理が厳密にO(H)だとは主張していない。
- 本番に旧/新の切替flagやfallbackが残っていない。
- productionとtestsのdiffを別々に監査し、`results/S03.md` に採用API、measurements、残課題、英語commit title/summaryを記録した。
