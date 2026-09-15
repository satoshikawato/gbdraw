# INSTRUCTION PROMPT — S07: Collinear cluster merge の評価と必要な改修

共通解析・cache改善後にCollinearのcluster mergeを再測定し、必要な場合に限って同値な改修を実装してください。実装を見送る場合も、判断を再現可能な測定とともにhandoffしてください。

## 必読と前提

- [総合計画書](MASTER_PLAN.md)、特に §3、§5.6、§7〜§9。
- `results/S01.md`、`results/S04.md`、`results/S05.md`。
- S06完了時は `results/S06.md`。S06判断待ちでもこのセッションは実行可能。
- cluster merge / max_conflictsの最新authorityと `gbdraw/analysis/collinearity.py`、対象tests。

## 担当範囲

- `_lossless_conflicts_between_clusters()` とmerge caller、直接必要なprivate endpoint/index処理。
- `tests/test_collinearity.py`、`tests/test_collinearity_units.py`、S01のmerge benchmark。
- `results/S07.md`。
- group inference、score、path/schema、source検索、cache方式は変更しません。

## 実行手順

1. Galleryと300/600/1200 anchorsのchain mergeで、conflict訪問数、sort回数、cluster copy量、時間とmemoryを測定する。
2. S01で固定した測定基準に照らし、現実の入力と合成stressの双方から改修の必要性を判断する。速度の小さな揺れだけを改善根拠にしない。
3. 必要ならquery orderを一度整列し、二分探索でstrictなquery interval内へ候補を絞る。subject intervalとcluster membershipの除外条件を維持する。
4. boolean merge判定だけでよいconsumerではmax_conflictsを超えた時点で終了する。exact countを要求するconsumerには正確な数を返す。
5. 順序済みendpointを再利用する。成長clusterのcopyを減らす場合、確定時のanchor順・derived属性・IDが一致するデータ構造にする。
6. initial clustering、merge候補順、singleton保存、最終sort/renumberの規則を変えない。
7. 高度なrange tree等は単純な索引後も測定上必要な場合だけ検討する。別algorithmのfallbackや小入力専用経路を増やさない。

## 必要な検証

- conflict=0、ちょうどmax、max超過、境界上のanchor、重複order、reverse orientation。
- 左右clusterに含まれるanchorの除外、間にあるsingletonの保存。
- chain merge、merge不可、複数pair、min_anchors、block IDと全anchor順。
- 既存resultとDataFrame/typed metadataの一致。
- Galleryの非profile時間、stressのscaling、sort/copy/visit counter、memory。
- focused collinearity / unit / protein tests。geometry差が疑われる場合はread-only SVG比較。

## 完了条件

- 改修した場合は、狭い区間へのcandidate絞り込みと同値性を示す。
- 区間queryがO(log N+k)になっても、merge全体が必ず線形とは報告しない。
- 見送った場合は現在の寄与、stressの限界、再着手条件を記録する。計測せず「優先度が低い」で終了しない。
- production/testsのdiffを別々に監査する。
- `results/S07.md` に実施/見送り、source revision、benchmarks、残る計算量、S08の検証対象、英語commit title/summaryを記録する。
