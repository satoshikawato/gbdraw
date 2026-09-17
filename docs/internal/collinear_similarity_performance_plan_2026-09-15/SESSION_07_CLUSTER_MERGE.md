# INSTRUCTION PROMPT — S07: Collinear cluster merge の評価と必要な改修

2026-09-16実行完了。性能評価はユーザーの明示指示により合格。実装、同値検証、時間・memory・counter、留保と引き継ぎは[results/S07.md](results/S07.md)を参照。S08は未実装。

S04→S06の成果を引き継ぎ、Collinearのcluster mergeを評価し、必要な場合に限って同値な改修を実装してください。実装を見送る場合も、判断を再現可能な測定とともにhandoffしてください。

## 必読と前提

- [総合計画書](MASTER_PLAN.md)、特に §3、§5.6、§7〜§9。
- `results/S01.md`、`results/S03.md`、`results/S04.md`、`results/S06.md`。
- S05は2026-09-15に却下済み。共通解析class、cache owner、Worker transaction、runtime変更を依存に含めない。
- S06の通常経路から全path列挙を除いた成果は了承済み。S06の測定留保は維持するが、Galleryの一律高速化や独立性能合格をS07開始条件にしない。
- S06未コミット成果をレビュー台帳と照合して一つのローカルcommitに固定し、最新origin/dev由来の専用worktreeへ未統合のS01〜S04/S06だけを適用する。
- cluster merge / max_conflictsの最新authorityと `gbdraw/analysis/collinearity.py`、対象tests。

## 担当範囲

- `_lossless_conflicts_between_clusters()` とmerge caller、直接必要なprivate endpoint/index処理。
- `tests/test_collinearity.py`、`tests/test_collinearity_units.py`、S01のmerge benchmark。
- `results/S07.md`。
- group inference、score、path/schema、source検索、cache方式は変更しません。

## 実行手順

1. 既存のGalleryとmerge-300/600/1200・merge-edgesの時間、profile、counter、memory、oracleを先に確認する。測定関数と依存・入力・設定・境界が一致する結果は再利用し、新しいisolated merge境界や不足counterだけを固定S07 baselineで追加測定する。無関係なファイルhash差だけで既存結果を無効にしない。
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

時間はwarmup 1回・7 samplesを基本とし、21 samplesのbaselineには21 samplesで対応する。中央値10%超の悪化はregression、MAD/median 5%超はinconclusive。timing、probe/counter、tracemalloc、browser memoryは別runとし、checkout/build/tests/別benchmarkを時間測定へ重ねない。S06の留保を維持し、S03/S04/S06の改善率をS07の成果へ流用しない。S07はS06とは別の一つの変更として引き継ぎ、S08以降を実装しない。
