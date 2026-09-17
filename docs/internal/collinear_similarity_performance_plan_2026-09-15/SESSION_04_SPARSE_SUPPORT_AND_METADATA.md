# INSTRUCTION PROMPT — S04: 疎な所属候補探索と metadata 逆引き

共通group inferenceの総当たりをevidenceに基づく探索へ変更し、同じ情報を繰り返し線形検索するmetadata処理を整理してください。所属・診断・順序・最終出力を保った実装と検証を完了してください。

## 必読と前提

- [総合計画書](MASTER_PLAN.md)、特に §3、§5.3、§7〜§9。
- `results/S01.md`、`results/S03.md`、採用runnerとsparse/dense fixtures。
- group構築、record-local競合、metadata変換、既存testsとauthority。

S03の変更が作業baseに存在することを確認してください。S02のProduct選択は不要です。

## 担当範囲

- `gbdraw/analysis/protein_colinearity.py` のsupport候補とmetadata lookup。
- `gbdraw/web_support/orthogroup_metadata.py` のRBH/member逆引き。
- `gbdraw/analysis/collinearity.py` の対応metadata lookupと不要display要求の整理に必要な範囲。
- 関連Python tests、S01 benchmark、`results/S04.md`。
- cache、公開path形式、block merge algorithmは変更しません。

## 実装手順

1. `_build_core_support_candidate()` と `_best_evidence_between_protein_and_members()` が消費する全evidenceとtie-breakを確認する。
2. `best_by_direction` からincoming/outgoing索引を一度作る。protein→groupから接続先group候補を求める。
3. 固定 `core_member_snapshot` に対する索引を使って未所属の追加を判定する。追加済みmemberをこのsnapshotへ混入させない。
4. group内でも実際に接続するevidenceからsame/cross supportとdiagnosticを求める。候補groupを絞った後に巨大group全memberを二回走査する経路を残さない。
5. evidenceの採点と選択順位は既存ownerに集約する。新索引用に別のscoring式を複製しない。
6. record-local競合判定では、所属拡張後の集合から一度索引を構築する。既存のlocal component判定と追加の段階順を保つ。
7. protein→RBH group、member情報、endpoint→edgeの逆引きを必要なconsumerへ適用する。多対多や重複時にdict上書きでfirst-match/出力順を変えない。
8. 最新selectorの `comparison_pairs=()` をCollinearの不要なSimilarity表示projectionの抑止に使えるか確認する。all_edges、group、必要なmetadataを保てる場合だけ使用する。新skip flagを追加する前に既存契約を使う。
9. 旧全group/全member走査と重複lookupを撤去する。索引寿命は一回の解析段階内に限定する。

## 必要な検証

- incoming-onlyとoutgoing-only、same/cross-record、domain-only高score、membership supportの優先。
- bestとsecond-best、完全tie、ambiguous、low/high confidence、related edge数と順序。
- 固定coreと拡張後集合の差が結果に現れるfixture。
- record-local競合あり/なし、local group追加、representative、names/descriptions。
- metadataのRBH IDs、path IDs、member counts、edge選択・順序の完全比較。
- group/unassigned=200/200、400/400、800/800のsparse caseでcandidate/edge訪問回数を計測。
- 少数巨大groupとdense evidenceでも結果一致・memory使用量を確認する。
- Gallery両modeで最終groups、display tables、CollinearityResult、typed metadataを比較する。

## 合格条件

- evidenceのないgroupを追加してもcandidate評価がU×Gで増えない。
- 関連するgroup内の仕事も接続evidenceへ絞られ、全member総当たりが残存していない。
- support/diagnostic、snapshot、record-local、ID/順序/metadataを維持する。
- 索引追加によるmemoryと構築時間を測り、dense caseの限界を明記する。
- focused protein/collinearity/metadata testsが通り、architectureのowner/pathと削除旧経路を記録する。
- `results/S04.md` にsource revision、実測、差分監査、後続cacheが扱うデータの可変性、英語commit title/summaryを記録する。
