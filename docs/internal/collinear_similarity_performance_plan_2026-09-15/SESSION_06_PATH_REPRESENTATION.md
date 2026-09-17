> 2026-09-16 task override: S05 is rejected; execute S06 from S04. PATH-B authority
> is integrated in origin/dev `5e0cb0fa` through PR #536. Reuse matching S04 final
> current-side measurements; do not rerun unchanged S04 benchmarks. This override
> replaces the rejected S05 dependency. Implement S06 only; do not publish it.

# INSTRUCTION PROMPT — S06: 選択済みの経路表現を実装

S02で具体化され、適用可能なbase authorityとして受理された経路契約を実装してください。推奨案と承認済みoutcomeを区別し、未選択の公開・保存挙動を実装しないでください。

## 必読・前提・開始条件

- [総合計画書](MASTER_PLAN.md)、特に §3、§5.5、§7〜§9。
- `results/S02.md`、存在する場合は `results/PATH_DECISION_PACK.md` と完全なProduct response / authority。
- `results/S04.md` と実在する解析・consumer境界。S05は却下済みで依存に含めない。
- 最新API、typed resource/session/catalog形式、compatibility履歴、Product/architecture ratchet。

S04を基準とし、選択されたoutcomeを適用可能なbase authorityが認可していることが開始条件です。authority-only候補とdependent runtimeを同じ変更で自己認可しないでください。既存authorityだけで完全に決まるとS02が証明した場合は、その根拠で進めてください。

条件不足なら、不足するdecision/統合を具体的に記録し、実装を推測して進めないでください。独立したconsumer確認や既存fixtureの検証は進められます。

## 担当範囲

- Python共通解析のpath/graph表現と必要なprivate関数。
- 必要なpublic API境界、`session_request_codec.py`、resource decode/encode、metadata/catalog/実consumer。
- 当該namespaceの認可済みreader/writer、focused tests、benchmark、正確な契約文書。
- `results/S06.md`。
- 新しい解析mode、推論規則、近似、未要求のpaging UIは追加しません。

## PATH-B / lossless-graph が選択された場合

1. S02で証明したnode/edge/重複/terminal規則から、lossless graphを共通の標準表現として実装する。
2. DAG性を構築boundaryで担保する。到達可能な循環がS02で未解決なら、edge削除や向き変更で迂回しない。
3. count/shared情報と、承認された既存path ID対応を全列挙なしで導出する。大整数の表現も選択済み契約に従う。
4. 旧APIのtuple等を維持すべきconsumerには明示的adapterを置く。二つの独立した推論algorithmを残さない。lazy型を黙って既存tupleへ代入しない。
5. ordinary Generate、typed resource転送、SVG metadata、catalog、Session save/load/replay、derived cacheが全経路をlist/JSONへ展開しないようにconsumerを接続する。
6. 全量取得を維持する契約なら、その明示的境界でだけ展開し、時間・出力量の制約を正確に説明する。上限で黙って切らない。
7. release/main履歴で必要と証明された旧readerを、既存compatibility ownerでcurrent modelへ一度変換する。branch-only artifactはcurrent writerで更新し、中間migration chainを削除する。
8. request/resource/session/cache/catalogのどのnamespaceを変更するか明示し、無関係なschemaをまとめて上げない。新しいCBが生じる場合はratchetの例外条件を満たす。
9. tests/reference_outputsを通常検証で更新しない。認可済みmetadata変更と意図しないgeometry差を分ける。public docs/figure更新が必要なら該当skillと生成ownerを利用する。

## PATH-A / exhaustive-current が選択された場合

現在の全経路配列・型・順序・ID・保存契約を保持してください。S02で同値と確認した不要中間copy等の改善だけを実装できます。通常の全量出力が指数コストを残すことを `results/S06.md` に明記し、グラフ化による通常処理の多項式化を達成したと報告しないでください。

## 必要な検証

- 小規模の全path、path順、edge IDs、最初のpath ID、sharedProteinIdsを旧oracleと一致させる。
- 複数source/sink、重複edge、coortholog、逆向きevidence、record-local、空group、single chain。
- PATH-BではR=8/12/16の既知count、R=24以上のcompact動作、大整数countを検証する。
- ordinary generation / serialize / saveで列挙関数が呼ばれないことをoperation counter等で確認する。公開全量取得のテストとは分ける。
- direct CLI/Python/Web、typed encode/decode、supported old session→current model→save→fresh load→regenerate。
- cache有/無、cancel/stale時、metadata popup、表示・block/path属性、source/record順変更。
- native/browserの時間、転送量、peak/retained memoryを測り、既存JSONコピーへの移動だけで終わっていないことを確認する。

## 合格条件と引き継ぎ

- 実装したoutcomeが受理済みauthorityと一致する。
- PATH-Bでは通常の解析・描画・保存でP/Lに比例した全量materializationがない。
- 全量consumerの不可避なコスト、互換経路、未達を明記する。
- 必要なfocused Python/Node/Browser testsが通り、owner/path/CBとrollbackを記録する。
- `results/S06.md` に実装契約、format/reader、計測、実際の限界、英語commit title/summaryを記録する。
