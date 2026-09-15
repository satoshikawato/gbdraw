# INSTRUCTION PROMPT — S05: 共通解析境界と容量制限付き中間 cache

両modeが共用する解析段階を整理し、Webで設定変更時に中間結果を再利用できるようにしてください。64件共有LRUの連続missを解消し、正しい入力検証とWorker/Result lifecycleを保ってください。

## 必読と前提

- [総合計画書](MASTER_PLAN.md)、特に §3〜§5.4、§7〜§9。
- `results/S01.md`、`results/S03.md`、`results/S04.md` とその実装。
- `gbdraw/web/CLAUDE.md`、最新OIPC・Product/architecture ratchet・raw/derived/session cache contracts。
- S02のconsumer inventoryがあれば読む。ただし未選択のgraph表現を実装しない。

## 担当範囲

- Python解析層の既存関数を、準備・推論・projectionとして必要な範囲で分解する。
- `gbdraw/web/js/app/python-helpers.js`、必要な `gbdraw/web_support/` 私有owner。
- `run-analysis.js`、diagram worker/protocol/servicesのtransportとlifecycle接続。
- 関連Python/Node/Browser tests、benchmark、`results/S05.md`。
- 検索schedulerの新設、別runtime、path/schema移行、optional inference UIの並行実装は対象外。

## 実装手順

1. filtered、member selection、normalized、inferred、block/displayの実際のread-setを調査し、各keyの必要項目を記録する。処理順とper-table score fitting等の範囲を変えない。
2. CLI/Python/Webと両modeが同じ解析関数へ収束する最小の境界を作る。Web helperの文字列内に第二の推論実装を置かない。
3. appearance/block設定変更で不要な再正規化・再推論を避ける。向き・record順・annotationはID/順位/命名を読むconsumerまで追い、証明していないkey項目を除去しない。
4. 64件共有LRUを、直近一解析の中間データを単一ownerがバイト予算内で管理する方式へ置換する。既存の適切なbounded cache ownerがあれば拡張し、並行ownerを作らない。
5. native/browser測定から予算と算定方法を決める。DataFrame/object/graphのpayloadと保持中の旧snapshot、転送copyの区別を記録する。保持量が予算を超える場合は保持せず、理由をdiagnosticsへ出す。
6. cache miss/over-budget時も同じ関数で結果を計算する。新cache hitでidentityや必要なvalidationを飛ばさない。
7. 不変な中間値を共有し、mutateする境界のみcopyする。Python converted JSONとJS derived cacheのconsumerを照合し、重複保持を整理する。cached drawings/catalog/SVG/final Resultを追加しない。
8. cache publishの確定点とerror/cancel/staleの破棄を既存Worker契約に合わせる。完了raw retryは既存ownerを呼び、pending derived結果と同一寿命へ統合しない。
9. Clear Cache、source置換、Session/History置換、Worker終了で適切に失効させる。遅れて完了したoperationが新ownerへcacheを入れない。
10. instrumentationは既存diagnosticsを拡張し、sequence/raw行/巨大manifestをlogへ出さない。

## 必要な検証

- 49/64/81方向tableを用い、一式が予算内ならwarmの色/block変更でparse missが0となる。
- 同じcaseで小さい予算を与え、保持上限・over-budget理由・出力一致を確認する。任意のtable数で必ずwarm hitすると主張しない。
- member/filter変更、向き、record順、metadata、source contents、runtime bindings、searchContext変更ごとの正しい段階失効。
- cache有/無、cold/warm、success/errorの同じ出力・例外・provenance。
- raw完了→downstream cancel→member limit変更→Generateでraw searchを再実行しない。
- raw limit変更、Clear Cache、Session/History置換、stale完了、Worker再作成時の挙動。
- frozen inputへの破壊的変更、旧snapshotと新snapshotの二重保持によるpeakを確認する。
- Node/helper testsと実browser/Pyodide経路を検証し、転送量・処理時間・メモリを別々に記録する。
- 共有Worker変更の影響としてCircularの既存生成とlazy initializationも確認する。

## 合格条件

- 段階ごとの依存表と実装keyが一致し、二つ以上の実consumerが同じ解析境界を使用する。
- 件数境界でのLRU thrashが消え、保持量が宣言した予算以内である。
- peak memoryの測定とcache予算を混同せず、oversize時の制約を報告する。
- existing raw retry、last successful Result、Session/History、cancellation/stale guardsを維持する。
- 旧共有LRU・不要なcompleted payload保持・移設元の推論経路を同じ変更で削除する。
- production/tests/generated artifactsを別々に監査し、`results/S05.md` に実在API、key/read-set、予算、確定点、検証、英語commit title/summaryを記録する。
