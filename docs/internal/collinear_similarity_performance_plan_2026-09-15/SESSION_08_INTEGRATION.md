# INSTRUCTION PROMPT — S08: 統合・性能・回帰検証と最終引き継ぎ

この計画の実装を統合された実経路で検証し、達成範囲と未達を明確にした最終handoffを作成してください。新機能や追加の最適化を際限なく広げず、検出したin-scopeの不具合を修正して必要なgateを完了してください。

## 必読と前提

- [総合計画書](MASTER_PLAN.md) 全体と `results/S01.md`〜`results/S07.md` の実在するhandoff。
- 最新baseのauthority、各実装revision、test inventory、generated browser wheel準備方法。
- repository `AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`、architecture/Product ratchet。

S03〜S05とS07の扱いが確定していることが前提です。S06が未決なら限定的検証は進められますが、計画全体の完了を宣言してはいけません。optional inference / limit記憶は `65f231af` の既存実装を回帰対象とし、必要な依存revisionを確認してください。

## 担当範囲

- 統合された解析、cache、path表現、block生成と各surfaceの回帰検証。
- 検出したin-scope不具合の最小修正とfocused regression tests。
- benchmark最終比較、既存契約文書の必要な更新、`results/S08.md`。
- deploy/publish/tag/push/merge、public showcaseの無関係な更新はこのpromptでは認可しません。

## 実行内容

1. 各handoffと実際のsourceを照合し、依存実装・Product authority・互換readerが同じ統合revisionで成立することを確認する。
2. S01で固定した入力・明示設定・runnerを使って変更前後を比較する。fresh defaultの違いや古いreportの秒数を改善率へ混ぜない。
3. HSP、sparse/dense、49/64/81 cache、path、mergeのoperation count・時間・転送量・memoryを記録する。改善率は対象stageとenvironmentを明示する。
4. CLI/Python/Webがそれぞれ公開する範囲で、Similarity、Collinear adjacent/all、inference ON/OFF、有限/無制限limit、multi-record sourcesを実経路で確認する。source job数とdirectional table数を別々に照合する。
5. warm Generate、色/block/member/filter変更、向き・順序・source変更を検証する。cache keys、実行段階、provenance、最終結果が整合することを確認する。
6. raw完了後cancel→member変更→retry、raw設定変更、Clear Cache、Session/History置換、stale完了、Worker再作成を実browserで実行する。last successful Resultが残ることを確認する。
7. supported old sessionのload/current save/fresh load/regenerate、typed resourceとmetadataを確認する。PATH-Bでは通常の全経路materializationがないことを実測する。
8. 共有Worker/cacheのCircular smoke、local assetsだけの実行、lazy initialization、メモリ保持とreleaseを確認する。
9. production / tests / docs / generated artifactsを別々にレビューする。旧全走査、並行scoring、共有LRU、不要payload保持、無根拠なcompatibility path、未使用flagが残っていないかを確認する。
10. 不具合を修正した場合はfocused testと影響するgateを再実行する。既に通った無関係な検証を理由なく繰り返さない。

## 必須 gate

```bash
pytest tests/test_protein_colinearity.py tests/test_collinearity.py tests/test_collinearity_units.py -v
pytest tests/test_web_feature_catalog.py tests/test_session_request_codec.py tests/test_session_compat.py -v
pytest tests/ -v -m "not slow"
pytest tests/test_output_comparison.py::TestOutputComparison -v
ruff check gbdraw/
```

最新のNode inventoryからcache、derived/raw identity、session、worker、architecture/Product testsを実行してください。代表対象とPlaywrightの代替手順は総合計画書 §8.3にあります。browser wheelが必要なら生成し、trackedな配布資産として扱わないでください。

tests/reference_outputsの再生成で意図しない差を隠さないでください。意図した契約変更のmetadata差は、geometryや科学的選択の差から分けて説明してください。

## 最終報告と合格条件

`results/S08.md` に以下をまとめてください。

- 比較baselineと最終revision、再現command、環境、fixture/hash。
- 各report指摘の状態: 既存修正の検証 / 今回修正 / 条件付き見送り / 判断待ち / 未達。
- native stage時間、browser end-to-end、転送量、retained/peak memory、operation countとnoise。
- 科学的結果・ID/順序・保存互換・lifecycleの検証結果と実行していないもの。
- owner/pathとcompatibilityの前後、削除した旧経路、Product選択の根拠。
- PATH-Aなら指数出力量の残存。PATH-Bでも明示的な全量取得の不可避なコスト。
- 外部依存が未統合なら残る具体的なgate。未実装OFFをmockだけでpassにしない。
- rollback、残課題、英語のproposed commit titleと短い英語summary。

全体完了は、選択済みoutcomeの実装と必要なgateが全て成立したときだけ宣言してください。計画が終わったことと、指数増加問題を解消できたことを混同しないでください。
