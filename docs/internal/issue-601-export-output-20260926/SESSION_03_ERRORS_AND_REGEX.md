# INSTRUCTION PROMPT — S03 原因保持・安全なDetails・regex構文案内

あなたはgbdraw Issue #601のerror処理とregex案内の実装担当です。失敗した操作に既知の原因と利用可能な次の操作を示し、任意Detailsへ許可済みの診断情報を表示してください。Color/LabelのPython regexとFeature SearchのJavaScript regexの評価経路を維持し、入力欄でそれぞれのsyntaxを案内します。

## branch取得と必読資料

SESSION_IDは`s03`です。[総合計画書](./MASTER_PLAN.md)の共通取得手順で**`fix/issue-601-export-output-20260926`**をremoteから取得し、専用cloneで最新push済みcommitから作業してください。AGENTS.md、CLAUDE.md、Web CLAUDE、総合計画、[Decision Pack 02](./DECISION_PACK_02_ERROR_DISCLOSURE.md)、S02結果、S01のerror契約とmerge状態を読みます。

通常はS02 push済みから開始し、shared fileの所有を引き継ぎます。PDF実現性が保留でもerror concernは独立して進められます。その場合は先行branch状態と対象fileの編集が停止していることを確認し、S01のerror契約がorigin/devに統合された後、実装branchへ通常mergeして作業します。未merge候補はruntime authorityではありません。

## 所有範囲

省略したJS pathは`gbdraw/web/js/`からの相対pathです。

- 比較identity例外のPython producer、既存Python render/helper adapter。
- `workers/diagram-generation-worker.js::serializeError`と`services/diagram-generation.js::deserializeWorkerError`。
- `services/error-normalization.js`: 利用者文言、safe params、Details、recoveryの一owner。
- `app/run-analysis.js`、`app/similarity-alignment.js`、rule/label import、`app-setup.js`: error result/causeと次の操作。
- `index.html`、Feature SearchとStandalone Interactive SVG検索: 実操作に合う表示とsyntax案内。
- 関連testsと`SESSION_03_RESULT.md`。

PDFのfont-selection ownerやassetはS02成果を使用し、別error storeや第二のregex evaluatorを作りません。authority、guard/CI、Session/request schemaをruntimeと同時に変更しません。

## 実施内容

1. 対象producer/callerの既知failureを一覧にし、現在利用者へ示している修正情報を記録します。比較identity不一致、Python regex syntax、Worker/font/library準備、unknown、cleanup、Circular track-slot等、移行する各failureのtyped code、actual stage、safe params、利用可能actionを定義します。
2. 比較identity不一致を識別可能な例外にし、既存ValueError捕捉を維持します。Pythonの実exception type/causeと`re.error`位置で分類し、raw messageの文字列heuristicを増やしません。自動補正やvalidation無視で成功扱いにしません。
3. Pyodideが例外をstringへ変える前に、既存render/helperの共通adapterでcode/stage/safe paramsを保存します。native CLI/APIの例外動作を変えません。serializer/deserializer双方を同時更新し、code/stageがWorker往復で失われないようにします。新しい永続protocolや旧schema readerは作りません。
4. `run-analysis.js`のerror returnを全監査し、causeが必要なcallerへ明示的なerror resultを返します。global errorLogの読み戻しでcauseを推測しません。Align reviewに原因を表示してdraft、record方向、Result、canonical request、Historyを維持し、retryできるようにします。
5. 一つのnormalizerでknown codeから利用者文言とactionを作ります。unknownはactual stage、原因不明、安全な次の操作を示します。既知のvalidation修正情報を全てunknownに潰さず、typed paramsへ移行します。二重正規化でcodeを失わず、prefixを重複させません。
6. Detailsを初期collapseの任意表示にし、code/stage/row/field/position/codepoint等のallowlist・型・長さを検査します。summary 1,000字、detail text 4,000字、section 8個の上限を維持します。raw exception/stack/traceback/stdout/stderr、入力sequence/全文/SVGを画面や通常consoleへ公開しません。cleanup二次失敗で主因を上書きしません。
7. 見出しとrole/nameを実操作へ合わせ、Exportを一律にGeneration Errorと呼ばないようにします。PDF recovery actionはS02のsame snapshotと利用可能性を使い、成功ResultやHistoryを変更しません。cancel/stale/supersededにerror UIを出しません。
8. Color/Label入力にPython regex、case-insensitive、`(?i)NADH`/`(?P<enzyme>NADH)`を案内します。Feature SearchとStandalone SVG検索はJavaScript regex、case-insensitiveの案内と実際のsyntax errorを一致させます。準備失敗をsyntax違反と断定せず、既存Python owner、atomic commit、History、staleを維持します。
9. 置き換えたraw summary/個別prefix、移行済みfailureのCircular特別extractorなどを同じ変更で削除します。表示先ごとのclassifierや二重ownerを残しません。

## 必須検証

source→adapter→Worker/client→caller→UIでcode/stage/paramsとactionを確認します。Pyodideでの型消失、末尾causeのあるtraceback、二重正規化、unknown、cleanup失敗、prepareとregexの区別を実際の経路で検査してください。

private sentinelを原例外・cause・stdout/stderr・uploaded contentに入れ、summary/Details/通常consoleへ出ないことをassertします。既知修正情報、上限、allowlist外fieldの非表示、cancel/stale非error、draft/Result/方向/History/canonical requestの保持、retryを自動検査します。

Python-only Color/Label、TSV/Session/preset、invalid rule atomic rejection、Feature SearchとStandalone SVGのJS拒否を保護します。BUG-19の画面・pattern・buildが未特定なら、その未確認を記録し、評価器を再実装して再現確認の代用にしません。

desktop/390 px、keyboard/focus、Details展開、利用可能actionをbrowserで検査します。新しいblocking受入を既存PR inventoryへ入れ、trusted-base checker、architecture contracts、関連Node/Python/Playwrightと既存PDFの回帰検査を実行します。

## 終了時のcommit・push

`SESSION_03_RESULT.md`にknown failure移行表、owner/path、Python例外互換、raw非公開とcause保持のcommands/results、UX、rollback/retry、残る再現調査を保存します。共通終了手順で**commitし、同名remote branchへpushしてください**。
Commit title例: `Preserve actionable error causes and clarify regex syntax`。
remote SHA一致を確認してS04へhandoffします。
