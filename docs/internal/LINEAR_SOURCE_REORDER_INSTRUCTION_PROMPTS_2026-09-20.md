# Linear入力File順序変更の復旧 — INSTRUCTION PROMPTS

このファイルは、gbdraw Web版Linearモードで`File 1`、`File 2`、`File 3`の順序を変更できない回帰を修正する実装セッションの開始指示である。前提となるチャット履歴はない。

問題、既存Product契約、目標アーキテクチャ、受入条件の管理先は[総合計画書](LINEAR_SOURCE_REORDER_MASTER_PLAN_2026-09-20.md)である。各プロンプトは、リポジトリルートをcurrent working directoryとする新しいセッションへコードブロック全体を渡して使用できる。

S0、S1、S2の順に使用する。同じ担当者が連続して実行してよい。各セッションの結果は総合計画書第16節の実施記録へ追記する。後続セッションは、総合計画書とリポジトリ内の証拠から状態を復元し、過去の会話を前提にしない。

## S0 — 最新baseのpreflightと実装境界の確定

```text
gbdraw Web版Linearモードの入力File順序変更回帰について、最新baseを監査し、
実装開始条件を確定してください。このセッションではruntimeを変更しません。

問題:
Input GenomesにはFile 1、File 2、File 3というsource単位カードが表示されるが、
Fileカードに順序変更操作がない。旧Up/Downは折りたたまれたrecord options内にあり、
linearSeqsの一recordだけを移動するため、multi-record Fileの移動には使えない。

最初に全て読むファイル:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/LINEAR_SOURCE_REORDER_MASTER_PLAN_2026-09-20.md
5. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
6. docs/internal/PRODUCT_IMPACT_RATCHET.md
7. docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.mdのPD-OI-018とOIC-015
8. docs/internal/WEB_CHANGE_POLICY.md
9. docs/REFERENCE/web-app.md

作業:
1. git status、現在branch、HEAD、origin/devとの差分を確認する。無関係な変更を保持する。
   branchやcommitを作る可能性がある場合はAGENTS.mdに従い、git fetch origin後の
   最新origin/devから追跡先なしの作業branchを作る。既に正しい作業branchがある場合は
   作り直さず、由来を確認する。
2. 総合計画書の監査基準eec715cと最新baseを比較する。少なくとも次を追跡する:
   - gbdraw/web/index.htmlのlinearSourceGroupsによるFileカードと旧Up/Downの位置
   - gbdraw/web/js/app/linear-sources.jsのsource grouping
   - gbdraw/web/js/app/app-setup.jsのapplyLinearSeqMutationとmoveLinearSeq系
   - gbdraw/web/js/services/session-request.jsのLinear records構築
   - tests/web/linear-multi-record.playwright.spec.jsのUI到達性coverage
3. Fileヘッダーにsource単位操作がないこと、旧actionがrecord一件を動かすこと、
   canonical requestがlinearSeqs順を使うことを最新baseで確認する。既に修正済みなら、
   重複実装せず総合計画書の受入条件で既存実装を評価する。
4. PD-OI-018 scenario revision 2とOIC-015が、one source = one File card、
   reordering後のrecord identity、placement、pair mapping、shared resources保持を
   authorityとしていることを確認する。新しいProduct Decisionや別authority storeを
   作らない。File moveで明示rowを再採番するなど、計画と異なる結果が必要なら、
   その差だけをPRODUCT_DECISION_REQUIREDとして止める。
5. 現在のテストを実行してbaselineを記録する:
   node --test tests/web/linear-sources.test.mjs
   npx playwright test tests/web/linear-multi-record.playwright.spec.js \
     --grep "one uploaded source stays one file card" --workers=1 --retries=0
   既存テストの成功をFile reorderの合格証拠としない。テスト本文が実際に何を操作し、
   何をassertしていないかを記録する。
6. NodeとPythonのPlaywright可用性、browser wheelの存在とsource対応を確認する。
   sandbox制約によるChromium失敗は、同じcheckを必要な権限で再実行して区別する。
7. 総合計画書第8節の予定ファイル、SR-01〜SR-14、owner/path設計が最新baseでも
   妥当か確認する。新規production module、Session field、watcher、別order stateが
   不要であることを確認する。

設計上の固定条件:
- linearSeqsを唯一の順序source of truthとする。
- source identityはgroupLinearSourceRecords()だけが決める。
- pure source-block transformはlinear-sources.jsが所有する。
- app-setup.jsは既存applyLinearSeqMutation()へ結果を渡すcoordinatorとする。
- record rowはUIDに付随したまま保持し、File moveで再採番しない。
- canonical request、Session schema、Worker、Pythonへ新経路を作らない。
- drag-and-dropや汎用sortable frameworkを導入しない。

このセッションの境界:
読み取り、baseline test、再現、authority/architecture確認までを完了する。
productionとtestを編集しない。総合計画書第16節へ実施記録を追記するが、受入条件を
結果に合わせて弱めない。push、PR、merge、tag、deployは行わない。

完了報告:
- base SHAとbranch
- 原因が残るか
- 既存testsの結果とcoverage gap
- Product Impact分類
- owner/pathのbefore/after予定
- S1を開始できるか
- 新しい判断、権限、外部操作が必要ならその正確な境界
を記載する。実装完了とは報告しない。
```

## S1 — source単位の順序owner、UI、focused testsの実装

```text
gbdraw Web版Linearモードで、Fileカードからsource全体を上または下へ移動できるように
実装し、focused testsと利用者向け文書を完了してください。

最初に全て読むファイル:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/LINEAR_SOURCE_REORDER_MASTER_PLAN_2026-09-20.mdの全体と第16節のS0実施記録
5. docs/internal/LINEAR_SOURCE_REORDER_INSTRUCTION_PROMPTS_2026-09-20.md
6. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
7. docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.mdのPD-OI-018とOIC-015
8. docs/REFERENCE/web-app.md

開始条件:
S0で記録したbase、作業branch、原因、テスト環境を確認する。最新sourceが変わった部分だけ
再監査する。無関係なdirty changesを保持する。Product結果が未解決なら依存箇所だけを
止め、独立して進められるテスト設計を終える。通常のローカル編集・testに再承認を求めない。

実装順序:
1. tests/web/linear-sources.test.mjsへ、source-block移動の失敗するdomain contractを追加する。
   single-record、2+3 multi-record、境界no-op、同名別File、非破壊性、source内順序、
   非連続な旧状態からの明示moveを検証する。実装の内部式をtestへ複製しない。
2. gbdraw/web/js/app/linear-sources.jsへ小さなpure transformを追加する。
   groupLinearSourceRecords()を再利用し、source groupを隣接交換してrecordsをflattenする。
   Vue、DOM、global state、File read、Session、cacheへ依存させない。引数を変更しない。
3. gbdraw/web/js/app/app-setup.jsへsource indexとdirectionを受けるactionと同じdomainに基づく
   availabilityを追加する。結果を既存applyLinearSeqMutation(next,
   { preserveLosatCacheInfo: true })へ一度だけ渡す。row、comparison、depth、cacheの
   並行reconcileを追加しない。pending discoveryのUID追跡を維持する。
4. gbdraw/web/index.htmlの各File NヘッダーへcompactなUp/Down操作を追加する。
   accessible nameはMove File N up/downのように対象と方向を含める。先頭/末尾をdisableし、
   keyboardで操作可能にする。File名、Remove、狭幅layoutを壊さない。
5. Record options内の旧moveLinearSeqUp/Down UIを削除する。production consumerがなくなる
   canMoveLinearSeqUp/Down、reorderLinearSeqs、moveLinearSeqUp/Downもtests専用に残さず削除する。
   recordの同一row内移動は既存moveLinearRecordWithinRow()だけを使う。
6. tests/web/linear-multi-record.playwright.spec.jsへ、実際のFile header buttonを操作する
   browser contractを追加する。window.__GBDRAW_APP__のmove action直接呼び出しだけで
   合格にしない。pointerとkeyboard、multi-record 2+3、DOM、linearSeqs、typed request、
   UID-row、selected endpoint、cache、collapsed disclosure、mobileを確認する。
7. 既存testsのmoveLinearSeq直接呼び出しを監査する。source移動を意図するtestは新しい
   source操作へ、同一row record移動を意図するtestはmoveLinearRecordWithinRow()へ移す。
   意味の違うtestを機械的に置換しない。
8. GFF3+FASTA pairとupload直後のasync discoveryをfocused browser testで確認する。
   File move後に全recordが正しいsource blockへ展開されることをassertする。
9. docs/REFERENCE/web-app.mdを更新し、File操作はsource全体、Record Layoutはrowと
   同一row内配置、明示rowはFile moveで再採番しないことを説明する。
10. focused checksを実行してin-scope failureを修正する:
    node --test tests/web/linear-sources.test.mjs
    node --test tests/web/session-request.test.mjs
    npx playwright test tests/web/linear-multi-record.playwright.spec.js \
      --grep "source.*order|File.*order" --workers=1 --retries=0
    git diff --check
11. production、tests、docsのdiffを別々に読む。dead action、重複owner、別state、
    不要なcompatibility branchがないことを確認する。

SOLID/KISS/DRY/YAGNI:
- source groupingとsource order transformを一ownerへ置く。
- UIはdomain groupingを再実装しない。
- reactive mutationは既存applyLinearSeqMutationを通す。
- fileOrder ref、watcher、Session field、drag framework、generic sorting abstractionを作らない。
- 旧record-level input reorder pathを新pathと並存させない。
- 将来必要かもしれない任意位置移動や複数選択を追加しない。

受入上の注意:
- File moveはrecord UIDに付随するrowを保持する。rowの上下はRecord Layoutが所有する。
- draft changeだけで既存Resultを変更しない。Generate成功後に新順序を反映する。
- 同名だが別uploadのsourceを結合しない。
- Load時の旧interleaved状態を自動migrationしない。明示move時だけsource blockへまとめる。
- source bytesとsemantically validなraw LOSAT cacheを保持し、numeric indexesを再解決する。

このセッションの境界:
production、focused tests、文書を整合した状態まで完成させる。失敗するtestだけを残して
終了しない。広いWeb suite、実Pyodide Generate、Save/fresh Loadの最終受入はS2でよい。
参照SVG、Gallery、public figureを必要なく更新しない。push、PR、merge、tag、deployは行わない。

引継ぎ:
- 総合計画書第16節へ以下を追記する
- base/branchと変更ファイル
- SR IDごとのfocused結果
- 正確なcommandとpass/fail
- 削除した旧owner/path
- 未実施のS2項目
- architecture concise evidenceの下書き
- 英語のproposed commit titleと短いsummary
を記録する。S2が残る間は全体完了とは報告しない。
```

## S2 — Generate、Session、accessibility、全体gateの最終受入

```text
gbdraw Web版Linear入力File順序変更の実装について、実ブラウザのGenerate、
Save/fresh Load、mobile/keyboard、広い回帰、architecture/policy gateを検証し、
見つかったin-scope不具合を修正して最終受入を完了してください。

最初に全て読むファイル:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/LINEAR_SOURCE_REORDER_MASTER_PLAN_2026-09-20.mdの全体と第16節
5. 第16節に記録されたS0/S1の対象branch、未commit差分、検証結果
6. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
7. docs/internal/WEB_CHANGE_POLICY.md
8. docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.mdのPD-OI-018とOIC-015

開始条件:
S1のsource、branch、production/test/docs diff、focused evidenceを確認する。別candidateの
証拠を無条件に流用しない。S1が未完了なら不足したローカル実装とfocused testsを先に直す。
無関係な変更を保持し、受入条件やtimeoutを結果に合わせて弱めない。

作業:
1. PlaywrightのNodeとPython経路を確認する。現在のPython sourceに対応するbrowser wheelを
   python tools/prepare_browser_wheel.pyで準備する。wheelはgenerated/gitignoredでありcommitしない。
2. 代表journeyを実ブラウザで実行する:
   - Linearを選ぶ
   - 2-record GenBank source、3-record GenBank source、single-record sourceを追加する
   - File headerのpointer操作でmiddle sourceを上/下へ動かす
   - keyboard TabとEnter/Spaceで別sourceを動かす
   - first Upとlast Downがdisabledであることを確認する
   - Generate Diagramを実行し、operation settlementを待つ
3. typed request、Result、Run Infoまたは安全なrecord metadataを使い、生成record identity/orderを
   確認する。古いResultを新しいGenerate結果として測定しない。Generate失敗時は既存Resultを
   保持する既存契約を壊していないことも確認する。
4. custom row、reverse、crop、definition、subtitle、depth、selected comparison endpointを設定した
   journeyでFileを動かし、UIDへの対応をbefore/afterで比較する。rowを暗黙に再採番しない。
   valid raw LOSAT cacheの再利用とderived numeric indexの更新を区別する。
5. Save Sessionを実行し、新しいpageでfresh Loadする。File card order、linearSeqs order、shared
   source resources、record UID、row、comparison pairを確認して再Generateする。新しいSession
   fieldやmigrationが追加されていないことをsaved JSONでも確認する。
6. GFF3+FASTA pair、同名別upload、upload直後moveからのasync multi-record discoveryを確認する。
7. viewportを390x844へ変更し、File名、Remove、Up/Down、collapsed record countが利用可能で
   横overflowや到達不能がないことを確認する。必要なら最終状態の使い捨てscreenshotを証拠にする。
8. focused evidenceを再利用しつつ、変更後に必要な広いchecksを実行する:
   node --test tests/web/architecture-contracts.test.mjs
   node --test tests/web/*.test.mjs
   npm run test:web:comparison-contracts
   node tools/check-web-change-budget.mjs --base origin/dev
   git diff --check
   Python側へ変更や失敗の波及がある場合は:
   pytest tests/ -v -m "not slow"
9. Web change policyの結果を読む。runtimeとguard/authorityを同じ変更で自己承認していないこと、
   dependency cycle、新production dependency、privileged surface expansionがないことを確認する。
10. architecture concise evidenceを完成する:
    - source grouping/order owner before/after
    - canonical path before/after
    - 削除した旧record-level input reorder path
    - user-visible verification
    - deterministic checks
    - rollback
    通常のnon-increasing changeを想定する。OE/PE/CBが増える、複数owner/pathを残す、
    compatibility pathを作る場合は自動的に完了扱いせず、ratchetの例外手続きを行う。
11. production、tests、docs、生成物のdiffを別々に最終確認する。Gallery tutorial/captureを
    実際に変更する場合だけリポジトリ指定のGallery保守skillと手順を使用する。

最終合格条件:
- 総合計画書SR-01〜SR-14が全てpassまたは明示した同等証拠で満たされる。
- File headerからsource全体を移動できる。
- source内record順、UID状態、rows、pairs、resources、Sessionが保持される。
- Input File cardsに旧record-level reorder ownerが残らない。
- linearSeqs以外のorder source of truth、別request path、compatibility branchがない。
- testsが内部actionだけでなく実ユーザー操作を通る。
- 必須gateに未解決failureがない。

外部操作:
push、PR作成・編集、merge、tag、deploy、外部投稿はこのプロンプトでは許可されない。
別途明示的に許可されている場合だけ、その対象branchと条件を確認して実行する。

完了報告:
- 結果を先に述べる
- 変更したowner/pathと利用者向け挙動
- 主要ファイルへのリンク
- SR-01〜SR-14の結果
- exact test commandsと結果
- browser/wheel source identity
- architecture concise evidenceとrollback
- 残る制約または未完了項目
- 英語のproposed commit titleと短いsummary
を記載する。必要な作業が残る場合は完了と報告しない。
```
