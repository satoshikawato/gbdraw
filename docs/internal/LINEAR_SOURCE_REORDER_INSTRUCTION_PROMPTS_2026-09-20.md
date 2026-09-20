# Linear File-row block修正 — 実装セッション用INSTRUCTION PROMPTS

このファイルは、過去のチャット履歴を持たない新規参加者へ渡す実行指示である。
管理文書は
[LINEAR_SOURCE_REORDER_MASTER_PLAN_2026-09-20.md](LINEAR_SOURCE_REORDER_MASTER_PLAN_2026-09-20.md)
である。各promptはリポジトリrootをcurrent working directoryとして使用する。

Product仕様は既に決定済みである。`PD-OI-018` scenario revision 3、
`LINEAR-FILE-ROW-BLOCK`を実装し、別のProduct結果を選ばない。authority-only PR `#549`は
`dev`へmerge済みで、merge commitは`ee0b4450`である。

S1、S2、S3の順に実行する。同じ担当者が連続実行してよい。後続担当者は会話ではなく、
repository、総合計画書第13節、Git history、test evidenceから状態を復元する。

## S1 — Domain、atomic action、UIを実装する

```text
gbdraw Web版LinearモードのFile moveを、PD-OI-018 scenario revision 3
LINEAR-FILE-ROW-BLOCKに従って実装してください。過去の会話は前提にしません。

最初に全て読む:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/LINEAR_SOURCE_REORDER_MASTER_PLAN_2026-09-20.md
5. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
6. docs/internal/PRODUCT_IMPACT_RATCHET.md
7. docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.mdのPD-OI-018 revision 3
8. docs/REFERENCE/web-app.md

開始確認:
- git status、branch、HEAD、upstreamを確認し、無関係な変更を保持する。
- runtime branchはauthority PR #549を含む最新origin/devから派生していなければならない。
- authority fileをruntime変更と同時に変更しない。
- 既存差分がある場合はproduction、tests、docsに分けて読み、正しい部分だけ再利用する。

実装する挙動:
- File headerのup/downは、一Fileの全recordsをsource blockとして移動する。
- defaultのArrange in rowsがONのnormal layoutでは、source順とvisual row順を同時に変える。
- normal layoutは、各Fileが一つのrowだけを使い、異なるFileがrowを共有しない状態である。
- normal moveでは現在のdistinct row番号を昇順にし、移動後のFile順へ割り当てる。
- File内record順とUIDに付随する全状態を保持する。
- Arrange in rowsがOFFでもnormal/customを判定する。normal moveでは、再度ONにしたときに
  File順とrow順が矛盾しないよう潜在rowも同時に再割当する。
- 一Fileが複数row、または複数Fileが共有rowのcustom layoutではFile moveをdisableする。
- custom時はRecord Layoutがvisual placementを所有することを説明し、Record Layoutへ案内する。
- blocked moveはFile順、rows、comparisons、cache metadata、Result、Historyを変更しない。
- successful moveは一つのatomic undoable draft transactionとする。
- current ResultはGenerate成功まで変更しない。

ownerと実装境界:
1. gbdraw/web/js/app/linear-sources.js
   - groupLinearSourceRecords()をsource identityの唯一ownerとして使う。
   - pure adjacent source-block transformを所有する。
   - Vue、DOM、History、rows、cache、Sessionへ依存させない。
2. gbdraw/web/js/app/linear-record-layout.js
   - sourceGroups、UID-row entries、sourceIndex、directionを受けるpure planを所有する。
   - availabilityとactionが同じplanを使えるallowed/reason/rowsを返す。
   - input arraysとrecordsを変更しない。
   - custom layoutを推測してnormal化しない。
3. gbdraw/web/js/app/app-setup.js
   - pure planとsource transformをcoordinatorとして結ぶ。
   - history.runUndoable('Move File', ...)の一transaction内でlinearSeqsとrow planを適用する。
   - applyLinearSeqMutation()の既存comparison/cache reconciliationを再利用する。
   - availabilityとactionに別の判定式を作らない。
4. gbdraw/web/index.html
   - File header control、disabled state、accessible name、custom説明だけを持つ。
   - titleまたは説明にcustom理由を示す。
   - Record options内の旧record-level File reorderを残さない。
   - Advanced Record LayoutのmoveLinearRecordWithinRow()は残す。

必ず保持するもの:
- record UID、source、selector、crop、reverse complement、definition、subtitle、depth、feature state
- File内record順
- explicit endpoint UIDとbiological endpoint
- compatible raw LOSAT evidence
- last successful Result
- GFF3とFASTAのpair binding
- pending discovery後のsource identity

追加してはいけないもの:
- fileOrderまたは別order state
- Session field、migration、request schema、Worker protocol、render path
- watcherによる後追い同期
- drag-and-drop、任意位置sort、汎用sortable framework
- custom layoutのsilent normalization
- 新しいproduction module。既存ownerで不可能なら実装を拡張せず理由を報告する。

focused unit tests:
- tests/web/linear-sources.test.mjs:
  single source block、2+3 records、boundary、same-name別upload、non-mutating、legacy interleave。
- tests/web/linear-record-layout.test.mjs:
  normal、非連続row slots、split File、shared row、layout OFF後の再ONを想定した潜在row、
  invalid boundary、non-mutating。

実行:
node --test tests/web/linear-sources.test.mjs tests/web/linear-record-layout.test.mjs
node --test tests/web/session-request.test.mjs
git diff --check

設計review:
- SOLID: source、row、coordination、UI、projectionのownerが分離されている。
- KISS: adjacent swapとnormal/custom規則だけである。
- DRY: source identity、layout classification、reconciliationが各一箇所である。
- YAGNI: state、schema、framework、future optionを追加していない。

終了条件:
- productionとpure testsが整合し、focused unit testsがpassする。
- 総合計画書第13節へbranch、変更owner、正確なcommands/results、未実施項目を記録する。
- browser、Session、full gateが残る場合は実装全体完了とは報告しない。
- push、PR、merge、tag、deployは明示的な許可がない限り行わない。
```

## S2 — Browser、History、comparison、Sessionを受け入れる

```text
LINEAR-FILE-ROW-BLOCK runtimeについて、実際のFile header操作からGenerate、History、
comparisons、Save/fresh Loadまでを検証し、見つかったin-scope不具合を修正してください。
過去の会話は前提にしません。

最初に全て読む:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/LINEAR_SOURCE_REORDER_MASTER_PLAN_2026-09-20.md
5. docs/internal/LINEAR_SOURCE_REORDER_INSTRUCTION_PROMPTS_2026-09-20.md
6. docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.mdのPD-OI-018 revision 3
7. S1のgit diffと総合計画書第13節

開始確認:
- runtime branchがmerged authority ee0b4450をancestorに持つことを確認する。
- production、tests、docsの既存差分を別々に読む。
- Node @playwright/testとPython Playwrightを確認する。
- browser wheelが必要ならpython tools/prepare_browser_wheel.pyで現在sourceから生成する。
  wheelはgitignoredでありcommitしない。

primary browser journey:
1. Linearで2-record GenBank Fileと3-record GenBank Fileをuploadする。
2. Arrange in rowsがdefaultでONであり、各Fileが一つのrowを持つことを確認する。
3. Generateし、Result、request、Adjacent pairs、raw cacheのbefore snapshotを取る。
4. 実File header buttonで下のFileを上へ動かす。内部action直接呼び出しだけで済ませない。
5. 次をassertする:
   - File cards、linearSourceGroups、linearSeqsが新順序である。
   - moved Fileの全recordsが上のrow slotへ移り、元のFileが下のrow slotへ移る。
   - File内record順とUIDが同じである。
   - Adjacent pairsが新occupied-row adjacencyから再導出される。
   - explicit endpoint testsではendpoint UIDを保持しnumeric indexだけ変わる。
   - incompatible derived comparison dataが残らない。
   - compatible raw cacheだけが残る。
   - Result bytesとdiagram run countはGenerate前に変わらない。
6. Historyが一件増え、Undoでsource順とrowsが同時に戻り、Redoで同時に進むことをassertする。
7. pointer、Space、Enterを使い、first Upとlast Downのdisabledを確認する。
8. Arrange in rowsをONのままGenerateする。OFFにしてprimary outcomeを証明してはならない。
9. request selectors、presentation.gridRow、SVG semantic record row/orderを確認する。
10. Saveし、新pageでfresh Loadする。File order、rows、UID、resources、pairsを確認して再Generateする。
11. saved sessionに新しいfileOrder fieldやschema migrationがないことを確認する。

custom layout journeys:
- 一Fileのrecordsを二rowへ分ける。
- 異なるFilesのrecordsを同じrowへ置く。
- それぞれで全File move buttonsがdisabled、説明がvisible、Record Layoutへ案内されることを確認する。
- Arrange in rowsをOFFにしてもcustom判定とdisabled説明が維持されることを確認する。
- blocked actionをdomain経路から試してもstate、History、cache、Resultが不変であることを確認する。
- Record Layoutをnormalへ戻した後、File moveが再び利用可能になることを確認する。

追加journeys:
- same-nameの別uploads
- paired GFF3+FASTA
- source move中に完了するmulti-record discovery
- selector、crop、reverse、definition、subtitle、depth、feature state
- failed/canceled/stale GenerateのResult preservation既存contracts
- viewport 390x844のoverflow、File名、Remove、move buttons、collapsed record count

主な対象test:
- tests/web/linear-multi-record.playwright.spec.js
- 必要な既存History/Result contract tests
- tests/ci/playwright-inventory.test.mjsはinventory数が実際に変わる場合だけ更新する。

focused実行:
npx playwright test tests/web/linear-multi-record.playwright.spec.js \
  --grep "File source order|source order|custom Record Layout" --workers=1 --retries=0
npm run test:web:comparison-contracts
node --test tests/web/session-request.test.mjs
git diff --check

失敗時:
- timeoutやassertionを弱めない。
- Arrange in rowsをOFFにして回避しない。
- custom placementを自動破棄しない。
- Chromium sandbox errorは同じcommandを必要な権限で再実行し、環境と実装を区別する。
- Result、request、SVGのどれを測ったかを明示し、古いResultを新しいGenerate証拠にしない。

終了条件:
- 総合計画書FR-01〜FR-19をbrowserまたは同等の自動証拠で満たす。
- production、tests、docs diffを別々にreviewする。
- 第13節へ正確なcommands/results、wheel identity、未解決事項を記録する。
- full Web/architecture gateが残る場合は全体完了とは報告しない。
- push、PR、merge、tag、deployは明示的な許可がない限り行わない。
```

## S3 — Full gates、architecture review、deliveryを完了する

```text
LINEAR-FILE-ROW-BLOCK runtimeの広い回帰、architecture fitness、Product Impact、
最終diffを検証し、許可された範囲でdeliveryを完了してください。過去の会話は前提にしません。

最初に全て読む:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/LINEAR_SOURCE_REORDER_MASTER_PLAN_2026-09-20.md
5. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
6. docs/internal/PRODUCT_IMPACT_RATCHET.md
7. docs/internal/WEB_CHANGE_POLICY.md
8. docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.mdのPD-OI-018 revision 3
9. S1/S2のgit diffと総合計画書第13節

full verification:
node --test tests/web/*.test.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/ci/*.test.mjs
npm run test:web:comparison-contracts
node tools/check-web-change-budget.mjs --base origin/dev --head HEAD
git diff --check

Python runtimeへ変更していない場合、Web/CI contractsとfocused browserを優先する。
Python側failureへ波及した場合だけpytest tests/ -v -m "not slow"を追加する。
long-running local testsは少なくとも30分を許容し、増分監視する。remote CIをpollする場合は
明示的な緊急指示がない限り10分間隔とし、10秒pollを行わない。

最終review:
1. production diff
   - source owner、row owner、atomic coordinator、UIだけか。
   - duplicate availability、watcher、fallback、dead actionがないか。
2. test diff
   - default Arrange in rows ONのprimary journeyか。
   - internal actionだけでなく実UIを通るか。
   - custom block、Undo/Redo、Adjacent、explicit endpoint、cache、Result、Sessionを区別するか。
3. docs diff
   - web-app reference、master plan、promptsがrev.3と一致するか。
   - 過去の会話を知らない読者に意味が通るか。
4. generated diff
   - wheel、screenshots、reference SVGを誤ってstageしていないか。

architecture evidence:
- owner before/after、canonical path、superseded record-level pathを記録する。
- expected ordinary non-increasing: OE 0->0、PE 0->0、CB 0->0。
- new source of truth、schema、compatibility path、dependency cycleがないことをgateで確認する。
- 増加があれば通常変更として押し切らず、ratchetの例外手続きへ移る。

Product evidence:
- runtime baseにPD-OI-018 revision 3が含まれることを示す。
- runtime PRでauthorityを変更しない。
- implementationがLINEAR-FILE-ROW-BLOCKの全Must preserve、custom rule、residual riskを満たす。

delivery:
- repository guidanceに従い、一session一commitとして英語のcommit titleとsummaryを用意する。
- commit前にbranchとupstreamを確認する。
- PR wordingを作る場合はwrite-clear-pull-request skillを読み、PR language checkerを一度通す。
- push、PR、mergeは現在のtaskで明示的に許可された場合だけ行う。
- remote mutationのretry前にremote stateを確認する。

完了報告:
- 結果を先に述べる。
- user-visible behaviorとnormal/custom ruleを短く説明する。
- 主要owner/pathと主要ファイルを示す。
- FR-01〜FR-20の結果を示す。
- exact commandsとpass/failを示す。
- architecture evidence、Product authority、rollbackを示す。
- 残るfailure、未実施、external blockerがあれば完了と表現しない。
- 総合計画書第13節を最終状態へ更新する。
```
