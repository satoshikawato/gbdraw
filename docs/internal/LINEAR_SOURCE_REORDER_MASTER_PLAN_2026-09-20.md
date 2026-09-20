# Linear File順序と描画row順序の整合 — 総合計画書

状態: runtime実装・local verification完了

対象base: `origin/dev`の`ee0b44502db14e7799b007ad07cb23b30db4ede8`

Product authority: `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`の
`PD-OI-018` scenario revision 3、選択肢`LINEAR-FILE-ROW-BLOCK`

Authority経路: PR `#549`、merge commit `ee0b4450`

実装セッションの開始指示は、別紙
[LINEAR_SOURCE_REORDER_INSTRUCTION_PROMPTS_2026-09-20.md](LINEAR_SOURCE_REORDER_INSTRUCTION_PROMPTS_2026-09-20.md)
にある。本書は、過去の会話を読んでいない新規参加者が、症状、原因、決定済み仕様、
設計、検証、完了条件を復元できるようにする管理文書である。

## 1. 対象と用語

gbdraw Web版のLinearモードは、複数のGenBank、またはGFF3とFASTAの組を
一つの図へ描画する。

| 用語 | 意味 |
| --- | --- |
| File / source | 一回のGenBankアップロード、または一組のGFF3とFASTA。Input Genomesでは`File 1`などのカードになる。 |
| record | File内の一つの配列。染色体、plasmid、contigなど。一つのFileに複数recordを含められる。 |
| `linearSeqs` | Linear recordの正準な順序付き配列。安定UID、source、selector、crop、reverse、表示情報、depthなどを保持する。 |
| `linearSourceGroups` | 同じsourceを参照するrecordsをFileカードへまとめる導出値。永続する別order stateではない。 |
| Record Layout | `Arrange in rows`とrecordごとのrow番号、同一row内順序を編集する高度設定。 |
| normal layout | 各Fileがちょうど一つのrowを占有し、異なるFile同士がrowを共有しない配置。 |
| custom layout | 一つのFileが複数rowに分かれる、または異なるFileが同じrowを共有する配置。 |
| current Result | 最後に成功したGenerateの不変成果物。draft編集だけでは置換しない。 |

## 2. 利用者から見える不具合

典型例では、`Vibrio harveyi`を追加した後、Fileカードの上矢印で`File 1`へ
移動しても、Generate後の図では下段に残る。

File番号と矢印は視覚順序の操作に見えるため、利用者は次を期待する。

1. Fileカードを上へ動かすと、そのFileに属する全recordが一つのblockとして上へ動く。
2. defaultの`Arrange in rows`がONでも、Fileカード順と描画row順が一致する。
3. 次のGenerateが成功するまで、現在表示中のResultは変わらない。

## 3. 原因監査

### 3.1 最初の回帰

multi-record対応で入力UIの反復単位がrecordからFileへ変更された際、Fileカードは
`linearSourceGroups`を使うようになった。一方、順序操作はrecord単位のまま折りたたみ内へ
残り、Fileヘッダーにsource単位操作が作られなかった。

これは「一つのFileはrecord数にかかわらず一つのFileカード」という正しい変更に、
同じ粒度の並べ替え操作を追加しなかったことによる状態・UI粒度の不一致である。

### 3.2 最初の修正候補が図へ反映されなかった原因

旧実装`53d2f3c6`（PR `#548`でmerge）は、File矢印からsource blockを
`linearSeqs`内で移動できるようにした。
しかしrecord UIDに対応する`linearRecordRows`の絶対row番号を保持した。

`Arrange in rows`がONのとき、描画順はrecord配列順だけでなくrow番号で決まる。このため、
Fileカードを`File 3`から`File 1`へ動かしても、recordsが以前の下段rowを保持し、図では
下段に残った。UI上のsource順と描画row順を別々に更新したことが直接原因である。

同候補の主acceptance testはGenerate直前に`Arrange in rows`をOFFにしていた。そのtestは
source配列順が変わることは証明したが、default-layoutでFile moveが図のrow順へ反映される
という主要journeyを検証していなかった。

### 3.3 バグを生みやすくした構造

- source orderは`linearSeqs`、row placementは`linearRecordRows`に分かれている。
- File move actionがsource orderだけを変更し、row ownerへ明示的な計画を要求しなかった。
- UI actionの可否と実行がnormal/custom layoutを共通判定していなかった。
- testがdefaultを変更してからGenerateし、利用者のdefault journeyを迂回した。

Python renderer、Worker、typed requestが順序を失ったことが原因ではない。Web draftから
canonical requestへ渡す前の状態更新が不完全だった。

## 4. 決定済みProduct仕様

`PD-OI-018` scenario revision 3は次を要求する。

### 4.1 normal layout

- File番号とup/downはsource順とvisual row順を同時に変更する。
- 一つのFileに属する全recordを一つのblockとして移動する。
- File内record順とrecord固有状態を保持する。
- Fileごとに一つのrowを維持し、File rowsはFileカード順に並ぶ。
- 操作は一つのatomicかつundoableなdraft transactionである。

File move前に使われていたFile-owned row番号の集合は保持してよい。例えばFile rowsが
`1, 3, 6`なら、その三つを新しいFile順へ割り当てる。重要なのは数値の連続性ではなく、
新しいFile順とrowの昇順が一致することである。

### 4.2 custom layout

一つのFileが複数rowに分かれている、または異なるFileが一つのrowを共有している場合、
File moveを無効にする。画面はcustom Record Layoutがvisual placementを所有していることを
説明し、Advanced comparison and layoutのRecord Layoutへ案内する。

blocked moveはFile順、row、comparison、cache metadata、current Resultのどれも変更しない。
custom placementを推測して書き換えるnormalizationは行わない。

### 4.3 comparisons、cache、Result

- explicit comparison endpointは安定record UIDに付随し、生物学的endpointを変えない。
- numeric record indexは新しい`linearSeqs`順で再解決する。
- Adjacent comparisonは移動後のoccupied-row adjacencyから導出し直す。
- derived comparison artifactsは無効化する。
- source、region、検索設定、endpointが互換なraw LOSAT evidenceだけを再利用する。
- current Resultは次のGenerate成功まで保持する。
- failed、canceled、superseded、stale Generateはcurrent Resultを置換しない。

### 4.4 persistenceとcompatibility

Save、fresh Load、regeneration、keyboard operation、Session replayでFile順とrow順を保持する。
新しいFile-order state、Session field、migration、request schema、Worker protocol、rendering
pathは追加しない。

## 5. 目標アーキテクチャ

```text
index.html: File up/downとcustom-layout説明
  -> app-setup.js: availabilityと一つのHistory transaction
       -> linear-sources.js: source blockの純粋な隣接交換
       -> linear-record-layout.js: normal/custom判定とrow再割当計画
       -> applyLinearSeqMutation(): 既存の一括reconciliation
            -> linearSeqs / linearRecordRows
            -> comparison plan reconciliation
            -> derived artifact invalidation / raw cache reindex
            -> 既存session-request / Worker / Python renderer
```

### 5.1 owner表

| 責務 | Owner | 境界 |
| --- | --- | --- |
| source identityとgrouping | `gbdraw/web/js/app/linear-sources.js` | File object/resource identityを既存規則でgroup化する。 |
| source block変換 | `linear-sources.js` | DOM、Vue、History、row、cacheへ依存しない。 |
| normal/custom判定とrow move plan | `gbdraw/web/js/app/linear-record-layout.js` | source groupとUID-rowだけを受ける純粋関数。引数を変更しない。 |
| atomic coordination | `gbdraw/web/js/app/app-setup.js` | 同じplanをavailabilityとactionで使い、既存mutation境界を一度呼ぶ。 |
| UI | `gbdraw/web/index.html` | event forwarding、disabled状態、説明だけを持つ。domain規則を複製しない。 |
| request/session projection | `gbdraw/web/js/services/session-request.js` | 変更しない。既存`linearSeqs`とrow projectionを使う。 |
| rendering | Worker/Python既存経路 | 変更しない。 |

### 5.2 状態規則

- `linearSeqs`を唯一のFile/record順序source of truthとする。
- `linearSourceGroups`は常に導出する。別のFile order配列は作らない。
- `linearRecordRows`は既存UID-row mapであり、新形式を作らない。
- normal moveでは現在のdistinct row slotsを昇順にし、移動後のFile順へ割り当てる。
- `Arrange in rows`がOFFでもrowはSessionに残り、再度ONにできる。そのためnormal/custom
  判定はON/OFFに依存せず、normal moveでは潜在rowも同時に再割当する。
- availabilityと実行は同じpure planを使い、enabledなのにno-opとなる分岐を作らない。

## 6. SOLID、KISS、DRY、YAGNI

### SOLID

- Single Responsibility: source grouping、row semantics、coordination、UI、projectionを既存ownerへ分離する。
- Open/Closed: schemaやrendererを変更せず、既存mutation pipelineへrow planを入力する。
- Liskov Substitution: 新しいclass階層を作らないため対象外とし、抽象化の口実にしない。
- Interface Segregation: pure plannerは必要なgroups、rows、index、directionだけを受ける。
- Dependency Inversion: UIはDOM上でgroupingを再計算せず、app actionとpure domain functionsへ依存する。

### KISS

- 隣接up/downだけを扱う。
- normalではrow slotを交換し、customでは明示的にblockする。
- drag-and-drop、任意位置移動、silent normalizationは実装しない。

### DRY

- source identityは`groupLinearSourceRecords()`一箇所で決める。
- normal/custom判定はrow planner一箇所で決め、UIとactionで複製しない。
- comparison/cache reconciliationは`applyLinearSeqMutation()`の既存経路を使う。

### YAGNI

- `fileOrder`、新watcher、Session migration、sortable frameworkを追加しない。
- custom layoutの自動変換や確認dialogを追加しない。
- 将来用の汎用transaction managerや新render pathを作らない。

## 7. 実装作業

1. `linear-sources.js`のpure source-block moveを維持する。
2. `linear-record-layout.js`へpure move planを追加する。
3. `app-setup.js`でplanをavailabilityとactionに共有する。
4. source順とrow planを一つの`history.runUndoable('Move File', ...)`で適用する。
5. `index.html`でFile header controlsを表示し、custom layoutではdisableと説明を出す。
6. Record options内の旧record-level input reorderを削除する。Record Layout内の
   `moveLinearRecordWithinRow()`は残す。
7. primary browser testは`Arrange in rows`をONのままGenerateする。
8. docs、unit、browser、session、history、comparison/cache testsを更新する。

## 8. 受入条件

| ID | 条件 |
| --- | --- |
| FR-01 | single-record File BをAの上へ動かすと、File cards、`linearSeqs`、rows、requestがB,A順になる。 |
| FR-02 | 2-record Fileと3-record Fileを入れ替えると、全recordsがblockで動き、File内順序とUIDが維持される。 |
| FR-03 | defaultの`Arrange in rows`をOFFにせずGenerateし、requestの`gridRow`とSVGのrow順がFile順に一致する。 |
| FR-04 | distinct row slotsが非連続でも、新しいFile順へ昇順で再割当される。Arrange in rowsがOFFでも潜在rowを同じように更新する。 |
| FR-05 | 一File複数rowではmove buttonsが無効になり、Record Layoutへの説明が表示される。 |
| FR-06 | 複数File共有rowでも同様にblockし、状態は完全なno-opである。 |
| FR-07 | blocked moveはHistory、comparison、cache、Resultを変更しない。 |
| FR-08 | successful moveはHistory一件であり、Undo/Redoがsource順とrowsを同時に戻す。 |
| FR-09 | explicit endpointsはUIDを維持し、numeric indexesだけ再解決する。 |
| FR-10 | Adjacent pairsは新row adjacencyから再導出され、旧derived artifactを使わない。 |
| FR-11 | compatible raw cacheは再利用され、incompatible entryは使われない。 |
| FR-12 | selector、crop、reverse、definition、subtitle、depth、feature stateがUIDに付随する。 |
| FR-13 | File moveだけではcurrent Resultとrun countが変わらず、Generate成功で新配置へ置換する。 |
| FR-14 | failed/canceled/stale Generateが旧Resultを保持する既存契約に回帰がない。 |
| FR-15 | Save、fresh Load、再GenerateでFile順、rows、resources、pairsを保持し、新Session fieldがない。 |
| FR-16 | same-name upload、GFF3+FASTA、pending multi-record discoveryがsource identityを失わない。 |
| FR-17 | pointer、Enter、Spaceで操作でき、boundary/custom disabledと一意なaccessible nameを持つ。 |
| FR-18 | 390px幅でFile名、Remove、move controls、record countに横overflowや到達不能がない。 |
| FR-19 | Input cardsの旧record-level reorder pathがなく、Advanced Record Layoutのrecord操作は残る。 |
| FR-20 | 新order state、schema、migration、Worker/render path、dependency、cycleを追加しない。 |

## 9. Test strategy

### 9.1 Pure contracts

- `tests/web/linear-sources.test.mjs`: source block、identity、same-name、boundary、legacy interleave。
- `tests/web/linear-record-layout.test.mjs`: normal plan、非連続row slots、split/shared custom、非破壊性、layout OFF。

### 9.2 Browser contracts

`tests/web/linear-multi-record.playwright.spec.js`で実File header buttonを操作する。

- primary 2+3-record journeyはdefault row layoutを維持する。
- before/afterのUID、rows、Adjacent pairs、cache、Result、Historyを検証する。
- Generate requestのselectorsと`presentation.gridRow`、SVG semantic row metadataを検証する。
- Save、fresh Load、再Generateを同じjourneyで検証する。
- custom split/shared layoutのdisable、説明、no-opを検証する。
- GFF3+FASTA、pending discovery、keyboard、mobileを検証する。

内部`moveLinearSource()`を直接呼ぶtestだけではUI acceptanceにしない。

### 9.3 Gates

最低限、次を実行する。

```bash
node --test tests/web/linear-sources.test.mjs tests/web/linear-record-layout.test.mjs
node --test tests/web/session-request.test.mjs
npx playwright test tests/web/linear-multi-record.playwright.spec.js --grep "File source order|source order|custom Record Layout" --workers=1 --retries=0
npm run test:web:comparison-contracts
node --test tests/web/*.test.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/ci/*.test.mjs
node tools/check-web-change-budget.mjs --base origin/dev --head HEAD
git diff --check
```

browser wheelが必要な場合は`python tools/prepare_browser_wheel.py`で現在sourceから作る。
wheelはgitignored生成物でありcommitしない。Chromium sandbox failureは同じcheckを必要な
権限で再実行して環境制約と実装失敗を区別する。

## 10. Architecture fitness evidence

想定はordinary non-increasing changeである。

- semantic owners: 既存`linear-sources.js`と`linear-record-layout.js`へ責務を置き、新ownerを作らない。
- canonical path: UI → app coordinator → existing mutation → existing request/Worker/render pathの一つ。
- superseded path: Input Record optionsのrecord-level File reorder。
- compatibility burden: 追加なし。
- expected changed-scope: `OE 0 -> 0`、`PE 0 -> 0`、`CB 0 -> 0`。

実測でowner/path/compatibilityが増える場合は完了扱いにせず、
`ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`の例外条件を適用する。

## 11. Workflowと停止条件

1. branchはmerged authorityを含む最新`origin/dev`から作る。
2. authorityとruntimeを同じPRにしない。authority PR `#549`は既にmerge済みである。
3. focused testsを先に通し、次にbrowser、広いunit、policy gateへ進む。
4. production、tests、docs、generated diffを別々にreviewする。
5. timeoutやassertionを弱めてpassさせない。
6. remote CIを監視する場合、明示的な緊急指示がない限りpoll間隔は10分とする。

次の場合だけ停止して判断を求める。

- normal/custom定義を変更する必要がある。
- custom layoutを自動変換する必要が生じる。
- Session/request/Worker schema変更が必要になる。
- current ResultをFile move時にlive updateする必要が生じる。
- signed decisionと両立しないuser-visible outcomeが残る。

## 12. Rollback

runtime rollbackは、File header controls、row planner、atomic coordinator、関連tests/docsを
一単位でrevertする。authority`PD-OI-018` revision 3はruntime rollbackだけでは取り消さない。
Session migration、resource変換、reference SVGのrollbackは不要である。

## 13. 実施記録

| 項目 | 記録 |
| --- | --- |
| authority | PR `#549`を`dev`へmerge。merge commit `ee0b4450`。 |
| runtime branch | `fix/linear-file-row-block-20260920`、`ee0b4450`から作成。 |
| prior runtime | PR `#548`の`53d2f3c6`がbaseにmerge済み。source-block UIは維持し、絶対row保持をrev.3 semanticsへ置換した。 |
| focused unit | `node --test tests/web/linear-sources.test.mjs tests/web/linear-record-layout.test.mjs tests/web/session-request.test.mjs`: 3 pass。 |
| primary browser | default Arrange in rows ONのFile move、History、Generate、Save/fresh Load、custom block: 1 pass。 |
| comparison browser | `npm run test:web:comparison-contracts`: 16 pass。 |
| Web unit | architectureを独立実行し137 pass、その他の`tests/web/*.test.mjs`は472 pass。 |
| CI contracts | `node --test tests/ci/*.test.mjs`: 59 pass。 |
| policy | `node tools/check-web-change-budget.mjs --base origin/dev`: Gate PASS、Review REQUIRED（pure planner exportとcomputed availabilityの可視化）。owner/path/schemaのblocking violationなし。 |
| architecture | 既存owners内のpure planningとcoordination。`OE 0 -> 0`、`PE 0 -> 0`、`CB 0 -> 0`。新state、schema、migration、Worker/render path、dependency cycleなし。 |
| browser wheel | current sourceで作成。SHA-256 `413b6f370b801ae90c863f9f5835a54cd7deb282c8b4546508fb44308614bab2`。gitignoredでcommit対象外。 |
| runner note | architectureを含む全Web unitの一括並列呼び出しはrunner競合で停止したため終了。architecture 137件と残り472件を分離し、同じsourceで全件passさせた。 |

実施結果はacceptanceを弱めず、この表を更新する。
