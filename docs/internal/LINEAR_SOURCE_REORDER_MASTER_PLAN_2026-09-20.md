# Linear入力File順序変更の復旧 — 総合計画書

状態: 実装前。2026-09-20時点の`dev`、`eec715c70323c1dc6173a36bda2f18e600f8384d`を監査基準とする。実装開始時には最新の`origin/dev`で原因と契約を再確認する。

実装セッションの開始指示は、別紙の[INSTRUCTION PROMPTS](LINEAR_SOURCE_REORDER_INSTRUCTION_PROMPTS_2026-09-20.md)にある。本書は、前提知識のない実装者が問題、既存契約、設計、作業範囲、受入条件を復元できるようにする管理文書である。

## 1. 背景と用語

gbdrawは、GenBankまたはGFF3とFASTAからゲノム図を生成するPythonソフトウェアである。Web版はVueの単一ページアプリであり、Linearモードでは複数の入力sourceと、それぞれに含まれる一つ以上の生物学的recordを一つの図へ配置できる。

| 用語 | 本書での意味 |
| --- | --- |
| source / File | 一回のGenBankアップロード、または一組のGFF3とFASTA。画面では`File 1`、`File 2`のカードとして表示する。 |
| record | source内の一つの配列レコード。染色体、プラスミド、contigなど。一つのFileに複数含まれ得る。 |
| `linearSeqs` | Linearで選択されたrecordの正準な順序付きリアクティブ配列。各要素は安定UID、source参照、selector、表示、depthなどを持つ。 |
| `linearSourceGroups` | 同じsourceを参照する`linearSeqs`要素をFileカードへまとめる導出ビュー。永続状態ではない。 |
| Record Layout | recordごとのrow番号と、同じrow内の左右順を編集する高度設定。File順序とは別の利用者所有状態。 |
| canonical request | Web状態から構築し、Worker、Python描画、Session保存で共有するschema 7の型付き描画要求。 |

## 2. 問題

### 2.1 利用者から見える症状

LinearのInput Genomesには`File 1`、`File 2`、`File 3`という順序付きカードが表示されるが、Fileカードのヘッダーには順序変更ボタンもドラッグ操作もない。そのため、アップロード後にFile単位で順序を変更できない。

record用の`Up`と`Down`は残っているが、各Fileの`Record options`内に隠れている。multi-record sourceでは、最初に`Number of records: N`、次に各recordの`Record options`を開かなければ見えない。さらに、この操作が移動するのはFile全体ではなく一つのrecordである。

### 2.2 影響範囲

- 対象はWeb版Linearモードの入力source順序である。
- Circular、CLI、Python APIの入力順序操作は直接の原因ではない。
- Python描画やWorkerが順序を失っているのではない。canonical requestは`linearSeqs`順にrecordsを構築する。
- sourceをアップロード時から正しい順に並べた場合、生成はその状態を使用できる。
- Record Layoutのrow編集は利用可能だが、Fileカードの順序変更を代替するものではない。

### 2.3 安全な暫定回避策

修正前は、sourceを最初から必要な順で追加する。生成図の上下rowまたは同じrow内の左右配置だけを変えたい場合は、Advanced comparison and layoutのRecord Layoutを使用する。multi-record Fileに対して、隠れたrecord単位のUp/DownをFile移動の代用にしてはならない。

## 3. 原因

### 3.1 回帰を導入した変更

2026-09-14の`77563b6bf40b93a9637ea6db1a605380bca96640`、`Restore complete multi-record comparisons and default layouts (#526)`に含まれる「Keep one file-input card per multi-record source」で、入力表示の反復単位がrecordからsourceへ変更された。

変更前は各`linearSeqs`要素のカード上部にrecord番号とUp/Downが常時表示されていた。変更後は`linearSourceGroups`ごとにFileカードを作り、旧Up/Downを各recordの折りたたみ領域へ移した。Fileカードのヘッダーには対応するsource単位操作が追加されなかった。

この変更自体の目的は正しい。一つのmulti-record sourceをrecord数だけ重複したファイルアップローダーとして表示しないためであり、`PD-OI-018`の「一つのsourceは一つのFileカード」という契約を実現した。回帰は、表示単位を変えたときに並べ替え操作の対象単位を変えなかったことにある。

### 3.2 状態と表示の粒度不一致

現在の責務は次のように分かれている。

| 層 | 現在の単位 | 問題 |
| --- | --- | --- |
| Fileカード | source | 順序操作がない。 |
| `linearSeqs` | record | 正準順序は保持できる。 |
| `moveLinearSeqUp/Down` | record一件 | source全体を移動できない。 |
| `linearSourceGroups` | sourceへ再集約した導出値 | 非連続になった同一sourceのrecordも一つのカードへ再集約する。 |
| canonical request | `linearSeqs`のrecord順 | 画面のカード順と内部record順が一致しない状態を表現し得る。 |

例えば`linearSeqs = [A1, A2, B1, B2]`でA2だけを下へ移すと、内部順序は`[A1, B1, A2, B2]`になる。一方、FileカードはA1とA2を再集約してFile Aへ表示する。この状態ではFileカードを上から読んだ順序と生成に使うrecord順が一致しない。

### 3.3 テストが検出しなかった理由

- `tests/web/linear-multi-record.playwright.spec.js`のsource-cardテストは、Fileカード数、record数、折りたたみ、置換、削除を検証するが、Fileヘッダーから順序を変えない。
- テスト名に`through record moves`を含むケースは、実際にはrow番号を変更しており、record配列の順序もFile順序も変更していない。
- 別のテストは`window.__GBDRAW_APP__.moveLinearSeqUp()`を直接呼ぶ。内部関数を呼べることは確認するが、利用者が到達できるFile単位の操作は確認しない。
- `tests/web/linear-sources.test.mjs`はsource groupingとLOSAT batchingを確認するが、source groupの並べ替え契約を持たない。

## 4. 既存のProduct契約と採用する挙動

### 4.1 Authority

適用する既存authorityは`docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`の`PD-OI-018: Complete Linear records, placement, and comparison scope`、scenario revision 2と`OIC-015`である。

この契約は次を要求する。

- File数とrecord数を区別する。
- 一つのアップロードsourceはrecord数にかかわらず一つのFileカードを持つ。
- record placementは独立して編集可能である。
- Save、fresh Load、regeneration、reordering、cache reuseでrecord identity、placement、pair mapping、shared source resourcesを維持する。

本修正は新しいProduct結果を選ぶ変更ではなく、既存の順序変更continuationをsource単位UIへ復旧する`IMPLEMENT_EXISTING_AUTHORITY`として扱う。

### 4.2 完成後の利用者向け契約

1. Fileカードのヘッダーに、常時発見可能なsource単位の上移動と下移動を置く。
2. sourceを移動すると、そのsourceに属する全recordが連続した一つのblockとして移動する。
3. source内のrecord順序は変えない。
4. 先頭sourceの上移動と末尾sourceの下移動は無効である。無効操作は状態を変更しない。
5. 操作後、`File 1`、`File 2`などの番号とカードのDOM順は新しいsource順を表す。
6. record UID、selector、crop、reverse complement、definition、subtitle、depth、明示的comparison endpoint、raw cache identity、source bytesを対応するrecordと共に保持する。
7. Record Layoutのrow番号はrecord UIDに付随した利用者所有状態として保持する。File移動はrow番号を暗黙に付け替えない。
8. Arrange in rowsが無効な場合、canonical record順は移動後のsource block順になる。Arrange in rowsが有効な場合、rowは明示配置を支配し、同じrow内の左右順では移動後のrecord順を使用する。
9. 入力設定はdraftである。既存ResultはFile移動だけでは書き換えず、次のGenerate成功時に新しい順序を反映する。
10. Save Sessionはdraftとcanonical requestの既存規則に従い、新規schemaや移行処理なしで順序を往復する。
11. 操作はpointerとkeyboardで利用でき、各ボタンには対象Fileと方向を含む一意なaccessible nameを付ける。

### 4.3 非目標

- drag-and-drop並べ替えを導入しない。
- `fileOrder`、source用UID配列、並行する順序stateを追加しない。
- Session schema、canonical request schema、Worker protocol、Python APIを変更しない。
- File移動に連動してRecord Layoutのrow番号を自動再採番しない。
- multi-record source内の任意recordを別sourceの間へ配置する新機能を作らない。
- History、Undo/Redo、preview live edit、comparison仕様をこの修正だけのために拡張しない。
- Fileカード全体、Input Genomes全体、Record Layout全体を再設計しない。
- 既存のsource groupingをファイル名比較へ変更しない。同名の別アップロードは別sourceのままにする。

## 5. 目標アーキテクチャ

### 5.1 正準経路

```text
index.htmlのFileヘッダーボタン
  -> app-setup.jsのsource move action
  -> linear-sources.jsの純粋なsource-block順序変換
  -> app-setup.jsの既存applyLinearSeqMutation()
  -> state.jsのlinearSeqs（唯一の順序source of truth）
  -> session-request.jsの既存canonical request構築
  -> Worker / Python render / Session save
```

新しい永続状態、watcher、生成経路を追加しない。

### 5.2 責務

| 責務 | 既存または目標owner | 方針 |
| --- | --- | --- |
| source identityとrecord grouping | `gbdraw/web/js/app/linear-sources.js` | 既存ownerを維持する。 |
| source blockの純粋な順序変換 | `gbdraw/web/js/app/linear-sources.js` | groupingと同じ意味ownerへ追加する。DOMやreactive stateへ依存させない。 |
| reactive mutationと下流reconciliation | `gbdraw/web/js/app/app-setup.js::applyLinearSeqMutation` | 既存経路を再利用する。comparison、row、cacheの別reconcileを追加しない。 |
| visible control | `gbdraw/web/index.html` | 表示とevent forwardingだけを持つ。 |
| record rowと同一row内配置 | `gbdraw/web/js/app/linear-record-layout.js`と既存app action | source移動ロジックへ混ぜない。 |
| canonical request / Session projection | `gbdraw/web/js/services/session-request.js` | 変更しない。新しいFile順序モデルを持たせない。 |

### 5.3 Architecture ratchet

通常のnon-increasing changeとして設計する。

- source grouping/order capabilityの意味ownerは`linear-sources.js`一つとする。
- `app-setup.js`はownerを複製せず、純粋な結果を既存mutation境界へ渡すcoordinatorとする。
- 旧record-level input-card reorder pathを同じ変更で削除する。
- canonical request pathは一つのまま維持する。
- compatibility pathを追加しない。
- 想定するchanged-scopeは`OE 0 -> 0`、`PE 0 -> 0`、`CB 0 -> 0`である。実装時に別owner、別path、互換分岐が必要になった場合は、この通常計画を拡張せず、Architecture Fitness Function Ratchetの例外手続きへ移る。

## 6. 実装設計

### 6.1 純粋なsource-block変換

`linear-sources.js`に小さな純粋関数を追加する。名前は実装時の既存命名へ合わせるが、契約は次のとおりとする。

```text
入力:
  ordered sequences
  source group index
  direction (-1 または +1)

処理:
  groupLinearSourceRecords(sequences)を一度呼ぶ
  対象groupと隣接groupを交換する
  各group.recordsのsequenceを元の順でflattenする

出力:
  新しいsequence配列
```

境界外、整数でないindex、不正なdirectionは変更なしとして扱う。関数は引数配列やsequence objectを直接変更せず、File API、Vue、DOM、cache、Sessionへ依存しない。

この関数にsource identity判定を再実装しない。必ず既存`groupLinearSourceRecords()`の結果を使う。

### 6.2 reactive action

`app-setup.js`にsource indexとdirectionを受ける単一actionを追加する。

- 可否判定と実行で同じsource-group domainを使用する。
- 実行結果を`applyLinearSeqMutation(next, { preserveLosatCacheInfo: true })`へ渡す。
- raw LOSAT cacheはsource内容と検索設定が同じ場合に再利用可能とし、derived comparison artifactsとnumeric indexesは既存reconciliationに任せる。
- pending record discoveryは安定UIDを維持する。source移動中にdiscoveryが完了しても、現在indexをUIDで再解決する既存処理を壊さない。
- 空の入力カードもそのUIDをsource identityとして移動できる。

`canMoveLinearSeqUp/Down`、`reorderLinearSeqs`、`moveLinearSeqUp/Down`はInput File順序のproduction consumerを失うため削除候補とする。テスト専用exportとして残さない。recordの同一row内移動は既存`moveLinearRecordWithinRow()`へ収束させる。

### 6.3 FileカードUI

`File N`ヘッダー右側へ二つのbuttonを追加する。

- desktopと狭幅sidebarでFile名やRemove操作を押し出さないcompactな配置にする。
- iconだけを使う場合も`aria-label="Move File 2 up"`のような名前とtooltip/titleを持たせる。
- 先頭のUp、末尾のDownを`disabled`にする。
- Fileが一つだけでもレイアウトを不自然に変えない。
- multi-record listの初期collapsed契約は維持する。

各recordの`Record options`内にある旧Up/Downは削除する。Record Layoutにあるrecord配置操作は残す。二つの異なる階層に同名のUp/Downを残さず、File順序とrecord placementの操作対象を明確にする。

### 6.4 Sessionと既存データ

新しいSession fieldやreader migrationは追加しない。`linearSeqs`順が既存のcanonical record順であり、現在のwriterとreaderをそのまま使用する。

現行バグの隠れたrecord操作によって、同一sourceのrecordが`linearSeqs`内で非連続になったSessionが存在する可能性がある。Loadだけではその順序を自動変更しない。利用者がFile移動を実行したときは、画面上で既に一つに見えているsourceを一つの連続blockとして移動する。これは明示操作の結果であり、互換migrationではない。

### 6.5 文書

`docs/REFERENCE/web-app.md`のLinear入力説明を、Fileカードとrecord controlsの階層に合わせる。次を短く明記する。

- Fileヘッダーの操作はsource全体を移動する。
- multi-record Fileでは全recordが一緒に移動する。
- Record Layoutはrowと同じrow内配置を担当し、File移動は明示rowを再採番しない。

新しい説明ページは作らない。Gallery tutorialやスクリーンショットがこの操作を説明している場合だけ、既存のGallery保守手順で更新する。

## 7. SOLID・KISS・DRY・YAGNIの適用

### SOLID

- Single Responsibility: source groupingとsource-block順序変換は同じdomain ownerに置き、reactive mutation、UI、record layout、request projectionを分離する。
- Open/Closed: canonical requestとSession schemaを変更せず、既存mutation境界へ一つの変換を追加する。
- Liskov Substitution: 継承やsubtypeを導入しないため対象外。SOLIDを不要なclass階層追加の理由にしない。
- Interface Segregation: 純粋関数はsequence配列、source index、directionだけを受け、app state全体を受け取らない。
- Dependency Inversion: UIはgrouping実装を再現せず、app actionへ依存する。coordinatorはDOMではなく純粋なdomain変換へ依存する。

### KISS

- Up/Downの隣接交換だけを実装する。
- drag lifecycle、drop target、pointer capture、animation、別order stateを導入しない。
- 既存の`applyLinearSeqMutation()`とUID reconciliationを使う。

### DRY

- source identity判定は`groupLinearSourceRecords()`一箇所に保つ。
- comparison、row、depth、cache reconciliationをsource move用に複製しない。
- record order操作をInput FileカードとAdvanced Record Layoutの両方で所有しない。

### YAGNI

- 任意位置drop、複数選択、一括sort、alphabetical order、永続`fileOrder`を追加しない。
- branch-onlyなSession形式の互換readerを作らない。
- 将来の一般sortable frameworkを先行導入しない。

## 8. 変更対象

| ファイル | 予定する責務 |
| --- | --- |
| `gbdraw/web/js/app/linear-sources.js` | source blockの純粋な隣接移動。 |
| `gbdraw/web/js/app/app-setup.js` | source move actionを既存mutation境界へ接続し、旧record-level input reorder actionを除去。 |
| `gbdraw/web/index.html` | Fileヘッダーの操作と旧Record options内操作の除去。 |
| `tests/web/linear-sources.test.mjs` | 純粋変換、不変条件、境界、同名別sourceの単体テスト。 |
| `tests/web/linear-multi-record.playwright.spec.js` | 実UI、request、Session、accessibilityの回帰テスト。 |
| `tests/web/session-request.test.mjs` | 既存coverageで不足する場合だけ、並べ替え後のprojectionを追加。 |
| `docs/REFERENCE/web-app.md` | File順序とRecord Layoutの利用者向け説明。 |

新規production moduleは作らない。実装調査で追加ファイルが必要になった場合は、既存ownerへ置けない理由とArchitecture ratchet上の分類を先に記録する。

## 9. 受入条件

| ID | シナリオ | 合格条件 |
| --- | --- | --- |
| SR-01 | single-record source A、B、CでBを上へ移動 | Fileカード、`linearSourceGroups`、`linearSeqs`、canonical requestがB、A、C順になる。 |
| SR-02 | Aが2 records、Bが3 recordsでBを上へ移動 | `B1,B2,B3,A1,A2`になり、source内順序とUIDを維持する。アップローダーは2個のままである。 |
| SR-03 | 同名の別ファイル二つ | 別sourceとして移動し、誤って一つへgroup化しない。 |
| SR-04 | 先頭Up、末尾Down、不正index | buttonは無効、domain actionはno-op、stateとcacheを変更しない。 |
| SR-05 | multi-record sourceの初期表示 | `Number of records: N`はcollapsedのまま。File moveは開閉を要求しない。 |
| SR-06 | record固有状態 | selector、crop、reverse、definition、subtitle、depth fileが対応UIDと共に移動する。 |
| SR-07 | Record Layoutとselected comparisons | UIDからrowへの対応と明示endpointを維持し、numeric indexだけ新しい順序へ再解決する。 |
| SR-08 | LOSAT cache | semantically equivalentなraw cacheだけ再利用し、derived pair/index表示は新しい順序と一致する。 |
| SR-09 | Generate | 利用者操作でFileを移動してGenerateすると、実際に送るtyped requestとResultのrecord identity/orderが一致する。移動前ResultはGenerate成功まで保持する。 |
| SR-10 | Save、fresh Load、再Generate | source順、shared source resources、UID、row、pair mappingを往復する。新規Session fieldはない。 |
| SR-11 | keyboardと狭幅画面 | File move buttonへTab移動でき、Enter/Spaceで動作する。accessible nameが一意で、390px幅でも操作可能である。 |
| SR-12 | 役割分離 | Input Fileカード内の旧record-level Up/Downがなく、Record Layoutのrecord配置操作は動作する。 |
| SR-13 | async discovery | upload直後のsource移動と、その後のmulti-record discoveryでsource block位置、UID、record全件を維持する。 |
| SR-14 | GFF3 + FASTA | pairを一sourceとして全record一緒に移動し、二つのfile bindingを維持する。 |

## 10. テスト戦略

### 10.1 Unit

`linear-sources.test.mjs`でdomain変換を直接検証する。

- single-record groupsの隣接交換。
- 2-recordと3-record groupのblock交換。
- source内順序とobject identityの保持。
- 入力配列の非破壊性。
- 先頭、末尾、不正direction/indexのno-op。
- 同じファイル名だが異なるFile objectは別group。
- Session resource descriptorを共有するrecordsは一group。
- 既存の非連続source recordsを明示移動したとき、各sourceを連続blockへする。

### 10.2 Browser contract

`linear-multi-record.playwright.spec.js`に、内部app methodを直接呼ばず、Fileヘッダーのbuttonをclickまたはkeyboard activateするテストを追加する。既存の2-record + 3-record fixtureを再利用する。

browser assertionは次を同じjourneyで結び付ける。

- visible File orderとfilename。
- `linearSeqs`のUID/selector順。
- source card数とrecord card数。
- canonical requestのrecordKeyとsource resource association。
- row、selected comparison、cache info。
- Save、page reload、fresh Load、Generate。
- mobile viewportとaccessible button state。

テスト名だけでcoverageを主張せず、実際にFile move操作を含める。既存の`moveLinearSeqUp/Down`直接呼び出しは、source移動なら新しいUIまたはsource actionへ、同一row record移動なら`moveLinearRecordWithinRow()`の契約へ置き換える。

### 10.3 Regression scope

focused checksの後、全Web unit tests、関係するbrowser comparison contracts、architecture contract、Web change policyを実行する。Python productionを変更しないためPython全体の再実行は失敗や影響分析が示す場合に広げるが、実Worker/Python Generateを少なくとも一つの代表fixtureで確認する。

## 11. 実装セッション

| セッション | 主な成果 | 完了条件 |
| --- | --- | --- |
| S0: preflight | 最新base、回帰、authority、owner、テスト環境を再確認する。 | 実装可能範囲と変更対象が確定し、新しいProduct判断やarchitecture例外がない。 |
| S1: implementation | domain変換、app action、File header UI、旧path削除、unit/browser focused tests、文書を一つの整合した変更として実装する。 | SR-01〜SR-08、SR-12〜SR-14の該当focused evidenceが通る。productionとtest差分を別々に監査済み。 |
| S2: acceptance | Generate、Save/fresh Load、mobile/keyboard、広い回帰、policy gateを検証してin-scope不具合を修正する。 | SR-01〜SR-14が成立し、必須gateが通り、未検証事項がない。 |

セッションごとの具体的な開始指示は別紙にある。S1とS2を同じ担当者が連続実行してよい。セッション分割は、通常の編集、build、テストに追加承認が必要という意味ではない。

## 12. 検証コマンド

実装時の実ファイル名とtest名に合わせてgrepは調整する。長時間testは途中出力を監視し、リポジトリ既定のtimeoutを短縮しない。

```bash
node --test tests/web/linear-sources.test.mjs
node --test tests/web/session-request.test.mjs
npx playwright test tests/web/linear-multi-record.playwright.spec.js \
  --grep "source.*order|File.*order" --workers=1 --retries=0
node --test tests/web/architecture-contracts.test.mjs
node --test tests/web/*.test.mjs
python tools/prepare_browser_wheel.py
node tools/check-web-change-budget.mjs --base origin/dev
git diff --check
```

必要に応じて、既存のcomparison contract suiteとfast Python suiteを追加する。

```bash
npm run test:web:comparison-contracts
pytest tests/ -v -m "not slow"
```

browser verificationの前に、NodeとPythonのPlaywright経路を確認する。Node runnerがなければPython Playwrightで同じ利用者操作を検証する。Chromiumがsandbox制約で失敗した場合は、同じローカルcheckを必要な権限で再実行し、利用不能とは結論しない。

## 13. リスクと対策

| リスク | 対策 |
| --- | --- |
| File moveが一recordだけを動かす | groupを交換してから全recordsをflattenするunit contractを置く。 |
| Fileカード順とcanonical順がずれる | browserでDOM、`linearSeqs`、requestを同時assertする。 |
| rowやcomparisonをindexで誤って付け替える | UID対応をbefore/afterでassertし、既存reconcileだけを使う。 |
| raw cacheを無条件に破棄または誤再利用する | `preserveLosatCacheInfo`と既存semantic keyを使い、derived index更新を検証する。 |
| 同名別uploadが結合される | file identity/descriptorを使う既存groupingを再利用し、同名別source testを置く。 |
| async discoveryが旧indexへ展開する | UID再解決を維持し、upload直後moveのbrowser testを置く。 |
| UI追加でsidebarが横overflowする | 390px viewportでFile名、Remove、move controlsを確認する。 |
| 旧Sessionを黙って並べ替える | Loadだけでは正規化せず、明示move時だけblock化する。 |
| テストが再び内部関数だけを呼ぶ | pointer/keyboardからrequestまでを一つのbrowser journeyで検証する。 |

## 14. Rollback

変更はSession schemaやPython出力形式を変更しないため、rollbackは局所的である。

- Fileヘッダーcontrols、source move action、純粋helper、対応tests/docsを同じ変更単位で戻す。
- Session migration、resource変換、reference SVGのrollbackは不要である。
- rollback後は既知のFile順序変更不能が再発するため、暫定回避策を明示する。
- 不具合回避のために旧record-level input reorderだけを復活させない。source/record粒度不一致を再導入する。

## 15. Definition of Done

- SR-01〜SR-14の結果が実際のテストまたは明示した同等証拠で確認されている。
- File順序の唯一のstate ownerが`linearSeqs`である。
- source identity判定は`groupLinearSourceRecords()`へ集約されている。
- Input Fileカードにsource単位操作があり、旧record単位操作は競合しない。
- record placement、comparison endpoint、cache、shared resource、SessionをUIDで維持する。
- canonical request、Worker、Pythonに並行経路を追加していない。
- production、tests、docs、生成物の差分を別々に確認している。
- architectureのowner/path証拠、user-visible behavior、checks、rollbackを記録している。
- 受入条件やtimeoutを失敗に合わせて弱めていない。
- push、PR、merge、tag、deployは別途明示的に許可された場合だけ実行する。
- 完了時に英語のproposed commit titleと短いsummaryを提示する。

## 16. 実施記録テンプレート

実装セッション間の引継ぎでは、チャット履歴ではなく次を記録する。

| 項目 | 記録 |
| --- | --- |
| base SHA / branch | 未開始 |
| audited cause still present | 未開始 |
| changed production files | 未開始 |
| changed test files | 未開始 |
| changed docs/generated files | 未開始 |
| acceptance IDs completed | 未開始 |
| exact commands and results | 未開始 |
| browser/wheel source identity | 未開始 |
| architecture owner/path evidence | 未開始 |
| remaining work or blockers | 未開始 |
| rollback note | 未開始 |

### S0実施記録

| 項目 | 記録 |
| --- | --- |
| base SHA / branch | `origin/dev` = `eec715c70323c1dc6173a36bda2f18e600f8384d`; `fix/linear-source-reorder-20260920`。最新`origin/dev`から追跡先なしで作成し、計画文書コミットを適用した。 |
| audited cause still present | あり。Fileカードは`linearSourceGroups`単位だがヘッダーに移動操作がなく、`moveLinearSeqUp/Down()`は`linearSeqs`の一recordだけを移動する。canonical requestは`linearSeqs`順からrecordsを構築する。 |
| baseline commands and results | `node --test tests/web/linear-sources.test.mjs`: 1 pass。`npx playwright test tests/web/linear-multi-record.playwright.spec.js --grep "one uploaded source stays one file card" --workers=1 --retries=0`: Chromium 1 pass。 |
| coverage gap | unitはsource groupingとLOSAT batchingだけを検証し、browserはFileカード数、collapsed record list、置換、削除を検証するが、Fileヘッダーからsource順を変更しない。既存passはSR-01〜SR-14のreorder証拠ではない。 |
| Product Impact | `IMPLEMENT_EXISTING_AUTHORITY`。`PD-OI-018` scenario revision 2と`OIC-015`が、one source = one File cardと、reordering後のidentity、placement、pair mapping、shared resources保持を選択済み。別のProduct判断はない。 |
| owner/path plan | source identityとblock変換は`linear-sources.js`、coordinationは`app-setup.js`から既存`applyLinearSeqMutation()`、表示は`index.html`、request projectionは既存`session-request.js`のまま。旧record-level input reorder pathを削除する。 |
| environment | Node Playwright `1.61.0`、Python Playwright import、Node `@playwright/test`を確認。`gbdraw/web/gbdraw-0.14.0-py3-none-any.whl`が存在する。S2で現在sourceから再生成する。 |
| S1 readiness / blockers | S1開始可能。新しいProduct判断、architecture例外、外部権限は不要。 |

### S1実施記録

| 項目 | 記録 |
| --- | --- |
| base / branch | S0と同じ`eec715c70323c1dc6173a36bda2f18e600f8384d` / `fix/linear-source-reorder-20260920`。 |
| changed production files | `gbdraw/web/js/app/linear-sources.js`: pure source-block move。`gbdraw/web/js/app/app-setup.js`: source actionを既存mutation境界へ接続。`gbdraw/web/index.html`: File header controlsと旧record controls削除。 |
| changed test files | `tests/web/linear-sources.test.mjs`: block、identity、boundary、same-name、Session descriptor、legacy interleave。`tests/web/linear-multi-record.playwright.spec.js`: pointer/keyboard、2+3 records、request、row/pair/cache、Result、Session、mobile、async discovery、GFF3+FASTA。`tests/ci/playwright-inventory.test.mjs`: comparison contract総数を追加後の16本へ同期。 |
| changed docs / generated files | `docs/REFERENCE/web-app.md`: File順序とRecord Layoutの責務を説明。generated fileは未変更。 |
| focused acceptance | SR-01〜SR-08、SR-11〜SR-14をfocused testsで確認。SR-09/SR-10のGenerateとfresh Loadも同じbrowser journeyで確認したが、S2の広いgateとwheel再生成は未実施。 |
| exact commands and results | `node --test tests/web/linear-sources.test.mjs`: pass。`node --test tests/web/session-request.test.mjs`: pass。`npx playwright test tests/web/linear-multi-record.playwright.spec.js --grep "source.*order\|File.*order" --workers=1 --retries=0`: Chromium 3 pass。`git diff --check`: pass。 |
| removed owner/path | Input File record-level `canMoveLinearSeqUp/Down`、`reorderLinearSeqs`、`moveLinearSeqUp/Down`とRecord optionsのUp/Downを削除。recordの同一row移動は`moveLinearRecordWithinRow()`へ収束。 |
| architecture concise evidence draft | source identity/order semanticsは`linear-sources.js`一owner。UI -> `moveLinearSource()` -> `moveLinearSourceGroup()` -> `applyLinearSeqMutation()` -> `linearSeqs` ->既存request/Worker経路。新state、watcher、schema、compatibility pathなし。想定`OE 0 -> 0`、`PE 0 -> 0`、`CB 0 -> 0`。 |
| remaining S2 work | browser wheel再生成、実LOSAT/Worker代表journey、全Web unit、comparison contracts、architecture contract、Web policy gate、production/test/docs/generated diff最終監査。 |
| proposed commit title / summary | `Restore Linear source reordering from File cards` — Add accessible File-level moves that keep multi-record sources, record-owned state, comparisons, cache, and Session replay aligned while removing the obsolete record-level input reorder path. |

### S2実施記録

| 項目 | 記録 |
| --- | --- |
| browser / wheel identity | Node Playwright `1.61.0`、Python Playwright、Node `@playwright/test`を確認。`python tools/prepare_browser_wheel.py`で現在sourceから`gbdraw-0.14.0-py3-none-any.whl`を生成。SHA-256 `8b5cb22f096f7a7e646bbaec006f0ea709085f1aa25bcf4f8ac36d62806e2b55`、wheel metadata version `0.14.0`。 |
| browser acceptance | source-order focused journey 3 pass。実LOSAT/Wasm offline journey 1 pass。comparison contract suite 16 pass。旧record順序testを正しいRecord Layout/source actionへ移した2 journeyもpass。pointer、Space、Enter、disabled boundary、Generate、fresh Load、cache/index reconciliationを確認。 |
| mobile / visual | 390x844でFile名、Remove、Up/Down、collapsed countを検証し、各File cardで`scrollWidth <= clientWidth`。使い捨て`file-source-order-mobile.png`を目視し、横overflowまたは到達不能なし。 |
| unit / architecture | `node --test tests/web/*.test.mjs`をPython子プロセス利用のためsandbox外で実行し606 pass。`node --test tests/web/architecture-contracts.test.mjs`は137 pass。`node --test tests/ci/*.test.mjs`は59 pass。sandbox内一覧runは既存`feature-color-actions.test.mjs`のPython childで停止したが、同じsuiteのsandbox外再実行で解消。 |
| policy / deterministic checks | `node tools/check-web-change-budget.mjs --base origin/dev --head HEAD`: Gate PASS、Review REQUIRED（public export追加）。production 3 files、gross churn 108、net +2、cycle 0 -> 0、authority/dependency/guard deltaなし。`git diff --check`はpass。 |
| acceptance IDs | SR-01〜SR-14 pass。source block、same-name別upload、boundary no-op、collapsed records、record state、UID-row、explicit pair/cache、draft Result、Generate、Session、keyboard/mobile、旧path削除、async discovery、GFF3+FASTAを自動testまたは同等browser evidenceで確認。 |
| architecture concise evidence | Before: source groupingは`linear-sources.js`、旧input reorder semanticsは`app-setup.js`のrecord-level path。After: source identityとblock変換は`linear-sources.js`一owner、UI -> `moveLinearSource()` -> `moveLinearSourceGroup()` ->既存`applyLinearSeqMutation()` -> `linearSeqs` ->既存request/Worker path。旧record-level pathを削除。新state、watcher、schema、compatibility pathなし。通常non-increasing changeで`OE 0 -> 0`、`PE 0 -> 0`、`CB 0 -> 0`。 |
| rollback | File header controls、source action、pure helper、tests/docsを同じ単位で戻す。Session migration、resource変換、reference SVGのrollbackは不要。旧record-level input reorderだけは復活させない。 |
| remaining work | ローカル実装・受入は完了。push、PR、CI確認、問題がなければmerge。 |
