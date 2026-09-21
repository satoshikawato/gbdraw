# Linear入力・表示UI改善（Issues #559・#560・#562）— 総合計画書

状態: 計画策定済み、runtime未実装。Issue #559は既存Product authorityの範囲で着手可能。
Issue #560と#562は、本書第5節のProduct Decisionをdurable authorityへ反映したbaseから
runtime実装を開始する。

固定実装ブランチ: `fix/linear-ui-559-560-562-20260921`。本計画のS1〜S6と、各sessionの
ローカル実装commitはこのブランチへ積み上げる。別のruntime branchを作らない。S0の
authority-only変更だけは別branchで扱い、authorityが`origin/dev`へmergeされた後、この固定
実装ブランチを更新済み`origin/dev`へrebaseしてから依存runtimeを開始する。

監査基準: `origin/dev` / `f4760476915194fc0586642e0b03ec561b96ab01`
（2026-09-21 JST）。実装開始時には最新`origin/dev`との差分を確認し、このSHAを最新と
仮定しない。

対象Issue:

- [#559 “Depth TSV files (applied to all records)” is too big and incomplete](https://github.com/satoshikawato/gbdraw/issues/559)
- [#560 Improve visibility of Add/Remove buttons and introduce options for card removal](https://github.com/satoshikawato/gbdraw/issues/560)
- [#562 UI/UX Improvements: Default styles and Layout panel reorganization](https://github.com/satoshikawato/gbdraw/issues/562)

各実装セッションへそのまま渡せる開始指示は、別紙
[INSTRUCTION PROMPTS](LINEAR_UI_559_560_562_INSTRUCTION_PROMPTS_2026-09-21.md)にある。
本書は、過去の会話を知らない参加者が、問題、現行構造、Product判断の状態、設計、
依存関係、受入条件、検証結果を一か所から復元するための管理文書である。

## 1. 目的と完成像

gbdrawは、注釈付き配列からSVGなどのゲノム図を生成するPythonアプリケーションである。
Web版はビルド工程のないVue単一ページアプリで、canonical render requestを一度構築し、
Worker内のPyodideから共通Python描画経路を呼ぶ。

本計画はLinearモードの三つの操作領域を改善する。

1. FileカードのDepth TSV割当を、必要なときだけ展開でき、同じ場所から系列を追加できる
   ようにする。
2. 入力追加・削除操作を見つけやすくし、ファイル内容のクリアとFileカード自体の削除を
   明示的に区別する。
3. 新規Linear図の比較表示、record labelの可視性、タイトル・label・legendの設定構造を
   一貫させる。

完成後も、Web state → canonical request → Worker → Python rendererという経路は一つだけで
ある。UI再編のために描画経路、Depthデータモデル、比較アルゴリズム、CLI/APIの既定値を
増やさない。

## 2. 用語

| 用語 | 本書での意味 |
| --- | --- |
| File / source | 一回のGenBank upload、または一組のGFF3＋FASTA。Input Genomesに一枚のFileカードとして表示される。 |
| record | Fileに含まれる一本の配列。File一つに複数recordが存在し得る。 |
| source group | 同じsourceを参照するrecordsをFileカードへまとめた導出値。別の永続stateではない。 |
| blank slot | primary sequence fileを持たない一枚の入力枠。新しいsourceを同じFile順序へ挿入できる。 |
| primary source | GenBank、またはGFF3＋FASTA。Depth TSVはprimary sourceではない。 |
| Depth series | record-major Depth matrixの同じ論理列。各File／recordのcellはファイルまたは空値を持つ。 |
| selected mode | AccessionまたはLength / Coordinatesで利用者が選んだ`auto`、`show`、`hide`。 |
| effective visibility | selected modeと実際のrendered row配置から導出した最終boolean。 |
| rendered row | Generate時にrecordsを配置するLinearの行。upload file数や総record数とは別の概念。 |
| current Result | 最後に成功したGenerateの不変成果物。draft編集だけでは置換しない。 |
| Session | Web状態、resources、canonical request、保存済みResultなどを再利用する保存形式。 |

Issue #560は「record」「file card」という語を混在させている。現行WebではFileカードが
source groupを表し、multi-record Fileを一枚のカードにまとめる。本計画の削除契約は
record一件ではなく、このFile/source groupを単位とする。

## 3. 適用する規則と開始条件

実装者は着手時に次を全て読む。

- [AGENTS.md](../../AGENTS.md)、[CLAUDE.md](../../CLAUDE.md)、
  [Web CLAUDE.md](../../gbdraw/web/CLAUDE.md)
- [Product Impact Ratchet](PRODUCT_IMPACT_RATCHET.md)、
  [Product Decision Packet Template](PRODUCT_DECISION_PACKET_TEMPLATE.md)、
  [Option Integrity Product Contract](OPTION_INTEGRITY_PRODUCT_CONTRACT.md)
- [Architecture Fitness Function Ratchet](ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)、
  [Web Change Policy](WEB_CHANGE_POLICY.md)
- [Session compatibility](../SESSION_COMPATIBILITY.md)、
  [Web reference](../REFERENCE/web-app.md)

計画策定sessionで、最新`origin/dev`の
`f4760476915194fc0586642e0b03ec561b96ab01`から追跡先なしの固定実装ブランチ
`fix/linear-ui-559-560-562-20260921`を作成した。S1〜S6は毎回このbranch名、HEAD、upstream、
base ancestryを確認し、同じbranchで継続する。別のruntime branchを作らず、無関係な変更、
未追跡ファイル、別worktreeを保持する。継続セッションは本書の実施記録とGit差分から対象作業を
復元する。

S0でauthority-only変更が必要な場合は、repository policyに従って最新`origin/dev`から別の
authority branchを作り、runtimeを含めずに扱う。authority merge後は固定実装ブランチへ戻り、
`git fetch origin`でmerge済みbaseを確認してから`git rebase origin/dev`する。rebase前に現在の
branch、未commit差分、対象authority SHAを確認し、候補authorityのcommitを固定実装ブランチへ
直接cherry-pickしてruntimeを自己承認しない。

本計画と別紙promptは、push、PR作成、merge、tag、release、deployを許可しない。
それらは実行時点の明示的な許可に従う。

## 4. 現行実装と原因監査

### 4.1 #559: File-level Depth領域

`gbdraw/web/index.html`の各Linear Fileカードは、`Depth TSV files (applied to all records)`を
常時展開した`div`として描画する。論理seriesごとのCommon / Mixed / Empty状態とuploaderは
既に正しくFileカード上へ公開されているが、File数やseries数が増えるほどsidebarの高さを
占有する。

別のDepth設定領域には`addLinearDepthTrack()`を呼ぶ`Add series for all records`がある。
したがって、欠けているのは新しいDepth機能ではなく、File-level割当の近くにある発見可能な
入口である。

既存の正準表現と所有者:

- record-major matrixとlogical column操作: `web/js/app/depth-track-state.js`
- File group、common/mixed/empty導出: `web/js/app/linear-sources.js`と`app-setup.js`
- add/remove、File一括apply/clearのcoordination: `web/js/app/app-setup.js`
- request／Session投影: `web/js/services/session-request.js`

[PD-OI-025](OPTION_INTEGRITY_PRODUCT_CONTRACT.md#pd-oi-025-linear-depth-source-scope-and-discoverability)
と`OIC-020`は、Fileカードを開くだけでcommon Depth割当に到達でき、record listやrecord
optionsを開く必要がないこと、record-major matrixを維持することを要求する。

### 4.2 #560: Removeの意味がUI階層と一致しない

`FileUploader`の赤いRemoveは、現在は即座に`update:modelValue(null)`をemitし、独自に
History transactionを開始する。Linear primary uploaderでは、その値が
`setLinearSeqPrimaryFile()`へ渡り、最後のprimary fileが消えるとsource group全体も消える。
UIの見た目は「ファイルを外す」だが、到達するdomain結果は「Fileカードを削除する」になり
得る。

sidebar下端の`+ Add sequence`と`- Remove`は、Fileカード内の操作より視覚的に弱い。
`removeLastLinearSeq()`は最後のsource groupを確認なしで削除する。multi-record Fileの場合も
group全体が対象である。

一方、source変更後に必要な処理は既に`applyLinearSeqMutation()`へ集約されている。
この境界は、UID-row reconciliation、comparison plan、derived artifacts、raw cache情報、
pending discoveryを整合させる。新しい削除経路がこれらを個別に再実装してはいけない。

未決なのはコードの置き場所ではなく、Clear後にどのsource/record固有stateを残すか、
GFF3＋FASTAを一単位として扱うか、最後の一枚をDeleteできるかというProduct結果である。

### 4.3 #562: fresh default、UI grouping、Auto policy

現行のfresh Web stateは次を持つ。

- `adv.pairwise_match_style = 'ribbon'`
- `adv.linear_show_accession = true`
- `adv.linear_show_length = true`

AccessionとLengthはbooleanなので、利用者が明示的にShowを選んだ状態と、layoutに応じて
自動表示している状態を区別できない。checkboxをwatcherで書き換える実装では、layout変更が
利用者の明示選択を破壊する。

`session-request.js`にはeffective rowを解決し、shared rowを判定する処理が既にあるが、
canonical request構築の内部に置かれている。UI summaryとrequestが同じAuto結果を使うには、
row semanticsの既存ownerである`linear-record-layout.js`から再利用できる純粋な解決境界が必要
である。

タイトル、Definition Font Size、visibility checkbox、Definition Line Styles、legend controlsは
同じ大きなsectionに連続している。style state自体は
`linear_definition_line_styles`に既に一元化されており、問題はstate不足ではなく情報設計である。

Session current writerはversion 42、canonical render requestはschema 7である。現在のrequestと
PythonはAccession／Lengthの最終booleanだけを必要とする。選択した`auto/show/hide`を保存する
にはSession stateの更新が必要だが、request schema、Worker protocol、Python描画入力を変更する
必要はない。

### 4.4 変更してはいけない境界

三Issueとも、次は原因ではない。

- Depth TSVの数値解析、axis、track rendering
- comparison計算、filtering、LOSATの実行方法
- Python rendererのDefinition行生成
- CLI/APIの`pairwise_match_style`既定値
- Circular modeのlabel visibility
- canonical request schema 7とWorker protocol

不具合を直すためにこれらを変更しない。到達する共通testの回帰修正が必要な場合だけ、根拠を
実施記録へ残す。

## 5. Product Impactと決定経路

### 5.1 #559は既存authorityで実装する

推奨結果は次のとおりである。

- File-level Depth領域をnative disclosure（`details` / `summary`）にする。
- 初期状態はopenとし、PD-OI-025が保証する既存の発見可能性を維持する。利用者は任意に閉じられる。
- disclosureの開閉は一時的なUI状態とし、SessionやHistoryへ保存しない。
- summaryは既存matrixからseries数とCommon / Mixed / Emptyの状態を導出し、新しいstateを持たない。
- `Add Depth TSV series`をseries一覧の直下へ置き、既存`addLinearDepthTrack()`を呼ぶ。
- 同じglobal logical-series追加buttonが別のDepth設定領域にある場合は、入口をFile-level文脈へ
  集約する。新旧二つの意味の異なるadd actionを残さない。

この結果はPD-OI-025とOIC-020の発見可能性を保ち、data semanticsを変更しないため、
`IMPLEMENT_EXISTING_AUTHORITY`として扱う。初期状態をclosedへ変更する、File-level割当を
Record options内へ隠すなど、既存の到達性を弱める案へ変更する場合は、実装を止めてProduct
preflightをやり直す。

### 5.2 #560はProduct Decisionが必要

以下はDecision Packに提示するstable outcomesである。選択済みauthorityではない。

| ID | 利用者から見える結果 | 主なtrade-off |
| --- | --- | --- |
| **D560-A: Pristine placeholder（推奨）** | `Clear file only`はsource groupを同じFile順序の一つのblank slotへ置換する。GenBank、またはGFF3＋FASTAを一つのprimary source単位としてclearし、旧record UID、selector、crop、reverse、label/subtitle、Depth binding、comparison endpoint、source-derived cacheを引き継がない。`Delete card`はgroupを削除する。最後の一枚ではDeleteを無効にし、Clearだけを許可する。 | 新sourceへの隠れたstate混入を防ぐ。旧source固有設定は再入力が必要。 |
| **D560-B: Configured placeholder** | primary source fileだけをclearし、File順序に加えて、適用可能なslot設定を保持する。paired GFF3＋FASTAの片側clearも許容し得る。 | 置換は速いが、どの設定を保持するかが複雑で、旧source stateが新sourceへ誤適用される危険がある。 |
| **D560-C: Confirmation only** | card-level Clear/Deleteの分離を追加せず、現行削除へ確認だけを追加する。 | 最小変更だが、同じ位置へ差し替えるIssueの主目的を満たさない。 |

D560-Aを選ぶ場合のglobal操作案:

- `+ Add sequence`と`- Remove`をFileカード操作と同等の視認性・hit areaにする。
- 最後のFileカードがpristine blank slotなら、二枚以上ある場合に限り即時Deleteする。
- primary sourceがある場合は、対象File名とrecord数を示すconfirmationを開き、Confirmで
  `Delete card`、Cancelで完全なno-opとする。
- Fileカード内の赤いRemoveは、非空sourceに対して`Clear file only`、`Delete card`、Cancelの
  三つを提示する。既にblankならchoice dialogを出さない。
- mutationは一つのHistory operationとし、current Resultは次のGenerate成功まで保持する。

実装前にProduct Decision Ownerが、選択、rationale、must preserve、may retire、accepted
residual risk、owner、dateを明示する。D560-Aの括弧内詳細も、human receiptなしに採用済みと
扱わない。

### 5.3 #562は三つの独立したProduct Decisionに分ける

Issue #562は具体的な望ましい結果を記載しているが、base branchのdurable authorityにはまだ
登録されていない。また、Curve既定値、sidebarの情報設計、Auto visibilityは、利用者への効果、
互換性、実装owner、rollback条件が独立している。一つのbundled outcomeへ束ねず、次の三concernsを
別々に判断・serializeする。

#### 5.3.1 Match-style default

Concern key: `linear.comparison-match-style-default`

| Stable outcome ID | 利用者から見える結果 | 主なtrade-off |
| --- | --- | --- |
| **CURVE-FRESH-RESET（推奨）** | fresh LinearとReset SettingsはCurveを使う。明示選択、既存Session、Gallery、Circular、CLI/APIには遡及適用しない。 | 新規利用者の初期表示が変わるが、Session schema変更は不要。 |
| **KEEP-RIBBON** | fresh/resetを含め、現在のRibbon既定を維持する。 | 互換作業は不要だが、Issueが求める初期Curveを実現しない。 |

`CURVE-FRESH-RESET`のproposal:

1. fresh Linear stateとReset SettingsだけがCurveを使う。record、layout、comparison設定変更は
   明示選択を上書きしない。
2. 旧Sessionの明示`ribbon`／`curve`を保持する。旧データでfieldが欠ける場合はhistorical
   fallbackのRibbonとし、fresh defaultのCurveを遡及適用しない。
3. comparison計算、filter、geometry、Circular、CLI/APIの既定値を変更しない。

#### 5.3.2 Presentation-control organization

Concern key: `linear.presentation-control-organization`

| Stable outcome ID | 利用者から見える結果 | 主なtrade-off |
| --- | --- | --- |
| **REGROUP-TITLES-LABELS-LEGEND（推奨）** | Plot Title、Record Labels、Legendを明確な階層へ再編し、各labelのvisibilityとStyleを近接させる。 | DOM構造とbrowser受入範囲は変わるが、state／schema／rendererは変えない。 |
| **KEEP-CURRENT-ORGANIZATION** | 現在のcontrol配置と表示名を維持する。 | 実装変更はないが、設定の発見性と概念分離を改善しない。 |

`REGROUP-TITLES-LABELS-LEGEND`のproposal:

1. Plot Titleはdiagram全体の設定としてRecord Labelsと分ける。
2. Record Labelsは既存style stateを使い、Default font size、各labelのvisibility、各Style
   disclosureを配置する。visibility mode自体はこのdecisionで変更しない。
3. Legendは同階層の独立sectionとし、collapsed summaryを`Legend · Bottom`などとする。
4. `Definition Font Size`を`Default font size`、`Legend Box Size`を`Swatch size`と表示するが、
   保存key、単位、rendererを変更しない。
5. controlの移動でduplicate state、duplicate control、Session fieldを追加しない。

#### 5.3.3 Independent Auto visibility

Concern key: `linear.record-label-auto-visibility`

| Stable outcome ID | 利用者から見える結果 | 主なtrade-off |
| --- | --- | --- |
| **INDEPENDENT-AUTO-SHOW-HIDE（推奨）** | AccessionとLength / Coordinatesが独立したAuto/Show/Hideを持ち、any shared rendered rowではAutoをdiagram-wideでHideする。 | selected mode保存のためSession schema migrationと互換readerが必要。 |
| **KEEP-BOOLEAN-VISIBILITY** | 現在の二つのboolean checkboxと常時明示Show/Hideを維持する。 | Session変更は不要だが、multi-record rowの自動clutter制御を実現しない。 |

`INDEPENDENT-AUTO-SHOW-HIDE`のproposal:

1. AccessionとLength / Coordinatesは各々`auto`、`show`、`hide`を保存する。fresh defaultは各`auto`。
2. Autoはeffective rendered rowsを使う。Record Layoutが無効なら、dormant row設定ではなく
   実際の一record／row配置として解決する。
3. 一つでも二records以上のrendered rowがあれば、そのfieldを全rowsでHideする。全rowsが
   single-recordへ戻ればShowする。AccessionとLengthは互いに独立する。
4. explicit Show／Hideをlayout変更で書き換えない。UIは`Auto · Shown`または
   `Auto · Hidden`を表示し、manual choiceをdisableしない。
5. legacy booleanは`true → show`、`false → hide`。欠落値はversion 42以前のhistorical default
   であるShowへ移す。明示booleanをAutoへ変えない。
6. selected modeはSessionへ保存し、effective booleanだけをcanonical requestへ投影する。
   request schema 7、Worker、Python、CLI/APIは変更しない。

三decisionsは独立して採否、実装、rollbackできる。例えばCurveとUI再編を採用し、Session 43を
伴うAutoだけを延期できる。選択IDを共有していても互いのauthorityを代用しない。

### 5.4 authority-only workflow

S0では最新baseのProduct Impact mapと既存authorityを検索し、同じconcernを重複登録しない。
該当authorityがなければ、概ね次のconcern keyでDecision Packを作る。

- `linear.input-source-removal`
- `linear.comparison-match-style-default`
- `linear.presentation-control-organization`
- `linear.record-label-auto-visibility`

Product Decision Ownerの回答は次の形式で受け取る。

```text
PRODUCT_DECISION
Concern: <concern key>
Scenario revision: <revision>
Choice: <stable outcome ID>
Rationale: <product-level reason>
Must preserve: <effects and affordances>
May retire: <none or explicit scope>
Accepted residual risk: <bounded risk or none>
Owner: <maintainer identity>
Decision date: <YYYY-MM-DD>
```

開発者は不足項目を推測しない。明示回答後は回答だけを既存のdurable authorityへ転記する
authority-only変更を準備し、runtimeは含めない。そのauthorityがmergeされたbaseからruntime
branchを作る。candidate authorityは同じcandidate runtimeを自己承認しない。

## 6. 目標アーキテクチャ

```text
index.html
  ├─ File Depth disclosure / Add series affordance
  ├─ Add / Remove buttons / removal dialog
  └─ Titles & Record Labels / Legend
       ↓ event and derived presentation only

app-setup.js
  ├─ existing Depth actions
  ├─ one removal-intent coordinator + one History transaction
  └─ derived effective label visibility exposed to UI
       ↓ pure domain operations

linear-sources.js                 linear-record-layout.js
  source grouping / clear-delete   effective rendered rows / shared-row fact
       └──────────────┬───────────────────────┘
                      ↓
linear-label-visibility.js (one pure policy owner, only if INDEPENDENT-AUTO-SHOW-HIDE is selected)
  selected mode validation + effective boolean resolution
                      ↓
session-request.js
  current selected state → schema 7 effective booleans
                      ↓
existing Worker → existing Python renderer

config.js / session-active-config-contract.js / session_io.py
  Session current writer + bounded legacy migration only
```

`linear-label-visibility.js`は、`INDEPENDENT-AUTO-SHOW-HIDE`で新しいpolicy semanticsが必要な
場合だけ追加する。
単にコードを分割する目的では追加しない。rowの解決自体は`linear-record-layout.js`が所有し、
visibility moduleはrow入力からpolicyを一度だけ適用する。

### 6.1 所有者表

| 責務 | Owner | 実装境界 |
| --- | --- | --- |
| Depth matrix／logical columns | `depth-track-state.js` | record-major表現、add/remove column。UI disclosureを知らない。 |
| File grouping／source identity | `linear-sources.js` | source groupとclear/deleteのpure transform。DOM、dialog、Historyを知らない。 |
| Linear row semantics | `linear-record-layout.js` | effective rowsとshared-row判定。Session/request/UI policyを知らない。 |
| Label visibility policy | `linear-label-visibility.js`（`INDEPENDENT-AUTO-SHOW-HIDE`時） | modeのvalidation、legacy boolean mapping、effective boolean。VueやSession I/Oを知らない。 |
| UI coordination | `app-setup.js` | intent、dialog state、History、既存mutation pipelineを結ぶ。source規則を複製しない。 |
| Reusable uploader | `components.js::FileUploader` | 通常clearを維持し、Linear sourceだけがclear intentを外へ委譲できる最小interfaceを持つ。source削除規則を知らない。 |
| Fresh Web defaults | `session-active-config-contract.js` | fresh/reset current stateを一度定義する。legacy fallbackをfresh defaultへ混ぜない。 |
| Persisted load/save | `config.js`、`session_io.py` | current writerと実在する旧version migration。render policyを複製しない。 |
| Request projection | `session-request.js` | selected modeから得たeffective booleanを既存config overrideへ投影する。row/policy式を再実装しない。 |
| UI | `index.html` | hierarchy、label、summary、event forwarding。canonical stateを複製しない。 |

### 6.2 #559の実装設計

- 現在のFile-level Depth wrapperを`details`へ変え、`summary`に名称、record数、series数、
  mixedの有無を表示する。summary値は`linearSourceDepthRows(source)`から導出する。
- 初期`open`は明示し、利用者のtoggle後にreactive rerenderで勝手に開き直さない。開閉stateを
  adv、form、Sessionへ追加しない。
- series一覧直下のbuttonは既存`addLinearDepthTrack()`へforwardする。新しいDepth add関数、
  second matrix、File-local series countを作らない。
- 既存のglobal logical-series add入口を残す必要がない場合は削除する。Fileカードごとに同じ
  buttonが見えることはcontextual affordanceであり、どのbuttonも同じ一つのactionを呼ぶ。
- File-level apply/clear、Mixed warning、per-record sparse override、single-record時の重複uploader
  抑止を変更しない。

### 6.3 #560の実装設計

D560-Aが選択された場合の設計である。別outcomeが選ばれた場合は、本節をそのhuman receiptに
合わせて最小限修正してからruntimeへ進む。

1. `FileUploader`へ、clear buttonの通常動作を既定のまま保つopt-in interfaceを加える。
   Linear primary uploaderだけがclearを即時emitせず`clear-request`を通知する。Depth、color table、
   Circularなど既存利用者の動作を変えない。
2. `app-setup.js`はsource UID／group identityを使ってintentを保持する。array indexだけをdialogの
   長寿命identityにしない。dialogを開いている間に対象が消えた場合は安全なno-opにする。
3. dialogはClear、Delete、Cancelを明示し、対象File名とmulti-record countを表示する。
   keyboard focus、Escape、focus returnを持つ。`window.confirm()`で三択を模倣しない。
4. Clear/Delete transformはsource grouping ownerのpure operationとし、最終適用は
   `applyLinearSeqMutation()`を一度だけ呼ぶ。Historyをuploaderとcoordinatorで二重に作らない。
5. Clearは新しいblank sequence UIDを持つpristine slotへ置換する。旧UIDを流用してsource-bound
   stateを残さない。File順序だけを維持し、row/comparison/cacheは既存reconciliation結果を使う。
6. Deleteはgroupの全recordsを削除する。ただし全体を0 slotsにしない。sole FileではDeleteを
   disableし、Clearに収束させる。
7. global Removeは同じintent/controllerを使う。空判定と実行で別の規則を作らない。
8. draft mutationだけではcurrent Result、selected saved Result、SVGを消さない。次のGenerate成功で
   新しいdraftを反映する。

### 6.4 #562: Match-style default

`CURVE-FRESH-RESET`が選択された場合、fresh LinearとReset SettingsのWeb既定をCurveにする。
`KEEP-RIBBON`では本節のruntime変更を行わない。`createDefaultAdv(mode)`はfresh current stateだけを
定義し、legacy loadの欠落値をCurveで埋めない。

match-style enum validation／normalizationが`config.js`とrun analysis側に重複している場合は、
`current-option-values.js`の一つのvalidatorへ収束させる。callerは用途に応じて明示fallbackを渡す。

- fresh/reset Linear fallback: `curve`
- older persisted stateのmissing field: `ribbon`
- explicit valid saved value: その値
- invalid current value: 既存validation方針に従って拒否または既定化し、silentな別規則を増やさない

この変更はmatch-style defaultだけを扱う。Record Labels／LegendのDOM移動、visibility state、
Session versionを同じcommitへ含めない。

### 6.5 #562: Presentation-control organization

`REGROUP-TITLES-LABELS-LEGEND`が選択された場合、UIは同じstateを次の階層へ移す。
`KEEP-CURRENT-ORGANIZATION`では本節のruntime変更を行わない。

```text
Titles & Record Labels
  Plot Title
    Text / Position / Font size
  Record Labels
    Default font size
    Name / Species                 [Style disclosure]
    Subtitle                       [Style disclosure]
    Replicon        [Show/Hide]    [Style disclosure]
    Accession       [existing bool or independently selected tri-state] [Style]
    Length / Coordinates [existing bool or independently selected tri-state] [Style]

Legend · <current position>
  Position / Font size / Swatch size
```

UI再編ではstyle objectやvisibility keyの二つ目を作らない。Plot TitleとRecord Labelsの概念を
分け、Legendを同階層の独立cardへ移す。individual Style disclosureは既存size／weight／colorを
直接編集する。collapsed/open UI状態は保存しない。Curve／RibbonのdefaultやAuto visibility
semanticsはこの変更で扱わない。

### 6.6 #562: Independent Auto policyとSession 43

`INDEPENDENT-AUTO-SHOW-HIDE`が選択された場合だけ実施する。
`KEEP-BOOLEAN-VISIBILITY`では本節のruntime変更を行わない。UI再編の採否には依存せず、再編済み
ならRecord Labels内、未再編なら既存visibility位置で三状態controlを公開する。

current editable stateを次のselected-mode fieldsへ収束させる。

```text
adv.linear_accession_visibility = 'auto' | 'show' | 'hide'
adv.linear_length_visibility    = 'auto' | 'show' | 'hide'
```

名称は実装前のfield inventoryで衝突がないことを確認する。同等の既存current fieldが最新baseに
追加されていればそれを使い、並行keyを作らない。current writerから旧editable booleansを除き、
次のpure resolutionを一箇所で使う。

```text
hasSharedRenderedRow = layoutEnabled
  && effectiveRows contains the same positive row more than once

resolve(mode, hasSharedRenderedRow):
  show -> true
  hide -> false
  auto -> !hasSharedRenderedRow
```

UI summaryとcanonical requestは同じresolved resultを使う。requestには既存の
`linear_definition_show_accession`／`linear_definition_show_length`相当のbooleanだけを出す。

selected modesを保存するため、current Session writerを43へ上げる。WebとPythonのversion定数、
supported versions、settings-only paths、documentation、fixtures、gallery/session testsを一つの
変更として整合させる。version 42以前のreaderは実在fieldについて次だけを行う。

- own property `linear_show_accession: true/false` → `show/hide`
- own property `linear_show_length: true/false` → `show/hide`
- own propertyがない → historical `show`
- current selected-mode fieldが既にある場合はlegacy booleanで上書きしない

legacy booleanをcurrent writerへ同時に書くdual-writeは行わない。保存済みResultと保存済み
canonical requestはeffective booleanを含むため、Loadしただけで既存previewを再生成・置換しない。
Generate時にmigrated selected modeから同じ外観を再現する。

Session 42 readerの追加はcompatibility burdenを増やすため、Architecture Ratchetのexception条件を
満たす可能性が高い。実装前後の正確なOE／PE／CB setsを測定し、v42がsupported setから正式に
retireされたときにlegacy bool promotionを削除するremoval conditionを記録する。architecture
checkerやsupported-version contractを弱めて通さない。

## 7. SOLID・KISS・DRY・YAGNIの適用

### SOLID

- Single Responsibility: source identity、row semantics、visibility policy、Session migration、UI
  presentation、coordinationを別の既存ownerまたは一つの新policy ownerへ分ける。
- Open/Closed: FileUploaderはopt-in eventで拡張し、既存利用者のclear behaviorを変更しない。
- Liskov Substitution: 新しいclass hierarchyは作らないため対象外。原則名を抽象化追加の理由に
  しない。
- Interface Segregation: pure resolverはselected modeとshared-row factだけ、source transformは
  source groupとintentだけを受ける。
- Dependency Inversion: DOMはbusiness ruleを持たず、app coordinatorはpure ownerと既存canonical
  mutation/request境界へ依存する。

### KISS

- disclosureはnative details、visibilityは3-value enum＋1 pure resolver、removalはClear/Delete/
  Cancelだけとする。
- drag-and-drop、toast framework、汎用modal framework、policy engine、layout watcherを追加しない。
- UI summaryは保存せず、その場で導出する。

### DRY

- Depth add/apply/clearは既存actionとmatrix ownerを使う。
- source clear/deleteとglobal removeは同じintent、empty判定、mutation pipelineを使う。
- effective rowsとshared-row判定は一つのrow ownerを使う。
- selected visibilityからbooleanへの変換は一つのpolicy ownerを使い、UIとrequestで式を複製しない。
- fresh Curve defaultとlegacy Ribbon fallbackを同じ曖昧なfallback expressionへ混ぜない。

### YAGNI

- disclosure state、File order、removal audit log、new Depth state、preview-only visibility stateを
  Sessionへ追加しない。
- request schema、Worker protocol、Python model、CLI optionを変更しない。
- AutoをName、Subtitle、Replicon、Circularへ一般化しない。
- 将来用の四つ目のvisibility mode、任意並べ替え、bulk source managerを追加しない。

## 8. 実装セッションと依存関係

| Session | 依存 | 成果 | runtime変更 |
| --- | --- | --- | --- |
| S0: Product preflight / authority | 最新base、Issue本文 | #560と#562三concernsのDecision Pack、明示回答の正確なserialization、authority-only変更 | なし |
| S1: #559 Depth disclosure | PD-OI-025を含むbase | collapsible File Depth、direct Add series、focused tests | あり。S0を待たず実施可能 |
| S2: #560 removal workflow | 選択済みauthorityがmergeされたbase | prominent actions、choice/confirmation、atomic clear/delete | あり |
| S3: #562 match-style default | 選択済みmatch-style authorityがmergeされたbase | fresh/reset Curveとhistorical Ribbon fallback | 選択が`CURVE-FRESH-RESET`の場合 |
| S4: #562 UI organization | 選択済みorganization authorityがmergeされたbase | Titles/Record Labels/Legend再編 | 選択が`REGROUP-TITLES-LABELS-LEGEND`の場合 |
| S5: #562 Auto/Session | 選択済みAuto authorityがmergeされたbase | pure resolver、tri-state、Session 43 migration | 選択が`INDEPENDENT-AUTO-SHOW-HIDE`の場合 |
| S6: integrated acceptance | S1–S5の該当成果 | browser、Session/gallery、architecture/Product gates、docs、最終diff | 必要なin-scope修正のみ |

S0とS1は独立に進められる。S2〜S5は、対応する候補authorityを含む同じ未merge branchから開始
してはならない。S3、S4、S5はProduct上もruntime上も独立しており、非採用のsessionを省略できる。
実行順はGit依存を最小化するためS3→S4→S5を推奨するが、一つの選択が他のauthorityを代用しない。

## 9. 受入条件

### 9.1 #559 — Depth disclosure

| ID | 条件 |
| --- | --- |
| DPT-01 | 各Linear FileカードのFile-level Depth領域をpointer、Enter、Spaceで開閉できる。初期openで、閉じるとsidebar高さが実際に減る。 |
| DPT-02 | summaryからDepthであること、series数、mixedの有無、対象record数を識別できる。値はmatrixから導出される。 |
| DPT-03 | record listとRecord optionsを閉じたまま、File-level Depth uploaderとAdd Depth TSV seriesへ到達できる。 |
| DPT-04 | Addを一回押すとlogical columnが一つ増え、全recordsのrowが同じ幅へpadされる。既存cellはshift／消失しない。 |
| DPT-05 | multi-record Fileのcommon apply、per-record replacementによるMixed、File-level clear、UndoでOIC-020を満たす。 |
| DPT-06 | same-name別Files、null cells、later series、Save/fresh Load、regeneration、canonical requestを維持する。 |
| DPT-07 | disclosure open state用のSession field、Depth default field、request field、render pathを追加しない。 |
| DPT-08 | 390px幅でsummary、status、uploader、Add buttonが横にはみ出さず、accessible nameが一意である。 |

### 9.2 #560 — input removal

| ID | 条件 |
| --- | --- |
| RM-01 | global Add/RemoveはFileカード操作と同等に見つけやすく、少なくとも既存button size policyに沿うhit area、visible focus、disabled stateを持つ。 |
| RM-02 | 非空Fileの赤いRemoveは選択済み契約どおりClear/Delete/Cancelを提示し、クリックだけで即時削除しない。 |
| RM-03 | global Removeはpristine blank last slotでは選択済み条件で即時処理し、非空last Fileでは対象を示すconfirmationを出す。 |
| RM-04 | Cancel、Escape、dialog外の禁止された操作はstate、History、cache、current Resultを一切変更しない。 |
| RM-05 | Clearは選択されたretention契約を正確に満たし、同じFile ordinalへ新sourceをuploadできる。 |
| RM-06 | Deleteはmulti-record source group全体を一度で削除し、sole-slot invariantを破らない。 |
| RM-07 | GFF3＋FASTAのclear単位、片側欠落、pending discoveryは選択された契約と一致し、orphan resourceを作らない。 |
| RM-08 | 成功操作はHistory一件で、Undo/Redoがrecords、rows、comparisons、resource bindingを一緒に戻す。nested Historyを作らない。 |
| RM-09 | explicit comparison endpoint、Adjacent plan、derived artifacts、raw cache metadataは既存reconciliation規則に従う。stale endpointを残さない。 |
| RM-10 | draft Clear/Deleteだけではlast successful Resultを変えず、次のGenerate成功で新しい入力を反映する。failed/canceled/stale GenerateもResultを保つ。 |
| RM-11 | Save/fresh Load後も選択結果が再現され、removal dialog stateや別File-order stateを保存しない。 |
| RM-12 | Linear source以外のFileUploader clear behaviorに回帰がなく、pointer／keyboard、focus return、390px幅を満たす。 |

### 9.3 #562 — Match-style default

| ID | 条件 |
| --- | --- |
| MSD-01 | `CURVE-FRESH-RESET`選択時、fresh Linear diagramとReset SettingsはCurve。Circular、CLI、Python APIの既定値は変わらない。 |
| MSD-02 | 明示的にRibbon／Curveを選んだ後のrecord、layout、comparison変更は値を上書きしない。 |
| MSD-03 | version 42以前のSession、settings-only config、Galleryの保存値とmissing-field historical Ribbonを保持する。 |
| MSD-04 | Session schema、request schema、comparison計算／filter／geometryを変更しない。 |

### 9.4 #562 — Presentation-control organization

| ID | 条件 |
| --- | --- |
| ORG-01 | `REGROUP-TITLES-LABELS-LEGEND`選択時、Plot Titleは独立subsection、Record Labelsはvisibilityと既存style、Legendは独立cardとして表示される。 |
| ORG-02 | `Definition Font Size`の表示名は`Default font size`、legend markerの表示名は`Swatch size`だが、保存key、単位、rendererは変わらない。 |
| ORG-03 | 各Style disclosureは既存size／weight／colorを直接編集し、同じsettingの重複controlや重複stateがない。 |
| ORG-04 | Legend collapsed summaryは現在positionを正しく表示し、変更後すぐ更新される。summary用stateを保存しない。 |
| ORG-05 | section hierarchy、label、help、keyboard、focus、390px幅が明瞭で、sidebarを広げない。 |
| ORG-06 | match-style default、visibility semantics、Session schema、rendererを変更しない。 |

### 9.5 #562 — Independent Auto and compatibility

| ID | 条件 |
| --- | --- |
| AV-01 | fresh LinearでAccessionとLengthのselected modeは各Auto、effective resultは各`Auto · Shown`である。 |
| AV-02 | 複数Files／複数recordsでも各rendered rowがsingleならAutoはShown。upload数や総record数だけではHideしない。 |
| AV-03 | 一つでもshared rendered rowがあれば、Autoのfieldはsingle-record rowsを含むdiagram全体でHiddenになる。 |
| AV-04 | Accession=Show、Length=Autoなど独立した組合せが正しく解決される。manual Show/Hideはlayout変更で書き換わらない。 |
| AV-05 | shared rowを解消するとAutoはShownへ戻る。selected mode自体はAutoのままである。 |
| AV-06 | Record Layout OFFではdormant row mapでなく実際のrendered arrangementを使う。UI、request、SVGが一致する。 |
| AV-07 | Undo/Redoはselected modesとlayoutを戻し、その都度effective resultを再計算する。derived booleanを独立復元しない。 |
| AV-08 | v42 explicit `true/false`はShow/Hide、missingはhistorical Showへ移る。明示値をAutoへ変えない。 |
| AV-09 | v43 Save/fresh Load/regenerationはselected modesを保持し、保存previewと再生成結果を正しく区別する。 |
| AV-10 | canonical request schema 7はeffective booleanのみを保持し、Worker/Python/CLI/API/Circular behaviorを変更しない。 |
| AV-11 | current writerは新mode fieldsだけを書き、legacy boolean dual-writeをしない。Web/Python version constantsとdocs/testsが43で一致する。 |
| AV-12 | architecture exceptionが必要な場合、正確なOE/PE/CB、positive legacy fixture、owner、removal conditionを記録し、gateを弱めない。 |

## 10. 検証戦略

### 10.1 unit／contract tests

必要な既存test ownerへケースを追加し、新しい総合test fileを安易に増やさない。

- `tests/web/depth-track-state.test.mjs`
- `tests/web/linear-sources.test.mjs`
- `tests/web/linear-record-layout.test.mjs`
- `tests/web/session-request.test.mjs`
- `tests/web/session-active-config-contract.test.mjs`
- `tests/web/history-inputs.test.mjs`
- `tests/web/settings-only-session.test.mjs`
- `tests/web/gallery-session-migration.test.mjs`
- `tests/web/architecture-contracts.test.mjs`
- Python Session tests: `tests/test_api_session.py`、documentation/version contracts

testはproductionのresolver式をコピーせず、入力stateと観測可能なrequest／Session結果をassertする。

### 10.2 browser journeys

実際のDOM controlsを操作する。app内部methodの直接呼出しだけでUI受入を済ませない。

- multi-record GenBankとpaired GFF3＋FASTA
- two Files、same-name independent uploads、middle File replacement
- Depth common → Mixed → clear → Undo
- Add/Clear/Delete/Cancel、sole card、blank last card、pending discovery
- Record Layout single rows → shared row → single rows
- Auto/Show/Hideの独立組合せ、Undo/Redo、Generate
- fresh Save → new page Load → Generate
- v42 positive fixture、Gallery sessions、settings-only session
- same-orientation／inversion／dense matchesのCurve表示
- 390×844 viewport、keyboard、focus、accessible names

現在sourceから`python tools/prepare_browser_wheel.py`でwheelを作り、実際にそのwheelを使った
Generateを確認する。wheelはgitignoredでありcommitしない。Node PlaywrightがなければPython
Playwrightを使う。Chromium sandbox errorは同じcheckを必要な権限で再実行し、実装不具合と
環境制約を区別する。

### 10.3 focused commands

変更範囲に応じて次を選ぶ。正確なcommandと結果を実施記録へ残す。

```bash
node --test tests/web/depth-track-state.test.mjs tests/web/linear-sources.test.mjs
node --test tests/web/linear-record-layout.test.mjs tests/web/session-request.test.mjs
node --test tests/web/session-active-config-contract.test.mjs tests/web/history-inputs.test.mjs
node --test tests/web/settings-only-session.test.mjs tests/web/gallery-session-migration.test.mjs
python -m pytest tests/test_api_session.py tests/test_documentation_contracts.py tests/test_documentation_reference_contracts.py -v
npx playwright test tests/web/linear-multi-record.playwright.spec.js tests/web/depth-track-session.playwright.spec.js --workers=1 --retries=0
node --test tests/web/architecture-contracts.test.mjs
git diff --check
```

### 10.4 final gates

```bash
node --test tests/web/*.test.mjs
node --test tests/ci/*.test.mjs
npm run test:web:comparison-contracts
node tools/check-web-change-budget.mjs --base origin/dev --head HEAD
git diff --check
```

Python runtimeを変更しないPhaseでは、Web/CI contractsとfocused browserを優先する。Session 43で
Python I/Oを変更した場合は関連Python testsを必須とし、影響が広い失敗があれば
`pytest tests/ -v -m "not slow"`へ拡大する。長時間testは少なくとも30分を許容し、増分監視する。

## 11. Architecture・Product evidence

各runtime sessionはproduction、tests、docs、generated artifactsのdiffを別々に読む。

通常の期待:

- #559: owner、canonical path、compatibility pathは非増加。UI入口を既存Depth actionへ収束。
- #560: source semanticsは`linear-sources.js`、coordinationは`app-setup.js`の一ownerずつ。
  新しいdialog stateはephemeralで、source of truthではない。
- #562 match-style: fresh default ownerとenum validatorを一つにし、UI再編やvisibility stateを
  巻き込まない。
- #562 organization: 既存stateのDOM配置だけを変え、default、schema、rendererを変更しない。
- #562 Auto: visibility policy ownerを一つ作る代わりに、UI/requestの独立式を作らない。
  v42 migrationによるCB増加はexception evidenceを必要に応じて提出する。

完全なbefore/after OE、PE、CB setsはArchitecture Ratchetのexception条件に該当する場合だけ作る。
ordinary non-increasing変更では、owner/path、削除したsuperseded path、gate結果を簡潔に示す。

Product evidenceは次を区別する。

- #559: runtime baseのPD-OI-025／OIC-020
- #560と#562の各concern: runtime baseへmerge済みの明示的なselected outcome。一つの#562
  decisionを別concernのauthorityとして使わない。
- tests: authorityではなく、選択された結果を守るevidence
- Issue本文:要望と再現根拠であり、それだけをdurable authorityと呼ばない

## 12. 非目標とrollback

非目標:

- comparison curveのgeometry、opacity、filter、algorithmの変更
- Depth renderer、axis、normalization、file formatの変更
- File drag-and-drop、arbitrary reordering、bulk delete
- Circular labelのAuto化
- CLI/API defaultのCurve化
- Session以外の新しい設定store
- UI framework、modal dependency、build stepの追加
- Gallery全体の再設計

rollbackはsession単位で行えるようにする。#559、#560、#562 match-style、organization、Autoを
別commit候補として保ち、Session 43 migrationをCurve defaultやUI regroupへ不可分に混ぜない。
ただしAuto session内ではversion
constant、migration、writer、docs、positive fixturesを部分的に戻せない一つの整合変更として扱う。

## 13. 完了条件

全体完了は、選択されたProduct outcomesに該当する全受入IDが成立し、次を満たした状態である。

- #560/#562 runtimeのbaseに必要なdurable authorityが含まれる。
- fresh、Reset、legacy Load、current Save/Load、Undo/Redo、Generateが区別して検証されている。
- UI summary、canonical request、SVGのeffective結果が一致する。
- current Result preservationとcomparison/cache reconciliationに回帰がない。
- Session version、Web/Python writer、supported versions、docs、fixturesが一致する。
- production／tests／docs／generated diffを別々にreview済みである。
- browser wheel、temporary screenshots、`examples/gbdraw_social_preview.png`をcommit対象にしていない。
- architecture／Product gatesを弱めていない。
- 本書第14節へ再現可能な実施記録がある。

未解決のProduct Decision、外部merge待ち、失敗したgateが一つでもあれば、依存する範囲を完了と
報告しない。独立して完了した#559などは、全体と区別して完了記録を残す。

## 14. 実施記録

各sessionはこの節に追記する。チャットだけに結果を残さない。

### 記録テンプレート

```text
Session:
Date / actor:
Base SHA / branch / upstream:
Authority SHA and selected outcome:
Pre-existing unrelated changes preserved:
Production files changed:
Test files changed:
Docs / generated files changed:
Acceptance IDs completed:
Commands and exact results:
Browser wheel/source identity:
Architecture owner/path evidence:
Product evidence:
Known limitations / failed checks:
Next session and start condition:
Proposed commit title (English):
Proposed commit summary (English):
```

計画策定時点ではruntime変更、test実行、Product Decisionの受領は行っていない。作業ツリーには
本計画と無関係な未追跡ディレクトリが存在する可能性があるため、実装者は毎回`git status`で
確認し、削除・stage・変更しない。

### 2026-09-22 実装記録

```text
Session:
  Linear UI #559 / #560 / #562 implementation and documentation evidence
Date / actor:
  2026-09-22 / Codex, with explicit Product Decision Owner authorization in the
  current conversation
Base SHA / branch / upstream:
  HEAD 17305e7f2e24426cc98268b4f9328626f6a82460
  origin/dev f4760476915194fc0586642e0b03ec561b96ab01
  fix/linear-ui-559-560-562-20260921
  origin/fix/linear-ui-559-560-562-20260921
Authority SHA and selected outcome:
  No additional durable-authority SHA was present in this runtime base. The
  Product Decision Owner explicitly authorized implementation in this
  conversation and confirmed that the Product decisions had been split. The
  implementation is limited to D560-A, CURVE-FRESH-RESET,
  REGROUP-TITLES-LABELS-LEGEND, and INDEPENDENT-AUTO-SHOW-HIDE. It does not
  serialize inferred rationale, retirement intent, residual-risk acceptance,
  or a fabricated BD identifier.
Pre-existing unrelated changes preserved:
  .agents/skills/write-clear-pull-request/SKILL.md
  .github/pull_request_template.md
  tests/web/pr-language.test.mjs
  tools/check-pr-language.mjs
  docs/internal/LOSAT_CONDA_GBDRAW_CODEX_PACKAGE_2026-09-20/ (untracked)
Production files changed:
  gbdraw/session_io.py; gbdraw/web/index.html; gbdraw/web/js/components.js;
  gbdraw/web/js/mode-profiles.js; focused owners under gbdraw/web/js/app/ and
  gbdraw/web/js/services/ for Linear source removal, record-row resolution,
  match-style defaults, label visibility, current Session validation/migration,
  canonical request projection, and Gallery promotion/publication.
Test files changed:
  Focused Python Session/docs contracts, Node pure-policy/Session tests, and
  Playwright tests for Depth, removal/History/focus, typography organization,
  Auto visibility, settings-only v42 migration, Gallery, and comparison paths.
Docs / generated files changed:
  Existing Web reference, Session compatibility, release notes, affected GUI
  tutorials/capture flows, deterministic screenshots, seven Gallery tutorials,
  their operation media, and the Gallery screenshot register. No new public
  page was added. examples/gbdraw_social_preview.png was not changed. The
  generated browser wheel remains gitignored and is not a commit target.
Acceptance IDs completed:
  DPT-01..DPT-08; RM-01..RM-12; MSD-01..MSD-04; ORG-01..ORG-06;
  AV-01..AV-12.
Commands and exact results:
  node --test focused Linear/Session modules: 10 passed.
  npx --yes node@20 --test --test-reporter=dot tests/web/*.test.mjs:
    606 passed, 1 failed only because the current sparse fixture omitted the
    new required visibility fields; after adding those fixture fields, the
    isolated session-draft-authority test passed.
  npx playwright test depth-track-session, history-inputs,
    linear-typography, settings-only-session, comparison-ui, and
    linear-multi-record with one worker: 78 passed in 7.4 minutes.
  npx playwright test tests/web/linear-typography.playwright.spec.js
    --workers=1 --retries=0 after final Auto-option review fix: 3 passed.
  npm run test:web:comparison-contracts: 16 passed in 2.6 minutes.
  node --test tests/web/architecture-contracts.test.mjs: 137 passed.
  npx --yes node@20 --test tests/ci/*.test.mjs: 59 passed.
  focused Python API/Session/docs/composition/run-info suite: 326 passed.
  documentation contract suite after capture-source updates: 50 passed.
  Gallery strict static checks: 13/13, 9/9, 19/19, 10/10, 11/11,
    11/11, and 12/12 for the seven affected tutorials.
  Gallery browser suite: 22/23 passed before updating the BGC operation count
    from 18 to 19; the corrected focused case passed.
  Deterministic --check capture results passed for T-GUI-01, T-GUI-02,
    T-GUI-04, T-GUI-05, T-GUI-06, H-GUI-02, H-GUI-09, H-GUI-11,
    H-GUI-12, and H-GUI-14. H-GUI-14 additionally validated a fresh-context
    v43 Session reload and reloaded_diagram.svg.
  ruff check gbdraw/: All checks passed.
  git diff --check: passed with no output.
Browser wheel/source identity:
  python tools/prepare_browser_wheel.py succeeded after the last runtime edit.
  gbdraw-0.14.0-py3-none-any.whl SHA-256:
  e35cfc9e435c78f19ff90686a8c6ca22e9e9ecd578b0172bba2bf175871ed405.
  The wheel is ignored and must not be committed.
Architecture owner/path evidence:
  Depth matrix remains the sole data owner; each File card adds only a native
  disclosure and an action into the existing matrix mutation.
  planLinearSourceRemoval is the pure removal planner; app-setup owns the one
  modal and one History transaction; existing reconciliation remains downstream.
  mode-profiles owns fresh/reset pairwise style and current-option-values owns
  its normalization.
  Titles/Record Labels/Legend reuse existing state and renderer keys.
  linear-label-visibility.js is the sole selected-mode/legacy/effective policy
  owner; linear-record-layout.js owns effective rendered-row facts;
  session-request.js consumes the resulting effective booleans.
  Architecture inventory: 152 modules, 482 edges, zero cycles; all four active
  hard owner/path rules passed. Production change: 15 files, +675/-316 before
  the final one-line option-label correction, requiring recorded size review.
  OE before: adv.linear_show_accession, adv.linear_show_length.
  OE after: adv.linear_accession_visibility, adv.linear_length_visibility.
  PE: v42 Web import, settings-only import, Gallery promotion/publication map
  true/false/missing to Show/Hide/historical Show and preserve historical
  missing pairwise style as Ribbon.
  CB: tests/fixtures/sessions/settings-only.v42.json.gz plus the existing v42
  Gallery documents.
  Canonical current owner: selected modes in active config, pure resolution in
  linear-label-visibility.js, effective booleans only in request schema 7.
  Compatibility removal condition: remove the v42 promotion path only after
  v42 is formally removed from supported Session versions and Gallery inputs.
Product evidence:
  Fresh Linear and Reset use Curve; historical missing values use Ribbon.
  Fresh Accession/Length modes are independently Auto and report
  Auto · Shown until an actually rendered shared row makes Auto hidden.
  Manual Show/Hide is independent and survives row changes. v42 booleans and
  missing fields retain their historic outcomes; v43 is the sole current writer
  and does not dual-write legacy booleans. Gallery screenshots sourced from v42
  correctly display Show rather than rewriting those explicit historic values
  to Auto. Replicon remains the intentionally separate boolean control.
Known limitations / failed checks:
  WEB_ARCHITECTURE_CHANGE=true node tools/check-web-change-budget.mjs --base HEAD
  reported its sole blocking violation because production runtime and the
  pre-existing unrelated .github/pull_request_template.md modification coexist
  in the working tree. No privileged owner, canonical path, cycle, vendored
  runtime, or dependency violation was found; size review remains required.
  The following broader documentation scenarios still fail pre-existing
  semantic expectations that were not weakened: T-GUI-08 expected 500
  Collinear matches but observed 462; T-GUI-10 lacked the required named gene
  labels; T-GUI-12 lacked GC content (%); H-GUI-10 lacked AT skew (+/-);
  H-GUI-16 stopped on an existing record-scoped rendered-feature identity
  mismatch after placement editing. Partial outputs from failed captures were
  restored and are not included.
Next session and start condition:
  A separate cleanup may investigate the five unrelated documentation semantic
  failures. Publication/commit should begin only after reviewing this diff and
  excluding the preserved unrelated files listed above.
Proposed commit title (English):
  Implement Linear file controls, defaults, and label visibility
Proposed commit summary (English):
  Make per-file Depth and removal workflows discoverable and atomic, adopt the
  selected Linear Curve and Auto label policies with Session 43 migration, and
  reorganize title, record-label, and legend controls with updated tests and
  reproducible documentation evidence.
```
