# Work package B: LOSAT・比較UIの実装計画

- Date: 2026-08-10
- Status: planned
- Baseline: `docs_renovation` / `3d1d9d287cc9`
- Source: `gbdraw v0.14.0 Release Roadmap and Codex Implementation Brief` の Work package B

この文書は、Work package B「Clean up the LOSAT and comparison UI」を実装可能な作業単位へ分解した計画である。文書作成時点ではproduction code、test、Gallery artifactを変更していない。

## 1. 結論

入力を最初に完了できる画面へ変える。比較は、明示的な実行意図、現在の設定概要、必要時だけ開く詳細設定の3層に分ける。

Linearでは既存の`linearComparisonPlan`と`linearComparisonResolution`を使い続ける。新しい比較state、request field、session versionは追加しない。

新規Linear sessionとReset Settingsの既定値は`No comparison`に変更する。2件目のsequenceを追加しただけでLOSAT intentが生じる現在の挙動を廃止する。保存済みWeb sessionは保存された比較意図を復元し、CLI-only replayはcommitted renderとeditable Web draftを分ける。

## 2. 目標と対象範囲

### 2.1 目標

- 1280 x 720で最初のsequence uploaderと上部の`Add sequence`をスクロールなしで操作できる。
- 比較を使わないGenerate/UpdateではLOSAT、比較file、cache、過去artifactに触れない。
- 選択したprogramとpresentationに関係するcontrolだけを表示する。
- pair編集、thread、cache、raw resultは必要な利用者だけが開ける。
- accordionを閉じても値を失わず、sessionの意味状態と開閉状態を混ぜない。
- fresh、Reset、session replayで同じLinear実行許可の規則を守る。

### 2.2 対象

- Linearのsequence入力、comparison summary、LOSAT settings、pair editor、record layoutの情報設計
- fresh stateとReset SettingsのLinear既定値
- program、LOSATP presentation、filter、appearance、runtime設定のprogressive disclosure
- session save/load、inactive draft、cache、canonical requestの回帰防止
- keyboard、focus、error表示、desktop/narrow viewport
- 変更した操作説明、Gallery tutorial、capture script、screenshot

### 2.3 対象外

- LOSATN、TLOSATX、LOSATPの検索アルゴリズム
- similarity-groupまたはcollinear-blockの生物学的意味と出力形式
- renderer geometry、SVG metadata、canonical render request schema 5
- Python API、CLI、worker protocol、public config key
- session versionの更新または新しいmigration chain
- Work package Eのpreview再描画、debounce、LOSAT rerun policy
- Circular `Pairwise Comparisons`の情報設計、明示opt-out、managed conservation slot、Linearとのstate統合
- 新しいUI framework、build step、CDN dependency

Work package Bのroadmap例は、Linearのadjacent/selected pair、LOSATP presentation、record timelineを対象としている。Circularは現在details内で初期closedであり、Linearと同じviewport問題を起こしていない。Circularへ`No comparison`を追加するには、`enabled`だけでなくmanaged conservation slot、user-uploaded TSV、生成済みLOSAT artifact、session replayのownerを分離する必要がある。このstate再設計をUI整理へ混ぜない。releaseでCircular parityを必須にする場合は、別計画を先に作る。

## 3. 現状の根拠

2026-08-10のplanning auditで次を確認した。Phase 0では同じ測定を実装開始時のworktreeで取り直す。

| 現状 | 根拠 | Work package Bへの影響 |
| --- | --- | --- |
| Linearのfresh planは`adjacent + losat` | `gbdraw/web/js/app/linear-comparisons.js`の`createDefaultLinearComparisonPlan()` | 2件目を追加すると、明示操作なしでLOSAT intentが生じる |
| Resetも`adjacent + losat`へ戻る | `gbdraw/web/js/services/reset.js`の`resetLinearComparisonPlan()` | Reset後の次回Generateで比較が再開し得る |
| 比較選択、LOSAT実行設定、program設定、cache操作、Record Layoutが最初のuploaderより前にある | `gbdraw/web/index.html`のLinear `Input Genomes` block | sequence入力が設定群に押し下げられる |
| 1280 x 720では初期uploader上端が約398 px、`Add Seq`は最初のrecord card終端約746 pxより下 | Chromiumでのplanning measurement | 追加操作が初期viewportに入らない |
| 2件目追加後はLOSAT panelが約281から639 pxを占め、最初のuploader上端が約772 pxになる | 同measurement | comparison settingsがprimary inputを完全に押し出す |
| global computed stateはselected planで`selected`を返すが、DOMにはnone/losat/uploadのradioしかない | `gbdraw/web/js/app/app-setup.js`と`gbdraw/web/index.html` | selectedまたはmixed planで全radioが未選択に見える |
| `none`またはactive edgeなしでは、canonical requestがdormant fileやartifactを見る前に`comparisons: []`を返す | `gbdraw/web/js/services/session-request.js` | Linearの実行許可境界は再利用できる |
| LOSAT runtime warmupとexecutorはresolved planの`hasLosatIntent`に依存する | `gbdraw/web/js/app/run-analysis.js` | UIから新しいruntime flagを作る必要はない |
| comparison plan、LOSAT settings、Circular conservation、record layoutはsession configへ保存される | `gbdraw/web/js/services/config.js` | 意味状態は既存schemaで保存できる |
| disclosure open/closedはsessionへ保存されない | 現行session UI projection | accordionは非永続のままにできる |
| Circularでは`enabled`をfile数から書くwatcher、managed conservation slot、uploaded/derived comparison fileが結合している | `gbdraw/web/js/app/app-setup.js`、`circular-track-slots.js`、`services/config.js` | Linearと同じ3 actionを見た目だけ移植できないため本計画から外す |
| DOM順と`Apply to all adjacent gaps`を固定するtest、tutorial、capture scriptがある | `tests/test_web_ux_profile.py`、`tests/web/linear-multi-record.playwright.spec.js`、`docs/capture`、`gbdraw/web/gallery/tutorials` | markup変更と同時にtestと再生成手順を更新する必要がある |

## 4. 実装前に固定する製品判断

### D1. fresh stateとResetは`No comparison`

`createDefaultLinearComparisonPlan()`は`mode: 'none'`を返す。Reset Settingsも`mode: 'none'`へ戻し、保持対象のupload fileとcustom raw filenameはinactive draftとして残す。

この判断は`docs/internal/IMPLEMENTATION_PLAN_316_LOSAT_OPT_OUT.md`にある「new sessionとResetは`adjacent + losat`」を置き換える。実装時に同文書へsuperseded noteを追加し、相反する計画を残さない。

既存Web sessionの比較意図は変えない。current-writer sessionに保存された`none | adjacent | selected`はそのまま復元する。supported pre-40 Web sessionは既存migration規則を維持する。disabledまたはabsent legacy layoutは従来のadjacent intent、enabled explicit listは`selected`、enabled empty listまたは明示`blastSource: none`は`none`へ移行する。

CLI-only replay sessionは例外とする。これは意図的にsynthetic Web comparison draftを持たないため、保存済みcanonical renderはそのままreplayするが、Webで次のGenerateを行うdraftはfresh defaultの`none`から始める。CLI requestを推測してWeb draftへ変換しない。

`normalizeMode()`のinvalid/missing value fallbackはこのUI既定値変更に巻き込まない。fresh stateは明示的なfactory、旧sessionはmigration、壊れた入力はnormalizerという3つの境界を別々にtestする。

### D2. top-levelの3択はradioではなくaction control

表示する選択肢はroadmapどおり次の3つとする。

```text
[ No comparison | Run LOSAT | Upload BLAST TSV ]
```

Linearには`selected`とmixed-source planがあるため、この3つはstate全体を列挙するradio groupではない。`Set all adjacent comparisons`というaccessible nameを持つcommand groupとして実装する。buttonのaccessible nameは`Set no comparison`、`Run LOSAT for all adjacent pairs`、`Use uploaded BLAST TSV for all adjacent pairs`とし、toggle radioを表す`aria-pressed`は使わない。

現在値はbuttonとは別のstatus行へ常に表示する。adjacent/noneでは`Current: No comparison`などを表示し、selected planでは`Current: Selected pairs (N; X LOSAT, Y upload)`を表示する。見た目上も`Custom` status badgeを付ける。これにより、3つのbulk commandが未選択に見える状態と、現在のcustom intentを区別する。

global actionを押したときだけ、active topologyを`none`または`adjacent`へ切り替える。edge draft、upload、custom filenameは削除しない。

### D3. search methodは既存fieldのUI projection

画面上は次のmethodを一つのselectorとして見せる。

| 表示値 | 既存stateへの写像 |
| --- | --- |
| `LOSATN` | `losatProgram = 'blastn'` |
| `TLOSATX` | `losatProgram = 'tblastx'` |
| `LOSATP Pairwise` | `losatProgram = 'blastp'`、`losat.blastp.mode = 'pairwise'` |
| `Similarity groups` | `losatProgram = 'blastp'`、`losat.blastp.mode = 'orthogroup'` |
| `Collinear blocks` | `losatProgram = 'blastp'`、`losat.blastp.mode = 'collinear'` |

新しいpersisted enumは作らない。selector actionは上記2 fieldだけを変更し、inactive methodの設定値を初期化しない。`Collinear blocks`への切替時に`adv.pairwise_match_style`を直接`curve`へ書き換える現行side effectは削除し、保存済みappearanceを維持する。

selected/mixed planでLOSAT edgeがある場合、表示するmethodは`LOSATN`、`TLOSATX`、`LOSATP Pairwise`に限る。`Similarity groups`と`Collinear blocks`はall-adjacent LOSATのpresentationであり、selected topologyでは既存resolverが拒否する。UIはこの2 optionをdisabledにして理由と`Use all adjacent LOSAT` actionを示す。actionを押した場合だけtopologyをadjacentへ変更する。自動変換しない。

### D4. disclosureは意味状態ではない

`Settings`、`Selected pairs (N)`、`Advanced and reproducibility`、各recordの`Record options`はnative `<details>/<summary>`を使う。open stateはconfig、session `ui`、undo snapshot、render requestへ追加しない。開閉だけではvalidation、cache probe、runtime warmupを実行しない。

### D5. Circularへ見た目だけのparityを作らない

本計画ではCircular `Pairwise Comparisons`を変更しない。Linear用buttonやsummary helperをCircular stateへ接続しない。Circularの明示opt-outを実装する別計画では、少なくとも次を決める。

- active intentとfile-derived watcherの整理
- managed conservation slotをinactiveにしたときのorder、axis、appearance保持
- user-uploaded TSVとgenerated LOSAT resultを別ownerにする方法
- last committed renderとeditable draftが異なるsessionのvalidation
- Upload、LOSAT、noneを往復した後のsave/load

### D6. 一つの設定に一つのvisible ownerを置く

現在のLOSATP blockと離れた`Pairwise Match` disclosureに分散しているfilterとappearanceを整理する。同じ`v-model`を複数箇所へ複製しない。markupを移動し、active modeで一つだけ表示する。

### D7. 自動open/focusはstructured comparison issueに限定する

`linearComparisonResolution.errors`は`code`、`edgeId`、`edgeKey`を持つ。この既存構造だけをfocus routingに使う。

- `missing-upload`は`Selected pairs`と対象edgeのuploaderへ対応させる。
- `selected-losat-requires-pairwise`は`Settings`とsearch methodへ対応させる。
- duplicate、missing UID、self、same-row、non-adjacentは`Selected pairs`と対象edge rowへ対応させる。

`run-analysis.js`のplain error stringを文面で分類しない。TLOSATX gencodeなどstructured targetを持たないerrorは現行のcentral error logを維持する。別のerror型をWork package Bだけのために追加しない。

## 5. 目標UI

### 5.1 Linear

```text
Input sequences                                      [ + Add sequence ]
[ GenBank | GFF3 + FASTA ]

Sequence 1
[ Choose GenBank file ]
Record options ▸

Comparison
[ No comparison | Run LOSAT | Upload BLAST TSV ]
LOSATN · 3 adjacent pairs · filters shown as values     [ Settings ▸ ]
Selected pairs (3) ▸

Basic figure settings
...

[ Generate / Update preview ]    (既存fixed action bar。DOM上はここ)

Advanced comparison and layout ▸
  Record Layout
  Runtime and reproducibility
```

実際のDOM順とtab orderもこの順に合わせる。既存Generate barはfixed表示を維持するが、DOM anchorをBasic figure settingsの後、Advanced controlsの前へ移す。buttonを複製しない。CSSだけで見た目を並べ替えない。`Add sequence`はsection headerに置き、長いrecord listの末尾にも既存buttonを残す。上部buttonは常に見えるが、最大件数など既存のdisable条件を共有する。

各record cardはuploaderを先頭に置く。organism/strain、subtitle、region、depth、genetic codeなどの補助fieldは`Record options`へ入れる。ただしactiveなTLOSATXが必須とするgenetic-code errorは閉じたsummaryにも表示する。

### 5.2 active modeごとの表示

| Active mode | 通常の`Settings`に表示するもの | 非表示またはadvancedへ移すもの |
| --- | --- | --- |
| `No comparison` | summaryのみ | 全LOSAT、upload、filter、appearance、runtime control |
| `LOSATN` | search method、task、共通result filter、comparison appearance | protein、gencode、grouping、collinearity control |
| `TLOSATX` | search method、共通result filter、active recordのgencode | blastn task、protein presentation control |
| `LOSATP Pairwise` | search method、max hits、共通result filter、appearance | grouping、collinearity control |
| `Similarity groups` | search method、groupingに必要なoption、共通result filter | pairwise max hits、collinearity control |
| `Collinear blocks` | search method、min anchors、max gap、evidence scope、color mode、共通result filter | pairwise max hits。diagonal driftとmerge conflictはadvanced |
| `Upload BLAST TSV` | active pairのfile readiness、共通result filter、appearance | LOSAT program、search implementation、LOSAT runtime control |
| selected/mixed | LOSATN、TLOSATX、LOSATP Pairwise、共通filter、source内訳。file assignmentはpair section | Similarity groups、Collinear blocks。global radio風のselected表示 |

compact summaryはactiveな値から求める。default判定をする場合はmode-profileの既存defaultを参照し、UI moduleへ数値を複製しない。共有default ownerがないfieldは`Default filters`と推測せず、`E-value <= ...`など実値を短く表示する。

### 5.3 Advanced and reproducibility

次を初期状態で閉じたsectionへ移す。

- execution mode
- total thread budget
- parallel runs
- threads per run
- job estimateとworker estimate
- cache status、clear、reuse説明
- raw result filename、download、retained artifactの状態
- collinear diagonal drift、merge conflictsなど通常変更しないsearch detail
- Record Layout

pair membership、source、batch action、custom endpointの編集は`Selected pairs (N)`へ置く。既存の`linearComparisonTimeline`を表示データとして再利用し、別のpair topologyを作らない。

## 6. 状態とsurfaceの契約

| Surface | 現在のowner | 既定値またはinactive値 | 保存 | validation / normalization | 最終consumer | 計画上の変更 |
| --- | --- | --- | --- | --- | --- | --- |
| Linear比較intent | `linearComparisonPlan` | fresh/resetを`none`へ変更 | session config | `normalizeLinearComparisonPlan()`とresolver | request作成、LOSAT job planner | ownerは維持。defaultだけ変更 |
| Linear active snapshot | `linearComparisonResolution` | edge 0ならLOSAT intentなし | renderごとのsnapshot | topology、source、file readiness | `run-analysis.js`、`session-request.js` | UI summaryの入力に使う。変更しない |
| pair表示 | `linearComparisonTimeline` | planとrecord orderから導出 | 保存しない | edge UIDとrow topology | pair editor | collapsed sectionへ移す |
| Linear UI summary | 新しいpure projection | plan、resolution、programから導出 | 保存しない | projection unit test | `index.html` | label、count、visibilityだけを返す |
| disclosure open state | DOMの`<details>` | closed | 保存しない | browserのnative動作 | UIのみ | stateへ追加しない |
| LOSAT method | `losatProgram`と`losat.blastp.mode` | 現行mode profile | session config | 既存normalizer | LOSAT job planner | 一つのsearch-method selectorへ投影 |
| filters / appearance | `adv`と`losat`の既存field | mode-profile default | session config | active intent/methodだけを既存normalizerへ渡す | requestとrenderer | visible ownerを一本化。noneではdormant値を変更しない |
| Linear upload/raw draft | comparison edgeとsession resource binding | inactive可 | config metadataとresources | fileActive、losatFilenameActive | active edgeだけrequestへ | noneやmode切替で削除しない |
| cache | 既存LOSAT cache | content/program/argument key | 既存session/cache契約 | cache validator | LOSAT executor | UIをadvancedへ移すだけ |
| canonical request | schema 5 `comparisons` | no intentなら`[]` | committed render/session | request builder | worker/Python renderer | schemaを変更しない |
| CLI / Python / renderer | 現行public contract | 現行default | 各surfaceの既存契約 | 既存validator | diagram assembly | 変更しない |

意味状態の流れは次を維持する。

```text
UI action
  -> editable comparison intent
  -> immutable resolution at Generate
  -> active file reads / LOSAT jobs only when permitted
  -> canonical schema-5 request
  -> committed render result
```

summary、accordion、buttonの押下表示はこの流れへ新しいauthorityを加えない。

## 7. 実装フェーズ

### Phase 0: baselineと変更境界を固定する

Status: pending

#### 作業

1. branch、HEAD、dirty fileを記録し、production、test、documentation、generated artifactの既存差分を分ける。
2. 1280 x 720と390 x 844で、初期表示と2件目追加後のDOM順、bounding box、scroll position、tab orderを記録する。
3. Linearについて、active controlがbindする既存fieldを一覧化する。同じfieldを編集するDOMが複数ある場合は移動対象を決める。
4. fresh、Reset、current Web session、supported legacy Web session、CLI-only replayの比較intentを別fixtureで記録する。
5. focused testを変更前に実行し、pre-existing failureを分ける。

#### Baseline command

```bash
git status --short --untracked-files=all
git diff --stat
git rev-parse --short=12 HEAD
TMPDIR=/tmp node --test tests/web/linear-comparisons.test.mjs tests/web/losat-settings.test.mjs
TMPDIR=/tmp node --test tests/web/session-request.test.mjs tests/web/session-draft-authority.test.mjs
TMPDIR=/tmp node --test tests/web/run-analysis-simple-path.test.mjs
pytest tests/test_web_ux_profile.py -q
```

#### 完了条件

- viewport measurementを再現するか、planning auditとの差を記録できる。
- fresh defaultとlegacy restoreの期待値を別々に固定できる。
- production change前のfailureを記録できる。

### Phase 1: 明示的な比較intentを固定する

Status: pending

#### 対象owner

- `gbdraw/web/js/app/linear-comparisons.js`
- `gbdraw/web/js/services/reset.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/gallery-session-migration.js`（既存migrationの確認。default couplingがある場合だけ変更）
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/services/session-request.js`
- 必要な既存unit test

#### Linear

1. fresh factoryを`mode: 'none'`へ変更する。`defaultSource`はinactive draftの次回明示選択に使えるため`losat`のまま保持してよいが、active intentにはしない。
2. Resetはretained edgeの`included`、`fileActive`、`losatFilenameActive`をfalseにし、plan modeを`none`へする。file objectとcustom filenameは保持する。
3. supported pre-40 Web sessionのmigration matrixを維持する。disabled/absent layoutはadjacent、enabled explicit listはselected、enabled empty listまたは明示noneはnoneとする。
4. current-writer sessionは保存済みplanを変更せず復元する。session versionは上げない。
5. CLI-only replayはcanonical renderを維持し、editable Web planを合成しない。次のWeb Generateではfresh noneを使うfixtureを追加する。
6. global actionは既存invalidation ownerを通す。deep watcherは追加しない。
7. `activeComparisonPlanSnapshot`をcomparison-only validationとnormalizationのpermission boundaryにも使う。`hasComparisonIntent`はPairwise Match Heightと共有appearance/filter、`hasLosatIntent`はLOSAT固有設定、LOSATP method/submodeは各protein設定をgateする。
8. noneでは`adv.comparison_height`をvalidateせず、canonical `comparisonHeight` overrideをemitしない。dormant filter、appearance、LOSAT method値をnormalizeまたは初期化しない。

#### 完了条件

- freshとResetは`none`になる。
- saved current Web sessionとsupported legacy Web sessionのactive intentが変わらない。
- CLI-only replayのcommitted renderは変わらず、次のWeb Generateはnoneになる。
- LinearのnoneではGenerate/Update中のruntime、cache、comparison file read、active comparison resourceが0になる。
- noneはinvalidなdormant Pairwise Match HeightでもGenerateでき、comparison-only draft値を変更しない。
- inactive draftはmodeを往復しても残る。

### Phase 2: presentation-only projectionを切り出す

Status: pending

#### 対象owner

- 新規 `gbdraw/web/js/app/comparison-ui.js`
- `gbdraw/web/js/app/app-setup.js`
- 新規 `tests/web/comparison-ui.test.mjs`

#### 方針

`comparison-ui.js`はpure helperを持つ。最初からsubfolderを増やさず、summaryとvisibilityの責務を分ける必要が生じた場合だけ`app/comparison-ui/`へ分割する。

projectionは次を返す。

- current intent key: `none | losat | upload | custom`
- topology labelとactive pair数
- LOSAT/upload/mixedのsource内訳
- active search method key
- compact filter summary
- retained inactive draft数
- Settings、Selected pairs、Advancedに表示するsection key
- error badgeと、開くべきdisclosure key

projectionはstateを変更せず、Fileを読み取らず、cache APIを呼ばない。labelのためにdomain resolverを再実装しない。pair数とsourceは`linearComparisonResolution`、draft数はnormalized plan、program表示は既存LOSAT stateから求める。focus targetはD7で列挙したresolver issue codeだけを扱う。

search-method actionは既存fieldをまとめて変更する薄いmutationにする。mode switch時に`losat` objectを作り直さない。collinearへの切替もappearance draftを上書きしない。

#### Unit test

- fresh none、adjacent LOSAT、adjacent upload、selected upload、selected mixedのsummary
- recordが1件のadjacent LOSATで、選択intentとactive job 0を区別する表示
- LOSATN、TLOSATX、LOSATP 3 presentationのvisibility matrix
- selected/mixedではSimilarity groupsとCollinear blocksをactive methodにできず、明示actionだけがadjacentへ変える
- defaultとcustom filter summary
- dormant draft countがactive pair countへ混ざらない
- projection呼出し前後で入力objectが同一内容である

#### 完了条件

- `app-setup.js`とtemplateにsummary規則が重複しない。
- projection outputはsession、request、undoへ保存されない。
- selected/mixed planが未選択radioに見える問題を解消できるmodelになる。

### Phase 3: DOM順とprogressive disclosureを変更する

Status: pending

#### 対象owner

- `gbdraw/web/index.html`
- `gbdraw/web/js/app/app-setup.js`
- 必要なfocused CSS classとcomponent prop
- `tests/test_web_ux_profile.py`

#### 作業

1. Linear `Input Genomes`のsection headerへ`Add sequence`を置き、その後をinput type、record listの順に組み直す。各record cardではuploaderを`Record options`より前に置く。
2. comparison cardをrecord listの後へ移し、3 action、summary、`Settings`を置く。
3. 現在record境界へ挿入しているpair UIをrecord timelineから外し、`Selected pairs (N)`へ移す。既存`linearComparisonTimeline`のrow、edge key、batch actionを使い続ける。
4. `Record Layout`をprimary uploadより前から外し、`Advanced comparison and layout`へ移す。
5. LOSAT runtime、thread、estimate、cache、raw filename、downloadを`Advanced and reproducibility`へ移す。
6. fixed Generate barのDOM anchorをBasic figure settingsの後、Advanced controlsの前へ移す。表示位置と既存actionは変えない。
7. active methodに関係するcontrolだけをDOMへ出す。CSSで見えなくするだけではなく、`v-if`でinactive controlをtab orderとvalidation対象から外す。
8. Pairwise Match sectionにある共有filterとappearanceをcomparison Settingsへ移し、旧markupを削除する。
9. primary labelは12から14 px相当とする。補助説明だけを10 px相当にする。既存sidebar幅でtruncateする内容には`title`または隣接summaryを付ける。

#### disclosure動作

- 初期状態はclosed。
- native keyboard操作を維持する。
- closeでは値、file、cache metadataを変更しない。
- D7のstructured comparison issueがclosed section内を指す場合、そのdetailsを開いて対象controlへfocusする。
- errorの短い説明と件数はclosed summaryにも表示する。

#### 完了条件

- 実DOMで最初のuploaderがcomparison settingsより前にある。
- 1280 x 720で最初のuploaderとheader `Add sequence`が見える。
- 2件目追加後もSettingsはclosedで、primary inputを押し下げない。
- 390 x 844で横scrollがなく、segment、summary、buttonが自然にwrapする。
- selected/mixed planはsummaryとpair countが明示される。
- DOMとtab orderはInput、Comparison、Basic、Generate、Advancedの順になる。

### Phase 4: state、session、requestの回帰testを追加する

Status: pending

#### 対象test

- `tests/web/linear-comparisons.test.mjs`
- `tests/web/comparison-ui.test.mjs`
- `tests/web/session-request.test.mjs`
- `tests/web/session-draft-authority.test.mjs`
- `tests/web/session-active-files.test.mjs`
- `tests/web/session-resources.test.mjs`
- `tests/web/run-analysis-simple-path.test.mjs`
- `tests/web/losat-settings.test.mjs`
- `tests/web/gallery-session-migration.test.mjs`
- `tests/web/mode-profiles.test.mjs`
- verification-onlyの`gbdraw/session_io.py`と`tests/test_session_io.py`

#### 必須scenario

| Scenario | Assertion |
| --- | --- |
| fresh Linear | plan modeはnone、comparison summaryはNo comparison |
| Reset with retained upload/raw name | planはnone、payloadはinactiveで残る、requestには出ない |
| load saved adjacent/selected plan | 保存済みmode、source、edge、program、filterを維持する |
| pre-40 Web, layout disabled/absent | 従来どおりadjacent intentへ移行する |
| pre-40 Web, enabled list/empty/explicit none | listはselected、emptyとnoneはnoneへ移行する |
| pre-40 CLI-only replay | committed renderは維持し、synthetic Web planを作らず、次のWeb Generateはnone |
| adjacent LOSATN -> Collinear -> LOSATN | 各methodのdraft値とpairwise appearanceを維持する |
| Generate/Update none with stale Linear file/cache/artifact | comparison file read、cache probe、Pyodide、LOSAT executor、active comparison resourceは0 |
| Generate none with invalid height/dormant LOSAT values | Generateは成功し、height overrideをemitせず、dormant値は変更されない |
| selected mixed plan | active edgeだけをrequestへ出し、summaryのsource内訳が一致する |
| selected mixed + grouping/collinear | methodは選択不能。明示bulk action後だけadjacentへ変わる |
| disclosure close/open | semantic stateとsession serializationが同一 |

session round-trip testでは、open/closedをJSONへ追加していないことも確認する。Save Sessionはinactive file bytesを読み、既存resource inventoryと`webFiles.bindings.linearComparisons`へ保存してよい。Generate/Updateのactive canonical requestはそれを読まず、`comparisons: []`のままとする。この2つの境界を別testで確認する。

#### 完了条件

- current writerとsupported readerのversionを変えずにround-tripが通る。
- request schema snapshotにfield追加がない。
- no-comparison testは単なるDOM非表示ではなく、Generate/Updateのruntime/file/active-resource counterが0であることを検証する。

### Phase 5: browser UXとaccessibilityを検証する

Status: pending

#### 対象test

- `tests/web/linear-multi-record.playwright.spec.js`
- 必要なら新規 `tests/web/comparison-ui.playwright.spec.js`
- `tests/test_web_ux_profile.py`

#### browser scenario

1. 1280 x 720でfresh Linearを開き、最初のuploaderとheader `Add sequence`のbounding boxが設定paneのvisible area内にあることを確認する。
2. 2件目を追加し、LOSAT Settingsがclosedのまま、最初のuploaderとheader `Add sequence`がvisible area内に残ることを確認する。2件目のrecord cardまでの間にcomparison detailが挿入されないことも確認する。
3. Tab、Space、Enterだけで3 bulk actionと全summaryを操作する。focus-visibleが消えないことを確認する。button groupとcurrent-status elementのrole、accessible name、stateを別々にassertする。
4. 各search methodで、表示対象controlと非表示対象controlのmatrixを確認する。
5. selected/mixed sessionをloadし、summary、pair数、source内訳、bulk command group、current statusを確認する。
6. Settingsを閉じたselected planでmissing uploadとnon-Pairwise LOSATPをGenerateし、error summary、details open、対象controlへのfocus移動を確認する。
7. noneへ切り替えた後にGenerateし、dormant file/cacheが読まれないことをbrowser側counterで確認する。
8. 390 x 844で横overflow、DOM順、tab order、fixed Generate barとの重なりを確認する。

#### Playwright経路

```bash
command -v playwright && playwright --version
node -e "console.log(require.resolve('@playwright/test'))"
python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
npx playwright test tests/web/linear-multi-record.playwright.spec.js
npx playwright test tests/web/comparison-ui.playwright.spec.js
```

新規specを既存specへ統合した場合は2つ目のcommandを省く。Node runnerがない環境ではPython Playwrightで同じtargeted checkを行い、browser sandbox failureは必要な権限で同じcheckを再実行する。

#### 完了条件

- roadmapのviewport、keyboard、program visibilityのacceptanceをbrowserで検証できる。
- selectorやtextの存在だけを調べるstatic testでviewport acceptanceを代用しない。

### Phase 6: user documentationとGalleryを更新する

Status: pending

#### 対象

- `docs/FAQ.md`
- `docs/REFERENCE/web-app.md`
- `docs/REFERENCE/comparison-programs-thresholds-and-results.md`
- Linear comparisonを扱う`docs/TUTORIALS/GUI/`
- `docs/capture/flows/`
- `gbdraw/web/gallery/tutorials/`
- `docs/internal/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md`
- 関連capture contract test
- `docs/internal/IMPLEMENTATION_PLAN_316_LOSAT_OPT_OUT.md`のsuperseded decision note

#### 方針

1. `Apply to all adjacent gaps`、`LOSAT Settings`、Linearの`Add Seq`に依存する手順とselectorを検索し、新しいlabelと操作順へ更新する。
2. `No comparison`がfresh defaultであり、LOSATを実行するには`Run LOSAT`を明示選択することを書く。
3. selected planでは3 buttonがbulk actionであり、現在のcustom pair stateはstatus、summary、pair sectionに表示されることを説明する。
4. screenshotを手編集しない。`web-gallery-screenshot-maintenance` skillを使い、各tutorialのcapture recipeから実GUIを再撮影する。
5. screenshotのalt textとcaptionは、開いたSettings、active method、見えているcontrolを正確に書く。
6. affected tutorial sessionを再生成する場合はgeneratorから行い、session migrationとmedia contractを再検証する。

Circular `Pairwise Comparisons`のlabelと操作は本計画で変えない。検索結果を一括置換せず、Linear操作だけを更新する。

Gallery作業では`web-gallery-screenshot-maintenance` skillに従う。まずtutorial JSONが参照するmediaだけをinventoryし、capture metadataへexact session、expected app state、visible control/textを記録する。data-dependent cropを他exampleから流用しない。real UIをdevice scale factor 2以上、dense formは3で撮影し、更新前後を同じ表示寸法で比較する。desktopとmobileのGallery表示も目視する。

#### Gallery verification

```bash
python -m json.tool gbdraw/web/gallery/tutorials/hepatoplasmataceae_orthogroup.json >/tmp/hepatoplasmataceae_orthogroup.json.check
python tools/capture_gallery_tutorial_screenshots.py --example hepatoplasmataceae_orthogroup --check
node --check gbdraw/web/gallery/gallery.js
npx playwright test tests/web/gallery-tutorial.playwright.spec.js --project=chromium
```

上記を基準に、Phase 0で特定した各affected Linear tutorialへJSON checkとcapture `--check`を実行する。bitmapを変更したtutorialはclean captureを実行し、caption、alt text、selected value、cropの読みやすさを目視確認する。

#### 完了条件

- UIに存在しないlabelやDOM位置へ誘導する文書がない。
- capture scriptをclean runし、更新した画像を目視確認できる。
- operation register、tutorial JSON、Markdown、capture contractが一致する。

### Phase 7: 最終verificationとdiff audit

Status: pending

#### Focused gate

```bash
TMPDIR=/tmp node --test tests/web/linear-comparisons.test.mjs tests/web/comparison-ui.test.mjs tests/web/losat-settings.test.mjs tests/web/mode-profiles.test.mjs
TMPDIR=/tmp node --test tests/web/session-request.test.mjs tests/web/session-draft-authority.test.mjs tests/web/session-active-files.test.mjs tests/web/session-resources.test.mjs tests/web/gallery-session-migration.test.mjs
TMPDIR=/tmp node --test tests/web/run-analysis-simple-path.test.mjs
pytest tests/test_web_ux_profile.py tests/test_gui_nucleotide_comparison_capture_contracts.py tests/test_gui_protein_comparison_capture_contracts.py tests/test_documentation_capture_contracts.py tests/test_gallery_capture_contracts.py -q
pytest tests/test_session_io.py -q
npx playwright test tests/web/linear-multi-record.playwright.spec.js
```

#### Broad gate

```bash
TMPDIR=/tmp node --test tests/web/*.test.mjs
pytest tests/ -v -m "not slow"
pytest tests/test_output_comparison.py::TestOutputComparison -v
```

Gallery screenshotを変更した場合は対象capture flowに加え、`tests/web/gallery-tutorial.playwright.spec.js`と`tests/web/gallery-session-regeneration.playwright.spec.js`も実行する。browser wheelを使うtestが失敗した場合だけ`python tools/prepare_browser_wheel.py`でgitignored wheelを準備する。cache-bust tokenはdeployable bundleを作る別作業まで更新しない。

#### Diff audit

1. production diffでstate owner、mutation、request境界が重複していないことを確認する。
2. test diffで既存assertionを削った箇所を確認し、同じcontractを新しいUI構造で検証していることを示す。
3. documentation diffとgenerated screenshot diffを分けて確認する。
4. `tests/reference_outputs/`が変更されていないことを確認する。UI変更だけでreference SVGを再生成しない。
5. session version、request schema、Python signature、CLI helpに差分がないことを確認する。

#### 完了条件

- focused gateとbroad gateが通る。
- screenshotは再生成証拠と目視確認を持つ。
- production、test、documentation、generated artifactの差分を別々に説明できる。
- pre-existing failureがある場合は、再現commandとWork package Bとの差を記録する。

## 8. Acceptance criteria

### 情報設計

- Linearのbase sequence uploaderはcomparison settingsより前にある。
- comparison cardは3 action、1行summary、初期closedのSettingsを持つ。
- pair editor、Record Layout、runtime、cache、raw resultはprimary inputを押し下げない。
- active methodに関係しないcontrolはDOMとtab orderから外れる。
- DOMとtab orderはInput、Comparison、Basic、Generate、Advancedの順になる。

### Responsive layout

- 1280 x 720で最初のuploaderとheader `Add sequence`がスクロールなしで見える。
- 2件目追加後もcomparison detailは自動展開しない。
- 390 x 844で横overflowとfixed action barの重なりがない。

### Stateと実行

- fresh LinearとResetは`No comparison`になる。
- current Web sessionとsupported legacy Web sessionの保存済み比較intentを維持する。
- CLI-only replayはcommitted renderを維持し、次のWeb Generateを`No comparison`から始める。
- mode switchとaccordion closeはinactive draftを消さない。
- LinearのnoneはGenerate/Update中にLOSAT初期化、comparison file read、cache probe、active comparison resource生成を行わない。Save Sessionはdormant binding保存のためfileを読んでよい。
- Linearのnoneはhidden comparison fieldをvalidate/normalizeせず、invalidなdormant heightで作図を止めない。
- selected/mixed planでは3 bulk commandと`Current: Selected pairs` statusを区別して表示する。
- selected/mixed planでSimilarity groupsまたはCollinear blocksを暗黙に選択しない。

### Accessibilityとerror

- 3 action、summary、file input、pair actionをkeyboardだけで操作できる。
- focus-visible、command group名、button名、current status、pair countが確認できる。
- D7のstructured comparison issueはsummaryに出て、Generate時に該当sectionが開き対象controlへfocusが移る。

### Persistenceとsurface parity

- session round-tripはplan、program、presentation、filter、appearance、inactive file/raw filenameを維持する。
- disclosure open stateはsessionへ保存されない。
- canonical request schema 5、session version、Python API、CLI、renderer outputは変わらない。
- tracked reference SVGに差分がない。

### Documentation

- Web reference、FAQ、GUI tutorial、Gallery tutorial、capture selectorが新UIと一致する。
- affected screenshotは実行可能なcapture scriptから再生成され、desktop/mobileで目視確認される。

## 9. Riskと中止条件

| Risk | 防止策 | 中止条件 |
| --- | --- | --- |
| fresh default変更がlegacy restoreにも波及する | fresh factoryとmigration fallbackを別testで固定する | supported sessionのactive intentが変わる |
| UI projectionが第2のdomain modelになる | resolver outputと既存stateだけを入力にし、mutationとfile accessを禁止する | summary helperがedge topologyを再実装する |
| hidden controlがvalidationを止める | resolved intentとactive methodだけをvalidate/normalizeする | noneがdormant comparison valueで失敗する、またはdormant値がGenerateで変わる |
| mode切替でdraftを失う | object再作成を避け、既存activation flagを使う | file、raw filename、method固有値が往復で消える |
| noneでもruntimeが起動する | resolved Linear planを最初のpermission boundaryに置きcounter testを使う | Generate/Update中にcomparison file/cache/runtimeへの呼出しが1回以上ある |
| Circularへ見た目だけのparityを足す | 本計画のscope stopを守り別planへ分ける | Circular stateまたはmanaged slotに変更が必要になる |
| 同じfilterを二重表示する | visible ownerを一つにし旧markupを削除する | 同じfieldへ複数のactive `v-model`が残る |
| CSS順とDOM順が食い違う | markup自体を並べ替える | screen readerの順序が見た目と異なる |
| screenshotだけ更新し再現手順が壊れる | capture scriptを先に更新して再撮影する | manual cropしか再現方法がない |

次の回避策は採用しない。

- comparison panelをCSSで小さくするだけでDOM順を残す。
- accordion closeを`No comparison`と同義にする。
- 設定を隠すためにinactive draftを削除する。
- Generate時にfile有無からcomparison intentを推測する。
- viewport testを通すためにprimary controlをabsolute positionへ移す。
- UI cleanupを理由にrequest schema、session version、rendererへfieldを追加する。

## 10. 計画review記録

2026-08-10にread-onlyの独立reviewを2回行った。初回で見つかったblocker/high findingは、次の製品判断と境界追加で解消した。最終判定はready to executeである。

| Finding | 計画への反映 |
| --- | --- |
| Circularのinactive slotとuploaded/derived artifact ownerが未定義 | Circular parityをscope外とし、別計画で必要なowner判断をD5へ記録 |
| selected/mixedとSimilarity groups/Collinear blocksが非互換 | selectedでは3つのpairwise-compatible methodだけを許可し、adjacent化は明示bulk actionに限定 |
| CLI-only replayにeditable Web intentがない | committed renderを維持し、次のWeb Generateはfresh noneとする例外をD1へ記録 |
| Generate zero-workとSave Sessionのdormant file保存が混在 | Generate/Updateのpermission boundaryとsession resource bindingを別testへ分離 |
| selected状態のbutton semanticsとerror focusが曖昧 | bulk command/current statusを分離し、structured resolver issueだけをfocus routing対象に限定 |
| hidden Pairwise Match Heightがnoneを失敗させる | `run-analysis.js`と`session-request.js`をownerへ加え、active intentだけをvalidate/projectするtestを追加 |

## 11. 実装ledger

完了はtest、browser measurement、capture log、diffのいずれかで証明する。コードを変更しただけではcompletedにしない。

| Phase | Status | 完了証拠 |
| --- | --- | --- |
| Phase 0: baselineと境界 | pending | 実装開始時に記録する |
| Phase 1: 比較intent | pending | unit、session、runtime counter test |
| Phase 2: UI projection | pending | pure helper unit test |
| Phase 3: DOMとdisclosure | pending | static contractとbrowser DOM audit |
| Phase 4: state/session/request | pending | Node test結果 |
| Phase 5: browser UX | pending | 1280 x 720、390 x 844のPlaywright結果 |
| Phase 6: docsとGallery | pending | capture log、media check、目視記録 |
| Phase 7: final gate | pending | focused/broad gateとdiff audit |

## 12. 完了時のhandoff

最終報告には次を含める。

- 変更した製品契約と、変更しなかったsurface
- fresh、Reset、saved session、legacy sessionの比較intent結果
- Linear noneのGenerate/Updateで実行処理が0だった証拠と、Save Sessionでdormant bindingが残る証拠
- 1280 x 720と390 x 844のviewport結果
- focused test、broad test、Playwright、captureのcommandと結果
- reference SVG、request schema、session versionが不変だった確認
- production、test、documentation、generated artifactのdiff summary

Proposed commit title: `web: simplify LOSAT and comparison controls`
