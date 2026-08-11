# Work package G: Annotation UX finish implementation plan

- Date: 2026-08-11
- Status: ready to execute; implementation not started
- Baseline: `docs_renovation` / `1a1780d3a9a7` と文書作成時点のworking tree
- Source: [`gbdraw_v0.14.0_codex_roadmap.md`](./gbdraw_v0.14.0_codex_roadmap.md) の Work package G
- Target release: `v0.14.0`

この文書は、Work package Gを実行可能な作業単位へ分解した計画である。本作業が変更したのはroadmapのG関連箇所と本計画だけであり、production code、test、public documentation、Gallery media、documentation screenshot、reference outputは変更していない。

shared working treeには本計画と無関係なroadmap変更が存在する。本作業ではそれらを保持した。実装開始時にはPhase 0でbaselineと差分を取り直し、同じownerへ後から入った変更を優先して本計画を調整する。本文中の行番号ではなく、symbol、DOM label、capture ownerを再確認してから編集する。

## 1. 結論

`v0.14.0`で追加するのは、Region Annotationsの現在のeditor draftを既存TSV encoderで直接保存する`Download TSV`だけである。

- top rowは`[Add set] [Import TSV] [Download TSV]`にする。
- Generate前でもofflineでも`annotations.tsv`を保存できる。
- 保存対象は現在編集中のannotation rowsであり、最後にcommitされたResultではない。
- annotation rowが0件なら、ボタンを理由付きでdisabledにする。
- round-tripはTSVで表現可能な有効な非空rowのeffective semanticsに限定する。
- Web encoder、Web parser、Python `read_annotation_table()`が同じtable contractを共有する。
- current codecで確認できるexplicit `fill: null`のdriftを、ボタン公開前に修正する。
- toolbar/popupの`Annotate selection` shortcutは追加しない。既存のset内`Selected features` actionは維持する。
- Python API、CLI、canonical request、session、worker、renderer、SVG、reference outputは変更しない。

UI、codec修正、focused tests、実browser download、user documentation、該当する2種類のscreenshotを一つの機能単位として完成させる。buttonだけを先にlandしない。

## 2. 目標と対象範囲

### 2.1 必須目標

- Region Annotationsを開けば`Download TSV`が常に見つかる。
- zero setsとempty setsのどちらでもheader-only downloadを作らない。
- 最初のannotation rowを追加またはimportするとdownloadが有効になる。
- 1 activationで1 fileだけを保存する。
- 一つのBlob URLを作成し、同じactivation内で一度だけrevokeする。
- filenameを`annotations.tsv`、MIME typeを`text/tab-separated-values;charset=utf-8`に固定する。
- 既存のcolumn orderとterminal newlineを維持する。
- current editor draftを正規化して保存する。
- Webでdownloadした実ファイルを既存のImport TSVへ戻せる。
- 同じ実ファイルをPython readerで読める。
- multi-recordのrecord IDと`#N` selectorを維持する。
- coordinate target、feature target、4 marks、label、lane、effective legend/styleを維持する。
- explicit no-fillをWeb default fillへ変えない。
- downloadによってrender、analysis、history、Result、canonical owner、selection、viewportを変更しない。
- runtime networkを遮断しても同じ操作が完了する。

### 2.2 成果物

1. 既存annotation codecのbounded semantic round-trip修正。
2. Region Annotations top rowのdownload actionとempty-state説明。
3. Web codec unit、Python reader、CI収集対象のbrowser acceptance。
4. Web app reference、TSV schema reference、annotated chloroplast GUI tutorialの更新。
5. exact Gallery sessionからのRegion Annotations operation crop再撮影。
6. documentation capture ownerによる`T-GUI-05/03-annotation-table.png`再撮影。
7. focused、browser、offline、broad、no-render-diffのevidence ledger。

### 2.3 対象外

- toolbarまたはfeature popupの新しい`Annotate selection` shortcut
- stable feature hash selectorや新しいrecord-binding path
- grouped selected-feature spanという新しいdomain behaviour
- annotation model、target kind、mark、TSV column、table versionの追加
- empty annotation setのTSV保存
- set defaultとitem overrideの所有構造の復元
- set/item metadataのTSV保存
- cell内tab/newlineとsurrounding whitespaceのbyte-for-byte保存
- drag-to-resize endpoint
- point callout
- nested style editor
- GenBank/GFF3 biological featureへの変換
- Python public API、CLI、typed request、canonical request、session versionの変更
- worker、render pipeline、track layout、SVG semantic hook、reference SVGの変更
- Work package Eのannotation direct-preview allowlist拡張
- Gallery session、source file、example inventory、thumbnail、final previewの再生成

## 3. 現状の根拠

| 現状 | Owner | 実装への影響 |
|---|---|---|
| Region Annotations top rowには`Add set`と`Import TSV`だけがある | [`gbdraw/web/index.html`](../../gbdraw/web/index.html) の Region Annotations card | 同じrowに第三buttonと0-row説明を置く |
| annotation actionはfeature moduleへ集約される | [`gbdraw/web/js/app/annotations.js`](../../gbdraw/web/js/app/annotations.js) の `createAnnotationEditor()` | download可否とactionをここに置き、templateへserializerを持ち込まない |
| app surfaceはeditor methodを公開する | [`gbdraw/web/js/app/app-setup.js`](../../gbdraw/web/js/app/app-setup.js) のannotation bindings | download methodとrow availabilityだけを公開する |
| Web TSV parser/encoderは一つのmoduleにある | [`gbdraw/web/js/app/annotations/table-codec.js`](../../gbdraw/web/js/app/annotations/table-codec.js) の `parseAnnotationTable()` / `encodeAnnotationTable()` | 新しいdownload serializerを作らない |
| Generate時のCLI helperも同じencoderを使う | [`gbdraw/web/js/app/run-analysis.js`](../../gbdraw/web/js/app/run-analysis.js) の`web_annotations.tsv` staging | direct downloadとGenerate helperのbytesを分岐させない |
| local text downloadはBlob URLを同期的に作成・破棄する | [`gbdraw/web/js/services/text-download.js`](../../gbdraw/web/js/services/text-download.js) の `downloadTextFile()` | filename/MIME/1 click 1 downloadを既存ownerで満たす |
| encoderは`annotation.style || set.defaultStyle`とlegend fallbackを各rowへ展開する | `table-codec.js` の `encodeAnnotationTable()` | TSV contractをeffective valueに限定する |
| parserはrow styleをitem styleへ置き、set default/legendはWeb defaultへ戻る | `table-codec.js` の `styleFromRow()` / `parseAnnotationTable()` | effective outputは保てるがinheritance ownershipは保てない |
| explicit `fill: null`と他のstyle fieldをencodeするとfill cellがblankになる | `encodeAnnotationTable()` のnullish serialization | Web parserがblankを`#94a3b8`へ戻すdriftを修正する |
| empty setはdata rowを持たない | row-oriented table schema | total row count 0ではdownloadを無効にし、empty set preservationを主張しない |
| Web default fillは`#94a3b8` | [`gbdraw/web/js/app/annotations/state.js`](../../gbdraw/web/js/app/annotations/state.js) | blank fillの意味をdefault colorと混同しない |
| Python table schemaは既存Web columnsを受ける | [`gbdraw/annotations/io.py`](../../gbdraw/annotations/io.py) の `ANNOTATION_TABLE_COLUMNS` / `read_annotation_table()` | Python production codeやschemaを変えず、実downloadをreaderへ渡してparityを証明する |
| Web unitは基本的なcoordinate/feature round-tripと一部record selectorだけを見る | [`tests/web/annotations.test.mjs`](../../tests/web/annotations.test.mjs) | full semantic projectionとno-fill regressionを追加する |
| Python typed request testはannotation modelを深く検証するがWeb TSV bytesではない | [`tests/test_session_request_codec.py`](../../tests/test_session_request_codec.py) | request/session evidenceをtable evidenceの代用にしない |
| multi-record browser testはrecord selector UIを確認するがdownload/re-importしない | [`tests/web/linear-multi-record.playwright.spec.js`](../../tests/web/linear-multi-record.playwright.spec.js) | new browser acceptanceで`#2` download境界を完結させる |
| CIのWeb jobはNode unitと`pytest -m "browser and not slow"`を実行する | [`.github/workflows/test.yml`](../../.github/workflows/test.yml) | browser acceptanceをpytest browser markerで既存gateへ入れる。workflow変更はcollection不足時だけ行う |
| tobacco GalleryにはRegion Annotationsのoperation cropがある | [`gbdraw/web/gallery/tutorials/tobacco-chloroplast.json`](../../gbdraw/web/gallery/tutorials/tobacco-chloroplast.json) と`manual-06-01-region-annotations.webp` | exact session、state/visible contractを補強し、このcropだけを再撮影する |
| T-GUI-05はannotation import直後を撮る | [`docs/capture/flows/tutorials/gui_annotated_chloroplast.py`](../capture/flows/tutorials/gui_annotated_chloroplast.py) の`03-annotation-table.png` owner | Download visible/enabled assertionを追加し、scenario ownerで再撮影する |

## 4. 固定する製品判断

### G1. download ownerはcurrent editor draft

入力は`state.annotationSets`の現在値である。最後にGenerateしたResult、`committedCanonicalSession`、CLI helper file、session snapshotから逆算しない。

この区別により、ユーザーは未生成のannotation draftも保存できる。downloadはauthoring-data exportであり、rendered-result exportではない。

### G2. serializerは既存encoder一つ

`encodeAnnotationTable(state.annotationSets)`を直接使う。column list、escaping、record-selector formattingをbutton handlerやtemplateへ複製しない。

GenerateがCLI helperとしてstageする`annotations.tsv`も同じencoder consumerのままにする。codec修正が必要なら、direct downloadとGenerate helperの両方に同じ修正が届くことをunitで固定する。

### G3. row countがavailabilityを決める

availabilityはset数ではなく、正規化後の全setに含まれるannotation row数の合計で決める。

| Draft | Download state | 理由 |
|---|---|---|
| zero sets | disabled | tableにdata rowがない |
| one empty set | disabled | setだけを表すrowがschemaにない |
| empty setとnon-empty set | enabled | non-empty rowsを保存し、empty setは表現しない |
| one or more valid rows | enabled | bounded contractの対象がある |

disabled buttonの近くに`Add an annotation before downloading.`相当の短い説明を表示し、accessible relationを付ける。handler側にも0-row guardを置き、DOMから直接呼ばれてもheader-only fileを作らない。

### G4. filenameとbytesを固定する

- filename: `annotations.tsv`
- MIME: `text/tab-separated-values;charset=utf-8`
- header: existing `COLUMNS` order
- newline: LF-separated rows and one terminal LF
- cell normalization: encode時のtab/CR/LF→spaceとimport時のsurrounding-whitespace trim

timestamp、output prefix、mode、Result nameをfilenameへ混ぜない。

### G5. round-tripはeffective semantics

set ownershipを含むobject deep equalityではなく、各有効rowを次のprojectionへ正規化して比較する。

```text
normalizeCell(text)  = trim(replaceTabCrLfWithSpaces(text))
effectiveStyle(row)  = annotation.style ?? set.defaultStyle
effectiveLegend(row) = blankToNull(normalizeCell(annotation.legendLabel ?? set.legendLabel))
```

encoderがset defaultをitem rowへmaterializeし、parserがitem overrideとして復元しても、このprojectionが一致すればtable contractを満たす。再import後にset defaultを変更したときの追従関係はTSV contractに含めない。

### G6. blank fillはexplicit no-fill

encoderがfull effective styleをrowへ書いた状態で`fill`だけblankなら、そのrowのeffective fillは`null`である。Web parserは他のstyle cellがあることを理由にWeb default fillを注入しない。

column自体を省略した外部tableと、既存encoderが出したblank fillの扱いを混同しない。Phase 1でheader presenceとrow-style detectionをtest-firstで固定し、Python readerの既存`None` semanticsと一致させる。互換tableの解釈を広く変える必要が生じた場合はstop ruleを適用する。

### G7. record selectorはtable textがauthority

- unique record IDはそのIDを保存する。
- duplicate IDなどindex bindingが必要な場合はone-based `#N`を保存する。
- blank recordはsingle-recordで許される既存意味を維持する。
- `#0`、負数、不正tokenは既存validationどおりrejectする。

downloadのためにrecord catalogを再解決したり、IDへ強制変換したりしない。

### G8. downloadは副作用を持たない

buttonにはhistory captureを避ける既存markerを付ける。actionはencoderとdownload serviceだけを呼ぶ。

次を呼ばない、更新しない。

- Generate / render worker
- LOSAT / comparison planning / cache
- Result list / selected Result
- feature catalogue / editor binding
- committed canonical request/resources
- session writer
- Undo/Redo stack
- selected feature/annotation
- zoom / pan
- external fetch

### G9. shortcutはdeferred

既存`featureTargetsFromSelection()`は`locus_tag`、`gene`、bare `#N`の順でtargetを作る。duplicate qualifierは複数featureへ解決され、fallbackはstable feature identityではない。duplicate record IDもcatalogue bindingなしでは曖昧になる。

したがって、別の入口だけを追加しない。将来実装する場合はstable hash selector、record catalog binding、destination set UI、duplicate qualifier/record ID/missing qualifier testsを先に設計する。

### G10. persisted/rendering contractを変えない

annotation draft自体は既存session ownerで保存されるが、download actionのための新しいpersisted fieldはない。request/session versionを上げず、renderer input、SVG、reference outputを変えない。

## 5. Bounded semantic equivalence

### 5.1 比較projection

各setをfirst appearance order、各annotationをrow orderで比較する。JavaScriptのnumberとPythonのint/floatは値として比較し、言語固有classのdeep equalityを要求しない。

record selectorは両言語のinternal objectを、blank、record ID、またはone-based `#N` textへ正規化して比較する。feature selectorsもordered `(key, value)` pairsへ正規化する。labelなどのfree textは`normalizeCell()`後、optional legendのblankは`null`として比較する。

| Domain | 比較する値 |
|---|---|
| Row identity | set ID、annotation ID、set order、row order |
| Coordinate target | record selector、start、end、coordinate space、wraps origin、out-of-bounds policy |
| Feature target | record selector、ordered selectors、envelope、circular path |
| Mark/content | `highlight`、`band`、`line`、`bracket`、label |
| Placement | lane `null`、`0`、representative non-zero lane |
| Legend | effective legend labelまたは`null` |
| Stroke | color、width、dash array、line cap |
| Fill | colorまたは`null`、opacity |
| Hatch | `null`またはangle、spacing、color、width、cross |
| Label style | color、font sizeまたは`null`、orientation、position、offset |

### 5.2 必須fixture matrix

一つのcompact matrixで次を少なくとも一度ずつ通す。

- two non-empty setsとmultiple rowsによるset/row order、set/item IDs
- coordinate targetとfeature target
- blank record、record ID `RecB`、record index `#2`
- four marks
- labelあり/なし
- lane `null`、`0`、`2`
- set-level legend fallback、item-level legend override、blank→`null`
- set default style fallbackとitem style override
- stroke width、non-empty dash array、`none`/`tick`/`arrow` line cap
- colored fillとexplicit `fill: null`
- fill opacity `0`とnon-zero value
- hatch absent、single hatch、cross hatchと全hatch subfields
- label font size null/non-null、orientationとpositionの各enum値をrows全体で少なくとも一度
- coordinate wrap/out-of-boundsとfeature envelope/circular path

all enum valuesを一行へ無理に詰めず、4 marksを軸に複数rowへ割り当てる。失敗時にどのfieldがdriftしたか分かるsemantic projection helperをtest ownerへ置く。

### 5.3 明示的に比較しない値

- empty setの存在と順序
- valueがset default由来かitem override由来か
- set/item metadata
- tab、CR、LF、surrounding whitespaceを含むoriginal cell bytes
- object prototype、mapping proxy、JavaScript/Python class identity
- rendererが導出したlane geometryやSVG

これらが必要になった場合、既存tableの「bug fix」として暗黙に追加しない。schema/model extensionとして別計画へ戻す。

## 6. Ownerと変更ファイル

### 6.1 Production owner

| File | Planned change | 変更しないもの |
|---|---|---|
| [`gbdraw/web/js/app/annotations/table-codec.js`](../../gbdraw/web/js/app/annotations/table-codec.js) | explicit no-fillとoptional-column presenceをbounded contractどおりにする最小修正 | columns、target model、record-selector syntax、separate encoder |
| [`gbdraw/web/js/app/annotations.js`](../../gbdraw/web/js/app/annotations.js) | row availabilityとguarded `downloadAnnotationTable()`を追加し、existing encoder/downloaderを呼ぶ | import replace semantics、selected-feature target generation |
| [`gbdraw/web/js/app/app-setup.js`](../../gbdraw/web/js/app/app-setup.js) | editor actionとavailabilityをtemplate surfaceへ公開 | state owner、session/request projection、worker orchestration |
| [`gbdraw/web/index.html`](../../gbdraw/web/index.html) | top-row button、disabled reason、accessible/history-ignore attributes | annotation nested editor、new shortcut、new style controls |

[`gbdraw/web/js/services/text-download.js`](../../gbdraw/web/js/services/text-download.js) と [`gbdraw/web/js/app/run-analysis.js`](../../gbdraw/web/js/app/run-analysis.js) はreuse対象であり、原則変更しない。変更が必要になった場合は既存consumer全体のregression scopeをPhase 0へ戻して再評価する。

### 6.2 Test owner

| File | Planned evidence |
|---|---|
| `tests/fixtures/annotation_download_full.tsv` | Web encoderのexact bytesとPython readerが共有するbounded table oracle |
| [`tests/web/annotations.test.mjs`](../../tests/web/annotations.test.mjs) | semantic projection、fixture matrix、header/order/newline、explicit no-fill、0-row/header-only codec behaviour |
| [`tests/test_annotations.py`](../../tests/test_annotations.py) | same table columnsをPython readerで読み、typed semantic projectionを確認 |
| `tests/test_annotation_download.py` | new pytest/Playwright acceptance。disabled/enabled、real download、Web re-import、Python reader、offline/no-side-effectを一つのownerで確認 |
| [`tests/web/session-request.test.mjs`](../../tests/web/session-request.test.mjs) | annotation style/lane/record selectorの既存session projectionが不足するときだけdeep assertionを追加。production session ownerは変更しない |
| [`tests/test_documentation_scenario_contracts.py`](../../tests/test_documentation_scenario_contracts.py) | `test_t_gui_05_annotation_download_capture_contract`を追加し、flowにDownload visible/enabled assertionがあることを固定 |
| [`tests/test_gallery_capture_contracts.py`](../../tests/test_gallery_capture_contracts.py) | `test_tobacco_annotation_download_capture_has_exact_contract`を追加し、operationがexact session/state/visible contractを宣言することを固定 |

### 6.3 Public documentationとgenerated media owner

| File/artifact | Planned change |
|---|---|
| [`docs/REFERENCE/web-app.md`](../REFERENCE/web-app.md) | current draft download、filename、zero-row state、Generate不要を説明 |
| [`docs/REFERENCE/input-formats-and-tsv-schemas.md`](../REFERENCE/input-formats-and-tsv-schemas.md) | bounded effective-semantic contractと表現しない値を明記 |
| [`docs/TUTORIALS/GUI/build-an-annotated-chloroplast-map.md`](../TUTORIALS/GUI/build-an-annotated-chloroplast-map.md) | import後にcurrent tableを保存できる短い手順を追加 |
| [`gbdraw/web/gallery/tutorials/tobacco-chloroplast.json`](../../gbdraw/web/gallery/tutorials/tobacco-chloroplast.json) | operation textとcapture metadataを新UIへ同期 |
| [`docs/internal/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md`](./WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md) | Work package Gのrecapture decision/evidenceを追記 |
| `gbdraw/web/gallery/media/tobacco-chloroplast/manual-06-01-region-annotations.webp` | exact tobacco sessionから対象operationだけを再撮影 |
| [`docs/capture/flows/tutorials/gui_annotated_chloroplast.py`](../capture/flows/tutorials/gui_annotated_chloroplast.py) | `Download TSV` visible/enabled assertionを追加 |
| `docs/images/t-gui-05/03-annotation-table.png` | T-GUI-05 owner flowから再撮影 |

final Gallery preview、Gallery session/source/thumbnail、他のT-GUI-05 screenshots、reference SVGはno-change対象である。

## 7. 実装phase

### Phase 0 — Baselineとcontract oracle

Status: not started

Depends on: none

Owned systems: planning evidence only

#### Work

1. current HEAD、branch、working treeを記録する。
2. annotation codec、editor、template binding、Python reader、CI browser collectionを再確認する。
3. current request schema/session versionとreference-output baselineを記録する。
4. Section 5のfixture matrixとsemantic projectionをtest名へ割り当てる。
5. current expected failureとしてexplicit no-fill driftを再現し、他のbaseline failureと区別する。
6. Gallery/tutorial/capture inputsのhashとownerを記録し、対象mediaを2点に固定する。

#### Gate

```bash
git status --short --untracked-files=all
git rev-parse --short=12 HEAD
git diff --stat
rg -n '^CANONICAL_REQUEST_SCHEMA|^CURRENT_SESSION_VERSION' gbdraw/session_request_codec.py gbdraw/session_io.py
rg -n '^export const SESSION_VERSION' gbdraw/web/js/services/config.js
TMPDIR=/tmp node --test tests/web/annotations.test.mjs
TMPDIR=/tmp node --test tests/web/session-request.test.mjs
TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest tests/test_annotations.py tests/test_session_request_codec.py -q
TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest --collect-only -m "browser and not slow" -q
```

#### Evidence

- SHA、dirty paths、existing test counts
- Web/Python table columns
- request/session version
- explicit no-fill reproduction
- browser collectionにnew testを入れられる既存CI path
- before hashes for the two recapture artifacts

#### Completion conditions

- Section 5の各fieldにautomated oracleが割り当てられている。
- unrelated dirty diffとin-scope diffが分離されている。
- schema expansionなしでbounded contractを満たせる見通しがある。

### Phase 1 — Codec semantic round-trip

Status: not started

Depends on: Phase 0

Primary files:

- `gbdraw/web/js/app/annotations/table-codec.js`
- `tests/web/annotations.test.mjs`
- `tests/test_annotations.py`

#### Work

1. test側にlanguage-neutralなsemantic projectionを定義する。
2. Section 5.2のmatrixで`parseAnnotationTable(encodeAnnotationTable(sets))`を比較する。
3. exact header、column order、one terminal newlineを固定する。
4. record IDと`#2`を別rowで確認する。
5. explicit `fill: null`をred testにし、Web parserの最小修正でgreenにする。
6. optional style columnがない外部tableの既存解釈をregression testに残す。
7. Web encoderのexact outputをshared `tests/fixtures/annotation_download_full.tsv`へ固定する。
8. Python `read_annotation_table()`へ同じfixtureを渡し、typed projectionを比較する。
9. empty setがheader-onlyになる現状を記録するが、empty set round-tripは追加しない。

#### Gate

```bash
TMPDIR=/tmp node --test tests/web/annotations.test.mjs
TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest tests/test_annotations.py -q
```

#### No-change gate

- `ANNOTATION_TABLE_COLUMNS`は不変。
- `gbdraw/annotations/models.py`と`gbdraw/annotations/io.py`のproduction diffは0。
- request/session codecとversionは不変。

#### Completion conditions

- Web→Webとtable→Pythonの全bounded semanticsがgreenである。
- no-fill driftが再現test付きで解消している。
- inheritance ownershipとempty setをlosslessと主張していない。

### Phase 2 — Visible download actionとzero-row UX

Status: not started

Depends on: Phase 1

Primary files:

- `gbdraw/web/js/app/annotations.js`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/index.html`

#### Work

1. annotation editorにtotal row availabilityを一箇所で追加する。
2. guarded download actionから既存encoderと`downloadTextFile()`を呼ぶ。
3. app setupからtemplateへ必要なmethodだけを公開する。
4. top rowを3 actionsへし、`Download TSV`を常時表示する。
5. buttonへstableなaccessible nameを付け、0-rowではdisabled、短いreason、accessible relationを表示する。
6. download buttonをhistory capture対象外にする。
7. filename/MIMEを固定し、output prefixやResult stateへの依存を持たせない。
8. 既存Import TSVとSelected features actionを変更しない。

#### Gate

```bash
node --check gbdraw/web/js/app/annotations.js
node --check gbdraw/web/js/app/app-setup.js
node --check gbdraw/web/js/app/annotations/table-codec.js
TMPDIR=/tmp node --test tests/web/annotations.test.mjs
```

#### Completion conditions

- zero setsとempty setでguardが効く。
- 1 rowでenabledになる。
- action ownerがannotation editor一つである。
- serializer、Blob URL owner、history ownerを複製していない。
- Phase 1とPhase 2が同じimplementation changeとしてreview可能である。

### Phase 3 — Browser download、multi-record、cross-reader、no side effects

Status: not started

Depends on: Phase 2

Primary files:

- `tests/test_annotation_download.py`
- `tests/web/session-request.test.mjs` only if an existing persistence assertion is insufficient

#### Work

一つのPython Playwright acceptanceで次を通す。pytestの`browser` markerを付け、既存Web CI jobへ自動的に収集させる。

1. fresh appでRegion Annotationsを開き、zero setsのvisible disabled reasonを確認する。
2. empty setを作り、なおdisabledであることを確認する。
3. coordinate/feature、four marks、lane、legend、full style、explicit no-fillを持つdraftを作る。
4. duplicate record IDを含むmulti-record catalogでone-based `#2` bindingを作る。
5. Generate前にdownload eventを捕捉し、1 clickで1 fileだけ、filename/header/newlineを確認する。
6. `URL.createObjectURL`と`URL.revokeObjectURL`をspyし、前者へ渡されたBlobの`type`がexpected MIMEであることと、同じ一つのURLが一度ずつ作成・破棄されることを確認する。保存後のpathからMIMEを推測しない。
7. downloaded pathをPython `read_annotation_table()`へ渡し、Section 5 projectionを比較する。
8. same pathをImport TSV inputへ渡し、Web app stateのprojectionを比較する。
9. current annotation draft、history depth、Results/canonical owner、selection、viewportのbefore/after snapshotを比較する。
10. worker/render/LOSAT entry pointをspyまたはfail-fast stubにし、call count 0を確認する。
11. external networkをabortし、download/re-importが成功することを確認する。
12. 必要ならSave Session→fresh page→Load Sessionを一度通し、既存annotation projectionが変わらず再download可能であることを確認する。新しいsession fieldは追加しない。

#### Gate

```bash
TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest tests/test_annotation_download.py -m browser -q
TMPDIR=/tmp node --test tests/web/session-request.test.mjs
TMPDIR=/tmp npx playwright test tests/web/linear-multi-record.playwright.spec.js --project=chromium --workers=1 --grep "Region annotations expose"
```

Chromiumがagent sandboxで`Operation not permitted`になる場合は、同じcommandをsandbox escalationで再実行する。browser unavailableとしてskipしない。

#### CI gate

```bash
TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest --collect-only tests/test_annotation_download.py -m "browser and not slow" -q
```

既存`.github/workflows/test.yml`のWeb jobへcollectionされる限りworkflowを変更しない。collectionされない、またはrequired browser evidenceが実行されない場合だけ、同jobへfocused commandを追加する。

#### Completion conditions

- actual browser download→Web import→Python readerが同じfileでgreenである。
- `#2`とexplicit no-fillをactual file境界で確認している。
- zero analysis/render/history/Result mutationがcounterとsnapshotで証明されている。
- no unexpected skip/retryがない。

### Phase 4 — User documentationとtargeted recapture

Status: not started

Depends on: Phase 3

Primary owners:

- `docs/REFERENCE/web-app.md`
- `docs/REFERENCE/input-formats-and-tsv-schemas.md`
- `docs/TUTORIALS/GUI/build-an-annotated-chloroplast-map.md`
- `gbdraw/web/gallery/tutorials/tobacco-chloroplast.json`
- `docs/internal/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md`
- `docs/capture/flows/tutorials/gui_annotated_chloroplast.py`
- the two media outputs listed in Section 6.3

#### Work

1. Web referenceにcurrent-draft、Generate不要、filename、zero-row stateを記載する。
2. TSV schema referenceにeffective semanticsと非表現領域を記載する。
3. chloroplast GUI tutorialにimport後のdownload/re-import確認を短く追加する。
4. tobacco Gallery operationを`dataDependent: true`にし、exact tobacco sessionを宣言する。
5. capture metadataで`plastome_regions`と`lsc,irb,ssc,ira`をassertする。`visibleText`で`Import TSV`と`Download TSV`がcrop内にあることを確認し、capture actionはdownload buttonが存在して`disabled === false`でなければthrowする。現行`visibleControls` schemaがdisabledを検証できるとは仮定しない。
6. alt textとcaptionを新しいtop rowに同期する。
7. operation screenshot registerへtarget、required state、status/evidenceを追記する。
8. filename filterで`manual-06-01-region-annotations.webp`だけをrecaptureする。
9. T-GUI-05 capture flowでDownload visible/enabledをassertしてowner scenarioを再実行する。
10. GalleryとT-GUI-05のfocused source contract testsを追加する。
11. scenario実行前にT-GUI-05の兄弟4画像のSHA-256とbackupを`/tmp/gbdraw-wp-g-tgui05-before/`へ保存する。実行後にhashを比較し、差があれば受け入れず、原因を調べてreview済みbefore copyへ戻す。
12. 両画像をreadable scaleで目視し、button label、disabled stateではないこと、4 rows、crop切れを確認する。

#### Gallery gate

```bash
/home/kawato/micromamba/bin/python -m json.tool gbdraw/web/gallery/tutorials/tobacco-chloroplast.json
TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest tests/test_gallery_capture_contracts.py -k tobacco_annotation_download -q
TMPDIR=/tmp /home/kawato/micromamba/bin/python tools/capture_gallery_tutorial_screenshots.py --example tobacco-chloroplast --operation manual-06-01-region-annotations.webp
TMPDIR=/tmp /home/kawato/micromamba/bin/python tools/capture_gallery_tutorial_screenshots.py --example tobacco-chloroplast --check --strict
TMPDIR=/tmp npx playwright test tests/web/gallery-tutorial.playwright.spec.js --project=chromium --workers=1 --grep "tobacco chloroplast"
```

#### Documentation capture gate

```bash
TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest tests/test_documentation_scenario_contracts.py -k t_gui_05_annotation_download -q
mkdir -p /tmp/gbdraw-wp-g-tgui05-before
cp docs/images/t-gui-05/01-input-ready.png docs/images/t-gui-05/02-first-diagram.png docs/images/t-gui-05/04-track-settings.png docs/images/t-gui-05/05-finished-diagram.png /tmp/gbdraw-wp-g-tgui05-before/
sha256sum /tmp/gbdraw-wp-g-tgui05-before/01-input-ready.png /tmp/gbdraw-wp-g-tgui05-before/02-first-diagram.png /tmp/gbdraw-wp-g-tgui05-before/04-track-settings.png /tmp/gbdraw-wp-g-tgui05-before/05-finished-diagram.png
TMPDIR=/tmp /home/kawato/micromamba/bin/python docs/capture/run_all.py --scenario T-GUI-05 --tier extended
TMPDIR=/tmp /home/kawato/micromamba/bin/python docs/capture/run_all.py --scenario T-GUI-05 --tier extended --check
diff --brief /tmp/gbdraw-wp-g-tgui05-before/01-input-ready.png docs/images/t-gui-05/01-input-ready.png
diff --brief /tmp/gbdraw-wp-g-tgui05-before/02-first-diagram.png docs/images/t-gui-05/02-first-diagram.png
diff --brief /tmp/gbdraw-wp-g-tgui05-before/04-track-settings.png docs/images/t-gui-05/04-track-settings.png
diff --brief /tmp/gbdraw-wp-g-tgui05-before/05-finished-diagram.png docs/images/t-gui-05/05-finished-diagram.png
```

#### No-change gate

- `manual-08-01-chloroplast-preview.webp`は不変。
- tobacco Gallery session、source、download fixtures、thumbnail、examples inventoryは不変。
- T-GUI-05の`01`、`02`、`04`、`05`は不変。
- H-GUI-10はcontextual coverageだけであり、再撮影しない。
- reference SVGは不変。

#### Completion conditions

- public textとactual UIが同じcontractを説明する。
- capture metadataがexact data/state/control contractを持つ。
- intended two imagesだけが変わり、visual review evidenceがある。

### Phase 5 — Final verificationとdiff audit

Status: not started

Depends on: Phases 1–4

Owned systems: evidence and closure only

#### Focused gates

```bash
TMPDIR=/tmp node --test tests/web/annotations.test.mjs
TMPDIR=/tmp node --test tests/web/session-request.test.mjs
TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest tests/test_annotations.py tests/test_gallery_capture_contracts.py tests/test_documentation_scenario_contracts.py -q
```

#### Browser and offline gates

```bash
TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest tests/test_annotation_download.py -m browser -q
TMPDIR=/tmp npx playwright test tests/web/linear-multi-record.playwright.spec.js tests/web/gallery-tutorial.playwright.spec.js --project=chromium --workers=1
TMPDIR=/tmp /home/kawato/micromamba/bin/python tools/verify_gui_offline.py
```

#### Broad gates

```bash
TMPDIR=/tmp node --test tests/web/*.test.mjs
TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest tests/ -v -m "not slow"
TMPDIR=/tmp /home/kawato/micromamba/bin/python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
/home/kawato/micromamba/bin/ruff check gbdraw/
```

long-running pytestは少なくとも30分を許容し、60秒以内の間隔で進捗を確認する。unrelated baseline failureがある場合はexact failing test、baseline再現、in-scope影響をledgerへ分けて記録する。

#### Diff audit

production、test、documentation、generated artifactを別々にreviewする。

```bash
git diff --check
git diff -- gbdraw/web/index.html gbdraw/web/js/app/annotations.js gbdraw/web/js/app/app-setup.js gbdraw/web/js/app/annotations/table-codec.js
git diff -- tests/fixtures/annotation_download_full.tsv tests/web/annotations.test.mjs tests/test_annotations.py tests/test_annotation_download.py tests/test_gallery_capture_contracts.py tests/test_documentation_scenario_contracts.py
git diff -- docs/REFERENCE docs/TUTORIALS/GUI/build-an-annotated-chloroplast-map.md docs/internal/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md gbdraw/web/gallery/tutorials/tobacco-chloroplast.json docs/capture/flows/tutorials/gui_annotated_chloroplast.py
git diff --numstat -- gbdraw/web/gallery/media/tobacco-chloroplast/manual-06-01-region-annotations.webp docs/images/t-gui-05/03-annotation-table.png
git diff --exit-code -- tests/reference_outputs gbdraw/annotations gbdraw/api gbdraw/render gbdraw/diagrams gbdraw/session_io.py gbdraw/session_request_codec.py gbdraw/web/js/services/config.js gbdraw/web/js/services/session-request.js gbdraw/web/js/app/run-analysis.js gbdraw/web/js/services/text-download.js
git diff --exit-code -- gbdraw/web/gallery/sessions gbdraw/web/gallery/sources gbdraw/web/gallery/thumbnails gbdraw/web/gallery/examples.json gbdraw/web/gallery/media/tobacco-chloroplast/manual-08-01-chloroplast-preview.webp
git diff --exit-code -- docs/images/t-gui-05/01-input-ready.png docs/images/t-gui-05/02-first-diagram.png docs/images/t-gui-05/04-track-settings.png docs/images/t-gui-05/05-finished-diagram.png
```

Phase 0でこれらのpathに既存差分が記録されていた場合は、raw `--exit-code`の代わりにPhase 0のpatch/hashと実装後を比較する。既存差分を消してgateを通してはならない。

#### Completion conditions

- all mandatory gates pass with exact counts/environment recorded。
- no unexpected skip、retry-only pass、unexplained binary diffがない。
- production/test/docs/generated diffsがそれぞれownerと一致する。
- reference output、schema/version、API/CLI/session/rendererのdiffが0である。
- Section 10のcompletion ruleをすべて満たす。

## 8. Acceptance matrix

| ID | Requirement | Primary oracle | Evidence gate |
|---|---|---|---|
| G-UI-01 | top rowに3 actionsが常時見える | browser role/label assertion | Phase 3 |
| G-UI-02 | zero setsとempty setは理由付きdisabled | browser state + visible text | Phase 3 |
| G-UI-03 | first rowでenabled | browser state | Phase 3 |
| G-DL-01 | 1 activation = 1 `annotations.tsv` | Playwright download event | Phase 3 |
| G-DL-02 | header、order、terminal LF | downloaded bytes assertion | Phases 1, 3 |
| G-DL-03 | current draft、Generate前、offline | browser state/network fail-fast | Phase 3 |
| G-DL-04 | Blob URLを一度作成し一度破棄 | browser lifecycle spies | Phase 3 |
| G-DL-05 | expected MIMEでBlobを作る | `createObjectURL`へ渡された`blob.type` | Phase 3 |
| G-RT-01 | coordinate/feature target parity | JS/Python semantic projection | Phase 1 |
| G-RT-02 | record IDと`#2` parity | unit + actual download | Phases 1, 3 |
| G-RT-03 | 4 marks、labels、lanes parity | JS/Python semantic projection | Phase 1 |
| G-RT-04 | effective legend/style parity | JS/Python semantic projection | Phase 1 |
| G-RT-05 | explicit no-fill remains null | regression + actual file | Phases 1, 3 |
| G-RT-06 | browser download re-import | same downloaded path | Phase 3 |
| G-STATE-01 | no render/analysis/history/Result mutation | fail-fast counters + before/after snapshot | Phase 3 |
| G-COMPAT-01 | no model/schema/API/CLI/session change | no-change diff + existing tests | Phase 5 |
| G-RENDER-01 | no SVG/reference-output change | output comparison + diff | Phase 5 |
| G-DOC-01 | reference/tutorial text matches UI | doc review + contract tests | Phase 4 |
| G-VIS-01 | only owned operation/docs screenshots change | strict capture + visual/diff audit | Phase 4 |
| G-SCOPE-01 | global selection shortcut absent | DOM/code search + diff review | Phase 5 |

## 9. Risks、stop rules、atomic landing

### R1. “lossless”がobject ownershipまで拡大する

set inheritance source、empty set、metadataの復元が必須になったら停止する。既存row schemaでは表現できないため、table version/model extensionを別途設計する。本計画へhidden columnsやcomment recordsを足さない。

### R2. blank fillの外部互換性が曖昧

explicit no-fill修正が既存のpartial external tableを広く再解釈する場合は停止する。header presence、full Web encoder output、Python reader semanticsを分けたcompatibility ruleを先に追加する。

### R3. empty setを黙って失う

total rows 0では保存させない。mixed empty/non-empty draftではnon-empty rowsだけがtable対象であることをdocsへ明記する。empty set preservation用のsynthetic rowを作らない。

### R4. current draftとcommitted Resultが混ざる

actionがrun-analysis、committed session、Results file stagingへ依存し始めたら停止する。current authoring stateからexisting encoderへの一方向だけに戻す。

### R5. download clickがhistoryを増やす

history stackが変化したらrelease blockerとする。UI markerとhandlerの両方を確認し、history側へdownload-specific例外を増やさない。

### R6. duplicate serializer

button用のcolumn list、TSV escaping、Blob helperが新しく現れたらreviewを止める。既存ownerを再利用する。

### R7. shortcut scope creep

toolbar/popup action、destination modal、stable selector changeがdiffへ入ったら本packageから外す。既存set内`Selected features`のbugを発見した場合はissue/別planへ記録する。

### R8. screenshot refreshが広がる

Gallery final preview、thumbnail、session/source、他scenario、reference SVGへ差が出たらrecaptureを続けず原因を調べる。UI top rowが写る2 artifact以外は自動的に受け入れない。

### R9. browser evidenceがCIに入らない

local-only JS Playwright passだけでcompleteにしない。pytest browser collectionへ入れるか、focused commandを既存Web CI jobへ明示追加する。

### Atomic landing rule

次を一つのreviewable changeとして揃える。

1. codec regression testとminimal fix
2. visible actionと0-row guard
3. actual browser download/cross-reader evidence
4. public docsとowned recaptures
5. no-change/output evidence

production buttonだけ、またはbinary screenshotだけを先にlandしない。実装中にphase単位のcommitを使う場合も、release branchへ統合する時点では全gateを満たす一つの機能境界にする。

## 10. Completion rule

Work package Gをcompleteにする条件は次のすべてである。

1. `Download TSV`がRegion Annotations top rowにあり、keyboardで操作できる。
2. zero rowsが理由付きdisabledで、one rowがenabledである。
3. current draftをGenerate前/offlineに`annotations.tsv`として一度だけ保存できる。
4. bounded semantic matrixがWeb encode→Web importとdownload bytes→Python readerで一致する。
5. `#2`、four marks、lanes、effective style/legend、explicit no-fillのevidenceがある。
6. download→Import TSVがactual browser fileで一致する。
7. render、analysis、history、Result、canonical owner、selection、viewportが変わらない。
8. selection shortcut、new schema/model/API/CLI/session/renderer pathを追加していない。
9. focused Node/Python、CI-owned browser、offline、broad non-slow、output comparisonがpassする。
10. exact Gallery operationとT-GUI-05 annotation screenshotだけが更新され、strict validationとvisual reviewを通る。
11. production/test/documentation/generated diffsが説明され、no-change ownersに差がない。
12. evidence ledgerが実数、environment、deviation、remaining riskで埋まっている。

production codeが動くだけ、unit testだけ、またはscreenshotが新しいだけではcompleteにしない。全条件成立後にだけ文書先頭のStatusを`completed`へ変える。

## 11. Evidence ledger

gateを実行するまでは`complete`と記録しない。計画作成時点の初期値はすべて`not started`である。

| Phase | Status | Behavior implemented | Evidence | Deviations | Remaining risk |
|---|---|---|---|---|---|
| Phase 0 — Baseline/contract | not started | — | — | — | contract drift at implementation start |
| Phase 1 — Codec parity | not started | — | — | — | blank-fill compatibility |
| Phase 2 — Download UI | not started | — | — | — | history or duplicate-owner coupling |
| Phase 3 — Browser/cross-reader | not started | — | — | — | CI collection and browser lifecycle |
| Phase 4 — Docs/captures | not started | — | — | — | unintended recapture spread |
| Phase 5 — Closure | not started | — | — | — | broad-suite baseline failures |

各entryは少なくとも次を含める。

```text
Date/time:
Commit/worktree:
Environment:
Command:
Exit status:
Passed/failed/skipped counts:
Fixture or artifact:
Semantic/visual oracle:
Observed result:
Deviation or retry:
Remaining risk:
```

## 12. Cross-package handoff

- C/D/E/Fのrenderer/session implementationはTSV downloadのruntime dependencyではない。Phase 0でcurrent annotation ownerを再確認した後、existing editor draftとcodecだけで実装する。
- roadmapのlanding orderどおりstate contractsが安定した後にGを統合し、同じ`index.html`/`app-setup.js`周辺のmerge churnを避ける。
- Work package Eのannotation direct-preview allowlistを広げない。downloadはrender effectではない。
- shortcut deferにより、Work package C/Dのstable feature identityやrecord provenanceをGへ持ち込まない。
- Work package A1はfinal candidateでpublic docsとcapturesを再検証するが、G-owned stale operation screenshotの修正をA1まで延期しない。
- Work package JへAcceptance matrix、exact commands、CI job、two changed image hashes、no-change evidenceを渡す。
- implementation handoffでは英語のcommit titleと短いsummaryを提案する。推奨titleは`feat(web): download annotation tables from the editor`である。
