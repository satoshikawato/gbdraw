# Issue #563 Feature popupからのrecord回転 — 総合実装計画書

- 状態: 実装前（計画と実行境界のみ確定）
- 作成日: 2026-09-22
- 対象Issue: [#563 Add feature-based rotation of circular records from the feature popup](https://github.com/satoshikawato/gbdraw/issues/563)
- 必須実装ブランチ: `issue-563-feature-popup-record-rotation-20260922`
- 作成時base: `origin/dev` @ `a9eaeadd105e0c26e46086626feaa695bdd33c94`

## 0. この文書の使い方

本書は、過去の会話を一切知らない実装者がIssue #563の目的、現状、設計、作業順、
検証、停止条件を復元するための唯一の総合計画書である。実装セッションごとの実行指示は
同じディレクトリの`SESSION_01_*.md`から`SESSION_07_*.md`に分けてある。

すべての実装セッションは、必ず次のブランチを使用する。

```text
issue-563-feature-popup-record-rotation-20260922
```

別のsession branchを作らない。`main`、`dev`、他Issueのbranchへ実装しない。開始時に
`git branch --show-current`、`git status --short --branch`、`git rev-parse HEAD`、upstreamを
確認する。このbranchは最新`origin/dev`から作成済みで、`dev`をupstreamにしない。
後から`origin/dev`が進んでも、明示的な依頼なしにrebaseやmergeをしない。

## 1. Product goal

Feature popupで開いている生物学的featureを基準に、そのfeatureが属する一つのcircular
recordだけを非破壊で回転し、必要ならrecord全体の表示方向を変更して再描画できるようにする。

利用者はpopup内の独立した`Record actions`から次を指定する。

- Anchor: `5′ end`、`Midpoint`、`3′ end`
- 符号付き整数offset（bp）
- 任意の`Orient this feature forward`
- 独立したpreset `Place this feature at the end`
- 適用前のrecord、feature、解決済みsource coordinate、orientation変化のpreview
- `Apply and regenerate`または`Cancel`

対象はdiagram modeではなくrecord topologyで判定する。Circular modeのcircular recordにも、
Linear modeでlinearized表示されているcircular recordにも適用できる。source sequence、annotation、
qualifier、biological identityは変更しない。

## 2. 対象外

次はIssue #563に含めない。

- featureの自動選択、複数recordのbatch回転
- similarity group memberへの自動伝播やalignment
- crop-and-rotate
- re-originated GenBank/FASTAの書き出し
- CLI向けfeature検索UIまたは新しいCLI option
- 新しいrenderer、座標変換engine、Worker protocol、SVG座標からの逆算
- popup actionと無関係なフォーム変更の暗黙的な確定

## 3. 現行実装の確認済み事実

### 3.1 再利用できるbackend機構

- `gbdraw/api/requests.py`の`RecordPresentation.reverse_complement`がrecord表示方向を持つ。
- 同ファイルの`RecordDisplayOptions.start_coordinate`が1-based source display originを持つ。
- `gbdraw/layout/record_coordinates.py`の`RecordDisplayTransform`がsource/display座標変換を一元化する。
- feature、label、tick、depth、GC、comparisonなどは既に同transformを利用する。
- したがってIssue #563用の回転engineやSVG geometry計算は不要である。

### 3.2 Web stateとrequest

- `gbdraw/web/js/app/record-display-options.js`がsource-boundなrecord表示draftを所有する。
- 現在のdraftは`scope/sourceUid/selector/recordId/topologyOverride/startCoordinate`を持つ。
- 既存shortcutはglobal selectionに依存し、5′とmidpointだけを解決する。
- `gbdraw/web/js/services/session-request.js`がWeb stateからcanonical request schema 7を作る唯一のownerである。
- `cardinality: all`のrecordは、record別display差があるとexact-one recordsへmaterializeされる。
- Session current versionは43、feature catalog current schemaは3である。
- popupのtargetは`state.clickedFeature`であり、`clickedFeature.feat`にsource-bound feature metadataを持つ。

### 3.3 Generate、Result、History

- `gbdraw/web/js/app/run-analysis.js`がcanonical request構築、comparison準備、Worker実行、
  result validation/admission、current Result置換を調停する。
- `gbdraw/web/js/services/history.js`の`runUndoableArtifactReplacement()`が生成済みartifactの
  atomicな置換とrollbackを所有する。
- current Resultは成功したcandidateだけで置換され、failed/canceled/stale candidateはadmitされない。

### 3.4 LOSATについての重要な現状認識

通常の`Generate Diagram`でも、display startまたはcomplete-record reverse complementだけを
変更した場合、互換なraw LOSAT evidenceは既に再利用され、LOSAT executorは再実行されない。
根拠は`tests/web/linear-multi-record.playwright.spec.js`の
`LOSATP source jobs are reused after display start, reverse complement, and fresh Load`である。

したがって、本計画のcommitted-base candidate経路を設ける主目的は「LOSATを止めること」ではない。
主目的は次の二つである。

1. popup actionと無関係な未適用フォーム変更をcanonical requestへ混入させない。
2. target recordのorigin/orientation変更と成功Resultを一つのUndo/Redo単位にする。

candidateは既存comparison cacheを使い、追加LOSAT executor jobが0であることを回帰testで固定する。
comparison planning/cache lookupを通ること自体は許容する。性能上の追加最適化は計測で必要性が
示されない限り行わない。

## 4. Product semantics

### 4.1 Identityとeligibility

targetはaccession文字列や現在のSVG fragment IDではなく、次の組で解決する。

```text
(source-bound record identity, biologicalFeatureId)
```

少なくとも`recordKey`、`sourceUid`、exact record selector、resource bindingを照合する。
同名recordを含む複数file、1 file内の複数record、同一rowの複数record、seamで分割された
featureを区別する。

actionを有効にする条件:

- target recordのeffective topologyがcircular
- complete source recordと正のrecord lengthが利用可能
- recordがcroppedでない
- popup featureを現在のsource resourceとrecordへ一意に再bindできる
- 選択したoperationに必要なlocation precision、part order、strandが確定している

source replacement、古いResult由来のpopup、resource mismatchはstaleとしてdisableする。

### 4.2 Coordinate conventions

- UIとcanonical requestのdisplay startは1-based source base coordinate。
- feature partは内部で0-based half-open interval `[start, end)`。
- 計算はsource feature metadataだけを使用する。現在のSVG位置や表示fragmentを使用しない。
- circular wrapは非負moduloで行う。
- 入力offsetはsafe integerだけを受け付け、小数を切り捨てない。

record lengthを`L`、1-based anchorを`a`、strand方向を`+1/-1`、offsetを`o`とすると、
stranded featureのdisplay startは次である。

```text
1 + mod(a - 1 + strandDirection * o, L)
```

unstranded featureのoffsetはsource coordinate正方向を基準とし、`1 + mod(a - 1 + o, L)`とする。
その場合はupstream/downstreamという説明を表示しない。

### 4.3 Anchor

source location partsを確定したtraversal orderで走査する。

- 5′: biological traversalの最初のincluded base
- midpoint: covered basesだけのoffset `floor((N - 1) / 2)`
- 3′: biological traversalの最後のincluded base

midpointはmin/maxの平均ではない。part間gapを数えず、偶数長ではtraversal上の早い中央baseを
選ぶ。known consistent minus strandでは各part内も逆向きに走査する。

unstranded featureではbiological 5′/3′をdisableする。単一part、またはsource metadataが
順序を明示できるexact pathではmidpointを許可する。mixed-strand、fuzzy、unordered、
意味が失われたcompound locationでは安全に解決できないoperationだけを理由付きでdisableする。

### 4.4 Orientation

`Orient this feature forward`はtoggleではなくabsolute requestである。

- OFF: 現在のeffective record orientationを保持する。
- ONかつsource feature strandが`+`: target reverse complementは`false`。
- ONかつsource feature strandが`-`: target reverse complementは`true`。
- unstrandedまたはambiguous: ONをdisableする。

同じrequestを繰り返してもorientationは蓄積せず、同じabsolute stateへ収束する。
offsetで解決したsource coordinateは、その後に選ぶorientationから独立する。

### 4.5 Place this feature at the end

これは3′ anchorのaliasではない。最終orientationを先に解決し、そのdisplay traversal方向で
最後に表示されるincluded baseの直後をoriginにする。

display source stepを`d`（forward=`+1`、reverse complement=`-1`）、最後のincluded source
baseを`b`とすると、originは`1 + mod(b - 1 + d, L)`である。known-strand featureでは、
display directionとbiological traversalが一致すると3′ base、逆なら5′ baseが`b`になる。
明確なsource-forward pathを持つunstranded featureではdisplay directionに従ってpath端を選ぶ。

preset選択時はoffsetを0へ戻す。その後offsetを変更した場合はcustom placementとして表示する。
outgoing boundaryを一意に決められないlocationではpresetをdisableする。

## 5. Target architecture

```text
feature popup (index.html)
  -> record action controller (explicit clickedFeature target)
       -> pure feature anchor resolver
            -> source feature anchor profile + location parts
       -> target-only committed-request projector
            -> last successful canonical request clone
            -> one record's display.startCoordinate / presentation.reverseComplement
       -> existing canonical candidate execution in run-analysis.js
            -> existing comparison cache planning
            -> existing Worker / typed response / SVG admission
       -> extended artifact replacement transaction
            -> target record intent before/after
            -> generated artifact before/after
  -> existing record display controls reflect committed effective state
  -> Session save/load preserves absolute transform and non-authoritative provenance
```

### 5.1 Owner matrix

| 責務 | Semantic owner | 方針 |
| --- | --- | --- |
| source locationのprecision/order facts | Python feature metadata/catalog | Biopython locationを失う前に最小のanchor profileへ正規化する。anchor座標自体は決めない。 |
| anchor、offset、end placement、orientation解決 | `gbdraw/web/js/app/record-display/feature-anchor.js`（新規private domain module） | DOM、Vue、History、Workerへ依存しないpure function。sidebarとpopupが共有する。 |
| per-record display draftとstale binding | `record-display-options.js` | topology/start/orientation/provenanceの一つのeffective-state ruleを持つ。 |
| canonical request projection | `session-request.js` | schema 7の唯一ownerを維持し、record materializationとtarget-only overlayをここだけで行う。 |
| candidate execution/admission | `run-analysis.js` | normal Generateとpopup candidateが同じWorker/result admission helperを使う。 |
| atomic history | `services/history.js` | artifact replacement transactionへ小さいintent checkpoint hookを追加する。別transaction engineを作らない。 |
| popup presentation | `index.html` + focused controller | input、preview、disable reason、Apply/Cancelのみ。domain式をtemplateへ書かない。 |
| backend geometry | 既存request/render/`RecordDisplayTransform` | 変更しない。 |

### 5.2 State model

`recordDisplayDrafts`をrecord単位表示intentの唯一のWeb ownerとして拡張する。目標形は次の意味を持つ。

```js
{
  scope,
  sourceUid,
  selector,
  recordId,
  topologyOverride,
  startCoordinate,
  reverseComplementOverride,
  anchorIntent
}
```

- `reverseComplementOverride`: complete recordのabsolute override。`null`は既存base値を継承する。
- cropped recordのreverseは引き続きregion bindingが所有し、popup action対象外。
- `anchorIntent`: popupを再表示するためのprovenanceであり、rendering authorityではない。
- rendering authorityは常にcanonical request内のabsolute `startCoordinate`と
  `reverseComplement`だけである。
- manual start/orientation変更は古い`anchorIntent`をclearする。
- provenanceから自動的に再計算してResultを変えない。source mismatch時はstaleと表示する。

`anchorIntent`の最小内容:

```js
{
  schema: 1,
  recordKey,
  biologicalFeatureId,
  placement: 'anchor' | 'feature-end',
  anchor: 'five-prime' | 'midpoint' | 'three-prime' | null,
  offsetBp,
  orientForward
}
```

resource digestやresolved coordinateを二重保存しない。source bindingとabsolute display stateは既存
canonical/session ownerから取得する。

### 5.3 Feature anchor profileとschema

現在のcatalog schema 3はexact/fuzzy、compound operator、順序の確実性を完全には保持しない。
安全に「推測しない」を実現するため、Python側で小さなanchor capability profileを生成し、
feature catalog schemaを4へ上げる。profileは少なくとも次を区別できればよい。

- exactかfuzzyか
- single/join/order/unknown
- part orderが`biological`、`source-forward`、`ambiguous`のどれか
- consistent strandが`+`、`-`、unstranded、mixedのどれか

profileはcapability factsのみを持ち、5′/midpoint/3′の座標やUI messageは持たない。
その判断はJS pure resolver一箇所へ残す。

catalog 4をSessionで保存するためSession versionは44へ上げる。計画作成時点の`main` current writerは
Session 42であり、Session 43は`origin/dev`だけの未release formatなので、v43専用readerは追加しない。
branch-owned v43 artifactが存在する場合はv44へ書き換え、mainに存在するv42/catalog 3のreaderが
保存済みResultを壊さず読み込む。既存partsから安全に証明できるsingle exact locationだけを
限定的に復元し、それ以外のpopup rotationは`Generate again to refresh feature location metadata`
相当の明示理由でdisableする。legacy dataを推測して完全対応に見せない。

canonical render requestは既存fieldだけで完全に表せるためschema 7のままとする。Worker protocol、
Python request codec、CLI schemaをIssue #563のために増やさない。

### 5.4 Exact-one materialization

一つのsource requestが`cardinality: all`で複数recordsを表す場合、targetだけのstartまたは
orientation差を表現するにはexact-one recordsへ展開する。現行のdisplay差materializationを拡張し、
次のどちらかがrecord間で異なるときに同じ一箇所で展開する。

- `record.display`
- complete-record `record.presentation.reverseComplement`

record order、grid row、tracks、depth source index、comparison endpoint、resource bindingを維持する。
target以外のrecord requestはdeep equalityで不変であることをtestする。

### 5.5 Committed-base candidate

popup Applyは現在のform全体からrequestを再構築しない。最後に成功したcanonical requestと
committed resource bindingsをbaseとして、target recordの次の二fieldだけをoverlayする。

```text
records[target].display.startCoordinate
records[target].presentation.reverseComplement
```

必要なexact-one materialization、validation、clone規則は`session-request.js`が所有する。
popup/controllerがJSONを直接編集しない。base requestのrecord順、他record、diagram options、
tracks、comparisons、resource IDsを保持する。

このcandidateを、normal Generateから抽出した共通のcanonical candidate execution/admission helperへ
渡す。新しいWorker constructor、SVG insertion、result validator、sanitizer、error fallbackを作らない。

### 5.6 Atomic transaction

Applyの状態遷移:

```text
validate fresh target + resolve preview
  -> capture before artifact + target draft intent
  -> build committed-base candidate
  -> execute/admit candidate
     success: commit target draft + candidate artifact as one history entry
     fail/cancel/stale/superseded: restore before artifact + before target draft
```

`runUndoableArtifactReplacement()`へ任意の`captureIntent`/`restoreIntent`相当の小さいhookを追加し、
既存artifact handleと同じentryでbefore/afterを保持する。通常Generateの呼び出しはhookを渡さず、
現行挙動を保つ。

一回のUndoでorigin、orientation、anchor provenance、Resultがすべてbeforeへ戻る。一回のRedoで
afterへ戻る。失敗時にtarget draftだけが先行して残らない。unrelated pending form stateはbeforeでも
afterでも変更しない。

## 6. UI workflow

popupの`Edit` surfaceに、feature color/placementと区別した`Record actions` sectionを置く。
rich/simple popupのどちらでも到達可能にし、popupのdrag/resize操作と競合しないよう`data-nodrag`
boundaryを使う。

初期値:

- Anchor: 5′ end
- Offset: 0
- Orient: OFF（現在のeffective orientation保持）

常に表示する情報:

- target record labelとsource-bound identityを人間向けに要約した値
- target feature label
- resolved 1-based source coordinate
- record orientation: unchanged / forward / reverse-complemented
- displayed feature strand before -> after（解決できる場合）
- `Coordinates refer to the original record.`

invalid/disabled時は、単なるdisabled buttonにせず具体的な理由を表示する。Cancelはpopup actionの
ephemeral inputだけを捨て、record draft、Result、Historyを変更しない。

成功後:

- feature search queryを保持する。
- `(recordKey, biologicalFeatureId)`でtargetを再同定する。
- popupを同featureへ再bindするか、少なくとも検索結果でtargetを発見可能にする。
- zoom/panは新Resultの表示安全性が保てる範囲で維持する。既存preview lifecycleを迂回しない。

## 7. SOLID / KISS / DRY / YAGNI

### SOLID

- Single Responsibility: source facts、anchor resolution、request projection、execution、history、UIを分ける。
- Open/Closed: 既存typed request/transformへabsolute valuesを渡し、rendererを分岐させない。
- Liskov Substitution: 新しいclass hierarchyを作らないため適用対象を増やさない。
- Interface Segregation: resolverはrecord length、parts/profile、current orientation、intentだけを受ける。
- Dependency Inversion: UIはDOM/SVG座標ではなくcontrollerとpure resolverへ依存する。

### KISS

- 回転は`startCoordinate`、方向は`reverseComplement`という既存表現を使う。
- targetはpopup feature一件、record一件に限定する。
- legacy metadataは推測せず、再Generate案内で安全側に倒す。
- cache preparationを通ることは許容し、executor 0 jobという意味のある契約だけを固定する。

### DRY

- sidebar shortcutとpopupは同じresolverを使う。
- normal Generateとpopup candidateは同じWorker/result admissionを使う。
- record effective orientation、stale binding、materializationを各一ownerで決める。
- Undo/Redoは既存artifact transactionを拡張し、二つ目のtransaction managerを作らない。

### YAGNI

- generic mutation framework、command bus、plugin systemを追加しない。
- batch、automatic gene lookup、similarity propagation、source exportを先回りしない。
- performance最適化は追加LOSAT executor jobまたは計測済みbottleneckがない限り行わない。
- future schema用の未使用fieldやaliasを追加しない。

## 8. Product Impact / architecture preflight

Issue #563はproduct outcomeを詳細に指定しているが、GitHub Issue本文だけでbase-branchのdurable
authority要件を満たすとは自動的に仮定しない。Session 01で
`docs/internal/PRODUCT_IMPACT_RATCHET.md`に従いpreflightする。

- 既存のmerged authorityまたはmaintainer instructionが完全なoutcomeを選択している場合は
  `IMPLEMENT_EXISTING_AUTHORITY`として根拠を記録する。
- base-branch authority receiptが必要な場合はauthority-only changeを先に準備し、runtime実装を
  stopする。実装者が別のproduct option、rationale、risk acceptanceを補完しない。
- Issue本文にないuser-visible choiceが見つかった場合はDecision Packetを作り、affected convergenceを
  stopする。独断で選ばない。

Architecture ratchetでは、次のowner/pathを増やさないことを基本判定とする。

- canonical request owner: `session-request.js`のまま1
- candidate execution/result admission path: `run-analysis.js`の既存pathを共有して1
- SVG admission: 既存ownerのまま1
- history transaction owner: `services/history.js`のまま1
- record coordinate transform owner: Pythonの既存`RecordDisplayTransform`のまま1

Session 44がmainに存在するv42/catalog 3 readerを新しいcompatibility pathにするため、Session 01で
policyに従うbefore/after CB evidence、namespace、positive fixture、retirement conditionを記録する。
Session 43が実装開始前にfirst-parent `main`またはrelease tagへ到達した場合だけ、その時点の証拠に
基づいてv43を同じcompatibility namespaceへ追加する。schema bumpを避けるために
曖昧なmetadataを推測したり、schema 3へ意味を黙って追加したりしない。

## 9. Work packages and session order

| Session | 成果 | 依存 |
| --- | --- | --- |
| 01 | baseline、authority/preflight、architecture ledger、red contract tests | なし |
| 02 | Python anchor profile、catalog 4、pure shared resolver | 01 |
| 03 | per-record orientation/provenance、request materialization、Session 44 compatibility | 02 |
| 04 | committed-base target projector、shared candidate execution、atomic history | 03 |
| 05 | popup Record actions controller/UI、stale handling、interaction continuity | 04 |
| 06 | end-to-end browser acceptance、LOSAT 0-job、tracks/comparison/identity isolation | 05 |
| 07 | full gates、docs、diff/architecture/product review、handoff | 06 |

順番を入れ替えない。後続sessionは前sessionのproduction/tests/docs差分を確認し、未完了なら新しい
parallel implementationを足さずに不足を補う。

## 10. Test strategy

### 10.1 Pure JS resolver

- plus/minus: 5′、midpoint、3′、positive/negative offset、origin wrap
- even/odd covered length、multipart gap exclusion、origin-spanning ordered parts
- orientation absolute/idempotent、preserve-current
- 3′ anchorとfeature-endの差、forward/reverse display
- unstranded single/ordered path、mixed strand、fuzzy、unordered、invalid length/input
- input non-mutationとdeterministic output

### 10.2 Python/catalog contract

- simple、compound join、origin-spanning、minus、unstranded、mixed、fuzzy、order/unknown
- biological/source traversal profileがsource transform後も安定
- catalog schema 4 validation/compaction/admission
- biologicalFeatureIdとrecordKeyが変わらない

### 10.3 State/request/session/history

- one target onlyのdraft patch、manual editでprovenance clear
- `cardinality: all`から必要時だけexact-one materialization
- duplicate record IDs、same file multi-record、same row multi-record
- target以外のrequest subtree、tracks、comparison、resource IDsが不変
- v44 save/fresh Load/reopen、released v42/catalog 3 safe migration、future/invalid rejection
- Apply一件、Undo一件、Redo一件
- render failure、cancel、stale/superseded completionでartifactとintentをrollback
- unrelated pending form editsがcandidateへ入らず、Apply後もpendingのまま

### 10.4 Browser journeys

Issueの10 acceptance scenariosを最低一つのpointer journeyで覆い、keyboard/mobile/disabled reasonを
補助testで確認する。内部controllerを直接呼ぶtestだけをacceptanceにしない。

特に次を固定する。

- search -> popup -> rotate -> regenerated Result（sidebar selection不要）
- Circular/Linearで同じsource coordinate
- chromosome I/IIを独立に操作
- plus/minusの5′ -100 bp
- initial forward/reverseと再適用idempotence
- 3′とfeature-endの差
- compound/origin-spanning/even/odd midpoint
- unstranded/ambiguous/linear/cropped/stale/replaced source
- duplicate IDs/split fragments
- tracks、labels、ticks、comparison geometryの一貫変換
- LOSAT executor追加job 0
- search query、stable target identity、実用的なzoom/pan continuity

### 10.5 Gates

各sessionはfocused testsを実行し、Session 07で少なくとも次を実行する。実際のfile名やgrepは実装時に
追加されたtestへ合わせる。

```bash
python -m pytest tests/test_web_feature_metadata.py tests/test_web_feature_catalog.py -v
node --test tests/web/record-display-options.test.mjs
node --test tests/web/feature-anchor.test.mjs
node --test tests/web/session-request.test.mjs
node --test tests/web/history.test.mjs tests/web/run-analysis-simple-path.test.mjs
npx playwright test tests/web/linear-multi-record.playwright.spec.js --grep "record rotation|LOSATP source jobs" --workers=1 --retries=0
npx playwright test tests/web/interactive-svg-v3.playwright.spec.js --grep "record rotation" --workers=1 --retries=0
node --test tests/web/*.test.mjs
python -m pytest tests/ -v -m "not slow"
ruff check gbdraw/
node tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs --base origin/dev --head HEAD
git diff --check
```

Node Playwrightがない場合はrepository guidanceに従ってPython Playwrightで同じ境界を確認する。
Chromium sandbox errorは同じcheckを必要な権限で再実行し、未検証扱いで終えない。browser wheelが
必要なら`python tools/prepare_browser_wheel.py`で生成するが、gitignored wheelをcommitしない。

## 11. Acceptance matrix

| ID | 完了条件 |
| --- | --- |
| AC-01 | popupだけでtarget featureのrecordを回転でき、別selectionを参照しない。 |
| AC-02 | effective circular recordはCircular/Linear両modeで同じsource requestを解決する。 |
| AC-03 | same-file multi-recordを含め、target以外のrecord状態とlayoutは変わらない。 |
| AC-04 | plus/minusのoffsetがfeature方向基準で解決され、wrapする。 |
| AC-05 | orientation requestがabsoluteかつidempotentで、OFFでは現在値を保持する。 |
| AC-06 | 3′ base anchorとfeature-end presetが区別される。 |
| AC-07 | multipart、origin-spanning、odd/even midpointがcovered traversalで正しい。 |
| AC-08 | unstranded/ambiguous/fuzzy/cropped/linear/staleにoperation別の理由が出る。 |
| AC-09 | duplicate record IDとsplit fragmentでもstable source-bound identityを維持する。 |
| AC-10 | originとorientationが一つのUndo/Redo entryでResultと共に戻る。 |
| AC-11 | Session save/fresh Load後もabsolute transformとprovenanceが復元される。 |
| AC-12 | failed/canceled/stale renderで以前のResultとtransformが残る。 |
| AC-13 | unrelated pending editsは適用も破棄もされずpendingのまま残る。 |
| AC-14 | LOSAT executorの追加jobは0で、既存raw evidenceを再利用する。 |
| AC-15 | feature/label/tick/depth/statistics/comparison geometryが同じtransformに従う。 |
| AC-16 | manual start/orientation変更で古いanchor provenanceがclearされる。 |
| AC-17 | Cancelはdraft、Result、Historyを完全に変更しない。 |
| AC-18 | search queryとstable target再同定を保持し、popup操作がkeyboard/mobileで到達可能。 |
| AC-19 | request schema 7、Worker protocol、renderer pathを増やさない。 |
| AC-20 | Product Impactとarchitecture ratchet evidenceがreview可能で全gateがpassする。 |

## 12. Failure and stop conditions

次の場合はaffected workを止め、独断のfallbackを追加しない。

- Issue本文とmerged Product authorityが異なるoutcomeを要求する。
- source locationからtraversal/outgoing boundaryを安全に表現できず、catalog contract変更にも
  authorityがない。
- committed artifact/resourceからtarget-only candidateを構築できず、unrelated pending formを
  commitしない契約を守れない。
- per-record orientationを表すためにrequest schema/Worker protocolを変える必要があるように見える。
  まず既存`presentation.reverseComplement`とexact-one materializationを再監査する。
- Undo/Redoがintentとartifactを別entryにする設計しか成立しない。既存history ownerの拡張を再検討する。
- LOSAT executorが動く。silent skipではなくcache key/compatibility regressionを診断する。

## 13. Review checklist

Production、tests、docs/generated artifactsを別々に一度reviewする。

### Production

- target identityとsource freshnessがApply直前にも検証されるか
- source/display、0/1-based、base/boundaryの意味が名前で区別されるか
- resolver、projection、history、UIの責務が混ざっていないか
- normal Generateとpopupが同じadmission pathか
- target外のrecord/comparison/resourceが変わらないか
- manual edit時にprovenanceがstaleにならないか

### Tests

- default pathをUIから通っているか
- failure/cancel/stale/supersededを含むか
- plus/minus、compound、duplicate identity、both modesを含むか
- job countがcache metadataではなくexecutor invocationを測っているか
- snapshotだけでなくcanonical requestとvisible Resultを検証しているか

### Docs and schemas

- schema 4/v44 migrationとretirement conditionが明記されるか
- request schema 7不変が確認されるか
- public manualがsource coordinate、offset、end preset、disabled reasonを説明するか
- generated wheel、`dist/`、`gbdraw.egg-info/`、reference SVGを意図せずcommitしていないか

## 14. Completion handoff

Session 07のhandoffには次を含める。

- 実装結果をuser-visible outcomeから始める短いsummary
- branch名とHEAD
- owner/path/compatibility evidence
- 実行したcommandとpass/fail/未実行理由
- Issue #563 acceptance matrixの結果
- remaining riskまたは`none`
- Englishのproposed commit titleとshort summary

runtime実装のpush、PR作成、mergeは、そのsessionで明示的に許可された場合だけ行う。
