# Work package D: Manual feature placementの実装計画

- Date: 2026-08-11
- Status: ready to execute; implementation not started
- Baseline: `docs_renovation` / `a2c69de1` と文書作成時点のworking tree
- Source: [`gbdraw_v0.14.0_codex_roadmap.md`](./gbdraw_v0.14.0_codex_roadmap.md) の Work package D
- Target release: `v0.14.0`

この文書は、Work package D「Manual feature placement for overlap resolution」を実装可能な作業単位へ分解した計画である。文書作成時点ではroadmap以外のproduction code、test、public documentation、generated artifactを変更していない。

working treeには本計画と無関係な広範囲の変更が存在する。実装開始時のPhase 0で現状を取り直し、関連しない差分を巻き戻さない。本文中の行番号ではなく、ownerとして記載したsymbolとfileを再確認してから変更する。

## 1. 結論

manual placementとautomatic overlap resolutionを別の状態として実装する。

- manual placementは`resolve_overlaps`の値にかかわらず適用する。
- `resolve_overlaps`はAuto featureだけを移動する。
- `Auto`はoverrideの不在として保存し、選び直したときはderived laneをすべて再計算する。
- fixed featureを先に配置し、Auto featureはその後で既存のdeterministic orderを使う。
- raw `feature_track_id`は永続化しない。
- release minimumはAuto/Mainに加えて、bidirectional combined-strand feature slotの両方向lane 1を含む。
- `feature_overlap_tolerance_bp`は独立したdefault-0設定とし、manual placementと同じcollision predicateを使う。

placementのcanonical ownerはtyped `DiagramOptions`、toleranceのownerは`GbdrawConfig.canvas`とする。Web、CLI、beginner Python APIは入力形を変換するだけで、placement policyを複製しない。

## 2. 目標と対象範囲

### 2.1 必須目標

- resolver OFFでも任意のforeground featureをMainまたは対応するlane 1へ固定できる。
- resolver ONではfixed featureを先に予約し、Auto featureを周囲へ移動できる。
- 自動的に移動されたfeatureをMainへ固定すると、以前MainだったAuto featureが移動する。
- requested placementとresolved laneを区別する。
- multipart featureは一つのplacement unitとして全blockとconnectorを同じlaneへ移動する。
- final laneからfeature geometry、label、leader line、track reservation、canvas boundsを再計算する。
- source identityはrotation、crop後も残るfeature、reverse complement、record reordering、session replayを通じて維持する。sourceには存在するが現在のcrop外にあるfeatureのoverrideはdormantにする。
- Web、CLI、package-root Python API、typed API、canonical request、session replayで同じ結果にする。
- overrideなし、tolerance 0のdefault outputを変更しない。

### 2.2 対応するplacement

| Resolved feature-slot geometry | Auto | Main | 外側/上 lane 1 | 内側/下 lane 1 |
|---|---:|---:|---:|---:|
| Circular `lane_direction=split`, combined strands | yes | yes | Outward | Inward |
| Linear feature slot `side=overlay`, combined strands | yes | yes | Above | Below |
| Circular `lane_direction=inside` / `outside` | yes | yes | no | no |
| Linear feature slot `side=above` / `below` | yes | yes | no | no |
| Separate strands | yes | yes | no | no |

Separate-strandのMainはglobalな整数lane 0ではなく、そのfeatureが属するstrand poolのnominal laneを意味する。

Built-inのCircular `middle` presetは通常`lane_direction=split`へ、Linear `middle` presetは通常feature slot `side=overlay`へ解決される。ただしcustom track stackがslot geometryを上書きできるため、preset名だけで方向対応を判定しない。構造上のtoken/level検証はdecoderが行い、方向の可否はtrack slot解決後にshared placement plannerが検証する。

### 2.3 対象外

- lane 2以降
- arbitrary pixel dragging
- per-blockまたはper-exonの別lane配置
- multipart collisionをblock-onlyへ変えること
- `occupation_ratio`によるautomatic priority変更
- percentage overlap threshold
- label、annotation、track stackのpixel collision tolerance
- hidden/underlay featureのforeground lane予約
- Web-only DOM/SVG patchでplacementを表現すること
- Work package E全体のdebounce/auto-update実装

`feature_overlap_tolerance_bp`はrelease targetだが、schedule上の最初のscope-cut候補である。削減する場合はroadmapと本計画のledgerに明記し、manual lane 1契約を代わりに削らない。

## 3. 現状の根拠

| 現状 | Owner | 計画への影響 |
|---|---|---|
| CircularとLinearは同じfeature allocatorを使う | [`gbdraw/features/factory.py`](../../gbdraw/features/factory.py) の `_build_feature_layers()` と [`gbdraw/features/tracks.py`](../../gbdraw/features/tracks.py) の `arrange_feature_tracks()` | fixed-firstとAuto allocationはこの共有境界に置く |
| `feature_track_id`は`FeatureObject`生成時に0で始まるtransient値 | [`gbdraw/features/objects.py`](../../gbdraw/features/objects.py) | persisted intentやWeb keyに使わない |
| collisionはhalf-open intervalの任意の正の交差を拒否する | [`gbdraw/features/tracks.py`](../../gbdraw/features/tracks.py) の `check_feature_overlap()` | overlap長を返す共有predicateへ拡張し、default 0を維持する |
| multipart collisionは最小startから最大endまでのouter envelopeを使う | 同fileの `get_feature_ends()` | intron内featureはtoleranceだけでは救済できない。manual placementで解決する |
| `occupation_ratio`は計算されるがsort keyでは未使用 | 同fileの `arrange_feature_tracks()` | v0.14ではpriorityを変更しない |
| Circularの両方向配置はresolved feature slotの`lane_direction=split`で決まる | [`gbdraw/diagrams/circular/assemble.py`](../../gbdraw/diagrams/circular/assemble.py) | preset名ではなくfinal slot geometryでdirection supportを決める |
| Circular radial layoutは実在するtrack IDだけを集めて密にpackする | [`gbdraw/diagrams/circular/radial_layout.py`](../../gbdraw/diagrams/circular/radial_layout.py) | manual lane 1だけが存在してもMainとlane 1の意味的なbandを保持する必要がある |
| Circular label側の一部はstrandと`resolve_overlaps`をlane sideのproxyにする | [`gbdraw/labels/circular.py`](../../gbdraw/labels/circular.py) | final physical sideをlabel plannerへ渡す |
| Linearの実効layoutはfeature slotの`side`がprofile presetを上書きする | [`gbdraw/layout/linear.py`](../../gbdraw/layout/linear.py)、[`gbdraw/labels/linear.py`](../../gbdraw/labels/linear.py) | `side=overlay`だけを両方向対応にし、featureとlabelが同じresolved geometryを使う |
| Circular CLIとWebはseparate strands時にresolverを無効化する | [`gbdraw/circular.py`](../../gbdraw/circular.py)、[`gbdraw/web/js/app/app-setup.js`](../../gbdraw/web/js/app/app-setup.js) | Auto/Mainをseparate strandsで検証するため、この禁止を削除してlayout/labelを完成させる |
| Webのstable override keyは`recordKey + biologicalFeatureId` | [`gbdraw/web/js/services/feature-catalog.js`](../../gbdraw/web/js/services/feature-catalog.js) | canonical identityも同じ二要素にする |
| duplicate biological hashはsource indexで`~index` disambiguationされる | [`gbdraw/web_support/feature_catalog.py`](../../gbdraw/web_support/feature_catalog.py) | identity生成をWeb専用moduleからshared coreへ移す |
| source indexはrecord transformを通してcopyされる | [`gbdraw/core/record_metadata.py`](../../gbdraw/core/record_metadata.py) | placement解決はtransform後も同じsource identityを使える |
| typed shared option ownerは`_ModeDiagramOptions` | [`gbdraw/api/options.py`](../../gbdraw/api/options.py) | placement tupleとtable inputをここへ追加する |
| `resolve_overlaps`は`CanvasConfig`にある | [`gbdraw/config/models/canvas.py`](../../gbdraw/config/models/canvas.py) | toleranceも同じglobal collision-policy ownerへ置く |
| Web canonical projectionは一つのschema-5 requestを作る | [`gbdraw/web/js/services/session-request.js`](../../gbdraw/web/js/services/session-request.js) | placement mapをここで一度だけsorted arrayへ射影する |
| current main writerはrequest schema 5 / session 40で、Work package C Phase 3がschema 6 / session 41 writerを所有する | [`gbdraw/session_request_codec.py`](../../gbdraw/session_request_codec.py)、[`gbdraw/session_io.py`](../../gbdraw/session_io.py)、[`WORK_PACKAGE_C_RECORD_COORDINATE_ANCHORING_IMPLEMENTATION_PLAN_2026-08-11.md`](./WORK_PACKAGE_C_RECORD_COORDINATE_ANCHORING_IMPLEMENTATION_PLAN_2026-08-11.md) | D persisted fieldsはCの原子的writer updateへfoldし、D-only 6/41を作らない |
| factoryは`SeqRecord`だけを受け、canonical `recordKey`を知らない | [`gbdraw/features/factory.py`](../../gbdraw/features/factory.py)、[`gbdraw/api/request_render.py`](../../gbdraw/api/request_render.py) | provenance順の`ResolvedPlacementInputs`をrequest planからdiagram/assembler/configurator/factoryへ明示的に運ぶ |
| sessionの`renderRequest`はcommitted renderのowner | [`docs/SESSION_COMPATIBILITY.md`](../SESSION_COMPATIBILITY.md) | draft mapとcommitted placementを区別し、`editorState`を第二のauthorityにしない |

## 4. 固定する製品判断

### D1. `resolve_overlaps`はAutoだけを制御する

| Requested placement | Resolver OFF | Resolver ON |
|---|---|---|
| Auto | nominal unresolved lane。Auto同士やfixedとのoverlapを許す | fixed laneを避けてexisting deterministic allocatorで配置 |
| Main | nominal main laneへ固定 | nominal main laneを先に予約 |
| Directional lane 1 | 指定laneへ固定 | 指定laneを先に予約 |

fixed-fixedが同じphysical laneを要求し、そのenvelope intersectionがtoleranceを超える場合は両resolver stateで`ValidationError`にする。一方がAutoでresolver OFFならoverlapを許す。

### D2. `Auto`はoverrideの削除

Auto rowをcanonical requestやsessionへ書かない。WebでAutoを選ぶ、Pythonでoverrideを除く、TSVで`auto`を指定する、の3経路は同じ空のoverrideへ正規化する。

再描画時に前回のresolved laneをseedしない。fresh `FeatureObject`からfixed、Autoの順に割り当て直す。

### D3. requested intentとresolved geometryを分ける

persistするのは次のrequested intentだけである。

```text
Main
Outward lane 1
Inward lane 1
Above lane 1
Below lane 1
```

renderer内部ではsemantic targetをphysical side、strand pool、levelへ解決してからtransient `feature_track_id`を導出する。fixed conflictはraw strand文字列ではなく、同じphysical pool/laneを共有するかで判定する。undefined strandは既存policyが割り当てるnominal poolに属し、独立した第三poolを暗黙に作らない。labelとcanvasはraw IDや`resolve_overlaps`ではなく、同じresolved assignmentを参照する。

### D4. unsupported mappingはerror

例:

- Circular resolved feature slot `lane_direction=outside`でInward
- Circular separate strandsでInward
- Linear resolved feature slot `side=above`でBelow
- Linear separate strandsでAbove
- level 2以上

これらをMain、native side、Autoへ黙って変換しない。Webはfinal slot geometryが分かるときに選択肢をdisabledにする。typed request/session decoderはtokenとshapeを検証し、config/track-slot解決後のshared plannerがlayout可否を必ず検証する。

### D5. source-knownだが非表示のoverrideはdormant

選択されたoriginal source recordにsource identityが存在する限りoverrideを保持する。hidden、underlay、または現在のcrop外にあるときはforeground occupancyへ入れない。visibility、shape、cropが後からfeatureをforegroundへ戻したときに同じoverrideを再適用する。

original source catalogに存在しないunknown/stale identityはdormant扱いにせずerrorにする。crop後も残るfeatureはsource indexから同じidentityへbindし、request-owned cropで完全に除外されたsource-known featureとunknown IDを混同しない。callerが既にcrop済みの`SeqRecord`だけを渡し、元sourceがrequest/resourceに存在しない場合は、見えない元featureを推測せずunknownとして扱う。

### D6. collision toleranceは一つのpredicate

`feature_overlap_tolerance_bp`はboolを拒否する非負整数、default 0とする。

coreは概念的に次を一箇所で所有する。

```python
feature_overlap_bp(a, b, *, separate_strands, genome_length) -> int
features_conflict(a, b, *, tolerance_bp, separate_strands, genome_length) -> bool
```

conflictは`feature_overlap_bp(...) > tolerance_bp`である。`IntervalIndex`はcandidate抽出だけを担当し、fixed validationと全Auto allocatorが同じfinal predicateを呼ぶ。Circular origin splitはdisjoint pieceの交差長を合計する。

### D7. multipartはatomic placement unit

全blockとconnectorを一緒に移動する。collision basisは現在のouter envelopeを維持する。block-onlyへ変えるとintron connectorが別feature glyphを横切る可能性があるため、v0.14のdefaultにはしない。

### D8. stable identityはshared coreが生成する

canonical identityは次のtupleである。

```text
(recordKey, biologicalFeatureId)
```

`biologicalFeatureId`は全location partを含むstable hashをbaseにし、同一recordでbaseが重複する場合はsource feature indexを加える。identity生成とvalidationを`web_support`へ依存させず、renderer、typed resolver、feature catalogが同じhelperを使う。

exact overrideはprovenance順の`ResolvedPlacementInputs`へmaterializeする。各record entryは`record_key`と、original source上で解決した`(biological_feature_id, source_feature_index, target)`を持つ。factoryはcrop/reverse complement後の座標からIDを再hashせず、transformを通じてcopyされたsource indexでruntime featureへbindする。`PreparedDiagramInputs`、request plan、`gbdraw/api/diagram.py`、mode assembler、`FeatureDrawingConfigurator`、factoryの順にこのrecord-aligned valueを明示的に渡す。`SeqRecord.id`やlist positionだけでrecordを再同定せず、duplicate record IDを持つmulti-recordでもcanonical `recordKey`で分離する。

### D9. request/session formatはCと一度だけ進める

planning時点の`main`にはcanonical request schema 5とsession 40があり、最新release tagは`0.13.0`である。CとDのfinal persisted fieldsを一つのformat changeとして扱う。

Work package C Phase 3をschema 6 / session 41の唯一の次期writer ownerとする。Dのpersisted fieldsは、Cのrecord `display`、format-era predicate、v40→v41 promotionと同じ原子的changeへfoldする。D単独で6/41 writerを作ってはならない。

統合時のtargetは次のとおりである。

- canonical request schema 6
- session version 41
- schema 5 / session 40をmain-era positive fixtureとして読み続ける
- schema 5 / session 40にplacementがない場合はempty overrides、tolerance missingは0

Phase 0でWork package Cの実装状態を再確認する。C Phase 3が未着手または未landならC+D統合phaseとして6/41を一度だけ書く。Cの6/41 writerがactive branch上だけにある場合は完成前にDをfoldし、さらにversionを増やさない。Cのformatがmainまたはrelease-backed historyへDなしで既に入っている場合はschema 6/session 41を黙って再定義せず、次のversion pairをcompatibility evidence付きで決める。

### D10. last successful renderを守る

Webのdraft placementは`config`側のeditable state、committed placementは`renderRequest.diagramOptions.featurePlacements`とする。両者が異なる状態を許し、stale indicatorを表示する。

invalid、cancelled、stale、failed candidateはpreview、feature catalog、selection、export resultを更新しない。placementとtoleranceはrender-only classificationとし、committed comparison evidenceを再利用してLOSATを起動しない。

## 5. 公開contract

### 5.1 Canonical request

`diagramOptions.featurePlacements`はnon-Auto rowだけを持つsorted arrayとする。

```json
{
  "featurePlacements": [
    {
      "recordKey": "record-1",
      "biologicalFeatureId": "f1234abcd",
      "placement": {"kind": "main"}
    },
    {
      "recordKey": "record-1",
      "biologicalFeatureId": "f5678ef00~12",
      "placement": {"kind": "lane", "side": "outward", "level": 1}
    }
  ]
}
```

規則:

- row orderは`recordKey`, `biologicalFeatureId`でsortする。
- duplicate identityを拒否する。
- `kind=main`には`side`と`level`を書かない。
- `kind=lane`にはmodeに合う`side`と`level=1`が必要。decoderはこの構造を検証し、resolved feature slotがそのsideを支えるかはshared plannerが後で検証する。
- `auto`、NUL結合Web key、SVG ID、raw track IDをwire formatへ書かない。
- unknown fieldは現行policyどおり拒否する。

### 5.2 Typed API

公開型を増やす場合は次の二つに限定する。

```text
FeaturePlacementTarget
FeaturePlacementOverride
```

`FeaturePlacementOverride`は`record_key`、`biological_feature_id`、`target`を持つfrozen dataclassとする。`FeaturePlacementTarget`はMainまたはmode-specific lane 1を表す。

`_ModeDiagramOptions`へ次を追加する。

```text
feature_placements
feature_placement_table
feature_placement_table_file
```

exact overridesとtable/DataFrame/pathを同時に指定した場合は`ValidationError`にする。`resolve_request()`はrecord selection/normalization後、cropやvisibility filterの前にoriginal source catalogを使ってselector rowをexact overrideへ変換し、materialized requestからtable/file fieldを除く。canonical encoderはmaterialized exact overridesだけを受け入れる。

package-root `FeatureOptions.placements`はpath、DataFrame、またはexact override sequenceを受け、`_base_options()`で一度だけtyped optionsへ射影する。新しいpublic型は`gbdraw.api`とpackage rootの必要な`__all__`、public contract fixtureへ同時に追加する。

`feature_overlap_tolerance_bp`は`CanvasConfig`が唯一のownerである。Pythonでは`config_overrides={"canvas.feature_overlap_tolerance_bp": 1}`を使い、beginner facadeへ重複fieldを追加しない。

### 5.3 Helper TSV

CLIとpackage-root Python APIは同じUTF-8 TSV readerを使う。

```tsv
record	feature_selector	placement	level
#1	locus_tag=exampleA	main
#1	locus_tag=exampleB	outward	1
```

| Column | Required | Contract |
|---|---|---|
| `record` | multi-record時required | existing `RecordSelector` syntax。single recordではblank可 |
| `feature_selector` | yes | existing `FeatureSelector` syntax。exactly one source featureへ解決すること |
| `placement` | yes | `auto`, `main`, `outward`, `inward`, `above`, `below` |
| `level` | directional時optional | blankは1。v0.14では1だけを受理 |

zero match、multiple match、duplicate resolved identity、unknown column、wrong-mode placementをdirect errorにする。qualifier selectorのambiguityをsource orderで決めない。

CLI optionは両mode共通で次を使う。

```text
--feature_placement_table TSV
--feature_overlap_tolerance_bp BP
```

CLIはtable policyを実装せず、pathをtyped requestへ渡す。session saveでhelper pathを認識するが、replay authorityはmaterialized canonical overridesであり、original TSVを第二の意味ownerにしない。

## 6. Surface matrix

| Surface | Placement owner/input | Tolerance owner/input | Final consumer | Persistence |
|---|---|---|---|---|
| Shared typed core | `_ModeDiagramOptions.feature_placements` + record-aligned `ResolvedPlacementInputs` | `GbdrawConfig.canvas.feature_overlap_tolerance_bp` | resolved-slot-aware placement planner | canonical request |
| Package-root Python | `FeatureOptions.placements` | `config_overrides` | typed options adapter | caller-owned unless session writer used |
| Typed Python API | `FeaturePlacementOverride`またはhelper table | `config_overrides` | `resolve_request()` | schema-6 request/session-41 target |
| CLI | `--feature_placement_table` | `--feature_overlap_tolerance_bp` | same typed request | materialized session request |
| Web | feature editor draft map | Advanced numeric input | `session-request.js` projection | draft config + committed renderRequest |
| Circular renderer | shared resolved assignment | render profile | radial layout、labels、canvas | derived only |
| Linear renderer | shared resolved assignment | render profile | lane geometry、labels、corridor | derived only |

## 7. 実装phase

各phaseはproduction diff、test diff、documentation diff、generated artifact diffを分けて確認する。前phaseのgateが通るまで次phaseをcompleteにしない。

### Phase 0: baselineとformat decisionを固定する

Status: not started

Owned systems:

- current worktree and branch baseline
- Work package Cとのschema/session coordination
- existing reference-output diff
- allocator、identity、sessionのcharacterization evidence

Work:

1. `git status --short --branch`とin-scope diffを記録する。
2. Work package Cの実装/plan状態、current request/session writer、first-parent/tag evidenceを確認する。
3. Work package C Phase 3と統合してschema 6 / session 41を一度だけ書けるか確認する。Cのwriterが既にmain/release-backed historyへ入っている場合は、D用の次version pairを決める。D-only 6/41は禁止する。
4. `tests/reference_outputs/`のpre-existing diffを記録し、Dのevidenceと混ぜない。
5. current allocator、duplicate identity、reverse-complement identityのfocused testsを実行する。
6. `NC_001879.gbk`などmultipart fixtureから、Mainを占めるsparse featureとintron内featureを含むQA caseを選定する。

Gate:

```bash
python -m pytest \
  tests/test_feature_tracks.py \
  tests/test_feature_catalog_source_identity.py \
  tests/test_reverse_complement_feature_identity.py -q
```

Evidence:

- exact command/result
- C Phase 3との統合点、または次version pairを必要とする履歴証拠を含むfinal schema/session decision
- selected fixtures and why they cover required geometry
- pre-existing reference-output status

### Phase 1: placement model、selector resolution、shared identity

Status: not started

Depends on: Phase 0 and Work package C Phase 2's shared record-aligned provenance. If C Phase 2 is not complete, establish that common provenance contract first; do not create a D-only parallel record identity path.

Primary files:

- new `gbdraw/features/placement.py`
- [`gbdraw/features/ids.py`](../../gbdraw/features/ids.py)
- [`gbdraw/core/record_metadata.py`](../../gbdraw/core/record_metadata.py)
- [`gbdraw/web_support/feature_catalog.py`](../../gbdraw/web_support/feature_catalog.py)
- [`gbdraw/api/options.py`](../../gbdraw/api/options.py)
- [`gbdraw/api/prepared.py`](../../gbdraw/api/prepared.py)
- [`gbdraw/api/request_render.py`](../../gbdraw/api/request_render.py)
- [`gbdraw/api/record_planning.py`](../../gbdraw/api/record_planning.py)
- [`gbdraw/api/io.py`](../../gbdraw/api/io.py)
- [`gbdraw/annotations/models.py`](../../gbdraw/annotations/models.py) and selector matching owner as needed
- new `tests/test_feature_placement.py`
- [`tests/test_feature_catalog_source_identity.py`](../../tests/test_feature_catalog_source_identity.py)

Work:

1. semantic target、exact override、resolved physical assignmentのimmutable internal modelsを定義する。
2. biological feature IDとduplicate disambiguationをshared core helperへ移し、Web catalogとplacement resolverが同じ実装を使う。
3. existing `FeatureSelector` matchingをshared helperへ抽出し、placementはexactly-one規則を追加する。annotationのmulti-match semanticsは変えない。
4. `_ModeDiagramOptions`へinternal typed fieldsとtable/file mutual exclusionを追加し、renderer phasesがprivate bypassなしでmanual intentを受け取れるようにする。
5. `record_planning.py`でrequest-owned cropを適用する前に、全source featureのimmutable identity/index catalogを一度生成し、Cのrecord-aligned provenanceまたはそれに隣接する一対一aligned valueへ保持する。resourceを再読込してcatalogを再構築しない。
6. table rowをrecord selection/normalization後、cropやvisibility filterがsource featureを除去する前のoriginal source catalogへ解決する。
7. Auto rowをdropし、exact overrideをidentity順へnormalizeする。
8. provenance順の`ResolvedPlacementInputs`を作り、各entryをcanonical `record_key`、original `biological_feature_id`、`source_feature_index`、targetと対応づけて`PreparedDiagramInputs`とrequest planに保持する。
9. hidden/underlay、source-known crop-excluded、unknown/stale、duplicate targetを別々にtestする。
10. duplicate `SeqRecord.id`を持つmulti-recordでも二つのrecord keyが混線しないことをtestする。

Gate:

```bash
python -m pytest \
  tests/test_feature_placement.py \
  tests/test_feature_catalog_source_identity.py \
  tests/test_reverse_complement_feature_identity.py \
  tests/test_feature_selector_values.py \
  tests/test_api_requests.py \
  tests/test_record_planning.py \
  tests/test_annotations.py -q
```

### Phase 2: shared collision lengthとtolerance

Status: not started

Depends on: Phase 1のidentity/model。manual placementより先にpredicateを固定する。

Primary files:

- [`gbdraw/features/tracks.py`](../../gbdraw/features/tracks.py)
- [`gbdraw/layout/spatial.py`](../../gbdraw/layout/spatial.py) if candidate-query support is required
- [`gbdraw/config/models/canvas.py`](../../gbdraw/config/models/canvas.py)
- [`gbdraw/config/models/render_profiles.py`](../../gbdraw/config/models/render_profiles.py)
- [`gbdraw/data/config.toml`](../../gbdraw/data/config.toml)
- [`tests/test_feature_tracks.py`](../../tests/test_feature_tracks.py)
- [`tests/test_spatial_index.py`](../../tests/test_spatial_index.py)

Work:

1. overlap lengthとboolean conflictを分ける。
2. current half-open、separate-strand、origin-spanning semanticsをtolerance 0で維持する。
3. tolerance 1/2 bp boundaryを追加する。
4. `CanvasConfig.from_dict()`はmissing valueを0としてsupported older configを読む。
5. bool、negative、fraction、nonnumericをrejectする。
6. `IntervalIndex`をcandidate prefilterのまま維持し、O(n²) fallbackを追加しない。

Gate:

```bash
python -m pytest \
  tests/test_feature_tracks.py \
  tests/test_spatial_index.py \
  tests/test_config_override_paths.py \
  tests/test_render_profiles.py -q
```

Scope-cut gate:

このphaseをreleaseから外す場合は、manual plannerがhard-coded tolerance 0で同じpredicateを使える状態を残す。roadmap、public docs、CLI/Web control、ledgerからtolerance claimを同じchangeで削除する。

### Phase 3: fixed-first shared allocator

Status: not started

Depends on: Phase 1、Phase 2またはdocumented tolerance scope cut

Primary files:

- [`gbdraw/features/objects.py`](../../gbdraw/features/objects.py)
- [`gbdraw/features/factory.py`](../../gbdraw/features/factory.py)
- [`gbdraw/features/tracks.py`](../../gbdraw/features/tracks.py)
- [`gbdraw/api/prepared.py`](../../gbdraw/api/prepared.py)
- [`gbdraw/api/request_render.py`](../../gbdraw/api/request_render.py)
- [`gbdraw/api/diagram.py`](../../gbdraw/api/diagram.py)
- [`gbdraw/configurators/features.py`](../../gbdraw/configurators/features.py)
- [`gbdraw/diagrams/circular/assemble.py`](../../gbdraw/diagrams/circular/assemble.py) for record-aligned transport only
- [`gbdraw/diagrams/linear/assemble.py`](../../gbdraw/diagrams/linear/assemble.py) for record-aligned transport only
- [`tests/test_feature_tracks.py`](../../tests/test_feature_tracks.py)
- new `tests/test_feature_placement.py`

Work:

1. `ResolvedPlacementInputs`をrequest planから`gbdraw/api/diagram.py`、mode assembler、`FeatureDrawingConfigurator`、factoryへrecord-aligned valueとして渡す。assemblerは検証済みprovenance orderからrecord entryを選び、shared configuratorへ一つのrecord mappingを誤って使い回さない。
2. factoryへcanonical `record_key`とそのrecordだけのresolved overridesを渡し、copy済みsource indexでsource identityをruntime featureへbindする。transform後の座標を再hashせず、`SeqRecord.id`をrecord-key substituteにしない。
3. semantic laneをresolved feature slot、mode、strand poolへ解決する。token/levelはdecoder、direction/layout compatibilityはこの境界で検証する。
4. explicit foreground featureをstable orderでseedし、fixed-fixed conflictを先に検証する。
5. resolver ONではexisting Auto sortとnearest-lane policyを維持してfixed occupancyの周囲へ配置する。
6. resolver OFFではAutoをnominal laneへ置き、fixedとのoverlapを解消しない。
7. final assignmentからtransient track IDを導出する。
8. Autoへの復帰が過去laneを再利用しないことを証明する。
9. 100-lane limitとerror typeを維持する。

Gate cases:

- OFF + manual directional
- ON + B MainでA Autoが移動
- 3-way overlap
- fixed-fixed conflictとtolerance boundary
- non-overlapping fixed featuresのlane共有
- hidden/underlay/source-known crop-excluded dormant
- duplicate record IDを持つmulti-recordのrecord-key isolation
- multipart envelope
- origin wrap
- deterministic repeated render

```bash
python -m pytest \
  tests/test_feature_placement.py \
  tests/test_feature_tracks.py \
  tests/test_feature_underlay_rendering.py -q
```

### Phase 4: Circular geometryとlabels

Status: not started

Depends on: Phase 3

Primary files:

- [`gbdraw/diagrams/circular/assemble.py`](../../gbdraw/diagrams/circular/assemble.py)
- [`gbdraw/diagrams/circular/radial_layout.py`](../../gbdraw/diagrams/circular/radial_layout.py)
- [`gbdraw/layout/circular.py`](../../gbdraw/layout/circular.py)
- [`gbdraw/labels/circular.py`](../../gbdraw/labels/circular.py)
- [`gbdraw/labels/circular_candidates.py`](../../gbdraw/labels/circular_candidates.py)
- [`gbdraw/circular.py`](../../gbdraw/circular.py) for separate-strand restriction removal
- relevant Circular tests

Work:

1. resolved feature slotが`lane_direction=split`でmanual laneがある場合は、resolver OFFでもMain、outer、inner bandを測定する。
2. lone lane 1がdense packingでMainへcollapseしないようsemantic levelをmaterializeする。
3. labelsとleader anchorをfinal physical sideから決める。
4. `resolve_overlaps`をlane-clearanceのproxyにしているbranchをresolved occupancyへ置き換える。
5. adjacent numeric/depth/annotation slotsとcenter reservationを再計算する。
6. Circular separate strands + resolverのPython/CLI禁止を外し、Auto/Mainの両stateをlayout/labelまでtestする。Web側のwatcherはPhase 7で同じcontractへ合わせる。
7. custom slotの`lane_direction=inside` / `outside`を含め、split以外のdirectional overrideをtyped `ValidationError`にする。

Gate:

```bash
python -m pytest \
  tests/test_circular_radial_layout.py \
  tests/test_circular_label_placement.py \
  tests/test_circular_radial_labels.py \
  tests/test_circular_label_track_policy.py \
  tests/test_circular_track_slots.py \
  tests/test_circular_composition.py -q
```

Visual evidence:

- `lane_direction=split`（built-in middleとcustom slot）、resolver OFF、Outward/Inward lane 1
- resolver ON、Main/Auto precedence
- multipart/intron case
- separate strandsのAuto/Main
- labels、leader lines、center text、adjacent tracks、canvas clippingをreadable scaleで確認

### Phase 5: Linear geometryとlabels

Status: not started

Depends on: Phase 3

Primary files:

- [`gbdraw/diagrams/linear/precalc.py`](../../gbdraw/diagrams/linear/precalc.py)
- [`gbdraw/diagrams/linear/positioning.py`](../../gbdraw/diagrams/linear/positioning.py)
- [`gbdraw/diagrams/linear/assemble.py`](../../gbdraw/diagrams/linear/assemble.py)
- [`gbdraw/layout/linear.py`](../../gbdraw/layout/linear.py)
- [`gbdraw/labels/linear.py`](../../gbdraw/labels/linear.py)
- relevant Linear tests

Work:

1. resolved feature slotが`side=overlay`でmanual above/below laneがある場合は、resolver OFFでもfeature bandへ含める。
2. combined `side=overlay`のAbove/Below lane 1をraw ID signに依存せず解決する。
3. label anchor、occupied label band、track stack、record spacing、comparison corridorを同じresolved geometryから計算する。
4. custom slotを含む`side=above` / `below`とseparate strandsのunsupported directionをrejectする。
5. Mainを各layout/strand poolのnominal laneとしてtestする。

Gate:

```bash
python -m pytest \
  tests/test_linear_vertical_layout.py \
  tests/test_linear_track_layout.py \
  tests/test_linear_label_placement.py \
  tests/test_linear_track_slots.py \
  tests/test_linear_composition.py \
  tests/test_linear_multi_record_layout.py -q
```

Visual evidence:

- combined `side=overlay`（built-in middleとcustom slot）、resolver OFF/ON、Above/Below lane 1
- separate strandsのAuto/Main
- multipart feature
- multi-record comparison corridorがfeature/label bandを横切らない

### Phase 6: public API、beginner API、CLI、C+D request/session integration

Status: not started

Phase 6A depends on: Phase 1、Phase 3、Phase 0のformat decision。

Phase 6B depends on: Phase 6A、Circular/Linear renderer gates、Work package C Phase 2。Phase 6BはWork package C Phase 3そのものであり、別々に実行・land・completeにしない。

Primary files:

- [`gbdraw/api/options.py`](../../gbdraw/api/options.py)
- [`gbdraw/api/prepared.py`](../../gbdraw/api/prepared.py)
- [`gbdraw/api/request_render.py`](../../gbdraw/api/request_render.py)
- [`gbdraw/api/__init__.py`](../../gbdraw/api/__init__.py)
- [`gbdraw/api/io.py`](../../gbdraw/api/io.py)
- [`gbdraw/interface.py`](../../gbdraw/interface.py)
- [`gbdraw/__init__.py`](../../gbdraw/__init__.py)
- [`gbdraw/cli_utils/common.py`](../../gbdraw/cli_utils/common.py)
- [`gbdraw/cli_utils/session.py`](../../gbdraw/cli_utils/session.py)
- [`gbdraw/circular.py`](../../gbdraw/circular.py)
- [`gbdraw/linear.py`](../../gbdraw/linear.py)
- [`gbdraw/session_request_codec.py`](../../gbdraw/session_request_codec.py)
- [`gbdraw/session_io.py`](../../gbdraw/session_io.py)
- [`gbdraw/session.py`](../../gbdraw/session.py)
- [`gbdraw/api/session_compat.py`](../../gbdraw/api/session_compat.py)
- [`gbdraw/web/js/state.js`](../../gbdraw/web/js/state.js)
- [`gbdraw/web/js/services/session-request.js`](../../gbdraw/web/js/services/session-request.js)
- [`gbdraw/web/js/services/config.js`](../../gbdraw/web/js/services/config.js)
- [`gbdraw/web/js/services/session-authority.js`](../../gbdraw/web/js/services/session-authority.js)
- contract/session/API tests and fixtures

Phase 6A — public and format preparation:

1. Phase 1のinternal model/typed fieldsをpublic dataclassesとexportsへ昇格し、public contractを固定する。
2. CLI flagsを両modeへ追加し、Phase 1のshared reader/materializerへ渡す。
3. canonical codecのexplicit nested placement encoder/decoderとWeb projection helperを準備する。generic dataclass JSON fallbackへ任せない。
4. public exportsと`tests/fixtures/public_contract.json`はintentional diffだけを更新する。
5. 6Aをschema 5/session 40 writerのまま単独landしない。runtime surfaceとpersisted surfaceのparityは6Bまで一つのchangesetとして保つ。

Phase 6B — Work package C Phase 3と同一のatomic writer update:

1. Webの`featurePlacementOverrides` draft mapとtolerance draft valueを追加し、`session-request.js`でnon-Auto placementをsorted canonical arrayへ一度だけ射影する。
2. canonical requestからWeb draftへ逆射影し、schema 5/session 40ではempty placementとtolerance 0を復元する。`editorState`を第二のplacement authorityにしない。
3. Cのrecord `display`、format-era predicate、v40→v41 promotionとD placement/toleranceを同じschema 6/session 41 writer changeとして実装する。
4. request schema/session versionとPython/Web constants、fixtures、compatibility docsをatomically更新する。D-only 6/41 writer、C-only 6/41 writer、または片方のWeb projectionを欠くfixtureを作らない。
5. schema 5/session 40 positive fixtureがdisplay default、empty placement、tolerance 0でreplayし、v40 load→saveが完全なschema 6/session 41へpromotionされることをtestする。
6. CLI session saveがhelper tableを扱い、replayはmaterialized placementsを使うことをtestする。
7. Phase 0でCの6/41がmain/release-backed historyへ既に入っていると判定した場合だけ、決定済みの次version pairへこのmatrixを移す。

Phase 6 completion rule:

- 6Aと6BのPython/Web testsが同じcommit candidateで通るまでPhase 6をcompleteにしない。
- Work package C Phase 3も同じevidenceを参照し、C-only writer evidenceではcompleteにしない。

Gate:

```bash
python -m pytest \
  tests/test_api_requests.py \
  tests/test_api_request_render.py \
  tests/test_api_library_usage.py \
  tests/test_dead_api_cleanup.py \
  tests/test_public_contract.py \
  tests/test_cli_tables.py \
  tests/test_session_request_codec.py \
  tests/test_session_io.py \
  tests/test_session_compat.py -q

node tests/web/feature-placement-request.test.mjs
node tests/web/session-request.test.mjs
node tests/web/session-file.test.mjs
node tests/web/session-draft-authority.test.mjs
```

### Phase 7: Web editor、history、render orchestration

Status: not started

Depends on: Phase 6、Circular/Linear renderer gates

Primary files:

- [`gbdraw/web/js/app/feature-editor.js`](../../gbdraw/web/js/app/feature-editor.js)
- new `gbdraw/web/js/app/feature-editor/placement-actions.js`
- [`gbdraw/web/js/app/app-setup.js`](../../gbdraw/web/js/app/app-setup.js)
- [`gbdraw/web/index.html`](../../gbdraw/web/index.html)
- [`gbdraw/web/js/services/history-snapshot.js`](../../gbdraw/web/js/services/history-snapshot.js)
- [`gbdraw/web/js/services/history.js`](../../gbdraw/web/js/services/history.js)
- [`gbdraw/web/js/app/run-analysis.js`](../../gbdraw/web/js/app/run-analysis.js)
- Web unit and browser tests

Work:

1. Phase 6Bの`featurePlacementOverrides` mapを既存biological identity keyで編集するcontrolをfeature popup/editorへ追加する。
2. Autoはmap entryを削除し、single/bulk actionを一つのundo itemにする。
3. committed requestとdraft mapを区別し、draft mismatchを既存stale-state contractへ接続する。
4. tolerance inputをAdvancedへ追加し、Phase 6Bのconfig projectionを更新する。
5. Circular separate strandsとresolverのmutual disable watcher/help textを削除し、supported combinationsを表示する。
6. invalid candidateはlast preview/catalog/selectionを残す。
7. placement/tolerance changeがcomparison analysis fingerprintを変えず、LOSAT artifactを再利用することをtestする。
8. Work package Eが未実装ならGenerate/Update actionを維持し、Dだけで新しいdebounce loopを作らない。

Unit gate:

```bash
node tests/web/feature-placement.test.mjs
node tests/web/session-request.test.mjs
node tests/web/feature-override-identity.test.mjs
node tests/web/history.test.mjs
node tests/web/session-draft-authority.test.mjs
```

Browser gate:

```bash
python tools/prepare_browser_wheel.py
npx playwright test \
  tests/web/feature-placement.playwright.spec.js \
  --project=chromium --workers=1
```

Browser cases:

- popupでOutward/Aboveを選択、resolver OFFでGenerate、geometry確認
- resolver ONでauto-displaced featureをMainへ戻す
- Auto、undo、redo
- session save/load/re-render
- unsupported directionのdraft errorとlast preview保持
- tolerance 0/1 boundary
- placement/tolerance変更でLOSATが起動しない

### Phase 8: public documentation、visual QA、broad gates

Status: not started

Depends on: Phase 1–7

Primary documentation:

- [`docs/REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md`](../REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md)
- [`docs/REFERENCE/input-formats-and-tsv-schemas.md`](../REFERENCE/input-formats-and-tsv-schemas.md)
- [`docs/REFERENCE/command-line.md`](../REFERENCE/command-line.md)
- [`docs/REFERENCE/python-api.md`](../REFERENCE/python-api.md)
- [`docs/REFERENCE/typed-requests.md`](../REFERENCE/typed-requests.md)
- [`docs/REFERENCE/session-and-request-compatibility.md`](../REFERENCE/session-and-request-compatibility.md)
- [`docs/SESSION_COMPATIBILITY.md`](../SESSION_COMPATIBILITY.md)
- generated [`docs/CLI_Reference.md`](../CLI_Reference.md)
- release notes and roadmap completion ledger

Work:

1. resolver truth table、supported layout matrix、identity、TSV、API、conflict、toleranceを正確に記載する。
2. CLI inventoryをgeneratorで更新し、hand-editしない。
3. package-root Python exampleとtyped request exampleをliteral executionする。
4. `NC_001879.gbk`などのmultipart fixtureでCircular/Linear QA artifactsを一時directoryへ生成し、readable scaleで目視する。
5. existing `tests/reference_outputs/`をcomparison-onlyで確認する。default diffが発生した場合は原因を直し、manual geometryを理由に既存一式を再生成しない。
6. 新しいpublic Gallery exampleは必須にしない。Work package Aのdocumentation regenerationがplacementを教えると決めた場合だけ、Gallery-quality sourceから生成する。
7. production、test、docs、generated artifact diffを別々にreviewする。

Documentation/CLI gate:

```bash
python tools/update_cli_reference_help.py
python tools/update_cli_reference_help.py --check
python -m pytest \
  tests/test_cli_reference.py \
  tests/test_documentation_reference_contracts.py \
  tests/test_documented_recipes.py -q
```

Default rendering gate:

```bash
python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
```

Broad gate:

```bash
ruff check gbdraw/
python -m pytest tests/ -v -m "not slow"
python -m build
```

## 8. Acceptance matrix

| Case | Required result | Primary evidence |
|---|---|---|
| Empty overrides, tolerance 0 | current assignments and reference SVGs unchanged | core + output comparison |
| Resolver OFF + fixed directional | fixed feature moves; Auto stays nominal and may overlap | core + Circular/Linear SVG |
| Resolver ON + B Main | B remains Main; previous Auto A moves | core + both modes |
| Main with resolver OFF | same nominal geometry; override survives and becomes active precedence when resolver turns on | request/session + core |
| Auto reset | override absent; fresh deterministic allocation | core + Web history |
| Fixed-fixed conflict | direct error names both identities and lane | core + Web failed-candidate flow |
| Touching intervals | no conflict at tolerance 0 | unit |
| 1/2 bp overlap | conflict only when intersection exceeds tolerance | unit + API/CLI/Web parity |
| Circular origin wrap | split intersections are summed correctly | unit + render |
| Multipart/intron | whole feature moves together; envelope policy remains | core + visual |
| Hidden/underlay | override persists but does not reserve foreground lane | core + session |
| Circular resolved `lane_direction=split` | Main/Outward/Inward lane 1, labels and canvas correct | built-in + custom-slot render/visual |
| Linear resolved feature slot `side=overlay` | Main/Above/Below lane 1, labels/corridor correct | built-in + custom-slot render/visual |
| One-sided/separate direction | explicit validation error | typed + CLI + Web |
| Separate strands Auto/Main | resolver OFF/ON both work | Circular/Linear tests |
| Rotation/reverse/crop-retained | exact same source feature receives override | identity + request replay |
| Source-known crop-excluded | override stays dormant and reactivates when crop includes it | identity + request replay |
| Duplicate base hash | source-index-disambiguated features independently targetable | catalog + placement tests |
| Multi-record | record-scoped override cannot leak to another record | core + request/session |
| Duplicate record IDs | distinct canonical record keys receive only their own overrides | request plan + render |
| CLI/Python/Web parity | equivalent request produces equivalent resolved lane geometry | integration |
| Save/load/re-render | sparse overrides and tolerance survive current writer round-trip | Python + Web session |
| Undo/redo | requested map and stale state restore atomically | Web unit/browser |
| Invalid draft | committed preview, catalog, selection, export remain | browser |
| Render-only | no LOSAT run or raw evidence invalidation | Web orchestration test |

## 9. Evidence ledger

future implementation sessionは、各phase完了時に次を更新する。gateが通る前に`complete`へ変更しない。

| Phase | Status | Behavior implemented | Evidence | Deviations | Remaining risk |
|---|---|---|---|---|---|
| 0 Baseline/schema | not started | — | — | — | C Phase 3 integration未確認 |
| 1 Model/identity | not started | — | — | — | shared identity未実装 |
| 2 Tolerance | not started | — | — | — | scope-cut可能 |
| 3 Allocator | not started | — | — | — | core precedence未実装 |
| 4 Circular | not started | — | — | — | label/layout未実装 |
| 5 Linear | not started | — | — | — | label/corridor未実装 |
| 6 API/CLI/session | not started | — | — | — | public/persisted contract未実装 |
| 7 Web | not started | — | — | — | editor/history/transaction未実装 |
| 8 Docs/QA | not started | — | — | — | broad/visual evidenceなし |

Evidence entryは最低限次の形にする。

```text
Phase:
Behavior implemented:
Evidence: exact command, result, and inspected artifact
Deviations from plan:
Remaining risk:
```

## 10. Completion rule

Work package Dをcompleteにできるのは次がすべて成立したときだけである。

1. roadmapのrelease-blocking manual placement contractを全surfaceで満たす。
2. requested intentとresolved geometryに複数のauthorityがない。
3. default outputとerror typeの回帰がない。
4. Circular/Linearのfeature、label、leader、track reservation、canvas/corridorを実物で確認している。
5. request/session versionとpositive compatibility fixtureがroadmap、Python、Web、docsで一致し、Work package Cのrecord-display fieldsとDのfieldsが同じcomplete writerに含まれる。
6. CLI/Python/Webの再現手段が同じsemantic overrideへ収束する。
7. invalid candidateとundo/redoがlast-successful-render contractを守る。
8. focused、browser、reference comparison、not-slow、lint、build gatesの結果をledgerへ記録している。
9. toleranceをscope-cutした場合は、roadmap、public claims、UI/CLI/API、testsから同じ変更で外し、未完了をcomplete扱いしていない。

実装完了時にはrepository guidanceに従い、English commit titleと短いsummaryを提案する。
