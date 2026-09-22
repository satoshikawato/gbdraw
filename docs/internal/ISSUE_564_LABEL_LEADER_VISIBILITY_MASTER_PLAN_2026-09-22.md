# Issue #564 Feature label・leader line可視性修正 — 総合計画書

状態: 計画策定済み、runtime未実装。

固定実装ブランチ: `fix/issue-564-label-leader-visibility-20260922`

計画commit時のbase: `origin/dev` / `11aae136694a4433cabc68c0dae31edf77222740`
（2026-09-22 JST）。本書と別紙promptはこの固定ブランチへcommitする。S0〜S3の実装も
同じブランチへ積み上げ、別のruntime branchを作らない。実装開始時には最新のremote状態を
確認するが、この記録済みSHAを無条件に最新とみなしたり、公開済みブランチを無断で
rebase／force-pushしたりしない。

対象Issue:

- [#564 Leader lines remain visible after hiding feature labels in linear views](https://github.com/satoshikawato/gbdraw/issues/564)

各セッションへそのまま渡せる開始指示は、別紙
[INSTRUCTION PROMPTS](ISSUE_564_LABEL_LEADER_VISIBILITY_INSTRUCTION_PROMPTS_2026-09-22.md)
にある。本書は、過去の会話や`/tmp`の一時ファイルを知らない参加者が、再現結果、原因、
設計判断、受入条件、実施状況を一か所から復元するための管理文書である。

## 1. 目的と完成像

gbdraw Web版のFeature editorでは、生成済みSVG上のfeature labelを個別に表示・非表示へ
切り替えられる。外側labelは、文字`<text>`と一本または複数本のleader `<line>`から成る。
現行処理は文字だけを隠すため、線が空の位置を指す「orphaned leader line」が残る。

完成後は、次の一つの論理単位を扱う。

```text
LabelVisualUnit = one label text + zero or more leader-line segments
```

`LabelVisualUnit`は設計上の概念であり、新しい公開classや永続modelを必須にするものではない。
rendererが既存のrecord-qualified feature identityを各構成要素へ投影し、Web editorはその
identityだけを使って単位全体の可視性を同期する。

完成条件は次のとおりである。

- Hide／Restoreは文字と全leader segmentへ同期して即時反映される。
- Auto reflowがオフでも正しく、他labelの座標を変えない。
- Auto reflowがオンでも、再配置前後で同じ可視性になる。
- embedded labelは従来どおりleaderを持たない。
- 同名label、複数record、同一hashのfeatureでも別単位を誤操作しない。
- current Result、Save／Load、再生成、exportで可視性の意味が一致する。
- LinearとCircularで、同じ一般契約を別々に実装しない。
- label配置、衝突解決、layer順、feature geometryは変更しない。

## 2. 用語

| 用語 | 本書での意味 |
| --- | --- |
| feature identity | 一つのrendered feature instanceを識別する、recordとduplicateを考慮した既存ID。表示文字列ではない。 |
| stable feature ID | 生物学的feature geometry等から得る基礎ID。複数recordやduplicateだけでは一意とは限らない。 |
| rendered feature ID | record instanceと必要なsource ordinalで修飾した、現在のSVG内のfeature identity。 |
| label text | feature labelを表示するSVGの`<text>`要素。 |
| leader segment | external labelとfeatureを結ぶSVGの`<line>`要素。一つのlabelが複数segmentを持ち得る。 |
| complete binding | label textと、そのlabelに属する全leader segmentが同じexact identityを持つことをrendererが保証した状態。 |
| direct visibility edit | geometry再計算を待たず、mounted SVGとcurrent Resultを同期更新する既存のlive-edit経路。 |
| authoritative rerender | canonical overrideを入力としてPython rendererからSVGを再構築する既存経路。 |
| saved Result | Session内に保存された、最後に成功した生成結果のSVG。Loadだけでは再生成されない。 |

## 3. 適用規則、ブランチ、変更境界

実装者は着手時に次を全て読む。

- [AGENTS.md](../../AGENTS.md)、[CLAUDE.md](../../CLAUDE.md)、
  [Web CLAUDE.md](../../gbdraw/web/CLAUDE.md)
- [Architecture Fitness Function Ratchet](ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)
- [Product Impact Ratchet](PRODUCT_IMPACT_RATCHET.md)と
  [Product Decision Packet Template](PRODUCT_DECISION_PACKET_TEMPLATE.md)
- [Web Change Policy](WEB_CHANGE_POLICY.md)
- [Session compatibility](../SESSION_COMPATIBILITY.md)と
  [current compatibility reference](../REFERENCE/session-and-request-compatibility.md)
- [SVG semantic hooks](../SVG_SEMANTIC_HOOKS.md)
- [multipart label remediation precedent](SESSION05A7_E_MULTIPART_LABEL_REMEDIATION.md)

固定実装ブランチは`fix/issue-564-label-leader-visibility-20260922`である。各セッションは最初に
branch名、HEAD、upstream、`origin/dev`とのancestry、worktreeを確認し、同じbranchで継続する。
無関係な変更と未追跡ファイルを保持し、本件のcommitへ混ぜない。

本計画の作成・commit・同名remote branchへのpushは依頼者が許可済みである。後続の実装、
追加push、PR作成、merge、tag、release、deployは、本書だけでは許可されない。各作業時点の
明示的な依頼範囲に従う。

## 4. 再現結果

### 4.1 監査対象

Issue本文は2026-09-21作成、2026-09-22の確認時点でopen、comments 0だった。本文は元報告の
正確なversion、Auto reflow設定、操作順が未確定であると明記している。このため、元利用者の
入力を再現したとは主張せず、現行baseで同じfailure mechanismを独立して再現した。

監査対象は計画commit時の`origin/dev`、
`11aae136694a4433cabc68c0dae31edf77222740`である。Chromium/PlaywrightからローカルWebを開き、
repositoryに追跡されたGallery Sessionを読み込んだ。実装後の合格証拠には、実装対象sourceから
準備したbrowser wheelを使う。

### 4.2 Linearでの再現

入力: `gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json`

手順:

1. SessionをLoadし、Feature editorのlabel一覧を同期する。
2. leader lineを持つexternal labelを一件選ぶ。
3. Auto reflowをオフにする。
4. 個別label visibilityをOffにしてApplyする。
5. mounted SVG、current Resultの文字と線を調べる。
6. visibilityをOnへ戻す。

確認した対象はrendered feature ID `f40a677cb`、表示文字
`hypothetical protein`だった。Hide前は文字と次のleaderが表示されていた。

```svg
<line stroke="gray" stroke-width="0.5"
      x1="1972.722774318585" x2="1972.722774318585"
      y1="-28.0" y2="-67.40546875" />
```

Hide後の実測結果:

- text: `display="none"`
- text: `data-gbdraw-label-visibility-preview="off"`
- leader line: 同じDOMへ接続されたまま、`display`なし
- canonical `labelVisibilityOverrides.f40a677cb`: `"off"`
- serialized current Result: textは非表示だが上記leader lineは表示状態

Restoreでは文字が再表示されたが、leaderはHide中も継続して表示されていた。このSessionには
external labelが73件あり、embedded labelだけの偶然ではない。exportはcurrent SVGをcloneする
ため、orphan lineもexportへ残る。

### 4.3 Circularでも確認した同じ一般不具合

入力: `gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json`

rendered feature ID `fb8ff22d9`、表示文字`tRNA-Phe`で同じ操作を行うと、textは非表示になり、
そのlabelの二本のleader segmentは表示されたままだった。

Issue titleはLinearを示すが、原因は共通のWeb visibility処理であり、Circular rendererもtextだけに
identityを付け、leaderには付けていない。Circularを同じ契約へ収束させることは機能拡張ではなく、
一つの不具合を二つのmodeで別実装しないためのDRY修正である。

### 4.4 Auto reflowに関する証拠の限界

Auto reflowオンの最終状態は、監査用wrapperの待機条件が実際の非同期完了と一致せず、独立した
合格証拠を得ていない。静的には、direct edit後の任意reflowがPython出力を置換するため、症状を
隠す場合がある。実装者は「監査時にAuto reflowで直った」と仮定せず、S3でimmediate stateと
post-reflow stateを別々に計測する。

## 5. 原因と既存の正準経路

### 5.1 Web側の直接原因

`gbdraw/web/js/app/feature-editor/label-actions.js`では、
`applyLabelVisibilityPreview(textEl, modeRaw)`が一つの`<text>`だけを更新する。
`applyStoredVisibilityOverridesToSvg()`もeditable textを列挙し、
`applyDirectVisibilityToCurrentSvg()`も`data-label-key`からtext一件だけを解決する。

canonical stateである`labelVisibilityOverrides`の更新は成功している。問題はstate不足ではなく、
direct SVG projectionがlogical unitの一部しか対象にしないことである。したがって新しいvisibility
stateやwatcherを追加してはいけない。

`queueLabelReflow()`はAuto reflowオン、またはforceされた場合だけ既存の再描画を要求する。
可視性同期を常にreflowへ依存させると、Auto reflowオフの意味とvisibility-only editの即時性を
壊す。

### 5.2 Linear rendererの欠落

- `gbdraw/labels/linear.py::prepare_label_list_linear()`は、配置に必要なlabel entryを作るが、
  exact rendered feature identityをentryへ保持しない。
- `gbdraw/render/groups/linear/seq_record.py`はleader `<line>`をfeatureより前に描き、textをfeatureより
  後に描く。これはlayer順を守るための意図された分離である。
- `gbdraw/render/drawers/linear/labels.py`はtextを描くがidentityを出力しない。
- `FeatureDrawer`と`build_linear_feature_dom_index()`には、record countとduplicate source ordinalを
  考慮したexact rendered identityの既存計算がある。

したがってtextとlineを一つのDOM groupへ移動する必要はない。既存のidentity計算結果をlabel
entryへ渡し、離れたlayerにある全partへ投影すればよい。

### 5.3 Circular rendererの欠落

- `gbdraw/labels/circular.py`は既にexact `feature_id`をlabel entryへ投影する。
- `gbdraw/render/drawers/circular/labels.py`はtextへ`data-label-feature-id`を出力する。
- `gbdraw/render/groups/circular/labels.py`は一本または二本のleaderを描くが、同じidentityを
  出力しない。
- Circular multi-record copyでは`gbdraw/api/diagram.py`が既存の
  `data-label-feature-id`をrecord-qualified identityへrebindする。この処理は同じ属性を持つlineにも
  適用できる。

### 5.4 sanitization、Result、export

- `gbdraw/web/js/services/svg-sanitization.js`は既に`data-label-feature-id`と
  `data-gbdraw-label-visibility-preview`をallowlistしている。
- `serializeCleanSvg()`はcurrent SVGをcloneしてResultへ保存し、該当lineを削除しない。
- exportもcurrent SVGのcloneを使うため、direct editの完全性がそのままexportの完全性になる。
- `tests/web/right-drawer.playwright.spec.js`の既存live-edit testはtextだけをassertしており、
  connectorを検出できない。
- `tests/test_output_comparison.py`は`data-label-feature-id`を非visual metadataとして除外してから
  geometry比較する。新しいmarkerも同じ意図で扱い、identity自体は独立testで必須化する。

## 6. Product Impactと互換性preflight

### 6.1 推奨classification

基本挙動は`IMPLEMENT_EXISTING_AUTHORITY`とする。

根拠:

- Issue #564は文字とleaderを一つの論理単位として扱い、Auto reflowから独立して同期する結果を
  明示している。
- [Web CLAUDE](../../gbdraw/web/CLAUDE.md)は、visibility editがcanonical overrideを先にcommitし、
  mounted targetとcurrent Resultをoptional geometry reflowより前に同期することを要求する。
- [SVG semantic hooks](../SVG_SEMANTIC_HOOKS.md)は、compound feature identityをDOM idやgeometryから
  推測しないことを要求する。
- 過去のCircular multipart label remediationは、rendererがidentityを投影し、Webがgeometryから
  推測しない同じ依存方向を採用している。

S0では最新baseのProduct Impact map、active BD/PD、Issue本文・commentsを再確認する。上記と
衝突せず、一つの完全な結果だけが残る場合は新しいProduct Decisionを作らない。saved Resultの
扱いに二つ以上のproduct-valid outcomeが残る場合だけ、影響箇所を止めてDecision Packを作る。

### 6.2 保存済みmetadata-free Result

Sessionは保存済みResultをLoad時に再生成しない。現行のtracked Gallery Session 42は、textと
leaderのcomplete binding markerを持たないpositive fixtureである。このSVGに対し、DOM adjacency、
label text、座標近接からleaderを推測してはならない。

推奨するbounded behavior:

1. canonical visibility overrideを先にcommitする。
2. 対象textにcomplete binding markerがなければ、textだけを部分的に隠さない。
3. 既存のforce rerender経路を一度要求し、canonical overrideを含む現行renderer出力でResultを
   atomicに置換する。
4. 成功後はcomplete bindingがあるため、以後のvisibility-only editはdirect pathを使う。
5. rerender失敗時は、旧SVGの文字とleaderを共に従来表示のまま保ち、canonical overrideを保持し、
   既存UIで失敗を報告する。orphanを作らない。

これは通常のeditごとにreflowを強制する設計ではない。新しいSVGでは常にdirect pathを使い、
Auto reflow設定はgeometry再配置だけを制御する。

ただし、metadata欠落を条件にしたrerenderがArchitecture Ratchet上の新しいpersisted compatibility
pathに該当するかはS0で正式に分類する。既存の「安全なdirect projectionが成立しなければ正規
rerender」というcanonical pathの再利用であり`delta(CB)=0`と証明できる場合は、ordinary
non-increasing evidenceを残す。新規compatibility pathで`delta(CB)>0`となる場合は、次を満たす
architecture exceptionなしに実装しない。

- first-parent `main`またはrelease tagで旧形式の存在を示す。
- tracked Session 42を代表positive fixtureとして固定する。
- stable compatibility IDを付ける。
- 完全なbefore/after OE・PE・CB setと算術を記録する。
- removal conditionを定める。
- exact proposed headについてmaintainer decisionを得る。

新しいSession version、canonical request schema、feature catalog schema、別のmigration storeは
原則不要である。SVG metadataはResult内のadditive internal rendering metadataであり、既存の
Circular identity追加と同じdispositionを第一候補にする。S0の実在履歴が反証した場合は記録を
更新する。

## 7. 修正アーキテクチャ

### 7.1 所有者と依存方向

| 責務 | 正準owner | 修正後の役割 |
| --- | --- | --- |
| rendered feature identity | `gbdraw/features/ids.py`と既存Linear identity index | stable ID、record qualification、duplicate instance qualificationを一度だけ解決する。 |
| label geometry | `gbdraw/labels/linear.py`、`gbdraw/labels/circular.py` | text・leaderの座標を決め、identityをopaque値としてlabel entryへ運ぶ。visibilityは決めない。 |
| SVG part projection | Linear/Circularのlabel drawer/group | textと全leaderへ同じ`data-label-feature-id`を出力し、textでcomplete binding versionを宣言する。 |
| canonical visibility intent | `labelVisibilityOverrides` | `on`／`off`の唯一の編集state。新しいmirrorを作らない。 |
| direct visibility projection | `feature-editor/label-actions.js` | exact identityで一つのLabelVisualUnitを解決し、全partへ一回の操作を適用する。 |
| geometry reflow | 既存`runLabelReflow()`経路 | label配置だけを再計算する。通常のvisibility editを成立させる前提にしない。 |
| Result／export | 既存serialization/export owner | 完全に更新済みのSVGをcloneする。label固有の第二export修正を作らない。 |

依存方向は、renderer-owned semantic identity → SVG metadata → Web projectionとする。WebからPython
geometryへ逆依存せず、rendererにWeb visibility stateを持ち込まない。

### 7.2 SVG binding contract

推奨する最小contract:

```svg
<line data-label-feature-id="feature-instance-id" ... />
<line data-label-feature-id="feature-instance-id" ... />
<text data-label-feature-id="feature-instance-id"
      data-gbdraw-label-binding-schema="1" ...>...</text>
```

- `data-label-feature-id`はtextと、そのlabelに属する全leader segmentで同一。
- `data-gbdraw-label-binding-schema="1"`はtextだけに置き、このlabel unitの全partへexact bindingを
  投影済みであることを示す。
- embedded labelもschema markerを持つが、leaderは0本のまま。
- schema markerは配置、style、visibilityを表さない。
- attribute名はS0でrepository-wide collisionと既存命名を検索し、同じ意味の既存contractがあれば
  それを使う。別名を選ぶ場合も上記semanticsを変えない。
- `data-label-feature-id`をpublic SVG semantic hookへ昇格する必要はない。公開する場合は別の
  compatibility判断を要する。本件ではinternal editor metadataを維持する。

document rootにcomplete markerを置く案は、個別label単位の証明より広く、部分的に古いsubtreeを
copyする経路で過剰な保証になりやすい。まずper-text markerを採用する。S0が全renderer boundaryで
root保証を厳密に証明できた場合だけ、より小さいchange setとの比較を記録して変更してよい。

### 7.3 Linear identity projection

Linearではfeature描画とlabel描画が同じexact rendered IDを使わなければならない。

1. `prepare_label_list_linear()`がlabel対象のsource identityを失わないよう、既存
   `source_feature_index`または同等のopaque lookup keyをlabel entryへ保持する。
2. `build_linear_feature_dom_index()`相当の既存解決結果をrendererへ渡すか、その解決をneutralな
   identity ownerへ抽出する。
3. record index/count、stable hash、同一hashのduplicate ordinalを再計算する第二式をlabel専用に
   書かない。
4. precalculated label pathとdirect label preparation pathが同じidentityになるようtestする。
5. `LabelDrawer`はtextへidentityとbinding schemaを投影する。
6. `SeqRecordGroup`は`leader_line`分岐と通常external分岐の両方で、作成する全lineへ同じidentityを
   投影する。
7. leader → feature → textという既存の描画順を変えない。

既存identity resolverを一般化するために小さなprivate helperを置くことはできるが、将来のlabel
framework、汎用DOM registry、公開classを作らない。現在二つ以上の実利用者が共有する最小interface
だけを抽出する。

### 7.4 Circular identity projection

Circularはlabel entryの`feature_id`を既に持つ。

1. `gbdraw/render/groups/circular/labels.py`が全leader segmentへ同じ
   `data-label-feature-id`を付ける。
2. `gbdraw/render/drawers/circular/labels.py`がtextへbinding schema markerを付ける。
3. multi-record copy/rebindingがtextとlineの両方を同じ値へ変換することを確認する。
4. horizontal／radial、embedded／external、multipartを同じcontractで扱う。

新しいCircular-specific Web処理は作らない。

### 7.5 Webでのatomic visibility projection

`label-actions.js`内の既存ownerへ、狭いprivate helperを置く。

概念的な処理:

```text
resolve complete unit from selected text
  -> require non-empty exact feature identity
  -> require supported complete-binding schema on text
  -> select every element with the exact escaped data-label-feature-id
  -> require selected text to be in that set
  -> apply one visibility operation to all parts
```

`off`では各partに`display="none"`と既存preview markerを付ける。`on`またはpreview付きのdefaultへ
戻す場合は、同じunitのmarkerと`display`を除く。別featureのpartには触れない。selectorは
`CSS.escape()`を使い、文字列連結でunescaped IDを解釈しない。

`applyDirectVisibilityToCurrentSvg()`と`applyStoredVisibilityOverridesToSvg()`は同じunit resolverと
mutation helperを使う。mounted SVGの変更後、既存ownerからcurrent Resultを一度だけserializeする。
mode別helperやexport時の修正を追加しない。

complete markerがない対象では、helperは「変更なし」ではなく「direct projection unavailable」を
呼出元へ明示できる小さな結果を返す。呼出元は第6.2節のpreflight結果に従い、既存force rerenderを
要求する。新しいreactive flag、別queue、retry loopは追加しない。

### 7.6 SOLID・KISS・DRY・YAGNI

- **SRP:** identity解決、label geometry、SVG投影、visibility intent、direct mutation、reflowを
  それぞれ既存ownerに保つ。
- **OCP:** Web editorはLinear/CircularのDOM配置ではなくsemantic identity contractへ依存する。
  leaderが一、二本、または0本でも同じhelperで扱える。
- **LSP:** embedded、horizontal、radial、precalculated labelは、同じbinding contractを満たす限り
  consumerを変えずに置換できる。
- **ISP:** helperはSVG root、target text、exact identity、visibility modeという必要最小限だけを
  受け取り、renderer configやSession全体を要求しない。
- **DIP:** Webはgeometryやtext matchingではなくrendererが出したsemantic contractへ依存する。
- **KISS:** 既存属性、既存override、既存rerender、既存serializationを使い、DOM regroupingをしない。
- **DRY:** identity resolverとvisibility mutationを一つずつにし、Linear/Circularやdirect/replayで
  式を複製しない。
- **YAGNI:** label layout algorithm、placement policy、public SVG API、新Session field、new watcher、
  new state store、常時reflow、汎用component frameworkを追加しない。

## 8. 実装セッション

依存順はS0 → S1 → S2 → S3である。各sessionは本書第13節の実施記録を更新する。入力、source、
環境、受入条件が変わっていないpassing evidenceは再利用し、変更・失敗・未解決箇所だけを再実行する。

### S0 — baseline、authority、compatibility、test contract

- 最新Issue、base、branch、dirty worktreeを確認する。
- LinearとCircularのreproducerを、tracked inputから再実行する。
- Product classificationを確定し、旧saved Resultの扱いを確認する。
- complete binding marker名・owner・sanitizer境界を確定する。
- Architecture Ratchetのordinary／exception分類を行う。
- 受入IDへ対応する失敗testとpositive legacy fixtureを配置する。
- runtime production codeは変更しない。

### S1 — renderer identity contract

- Linearのexact identityをlabel entry、text、全leaderへ投影する。
- Circularの既存identityを全leaderへ投影する。
- textへcomplete binding schema markerを付ける。
- multi-record、duplicate、multipart、embedded、external、precalculated pathをPython testで守る。
- geometry、element order、strict serializationが変わらないことを検証する。

### S2 — Web atomic visibilityと旧Result経路

- sanitizerへbinding markerをallowlistする。
- exact unit resolverと一つのvisibility mutation helperを実装する。
- direct edit、stored override replay、Result serializationを同じhelperへ収束させる。
- S0で認められた場合だけ、metadata-free saved Resultを既存force rerenderへ送る。
- unit、security、Result lifecycle、failure testを完了する。

### S3 — real browser、artifact、final gates

- current sourceからbrowser wheelを準備し、Linear／Circularの実journeyを検証する。
- Auto reflow off/on、Restore、embedded、multi-record、same text、legacy Session、Save/Load、
  regenerate、export、failureを確認する。
- additive metadataで変わるtracked Gallery/session/recipe artifactだけをowner toolから再生成する。
- production、tests、docs、generated artifactsを別々にreviewし、全gateを完了する。

## 9. 受入条件

| ID | 条件 | 主な証拠 |
| --- | --- | --- |
| LV-01 | fresh Linear external labelをAuto reflow offでHideすると、textと全leader segmentが同じtask内で非表示になる。 | browser DOM、current Result |
| LV-02 | 同labelをRestoreすると、textと全leader segmentが戻り、preview markerが残らない。 | browser DOM、serialized SVG |
| LV-03 | Auto reflow offのvisibility-only editで、非対象labelの座標・transform・leader geometryが変わらず、rerenderを起動しない。 | spy／before-after geometry |
| LV-04 | Auto reflow onでimmediate stateとpost-reflow stateが共に完全で、hidden intentを失わない。 | browser two-checkpoint evidence |
| LV-05 | embedded labelのHide／Restoreはtextだけを扱い、leaderを新規作成しない。 | renderer＋browser control |
| LV-06 | 同名label、複数record、同一stable hashのduplicateでも、選択したrendered instanceだけを更新する。 | Python identity＋browser |
| LV-07 | Circularの一本／二本leaderでも同じunit contractが成立する。 | Python＋browser |
| LV-08 | reflowとfull regenerationがcanonical overrideを再適用し、hidden labelのpartを全て非表示にする。 | candidate render／browser |
| LV-09 | current Resultとexported SVGに、hidden label由来のvisible leaderが存在しない。 | serialization／downloaded SVG |
| LV-10 | supported metadata-free saved Resultは推測やpartial hideを行わず、S0で認めた一回のcanonical refreshを行う。失敗時もorphanを作らずerrorを示す。 | legacy positive fixture＋failure injection |
| LV-11 | sanitization、Save／fresh Loadでidentity、binding marker、visibility overrideが必要な範囲で保持される。 | sanitizer／session browser test |
| LV-12 | label配置、混在embedded/external policy、layer順、feature geometry、canvas boundsは変わらない。 | reference comparison＋visual review |
| LV-13 | svgwrite strict validationを破らず、custom attributeを持つ全partをserializeできる。 | focused Python test |
| LV-14 | production codeにDOM adjacency、label text、coordinate proximityによるleader associationがない。 | source review＋negative same-text test |

## 10. 検証計画

### 10.1 Python renderer tests

既存の最も近いtest ownerを優先し、必要なら
`tests/test_linear_label_identity.py`を追加する。

- Linear single-record external: textとlineのidentity一致。
- Linear multi-record: record-qualified ID一致。
- equal stable hash duplicate: source ordinalで分離。
- multipart feature:全feature part、label text、leaderが同じrendered instanceへ解決。
- embedded: schema marker付きtext一件、leader 0件。
- explicit leader branchと通常external branchの両方。
- precalculated labelsとdirect preparationの同一性。
- Circular horizontal／radial、single／two-segment、embedded／external、copied multi-record。
- element orderと全geometry attribute不変。
- output comparison normalizerが新しい非visual markerだけを除外し、identity contract testは属性を
  実際に要求する。
- exact replay hashが影響する場合、full candidate hashと「追加属性だけを除くと既存hashへ戻る」
  dual oracleを使う。既存hashを無根拠に置換しない。

### 10.2 JavaScript contract tests

- exact identityのunitだけを列挙する。
- off／on／defaultで全partをatomicに更新する。
- repeated label textはassociationに使わない。
- unrelated ID、feature path、legend、annotationに触れない。
- invalid／empty ID、unsupported marker、missing targetはfail closed。
- stored override replayも同じhelperを通る。
- serializationは一actionにつき一回。
- sanitizerがidentity、binding marker、preview markerを保持し、危険なattributeは従来どおり除く。
- incomplete old SVGはpartial mutationをしない。

### 10.3 Browser journeys

実装対象sourceから`python tools/prepare_browser_wheel.py`でwheelを用意する。生成wheelはcommitしない。

1. Lambda Linearのfresh generationまたは現行metadata付きSessionでLV-01〜LV-05、LV-08、LV-09。
2. 二record以上、同じlabel textを持つfixtureでLV-06。
3. HmmtDNA Circularの二本leaderでLV-07。
4. tracked metadata-free Session 42でLV-10。rerender invocation countも確認する。
5. Save → 新しいpageでLoad → Hide/Restore → regenerate → exportでLV-08〜LV-11。
6. rerender failureをinjectし、旧visualが部分状態にならず、overrideとerrorが残ることを確認する。
7. DOMの`display`だけでなく、Result文字列と実download SVGをparseする。
8. readable scaleで、hidden positionにorphan lineがなく、unrelated labelが動いていないことを目視する。

Nodeの`@playwright/test`が使えなければPython Playwrightで等価journeyを実行する。Chromiumが
`sandbox_host_linux.cc ... Operation not permitted`で失敗した場合は、同じcheckを必要なsandbox
権限で再実行し、Playwright unavailableとは結論しない。

### 10.4 focused commandsとgates

正確なtest file名はS0で既存ownerに合わせ、実施記録へ保存する。少なくとも次を含む。

```bash
pytest tests/test_circular_label_identity.py -v
pytest tests/test_output_comparison.py::TestOutputComparison -v
node --test <affected label/editor/sanitization contract tests>
npx playwright test <focused Issue-564 specs> --workers=1 --retries=0
node --test tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs
ruff check gbdraw/
git diff --check
```

必要な選択testを追加し、test timeoutや受入assertionを弱めない。通常testで
`tests/reference_outputs/`を更新しない。意図したmetadata-only差分がreference更新を必要とする場合も、
正式な`--update-reference-outputs`手順で生成し、SVG tree／geometry差分をreviewしてからcomparisonを
再実行する。

Gallery Session、recipe、manifestなどに保存されたSVGへadditive metadataを反映する必要がある場合、
既存のowner toolを使う。JSONやSVGを手編集しない。`examples/gbdraw_social_preview.png`には触れない。
Gallery screenshotsやtutorial contentを変更する場合だけ、repository指定のGallery screenshot skillを
読み、実画像を再取得する。本件で不要なら画像文書を増やさない。

## 11. Architecture Fitness evidence

第一候補はordinary non-increasing correctionである。

- Owner before/after: canonical visibility intentは`labelVisibilityOverrides`の一ownerのまま。
- Path before: canonical override → text-only direct projection → Result clone。
- Path after: canonical override → exact LabelVisualUnit direct projection → Result clone。
- Superseded behavior: text-only visibility mutationを同じ変更で削除する。
- Renderer identity owner:既存feature identityを再利用し、label-specific identity ownerを増やさない。
- New module:原則なし。neutral helperの抽出が必要ならprivate decompositionで、state、validation、
  lifecycle、public exportを所有しないことを示す。
- Rollback: runtime/test/artifact変更を一つのcommit seriesとしてrevertし、Session/request schemaを戻す
  migrationは不要。

S0でmetadata-free saved Result分岐が新compatibility pathと分類された場合はordinaryとして進めず、
第6.2節のexception evidenceを作る。architecture checkerやallowlistを緩めて回避しない。

## 12. 非目標、リスク、rollback

### 非目標

- Auto placementで全labelをexternalへ揃えること。
- `External only`の別のplacement不具合を同時に直すこと。
- label collision／packing algorithmの変更。
- leaderの色、太さ、anchor、形状の変更。
- feature ID全体、SVG public API、Session schemaの再設計。
- right drawerやreflow architectureの全面リファクタ。

### 主なリスクと対策

| リスク | 対策 |
| --- | --- |
| duplicate featureが同じlabel unitに混ざる | exact rendered IDを既存resolverから再利用し、source ordinalとmulti-record testを置く。 |
| textとlineをgroup化してlayer順が変わる | group化せず、離れたlayerへ同じidentityを投影する。 |
| old SVGでtextだけを隠す | complete markerがない場合はfail closedし、認められたcanonical refreshだけを使う。 |
| Auto reflowを毎回起動して性能・操作感が悪化 | complete current SVGは常にdirect path。legacy refresh invocation countをtestする。 |
| sanitizerが新markerを落とす | allowlistとsanitize round-trip testを同じsessionで追加する。 |
| nonvisual属性でreference/hashが変わる | geometry comparisonとidentity contractを分け、必要ならdual-hash oracleを使う。 |
| rerender失敗でintentを失う | overrideを先にcommitし、旧complete visualを保持してerrorを報告する。 |

rollbackは、renderer metadata、Web unit projection、関連test/artifactを一緒にrevertする。新しい
Session/request schemaや永続stateを作らないため、data migration rollbackは発生させない。

## 13. 実施記録

各セッションは終了前に一行を更新し、その下へ証拠を追記する。`未着手`を結果で上書きし、
command、pass/fail/skip件数、browser/wheel identity、未解決事項を省略しない。

| Session | 状態 | 対象HEAD／branch | 変更・判断 | 検証 | 次の開始条件 |
| --- | --- | --- | --- | --- | --- |
| Plan | 完了 | `11aae136694a4433cabc68c0dae31edf77222740` / `fix/issue-564-label-leader-visibility-20260922` | 現行Linear/Circularで再現。master planとsession promptsを作成。 | Issue/API、source audit、Chromium独立再現。runtime testは未実施。 | S0で最新base、authority、compatibility分類を再確認。 |
| S0 | 完了 | `3dea935372da3444e394973b410c7c9abf028219` / `fix/issue-564-label-leader-visibility-20260922` | 最新base・Issue・authority・互換性を再確認。per-text schema `data-gbdraw-label-binding-schema="1"`を確定し、renderer failing contractを追加。Productは`IMPLEMENT_EXISTING_AUTHORITY`、architectureはordinary non-increasing、`delta(CB)=0`。 | Chromium再現: Linear 文字hidden／leader 1本visible、Circular 文字hidden／leader 2本visible。current Resultも同数visible、overrideはoff、force rerender 0。Python contractは期待どおり14 fail/1 pass（うちCircular marker 12 fail、Linear contract 2 fail）。 | Product／architecture未解決なし。S1開始可。 |
| S1 | 完了 | `3dea935372da3444e394973b410c7c9abf028219`＋実装candidate / `fix/issue-564-label-leader-visibility-20260922` | Linear/Circularのlabel textと全leader partへexact rendered feature IDを投影し、textへbinding schema markerを追加した。Linear identity ownerはpairwise専用moduleからneutral ownerへ移した。 | renderer focused 123 pass。tracked SVGは新binding attributeを除いてbyte一致し、interactive metadataを持つ一例も既存metadataと新attributeを除いてbyte一致。exact replayとarrow examplesもpass。 | S2開始条件を充足。 |
| S2 | 完了 | 同上 | Webのvisibility投影をcomplete binding単位のatomic helperへ集約。`off`は全partを隠し、`on`/`default`はpreview所有のstyleだけを解除する。incomplete legacy Resultは既存canonical force reflowへ一回だけ送り、失敗時は旧visualを保持する。 | affected Node contract 2 pass、right-drawer browser smoke 1 pass、sanitizer round-trip pass。新state、schema、watcher、queue、retry、lifecycleなし。 | S3開始条件を充足。 |
| S3 | 完了 | 同上 | offline browser wheel、real-browser acceptance、owner artifact、全体gateを完了した。Gallery Sessionはschemaを維持し、binding metadataだけを追加した。 | Issue専用Playwright 5 pass、fast suite 6117 pass/17 skip/11 deselect、architecture 137 pass、`ruff` pass、build成功、Web change budget PASS、`git diff --check` pass。 | 未解決事項なし。commit、push、PR、merge可。 |

### 計画作成時の再現要約

- Linear: `lambda_basic_linear.gbdraw-session.json`、`f40a677cb`でtextのみhidden、leader visible。
- Circular: `HmmtDNA_basic_circular.gbdraw-session.json`、`fb8ff22d9`でtextのみhidden、二本のleader visible。
- current Resultとexportがlive SVGをcloneするため、visible orphanもpersist/exportされる。
- Auto reflowオンの独立完了証拠は未取得。S3の必須checkとする。

### S0実施記録

- `git fetch origin --prune`後も`origin/dev`は
  `11aae136694a4433cabc68c0dae31edf77222740`で、実装branchは計画commitだけを一件積んだ
  `3dea935372da3444e394973b410c7c9abf028219`だった。merge-baseは`origin/dev`と一致し、
  ahead/behindは`0/1`、upstreamは同名remote branch、worktreeは開始時cleanだった。
- Issue #564は`open`、comments 0、`updated_at=2026-09-21T09:41:40Z`で、本文と計画時の
  authority前提に差はなかった。Product Impact mapのactive concernはcanonical render requestと
  saved Session regenerationだけで、本変更がregistered subjectを変更する箇所はない。active
  `BD-###`もない。Issue本文、Web live-edit contract、SVG semantic hookが同じ完全な結果を
  指定するため、procedural classificationは`IMPLEMENT_EXISTING_AUTHORITY`である。
- 最新sourceのtracked Gallery SessionをChromiumでLoadし、Editor同期後にAuto reflowをoffへ固定して
  再現した。Linear `f40a677cb`はtextとResult textがhidden、leader 1本はmounted SVGとcurrent
  Resultの双方でvisibleだった。Circular `fb8ff22d9`もtextだけがhiddenで、leader 2本は双方で
  visibleだった。どちらもcanonical overrideは`off`、force rerender requestは0だった。
- repository-wide collision検索では、計画書以外に`data-gbdraw-label-binding-schema`は存在しなかった。
  contractは全partの`data-label-feature-id`一致と、text単位の
  `data-gbdraw-label-binding-schema="1"`に確定した。document-root markerは採用しない。
- Linearは既存`build_linear_feature_dom_index()`のrecord qualificationとduplicate source ordinalを
  再利用し、label geometryはsource feature ordinalだけをopaqueに運ぶ。Circularは既存
  `feature_id`をlineへ投影し、multi-record copyの既存attribute rebindをそのまま利用する。
- `origin/main`のfirst-parent `4556e04e929a4a85ad28d1833ce7304bd764881c`にtracked
  `lambda_basic_linear` Session 42が実在し、そのsaved Resultにはcomplete markerも
  `data-label-feature-id`もない。これはpositive legacy fixtureとして使用できる。
- metadata-free Resultでdirect projectionを拒否して既存force rerenderへ一度送る処理は、既存の
  「mounted targetへ安全に直接適用できなければauthoritative rerender」というcanonical pathの
  availability条件を厳密化するものである。新しいreader、migrator、schema、state、queue、retry、
  lifecycleは増えず、owner/pathは一つのままなので`delta(OE)=0`、`delta(PE)=0`、
  `delta(CB)=0`のordinary non-increasing correctionと分類した。
- failing contractは`tests/test_linear_label_identity.py`と
  `tests/test_circular_label_identity.py`へ置いた。LV-05〜LV-07、LV-12〜LV-14をrendererで、
  LV-01〜LV-04とLV-08〜LV-11をWeb contract／browser ownerで守る。現行productionに対する
  `pytest tests/test_linear_label_identity.py tests/test_circular_label_identity.py -q`は、test自身の
  setup誤りを修正後、Linearのmarker欠落2件とCircularのmarker欠落12件が期待どおり失敗し、
  unrelated repeated-table test 1件はpassした。tracked reference／Gallery artifactは変更していない。

### S1実施記録

- `gbdraw/render/label_binding.py`をbinding属性の唯一のwriterとし、label textへ
  `data-label-feature-id`と`data-gbdraw-label-binding-schema="1"`、leader lineの全segmentへ同じ
  exact IDを投影した。既存layer順を変えず、textとlineを新しいSVG groupへまとめていない。
- Linearのexact rendered ID ownerを`pairwise_match.py`からneutralな
  `gbdraw/render/groups/linear/feature_identity.py`へ移し、既存のrecord qualification、duplicate source
  ordinal、embedded table ID contractをそのまま再利用した。label geometryからはsource ordinalだけを
  運び、描画時に同じindexで解決する。Circularは既存`feature_id`を全leader segmentへ投影した。
- `python -m pytest tests/test_linear_label_identity.py tests/test_circular_label_identity.py
  tests/test_linear_label_placement.py tests/test_circular_label_placement.py
  tests/test_output_comparison.py::TestOutputComparison -q`は123 passだった。multi-record duplicate、
  embedded tableの0-leader、Circularの2-segment leader、repeated recordを含む。
- owner recipeで更新したstatic SVGは、新binding属性だけを除くとtracked版とbyte一致した。
  interactive BGC Gallery SVGも、既存のembedded interactive metadataと新binding属性を除くとbyte一致し、
  geometry／paint／layer順の差はない。生成中に検出した本件外の5件の既存geometry driftは採用せず、
  対象artifactへ戻した。
- `tests/test_run_info_exact_replay.py`はcomplete bindingを含む新hash、bindingだけを除く旧hash、全label
  identityを除くlegacy hashの三つを明示し、1 passだった。arrow example contractは2 passだった。

### S2実施記録

- `gbdraw/web/js/app/feature-editor/label-actions.js`へcomplete binding resolverとatomic visibility helperを
  置いた。resolverは非空のexact IDとschema marker `1`を要求し、同じexact IDを持つ全partを列挙する。
  IDをCSS selectorへ補間しないため、similar IDや特殊文字で別featureへ波及しない。
- `off`はlabel textと0/1/2本の全leader partへpreview markerと`display: none`を付ける。`on`と
  `default`はpreviewが所有したstyleだけを除去するため、author-authored `display: none`を復元しない。
  direct editとsaved override replayは同じhelperを使い、text＋visibilityの複合editもserializeは一回だけである。
- incomplete current／saved Resultでは部分的に隠さずfail closedとし、S0で認めた既存force reflowを
  一回だけqueueする。force failure時もoverride intentと旧complete visualを保持し、既存error surfaceへ
  通知する。新しいpublic state、Session/request schema、watcher、queue、retry、lifecycleは追加していない。
- sanitizerへschema markerをallowlistし、round-trip contractを追加した。
  `node --test tests/web/feature-label-visual-unit.test.mjs tests/web/svg-sanitization.test.mjs`は2 pass、
  `npx playwright test tests/web/right-drawer.playwright.spec.js`の対象smokeは1 passだった。

### S3実施記録

- `python tools/prepare_browser_wheel.py`でcurrent sourceのoffline wheelを準備した。wheelはgitignoredの
  generated assetでありcommit対象外である。
- `tests/web/issue-564-label-leader-visibility.playwright.spec.js`は5 passだった。fresh Linearで
  off、on/default、Auto reflow on/off、save/load、regeneration、exportを確認し、fresh Circularで
  2-segment leaderを確認した。main first-parent由来のmetadata-free Session 42はcanonical refreshが
  一回だけ発生し、failure injectionでは旧visual、override intent、errorが保持された。同名labelを持つ
  2 recordとembedded tableもexact IDで分離された。
- Linear hidden状態とCircular restored状態のscreenshotを可読scaleで目視し、前者にorphan leaderがなく、
  後者ではtextと2本のleaderが一緒に復元されることを確認した。
- Gallery 3例は旧Session schemaを保持したままcurrent rendererでsource SVGを再生し、owner toolで
  interactive assetとmanifestを更新した。recursive normalized JSON比較ではbinding metadata以外の
  Session差はなく、`python tools/gallery_artifact_manifest.py`はpassした。documented owner recipeの
  Session artifactに含まれる`createdAt`／file `lastModified`は生成timestampとして保持した。
- `python -m pytest tests/ -q -m "not slow"`は6117 pass、17 skip、11 deselectだった。
  `node tests/web/architecture-contracts.test.mjs`相当のarchitecture contractは137 pass、
  `ruff check gbdraw/`はpass、`python -m build`は成功、`node tools/check-web-change-budget.mjs`はPASSだった。
  change budgetのreview表示はregistered Gallery Session 3件の更新によるもので、上記のschema-preserving
  binding-only差分として監査済みである。`git diff --check`もpassし、Product／architectureの未解決事項はない。
- 専用5-case browser specはfull functional acceptanceの所有とし、既存の
  `right-drawer.playwright.spec.js` PR smokeがatomic label/leader投影を継続して守る。CI inventoryは
  上限を変更せず13 PR casesに収まり、`node --test tests/ci/*.test.mjs`は59 passだった。
- metadata-free Gallery Sessionのcommitted SVGとcurrent renderer出力を比較するpublication parityでは、
  この修正の非visual binding 2属性だけを明示的に除外する。比較utilityの既定は変更せず、geometry、
  paint、構造、その他のmetadata差分は従来どおり失敗させる。
