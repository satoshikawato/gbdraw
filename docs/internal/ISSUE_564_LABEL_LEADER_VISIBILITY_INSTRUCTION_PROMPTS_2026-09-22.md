# Issue #564 Feature label・leader line可視性修正 — INSTRUCTION PROMPTS

このファイルは、gbdraw Issue #564を修正する各実装セッションの開始指示である。読者は過去の
会話を知らないものとする。

問題、独立再現、設計、Product／architectureの扱い、受入条件、実施記録の管理先は
[総合計画書](ISSUE_564_LABEL_LEADER_VISIBILITY_MASTER_PLAN_2026-09-22.md)である。
S0 → S1 → S2 → S3の順に使用する。各code block全体を新しい担当者へのinstructionとして
そのまま渡せる。

全sessionの固定実装ブランチは
`fix/issue-564-label-leader-visibility-20260922`である。計画作成時のbaseは
`origin/dev`の`11aae136694a4433cabc68c0dae31edf77222740`。別のruntime branchを作らず、各sessionを
同じbranchへ積み上げる。無関係な変更と未追跡fileは保持し、本件へ含めない。

## S0 — baseline、Product／architecture preflight、test contract

```text
gbdraw Issue #564の実装前baseline、Product Impact、persisted Result互換性、architecture
classificationを確定し、失敗test contractを準備してください。runtime production codeは
変更しません。過去の会話や/tmpの一時ファイルは前提にしません。

対象Issue:
https://github.com/satoshikawato/gbdraw/issues/564

問題:
Feature editorでexternal labelを個別にHideすると、文字だけがdisplay=noneになり、別SVG
elementとして描かれたleader lineが見えるまま残る。文字と全leader segmentは一つの
logical LabelVisualUnitとして動く必要がある。Auto reflowのオン／オフ、同名label、
multi-record、Result、exportでもassociationを失ってはならない。

最初に全て読む:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/ISSUE_564_LABEL_LEADER_VISIBILITY_MASTER_PLAN_2026-09-22.md
5. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
6. docs/internal/PRODUCT_IMPACT_RATCHET.md
7. docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md
8. docs/internal/WEB_CHANGE_POLICY.md
9. docs/SESSION_COMPATIBILITY.md
10. docs/REFERENCE/session-and-request-compatibility.md
11. docs/SVG_SEMANTIC_HOOKS.md
12. docs/internal/SESSION05A7_E_MULTIPART_LABEL_REMEDIATION.md

開始確認:
- git fetch origin後、branchがfix/issue-564-label-leader-visibility-20260922であること、HEAD、
  upstream、origin/devとのancestry、worktreeを確認する。別のruntime branchを作らない。
- 計画作成時のbaseは11aae136694a4433cabc68c0dae31edf77222740だが、最新と仮定しない。
- 無関係な変更、未追跡file、別worktreeを保持する。特に計画書に記録された既存dirty fileを
  stage、restore、書換えしない。
- GitHub Issue #564の本文、comments、state、updated_atを取得し、総合計画書の前提との差を読む。

既知の独立再現:
- Lambda Linear Sessionでrendered feature f40a677cbのhypothetical proteinをAuto reflow offで
  Hideすると、textはdisplay=noneとpreview markerを持つがleader lineはvisibleのまま。
- HmmtDNA Circular Sessionでfb8ff22d9のtRNA-PheをHideすると、textだけが消え二本のleaderが残る。
- canonical labelVisibilityOverridesはoffへ更新され、壊れた部分状態がcurrent Resultとexportへ
  serializeされる。
この記録を修正後の合格証拠にせず、最新sourceで再現する。

作業A — baselineとauthority:
1. tracked Gallery Sessionを使い、LinearとCircularのHide/Restoreを実ブラウザで再現する。
   text、全leader、override map、current Resultを別々に記録する。Auto reflow offを必須とする。
2. Auto reflow onはimmediate stateとrerender完了後を別checkpointで測定する。完了待機条件を
   DOM mutationやResult generation sequenceで明示し、古いResultを新しい結果と誤認しない。
3. latest baseのProduct Impact map、active BD/PD、Web live-edit contract、Issue本文を検索する。
   desired outcomeが一つに決まる場合はIMPLEMENT_EXISTING_AUTHORITYとし、新しいProduct Decisionを
   作らない。他のproduct-valid outcomeが残る場合だけ影響箇所を止め、Decision Packを作る。
4. renderer、label-actions、candidate replay、sanitizer、serialization、export、session Loadの
   actual call pathを確認し、総合計画書第5節を最新化する。

作業B — binding contract:
1. repository-wideでdata-gbdraw-label-binding-schemaまたは同義attributeの衝突を検索する。
2. 一つのtextと0..N leader segmentsが同じdata-label-feature-idを持ち、textがcomplete binding
   schema=1を宣言するcontractを、最小のownerで実現できるか確認する。
3. Linearでは既存build_linear_feature_dom_indexまたはneutralな既存identity ownerを再利用し、
   record qualificationとduplicate source ordinalをlabel専用に再実装しない方針を確定する。
4. Circularの既存feature_idとmulti-record rebindがlineにも適用できることを確認する。
5. root-wide markerを使う必要がない限りper-text markerを選ぶ。別案を選ぶ場合は、部分的な
   copied subtreeでも完全性を証明でき、change pointが小さい理由を実施記録へ書く。

作業C — persisted metadata-free Result:
1. main first-parentまたはrelease tagでmetadata-free saved Resultが実在したか確認する。
2. tracked Session 42をpositive fixtureとして固定できるか確認する。
3. incomplete SVGでtextだけを隠すこと、DOM adjacency、label text、座標近接でlineを推測することを
   不採用にする。
4. canonical overrideを先にcommitし、direct projectionを行わず、既存force rerenderへ一度だけ
   送る推奨案をWeb contractと照合する。failureでは旧text/lineを共に保持しerrorを報告する。
5. この分岐が既存canonical fallbackの再利用でdelta(CB)=0か、新しいpersisted compatibility
   pathでdelta(CB)>0かをArchitecture Ratchetに従って分類する。
6. delta(CB)>0ならstable ID、完全なOE/PE/CB setsと算術、positive fixture、first-parent/release
   evidence、removal condition、maintainer decisionを用意するまでdependent runtimeを開始しない。
   checkerや受入条件を緩めない。

test-first deliverables:
- 総合計画書LV-01〜LV-14を既存test ownerへmappingする。
- current codeで失敗する最小Linear external caseとCircular multi-segment caseを追加する。
- equal label text、multi-record／duplicate、embedded control、metadata-free Result、rerender failureの
  fixture strategyを確定する。
- implementation式をtestへコピーせず、rendered SVGとuser journeyのobservable resultをassertする。
- tracked referenceやGallery artifactを通常testで上書きしない。

設計原則:
- SOLID: Product outcome、identity、geometry、SVG projection、Web visibility、reflowのownerを分ける。
- KISS: 既存feature identity、override、rerender、serializationを使う。
- DRY: mode別identityやvisibility helperを計画しない。
- YAGNI: new Session field、watcher、state store、public SVG API、layout redesignを追加しない。

このsessionの境界:
baseline、classification、contract、失敗test、実施記録まで。renderer、Web runtimeのproduction
修正は行わない。未解決のProductまたはarchitecture decisionがある場合は、その依存部分だけを
止め、独立したtest準備と証拠収集を終える。

終了条件:
- 総合計画書第13節S0へbase/branch、Issue state、reproduction、classification、binding contract、
  compatibility disposition、追加test、command/result、S1開始可否を記録する。
- production、tests、docsのdiffを別々にreviewする。
- repository guidanceに従う英語のproposed commit titleと短いsummaryを示す。
- push、PR、merge、tag、deployは、その時点で明示的に許可されない限り行わない。
```

## S1 — renderer-owned complete label binding

```text
gbdraw Issue #564について、LinearとCircular rendererがlabel textと全leader segmentへexact
feature identityを投影するcontractを実装し、focused Python verificationを完了してください。
過去の会話は前提にしません。Webのvisibility runtimeはこのsessionでは変更しません。

最初に全て読む:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/ISSUE_564_LABEL_LEADER_VISIBILITY_MASTER_PLAN_2026-09-22.md
5. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
6. docs/SVG_SEMANTIC_HOOKS.md
7. 総合計画書第13節のS0実施記録と、S0が追加した失敗test

厳格な開始条件:
- branchがfix/issue-564-label-leader-visibility-20260922であること、HEAD、upstream、base ancestry、
  worktreeを確認する。別branchを作らず、無関係な変更を保持する。
- S0でcomplete binding attribute名、Product classification、architecture routeが確定していること。
- persisted compatibilityにmaintainer decisionが必要でも、rendererのfresh-output contractが独立して
  実装可能ならそこだけ進める。未承認のWeb compatibility branchは実装しない。

実装するcontract:
- LabelVisualUnitは一つのtextと0..N leader segmentsから成る。
- 全partが同じexact data-label-feature-idを持つ。
- textがcomplete binding schema version 1を宣言する。
- embedded labelはtextだけで完全。存在しないleaderを作らない。
- IDはlabel text、DOM adjacency、座標から推測せず、renderer-owned feature identityから得る。

Linear実装:
1. prepare_label_list_linear()がsource identity lookupに必要な既存source_feature_indexまたは同等の
   opaque keyをlabel entryへ保持する。placement計算とは混ぜない。
2. build_linear_feature_dom_index()相当の既存exact resolutionを再利用する。record_index/count、
   stable hash、duplicate ordinalの第二式をlabel moduleへ書かない。
3. 現在のidentity ownerがpairwise-match固有fileにあるため再利用が不自然なら、二つ以上の現実の
   consumerが使う最小private resolverだけをneutral ownerへ抽出する。新しいpublic registryや
   mutable cacheを作らず、旧計算を同じ変更で削除する。
4. precalculated label pathとdirect pathが同じIDを受け取るよう、assemble→builder→SeqRecordGroupの
   既存data flowを最小限延長する。
5. Linear LabelDrawerはtextへidentityとcomplete markerを付ける。
6. SeqRecordGroupはexplicit leader_line分岐と通常external分岐の全lineへ同じidentityを付ける。
7. leader→feature→textのelement order、line coordinates、stroke、text geometryを変えない。

Circular実装:
1. labels.pyが既に持つfeature_idを使い、Circular LabelsGroupの全leader segmentへ
   data-label-feature-idを付ける。
2. Circular LabelDrawerのtextへcomplete binding markerを付ける。
3. copied multi-record canvasの既存rebindがtextとlineを同じrecord-qualified IDへ変えることを守る。
4. horizontal/radial、multipart、embedded/externalで同じcontractにする。Circular専用Web処理や
   identity resolverを増やさない。

svgwriteとoutput:
- custom attributeを付けるelementだけ、既存precedentに従って必要最小限のdebug=False/Parameterを
  使用する。global validationを無効化しない。
- tests/test_output_comparison.pyのmetadata normalizerへcomplete markerを加える必要がある場合、
  非visual attributeだけを除外する。geometry/style/text差分を隠す正規表現にしない。
- exact replay hashへ到達する場合はfull current hashと、追加属性だけを除いたhistorical hashの
  dual oracleを使う。tracked referenceを手編集しない。

必須test:
- Linear single/multi-record、equal stable hash duplicate、source ordinal、multipart。
- Linear embedded、one-line external、explicit leader_line branch。
- precalculated/direct label preparation parity。
- Circular horizontal/radial、one/two segment、embedded/external、copied multi-record。
- textと全lineのID一致、feature partsへのexact resolution、unrelated IDの分離。
- element orderとgeometry attribute不変、strict serialization成功。
- existing tests/test_circular_label_identity.pyを拡張し、必要ならfocused Linear fileを追加する。

focused verification:
pytest <new-or-existing-linear-label-identity-tests> -v
pytest tests/test_circular_label_identity.py -v
pytest tests/test_output_comparison.py::TestOutputComparison -v
ruff check gbdraw/
git diff --check

設計review:
- SOLID: identity resolverはgeometryやWeb visibilityを所有していない。
- KISS: additive metadataだけでlayer topologyを維持する。
- DRY: featureとlabelが同じexact identity resolutionを使う。
- YAGNI: layout algorithm、new schema store、public API、future label frameworkを追加しない。

このsessionの境界:
renderer contractとPython testsまで。label-actions.jsのbehavior、legacy saved Result fallback、
browser artifact regenerationはS2/S3へ残す。in-scope test failureは修正し、timeoutやassertionを
弱めない。

終了条件:
- 総合計画書LV-05〜LV-07、LV-12〜LV-14のrenderer側証拠を満たす。
- 第13節S1へHEAD、files、identity data flow、owner/path evidence、commands/results、未実施項目を
  記録する。
- productionとtestsのdiffを別々にreviewする。
- 英語のproposed commit titleと短いsummaryを示す。
- push、PR、merge、tag、deployは明示的な許可がない限り行わない。
```

## S2 — Web atomic visibility、Result、metadata-free saved Result

```text
gbdraw Issue #564について、rendererのcomplete binding contractを使ってFeature editorの
visibility direct editをLabelVisualUnit全体へ適用し、Result replayと認められたlegacy pathを
実装してください。過去の会話は前提にしません。

最初に全て読む:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/ISSUE_564_LABEL_LEADER_VISIBILITY_MASTER_PLAN_2026-09-22.md
5. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
6. docs/internal/PRODUCT_IMPACT_RATCHET.md
7. docs/internal/WEB_CHANGE_POLICY.md
8. docs/SESSION_COMPATIBILITY.md
9. 総合計画書第13節のS0/S1記録と現在のgit diff

開始条件:
- branchがfix/issue-564-label-leader-visibility-20260922であること、HEAD、upstream、base ancestry、
  worktreeを確認し、別branchを作らない。無関係な変更を保持する。
- S1のrenderer outputでtextと全leaderが同じexact IDを持ち、textのcomplete binding markerが
  Python testsでpassingであること。
- S0のProduct/architecture classificationを読む。metadata-free saved Resultが新compatibility
  pathと分類され、必要なexception/maintainer decisionがない場合、その分岐は実装せず停止する。
  current complete SVGのdirect pathと独立testは進める。

実装:
1. svg-sanitization.jsへS0で確定したcomplete binding markerをallowlistし、既存identityとpreview
   markerを保持する。sanitize policy全体を緩めない。
2. feature-editor/label-actions.jsの既存owner内へ、小さなprivate resolverを一つ置く。
   target textのnon-empty data-label-feature-idとsupported binding schemaを検証し、SVG root内の
   exact escaped identityを持つ全partを返す。target textが集合に含まれることも確認する。
3. 一つのvisibility mutation helperでoff/on/defaultを全partへ適用する。offはdisplay=noneと
   既存preview marker、restoreはこのpreviewが所有するdisplay/markerだけを戻す。
4. applyDirectVisibilityToCurrentSvg()とapplyStoredVisibilityOverridesToSvg()を同じresolver/helperへ
   収束させる。Linear/Circular、text/line、direct/replayの別実装を作らない。
5. canonical labelVisibilityOverridesをaction内で先にcommitする既存順序を維持する。mounted SVGを
   完全更新した後、current Resultを一度だけserializeする。
6. complete current SVGではvisibility-only editのためにforce rerenderしない。Auto reflowオン時の
   optional geometry reflowは従来どおり別責務とする。
7. S0で許可された場合だけ、complete markerがないsaved Resultではpartial hideを行わず、既存
   force rerenderを一回要求する。new queue、retry loop、reactive compatibility flagを追加しない。
8. legacy rerender failureでは旧SVGのtext/leaderを共に保持し、canonical overrideを失わず、既存
   labelReflowLastError経路で報告する。成功後は新complete SVGにoverrideが適用される。
9. serializeCleanSvg/exportへlabel固有の後処理を追加しない。正しいlive SVGを既存clone経路へ渡す。

禁止事項:
- DOM sibling/parent順、label text、x/y/d、線の近接によるassociation。
- 全visibility editでreflowをforceすること。
- Auto reflow設定を黙って変更すること。
- labelVisibilityOverridesのmirror state、mode別override、新Session field。
- textとlineを一groupへ移してrenderer layer順を変えること。
- binding不明時にtextだけを隠すbest-effort fallback。

必須contract tests:
- exact IDのtextと0/1/2 lineをoff/on/defaultでatomic更新。
- 同じ文字列・別ID、似たID、特殊文字を含むID、multi-record IDの隔離。
- unrelated feature path/legend/annotationへ非干渉。
- missing/unsupported marker、empty ID、missing targetがfail closed。
- stored override replayとdirect editが同じhelperを使う。
- Result serialization一回、preview marker cleanup、sanitization round-trip。
- complete SVGではforce rerender invocation 0。
- legacy positive fixtureでは認められたrefresh invocation 1、partial DOM mutation 0。
- injected rerender failureでold visual complete、override retained、error reported。
- candidate rerenderへoverrideを再適用し、hidden unitの全partがhidden。

既存browser testの修正:
tests/web/right-drawer.playwright.spec.jsのindividual Feature, Label, and Legend edit testは、textの
displayだけでなく、選択labelの全leader、Result content、unrelated labelをassertする。既存の
Feature/Legend assertionsを削らない。より小さい専用specが責務上明確なら追加し、同じjourneyを
二重所有しない。

focused verification:
node --test <affected label/editor/sanitization/result contract tests>
npx playwright test <small focused direct-edit spec> --workers=1 --retries=0
node --test tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs
git diff --check

Node @playwright/testが使えなければPython Playwrightで等価の実操作を行う。Chromium sandbox
errorは必要な権限で同じcheckを再実行する。

設計review:
- SOLID: Web helperはsemantic bindingを消費するだけでidentity/geometryを再計算しない。
- KISS: 一resolver、一mutation helper、既存serialization/reflowだけを使う。
- DRY: direct/replay、Linear/Circularで同じunit operationを使う。
- YAGNI: state、watcher、export patch、migration framework、always-reflowを追加しない。

終了条件:
- 総合計画書LV-01〜LV-04、LV-08〜LV-11、LV-14のWeb側focused evidenceを満たす。
- architecture ordinary evidence、または承認済みexception evidenceを正確に記録する。
- 第13節S2へHEAD、files、commands/results、serialization count、legacy/failure結果、S3の残作業を
  記載する。
- productionとtestsのdiffを別々にreviewする。
- 英語のproposed commit titleと短いsummaryを示す。
- push、PR、merge、tag、deployは明示的な許可がない限り行わない。
```

## S3 — real browser、保存・export、artifact、最終受入

```text
gbdraw Issue #564の修正について、現在sourceから作ったbrowser wheelで実browser journey、
Save/Load、reflow/regeneration、export、artifact差分、最終gateを完了してください。過去の会話は
前提にしません。

最初に全て読む:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/ISSUE_564_LABEL_LEADER_VISIBILITY_MASTER_PLAN_2026-09-22.md
5. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
6. docs/internal/PRODUCT_IMPACT_RATCHET.md
7. docs/internal/WEB_CHANGE_POLICY.md
8. docs/SESSION_COMPATIBILITY.md
9. docs/REFERENCE/session-and-request-compatibility.md
10. 総合計画書第13節のS0〜S2記録と現在の全diff

開始条件:
- branchがfix/issue-564-label-leader-visibility-20260922であること、HEAD、upstream、base ancestry、
  worktreeを確認し、別branchを作らない。無関係な変更を保持する。
- S1 renderer identity testsとS2 Web contract testsがpassingであること。失敗や未承認の
  compatibility pathがあれば、依存journeyを合格扱いにしない。
- current sourceからpython tools/prepare_browser_wheel.pyでwheelを準備する。wheelとcandidate source
  の一致を記録し、generated wheelをcommitしない。

real-browser acceptance:
1. Linear fresh/current output:
   - external label、Auto reflow off、Hide直後にtextと全leaderがhidden。
   - current Resultも同じで、非対象labelの座標/transform/leader geometryはbyte-equivalent。
   - reflow invocationは0。
   - Restoreで全partとpreview markerが正しく戻る。
2. Linear Auto reflow on:
   - action直後のdirect stateをrerender開始前に確認する。
   - rerender完了後も同じoverrideで全partがhidden。
   - stale completion、失敗、retryでhidden intentを失わない。
3. embedded control:
   - Hide/Restoreでき、leaderは0のまま。
4. identity isolation:
   - repeated label names、multi-record、可能ならequal stable hash duplicateで対象instanceだけ変わる。
5. Circular:
   - HmmtDNA等の二本leader labelでHide/Restore、Result、post-reflowを確認する。
6. persisted metadata-free Result:
   - S0で承認されたcontractどおり、tracked Session 42をLoadして最初のeditがpartial orphanを作らず、
     canonical refreshは一回だけで、成功後はdirect pathへ移る。
   - failure injectionではold visualがcomplete、overrideがretained、errorがvisible。
7. lifecycle:
   - Hide後にSave Sessionし、新しいpageでLoadする。
   - reflow、full Generate、Historyの対象操作、exportを行う。
   - mounted SVG、current Result、saved Result、downloaded SVGをparseし、visible orphanが0である。

browser evidence:
- textだけでなく、対象identityを持つ全elementsを数える。
- display/preview marker、Result generation、override map、reflow invocation countを記録する。
- unrelated labelのgeometry snapshotを比較する。
- 実download SVGをDOM文字列とは別にparseする。
- readable scaleのLinearとCircularを目視し、orphan、layout shift、layer順、clippingがないことを確認する。
- Node PlaywrightがなければPython Playwrightを使う。Chromium Operation not permittedは必要な
  sandbox権限で同じcheckを再実行する。

artifact handling:
1. additive SVG metadataでtracked Gallery/session/recipe artifactsが変わる範囲をproducerから特定する。
2. 必要なものだけ既存owner toolで再生成する。JSON/SVGを手編集しない。
3. 新旧SVG treeから許可したidentity/binding attributesだけを除いた比較を行い、geometry、style、
   text、element orderが同じことを証明する。
4. reference outputsを更新する必要がある場合は正式な--update-reference-outputs手順を使い、差分を
   reviewし、TestOutputComparisonを再実行する。
5. examples/gbdraw_social_preview.pngには触れない。
6. Gallery tutorial screenshotや説明を変更する場合だけweb-gallery-screenshot-maintenance skillを読み、
   実browserから再取得する。不要なら画像を変更しない。

minimum gates:
pytest <all affected label identity/renderer tests> -v
pytest tests/test_output_comparison.py::TestOutputComparison -v
node --test <all affected Web contract tests>
npx playwright test <all focused Issue-564 specs> --workers=1 --retries=0
node --test tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs
ruff check gbdraw/
git diff --check

変更の到達範囲に応じてrepositoryのfast test selection、recipe/Gallery checks、exact replay checksを
追加する。30分未満でlong testをtimeout扱いせず、進行を監視する。valid evidenceはsource/input/
environment/acceptanceが変わっていない限り再利用し、最終変更後に無効になったcheckだけを再実行する。

最終review:
- production diff: owner、identity flow、visibility path、superseded text-only pathを確認。
- tests diff: user-visible assertionsが実装式のcopyでないこと、negative controlsを確認。
- docs diff: current behaviorとcompatibility dispositionだけを記述していることを確認。
- generated diff: metadata-onlyかをtree/geometryと目視の両方で確認。
- SOLID/KISS/DRY/YAGNI: 一owner、一path、一state、最小helper、不要なschema/module/fallbackなしを確認。
- Architecture Ratchet: ordinary concise evidence、または承認済みexception packetを最終HEADに合わせる。
- Product Impact: selected outcome、journey、checkpoints、failure continuationを最終結果と照合する。

完了条件:
- 総合計画書LV-01〜LV-14が全てpass、または未達理由が明記されている。
- 第13節S3へfinal HEAD、wheel identity、全commands/results、browser evidence、artifact review、
  known limitations、rollbackを記録する。
- in-scope failureを残して完了報告せず、受入条件やcheckerを弱めない。
- 全条件成立時だけ修正完了と報告し、英語のproposed commit titleと短いsummaryを示す。
- push、PR、merge、tag、release、deployは、その時点で明示的に許可されない限り行わない。
```
