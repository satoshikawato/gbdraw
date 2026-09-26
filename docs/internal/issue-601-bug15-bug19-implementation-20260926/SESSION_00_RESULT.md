# S00 — エラー・regex実装の開始条件

実施日: 2026-09-26。対象は[Issue #601](https://github.com/satoshikawato/gbdraw/issues/601)のBUG-15/19。
runtime、tests、Product authorityを変更せず、S01以降は開始していない。
本結果は調査時点の記録であり、製品契約、owner移管、実装合格を新たに宣言するものではない。
S00の終了には、この資料を含む検証済みHEADの作業ブランチ・devへのpushとremote一致確認が必要。
最終commit SHAとpush結果はコミット後のhandoffで報告する。

## 1. 基準状態とsource inventory

専用clone: `/tmp/gbdraw-issue601-s00-vIoB52/repo`。
共有作業ツリー、他sessionのindex・未コミット変更・環境は操作していない。
開始時の中断された自分のcheckoutだけを復旧し、cleanな以下の状態から監査した。

| 対象 | 確認値 |
| --- | --- |
| 作業branch | fix/issue-601-bug15-bug19 |
| upstream | origin/fix/issue-601-bug15-bug19 |
| 最新remoteからの開始SHA | be4000ede35ea5692d7a88d46c0130ffec992f5f |
| 確認したorigin/dev | af5d942af60353dda199aa487da9152a3576b3fe |
| 先行調査SHA | 2edc00aebc74e01003da643dfc957b513d5dcfe5 |
| origin/main（変更なし） | 4556e04e929a4a85ad28d1833ce7304bd764881c |
| 開始branch対dev | 計画directory内の15追加fileのみ。runtime/test/tool/workflow差分なし |

[inventory生成script](./evidence/s00-inventory.py)と[JSON](./evidence/s00-inventory.json)に、
27 source file、6 authority fileのSHA-256、branch/result SHA、receipt比較を保存した。
gbdraw/tests/tools/.githubは三地点でtree objectがすべて同一。
以後source line番号は開始SHAに対するもの。最新devを取り込む後続sessionは再確認する。

環境はNode 26.8.2、Python 3.13.3、Biopython 1.85、pandas 2.3.0、svgwrite 1.4.3、
pytest 9.0.2、Python Playwright 1.61.0。
Playwright CLIとPython sync APIは利用可能。専用cloneのNode `@playwright/test`は未導入。
S00はbrowser bundleを変更・再生成していない。Python 3.10–3.12、Firefox/WebKit、公開buildは未検証。

## 2. エラーの入口・伝達・表示

| 境界 / 現在のowner | 入口から出口、残る欠落 |
| --- | --- |
| [Python render wrapper](../../../gbdraw/web/js/app/python-helpers.js) run_canonical_request_wrapper（112–199） | native render例外をtype/message/tracebackのJSONにする。文字列化前のcause、code/operation/stage/contextを保持しない |
| [helper producers](../../../gbdraw/web_support/rule_matching.py)、[feature metadata](../../../gbdraw/web_support/feature_metadata.py) | rule helperはnative例外がPyodideへ伝播。一方feature_metadataの618/646等はstr(exc)だけのerrorを返す。helperごとのcatch/文字列化が統一されていない |
| [Worker](../../../gbdraw/web/js/workers/diagram-generation-worker.js) callJsonHelper（479付近） | 一つのlazy Workerからhelper JSONを返す。helper.destroy()のfinally例外が主因を置換し得る。renderのresults.errorとtransportのok:falseは別の出口 |
| 同Worker serializeError（135付近） | 先にnormalizeし、name/message/details/notes/stackを返す。code/operation/stage/contextを落とす |
| 同Worker cleanup（169、308付近） | generation cleanupはJS主因→Python error→destroy→workspaceの優先。workspaceは主因を保護するが自由なcleanup notesを添付する。副因も許可済みの有界識別へ移行が必要 |
| [client](../../../gbdraw/web/js/services/diagram-generation.js) deserializeWorkerError（170付近） | Errorを復元するがstructured fieldを保持しない。resource staging/init/protocol/helper失敗をregex syntaxと分ける |
| [normalizer](../../../gbdraw/web/js/services/error-normalization.js) | summary/detailsのみ。summary 1,000字、detail 4,000字、section 8個。特定data keyやstackを除外しても自由string/notes/stdout/stderr全体の安全性は保証しない。traceback末尾causeが除去され、二度の正規化で型prefixが重なる場合もある |
| [run-analysis](../../../gbdraw/web/js/app/run-analysis.js) | executeCanonicalRenderCandidateはengineErrorを持つが、下記のcommitted経路で返却が落ちる。Circular slotの特別extractorはmessage/stdout/stderr/tracebackを解析して必要情報を救済している |
| [Align Apply](../../../gbdraw/web/js/app/similarity-alignment.js) applyPlan（787–829） | causeがあってもsummaryのみをlocal errorへ渡しDetailsを捨てる。objectをStringへ変換する出口とcauseなしgeneric fallbackがある。action/artifactの現在性、retry review、cancel/staleは維持する |
| [rule action](../../../gbdraw/web/js/app/feature-editor/rule-actions.js) commitPrepared / setSpecificRuleField（49–85） | preparation/History/commitの任意rejectをInvalid ruleと呼ぶ。finallyでaccepted fieldへ無条件復帰。synthetic非syntax失敗でも両挙動を確認 |
| [Label import](../../../gbdraw/web/js/app/feature-editor/label-actions.js)（1000–1108） | parse→Python evaluation→revision/intent/file/catalog/mode/Result/sourceSVG確認→commit。失敗はraw consoleとalert。validation failure時は旧overrideを保護 |
| [Color等のfile watcher](../../../gbdraw/web/js/app/watchers.js)（450–570） | file importを直列化。specific colorはparse/preparation/現在性確認後にrules/legendを更新し、失敗時に前fileへ戻す。catchはraw console/alert。default color、priority、whitelistは各別owner |
| [Export](../../../gbdraw/web/js/services/export.js) / [composition](../../../gbdraw/web/js/app/app-setup.js)（3397付近） | snapshot→SVG/PNG/PDF処理→runExportAction catchでnormalize→Export errorとformat prefixを追加。statusのみ返す。PDF finallyのcleanupが主因を置換し得る |
| [alert UI](../../../gbdraw/web/index.html)（4596付近） | exportにもGeneration Error / Error Occurredの固定見出し。任意Detailsは存在するがCopy diagnosticsは未実装。operationを保つ一つの表示modelへ配線する |

### 全run-analysis error return

9 return文、二つの条件分岐を含む11経路を確認した。返却原因をglobal errorLogから推測しない。

| 行 / 経路 | 現在のoutcome |
| --- | --- |
| 1999 Linear record準備 | reflowはstatusのみ、manualはerrorLogをerrorへ返す |
| 2035 Circular discovery | manualへerrorを返す |
| 2050 Depth input検査 | reflowはstatusのみ、manualはerrorを返す |
| 4665 engine-error reflow | statusのみ。labelReflowLastErrorにsummaryを保持 |
| 4669 engine-error manual | errorを返し、旧transaction状態を復元 |
| 4978 catch reflow | statusのみ |
| 4982 catch manual | 復元後にerrorを返す |
| 5169 committed-candidate engine-error | **errorLogを設定するがerrorを返さない** |
| 5282 committed-candidate catch | **復元後にerrorLogを設定するがerrorを返さない** |

runAnalysis / committed wrapper（5030/5293付近）はgeneratedArtifactCandidateを除く際に
`{status}`へ再構築する。現在のsuccess用途だけで将来のerror保持まで保証しない。
record discoveryやHistory finalization/rollbackにmain try外のrejectもある。
rollbackが失敗した場合まで「旧Resultは復元された」と案内してはならない。

通常のAlign Applyはapp-setup（2775付近）がmanual runAnalysisを注入している。
devに含まれる既存commit 5bcdd4f5で通常経路の原因返却は修正済み。
残存committed-candidate経路はrecord rotation（app-setup 2240付近）等から使われる。
従って5169/5282だけを根拠に「現行の通常Alignで原報告を再現した」と判断しない。
summary-only表示やtransport欠落は残り、原因をstubするAlign unitだけでは実接続の受入にならない。

### known failure移行表

これは識別すべきfailure familyのinventory。新code registryや新Product判断ではない。
safe contextは有限field ID、ordinal、必要な数値・位置・codepoint等を型/長さ検査して扱う。
自由なpattern、file/record名、path、SVG、例外文をcontextへ移してはならない。

| 原因 / producer・実際のstage | 保持する修正情報・safe context | 既存ownerでの次action |
| --- | --- | --- |
| GenBank欠落 / GFF3+FASTA不足、run-analysis入力検査 | 必要な入力組合せ、sequence ordinal、固定field ID | 該当uploadを補いGenerate |
| recordsなし、selector範囲外/曖昧、record planning/discovery | 件数・ordinal、曖昧selectorには既存#index案内。selector原文は公開しない | recordを選び直し、必要なら#indexで区別 |
| Region/cropの片側欠落、非数値/非整数/1未満/start≥end、planning | 両端指定条件、座標条件、固定start/end field・sequence ordinal | Region欄を修正して再生成 |
| annotation target/comparison plan不成立、request準備 | 不成立種別と対象ordinal。record名・input本文は除外 | annotation/comparison対象を修正 |
| BLAST/比較FASTA不足・空sequence、pairwise input/LOSAT準備 | BLAST outfmt 6/7、FASTAが必要という既存案内とordinal | 比較入力を追加/交換して再試行 |
| Pairwise Match Height不正、input検査 | Autoまたは正の有限値 | Heightを修正 |
| Depthなし/series不整合/record index無効、depth validation | uploadかtrack無効化、series/record ordinal等の既存案内 | Depth input/対応recordを修正、またはtrackを無効化 |
| Depth min>max、window/step/tick/font≤0、run-analysis検査 | min≤max・正値条件、該当固定field | Depth設定を修正 |
| GC min/max非有限・逆転、tick/font≤0、run-analysis検査 | finite/min≤max・正値条件、該当固定field | GC設定を修正 |
| Circular slotが内側に収まらない、Python layout | cannot fit insideの有用な修正情報、slot ordinal/固定geometry field。custom slot名を流さない | track設定・配置を見直して再生成 |
| Python regex syntax、colors/labels/whitelist/visibility compile | re.errorまたは明示causeの有限reason、TSV row/固定field、Python character位置（不明は不明） | pattern/TSVを修正。既存Color field draftには後続Retry/Revert |
| TSV列/必須値/color/action不正、file/label/visibility parsers | 行・必要column・固定field/有限reason。syntaxとは区別 | TSVを修正して再import |
| source/view comparison identity不一致、pairwise_match.py 254–300 | 識別可能な整合性failureと実render stage。identity原文は除外 | 比較入力/操作を見直し、必要なら明示保存したSessionで別途再現調査。identityを自動補正しない |
| Worker asset/init/resource staging/protocol/helper失敗 | 実際に失敗した段階、finite kind。regex原因と断定しない | 同じ操作を再試行、必要ならSession保存後に再起動等の実在する回復 |
| Result list/catalog/sanitizer/admission不正、candidate admission | 固定failure種別とactual stage | 旧Resultが保護された事実を確認して再試行 |
| History finalize/rollback失敗、transaction owner | finalizeとrollbackの区別、復元成否。復元未確認を成功と表示しない | ownerが確認できる状態から復旧、必要ならSession保存 |
| preset fetch/import read失敗、rule/file preparation | import/preset準備stageとfinite reason。任意準備失敗をInvalid ruleと呼ばない | 再import/再試行、旧canonical維持 |
| SVGなし/空content/catalogなし、export capture | 現在Resultの有無、SVG export/interactive準備stage | 図生成または利用可能なexportを選択 |
| PNG dimension/DPI/parse/conversion/load/blob、export実stage | dimension/DPI修正案内、該当固定field、実conversion stage | DPI修正またはSVG export等の既存操作 |
| PDF library/font/glyph不対応、export/pdf-fonts | bounded codepoint、既存Use SVG to retain this text.の代替案内 | SVGへ切替。日本語/中国語font実装はBUG-07 owner |
| cleanup secondary、Worker/export cleanup | 主因を維持し、bounded secondary種別・stageのみ | 主因の回復を妨げない |
| unknown | 分かるoperation/stageとstable code、不明cause/位置を捏造しない | 実operationのretry/編集/明示保存へ接続 |

cancel/stale/supersededはこの表のerrorへ変換しない。
raw consoleの移行対象にはrule preset、Color/Label/whitelist import、run-analysisのwarmup/file/CDS関連catchもある。
同じfailureの旧prefix/serializer/Circular slot extractor（1329付近）は、その移行時に除去する。
新旧classifierを併存させず、必要なCircular修正情報をunknownへ潰さない。

## 3. Regex入口・evaluator・commit

| 入口 | evaluator / 現在のcommit境界 | S00所見 |
| --- | --- | --- |
| Color既存field/manual新規 | rule-actions→createRulePreparation→EVALUATE_RULES→web_support/rule_matching→native colors re.I→現在性確認→History/atomic live commit | JS syntax compileなし。新規の長さ/ReDoS heuristic confirmは構文compilerではない。既存fieldのfailure finallyがdraftを消す |
| Color specific TSV / preset | file-imports parse→watchersまたはpreset action→同じPython preparation→現在性確認→rules/legend commit | parse failureとruntime/syntaxを分ける。preset fetch失敗にraw console。default color/priorityはregexではない |
| Label TSV | 5-column parse→同じpreparation kind=label→native _build_label_override_rules re.I→intent/revision/Result/sourceSVG等を確認→override commit | 空catalog/SVGなしでもcompile評価。Label re.errorはParseErrorの明示causeになる |
| Label manual text | 既存label action→escaped exact literal selector/override→Generate native label owner | arbitrary regex入力欄として扱わない。Label arbitrary-regex preset入口は見つからない |
| whitelist manual / TSV | index.html 4047付近のmanualWhitelist、watcher parseWhitelistRules（3 columns）→canonicalへ追加→Generateでnative _build_whitelist_map re.I | import時はPython preparationなし。構文検査はGenerateへ遅延。Color live atomic validation/draftと同じ契約に広げない |
| visibility manual rule | visibility-actions 430付近→normalized manual rulesへ直接mutation→Generateのvisibility TSV staging→native compile_feature_visibility_rules re.I | manual mutation時の構文preparationなし。native compileはfeatures走査前に実行 |
| visibility TSV | feature-visibility parser→session-request 4032付近のadmitted typed resource復元→canonical→同じnative Generate | parserは構造検査のみ。独立したupload watcher/UIは見つからないため想定しない |
| visibility live hash fallback | direct override→selector cache→cache missかつwildcard record/type・hash qualifierに限りfeature-visibility.js 429付近 new RegExp(value) | **Python-domain patternのJS再compileあり。i flagなし。例外はdefaultへ黙って戻す**。Color/Label事前検査でもSearchでもない |
| Feature Search Regex | search-core.js 409付近 new RegExp(query,'i')→search targetsのみ | 仕様としてJS。Invalid regexは方言不明。Word/IUPAC検索と対象集合は維持 |
| standalone Regex | standalone-interactivity-assets.jsのcompileSearchMatcher（2676付近）→export内検索targets | 同じJS仕様。Pythonを埋め込まない。download bytesと実document起動は後続browser受入 |

[rule-matching owner](../../../gbdraw/web/js/app/rule-matching.js)はWeakMapの一時prepared cacheを使う。
pendingはnullでありvalid/non-matchと同じではない。revision/catalog/Result/mode/rules snapshotのguardと
一Worker・一評価・prepared reuse・優先順を保つ。
visibilityのJS fallbackは別の到達可能な境界として記録した。
Color field回復のreceiptだけでvisibility semantics変更まで許可されたと解釈せず、後続ownerが
supported behaviorとProduct preflightを確認する。S00で修正方針を選択していない。

### 実測matrix

[Python script](./evidence/s00-regex-probe.py) / [結果](./evidence/s00-regex-probe.json):
8 patterns × corpus/unrelated/empty 3 catalogs × 6 owner calls = 144 calls。
Color/Label helperとnative Color/Labelのacceptance・targets同値をassertした。
whitelist/visibilityも実native compilerとmatcherを呼ぶが、Worker/Generate/UIを起動しない。

corpus ordinal: 0 NADH、1 nadh、2 β-lactamase、3 ı、4 i、5 İ、6 unrelated。

| pattern | Pythonの全6 owner / corpus targets | JS Search構文 / standalone corpus targets |
| --- | --- | --- |
| (?i)NADH | 受理、0/1 | 拒否 / なし |
| (?P<enzyme>NADH) | 受理、0/1 | 拒否 / なし |
| [ | 拒否、position 0 | 拒否 / なし |
| (?<enzyme>NADH) | 拒否、position 1 | 受理 / 0/1 |
| NADH\Z | 受理、0/1 | 受理 / なし（同じ意味ではない） |
| \bβ | 受理、2 | 受理 / なし（word boundary差） |
| i | 受理、3/4/5 | 受理 / 4（Unicode casefold差） |
| β | 受理、2 | 受理 / 2 |

unrelated/emptyでも有効式は受理してtargetsなし、不正二式は全6 ownerで拒否された。
Python 3.13ではre.errorの実type名はPatternError。Label/whitelist/visibilityはParseErrorのcauseに保持。
synthetic patternのnative logger出力も捕捉した。将来consoleへ自動公開しない検査が必要。
[JS script](./evidence/s00-js-probe.mjs) / [結果](./evidence/s00-js-probe.json)は
実Searchをempty catalogで呼び、standaloneの正確なregex関数branchを固定anchorで抽出して実行する。
standalone全文のdocument初期化やdownload integrationの代用ではない。
visibility hash fallbackはPython専用式と不正式をdefaultへ戻し、JS named groupをoffへ適用し、
NADH hashに対する^nadh$をonとした。native Pythonとの意味の差を確認した。
非syntax準備rejectを注入した実field actionはInvalid rule: Synthetic runtime unavailableを表示し、
inputをNADHへ戻した。History/SVG/preparationはstubであり、実Worker失敗の証拠とは区別する。

### 原監査の限界

Issue本文とcommentsを取得した（updatedAt 2026-09-26T04:24:45Z、commentsなし）。
BUG-19はJS preliminary validationという報告だが、実pattern・入力画面・公開buildは記載されていない。
参照されたGUI_AUDIT_DEV_20260926.mdは開始/dev/関連branch treeおよび取得済み全refの当該path履歴にない。
原監査artifact・buildは取得できず、未確認のまま。
現在のColor/Label Python成功や発見したvisibility差分を原報告と同一の原因と断定しない。

## 4. 他作業との重複とwriter

related refs/resultsの完全なfile SHA・lastCommitはinventory JSON参照。
branch対devの**三点diff**とancestorを調べ、古いauthorityを持つbranchの二点diffを
新しいdev契約の退役意図と取り違えない。

| branch / 監査head | 実装・dev統合状態 / 引継ぎ |
| --- | --- |
| fix/issue-601-export-output-20260926 / 43836eb77798924bc826d0064d3cc5e219379d91 | 未統合は計画資料のみ、runtime差分なし、session resultなし。MASTER_PLAN/APPROVED_DECISIONS/DECISION_PACK_02/SESSION_03_ERRORS_REGEXがBUG-15/19を依然担当。**owner移管は確認できない** |
| fix/issue-598-alignment-direction-reset-20260926 / d5edc4bb00b0d793d94361202c31a270f78a4aa0 | evidence S00/S01を確認。S01 domain projection実装済みだがruntimeはdev未統合。base 302dfa、通常merge fe6ec7f1。Worker/client/Alignは共有owner。S02/S03以降を完了扱いにしない |
| fix/issue-602-linear-live-edit-20260926 / 94faf8a98eddf823f5ebfed859d324e258daa689 | results S00–S05を確認。S01 disclosure、S02 default、S03 status data、S04 status UI、S05 compact Editorがbranchにある。base af5d942aをmerge edd3ce05で取り込み済み。runtimeはdev未統合。S06 compact review/S07は未完了 |

#598のruntime差分はapi/record_planning.py、diagrams/linear/assemble.py、similarity-alignment.js、
diagram-generation.js、diagram-generation-worker.js、web_support/similarity_alignment.pyの6 files。
resource stagingとdomain projectionの契約をS02 transport変更やS03 Align変更で上書きしない。
S00で記録された古いPD-OI-039のMatch競合は最新devのrevision 2（PR #609）を基準に再評価する。

#602のruntime差分はindex.html、app-setup.js、generation-status.js、svg-styles.js、config.js、
session-active-config-contract.js、session-request.jsの7 files。
derived statusのpending/live failure、compact Editor/canvasの契約をS03/S04で保護する。
S05は独立したEditor phase。metadata-free Sessionの#564由来staging失敗と、
Legend edit→retryのSanitized SVG content is missing a Legend binding.は別の既知問題として記載されている。
当該branchの受入をS00や全dev staging合格へ流用しない。

dev統合済みのauthorityは#602 PR #604、#598 PR #608と#609。両branchのruntimeとは別。
通常Alignの原因返却5bcdd4f5とPython rule owner化#538も既存devにある。

専用cloneでこのsessionは一writer。process照合では別のissue601 writerを特定できなかったが、
remote writerのlease・担当者の作業継続/停止・分散排他は確認できない。
他sessionのtree/index/envは調査・変更していない。process不在はowner移管の証拠ではない。
runtime開始前にexport側BUG-15/19担当の移管と、#598/#602共有fileの作業順・基点SHAを明示する。

関連resultのimmutable参照:
[598 S01](https://github.com/satoshikawato/gbdraw/blob/d5edc4bb00b0d793d94361202c31a270f78a4aa0/docs/internal/issue-598-alignment-direction-reset-20260926/evidence/S01.md)、
[602 S05](https://github.com/satoshikawato/gbdraw/blob/94faf8a98eddf823f5ebfed859d324e258daa689/docs/internal/issue-602-proposal-20260926/results/S05.md)、
[export計画](https://github.com/satoshikawato/gbdraw/blob/43836eb77798924bc826d0064d3cc5e219379d91/docs/internal/issue-601-export-output-20260926/MASTER_PLAN.md)。

## 5. Product authorityと二receipt

最新dev OIPC contract revisionは21、SHA-256は
62e3a9c08ceb64acc349ace81a97a9c87dc45db4d18187a3b78e6e131b2dc595。
正式な[OIPC](../OPTION_INTEGRITY_PRODUCT_CONTRACT.md)、[Product map](../../../tools/web-product-impact-map.json)、
[BD store](../../../tools/web-product-decisions.json)を照合した。
BD storeはschema 1、maintainerLoginsにsatoshikawato、decisionsは空。
mapにはcanonical render requestとsaved Session continuityの二concern、および4 hard complete architecture coverage。
今回の二concernは未登録。存在しないBD/PD番号は予約していない。

| 既存authority | 今回も独立に維持する条件 |
| --- | --- |
| OIPC-C07、PD-OI-016 rev1 | failed/canceled/stale/supersededで成功Result/requestを置換しない。draft/修正案内/最後の成功identity/corrected retry |
| PD-OI-027 rev5、029 rev3、031 rev5、034 rev5 | record-owned directionとReset、local no-Worker review、最終一batch validation、actionable cause、失敗retry、atomic Apply/旧artifact保持 |
| PD-OI-035 rev3、039 rev2 | canvas identity/keyboard/focus/overlay、exclusive directions。039の現選択はA / EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH。旧Match維持を現条件としない |
| PD-OI-037 rev1、038 rev1 | application statusをderivedに保ち、compact Editorで閉じるowner/選択tabを保つ |
| PD-OI-044/045各rev1 | truthful discoveryとexclusive Session operation。error表示がsource replacement/continuationを変更しない |

390×844/740でcanvas可視高200px以上等のcompact条件は独立要求。
Details/draftの配置だけで既存canvas/Editorの意味を再選択しない。

[診断receipt](./decisions/DECISION_01_ERROR_DISCLOSURE.md)と
[field回復receipt](./decisions/DECISION_02_REGEX_EDIT_RECOVERY.md)は各九項目の本文とJSONが完全一致。
scenarioRevision 1、owner/date/理由/維持/退役/残余リスクも保持され、各digestは次のとおり。

- web.errors.diagnostic-disclosure / A / GUIDANCE_WITH_BOUNDED_DIAGNOSTICS:
  8617888cb1838521f78640db4653e2dcf88dd427cd466bbbc33d59efacaceae8
- web.rules.rejected-pattern-edit-recovery / A / KEEP_REJECTED_PATTERN_DRAFT:
  9cc66bc30f97cc995b9b0002b094ba74cfa00a8e77b163c33074689f0e73b496

export計画の別candidateはweb.errors.user-facing-diagnostic-disclosure / rev1 /
A / ACTIONABLE_SUMMARY_WITH_SAFE_DETAILS。devのactive contractには存在しない。
既知cause・actual stage・optional collapsed bounded Details・入力privacy・cancel・draft・Result/History・
keyboard/390pxは同じ意味を持つ。今回receiptはさらにoperation、手動Copyの安全な可視情報限定、
clipboard不可の選択コピー、自由ログ等の明示除外、初回/旧Result区別を独立に要求する。
別candidateのoption IDが異なることだけを製品矛盾としないが、当該receiptは九項目の完全一致ではなく、
Copy等の独立要求を満たす正式authorityの代用にならない。
既存active判断と今回選択との実質的な製品矛盾は監査時点では見つからなかった。
同義の二active concernを二重登録する収束は止め、S01で担当間のcandidate調整と全要求coverageを確認する。
未承認のkey変更、既存receipt退役、rationale/残余リスク追加を推測しない。
KEEP_REJECTED_PATTERN_DRAFTは既存Color patternだけ。TSV/new/preset/Label/whitelist/visibilityへ拡張しない。

## 6. S01–S05の直列owner/path引継ぎ

| session | 所有path / 次ownerへ渡す契約と検査 |
| --- | --- |
| S01 authority-only | 専用product branchの既存OIPCだけ。二つの九項目receiptとsame-semantic candidateのcoverageを確認。PD番号は最新baseで決定。runtimeを含めず、正式dev統合SHA/ancestorをS02へ渡す |
| S02 cause/transport | exceptions.py、pairwise_match.py、限定web_support adapter、python-helpers.js、Worker/client、error-normalization。#598共有変更の受入基点を確定してから一writerで更新。ValueError捕捉互換、typed re.error/cause、actual render/helper往復、secondary/unknown/bounds/sentinel/double normalizationを検査し、finite modelをS03へ渡す |
| S03 callers/UI | run-analysis、Align、rule/import/export caller、app-setup/index、Search/standalone。#598 Align・#602 composition/status/compactを保持。全error return、actual generation→caller→UI、retry/rollback/cancel、Summary/Details/Copy/consoleのsentinel、clipboard不可、390px/keyboardを確認。field error modelをS04へ渡す |
| S04 Color field draft | rule-actionsと既存composition/index、必要なfeature-editor内focused helper。S03のmodelを再利用し、同じfileを同時編集しない。row/revision/documentに一owner、失敗text保持、Retry/Revert、History/Session/reset/row削除/一時mode/drawer/stale全transitionを確認 |
| S05 integrated acceptance | tests・既存reference docsと原因ownerの限定fix。S00–S04契約/receipt条件をAND-of-ORで覆い、実render/helper→transport→UI、actual Align、Python parity・JS download bytes、Save/fresh Load/Generate/Export、25,000-feature reuse/Worker数、科学出力/native互換を確認 |

各受渡しで最新remote、owner、受入HEAD、共有fileの差分、superseded pathの除去、必要なregressionを記録する。
意味ownerはproducer / transport / normalizer / operation owner / draft ownerの既存境界へ収束させる。
汎用error dispatcher、第二regex evaluator、汎用draft framework、新永続schemaを前提にしない。
本S00のproduction owner/path差分は空で、OE/PE/CBは変化しない。

既存tests/web/similarity-alignment-actions.test.mjsはrunAnalysis/Historyをstubするため、
callerのcause保持は測れても実producerで原因が落ちないことは測れない。
S00のprobeも関数/native境界まで。S02はreal helper/render→Worker/client、S03/S05は実orchestrationと
browser操作→UI/retryの証拠が別途必要。
既存Python rule testsとChromium二シナリオはdialect/canonical保持のbaselineであり、
未実装draft/Copy・将来のproduction変更合格には流用しない。

## 7. 実行command・証拠の再利用・検証

以下は専用clone rootで実行。logs/temp/cacheは当該sessionの/tmp配下へ隔離。

~~~bash
git fetch origin
git switch --track -c fix/issue-601-bug15-bug19 origin/fix/issue-601-bug15-bug19
git pull --ff-only
git status --short --branch
git rev-parse HEAD origin/dev
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
gh issue view 601 --repo satoshikawato/gbdraw --json number,title,body,comments,updatedAt,url
PYTHONDONTWRITEBYTECODE=1 python docs/internal/issue-601-bug15-bug19-implementation-20260926/evidence/s00-inventory.py
PYTHONDONTWRITEBYTECODE=1 python docs/internal/issue-601-bug15-bug19-implementation-20260926/evidence/s00-regex-probe.py
node docs/internal/issue-601-bug15-bug19-implementation-20260926/evidence/s00-js-probe.mjs
node docs/internal/issue-601-bug15-bug19-implementation-20260926/evidence/probe.mjs
~~~

inventoryのtree/receipt/9 returns assertions、144 Python calls、JS probeのassertionsは成功。
既存probe出力は保存済みprobe.jsonとbyte一致。
Issue #598/#602、関連branchの計画とresults、PR #604/#608/#609、dev source logもread-onlyで確認した。
参照先・receipt/matrix/authority整合性・変更scope・whitespaceはコミット前と統合後に検査する。

| 先行証拠 | 再利用可能な範囲 / 根拠 |
| --- | --- |
| 41 Node / 9 Python | 保存済みbaseline。source/tests/input tree同一、Node/Python version同一。正常rule評価とstubbed Alignなど実際に測った範囲だけ。新受入にはしない |
| 2 Chromium | 保存済みbrowser baselineとして参照。source/inputは同一だが元browser binary・bundleの同一性は未確認。今回再実行・現公開build成功とは表現しない |
| initial synthetic probe | 同source/inputで再実行しJSON byte一致。transportの情報欠落等の限定再現に使う |

docs/配下の全変更はtrusted-base CI policyでdocumentationに分類され、candidate admissionの必要jobは
recipes-standard。[SELECTIVE_CI](../SELECTIVE_CI.md)のDOCUMENTATION_ONLY_PRはbase stagingを継承しない。
別にtrusted-base Web change budgetも実行する。candidate自身のcheckerをauthorityとして使わない。

~~~bash
git archive af5d942af60353dda199aa487da9152a3576b3fe tools | tar -x -C /tmp/gbdraw-issue601-s00-vIoB52/trusted-tools
node /tmp/gbdraw-issue601-s00-vIoB52/trusted-tools/tools/check-web-change-budget.mjs --base af5d942af60353dda199aa487da9152a3576b3fe
PYTHONDONTWRITEBYTECODE=1 python -m pytest tests/ -m 'recipe and not slow' --durations=30 --basetemp /tmp/gbdraw-issue601-s00-vIoB52/pytest -o cache_dir=/tmp/gbdraw-issue601-s00-vIoB52/pytest-cache
~~~

preflight Web gateはGate PASS / Review CLEAR、4 hard rules CONFORMING。
recipe検査は189 passed、6090 deselected、149.53s（Python 3.13）。
83 Markdown参照、27 source/6 authority hash、144 Python owner結果とJS field evidenceの整合性検査は合格。
recipe実行が自分のclone内のLOSAT executable bitを変更したため、終了後に元の100644へ戻した。
production/test/generatedの最終差分は空とし、資料diffと別に確認する。
最終full-SHA Web gateと統合後scope/whitespaceはコミット後に確認してhandoffで報告する。
Python 3.13でのlocal recipeはCIのPython 3.11を実行した証拠ではない。

remote baselineの[Tests run 36244964256](https://github.com/satoshikawato/gbdraw/actions/runs/36244964256)は
af5d942aに対してfailure。CI impact planがDOCUMENTATION_BASE_EVIDENCE_UNAVAILABLE /
RUN_NOT_SUCCESSFULで停止し、Dev staging / gateもfailure、leaf jobsはskipだった。
このdirect-parentのstaging成功は再利用できない。docs-only candidate admissionの免除と
dev integrated stagingは別境界であり、S00のlocal合格やdev pushをstaging成功と報告しない。
CI guard修正・手動full workflow・release/deployをこのsessionで追加実行していない。

## 8. S01開始条件と未確認事項

| 条件 | 状態 / 次sessionが確認する境界 |
| --- | --- |
| S00 result・計画のwork branch/dev公開 | commit後に両remote一致・cleanをhandoffで確認する。dev push未完了ならS00未完了 |
| 二つの完全なA receipt | 充足。九項目・digest一致。正式authority/実装完了とは別 |
| 最新dev authority inventoryと製品条件 | 監査時点で充足。S01で最新baseと新候補を再確認 |
| same-semantic disclosure candidateとの調整 | 未充足。export candidateが残る。二重active化・key変更・退役を推測しない |
| BUG-15/19 owner移管・remote writerの直列化 | 未確認。export計画に残存。S01 authority writerの調整、特にS02以降のruntime開始前に明示する |
| authority候補のdev統合 | 未充足。S01未開始。このS00 push許可はauthority/runtimeのdev統合許可へ拡張しない |
| #598/#602の共有runtime引継ぎ | 未充足。監査headではdev未統合。最新result/基点SHAと共有file順序をruntime前に確定 |
| 原監査pattern・画面・公開build | 未確認。BUG-19 closure証拠ではない |

S00公開後、S01のreceipt照合とauthority-only準備の前提は揃う。
ただしsame-semantic candidate/writerの調整なしに無条件のS01収束完了・push可能と宣言しない。
二つの正式authorityがdevに統合され、ancestor確認とowner移管/共有file引継ぎが満たされるまで
S02以降の重複runtime編集は開始不可。未充足を解消するためにS00の変更scopeを広げない。
