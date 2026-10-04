# Linear入力・表示UI改善 — 実装セッション用INSTRUCTION PROMPTS

このファイルは、gbdrawのIssues #559、#560、#562を扱う各セッションの開始指示である。
各promptの読者は過去の会話を知らない。リポジトリrootをcurrent working directoryとして、
指定したcode block全体を一つのセッション指示として使用する。

問題、監査結果、設計、Product判断の状態、受入条件、実施記録の管理先は
[総合計画書](LINEAR_UI_559_560_562_MASTER_PLAN_2026-09-21.md)である。

対象Issue:

- [#559 Depth TSV panel](https://github.com/satoshikawato/gbdraw/issues/559)
- [#560 Add/Remove visibility and removal choices](https://github.com/satoshikawato/gbdraw/issues/560)
- [#562 Linear defaults, Record Labels, Legend, Auto visibility](https://github.com/satoshikawato/gbdraw/issues/562)

## 使用順序

| Prompt | いつ使うか | 依存 |
| --- | --- | --- |
| S0 | #560および#562三concernsのruntime着手前 | Product Decisionとauthorityだけを扱う |
| S1 | #559を実装するとき | PD-OI-025を含むbase。S0を待たなくてよい |
| S2 | #560を実装するとき | 選択済み#560 authorityがruntime baseへmerge済み |
| S3 | #562のCurve defaultを実装するとき | 選択済みmatch-style authorityがruntime baseへmerge済み |
| S4 | #562のsidebar再編を実装するとき | 選択済みorganization authorityがruntime baseへmerge済み |
| S5 | #562のAutoとSession互換を実装するとき | 選択済みAuto authorityがruntime baseへmerge済み |
| S6 | 選択された全runtime変更を最終受入するとき | 該当するS1–S5完了 |

S0とS1は独立に実行できる。S2–S5は対応する候補authorityとruntimeを同じ未merge変更に含めない。
S3、S4、S5は互いに独立して採否・実行でき、一つのdecisionを別sessionのauthorityに使わない。
各セッションの終了時に総合計画書第14節へ記録し、後続担当者は会話ではなく、その記録、
Git history、worktree、test evidenceから状態を復元する。

## 固定実装ブランチ

S1〜S6のruntime、tests、docs、検証記録は、計画策定sessionで最新`origin/dev`から作成した
`fix/linear-ui-559-560-562-20260921`へ積み上げる。各sessionは最初にこのbranch名、HEAD、
upstream、worktree、base ancestryを確認する。別のruntime branchを新規作成せず、現在branchが
異なる場合は未commit差分と無関係な変更を確認してからこのbranchへ戻る。各sessionを一commit
候補として扱い、main/devへ直接commitしない。

S0のauthority-only変更だけは例外である。必要なauthority branchを最新`origin/dev`から別に作り、
runtimeを含めない。authorityが`origin/dev`へmergeされた後、固定実装ブランチへ戻り、fetchした
merge済み`origin/dev`へrebaseしてからS2〜S5の該当runtimeを開始する。候補authority commitを
固定実装ブランチへ直接cherry-pickし、同じcandidateのruntimeを自己承認してはならない。

## S0 — Product preflightとauthority-only変更

```text
gbdraw Issue #560と、Issue #562に含まれる三つの独立concernsについてProduct Impact
preflightを行い、必要なhuman
Product Decisionを既存のdurable authority経路へ反映する準備をしてください。
runtime codeは変更しません。過去の会話は前提にしません。

対象:
- #560: Linear Fileカードで、ファイル内容のClearとカードDeleteを分け、global Removeの
  誤操作を防ぎ、Add/Removeを見つけやすくする。
- #562 match-style: fresh/reset Linearのmatch styleを変更する。
- #562 organization: Plot Title・Record Labels・Legendのsidebar hierarchyを変更する。
- #562 Auto: Accession/Lengthを独立したAuto/Show/Hideへ変更する。
- #559: 既存PD-OI-025で実装可能かだけ確認する。新しい判断を作らない。

最初に全て読む:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/LINEAR_UI_559_560_562_MASTER_PLAN_2026-09-21.md
5. docs/internal/PRODUCT_IMPACT_RATCHET.md
6. docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md
7. docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md
8. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
9. docs/internal/WEB_CHANGE_POLICY.md
10. docs/SESSION_COMPATIBILITY.md

開始確認:
- 固定実装ブランチfix/linear-ui-559-560-562-20260921で計画と現在の実施記録を読む。
  git status、HEAD、upstream、origin/devとの差を確認し、無関係な変更を保持する。
- authority変更が必要な場合だけ、git fetch origin後、最新origin/devから追跡先なしの別の
  authority branchを作る。main/devへ直接commitせず、runtimeをauthority branchへ含めない。
- GitHub Issues #559/#560/#562の現在の本文とcommentsを読み、更新日時を記録する。
- 最新baseのProduct Impact map、accepted PD、BD recordを検索する。同じconcern、同じ
  scenarioを表すauthorityがある場合は並行したauthorityを作らない。

監査済みの出発点:
- 計画策定時のbaseはf4760476915194fc0586642e0b03ec561b96ab01だった。
- #559はPD-OI-025 FILE-BULK-WITH-RECORD-OVERRIDESとOIC-020の範囲で、
  File-level Depth disclosureを初期openのままclosableにし、既存Add series actionへの
  direct affordanceを置く提案である。
- #560と#562には、計画策定時点でbaseへ反映済みのdurable decisionがなかった。
この情報を最新base確認の代わりにしない。

作業A — developer preflight:
1. 各提案で到達可能なuser-visible differencesとcontinuationsを列挙する。
   implementation file名ではなく、Clear後に何が残るか、legacy Sessionが何を表示するか、
   Autoがどのrow条件でHideするかなどの結果を比較する。
2. #559は総合計画書第5.1節の結果がPD-OI-025/OIC-020を維持することを確認する。
   初期closed、Record options内への移動など到達性を弱める変更は含めない。
3. #560はconcern候補linear.input-source-removalについて、総合計画書第5.2節の
   D560-A/B/CをDecision Packにする。推奨はD560-Aだが、開発者が選択済みにしない。
   GenBankとGFF3+FASTAのclear単位、保持state、sole-card、global Remove、Undo/Resultを
   optionごとに明記する。
4. #562は一つのDecision Packへ束ねず、次の三concernsを別々に作る。
   a. linear.comparison-match-style-default:
      CURVE-FRESH-RESET（推奨） / KEEP-RIBBON。
      fresh/reset、legacy missing Ribbon、explicit saved value、Circular/CLI/API非変更を扱う。
   b. linear.presentation-control-organization:
      REGROUP-TITLES-LABELS-LEGEND（推奨） / KEEP-CURRENT-ORGANIZATION。
      state semanticsを変えないsidebar情報設計だけを扱う。
   c. linear.record-label-auto-visibility:
      INDEPENDENT-AUTO-SHOW-HIDE（推奨） / KEEP-BOOLEAN-VISIBILITY。
      effective row rule、selected modes、Session 43、legacy boolean migrationを扱う。
   三つをall-or-nothing optionへ再結合せず、開発者が選択しない。
5. 各concern/choiceについてIMPLEMENT_EXISTING_AUTHORITY、EVIDENCE_REQUIRED、
   PRODUCT_DECISION_REQUIRED、NOT_ALLOWEDのrouteと理由を示す。
6. Product Decision OwnerがコピーできるPRODUCT_DECISION response templateを出す。
   rationale、must preserve、may retire、accepted residual risk、owner、dateを空のまま示し、
   開発者が補完しない。

作業B — 明示回答が既に存在する場合だけ:
1. human responseのconcern、revision、choice、rationale、must preserve、may retire、risk、
   owner、dateが全て明示されていることを確認する。不足項目を推測しない。
2. 選択された結果だけを、最新baseで使われている既存のdurable authorityへ正確に転記する。
   新しいdecision store、runtime flag、parallel JSON schemaを作らない。
3. authority testsとdocumentation contractsを実行し、生成したmachine representationを
   human responseと突き合わせる。
4. authority-only diffへruntime、UI、Session version変更を含めない。
5. merge前candidateはdependent runtimeを許可しない。外部push、PR、mergeは、このtaskで
   明示的に許可されている場合だけ行う。

human responseがない場合:
- Decision Pack、選択肢、response template、独立調査を完成させる。
- authority未反映の#560または#562 concernに依存するruntimeを開始しない。
- 判断待ちを「実装困難」や「全作業blocked」と表現せず、S1 #559は独立して開始可能と記録する。

設計原則:
- SOLID: Product outcomeとimplementation ownerを混同しない。
- KISS: #560一concernと#562三concernsに限定し、一般的なUX policyを作らない。
- DRY: 既存Option Integrity contractとProduct Ratchetを唯一のauthority経路として使う。
- YAGNI: 決定UI、bot、汎用schema、将来optionを追加しない。

終了条件:
- 総合計画書第14節へbase SHA、branch、Issue更新日時、classification、evidence、
  Decision Pack、human responseの有無、authority diff/test、各runtime sessionの開始可否を記録する。
- 対応するauthorityがmergeされていない場合はS2〜S5の該当sessionを開始可能と報告しない。
- runtime修正完了とは報告しない。
- repository guidanceに従う英語のproposed commit titleと短いsummaryを示す。
```

## S1 — Issue #559: File-level Depth disclosureとAdd series

```text
gbdraw Web版LinearモードのIssue #559を、既存PD-OI-025とOIC-020を維持して実装し、
focused browser verificationまで完了してください。過去の会話は前提にしません。

問題:
各Fileカードの"Depth TSV files (applied to all records)"が常時展開され、series数や
File数が増えるとsidebarを大きく占有する。File-level uploaderの直下には新しいDepth
seriesを追加する明確なbuttonがない。

最初に全て読む:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/LINEAR_UI_559_560_562_MASTER_PLAN_2026-09-21.md
5. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
6. docs/internal/PRODUCT_IMPACT_RATCHET.md
7. docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.mdのPD-OI-025とOIC-020
8. docs/REFERENCE/web-app.md

開始確認:
- branchがfix/linear-ui-559-560-562-20260921であること、HEAD、upstream、base ancestryを確認し、
  無関係な変更を保持する。別のruntime branchを作らず、main/devへ直接commitしない。
- この固定branchのbaseがPD-OI-025を含むことを確認する。
- 最新sourceでindex.html、app-setup.js、linear-sources.js、depth-track-state.js、
  session-request.jsを読み、既に直った部分を重複実装しない。
- 総合計画書第14節にS0または他sessionの記録があれば、同じworktreeの差分を識別する。
  #559は#560/#562の未決Product Decisionを待たない。

実装する結果:
1. 各Linear FileカードのFile-level Depth領域をnative details/summary disclosureにする。
2. 初期状態はopenとし、利用者が閉じてsidebar高さを減らせるようにする。
   disclosure open stateはform、adv、Session、Historyへ保存しない。
3. summaryにDepth TSVであること、対象record数、logical series数、Mixedの有無を簡潔に示す。
   matrixから導出し、summary専用のmutable stateを作らない。
4. record listとRecord optionsを閉じたまま、File-level upload/clearとAdd seriesへ到達できる
   hierarchyを保つ。
5. series cardsの直下に"Add Depth TSV series"相当の明確なfull-width actionを置き、既存の
   addLinearDepthTrack()を呼ぶ。一回の操作でlogical columnを一つだけ追加する。
6. 別のDepth settings領域に同じglobal series追加buttonがあり、二つの概念的入口を残す必要が
   なければFile-level文脈へ集約する。action実装は複製しない。
7. Common/Mixed/Empty表示、mixed warning、File-level apply/clear、per-record sparse override、
   single-record Fileの重複uploader抑止を維持する。

変更してはいけないもの:
- record-major Depth matrix、logical indexes、null cells
- Depth Session field、canonical request semantics、Worker/Python renderer
- track axis、normalization、color、heightなどの描画設定
- Circular Depth UI
- current Result lifecycle

test-first evidence:
- multi-record Fileをrecord list closedのまま操作するcaseを既存test ownerへ追加する。
- initial open、toggle close/open、summary、Add後のmatrix width、既存cell非移動をassertする。
- common → per-record Mixed → File-level clear → one Undoで完全復元するOIC-020 caseを維持する。
- same-name independent Files、later columns、null cells、Save/fresh Load、regenerationを確認する。
- 390x844、pointer、Enter、Space、unique accessible namesを実ブラウザで確認する。

主な対象:
- gbdraw/web/index.html
- gbdraw/web/js/app/app-setup.js（既存action/derived summaryの最小変更だけ）
- tests/web/depth-track-state.test.mjs
- tests/web/depth-track-session.playwright.spec.js
- tests/web/linear-multi-record.playwright.spec.jsまたは既存の最も近いjourney
- docs/REFERENCE/web-app.md（利用者向け挙動が変わる箇所だけ）

focused commands:
node --test tests/web/depth-track-state.test.mjs tests/web/linear-sources.test.mjs
node --test tests/web/session-request.test.mjs
npx playwright test tests/web/depth-track-session.playwright.spec.js --workers=1 --retries=0
node --test tests/web/architecture-contracts.test.mjs
git diff --check

browser wheelが必要なら現在sourceからpython tools/prepare_browser_wheel.pyで作る。wheelは
commitしない。Node @playwright/testが使えなければPython Playwrightで同じDOM journeyを
行う。Chromium sandbox errorは必要な権限で同じcheckを再実行する。

設計review:
- SOLID: UI disclosureはDepth matrix semanticsを所有しない。
- KISS: native disclosureと既存add actionだけである。
- DRY: add/apply/clearのproduction actionを複製していない。
- YAGNI: disclosure persistence、新Depth field、新module、新render pathを追加していない。

終了条件:
- 総合計画書DPT-01〜DPT-08を満たす。
- production、tests、docs、generated diffを別々にreviewする。
- 第14節へbase/branch、変更file、受入ID、正確なcommands/results、browser wheel identity、
  未実施項目を記録する。
- 英語のproposed commit titleと短いsummaryを示す。
- push、PR、merge、tag、deployは明示的な許可がない限り行わない。
```

## S2 — Issue #560: Add/Removeとatomic removal workflow

```text
gbdraw Web版Linear入力のIssue #560を、runtime baseへmerge済みのProduct Decisionどおり
実装し、Clear/Delete/Cancel、History、comparison、Sessionを検証してください。
過去の会話は前提にしません。

問題:
- sidebar下端のglobal Add/Removeが小さく見つけにくい。
- global Removeが非空last Fileを確認なしで削除する。
- FileUploaderの赤いRemoveはファイルclearに見えるが、Linear primary sourceではFileカード
  全体を消し得る。利用者は同じ位置にblank slotを残すClearとcard Deleteを選べない。

最初に全て読む:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/LINEAR_UI_559_560_562_MASTER_PLAN_2026-09-21.md
5. docs/internal/PRODUCT_IMPACT_RATCHET.md
6. docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.mdの、このissue用にmerge済みのauthority
7. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
8. docs/SESSION_COMPATIBILITY.mdとdocs/REFERENCE/web-app.md
9. 総合計画書第14節のS0/S1記録と現在のgit diff

厳格な開始条件:
- 最新baseのauthorityからconcern、scenario revision、selected outcome、must preserve、
  may retire、accepted residual riskを引用できること。
- 固定実装ブランチfix/linear-ui-559-560-562-20260921が、そのauthority merge commitを
  ancestorに持つこと。持たない場合はruntimeへ進まず、authority merge後のorigin/devへ
  同じ固定branchをrebaseする。
- authority candidateとruntimeを同じ未merge branchに置かないこと。
- authorityがない、または総合計画書の推奨D560-Aと別outcomeが選ばれている場合は、
  D560-Aを黙って実装しない。選択結果に合わせて受入表と小さな設計差分を先に記録する。
- branch名、HEAD、upstream、base ancestryを確認し、無関係な変更を保持する。別のruntime
  branchを作らない。

D560-A Pristine placeholderが選択された場合のnormative result:
- Clear file onlyはsource groupを同じFile順序の一つのpristine blank slotへ置換する。
- GenBank、またはGFF3+FASTAをprimary source単位でclearする。
- 旧record UID、selector/crop/reverse、label/subtitle、Depth binding、comparison endpoint、
  source-derived cacheを新slotへ引き継がない。
- Delete cardはsource group全体を削除する。sole FileではDeleteを無効にし、0 slotsにしない。
- global Removeはpristine blank last slotなら二枚以上で即時Delete、非空last Fileなら対象を
  示すconfirmation後にDeleteする。
- Cancel/Escapeは完全なno-op。draft変更だけではcurrent Resultを置換しない。

実装境界:
1. components.js::FileUploader
   - 既定のclear behaviorとHistoryを全既存利用者で維持する。
   - Linear primary sourceだけがopt-inでclear intentを外へ通知できる最小prop/eventを追加する。
   - source grouping、Delete、dialog、comparison規則をcomponentへ入れない。
2. linear-sources.js
   - groupLinearSourceRecords()をsource identityの唯一ownerとして使う。
   - 必要ならclear/delete対象を決めるpure planを追加する。DOM、Vue、History、cacheに依存させない。
   - createLinearSeq()がstate ownerにある場合、pure planへfactoryやVue stateを持ち込まず、
     coordinatorがplan結果からblankを作る。
3. app-setup.js
   - source UID/group identityでpending intentを持ち、index driftに耐える。
   - Clear/Delete/global Removeを同じempty判定と一つのcoordinatorへ収束させる。
   - 一操作につきhistory.runUndoable()を一回だけ呼び、applyLinearSeqMutation()を一回通す。
   - layout、comparison plan、derived artifact、raw cache、pending discoveryの既存reconciliationを
     個別に再実装しない。
4. index.html
   - Add/Removeを視認しやすいsize/hierarchyにする。
   - 非空cardの赤いRemoveからaccessible choice dialogを開く。
   - File名、record count、Clear/Deleteの差を表示し、Cancel、Escape、focus returnを持つ。
   - window.confirmで三択を模倣せず、新しい汎用modal frameworkも追加しない。

atomic state rules:
- dialogを開いただけではHistoryを増やさない。
- 対象sourceがdialog中に消えた場合は安全にcancelする。
- nested uploader Historyとcoordinator Historyを作らない。
- successful Clear/Deleteはrecords、rows、comparisons、resource bindingsを一つのUndoで戻す。
- explicit endpointはstale UIDを残さず、Adjacent/derived artifacts/cacheは既存contractに従う。
- selected/current Resultは次のGenerate成功まで保持する。
- dialog state、source order、removal intentをSessionへ保存しない。

tests:
- pure: non-empty GenBank、multi-record group、GFF3+FASTA、blank、sole source、missing target、
  no mutation、same-name independent upload。
- component: opt-in Linear sourceはclear-request、通常FileUploaderは従来update emission。
- History: Clear/Delete各一件、Cancelゼロ、Undo/Redo完全復元、nested operationなし。
- browser: first/middle/last File、global Remove、card Remove、Cancel/Escape、re-upload同じordinal、
  pending multi-record discovery、390x844、pointer/Enter/Space/focus。
- lifecycle: request、Save/fresh Load、current Result、failed/canceled/stale Generate。

主な対象test:
- tests/web/linear-sources.test.mjs
- tests/web/history-inputs.test.mjs
- tests/web/history-inputs.playwright.spec.js
- tests/web/linear-multi-record.playwright.spec.js
- tests/web/session-request.test.mjs
- 必要な既存comparison/Result contract tests

focused commands:
node --test tests/web/linear-sources.test.mjs tests/web/history-inputs.test.mjs
node --test tests/web/session-request.test.mjs tests/web/linear-record-layout.test.mjs
npx playwright test tests/web/history-inputs.playwright.spec.js tests/web/linear-multi-record.playwright.spec.js --workers=1 --retries=0
npm run test:web:comparison-contracts
node --test tests/web/architecture-contracts.test.mjs
git diff --check

設計review:
- SOLID: uploader、source semantics、coordination、UI、requestのownerが分かれている。
- KISS: Clear/Delete/Cancelとglobal confirmationだけである。
- DRY: card/globalが同じsource identity、empty rule、mutation pipelineを使う。
- YAGNI: removal log、new order state、modal dependency、watcher、schemaを追加していない。

終了条件:
- 選択済みoutcomeと総合計画書RM-01〜RM-12を満たす。
- production、tests、docs、generated diffを別々にreviewする。
- 第14節へauthority SHA、base/branch、files、受入ID、commands/results、browser evidence、
  architecture owner/path evidence、未実施項目を記録する。
- 英語のproposed commit titleと短いsummaryを示す。
- push、PR、merge、tag、deployは明示的な許可がない限り行わない。
```

## S3 — Issue #562: Curve default

```text
gbdraw Issue #562のmatch-style concernについて、`CURVE-FRESH-RESET`を選択したdurable
Product authorityの範囲だけでfresh/reset defaultを実装し、既存Session/Galleryの外観保持を
検証してください。sidebar再編とAuto/Show/Hideはこのsessionでは実装しません。
過去の会話は前提にしません。

最初に全て読む:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/LINEAR_UI_559_560_562_MASTER_PLAN_2026-09-21.md
5. docs/internal/PRODUCT_IMPACT_RATCHET.md
6. docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.mdの
   linear.comparison-match-style-default authority
7. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
8. docs/SESSION_COMPATIBILITY.md、docs/REFERENCE/web-app.md
9. 総合計画書第14節と現在のproduction/tests/docs diff

厳格な開始条件:
- authorityのselected outcomeがCURVE-FRESH-RESETであることを確認する。
  KEEP-RIBBONならruntime変更を行わず、このsessionをnot applicableとして記録する。
- 固定実装ブランチfix/linear-ui-559-560-562-20260921はauthority merge commitをancestorに持つ。
  持たない場合はruntimeへ進まず、authority merge後のorigin/devへ同じ固定branchをrebaseする。
- Product authorityとruntimeを同じ未merge candidateにしない。
- branch名、HEAD、upstream、base ancestryを確認し、無関係な変更を保持する。別のruntime
  branchを作らない。

実装する結果:
1. fresh Linear stateと明示的なReset SettingsだけをCurveにする。
2. Circular、CLI、Python APIのdefaultを変更しない。
3. explicit Ribbon/Curveはrecord、layout、comparison option変更で上書きしない。
4. v42以前のSession、settings-only config、Galleryで明示保存された値を保持する。
5. old persisted fieldがmissingの場合はhistorical Ribbonを使い、fresh Curveを遡及適用しない。
6. createDefaultAdv(mode)のfresh defaultとconfig readerのhistorical fallbackを明示的に分ける。
7. pairwise style normalization/validationがconfig.jsとrun pathで重複していれば、
   current-option-values.jsの一つのvalidatorへ収束させる。callerごとのfallbackを隠さない。
8. comparison computation、filter、match geometryそのものは変更しない。

tests:
- createDefaultAdv('linear')、Reset、explicit match-style preservation、legacy missing field Ribbon。
- settings-only Save/Load、v42 fixture、Gallery sessionsのmatch style。
- same-orientation、inversion、dense matchesをCurveでGenerateし、計算結果が変わらないこと。
- Circular UI/defaultに回帰がないこと。
- index.htmlのTitles/Labels/Legend構造、visibility state、Session versionが変わっていないこと。

主な対象:
- gbdraw/web/js/services/session-active-config-contract.js
- gbdraw/web/js/app/current-option-values.js
- gbdraw/web/js/services/config.js
- gbdraw/web/js/services/reset.js
- tests/web/session-active-config-contract.test.mjs
- tests/web/settings-only-session.test.mjs
- tests/web/gallery-session-migration.test.mjs
- 適切な既存comparison browser test

focused commands:
node --test tests/web/session-active-config-contract.test.mjs tests/web/settings-only-session.test.mjs
node --test tests/web/gallery-session-migration.test.mjs tests/web/session-request.test.mjs
npx playwright test tests/web/settings-only-session.playwright.spec.js tests/web/comparison-ui.playwright.spec.js --workers=1 --retries=0
node --test tests/web/architecture-contracts.test.mjs
git diff --check

SOLID/KISS/DRY/YAGNI review:
- fresh default、legacy fallback、enum validationのownerが明示されている。
- sidebar DOM、visibility state、Session schema、render pathを変更していない。
- Auto tri-stateや汎用style componentを先行実装していない。

終了条件:
- 総合計画書MSD-01〜MSD-04を満たす。
- production、tests、docs、generated diffを別々にreviewする。
- 第14節へauthority/base、branch、変更file、受入ID、commands/results、legacy/Gallery/browser
  evidenceを記録する。
- organization/Auto concernsの採否や実装状態をCurve完了条件へ混ぜない。
- 英語のproposed commit titleと短いsummaryを示す。
- push、PR、merge、tag、deployは明示的な許可がない限り行わない。
```

## S4 — Issue #562: Titles、Record Labels、Legendの再編

```text
gbdraw Issue #562のpresentation-control organization concernについて、
`REGROUP-TITLES-LABELS-LEGEND`を選択したdurable Product authorityの範囲だけでsidebarを
再編し、同じstate、保存key、rendererが維持されることを検証してください。Curve defaultと
Auto/Show/Hide semanticsはこのsessionでは実装しません。過去の会話は前提にしません。

最初に全て読む:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/LINEAR_UI_559_560_562_MASTER_PLAN_2026-09-21.md
5. docs/internal/PRODUCT_IMPACT_RATCHET.md
6. docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.mdの
   linear.presentation-control-organization authority
7. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
8. docs/REFERENCE/web-app.md
9. 総合計画書第14節と現在のproduction/tests/docs diff

厳格な開始条件:
- authorityのselected outcomeがREGROUP-TITLES-LABELS-LEGENDであることを確認する。
  KEEP-CURRENT-ORGANIZATIONならruntime変更を行わず、not applicableとして記録する。
- 固定実装ブランチfix/linear-ui-559-560-562-20260921が、そのauthority merge commitを
  ancestorに持つ。持たない場合はruntimeへ進まず、authority merge後のorigin/devへ
  同じ固定branchをrebaseする。
- Product authorityとruntimeを同じ未merge candidateにしない。
- match-style/Auto decisionsの採否を、このsessionの開始条件にしない。
- branch名、HEAD、upstream、base ancestryを確認し、無関係な変更を保持する。別のruntime
  branchを作らない。

実装する結果:
- Titles & Record Labels card内でPlot TitleとRecord Labelsを別subsectionにする。
- Plot TitleはText、Position、Font sizeを持ち、Record Labelsへ混ぜない。
- Definition Font Sizeの表示名をDefault font sizeにする。既存adv.def_font_sizeとauto sizingを使う。
- Record LabelsはName / Species、Subtitle、Replicon、Accession、Length / Coordinatesを並べる。
- visibility controlは現在baseにあるbooleanまたはtri-stateをそのまま使う。mode semanticsを
  このsessionで変えない。
- 各labelの既存size/weight/colorをStyle disclosure内へ移す。
  linear_definition_line_stylesを直接編集し、style objectやvisibilityのcopyを作らない。
- main visibility controlはlabel名の近くに置き、Style disclosure内へ隠さない。
- LegendをColorsやFeaturesと同じhierarchyの独立cardへ移す。
- collapsed summaryをLegend · <current position>とし、form.legendから導出する。
- Legend Box Sizeの表示名をSwatch sizeへ変える。adv.legend_box_size、unit、rendererは変えない。
- disclosure open stateとsummary textをSessionへ保存しない。

変更してはいけないもの:
- fresh/reset Ribbon/Curve default
- Accession/Lengthのvisibility mode、effective result、legacy migration
- Session version、request schema、Worker/Python renderer
- Plot Title、label style、Legendの保存keyと描画結果
- Circular固有のposition optionsとrendering behavior

accessibility/layout:
- native details/summaryまたは既存card patternを使う。
- summaryとinputsに一意なaccessible name、keyboard操作、visible focusを持たせる。
- main visibility controlをStyleの中へ隠さない。
- 390x844でhorizontal overflow、切れたselect、到達不能buttonがないことを確認する。
- stronger colorや装飾だけでsection境界を表さず、headingとspacingを使う。

tests:
- Record Labelsに各visibility/style controlが一つだけ存在すること。
- style編集が同じstate、canonical request、generated SVGへ到達すること。
- Legend position変更がcollapsed summaryへ即時反映し、Swatch sizeが同じkeyを更新すること。
- Plot Title text/position/font、Default font size、line-specific overrideのSave/LoadとGenerate。
- current baseのRibbon/Curve defaultとboolean/tri-state semanticsが変わらないこと。
- Linear/Circular、desktop/390x844、pointer/keyboardを実ブラウザで確認すること。

主な対象:
- gbdraw/web/index.html
- gbdraw/web/js/app/app-setup.js（derived summaryが必要な場合の最小変更だけ）
- 既存のdefinition-line-style-state.js（state semanticsを変える必要がない限り編集しない）
- tests/web/linear-typography.test.mjs / .playwright.spec.js
- tests/web/session-request.test.mjs
- 適切な既存visual-state／accessibility browser test

focused commands:
node --test tests/web/linear-typography.test.mjs tests/web/session-request.test.mjs
node --test tests/web/session-active-config-contract.test.mjs
npx playwright test tests/web/linear-typography.playwright.spec.js tests/web/visual-state-regressions.playwright.spec.js --workers=1 --retries=0
node --test tests/web/architecture-contracts.test.mjs
git diff --check

SOLID/KISS/DRY/YAGNI review:
- organizationはDOM hierarchyだけを所有し、defaultやvisibility policyを決めない。
- existing stateを一か所から編集し、duplicate controls/stateがない。
- summaryは導出値であり、Session/Historyへ保存しない。
- 新UI framework、汎用style schema、renderer pathを追加していない。

終了条件:
- 総合計画書ORG-01〜ORG-06を満たす。
- production、tests、docs、generated diffを別々にreviewする。
- 第14節へauthority/base、branch、変更file、受入ID、commands/results、browser evidenceを記録する。
- match-style/Auto concernsの採否や未完了をorganization完了条件と混同しない。
- 英語のproposed commit titleと短いsummaryを示す。
- push、PR、merge、tag、deployは明示的な許可がない限り行わない。
```

## S5 — Issue #562: Independent Auto visibilityとSession 43

```text
gbdraw Issue #562のrecord-label-auto-visibility concernについて、Accessionと
Length / Coordinatesの独立した
Auto/Show/Hide、effective rendered-row resolution、Session 43 compatibilityを実装し、
request、SVG、Undo/Redo、legacy Loadを検証してください。過去の会話は前提にしません。

このsessionは`INDEPENDENT-AUTO-SHOW-HIDE`がdurable authorityとしてruntime baseへmerge済みの
場合だけ実行する。`KEEP-BOOLEAN-VISIBILITY`では実行しない。

最初に全て読む:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/LINEAR_UI_559_560_562_MASTER_PLAN_2026-09-21.md
5. docs/internal/PRODUCT_IMPACT_RATCHET.md
6. docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.mdの
   linear.record-label-auto-visibility authority
7. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
8. docs/internal/WEB_CHANGE_POLICY.md
9. docs/SESSION_COMPATIBILITY.mdとdocs/REFERENCE/session-and-request-compatibility.md
10. 総合計画書第14節の実施記録、現在のgit diff、Session version inventory

開始条件:
- 固定実装ブランチfix/linear-ui-559-560-562-20260921が、INDEPENDENT-AUTO-SHOW-HIDE authority
  merge commitをancestorに持つ。持たない場合はruntimeへ進まず、authority merge後の
  origin/devへ同じ固定branchをrebaseする。
- match-style defaultとorganization decisionsの採否を、このsessionの開始条件にしない。
- organizationが実装済みならRecord Labels内、未実装なら既存visibility位置でtri-stateを公開する。
- branch名、HEAD、upstream、base ancestryを確認し、無関係な変更を保持する。別のruntime
  branchを作らない。
- rgでCURRENT_SESSION_VERSION、SESSION_VERSION、hardcoded 42、supported versions、
  gallery publication、capture scripts、docs、fixturesを列挙し、更新checklistを実施記録へ置く。
- main first-parentまたはrelease tagでv42 booleansの実在を確認し、positive legacy fixtureを選ぶ。

normative behavior:
- fresh LinearのAccession selected modeとLength selected modeは各auto。
- showは常にtrue、hideは常にfalse、autoはshared rendered rowがなければtrue、あればfalse。
- shared rendered rowはeffective Generate layoutの同じrowに二records以上ある状態である。
- upload File数、File内record数、総record数だけで判定しない。
- Record Layout OFFではdormant row mapを無視し、実際の一record/row layoutとして扱う。
- any shared rowがあればAuto fieldをdiagram-wideで隠す。二fieldsは独立する。
- layout変更はselected modeを書き換えない。Show/Hideを自動上書きしない。
- UIはAuto · Shown / Auto · Hiddenを表示し、manual optionsをdisableしない。
- Undo/Redoはselected modeとlayoutを戻し、effective booleanを再計算する。
- Auto policyはLinearだけに適用する。

state and owner design:
1. current editable stateを次の一組へ収束させる。最新baseに同等keyがある場合はそれを使う。
   adv.linear_accession_visibility = auto|show|hide
   adv.linear_length_visibility = auto|show|hide
2. current writerからlinear_show_accession/linear_show_lengthをeditable stateとして除く。
   legacy reader以外で二つの表現を同期するwatcherやdual-writeを作らない。
3. effective rowsとshared-row factをlinear-record-layout.jsへ移す/公開する。
   session-request.js内部の重複row解決を削除する。layout inputsを変更しないpure functionにする。
4. selected mode validation、legacy bool mapping、effective booleanを一つのpure policy owner
   （計画ではlinear-label-visibility.js）へ置く。UIとrequestが同じresolverを呼ぶ。
5. session-request.jsは既存schema 7のboolean config overrideだけを出す。
   request schema、Worker protocol、Python render modelを変更しない。

Session compatibility:
- Web/Python current writerを42から43へ上げる。
- supported version setには既存versionsを保持し、43をcurrentとして追加する。
- v42以前のown boolean true -> show、false -> hide。
- boolean missing -> historical Show。明示booleanをAutoへ変えない。
- current selected-mode fieldがある場合はlegacy boolで上書きしない。
- v43 current writerはselected modeだけを書き、legacy booleansをconfig stateへdual-writeしない。
- 保存済みcanonical request／Resultのeffective booleanはそのまま保持し、Loadだけでpreviewを
  regenerate/replaceしない。次のGenerateでmigrated modeから同じ外観を再現する。
- settings-only、full Session、Python read/write、Web publication、Gallery migration、capture scripts、
  docsのversionを一致させる。

最低限確認するversion owner/consumer:
- gbdraw/web/js/services/config.js
- gbdraw/session_io.py
- gbdraw/web/js/services/gallery-session-publication.js
- docs/SESSION_COMPATIBILITY.md
- docs/REFERENCE/session-and-request-compatibility.md
- docs/capture/flows/how_to/interactive_sessions.py
- tests/test_api_session.py、tests/test_session_io.py
- tests/test_documentation_contracts.py、tests/test_documentation_reference_contracts.py
- tests/test_composition_surface_contracts.py、tests/test_gui_interactive_capture_contracts.py
- tests/web/gallery-session-publication.test.mjs、settings-only/session request tests
最新rg inventoryを優先し、このlistだけで完了としない。

architecture exception:
- legacy boolean promotionはcompatibility pathである。Architecture Ratchetのexception条件を確認する。
- 必要なら正確なbefore/after OE、PE、CB sets、positive v42 fixture、canonical current owner、
  removal owner、removal conditionを記録する。
- removal conditionは、v42がsupported version setから正式にretireされた後にlegacy promotionを
  削除することを含める。
- checker、supported version contract、fixture assertionsを弱めて通さない。

unit acceptance matrix:
- 1/1/1 records on separate rows: Auto true。
- four records on four rows: Auto true。
- rows 1,1,2: Auto false for Auto field on all records。
- Accession=show + Length=auto: true/false in shared layout。
- Accession=hide + Length=show: false/true in any layout。
- shared -> separate -> shared: selected modes unchanged, derived result toggles。
- layout disabled with dormant shared rows: Auto true。
- invalid selected mode: current validation policyに従う明示error、silent Auto化なし。
- legacy true/false/missing/current-mode precedence。

browser/session acceptance:
- real controlsでAuto/Show/Hideを選び、summary、canonical request、generated SVG textを比較する。
- Record Layoutでshared rowを作り、Generate後のAccession/Lengthの有無をsemantic SVGで確認する。
- layoutを戻し、Autoだけ復帰しmanual Show/Hideが維持されることを確認する。
- History Undo/Redo、Save、新page fresh Load、再Generateを行う。
- v42 positive fixture、settings-only、Gallery examplesが元のRibbon/visibilityを維持する。
- Circular、CLI/API、request schema 7に変更がないことを確認する。

focused commands:
node --test tests/web/linear-record-layout.test.mjs tests/web/session-request.test.mjs
node --test tests/web/session-active-config-contract.test.mjs tests/web/settings-only-session.test.mjs
node --test tests/web/gallery-session-migration.test.mjs tests/web/gallery-session-publication.test.mjs
python -m pytest tests/test_api_session.py tests/test_session_io.py tests/test_documentation_contracts.py tests/test_documentation_reference_contracts.py -v
npx playwright test tests/web/linear-typography.playwright.spec.js tests/web/settings-only-session.playwright.spec.js --workers=1 --retries=0
node --test tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs --base origin/dev --head HEAD
git diff --check

SOLID/KISS/DRY/YAGNI review:
- row ownerはrows、visibility ownerはpolicy、Session ownerはmigration、request ownerはprojectionだけを持つ。
- enum二つとpure resolver一つであり、watcher/state mirrorがない。
- UI/request/SVGが同じeffective resultを使う。
- AutoをCircularや他labelへ一般化せず、request/Python schemaを増やしていない。
- legacy branchは実在v42だけを扱い、dual-writeや推測migrationがない。

終了条件:
- 総合計画書AV-01〜AV-12を満たす。
- Web/Python version constants、supported versions、docs、tests、capture codeが43で一致する。
- production、tests、docs、generated diffを別々にreviewする。
- 第14節へauthority/base、branch、version inventory、positive fixture、変更file、受入ID、
  commands/results、browser wheel identity、architecture exception evidence、未実施項目を記録する。
- 英語のproposed commit titleと短いsummaryを示す。
- push、PR、merge、tag、deployは明示的な許可がない限り行わない。
```

## S6 — 統合browser受入、full gates、最終handoff

```text
gbdraw Issues #559/#560/#562について、選択・実装された範囲の統合受入、full gates、
architecture/Product evidence、最終diff reviewを完了し、in-scope不具合を修正してください。
過去の会話は前提にしません。

最初に全て読む:
1. AGENTS.md
2. CLAUDE.md
3. gbdraw/web/CLAUDE.md
4. docs/internal/LINEAR_UI_559_560_562_MASTER_PLAN_2026-09-21.md
5. docs/internal/LINEAR_UI_559_560_562_INSTRUCTION_PROMPTS_2026-09-21.md
6. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
7. docs/internal/PRODUCT_IMPACT_RATCHET.md
8. docs/internal/WEB_CHANGE_POLICY.md
9. runtime baseにある#559/#560/#562のapplicable Product authority
10. 総合計画書第14節、git history、現在のproduction/tests/docs/generated diff

開始確認:
- branchがfix/linear-ui-559-560-562-20260921であること、HEAD、upstream、base ancestryを確認し、
  無関係な変更を保持する。別のruntime branchを作らない。
- #559はPD-OI-025、#560と#562三concernsは実際のselected outcomeとauthority SHAを個別に特定する。
- CURVE-FRESH-RESETならMSD、REGROUP-TITLES-LABELS-LEGENDならORG、
  INDEPENDENT-AUTO-SHOW-HIDEならSession 43とAV受入を必須にする。非採用concernのruntime受入を
  要求しない。
- 各sessionの未実施、失敗、temporary evidence、browser wheel identityを一覧化する。
- source/input/environment/acceptanceが変わっていない有効なevidenceは再利用し、変更・失敗・
  未解決箇所だけを無条件に重複実行しない。

integrated primary journey:
1. 現在sourceからpython tools/prepare_browser_wheel.pyでwheelを準備し、sourceとの対応を記録する。
2. fresh Linear pageでtwo Filesを追加する。一方はmulti-record GenBank、必要なら他方はpaired
   GFF3+FASTAとし、same-name independent source caseも別に持つ。
3. File-level Depth disclosureをclose/openし、record list closedのままDepth seriesを追加する。
   common apply、per-record Mixed、File-level clear、Undoを行う。
4. global Add/Removeとcard Removeをpointer/keyboardで操作する。選択済み#560契約に従って
   Clear/Delete/Cancel、sole slot、blank last slot、middle replacementを確認する。
5. match-styleとorganizationは、それぞれのselected outcomeを独立に確認する。
   CURVE-FRESH-RESETならfresh/reset Curveとhistorical Ribbonを確認する。
   REGROUP-TITLES-LABELS-LEGENDならTitles/Record Labels、Style disclosures、Legend summary/
   Swatch sizeを確認する。KEEP outcomeなら現在の対応する結果を維持する。
6. INDEPENDENT-AUTO-SHOW-HIDEが選択されている場合はAccession/LengthのAuto/Show/Hideを
   独立設定し、separate rows、shared row、layout disabled、Undo/Redoを操作する。
7. Generateごとに完了を待ち、canonical request、current Result identity、SVG semantic textを
   対応付ける。古いResultを新しいGenerateの証拠にしない。
8. Saveし、新しいpageでfresh Loadする。selected modes、source order、blank/deleted result、Depth
   matrix、match style、label styles、Legend、保存Resultを確認して再Generateする。
9. v42 positive Sessionと代表Gallery sessionをloadし、historical Ribbonと明示visibilityを確認する。
10. failed/canceled/stale Generateでlast successful Resultが保持される既存contractを確認する。

visual/accessibility review:
- viewport 390x844とdesktopでsidebarの横overflow、control clipping、scroll到達性を確認する。
- 実装されたdisclosure、dialog、buttons、tri-state selects、Style、Legendをpointer、Enter、Space、
  Escape、Tabで操作し、visible focus、focus order、focus return、unique accessible namesを確認する。
- CURVE-FRESH-RESETの場合はCurveをsame-orientation、inversion、dense match setで目視し、
  filter/computationが変わっていないことをrequest/result evidenceで確認する。
- REGROUP-TITLES-LABELS-LEGENDの場合はRecord Labelsのvisibility/styleとPlot Titleが概念的にも
  DOM hierarchyでも混ざらず、Legendが独立sectionであることを確認する。
- screenshot差分だけでbehaviorを合格にせず、request/Session/SVGも確認する。

Gallery/docs:
- 既存Gallery sessionsは保存外観を維持する。fresh defaultを適用するためだけにsessionを再保存しない。
- Gallery tutorial text/screenshotsを変更する場合だけ、repositoryの
  web-gallery-screenshot-maintenance skillを読み、その再現手順に従う。
- 公開figureが変わる場合はrealistic Gallery-quality recipeから生成し、readable scaleで目視する。
- examples/gbdraw_social_preview.pngは変更しない。
- browser wheel、temporary screenshots、test artifactsをcommitしない。

full verification:
node --test tests/web/*.test.mjs
node --test tests/web/architecture-contracts.test.mjs
node --test tests/ci/*.test.mjs
npm run test:web:comparison-contracts
node tools/check-web-change-budget.mjs --base origin/dev --head HEAD
git diff --check

Session 43/Python I/Oを変更した場合:
python -m pytest tests/test_api_session.py tests/test_session_io.py tests/test_documentation_contracts.py tests/test_documentation_reference_contracts.py -v

Python側failureやimpact policyが広いgateを要求する場合:
pytest tests/ -v -m "not slow"

長時間testは少なくとも30分を許容し、増分監視する。timeout/assertionを弱めない。
remote CIをmonitorする場合はrepositoryのpoll間隔を守る。

最終diff review:
1. production
   - Depthは既存matrix/action、removalはsource owner+one coordinator、Autoはrow owner+one policy、
     Sessionはone writer/migration、requestはone projectionになっているか。
   - watcher、dual-write、parallel state、fallback、dead action、duplicate controlがないか。
2. tests
   - internal methodだけでなくreal DOM journeyを通るか。
   - fresh、Reset、legacy Load、current Save/Load、Undo/Redo、Generateを区別するか。
   - request、Result、SVG、Sessionのどれを測るtestか明確か。
3. docs
   - web reference、Session compatibility、master plan、authorityが実装結果と一致するか。
   - 過去の会話を知らない読者に意味が通るか。
4. generated artifacts
   - wheel、screenshots、reference SVG、egg-info、distに意図しない変更がないか。

architecture/Product evidence:
- #559はPD-OI-025/OIC-020へtraceする。
- #560/#562はruntime baseのselected outcomeへtraceし、runtime diffでauthorityを変更しない。
- ordinary non-increasing変更はowner/pathとsuperseded pathを簡潔に示す。
- Session 43 compatibilityなどexception条件は正確なOE/PE/CB、positive fixture、removal conditionを
  示す。architecture checkerを緩めない。
- testsやIssue本文をProduct authorityと呼ばない。

終了条件:
- 総合計画書第13節と、選択された全DPT/RM/MSD/ORG/AV IDが成立するか、未成立の境界を
  具体的に示す。
- browser、Session/Gallery、full gates、architecture/Product evidence、final diff reviewを完了する。
- 検証で見つかったin-scope不具合を修正し、その変更で無効になったevidenceを更新する。
- 第14節へ最終base/branch、authority、source/wheel、commands/results、acceptance matrix、limitations、
  rollback、変更fileを記録する。
- 全条件成立時だけ実装完了と報告する。
- repository guidanceに従い、英語のproposed commit titleと短いsummaryを示す。
- commit前にbranch/upstreamを確認する。push、PR、merge、tag、deployは明示的な許可がない限り
  行わない。PR wordingを作る場合はwrite-clear-pull-request skillとlanguage checkerを使う。
```
