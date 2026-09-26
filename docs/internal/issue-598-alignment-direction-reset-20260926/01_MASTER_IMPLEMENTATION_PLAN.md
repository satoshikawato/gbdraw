# Issue #598 修正実装総合計画書

- Author / Product Decision Owner: `satoshikawato`
- 作成日: 2026-09-26
- 対象: https://github.com/satoshikawato/gbdraw/issues/598
- 実装ブランチ: `fix/issue-598-alignment-direction-reset-20260926`
- 初期ベース: 作成時の最新 `origin/dev`、`d457b7189b137185a8dec800819a312c30b969fa`
- 決定更新ブランチ: `product/issue-598-decisions-20260926`
- 決定更新コミット: `eba01dd518ff6edd4effb4b00d519cbd78ff8623`（Contract単独変更）
- 状態: Productの選択は承認済み。本文は実装指示であり、実装完了の報告ではない。

## 1. 問題と修正範囲

gbdraw Webは複数のゲノムレコードを線形図に配置し、Similarity Group内の選択フィーチャーの中心を水平方向に揃えられる。referenceは位置合わせの基準として選ばれた正確なフィーチャーであり、グループ内で最も多い方向を代表するとは限らない。

現行のレビューは `Match reference direction` という単一のチェックボックスを持つ。適用に伴ってレコードが反転しても、`Reset Align` は位置だけを戻す。Issue #598のBUG-03とBUG-04に対応し、表示方向を明示的に選べるレビューと、方向の復元範囲を選べるResetを実装する。

**BUG-17は対象外。線形フィーチャーを表示左端に置く新しいクロップ操作、線形回転の代替機能、関連する新規テスト・例・チュートリアルを追加しない。** 既存の線形region編集と円形Rotateの動作は維持する。

## 2. 決定と実装の順序

正式なProduct authorityは [Option Integrity Product Contract](../OPTION_INTEGRITY_PRODUCT_CONTRACT.md) である。[承認文保存記録](00_ACCEPTED_PRODUCT_DECISIONS.md) は、関心別に分けた4本の承認済みDecision Packsへの索引である。各Packは対応する完全なowner receipt、責任、受入条件、実装routeを保存する補助資料であり、別のdecision storeではない。

| 決定 | scenario revision | 承認済みchoice | 責任 |
| --- | ---: | --- | --- |
| PD-OI-027 | 5 | A / EXPLICIT_DISPLAY_DIRECTION_MODES | 表示方向と幾何変換 |
| PD-OI-029 | 3 | A / SELECTABLE_RESET_WITH_ALIGNMENT_DIRECTION_RESTORE | ResetとHistory |
| PD-OI-031 | 5 | A / AUTO_APPLY_WITH_EXPLICIT_DIRECTION_REVIEW | 操作面と自動Align |
| PD-OI-034 | 5 | A / LOCAL_EXCLUSIVE_DIRECTION_REVIEW_WITH_RETRY | draft、検証、再試行 |

Product Contract revision 17の更新は決定更新ブランチで独立して提出する。[Product Impact Ratchet](../PRODUCT_IMPACT_RATCHET.md) と [Web Change Policy](../WEB_CHANGE_POLICY.md) により、candidate authorityは同じcandidate runtimeを認可できない。**4件の正確な決定が `origin/dev` に取り込まれてから、実装ブランチへそのdevをmergeし、S01以降を開始する。** 承認済みテキストがwork branchにあるだけでは開始条件を満たさない。実装セッションは決定更新ブランチを実装ブランチへ直接mergeしない。devへのmergeやPR公開は、この計画によって自動的には許可されない。

PD-OI-026/028/030/032/033/035を含む他の決定は維持する。特に、正確なreference選択、Pythonによる候補適格性判定、決定的なアンカー解決、計画の再生成、CLIの曖昧性拒否、既存の狭幅palette coverageの契約を変更しない。

## 3. 操作仕様

### 3.1 Alignレビュー

方向は単一のradio groupで選ぶ。複数モードは同時に成立しない。

| モード | 結果 |
| --- | --- |
| Keep current directions（既定） | すべて現在の方向を維持 |
| All selected features right-facing → | 選択アンカーが既知strandのレコードを、アンカーが表示上右向きになる方向へ変更 |
| All selected features left-facing ← | 同じ対象を、アンカーが表示上左向きになる方向へ変更 |
| Custom | referenceを含む各対象でKeep / right-facing / left-facingを選択 |

「all」の範囲は正確なreferenceとAlignに参加する選択済みtargetアンカーであり、入力中の全レコードやレコード内の全遺伝子を同じstrandに書き換える意味ではない。変えるのはレコード全体の表示方向であり、生物学的なsource strandは編集しない。strand不明・Skip・missing・unusableは変えず、理由を表示する。不明strandを多数派やreferenceから推測しない。

referenceが一つだけ左向きで、他が右向きなら、all-rightでreferenceレコードだけが反転する。特定レコードだけ変える一般的な操作にはCustomを使う。別の「少数派だけ反転」モードは作らない。

自動解決できる通常のAlignはKeepで自動適用し、レビューを強制しない。曖昧性がある場合、または利用者がレビューを開いた場合に選択肢を示す。各行にrecord名、選択アンカー、現在→結果の矢印、適用対象外の理由を出す。referenceの反転によるレコード左端の移動も予告する。

### 3.2 位置の意味

Alignは、実行直前のreferenceフィーチャー中心の論理canvas xを保ち、target中心をそこへ揃える。referenceレコード左端のx固定ではない。referenceを反転すると中心がレコード内で移るので、その分レコード配置を補正する。すべてのレコードのyは保つ。DOMのscreen座標、zoom、CSSサイズ、fit-to-canvas後のviewBox余白を意味の基準にしない。

通常のReverseは独立したレコード操作であり、active planを方向policyとして解釈しない。再生成時の整列中心は現在のrecord方向から計算する。上記の「直前のx固定」は各Align操作の確定時の契約であり、以後の表示倍率や異なるlayout設定まで固定する機能ではない。

### 3.3 Reset

`Reset alignment…` で次のexclusiveなscopeを選択する。通常scopeはpositions-onlyとする。

| scope | 位置 | 方向 |
| --- | --- | --- |
| Reset positions | 最新の成功したAlignの直前配置へ戻す | すべて現在の方向を維持 |
| Reset positions and alignment direction changes | 同じ直前配置へ戻す | そのAlignで実際に反転したレコードだけ、直前の絶対方向へ戻す |

combined scopeは対象名・件数・現在→復元方向と「対象でAlign後に行った手動方向編集も置き換える」ことを操作前に示す。Alignで方向が変わらなかったレコードの後続手動Reverseは保持する。referenceはそのAlignで実際に反転した場合だけ対象となる。反転回数を数えてtoggleし直す方法は採らない。

両scopeともplanと復元receiptを消費する。positions-only実行後にcombinedへ切り替えるには、そのResetをUndoする。最新Alignの反転差分が空ならcombinedは無効で理由を示す。古いSessionに信頼できる方向の履歴がなければpositions-onlyを使え、combinedは理由付きで無効にする。source方向や最初のAlignから推測しない。

Align Aの後にAlign Bを実行した場合、Bの直前配置を記録し、receiptはBのものに置き換える。Reset BはB直前の配置に戻してactive planを消す。Aのplanを再有効化する操作ではない。Undo Bは完全な前artifactを復元するため、Aのplanも復元できる。

## 4. 実装アーキテクチャ

### 4.1 所有者と経路

| 意味・責任 | 既存所有者／拡張箇所 | 方針 |
| --- | --- | --- |
| source identity、候補の適格性、strand・centerの確定 | `gbdraw/layout/similarity_alignment.py`、`gbdraw/api/record_planning.py`、既存Web adapter | Pythonを唯一のdomain判定元とする |
| レビューの方向intentとpreview投影 | `gbdraw/web/js/app/similarity-alignment.js` と必要最小限のprivate helper | 単一resolverでPython factsから結果を投影 |
| 現在方向 | 既存record region/presentation状態とrecord-display owner | absolute directionを同じsetterで適用 |
| plan | typed `SimilarityAlignmentPlan` と canonical request layout | 正確なanchorsと配置関係のみ。方向policyを持たない |
| base translationと描画配置 | canonical `recordTranslations`、`linear_multi_record.py`、`diagrams/linear/assemble.py` | 既存placement計算を再利用。別render engineを作らない |
| rendered配置からtranslationへの橋渡し | `app/legend-layout/composition-actions.js::materializeRecordTranslations` | keyed metadataの既存bridgeのみを使う |
| request構築・admission・Worker | `services/session-request.js`、既存candidate execution/admission | 既存canonical pipelineを使う |
| Apply / Resetの確定 | `app/run-analysis.js::runCommittedCanonicalCandidate` と生成artifact transaction | 確定済みcanonicalへの限定変更として原子的に確定 |
| Session / History | `services/config.js`、`history-snapshot.js`、`history.js` と既存artifact capture/restore | 復元receiptを既存capture対象へ追加。第二Historyを作らない |

`runAnalysis` にlive formを再構築させる経路をAlign/Resetの専用確定経路として残さない。未適用form編集、target外設定、比較結果を保持するため、最後に確定したcanonical artifactを起点に、対象plan・方向・translationだけを変更する。candidate実行、SVG sanitizer/admission、preview readiness、History finalizationは既存の共通所有者へ委譲する。既存経路が欠く最小入力はその所有者で拡張し、別のApply/Reset pipelineを作らない。

### 4.2 一つの方向draft

概念型は次のtagged unionとする。これはtransientなreview stateであり、Sessionやplanに保存しない。

```text
DirectionIntent =
  { mode: 'keep' }
  | { mode: 'right' }
  | { mode: 'left' }
  | { mode: 'custom', byRecordKey: Map<RecordKey, 'keep' | 'right' | 'left'> }
```

Custom以外のbranchに行ごとのoverrideを持たせない。candidateやSkipを変えると同じresolverで対象とpreviewを再計算し、除外レコードの残ったrow値を適用しない。known displayed strandが希望する矢印と異なる時だけ、現在のabsolute reverse値を反転した値に投影する。region reverseとpresentation reverseの解釈を第二の式として複製しない。

方向・candidate・Skipの編集はローカルで行い、Workerを呼ばない。Applyの試行ごとに既存Python batch validationを一度行う。最終factsで方向またはreference中心補正がpreviewと変われば、previewを更新してレビューに戻り、次のApplyを待つ。黙って新しい出力を確定しない。自動解決直後に既に得たfinal validation結果は同一bindingの間で再利用し、冗長な同一validationを追加しない。

### 4.3 reference中心補正と位置baseline

既存planはreferenceの現在中心を用いてtargetを揃える。方向反転前の中心を固定するため、次の補正をcanonical record translationへ一度だけ適用する。

```text
x_before = placement_before.x_for_position(center_before) + base_before.x
base_after.x = x_before - placement_after.x_for_position(center_after)
reference_position_delta = base_after.x - base_before.x
```

ここでplacementは既存linear layoutの論理canvas geometryである。centerはPythonのsource-bound projectionから得る。単純な `length - center` をcrop・複合feature・配置の代わりにJSで実装しない。必要なfacts/geometryは既存projection/placement ownerから最小限取り出す。previewと最終candidateが同じ計算を使うことを検証する。

base translationは現在のcanonical ownerに一つだけ置く。Reset receiptにはreference補正が実際に生じた時の有限なsparse Δxだけを補助履歴として保存し、positions Resetはそれを取り消してbase_beforeを復元する。全レコードtranslationの第二snapshotは作らない。targetの整列offsetは既存planによって計算されるため、planを消すとその直前baseへ戻る。再Align時は既存materialization bridgeで現在の実配置を新baseにし、旧Δxを引き継がない。

この表現のbinding/admission/Save・Load・Undoの一貫性はS00/S01で実証する。既存contractでは正しく保持できない場合、独自のfallbackやpolicy fieldをplanへ足す前に、既存translation owner内で同じ結果を実現する最小の表現変更を検討する。ユーザーの復元位置を変える妥協は認めない。

### 4.4 復元receiptとライフサイクル

概念上の最小receiptは、active plan/source binding、実際に変わったrecordKeyごとの絶対方向before/after、および必要なreference placement Δxから成る。フィールド名はS02で既存Session capture形に合わせて固定する。preview時には保存せず、成功して確定するartifactの実際のbefore/afterから作る。生成結果の成功だけでなく、preview readinessとHistory finalizationまで成功した状態を確定とする。

receiptはdirection selectionを描画するデータではなく、Resetが明示的に変更するための履歴情報である。Rendererはreceiptを方向指定として読まない。plan/requestにmode、Custom値、effective orientationを追加しない。

| 事象 | plan / receiptの扱い |
| --- | --- |
| 新しいAlign成功 | 一緒に置換。beforeはそのAlignの直前 |
| style再生成、identityを保つreorder | source/anchor bindingを検証して保持 |
| 普通のReverse | 現在方向を更新。receiptの絶対before方向は書き換えない |
| source/region/selector変更、手動位置移動などのplan invalidation | 既存invalidation条件に従い一緒に消す。必要なtranslation materializationは既存ownerで実行 |
| Reset成功、明示的clear | 一緒に消す |
| Undo / Redo | SVG、canonical、位置、方向、plan、receipt、関連資源を一緒に復元／再適用 |
| validation/render/readiness/finalization失敗、Cancel、stale、superseded | prior artifactとreceiptを保持し、Historyを追加しない |
| Session Save / fresh Load | 同じbindingと値を保持。壊れたreceiptはadmissionで拒否し、部分的に消して成功扱いにしない |

「古いSessionにreceiptがない」と「新形式receiptが壊れている」を区別する。前者は既存の位置baselineを使い、方向復元不可を明示する。後者はfresh Loadの候補を拒否し、現在のartifactを維持する。新しいreference Δxを持つartifactでは、その履歴だけを落としてpositions Resetが誤る状態を受け入れない。compact bindingには既存のsource fingerprintと正確なplan/anchor identityを再利用し、日時やrow indexをidentityにしない。

### 4.5 compatibilityの制約

調査ベースではSession writer=44、canonical request=8、plan schema=2、catalog=4であり、first-parent mainでwriter=44の公開証拠は見つかっていない。これは2026-09-26の初期観測であり、未公開の最終判定ではない。S00は実際の最新main/release tagsとfixturesを調べる。

公開済みの形だけに、既存admission ownerで必要最小限のcompatibilityを提供する。devだけの中間形式にmigrationを追加せず、fixtures/内部docsを現行表現へ書き換える。Session versionを上げる必要があるかは、この公開証拠と既存unknown-field/admission契約を根拠に決める。独立したreceipt schema/version registryは追加しない。request=8とplan=2の方向非依存構造は維持し、方向policyのためのversion更新をしない。

## 5. SOLID / KISS / DRY / YAGNIの実行基準

- **SRP**: Python domain判定、review intent投影、record transform、artifact確定、Session履歴の責務を分ける。ファイル分割自体を目的にしない。
- **OCP・LSP**: 既存のtyped factsとcanonical候補契約を最小限拡張し、Keep/通常Reverse/CLIの既存契約を保つ。名前だけ異なる別render pathを作らない。
- **ISP・DIP**: controllerは必要なvalidationとcanonical transactionのインターフェースを受け、Worker、DOM admission、History内部を再実装しない。将来用の汎用strategy/plugin frameworkは追加しない。
- **KISS**: 方向は4モードの一択、Resetは2scope、Applyは一つの確定経路。多数派推測や重複するbulk/override組み合わせを作らない。
- **DRY**: previewとApplyの方向・位置計算、ApplyとResetのartifact確定、SessionとHistoryのreceipt検証を既存ownerで共有する。廃止するMatch経路・state・tests/docsを同じ機能変更で除去する。
- **YAGNI**: BUG-17、CLI方向flag、新しいorientation/render/History owner、decision store、review policy永続化、dev-only migrationsを追加しない。

[Architecture Ratchet](../ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md) に従い、通常の非増加変更はowner/pathのbefore/afterと廃止経路を簡潔に示す。例外条件に当たる時だけ完全なOE/PE/CB集合と承認経路を提示する。候補を通すためにchecker、authority、workflowを弱めない。

## 6. セッション分割と実行順序

同じremote branchへの実装writerは一度に一セッションとする。順序は固定し、各セッションのコミット・プッシュ・remote SHA確認後に次を開始する。途中HEADは作業用であり、S04完了まで配布しない。

| セッション | 指示ファイル | 所有責務 | 完了条件 |
| --- | --- | --- | --- |
| S00 | [00_PREFLIGHT.md](sessions/00_PREFLIGHT.md) | authorityのdev反映確認、隔離、公開format証拠、baseline、owner/path intake | 実装開始条件と証拠を保存してpush |
| S01 | [01_DOMAIN_PROJECTION.md](sessions/01_DOMAIN_PROJECTION.md) | Python facts/placementを使う方向・中心補正のpure projection | 式、bindings、幾何の意味あるtestsをcommit/push |
| S02 | [02_ATOMIC_ALIGNMENT_AND_RESET.md](sessions/02_ATOMIC_ALIGNMENT_AND_RESET.md) | Align/Reset、receipt、Session/History、機能する最小UIの一括統合 | 両機能と保存復元を一経路で通してcommit/push |
| S03 | [03_REVIEW_UI_AND_ACCESSIBILITY.md](sessions/03_REVIEW_UI_AND_ACCESSIBILITY.md) | review/reset UI、失敗再試行、desktop/390px/keyboard実証 | 操作仕様と既存palette契約の両立をcommit/push |
| S04 | [04_DOCUMENTATION_AND_FINAL_VERIFICATION.md](sessions/04_DOCUMENTATION_AND_FINAL_VERIFICATION.md) | public docs、必要なGallery更新、全体検証、最終review | 最終コード・docs・成果物の証拠とhandoffをcommit/push |

S00がauthority未mergeを検出したら、独立したbaseline/format調査と文書のみ完了し、commit/pushして開始条件未達を報告する。S01以降を進めるためにcandidateを自己承認したり、未承認のdevへのmergeを行ったりしない。新たなmaterial Product outcomeが必要になった場合は影響範囲だけを停止し、関心ごとにDecision Packを分ける。実装上のroutineな選択について再承認を求めない。

## 7. 他セッションに干渉しない手順

各INSTRUCTION PROMPTに専用clone作成とbranch取得手順を含める。作業中の共有checkoutのbranchをswitch/resetしない。既存のuncommitted changesを取り込まない。実装ブランチを最新devから新規に作り直すのではなく、**この計画を含む指定のremote実装ブランチをfetchして使う**。authority反映後のdev取り込みはS00が通常mergeで行い、履歴をforce-pushしない。

専用cloneまたはbranch競合のない専用worktree、checkout内のPython environment、個別server portを使う。共有editable install、別セッションのwheel/fixtures/server/processを変更しない。NodeのPlaywrightとPythonのPlaywrightを両方確認し、Node packageがなければPythonで必要なブラウザ検証を行う。Chromium sandbox障害は承認されたsandbox escalationで同じ検証を再実行する。

各セッション末尾で、実装・tests・docs・生成物のdiffを別々に確認し、対象pathだけstageして英語messageでcommitする。同名remote branchだけへpushし、HEADとremote SHAの一致を確認する。remoteが進んでいたら他writerの終了状況を確認して正規に取り込み、force-pushしない。dev/mainへ直接pushしない。

## 8. 受け入れ条件と証拠

以下はすべて必要な独立条件である。同じchoice IDであることをもって、一部条件の不足を代替しない。`evidence/` に再実行可能なcommand、環境、fixture、結果、commit SHAを記録し、過去の実行を今回のHEADの結果として表示しない。

| ID | 条件 | 主担当 |
| --- | --- | --- |
| A01 | Keep自動Alignとexact reference選択が既存どおり。CLI曖昧性拒否とdefaults維持 | S01/S02 |
| A02 | right/left/Customが一択。reference少数派、referenceだけ反転、targetだけ反転、sourceの+/-の両方 | S01/S03 |
| A03 | unknown/Skip/missing/unusableは方向不変、理由とbefore/after矢印が一致 | S01/S03 |
| A04 | source bytes/identity/strand不変、record全体のfeature/label/annotation/ribbon整合、text可読 | S01/S04 |
| A05 | reference中心x固定、全y固定、左端補正、compound/crop/逆方向/不等長/非ゼロbase配置 | S01/S02 |
| A06 | local draft編集はWorkerゼロ。Applyのfinal batchは一回、変更factsは再previewと別Apply | S02/S03 |
| A07 | positions Resetは直前base、現在方向維持。combinedは実際のAlign反転対象だけ絶対beforeへ | S02 |
| A08 | affected/unaffectedの後続手動Reverseを区別。reference含有、no-op差分、Align A→B、両Reset消費 | S02/S03 |
| A09 | Save→fresh Load→両Reset、style regeneration、reorder、binding破損拒否、公開済み旧Sessionの理由付き制限 | S02/S04 |
| A10 | Undo/Redoはplan/receipt/方向/位置/SVG/資源を一緒に復元。pending formとtarget外設定維持 | S02/S04 |
| A11 | validation/render/preview readiness/History finalization失敗、Cancel/stale/supersededでartifact不変・History増加なし、retry editable | S02/S03 |
| A12 | Resetの追加LOSAT jobはゼロ。比較結果の安全な再利用と方向変更後のribbon geometry | S02/S04 |
| A13 | desktopと390pxでradio/custom/resetが操作可能。keyboard、focus、名前・理由・scope読み上げ、既存palette契約維持 | S03 |
| A14 | Match flag/旧経路/古い文言の除去、単一owner/path、architecture/Product gates合格 | 全セッション |
| A15 | BUG-17機能・資料は追加なし。既存円形Rotate/線形region操作は回帰なし | S04 |
| A16 | 実操作で再生成できるdocs/必要なGallery例。social previewは変更なし | S04 |

過去のテスト結果を実装後の合格証拠に流用しない。S00で現環境と最新ベースのbaselineを記録する。既存候補testsは `tests/web/similarity-alignment-actions.test.mjs`、`session-request.test.mjs`、`record-display-options.test.mjs`、`feature-record-rotation.test.mjs`、History群、`similarity-alignment-ui.playwright.spec.js` と `tests/test_similarity_alignment_web_adapter.py`。変更箇所に応じて関連typed/layout testsを選ぶ。新しいtestsは振る舞い境界を保護し、単なる式や実装の複写は避ける。

最終ゲートは関連Node/ブラウザtests、`pytest tests/ -v -m "not slow"`、`ruff check gbdraw/`、必要なoutput comparison、`node tools/check-web-change-budget.mjs --base origin/dev --head HEAD`。architecture-bearing変更はpolicyの適切なprofileとevidenceも実行する。checkerはGate/Review本文を必ず確認し、shell exitだけで合格としない。slow境界や全buildは未解決リスクがある場合に追加し、根拠なく全チェックを反復しない。

テストは30分以上の実行を許容して増分監視する。CI pollは5分以上間隔。`tests/reference_outputs/` は通常read-onlyとし、意図したgeometry差分だけreview後に専用commandで更新する。`dist/`・egg-infoを手編集しない。ブラウザwheelは必要な場合だけ専用checkoutでprepareしgitに入れない。deploy準備でない作業にcache-bust刷新を加えない。

## 9. 最終成果物

実装ブランチに完成コード、適切なtests、再現できるdocs、セッションごとの証拠、最終handoffを保存し、各セッションでcommit/pushする。handoffは実装結果、未解決事項、適用された4決定、検証と限界、remote SHAを自立した文章で記す。

Proposed final commit title: `Fix alignment direction choices and reset restoration`

Proposed final summary: `Add exclusive display-direction review and selectable reset scope through the canonical artifact transaction. Preserve source annotations, reference anchor placement, Session/History continuity, and unrelated edits.`

この英語タイトル・要約は最終handoff用であり、複数セッションの履歴を無断squashする指示ではない。PRの作成・本文変更が別途許可された場合は `.agents/skills/write-clear-pull-request/SKILL.md` を読み、同一本文でlanguage checkerを実行する。
