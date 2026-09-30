<!-- Raw design report of workstream W1a (sub-agent output, 2026-09-30, base origin/dev 4c89bab1). Evidence for 01_REMEDIATION_PROPOSAL.md, not Product or architecture authority. -->

# W1a 修正提案: 複数 Result のエディタ状態と label/visibility 編集（dev 4c89bab1）

対象は FE-01、FE-02（PV-09 を含む）、FE-03、FE-04、FE-10 の 5 件。パスはリポジトリ直下から書き、行番号は 4c89bab1 のもの。コードの変更やファイルの作成はしていない。

## PR 計画（要約）

| PR | 対象 | 分類 | 本番ファイル数 / 変更量 | プロファイルと Review |
|---|---|---|---|---|
| PR-1 | FE-01 | IMPLEMENT_EXISTING_AUTHORITY | 8 / gross 70–100、net −25〜−40 | Ordinary。Session writer の field を削除するので Review REQUIRED（compatibility-path change） |
| PR-2 | FE-04（FE-10 を同梱してよい） | IMPLEMENT_EXISTING_AUTHORITY | 2（+3）/ 約 60 | Ordinary。新 export 1 件が review signal になる |
| PR-3 | FE-03 | IMPLEMENT_EXISTING_AUTHORITY | 3 / 30–45 | Ordinary。reactive 宣言 1 件が review signal になる |
| PR-4 | FE-02 と PV-09 | PRODUCT_DECISION_REQUIRED（決定後に着手） | 5–7 / 300–600 | Architecture。`architecture-change` label を付け、Review REQUIRED |

順序: PR-1 は他に依存しないので最初に出す。PR-2、PR-3 を先に入れてから PD を決め、PR-4 に進む。

---

## FE-01（P1）label の文字と表示の override が黙って消える

**1. 分類: IMPLEMENT_EXISTING_AUTHORITY**
- 根拠:
  - OIPC-C06: 置換と削除は明示操作だけで起きる。inactive な mode や生成の失敗は削除を意味しない。
  - OIPC-C05: 有効な値は、表示中の画面に editor が無いという理由で削除しない。
  - web CLAUDE.md の live-edit 条項: 「Route the mutation through the owning editor action so History records the same operation」。今回の消去はこれに当たらない。
- 補強（authority ではない）:
  - コードのコメント `gbdraw/web/js/app/feature-editor/label-actions.js:683-684`「Keep its label intent dormant」。
  - 内部文書 `docs/internal/SESSION05A8_R_SOURCE_RECONCILIATION.md:69-72`: feature に結び付いた label の map は「preserve its existing owner」。
- source 置換（別の genome を読み込んで Generate）のときの扱いを決めた authority は無い。PR-1 では現在の結果を変えず、発火する場所だけを移す。細かい扱いを変えるなら任意の PD にする（4. を参照）。

**2. 根本原因の再確認**
- 監査の指摘どおり。`syncLabelEditor`（`label-actions.js:689-703`）は、mounted SVG の feature-ID の集合（`buildContextKey` `:342-351`）が変わったことを「genome の置換」の代わりの指標として使い、`clearOverrides`（`:396-402`）を呼ぶ。
- 監査への補足:
  - `labelTextBulkOverrides` と `labelTextFeatureOverrideSources` も消える。bulk override しか無い場合、`hasFeatureScopedOverrideInSvg`（`:353-364`）は feature 側の map しか見ないため、Result を切り替えるだけで必ず消える。
  - 経路がもう一つある。live で非表示にした後の自動 label reflow で、reflow の SVG から hidden feature が消える。コードから判断したもので、再現はしていない。
  - record selector の往復は、コード経路としては確認できる。専用の spec の出力は無い。
- source 置換を判定する明示的な述語は既にある。`sourceReplaced`（`gbdraw/web/js/app/run-analysis.js:1952-1957`、#510/#511）。

**3. 修正案（コード）**
- owner は `label-actions.js` の `syncLabelEditor`。
- context-key による判定を削除する。代わりに `syncLabelEditor({…, sourceReplaced = false})` とし、`sourceReplaced && !hasFeatureScopedOverrideInSvg(svg, …)` のときだけ `clearOverrides()` とダイアログを閉じる処理を行う。
  - 既存の述語をそのまま使い、消去の順序も「消去してから投影する」のまま変えない。発火するのは、成功した source 置換 Generate の bind のときだけになる。
- `run-analysis.js`:
  - `:1952` の `hasSourceBoundEditorIntent` に label の map（text、bulk、visibility）を加える。
  - `:4913` の `bindingOptions` に `sourceReplaced` を加える。
  - reflow（`isReflow`）と target-record-transform（`:5347`）では false のままにする。
- `gbdraw/web/js/app/app-setup.js:2284-2291` から `context.bindingOptions.sourceReplaced` を渡す。
- `labelOverrideContextKey` の経路をすべて削除する（CW-05: consumer を消したら、producer と保存状態も消す）。
  - `gbdraw/web/js/state.js:464,965`
  - `app-setup.js:369,4798`
  - `gbdraw/web/js/services/history-snapshot.js:47,94,368,388,973,1021-1023,1366`
  - `gbdraw/web/js/services/config.js:3715,3883,3913,4161`
  - `gbdraw/web/js/services/session-authority.js:80`
  - `gbdraw/web/js/services/reset.js:111`
  - `label-actions.js:385,992`
  - 旧 Session はそのまま読める。reader の `copyFields`（`session-authority.js:497`）が未知の field を落とす。
- 再利用するもの:
  - `clearOverrides`、`hasFeatureScopedOverrideInSvg`、`sourceReplaced`。
  - Generate の pre-state handle。label の map は `history-snapshot.js:1011-1020` に含まれるので、Undo と rollback で復元できる。
- 追加しないもの: fingerprint、新しい state や Session field、watcher、Result ごとの override map。

**4. 上位の対応**
- web CLAUDE.md に invariant を 1 文加える（横断 R1）。
- `docs/SESSION_COMPATIBILITY.md` に次の 1 行を加える。「current writer は `features.labelOverrideContextKey` を出力しない。reader は無視する。」
- 任意の PD（source 置換時の細部）:
  - **A（推奨）**: 対象単位で prune する。`reconcileFeatureVisibilityOverrides` と同じ述語を使い、bulk override は matcher として残す（05A8 の visibility rules と同じ扱い）。
  - **B**: 現在の all-or-nothing を維持する。PR-1 の結果はこれになる。
  - **C**: 明示的な Reset まで全部残す。

**5. テスト**
- Node（`tests/web/feature-label-visual-unit.test.mjs`）:
  - feature 集合が異なる SVG に順に `syncLabelEditor` しても、3 つの map が変わらない。監査の不具合はこれで検出できた。
  - `sourceReplaced=true` で target が無ければ全部消える。
  - `sourceReplaced=true` で target が有れば残る。
- Node（`run-analysis-simple-path.test.mjs`）: label の map しか無いときでも、source 置換で `sourceReplaced` が `bindingOptions` に載る。
- Playwright（`tests/web/gui-audit-regressions.playwright.spec.js`）: 監査 spec の batch-label-loss と label-context を assert 化する。
  - Result の往復後も `labelTextFeatureOverrides` が `{[id]:'EDITED_B_TRNA'}` のまま。
  - Generate 後の SVG の label が `EDITED_B_TRNA`。
  - 非表示 → Generate → 表示 → Generate で `AUDIT_X` が戻る。
  - Session 保存に override が含まれる。
- 既存テストの更新: `history.test.mjs`、`session-draft-authority.test.mjs:731`、`helpers/mode-transition.cjs:91`、`feature-label-visual-unit.test.mjs:72`。

**6. 規模と PR**
- 本番ファイル 8（label-actions、run-analysis、app-setup、state、history-snapshot、config、session-authority、reset）。Ordinary の上限ちょうど。net はマイナス。
- Review REQUIRED（compatibility-path）。FE-01 だけの PR にする。

**7. 依存とリスク**
- 他に依存しない。
- 使われない override が label table に `*\t*\thash\t^f…$` の行として残る。Python では一致しないので出力は変わらず、Session が少し大きくなる。
- batch で「Reset all label text overrides」を押すと、全 Result の override が消える。明示操作なので許容できるが、文言は確認する。
- batch の Label TSV 取り込みは、mounted Result の label だけで評価したうえで `clearOverrides` を呼ぶ（`label-actions.js:1036-1077`、コードからの判断）。FE-02 と同じ系統の残課題として記録しておく。

---

## FE-02（P2）と PV-09: batch で範囲指定の色・非表示・凡例の編集が表示中の Result にしか反映されない

**1. 分類: PRODUCT_DECISION_REQUIRED**
- authority の検索結果:
  - web CLAUDE.md の live-edit 条項は「mounted target と current Result を同期して更新する」としか書いておらず、batch で表示していない Result の扱いを定めていない。
  - PD-OI-037 と OIC-024 は、操作ごとの Live edit の分類が事実どおりであることを求める。ところが scope dialog は他の Result を含む件数を示している（3 件、`out-batch-scope.json`）。現状は事実どおりではない。
  - PD-OI-052 の「batch全出力を対応identityへ」は Generate の装飾 delta についての決定で、live edit には及ばない。
- 結果の候補:
  - **A: 編集時に伝播する。** target を含む、表示していない Result をすべて parse → 投影 → serialize する。Session と export は常に一致する。
    - ただし label の DOM identity は mounted binder でしか作れない（`gbdraw/web/js/services/svg-result-ingestion.js:511`）。未 mount の Result の label は遅れて反映されることになり、方式が混在する。
    - 編集ごとに O(Result 数) の parse と serialize がかかる。
  - **B（推奨）: 表示時に正本から再投影する。** Result を mount するときに、正本の intent を投影してから bind する。対象は palette と rules、visibility の override と editor qualifier rule、stroke、凡例の削除・改名・色、label。
    - 観察できる出力（preview と、その Result の export）は常に一致する。
    - 残るリスク: 一度も mount していない Result は、Session に保存される SVG の bytes が次の mount か Generate まで古い。Load 後にその Result を選べば再投影される。
  - **C: 現状を維持して明示する。** 即時に反映するのは表示中の Result だけで、他は Generate で反映される、と明示する。dialog の件数も mounted Result に限る。
    - LSP に反する（single と batch で同じ操作が別の結果になる）。表示していない Result の export は古いまま。
- B を推奨する理由:
  1. label は構造上 B しか取れないので、全 domain を B に揃えれば方式が 1 つで済む。
  2. History の after-apply hook（`app-setup.js:2643-2667`）が、正本から mounted Result への投影をすでにひととおり持っている。mount のときにも同じ関数を呼べば、2 つの経路が 1 つにまとまる。
  3. 編集ごとの N 件の parse が無い。
  4. export は mounted Result 単位（`app-setup.js:3621-3646`）なので、出力の正しさは B で保証できる。

**2. 根本原因の再確認**
- 監査の指摘どおり。投影の仕組みが 2 つある。
  - E1（Generate と Session の取り込み、全 Result 対象）: `gbdraw/web/js/app/candidate-render.js:103-315` → `svg-result-ingestion.js:360-470`。
  - E2（mounted DOM のみ）: `gbdraw/web/js/app/svg-styles.js:370-411`、`gbdraw/web/js/app/feature-editor/visibility-actions.js:210-222`、`gbdraw/web/js/app/legend/entry-actions.js:543-660`、`label-actions.js:541-618`。
  - Generate と Generate の間に表示していない Result を投影する仕組みは無い。binder の同期（`app-setup.js:2237-2307`）は label だけを扱う（`:2284-2291`）。
- PV-09 で Result 2 に切り替えると Legend editor に `GC content` が再び出る（editorOn2）のは、`adoptLegend`（`:2238-2250`）が mounted SVG から `legendEntries` を作り直すため。
- 次の Generate で削除が全 Result に効くのは、`candidate-render.js:292-296`（allResultIndexes）が diagram 全体に適用する正しい意味論で、バグではない。
- PV-09 の「rename lost」は PV-02（`candidate-render.js:228-242`）の重複で、batch 固有ではない。PV-02 は別のワークストリームが担当する。

**3. 修正案（B を採用した場合。決定前は実装しない）**
- owner は `app-setup.js` の binder の構成。
- History の after-apply の処理本体を `reconcileMountedEditorIntent({domains})` として 1 つにまとめる。`setAfterApplyHistoryIntent` と、binder の `phase === 'result-selection'` の両方から呼ぶ。label だけを特別扱いしている同期はこの共通呼び出しに吸収する。
- Legend: result-selection のときは、inventory を抽出する前に正本の凡例 intent を適用する。
  - 新しい投影は作らない。E1 の凡例操作（`compileDirectEditorMutationPlan` の legendDeletes/Renames/Fills/Strokes と `applyLegendOperations`）を、mounted root に対して idempotent（allowMissing）で再利用する。
  - このために `svg-result-ingestion.js` から export を 1 つ出す。
- FE-04 の修正後の `reconcileFeatureVisibility` を使う。
- 全 intent が空なら何もしない（CW-03 の zero fast path）。`recordStructuralMetric` で計測する。
- 削除するもの: binder の label 専用の同期条件。History hook 内に並んでいる個別の呼び出し（共通関数へ移す）。
- 追加しないもの: Result ごとの override の複製、編集時の全 Result の parse、新しい schema、Worker、Session field。

**4. 上位の対応**
- 静的な Product Contract に PD-OI-056（concern の例: `web.batch-live-edit-projection`）と OIC-027 を加える。co-change route でよい。
  - OIC-027 の内容: batch のすべての Result が、観察時に live edit を反映する。export、Save → Load → 選択、Undo/Redo、zero fast path、stale と cancel のときの旧 Result 保持を含む。
- web CLAUDE.md の live-edit 条項に「表示していない Result は、mount のときに同じ reconcile で投影する」を加える。
- Architecture ratchet 向けに、2 経路（History apply と Result mount）をまとめたことの owner/path の根拠を記録する。

**5. テスト**
- Playwright（監査 spec の batch-scope と batch-legend を assert 化）:
  - Result 1 で、label 範囲の色、Exact product での非表示、凡例の削除を行う。
  - Result 2 を選んだ後、mounted SVG、`results[1]`、`downloadSVG` の出力のすべてで次を確認する。TESTB_0001/0002 が `#ff0000`、TESTB_0004 が hidden、`GC content` が無い。`legendEntries` にも `GC content` が無い。
  - Save → Load → Result 2 の選択でも同じになる。
  - Undo で両方の Result が元に戻る。
- Node: History apply と result-selection が同じ関数を呼ぶこと。override が空なら投影の metric が 0 であること。

**6. 規模と PR**
- 5–7 ファイル（app-setup、legend/entry-actions、svg-result-ingestion、candidate-render。必要なら preview-runtime と svg-styles）。
- Architecture profile。Review REQUIRED（責務の移動、新 export、Product Contract の co-change）。
- FE-02 と PV-09 は同じ PR にする。

**7. 依存とリスク**
- PD が決まるまで runtime には着手しない。
- 先に入れておくもの:
  - FE-04: reconcile を正しくするため。
  - FE-03: Result 2 で編集できる状態で検証するため。
  - FE-01: label が残らないと、B の label 部分が成り立たない。
- 凡例の owner は PV-02 と PV-03 と重なるので、PV のワークストリームと順序を合わせる。PV-02 を先に直すと、E1 の凡例操作を再利用しやすくなる。
- リスク:
  - 大きな batch では、mount のたびの投影コストがかかる。
  - `rulePreparation.prepare()` が非同期なので、mount が遅れることがある。
  - 凡例の rename が idempotent であること。

---

## FE-03（P2）batch で Features drawer が record 1 を並べ、Edit を押しても何も起きない

**1. 分類: IMPLEMENT_EXISTING_AUTHORITY**
- 根拠（web CLAUDE.md）:
  - 「one canonical selected value / no mirrors」。
  - 「same availability predicate for enabled rendering and action resolution… must not silently return」。
- 先例: Preview 検索は mounted SVG に描画された feature に限っている（`gbdraw/web/js/app/feature-search/preview-actions.js:161-166`）。drawer の文言も「in current diagram」（`gbdraw/web/js/app/feature-editor/svg-actions.js:400`）。
- 「drawer の Edit で該当する Result に切り替える」という案は新しい導線になり、この修正には必要ない。望むなら別に PD を立てる。

**2. 根本原因の再確認**
- 監査の指摘どおり。batch で「どの record を表示しているか」を決める正本は `selectedResultIndex` で、`selectedFeatureRecordIdx` はその写しになっている。
  - Generate の後、record 側は 0 に戻る（`run-analysis.js:4929`）が、Result の選択は保たれる。
- 一覧（`state.js:712-735`）と Edit（`svg-actions.js:404-413`、mounted SVG に要素があるか）は別の述語を使っている。

**3. 修正案（コード）**
- owner は `state.js` の `filteredFeatures`。
- 一覧の母集合を、mounted Result に描画された feature に限る。
  - 既存の Result ごとの metadata を使う: `getCommittedSvgResultMetadata(results.value[selectedResultIndex.value])?.renderedFeatureIdentities?.renderedIds`（`svg-result-ingestion.js:177-179`）。
  - この metadata は symbol property なので、spread で copy しても保たれる。
  - metadata が無い場合（legacy）は、従来どおり全件を並べる。
- record picker は「mounted Result に record が 2 つ以上あるとき」だけ表示し、適用する。computed を 1 つ作り、`filteredFeatures` と `gbdraw/web/index.html:6590` の v-if で共用する（`featureRecordIds.length > 1` を置き換える）。
  - single、grid、Linear は 1 つの Result に全 record があるので挙動は変わらない。
- `app-setup.js:1373-1376` の既存 watcher の source に `selectedResultIndex` を加える。新しい watcher は作らない。
- 追加しないもの: record と Result の選択を同期する watcher、`result_index` field、`svg-actions.js` の新しいメッセージ（述語が揃えば到達しなくなる）。
- `svg-result-ingestion.js` から `state.js` への import cycle は無いことを確認した。

**4. 上位の対応**: 不要。

**5. テスト**
- Playwright（監査の batch-drawer を assert 化）:
  - Result 2 のとき、行がすべて `TESTB:*` で picker が出ず、Edit で TESTB の `clickedFeature` が開く。
  - grid では picker が従来どおり動く。
- Node: `session-draft-authority.test.mjs:44` と同じ形で `state.js` を import し、2 つの Result を legacy admission と fake-svg-dom で用意できれば、`filteredFeatures` を検証する。

**6. 規模と PR**
- 3 ファイル（state.js、index.html、app-setup.js）。Ordinary。
- 単独の PR にし、FE-02 より前に出す。

**7. 依存とリスク**
- 他に依存しない。
- mounted Result の record 数は `record_key` で数え、Circular と Linear で共通にする。

---

## FE-04（P2）Redo や無関係な Undo のあと、Exact product/protein ID で隠した feature が再び表示される

**1. 分類: IMPLEMENT_EXISTING_AUTHORITY**
- 根拠: web CLAUDE.md（History は owning action と同じ操作を記録する）と OIPC-C04。
- Redo で rule は戻るが、`results[selected]` は visible のまま flush されてしまう。正しい結果は「action の直後と同じ表示」の 1 つしかない。

**2. 根本原因の再確認**
- 監査の指摘どおり。reconcile（`visibility-actions.js:491-500`）は override だけを見る `getFeatureVisibility`（`:434-438`）を使う。一方、live action は scope の一致判定（`getMatchingQualifierFeatures` `:123-131`）で決める。同じ実効 visibility を、2 つの別々の評価器が出している。
- `resolveEffectiveFeatureVisibility`（`gbdraw/web/js/app/feature-visibility.js:412-437`）は hash rule しか評価しない。
- 発火するのは History apply の `features` domain（`app-setup.js:2657-2660`）だけ。label 編集の Undo などが当たる。config の Undo では起きない。

**3. 修正案（コード）**
- owner は `feature-visibility.js`。
- `resolveEffectiveFeatureVisibility` が feature を受け取れるようにし、次の順に評価する。
  1. override
  2. manualRules を順に first-match。既存の hash rule に加え、`isEditorExactQualifierRule`（`:451-461`）に当たる rule を評価する。type と qualifier の全値に対して `new RegExp(value,'i')` を使い、Python（`gbdraw/features/visibility.py:216` の IGNORECASE search）と揃える。
  3. どれにも当たらなければ `'on'`。
- matcher を 1 つ export し、`getMatchingQualifierFeatures` をそれで置き換える。action と reconcile が同じ述語を使うようにする。
- `reconcileFeatureVisibility` と `applyVisibilityPreviewForScope` の default 分岐は、feature を渡して resolver を呼ぶ。
- 追加しないもの: 表から作った rule の live 評価（Applies on Generate のまま）、History の変更。

**4. 上位の対応**: 不要（横断 R2 で予防する）。

**5. テスト**
- Node（`feature-visibility.test.mjs`）: editor の product rule で `'off'` になる。大文字小文字を区別しない。override が優先される。editor 以外の rule は無視する。
- Node（`feature-visibility-actions.test.mjs`）: `handleFeatureVisibilityScopeChoice('product')` の後に `reconcileFeatureVisibility()` を呼ぶと、`mode:'off'` の change が出る。現状は `'default'` になるので、これで検出できた。
- Playwright（監査の visibility-redo を assert 化）:
  - Exact product で非表示 → Undo → Redo のあと、mounted SVG と Result の ND1 が `display=none`。
  - 別の label 編集を Undo した後も hidden のまま。

**6. 規模と PR**
- 2 ファイル、gross 40–60。Ordinary。FE-10 と同梱してよい。

**7. 依存とリスク**
- FE-02 の B の前提になる。
- 値を複数持つ qualifier では、live の対象が Python 側に揃う方向に変わる。改善だが、テストで明示しておく。

---

## FE-10（P3）キャンセルした Reset fill ダイアログの既定色が、次の Reset で使われる

**1. 分類: IMPLEMENT_EXISTING_AUTHORITY**
- Reset fill の正解は、その type の palette 既定色（`color-actions.js:1284`）の 1 つだけ。
- web CLAUDE.md の「no mirrors」と「Close は可視性だけを変える」に当たる。

**2. 根本原因の再確認**
- 監査より範囲が広い。キャンセルしたときだけでなく、dialog を一度でも完了した後に、sibling の無い feature を Reset すると古い色が使われる。この Reset は dialog を経ずに `doResetFillColor('this')` を直接呼ぶ（`:1302`）。
- `resetColorDialog.defaultColor` は閉じても消えない写し。`:1299` で設定され、`:1309` とテンプレートの Cancel では消えず、`:1317` で優先して使われる。

**3. 修正案（コード）**
- `:1317` は `appliedPaletteColors.value[feature.type]` だけを使う。
- `:1299` と、`state.js:526-531` の `defaultColor` を削除する。
- `index.html:6366` の Cancel は `handleResetColorChoice('cancel')`（owner の遷移）を呼ぶ。

**4. 上位の対応**: 不要。

**5. テスト**
- Node（`feature-color-actions.test.mjs`）: sibling のある tRNA で dialog を開き、cancel する（'this' で完了する場合も同じ）。その後 sibling の無い rRNA を Reset すると、確定した rule の色が `palette.rRNA` になる。
- `linear-multi-record.playwright.spec.js:2665` の `defaultColor` への代入を削除する。

**6. 規模と PR**
- 3 ファイル、約 8 行。Ordinary で、Review は CLEAR の見込み。

**7. 依存とリスク**: なし。

---

## 横断的な提案

- **共通する原因:** diagram 全体のエディタ intent（安定した identity、caption、matcher で key する）と、1 つしか mount されない Result の view 状態が混同されている。現れ方は次の 5 つ。
  - 表示の変化を代わりの指標にして intent を消す（FE-01）。
  - mounted Result にしか投影しない（FE-02、PV-09）。
  - 一覧の母集合が Result の選択から導出されない（FE-03）。
  - reconcile が action と別の評価器を使う（FE-04）。
  - dialog が導出値を持ち続ける（FE-10）。
- **R1（intent の寿命）:** 次の文を web CLAUDE.md の live-edit invariants に加える（OIPC-C05/C06 の Web 版）。「editor override は、view の変化（Result 選択、mount、record 選択、mode、非表示、reflow）では作成・剪定・削除されない。削除は明示操作（Reset、Import、Undo/Redo、Session 置換）と、成功した source 置換 Generate の中での owner reconcile に限る。」
  - 強制手段 1: Node contract。feature 集合が重ならない SVG を順に bind しても、override の map が変わらないことを確認する。
  - 強制手段 2: `architecture-contracts.test.mjs` に静的 assertion を加える。override の map を delete/clear してよいのは、許可リストの owner のパスだけにする。
- **R2（domain ごとに投影は 1 つ）:** live action、History apply、Result mount は同じ投影関数を呼ぶ。
  - 強制手段: 各 action の直後に reconcile を呼んでも DOM が変わらないことを確かめる parity test。FE-04 型の退行を検出できる。
- **R3:** drawer、検索、editor の対象は、mounted Result の committed metadata から導出する。同じ概念の選択 ref を 2 つ持たない。
- **R4:** dialog の reactive object は表示用の値だけを持つ。確定に使う値は、確定の時点で正本から導出する。
- **LSP:** 2 record の batch fixture（監査の `multi.gb` 相当）を `tests/web/helpers` に置き、editor 系の回帰テストを single と batch の両方で実行する。
- **CW-05:** consumer を消したら、Session writer の field も消す（`labelOverrideContextKey` がその例）。
- FE-02 の PD を待っている間に、FE-01、FE-03、FE-04、FE-10 を先に出せる。
- 他ワークストリームとの境界:
  - 凡例の owner は PV のワークストリームと共有している（PV-02、PV-03）。
  - FE-04 の修正は History 側（SE のワークストリーム）に手を入れない。
  - X-01 と X-02 には依存しない。

## 監査の分類で見直すべき点、バグではない点

- **PV-09**: 「次の Generate で削除が全 Result に適用される」のは diagram 全体の意味論として正しい。問題は Generate までの乖離だけ。「rename lost」は PV-02 の重複。
- **FE-10**: キャンセルしたときに限らず、dialog を一度使った後なら常に起きる。
- **FE-01**: 監査に書かれていない影響が 2 つある。bulk override と sources の map も消えること、reflow の経路もあること。record selector の経路はコードから判断したもの。
- **FE-04**: 「無関係な Undo」は、`features` domain の Undo に限られる。
- **FE-03**: 監査の見立てどおり。picker で record を選び直す回避策があるので、P2 は妥当。
