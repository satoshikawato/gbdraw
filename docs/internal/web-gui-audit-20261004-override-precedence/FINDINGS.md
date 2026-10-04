# Web GUI 監査: 個別 override と全体設定の食い違い（2026-10-04）

対象は PR #757（`fix/label-override-respects-label-scope`、head `c071e477`）を `dev`（`e73df62b`）に載せたコード。
#757 は、Linear の Show Labels = First Record Only で popup の **Label visibility: On** が効かず、次の Generate が
`UNKNOWN` で失敗した不具合を直した。この監査は、同じ種類の不具合が GUI の他の場所に残っていないかを調べた。

## 方法と表記

- クラスは依頼のとおり: **A** 個別 override と全体設定の優先順位、**B** Web の ID と Python の selector の食い違い、
  **C** 何も起きない live edit、**D** 生成後の検証が診断情報なしのエラーになる、**E** 個別編集のあとに全体設定を変える、
  **F** label reflow の `reason` の受け渡し。
- すべての不具合は `c071e477` の Web app を Chromium（Playwright）で動かして再現した。B の一部は、保存した Session を
  CLI（`gbdraw ... --session`）で描き直して Python が行を照合するかを確かめた。`UNKNOWN` の生の例外は、
  `failOperation`（`app/run-analysis.js`）に一時的に入れた `console.error` で取った（commit していない）。
- 再現 probe とログはリポジトリに入れていない（`tests/web/_audit_probe/`、手元の写しは
  `/home/kawato/gbdraw-baselines/override-precedence-audit-20261004/probes/`）。修正 PR が失敗するテストを足す。
- 重大度: **図の誤り** > **Generate の失敗** > **無言の no-op** > **文言の不一致**。
- 判断の基準: Owner 方針（2026-10-04）「popup で明示した個別設定は全体設定に勝つ。個別編集の副作用で全体設定を書き換えない。
  編集が見えないままになるときは具体的な選択肢のダイアログで確認する」、OIPC-C03（受け付けた値は必ず使われるか、
  実行前に拒否される）、PD-OI-069（Web が作る色ルールは Python が照合できる値にする。label の instance 単位の編集は維持する）、
  R1–R12（`gbdraw/web/CLAUDE.md`）。
- JS のパスは `gbdraw/web/js/` からの相対パス。

## 要約

不具合は 13 件（図の誤り 4、Generate の失敗 5、無言の no-op 2、文言・状態の不一致 1、検査の抜け 1）。
このほか、#757 自身のテスト helper の取り残し 1 件をこのセッションで直した（OV-14）。

| ID | クラス | 重大度 | 症状 | 根本原因 |
|---|---|---|---|---|
| OV-01 | B | 図の誤り | crop / 逆相補した record で、Feature visibility と label 文字の編集が Python に届かないか、別の feature に当たる | RC-1 |
| OV-02 | B, C | 図の誤り | crop した record の `record_location` 色ルールが、live と Generate で別の feature を塗る | RC-1 |
| OV-03 | B | 図の誤り | 選択から作った annotation が crop した record で別の feature に付くか、黙って捨てられる | RC-1 |
| OV-04 | B | 図の誤り | 同じ record を 2 回読み込むと、片方だけへの label・visibility の編集が Generate で両方に当たる | RC-2 |
| OV-05 | B, D | Generate の失敗 | 描画 ID に `__instance_` が付く feature（同じ座標の重複、重複 record の Circular canvas）の編集が Python に届かない。Label On は `UNKNOWN` | RC-2, RC-3 |
| OV-06 | A, D | Generate の失敗 | underlay で描く feature に Label On → live は何も起きず、Generate は `UNKNOWN`。以後の Generate と Session 読み込み後も失敗し続ける | RC-3 |
| OV-07 | A, D | Generate の失敗 | Label Rendering = Embedded Only で収まらない label に On → OV-06 と同じ | RC-3 |
| OV-08 | E, D | Generate の失敗 | Circular で Feature placement を lane にしたまま Linear に切り替えると、Linear の Generate がすべて `UNKNOWN` | RC-4 |
| OV-09 | E, D | Generate の失敗 | lane の Feature placement のあとで Track type や Separate strands を変えると、`RENDER_FAILED`（ValidationError）で、原因の feature も次の行動も示さない | RC-4 |
| OV-10 | A, C | 無言の no-op | 隠した feature に Label On → 何も起きない。popup は "Choose On to show it." と案内し続ける | RC-3 |
| OV-11 | E, C | 無言の no-op | 個別編集のあとで crop や逆相補を変えると、その編集が黙って効かなくなる | RC-5 |
| OV-12 | B | 状態の不一致 | Linear の record の並べ替えのあと、label を隠した feature の popup が Label visibility を Default と表示する | RC-5 |
| OV-13 | D | 検査の抜け | `tests/web/error-producer-coverage.test.mjs` が OV-05〜07 の throw を含む 4 つの owner を数えていない | RC-6 |

## 根本原因

- **RC-1: Web が作る selector 値の座標系が、Python の照合と違う。** Web の feature catalog の `selector.hash`、
  `location`、`record_location` は元の座標（crop・逆相補の前）で作られる（`gbdraw/web_support/feature_metadata.py::_biological_selector_values`）。
  Python の visibility、label、色、annotation の照合は、描画した record（crop・逆相補の後）の値と比べる
  （`gbdraw/features/selector_values.py`）。`hash` 行については、描画した feature の hash を使う契約が既にある
  （PD-OI-069、`app/rule-matching.js::ruleFeaturePayload` のコメント）。この契約に従っていない builder:
  `app/feature-selector.js::normalizeFeatureSelectorMetadata`（`stableFeatureId`、`recordLocation`、`position`）を通る
  `app/feature-visibility.js::buildFeatureVisibilitySelectorCache`、`app/feature-editor/label-override-table.js::selectStableFeatureKey`
  （`record_location`）、`app/rule-matching.js::ruleFeaturePayload`（`location`、`record_location`）、
  `app/annotations/target-actions.js::featureTargetsFromSelection`（`selector.hash`）。
- **RC-2: 描画 ID の接尾辞と複製の区別。** `app/feature-utils.js::getFeatureHashCandidates` は `_record_<n>` だけを外し、
  `__instance_<…>` を外さない（Circular は record ID の重複や表示開始位置の指定で、両 mode は同じ座標の重複で付ける。
  `gbdraw/render/interactive_context.py`、`gbdraw/labels/circular.py::instance_svg_id`）。また Python の label・visibility の行は
  `record_id` と hash で照合するので、同じ record の複製や同じ座標の重複のうち 1 つだけを指定できない。
- **RC-3: Python が描けない Label On を、Web が必須にして素の Error で失敗させる。** Python は underlay の feature
  （`gbdraw/features/factory.py`、`include_label = ... and rendering != "underlay"`）、Embedded Only で収まらない label
  （`gbdraw/labels/circular.py`、`gbdraw/labels/linear.py` の `embedded_only` 分岐）、隠した feature に label を描かない。
  `app/run-analysis.js` は Label visibility On の feature をすべて必須の binding にし、
  `app/feature-editor/label-actions.js::requireUniqueEditableLabelBindings` が素の `Error` を投げ、
  `normalizeUserFacingError` が `UNKNOWN` にする（R6 違反）。live の reflow はエラーも出さずに何も描かない。
- **RC-4: Feature placement の override が mode と slot に追従しない。** `services/session-request.js` は
  `canonicalFeaturePlacements(state.featurePlacementOverrides, mode)` を全行に適用する。record key は mode ごとに違い
  （Linear は `seq.uid`、Circular は `circularRecordKey`）、`services/feature-placement.js::canonicalFeaturePlacements` は
  他 mode の side を素の `Error` で拒否する。slot の向きが変わっても override は残り、Python の
  `gbdraw/features/placement.py::FeaturePlacementSlot.validate_target` が `diagnostic=` のない `ValidationError` を投げる。
- **RC-5: 個別 override の key が描画 ID に依存する。** visibility と label の override は描画 ID（`<hash>_record_<n>`）、
  色は描画 hash を key にする。crop・逆相補・並べ替えで描画 ID が変わると、override は残ったまま当たらなくなる。
  `app/feature-visibility.js::pruneUnmatchedFeatureOverrides` は元の座標の hash で生存を判定するので、消しもしない。
  Feature placement だけは (`record_key`, `biological_feature_id`) を key にし、Python が元の catalog から解決するので、
  crop・逆相補・複製のどれでも正しく当たる（`gbdraw/features/placement.py::resolve_placement_inputs`）。
- **RC-6: 検査の抜け。** `tests/web/error-producer-coverage.test.mjs` は `app/feature-editor/label-actions.js`
  （未分類 1）、`services/svg-result-ingestion.js`（22）、`app/candidate-render.js`（3）、`app/preview-runtime.js`（7）を数えない。

## 所見の詳細

### OV-01: crop・逆相補した record で visibility と label 文字の編集が Python に届かないか、別の feature に当たる

- 再現 A: Linear、`tests/fixtures/web_batch_two_records.gb` の TESTA を 201..3800 で crop、TESTB を逆相補、misc_feature を表示。
  TESTA の misc_feature と TESTB の tRNA を Feature visibility Off、TESTA の tRNA の label 文字を `PROBE_A_TRNA` にする。
  保存した Session の行は `TESTA misc_feature hash ^fb5977f81$ off`（元座標の hash。描画 hash は `f3b928d8c`）、
  `TESTA tRNA record_location ^TESTA:2549\.\.2620:-$ PROBE_A_TRNA`（元座標）。CLI で Session を描き直すと、
  隠したはずの 2 feature が描かれ、`PROBE_A_TRNA` は出ない。
- 再現 B: 同じ crop で、misc_feature 2601..2700（Y）と tRNA complement(2750..2820)（TY）を足した record。
  X（misc_feature 2401..2500）を Off、TX（tRNA complement(2550..2620)）の label を `EDITED_X` にして Generate すると、
  Web の Result から **Y が消え、TY に `EDITED_X` が付く**（Y の描画 hash が X の元座標 hash と一致する）。
- 期待: 編集した feature だけが変わる。Web の Result と、Session を CLI で描き直した図が一致する（OIPC-C04）。
- 補足: Web は Generate 後に live の `display="none"` と label 文字を Result に重ねるので、Python が行を照合しなくても
  Web の表示は正しく見えることがある。そのため Web 上では気づきにくい。
- 修正案: 根本原因 RC-1 の builder が、Python と同じ描画座標の値（`hash` は `getFeatureGenerationHash`）を使う。
  分類は IMPLEMENT_EXISTING_AUTHORITY（PD-OI-069 の契約と OIPC-C03）。ただし OV-11 の扱い（Q4）で key の設計が変わるので、
  Q4 の回答を待って実装する。

### OV-02: crop した record の `record_location` 色ルールが live と Generate で別の feature を塗る

- 再現: OV-01 再現 B の record に、色ルール `misc_feature record_location ^TESTA:2400\.\.2500:\+$` を足す。
  live は X を塗り、Generate は Y を塗って X は既定色に戻る。Label TSV の読み込み
  （`app/feature-editor/label-actions.js::loadLabelOverrideTable`）も同じ payload を使う。
- 原因: `app/rule-matching.js::ruleFeaturePayload` が `hash` だけ描画座標で、`location` と `record_location` は元座標で送る。
- 修正案: OV-01 と同じ（RC-1）。

### OV-03: 選択から作った annotation が crop した record で別の feature に付くか、捨てられる

- 再現: OV-01 再現 B で X を選択し、選択から annotation を作る。対象は `hash=fb5977f81`（元座標）。Generate 後の annotation は
  Y の位置に付く。重なる feature がなければ Python は `feature_selector_unmatched` で annotation を捨てる。
- 原因: `app/annotations/target-actions.js::featureTargetsFromSelection` が `selector.hash` を使い、
  `gbdraw/annotations/resolve.py::_feature_matches` は描画 hash と比べる。
- 修正案: OV-01 と同じ（RC-1）。

### OV-04: 同じ record を 2 回読み込むと、片方だけの編集が Generate で両方に当たる

- 再現（Linear、TESTA を 2 回）: 1 つ目の複製の feature に Feature visibility Off、Label visibility Off、label 文字の変更を
  それぞれ行う。live は 1 つ目だけが変わり、Generate で両方が変わる。Circular の record ID が重複する canvas でも、
  visibility Off が両方の複製から feature を消す。
- 期待: PD-OI-069 は「label の instance 単位の編集」を維持すべきものとしている。visibility は決定がない。
  「This feature only」の色が両方に付くのは PD-OI-069 が受け入れた残余リスクで、不具合ではない。
- 原因: RC-2。Python の行は `record_id` と hash で照合し、複製はどちらも同じ。live の投影は描画 ID ごと。
- 修正案: Q4（個別 override の identity）。

### OV-05: `__instance_` が付く feature の編集が Python に届かず、Label On は `UNKNOWN`

- 再現 1（同じ座標の重複）: 同じ type・同じ座標（3001..3600 +）の CDS 2 つを含む record。Show Labels None で Generate し、
  片方に Label On を Apply する。描画 ID は `f1aff4c2b__instance_4_…` と `f1aff4c2b__instance_5_…`、送る行は
  `hash ^f1aff4c2b__instance_5_…$` で、Python の hash `f1aff4c2b` と一致しない。live は何も描かず、Generate は `UNKNOWN`。
- 再現 2（Circular、同じ record を 2 回並べた multi-record canvas）: 描画 ID は `ffa1f4c4a__instance_record_1_…`。
  「This feature only」の色は live では付き、Generate で既定色に戻る。label 文字・表示の行も照合されない
  （Session を CLI で描き直すと `COPY1_ONLY` は 0 回）。
- 原因: RC-2。label 文字だけの編集（Keep hidden を選んだもの）も `hash` 行になり、Python は空でない `hash` 行を
  On と扱う（`gbdraw/labels/filtering.py`、docs の Feature presentation 節）。行が照合されるようになると、Keep hidden を選んでも
  label が出る（今は行が照合されないので出ない）。
- 修正案: 色は `__instance_` も外した描画 hash を使う（PD-OI-069 の A / STABLE-HASH-ONLY。複製は色を共有する）。
  label を 1 つの instance に当てるには Python の selector の拡張が要る（Q4）。文字だけの行は On にしないこと。

### OV-06: underlay の feature に Label On → live は何も起きず、Generate は `UNKNOWN` のまま続く

- 再現: `repeat_region`（既定で underlay）を含む record、Show Labels Out / All Records。popup で Label On、文字 `RPT_FORCED`、
  Apply Label。popup は "This feature has no label in the current Result. Choose On to show it." と案内する。reflow 後も
  label は無く、`labelReflowLastError` は null。Generate は `Code: UNKNOWN / Stage: render`。生の例外は
  `Sanitized SVG content is missing or ambiguously binds an editable Label.`（`requireUniqueEditableLabelBindings`）。
  popup を Default に戻すまで、以後の Generate も、Session を保存して読み直したあとの Generate も失敗する。
- 原因: RC-3。
- 修正案: Q1。どの選択でも、必須 binding の失敗は分類したエラーにする（R6、PR-2）。

### OV-07: Label Rendering = Embedded Only で収まらない label に On → OV-06 と同じ

- 再現: `adv.label_rendering = 'embedded_only'` で Generate し、90 bp の CDS（長い product 名）に Label On。
  reflow は何も描かず、Generate は `UNKNOWN`。
- 原因: RC-3（`embedded_only` の分岐が強制 label も捨てる）。
- 修正案: Q1。

### OV-08: Circular の lane placement が Linear の Generate をすべて失敗させる

- 再現: Gallery Session `HmmtDNA_basic_circular`、Generate。CDS の popup で Feature placement を Outward lane 1 にして Generate
  （成功）。Linear に切り替え、`tests/test_inputs/HmmtDNA.gbk` を入れて Generate すると
  `Code: UNKNOWN / Operation: generate / Stage: request-validation`。Linear の入力が別のファイルでも同じ。
- 期待: R2 により mode の切り替えは override を消さない。Circular の placement は Circular の record にだけ当たり、
  Linear の Generate は影響を受けない。
- 原因: RC-4（他 mode の行も request に載せ、素の `Error` で拒否する）。
- 修正案: request には現在の mode の record に属する行だけを載せる（PR-1、IMPLEMENT_EXISTING_AUTHORITY: R2 と mode ごとの record key）。

### OV-09: lane placement のあとで slot の向きを変えると、Generate が原因を示さずに失敗する

- 再現 1（Circular）: `HmmtDNA_basic_circular`（Track type middle）で CDS を Outward lane 1 にして Generate（成功）。
  Track type を spreadout または tuckin にして Generate すると `RENDER_FAILED`（"Python exception: ValidationError"）。
- 再現 2（Linear）: `lambda_basic_linear` で Separate strands を外して Generate、CDS を Above lane 1 にして Generate（成功）。
  Separate strands をオンにして Generate すると同じ `RENDER_FAILED`。popup の Feature placement は `above` のまま、
  その選択肢は "Unavailable in the current draft feature slot." で無効と表示される。
- 期待: 失敗するなら、原因の feature と次の行動（placement を Auto か Main に戻す、または slot を戻す）を示す（R6）。
  失敗させずに扱う方法は Q3。
- 原因: RC-4（slot の変更で override を照合しない。Python の例外に `diagnostic=` がない）。
- 修正案: Python が分類した診断を出す（PR-3、IMPLEMENT_EXISTING_AUTHORITY: R6）。ダイアログか自動の扱いかは Q3。

### OV-10: 隠した feature に Label On → 何も起きず、popup の案内も誤り

- 再現: Feature visibility Off の feature に Label On（順序を逆にしても同じ）。label は出ず、ダイアログもなく、Generate は成功する。
  popup は "Choose On to show it." と案内し続ける。feature type を Features から外した場合も同じで、override は残る。
- 原因: RC-3（隠した feature には label を描かない。Web はそれを知らせない）。
- 修正案: Q2。

### OV-11: 個別編集のあとで crop や逆相補を変えると、編集が黙って効かなくなる

- 再現: OV-01 再現 A の 4 つの編集のあと、TESTA の crop 開始を 201 → 101、TESTB の逆相補を外して Generate。
  4 つとも図から消える（feature は表示、label は既定、色は既定）。override は古い描画 ID の key のまま残り、
  `featureColorOverrides` には色が残っているのに図は既定色。`pruneUnmatchedFeatureOverrides` は元座標の hash が残っているので消さない。
- 原因: RC-5。
- 修正案: Q4。

### OV-12: record の並べ替えのあと、popup が古い key の override を読めない

- 再現: Linear で TESTB の tRNA の label を Off にし、record を並べ替える。図の編集は保たれるが、tRNA の popup は
  Label visibility を Default と表示する。override の key は `f88047061_record_2` のままで、feature は `_record_1` になった。
  popup から Default に戻しても古い key に届かない。
- 原因: RC-5。
- 修正案: Q4。

### OV-13: error-producer coverage が 4 つの owner を数えていない

- 内容: RC-6。OV-05〜07 の throw は `label-actions.js` にあり、shrink-only の baseline に入っていなかった。
  svg-result-ingestion の 22、candidate-render の 3、preview-runtime の 7 は、コードを読んだ限りすべて内部の不変条件で、
  利用者の操作だけでは届かない（41 行目は分類済みの `RESULT_INVALID`）。
- 修正案: 4 つの owner を baseline に加え、`requireUniqueEditableLabelBindings` を `diagnosticError` にする（PR-2）。

### OV-14（#757 で修正済み）: テスト helper が撤去した Enable Labels を待っていた

- #757 の CI で `tests/web/visual-state-regressions.playwright.spec.js` の B-06 が 2 回とも timeout した。
  `tests/web/helpers/visual-state.cjs::label` が Enable Labels のダイアログを待ち、新しい Label Not Shown が開いたまま
  次の click を遮った。#757 に `e1b3d701` を足し、Show this label を選ぶようにした（B-06 と
  `mobile-feature-popup` をローカルで確認）。ボタンの accessible name にはアイコン文字が入るので、正規表現で探す。

## 再現しなかった候補

- **Feature visibility On と feature type の除外（A）**: tRNA を Features から外しても、On にした tRNA は描かれる
  （`gbdraw/features/visibility.py::should_render_feature` で show が勝つ）。
- **Label On と Blacklist（A）**: Blacklist に当たる feature の Label On は live と Generate の両方で表示される。文字だけの編集は
  Label Not Shown を開き、Keep hidden で隠れたまま。Blacklist の語を含む文字に変えても表示される（filter は元の文字を見る）。
- **「This feature only」の fill と palette（A, E）**: palette を変えても個別の色は live と Generate で残る。live の色ルール照合は
  Python の優先順位を使う（`app/rule-matching.js::firstMatchingRule`）。
- **mode の往復（E）**: Circular → Linear → Circular で Label On、Feature visibility Off、色の override が残る。
- **busy 中の操作（C）**: reflow 中は popup の操作が `semanticMutationAvailable` で無効になるので、利用者の操作では無言の no-op にならない。
- **batch Results（B）**: 片方の Result の visibility と label 文字の編集は、もう片方に漏れない。
- **複数 record の Linear（crop なし）（B）**: 描画 hash と元座標の hash が等しく、#757 の label visibility の `hash` 行は逆相補でも照合される。
- **Feature placement の identity（B）**: (`record_key`, `biological_feature_id`) を元の catalog から解決するので、crop・逆相補・複製で正しく当たる（コードを読んだ確認）。
- **TSV の書き出しと読み込み（B）**: Download Label TSV などは Generate と同じ builder を使い、読み込み直すと Web では同じ結果になる。
  行の誤りは OV-01、OV-05 と同じ原因で、別の不具合ではない。
- **重複 feature で Keep hidden でも label が出る（D の既知候補）**: 現在は出ない。OV-05 の行が照合されないためで、OV-05 を直すと出うる。
- **生成後の他の `throw new Error`（D）**: `preview-runtime.js`、`candidate-render.js`、`svg-result-ingestion.js` の残りは内部の不変条件。
  underlay の feature の fill と visibility Off は Generate に成功する。
- **未使用のコード**: `featureVisibilityRulesFromOverrideCache` と `buildExactHashFeatureVisibilityRule` は `_record_<n>` 付きの値を作るが、実行時にはどこからも呼ばれない。

## F: label reflow の `reason` の受け渡し

- 現状: `label-actions.js::queueLabelReflow` が `labelReflowRequestReason` / `labelReflowForceRequestReason`（`state.js` の ref）に
  文字列を書き、`watchers.js` がそれを `runLabelReflow(reason)` に渡し、`pendingReflowReason` を経て
  `runLabelReflowCandidate({ reason })` に届く。受け取った側は使わない。分岐もログもない。
- 読むのはテストだけ: `tests/web/feature-label-visual-unit.test.mjs`（どの reflow を積んだかを reason の文字列で確かめる）、
  `tests/web/error-boundary.playwright.spec.js`、`tests/web/session-operation-consistency.playwright.spec.js`（seq の前に reason を書く）。
- 判断: **削除する。** 診断に使われておらず、残すと「reason で挙動が変わる」と読み違えやすい。テストは seq の増分と force の別で確かめれば足りる。
  振る舞いは変わらない私的な整理なので、Product の判断は要らない（IMPLEMENT_EXISTING_AUTHORITY）。単独の小さな PR にする（PR-6）。

## 修正の計画

| PR | 対象 | 分類 | 状態 |
|---|---|---|---|
| PR-1 | OV-08: request には現在の mode の record の placement だけを載せる | IMPLEMENT_EXISTING_AUTHORITY（R2） | 実装する |
| PR-2 | OV-05〜07 の失敗の分類と OV-13: 必須 binding の失敗を `diagnosticError` にし、4 owner を coverage に加える | IMPLEMENT_EXISTING_AUTHORITY（R6） | 実装する |
| PR-3 | OV-09: 対応しない lane placement の診断（feature と次の行動） | IMPLEMENT_EXISTING_AUTHORITY（R6） | 実装する |
| PR-4 | OV-05 の色: 「This feature only」の行から `__instance_` を外す | IMPLEMENT_EXISTING_AUTHORITY（PD-OI-069） | 実装する |
| PR-6 | F: `reason` の受け渡しを削除 | IMPLEMENT_EXISTING_AUTHORITY（私的な整理） | 実装する |
| PR-7 | OV-06、OV-07、OV-10: Label On を Apply する時点のダイアログ（Q1、Q2）。描けない On で Generate を失敗させない | Owner の決定 Q1、Q2 | PR-2 のあと |
| PR-8 | OV-09: lane placement を使えなくする設定変更の時点のダイアログ（Q3） | Owner の決定 Q3 | PR-3 のあと |
| DOC-1 → PR-9 以降 | OV-01〜04、OV-05 の label、OV-11、OV-12: 個別 override を (`record_key`, 元の feature の ID と instance) で指す設計文書、続いて実装（PR-5 を置き換える） | Owner の決定 Q4 | 設計文書から |

PR-1 と PR-3 は同じ根本原因（RC-4）なので 1 つの PR にまとめる。

## Owner の決定（2026-10-04）

Owner（satoshikawato）が 2026-10-04 に、次節の選択肢から選んだ。依存する実装はこの決定を引用する。

| 質問 | 選択 | 選択肢の文言（提示したまま） |
|---|---|---|
| Q1 | B: Apply 時に確認 | 描けない場合は Apply の時点でダイアログを出す（例: underlay には label を描けない → Keep without label / Cancel）。On は描けるときだけ効く。 |
| Q2 | A: ダイアログで選ぶ | "Show feature and label" / "Keep feature hidden" を選ばせる。後者は On を保存し、feature を表示したときに効く。popup の案内も直す。 |
| Q3 | B: 変更時にダイアログ | Track type などを変える時点で "Reset N placements to Auto" / 変更を取り消す を選ばせる。 |
| Q4 | A: 元の feature で指す | Feature placement と同じく (record_key, 元の feature の ID と instance) を key にし、Python が元の catalog から解決する。crop・逆相補・並べ替え・複製でも同じ feature に当たり続ける。Session と生成表の形式が変わるので、設計文書を先に作ってから実装する。色の This feature only は PD-OI-069 のまま。 |

この決定で、修正の計画の PR-5 は Q4 の設計文書と実装に置き換える。OV-06、OV-07、OV-10 は Q1 と Q2 のダイアログ、
OV-09 の後半は Q3 のダイアログで扱う。PR-2 と PR-3 の分類したエラーは、ダイアログを通らない経路（Session の読み込み、
CLI、想定外の不一致）の安全網として残す。

## Owner に確認したこと（提示した選択肢）

### Q1: Python が描けない Label On（OV-06 underlay、OV-07 Embedded Only）

- A（推奨）: 明示した On を優先し、Python が描く。underlay の feature には通常の feature と同じ外側の label を、
  Embedded Only で収まらない label には外側の label を描く。全体設定は変えない。
- B: 描けない場合は Apply の時点でダイアログを出す（例: "Labels are not drawn for underlay features" と、
  Keep without label / Cancel）。On は描けるときだけ選べる。
- C: Generate は成功させ、描けなかった label を通知に列挙する。

### Q2: 隠した feature の Label On（OV-10）

- A（推奨）: Label On を Apply したとき feature が隠れていれば、ダイアログで "Show feature and label" /
  "Keep feature hidden" を選ばせる。後者は label の On を保存し、feature を表示したときに効く。popup の案内も直す。
- B: feature が隠れている間は Label visibility を無効にし、理由を表示する。

### Q3: slot の変更で使えなくなる lane placement（OV-09）

- A（推奨）: PR-3 の分類した Generate エラーに、該当する placement を Auto に戻す操作を付ける。
- B: Track type、Separate strands、slot を変える時点でダイアログを出す（"Reset N placements to Auto" / 変更を取り消す）。
- C: Python が使えない lane を Main で描き、通知する。

### Q4: 個別 override の identity（OV-04、OV-05 の label、OV-11、OV-12、PR-5 の形）

- A（推奨）: Feature visibility、label の文字と表示、annotation の対象を、Feature placement と同じ
  (`record_key`, `biological_feature_id`) を key にした行で Python に渡し、Python が元の catalog から解決する。複製と同じ座標の重複は
  instance（record_key と source feature の順番）で区別する。crop・逆相補・並べ替えのあとも編集は同じ feature に当たり続ける。
  Session と生成表の形式が変わるので、別の設計文書を作ってから実装する。色の「This feature only」は PD-OI-069 のまま。
- B: 今の描画 hash の契約を保つ（PR-5 で Web を Python に合わせる）。crop・逆相補・並べ替えで描画 ID が変わるときは、
  Web が override の key を新しい描画 ID に付け替える。複製の区別は Python に instance 番号の照合を足す。
- C: 今の契約を保ち、複製は編集を共有し、crop などを変えると編集が外れることを文書にして、外れたときに通知する。
