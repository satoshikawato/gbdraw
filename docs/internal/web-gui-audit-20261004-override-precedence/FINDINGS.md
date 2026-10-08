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
- OV-17 以降（2026-10-05 の残る課題の修正中に見つけたもの）は、クラスに **G**（Preview と standalone SVG、Web と Python の読み込み・描画の食い違い）、
  **H**（遷移と key の抜け）、**I**（アクセシビリティ）を足した。
- 日付は日本時間。PR の merge 日は `gh pr view` の `mergedAt`（UTC）を日本時間にしたもの。
- 重大度: **図の誤り** > **Generate の失敗** > **操作の失敗**（popup の操作が失敗する。OV-17 で足した）> **無言の no-op** > **文言の不一致**。
- 判断の基準: Owner 方針（2026-10-04）「popup で明示した個別設定は全体設定に勝つ。個別編集の副作用で全体設定を書き換えない。
  編集が見えないままになるときは具体的な選択肢のダイアログで確認する」、OIPC-C03（受け付けた値は必ず使われるか、
  実行前に拒否される）、PD-OI-069（Web が作る色ルールは Python が照合できる値にする。label の instance 単位の編集は維持する）、
  R1–R12（`gbdraw/web/CLAUDE.md`）。
- JS のパスは `gbdraw/web/js/` からの相対パス。

## 要約

不具合は 14 件（図の誤り 5、Generate の失敗 5、無言の no-op 2、文言・状態の不一致 1、検査の抜け 1）。
OV-15 は修正 PR のレビュー中に見つけた、`dev` に前からある不具合。OV-16 は修正中に見つけた error boundary の潜在的な欠陥で、
利用者の操作で起きる経路は見つからなかった（件数に含めない）。
このほか、#757 自身のテスト helper の取り残し 1 件をこのセッションで直した（OV-14）。
各 OV の修正 PR と残る課題は「修正の状況（2026-10-05）」にまとめた。

2026-10-05 の「残る課題」の修正（R-1〜R-7）で、さらに 18 件（OV-17〜OV-34）を見つけた。OV-17 は Owner が gbdraw.app の Gallery で報告した。
OV-29、OV-31、OV-32 は見つけた PR の中で直した。残る課題そのものではなく、同じ種類の不具合（個別の設定と全体の仕組みの食い違い）を
直す途中で見つかった。

2026-10-06 の Phase E（owner-coupling、#813）と、その途中で行った調査で、さらに OV-35〜OV-48 を見つけた（OV-37 は不具合ではなかった。
OV-41 は見た目だけ）。OV-35 は #810 の修正を #824 に載せ直したときに見つけた。OV-36 は E5（#832）、OV-38〜OV-40 は残る観察の調査、
OV-42〜OV-45 は live と Generate の食い違いを調べるテスト（#839）、OV-46 は OV-45 の修正中に、OV-47 は OV-42〜OV-44 の修正（#857）の確認中に、OV-48 は #871 の CI の失敗を調べて見つけた。

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
| OV-15 | A, B | 図の誤り | feature_type が `*` の色ルールか Feature visibility の行があると、GFF3 の入れ子の feature（gene の下の CDS など）が Generate で描かれない | RC-7 |
| OV-16 | D | 潜在（Generate の失敗） | Generate が `runAnalysisInternal` の `try` の前で失敗すると、元の Result に戻す処理が終わらず、「Generating Diagram...」のまま止まる | RC-8 |
| OV-17 | G | 操作の失敗 | standalone SVG の collinearity block（cluster）の popup で span の配列をコピーできない（"Match feature endpoint identity is invalid."）。Gallery の collinear 例で起きる | RC-9 |
| OV-18 | G | 状態の不一致 | standalone SVG の block popup が、group の member でない anchor protein まで各 group の下に並べ、複数 anchor の Query/Subject を 1 行にまとめる。Web app の popup と内容が違う | RC-9 |
| OV-19 | C | 無言の no-op | Features → Feature Visibility の規則の編集（追加、変更、移動、削除）が Result に反映されず、"Applies on Generate" の表示もない | RC-10 |
| OV-20 | G | 図の誤り | 表の値の `#` で行が切られる。visibility の行は Generate が `Missing values` で失敗し、label の文字 `Gene #1` は `Gene ` になる | RC-11 |
| OV-21 | H | 図の誤り | 選択から作った annotation の対象が record key だけで選ばれ、同じ record key を持つ別 mode の Generate にも描かれる | RC-12 |
| OV-22 | G | Generate の失敗 | styling と override の表を pandas の CSV 引用で読むので、`"` で始まる値（`"lead`）が ParseError になるか、`"quoted"` の引用符が消える。Web は値をそのまま書く | RC-11 |
| OV-23 | G | 図の誤り（軽微） | Web の whitelist の書き出しが keyword 中の tab を直さず、列がずれる（visibility と override の書き出しは直す） | RC-11 |
| OV-24 | G | Generate の失敗 | annotation の表を Python が `csv.reader(strict=True)` で読むので、`"lead` の label を含む表が失敗し、引用した label は Web と Python で値が違う | RC-11 |
| OV-25 | G | 図の誤り | styling の表の行が列を余分に持つと、pandas の index のずれで黙って別の列として読まれる（whitelist `CDS\tproduct\ttwo\twords` の feature_type が `product` になる） | RC-11 |
| OV-26 | G | Generate の失敗 | annotation の表の `#` で始まる行を、Web の Import は飛ばし、Python は行として読んで列数エラーにする | RC-11 |
| OV-27 | G | Generate の失敗 | Specific colors と Default colors の表の `#` で始まる行を Python が行として読み（`Missing values`）、Web は飛ばす | RC-11 |
| OV-28 | G | 状態の不一致 | Web の whitelist、priority、色の表の読み込みが、列の足りない行を飛ばし余分な列を無視する（CLI は OV-25 の後は拒否する） | RC-11 |
| OV-29 | H | 図の誤り | Circular の depth track を削除しても `circular_track_slots_axis_index` が下がらず、次の行が Axis をまたぐ | RC-13 |
| OV-30 | I | アクセシビリティ | アイコンだけのボタンの accessible name が Phosphor の glyph（私用領域の文字）になる。`aria-hidden` のないアイコンが 245 のうち 224 ある | RC-14 |
| OV-31 | H | 状態の不一致 | 自動 rerender が catalog を捨てるので、On にした feature が catalog の行なしで描かれる。Undo しても描かれたまま | RC-13 |
| OV-32 | H | 状態の不一致 | 個別の visibility を持つ描かれない source feature と、type を絞った GFF3 の Off の行の feature が、biological catalog から落ちる | RC-13 |
| OV-33 | G | Generate の失敗 | crop した record、同じ座標の重複 CDS、別の record があると、Web は `RENDER_FAILED`、CLI は "Feature metadata identity does not agree with rendered SVG ID" で失敗する | RC-2 |
| OV-34 | H | 無言の no-op | Auto Reflow を切ると、規則の削除など「描かれていない feature を描かせる」編集が、Generate まで何も表示されない | RC-10 |
| OV-35 | H | 図の誤り | Auto Reflow を切ると、Feature Visibility の規則で feature を隠しても、その label が Result に描かれたまま残る | RC-10 |
| OV-36 | H | 状態の不一致 | live の編集の label の rerender が失敗したあと、Generate が成功しても "Live edit failed: …" の note が Result の上に残る | RC-13 |
| OV-38 | G | 操作の失敗 | 0.13.0 の Gallery Session が Web app で読み込めず、"The operation failed without recognized diagnostic information." と出る | RC-16 |
| OV-39 | G | 図の誤り | Worker が起動できないとき、feature の編集を持つ Session 31〜33 の Load が、その編集をすべて「対応する feature がない」として捨てて成功する。次の Save で編集が失われる | RC-16 |
| OV-40 | G | 操作の失敗 | Session 31〜39 の Default colors、Label whitelist、Qualifier priority の行に余分な cell（tab を含む値）があると、Load が表も Session も名指さないエラーで失敗する（OV-28 の副作用） | RC-16 |
| OV-41 | H | 文言の不一致（軽微） | Session を読み込むと、検査する前の Record が "(not found)" と表示される（Records not inspected の下） | — |
| OV-42 | H | 状態の不一致 | Auto Reflow を切ると、feature の visibility の編集で legend が変わらない（Generate は型の行の順を変え、描かれない型の行を落とす） | RC-15 |
| OV-43 | H | 状態の不一致 | live の specific color 規則が既定の caption（`CDS`）を残して規則の行を末尾に足す。Generate は規則の行を先頭にし、残りを `other proteins` にする | RC-15 |
| OV-44 | H | 状態の不一致 | batch で、1 つの Result で編集した規則の legend の行が、もう一方の Result を表示したときに届かない | RC-15 |
| OV-45 | H | Generate の失敗 | Circular の batch で、前回の Generate で塗った feature に 2 回目の「This feature only」の色を付けると、次の Generate が `UNKNOWN` で失敗する | RC-15 |
| OV-46 | H | Generate の失敗 | 1 つの Result だけが描く generated の legend の行に Legend だけの色を付けると、Generate が `RESULT_INVALID` で失敗する（OV-44 と同じ種類） | RC-15 |
| OV-47 | H | 状態の不一致 | batch で、Result 2 を表示したまま Generate（または #857 の自動 rerender）をすると、Result 1 の Undo と "Sort by default" が Python の legend の順に戻らない | RC-15 |
| OV-48 | H | Generate の失敗 | 自動 rerender が Generate の実行中に始まると、Generate が Result もメッセージも出さずに終わる | RC-17 |

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
- **RC-7: GFF3 の「すべての type を読む」が入れ子を平らにしない。** `gbdraw/features/visibility.py::resolve_candidate_feature_types` は、
  色表か visibility 表に feature_type `*` の行が 1 つでもあると（action に関係なく）すべての type を読むと決める。
  `gbdraw/io/genome.py::load_gff_fasta` はそのとき BCBio の記録をそのまま返し、Parent でつながる子（gene の下の CDS）は
  `sub_features` に残る。type で絞る経路だけが `scan_features_recursive` で平らにする。描画は `record.features` しか見ない。
- **RC-8: 元の Result に戻す処理が、再 mount を待つだけで自分では bind しない。** `app/preview-runtime.js::restorePreviousSelectedResult` は
  readiness の receipt を待ち、それを作る bind は `app/watchers.js` の mount watcher だけが始める。候補を有効にする前の失敗では、
  戻す Result が表示中のものと同じなので Vue は再 mount せず、bind が始まらない（R10: watcher の実行を不変条件の仕組みにしない）。

- **RC-9: standalone SVG の popup が、Web app の popup の block の規則に従っていない。** standalone の runtime
  （`services/standalone-interactivity-assets.js`）は、collinearity block の cluster が複数の endpoint を claim すること
  （`queryFeatureReferences`）と、group の member でない anchor を並べないことを、Web app の popup と別に実装していた。
- **RC-10: Feature Visibility の規則の編集に、投影の遷移がない。** `app/feature-editor/visibility-actions.js` の 7 か所が
  `featureVisibilityManualRules` を書き、どれも Python に規則の照合を頼む処理（`prepareDrawn`）と投影を呼ばなかった。
  Auto Reflow を切ると、描かれていない feature を描かせる編集は、形（geometry）を作る自動 rerender を呼ばない限り表示できない。
- **RC-11: 表を読む側と書く側が、表ごとに別の規則を持つ。** pandas の `comment="#"`、既定の CSV 引用、`csv.reader(strict=True)` は
  読む表ごとに違い、Web の書き出しと読み込み（`file-imports.js`、`run-analysis.js`）も別の規則だった。コメント行、値中の `#`、引用符、
  列数、tab の扱いが Web と Python でそろっていない。
- **RC-12: 個別の編集の key が、mode を持たない。** R-1 と同じ。annotation の対象の `featureIdentity` が record key だけで選ぶ
  （`app/annotations/state.js::annotationOptionsPayload`）。
- **RC-13: 遷移が、すべての書き手を通らない。** slot を変える編集、自動 rerender、catalog の admission が、それぞれ自分の状態だけを
  更新する（`removeCircularDepthTrack`、`runLabelReflowCandidate`、`gbdraw/web_support/feature_metadata.py` の source の skip）。
- **RC-14: アイコンの markup に規則がない。** Phosphor のアイコンは CSS `::before` の glyph で、`aria-hidden` と `aria-label` の付け方を
  検査する仕組みがなかった。
- **RC-15: legend の行を、Result ごとの、図から導く状態として扱っていない。** Python は legend の行を描いた feature から導く
  （型の行の順、描かれない型の行の削除、残りの caption）。live の遷移は feature を隠す・色を付けるだけで、導く行を作り直さない
  （`projectFeatureVisibility`、`commitSpecificRules`）。batch では、規則の legend の行は表示中の legend にだけ書かれ、
  Generate の compiler（`app/candidate-render.js`）は行をすべての Result に必須とする。
- **RC-16: 古い Session の読み込みが、現在の厳密な検査をそのまま当て、失敗の種類を区別しない。** `validateImportedCircularTrackSlots` は
  古い writer が書いた null の欄を拒否し（OV-38）、`extractSessionSourceFeatures` は runtime の失敗とデータの誤りを同じ `null` にし（OV-39）、
  OV-28 の厳密な reader は古い writer の tab を含む値を拒否する（OV-40）。どれも、失敗の意味を作る側が分類していない（R6）。
- **RC-17: Generate と自動 rerender が 1 つの世代番号を共有し、互いの実行を確かめない。** `app/run-analysis.js` の Generate と
  自動 rerender は、どちらも `latestGenerationToken` を進め、番号が進んだ操作は `stale` として何も言わずに終わる。
  rerender は `processing` を見ないので、Generate の後に始まった rerender が Generate を捨てる。

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

### OV-15: feature_type `*` の行で、GFF3 の入れ子の feature が描かれない

- 再現（CLI、`dev` `e40e058e`）: `tests/test_inputs/NC_013668.gff3` と `.fasta` を Circular で描く。表がなければ 136 の feature を描く。
  `--feature_visibility_table` に 1 行（record_id `*`、feature_type `*`、qualifier `product`、value `.*`、action `show`）を足すと
  2 になり、CDS がすべて消える。警告は出ない。feature_type を `CDS` にした行なら変わらない。色表（`-t`）の feature_type `*` の行でも、
  Linear でも同じ。Python API の `read_gff()` を `features` なしで呼んだ場合も 2 になる。
- 再現（Web）: 同じ GFF3 と FASTA で Generate（136）。手入力の Feature visibility 規則を既定の値（record `*`、type `*`、`product`、Off）のまま、
  どれにも当たらない value で足して Generate すると 2 になる。popup の catalog は入れ子の CDS を持つ（134）ので、popup と図が食い違う。
  popup の Feature visibility Off は feature 自身の type と record の行を作るので、この不具合に当たらない（136 → 135）。
  Web の色規則の画面は具体的な type しか選べないが、読み込んだ色の TSV は `*` を持てる。
- 期待: 受け付けた行はその行の意味どおりに使われ、関係のない feature を消さない（OIPC-C03）。
- 原因: RC-7。#771 は identity の行についてだけ、すべての type ではなく必要な type を足して読み直すことで避けていた
  （`gbdraw/api/request_render.py::_load_request_records`）。
- 修正: すべての type を読むときも、type で絞るときと同じく入れ子を平らにする（#781、IMPLEMENT_EXISTING_AUTHORITY: OIPC-C03）。
  逆相補した GFF3 の record では popup の catalog の並び（と位置の `id`）が描画の開始位置の順になり、GenBank の record と同じになる。
  `feature_index`、`svg_id`、selector、hash は変わらない（Owner に委ねられた選択として #781 に記録）。

### OV-16: Generate が早い段階で失敗すると、元の Result に戻す処理が終わらない（潜在）

- 見つけた経緯: #784 の修正中に、古い metadata（`record_key` がない）で `isCurrentFeature` が投げ、Generate が止まった。#784 はその原因を直した。
- 再現（`dev` `4ae723bd`、fault injection）: 既存の `__GBDRAW_TEST_HOOKS__.onSessionLifecycleEvent` で `generation-input-resolution-start` に
  1 回だけ例外を投げる。Session の読み込み後でも、普通に Generate した後でも、「Generating Diagram...」のまま止まる。
  `processing` と History の `restoring` が true のまま、Undo と Save Session は無効、エラーは出ない。Cancel でも終わらない
  （`cancelRunAnalysis` は Generate 自身の token の expectation だけを拒否する）。trace は `artifact.rollback-started` →
  `preview.restore-bind-started` で止まる。
- 利用者の操作で起きる経路: 見つからなかった。`try` より前で投げうる手順（`refreshCircularRecordOrder`、`prepareLinearRecordCatalog`、
  `validateDepthInputPresence`、`resolveLinearComparisonPlan`、`isCurrentFeature`）を調べ、古い Session の fixture でも試した。
  "Rotate record to feature" の失敗時の戻しにも同じ欠陥がある。
- 関連（`dev`、#784 で解消）: v40 より前の Session を読み込むと feature に `record_key` がなく、`isCurrentFeature` が false を返すので、
  source に結びついた編集があると `sourceReplaced` が true になる。
- 原因: RC-8。
- 修正: 戻した後も同じ root が mount されたままで、active な runtime が戻した Result のものなら、戻す処理が自分でその root を bind する
  （#785、IMPLEMENT_EXISTING_AUTHORITY: R6、R10、PD-OI-055）。失敗は分類したエラーとして表示し、元の Result と History を保ち、再試行できる。

### OV-17: collinearity block の popup で span の配列をコピーできない（standalone SVG）

- 再現（`dev` `cedbe18d`、Owner が gbdraw.app の Gallery `hepatoplasmataceae_collinear` と `vibrio-harveyi-group-collinear` で報告）:
  standalone SVG の block の popup から span の配列をコピーする。cluster の block は失敗し、singleton の block は成功する。
- 期待: どちらの block でもコピーできる。
- 実際: "Match feature endpoint identity is invalid."。
- 原因: RC-9。cluster の block は複数の endpoint を claim する（Gallery の block 0 は 15）が、`resolveEmbeddedMatchSource` の guard が
  endpoint を 1 つだけ解決する `resolvedCatalogMatchFeature` を使っていた。Web app の popup は `data-*` 属性から span を作るので影響を受けない。
- 修正: #794（IMPLEMENT_EXISTING_AUTHORITY）。`catalogMatchEndpointResolved` が、複数の claim では catalog の展開が全部解決したか
  どうかの結果を信頼する。`docs/images` の h-cli-13 と t-py-08 の interactive SVG を作り直した。
- Gallery: Gallery の 10 個の例 SVG はすべて runtime を埋め込む（collinear の 2 つは 15〜16 MB）。#794 では作り直さず、
  runtime の PR がそろった後の 1 回の Gallery の更新（#851）に回した。gbdraw.app で直るのは、`main` への昇格（Owner が決める）と
  Gallery の更新の後。

### OV-18: standalone の block popup が、member でない anchor と複数 anchor の行を Web app と違う形で出す

- 再現（OV-17 の修正中に見つけた。standalone のみ）: block の popup の "Similarity groups covered" が、各 group の Query/Subject の member の下に
  anchor の protein をすべて並べる。複数 anchor の Query/Subject の feature 行が `a0;a1;a2` の 1 行になる。
- 期待: Web app の popup と同じ内容。member でない anchor は並べず、anchor ごとに 1 行。
- 原因: RC-9。`materializedBlockMemberLabels` が member でない feature の描画 protein ID にフォールバックする
  （Web の `resolveBlockMemberLabels` は member でないものを飛ばす）。`materializedMatchFeatureRows` が複数の claim を 1 行にまとめる。
- 修正: #797。Preview の binder と standalone SVG を同じ fixture に mount して、同じ popup の内容を要求する Playwright のテストを足した。
  `docs/images` の h-cli-13 と t-py-08 を作り直した。Gallery の SVG は presentationScope のない古い metadata の group を持つので、
  Gallery の更新（#851）が要る。

### OV-19: Feature Visibility の規則の編集が Result に反映されない

- 再現（`dev` `8e5234ac`）: `tests/fixtures/forced_label_underlay.gb`（Circular）を Generate し、規則 `FORCEDLBL` / `CDS` / `locus_tag` / `^fl1$` / Off を足す。
- 期待: FL1 が消える（Generate が描く図と同じ。PD-OI-066、OIC-027、R1、R3、R10）。
- 実際: Generate、Undo、Result の切り替え、label の rerender のどれかまで FL1 が描かれたままで、パネルは "Applies on Generate" と表示しない。
- 原因: RC-10。`setFeatureVisibilityRuleField`、`addFeatureVisibilityRule`、`moveFeatureVisibilityRule*`、`removeFeatureVisibilityRule` が規則を書くだけで、
  どの遷移も照合の準備と投影を呼ばなかった。
- 修正: #810 は merge せずに 2026-10-06 に閉じ、#824（Phase E の E1「visibility⇄label」）が同じ挙動を #807 の rerender の上に載せ直して入れた。
  `editFeatureVisibilityRules` を規則の編集の遷移にし、投影を `projectFeatureVisibility` 1 つにまとめた。Generate が拒否する regex は
  draft に残し、Result は変えず、Generate と同じ `REGEX_SYNTAX` のエラーをパネルに出す。ガードのテスト（`feature-visibility-rule-writers.test.mjs`）が、
  規則の書き手の一覧を固定する。
- 同じ PR で直した軽微な不具合: "Load Feature Edits TSV" が、照合の準備なしに visibility を投影していた（ある規則で隠れる feature の個別編集を消す TSV を読むと、
  Generate は隠すのに live は描いたままだった）。
- 未確認だった点: Auto Reflow を切って規則の編集で feature を隠すと、弧は消えるが label の文字が SVG に残る。→ 確かめて OV-35 とし、#824 で直した。

### OV-20: 表の値の `#` で行が切られる

- 再現（CLI）: visibility の表に値が `foo#bar` の行を入れると `ValidationError: Missing values`。label override の `Gene #1` は `Gene ` になって描かれる。
  Web の Generate は `TABLE_INVALID`。
- 期待: `#` は、行頭（最初の空白でない文字）にあるときだけコメント。値の中の `#` は値の一部。
- 原因: RC-11。`gbdraw/features/visibility.py::read_feature_visibility_file` などが pandas の `comment="#"` を使い、行の途中から切る。
- 修正: #798（R-2 の関連）。`gbdraw/io/table_text.py` を新設し、行頭のコメントだけを飛ばす。`utf-8-sig` も扱う。
  `gbdraw/` で pandas の `comment=` を使うことを禁じるテストを足した。

### OV-21: 選択から作った annotation が、同じ record key の別 mode にも描かれる

- 再現（`dev`）: `lambda_basic_linear` を読み込み、Circular の `NC_001416` を足した multi-record の canvas（両方の mode で record-1）で、
  Circular で選択から annotation を作る。Linear の request にも届き、描かれる。
- 期待: 対象は、それを描く record と mode にとどまる（設計 Q4 3.2、R2）。
- 原因: RC-12。`annotationOptionsPayload`（`app/annotations/state.js`）が record key だけで選ぶ。
- 修正: #806（R-1 の annotation 側）。draft の `featureIdentity` の対象が `scope` を持ち、`annotationOptionsPayload(sets, mode, records)` が
  絞って外す（request schema 9 と CLI の表は変えない）。`scope` のない対象は `INPUT_INVALID`。

### OV-22: styling と override の表の引用符

- 再現: label の表で値 `"quoted"` が `quoted` になる。先頭が `"` で閉じない値は ParseError。Web の GUI で規則の値 `"lead` を書くと、Web は値をそのまま書くので Generate が失敗する。
- Owner の決定（2026-10-05）: A（styling の表は引用を切って読む。`"` は値の一部。文書に書く）。
- 修正: #800。`read_literal_table`（`quoting=csv.QUOTE_NONE`）を、visibility、label override、label filter、qualifier priority、specific colors、
  default colors に使う。feature override と placement の表（`label_text` の CSV 引用を文書化済みで、Web の書き出しが引用する）、
  cli_tables の manifest、canonical な request の表、BLAST（もとから `QUOTE_NONE`）、annotation の表（OV-24）は変えない。

### OV-23: Web の whitelist の書き出しが keyword 中の tab を直さない

- 再現: Web の whitelist の書き出し（`manual_wl.tsv`、`run-analysis.js`）は、keyword の tab を直さずに書く。visibility と override の書き出しは tab を直すが、whitelist と qualifier priority は直さなかった。tab を含む keyword は列数がずれる。
- 原因: RC-11。
- 修正: #801。`normalizeTsvCell`（`utils/tsv-cell.js`）と `file-imports.js` の serializer（whitelist と qualifier priority の 4 か所の書き出し）に集め、ガードのテストを足した。

### OV-24: annotation の表の引用符

- 再現: Web の Region Annotations の download は label をそのまま書き、Web の Import は分割し、Python の `read_annotation_table` は `csv.reader(strict=True)` で読む。
  label `"lead` と別の行を含む表は Python で `ValidationError`。引用した label は Web と Python で値が違う。
  `test_annotation_file_preserves_quoted_cells_and_blank_lines` が CSV 引用を固定していた。文書は引用について何も書いていない。
- Owner の決定（2026-10-05）: A'（annotation の表の引用符は値の一部。OV-22 と同じ。そのテストを書き直す）。
- 修正: #802。`csv.QUOTE_NONE`。共有ベクトル `tests/fixtures/annotations/tsv-import-cases.json` を Python と Web の両方が通す。

### OV-25: styling の表の列を余分に持つ行が、黙って別の列として読まれる

- 再現: whitelist の行 `CDS\tproduct\ttwo\twords` の feature_type が `product` になる。pandas の index のずれ。
- 期待: 受け付けた値は使われるか、実行前に拒否される（OIPC-C03）。
- 修正: #803。`read_literal_table(label, filepath)` が行の列数を調べ、余分な列に行番号つきの `ParseError`。重複していた 2 つの事前検査を消した。

### OV-26: annotation の表の `#` の行

- 再現: Web の Import は `#` で始まる行を飛ばし、Python は行として読んで列数エラー（または header の欠落）にする。文書はコメント行を約束していない。
- Owner の決定（2026-10-05）: 推奨どおり、Python が行頭の `#` の行を飛ばす（`read_table_lines`）。文書に書く。
- 修正: #812（2026-10-05 に merge）。

### OV-27: Specific colors と Default colors の表の `#` の行

- 再現: 両方の reader が `#` の行をデータとして読み、`Missing values`。Web の `parseColorTable` と `parseSpecificRules` は飛ばす。文書がコメントを約束するのは whitelist、label override、feature visibility だけ。
- Owner の決定（2026-10-05）: OV-26 と同じ。Qualifier priority も同じ規則にする（確認済み）。
- 修正: #812（2026-10-05 に merge）。`read_literal_table` が path を `read_table_lines` 経由で読む。`docs/REFERENCE/input-formats-and-tsv-schemas.md` を直した。

### OV-28: Web の表の読み込みが、列の足りない行と余分な列を許す

- 再現: `parseWhitelistRules`、`parsePriorityRules`、`parseColorTable` が、短い行を飛ばし、余分な列を無視する。CLI は OV-25 の後は拒否する。
- 修正: #804。ちょうど 3、2、2 列を求め、`TABLE_INVALID` の reason `FIELDS` を使う。注意: Session の復元も同じ parser を使うので、保存された列数の合わない行（OV-23 の前の whitelist の tab を含む keyword など）は、この診断で復元に失敗する。

### OV-29: depth track を削除すると、次の行が Axis をまたぐ（#805 で修正）

- 原因: `removeCircularDepthTrack` が `circular_track_slots_axis_index` を下げなかった。
- 修正: #805（R-3 の中で）。

### OV-30: アイコンだけのボタンの accessible name が glyph になる

- 再現（R-3 の agent が見つけた）: 「Move outside Axis」のようなアイコンだけのボタンの name が、Phosphor の glyph（私用領域の文字）になり、`getByRole('button', { name: 'Move outside Axis' })` が見つからない。`aria-hidden` のないアイコンが 245 のうち 224。
- 原因: RC-14。
- 修正: #809（R12。2026-10-05 に merge）。224 のアイコンに `aria-hidden`、51 のアイコンだけのボタンに `aria-label`（静的な title と同じ。title が状態で変わるときは安定した動作名）。Editor の toggle の name が状態で変わっていた
  （"Open editor" / "Close editor"。R12 違反）ので "Editor" と `aria-expanded` にした。ガード: `tests/web/icon-accessible-names.test.mjs` と accessibility の spec 4 件。

### OV-31: 自動 rerender が catalog を捨てる（#807 で修正）

- 原因: `runLabelReflowCandidate`（`app/run-analysis.js`）が `results.value` だけを入れ替え、catalog の admission を捨てた。
  On にした feature が catalog の行なしで描かれ、popup、label の binding、一覧、Undo の投影が届かなかった。Undo しても描かれたままで、Auto Reflow を切ると描かれない feature の On は Generate まで待った。
- 修正: #807（2026-10-05 に merge）。自動 rerender が catalog と、catalog から作る projection（`extractedFeatures`、`biologicalFeatures`）を Result と一緒に入れる。

### OV-32: 個別の visibility を持つ描かれない feature が catalog から落ちる（#807 で修正）

- 原因: `gbdraw/web_support/feature_metadata.py` が、描かれない `source` feature を biological catalog から外した。GFF3 の type の読み込みが On の行しか数えず、type を絞ったときの Off の行の feature を読まなかった。
- 修正: #807（2026-10-05 に merge）。個別の設定を持つ `source` は残し、GFF3 は identity の行が名指す type を On と Off の両方で読む（`_gff_types_named_by_visibility_rows`）。描画は変えない。

### OV-33: crop した record の同じ座標の重複 CDS で、別の record があると失敗する

- 再現: crop した record、同じ座標の重複 CDS、別の record の組み合わせ。Web は result の admission で `RENDER_FAILED`、CLI は "Feature metadata identity does not agree with rendered SVG ID `..._record_1__instance_1_...`"。
- 原因: `gbdraw/render/interactive_svg.py::_rendered_space_stable_id_candidates` が、`_record_<n>` を ID の末尾からしか外さなかった（RC-2 と同じ種類）。
- 修正: #808。上位集合を作る 1 行の修正。

### OV-34: Auto Reflow を切ると、規則の編集で描かれる feature が Generate まで出ない

- 再現: Auto Reflow を切り、前回の Generate で規則に隠された feature を、規則の削除で描かせる。形がないので、Generate まで何も出ない。
  #810 の本文は「同じことが個別の visibility、Undo、Result の表示にも当たる」と書いた（#807 が個別の visibility、History、Load Feature Edits TSV には直した）。
- 原因: RC-10。規則の編集の遷移が、#807 の「描かれていない feature が描かれるようになるとき自動 rerender を入れる」処理を呼ばない。
- 修正: #824（2026-10-06 に merge。#810 を #807 の後に載せ直した E1）。規則の編集の遷移が、結果として描かれていない feature を描く編集のとき、他の visibility の編集（#807）と同じく
  自動 rerender を要求する。Auto Reflow を切っていても、Generate で規則に隠された feature を規則の削除で描かせると、feature と label が live で描かれる。

### OV-35: Auto Reflow を切ると、規則で隠した feature の label が Result に残る

- 再現（E1 の途中、`6d29c934`。Auto Reflow を切る）: `tests/fixtures/forced_label_underlay.gb`（Circular、labels Out）を Generate し、
  規則 `FORCEDLBL` / `CDS` / `locus_tag` / `^fl1$` / Off を足す。
- 期待: FL1 の弧と label が一緒に消える（Generate が描く図と同じ。PD-OI-066）。
- 実際: 弧は `display="none"` になるが、label "alpha protein" は描かれたまま（label は未 bind）。Undo、Redo、表示中の別の batch Result でも同じ。
- 原因: RC-10。label の visibility の投影（`applyStoredVisibilityOverridesToSvg`）が、editor が bind した label だけを読む。editor は drawer、popup、label の編集が
  要るときだけ label を bind するので、規則の編集、History の apply、Result の表示、TSV の読み込みのあとに、隠れた feature の label が残った。
  popup を開くと先に bind するので、個別の編集では起きなかった。
- 修正: #824（2026-10-06 に merge）。投影が、Python が feature に結んだ label（`data-gbdraw-label-binding-schema="1"` と `data-label-feature-id`）も読む。Python は描く feature の label だけを描き、
  resolver が「描かれない」と答えたときだけ隠す（不明なら残す）ので、label は feature と一緒にだけ変わる。Phase E の E1 の中で直した。

### OV-36: live の編集の失敗の note が、Generate が成功しても残る

- 再現（E5 の作業中に見つけた。#838）: `tests/fixtures/web_batch_two_records.gb` を Circular の batch で Generate。`__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse` が投げるようにして
  label の rerender を強制し（失敗して note が出る）、hook を外して Generate する。
- 期待: 新しい Result がすべての編集を描くので、note は消える。
- 実際: "Live edit failed: direct edits already applied are kept; geometry may still need updating. … Retry the live edit or use Generate."
  （`labelReflowLastError`、`[data-live-application-feedback]`）が Result の上に残る。
- 原因: RC-13。Generate の開始は label の build の警告だけを消し、rerender の失敗を消す Generate の経路がなかった。
- 修正: #838（2026-10-06 に merge）。Generate が History の step を commit したあと、label owner の port `clearLabelBuildNotices({ rerender: true })` を呼ぶ。
  失敗、cancel、supersede された Generate は以前の Result を残すので、note も残す。

### OV-37: 調べた。不具合ではない

- 疑い（E4 の agent、未確認）: live の specific color 規則の caption が正規化されると、feature の色 override が rebind されず（Generate の `rebindRuleColorOverrides` だけが行う）、
  古い caption の key の override が次の Generate で当たらなくなる。
- 結果: 当たらない。`applyRulePreview` → `refreshFeatureOverrides` が、描かれた feature の override を改名後の規則から書き直す。
  Generate が残すのも描かれた feature の override だけ。

### OV-38: 0.13.0 の Gallery Session が Web app で読み込めず、原因が出ない

- 再現（2026-10-06 の調査。#834）: `tests/fixtures/sessions/BGC0000708-BGC0000713.v30.gbdraw-session.json.gz`（0.13.0 の Gallery Session、Session 30）を fresh な app で Load。
- 期待: 0.13.0 の GUI Session を読む、または原因を具体的に示して拒否する。`docs/SESSION_COMPATIBILITY.md` は Session 27〜30 を CLI replay だけとしていたが、
  Web app は 0.13.0 の CLI Session を読む（B15）。
- 実際: "The operation failed without recognized diagnostic information."（`UNKNOWN`）。実際の例外は "Custom Track Slots use obsolete field 'circular_track_slots[0].spacing'"。
- 原因: RC-16、RC-6。`validateImportedCircularTrackSlots` が、退役した欄 `spacing` を値と Custom Track Slots の有無にかかわらず拒否する（`main` で `94c158e9`、2026-08-24 から）。
  0.13.0 は slots が off でも全 row に `spacing: null` を保存した。投げるのは分類されない素の `Error`。`session-active-config-contract.js` は R6 の走査対象にないので、
  未分類の throw が coverage に出なかった。
- 修正: #834（2026-10-06 に merge。Owner の決定 2026-10-06）。Session 27〜33 の reader（`migrateSessionDataToCurrent`）が、slots が off のとき schema 4 の row の null の `spacing` を落とす。
  他の退役した欄、非 null の値、slots が on の Session は拒否し、`diagnosticError('TRACK_INVALID', …, reason: 'OBSOLETE_TRACK_FIELD')` で欄と track の row を名指す。
  `docs/SESSION_COMPATIBILITY.md` と `docs/REFERENCE/session-and-request-compatibility.md` の「CLI replay だけ」を、Web が保存した設定から読む実態に直した。

### OV-39: Session の source を読み直せないとき、編集をすべて捨てて Load が成功する

- 再現（2026-10-06 の調査。#833）: catalog を持たない Session（31〜33）で、crop・逆相補した record の feature の編集を持つものを、
  `pyodide.asm.wasm` と `pyodide.asm.js` の request を abort して Load する（`tests/web/feature-identity-overrides.playwright.spec.js` の Session 33 の変形）。
  16 workers で WORKER_INIT のエラーが出たときに、このテストが 1 件ではなく 5 件の脱落を数えたことが手がかりだった（「未確認の観察」の 3 つ目）。
- 期待: runtime が起動しないなら Load は失敗し、以前の Session が残る（PD-OI-045: 省略された data は受け入れた残余リスクではない）。
- 実際: "Session loaded successfully! 5 feature edit(s) from an older Session could not be matched to a feature of its saved diagram and were dropped."。
  編集の row は 0。次の Save で編集が失われる。
- 原因: RC-16。`extractSessionSourceFeatures`（`app/session-feature-metadata.js`）が `catch { return null; }` で runtime の失敗も握りつぶし、
  保存した metadata（描いた hash だけ）にフォールバックして、crop・逆相補した record の編集をすべて不一致と数えた。
- 修正: #833（2026-10-06 に merge。Owner の決定 2026-10-06: 推奨）。再試行で成功しうる失敗（`WORKER_INIT`、`RESULT_INVALID`、`UNKNOWN`、Worker の cancel など。
  `liveEditFailure` と同じ分類を `retryCanSucceed` に切り出した）は rethrow し、Load を `WORKER_INIT` の診断で失敗させて以前の Session を残す。
  source がない場合とデータの誤りは今までどおり `null`。runtime が起動できるときに同じ file を読み直すと、編集を migrate する。

### OV-40: Session 31〜39 の表に余分な cell があると、Load が表も Session も名指さずに失敗する（OV-28 の副作用）

- 再現（2026-10-06 の調査。#835）: 新しい fixture `tests/fixtures/sessions/whitelist-tab-keyword.v39.gbdraw-session.json.gz`
  （`17e2c9de` の app で Save Session した download。HmmtDNA.gbk、Circular、Labels Out、Generate、Label filtering Whitelist、Pattern を `cytochrome<TAB>c oxidase`）を Load。
  保存された row は `CDS<TAB>product<TAB>cytochrome<TAB>c oxidase`。
- 期待: 古い writer が書いた Session は、その writer が書いたとおりに読める。失敗するなら、どの表のどの Session の問題かを示す。
- 実際: `TABLE_INVALID`（reason `FIELDS`）で "The table is invalid. Line 1. Check the required fields. Required columns: 3."。表も Session も名指さず、
  edit-table と retry の操作は当たらない。Session 40〜45 は読める。
- 原因: RC-16。OV-28（#804）の厳密な reader（2、3、2 列）が、Session 31〜39 の projection（`projectCanonicalSessionRequest`）でも Default colors、Label whitelist、Qualifier priority に走る。
  この世代の writer は cell を normalize せずに結合した（normalize は OV-23 で入った）。
- 修正: #835（2026-10-06 に merge。Owner の決定 2026-10-06: A）。Session 31〜39 の 3 つの表だけ、現在の writer の書き方で読む。余分な cell は最後の欄に空白 1 つで結合し、必須の欄がない row は捨て、
  Load の通知が表名と行番号を報告する。残る表の失敗は表名と Session を名指し、操作は `select-input` だけ。Session 40〜45、表の file import（厳密のまま）、CLI replay（拒否のまま）は変えない。

### OV-41: Session を読み込むと、Record が "(not found)" と表示される

- 再現（OV-38/OV-40 の調査中に見つけた。#852 で Session 45 でも再現した）: Circular、Multi-Record Canvas を切る。
  `tests/fixtures/sessions/cli-web-mito.gb` を読み込み、record `NC_012920.1` を選んで Generate し、Session を保存する。
  page を読み直してその Session を Load する。Session の版には依らない。
- 期待: まだ検査していない record は、存在しないとは出さない。
- 実際: Record が "NC_012920.1 (not found)" と、"Records not inspected" の下に出る。record は存在する。
- 原因: Session の Load は source の record を検査しないので、discovery の状態は `deferred` のままで、catalog は空。
  保存された selector は `missing` に解決され、`app/app-setup.js` の `circularRecordPresentationOptions` は状態を見ずに
  "(not found)" を付けていた。隣の `circularRecordPresentationError` は `ready` を待っていた。
- 修正（#852）: 状態が `ready` でないあいだは "(not inspected)" と出す。"Inspect source records" か Generate の後は
  本当の行（"NC_012920.1 (16,569 bp)"）に替わる。本当にない selector は、検査の後に "(not found)" と出る。
  Load のたびに自動で検査することはしなかった（Worker の読み込みが毎回かかる。別の判断）。

### OV-42: Auto Reflow を切ると、feature の visibility の編集で legend が変わらない

- 再現（live と Generate の parity の調査、`b7e0426e`。#839）: `tests/fixtures/forced_label_underlay.gb`、Circular、labels out、Auto Reflow を切って Generate（legend は CDS、repeat_region、GC …）。
  規則 `FORCEDLBL` / `CDS` / `locus_tag` / `^fl1$` / Off を足す。
- 期待: Generate と CLI の順 repeat_region、CDS（Python は型の行を最初に描いた feature の順に並べる。FL1 が最初の CDS だった）。
- 実際: live の legend は変わらない。すべての CDS を隠すと、Generate は CDS の行を落とし、live は残す。feature 個別の Off と Redo、Circular と Linear でも同じ。
  Auto Reflow を入れると正しい（rerender が legend を導き直す）。
- 原因: RC-15。`projectFeatureVisibility` と `editFeatureVisibilityRules` は feature と label を隠すだけで、legend は変えない。
- Owner の決定（2026-10-06）: A。提示した選択肢は、A 自動 rerender を要求する（推奨）、B JavaScript で導く、C Python の legend だけの操作。
- 修正: #857。edit が legend の導く内容を変えるとき、Auto Reflow を切っていても自動 rerender を要求する。legend の導出は Python だけ。
  #839 の `test.fail` の case がこの修正で外れる。

### OV-43: live の specific color 規則が、既定の caption を残して規則の行を末尾に足す

- 再現（#839）: 規則 `CDS` / `locus_tag` / `^FL1$` / `#e63946` / `alpha` を足す（popup の "This feature only" でも同じ）。
- 期待（Generate）: alpha、other proteins、repeat_region …。
- 実際（live）: CDS、repeat_region、GC …、alpha。使った規則を live で削除すると、Generate は "CDS" に戻すが live は "other proteins" のまま。Auto Reflow を入れても変わらない。
- 原因: RC-15。`commitSpecificRules`（`prepareFileLegendEntries`、`legend.apply`）が規則の行を足し、消すだけで、generated の行（残りの caption、Python の順）を導き直さない。
  reflow は live の legend の状態を再適用する。live equals Generate の別のテスト（G-A）が "This feature only" の fill で legend を飛ばしている（`{ legend: false }`）のが、この差を決定なしに隠している。
- Owner の決定（2026-10-06）: A（OV-42 と同じ）。
- 修正: #857。

### OV-44: batch の規則の legend の行が、もう一方の Result に届かない

- 再現（#839）: `tests/fixtures/web_batch_two_records.gb` の Circular の batch を Generate。TESTA_0001 の popup で fill `#c83366`、同じ product の範囲
  （規則 `CDS` / `product` / `^duplicate protein$`）。
- 期待: もう一方の Result（TESTB の feature も `#c83366` になる）にも "duplicate protein" の行（Generate と同じ）。
- 実際: Result 2 を表示すると "duplicate protein" の行がない。色を変えたあとは swatch が古い色のまま、削除したあとは行が残る。
- 原因: RC-15（R3）。規則の legend の行は、mount された legend にだけ書かれる。`applyEditorOperationsToMountedSvg` の `legendAdds` は直接の追加だけを扱い、規則の行を扱わない。
- Owner の決定（2026-10-06）: A（OV-42 と同じ）。
- 修正: #857。

### OV-45: batch の Generate が、2 回目の "This feature only" の色のあと `UNKNOWN` で失敗する

- 再現（#839 で見つけ、#845 で直した）: `tests/fixtures/web_batch_two_records.gb` を Circular、Record は All records (separate diagrams) で Generate。
  TESTA_0005 を fill `#c83366`、"This feature only"（規則 CDS / hash / `#c83366` / "codon start two" ができる）、Generate（成功）。
  もう一度 TESTA_0005 を fill `#123456`、"This feature only"、Generate。
- 期待: 成功し、TESTA が "codon start two" を `#123456` で描き、TESTB は変わらない。
- 実際: `UNKNOWN`（stage render）。前の Result は残る。
- 原因: RC-15。2 回目の fill が "codon start two" のただ 1 つの寄与を置き換えるので `legendColorOverrides["codon start two"]` も書く。
  `compilePlanBundle`（`app/candidate-render.js`）の `styledCaptions` の loop が、その fill を batch のすべての Result に `allowMissing: false` で足し、
  行を描かない TESTB で Result の admission が素の `Error`（"Sanitized SVG content is missing a Legend binding."）を投げ、`UNKNOWN` になった。
  stroke の枝は `renderedResultIndexes` で絞っていたが、同じ `allowMissing` を使っていた。
- 修正: #845（2026-10-06 に merge）。fill と stroke が `allowMissing` を Result ごとに決める。行の feature を 1 つも描かない Result は行がなくてよい。
  残る欠落の失敗は producer（`requireLegendEntries`）が `RESULT_INVALID`（stage `result-admission`）で出す。

### OV-46: 1 つの Result だけが描く generated の legend の行に Legend だけの色を付けると、Generate が失敗する

- 再現（OV-45 の修正中に見つけた。#845 の本文に記録）: Circular の batch（`web_batch_two_records.gb`）で、TESTA_0005 を fill `#c83366`、"This feature only"、Generate
  （TESTA の legend は "other proteins"、TESTB は "CDS"）。Legend editor で "other proteins" の色を `#00aa00` にして Generate。
- 期待: Generate が成功する。
- 実際: `RESULT_INVALID`（result-admission。OV-45 の前は `UNKNOWN`）。
- 原因: RC-15。この行は `featureIds` を持たず、どの direct fill もその caption を名指さないので、compiler がどの Result が描くかを判断できず、すべての Result に必須とする。
  直すには設計の選択が要る（compiler が Legend entry の元の Result を知る、たとえば `prepareCommitInput` に `selectedResultIndex` を足す、または catalog が Result ごとの Legend の行を持つ）。
  OV-44（legend の状態が Result ごと）と同じ種類。
- 状態: A（#857 の自動 rerender）では直らない（#857 の probe で確かめた）。修正は #871（Owner に委ねられた選択 A）:
  どの Result が描くか compile のときに分からない行は、各 Result では欠けてよいとし、Result の受け入れで、その Generate のどれか 1 つの Result が描くことを求める。
  他の選択肢は B（legend の edit の元の Result を記録し、そこで必須にする。永続の欄が増える）、C（必須にしない。どの Result も描かない古い override を検出できなくなる）。

### OV-47: batch で、Undo と "Sort by default" が Result ごとの Python の legend の順に戻らない

- 再現（#857 の確認中に見つけた。`dev` `9a83e8f8` で Generate の経路でも起きる）: Circular の batch（`web_batch_two_records.gb`）。Result 2 を表示して
  TESTB_0001 を fill、"This feature only"。Result 2 を表示したまま Generate。Result 2 で Legend の Sort Z-A。Result 1 を表示して Undo、または "Sort by default"。
- 期待: Result 1 は Python の順（`CDS, tRNA, GC content, GC skew (+), GC skew (-)`）、Result 2 も Python の順に戻る。
- 実際: Result 1 は `CDS` が最後になる。Result 2 は Undo の後も "Sort by default" の後も Python の順にならない。
- 原因: RC-15。既定の順（`originalLegendOrder`）は 1 つの一覧で、Generate や rerender のときに表示している Result から取る。
  Result 2 を表示しているときに作ると、Result 1 だけが描く行（`CDS`）が一覧から落ちる。Result ごとの Python の順は 1 つの一覧では表せない。
  #857 では、自動 rerender が Result 2 を表示したまま起きるので、Generate なしでもこの状態になる。
  `tests/web/multi-result-edit-matrix.playwright.spec.js` の B19 と "Sort by default" の case は、正しい期待のまま `test.fail`（OV-47）にした（#857）。#870 がその印を外す。
- 修正: #870（Result ごとの既定の順。Generate や rerender のとき、または Result を初めて表示したときに、編集した順を当てる前の Python の描画から取る。
  Session の形は変えない）。

### OV-48: 自動 rerender が Generate の実行中に始まると、Generate が何も出さずに終わる

- 再現（#871 と #877 の CI、`tests/web/gui-audit-20260930-editor.playwright.spec.js` の "a Feature Visibility rule that hides a feature hides its label too"）:
  Auto Reflow を切り、Feature Visibility の規則を続けて 4 回編集（Undo、Redo、show、off）してから Generate。
- 期待: Generate が新しい Result を出す。
- 実際: 180 秒待っても `resultGenerationKey` が増えず、`processing` は false、エラーもない（タイミングによって起きる）。
- 原因: RC-17。規則の編集が要求した rerender が順に回り、Generate が番号を取った後に始まった回が番号を進め、Generate が `stale` で終わる。
  #857 で、規則と visibility の編集が Auto Reflow を切っていても rerender を要求するようになり、起きやすくなった。
- 修正: #879。`processing` の間は rerender を始めず、要求は Generate が終わった後に 1 回だけ実行する。Generate が実行中の rerender を追い越すのは今までどおり。
  同じ spec の `:888`（Label On の後に編集が見えない）はテストが rerender の終わる前に読んでいただけで、製品の不具合ではなかった（OV-48b、#879 でテストを直した）。

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
- **live の色規則の caption の rebind（OV-37 の候補）**: 不具合ではなかった。詳細は OV-37。
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

## 修正の状況（2026-10-05〜06）

Q4 の設計文書（`DESIGN-Q4-FEATURE-IDENTITY.md`、#763）は、PR-5 と DOC-1 の後の実装を PR-Q4-1〜5 に分けた。

| ID | 状態 | PR |
|---|---|---|
| OV-01 | 修正。CLI と Source recipe は `--feature_override_table`、Web は identity の行（request schema 9、Session 45） | #766、#771、#779、#784 |
| OV-02 | 修正。live の照合を描画した座標の selector（feature catalog 5）で行う | #784 |
| OV-03 | 修正。選択から作る annotation の対象を identity にする（PR-Q4-5） | #788 |
| OV-04 | 修正。label と visibility は複製ごと（設計の決定 Q1 = A）。色の「This feature only」は PD-OI-069 のとおり複製で共有 | #771、#784 |
| OV-05 | 修正。色は #760、失敗の分類は #761、label と visibility は identity の行 | #760、#761、#784 |
| OV-06、OV-07 | 修正。Apply の時点のダイアログ（Q1）。ダイアログを通らない経路は `LABEL_NOT_DRAWN` | #761、#764 |
| OV-08 | 修正。request には現在の mode の record の placement だけを載せる | #762 |
| OV-09 | 修正。Python が `FEATURE_PLACEMENT` の診断を出し、設定を変える時点でダイアログ（Q3） | #762、#768 |
| OV-10 | 修正。Apply の時点のダイアログ（Q2）と popup の案内 | #764 |
| OV-11 | 修正。CLI と Source recipe は #779、Web は #784 | #779、#784 |
| OV-12 | 修正 | #784 |
| OV-13 | 修正。4 owner を coverage に加えた | #761 |
| OV-14 | 修正 | #757 |
| OV-15 | 修正 | #781 |
| OV-16 | 修正（潜在の欠陥） | #785 |
| F | `reason` の受け渡しを削除した | #759 |

規則の変更: #771 の承認の手順をきっかけに、Owner が ARCHITECTURE_EXCEPTION の承認を普通の approval にした（#778）。

### 2026-10-05: 残る課題

R-1〜R-7 の状態。原因と修正の一行は「残る課題」に書いた。

| ID | 状態 | PR |
|---|---|---|
| R-1 | 修正。draft の行が `scope`（mode）を持つ。2026-10-08 の #945 で mode ごとの drawing に置き換えた（`gbdraw/web/CLAUDE.md` R2、R15） | #796（annotation 側は OV-21: #806）、#945 |
| R-2 | 修正。「描かれるか」を Generate と同じ規則で答える | #795 |
| R-3 | 修正。slot を変える編集を 1 つの遷移に通す | #805 |
| R-4 | 修正。ダイアログの文を理由の表から作る | #799 |
| R-5 | 修正。隠れた feature を Features 一覧と Search から届くようにする | #807 |
| R-6 | 修正。GFF3 を 1 回だけ解析する | #793 |
| R-7 | 修正。変形なしの record の `hash=` だけ移す | #811 |

新しく見つけた OV-17〜OV-48 の状態（OV-35 以降は 2026-10-06）。

| ID | 状態 | PR |
|---|---|---|
| OV-17 | 修正（runtime）。Gallery の更新は #851。昇格が残る | #794、#851 |
| OV-18 | 修正（runtime）。Gallery の更新は #851。昇格が残る | #797、#851 |
| OV-19 | 修正。#810 は merge せずに閉じ、#824（E1）に載せ直した | #824（#810 は閉じた） |
| OV-20 | 修正 | #798 |
| OV-21 | 修正 | #806 |
| OV-22 | 修正（Owner の決定 A） | #800 |
| OV-23 | 修正 | #801 |
| OV-24 | 修正（Owner の決定 A'） | #802 |
| OV-25 | 修正 | #803 |
| OV-26 | 修正（Owner の決定: 推奨） | #812 |
| OV-27 | 修正（Owner の決定: 推奨） | #812 |
| OV-28 | 修正 | #804 |
| OV-29 | 修正（#805 の中） | #805 |
| OV-30 | 修正 | #809 |
| OV-31 | 修正（#807 の中） | #807 |
| OV-32 | 修正（#807 の中） | #807 |
| OV-33 | 修正 | #808 |
| OV-34 | 修正。#824 が、描かれていない feature を描かせる規則の編集で自動 rerender を要求する | #824 |
| OV-35 | 修正（#824 の中）。Python が結んだ label も投影する | #824 |
| OV-36 | 修正。Generate の commit で note を消す | #838 |
| OV-37 | 調べた。不具合ではない | なし |
| OV-38 | 修正（Owner の決定: 推奨）。slots が off の null の `spacing` を読み、他は欄と row を名指して拒否する | #834 |
| OV-39 | 修正（Owner の決定: 推奨）。runtime の失敗は Load を失敗させる | #833 |
| OV-40 | 修正（Owner の決定: A）。Session 31〜39 の 3 つの表を現在の writer の書き方で読む | #835 |
| OV-41 | 修正。検査する前は "(not inspected)" と出す | #852 |
| OV-42 | 修正（Owner の決定 A）。legend の導く内容が変わる edit で自動 rerender を要求する | #857 |
| OV-43 | 修正（Owner の決定 A）。legend の導く内容が変わる edit で自動 rerender を要求する | #857 |
| OV-44 | 修正（Owner の決定 A）。legend の導く内容が変わる edit で自動 rerender を要求する | #857 |
| OV-45 | 修正。per-Result の `allowMissing` と `RESULT_INVALID` | #845（見つけたテストは #839） |
| OV-46 | 修正 PR。描く Result が分からない行は、どれか 1 つの Result が描けばよい（選択 A） | #871 |
| OV-47 | 修正。Result ごとの既定の順 | #870（`test.fail` は #857） |
| OV-48 | 修正。rerender は Generate の実行中に始まらず、終わった後に 1 回だけ走る | #879 |

OV-42〜OV-45 を見つけた live と Generate の parity のテスト（17 の case と 2 の表の確認）は #839 で入った。OV-42、OV-43、OV-44 の case は修正まで `test.fail`、
OV-45 の case は #845 で外れた。

未確認の観察（原因を調べていない）:

- `docs/images/t-cli-11/interactive_human_mitochondrion.interactive.svg` は前から古い（`run_cli_scenarios.py --scenario T-CLI-11 --check` が `cedbe18d` で失敗する。どのテストも見ていない）。Gallery の更新（#851）で作り直した。
- `tests/web/right-drawer.playwright.spec.js` の "preview similarity-group copy actions report isolated accessible outcomes" が、この機械で `page.clock.install()` の後の feature popup を開く行で timeout する。`origin/dev` でも同じ。CI は通る。
  → 調べた（2026-10-06）: テストの不具合で、clock ではない。Gallery の Session を読み込んだあとの最初の popup が diagram の Worker で Pyodide を起動し（9.7 秒のうち約 5.4 秒）、
  既定の 30 秒を並列の負荷で使い切った。#829 が test の timeout を延ばした。`:1215` と `:1284` は既定のまま。
- `feature-identity-overrides` の Session 33 のテストが、16 workers で `WORKER_INIT` のエラーが出ると 1 件ではなく 5 件の編集の脱落を数える。3 workers では通る。負荷の問題と見ている。
  → 調べた（2026-10-06）: 製品側の抜けで、数え方の競合ではない。OV-39（#833）。
- `tests/web/feature-color-actions.test.mjs` の 699 行目の 1 回だけの失敗 → 調べた（2026-10-06）: 環境の問題。負荷で Python の helper の subprocess が 1 つ失敗した。

### Owner に委ねられた選択（2026-10-05）

標準の指示（2026-09-29）で、推奨の選択肢を採った。PR ごとに 1 行。

- #793（R-6）: 入れ子のまま解析して既存の `filter_features_by_type` をメモリ上でかける（すべての type を平らに読む案は、よくある場合で約 9% 遅かった）。Web の解析結果の cache の key は type を含まない。
- #794（OV-17）: catalog の producer ではなく runtime の endpoint の確認を直す。複数の claim は、catalog の展開が全部解決していれば信頼する。Gallery の SVG はこの PR では作り直さない。
- #795（R-2）: visibility の規則の照合は、JavaScript に regex を移さず、既存の `evaluateRules` の新しい kind で Python に頼む。catalog 5 は変えず、`hash`、`location`、`record_location` の規則を描かれていない feature に当てるときの答えは「不明」（ダイアログを出さない）。
- #796（R-1）: 行に `scope` を持たせる。record key を mode ごとに別にしない（Session の保存形式と Gallery Session の変更が要る）。Session 41〜44 の Main の placement の行は両方の mode に残す。
- #797（OV-18）: 複数 anchor の行の分け方も同じ PR で直す。ガードは、standalone の runtime が DOM を要るので、node ではなく 1 つの browser のテストで Preview と比べる。
- #798（OV-20）: 共有の helper を `gbdraw/io/` に 1 つ置く。インデントした `#` の行もコメント。引用符の扱いはこの PR では変えない（OV-22）。
- #799（R-4）: scope と filter の理由は Default の label visibility にだけ当てる（On は両方に勝つ）。`scope_orthogroup_top` は設定の名前を挙げ、この feature が外れたとは言わない。
- #800（OV-22）: `read_literal_table` を 1 つの helper にする。annotation の表は決定を待つので変えない。
- #801（OV-23）: Qualifier priority の書き出しも同じ PR で直す。空の whitelist は stage しない。
- #802（OV-24）: pandas の helper ではなく標準の `csv.reader`（`QUOTE_NONE`）にして、#800 に依存しない。
- #803（OV-25）: 表の名前を `label` で渡し、path と stream の両方に同じ列数の検査をかける。コメントの扱いはこの PR では変えない。
- #804（OV-28）: 新しい reason を足さず、`TABLE_INVALID` の reason `FIELDS` を使う。短い行も拒否する。
- #805（R-3）: Undo、Redo、Session の読み込み、Import、Reset Settings は確認しない。Depth と Comparison の行などの reconcile も確認しない。他の mode の lane を失うときは確認し、その mode を名指す。
- #806（OV-21）: `scope` のない `featureIdentity` の対象は、保持せず `INPUT_INVALID` で拒否する。request だけの Session の対象は request の mode を使う。
- #807（R-5）: チェックボックスは効いている表示状態を示し、操作は常に On か Off を書く（Default に戻すのは popup）。type は committed request のもの。規則の照合が不明なときは Result の catalog が描くかどうかを使う（catalog 5 は変えない）。描かれていない feature を表示にすると、Auto Reflow を切っていても自動 rerender を入れる（右の editor の live-edit の規則）。
- #809（OV-30）: Editor の toggle の安定した name は "Editor"（状態で変わる name は R12 違反）。`aria-label` は静的な title と同じ。
- #810（OV-19）: Generate が拒否する regex は draft に残し Result は変えず、Generate と同じエラーを出す。規則を書いてから、Python の照合が届いて投影する（History の apply と同じ）。Reset Settings は遷移の外のまま。#810 は閉じ、#824 が同じ挙動を引き継いだ。
- #811（R-7）: 移した対象の `scope` は Session の mode。Session 31〜39 は catalog がないので移さない。rotation も crop、逆相補と同じく変形として扱う。record を確実に結べない対象は移さない。
- #812（OV-26/27）: Qualifier priority にも同じコメントの規則を適用する（同じ helper で、Web は飛ばしている。Owner が確認した）。annotation の表も `read_table_lines` を `csv.reader` の前にかける。

### Owner に委ねられた選択（2026-10-06）

標準の指示（2026-09-29）で、推奨の選択肢を採った。PR 本文に書かれたものを 1 行ずつ。

- #824（OV-35、OV-34、OV-19）: #810 の挙動を #807 の rerender の上に載せて引き継ぐ。Generate が拒否する規則の表がある間は、History の apply と Result の表示でも何も投影しない（#810 から）。
  label の投影は、Python が結んだ label も読む（Python が描く feature の label だけが変わる。resolver が不明と答えたら残す）。Load Feature Edits TSV の label の配置は `reflow` の option で保つ。
- #838（OV-36）: 成功した Generate の commit だけ note を消し、失敗、cancel、supersede された Generate は残す。note は History の状態にしない（Undo で戻らない）。
  record の回転と Similarity alignment の実行は note を残したまま（Generate ではないので変えない）。
- #834（OV-38）: slots が off で null の `spacing` だけ読み、他の退役した欄は拒否して欄と row を名指す。`session-active-config-contract.js` は R6 の走査対象に加えない（baseline の追加になる）。
- #833（OV-39）: 再試行で成功しうる失敗だけ Load を失敗させ、データの誤りは `null` のまま。分類は `liveEditFailure` と同じ判定を `retryCanSucceed` として export して共有する。
  Session を読み込んで「migrate しなかった」と報告する案は採らなかった。
- #835（OV-40）: 余分な cell は最後の欄に連結し、Default colors では色の欄に連結する（key に tab を持つ行は色が無効で捨てて報告する）。Web の Load と CLI replay が食い違う
  （CLI は今までどおり拒否する）。Owner の決定は Web の Load だけ。
- #839（parity のテスト）: 許す差は 1 つだけ（Auto Reflow を切ったときの、両側に描く label の leader の位置は reflow を待ってよい。Owner の OV-35 の指示）。
  OV-43 の legend の差は決定がないので許さない。
- #845（OV-45）: 行を必須にするのは、その行の feature を描く Result だけにする。executor は変えない。OV-46 は設計の選択が要るので直さず、別の finding にした。

## 残る課題

修正 PR のレビューと検証で見つけ、それぞれの PR の範囲外にしたもの。2026-10-05 に R-1〜R-7 を直した（状態は「修正の状況」の表）。

- **R-1 mode をまたぐ record key（#762）→ 解決、#796**: Circular と Linear が同じ record key（Gallery Session の `record-1` など）と同じ feature を持つと、
  Main の placement の行が両方の mode に当たる。原因: draft の key が mode を持たない。placement だけでなく feature override も同じだった。
  行に `scope` を持たせ、key を [scope, recordKey, biologicalFeatureId] にした。annotation の対象は OV-21（#806）。
- **R-4 Label Not Shown の最初の文（#764）→ 解決、#799**: Show Labels と filter のことしか書いていなかった。
  原因: 文が判定と別に書かれていた。`labelAbsenceReason` と理由の表 `LABEL_ABSENCE_REASONS` から、ダイアログと popup の文を作る。
- **R-2 TSV の visibility 規則で隠した feature（#764）→ 解決、#795**: 読み込んだ規則で隠れた feature の Label On が、Q2 のダイアログを出さなかった。
  原因: 「描かれるか」の判定が 2 か所にあった。`resolveFeatureDrawn` が `should_render_feature` と同じ規則（規則の照合は Python）で答える。
- **R-5 Keep feature hidden の後（#764、前から）→ 解決、#807**: 隠れた feature の popup は開けない。
  原因: 一覧が catalog の描かれた feature（`extractedFeatures`）だけを並べた。Owner の決定（2026-10-05）どおり、右の Features 一覧と Search で隠れた feature に届く。
- **R-3 ダイアログを通らない設定変更（#768）→ 解決、#805**: Use custom stack、preset の Reset、feature の行の削除・無効化・並べ替え、
  Circular の panel の Separate Strands は Q3 のダイアログを出さなかった。原因: 各 control が slot を別々に書いた。
  `changeTrackLayout(apply, control)` の 1 つの遷移で適用し、lane を失うときにダイアログを出す。
- **Python の Session と crop（#771、#786 で解消）**: Python の Session の writer は crop した record を crop 後の GenBank として保存していたので、
  crop した Python の request から保存した identity の行と placement は、描き直すときに解決しなかった。別のセッションの #786 が、
  入力ファイルと crop・向きを保存するようにした。`dev` `b7911c0f` で確認: `tests/fixtures/b_collide.gb` を Linear で
  `--region TESTA:201-3800 --reverse_complement 1` と `--feature_override_table`（Off 1 行、label 文字 1 行）で描き、`--session_output` の
  Session を CLI で描き直すと、SVG が一致し、2 つの編集が残る。
- **R-6 GFF3 の 2 回目の読み込み（#771）→ 解決、#793**: identity の行が読み込んでいない type の feature を表示にすると、GFF3 をもう一度読んだ。
  原因: `_load_request_records` が type の絞り込みと identity の解決で別々に解析した。入れ子のまま 1 回解析し、type の絞り込みをメモリ上で行う（解析 2 回 → 1 回）。
- **Q4-4 の後に続くもの → 解決、#792**: Web の "Export / Load feature edits TSV"（`--feature_override_table` の形式）は #789 で足した。
  Gallery Session は #790 で Session 45 に作り直した。#784 の後、`tools/refresh_gallery_sessions.py` が catalog 5 の `drawnSelector` を
  生物学的な欄の重複として拒否していたので、#790 で直した。architecture ratchet の CW scope に残っていた
  "feature-selector metadata and uniqueness index" の行は、#792（2026-10-05）が消した（解決）。
- **R-7 古い Session の annotation（#788）→ 解決、#811**: v45 より前の Session の `hash=` で指す annotation は identity に移さなかった。選択から作った行と
  手で書いた行を区別できないため。Owner の決定（2026-10-05）どおり、crop・逆相補・rotation のない record にあり、catalog のちょうど 1 つの feature に
  当たる行だけを Load の時点で `featureIdentity` に移す（読み込んだ図は変わらず、移した件数を通知する）。それ以外は今までどおり残る。

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

### 2026-10-05

Owner が 2026-10-05 に、残る課題の修正の途中で決めた。

| 課題 | 決定 | 内容 |
|---|---|---|
| R-5 | 隠れた feature は Features 一覧と Search から届く | 選んだ type の feature は、個別の Off や規則で隠れていても一覧に残し、Visibility のチェックボックスで効いている表示状態を示す（On にすると個別の On、外すと個別の Off。個別の On は type と規則に勝つ）。個別の設定を持つ feature（On から Off に戻したものも）と、描かれている feature も載せる。Search は隠れた feature も見つけ、popup は「この feature は隠れている」と示す。新しい Hidden features の一覧は足さない。 |
| R-7 | 変形なしの record だけ移す | v45 より前の Session の `hash=` の featureSpan の対象のうち、crop も逆相補もしていない record にあり、source catalog のちょうど 1 つの feature に当たる行だけを `featureIdentity` に移す。読み込んだ図は変えず、移した件数を通知する。それ以外は今のまま。 |
| OV-22 | A | styling の表は引用を切って読む。`"` は値の一部。文書に書く。 |
| OV-24 | A' | annotation の表の引用符も値の一部（OV-22 と同じ。CSV 引用を固定していたテストを書き直す）。 |
| OV-26、OV-27 | 推奨 | Python が annotation、Specific colors、Default colors の表の行頭の `#` の行を飛ばす（`read_table_lines`）。文書に書く。Qualifier priority も同じ規則に含める（確認済み）。 |

### 2026-10-06

Owner が 2026-10-06 に、Phase E と残る観察の調査の途中で決めた。

| 課題 | 決定 | 内容 |
|---|---|---|
| OV-38 | 推奨どおり修正 | Web app は 0.13.0 の GUI Session を、Custom Track Slots が off で退役した欄が null のとき（損失なしに）読む。非 null の値と slots が on の Session は、欄と row を名指して拒否する（#834）。 |
| OV-39 | 推奨どおり修正 | `null` を返すのは source がないときとデータの誤りだけ。runtime の失敗は Load を `WORKER_INIT` の診断で失敗させ、以前の Session を残す（#833）。 |
| OV-40 | A | Session 31〜39 だけ、3 つの表を現在の writer の書き方で読む。余分な cell は最後の欄に空白 1 つで連結し、必須の欄がない row は捨て、Load の通知に表名と行番号を書く。残る表のエラーは表名を名指す（#835）。他の選択肢は B（拒否したままで表名と Session を名指す）、C（古い lenient な reader を restore だけに使う。keyword を黙って切り、reader が 2 つになる）。 |
| OV-42、OV-43、OV-44 | A | edit が legend の導く内容を変えるときは、Auto Reflow を切っていても自動 rerender を要求する。legend の導出は Python だけ。選択肢は B（JavaScript で導く）、C（Python の legend だけの操作）。修正は #857。 |

OV-46 は OV-44 と同じ種類なので A に照らして確かめた。A では直らない（#857）。修正は、Owner に委ねられた選択 A（どれか 1 つの Result が描けばよい。#871）。

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
