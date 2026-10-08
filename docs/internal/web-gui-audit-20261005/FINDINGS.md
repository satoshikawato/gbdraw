# Web GUI 再監査 — dev 8e5234ac（2026-10-05）

対象は `origin/dev` の `8e5234ac`（PR #796 のマージ時点）。読み取り専用の worktree
`/home/kawato/gbdraw-work/.worktrees/gui-audit-20261005` で監査した。コードの変更、commit、issue 作成はしていない。

**固有の不具合は 59 件（P1 1、P2 14、P3 44）。** すべて Chromium（Python Playwright）で dev を動かして再現した。
以前の監査で直した項目（2026-09-30 の 65 件、OV-01〜OV-18、OV-20）に退行はない。ただし IN-04 と CO-10 は修正が一部にとどまっている（CI-03、CI-05）。

## 方法と表記

- 5 領域を並列に監査した: TK カスタムトラック、UJ ユーザージャーニー、UI UI/UX/アクセシビリティ、CI 入力と Linear の比較、FL 色・検索・注釈・凡例・出力。
  各領域の報告は `reports/<area>.md`、probe・ログ・スクリーンショットは `probes/<area>/` にある。
- 実行中の override-residuals セッション（gbdraw-51）が扱っている項目は対象から外した:
  - R-3、R-5、R-7
  - OV-19、OV-21、OV-22、OV-24、OV-25
  - PR #799〜#801
- 次の P1/P2 は、担当とは別に orchestrator が再確認した:
  - TK-01: CLI で再現
  - UJ-01、FL-06: probe を再実行して再現
  - UI-01: CSSOM とスクリーンショットで確認
  - TK-04、CI-02: コードで確認
- 重大度は 2026-09-30 の監査と同じ。
  - **P1**: 誤った科学的出力、編集の消失、有効な入力で Generate できない、CLI も同じく誤る。
  - **P2**: 機能が壊れている、または Generate まで古い状態が残る。回避策はある。
  - **P3**: 軽微、端のケース、表示、アクセシビリティ。
- JS のパスは `gbdraw/web/js/` からの相対パス。それ以外はリポジトリ直下から。

## P1（1 件）

| ID | 症状 | 原因 |
|---|---|---|
| TK-01 | Circular の custom stack で、gc_content や gc_skew の行に **Radius** を入れると、Generate が TRACK_LAYOUT で失敗する。入れた値がその行の "(auto)" と同じ値（0.6 R）でも、preset のどの track type でも、Gallery の `HmmtDNA_ATskew` Session でも失敗する（例: "Available band: 275–193 px"）。CLI も `--circular_track_slot gc_content:dinucleotide_content@side=inside,r=0.6`（features と ticks の行も指定）で "slot 'features' cannot fit inside between 274.9px and 193.0px" となる。Radius を外せば成功する。結果として、GC の行の半径は指定できない | `gbdraw/diagrams/circular/radial_layout.py` の `_resolve_circular_radial_layout`。`inside_max_outer` が、順序に関係なくすべての pinned inside slot の内縁で制限される。さらに `_inside_placement_window` が、次の pinned slot から内側の限界を取る。そのため、上にある unpinned の行に割り当てる範囲が誤る |

## P2（14 件）

### カスタムトラック

| ID | 症状 | 原因 |
|---|---|---|
| TK-02 | radius の重なり（gc_content と gc_skew をどちらも 0.5 にする、中央の definition と重なる 0.35 にする）が、"diagram engine failed … Python exception: ValidationError" と Save Session の案内になる。CLI は重なる slot の名前を出す | `radial_layout.py` が `diagnostic=` を付けずに `ValidationError` を投げる（約 830、1594/1600、1920 行）。そのため `gbdraw/web_support/error_adapter.py` の `_classify_native` が RENDER_FAILED にする |
| TK-03 | Depth の track index を一時的に範囲外にすると（↑ を 1 回押して ↓ で戻すだけでよい）、見えない Depth series が残る。その後の Generate と Save Session が両方 UNKNOWN で失敗する。2 回 Undo すると戻る | `app/app-setup.js` で、描画中の getter（`getDepthTrackLegendLabelForSlot` → `ensureDepthTrackEditableConfigCount`）が `adv.depth_tracks` を伸ばす。さらに `services/session-request.js` `buildDepthResources` の "has no source" という文言が、`services/error-normalization.js` の正規表現（"no TSV source" を探す）に一致しない |
| TK-04 | 行の renderer を変えると、前の renderer の parameter（`tick_label_layout`、`nt`、色）が残る。Generate は "retains '…' from another renderer" で止まり、元に戻しても解除されない。parameter を消す UI はなく、行を削除するしかない。Linear では GC の行が Depth series の legend label まで引き継ぐ | `app/circular-track-slots.js` の `updateCircularTrackSlotRenderer` と `app/linear-track-slots.js` の `updateLinearTrackSlotRenderer` が `slot.params` をそのまま残し、`app/track-slot-validation.js` の `PARAM_OWNER` がそれを拒否する（コードで確認） |
| TK-05 | Axis の外側へ移した Depth の行がある状態で、Depth series の **Remove** を押すと、Features の行が axis の外側へ移る。次の Generate で、feature が目盛の外に描かれる | `app/app-setup.js` の `removeCircularDepthTrack` が `circular_track_slots_axis_index` を直さない。アップローダーの Remove と Linear の経路は直している |

### ユーザージャーニー、凡例、色

| ID | 症状 | 原因 |
|---|---|---|
| UJ-01 | Show Labels None（Circular の既定）で、ND1 にラベル文字を入れて **Show this label** を選び、Ctrl+Z を 2 回押す。state は編集前とまったく同じになるが、preview、`results[].content`、SVG の出力には、編集前の図になかった既定のラベル "NADH dehydrogenase subunit 1" が残る。次の Generate で消える | `app/feature-editor/label-actions.js` の `projectLabelIntent` → `applyLabelVisibilityPreview`（History の restore から `reconcileLabelOverrides` 経由で呼ばれる）。visibility が `default` のとき preview の属性を外すだけで、`on` の reflow が足したラベルを消さず、reflow も予約しない。Label Mode が Out なら正しく戻る |
| FL-02 | tRNA-Phe の Fill を変えて **Apply to all tRNA (22)** を選び、SPECIFIC RULES の **Clear All** を押して Generate する。feature は既定の金色に戻るが、凡例と出力の tRNA は青 `#0000ff` のまま。図と凡例が食い違い、この状態は Session にも保存される | `app/feature-editor/color-actions.js` が書く `legendColorOverrides` を、`app/feature-editor/rule-actions.js` の `clearAllSpecificRules` と `removeSpecificRule` が片付けない。`app/candidate-render.js` が Generate のたびにそれを当て直す |
| FL-03 | FL-02 の状態で、caption のある tRNA ルール（例 `tRNA product tRNA-Ser #ffee00 Ser`）を足すと、Generate が UNKNOWN で失敗する。新しい Session でも CLI でも、同じ表は描ける | 推定: 古い `tRNA` の凡例 override の当て先がなく、`candidate-render.js` の `legendFills`（allowMissing:false）が例外を投げる。その例外が分類されない |
| FL-06 | GC content または GC skew の凡例名を変え、Legend の Position を None にすると、Generate が UNKNOWN で失敗する。feature の項目の改名、削除、移動、Sort Z-A では起きない。Position を left に戻すと成功する | 推定: 凡例のない図に、GC の行の改名を当てようとしている（`candidate-render.js` の `legendRenames`） |
| FL-07 | GC content や GC skew の凡例名を変えて Generate すると、その行だけ 12 px 下にずれて描かれる（行の間隔が 53 px と 29 px になる）。出力にもそのまま残る。live では全行が +12 px ずれる | 推定: `app/legend/entry-actions.js` の `moveLegendEntryToAnchor` と `legendEntryAnchor` |

### 入力、比較、再現

| ID | 症状 | 原因 |
|---|---|---|
| CI-02 | Circular で Up/Down で変えた record の並び順が、Run Info の Source recipe に入らない。CLI は元の順に描く。`--session` の replay と Linear は正しい | `app/run-info.js` の records 表の builder が `order: index + 1` と書き、`column` も空にする。records 表を使うときは `--multi_record_position` も出さない（コードで確認） |
| CI-03 | GenBank ファイルを置き換えたり Remove したりしても、Circular の単一 record の crop、逆相補、Record label、Subtitle が残る。Generate は "record selection is invalid"（"REC_C.1 (not found)"）で失敗する。別の record を選ぶと、古い crop、逆相補、タイトルで黙って描かれる。Linear でも、単一 record のファイルを単一 record のファイルで置き換えると、region、逆相補、Definition が残る（IN-04 の修正は一部だけ） | Circular: `app/app-setup.js` の `setCircularRecordPresentationSelector` と source watcher が `form.circular_region_*`、`circular_reverse`、`circular_record_label/subtitle` を初期化しない。Linear: `setLinearSeqPrimaryFile` が `group.records.length > 1` のときだけ初期化する |
| CI-04 | query と subject を入れ替えた BLAST 表で、`COMPARISON_IDENTITY` の汎用文になる（context は空。pair も ID も出ない）。CLI は "Swap the query and subject columns…" と原因を示す。表のほかの誤りも、Circular conservation 用の "Supply a comparison sequence file (FASTA, GenBank, or DDBJ)…" の文になり、pair、ファイル名、値を出さない | `gbdraw/web_support/error_adapter.py` の `_classify_native` と、`services/error-normalization.js` の COMPARISON_IDENTITY が pair と ID を落とす |
| CI-06 | pair に割り当てた BLAST ファイルを Selected pairs の Remove で外して Generate すると、`UNKNOWN`（stage request-validation）になる | plan の edge が `{source:"upload", file:null}` のまま残り、検証の例外が分類されない（原因の関数は未特定） |

### 見た目

| ID | 症状 | 原因 |
|---|---|---|
| UI-01 | `index.html` の `@apply` のルールがすべて効いていない。対象は `.card`、`.input-label`、`.form-input`、`.form-checkbox`、`.btn*`、`.upload-zone`、`details > summary`。入力欄 213 個に枠も角丸もない。`btn-secondary` 約 100 個（"Retry Generate" など）は地の文に見え、アップロード欄の点線も出ない。CSSOM では `.card { }` のように空になっている（`probes/ui/81-as-served.png` と `82-with-typed-style-experiment.png` を比較） | `gbdraw/web/index.html:26` は普通の `<style>` だが、Tailwind Play は `style[type="text/tailwindcss"]` しか処理しない。`@apply` は最初の commit `1ff0c843`（2025-12-17、当時は `cdn.tailwindcss.com`）からあり、一度も効いていない。**直すと全画面の見た目と幅が変わるので、Owner の判断が要る**（下記） |

## P3（44 件）

### カスタムトラック（TK）

| ID | 症状 | 原因 |
|---|---|---|
| TK-06 | 先頭の Depth series のファイルを外すと "Supply the required value" になり、"or remove the series" という回復策が消える | `services/error-normalization.js`（約 361 行）が REQUIRED に分類する |
| TK-07 | Circular の Depth track index に `1.5` や `-1` を入れても、黙って `1` に戻る（X-02 の類） | `circular-track-slots.js` の `normalizeTrackIndex`、`index.html` の `v-model.number` |
| TK-08 | slot の Depth legend label の help は「この slot の文字」と言うが、実際は series の名前が変わる | `app-setup.js` の `setDepthTrackLegendLabelForSlot` → `setDepthTrackLabel` |
| TK-09 | Hide GC Skew が AT skew の行も隠し、ダイアログでそれを "GC skew tracks" と数える。文法の誤り "those custom track slot" もある | `SUPPRESS_RENDERER_BY_KEY` が `nt` を見ず、renderer で選ぶ |
| TK-10 | Show Depth が OFF だと、stack の Reset が Depth の行を落とす。docs には "the loaded Depth series" から作り直すとある | `resetCircularTrackSlotsToPreset` |
| TK-11 | 行の Inner gap の入力が Outer gap の下にもぐり、Outer gap と legend label が行のカードからはみ出して切れる | `.track-slot-geometry-input` に `width:100%` と `min-width:0` がない（UI-01 を直した後に確認し直す） |
| TK-12 | Width/Radius に `0x10` を入れると 16 として受け付ける（CLI は拒否する）。`0` では同じ意味のエラーが専門用語で 2 つ出る | `track-slot-validation.js` の `parseOptionalCircularScalar` が `Number()` を使う |
| TK-13 | Features の行のない stack を、"feature underlays" を理由に拒否する。record に underlay の feature がなくても拒否し、CLI は同じ stack を描ける | `track-slot-validation.js`（約 1056 行）が、record ではなく設定した type を見る |
| TK-14 | docs には Linear の Depth TSV disclosure が開いた状態で始まるとあるが、閉じている（`25dd141f` で `open` が外れた） | `index.html` または `docs/REFERENCE/web-app.md:79` |
| TK-15 | Generate 前の "(auto)" は概算なのに、確定値と同じ形で出る。値も描画と違う（ticks が 0 px / 0.98 R と出るが、描画は 7.8 px / 0.77 R） | `track-slot-display.js` と `estimateCircularSlotGeometry` |

### ユーザージャーニー（UJ）

| ID | 症状 | 原因 |
|---|---|---|
| UJ-02 | mode を切り替えると、preview に別の mode の Result が残る。feature をクリックしても何も起きず、検索は 0/0、SVG ボタンはその古い図を出力する。sidebar には "…has changes pending Generate" と誤って出る | `record-display-options.js` の `hasPendingChanges` が、mode が違えば常に true を返す |
| UJ-03（+FL-01） | "This feature only" の色や凡例名を変えると、live の凡例と Generate の凡例が食い違う。live は項目が末尾に足され、残りの CDS は `CDS` のまま。Generate は `<新しい名前>, other proteins` になる。Reset fill の後も同じ（`other tRNAs` と `tRNA`、並び順）。出力される凡例が、編集の経路で変わる | live の凡例の投影と Python の分割が一致しない。docs は "on regeneration" としているので、仕様として受け入れるかを判断する |
| UJ-04 | Gallery の Session 10 件のうち 7 件で Output Prefix が `out` になっている。Generate 後のファイル名は `out.*` で、カードのコマンド（`-o <id>`）と合わない | Gallery の publication / refresh の経路（`tools/refresh_gallery_sessions.py` など）。CLI で書いた Session は prefix を保つ |
| UJ-05 | Linear で pair が 1 つのとき、**Upload BLAST TSV** で表をアップロードすると、ボタンが押されていない状態に戻り、CUSTOM のバッジが出る | `comparison-ui.js` の `intentKeyForPlan` が SELECTED を custom に対応させる |
| UJ-06 | 初めて開いたとき、例のデータへの道がない。Gallery は Session のダウンロードしか提供せず、Load Session を使う必要がある | UX の設計（Open in app がない） |
| UJ-07 | GenBank 欄に FASTA、テキスト、SVG、空のファイルを入れると、どれも "No records were found." になる。GenBank 形式でないことも、FASTA なら GFF3 + FASTA を使うことも案内しない | record の検査が形式の誤りを区別しない |
| UJ-08 | 正しいファイルに置き換えても、Generation Error の帯が次の Generate まで残る | ファイルの置き換えで `errorLog` を消さない |
| UJ-09 | Load Session は確認なしで未保存の作業を捨て、History も消す（9 → 0）。誤って Load しても Undo できない | 設計かどうかの判断が要る |
| UJ-10 | 変更のないまま Generate を Cancel しても "Current changes were not applied." と出る | cancel の通知の文言 |

### UI/UX/アクセシビリティ（UI）

| ID | 症状 | 原因 |
|---|---|---|
| UI-02 | 選択肢を出す modal 6 種（Label Not Shown、Color/Stroke Change Scope、Reset Fill Color、Feature Visibility Scope、Label Text Scope、Legend Rename）に role と名前がない。focus は背後に残り、Tab で外へ出られ、Escape でも閉じない。しかも Escape は背後の feature popup を閉じ、保留中の編集が宙に浮く。Label Not Shown には Cancel がない | 良い例（`trapDialogFocus`、Linear source の削除ダイアログ）がすでにあるが、使っていない |
| UI-03 | feature popup をキーボードで開いても focus が移らず、閉じるボタンまで Tab を約 14 回押す必要がある。Escape で閉じると focus が body に落ちる（popup での PV-11 と同類）。SVG の feature に tabindex がない | popup の focus 管理 |
| UI-04（+TK-16） | 名前のない操作がある: 「+」ボタン 2 個、feature を外す X ボタン 7 個（16 px）、凡例の caption 入力 6 個、similarity group の検索欄と並び替え。drawer の tab に `role=tab` と `aria-selected` がない。help tip 66 個がすべて "Help" という名前。行ごとの入力が同じ名前になる（Inner gap ×5、Outer gap ×5 など）。PNG、PDF、Interactive SVG のボタン名にアイコンの文字が混じる。batch の Result の選択に名前がない | `index.html` |
| UI-05 | 英語だけの UI なのに `<html lang="ja">` | `index.html:2` |
| UI-06 | slate-400 の文字のコントラストが 2.45:1 しかない（"Click to Browse"、"(auto)" 17 個など） | 配色 |
| UI-07 | 入力なしで Generate すると、Pyodide を先に起動する（Circular 26 秒、Linear 42 秒）。その後で "Supply GenBank input…" と出る。"Retry Generate" を押しても解決しない | 入力の検査が Worker の起動の後にある |
| UI-08 | 中身が GenBank でないファイル（空、FASTA、名前を変えた SVG）でも、アップロード欄には緑のチェックが出る。その下には "No records were found" と出る | upload-zone の表示が検査の結果を見ていない |
| UI-09 | "Single-record crop, orientation and titles" の全欄、Save Raw LOSAT TSV、Clear Cache が、理由を示さずに無効になる | 無効の理由を表示しない |
| UI-10 | 電話の幅（390 px）では、preview がページの約 1850 px 下から始まり、高さは 254 px しかない。Editor の drawer は 203 px で、実質使えない | レイアウト |
| UI-11 | fit-to-window がない。1280×800 では 996×817 の図が 896×386 の枠に入り、1920 でも切れる | preview |
| UI-12 | Generate Diagram のボタンは画面下に固定されているが、Tab 順では設定の途中（Basic の後）に来る | DOM の順序 |
| UI-13 | 1280 と 1920 で、設定 panel の中に横スクロールが出る（入力が 3 px はみ出す）。Generate の上に約 70 px の空きがある。"Custom Track Slots" の見出しが "Custom ..." に切れる。読み込むたびに Tailwind Play の警告が出る | UI-01 を直した後に確認し直す |

### 入力と比較（CI）

| ID | 症状 | 原因 |
|---|---|---|
| CI-01 | Linear で、record の長さを完全に超える region（60 bp の record に 1000–2000）を指定すると、説明のない RENDER_FAILED になる。CLI の文 "Start position (1000) must be less than end position (2000)" も誤り | `gbdraw/crop_genbank.py` の `check_start_end_coords` が end を切り詰めてから比べる。Web の Linear 経路に事前の検査がない |
| CI-05 | Match popup の "Query span / Subject span" の行は検索の座標系の値を出し、同じ popup の Interval の行や FASTA のヘッダと食い違う（crop、逆相補、rotate のとき）。CO-10 の残り | `pairwise-match-popup.js` の `buildMatchSpans` が、生の `data-qstart` などを使う |
| CI-07 | 比較が何も描かれないのに通知がない場合がある: (a) 1 ファイルに 2 record で Run LOSAT。両方が row 1 になり、pair がない。(b) Upload BLAST TSV を選んだがファイルがない。(c) 3 ゲノムで最初の pair に表を付けると plan が Selected に変わり、#2–#3 の比較が消える。(d) outfmt 7 の `# Fields:` の順が標準と違っても位置で読まれ、0 本になる（CLI も同じ） | 警告がない |
| CI-08 | gzip などのバイナリを GenBank 欄に入れると、`VALIDATION_UNCLASSIFIED` "Input validation failed." になる | 分類がない（CLI は utf-8 の decode エラーを出す） |

### 色・検索・注釈・凡例・出力（FL）

| ID | 症状 | 原因 |
|---|---|---|
| FL-04 | 凡例の改名ダイアログで "Merge into existing GC skew (-)" を選ぶと、live では統合されるが、Generate で失われる。同じ名前の行が色違いで 2 行になる | 凡例の改名の永続化 |
| FL-05 | 削除した凡例の項目を UI から戻せない。`restoreDeletedLegendEntries` は export されているが、どこからも呼ばれない | 配線がない |
| FL-08 | Interactive SVG の検索欄で Enter を押しても検索が走らない（Search ボタンだけが効く） | `standalone-interactivity-assets.js`（約 3451 行） |
| FL-09 | Interactive SVG の Qualifier key 欄の幅が 18 px しかない | 同上 |
| FL-10 | Output Prefix の誤り（使えない文字、`CON`、`nul.txt`、末尾の `.`）は、項目名のない `VALIDATION_UNCLASSIFIED` になる。300 文字では `CLEANUP_FAILED`（"save a Session and reload"）が出る。`../../x` は受け付けられ、Result の名前とダウンロード名が食い違う | prefix の検証と分類 |
| FL-11 | Bakta preset が、黙って CDS の既定色を `#cccccc` にし、凡例のフォントを 12 にする。panel の説明と合わず、palette も変わらない | preset の定義と説明 |
| FL-12 | qualifier の違うルールの間では、SPECIFIC RULES の Move up/down が効かない（描画は qualifier key の順に並べる。CLI も同じ）。UI はそれを説明しない | `gbdraw/features/selector_values.py` の順位付け |
| FL-13 | Specific colors の TSV で、Web は `#` のコメント行と空行を受け付けるが、CLI は "Missing values" で拒否する | Web と Python の reader の違い（#800 の `read_literal_table` で解消する可能性がある。merge 後に確かめる） |
| FL-14 | Region Annotations が黙って値を変える: 重複した id → `region_1_2`、空の id → `region`、重複した set id → `annotations_2`。空の TSV を Import すると、確認なしで全 set が消える。Import のエラーは native の `alert()` だけで、panel の notice は空のまま | annotations の入力の正規化 |

## 修正の状況（2026-10-08）

59 件の状況: 修正済み 45、別の PR で修正済み 11、説明を足した（挙動は維持） 2、0.15.0 1（合計 59）。修正した PR は #939（A-out）、#941（A-ui）と Lane B（`fix/web-gui-audit-b`）。FL-10 は Python 側が #939、Web 側が #945 なので「修正済み」に数えた。UJ-03 の行は FL-01 を含む。

| ID | 優先 | 状況 | PR / 確かめ方 |
|---|---|---|---|
| TK-01 | P1 | 修正済み | #939 |
| TK-02 | P2 | 修正済み | #939 |
| TK-03 | P2 | 別の PR で修正済み | #945（PER-MODE）。`depth-track-session` の TK-03 case と `depth-track-state.test.mjs` が、`29c631e2` で通る |
| TK-04 | P2 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-tracks）。`track-slot-row-edits.test.mjs` |
| TK-05 | P2 | 別の PR で修正済み | #805 |
| TK-06 | P3 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-depth）。`DEPTH_INVALID` / `DEPTH_SERIES_SOURCE` |
| TK-07 | P3 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-tracks）。欄に誤りを出し、行は直前の有効な index を保つ |
| TK-08 | P3 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-depth）。help の文を直した |
| TK-09 | P3 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-tracks）。Hide GC は **Dinucleotide** の行だけに届く |
| TK-10 | P3 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-tracks）。Reset は Show Depth に関わらず Depth の行を作る |
| TK-11 | P3 | 修正済み | #941（UI-01 の刷新で解消。guard で固定） |
| TK-12 | P3 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-depth）。小数の parser を一つにした |
| TK-13 | P3 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-depth）。underlay の検査を Python に一本化 |
| TK-14 | P3 | 修正済み | #941（docs を直した。disclosure は閉じたまま） |
| TK-15 | P3 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-depth）。Generate 前は「≈ … (estimate)」 |
| UI-01 | P2 | 修正済み | #941 |
| UI-02 | P3 | 修正済み | #941 |
| UI-03 | P3 | 修正済み | #941 |
| UI-04 | P3 | 修正済み | #941（TK-16 を含む） |
| UI-05 | P3 | 修正済み | #941 |
| UI-06 | P3 | 修正済み | #941 |
| UI-07 | P3 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-run）。Worker を起こす前に入力を検査する |
| UI-08 | P3 | 修正済み | #941（HANDOFF は A-out としていたが、#941 が直した） |
| UI-09 | P3 | 修正済み | #941 |
| UI-10 | P3 | 0.15.0 | Owner の決定（2026-10-07）。multi-drawing（v0.15.0）で扱う |
| UI-11 | P3 | 修正済み | #941（Fit ボタン） |
| UI-12 | P3 | 修正済み | #941 |
| UI-13 | P3 | 修正済み | #941 |
| CI-01 | P3 | 修正済み | #939 |
| CI-02 | P2 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-run）。`run-info.test.mjs` |
| CI-03 | P2 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-run）。単一 record から単一 record への置き換えは D-32 のとおり crop と Definition を残す |
| CI-04 | P2 | 修正済み | #939 |
| CI-05 | P3 | 修正済み | #939 |
| CI-06 | P2 | 修正済み | #939 |
| CI-07 | P3 | 修正済み | #939（a、b、d を直した。c は一行の説明を足し、挙動は維持。Owner-delegated） |
| CI-08 | P3 | 修正済み | #939 |
| UJ-01 | P2 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-labels）。`live-generate-parity` の `label Default` 5 case |
| UJ-02 | P3 | 別の PR で修正済み | #937（E1）。`per-mode-results` の UJ-02 case が `29c631e2` で通る |
| UJ-03 | P3 | 別の PR で修正済み | #857、#870（FL-01 を含む） |
| UJ-04 | P3 | 修正済み | #939 |
| UJ-05 | P3 | 修正済み | #939 |
| UJ-06 | P3 | 修正済み | Lane B（`fix/web-gui-audit-b`）（**Load an example**。例の Session は pip package にも入れる。D-B05） |
| UJ-07 | P3 | 修正済み | #939 |
| UJ-08 | P3 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-run）。入力のエラーだけを、検出が成功した時に消す |
| UJ-09 | P3 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-session）。**Load Session** の前に確認を出す |
| UJ-10 | P3 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-run）。変更がなければ「It matches the current settings.」 |
| FL-02 | P2 | 別の PR で修正済み | #940（OV-152）。`live-generate-parity` の case が `29c631e2` で通る |
| FL-03 | P2 | 別の PR で修正済み | #892（OV-63） |
| FL-04 | P3 | 別の PR で修正済み | #894（OV-62） |
| FL-05 | P3 | 別の PR で修正済み | #940（OV-154）。`source-legend-reconciliation` M14 が `29c631e2` で通る |
| FL-06 | P2 | 別の PR で修正済み | #892（OV-63） |
| FL-07 | P2 | 別の PR で修正済み | #940（OV-127）。`source-legend-reconciliation` M10 が `29c631e2` で通る |
| FL-08 | P3 | 修正済み | #939 |
| FL-09 | P3 | 修正済み | #939 |
| FL-10 | P3 | 修正済み | Python 側は #939、Web 側は #945（PER-MODE）。`gui-audit-20260930-options` の FL-10 case と `output-prefix.test.mjs` が `29c631e2` で通る |
| FL-11 | P3 | 説明を足した（挙動は維持） | #941（Bakta preset の panel の文と help を直した。Owner-delegated） |
| FL-12 | P3 | 説明を足した（挙動は維持） | #941（rule の順序の一行と owner の docs。Owner-delegated） |
| FL-13 | P3 | 別の PR で修正済み | #812（OV-27）。tab だけの行は GX-07 で #939 が直した |
| FL-14 | P3 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-labels）。panel の通知と、空の TSV の Import の確認 |

### その他の観察

| 観察 | 状況 | PR / 確かめ方 |
|---|---|---|
| feature 検索の query が mode を越えて残る | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-session）。検索は mode ごとに持つ |
| Linear の Auto の通知が Generate の後も未来形 | 修正済み | Lane B（`fix/web-gui-audit-b`）（b-session）。表示中の Result に合わせて現在形にする |
| Reset alignment のダイアログ | 修正済み | #941 |
| 1600 px の右 drawer | 説明を足した（挙動は維持） | #941。overlay は契約 `PD-OI-054` による |
| preview の popup が footer で切れる | 修正済み | #941、さらに Lane B（`fix/web-gui-audit-b`）（b-session）で全 popup の下限を一つにした |
| region と逆相補の順序 | 修正済み | #939（`docs/CLI_Reference.md` と入力形式の docs） |
| Interactive SVG の popup の題 | 修正済み | #939 |
| 凡例の改名後の移動、Feature Edits TSV、scope の選択で popup が閉じる | 本節の対象外 | `REVERIFICATION.md` の記録のまま |

### 新しく見つけた不具合（GX）

| ID | 扱い |
|---|---|
| GX-01 | 修正済み。#941（`index.html` の bindings）と Lane B（`fix/web-gui-audit-b`）（b-tracks。Duplicate の predicate） |
| GX-02 | 修正済み。#941（zoom の上限と刻み） |
| GX-03 | 修正済み。#939（Gallery の refresh tool） |
| GX-04 | 修正済み。Lane B（`fix/web-gui-audit-b`）（gx0405）。隣り合う二行は向き合う gap の大きい方を保つ。基準 SVG は変わらない |
| GX-05 | 修正済み。Lane B（`fix/web-gui-audit-b`）（gx0405）。GX-17 の規則で、例のコマンドが描ける。resolver への別の変更はない |
| GX-06 | 修正済み。#939（`TRACK_LAYOUT` / `CANNOT_FIT`） |
| GX-07 | 修正済み。#939（`read_table_lines`） |
| GX-08 | 修正済み。#939（Protein ID の行） |
| GX-09 | 修正済み。#939（`REGION_INVALID`） |
| GX-10 | 修正済み。Lane B（`fix/web-gui-audit-b`）（gx10、b-labels）。caption は 1 始まりで、strand は付けない（D-B01） |
| GX-11 | 修正済み。#939（Interactive SVG の検索欄） |
| GX-12 | 修正済み。#941（設定欄の `:disabled` と見た目） |
| GX-13 | 修正済み。#941（暗くするのは操作部品だけ） |
| GX-14 | 修正済み。#941（tutorial の caption と alt） |
| GX-15 | 修正済み。#941（match popup の文字色） |
| GX-16 | 別の PR で修正済み。#945（PER-MODE、OV-161）。`feature-fill-scope` の case が `29c631e2` で通る |
| GX-17 | 修正済み。Lane B（`fix/web-gui-audit-b`）（gx0405）。Auto の幅の行は、radius を指定しても Auto と同じく縮む（Owner の決定、2026-10-08） |
| GX-18 | 修正済み。Lane B（`fix/web-gui-audit-b`）（b-depth）。ticks の行の Radius の注は tick の anchor を示す |
| GX-19 | 修正済み。Lane B（`fix/web-gui-audit-b`）（gx0405、D-B02）。**Use custom stack** を入れただけでは図が変わらない |
| GX-20 | 修正済み。Lane B（`fix/web-gui-audit-b`）（b-depth）。Legend Name Scope の Cancel はすぐ閉じる |
| GX-21 | 修正済み。Lane B（`fix/web-gui-audit-b`）（b-labels）。Auto Reflow が on の時の label Off の Undo |
| GX-22 | 修正済み。Lane B（`fix/web-gui-audit-b`）（b-session）。popup のドラッグ後の click で閉じない |
| GX-23 | 修正済み。Lane B（`fix/web-gui-audit-b`）（b-session）。`test_web_packaging.py` の path |

## その他の観察（未分類、軽微）

- feature 検索の query が mode を切り替えても残り、すべての feature が薄く表示される（UJ）。
- Linear では、Generate の直後でも "Auto will show these fields… on the next successful Generate" が出続ける（UJ）。
- Reset alignment のダイアログが、右上の header や Session ツールバーに重なるカードとして出る（CI）。
- 1600 px の幅では、右の drawer が preview に重なる。preview の下端の popup が footer で切れる（CI）。
- Interactive SVG の popup は tRNA を gene（`TRNF`）で題し、アプリの popup は product（`tRNA-Phe`）で題す（FL）。
- 凡例の feature の項目を改名すると末尾へ移る。Feature edits TSV に fill と stroke の編集が入らない。scope を選ぶと popup が閉じる（FL。仕様かどうかは未確認）。
- 未確認の疑い: 一部の track の行の操作に `semanticMutationAvailable` の無効化がない（TK）。region と逆相補の順序（`500-3000:rc` と、`reverse_complement=1` の表の行）で結果が違うが、docs に順序の記述がない（CI）。

## 以前の修正の再確認（退行なし）

- 入力: IN-01、IN-02、IN-03、IN-06、IN-08（今は明確な `INPUT_REQUIRED`）
- Session と History: SE-01、SE-02、SE-03、SE-04、SE-05、GE-06
- Feature 編集: FE-01、FE-06、FE-07、FE-11、FE-12
- 比較: CO-02、CO-03（ジオメトリ）、CO-04、CO-05、CO-06（ribbon。文言は CI-04）、CO-07
- 凡例と出力: PV-02、PV-03、PV-04、PV-05、PV-10、PV-12
- トラック: TR-01、TR-02、TR-03、TR-06、TR-07、TR-08、TR-09、TR-10、TR-11
- 一部だけ: IN-04（CI-03）、CO-10（CI-05）、PV-11（Editor の drawer は直った。feature popup は UI-03）

## 確認して問題がなかったもの（要約）

- **ジャーニー**:
  - 初回の利用 → popup の編集 → SVG、PNG、PDF、Interactive SVG の出力
  - 複数 record、crop、逆相補 → Session の保存、再読み込み、Load → 同じ SVG
  - Linear と比較: BLAST 表、LOSATN、TLOSATX、LOSATP（similarity groups、collinear）
  - mode の往復、Undo/Redo（Generate、Reset Settings、mode 切替を含む）、batch の Result
  - Cancel（runtime の起動中、描画中、LOSAT の実行中）、Generate の二度押し
  - Gallery のリンク 247 本と Session 5 件の再現
- **トラック**:
  - 測定値の入力の検証（px と %、0、負、Infinity）
  - stack の Undo/Redo、Session の往復
  - Source recipe の CLI 再現（既定、編集後、Depth 1–2 series、Multi-Record Canvas、Linear Depth）
- **色と凡例**:
  - 範囲の違う fill と stroke、`-d` と `-t` を CLI と比べて一致
  - 凡例の改名と衝突ダイアログ、並べ替え、全 Position、ドラッグ
  - title、definition、Qualifier Priority を CLI と比べて一致
  - PDF の寸法が CLI と同じ（921.56×679.37 pt）
- **入力**: Prokka、DDBJ、CRLF、重複 ID、GFF3 の入れ子、BLAST の filter と不正値のメッセージ（項目名つきの INPUT_INVALID）
- **ページ全体**: 全体の横スクロールなし、重複 ID なし、ダークモードは未対応（明るいまま読める）

## 共通する原因

1. **エラーの分類が足りず、UNKNOWN や汎用文になる**: TK-02、TK-03、FL-03、FL-06、CI-04、CI-06、CI-08、FL-10、TK-06
2. **live の投影や後処理が、Generate や History と一致しない**: UJ-01、UJ-03/FL-01、FL-02、FL-04、FL-07
3. **状態の掃除が漏れる**（置き換え、削除、renderer の変更）: TK-04、TK-05、CI-03、FL-02、UJ-08
4. **再現性**（Source recipe、Gallery Session）: CI-02、UJ-04
5. **ダイアログと focus の管理**: UI-02、UI-03
6. **CSS の土台**: UI-01（TK-11、UI-13 にも影響する）

## Owner の決定（2026-10-05）

この session は監査だけなので、どれもまだ実装していない。

| 項目 | 決定 |
|---|---|
| UJ-03/FL-01 | live の凡例を Generate に合わせる（Generate と同じ分割と並び順） |
| UJ-09 | Load Session の前に確認を出す |
| TK-13 | CLI に合わせる（record に underlay の feature がなければ、Features の行のない stack を受け付ける） |
| UI-01 | **「提案」を採用する**（案 A に `probes/ui01_proposal.css` の内容を加えたもの。比較ページ Version 2: https://claude.ai/artifact/2ub5cMWzCBxURVGzXy6ckE）。次のセッションで、バグの修正と一緒に実装する |
| CI-07c、FL-12 | 未決定。次のセッションでは推奨案（今の挙動のまま、UI に説明を足す）で進め、Owner-delegated として記録する |

計画と次のセッションのプロンプト: `PLAN.md`、`NEXT-SESSION-PROMPT.md`（このディレクトリ）。

別セッション（override-residuals）と重なるもの（2026-10-05 夕方の時点）:
- TK-05 は OV-29 として #805 で修正済み。確認だけする。
- FL-13 は OV-27 の #812 で直る見込み。
- UI-04 の一部（アイコンだけのボタンの名前）は OV-30 の #809 で直る。

## Owner の判断が要るもの（提示した選択肢）

- **UI-01**: 次のどちらにするか。
  - A（推奨）: style を `type="text/tailwindcss"` で効かせ、設計どおりの見た目にする。全画面の幅と見た目が変わるので、レイアウトの test と docs のスクリーンショットを撮り直す。
  - B: 今の見た目を正とし、死んでいる `@apply` を消す。この場合も、入力欄の枠（WCAG 1.4.11）とボタンの見た目は別に直す必要がある。
- **UJ-03/FL-01**: live の凡例を Generate と同じ分割にするか、"on regeneration" のまま受け入れるか。
- **UJ-09**: Load Session の前に確認を出すか。History を残すか。
- **TK-13**: underlay の検査を record の feature で行い、CLI に合わせるか。
- **CI-07c と FL-12**: 今の挙動のまま、UI に説明を足すか。

### UI-01 案 A の実装（スクリーンショットで確かめた形）

- `<style>` の `type` を変えるだけでは足りない。Play は型付きの style を utility の後ろに出力する。そのため `.form-input` の `text-sm p-2` が、要素ごとの `px-1 text-xs` などを上書きする（`.btn` の `text-base` が `text-xs` を上書きするのも同じ）。
- `:where()` で詳細度を 0 にする方法も使えない。preflight（`summary { display: list-item }`、`button { background-color: transparent }`、`input { padding: 0 }`）に負けて、開閉の矢印、主ボタンの色、入力欄の余白が消える。
- 採った形:
  - 18 個の `@apply` ルールを、`<style type="text/tailwindcss">` の `@layer components` に置く。
  - `.form-input` を上書きする普通の CSS 3 つ（`.form-input-compact`、`.adv-options .form-input-compact`、`.auto-value-input`）も、同じ層の後ろに置く。
  - そうしないと compact の入力欄 58 個の文字が欠け、"(auto)" の表示 23 個が隠れる。
- この形で文字が欠ける入力欄は 1 個（Preset scheme の select、`h-7` と `p-2`）。コンソールエラーはない。
- 差し替えの処理は `probes/ui01_transform.py`、撮影は `probes/ui01_screens.py`。
- 案 A を採るなら、ほかにも直す:
  - "PRIORITY FILE (TSV)" が大きな見出しになる。
  - "CENTER RESERVED RADIUS" が 3 行に折り返す。
  - "Custom Track Slots" の見出しが縮む。
  - Retry Generate が小さいまま。
  - checkbox が 20 px になる。
  - 設定パネルが縦に長くなる。
  - レイアウトの test と docs のスクリーンショットを撮り直す。

### UI-01 の改善の提案（2026-10-05、Owner の依頼で作った見本）

Owner の指摘: 案 A のアップロード欄は縦に大きすぎる（60 px）。これを受けて、案 A の上に重ねる CSS を作った（`probes/ui01_proposal.css`）。
スクリーンショットは `screenshots/ui-01/raw/*-proposal.png`、3 案を並べた画像は `triples/`。比較ページは同じ URL（Version 2）。

- **アップロード欄**: 1 行で 32 px。1 px の点線にする。未選択は灰色、ホバーで青、選択後は薄い緑の実線。黄色の地はやめる。
- **文字の階層**: 4 段にする。
  - カードの見出し: 14 px、semibold。
  - グループ見出し: 11 px、semibold、大文字。大文字はこの段だけにする。
  - 項目名: 12 px、medium、普通の書き方。
  - 補足: 11 px 以上、slate-500。
  - 結果: 設定パネルの太字（700）が 124 個から 40 個に減る。
- **色**: アクセントは青 1 色にする。緑は成功と選択済み、赤は削除とエラーだけに使う。
  - Result Preview の緑の枠、Input Genomes の青い左帯、Reset Settings の黄色をやめる。
  - 出力ボタンは PNG だけが primary なので、4 つを同格にする。
  - ページの地を #eef2f6 にする。
- **大きさ**: checkbox を 16 px にする。
  - `.btn` と `.input-label` の `text-base` は、24 px の行の高さも持ち込む。要素の `text-[10px]` は文字の大きさしか変えないので、行の高さは 24 px のまま残る。これを 1.25〜1.3 にする。
  - Retry Generate に `.btn` の余白を付ける。
- **マークアップの変更が要るもの**:
  - Custom Track Slots の見出し行の混雑。
  - "Reset to Tuckin / Middle / Spreadout" の分け方。
  - 補足の `text-[9px]`（87 か所）と `text-[10px]`（423 か所）を、1 つの class にまとめる。

## 推奨する着手順

1. TK-01（P1。CLI も直る）。
2. UNKNOWN になる P2: TK-03、FL-03、FL-06、CI-06。原因が小さく、1 つの PR にまとめられる可能性がある。
3. 黙って誤る P2: FL-02、UJ-01、CI-03、TK-05、TK-04、CI-02。
4. UI-01 は Owner の判断の後で、P3 の見た目の項目（TK-11、UI-13）と一緒に直す。
5. 軽微な P3（UI-05 の `lang`、TK-14 の docs、UJ-10 の文言、TK-09 の文法、UJ-04 の Gallery の prefix など）は、まとめて小さな PR にできる。

注意: 実行中の override-residuals セッションと重なる箇所がある。UJ-01 は label-actions（R-4 #799 の近く）、FL-13 は #800。着手の前に、そちらの merge を待つか調整する。
