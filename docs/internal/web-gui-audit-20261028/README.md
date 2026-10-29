# Web GUI リリース前監査 — dev 33d9e234c、e7b3c3a06 と 639de764a（2026-10-28）

対象は 0.14.3 の promotion 前の `origin/dev` 3 点。release tier の `Tests`（`workflow_dispatch`、`tier=release`）run 40634318370（`33d9e234c`、seed 40634318370）、run 40648501950（`e7b3c3a06`、seed 40648501950）と run 40655413787（`639de764a`、seed 40655413787）の結果をまとめる。seed は run ID。

**件数**: 固有の問題 3 件（P2 3: OV-1985、OV-1987、OV-1988）。#1433 の review で同じ経路の OV-1986 を追加（再現せず）。どの run でも `Promotion audit / random-walk` だけが失敗し、その原因は run 40634318370 では OV-1985、run 40648501950 では OV-1987、run 40655413787 では OV-1988 である。修正は OV-1985 が #1433（`fix/ov1985-legend-stroke-caption`）、OV-1987 が #1438（`fix/ov1987-legend-stroke-leave`）、OV-1988 が #1441（`fix/ov1988-legend-fill-palette`）。

追加の乱歩（dev `d08c1878e`、seed 20261101 から 20261108、fixture ごとに 30 step）で 3 件（P2 3: OV-2187、OV-2188、OV-2189）。8 seed のうち 4 seed が失敗し、すべて live と Generate の Legend の差。修正は 3 件とも #1445（`fix/rw-legend-parity-0143`）。OV-2189 は Owner の判断（2026-10-29、A）に従い、Python の rule 合成が、置き換える rule と同じ条件・同じ色で名前だけを変える rule をその rule の表の位置に書く形にした。3 件とも 0.14.2 からの回帰ではない。OV-2187 と OV-2189 は v0.14.2 で同じ seed を再生し、同じ step で失敗した（seed 20261101、20261108）。OV-2188 は v0.14.2 で再生しておらず、関係するコード（追加する行の複製元を選ぶ行）が v0.14.2 と同一であることから判断した。

## 方法と表記

- 手順は `docs/internal/WEB_PERIODIC_AUDIT.md` に従う。`random-walk` は seed 付きの乱歩（`tests/web/live-generate-random-walk.promotion.spec.js`）で、各 step で live 表示と Generate 結果を比べる。
- 重大度の基準は 2026-09-30 の監査と同じ（P1 誤った出力・編集の消失・有効入力で Generate 不能、P2 機能の不具合・古い状態の残存、P3 軽微）。
- JS のパスは `gbdraw/web/js/` からの相対パス、それ以外はリポジトリ直下からのパス。行番号は書かない。

## 結果（ジョブ別）

| ジョブ | 結果 | 項目 |
|---|---|---|
| random-walk（run 40634318370、dev `33d9e234c`） | 失敗（fixture "Linear, one Result" の step 6 と 7） | OV-1985 |
| ほかの `Promotion audit / ...` と release ジョブ | 成功 | なし |
| random-walk（run 40648501950、dev `e7b3c3a06`） | 失敗（同じ fixture の step 5 から 8） | OV-1987 |
| random-walk（run 40655413787、dev `639de764a`） | 失敗（fixture "Circular, two-Result batch" の step 12） | OV-1988 |

## 見つかったもの

### P2（3 件）

| ID | 症状 | 原因 |
|---|---|---|
| OV-1985 | 同じ caption の color rule が 2 つの色で描かれている Legend 行（Python の key は `CDS [#54bcf8]`）に、feature popup の "Legend item" scope で stroke を付ける。その後、もう一方の rule の caption を変えると（`CDS` → `transfer RNA`）、Legend の CDS 行の枠線が live では既定の灰色（`#808080`）に戻り、Generate では指定色（`#123456`）で描かれる。Generate の後で caption を `CDS` に戻したときも同じく live だけ灰色になる。Generate すれば正しく描かれる | rule commit（`commitOnce` in `app/feature-editor/rule-actions.js`）は、key の変わる行の記録を新しい key に移す（`moveLegendRecords` in `app/legend/entry-actions.js`、OV-660）。stroke の編集（`legendStrokeOverrides`）も `CDS [#54bcf8]` から `CDS` に移るので、editor intent と Generate は正しい。ところが commit の live 表示（`projectPaletteAndRules` in `app/app-setup.js`）は fill と Legend 構造の domain だけを `showEditorIntent` に渡し、stroke の domain を compile しなかった。そのため stroke は隠れた古い行に残り、表示される行には付かなかった。#1042（v0.14.0）から存在し、v0.14.2 にもある（関係するコードは v0.14.2 と同一）。0.14.2 からの回帰ではない |
| OV-1987 | Legend の `CDS` 行に stroke を付けた後、CDS の feature を別の caption（`transfer RNA`）で描く color rule を足すと、Python は CDS の palette 行を "other proteins" と名付け、`CDS` の行は描かれない。その後の rule の色の変更（または Legend 行の削除）で、live では元の CDS 色の feature 3 つに stroke（`#7b2cbf`）が戻り、Generate では既定の灰色（`#808080`）で描かれる | Legend 行の stroke が届く feature は `legendRowFeatureIds`（`services/legend-svg.js`）が決める。live の compile は draft の行の色（`draftColor`、ここでは key `CDS` の palette 色 `#54bcf8`）を渡し、関数は Result にその行がないとき（Python の行の色が null）も、その色の feature に stroke を届けた。Generate の compile は draft の色を渡さないので、行がなければ何にも届かない。v0.14.2 から存在し（同じコード）、0.14.2 では Legend 行の削除など stroke を表示する step で見える。#1433 は fill の表示でも stroke を表示するので、rule の編集でも見えるようになった（bisect: `33d9e234c` は step 7 のみ、`6f5e2c430` と `e5ed39352` は step 5 から 8） |
| OV-1988 | visibility rule で隠した CDS の feature に、feature popup の "Use existing color" で色を付ける。これは caption `CDS`、色は CDS の palette 色の rule を書く。rule の色を塗る feature は描かれないので、Python はこの rule の行を描かず、`CDS` 行は palette の行のまま。その後 palette を切り替えると（`herbal_garden` → `earthtones`）、Legend の `CDS` 行の塗りが live では rule の色（古い palette の `#addc7f`）のまま、Generate では新しい palette の色（`#8a7f8d`）になる。2 つの行が同じ名前 `CDS` を見せる状態ではなく、比較した行は両側とも key `CDS` の同じ行 | live の compile が Legend 行に与える draft の色（`draftLegendRowColors` in `app/candidate-render.js`）は、caption を持つ rule をすべてその行の key（ここでは `CDS`）に割り当て、rule の色を palette の色より優先した。描かれる feature を 1 つも塗らない rule も含めた。Generate では Python がその行の色を決めるので、rule の行はなく palette の色になる。関係するコード（`candidate-render.js`、`services/specific-color-rules.js`）は v0.14.2 と同一で、0.14.2 からの回帰ではない |

**修正（OV-1985、#1433）**: drawing が Legend 行の stroke を 1 つでも持つ間は、feature の fill を表示するたびに（rule commit、palette の変更、Undo と Redo の色の step など）同じ compile で stroke の domain も表示する（`fillDomainsOf` in `app/app-setup.js`）。Legend 行の stroke は行の色で描かれた feature に届くので、fill の変更は行の key だけでなく、行が届く feature も変えうる。Legend 行の stroke がない drawing の作業は変わらない。

**関連（OV-1986、#1433 の review で指摘）**: rule commit が stroke 付きの Legend 行の届く feature を変えるとき（rule を消して行がなくなる、pattern を変えて feature が別の行に移る）、live の stroke が古いまま残るおそれ。review が示した 2 つの手順（rule `Special` の削除、rule `Other` の pattern 変更とその Undo と Redo）は dev `33d9e234c` でも live = Generate だった。どちらの commit も自動の再描画を要求し、それが stroke も描き直すため。上の修正は commit 自身の表示でも stroke を描くので、再描画を待たずに正しく、再描画しない経路でも正しい。再描画のない経路で差の出る手順は見つかっていない。

**修正（OV-1987、#1438）**: Result が描かない行の stroke は、draft がその key にどの色を与えても色では feature に届かない（`legendRowFeatureIds`）。color の編集で行に名指しされた feature や古い Session が行に列挙した feature には、live でも Generate でも従来どおり届く。Generate の結果は変わらない。

**修正（OV-1988、#1441）**: palette の切り替えでも Default colors の変更でも、rule が行に色を与えるのは、表示中の Result が描く feature のうち、その rule の group（caption と色）を最初の rule として取るものがあるときだけ。描かれる feature が最初に取る rule と、すべての match が分かっているかは 1 つの owner（`drawnRuleMatches` in `services/specific-color-rules.js`）が 1 回の走査で決め、`draftLegendRowColors` と Python が描かない行を決める `rowsLeftOut` がその結果を読む（live の compile では 1 回だけ計算して両方に渡す）。描かれる feature の match がまだ分からない間は、従来どおりすべての rule が行に色を与える（`rowsLeftOut` は従来どおり、分からない match を「塗らない」と読む）。Legend パネルの swatch 色（`draftLegendPanelColors`）と行の stroke が届く feature も同じ色を読む。Generate の結果は変わらない。

## guard の判断

- OV-1985: guard を追加した。`tests/web/live-generate-parity-legend.playwright.spec.js` の 2 件（caption を変えたとき、Generate の後で caption を戻したとき。いずれも base で失敗）を PR の Playwright functional shard が実行する。
- OV-1986: 同じ spec に 3 件（rule の削除、pattern の変更、その Undo と Redo）を追加した。base でも成功するので結果を固定する guard であり、失敗から始まる test ではない。
- OV-1987: 同じ spec に 1 件（base `e5ed39352` で失敗）と、`tests/web/candidate-render.test.mjs` に `legendRowFeatureIds` の 1 件を追加した。
- OV-1988: 同じ spec に 1 件（base `6ea8b36ca` で失敗）と、`tests/web/palette-default-colors.test.mjs` に `draftLegendRowColors` と `draftLegendPanelColors` の 1 件（base で失敗）を追加した。乱歩の seed 40655413787 も修正後に成功する。work guard の行は変わらない。
- work guard（`tests/web/helpers/live-work-allowlist.mjs` の `WORK_ALLOWLIST`）: setup が CDS 行に stroke を付けた後の fill を表示する 12 行（color rule commit、palette の変更、Apply to all、Save、Legend rename など）に `strokes` の stage を記録した。

## 追加の乱歩（dev `d08c1878e`、30 step）

| ジョブ | 結果 | 項目 |
|---|---|---|
| seed 20261101（fixture "Circular, two-Result batch" の step 8） | 失敗 | OV-2187 |
| seed 20261102（fixture "Circular, one Result" の step 21） | 失敗 | OV-2188 |
| seed 20261106（fixture "Circular, one Result" の step 12） | 失敗 | OV-2188 |
| seed 20261108（fixture "Linear, one Result" の step 18 と 19） | 失敗 | OV-2189 |
| seed 20261103、20261104、20261105、20261107 | 成功 | なし |

| ID | 重大度 | 症状 | 原因 |
|---|---|---|---|
| OV-2187 | P2 | batch の 1 つ目の Result の CDS に feature popup の "Use existing color" で色を付ける（caption `CDS`、CDS の palette 色の rule）。この rule は 2 つ目の Result の feature を塗らないので、2 つ目の Result では Python は palette の `CDS` 行を描く。2 つ目の Result を表示したまま rule を編集すると（caption の再入力、色の変更）、live ではその `CDS` 行が消え、Generate では描かれる。1 つの Result でも、rule が Feature visibility で隠した feature だけを塗るときに同じことが起きる | rule commit が退かせる行（`candidateLegendIntents` in `app/feature-editor/rule-actions.js` の `previousIntents`）に、caption を持つ rule をすべて入れた。表示中の Result で描かれる feature を 1 つも塗らない rule も含めた。Legend の準備（`prepareFileLegendEntries` in `app/legend/entry-actions.js`）は、名前と色が一致する行を rule の行と見なすので、palette の `CDS` 行を退かせた。関係するコードは v0.14.2 と同一で、v0.14.2 でも同じ seed が同じ step で失敗する |
| OV-2188 | P2 | Legend の並びを変えて GC content、GC skew、Depth のいずれかの行を先頭にした後（Sort、行の移動）、feature popup の Legend name や Legend の rename で rule の行を作ると、新しい行の枠線が live では無し（`none`）、Generate では block stroke の灰色（`#808080`）になる | rule の行をすぐ表示するとき、executor（`applyLegendOperations` in `services/svg-result-ingestion.js`）は表示中の並びで先頭の行を複製した。並びを変えると先頭は GC content の行（枠線なし）になる。Python は rule の行を block stroke で描き、Generate の executor は Python の並びの先頭（feature の行）を複製する。関係するコードは v0.14.2 と同一（v0.14.2 での再生はしていない） |
| OV-2189 | P2 | product の rule 2 つ（`alpha`、`beta`）の後、`beta` 行の唯一の feature に feature popup で Legend name `gamma` を付けると、live では `gamma` が `beta` の位置に残り（OV-158）、Generate では先頭に描かれる。乱歩の Sort（step 12）は並びを変えておらず、関係しない | 色の編集の rule 表は Python が作り（OV-244）、この編集ではその feature の rule を表の先頭に置く。Python は rule の行を表の順に描く。live の commit は OV-158（Owner 2026-10-07: 名前を変えた行はその位置に残る）に従い、行と生成時の並びの記録の両方を同じ位置で新しい key に移すので、並びは編集されたものと見なされず、Generate は Python の並びのまま描く。v0.14.2 でも同じ seed が同じ step で失敗する。#1445 で修正: Python の rule 合成は、feature がこれまで取っていた rule と同じ条件・同じ色で名前だけを変える rule を、その rule の表の位置に書く。それ以外の合成（hash rule の追加など）は変えない |

**修正（OV-2188、#1445）**: 追加する行は Python の並びの最初の行の複製にする。「feature Legend group の Python の並び」の owner を `pythonOrderedLegendRows`（`services/legend-svg.js`）1 つにした。editor の行（`direct-editor`）を除く行を、group が記録した Python の並び（`LEGEND_ORDER_RECORD`）の順に、記録がなければ描かれた順に返す。rule commit の行（`specific-color-file`）は、Python が次に描く rule の行の代わりとして含める（並びの記録はもともとそう数えている）。executor の行の追加（複製元）、並びの記録の書き込み（`orderLegendGroup`）、`legendRowTakesPlace` がこれを読む。`restoreLegendStructure` は記録をそのまま読む（`recordedGroupOrder`）。Generate の結果は変わらない。

**修正（OV-2187、#1445）**: rule commit が退かせる行は、表示中の Result で描かれる feature が最初の rule として取る rule の行だけにする（#1441 の `drawnRuleMatches` の `groupTaken`）。match が分からない間は従来どおりすべての rule の行。Legend editor が追加した行の記録（`addedLegendCaptions`）から rule の caption を外す処理は、従来どおりすべての rule の行で行う（`ownedCaptions`）。Generate の結果は変わらない。

**OV-2189（未修正）**: Product Contract に該当する行がなく、Owner の決定 OV-158（`tests/web/source-legend-reconciliation.playwright.spec.js` の M15: 行は live でも Generate でも位置を保つ）は live 側を意図した結果としている。直すには Generate の結果（または rule 表の順）を変える必要があり、Owner の判断を待つ。選択肢と推奨は handoff `/home/kawato/gbdraw-baselines/release-0143-20261027/handoffs/RW-LEGEND.md` にある。

guard: `tests/web/live-generate-parity-legend.playwright.spec.js` に OV-2188 の 2 件（popup の Legend name、Legend の rename。base で失敗）、OV-2187 の 2 件（caption の再入力、色の変更。base で失敗）、OV-2189 の `test.fail` 1 件（前提の rule の順と live の並びを case の中で確かめる。`test.fail` を外した実行では並びの差だけで失敗する）を追加した。`tests/web/specific-color-actions.test.mjs` に OV-2187 の 1 件（rule が塗らない Result では行を持たず、追加した行の記録からは caption を外す。base と前の head で失敗）を追加した。seed 20261101、20261102、20261106 は修正後に成功する。

## 証拠と再現方法

- run 40634318370（dev `33d9e234c`、2026-10-28）。失敗ログは `/home/kawato/gbdraw-baselines/release-0143-20261027/logs/tier-40634318370-failed.log`。
- 再現: `GBDRAW_RANDOM_WALK_SEED=40634318370 npx playwright test --config=tests/web/playwright/promotion.config.js --workers=1 -g "Linear, one Result"`。
- run 40648501950（dev `e7b3c3a06`、2026-10-28）。失敗ログは `/home/kawato/gbdraw-baselines/release-0143-20261027/logs/tier-40648501950-failed.log`。再現は上と同じコマンドで seed 40648501950。
- run 40655413787（dev `639de764a`、2026-10-28）。失敗ログは `/home/kawato/gbdraw-baselines/release-0143-20261027/logs/tier-40655413787-failed.log`。再現は seed 40655413787 で `-g "Circular, two-Result batch"`。
- 追加の乱歩（dev `d08c1878e`）: ログは `/home/kawato/gbdraw-baselines/release-0143-20261027/logs/rw/seed-<seed>.log`。再現: `GBDRAW_RANDOM_WALK_STEPS=30 GBDRAW_RANDOM_WALK_SEED=<seed> npx playwright test --config=tests/web/playwright/promotion.config.js --workers=1 -g "<fixture>"`。
- 調査メモ: `/home/kawato/gbdraw-baselines/release-0143-20261027/handoffs/OV-1985.md`、`/home/kawato/gbdraw-baselines/release-0143-20261027/handoffs/OV-1987.md`、`/home/kawato/gbdraw-baselines/release-0143-20261027/handoffs/OV-1988.md`、`/home/kawato/gbdraw-baselines/release-0143-20261027/handoffs/RW-LEGEND.md`。

## GUI 監査（dev `d08c1878e`、2026-10-28）: ヘッダー

0.14.3 の候補 `d08c1878e` を実ブラウザ（Chromium）で調べた（証拠は `/home/kawato/gbdraw-baselines/release-0143-20261027/evidence/gui-audit/FINDINGS.md` と `evidence/header-0143/`）。ヘッダーの 3 件はどれも v0.14.2 にもある（v0.14.2 の Web bundle を同じ probe で測った。v0.14.2 では en のはみ出しが 769〜778 px で最大 10 px、en の Session 名は 769〜910 px で 64 px 未満）。#1409 の Gallery & tutorials ボタンが原因ではない。

| ID | 重大度 | 症状 | 原因 | 修正と guard |
|---|---|---|---|---|
| OV-2177 | P3 | 1 段のヘッダーが 769〜800 px で window からはみ出す（fr 769〜800 px で最大 32 px、en と pt 769〜776 px で最大 8 px、id 769〜775 px で最大 7 px、es 769〜774 px で最大 6 px）。右端の翻訳 feedback リンクが切れ、ページが横にスクロールする | header の状態を選ぶ `chooseHeaderFit`（`app/header-fit.js`）は、最後の状態（アイコンだけ、「Session」の小見出しなし）でも足りないときに次の手がなかった。Session 名のボタンが 10 px まで縮んだ後は toolbar がはみ出す | 同じ owner に 2 つの状態を足した。Session 名がまだ縮むとき、閉じた Language メニューを表示中の言語名の幅にする（CSS `field-sizing: content`）。それでも Session 名が `SESSION_TITLE_MIN_PX`（80 px、既定の 16 px の root で 5rem）より狭いとき、768 px 以下と同じ 2 段のレイアウトにする。2 段のレイアウトの規則は `data-header-wrapped` の 1 か所にまとめ、768 px 以下では stylesheet の目印（`--app-header-two-row`）を header-fit が読んで常に 2 段にする。Circular と Linear の名前は 1 行に保ち（vi と th は 769〜930 px で 2 行に割れていた）、mode ボタンの transition は色だけにした（大きさの変化が測定を狂わせる）。header-fit は header、tagline の領域、Session 名のボタンの大きさの変化で測り直す（狭めたときも 2 段に移る）。guard: `tests/web/ui-language.playwright.spec.js` の feedback リンクの case が 11 言語で 769、800、850、800 px（最後は狭める手順）を確かめる（はみ出しがなく Session 名が 80 px 以上）。base では en 769 px が、#1446 の最初の head（`567f4ccf4`）では en の 850 → 800 px が失敗する。`tests/web/header-fit.test.mjs` が状態の順と 80 px の境を確かめる |
| OV-2178 | P3 | Session 名の欄が 810 px 付近から 769 px まで 0 px になる（en で 900 px 107 px、850 px 57 px、810 px 以下 10 px のボタンだけが残る）。名前が見えず編集もしにくい | OV-2177 と同じ。release notes の「小見出しを隠すので Session 名が幅を保つ」は 1 段のヘッダーの一部の幅でしか正しくなかった | OV-2177 と同じ修正。11 言語の 769〜1000 px（上りと下りの両方）で、1 段の Session 名は 80 px 以上、2 段では 160 px。1000〜1440 px と 768 px 以下（360〜768 px の高さと Language メニューの幅）は base と同じ。release notes の主張を「小見出しを隠すので Session 名が広くなる」に直した |
| OV-2186 | P3（未修正） | 768 px 以下の 2 段のヘッダーで、閉じた Language メニューが矢印だけになる（ja 480〜496 px、en 504〜520 px と 640〜660 px など、toolbar の行が折り返す幅のすぐ下。11 言語すべて） | `.app-language` の flex-basis 4.5rem が globe と feedback リンクを含むため、行の詰まる幅ではメニューに 28〜32 px しか残らない。OV-368 が受け入れたのは 390 px で English より長い名前が切れることで、矢印だけになることではない。v0.14.2 も同じ | 製品の判断が要るので直していない。選択肢: (A) このまま（0.14.2 からの回帰ではない）。(B) basis を 7.25rem にする。メニューは "English" の幅を保つが、ヘッダーが 1 行増える幅がある（en 360〜376 px、fr 360〜392 px で 135 → 173 px、ja 616〜648 px で 96 → 135 px）。(C) より小さい最小幅（例: メニュー 3.5rem）で数文字を見せ、1 行増える幅を狭める |

Show Definition と Default font size（Record Labels）が Generate まで図を変えないことは、R1（`gbdraw/web/CLAUDE.md`）と OIC-024 のとおり。両者は `data-setting="titles"`（`Applies on Generate`）の scope にあり、その card の先頭に注記がある。`tests/web/setting-application.test.mjs` は Generate の scope ごとに注記を 1 つ求め、control ごとには求めない。行ごとの注記を求める規則はない。
