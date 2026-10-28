# Web GUI リリース前監査 — dev 33d9e234c、e7b3c3a06 と 639de764a（2026-10-28）

対象は 0.14.3 の promotion 前の `origin/dev` 3 点。release tier の `Tests`（`workflow_dispatch`、`tier=release`）run 40634318370（`33d9e234c`、seed 40634318370）、run 40648501950（`e7b3c3a06`、seed 40648501950）と run 40655413787（`639de764a`、seed 40655413787）の結果をまとめる。seed は run ID。

**件数**: 固有の問題 3 件（P2 3: OV-1985、OV-1987、OV-1988）。#1433 の review で同じ経路の OV-1986 を追加（再現せず）。どの run でも `Promotion audit / random-walk` だけが失敗し、その原因は run 40634318370 では OV-1985、run 40648501950 では OV-1987、run 40655413787 では OV-1988 である。修正は OV-1985 が #1433（`fix/ov1985-legend-stroke-caption`）、OV-1987 が #1438（`fix/ov1987-legend-stroke-leave`）、OV-1988 が #PRNUM（`fix/ov1988-legend-fill-palette`）。

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

**修正（OV-1988、#PRNUM）**: rule が行に色を与えるのは、表示中の Result が描く feature のうち、その rule の group（caption と色）を最初の rule として取るものがあるときだけ（`ruleGroupDrawn` in `services/specific-color-rules.js`、`draftLegendRowColors` が読む）。描かれる feature の match がまだ分からない間は、従来どおりすべての rule が行に色を与える。Legend パネルの swatch 色（`draftLegendPanelColors`）と行の stroke が届く feature も同じ色を読む。Generate の結果は変わらない。

## guard の判断

- OV-1985: guard を追加した。`tests/web/live-generate-parity-legend.playwright.spec.js` の 2 件（caption を変えたとき、Generate の後で caption を戻したとき。いずれも base で失敗）を PR の Playwright functional shard が実行する。
- OV-1986: 同じ spec に 3 件（rule の削除、pattern の変更、その Undo と Redo）を追加した。base でも成功するので結果を固定する guard であり、失敗から始まる test ではない。
- OV-1987: 同じ spec に 1 件（base `e5ed39352` で失敗）と、`tests/web/candidate-render.test.mjs` に `legendRowFeatureIds` の 1 件を追加した。
- OV-1988: 同じ spec に 1 件（base `6ea8b36ca` で失敗）と、`tests/web/palette-default-colors.test.mjs` に `draftLegendRowColors` と `draftLegendPanelColors` の 1 件（base で失敗）を追加した。乱歩の seed 40655413787 も修正後に成功する。work guard の行は変わらない。
- work guard（`tests/web/helpers/live-work-allowlist.mjs` の `WORK_ALLOWLIST`）: setup が CDS 行に stroke を付けた後の fill を表示する 12 行（color rule commit、palette の変更、Apply to all、Save、Legend rename など）に `strokes` の stage を記録した。

## 証拠と再現方法

- run 40634318370（dev `33d9e234c`、2026-10-28）。失敗ログは `/home/kawato/gbdraw-baselines/release-0143-20261027/logs/tier-40634318370-failed.log`。
- 再現: `GBDRAW_RANDOM_WALK_SEED=40634318370 npx playwright test --config=tests/web/playwright/promotion.config.js --workers=1 -g "Linear, one Result"`。
- run 40648501950（dev `e7b3c3a06`、2026-10-28）。失敗ログは `/home/kawato/gbdraw-baselines/release-0143-20261027/logs/tier-40648501950-failed.log`。再現は上と同じコマンドで seed 40648501950。
- run 40655413787（dev `639de764a`、2026-10-28）。失敗ログは `/home/kawato/gbdraw-baselines/release-0143-20261027/logs/tier-40655413787-failed.log`。再現は seed 40655413787 で `-g "Circular, two-Result batch"`。
- 調査メモ: `/home/kawato/gbdraw-baselines/release-0143-20261027/handoffs/OV-1985.md`、`/home/kawato/gbdraw-baselines/release-0143-20261027/handoffs/OV-1987.md`、`/home/kawato/gbdraw-baselines/release-0143-20261027/handoffs/OV-1988.md`。
