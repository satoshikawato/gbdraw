# Web GUI リリース前監査 — dev 33d9e234c（2026-10-28）

対象は `origin/dev` の `33d9e234c`（0.14.3 の promotion head）。release tier の `Tests`（`workflow_dispatch`、`tier=release`）run 40634318370 の結果をまとめる。seed は run ID（40634318370）。

**件数**: 固有の問題 1 件（P2 1）。`Promotion audit / random-walk` だけが失敗し、その原因が OV-1985 である。修正は #1433（`fix/ov1985-legend-stroke-caption`）。

## 方法と表記

- 手順は `docs/internal/WEB_PERIODIC_AUDIT.md` に従う。`random-walk` は seed 付きの乱歩（`tests/web/live-generate-random-walk.promotion.spec.js`）で、各 step で live 表示と Generate 結果を比べる。
- 重大度の基準は 2026-09-30 の監査と同じ（P1 誤った出力・編集の消失・有効入力で Generate 不能、P2 機能の不具合・古い状態の残存、P3 軽微）。
- JS のパスは `gbdraw/web/js/` からの相対パス、それ以外はリポジトリ直下からのパス。行番号は書かない。

## 結果（ジョブ別）

| ジョブ | 結果 | 項目 |
|---|---|---|
| random-walk | 失敗（fixture "Linear, one Result" の step 6 と 7） | OV-1985 |
| ほかの `Promotion audit / ...` と release ジョブ | 成功 | なし |

## 見つかったもの

### P2（1 件）

| ID | 症状 | 原因 |
|---|---|---|
| OV-1985 | 同じ caption の color rule が 2 つの色で描かれている Legend 行（Python の key は `CDS [#54bcf8]`）に、feature popup の "Legend item" scope で stroke を付ける。その後、もう一方の rule の caption を変えると（`CDS` → `transfer RNA`）、Legend の CDS 行の枠線が live では既定の灰色（`#808080`）に戻り、Generate では指定色（`#123456`）で描かれる。Generate の後で caption を `CDS` に戻したときも同じく live だけ灰色になる。Generate すれば正しく描かれる | rule commit（`commitOnce` in `app/feature-editor/rule-actions.js`）は、key の変わる行の記録を新しい key に移す（`moveLegendRecords` in `app/legend/entry-actions.js`、OV-660）。stroke の編集（`legendStrokeOverrides`）も `CDS [#54bcf8]` から `CDS` に移るので、editor intent と Generate は正しい。ところが commit の live 表示（`projectPaletteAndRules` in `app/app-setup.js`）は fill と Legend 構造の domain だけを `showEditorIntent` に渡し、stroke の domain を compile しなかった。そのため stroke は隠れた古い行に残り、表示される行には付かなかった。#1042（v0.14.0）から存在し、v0.14.2 にもある（関係するコードは v0.14.2 と同一）。0.14.2 からの回帰ではない |

**修正**: commit が Legend 行の stroke の key を変えたときは、同じ compile で stroke の domain も表示する（`projectPaletteAndRules` の `strokes`）。stroke の key が変わらない commit の作業は変わらない。

## guard の判断

- OV-1985: guard を追加した。`tests/web/live-generate-parity-legend.playwright.spec.js` の 2 件（caption を変えたとき、Generate の後で caption を戻したとき。いずれも base で失敗）を PR の Playwright functional shard が実行する。work guard（`tests/web/helpers/live-work-allowlist.mjs` の `WORK_ALLOWLIST`）の 2 行（"color rule commit"、"Legend rename, rule row both Results draw"）は、stroke 付きの行の key を変える commit なので `strokes` の stage を記録した。

## 証拠と再現方法

- run 40634318370（dev `33d9e234c`、2026-10-28）。失敗ログは `/home/kawato/gbdraw-baselines/release-0143-20261027/logs/tier-40634318370-failed.log`。
- 再現: `GBDRAW_RANDOM_WALK_SEED=40634318370 npx playwright test --config=tests/web/playwright/promotion.config.js --workers=1 -g "Linear, one Result"`。
- 調査メモ: `/home/kawato/gbdraw-baselines/release-0143-20261027/handoffs/OV-1985.md`。
