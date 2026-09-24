# Issue #581 — S03 Synteny guide と図上選択 INSTRUCTION PROMPT

以下を一つの独立した作業指示として使用する。過去の会話は前提にしない。

## 目的と開始条件

gbdraw Issue #581 の floating alignment palette に対して、参照 feature に対応する
縦ガイド、番号付き候補、palette と図の双方向 hover／選択を追加する。利用者は元の
Linear 比較図を見たまま exact candidate を選べる。製品意味、対象外、受入条件は
`docs/internal/issue-581-feature-popup-alignment-ux/IMPLEMENTATION_PLAN.md` にある。
S01 の local draft と S02 の non-modal palette、必要な Product authority が完了してから
着手する。

## ブランチと資料

- runtime 作業は **`issue-581-feature-popup-alignment-ux-20260924`** に積む。
  branch、HEAD、upstream、worktree、総合計画第10節を確認する。
- `AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`、
  `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` の PD-OI-026〜029、031、
  `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md` と S00 の判断を読む。
- `gbdraw/web/js/app/feature-editor/svg-actions.js` の identity lookup、hover、
  preview transform interaction、`gbdraw/web/js/app/app-setup.js` の wiring、
  `gbdraw/web/js/app/similarity-alignment.js` の Select／Skip、
  `gbdraw/web/index.html` の Result preview、現行 browser tests を確認する。

## 実装

1. Resolver の exact anchor identity を既存の rendered feature lookup で照合する。
   重複 ID、分割 feature、非表示、crop、未描画を想定し、一意に対応付けられる
   rendered geometry だけを overlay 対象にする。描画対象がなくても Python が usable と
   判定した候補は palette で選択できる。
2. reference feature の実描画中心 x を通る薄い縦ガイドを preview 専用 overlay に描き、
   各候補へ resolver 順序に基づく安定した番号 badge を付ける。overlay は sanitized
   SVG Result、download、Session に入れない。ガイドは目視の補助であり、選択順位や
   Python Resolver の適格性判定には使わない。
3. pan／zoom／scroll／resize／Result 置換で overlay を更新または破棄する。
   毎フレームの全 SVG 再走査や別の geometry engine を作らず、既存 preview transform
   通知・feature lookup の範囲で更新する。palette close／Cancel／Apply／generation
   replacement では event listener、hover style、badge を確実に掃除する。
4. palette hover/focus は候補 feature を強調し、図上 hover は対応する palette 行を
   強調する。強調は選択 state ではない。前の通常 hover は後始末時に復元する。
5. badge または候補 feature の click は S01 の同じ
   `selectCandidate(recordKey, anchor)` に送る。候補 click で通常の feature popup が
   重複して開かないよう既存イベント経路で処理する。非候補 click は既存動作を保つ。
   palette の radio は常に keyboard 代替となる。

## 検証と終了

- 実ブラウザで、同一 record に複数候補がある realistic fixture を使い、badge と
  palette の番号、直接 click、radio、hover の両方向同期を確認する。
- pan／zoom／resize 後のガイドと badge 位置、見えない候補の palette 選択、
  popup click との競合、cancel／apply／Result 置換後の cleanup を確認する。
- fixture の位置を使った自動判定や候補 ranking が増えていないこと、各選択が Worker を
  呼ばないこと、export と saved Session に overlay がないことを確認する。
- focused tests、real-browser visual review、必要な lint を実行する。Node Playwright が
  なければ Python Playwright を使い、Chromium sandbox エラーなら必要な escalation で
  同じ確認を再実行する。
- production/test/visual diff と architecture owner/path を監査し、総合計画第10節へ
  変更と証拠、S04 開始可否を記す。

SOLID: SVG interaction owner が geometry と event lifecycle を持ち、controller が
choice を持つ。KISS/DRY: lookup と hover を再利用し、radio と canvas は一つの
action に収束する。YAGNI: 汎用 annotation engine や scoring を追加しない。
英語の proposed commit title と短い summary を示す。push／PR はこのセッションで
明示的に許可された場合だけ行う。
