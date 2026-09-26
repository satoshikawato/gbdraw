# 計画準備の検証記録

日付: 2026-09-26。対象は Issue #597 の実装計画と引継ぎ資料であり、runtime の修正は未実施。

- Branch: `fix/issue-597-input-session-20260926`。
- 作成時 base: `origin/dev` = `d457b7189b137185a8dec800819a312c30b969fa`。
- shared checkout を切り替えず、独立 worktree で作成した。
- 総合計画、共通作業規約、S00〜S08 の担当範囲・前提・受け入れ条件・commit/push 指示を確認した。
- 相対 Markdown links の参照先、code fences、trailing whitespace を検査した。
- 3件の承認済み human receipts と inert JSON の9項目の一致を検査した。署名は生成していない。
- 全 bash code blocks を `bash -n`、基準観察 script を Python AST で構文検査した。
- `git diff --cached --check` を実行した。
- `node tools/check-web-change-budget.mjs --base origin/dev` は Gate PASS / Review CLEAR。production、guard、active authority の差分は0。

基準観察の再利用可能な範囲は [evidence/README.md](../evidence/README.md) に記録した。
大規模性能、source collection、import Worker の実装後検証は未実施で、各セッションの実施対象。
この準備記録は S00〜S08 の検証完了を意味しない。
