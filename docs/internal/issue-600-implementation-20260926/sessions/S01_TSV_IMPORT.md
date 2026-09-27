# INSTRUCTION PROMPT — S01: Annotation TSV の付加列受入

あなたは gbdraw issue #600 の **S01** 担当です。この prompt の担当範囲を実装・検証し、
指定ブランチへ commit/push するところまで完了してください。他の session の作業は開始しない。
会話ログや未公開の提案を前提にせず、下記の committed files を入力に使ってください。

## Checkout と前提

他セッションの shared checkout を切り替えず、独立 clone に実装ブランチを取得する。
remote writer はこの session だけであることを確認し、前 session の pushed report があれば読む。

```bash
SESSION_DIR="$(mktemp -d /tmp/gbdraw-issue600-s01-XXXXXX)"
git clone --no-checkout https://github.com/satoshikawato/gbdraw.git "$SESSION_DIR"
cd "$SESSION_DIR"
git fetch origin
git switch --no-track -c fix/issue-600-annotations-styles-20260926 origin/fix/issue-600-annotations-styles-20260926
git status --short --branch
git rev-parse HEAD
```

この clone の docs を読み、次の内容を開始前に確認する。

- [総合計画](../MASTER_PLAN.md): scope、architecture、session 順、Git 運用、acceptance。
- [承認済み Product Decisions](../APPROVED_PRODUCT_DECISIONS.md): 完全な採用本文と署名者。
- `CLAUDE.md`、`gbdraw/web/CLAUDE.md`、`docs/internal/PRODUCT_IMPACT_RATCHET.md`、
  `ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`、`WEB_CHANGE_POLICY.md`。
- 当該 clone にある nearest `AGENTS.md` があれば、それも読む。

ブランチは最新 dev から作成された `fix/issue-600-annotations-styles-20260926` を継続利用する。
新たに dev から作り直して前 session の成果を捨てない。承認済み４ outcome は選び直さない。
runtime に新しい Product outcome が必要になった場合は、その範囲だけ evidence/Decision
を用意し、選択を推測しない。無関係な編集・process・server・tmp files に触れない。

## 開始条件・所有範囲

S00 report を読み、４ outcome の durable records が `origin/dev` にあることを確認する。
work branch の未 merge candidate だけで runtime を開始しない。authority 済み dev の必要な
commits を work branch に取り込み、その SHA を記録する。

担当: `gbdraw/annotations/io.py`、`gbdraw/web/js/app/annotations/table-codec.js`、
`annotations.js`、import notice の最小 markup、対応 tests と S01 report。
ほかの３修正を先に実装しない。総合計画 §4 と **TSV-01～03** に従う。

## 実装

1. BOM/trim 後の required/duplicate header を検証し、既知 annotation 列だけを射影する。
   pandas が header を mangle する前に duplicate を検出。付加列と row-width mismatch を
   区別する。known enum、座標、style、target、ID の既存検証を弱めない。
2. unknown header は捨て、１表につき集約した列名＋Session/TSV に保存しない旨を知らせる。
   cell contents はログに出さない。optional typo も通知し、fuzzy repair はしない。
3. 既存 JS array/Python tuple 戻り値を維持し、codec 内の１処理を wrapper/通知結果が共有する。
   notification のために Worker を起動しない。unknown metadata を schema に追加しない。
4. Web は complete validation 後に owner の `replaceSets` を１回実行する。known-invalid、
   file-read failure/stale のときは直前 draft/Result を保持する。
5. focused tests と実際の Web file import、通知の keyboard/status/390px、TSV download と
   Session round trip を確認する。通知と draft が次の unrelated import に混ざらないこと。

## Focused verification

```bash
pytest tests/test_annotations.py -q
node --test tests/web/annotations.test.mjs
```

notes/gene_desc/pmid、unknown optional typo、missing required、duplicate/BOM/trim、
malformed row、known-invalid の JS/Python parity を実証する。既存
`tests/web/annotation-download.playwright.spec.js` 等へ必要な browser assertions を統合。
未知列だけで strict rejection を期待する旧 tests を採用 authority に沿って更新し、
known-invalid の拒否証拠は残す。import notification と resolved Result warning を
１つの global owner にまとめない。

## 完了と publication

1. 担当 acceptance の実装・focused checks・必要な修正を完了し、scope の旧経路を除去する。
2. `docs/internal/issue-600-implementation-20260926/SESSION_RESULTS/S01.md` に start
   work SHA、authority base SHA、scope、owner/path before→after、実行 command/exit/result、
   artifacts、未実行範囲、remaining boundary を記す。未完了を completion と呼ばない。
3. production/test/docs/generated diff を別々に review し、`git diff --check` を実行する。
   command、事実、結論を区別し、過去 session の証拠は変更のない条件で再利用する。
4. `git fetch origin` と remote work SHA を確認。他 writer が進めていたら remote の実状態を
   調べてから統合方針を決める。force-push で他人の commits を消さない。
5. 現 branch と upstream を確認し、担当 path だけを `git add -- <明示 path>` で stage。
   session を１ commit とし、英語 title と短い summary を使う。**commit と push を実施して
   から session を終了する。**

```bash
git branch --show-current
git for-each-ref --format='%(upstream)' refs/heads/fix/issue-600-annotations-styles-20260926
git diff --cached --stat
git commit -m "fix: allow auxiliary columns in annotation tables"
git push --set-upstream origin HEAD:refs/heads/fix/issue-600-annotations-styles-20260926
git rev-parse HEAD
git ls-remote --heads origin refs/heads/fix/issue-600-annotations-styles-20260926
```

branch は指定名、upstream は未設定か同名 remote に限定。local HEAD と remote SHA の
一致を確認し、最終 handoff に pushed SHA、summary、checks、remaining boundary を記す。
push error が出たら先に実際の remote state を調べ、成功済み mutation を再実行しない。
main/dev への直接 push、force-push、他 branch への push、未承認の PR publication/merge/
deploy/tag は行わない。次 session の実装には着手しない。
