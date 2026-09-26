# INSTRUCTION PROMPT — S04: 純 pixel parsing の統一

あなたは gbdraw issue #600 の **S04** 担当です。この prompt の担当範囲を実装・検証し、
指定ブランチへ commit/push するところまで完了してください。他の session の作業は開始しない。
会話ログや未公開の提案を前提にせず、下記の committed files を入力に使ってください。

## Checkout と前提

他セッションの shared checkout を切り替えず、独立 clone に実装ブランチを取得する。
remote writer はこの session だけであることを確認し、前 session の pushed report があれば読む。

```bash
SESSION_DIR="$(mktemp -d /tmp/gbdraw-issue600-s04-XXXXXX)"
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

S00 authority と S01～03 の pushed completion を確認する。
総合計画 §7 と **PX-01～03** を満たす。
担当: Web `track-slot-validation.js`、`linear-track-slots.js`、`circular-track-slots.js`、
Python `tracks/scalars.py` の純 pixel helper と Circular/Linear text adapters、tests。

## 実装

1. 対象 text fields は Linear height/spacing、Circular inner_gap_px/outer_gap_px。
   finite decimal/exponent＋optional case-insensitive px（前後空白可）を受理。
   omitted/null/trim空欄は auto、height>0、gap/spacing>=0。pure-px 文法だけを所有する。
2. `px`、0x10、bool/array/object、NaN/Infinity、10%、10em、負の geometry を拒否。
   編集中 invalid text は draft に保持して row error、submission は fatal。
   normalizer が invalid→null/0 を返して error を隠す経路を除去する。
3. shared Web helper は validator/normalizer/payload が同じ関数を使う。既存 import 方向を
   保ち、cycle を作らない。`parseOptionalLinearPx`/`parseOptionalCircularGap` と同義
   normalizer/payload parsing を同じ変更で置換する。
4. Python text adapter は同じ grammar fixtures を使い、typed boundary の手前で変換。
   Circular canonical gap は数値、Linear は既存 {value:number,unit:'px'}。
   typed JSON reader を任意 text 許可にせず、Circular radius/width の factor/% を保つ。
5. disabled draft/Save/Load と現行 schema、auto/reserved geometry を維持し、retired alias
   reader を増やさない。font-size 等の全 dimension framework に scope を広げない。

## Focused verification

```bash
pytest tests/test_circular_track_slots.py tests/test_linear_track_slots.py -q
node --test tests/web/track-slot-validation.test.mjs tests/web/circular-track-slots.test.mjs tests/web/session-request.test.mjs
```

10 / 10px / 10PX / 10 px / .5px / 1e1px の同値を raw validation、draft、payload、
CLI slot/TSV で確認。0、auto、負数、unit-only、invalid types、factor fields を検証。
browser で両 mode の対象欄へ入力し、error/Generate/Save/Load/download の実際の値を
確認する。unchanged valid geometry の SVG は不変。

## 完了と publication

1. 担当 acceptance の実装・focused checks・必要な修正を完了し、scope の旧経路を除去する。
2. `docs/internal/issue-600-implementation-20260926/SESSION_RESULTS/S04.md` に start
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
git commit -m "fix: unify pixel parsing for track slots"
git push --set-upstream origin HEAD:refs/heads/fix/issue-600-annotations-styles-20260926
git rev-parse HEAD
git ls-remote --heads origin refs/heads/fix/issue-600-annotations-styles-20260926
```

branch は指定名、upstream は未設定か同名 remote に限定。local HEAD と remote SHA の
一致を確認し、最終 handoff に pushed SHA、summary、checks、remaining boundary を記す。
push error が出たら先に実際の remote state を調べ、成功済み mutation を再実行しない。
main/dev への直接 push、force-push、他 branch への push、未承認の PR publication/merge/
deploy/tag は行わない。次 session の実装には着手しない。
