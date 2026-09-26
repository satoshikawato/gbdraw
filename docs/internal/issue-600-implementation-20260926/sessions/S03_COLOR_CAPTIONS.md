# INSTRUCTION PROMPT — S03: 同名異色 caption の共有正規化

あなたは gbdraw issue #600 の **S03** 担当です。この prompt の担当範囲を実装・検証し、
指定ブランチへ commit/push するところまで完了してください。他の session の作業は開始しない。
会話ログや未公開の提案を前提にせず、下記の committed files を入力に使ってください。

## Checkout と前提

他セッションの shared checkout を切り替えず、独立 clone に実装ブランチを取得する。
remote writer はこの session だけであることを確認し、前 session の pushed report があれば読む。

```bash
SESSION_DIR="$(mktemp -d /tmp/gbdraw-issue600-s03-XXXXXX)"
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

S00 authority と S01/S02 の pushed completion を確認する。
総合計画 §6 と **CLR-01～05** を満たす。
担当: Python `features/colors.py`、`api/prepared.py::resolve_feature_inputs`、既存 Web
Python helper/dispatch、Web `specific-color-rules.js`、`rule-matching.js`、`watchers.js`、
`feature-editor/rule-actions.js`/`color-actions.js` の準備・commit 入口、必要な既存
legend owner と tests。state/template 全体の再設計を目的にしない。

## 実装

1. １つの I/O-free Python normalizer を compile の前に使い、table/compiled rules/legend/
   actual request/provenance が同じ canonical caption を使う。Web typed helper も同じ関数。
2. 同 caption・異色は全色へ `[lowercase-normalized-hex]` を付け、同名同色は共有する。
   空 caption は凡例なし。raw caption を予約して literal suffix 衝突を避け、allocation は
   deterministic/idempotent。rule 行順・regex/precedence/visibility を変えない。
3. 同じ色表現（named color と hex 等）は dedupe。既存 `_unique_legend_key` 相当の
   allocator を再利用し、JS 側に独立した suffix generator を作らない。新生成 caption の
   衝突で無関係な manual legend を置換しない。従来の type/numeric legend 命名体系全体は
   対象外。canonical label を採用した source rule を追い、suffix 文字列の逆解析をしない。
4. file rules と retained manual rules の full candidate を preparation し、正規化した
   caption から既存 solid legend intents を作る。`reject`/`last-wins` の旧分岐と native の
   raw-caption 上書き経路を除去する。matching cache hit でも caption processing は skip しない。
5. import/manual edit/feature recolor/legend edit は同じ preparation/action を通し、
   snapshot/stale/cancel/error を検証してから canonical rules/provenance/legend/Result を
   既存 History の１ transaction で commit。非同期処理前に legend だけ mutate しない。
6. Session/TSV は admitted caption を普通の string として保存。原 file bytes は保つ。
   Load は saved preview/draft を保持して自動 Generate せず、曖昧 draft は次の edit/Generate
   の normal preparation で通知付き正規化。native request/provenance も actual rules に揃える。
7. rename/recolor/remove/sort/再 import/Undo/Redo を支持する。１ caption に複数 swatch を
   持つ renderer/editor/schema や、新 compatibility migrator を追加しない。

## Focused verification

```bash
pytest tests/test_feature_visibility.py tests/test_web_rule_matching.py tests/test_legend_measurement.py -q
node --test tests/web/file-imports.test.mjs tests/web/rule-matching.test.mjs tests/web/legend-sync.test.mjs tests/web/history.test.mjs
```

same-caption/same-color、same-caption/different-color、named/hex、blank、literal suffix/
new-generated-vs-existing legend collision、idempotency、reorder、unused rules を検証。
`source-legend-reconciliation.playwright.spec.js`、`python-rule-parity.playwright.spec.js` 等の
既存 journey を使い、native fresh と Web live/manual/fresh、dual legends、Save/Load/
Generate、stale rollback の parity を実測する。representative 図の feature と各 swatch
の色・caption・extent を視覚確認し、metadata-only テストだけで完了としない。

## 完了と publication

1. 担当 acceptance の実装・focused checks・必要な修正を完了し、scope の旧経路を除去する。
2. `docs/internal/issue-600-implementation-20260926/SESSION_RESULTS/S03.md` に start
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
git commit -m "fix: disambiguate multicolor legend captions"
git push --set-upstream origin HEAD:refs/heads/fix/issue-600-annotations-styles-20260926
git rev-parse HEAD
git ls-remote --heads origin refs/heads/fix/issue-600-annotations-styles-20260926
```

branch は指定名、upstream は未設定か同名 remote に限定。local HEAD と remote SHA の
一致を確認し、最終 handoff に pushed SHA、summary、checks、remaining boundary を記す。
push error が出たら先に実際の remote state を調べ、成功済み mutation を再実行しない。
main/dev への直接 push、force-push、他 branch への push、未承認の PR publication/merge/
deploy/tag は行わない。次 session の実装には着手しない。
