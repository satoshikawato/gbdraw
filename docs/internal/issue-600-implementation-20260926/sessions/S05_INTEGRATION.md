# INSTRUCTION PROMPT — S05: 統合検証と修正完了

あなたは gbdraw issue #600 の **S05** 担当です。この prompt の担当範囲を実装・検証し、
指定ブランチへ commit/push するところまで完了してください。他の session の作業は開始しない。
会話ログや未公開の提案を前提にせず、下記の committed files を入力に使ってください。

## Checkout と前提

他セッションの shared checkout を切り替えず、独立 clone に実装ブランチを取得する。
remote writer はこの session だけであることを確認し、前 session の pushed report があれば読む。

```bash
SESSION_DIR="$(mktemp -d /tmp/gbdraw-issue600-s05-XXXXXX)"
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

S00 authority が origin/dev にあり、S01～04 の pushed completion が当該 work branch に
あることを確認。未完了 acceptance を completion と見なさない。
担当は総合計画の全 acceptance の統合確認と、その failure に必要な最小 fix、
既存 public docs の適切な owner 更新、S05 report。無関係な refactor/feature を足さない。

## 実施内容

1. S00 の complete records と４ runtime outcomes を照合。canonical request / Result
   admission の mapped requirements を AND-of-OR で維持し、option ID の一致だけで
   coverage を判断しない。各 session の owner/path、superseded paths の削除を確認。
2. focused evidence を source/input/environment/criterion が不変なら再利用し、未検証の
   interaction だけを追加する。４ outcome を１ realistic diagram/session で組み合わせて
   Web Generate → live edit → Undo/Redo → Save/Load → Generate → export を確認する。
3. Circular single/grid/batch、Linear、CLI/Python API、worker helper、warning の成功/fail/
   cancel/stale の帰属を検証。privacy/sanitizer/local processing と原 file bytes を維持。
4. raw table → canonical values → actual SVG に effect が届くことを確認。未知列を保持した
   と誤認させず、selector miss を partial geometry にせず、各使用色を説明し、pixel invalid
   を auto/zero にしない。390px/keyboard/status と readable-scale 図を確認する。
5. 既存 valid input の SVG を read-only reference comparison で検証。期待されない幾何差の
   原因を修正し、reference を先に更新しない。public showcase を smoke diagram に置換しない。
6. static behavior説明は既存 CLI/technical owner に統合し、規約と矛盾を残さない。
   public procedural docs/Gallery が本当に変わる場合だけ該当 skill/owner command を適用。
   不要な新 docs page、wheel commit、cache-bust、social-preview 更新はしない。

## Required gates

総合計画 §10 の Ruff、fast Python、read-only OutputComparison、Web architecture/
Product contracts を実行。`.github/workflows/test.yml` が指定する fast Web contracts と
対象 Playwright journeys を確認し、base checker の policy gate を actual accepted-base /
runtime-head commits に対して実行する。candidate checker/authority に自己承認させない。

shared change に適用される required regression/CI impact jobs を満たす。test run は30分
未満で timeout と判断せず増分監視、remote CI は５分間隔以上。fix 後は影響 checks と
必要な regression のみを再実行。期待値/threshold を緩めて failure を隠さない。

## Integration report

`SESSION_RESULTS/S05.md` には、全 acceptance ID→test/artifact の対応、commands/results、
base/work SHAs、owner/path evidence、remaining limits、rollback をまとめる。
物理ファイル/CI gate/画面で未実行の内容を PASS と呼ばない。後続 reader が同 branch の
report だけで release/integration 可否を判断できる内容とする。
必要な修正・docs・report を１ commit にして push。PR の作成/変更が別途許可された場合は
write-clear-pull-request skill と language check を適用する。deploy/merge/tag は実行しない。

## 完了と publication

1. 担当 acceptance の実装・focused checks・必要な修正を完了し、scope の旧経路を除去する。
2. `docs/internal/issue-600-implementation-20260926/SESSION_RESULTS/S05.md` に start
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
git commit -m "test: verify annotation and style input fixes"
git push --set-upstream origin HEAD:refs/heads/fix/issue-600-annotations-styles-20260926
git rev-parse HEAD
git ls-remote --heads origin refs/heads/fix/issue-600-annotations-styles-20260926
```

branch は指定名、upstream は未設定か同名 remote に限定。local HEAD と remote SHA の
一致を確認し、最終 handoff に pushed SHA、summary、checks、remaining boundary を記す。
push error が出たら先に実際の remote state を調べ、成功済み mutation を再実行しない。
main/dev への直接 push、force-push、他 branch への push、未承認の PR publication/merge/
deploy/tag は行わない。次 session の実装には着手しない。
