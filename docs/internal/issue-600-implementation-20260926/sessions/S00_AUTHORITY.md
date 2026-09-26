# INSTRUCTION PROMPT — S00: 承認済み仕様の durable registration

この prompt は実施済み S00 の履歴として保持する。以後の admission 作業には
使用しない。全候補の隔離 failure、docs/evidence → contract-only の統合順序と
publication 境界は [MASTER_PLAN](../MASTER_PLAN.md) と
[S00 admission report](../SESSION_RESULTS/S00_ADMISSION.md) に従う。
元 S00 の push 指示は新しい候補・branch の publication 許可ではない。

あなたは gbdraw issue #600 の **S00** 担当です。この prompt の担当範囲を実装・検証し、
指定ブランチへ commit/push するところまで完了してください。他の session の作業は開始しない。
会話ログや未公開の提案を前提にせず、下記の committed files を入力に使ってください。

## Checkout と前提

他セッションの shared checkout を切り替えず、独立 clone に実装ブランチを取得する。
remote writer はこの session だけであることを確認し、前 session の pushed report があれば読む。

```bash
SESSION_DIR="$(mktemp -d /tmp/gbdraw-issue600-s00-XXXXXX)"
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

## この session の開始条件と担当

４つの Product outcome は [承認記録](../APPROVED_PRODUCT_DECISIONS.md) の Choice A として
`satoshikawato` が 2026-09-26 に承認済み。signature/choice/rationale/risk を再質問しない。
この session の担当は authority registration とその evidence。runtime と既存 mapped
hard behavior contracts の実装・置換は対象外。

編集 owner は `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`、現在の static
supported behavior と矛盾する説明を持つ `docs/CLI_Reference.md` の該当 scope、
この計画パッケージの S00 report。未登録 concern を新 evaluator/JSON store に登録しない。

## 実施内容

1. fetch 後の `origin/dev` と work branch で、４ concern が既に登録されたかを確認する。
   同じ approval が完全に登録済みなら重複 record を作らず、その commit と一致する fields
   を report に記す。異なる approved authority があれば conflict を具体的に報告する。
2. 未登録の concern は static Product Contract に **４つの独立 record** として追記する。
   contract revision と record ID は実際の現行 inventory を読む。番号を先取りして他の
   authority session と衝突させない。名前・意味の異なる既存 concern を supersede しない。
3. 各 record に complete normative outcome、rationale、must preserve、may retire、
   accepted residual risk、owner/date、receipt、source、acceptance IDs を含める。
   [承認記録](../APPROVED_PRODUCT_DECISIONS.md) の inert JSON と全 field を照合し、
   S00 report に実際に登録した machine representation と採番対応を示す。
4. Annotation TSV の unknown rejection は Annotation にだけ例外を設ける。records/
   track 等の別 TSV の規則を広げない。Circular gaps の without-a-unit 文を pure pixel
   text grammar に限定して改訂し、typed numeric/factor semantics を維持する。
5. candidate authority で自分の runtime を承認しない。inert data/schema/link/active-record
   一意性と whitespace を検証。code を変えて test が通ったという証拠を作らない。
   新 acceptance test は実装 session が作成する予定であり、現在 PASS と称さない。
6. 適用される base policy の schema/authority check を **base checker** で実行する。
   この authority-only scope だけで要求されない runtime regression を大量に再実行しない。
   必要な public docs skill は actual workflow に適用し、対象変更を最小限にする。
7. authority-only commit と S00 report を指定 work branch へ push する。

## この session の停止境界

authority-only 候補の push を確認したら終了する。runtime は開始しない。
maintainer がこの authority-only prefix を dev へ統合した後、S01 は fetch して
`origin/dev` の４ records が receipt と完全一致することを確認する。PR/merge の実行は
この prompt だけから追加許可を推測しない。統合済みの場合は dev SHA を report に記す。
未統合の場合は「candidate pushed、base authority integration pending」と明記する。
この待機境界は outcome の再承認ではなく trusted-base admission の順序である。

## 完了と publication

1. authority registration の照合・focused validation・必要な修正を完了する。runtime を変更しない。
2. `docs/internal/issue-600-implementation-20260926/SESSION_RESULTS/S00.md` に start
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
git commit -m "docs: register approved behavior for issue 600"
git push --set-upstream origin HEAD:refs/heads/fix/issue-600-annotations-styles-20260926
git rev-parse HEAD
git ls-remote --heads origin refs/heads/fix/issue-600-annotations-styles-20260926
```

branch は指定名、upstream は未設定か同名 remote に限定。local HEAD と remote SHA の
一致を確認し、最終 handoff に pushed SHA、summary、checks、remaining boundary を記す。
push error が出たら先に実際の remote state を調べ、成功済み mutation を再実行しない。
main/dev への直接 push、force-push、他 branch への push、未承認の PR publication/merge/
deploy/tag は行わない。次 session の実装には着手しない。
