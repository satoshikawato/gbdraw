# INSTRUCTION PROMPT — S01 — 二つの承認outcomeの正式契約化

あなたは二つの承認済み製品動作を正式なbase authorityへ記録する担当です。
診断公開とColor rule field回復のreceiptを忠実にserializeし、runtimeを変更しないauthority-only候補を作成してください。

## 必ず当該ブランチを取得する

SESSION_CODEはs01です。
[総合計画書](../MASTER_PLAN.md)の共通取得手順で**origin/fix/issue-601-bug15-bug19を取得し、
専用cloneのfix/issue-601-bug15-bug19で作業してください**。
共有treeや他sessionのbranch/index/環境を使わず、同じremote branchのwriterを一つにします。
AGENTS.md、CLAUDE.md、Web CLAUDE、総合計画、
[診断公開の承認](../decisions/DECISION_01_ERROR_DISCLOSURE.md)、
[field回復の承認](../decisions/DECISION_02_REGEX_EDIT_RECOVERY.md)、先行SESSION_XX_RESULT.mdを読んでください。

## 開始条件とbranch例外

S00のpush済みresultと二つの承認記録を読んでください。
まず実装branchを取得してreceiptを確認し、authority/runtimeを分離するため専用authority branchへ移ります。
この例外は製品判断・実装scopeを変更するものではありません。

~~~bash
git fetch origin
git switch --no-track -c product/issue-601-errors-regex-authority-20260926 origin/dev
~~~

同名branchが既に存在する場合は作成を繰り返さず、remote内容と担当者の状態を確認する。
自分が確認した既存候補を継続するか、未使用branch名で進めるかを記録し、同名remoteだけへpushする。
実装branchのreceiptは以下で取得できる。

~~~bash
git show origin/fix/issue-601-bug15-bug19:docs/internal/issue-601-bug15-bug19-implementation-20260926/decisions/DECISION_01_ERROR_DISCLOSURE.md
git show origin/fix/issue-601-bug15-bug19:docs/internal/issue-601-bug15-bug19-implementation-20260926/decisions/DECISION_02_REGEX_EDIT_RECOVERY.md
~~~

OIPCへ写す本文の九項目は承認ファイルのJSONと文字列単位で比較する。

## 所有範囲

- authority branch: docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.mdのみ。
- 実装branchに戻った後: 本計画directoryのSESSION_01_RESULT.mdのみ。
- runtime、tests、detectors/rules、CI guard、map/BD storeへ同時変更を加えない。

## 作業

1. 最新devのmaintainer allowlist、OIPC revision/番号、同じsemantic concern、
   web.errors.user-facing-diagnostic-disclosureの別候補/正式契約を確認する。
   未登録Lane Bのconcernを新しいmachine-readable storeへ追加しない。
2. web.errors.diagnostic-disclosureとweb.rules.rejected-pattern-edit-recoveryのA receiptを
   個別に記録する。rationale、mustPreserve、mayRetire、acceptedResidualRisk、
   owner/dateを翻訳・要約・追加せず、scenarioRevisionとchoiceも保持する。
3. 同じ診断公開を異なるキーの競合active contractに複製しない。
   既存recordで同じoutcomeが完全に満たされていればそのauthorityを再利用し、差分を示す。
   Copy等の独立要求が未実現なら消さず、元承認にないキー変更や退役を推測しない。
   実質矛盾があればそのconcernのauthority収束を止め、具体的な差分を提示する。
4. OIPCの既存metadata形式に合わせてrevisionと利用可能なPD番号を管理する。
   他のPD、原監査の再現状態、科学的制約を変えない。未実装outcomeを実装済みとしない。
5. 自己candidate runtimeを許可せず、authority-only diffと既存checkerの要求を確認する。
   receiptのJSONと候補本文が同じ意味・全項目で一致することを検査する。
6. English title/summaryでauthority branchへcommitし、
   git push -u origin HEAD:refs/heads/product/issue-601-errors-regex-authority-20260926でpushする。
   使用名を変更した場合はその同名remoteへpushし、remote SHA一致を確認する。
7. PR作成/dev mergeは別の明示的許可が必要。許可があればPR wording skillを適用し、
   authority-onlyのレビュー/mergeを行う。直接dev pushは使わない。
   merge未完了ならAUTHORITY_PENDINGと記録し、S02開始を許可しない。
8. 元のfix/issue-601-bug15-bug19へ戻り、fetch/pull --ff-onlyしてから
   SESSION_01_RESULT.mdにauthority branch/SHA、記録PD/concern、元receipt一致、
   Gate/Review、dev統合SHAまたはpending、S02開始条件を保存する。
   この非規範handoffを別commitとして実装branchへpushする。authorityとruntimeを混ぜない。

## 受入

選択Aだけをserializeし、各receiptの九項目が欠落・改変なく保持されること。
B案や新しい残余リスク、既存権限緩和を含めない。
候補がpushされた状態と、正式authorityがdevに統合された状態を区別する。

Authority commit title: Record approved error disclosure and regex recovery behavior
Summary: Preserve the two complete product receipts in the existing static contract.
Handoff commit title: Record authority integration status for issue 601
Handoff summary: Record the authority candidate and the prerequisite for runtime work.

## セッション終了時のコミット・プッシュ

総合計画書の共通終了手順に従い、authority branchと非規範handoffをそれぞれ
**検証後にコミットし、当該同名remote branchへプッシュしてください**。
通常のtargetはfix/issue-601-bug15-bug19、authorityは上記専用branch。
branch/upstream、staged scope、remote実状態を確認し、force-pushやmain/dev直接pushは行いません。
remote/local SHA一致、result file、完了/未完了の開始条件を次sessionへ引き継ぎます。
