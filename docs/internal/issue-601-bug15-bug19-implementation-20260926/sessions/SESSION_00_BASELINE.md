# INSTRUCTION PROMPT — S00 — 基準状態・入力入口・重複作業の確認

あなたはIssue #601のBUG-15/19実装の開始条件を確認する担当です。
runtimeを編集せず、失敗経路とregex入口、他branchとの重複、authorityの現状を調査し、実装者が使うinventoryを作成してください。

## 必ず当該ブランチを取得する

SESSION_CODEはs00です。
[総合計画書](../MASTER_PLAN.md)の共通取得手順で**origin/fix/issue-601-bug15-bug19を取得し、
専用cloneのfix/issue-601-bug15-bug19で作業してください**。
共有treeや他sessionのbranch/index/環境を使わず、同じremote branchのwriterを一つにします。
AGENTS.md、CLAUDE.md、Web CLAUDE、総合計画、
[診断公開の承認](../decisions/DECISION_01_ERROR_DISCLOSURE.md)、
[field回復の承認](../decisions/DECISION_02_REGEX_EDIT_RECOVERY.md)、先行SESSION_XX_RESULT.mdを読んでください。

## 所有範囲

- 本計画directoryのBASELINE_EVIDENCE.md、必要なsynthetic evidence、SESSION_00_RESULT.md。
- production、tests、Product authority、別Issueの計画・branchは変更しない。

## 作業

1. 取得した実装branchの開始SHA、最新origin/dev SHA、対象sourceの差分を記録する。
   既存証拠が同じcode/input/環境を保っているなら再利用し、無条件に全suiteを繰り返さない。
2. error-normalization、Python render/helper wrapper、Worker/client、run-analysisの全error return、
   Align Apply、rule action、label/color import、exportの失敗を入口から表示まで追う。
   known cause、提供中の修正情報、実際のstage、安全なcontext、次actionを一覧にする。
   errorLogを設定して原因を返さない箇所と、任意failureをInvalid ruleと呼ぶ箇所を特定する。
3. Color/Labelのmanual/TSV/preset/whitelist/visibility、Feature Search、standaloneの各regex入口を
   evaluatorとcommitまで追う。Python構文をJSでcompileする入口があれば特定し、
   検索の正しいJS compileと混同しない。
4. (?i)NADH、(?P<enzyme>NADH)、不正[、JS named group、Unicode/empty catalogを確認する。
   原監査の画面・pattern・buildを取得できなければ未特定と記録する。
   Python Color/Labelの成功を原報告のclosure evidenceにしない。
5. origin/fix/issue-601-export-output-20260926とIssue #598/602の関連branch/resultを読み、
   本計画と共有するfile、既に実装・dev統合された差分、作業中writerを確認する。
   BUG-15/19のowner移管が確認できなければ重複runtime編集を開始可能としない。
   他branchや他セッションの未コミットfileを変更しない。
6. 最新devのOIPC、Product map/BD store、診断公開の同じ意味を持つ別concernを比較する。
   二つのA receiptと既存条件の一致/差分を示し、存在しないPD/BD番号を予約しない。
7. S01–S05で共有fileを直列に引き渡すowner/path inventoryと、必要な回帰検査を記録する。
   testで生成原因をstubしている場所は実経路の検証とは区別する。

## 受入と終了成果物

SESSION_00_RESULT.mdに、開始SHA、source/authority inventory、重複ownerの引継ぎ状態、
known failure移行表、regex入口表、実行command/result、未確認原監査、S01開始条件を保存する。
sourceが変わらない既存41 Node/9 Python/2 Chromiumの結果はbaselineであり、実装後合格ではない。

English commit title: Record error and regex implementation prerequisites
Summary: Inventory failure paths, regex owners, authority, and overlapping work.

## セッション終了時のコミット・作業ブランチとdevへのプッシュ

総合計画書の共通終了手順と「S00の追加終了条件」に従い、SESSION_00_RESULT.mdと対象変更を
**検証後にコミットし、fix/issue-601-bug15-bug19へプッシュしたうえで、
最新devへ統合してorigin/devまでプッシュしてください**。
S00の計画資料・調査結果のdev統合とpushは利用者の明示的な許可済みです。
S01以降のauthority/runtimeを同時に統合する許可には拡張しません。

専用cloneの作業ブランチで最新origin/devを通常mergeし、devとの差分が本計画directoryの資料だけであること、
branch/upstream、staged scope、remote実状態を確認してください。
統合後の参照・整合性・whitespaceと適用されるrequired gateが合格したHEADを、次の順にpushします。

~~~bash
git push origin HEAD:refs/heads/fix/issue-601-bug15-bug19
git push origin HEAD:refs/heads/dev
git ls-remote --heads origin fix/issue-601-bug15-bug19 dev
git rev-parse HEAD
~~~

両remote headとlocal HEADの一致を確認してください。
dev pushが競合した場合はremote状態を確認し、最新devの通常merge・必要な再検証・作業ブランチpushを
行ってから再試行します。force-push、mainへのpush、対象外変更の統合は行いません。
commit SHA、result file、検証結果、dev統合完了、完了/未完了のS01開始条件を次sessionへ引き継ぎます。
devへのpushが未完了なら、S00を完了扱いにしないでください。
