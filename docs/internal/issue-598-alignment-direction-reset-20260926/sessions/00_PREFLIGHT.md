# INSTRUCTION PROMPT: S00 — 実装開始条件とbaselineの確定

Author: `satoshikawato`

Issue #598のBUG-03（方向を戻せないReset）とBUG-04（単一のreference方向checkbox）を修正する実装セッションです。指定の責務を最後まで実装・検証し、同名実装ブランチへcommit/pushしてください。

## 開始手順：指定ブランチを専用checkoutで取得する

他セッションの共有checkoutは操作しない。同じremote実装ブランチのwriterは一度に一セッションとし、前セッションのremote SHAと完了証拠を確認してから着手する。このセッションの成果は指定ブランチへ必ずcommit/pushする。

```bash
set -e
task_checkout=$(mktemp -d /tmp/gbdraw-issue598-S00-XXXXXX)
git clone --branch fix/issue-598-alignment-direction-reset-20260926 --single-branch \
  https://github.com/satoshikawato/gbdraw.git "$task_checkout"
cd "$task_checkout"
git fetch origin \
  refs/heads/dev:refs/remotes/origin/dev \
  refs/heads/main:refs/remotes/origin/main \
  refs/heads/fix/issue-598-alignment-direction-reset-20260926:refs/remotes/origin/fix/issue-598-alignment-direction-reset-20260926
git status --short
git branch --show-current
git rev-parse HEAD origin/fix/issue-598-alignment-direction-reset-20260926
```

branch名と両SHAの一致、clean treeを確認する。古いdevから別のfix branchを作り直さない。既存の変更や他セッションのプロセスを破棄しない。Python environment・wheel・server portを専用checkoutに隔離し、共有editable installをしない。

`AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`、`docs/internal/PRODUCT_IMPACT_RATCHET.md`、`docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`、`docs/internal/WEB_CHANGE_POLICY.md` を読む。続いて同ディレクトリ上位の [総合計画書](../01_MASTER_IMPLEMENTATION_PLAN.md) と [承認文](../00_ACCEPTED_PRODUCT_DECISIONS.md) を読む。

`origin/dev:docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` のPD-OI-027/029/031/034がそれぞれrevision 5/3/5/5で、承認文の9フィールドすべてと一致することを確認する。candidate branch上だけの決定をauthorityとしない。S01以降は、S00が取り込んだauthority入りdev commitがHEADのancestorであることも確認する。最新devで当該決定が再更新された場合は契約を照合し、影響する実装を止めて新しいauthorityを解決する。条件未達時は独立した文書・調査だけを進め、依存するruntime変更を止める。決定更新ブランチを直接取り込んで代用しない。

BUG-17の新しい線形クロップ／回転代替は対象外。source strandを変更しない。未知方向の推測、多数派自動反転、第二render/History/decision owner、方向policyの永続化、dev-only migrationは追加しない。新しいmaterial Product outcomeだけは関心ごとのDecision Packへ分離し、routineな実装選択は既存authorityの下で進める。

## 所有範囲と目的

このセッションはruntimeを変更しない。`docs/internal/issue-598-alignment-direction-reset-20260926/evidence/S00.md` と、調査による必要最小限の計画補足を所有する。Product Contract、checker、workflow、runtime、reference outputsを編集しない。

1. 最新origin/devに4件の正確な決定があるかを確認する。なければbaseline/format調査を保存してcommit/pushし、authority先行mergeが未達と報告する。この未達はProduct Impact Ratchetのcandidate self-authorization禁止に由来する。runtime着手のためにapprovalを推定しない。
2. 決定が揃っていれば、指定fix branchで `git merge --no-edit origin/dev` を行い、承認済みauthorityを含むdevをancestorにする。競合は意味を保って解決し、他のdev変更を削除しない。authority branchを直接mergeしない。branchの作成ベースをmainに変更しない。
3. 最新mainとrelease tagsのfirst-parent公開証拠を調べ、Session44 / request8 / plan2 / catalog4の公開状態と旧Session fixturesを列挙する。tag名だけから公開を決めず、writer・admission実装・artifactを確認する。証拠のある公開済み形だけcompatibility対象にする。
4. source/display strand facts、candidate validation、linear placement、canonical target transaction、Session/History capture、plan invalidationの実際のowner/pathを追跡する。既存pure helperの再利用可否とreference Δx receiptの保存位置を確定する。record translationと第二baseline snapshotを重複させない。
5. 総合計画書A01〜A16をjourney/checkpointごとに割り当てる。AND-of-OR要件を使い、shared choice IDを証拠の代わりにしない。material outcomeの未決がなければ既存authority実装として自動進行する。
6. 関連Node tests、`node tests/web/session-request.test.mjs`、`python -m pytest tests/test_similarity_alignment_web_adapter.py -q` をbaselineとして実行する。可能なら既存alignment UI specも確認する。準備wheelは専用checkout内で必要時だけ生成する。
7. `node tools/check-web-change-budget.mjs --base origin/dev --head HEAD` のGate/Reviewを記録する。新規runtimeをまだ持たない状態の結果と明記する。環境、base/main/branch SHA、公開証拠、前提未達、S01 API指示を保存する。

独立調査はauthority未mergeでも完了できる。条件未達を隠してS01を開始可能と書かない。authority反映済みなら完了証拠とremote HEADを次セッションへ渡す。

## 完了手順：検証、コミット、プッシュ

総合計画書の当該A条件を検証し、`evidence/S00.md` にowner/path差分、削除した旧経路、fixture、再実行command、環境、結果、未解決事項を記す。本文はこの会話を知らない読者に通じるように書く。未実行をPASSとしない。通常の非増加変更は簡潔なowner/path証拠を使い、policyで求められる例外だけOE/PE/CB集合を作る。

適切なfocused testsとpolicy gatesを完了し、production・tests・docs・generated diffsを分けてreviewする。ブラウザ検証が必要ならNodeの `@playwright/test` とPython Playwrightの両方を確認する。Node packageがなければPythonを使う。sandboxによるChromium起動失敗は同じ検証を適切なescalationで再実行する。テストを30分未満でtimeout扱いにしない。実行中は増分監視し、remote CI pollは5分以上間隔とする。

```bash
set -e
git diff --check
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{u}'
# reviewした対象pathを明示してstageする。git add -Aは使わない。
git add <reviewed-paths>
git diff --cached --check
git commit --author='satoshikawato <kawato@kaiyodai.ac.jp>' -m "Record issue 598 implementation preflight and baseline"
node tools/check-web-change-budget.mjs --base origin/dev --head HEAD
# Gate本文がPASSであることとReview要件を確認してから、次へ進む。
git fetch origin \
  refs/heads/fix/issue-598-alignment-direction-reset-20260926:refs/remotes/origin/fix/issue-598-alignment-direction-reset-20260926
git merge-base --is-ancestor origin/fix/issue-598-alignment-direction-reset-20260926 HEAD
git push origin HEAD:refs/heads/fix/issue-598-alignment-direction-reset-20260926
git rev-parse HEAD
git ls-remote --heads origin fix/issue-598-alignment-direction-reset-20260926
```

branchとupstreamは同名実装branchであること。remoteが予期せず進んだ場合はpushを止めて状態を確認し、他writerの完了後に通常merge等で正しく統合する。force-push、dev/mainへのpush、他人の変更のrevertは行わない。GateがFAILならpushせず原因を修正して再検証する。Review REQUIREDはGate PASSと区別して報告し、承認を推定しない。remote SHA=HEADを確認してセッションを完了する。push失敗時は実際のremote状態を確認してから再試行する。

最終回答に実施内容、tests結果と限界、commit SHA、同名remoteへのpush確認、次セッションの開始条件、英語proposed commit title/summaryを記す。PR公開・merge・deploy・tagはこのセッションの自動実行対象ではない。
