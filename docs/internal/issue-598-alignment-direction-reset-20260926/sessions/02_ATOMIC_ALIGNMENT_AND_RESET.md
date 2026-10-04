# INSTRUCTION PROMPT: S02 — Align・Reset・receiptの原子的統合

Author: `satoshikawato`

Issue #598のBUG-03（方向を戻せないReset）とBUG-04（単一のreference方向checkbox）を修正する実装セッションです。指定の責務を最後まで実装・検証し、同名実装ブランチへcommit/pushしてください。

## 開始手順：指定ブランチを専用checkoutで取得する

他セッションの共有checkoutは操作しない。同じremote実装ブランチのwriterは一度に一セッションとし、前セッションのremote SHAと完了証拠を確認してから着手する。このセッションの成果は指定ブランチへ必ずcommit/pushする。

```bash
set -e
task_checkout=$(mktemp -d /tmp/gbdraw-issue598-S02-XXXXXX)
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

## 前提と所有範囲

S00/S01の完了remote HEADを取得する。このセッションは機能として不可分なApply / Reset / receipt / Save・Load / History統合を所有する。主な箇所は `app/similarity-alignment.js`、`app/run-analysis.js`、既存artifact capture/restore、`services/config.js`、`services/history-snapshot.js`、`state.js`、`app/app-setup.js`、必要なcanonical serializer/admissionである。`index.html` に機能するdirection/reset選択UIを同時に実装し、古いMatch checkboxを壊れた状態で残さない。S03はそのUIを実操作・accessibilityで仕上げる。

1. 確定済みcanonical artifactを起点にS01 resolverのorientations/translationとtyped planを適用し、既存 `runCommittedCanonicalCandidate` へ渡す。未適用formを採用する専用runAnalysis pathを除去する。candidate execution、比較再利用、SVG admission、preview readiness、History rollbackを再実装しない。
2. Apply試行のPython batch validationは一度。local mode/candidate/Skip編集はWorkerゼロ。final結果で方向・reference中心補正が変われば、新previewを示して別Applyを待つ。検証済み結果の同一binding再利用も守る。
3. 成功artifactの実際のbefore/afterからcompact reset receiptを作る。active plan/source binding、反転したrecordの絶対before/after、必要なreference Δxのみ。base translationsの全量コピーやdirection policyは保存しない。plan・receipt・方向・位置・Result・Historyを一つのtransactionで確定する。
4. `Reset positions` と `Reset positions and alignment direction changes` を同じcanonical transactionに通す。positions-onlyは現在方向を保持してΔxを取り消しplanをclear。combinedは実際に反転した対象だけ絶対beforeへ戻す。いずれもreceiptを消費する。resetで新しいLOSAT jobを起動しない。
5. combined previewにnames/count/current→restored/later-manual-edit置換を表示する。empty modern deltaとmissing historical evidenceを区別し、それぞれ理由付きでcombinedを無効化する。positions-onlyは使える。壊れた新receiptを「古いreceipt欠落」扱いにしない。
6. receiptの保存・fresh Load・History capture/restoreを既存ownerへ追加する。公開format証拠に基づき必要なら既存Session versionで分岐を一度設ける。dev-only形にはmigrationを増やさない。request8/plan2へpolicy fieldを追加しない。
7. style/reorder/通常Reverse/plan invalidation/clear/新Align/reset/failureで総合計画書4.4のライフサイクルを実装する。manual位置操作によるplan clear時も既存materialization bridgeで必要な座標を保持する。source/selector binding不一致は明示的に拒否し、一部復元で成功扱いにしない。
8. 最小限の完成したradio/custom UI、Reset scope UIを接続する。Match flag・旧projection・古いtest期待値・旧文言を当該経路から同じ変更で除去する。自動Keep・ambiguity review・underlying error/retryが動く状態をcommitする。
9. A06〜A12とA14をfocused testsで実証する。Align A→B、reference反転、affected/unaffectedの後続Reverse、old/current Session fresh Load、empty delta、positions-only消費後Undo、style/reorder、pending form、全failure境界とHistory-entry数、Reset LOSAT数、ribbon geometryを含める。

Session/Historyが使うreceipt検証は同じ既存admission境界で共有する。Rendererにreceipt解釈を追加しない。commitIntentとrollback checkpointを既存transaction APIで扱い、成功後の別state書き込みで部分確定させない。

## 完了手順：検証、コミット、プッシュ

総合計画書の当該A条件を検証し、`evidence/S02.md` にowner/path差分、削除した旧経路、fixture、再実行command、環境、結果、未解決事項を記す。本文はこの会話を知らない読者に通じるように書く。未実行をPASSとしない。通常の非増加変更は簡潔なowner/path証拠を使い、policyで求められる例外だけOE/PE/CB集合を作る。

適切なfocused testsとpolicy gatesを完了し、production・tests・docs・generated diffsを分けてreviewする。ブラウザ検証が必要ならNodeの `@playwright/test` とPython Playwrightの両方を確認する。Node packageがなければPythonを使う。sandboxによるChromium起動失敗は同じ検証を適切なescalationで再実行する。テストを30分未満でtimeout扱いにしない。実行中は増分監視し、remote CI pollは5分以上間隔とする。

```bash
set -e
git diff --check
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{u}'
# reviewした対象pathを明示してstageする。git add -Aは使わない。
git add <reviewed-paths>
git diff --cached --check
git commit --author='satoshikawato <kawato@kaiyodai.ac.jp>' -m "Apply alignment directions and reset through canonical artifact transactions"
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
