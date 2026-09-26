# INSTRUCTION PROMPT: S03 — レビューUIと操作性の実証

Author: `satoshikawato`

Issue #598のBUG-03（方向を戻せないReset）とBUG-04（単一のreference方向checkbox）を修正する実装セッションです。指定の責務を最後まで実装・検証し、同名実装ブランチへcommit/pushしてください。

## 開始手順：指定ブランチを専用checkoutで取得する

他セッションの共有checkoutは操作しない。同じremote実装ブランチのwriterは一度に一セッションとし、前セッションのremote SHAと完了証拠を確認してから着手する。このセッションの成果は指定ブランチへ必ずcommit/pushする。

```bash
set -e
task_checkout=$(mktemp -d /tmp/gbdraw-issue598-S03-XXXXXX)
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

S02の機能統合commitとevidenceを取得する。`gbdraw/web/index.html`、`app/similarity-alignment.js` のview modelと必要な既存UI wiring、関連browser testsを所有する。新たなdirection resolver、domain資格判定、History経路は作らない。S02の完成機能をdesktop/390px/keyboardで仕上げる。

1. 一つのradio groupにKeep / all-right → / all-left ← / Customを配置する。scopeは選択anchorsとreferenceを含むknown方向だけと説明する。CustomでのみrowごとのKeep/right/leftを出す。unknown/Skip対象は理由付きで適用不可とし、previewと実際の結果を一致させる。
2. referenceは正確なanchor identityとrecord名で示す。minority referenceのみの反転をall-right/all-leftで実演し、record全体が反転することとfeature center固定・record左端移動を明示する。source +/-を書き換える印象を与えない。
3. Reset scopeは一択、既定positions-only。combinedの対象names/count/current→restoredと後続manual-edit置換を表示する。receipt欠落とempty差分の理由を区別する。両scopeのconsume後にcombinedへ切替可能に見せない。
4. local candidate/Skip/mode/Custom変更でWorker呼び出しゼロを確認し、最終facts変化時の再previewと別Apply、failure後の保持draftとunderlying error、Cancel/stale回復を実操作で確認する。自動resolved Keepにreviewを強制しない。
5. semantic HTML、radio label、矢印のtext代替、focus順・復帰、keyboard選択とApply/Cancel、390pxのCustom行とReset対象一覧を検証する。既存PD-OI-035のpalette coverageを保持し、狭幅レイアウト全体のredesignを混ぜない。
6. `tests/web/similarity-alignment-ui.playwright.spec.js` と関連history/palette browser testsでA02/A03/A06/A08/A11/A13を実証する。Node package unavailableならPython Playwright scriptで同じ受入条件を検証し、再実行scriptを保存する。
7. 実際のブラウザ画像を読める倍率で目視する。internal QA画像はpublic showcaseとして代用しない。機能が完成していない部分を説明文で隠さない。evidence/S03.mdには環境・viewport・操作・結果・現行契約上の限界を書く。

軽微な文言/layout変更だけを写す無意味なunit testを増やさず、選択の排他性・scope・keyboard・worker job数など意味のある境界を保護する。

## 完了手順：検証、コミット、プッシュ

総合計画書の当該A条件を検証し、`evidence/S03.md` にowner/path差分、削除した旧経路、fixture、再実行command、環境、結果、未解決事項を記す。本文はこの会話を知らない読者に通じるように書く。未実行をPASSとしない。通常の非増加変更は簡潔なowner/path証拠を使い、policyで求められる例外だけOE/PE/CB集合を作る。

適切なfocused testsとpolicy gatesを完了し、production・tests・docs・generated diffsを分けてreviewする。ブラウザ検証が必要ならNodeの `@playwright/test` とPython Playwrightの両方を確認する。Node packageがなければPythonを使う。sandboxによるChromium起動失敗は同じ検証を適切なescalationで再実行する。テストを30分未満でtimeout扱いにしない。実行中は増分監視し、remote CI pollは5分以上間隔とする。

```bash
set -e
git diff --check
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{u}'
# reviewした対象pathを明示してstageする。git add -Aは使わない。
git add <reviewed-paths>
git diff --cached --check
git commit --author='satoshikawato <kawato@kaiyodai.ac.jp>' -m "Make alignment direction and reset choices accessible"
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
