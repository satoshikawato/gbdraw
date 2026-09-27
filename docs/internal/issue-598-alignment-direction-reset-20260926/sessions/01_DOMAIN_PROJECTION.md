# INSTRUCTION PROMPT: S01 — 表示方向とreference中心の投影

Author: `satoshikawato`

Issue #598のBUG-03（方向を戻せないReset）とBUG-04（単一のreference方向checkbox）を修正する実装セッションです。指定の責務を最後まで実装・検証し、同名実装ブランチへcommit/pushしてください。

## 開始手順：指定ブランチを専用checkoutで取得する

他セッションの共有checkoutは操作しない。同じremote実装ブランチのwriterは一度に一セッションとし、前セッションのremote SHAと完了証拠を確認してから着手する。このセッションの成果は指定ブランチへ必ずcommit/pushする。

```bash
set -e
task_checkout=$(mktemp -d /tmp/gbdraw-issue598-S01-XXXXXX)
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

S00の完了証拠と公開format判断、4決定を含むdev ancestorが必須。Pythonの `layout/similarity_alignment.py`、`api/record_planning.py`、既存Web adapter、`layout/linear_multi_record.py` / `diagrams/linear/assemble.py` の必要最小限のfacts抽出を所有する。JS側は `app/similarity-alignment.js` のpure projection、必要なら同ownerのprivate helperを扱う。controller確定経路・Session receipt保存・操作UIはS02で統合する。Python domain判定をJSへ複製しない。

1. exact referenceを含む選択済みknown-strand対象から、tagged-union intentを現在absolute orientationへ投影する単一pure resolverを設計・実装する。Python factsはsource/selector/anchor bindingを伴う。Keep / all-right / all-left / Customをpreviewと最終Applyに共通の結果shapeへ投影する。
2. reference少数派、sourceの+/-、既にreverseのrecord、unknown/Skip/missing/unusableを扱う。bulkは参加anchorsだけに作用する。Custom以外のrow overrideと多数派推測を作らない。
3. 総合計画書4.3の中心固定式を既存placement/source projectionで実装する。現在geometryのkeyed metadataを使う場合は既存composition bridgeを通す。screen座標、DOM row順、独自crop数学をauthorityとしない。
4. reference中心のbefore/after logical canvas x、actual Δx、absolute orientations、理由付きcoverageのpure結果を出す。planは方向非依存のまま。source bytesやfeature strandを変更しない。
5. previewと最終validationで比較する出力signatureを必要な方向・対象・reference中心補正に限定する。見た目だけの差で無限に再Applyを要求しない。binding変化は明示的staleにする。
6. 意味あるunit/typed/geometry testsでA01〜A05を検証する。不等長、非ゼロbase、region/crop、複合feature、反転reference、unknown、通常Reverse、idempotent整列を含める。単に実装式を複写したexpectedではなく、期待する矢印・中心・source不変を検証する。
7. S02が呼べる一つのAPIと結果shapeをevidence/S01.mdに記す。機能途中のUIや二つ目のvalidation pipelineを公開しない。既存Matchを利用するshipping経路の置換はS02の原子的統合で完了させる。

必要な事実を既存Worker responseへ追加する場合は同じtyped decoder/admissionを拡張し、別Worker operationを作らない。最小のowner/pathを保ち、汎用direction frameworkを追加しない。

## 完了手順：検証、コミット、プッシュ

総合計画書の当該A条件を検証し、`evidence/S01.md` にowner/path差分、削除した旧経路、fixture、再実行command、環境、結果、未解決事項を記す。本文はこの会話を知らない読者に通じるように書く。未実行をPASSとしない。通常の非増加変更は簡潔なowner/path証拠を使い、policyで求められる例外だけOE/PE/CB集合を作る。

適切なfocused testsとpolicy gatesを完了し、production・tests・docs・generated diffsを分けてreviewする。ブラウザ検証が必要ならNodeの `@playwright/test` とPython Playwrightの両方を確認する。Node packageがなければPythonを使う。sandboxによるChromium起動失敗は同じ検証を適切なescalationで再実行する。テストを30分未満でtimeout扱いにしない。実行中は増分監視し、remote CI pollは5分以上間隔とする。

```bash
set -e
git diff --check
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{u}'
# reviewした対象pathを明示してstageする。git add -Aは使わない。
git add <reviewed-paths>
git diff --cached --check
git commit --author='satoshikawato <kawato@kaiyodai.ac.jp>' -m "Define alignment display-direction and reference-center projection"
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
