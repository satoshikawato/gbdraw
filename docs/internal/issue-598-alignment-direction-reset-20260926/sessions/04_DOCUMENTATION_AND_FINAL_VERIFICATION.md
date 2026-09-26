# INSTRUCTION PROMPT: S04 — ドキュメントと最終検証

Author: `satoshikawato`

Issue #598のBUG-03（方向を戻せないReset）とBUG-04（単一のreference方向checkbox）を修正する実装セッションです。指定の責務を最後まで実装・検証し、同名実装ブランチへcommit/pushしてください。

## 開始手順：指定ブランチを専用checkoutで取得する

他セッションの共有checkoutは操作しない。同じremote実装ブランチのwriterは一度に一セッションとし、前セッションのremote SHAと完了証拠を確認してから着手する。このセッションの成果は指定ブランチへ必ずcommit/pushする。

```bash
set -e
task_checkout=$(mktemp -d /tmp/gbdraw-issue598-S04-XXXXXX)
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

S03までの完成remote HEADと全セッションevidenceを取得する。ユーザードキュメント、必要なGallery説明・操作画像、全体のfocused regression修正、最終handoffを所有する。既存policy/authority/checkerを機能に合わせて弱めない。BUG-17の新機能・資料は追加しない。

1. public user manual/procedural docsを変更する場合は `.agents/skills/love-me-love-my-docs/SKILL.md` を読み、その使用を通知する。Gallery screenshots/tutorial JSON/operation registerを変更する場合は `.agents/skills/web-gallery-screenshot-maintenance/SKILL.md` も読み、実GUIから再現可能に生成する。internal planの存在だけでskillsを適用しない。
2. Match reference direction / Reset Alignの説明を完成仕様へ更新する。minority reference、all-right/all-left、Custom、positions-only/combined、後続manual方向編集、receiptのない旧Session、Undo、source不変を利用者の操作として説明する。会話を知らない読者に通じない「前案」「それ」等の文章を残さない。
3. 既存の現実的なGallery-quality recipe/sessionを起点に、必要な例だけ更新する。labels、legend、colors、metadata、tracks、comparison contextを保つ。minimal smoke図をpublic例に使わない。最終図を目視し、記載recipeから再生成する。`examples/gbdraw_social_preview.png` は編集しない。
4. A01〜A16の実証を照合し、欠けた条件だけ補う。Session fresh Load、History/readiness/finalization rollback、pending form、Reset LOSATゼロ、reference center、ribbon geometry、円形Rotate/既存linear regionを最終HEADで確認する。
5. 関連Node/Playwright tests、`pytest tests/ -v -m "not slow"`、`ruff check gbdraw/`、必要なoutput comparison、適切なarchitecture/Product Gateを実行する。`node tools/check-web-change-budget.mjs --base origin/dev --head HEAD` のGate/Review本文も確認する。policy profileの追加条件があればその根拠と証拠を保存する。
6. `tests/reference_outputs/` はread-onlyのcomparisonを先に行う。intentional geometry差分だけreview後 `pytest tests/test_output_comparison.py::TestGenerateReferences --update-reference-outputs -v` で更新し、SVG diff目視とcomparison再実行を行う。既存Keep出力が無意図に変わる場合はbaseline更新で隠さず修正する。
7. 最終production/tests/docs/generated diffsを分けて監査する。旧Match state/path/docs、第二owner、dev-only migration、未使用helper、不要version/registry、誤ったscopeを除去する。関連required checksが通った後は新たな根拠なく反復しない。
8. evidence/S04.mdと `02_FINAL_HANDOFF.md` を保存する。handoffに完成behavior、4 accepted decisions、runtime基底authority SHA、全A条件の証拠、実行command/環境、未解決リスク、ブランチと最終commitの確認手順、英語commit title/summaryを記す。実際に未実施の条件を実施済みにしない。

最終handoffをcommit/pushしてremote SHAを確認する。PRを別途作成・更新する場合だけ `.agents/skills/write-clear-pull-request/SKILL.md` とlanguage gateを適用する。この指示自体はPR公開・merge・releaseを許可しない。

## 完了手順：検証、コミット、プッシュ

総合計画書の当該A条件を検証し、`evidence/S04.md` にowner/path差分、削除した旧経路、fixture、再実行command、環境、結果、未解決事項を記す。本文はこの会話を知らない読者に通じるように書く。未実行をPASSとしない。通常の非増加変更は簡潔なowner/path証拠を使い、policyで求められる例外だけOE/PE/CB集合を作る。

適切なfocused testsとpolicy gatesを完了し、production・tests・docs・generated diffsを分けてreviewする。ブラウザ検証が必要ならNodeの `@playwright/test` とPython Playwrightの両方を確認する。Node packageがなければPythonを使う。sandboxによるChromium起動失敗は同じ検証を適切なescalationで再実行する。テストを30分未満でtimeout扱いにしない。実行中は増分監視し、remote CI pollは5分以上間隔とする。

```bash
set -e
git diff --check
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{u}'
# reviewした対象pathを明示してstageする。git add -Aは使わない。
git add <reviewed-paths>
git diff --cached --check
git commit --author='satoshikawato <kawato@kaiyodai.ac.jp>' -m "Document and verify alignment direction and reset behavior"
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
