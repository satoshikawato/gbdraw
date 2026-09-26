# INSTRUCTION PROMPT S01 — 同じ図の装飾配置を Generate に継承

あなたは gbdraw の [Issue #599](https://github.com/satoshikawato/gbdraw/issues/599) BUG-05 の実装担当者である。legend、plot title、Linear scale bar を手動移動した後の Generate が位置差分を消す不具合を修正する。承認された結果は「同一対象の新しい自動配置へ手動差分を1回加算、対応不能な非ゼロは旧 Result を保持して候補公開前停止」である。

## 取得・開始条件

同じ remote branch の writer が自分だけであること、S00 の公開完了を確認する。共有作業ツリーを切り替えない。必ず次の実装ブランチを専用 clone で取得する。

```bash
task_dir=$(mktemp -d /tmp/gbdraw-issue599-s01.XXXXXX)
git clone --branch fix/issue-599-preview-layout-20260926 \
  https://github.com/satoshikawato/gbdraw.git "$task_dir/repo"
cd "$task_dir/repo"
git fetch origin
git pull --ff-only origin fix/issue-599-preview-layout-20260926
test "$(git branch --show-current)" = "fix/issue-599-preview-layout-20260926"
test "$(git rev-parse --abbrev-ref '@{upstream}')" = "origin/fix/issue-599-preview-layout-20260926"
git status --short
```

repository AGENTS/CLAUDE、Web CLAUDE、[総合計画](MASTER_PLAN.md)、[承認一覧](00_APPROVED_PRODUCT_DECISIONS.md)、[Pack 01](DECISION_01_COMPOSITION_CONTINUITY.md)、[進捗](SESSION_STATUS.md)、S00 の結果を読む。最新の Product/Architecture/Web policy を適用する。

**origin/dev に必要な承認 outcome が恒久記録されていることを、candidate 文書ではなく `git show origin/dev:<authority path>` で確認する。** 未 merge なら runtime を開始しない。入力 SHA・authority base SHA・record ID/revision を記録する。clean 実装ブランチへ origin/dev を取り込む。ff-only が ancestry の都合で不可能なら差分を確認して通常 merge を使う。公開済み branch を rebase/reset/force push しない。

## 所有する責任

`gbdraw/web/js/app/legend-layout/composition-actions.js` と必要な同責任 helper、`legend-layout.js`、`app/candidate-render.js`、`app/run-analysis.js`、`app/app-setup.js` の最小 wiring、および既存 composition/candidate/run-analysis の検査が主な所有範囲。候補の変換 hook は `services/svg-result-ingestion.js` の既存 seam を利用する。

CSS chrome、検索ドラッグ退役、Layout edit 説明は S02/S03 の責任として残す。他 issue の変更を revert しない。同じ Generate 境界に追加機能がある場合は既存 commit path に統合し、装飾と record の delta を二重適用しない。

## 実装タスク

1. Generate 前の committed request と旧 Result から、legend/title/Linear scale の小さい snapshot を取る。旧 clean SVG を正本にし、UI ref や latent delta map を新しい永続 state にしない。mounted root を使い、選択外 batch は必要な旧 SVG だけを読む。
2. mode/grouping、validated source identity/content、biological region、record identity 集合で旧新を照合する。filename、resource ID の文字列、prefix、配列順を identity にしない。既存 token/digest/catalog を再利用し、追加 genome read/hash をしない。証明できない identity は対応不能。
3. 既存 `compositionUserDeltas()` / `applyCompositionUserDeltas()` を用い、new automaticTranslation + delta を1回だけ適用する。scale は `#length_bar` の役割で特定し、primary 配列を丸ごと引き継がない。recordTranslations と alignment は既存意味を維持。
4. 小さい snapshot/transform provider を DI で接続する。既存 `prepareCandidateRenderCommit()` / `prepareReflowResultCommit()` と `transformSvg` / `callerTransforms` の sanitizer/editor mutation 後、commit 前へ統合する。通常 Generate、committed-candidate render、automatic reflow に同じ境界を適用する。
5. batch 全候補の変換を公開前に完了させる。非ゼロの target 消失、identity 変更・不明・曖昧は明示エラーと対象 Reset/Reset Layout または設定修正へ案内し、旧 Result/request/History を保持する。zero/fresh では追加確認や空 transform の MUTATING parse を追加しない。
6. async 前の snapshot と既存 lock/token/Cancel/stale/supersession guard を維持する。render/transform/bind 失敗は既存 failure isolation へ戻す。ready 後に reconcile する第二公開経路を作らない。
7. History/Session は保存済み Result の復元だけとし、差分を再加算しない。Reset は対象だけをゼロ化する。clipping/overlap を自動補正せず padding/Reset を使う。

SRP は composition/transaction/request/保存の既存 owner を維持する。OCP は既存 transform seam、LSP は成功・失敗・ready 契約、ISP/DIP は最小 provider で守る。KISS/DRY/YAGNI に従い、新 schema、汎用位置 store、drag engine、strategy registry、追加 Worker、全 History clone を作らない。

## 検証

総合計画 C01–C08、R01 の担当範囲を結果へ対応付ける。特に2回 Generate の非倍増、色/font/legend side、primary順、batch選択外・prefix・Save/Load、未知source/crop/mode、target消失、render/transform/bind失敗、Cancel/stale、zero fast path を意味のある regression test にする。C01/C03/C05/C07 の新保証は unchanged base で失敗することを disposable clone または有効な base 証拠で示す。

対象に応じて既存検査を実行・拡張する:

```bash
node --test tests/web/composition-layout.test.mjs \
  tests/web/candidate-render.test.mjs \
  tests/web/run-analysis-simple-path.test.mjs
node --test tests/web/responsiveness-guardrail.test.mjs
```

実 Circular/Linear、既存 Gallery session、実 pointer drag → Generate 2回の browser 検査を行う。Node/Python Playwright を確認し、必要な source一致 wheel を準備する。保存/復元と実current exportの座標も確認する。`observe_base.py` は不具合存在を検査する archived before 観測であり、after 合格試験に流用しない。全契約の integration の不足分は明示して S04 へ渡す。

## 完了・commit/push

`results/S01.md` に入力/検査/authority SHA、受入対応、コマンド・結果、artifact、未達項目、owner/path の before/after と旧経路撤去を記録し、進捗を更新する。production と tests と docs を別々に review する。architecture exception が発生すれば既存 ratchet の完全な evidence と必要判断を用意し、policy を緩めない。

branch/upstream を再確認して fetch。予期しない remote 更新を調べ、他者変更を上書きしない。対象パスだけを列挙して stage、`git diff --cached --check` を実行する。**実装終了後は必ず commit と push を行う。**

```bash
git commit -m "fix(web): preserve decoration offsets across regeneration"
git push origin HEAD:refs/heads/fix/issue-599-preview-layout-20260926
git rev-parse HEAD
git ls-remote --heads origin refs/heads/fix/issue-599-preview-layout-20260926
git status --short
```

local/remote SHA の一致を確認。エラー時は remote state を確認してから再試行する。main/dev直接 push、force push、PR作成・merge・deploy を自動実行しない。終了報告に公開 SHA、検証結果、未達条件、S02 開始条件、英語 commit title と short summary を含める。
