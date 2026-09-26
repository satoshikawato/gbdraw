# INSTRUCTION PROMPT S02 — Layout edit の対象と有効化方法を示す

あなたは gbdraw の [Issue #599](https://github.com/satoshikawato/gbdraw/issues/599) BUG-06 の実装担当者である。Preview の legend/title/Linear scale を動かすには Layout edit を ON にする必要があるが、OFF の cursor と説明からその操作を理解しにくい。明示 mode を維持し、対象と有効化方法を mouse/keyboard/touch で発見できるようにする。

## 取得・開始条件

S01 の commit/push 完了と、この remote branch の writer が自分だけであることを確認する。共有ツリーを切り替えず、必ず指定実装ブランチを専用 clone に取得する。

```bash
task_dir=$(mktemp -d /tmp/gbdraw-issue599-s02.XXXXXX)
git clone --branch fix/issue-599-preview-layout-20260926 \
  https://github.com/satoshikawato/gbdraw.git "$task_dir/repo"
cd "$task_dir/repo"
git fetch origin
git pull --ff-only origin fix/issue-599-preview-layout-20260926
test "$(git branch --show-current)" = "fix/issue-599-preview-layout-20260926"
test "$(git rev-parse --abbrev-ref '@{upstream}')" = "origin/fix/issue-599-preview-layout-20260926"
git status --short
```

repository AGENTS/CLAUDE、Web CLAUDE、[総合計画](MASTER_PLAN.md)、[Pack 02](DECISION_02_LAYOUT_AFFORDANCE.md)、[承認一覧](00_APPROVED_PRODUCT_DECISIONS.md)、[進捗](SESSION_STATUS.md)、S00/S01 の結果を読む。最新 origin/dev に必要 authority が存在することを確認し、入力 SHA・authority base/record と S01 の公開 SHA を記録する。未 merge の candidate authority で runtime を承認しない。

## 承認済み動作と所有範囲

選択 `A / EXPLICIT-MODE-WITH-DISCOVERABLE-TARGETS` は2026-09-26に `satoshikawato` が完全な本文で承認済み。OFF は従来 canvas pan と明示 toggle、ON は target drag。対象への pointer 操作で mode を自動 ON にしない。

主な所有範囲は既存 legend/composition drag owner、`gbdraw/web/js/app/legend-layout/`、`app/legend/drag-actions.js`、`gbdraw/web/index.html` の Preview hint と toolbar 説明、必要な `services/svg-serialization.js` transient inventory、既存 interaction/serialization tests。Generate 継承と候補 transaction は S01 の経路を維持し、検索の row/docking は S03 が担当する。

## 実装タスク

1. 既存 supported target の eligibility を共有して legend/title/Linear scale のみに説明を出す。動かせない要素へ移動可能な説明を出さない。hover ごとの全 SVG 走査を追加しない。
2. OFF は `help` cursor、控えめな hover outline、`Turn on Layout edit to move this item`。ON は `grab`、実 drag 中は `grabbing` とする。pan cursor の親継承と装飾移動の意味を区別する。
3. toolbar に常設説明を置き、toggle の `aria-pressed` と説明の関連付けを確認する。keyboard focus と touch でも同じ操作説明を読めるようにする。多数の SVG tab stop や touch専用の別 mode owner は追加しない。
4. feature click、label edit、legend 個別編集、Shift/Ctrl、background pan、record/alignment の既存優先順位を維持する。既存 drag owner を使用し、第二 pointer router、global event bus、新 drag engine を作らない。
5. hint は Preview 専用派生表示とする。HTML wrapper の hint/CSS を優先する。SVG transient class が必要なら `stripTransientPreviewState()` の明示 inventory へ追加して検査する。prefix が自動除去されるとは仮定しない。saved Result、History、plain/interactive SVG、PNG/PDF に hint/outline/cursor を混入させない。
6. Result/load/History の既存 bind によって再表示し、listener を重複させない。hint/focus/hover だけでは canonical state、History transaction、Worker 操作を増やさない。実 drag は既存1 History操作。

SRP/DRY は既存 target 判定と gesture owner の再利用、OCP/LSP は既存 bind/serialization 契約の維持で守る。KISS/YAGNI に従い明示 mode と派生 hint に限定し、新しい state store、gesture default、保存 schema、dependency を追加しない。

## 検証

総合計画 A01/A02 と C08 の hint/export 部分を担当する。実ブラウザーで OFF pan、help/hover、ON drag/active cursor、keyboard toggle、touch 説明、Shift/Ctrl、feature/label/legend個別編集、record/alignment を確認する。Undo/Redo・Session load・Generate で hint が正しく再bindし、二重 listener・History entry がないことを検査する。

対象に応じて既存検査を拡張する:

```bash
node --test tests/web/composition-layout.test.mjs \
  tests/web/legend-layout-actions.test.mjs
```

serialization/export の既存検査を `rg` で特定して実行する。clean Resultと保存Session、plain/interactive SVG に transient 内容がないことを検査し、PNG/PDF は実 export を使い対応する shared clean-source 経路と artifact を確認する。hover screenshot だけで保存互換性を合格にしない。

Node/Python Playwright を確認して既存 composition/preview navigation browser tests を拡張する。必要なら source一致 wheel を準備する。source/session/生物学的値を変えず、公開手順・screenshot は S04 の既存 docs owner へ引き継ぐ。

## 完了・commit/push

`results/S02.md` に入力/authority/検査 SHA、受入対応、browser操作、serialization/export証拠、コマンド・結果、未達項目を保存し、進捗を更新する。production/tests/docs/generated の担当 diff を別々に review。owner/path 維持と追加 transient の除去を確認する。

branch/upstream を再確認し fetch。予期しない remote 更新は差分と担当を確認してから統合する。他者変更を revert せず、対象パスのみ stage、`git diff --cached --check` を実行する。**終了時は必ず commit/push。**

```bash
git commit -m "fix(web): explain how to enable preview layout editing"
git push origin HEAD:refs/heads/fix/issue-599-preview-layout-20260926
git rev-parse HEAD
git ls-remote --heads origin refs/heads/fix/issue-599-preview-layout-20260926
git status --short
```

local/remote SHA が一致することを確認する。push エラー時は実 remote state を確認して再試行。main/dev直接 push・force push・無承認のPR作成/merge/deployは行わない。公開 SHA、検証、未達条件、S03開始条件、英語 commit title と short summary を報告する。
