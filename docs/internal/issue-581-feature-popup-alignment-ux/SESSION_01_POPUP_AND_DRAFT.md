# Issue #581 — S01 Popup、候補表示、ローカル選択 INSTRUCTION PROMPT

以下を一つの独立した作業指示として使用する。過去の会話は前提にしない。

## 目的と開始条件

gbdraw Issue #581 の第一段階として、feature popup の Record actions の高さと Cancel の
問題を直し、Similarity Alignment の候補を生物学的に読めるようにし、候補 click のたびに
Python Resolver が走る待機をなくす。詳細な製品結果と受入条件は
`docs/internal/issue-581-feature-popup-alignment-ux/IMPLEMENTATION_PLAN.md` にある。
Product Decision Owner は 2026-09-24 に `A / EDIT_DISCLOSURE` と
`A / LOCAL_BATCH_RETRY` を承認した。PD-OI-033、034 は PR #582 で
`origin/dev` @ `e52e3ea9` にマージされ、固定実装ブランチへ取り込まれた。
開始時にこの祖先関係と作業ツリーを再確認する。

## ブランチと資料

- Issue #581 runtime は **`issue-581-feature-popup-alignment-ux-20260924`** にだけ実装する。
  別の runtime ブランチを作らない。開始時に branch、HEAD、upstream、worktree、
  最新 `origin/dev` との差と総合計画第10節を確認する。
- `AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`、
  `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` の PD-OI-026〜029、031〜034、
  `docs/internal/PRODUCT_IMPACT_RATCHET.md`、
  `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md` を読む。
- `gbdraw/web/js/app/similarity-alignment.js`、
  `gbdraw/web/js/app/orthogroups.js`、`gbdraw/web/js/app/feature-utils.js`、
  `gbdraw/web/js/app/app-setup.js`、
  `gbdraw/web/js/app/record-display/feature-record-rotation.js`、
  `gbdraw/web/index.html`、既存の関連単体・ブラウザテストを読む。

## 実装

1. popup の rotation domain binding と表示を分ける。承認済みの placement に従い、
   初期状態でフォームを閉じ、利用者が開ける入口を主操作の近くに置く。
   `recordActionsExpanded` は transient な一つの UI state とし、feature 切替、
   popup close、欄内 Cancel で reset する。欄内 Cancel は rotation draft だけを
   cancel し、popup は残す。rich／simple のフォーム処理を複製しない。
2. `candidateView` が enriched member と feature catalog の既存 metadata を照合し、
   gene、locus tag、product、type、source 座標、display strand を見せる。
   record heading は現在の表示 record に対応する definition／accession を優先する。
   ID の一致は exact identity で行う。内部 hash を見出しに使わず、必要な debug facts は
   disclosure に置く。代表 status や edge 数による順位は作らない。
3. 初回 `resolveRequest()` は維持する。曖昧 record と候補はその応答で固定し、
   recordKey ごとの選択を一つのローカル draft owner で持つ。radio の Select／Skip は
   同期的に draft を更新し、Worker を呼ばない。全曖昧 record の choice がある時だけ
   Apply を可能にする。「plan ready」という表現は Python 検証前に使わない。
4. Apply は全 choice を一回の request にまとめ、既存 `validateResolution()` で
   Python 応答を検証する。plan が確定した場合だけ既存 `applyPlan()` と生成経路を使う。
   validation／generation 失敗では choices と retry に必要な baseline を保持する。
   Cancel、成功、外部 mutation、古い応答の扱いを区別する。二重 Apply を防ぐ。
5. 曖昧性なしの自動適用、exact reference、Select／Skip の意味、Result／History／
   Session の契約を維持する。Python Resolver、Worker protocol、保存形式は変えない。
   旧 per-click Resolve 経路は残さない。

## 検証と終了

- `tests/web/similarity-alignment-actions.test.mjs` で、複数 record、選択変更、Skip、
  incomplete Apply、初回 1 回・各選択 0 回・Apply 1 回の helper 呼び出し、失敗後の
  retry、Cancel、stale 応答を確認する。実装のコピーになる test は作らない。
- `tests/web/feature-popup-record-rotation.test.mjs` と関連ブラウザテストで、
  feature open 時の閉鎖状態、開閉、欄内 Cancel 後の popup／Result 保全を確認する。
- focused tests、必要な実ブラウザ確認、Ruff 等を行う。Node Playwright がなければ
  Python Playwright を確認し、ブラウザ確認を省略しない。
- production diff と test diff を別々に監査し、総合計画第10節へ HEAD、変更、
  test、未決事項、S02 開始可否を記す。

SOLID のため owner を増やさず、KISS のため existing workflow を使い、DRY のため
選択と Apply の一経路に収束させ、YAGNI のため ranking、debounce、保存 state、
第二の Resolver を作らない。英語の proposed commit title と短い summary を示す。

## コミット・push と次セッションへの引き継ぎ

本依頼は担当分のコミットと同名 remote work branch への push を明示的に許可する。
このセッションの作業は固定ブランチ `issue-581-feature-popup-alignment-ux-20260924` 内で完了する。
検証と差分監査の後、担当分を英語の題名で１コミットにまとめ、同名の remote work
branch へ `git push origin HEAD:refs/heads/issue-581-feature-popup-alignment-ux-20260924` で push する。
開始前と push 前に branch、upstream、作業ツリー、remote の状態を確認する。

回答の最後に、次の S02 `SESSION_02_FLOATING_PALETTE.md` を新規参加者が単独で実行できる
完全な INSTRUCTION PROMPT として提示する。そのプロンプトにも、同じ固定ブランチで
作業・検証・コミット・push まで行い、さらに次セッション用の完全なプロンプトを
回答末尾に提示する指示を含める。S04 が完了した場合は追加セッションを作らない。
