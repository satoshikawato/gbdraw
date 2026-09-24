# Issue #581 — S04 最終受入と文書更新 INSTRUCTION PROMPT

以下を一つの独立した作業指示として使用する。過去の会話は前提にしない。

## 目的と開始条件

gbdraw Issue #581 の feature popup 改善と Similarity Alignment palette を一つの
利用者ワークフローとして検証し、必要な既存文書を更新する。Issue 本文は
https://github.com/satoshikawato/gbdraw/issues/581、設計提案は
https://github.com/satoshikawato/gbdraw/issues/581#issuecomment-5811037046、
全受入条件は `docs/internal/issue-581-feature-popup-alignment-ux/IMPLEMENTATION_PLAN.md`
にある。2026-09-24 承認の PD-OI-033〜035 が `origin/dev` にマージされたことを確認し、
S00〜S03 の実施記録、Product authority、未完了事項を読んでから始める。
未決 Product outcome があれば、その依存部分を完了と報告しない。

## ブランチと必読資料

- 最終修正・test・docs は **`issue-581-feature-popup-alignment-ux-20260924`** に置く。
  branch、HEAD、upstream、作業ツリー、`origin/dev` との関係を確認する。
- `AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`、
  `docs/internal/PRODUCT_IMPACT_RATCHET.md`、
  `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`、
  関連する accepted Product authority、総合計画、S01〜S03 diff を読む。
- Browser 手順は `playwright.config.js`、`package.json`、既存
  `tests/web/*.playwright.spec.js` を確認する。公開 Gallery の文書や screenshot を
  編集する場合だけ `.agents/skills/web-gallery-screenshot-maintenance/SKILL.md` を読む。

## 統合受入

1. Record actions: rich／simple popup で feature を開いた初期状態、明示展開、
   別 feature への切替、フォーム内 Cancel、Apply、既存 popup 操作を確認する。
2. Alignment: exact reference、曖昧性なしの自動経路、複数 record の候補・Skip、
   生物学的 label の欠損、Apply 一回の batch 検証、生成成功の plan／summary／History、
   validation／generation error 後の retry、Cancel、stale／superseded 応答を確認する。
3. Canvas/palette: desktop と 390 px で図の比較文脈、drag、focus、keyboard、
   pan／zoom／resize、番号と hover/click 同期、隠れた候補、overlay cleanup、
   Result・download・Session への非混入を確認する。
4. Fresh Generate、Session save/load、Undo/Redo、Reset と既存 active plan の意味を
   回帰確認する。backend の候補決定規則を変える必要が生じたら S00 の authority と
   計画を再評価し、UI の都合で Python Resolver を変えない。

## 実行と文書

- focused JS tests: `node --test tests/web/similarity-alignment-actions.test.mjs tests/web/feature-popup-record-rotation.test.mjs`。
- focused Python tests: `pytest tests/test_similarity_alignment.py tests/test_similarity_alignment_web_adapter.py -v`。
- browser: `node -e "console.log(require.resolve('@playwright/test'))"` を確認し、
  利用可能な既存 Playwright spec を実行する。Node が使えなければ Python Playwright で
  対応する実ブラウザ操作を行う。Chromium sandbox 問題は必要な escalation で再実行する。
- 変更内容に応じて `pytest tests/ -v -m "not slow"`、`ruff check gbdraw/`、
  Web policy／architecture gate を実行する。pytest は30分以上を見込み、途中経過を
  観察する。図の geometry を意図的に変更していないなら reference SVG を更新しない。
- 利用者に説明が必要な変更は既存 `docs/REFERENCE/web-app.md` 等の適切な owner へ
  追記する。新規 public page、Gallery screenshot、宣伝用 figure は必要性と再現可能な
  完成図を確認した場合だけ作る。

## 最終監査と報告

本依頼は担当分のコミットと同名 remote work branch への push を明示的に許可する。
production、test、docs、generated diff を別々に一度ずつ監査する。Product Impact の
各効果と AND-of-OR requirement、architecture owner/path の変化を確認し、superseded
path を残さない。総合計画第10節を更新し、実際の command、結果、残る制限、
最終 HEAD を記録する。固定ブランチ `issue-581-feature-popup-alignment-ux-20260924` 内で担当分を英語の題名で
１コミットにまとめ、同名の remote work branch へ
`git push origin HEAD:refs/heads/issue-581-feature-popup-alignment-ux-20260924` で push する。
開始前と push 前に branch、upstream、作業ツリー、remote の状態を確認する。
PR と merge はその時点の明示許可に従う。

一つの Resolver、一つの選択 owner、一つの Result 経路を保つことを
SOLID／KISS／DRY の受入条件とする。未使用の extension point、互換 branch、
保存形式を残さないことを YAGNI の受入条件とする。
