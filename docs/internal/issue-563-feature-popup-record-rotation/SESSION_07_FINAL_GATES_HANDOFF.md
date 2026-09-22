# Session 07 Instruction Prompt — Documentation、full gates、final handoff

## Prompt

あなたはgbdraw Issue #563のSession 07を担当する。過去の会話は前提にしない。

### Branchとprecondition

必ず`issue-563-feature-popup-record-rotation-20260922`を使用し、別branchを作らない。
Session 01のpreflightがruntime実装を許可し、Session 02–06のfocused/acceptance testsがpassしていることを
確認する。branch、HEAD、upstream、statusを記録し、無関係な差分を保持する。

### 最初に全文を読む

1. `AGENTS.md`
2. `CLAUDE.md`
3. `gbdraw/web/CLAUDE.md`
4. `docs/internal/issue-563-feature-popup-record-rotation/IMPLEMENTATION_PLAN.md`
5. `BASELINE_AND_PREFLIGHT.md`
6. Session 02–06の全production/test/docs diff
7. `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`
8. `docs/internal/PRODUCT_IMPACT_RATCHET.md`
9. `docs/REFERENCE/web-app.md`

### Session goal

新規参加者と利用者向けdocumentationを完成させ、production、tests、docs/generated artifactsを別々に
reviewし、required gatesを実行してIssue #563をhandoff可能にする。新機能を拡張しない。

### Documentation work

1. `docs/REFERENCE/web-app.md`または現行canonical Web manualへ次を追加する。
   - popupへの到達方法
   - 5′/midpoint/3′、signed offset、feature-endの意味
   - source coordinateとdisplay orientationの区別
   - effective circular topology、crop/ambiguous/stale時の制約
   - Apply/Cancel、Undo/Redo、Session persistence
   - unrelated pending settingsはApplyされないこと
2. developer-facing docsへ次を記録する。
   - owner matrixとcanonical path
   - feature catalog 4 / Session 44 compatibility
   - request schema 7不変
   - LOSAT executorを通常Generateも既にcache reuseしている事実
   - target-only candidateの目的はpending-edit isolationとatomic historyであること
3. `BASELINE_AND_PREFLIGHT.md`のstatus/evidenceをactual implementationへ更新する。Product authority wordingを
   実装者の推測で変更しない。
4. public screenshotが本当に必要な既存documentation workflowで要求される場合だけ、repositoryの
   screenshot guidance/skillに従う。owner-maintained social previewを変更しない。

### Separate reviews

Production review:

- target identity/freshnessをApply直前に再検証
- source/display、0/1-based、base/boundary naming
- one resolver、one request owner、one execution/admission path、one history owner
- target-only mutationと他record/comparison/resource不変
- failure/cancel/stale/superseded rollback
- manual edit provenance clear
- no new dependency/build step/CDN/CSP delta unless explicitly justified

Test review:

- Issue 10 scenariosとAC-01–AC-20にevidenceがある
- primary journeyは実popup DOMを操作
- default mode/settingsをtestが迂回していない
- LOSAT executor invocationを測定
- fresh page Session Load、Undo/Redo、failure、unrelated pendingを含む
- flaky wait、過剰timeout、weakened assertionがない

Docs/generated review:

- UI labelとdocsが一致
- Session/catalog versionとmigration説明が一致
- generated wheel、`dist/`、`gbdraw.egg-info/`、reference outputsがdiffに入っていない
- plan/result文書に過去会話依存の表現がない

### Full gates

環境で利用可能なrequired gatesをすべて実行する。長いtestはrepository guidanceどおり30分以上を許容し、
60秒以内に進捗を共有する。失敗は診断して修正後にaffected gateを再実行する。

```bash
python -m pytest tests/test_web_feature_metadata.py tests/test_web_feature_catalog.py tests/test_session_io.py -v
node --test tests/web/feature-anchor.test.mjs tests/web/record-display-options.test.mjs
node --test tests/web/session-request.test.mjs tests/web/history.test.mjs tests/web/run-analysis-simple-path.test.mjs
npx playwright test tests/web/linear-multi-record.playwright.spec.js --grep "record rotation|LOSATP source jobs" --workers=1 --retries=0
npx playwright test tests/web/interactive-svg-v3.playwright.spec.js --grep "record rotation" --workers=1 --retries=0
node --test tests/web/*.test.mjs
python -m pytest tests/ -v -m "not slow"
ruff check gbdraw/
node tests/web/architecture-contracts.test.mjs
node --test tests/ci/*.test.mjs
node tools/check-web-change-budget.mjs --base origin/dev --head HEAD
git diff --check
```

command名がrepository baseで変更されていたら現行equivalentを使い、置換理由を記録する。
browser sandbox failureは権限付きで同じcheckを再実行する。test reference outputは更新しない。

### Architecture/Product closure

- architecture ledgerをactual before/after OE、PE、CB evidenceで閉じる。
- owner/pathが非増加であること、または承認済み例外があることを確認する。
- Product Impact classificationとactual user effects/checkpointsが一致することを確認する。
- unresolved Product decisionが残ればruntime完了を宣言せず、そのaffected convergenceをblockする。
- Issue本文、docs、tests、runtimeでproduct outcomeが一致することを確認する。

### Final acceptance report

`docs/internal/issue-563-feature-popup-record-rotation/FINAL_ACCEPTANCE.md`を作り、次を含める。

- branchとHEAD/base
- user-visible outcome
- owner/path/compatibility summary
- Issue scenariosとAC-01–AC-20のpass evidence
- command、result、duration/環境上の未実行理由
- schema/version changes
- LOSAT additional executor job count
- known residual riskまたは`none`
- generated/untracked artifact disposition

### Final handoff

回答はoutcomeから始め、変更file群、重要な設計判断、verification、remaining riskを簡潔に示す。
repository guidanceに従いEnglishのproposed commit titleとshort summaryを付ける。

このprompt自体はcommit、push、PR作成、mergeを許可しない。呼出時にuserが明示的に許可した範囲だけ
実行する。push前はbranch名とupstreamを確認し、同名remote work branch以外へpushしない。
