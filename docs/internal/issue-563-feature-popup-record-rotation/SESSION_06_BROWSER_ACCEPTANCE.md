# Session 06 Instruction Prompt — End-to-end acceptanceとregression closure

## Prompt

あなたはgbdraw Issue #563のSession 06を担当する。過去の会話は前提にしない。

### Branchとprecondition

必ず`issue-563-feature-popup-record-rotation-20260922`を使用し、別branchを作らない。
Session 01のpreflightがruntime実装を許可し、Session 02–05のfocused testsがpassしていることを
確認する。branch、HEAD、upstream、statusを記録し、無関係な差分を保持する。

### 最初に全文を読む

1. `AGENTS.md`
2. `CLAUDE.md`
3. `gbdraw/web/CLAUDE.md`
4. `docs/internal/issue-563-feature-popup-record-rotation/IMPLEMENTATION_PLAN.md`
5. `BASELINE_AND_PREFLIGHT.md`
6. Session 02–05のproduction/tests
7. Issue #563最新本文の10 acceptance scenarios

### Session goal

実DOM、実file binding、canonical request、Worker/result admission、History、Session save/fresh Loadを通る
end-to-end testsで全acceptanceを閉じる。失敗した契約はownerで修正し、test-only bypassやsilent fallbackを
追加しない。

### Acceptance journeys

次を最低限coverする。既存fixtureを優先し、最小の新fixtureだけを追加する。

1. Popup-only:
   feature search -> feature popup -> record action -> Apply -> fresh Result。sidebar selection不要。
2. Both modes:
   同じcircular source recordとanchor intentがCircularの12時、Linearの左端に同じsource coordinateを置く。
3. Independent chromosomes:
   same GenBank fileのchromosome Iを`dnaA`、IIを`parB`で別々に回転し、他方を変えない。
4. Strand-aware offsets:
   plus/minusの5′ -100 bpがそれぞれupstreamへ解決する。
5. Stable orientation:
   initial forward/reverse complement、orient-forward OFF/ON、同じaction再適用がtoggleしない。
6. Endpoint distinction:
   3′ base anchorとfeature-endがforward/reverse displayで異なる。
7. Circular/compound:
   origin-spanning、multipart gap、origin crossing offset、odd/even midpoint。
8. Safe eligibility:
   unstranded、mixed/fuzzy/unordered、linear topology、cropped、stale popup、source replacement。
9. Identity/isolation:
   duplicate record IDs across files、split SVG fragments、same row multi-record。
10. State/rendering:
   Undo/Redo、v44 save/fresh Load、sidebar sync、render failure rollback、unrelated pending edits、
   feature/label/tick/depth/statistics/comparison geometry、LOSAT executor 0 additional jobs。

### Evidence requirements

- UI textだけでなく、last canonical requestのtarget absolute valuesを検証する。
- target外record request subtree、order、gridRow、tracks、comparison endpoint/resource IDの不変を検証する。
- generated SVGはtarget geometryが変わり、他record identityとbiological feature countが保たれることを検証する。
- split fragmentsが同じ`biologicalFeatureId`を維持することを検証する。
- LOSATはcache entry数ではなくinstrumented executor invocation/job countで0追加を証明する。
- failed/canceled/stale run後にold Result content、target intent、history depthが保持されることを検証する。
- fresh Loadはpage reload後にSession fileを読み直し、in-memory stateの継続だけで済ませない。
- unrelated pending editはApply後もUIとstateに残り、committed candidateには含まれないことを検証する。

### Test organization

- pure calculationはSession 02 testsへ置く。browser specでformulaを再実装しない。
- request/session/history contractは既存focused unit testへ置く。
- browser specは利用者journeyとboundary integrationだけを所有する。
- enormous one-testを避けるが、save/loadやrollbackのatomic journeyを不自然に分割しない。
- test-only app APIでactionを直接呼ぶだけのtestをprimary acceptanceにしない。実popup controlを操作する。

### Focused verification

新規test titleに安定した`record rotation` markerを付け、少なくとも次を実行する。

```bash
python -m pytest tests/test_web_feature_metadata.py tests/test_web_feature_catalog.py -v
node --test tests/web/feature-anchor.test.mjs tests/web/record-display-options.test.mjs
node --test tests/web/session-request.test.mjs tests/web/history.test.mjs tests/web/run-analysis-simple-path.test.mjs
npx playwright test tests/web/linear-multi-record.playwright.spec.js --grep "record rotation|LOSATP source jobs" --workers=1 --retries=0
npx playwright test tests/web/interactive-svg-v3.playwright.spec.js --grep "record rotation" --workers=1 --retries=0
node tests/web/architecture-contracts.test.mjs
git diff --check
```

必要なbrowser wheelは`python tools/prepare_browser_wheel.py`でsourceから作る。wheelはcommitしない。
Node PlaywrightがなければPython Playwrightで同等境界を検証する。sandbox failureは権限付き再実行する。

### 修正方針

- test failureはfirst failing boundaryまで絞り、semantic ownerで直す。
- UI、controller、projectorに同じvalidationを複製しない。
- timeout延長、retry増加、assertion削除で隠さない。
- LOSAT cache missを強制reuseで隠さず、key/evidence compatibilityを診断する。
- geometry expectation更新は意図したtransformだけかvisual/source identityをreviewする。
- `tests/reference_outputs/`は通常read-only。意図したrenderer geometry変更がないため更新しない。

### Handoff

Issueの10 scenarioと`AC-01`–`AC-20`のevidence matrix、全command/result、環境制約、修正したroot
cause、残課題を報告する。commit/push/PRは呼出時のuser指示がある場合だけ行う。
