# Session 04 Instruction Prompt — Target-only candidateとatomic render transaction

## Prompt

あなたはgbdraw Issue #563のSession 04を担当する。過去の会話は前提にしない。

### Branchとprecondition

必ず`issue-563-feature-popup-record-rotation-20260922`を使用し、別branchを作らない。
Session 01のpreflightがruntime実装を許可し、Session 02/03のfocused testsがpassしていることを
確認する。branch、HEAD、upstream、statusを記録し、無関係な変更を保持する。

### 最初に全文を読む

1. `AGENTS.md`
2. `CLAUDE.md`
3. `gbdraw/web/CLAUDE.md`
4. `docs/internal/issue-563-feature-popup-record-rotation/IMPLEMENTATION_PLAN.md`
5. `BASELINE_AND_PREFLIGHT.md`
6. Session 02/03で変更されたproduction/tests
7. `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`
8. `docs/internal/PRODUCT_IMPACT_RATCHET.md`

### Session goal

最後に成功したcanonical requestをbaseにtarget recordだけを変更したcandidateを、通常Generateと同じ
Worker/result admission pathで実行する。target display intentとgenerated artifactを一つのUndo/Redo
entryとしてcommitし、failure/cancel/stale/superseded時は両方をrollbackする。UIはまだ追加しない。

### 重要な前提

通常Generateでもdisplay start/reverse complementだけの変更後は、互換なLOSAT evidenceをcacheから
再利用し、executorを再実行しない。専用candidateの目的はLOSAT回避そのものではなく、unrelated
pending form changesの隔離とatomic historyである。cache lookup/planningを迂回するための別経路を
作らない。

### Owner boundaries

- `session-request.js`: committed request clone、target identity解決、必要なexact-one materialization、
  two-field overlay、canonical validationの唯一owner。
- `run-analysis.js`: canonical candidate execution、comparison cache、Worker、typed response、result admissionのowner。
- `services/history.js`: artifactとintentのatomic replacement/rollback owner。
- feature action controller: 上記ownerを順に呼ぶだけ。JSON編集、Worker生成、SVG挿入をしない。

### 実装作業

1. `session-request.js`へpure target mutation/projectorを追加する。
   inputはlast successful canonical request、target source-bound identity、resolved absolute
   `startCoordinate/reverseComplement`。outputはvalidated canonical candidateとtarget resolution receipt。
2. projectorはtarget recordがちょうど一つであることを要求する。0件、複数件、stale resource、cropped、
   non-circular、length mismatchを具体的errorにする。
3. 変更可能fieldはtargetの`display.startCoordinate`と
   `presentation.reverseComplement`だけ。必要なexact-one materializationはSession 03 ownerを使う。
4. request order、他records、diagram options、tracks、comparisons、resource IDsを保持し、inputを変更しない。
5. `run-analysis.js`から、既にcanonicalなrequestを受け取るcandidate execution/admission helperを抽出する。
   normal Generateもpopup actionも同helperを使う。既存cancel/stale token、loading/error state、resource staging、
   comparison cache、Worker lifecycle、result validation、SVG sanitization/admissionを維持する。
6. `runAnalysis()`の現行挙動を変更しない。normal pathはcurrent formからcanonical requestを構築後、同helperへ渡す。
7. `runUndoableArtifactReplacement()`へ任意のintent checkpoint hooksを追加する。
   - before/after target draftだけをcaptureする。
   - Undo/Redoはartifact handleとintentを同じentryでrestoreする。
   - hookなしの既存callerは完全に同じ挙動。
   - byte/file retention、checkpoint、rollback diagnosticsを壊さない。
8. popupが後で使うfocused controller actionを追加する。explicit clicked feature identityとnormalized intentを受け、
   resolver -> projector -> execution -> commitを調停する。Vue templateへdomain logicを漏らさない。
9. candidate success時だけrecord draftのabsolute stateとanchorIntentをcommitする。candidate execution前の一時patchを
   global stateへ残さない。

### Required failure semantics

- resolver/projector validation失敗: Workerを起動せず、intent/Result/History no-op
- user cancel: no-op
- Worker/render/admission失敗: before intentとbefore Resultを維持
- stale/superseded completion: current stateを置換しない
- Undo: before origin/orientation/provenance/Resultへ一回で戻る
- Redo: afterへ一回で戻る
- unrelated pending form edits: candidateへ入らず、成功後もworking draftに残る

### Required tests

Pure/projector:

- one target two-field overlay
- target以外のcanonical subtree deep equality
- all-to-exact materialization reuse
- stale/duplicate/missing/non-circular/cropped rejection
- input non-mutation

Execution/history:

- normal Generateとtarget candidateが同admission helperを呼ぶ
- successでhistory entry一件
- Undo/Redoでintent+artifact一件
- failure/cancel/stale/superseded rollback
- existing artifact-only callerの回帰なし
- unrelated pending form fieldがcommitted requestに入らず、working stateに残る
- compatible LOSAT cacheでexecutor invocation 0 additional jobs

### Verification

```bash
node --test tests/web/session-request.test.mjs
node --test tests/web/history.test.mjs tests/web/run-analysis-simple-path.test.mjs
node --test tests/web/run-analysis-derived-cache.test.mjs tests/web/losat-cache.test.mjs
node tests/web/architecture-contracts.test.mjs
git diff --check
```

必要なら既存browser testのfocused contractを追加実行するが、popup journeyはSession 05/06で行う。

### 禁止事項

- `runAnalysis()`をそのままpopup actionから呼び、current form全体をcommitすること
- LOSATを強制skipするflag、cache miss時のsilent stale reuse
- 新しいWorker、SVG admission、sanitizer、result store、history stack
- intentとartifactを別々のhistory entryにすること
- failureを成功扱いにするfallback

### Handoff

共有execution path、projector invariant、history before/after、failure matrix、LOSAT invocation evidence、test結果を
報告する。commit/push/PRは呼出時のuser指示がある場合だけ行う。
