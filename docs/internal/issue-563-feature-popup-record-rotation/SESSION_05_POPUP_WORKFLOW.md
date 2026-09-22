# Session 05 Instruction Prompt — Feature popup Record actions workflow

## Prompt

あなたはgbdraw Issue #563のSession 05を担当する。過去の会話は前提にしない。

### Branchとprecondition

必ず`issue-563-feature-popup-record-rotation-20260922`を使用し、別branchを作らない。
Session 01のpreflightがruntime実装を許可し、Session 02–04のfocused testsがpassしていることを
確認する。branch、HEAD、upstream、statusを記録し、無関係な差分を保持する。

### 最初に全文を読む

1. `AGENTS.md`
2. `CLAUDE.md`
3. `gbdraw/web/CLAUDE.md`
4. `docs/internal/issue-563-feature-popup-record-rotation/IMPLEMENTATION_PLAN.md`
5. `BASELINE_AND_PREFLIGHT.md`
6. Session 02–04で変更されたproduction/tests
7. `docs/REFERENCE/web-app.md`
8. `docs/internal/PRODUCT_IMPACT_RATCHET.md`

### Session goal

open feature popupの生物学的featureをexplicit targetとして、`Record actions`のpreview、Apply、Cancel、
disabled reasons、成功後continuityを実装する。domain式、request mutation、Worker操作をtemplateへ書かない。

### UI requirements

feature popupのfeature-wide controlsと混ぜず、独立した`Record actions` sectionを設ける。
rich popup/simple popupの両方で利用可能にする。

表示/入力:

```text
Rotate record using this feature…
Record: <target record>
Feature: <target feature>
Anchor: 5′ end | Midpoint | 3′ end
Offset: <signed integer> bp
Orient this feature forward
Place this feature at the end
New display start: <1-based source coordinate>
Record orientation: <unchanged / forward / reverse-complemented>
Displayed feature strand: <before -> after when known>
Coordinates refer to the original record.
Apply and regenerate | Cancel
```

defaultsは5′、0 bp、orientation保持。popup targetは必ず`clickedFeature.feat`の
`(recordKey, biologicalFeatureId)`から解決し、global selectedFeaturesや以前のselectionを参照しない。

### Controller behavior

1. popup open/target change時にephemeral action draftを初期化する。Record stateはまだ変更しない。
2. anchor/offset/orientation/preset変更ごとにSession 02 resolverでpreviewを導出する。
3. `Place this feature at the end`選択時はoffsetを0へ戻す。後でoffsetを変更したらcustom placementとして
   表示し、「exactly at end」を示し続けない。
4. unstranded featureではmidpointと安全なend placementだけをcapabilityに応じて許可し、
   upstream/downstream wordingを使わない。
5. mixed/fuzzy/unordered/linear topology/cropped/stale/replaced sourceはoperation別の具体的reasonを表示する。
6. Apply直前にもsource freshness、record length/topology、feature identityを再検証する。
7. ApplyはSession 04のtransactional actionを一回だけ呼ぶ。double click中は再送しない。
8. Cancel/closeはephemeral draftだけを捨て、record state、Result、Historyを変更しない。
9. success後にsidebarのeffective start/orientation表示を更新する。
10. feature search queryを保持し、stable identityでfeatureを再同定する。
11. 既存preview lifecycleが許す範囲でzoom/panを保持し、不必要なfull-view resetを行わない。

### Accessibility/responsive requirements

- popup内controlに一意なaccessible nameを付ける。
- disabled controlとApplyにはreasonを可視textで関連付ける。
- keyboardでselect/input/checkbox/buttonsへ到達し、Enter/Spaceで操作できる。
- drag handleとform interactionを`data-nodrag`で分離する。
- 390px幅でも横overflow、到達不能、切れたApply/Cancelを作らない。
- pending/rendering/failed stateをaria-liveまたは既存status surfaceで伝える。
- signed offsetのfractional/empty/out-of-range errorをfield近傍に出し、silent truncateしない。

### Minimal implementation shape

- popup-specific ephemeral state/controllerはfocused moduleに置き、`index.html`はbindingとmarkupに限定する。
- 既存`clickedFeature` objectへ長期persistent fieldsを混在させない。
- availabilityとApplyが同じresolver resultを使う。別々のeligibility式を作らない。
- sidebarとpopupは同じrecord display controls/effective stateを読む。
- popup close/openでdefaultを復元するが、committed `anchorIntent`が同じfresh targetにある場合は再表示に
  必要な値を復元してよい。provenanceをrender authorityとして使わない。

### Required tests

Component/controller:

- popup targetとunrelated global selectionの分離
- default values、preview update、offset validation
- feature-end preset/reset/custom transition
- operation別disabled reason
- stale source recheck before Apply
- double-submit prevention、Cancel no-op
- success後のsidebar syncとstable target rebind

Browser:

- search -> open popup -> Apply and regenerateのpointer journey
- keyboard operation
- rich/simple popup
- 390px viewport
- pending/failed status accessibility
- unrelated pending editがUI上も残る

### Verification

```bash
node --test tests/web/feature-anchor.test.mjs tests/web/record-display-options.test.mjs
node --test tests/web/*feature*popup*.test.mjs
npx playwright test tests/web/interactive-svg-v3.playwright.spec.js --grep "record rotation" --workers=1 --retries=0
node tests/web/architecture-contracts.test.mjs
git diff --check
```

実際の新規test file名が異なる場合は同等のfocused commandを記録する。Chromium sandbox failureは
repository guidanceどおり権限付きで再実行する。

### 禁止事項

- HTML template内のcoordinate formula
- global selection fallback
- SVG fragment IDだけによるtarget解決
- Apply前のpersistent state mutation
- generic modal framework、new dependency、build step
- popup操作からnormal Generate buttonをclickさせる実装

### Handoff

visible workflow、disabled reason一覧、accessibility/mobile evidence、controller boundary、tests、remaining browser
coverageを報告する。commit/push/PRは呼出時のuser指示がある場合だけ行う。
