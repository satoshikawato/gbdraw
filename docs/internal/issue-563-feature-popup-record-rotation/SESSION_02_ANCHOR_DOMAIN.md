# Session 02 Instruction Prompt — Source anchor factsとpure resolver

## Prompt

あなたはgbdraw Issue #563のSession 02を担当する。過去の会話は前提にしない。

### Branchとprecondition

必ず`issue-563-feature-popup-record-rotation-20260922`を使用し、別branchを作らない。
Session 01の`BASELINE_AND_PREFLIGHT.md`が存在し、runtime実装が`BLOCKED`でないことを確認する。
blockedならruntime fileを変更せず理由を報告する。

開始時にbranch、HEAD、upstream、statusを確認し、無関係な差分を保持する。`dev/main`へ直接commit、
rebase、reset、stash、pushをしない。

### 最初に全文を読む

1. `AGENTS.md`
2. `CLAUDE.md`
3. `gbdraw/web/CLAUDE.md`
4. `docs/internal/issue-563-feature-popup-record-rotation/IMPLEMENTATION_PLAN.md`
5. `docs/internal/issue-563-feature-popup-record-rotation/BASELINE_AND_PREFLIGHT.md`
6. `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`
7. `docs/internal/PRODUCT_IMPACT_RATCHET.md`

### Session goal

DOM、Vue、History、Workerから独立した一つのfeature-anchor domain resolverと、そのresolverが安全性を
判定するためのsource location factsを実装する。catalog wire/schemaへの接続、Session migration、
popup UI、render実行はまだ行わない。

### Owner

- PythonはBiopython location semanticsを失う前に、precision/order/strand capability factsだけを
  正規化する。anchor座標やUI wordingを決めない。
- JSの新規private domain module
  `gbdraw/web/js/app/record-display/feature-anchor.js`が5′/midpoint/3′、offset、orientation、
  feature-endを解決する唯一のownerになる。
- `record-display-options.js`の既存`selectedFeatureDisplayStart()`は最終的にこのresolverのadapterにする。
  同じ計算式を残さない。

### 実装要求

1. Python側に最小のsource anchor profile helperを追加する。
   - exact/fuzzy
   - single/join/order/unknown
   - biological/source-forward/ambiguous part order
   - consistent `+`/`-`、unstranded、mixed
   - source transformやreverse display後もoriginal source semanticsを表す
2. このsessionではhelperをcatalog schemaへまだ公開しない。private helperとしてunit testし、Session 03で
   catalog 4へ一度だけ接続できる形にする。
3. JS resolverはplain dataを受け、inputを変更せず、deterministicなresultを返す。
4. intent:
   - `placement: anchor | feature-end`
   - `anchor: five-prime | midpoint | three-prime`
   - signed safe-integer `offsetBp`
   - `orientForward`
5. output:
   - eligibility/capabilitiesとoperation別reason code/message material
   - 1-based `startCoordinate`
   - absolute `reverseComplement`
   - source anchor coordinate/outgoing boundary
   - before/after displayed strand
   - normalized non-authoritative provenance
6. stranded offsetは`1 + mod(a - 1 + s * o, L)`、unstrandedはsource正方向で解決する。
7. midpointはcovered partsだけをtraversal順に数え、`floor((N - 1) / 2)`を使う。
8. orientationはabsoluteかつidempotentとする。known `+`ならRC false、known `-`ならRC true。
9. feature-endは最終display directionの最後のincluded base直後であり、3′ anchorと区別する。
10. unsafe locationはoperation単位でdisableし、guessしない。

### Required tests

Python unit tests:

- plus/minus simple location
- same-strand compound joinとorigin-spanning order
- exact unstranded single/ordered path
- mixed strand、fuzzy、`order`/unknownのconservative classification
- source coordinate mapping後もfactsがoriginal sourceを表す

JS unit tests:

- 5′/midpoint/3′ for plus/minus
- positive/negative offsetsと両方向のwrap
- multipart gap exclusion、odd/even midpoint
- preserve orientation、orient-forward、repeat idempotence
- 3′とfeature-endの1-base差をforward/reverseで検証
- unstranded offset wording capability、ambiguous disable reason
- invalid/fractional offset、zero/unknown length、non-circular/cropped
- input deep-freeze/non-mutation

test fileは既存命名規則に従い、JSは原則`tests/web/feature-anchor.test.mjs`、Pythonは既存
feature metadata/catalog testへfocused casesを追加する。

### Verification

```bash
python -m pytest tests/test_web_feature_metadata.py tests/test_web_feature_catalog.py -v
node --test tests/web/feature-anchor.test.mjs tests/web/record-display-options.test.mjs
ruff check gbdraw/web_support
git diff --check
```

### Design constraints

- SVG coordinates、DOM、Vue ref、global selectedFeaturesをresolverへ渡さない。
- resolver内でresource lookup、Session migration、request clone、renderをしない。
- generic coordinate frameworkやclass hierarchyを作らない。
- existing backend `RecordDisplayTransform`を変更しない。
- profileとlocation partsの重複表現を最小化する。
- 既存5′/midpoint計算を二重に残さない。

### Handoff

変更したowner、pure interface、edge-case table、test command/result、Session 03が接続すべきprofile shapeを
報告する。commit/push/PRは呼出時のuser指示がある場合だけ行う。
