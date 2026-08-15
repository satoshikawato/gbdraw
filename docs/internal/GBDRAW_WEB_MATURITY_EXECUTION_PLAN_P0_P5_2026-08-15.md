# gbdraw Web 成熟化実行計画 P0–P5

**文書種別:** 実装ロードマップ  
**対象:** `satoshikawato/gbdraw` Web application  
**基準日:** 2026-08-15  
**基準ブランチ:** `main`  
**基準SHA:** `3a098d8f25bb308833b7c0d82c6edd27e8e9a485`  
**基準イベント:** PR #336 `Preserve loaded session intent on regeneration` merge後  
**計画の目的:** この文書を各Codex実装プロンプトの上位仕様として使い、P0からP5までを一つずつ独立PRで進める。

---

## 1. 現在地

gbdraw Webは、機能面・回帰検出面ではすでに本格的なアプリケーションである。巨大session、persistent Worker、Pyodide、LOSAT、直接SVG編集、Save/Load、History、Undo/Redo、current/legacy session、実Vibrio fixtureまで含むcontractが存在する。

一方、現在の主要な弱点は、機能不足ではなく次の三点にある。

1. **一つの操作に関する責務が複数moduleへ分散している。**
2. **runtime testの強さに対して、JavaScriptの静的境界検証が弱い。**
3. **回帰を防ぐtest自体が巨大化し、次の変更の局所性を落とし始めている。**

したがってP0–P5では、追加機能を増やすよりも、以下の順序で構造を成熟させる。

```text
P0  基準点を固定し、既知のuser regressionを閉じる
P1  live SVG editのcommit/replay ownershipを一本化する
P2  current-session active config schemaをconfig.jsから分離する
P3  JavaScript境界moduleへ段階的な静的型検証を導入する
P4  browser contract harnessを分解し、test変更局所性を回復する
P5  complete generation identityにより完全同一GenerateをWorker前で終了する
```

---

## 2. この計画における「成熟」の定義

成熟度は機能数ではなく、以下の性質で評価する。

### 2.1 Authorityが一つである

各state domainについて、次を答えられること。

- 誰が書くのか。
- 誰がvalidateするのか。
- 誰がpersistするのか。
- 誰がrollbackするのか。
- 誰がfresh artifactへreplayするのか。

### 2.2 状態遷移が明示されている

少なくとも以下が同じowner modelで説明できること。

```text
success
no-op
failure
Cancel
stale result
superseded result
Undo
Redo
Save
Load
legacy migration
```

### 2.3 変更局所性が高い

一つの概念変更について、関係のないownerを多数編集しなくて済むこと。

目安:

- 一つのlive-edit domain追加で、production変更がowner module＋UI action module程度に収まる。
- session field追加で、`config.js`、request builder、state、testsの複数箇所に独立したfield listを増やさない。
- test helperを各specで複製しない。

### 2.4 契約が実行可能である

「そうなるはず」ではなく、structural metricまたはowner contractで検証する。

例:

```text
direct edit Worker run count = 0
mounted SVG Result serialization count = 1
current-session active config canonical overwrite count = 0
unchanged Generate Worker run count = 0
```

### 2.5 境界のshapeが静的に検証される

特に次の境界はJSDoc＋`checkJs`等で検証されること。

```text
session active config
preview runtime
editor SVG projection
diagram Worker protocol
GeneratedArtifactHandle
History APIs
generation identity
```

### 2.6 testが保守可能である

test coverageを維持しつつ、scenarioの意図が読めること。

- 一つのspecが全helper、全fixture、全comparisonを所有しない。
- 共通helperは一つのownerへ置く。
- scenario本体は入力、操作、期待結果を中心に記述する。

---

## 3. 全phase共通の実行規則

### 3.1 一phase一PR

P1からP5は原則として一つずつ独立PRにする。

```text
P1を実装中にP2へ着手しない
P2を実装中にP3を導入しない
P3でproduction behaviorを変更しない
P4でproduction refactorを混ぜない
P5でsession schema変更を混ぜない
```

一つのphaseが大きすぎる場合は、同じphase内で`P1a/P1b`のように分割できる。ただし、次phaseへ進んではならない。

### 3.2 PR開始条件

各PR開始時に必ず確認する。

```bash
git status --short --untracked-files=all
git rev-parse HEAD
git diff --stat
```

さらに以下を満たすこと。

- 直前phaseがmainへmerge済み。
- required CIがgreen。
- 直前phaseのremaining riskにcorrectness blockerがない。
- 計画書のbaseline SHAとstatusを更新した。

### 3.3 PR完了条件

各PRは以下を返す。

1. starting SHA;
2. root causeまたはstructural deficiency;
3. owner model;
4. files changed;
5. production/test line counts;
6. structural assertions;
7. focused tests;
8. broad tests;
9. browser contracts;
10. schema/artifact compatibility;
11. remaining risks;
12. proposed PR title/body.

### 3.4 禁止事項

P0–P5を通して次を避ける。

```text
repository-wide auditだけで終了する
measurement-only PR
test timeout引き上げで問題を隠す
session/reference SVGを安易に更新する
runtime network dependency追加
main-thread Pyodide復活
別render path追加
final SVG cacheを一般解として導入
Historyのfull checkpointへ戻す
untrusted input validationを省略
```

### 3.5 Evidenceの扱い

大きなevidence JSONやSVGは原則として`test-results/`へ出し、tracked artifactにしない。

tracked contractには以下を残す。

- structural count
- deterministic digest
- owner identity
- bounded broad timing ceiling
- exact semantic assertions

wall-clock値はlocal evidenceとし、tight CI thresholdにはしない。

---

## 4. Phase一覧

| Phase | 主題 | 依存 | 主な成果 | 状態 |
|---|---|---|---|---|
| P0 | Post-PR #336 baseline closure | なし | exact user session確認、基準点固定 | **実装済み／manual確認残り** |
| P1 | Live SVG edit commit/replay ownership | P0 | mounted edit ownerとfresh-SVG replay owner | 未着手 |
| P2 | Active config schema分離 | P1 | pure current-writer schema/validator/restore module | 未着手 |
| P3 | JS boundary static checking | P2 | JSDoc＋`tsc --checkJs`のCI gate | 未着手 |
| P4 | Browser contract harness分解 | P3 | 再利用可能helper、scenario spec小型化 | 未着手 |
| P5 | Pre-Worker unchanged-Generate no-op | P1–P4 | complete generation identity、即時no-op | 未着手 |

---

# P0. Post-PR #336 baseline closure

## 目的

PR #336で修正した以下のcontractを、実際にregressionを起こしたprivate sessionでも確認し、P1以降の基準点を固定する。

```text
Load
→ Generateなしで直接編集
→ DOMとselected Resultが同期
→ Save
→ fresh Load
→ editsが即時表示
→ Generate
→ plain diagramへ戻らない
→ Undo
→ Redo
```

## 現在の状態

PR #336はmainへmerge済みである。Gallery fixtureを用いたautomated contractsは追加済み。

残っているのは、private user sessionによるmanual acceptanceである。

## 実施項目

### P0-1. Exact user session acceptance

private fileをrepositoryへcommitせず、以下を実施する。

1. exact sessionをLoad。
2. Generate前にFeature fillを変更。
3. stroke、Feature visibility、label textまたはvisibilityを変更。
4. Worker/Python loadingが発生しないことを確認。
5. Save Session。
6. fresh pageでLoad。
7. editsが即時表示されることを確認。
8. Generate Diagramを1回実行。
9. editsが維持され、plain diagramへ戻らないことを確認。
10. Undo/Redo。
11. 再Save、再Load、再Generateを1回実行。

### P0-2. Baseline record

以下を記録する。

```text
main SHA
SESSION_VERSION
canonical request schema
feature catalog schema
Vibrio Generate 1 / Generate 2の最新既知値
direct-edit Worker/Python delta
session-load preflight
History baseline
prepared cache retained-byte estimate
```

高コストbenchmarkを再実行する必要はない。直近のverified evidenceをbaselineとして採用する。

### P0-3. Failure時の処理

exact fileで失敗した場合はP1へ進まない。

次を作る。

- private biological dataを除いたminimal sanitized fixture
- failing domainの一覧
- loaded artifact／active intent／direct edit／regenerationのどこで失われたか
- focused correctness prompt

## Exit criteria

```text
exact user session acceptance passed
main required CI green
known correctness blocker = 0
P1 baseline SHA recorded
```

## 成果物

- 本計画書のP0 status更新
- private acceptance結果の短い記録
- failureがあればsanitized reproducerとcorrectness prompt

## P0用Codex promptの扱い

exact fileが通る場合、P0用の実装promptは不要。

失敗した場合のみ、failure domainに限定したpromptを作る。

---

# P1. Centralize live SVG edit commit and replay

## 目的

一つのlive editを構成する以下の処理を、明示的な二ownerへ集約する。

```text
A. mounted SVGを編集しselected Resultへcommitするowner
B. fresh generated/reflowed SVGへcanonical editor stateをreplayするowner
```

## 背景

PR #336では、同じ概念に属する以下のregressionが見つかった。

```text
Feature visibilityはDOMだけ変え、Resultをflushしなかった
direct label visibilityが不要なWorkerを起動した
generated SVG admissionでlabel overridesが消えた
legend strokeはreplayされたがlegend colorは欠落した
```

これらは個別バグというより、live-edit protocolの責務分散である。

## Target owner model

### Mounted edit owner

候補:

```text
app/preview-runtime.js
```

または隣接する単一purpose module。

責務:

```text
active SVG取得
DOM mutation
no-op判定
index invalidation
dirty管理
SVG serialize
selected Result replacement
structured outcome返却
```

### Fresh SVG projection owner

候補:

```text
app/editor-svg-projection.js
```

責務:

```text
Feature fill
Feature stroke
Feature visibility
Label text
Label visibility
Legend fill
Legend stroke
specific-rule-derived fill
fresh SVGに必要なsuppression
```

action moduleはscope、validation、canonical override mutationのみを所有する。

## In scope

```text
single/bulk Feature fill
single/bulk Feature stroke
single/bulk Feature visibility
label text
label visibility
legend fill
legend stroke
legend renameのmounted SVG変更部分
normal Generate admission
label reflow admission
session preview restoration上のoverride replay
```

## Out of scope

```text
generation key
unchanged Generate short-circuit
session active-config schema分離
checkJs導入
test harness全面分解
History redesign
```

## Structural success metrics

```text
Feature/Label/Legend action moduleからserializeCleanSvg direct call = 0
Feature/Label/Legend action moduleからResult content direct replacement = 0
mounted edit commit production owner = 1
fresh SVG complete editor projection owner = 1
single edit Result serialization = 1
single edit Result replacement = 1
same value no-op serialization = 0
bulk edit serialization = 1
direct edit Worker construction = 0
direct edit Worker initialization = 0
direct edit Worker run = 0
Python render/helper = 0
```

## Required behavior matrix

各domainで以下を検証する。

| Domain | Canonical override | DOM即時変更 | Result同期 | Save/Load | Generate replay | Undo/Redo |
|---|---:|---:|---:|---:|---:|---:|
| Feature fill | 必須 | 必須 | 必須 | 必須 | 必須 | 必須 |
| Feature stroke | 必須 | 必須 | 必須 | 必須 | 必須 | 必須 |
| Feature visibility | 必須 | 必須 | 必須 | 必須 | 必須 | 必須 |
| Label text | 必須 | 必須 | 必須 | 必須 | 必須 | 必須 |
| Label visibility | 必須 | 必須 | 必須 | 必須 | 必須 | 必須 |
| Legend fill | 必須 | 必須 | 必須 | 必須 | 必須 | 必須 |
| Legend stroke | 必須 | 必須 | 必須 | 必須 | 必須 | 必須 |

## Label-specific requirements

```text
same stable targetがfresh SVGに存在する → overrideを保持
targetがbiological/filter changeで消失 → existing semanticsに従いunresolvedまたは除去
ambiguous target → fail closed
optional reflow failure → direct editを保持
```

## Verification gates

Focused:

```bash
node --test \
  tests/web/preview-runtime.test.mjs \
  tests/web/feature-color-actions.test.mjs \
  tests/web/feature-visibility-actions.test.mjs \
  tests/web/history.test.mjs \
  tests/web/architecture-contracts.test.mjs
```

Browser:

```bash
npx playwright test \
  tests/web/contracts/session-regenerate-intent.playwright.spec.js \
  --workers=1 \
  --retries=0

npx playwright test \
  tests/web/contracts/current-session-lazy-materialization.playwright.spec.js \
  --workers=1 \
  --retries=0
```

Broad:

```bash
node --test tests/web/*.test.mjs
npm run test:web:functional-smoke
ruff check gbdraw/
git diff --check
npm run test:web:vibrio-generate
```

## Exit criteria

```text
live edit commit owner = 1
fresh SVG projection owner = 1
all seven direct-edit domains pass the matrix
no hidden Generate
no Worker/Python on direct edit
Save/Load/Generate/Undo/Redo contracts preserved
no schema/artifact changes
```

## Expected PR title

```text
Centralize live SVG edit commit and replay
```

## Prompt source

既に作成したP1 promptを、この計画書と最新main SHAに合わせて更新して使用する。

---

# P2. Extract current-session active config schema

## 目的

`config.js`に集積したcurrent-writer active configのschema、validation、restore logicをpure ownerへ分離する。

候補:

```text
gbdraw/web/js/services/session-active-config.js
```

## 背景

PR #336で以下が明示化された。

```text
renderRequest / Results / catalog = committed artifact authority
config / UI controls / editor overrides = next Generate authority
```

しかし現在のactive config domain inventory、shape validation、compat alias、restore logicは`config.js`内にあり、save/load coordinationとschema ownershipが同居している。

## Target owner model

### `session-active-config.js`

所有するもの:

```text
CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS
domain kind
form field inventory
adv field inventory
collection entry schema
compat aliases
prototype-pollution rejection
supported value validation
restoreCurrentWriterActiveConfig()
```

所有しないもの:

```text
file input
session gzip
resource materialization
preview admission
Vue reactive state mutation
alert/download
legacy whole-session orchestration
```

### `config.js`

責務を次へ縮小する。

```text
session file read
version dispatch
migration coordination
canonical projection coordination
active config validator呼び出し
state application
preview commit
rollback
session export
```

## Purity requirement

active-config schema ownerは次を満たす。

```text
window依存なし
Vue依存なし
DOM依存なし
File/Blob依存なし
alert/downloadなし
可能ならstate.js direct importなし
```

field inventoryをreactive objectの現在のpropertiesから暗黙生成しない。current writerのschemaは明示的なcontractとして管理する。

## In scope

```text
domain inventory移動
shape validator移動
compat alias normalization
restore/fallback policy
buildConfigDataとのdrift test
current/legacy session tests
config.jsのactive-config code削減
```

## Out of scope

```text
SESSION_VERSION変更
canonical request schema変更
legacy session policy変更
P1 editor protocol変更
TypeScript/checkJs導入
test harness分解
generation key
```

## Structural success metrics

```text
config.js内のactive config domain list = 0
config.js内のform/adv field schema list = 0
active config pure owner = 1
active config module state.js import = 0（原則）
active config module DOM/window usage = 0
buildConfigData domain drift test = 1
```

## Required tests

### Current writer

```text
all supplied domains restore
omitted domains use documented canonical fallback
unknown domain rejects
unknown form/adv field rejects
unsafe prototype key rejects
invalid option value rejects
invalid collection shape rejects
invalid input fails before document replacement
```

### Compatibility

```text
colorsAreOverrides accepted only as documented alias
adv.losatProgram accepted only as documented alias
aliases stripped before runtime application
v39 migration unchanged
bare legacy config unchanged
```

### Drift protection

`buildConfigData()`が新domainを追加した場合、domain inventory testが失敗する。

schema module側のfield追加なしに新fieldがpersistされない。

## Verification gates

```bash
node --test \
  tests/web/session-draft-authority.test.mjs \
  tests/web/session-export-validation.test.mjs \
  tests/web/session-active-files.test.mjs \
  tests/web/session-definition-rehydration.test.mjs \
  tests/web/session-resource-backing.test.mjs

npx playwright test \
  tests/web/contracts/session-regenerate-intent.playwright.spec.js \
  --workers=1 \
  --retries=0

node --test tests/web/*.test.mjs
npm run test:web:functional-smoke
ruff check gbdraw/
git diff --check
```

Vibrio contractは最後に一度実行する。

## Exit criteria

```text
active config schema/validator/restoreのpure ownerが存在
config.jsはcoordination ownerへ縮小
current/legacy/divergent session contractsが維持
field/domain driftがtestで検出される
no schema version change
```

## Expected PR title

```text
Extract the current-session active config boundary
```

---

# P3. Add staged JavaScript boundary checking

## 目的

全体TypeScript化を行わず、壊れた場合の影響が大きいboundary moduleへJSDocと`tsc --checkJs`を段階導入する。

## 背景

runtime contractsは強いが、JavaScript object shape、optional method、read-only global、structured resultの誤用は実行時まで検出されにくい。

P1とP2で境界を明確化した後、その境界を静的に固定する。

## Initial checkJs scope

最低限:

```text
preview-runtime.js
editor-svg-projection.js
session-active-config.js
diagram-worker-protocol.js
diagram-generation.js
history.js
history-snapshot.jsのGeneratedArtifactHandle boundary
candidate-render.js
```

必要に応じて小さなshared typedef fileを使う。

## Tooling

候補:

```text
jsconfig.boundaries.json
typescript devDependency
npm script: test:web:types
tsc --noEmit --allowJs --checkJs
```

Web appは引き続きno-build SPAである。TypeScriptはCIの静的検証にのみ使用する。

## Required contracts

JSDocで最低限定義する。

```text
MountedDomEditOutcome
PreviewRuntime
EditorSvgProjectionInput
CurrentWriterActiveConfig
GeneratedArtifactHandle
DiagramGenerationRequest/Response
History artifact replacement entry
```

## Rules

```text
target filesに@ts-nocheckを追加しない
大量のanyで通さない
production runtime codeをtype checker都合だけで複製しない
generated declaration buildを導入しない
browser runtimeにTypeScriptを追加しない
```

外部libraryやVue reactive shapeで厳密化困難な部分は、狭いadapter boundaryで型を限定する。

## CI integration

既存Browserまたは独立した軽量jobで実行する。

```bash
npm run test:web:types
```

期待所要時間は短く保つ。

## In scope

```text
TypeScript devDependency
boundary jsconfig
JSDoc typedef
targeted type fixes
CI gate
```

## Out of scope

```text
.js→.ts rename
bundle/build step
全Web module checkJs
Vue template type checking
production architecture変更
session schema変更
```

## Structural success metrics

```text
boundary files checkJs coverage >= 指定一覧
target file @ts-nocheck = 0
type gate CI = green
runtime artifact change = 0
reference SVG change = 0
```

## Required tests

```bash
npm run test:web:types
node --test tests/web/*.test.mjs
npm run test:web:functional-smoke
ruff check gbdraw/
git diff --check
```

production behaviorが変わった場合はscope violationとして扱う。

## Exit criteria

```text
boundary type gateがCI required checkとして存在
P1/P2 owner APIsのshape driftが静的に検出される
no build step
no runtime dependency
no public Python/CLI change
```

## Expected PR title

```text
Type-check critical Web ownership boundaries
```

---

# P4. Decompose browser contract harness

## 目的

強いbrowser contractを維持しながら、巨大specとhelper duplicationを分解し、test側の変更局所性を回復する。

## 背景

`session-regenerate-intent.playwright.spec.js`は、session load/save、SVG comparison、raster evidence、direct edit、legacy config、History、Worker metricsを一つのfileで所有している。

coverageは強いが、新scenario追加のたびに巨大fileを編集する構造は成熟していない。

## Target structure

候補:

```text
tests/web/helpers/session-contract.cjs
tests/web/helpers/svg-evidence.cjs
tests/web/helpers/worker-evidence.cjs
tests/web/helpers/direct-edit-fixture.cjs

tests/web/contracts/session-active-intent.playwright.spec.js
tests/web/contracts/session-direct-edit-roundtrip.playwright.spec.js
tests/web/contracts/session-legacy-config.playwright.spec.js
```

必ずしもこの名前に固定しない。

## Helper ownership

### session helper

```text
load session
wait for interactiveReady
save session
fresh page load
current version assertions
```

### SVG evidence helper

```text
canonical SVG
semantic projection
raster digest
catalog digest
equivalence assertion
```

### Worker evidence helper

```text
construction
initialization
runs
resource transfer
Python/helper activity
```

### direct edit fixture

```text
deterministic Feature/Label/Legend target selection
actual UI/public action execution
canonical override evidence
```

## Constraints

```text
generic test frameworkを作らない
production helperをtestのために公開しない
assertionをhelper内に隠しすぎない
scenario固有の期待値はspecに残す
大きなevidence captureを複製しない
```

## Quantitative targets

目安:

```text
largest scenario spec <= 900–1200 lines
one helper file <= 700 lines
duplicate evidence-capture implementation = 0
functional CI duration regression <= 15%
test countとsemantic assertionsを減らさない
```

line数はabsolute gateではなく、読みやすさとowner分離の補助指標とする。

## In scope

```text
test-only helper extraction
scenario split
pure comparison helper unit tests
Playwright config/tag調整が必要なら最小限
CI command維持
```

## Out of scope

```text
production code変更
behavior change
performance optimization
schema change
new public fixture
reference SVG更新
P5 generation identity
```

## Required verification

```bash
node --test tests/web/*.test.mjs
npm run test:web:functional-smoke
npm run test:web:perf-smoke
npm run test:web:vibrio-generate
python tests/run_losat_cache_browser_acceptance.py
ruff check gbdraw/
git diff --check
```

各split specを単独でも実行する。

## Exit criteria

```text
production diff = 0
common session/SVG/Worker/direct-edit helper ownerが明確
scenario specsが入力・操作・期待結果中心になる
assertion coverage不変
CI時間が許容範囲
```

## Expected PR title

```text
Decompose Web session contract harnesses
```

---

# P5. Pre-Worker unchanged-Generate no-op

## 目的

完全に同じ最終artifactを既に表示している場合、Generateボタンを再度押してもWorker/Pythonを実行せず、即時に`unchanged`として終了する。

## 依存

P5はP1–P4の完了後にのみ着手する。

理由:

```text
P1: editor replay identityを一ownerへ集約する必要がある
P2: active configのcomplete schema inventoryが必要
P3: generation identity boundaryを静的に守る必要がある
P4: false-positive/false-negativeを検証するtest harnessが必要
```

## 設計原則

一つの巨大hashではなく、意味の異なる二componentを持つ。

```text
RenderCoreIdentity
EditorProjectionIdentity
```

### RenderCoreIdentity

Python/Workerが生成するbase diagramを決める。

含むもの:

```text
canonical render semantics
mode/topology
render/layout settings
biological input selection
resource identity: resourceId + cacheToken + size
LOSAT/orthogroup/annotation analysis identity
renderer/build compatibility version
```

含めないもの:

```text
full file bytes再hash
workspace path
filename alone
final SVG
feature catalog
transient UI
```

### EditorProjectionIdentity

P1のfresh-SVG projection ownerが適用するpost-render stateを決める。

含むもの:

```text
Feature fill
Feature stroke
Feature visibility
Label text
Label visibility
Legend fill
Legend stroke
specific rules
relevant suppression flags
```

### CommittedGenerationIdentity

```text
{
  renderCoreIdentity,
  editorProjectionIdentity,
  schema
}
```

成功したGenerateまたは安全にdirect-commitされたartifactに対応してruntimeで保持する。

## Direct editとの関係

P1でdirect editがmounted SVGとselected Resultへ完全commitされ、geometry reflow pendingでなければ、EditorProjectionIdentityを更新できる。

その場合:

```text
base render coreは同じ
editor projectionはcurrent artifactへ適用済み
pending required reflowなし
```

であればGenerateをno-opにできる。

ただし以下では必ずrenderする。

```text
geometry/layout setting変更
pending required label reflow
biological input変更
resource token変更
LOSAT/orthogroup/annotation identity変更
session load後でtrusted runtime identityがない
前回Generate failure/Cancel
identity construction failure
```

## Imported session policy

P5ではsession schemaを変更しない。

したがって:

```text
Load直後 → trusted runtime generation identityなし
最初のGenerate → full render
成功後 → identityをruntimeで保持
2回目の完全同一Generate → no-op可能
```

将来sessionへidentityをpersistするかはP5のscope外。

## Fail-open rule

identityが不完全、unknown、ambiguous、validation errorの場合は必ずrenderする。

```text
false negative: 許容（不要なrenderが起こる）
false positive: 不可（必要なrenderをskipする）
```

## Structural success metrics

完全同一warm Generate:

```text
generation identity comparisons = 1
Worker construction delta = 0
Worker initialization delta = 0
Worker run delta = 0
Python invocation delta = 0
resource materialization delta = 0
resource transfer delta = 0
Result admission delta = 0
History entry delta = 0
current committed owner unchanged
```

## Required scenarios

### A. 完全同一

```text
Generate A
→ Generate再押下
→ pre-Worker no-op
```

### B. Render-only change

```text
Generate A
→ show_scale/font/geometry変更
→ Generate B
→ Worker/Python実行
```

prepared biological cacheはhitしてよいが、drawingは実行する。

### C. Direct edit already applied

```text
Generate A
→ Feature fill/stroke等をdirect edit
→ pending reflowなし
→ Generate
```

P1のcontract上、current artifactがcompleteであればno-op可能。

### D. Direct edit requiring reflow

```text
Generate A
→ geometry creation/removalが必要なlabel/visibility変更
→ pending required reflow
→ Generateまたはauto reflow実行
```

skipしてはならない。

### E. Biological change

```text
source/selector/region/reverse/feature visibility/filter/annotation/orthogroup変更
→ identity miss
→ render
```

注: Feature visibilityのうち、post-render direct overrideとbiological filteringを区別する。

### F. Divergent loaded session

```text
loaded preview A
active intent B
→ first Generate
```

必ずBをrenderする。

### G. Imported no-draft session

```text
load A
→ first Generate full
→ second unchanged Generate no-op
```

### H. Failure/Cancel

failureまたはCancelでcandidate identityをcommittedにしてはならない。

## Performance target

構造的contractを第一とする。

local preferred target:

```text
unchanged warm Generate button settlement < 500 ms
preferred < 250 ms
Worker/Python = 0
```

canonical state snapshot／identity構築が500 msを超える場合は、そのowner内でcompact identity constructionを改善する。final SVGやcatalogをserializeしてkeyを作ってはならない。

## In scope

```text
runtime generation identity
identity schema
compact stable serialization/hash
successful commit時のidentity更新
direct-edit commit時のeditor identity更新
pre-Worker comparison
no-op structured result
History/no-op integration
structural metrics
```

## Out of scope

```text
session schemaへのidentity保存
final SVG cache
Worker result cache
prepared input cache redesign
canonical file serialization一般最適化
new UI
render algorithm変更
```

## Required verification

Focused:

```text
identity stable ordering
field inventory drift
resource token invalidation
editor projection invalidation
pending reflow blocking
failure/Cancel non-publication
```

Browser:

```text
unchanged warm Generate
render-only changed Generate
direct-edit no-op
direct-edit requiring reflow
divergent loaded session
imported first/full second/no-op
```

Broad:

```bash
node --test tests/web/*.test.mjs
npm run test:web:functional-smoke
npm run test:web:perf-smoke
npm run test:web:vibrio-generate
python tests/run_losat_cache_browser_acceptance.py
ruff check gbdraw/
git diff --check
```

## Exit criteria

```text
complete unchanged GenerateがWorker前で終了
false-positive skipを防ぐfield inventory contract
render-only/biological/divergent casesはrender
direct editsとpending reflowを正しく区別
failure/Cancelでidentity未commit
no schema/artifact changes
```

## Expected PR title

```text
Skip unchanged Generate before Worker dispatch
```

---

## 5. Phase間のdependency graph

```text
P0
 └─ P1: direct editとfresh replayを一本化
     └─ P2: active intent schemaを一本化
         └─ P3: P1/P2境界を静的検証
             └─ P4: browser evidence harnessを保守可能にする
                 └─ P5: complete identityで安全にskipする
```

P5は、単に「高速化したいから」先行してはならない。

P1/P2が未完の状態ではgeneration identityに含めるべきfieldが分散しており、必要なGenerateを誤ってskipする危険がある。

---

## 6. 各Codex promptの生成手順

各phaseのpromptは、この計画書をそのままコピーするのではなく、直前phaseの実績を反映して作る。

### Step 1. 現在地を確認

```text
latest main SHA
直前PR merge status
required CI
remaining risks
actual changed owner/module
最新structural metrics
```

### Step 2. Plan statusを更新

本計画書の以下を更新する。

```text
phase status
baseline SHA
completed exit criteria
remaining risk
次phaseへの設計変更
```

### Step 3. Promptを一phaseに限定

promptは必ず以下を含む。

1. prerequisite;
2. mission;
3. existing evidence;
4. explicit scope;
5. explicit non-goals;
6. required owner model;
7. structural contracts;
8. focused tests;
9. broad gates;
10. completion report.

次phaseの実装を明示的に禁止する。

### Step 4. Prompt中の数値を更新

過去の古いperformance値を使い続けない。

```text
最新verified warm Generate
最新session load
最新owner count
最新file/module構造
```

ただし「遅いことを再確認するだけのbenchmark」は要求しない。

### Step 5. 実装後にPlanと照合

Codexのcompletion reportをExit criteriaと照合する。

次の場合はphase未完了とする。

```text
structural ownerが複数残った
testだけ追加してduplicate production pathが残った
scope外変更が混ざった
required CI未完
remaining correctness blockerあり
```

---

## 7. Program-level metrics

P0時点とP5完了時点で比較する。

| 指標 | P0 baseline | P5 target |
|---|---:|---:|
| mounted SVG live-edit commit owner | 複数call-site判断 | 1 |
| fresh SVG complete editor replay owner | 分散 | 1 |
| active config schema owner | `config.js`内 | pure module 1 |
| critical JS boundary static check | なし／限定的 | CI required gate |
| session regeneration largest spec | 2,000行超 | scenario/helper分離 |
| exact unchanged warm Generate Worker run | 1 | 0 |
| exact unchanged warm Generate Python invocation | 1 | 0 |
| direct edit Worker/Python | 0をcontract化済み | 0をsingle ownerで保証 |
| session active intent canonical overwrite | 0 contract | 0をschema ownerで保証 |
| user-visible regression detection | fixture依存 | domain matrix＋private acceptance |

---

## 8. Risk register

### R1. P1が広がりすぎる

**Risk:** 全editor機能refactorへ拡大する。  
**Mitigation:** mounted commitとfresh replay以外は変更しない。scope外domainはinventoryへ記録して別途扱う。  
**Split condition:** production変更が15 filesまたはnet +600 linesを大きく超える場合、P1a mounted commit／P1b fresh replayへ分割する。

### R2. P2 validatorが過度にstrictになる

**Risk:** 実在するcurrent/legacy sessionを拒否する。  
**Mitigation:** current writer fixture、Gallery、v39、bare legacyをpositive fixturesとして維持。compat aliasesは明示的に列挙する。  
**Abort condition:** unknown fieldを受け入れるためにvalidatorを実質無効化する必要が生じた場合、writer schema inventoryを再検討する。

### R3. P3が全体TypeScript化へ膨張する

**Risk:** 数百の既存errorを直す大規模PRになる。  
**Mitigation:** boundary-specific jsconfigを使う。target module以外をincludeしない。  
**Abort condition:** target boundaryを型付けするために大量のunrelated module変更が必要なら、adapter typedefを導入してscopeを戻す。

### R4. P4でtest semanticsが弱くなる

**Risk:** helper抽出時にassertionが失われる。  
**Mitigation:** test count、scenario list、attachment/evidence keyをbefore/after比較する。production diffは0。  
**Abort condition:** CI時間短縮やline数削減のためにsemantic/raster/owner assertionを削る必要がある場合は実施しない。

### R5. P5のfalse-positive no-op

**Risk:** 必要なGenerateをskipし、古い図を表示する。  
**Mitigation:** fail-open、explicit field inventory、P1/P2 ownerからidentityを構築、unknown時はrender。  
**Release blocker:** render-only、biological change、divergent session、pending reflowのいずれかが誤ってskipされる。

### R6. Memory pressure

Prepared biological cacheのretained estimateは大きい。P0–P5では一般memory redesignを行わないが、各PRで新たなpersistent large ownerを増やさない。

P5 identityはcompactでなければならず、SVG、catalog、sequence、full requestのcopyを保持してはならない。

---

## 9. P0–P5で扱わないもの

以下は本計画の後に再評価する。

```text
outer session container streaming/chunking
OPFS project persistence
feature sequence on-demand retrieval
prepared cache retained-memory精密化
canonical file serialization一般cache
Canvas/WebGL renderer
render queue/latest-wins redesign
Web全体TypeScript化
UI internationalization
new biological analyses
```

P5完了後に最新phase timings、memory evidence、user feedbackを用いて次のroadmapを作る。

---

## 10. 全programのDefinition of Done

P0–P5は以下をすべて満たしたとき完了とする。

```text
1. exact user session regressionが閉じている。
2. live editはsingle mounted commit ownerを通る。
3. fresh SVG editor replayはsingle projection ownerを通る。
4. active current-session configはpure schema ownerを持つ。
5. critical JS boundariesがCIで静的に検証される。
6. browser contract harnessがscenario中心に読める。
7. 完全同一GenerateはWorker/Pythonを起動しない。
8. changed Generate、divergent session、pending reflowは誤ってskipされない。
9. Save/Load/Generate/Cancel/Undo/Redo/v39/current contractsが維持される。
10. session/request/catalog schemaとreference artifactsに意図しない変更がない。
```

---

## 11. 次のアクション

現在の次アクションはP1である。

使用するprompt:

```text
Centralize live SVG edit commit and replay
```

P1完了後、completion reportとmerged mainをこの計画書へ反映し、P2 promptを新規作成する。P2以降のpromptは前phaseの実装結果に依存するため、現時点では詳細実装を固定しすぎない。
