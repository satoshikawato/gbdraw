# Webアプリ性能改修計画

Date: 2026-08-10  
Status: active; Phases 0–3 completed/superseded by PR #320 (`6811f42178823cac614aa2f2d37748b22331b18d`); Phase 4 is next; Phase 5 remains measurement-conditional.<br>
Baseline: `docs_renovation` / `48973a0376f8` と、その時点の未コミット変更を含むcurrent worktree  
Scope: `gbdraw/web` の起動、session復元、Undo/Redo、SVG preview、補助Pyodide、色ルール、export読込

関連文書:

- [Web application maintenance](../../gbdraw/web/CLAUDE.md)
- [次セッション改善計画](PONYTAIL_NEXT_SESSION_IMPROVEMENT_PLAN_2026-08-09.md)
- [Session compatibility](../SESSION_COMPATIBILITY.md)

この文書は性能監査結果を実装可能な順序へ落とした計画である。文書作成時点ではproduction code、test、generated artifactを変更しない。

## 1. 結論と目的

現在の主要な待ち時間は、Pythonのdiagram assemblyではなく、ブラウザのmain threadで次を繰り返すことから生じている。

1. Undo履歴のために、設定変更ごとにSVGやfeature catalogを含む全状態をclone・`JSON.stringify()`する。
2. 大きなSVGをsanitize、parse、serializeした後、表示境界でもう一度sanitizeしてDOMを作り直す。
3. 既存凡例の復元で変更が不要な場合にも、main-thread Pyodideを初期化する。
4. Tailwind Play runtime、Vue development build、export専用コードを初期表示から読み込む。
5. 色ルール適用でfeatureとruleを複数回走査し、同じ操作中にSVGを複数回serializeする。

本計画の目的は、機能や互換性を減らさず、不要なコピー、再計算、runtime初期化を削除することである。新しいstate framework、汎用cache layer、常駐workerは導入しない。

## 2. 監査ベースライン

以下は2026-08-09〜10に同じlocal workspaceとChromiumで取得した診断値である。current worktreeが変わり得るため、実装開始時にPhase 0で再測定する。

| シナリオ | 観測値 |
| --- | --- |
| 空のSPA初期表示、Tailwind Playあり | 3回の中央値 `571.3 ms` |
| 空のSPA初期表示、Tailwind Playなし | 3回の中央値 `433.9 ms` |
| `WSSV_genome_comparison.gbdraw-session.json` 復元 | 約 `7.0 s` |
| 同sessionの最大main-thread long task | 約 `2.89 s` |
| 同sessionのSVG結果文字数 | `5,879,614` characters |
| 同sessionのDOMPurify | 約 `291.6 ms` と `326.9 ms` の2回 |
| 同session復元後の通常設定変更 `history.begin()` | `98.8–102.9 ms` |
| 同設定変更 `history.commit()` | `143.0–189.6 ms` |
| 0.87 MiB synthetic snapshot、12回編集 | commitが `31.6 ms` から `129.3 ms` へ増加 |

大規模sessionではartifactがさらに大きい。`vibrio-harveyi-group-collinear` の非圧縮sessionは約285 MiBで、history snapshot対象のresultsとeditor stateだけでも約91 MiBになる。200 MiBの履歴上限は保護になっているが、上限判定のために巨大snapshotを構築・JSON化してから破棄しており、計算自体を防げていない。

## 3. 守る契約

性能改善のために次を変更しない。

- diagram generationとgeneration-time feature extractionはmodule worker内のPyodideで実行する。main threadへ戻さない。
- Webはno-build SPA、same-origin、offline、browser-local data processingを維持する。
- session reader/writer、canonical schema-5 render request、resource authorityを変更しない。
- SVG sanitizerを無効化しない。session、worker result、importを含む全ての非信頼SVGを、live DOMへ挿入する前に共有profileでsanitizeする。
- Undo/Redoの表示、順序、最大action数、file retention、失敗時のstate保全を維持する。
- feature、label、legend、orthogroup、manual layoutの編集結果をsession saveとexportへ反映する。
- Circular/Linear、single/grid/batchのrender geometryとtracked reference outputを変更しない。
- generated browser wheel、Gallery session、reference SVGを手編集しない。
- Tailwind runtimeを外すために外部CDNやruntime network dependencyを追加しない。

## 4. 測定方法と成功判定

### 固定fixture

1. 空のappを新しいbrowser contextで開くcold start。
2. `gbdraw/web/gallery/sessions/WSSV_genome_comparison.gbdraw-session.json` を復元するlarge SVG path。
3. WSSV復元後、同じnumber/select/text controlを10回変更してUndo/Redoするhistory path。
4. 保存済み凡例がすでに一致するWSSV sessionを復元するmain-Pyodide no-op path。
5. 既存の25,000-feature interactive SVG search fixtureをsmall-change regression guardとして使う。

### 記録値

- navigationからVue mount完了まで
- session input開始からpreview SVG、history baseline、editor stateが安定するまで
- `PerformanceObserver`のlong task件数と最大時間
- result identityごとのsanitize回数と合計時間
- history begin/commit/undo/redoの各時間とsnapshot byte estimate
- worker Pyodideとmain-thread Pyodideの初期化回数
- 一つの色ルール操作あたりのfeature scan、serialize、result更新回数

同じbrowser、同じserver、同じfixtureでcold/warmを分け、各シナリオを最低3回実行して中央値を使う。異なるmachine間の絶対値をCI hard gateにしない。構造的な回数制約は自動testにし、時間値はbaseline host上のbefore/after evidenceとして保存する。

### 最終成功条件

| 対象 | 必須条件 |
| --- | --- |
| 通常設定のhistory | WSSVでbegin+commitの各p95が`50 ms`未満。履歴深度1〜10で後半中央値が前半中央値より20%以上悪化しない |
| WSSV session復元 | 同じartifactのfull-string sanitizeは1回以下。中央値をbaselineより20%以上短縮し、最大long taskを40%以上短縮 |
| main-thread Pyodide | 既存凡例が一致するsession復元では初期化0回。diagram workerの初期化は維持 |
| SPA初期表示 | Phase 4前から悪化させず、Phase 4完了後はbaseline中央値より15%以上短縮 |
| 色ルール | 一つのuser actionでresult serialize/updateは最大1回 |
| 機能契約 | focused Node/Playwright、offline packaging、session round-trip、SVG sanitizationが全て通る |

時間目標に届かない場合、閾値を緩めて完了にしない。profileを取り直し、残った最大のmain-thread taskを記録して計画を更新する。

## 5. Phase 0: 再現手順とbaselineを固定する

Status: completed/superseded by PR #320 (`6811f42178823cac614aa2f2d37748b22331b18d`)

### 作業

1. 実装開始時のbranch、HEAD、dirty filesを記録する。

   ```bash
   git status --short --untracked-files=all
   git diff --stat
   git rev-parse --short=12 HEAD
   ```

2. NodeとPythonのPlaywright経路を両方確認する。

   ```bash
   command -v playwright && playwright --version
   node -e "console.log(require.resolve('@playwright/test'))"
   python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
   ```

3. 既存testを変更前に実行し、pre-existing failureを分ける。

   ```bash
   node --test tests/web/history.test.mjs tests/web/history-inputs.test.mjs
   node --test tests/web/svg-sanitization.test.mjs tests/web/pyodide-startup.test.mjs
   pytest tests/test_web_packaging.py -q
   ```

4. 一つのPlaywright performance specに上記固定fixtureの計時とcounterをまとめる。既存の`interactive-svg-search-performance.playwright.spec.js`へ無関係なsession測定を混ぜず、必要なら`tests/web/webapp-performance.playwright.spec.js`を追加する。
5. 診断用counterはtestから注入する。productionへ常駐telemetry、analytics、console spamを追加しない。

### 完了条件

- 監査値をcurrent worktreeで再現するか、差分理由を説明できる。
- 測定開始・完了条件が固定され、任意の`setTimeout`だけでreadyを判定していない。
- browser sandbox failureは必要な権限で同じcheckを再実行し、未検証扱いにしない。

## 6. Phase 1: Historyを全artifact snapshotから分離する（P0）

Status: completed/superseded by PR #320 (`6811f42178823cac614aa2f2d37748b22331b18d`)

### 対象owner

- `gbdraw/web/js/app/history-inputs.js`
- `gbdraw/web/js/services/history.js`
- `gbdraw/web/js/services/history-snapshot.js`
- `gbdraw/web/js/services/history-files.js`
- historyを明示管理するfeature/label/legend/layout action

### 方針

既存の`runUndoableCommand()`を再利用する。全control用の新しいcommand frameworkは作らない。

1. generic form transactionは、変更可能なconfig/draft/UI値だけをbefore/afterとして保持する。SVG、results、feature catalog、sequence source、LOSAT cache、manifestを通常のcontrol focus/changeごとにdeep cloneしない。
2. generated artifactはgeneration/result identityが変わる境界で一度だけcheckpointし、その後のsetting historyからは安定した参照を使う。同一artifactをhistory entryごとに複製しない。
3. feature visibilityなど既にcommand化されている編集は、そのcommandのapply/revertと小さなmetadataを使い続ける。command actionを再びglobal snapshotへ戻さない。
4. Generate、session load、Reset、result switchなど複数domainを置換する操作だけがartifact checkpointを作る。通常のinput adapterと区別する。
5. snapshot signatureとbyte sizeはentry作成時に一度計算して保持する。`enforceLimits()`で全snapshotを再`JSON.stringify()`しない。
6. push/pop/evict時に合計bytesを増減し、file storeの参照解放は現在と同じ時点で行う。
7. `run-analysis.js`のcancellation rollbackをhistory snapshotへ統合しない。async cancellation、stale-result guard、worker lifecycleは別契約のままにする。

### テスト

- 30-action limitと200 MiB limitの既存挙動
- no-op transactionを履歴へ追加しない
- undo失敗時にstackを移動しない
- config edit → undo/redoで設定値だけが戻り、preview artifact identityは不変
- Generate → config edit → undo → result switch → undoの順序
- session load → editor command → undo/redo → session saveでartifactとeditor stateが一致
- file entry eviction後だけfile storeが解放される
- 10回の通常設定変更でsnapshot build、signature、byte estimateが履歴深度に比例して増えない

### 検証

```bash
node --test tests/web/history.test.mjs tests/web/history-inputs.test.mjs
node --test tests/web/session-authority.test.mjs tests/web/session-resources.test.mjs
npx playwright test tests/web/webapp-performance.playwright.spec.js
```

### 中止条件

- Undo後に別resultのSVG、feature catalog、resource bindingが混ざる場合は、参照共有を続けずartifact identity ownerを先に修正する。
- mutable objectを複数history entryで共有し、後の編集が過去entryを書き換える場合は完了扱いにしない。
- 履歴制限を小さくする、Undo対象を黙って減らす、巨大sessionだけhistoryを無効化する回避策は採用しない。

## 7. Phase 2: SVGを一度だけsanitizeして一度だけcommitする（P0）

Status: completed/superseded by PR #320 (`6811f42178823cac614aa2f2d37748b22331b18d`)

### 対象owner

- `gbdraw/web/js/services/svg-sanitization.js`
- `gbdraw/web/js/services/result-normalization.js`
- `gbdraw/web/js/app/candidate-render.js`
- `gbdraw/web/js/state.js`
- `gbdraw/web/js/app/watchers.js`
- `gbdraw/web/js/services/svg-serialization.js`
- `gbdraw/web/js/app/svg-styles.js`

### 方針

1. worker result、session result、imported resultがstateへ入る一つのingestion boundaryを決める。
2. non-live DOMで必要なoverrideを適用し、最終文字列を共有profileで一度sanitizeしてからcommitted resultへ保存する。
3. committed resultは「sanitize済み」という内部契約にし、`state.svgContent`表示時の同一文字列への再sanitizeを削除する。session loadを含む全ingressがboundaryを通るtestを先に置く。
4. overrideがなくresultが既にcurrent identityなら、candidate parse/serializeを行わない。
5. result commit後のlegend extraction、composition binding、drag、feature handler indexは一回のpost-commit処理にまとめる。`extractedFeatures` watcherから同じhandler indexを直後に再構築しない。
6. 一つのuser action内でpalette、specific rule、stroke処理をまとめ、最後に最大一回だけ`serializeCleanSvg()`とresult更新を行う。
7. このphaseではlive DOMを長期の唯一の正本にしない。重複sanitizeと重複serializeを除いても目標に届かない場合だけ、dirty DOMをexport/session save時にserializeする次段階を別計画にする。

### セキュリティテスト

次の入力をworker result、session result、candidate resultの各入口で確認する。

- `<script>`、event handler attribute、`foreignObject`
- `javascript:`や外部URLなど禁止scheme
- safe SVG geometry、semantic IDs、metadata、allowed styles
- sanitizerを通らずlive previewへ到達する経路がないこと

### 機能テスト

- Circular/Linear session load、result switch、Generate、Reset
- feature color、visibility、label、legend、manual dragの編集
- clean SVG、interactive SVG、PNG、PDF export
- undo/redoとsession save/load後の編集保持
- 一つのaccepted result identityにつきfull-string sanitize一回以下
- 一つの編集actionにつきserialize/result update一回以下

### 検証

```bash
node --test tests/web/svg-sanitization.test.mjs tests/web/svg-styles-track-groups.test.mjs
node --test tests/web/feature-color-actions.test.mjs tests/web/legend-sync.test.mjs
npx playwright test tests/web/linear-multi-record.playwright.spec.js
npx playwright test tests/web/interactive-svg-v3.playwright.spec.js
npx playwright test tests/web/webapp-performance.playwright.spec.js
```

### 中止条件

- sanitizer bypassを`trusted=true`のcaller任せにしない。信頼フラグが必要なら、その発行ownerと偽造不能なresult identityを先に定義する。
- DOMPurifyを軽量な正規表現置換へ変えない。
- performance testを通すためにlarge SVG機能、metadata、interactive hooksを削らない。

## 8. Phase 3: No-op session復元でmain-thread Pyodideを起動しない（P0）

Status: completed/superseded by PR #320 (`6811f42178823cac614aa2f2d37748b22331b18d`)

### 対象owner

- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/js/app/legend/entry-actions.js`
- `gbdraw/web/js/app/pyodide.js`
- 必要に応じてrecord discovery caller

### 方針

1. 保存済み凡例のkey、legacy label、color一致判定をJavaScriptだけで先に行う。
2. 追加・再計測が本当に必要なentryが0件なら`ensurePyodide()`を呼ばない。
3. session restore側からの先行`initializePyodide()`を削除し、必要なhelper action自身だけがlazy initializationを所有する。
4. 最初の修正ではPython helper全体をworkerへ移さない。no-op guardだけで測定し、main runtimeがなお大きな待ち時間を占める場合に限り、既存diagram workerへの小さなhelper messageまたはbrowser text metricsを別途比較する。
5. diagram-generation workerのpreinitializeとmain runtimeの状態表示を混同しない。

### テスト

- 既存凡例が一致するsession: main `loadPyodide()` 0回、凡例変更なし
- 一件だけ新しい凡例が必要: main runtime 1回、重複entryなし
- helper失敗: no-op pathは影響を受けず、required pathは明示error
- session restore中のdiagram worker ready/cancel/stale-result behaviorを維持

### 検証

```bash
node --test tests/web/pyodide-startup.test.mjs tests/web/legend-sync.test.mjs
npx playwright test tests/web/depth-track-session.playwright.spec.js
npx playwright test tests/web/webapp-performance.playwright.spec.js
```

## 9. Phase 4: 開発用runtimeと未使用export payloadを初期経路から外す（P1）

Status: next (current scope: 4C / P1-1a lazy export payload loading only)

### 4A. Vue production build

1. 同一versionの`vue.global.prod.js`をsame-origin vendor assetとして追加する。
2. `index.html`の参照とpackaging manifest/testを更新する。
3. development warning、devtools hooks、template behavior、CSP、offline起動を確認する。

### 4B. Tailwind Play runtime

1. static classとdynamic classの全sourceをinventoryし、current UIに必要なCSSをmaintenance-timeに生成する。
2. 生成したminified CSSをversioned same-origin assetとしてtrackする。runtime Tailwind compilerを配布しない。
3. 再生成手順は既存のrepository toolまたは一つの狭いscriptに置く。新しいfrontend build system、runtime npm install、CDNを追加しない。
4. arbitrary/dynamic classが静的生成から落ちる場合は、狭いsafelistまたはplain CSSへ置き換える。巨大なcatch-all safelistは使わない。
5. desktop/mobile、dark/light、popup、drag state、feature search highlightをvisual checkする。

### 4C. Export payloadのlazy load

1. `app-setup.js`から`services/export.js`のeager importを外し、最初のdownload actionでdynamic importする。
2. `standalone-interactivity-assets.js`はinteractive SVG export時だけ読む。
3. jsPDFとsvg2pdfのvendored UMD scriptはPDF export時に一度だけ読み込む。既存loaderがなければ、同じURLのPromiseを共有する最小helperだけを追加する。
4. SVG/PNG exportにPDF libraryを読み込まない。load failureはdownload errorとして表示し、silent fallbackしない。

### 検証

```bash
pytest tests/test_web_packaging.py -q
npx playwright test tests/web/interactive-svg-v3.playwright.spec.js
npx playwright test tests/web/webapp-performance.playwright.spec.js
python tools/verify_gui_offline.py
```

### 完了条件

- 初期pageにTailwind compiler、Vue development build、jsPDF、svg2pdf、standalone interactivity assetがない。
- runtime requestは全てsame-originで、offline起動と全export形式が成功する。
- Tailwindを無効化しただけではWSSV loadの主因は直らないため、Phase 1〜3を飛ばしてPhase 4だけで完了にしない。

## 10. Phase 5: 色ルール走査を一回にまとめる（P1、測定条件付き）

Status: measurement-conditional

Phase 2完了後、large feature fixtureで色ルール処理が残りのlong task上位にある場合だけ実行する。

### 方針

1. 一つのactionでmanual/hash/non-hash ruleを一度compileし、featureごとの最終色を一回だけ決める。
2. palette適用とspecific rule適用が同じfeature setを別々に走査しないよう、既存のmatching helperを一つのcallerから使う。
3. DOM差分を一回適用し、serialize/result updateを最後に一回だけ行う。
4. regex cache、rule engine、worker並列化はprofileで必要性が示されるまで追加しない。

### テスト

- rule precedence、hash/manual override、regex、disabled rule、track group
- 同じ入力から同じ色とsession stateが得られる
- feature×rule match回数とserialize回数のcounter
- 25,000-feature fixtureでbaseline比の改善、search性能の非回帰

### 検証

```bash
node --test tests/web/feature-color-actions.test.mjs tests/web/svg-styles-track-groups.test.mjs
npx playwright test tests/web/interactive-svg-search-performance.playwright.spec.js
npx playwright test tests/web/webapp-performance.playwright.spec.js
```

### スキップ条件

Phase 1〜4後のprofileで色ルールがlong task上位でなく、一action一serializeを既に達成している場合は、このphaseを`not needed`として測定値を記録し、追加実装しない。

## 11. Phase 6: 統合検証と差分監査

Status: pending

### Focused gate

```bash
node --test tests/web/*.test.mjs
npx playwright test tests/web/webapp-performance.playwright.spec.js
npx playwright test tests/web/interactive-svg-search-performance.playwright.spec.js
pytest tests/test_web_packaging.py tests/test_gallery_session_semantics.py -q
python tools/verify_gui_offline.py
```

Nodeの`@playwright/test`がない場合は、同じflowをPython Playwrightで実行する。browser launchがsandboxで失敗した場合は同じcheckを権限付きで再実行する。

### Broader gate

Webのsession/request ownerまたはPython wheelを変更した場合だけ、次を追加する。

```bash
python tools/prepare_browser_wheel.py
pytest tests/ -v -m "not slow"
python -m build
```

cache-bust tokenはdeployable bundleを作る明示的な作業でのみrefreshする。

### 差分監査

1. production、test、documentation、generated artifactのdiffを分けて確認する。
2. `tests/reference_outputs/`、Gallery owner assets、browser wheelに意図しない変更がないことを確認する。
3. sanitizer profile、CSP、worker lifecycle、session schema、canonical requestに意図しない差分がないことを確認する。
4. Phase 0と同じ条件でbefore/after表を更新する。
5. 未達の性能目標、残った最大long task、意図的に見送ったphaseをStatus欄へ記録する。

## 12. 実行順と独立性

実装は次の順序で進め、各phaseのfocused gateが通るまで次へ進まない。

```text
Phase 0 baseline
  → Phase 1 history
  → Phase 2 SVG ingestion/commit
  → Phase 3 Pyodide no-op
  → Phase 4 production assets/lazy export
  → Phase 5 color rules（profileで必要な場合のみ）
  → Phase 6 integration evidence
```

Phase 1とPhase 2は最も大きな体感改善を狙う。Phase 3は小さい変更で大きなsession復元コストを避ける。Phase 4はstartupとproduction hygieneの改善であり、large SVGの2.9秒taskの代替修正ではない。

## 13. 完了の定義

次を全て満たしたときだけ本計画をcompleteにする。

- Phase 4までの必須phaseが実装され、Phase 5は実装済みまたは測定根拠付きで`not needed`になっている。
- Section 4の性能条件を同一fixture・同一browser条件で満たす。
- history、session、SVG sanitization、main/worker Pyodide、全export形式、offline起動の回帰testが通る。
- Circular/Linearとinteractive editorのvisible behaviorが維持される。
- reference outputとGallery assetに意図しない変更がない。
- before/after evidence、test command、結果、残課題がこの文書へ追記されている。

## 14. やらないこと

- WebをReact等へ移行する、state managerを追加する、JavaScript build systemを全面導入する。
- main threadへdiagram generationを戻す。
- sanitizerを削除する、信頼境界を曖昧なbooleanだけで迂回する。
- Undo回数を減らす、large sessionではUndoを無効化する。
- 全SVGを常時canvasへ変換する、virtual DOMを独自実装する。
- performance telemetryやuser dataを外部へ送る。
- 性能改善と無関係なsession schema、Python public API、render geometryを同時に変更する。
- benchmarkを速く見せるためfixture、feature metadata、interactive behaviorを削る。

この計画では、追加機構より削除を優先する。履歴の全量コピー、SVGの二重境界、不要なruntime起動、初期経路の未使用payloadを取り除き、それでも残った処理だけを次のprofile対象にする。
