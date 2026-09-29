# S03 — Truthful discovery and Circular transform disclosure

日付: 2026-09-27。**S03 の discovery/presentation 実装・対象検証完了。**
BUG-02 / BUG-20 の S03 のみ。BUG-01 / S02 は再導入していない。
S04 以降、authority PR、deploy/tag は実施していない。

## Checkout / authority / measurement identity

専用 clone `/tmp/gbdraw-issue597-S03.Nmsmci`。
Branch `fix/issue-597-input-session-20260926`、upstream
`origin/fix/issue-597-input-session-20260926`。共有 checkout、他 session の環境・server は変更していない。
開始時の取得 target / S01 commit は `15bcbca89392cbd63a88fe80dc44b71ad4868061`。
追加 target commit はなかった。S01 は作業 HEAD の ancestor。
取得 dev `9d967f1c72f730b205420d62d383b73127c2a9b1`、main
`4556e04e929a4a85ad28d1833ce7304bd764881c`。
共通規約に従い clean tree で最新 dev を取り込み、integration merge
`ff02abbae948aaa870f18123805d7297a6c3495d` を作成した。
この dev 取り込みと、後続の S03 実装一 commit を区別する。

Product authority merge `af5d942af60353dda199aa487da9152a3576b3fe` は取得 dev / 作業 HEAD の ancestor。
`python tools/inspect_issue597_s01_contracts.py` は revision 21、PD-OI-044/045 の
承認記録各9項目（計18項目）、active Contract と取得 dev の byte 一致を確認した。
既存選択の再承認、authority / defaults / checker の変更はない。
PD-OI-044 の既存承認を実装する developer preflight として扱い、PD-OI-045 の
全 operation policy 実装は S04 に残す。import Worker permission の別 dev 統合待ちは
依存しない S03 の停止理由にしていない。

実測時 HEAD は上記 integration merge。実測 runtime はその HEAD **と未 commit の S03 差分**であり、
HEAD だけを runtime bytes の証明として使わない。
[validation manifest](../evidence/S03-validation.json) に個別 source / input / environment / artifact SHA-256 を保存した。
manifest hash は各対応 object の sorted-key compact JSON の SHA-256。

| 対象 | SHA-256 |
| --- | --- |
| Source files manifest | `fa0fb817cd7f3cec8d54e2d9130c7b6a17020f7075a5844a2818047a6b68c537` |
| Inputs / deterministic recipes manifest | `2ef8ce2af2aafd574bb67fb128469c024eb4170287d36cc34efede5f7babbadf` |
| Environment object | `c0e6e8f842fe276f5115db5a98b4ebdc634fac80c1ce050bd68271b3099b6e02` |
| Generated browser wheel | `0ef64bf6bfbe3f52be3fa31de57e9bf431e4fb22c88979d297e0e9eff25197f8` |

Python 3.13.3、Node 26.8.2、Node/Python Playwright 1.61.0、Chromium 149.0.7827.55。
通常 sandbox は bubblewrap の `/mnt/wslg/distro` mount エラーで起動不能だったため、
同じ専用 checkout の commands を escalation で実行した。Node Playwright は専用 clone のみに
`npm install --no-save --package-lock=false --ignore-scripts --no-audit --no-fund @playwright/test@1.61.0`
で準備した。tracked dependency manifest / lock は不変。生成 wheel は ignored、cache bust は不変。

## Resulting behavior and owners

| Owner / path | S03 の変更と除去した旧経路 |
| --- | --- |
| `app/record-discovery.js` | 既存 owner に純粋な `circularDiscoveryForInput` を追加。input type と primary/pair の File identity が一致する ready metadata のみ公開。欠けた pair は idle、未確認の完備 source は deferred。 |
| `app/run-analysis.js` / `app/watchers.js` | 重複 source predicate を上記 owner に統合。native upload/replace は自動探索。saved deferred は mode 往復で自動再探索しない。同じ ready/error source の自動再実行を防ぎ、明示 Inspect/Retry は既存 refresh を使う。既存 version/mode/type/File settlement guard と helper single-flight は維持。 |
| `services/config.js` | canonical Circular source 復元直後の偽 loading を deferred / idle に置換。schema、reader、migration、transport、rollback、lock は変更なし。 |
| `app/app-setup.js` / `app/annotations/record-catalog.js` | committed catalog を別 draft の検証済み metadata として流用しない。annotation catalog と表示が共有 predicate を使う。Load 後に開いていた details による自動 refresh を除去。未探索 annotation は Inspect を案内する。 |
| `app/app-setup.js` / `index.html` | Source records の idle/deferred/loading/ready/error と Inspect/Retry を表示。Inspect の enabled 表示と action guard、selector の enabled 表示と setter が各同一 predicate を使う。missing/ambiguous selector、grid、明示 batch は single transforms を許可しない。 |
| `app/app-setup.js` / `index.html` | 既存 native details を `Single-record crop, orientation and titles` に改称。applicable になる source/selector transition で展開。無関係更新で manual close を覆さない。grid/batch 理由と selector は details 外。grid 設定への導線は既存 control を開いて focus する。 |

Inspect は待機中に aria-disabled と同じ action predicate を使うため、操作元の focus を失わない。
source/selector の focus と desktop settings pane / mobile document の scroll anchor を維持する。
selector を rotation rows より前に置き、rows の増減による mobile anchor の移動を防ぐ。
Generate 前の native discovery は手動 refresh を介さない。
探索失敗は filename と Retry/Replace/Remove を示し、入力と既存 Result を保持する。

自動展開は既存 DOM details の open だけを変える。永続 disclosure state、History entry、Generate trigger は追加しない。
selector 操作が既存 grouping intent を変更する以外、source discovery や展開が
selection/grouping/topology/crop draft を自動変更することはない。fresh grid default も不変。

A-01 の ordinary non-increasing owner/path review:
既存 discovery lifecycle と parser/helper は一経路、共有 predicate が二つの重複 identity 判断を置換した。
presentation は既存 setup + native details、canonical request / Result admission は従来 owner のまま。
新 module、parser、schema、renderer、bindings、Circular 複数ファイル入力、compatibility reader、privileged owner はない。
Session writer 44 / request 8 / bindings 2 は不変。
新 export 一件、computed 三件、watch 一件は上記表示・availability・既存 details に限る。
full OE/PE/CB exception sets を要する例外や新 compatibility delivery は導入していない。

## Verification and acceptance

全 commands は専用 checkout で実行。raw logs / captures は
`/tmp/issue597-S03-evidence-Nmsmci`、個別 digest は validation manifest に保存した。
表の出力先はこの directory 以下。server は専用 port か OS 選択 port を使い、自分の server のみ終了した。

| Command / check | Result / log |
| --- | --- |
| `python tools/inspect_issue597_s01_contracts.py` | exit 0。ancestry / revision / 18 receipt fields / fetched dev bytes 一致。`authority.json`。 |
| `python tools/prepare_browser_wheel.py --no-build-isolation` | exit 0。生成 wheel 3,238,664 bytes。`wheel.log`。 |
| `node --test tests/web/session-file.test.mjs tests/web/session-request.test.mjs tests/web/session-active-files.test.mjs tests/web/record-display-options.test.mjs tests/web/record-metadata-inference.test.mjs tests/web/run-analysis-simple-path.test.mjs tests/web/session-draft-authority.test.mjs tests/web/record-selector.test.mjs tests/web/annotations.test.mjs` | exit 0、49 PASS / 0 skip。`node-final.log`。後の unit 差分は空白整理のみ。 |
| `GBDRAW_WEB_TEST_PORT=42974 node node_modules/@playwright/test/cli.js test --config=playwright.functional.config.js --workers=1 tests/web/record-display-discovery.playwright.spec.js tests/web/circular-record-presentation.playwright.spec.js --output=/tmp/issue597-S03-evidence-Nmsmci/browser-acceptance` | exit 0、14 PASS。最終 production bytes。`browser-acceptance.log`。 |
| 同 discovery spec の `--grep 'saved preview remains\|existing helper stays' --output=/tmp/issue597-S03-evidence-Nmsmci/browser-availability` | exit 0、2 PASS。直前 full run 後に加えた History / unavailable action assertions を含む最終 test bytes。`browser-availability.log`。production/input/environment は同一のため他12 tests の証拠を再利用。 |
| `GBDRAW_WEB_TEST_PORT=42976 node node_modules/@playwright/test/cli.js test --config=playwright.functional.config.js --workers=1 tests/web/composite-session-resources.playwright.spec.js --grep 'survive minimal' --output=/tmp/issue597-S03-evidence-Nmsmci/browser-cli-composite` | exit 0、1 PASS。既存三 Vibrio plasmids の CLI writer → Save → fresh Load → Generate、既存 composite source identity / content の保持。`browser-cli-composite.log`。 |
| `pytest tests/ -v -m 'not slow'` | **初回 exit 1**。6,266 PASS / 2 FAIL / 17 skip / 11 deselected。両 FAIL は comparison wrapper の test webServer startup exit 1。`pytest-fast.log`。全体 command を PASS と呼ばない。 |
| `GBDRAW_WEB_TEST_PORT=42975 pytest tests/test_linear_comparison_browser_contracts.py -v -x` | exit 0、上記両 shard を同じ assertions で再実行し 2 PASS。`pytest-comparison-diagnostic.log`。不変の残り6,266 tests の証拠を再利用。共有 server の停止、assertion / timeout の緩和なし。 |
| `ruff check gbdraw/` | exit 0。`ruff-final.log`。 |
| `node tools/check-web-change-budget.mjs --base origin/dev` | local trusted-base Gate **PASS**、Review **REQUIRED**。hard violations 0。`policy-final.log`。checker / detectors / authority は取得 dev と byte 一致。CI PASS や review waiver として扱わない。 |
| `git diff --check` / staged diff review | PASS。production/tests/docs/generated を分けて確認。commit 後も同じ checker の `--base origin/dev --head HEAD` を実行する。 |

| ID | Evidence and boundary |
| --- | --- |
| D-01 | native one/two/duplicate GenBank、DDBJ accession、GFF+FASTA は upload だけで ready、Python Worker 0。既存 topology parser/packaged Worker 比較二 tests も保持。forced lightweight reader miss は既存 helper 実処理を使用し、reply を保持して loading を確認後 settlement。Worker 一件 / helper 二件 / settled 二件。 |
| D-02 | incomplete pair idle、invalid error、Retry/Replace/Remove、rapid replacement/removal、mode/type change、History undo/redo の source restore と stale completion rejection。別 File/type/pair の metadata を公開しない unit assertion。missing/ambiguous selector は disabled、明示 #2 は有効、crop draft は保持。 |
| D-03 | **測定した Circular preview の受入 PASS**。既存 human mitochondrial JSON Load は deferred / 0 records / Python Worker 0、mode 往復・summary 操作でも再探索しない。Inspect は ready / Worker 0 / Result 不変 / History 不増、Generate は inspection 後に一 run。native replacement は自動 ready、旧 catalog を転用せず Result 保持。既存 CLI 三 plasmid gzip journey も PASS。S01 real large Linear preview の Worker-free FAIL は下記のまま残る。 |
| D-04 | applicable single の source/selector transition 展開、manual close 保持、明示 grid/batch の外部理由と導線。1280 / 390px、Enter/Space、upload/selector focus と2px以内の anchor、390px overflow なし。details/Inspect は History を増やさず、自動 grouping/selection/topology/crop 変更なし。 |
| W-01 | latest 同名 target の独立 clone、dev integration と S03 commit を区別、担当 files を明示 stage。同名 remote へ non-force push 後、local/remote SHA 一致を確認する。 |

初期 focused failures は修正してから上記最終 run を取得した。
setup initialization の getter 順序、既存 error prefix の互換、mobile rows 増減の anchor、
承認された自動展開に対応する従来 summary assertion を修正した。
非同期 `page.waitForFunction` に Promise を返してしまう新 test の premature wait は
`expect.poll(() => page.evaluate(async ...))` の実 settlement assertion に変更した。
偽の loading 再開始や、単に待ち時間を増やした成功として扱っていない。

## Visual review / generated files / S07 handoff

`browser-acceptance/**/disclosure-1280.png` / `disclosure-390.png` を読みやすい原寸で確認。
summary、外側 selector/reason、crop/title/reverse controls と keyboard focus は判読可能。
既存 round-trip test が実際に保存した `circular-record-presentation.svg` を Chromium で
1000px幅に描画した `circular-record-presentation.browser.png` も確認した。
BGC0000709 の選択・crop/reverse intent・日本語 title/subtitle は canonical request assertions と一致。
これは internal regression artifact であり公開 showcase の置換には使わない。
Cairo の補助 render は host font fallback により日本語が欠けたため、visual acceptance に使用していない。

focused command で SVG/UI captures を再生成できる。全文 SVG の補助 capture は、同じ app URL で
`window.__GBDRAW_APP__` を待ち、下記を実行後 `page.locator('svg').screenshot(...)` で再現できる。
app を先に開くことで packaged browser fonts を使う。server はこの clone の専用 port で起動する。

```javascript
await page.evaluate(async svg => {
  const root = new DOMParser().parseFromString(svg, 'image/svg+xml').documentElement;
  const width = parseFloat(root.getAttribute('width'));
  const height = parseFloat(root.getAttribute('height'));
  root.setAttribute('viewBox', `0 0 ${width} ${height}`);
  root.setAttribute('width', '1000');
  root.setAttribute('height', String(1000 * height / width));
  document.body.replaceChildren(root);
  document.body.style.cssText = 'background:white;margin:0;padding:0;overflow:visible';
  await document.fonts.ready;
}, savedSvgText);
```

S07 は `docs/REFERENCE/web-app.md` の upload/record-selection 説明と
`docs/TUTORIALS/GUI/highlight-mitochondrial-features.md` の single-record 操作を更新する。
公開 capture 対象は `docs/images/h-gui-12/presentation-settings.png`、
`docs/images/t-gui-10/presentation-settings.png`、必要なら record rows の位置が変わった
`docs/images/h-gui-16/01-record-start.png`。対応する result images は geometry が変わっていないため
再生成の要否を S07 で既存 recipe から判断する。
Gallery tutorials の `.gbdraw-session.json` は不変。S07 は既存
`gbdraw/web/gallery/tutorials/HmmtDNA_basic_circular.json` などの source records / rotation の
操作 crops、caption と新 section label の整合を確認する。public capture の列挙・撮影は S07 の技能/再現規約で実施する。
S03 は public docs/images/Gallery/`examples/gbdraw_social_preview.png` を変更していない。

production 七 files、tests 三 files、計画内 evidence/result 二 files を個別に review。
生成 wheel / screenshots / SVG / CLI session は repo 外または ignored。`tests/reference_outputs` は read-only。
pytest が変更した bundled LOSAT の executable mode はテスト終了後元に戻し、生成 binary の diff は残していない。
active maps、rules、detectors、guard files、dependency manifests、dist/egg-info は不変。

## Remaining boundaries and S04 entry

**Issue #597 全体完了ではない。** S01 の whole-object reply は選定済みであり、S03 は import transport を実装していない。
S01 の full Save/Load performance は FAIL、real large Linear preview の既存 Python/helper Worker construction も
未解消。S03 では再測定していない。小さい Circular preview / 三 plasmid の成功をその代用にしない。
full memory/heartbeat/copy metrics、全 historical browser matrix、六 replicon grid/batch/CLI sidecar journey、
offline runtime audit は今回未測定。S01 performance evidence を変更 runtime の性能証明として再利用していない。

次は **S04**。新しい独立 clone へ同名 remote の最新 HEAD を取得し、
[SESSION_WORKFLOW.md](../SESSION_WORKFLOW.md)、この結果と manifest、
[03_SESSION_OPERATIONS.md](../decisions/03_SESSION_OPERATIONS.md)、
[S04_INSTRUCTION_PROMPT.md](../sessions/S04_INSTRUCTION_PROMPT.md) を読む。
`python tools/inspect_issue597_s01_contracts.py` で既存 authority を再確認し、
S04 の既存 operation owner / availability から semantic-operation exclusivity を実装する。
S03 の Inspect/selector predicate を別 owner に複製せず、共通の operation availability と接続する。
S05 は別 import Worker permission の **dev 統合後**に一方式の transport/lifecycle を実装し、
S06 が full-pipeline 性能、S07 が public docs/screenshots を担当する。
この session は S04 を開始しない。

English commit title: **Expose applicable Circular transforms and truthful discovery states**

English summary: Keep saved Circular previews uninspected until Inspect or Generate, preserve native automatic discovery and source identity, and reveal applicable single-record transforms without changing focus, scroll position, or History.
