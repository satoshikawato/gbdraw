# gbdraw Ponytail Audit 2026-08-09

## 結論

現行 worktree には、約 11,500 行、約 13 MiB、条件つきで直接 npm 依存最大 1 個、CI ジョブ 1 個の削減余地がある。削減対象は、廃止済みアルゴリズム、同じ状態を復元する複数の経路、テスト専用の production API、実装文字列を固定するテスト、Git 履歴で足りる保存物である。

最初に扱うべき P0 は次の 6 件である。

1. Web の実装文字列テストと Python 内生成 Node テストを、既存の直接 JavaScript テストへ寄せる。
2. Python と standalone JavaScript に残る未到達の interactive match payload v1 を削除する。
3. 1 種類しかない orthogroup membership mode の背後に残る旧推論実装を削除する。
4. 利用面を持たない collinearity TSV と旧変換 API を削除する。
5. CLI セッション writer の raw argv 再投影を廃止し、canonical request を唯一の描画設定とする。
6. ドキュメント capture の scenario と fixture 定義を既存 manifest に集約する。

この文書は改善候補と実施方針を記録する。production code、テスト、生成物は変更していない。

## 監査条件

- 対象は `docs_renovation` ブランチの commit `ad38970e` を基点とした、未コミット変更を含む現在の worktree である。既存変更は利用者の作業として保持した。
- Python production 202 ファイル、92,744 行、vendor を除く Web JavaScript 128 ファイル、68,320 行、Python・JavaScript テスト 199 ファイル、113,401 行、`tools/` と documentation capture/recipe の Python 52 ファイル、30,479 行を走査した。行数は空行とコメントを含む物理行数である。
- static import、名前参照、AST の function range、entry point、manifest、workflow、package metadata、公開 API snapshot、session compatibility 文書を照合した。`ruff check gbdraw --select ARG,F401,F841` も補助証拠に使った。個別 finding の行数は主に function range、全体 inventory は物理行なので、合計は重複を控除して丸めた。
- `delete`、`stdlib`、`native`、`yagni`、`shrink` の観点だけを扱う。正しさ、security、performance の一般監査ではない。
- `gbdraw/web/vendor/` と生成済み reference output は通常のコード行集計から外した。ただし配布物に残る未使用 vendor fallback は容量監査の対象にした。
- [examples/gbdraw_social_preview.png](../../examples/gbdraw_social_preview.png) は repository rule に従って対象外とした。

優先度は次の意味で使う。

| 優先度 | 意味 |
| --- | --- |
| P0 | 未到達または二重所有が確認済み。最初の削減対象。 |
| P1 | 削減根拠は強い。対象 surface の focused test を通して実施する。 |
| P2 | 小規模または repository 運用上の整理。P0/P1 の後にまとめる。 |
| P3 | 機会的に削る小さな項目。専用 abstraction は作らない。 |

## 優先順位つき findings

各 finding は 1 行で、削るもの、根拠、置換方針、保守的な削減量を示す。

### P0

| 順位 | 判定 | Finding と方針 | 削減目安 |
| ---: | --- | --- | ---: |
| 1 | `shrink` | [tests/test_web_packaging.py:389](../../tests/test_web_packaging.py#L389) は 5,979 行の中に source-text contract 40 件、1,368 行、Python が production module を複製して inline JavaScript を生成するテスト 8 件、1,791 行、pytest から Node を呼ぶ同型 wrapper 18 件、144 行を持つ。wheel、sdist、CSP、bundle の検査だけを残し、未確認の振る舞いだけを既存の 63 `tests/web/*.test.mjs` と 10 Playwright spec に移し、Node runner を 1 個に parameterize する。 | 1,600〜2,000 行 |
| 2 | `delete` | [gbdraw/render/interactive_svg.py:1950](../../gbdraw/render/interactive_svg.py#L1950) の match payload v1 と private helper 1,073 行、および [standalone-interactivity.js:1349](../../gbdraw/web/js/services/standalone-interactivity.js#L1349) の caller がない V1 island 816 行と別の dead helper 13 行を削除する。両 runtime とも live dispatch は v2 だけを呼ぶため、新しい互換 wrapper は作らない。v1 範囲内でも live v2 が共有する `standaloneMatchKind` は残す。 | 約 1,902 行 |
| 3 | `delete` | [protein_colinearity.py:117](../../gbdraw/analysis/protein_colinearity.py#L117) の membership mode は `anchor_core_v1` だけで、legacy alias はすべて同じ値へ正規化され、[protein_colinearity.py:8066](../../gbdraw/analysis/protein_colinearity.py#L8066) 以下の旧推論 branch は到達不能である。旧 helper island を削除し active path から anchor-core を直接呼ぶ。persisted reader に加え、公開済み normalizer/type shim は API break まで残し、field を外す migration では [api/request_render.py:1451](../../gbdraw/api/request_render.py#L1451) の derived-cache descriptor も更新する。 | 約 1,450 行 |
| 4 | `delete` | 内部 package からも公開されず tests だけが使う [gbdraw/io/collinearity.py:1](../../gbdraw/io/collinearity.py#L1) 463 行、定義と module-local `__all__` 以外に利用がない [analysis/collinearity.py:247](../../gbdraw/analysis/collinearity.py#L247)、[同:320](../../gbdraw/analysis/collinearity.py#L320)、[同:430](../../gbdraw/analysis/collinearity.py#L430)、[同:1299](../../gbdraw/analysis/collinearity.py#L1299) から始まる不連続な変換 root/private helper 406 行、[protein_colinearity.py:434](../../gbdraw/analysis/protein_colinearity.py#L434) の evidence/forwarding helper 102 行を削除し、live block builder だけを残す。native TSV が必要になった時点で supported surface から実装する。 | 約 971 行と専用テスト |
| 5 | `delete` | [session_io.py:2227](../../gbdraw/session_io.py#L2227) は canonical request を必須にしながら、[session_io.py:3743](../../gbdraw/session_io.py#L3743) 以下で raw argv を Web config に再解釈している。Web 側はすでに [session-request.js:2652](../../gbdraw/web/js/services/session-request.js#L2652) で canonical request を投影するため、writer の exclusive projector graph を削除し、raw argv は provenance と legacy replay にだけ残す。 | 約 860 行 |
| 6 | `shrink` | [docs/capture/config.py:47](../capture/config.py#L47) の scenario/screenshot 定義と [docs/capture/config.py:282](../capture/config.py#L282) の fixture metadata は、[docs/scenarios/manifest.json:184](../scenarios/manifest.json#L184) と [tutorial-data/manifest.json:452](../../gbdraw/web/tutorial-data/manifest.json#L452) を再所有している。generic ID lookup と `scenario ID -> capture callable` だけを残し、[docs/capture/run_all.py:126](../capture/run_all.py#L126) の ID/tier 列挙も manifest から読む。 | 550〜700 行 |

### P1

| 順位 | 判定 | Finding と方針 | 削減目安 |
| ---: | --- | --- | ---: |
| 7 | `delete` | 完了済みで参照がない [IMPLEMENTATION_PLAN_WEBAPP_EXPLORATORY_QA_REMEDIATION.md:5](./IMPLEMENTATION_PLAN_WEBAPP_EXPLORATORY_QA_REMEDIATION.md#L5)、[IMPLEMENTATION_PLAN_LINEAR_COMPARISON_GAP_UI.md:5](./IMPLEMENTATION_PLAN_LINEAR_COMPARISON_GAP_UI.md#L5)、[IMPLEMENTATION_PLAN_310_SEVEN_VERTEX_ARROW.md:9](./IMPLEMENTATION_PLAN_310_SEVEN_VERTEX_ARROW.md#L9)、[WEB_GALLERY_SCALE_VISIBILITY_CAPTURE_PLAN.md:4](./WEB_GALLERY_SCALE_VISIBILITY_CAPTURE_PLAN.md#L4) を削除する。完了記録は Git history で保持し、manifest や contributor 文書が参照する記録だけを current tree に残す。 | 2,343 行 |
| 8 | `delete` | production consumer がない Web API を tests のために残さない。[depth-track-state.js:255](../../gbdraw/web/js/app/depth-track-state.js#L255) の旧 remove API と、`search-core.js`、`feature-utils.js`、`results.js`、`file-content-cache.js` などで caller がない export 約 150 行、[depth-file-codec.js:10](../../gbdraw/web/js/services/depth-file-codec.js#L10) の encoder-only path 112 行を対応テストごと削除し、live column removal と decoder は残す。feature visibility rule/cache chain は production consumer があるため対象外とする。 | 約 262 行 |
| 9 | `native` | 8 Playwright spec が独自 static server を持つ。[playwright.config.js:4](../../playwright.config.js#L4) の `webServer` または shared fixture を 1 owner にし、Range や特殊 header が必要な spec だけ専用 handler を残す。代表例は [composition-layout.playwright.spec.js:1](../../tests/web/composition-layout.playwright.spec.js#L1) と [gallery-tutorial.playwright.spec.js:53](../../tests/web/gallery-tutorial.playwright.spec.js#L53) である。 | 250〜320 行 |
| 10 | `delete` | supported API に再公開されず repository caller もない Python helper を削除する。対象は [definition_line_styles.py:193](../../gbdraw/definition_line_styles.py#L193)、[cli_utils/common.py:630](../../gbdraw/cli_utils/common.py#L630)、[svg/circular_ticks.py:281](../../gbdraw/svg/circular_ticks.py#L281)、[svg/linear_features.py:208](../../gbdraw/svg/linear_features.py#L208)、[diagrams/linear/builders.py:374](../../gbdraw/diagrams/linear/builders.py#L374) などで、direct internal import を守る shim は追加しない。 | 約 282 行 |
| 11 | `shrink` | test/capture の exact clone を既存 owner に寄せる。12 個の render stub、[run_cli_scenarios.py:1228](../recipes/run_cli_scenarios.py#L1228) と [run_python_scenarios.py:240](../recipes/run_python_scenarios.py#L240) の 45 行 validator、3 個の `_linear_pair`、6 個の manifest `_chapter` reader が対象である。共通 helper は既存 fixture/module に置き、新しい utility package は作らない。 | 190〜230 行 |
| 12 | `delete` | [examples/README.md:29](../../examples/README.md#L29) が「commit history のため」と明記する未参照 SVG 8 個と専用 recipe を current tree から削除する。対象は `NC_000921_spring.svg`、`NC_000962_psyche.svg`、`NC_001416.svg`、`NC_005042_pine_reflection.svg`、`NC_007205_oceanic_voyage.svg`、`NC_012920_middle_qualifier_priority_inner_axis5_def28_italic.svg`、`NC_016510_mint.svg`、`NZ_CP010822_orange.svg` である。 | 約 187 行、7,455,306 bytes |
| 13 | `shrink` | CLI test runner と input discovery が [tests/conftest.py:121](../../tests/conftest.py#L121)、[test_regression.py:26](../../tests/test_regression.py#L26)、[test_output_comparison.py:89](../../tests/test_output_comparison.py#L89) に重複する。subcommand と任意 BLAST file を受ける runner 1 個と search-path fixture 1 個にする。 | 170〜190 行 |
| 14 | `delete` | [tests/conftest.py:39](../../tests/conftest.py#L39) の marker hook は `pyproject.toml` の marker 定義と重複し、`project_root`、reference path、test-case tables、path-finder fixture は caller がない。`temp_output_dir` は built-in `tmp_path` に置換する。 | 約 150 行 |
| 15 | `native` | [color-utils.js:65](../../gbdraw/web/js/app/color-utils.js#L65) は CSS named color 148 個を手書きで持つ。browser で遅延生成する canvas parser に寄せる場合は、Node import 時に DOM を要求せず、英字名だけを対象にし、無効値を入力のまま返し、既存の uppercase hex 出力を保つ。named-color の直接 Node test は browser test に移す。 | 約 140 行 |
| 16 | `shrink` | [run-analysis.js:1606](../../gbdraw/web/js/app/run-analysis.js#L1606) の 202 行 cancellation snapshot は、[config.js:950](../../gbdraw/web/js/services/config.js#L950) と [history-snapshot.js:425](../../gbdraw/web/js/services/history-snapshot.js#L425) の domain snapshot を再所有する。domain builder/apply を使う generated-artifact transaction を 1 owner にし、result identity、CLI helper local、telemetry、cache/manifest の参照保持だけを cancellation 固有に残す。 | 90〜120 行 |
| 17 | `delete` | [pyproject.toml:37](../../pyproject.toml#L37) の `dev` extra は CairoSVG を含み、matrix job は `.[dev]` を導入済みである。[.github/workflows/test.yml:59](../../.github/workflows/test.yml#L59) の `test-with-cairosvg` は同じ fast suite を再実行するため削除する。CairoSVG 自体は残す。 | 29 行、CI job 1 個 |
| 18 | `yagni` | [phosphor-icons/regular/style.css:3](../../gbdraw/web/vendor/phosphor-icons/regular/style.css#L3) の WOFF、TTF、SVG は現在 fallback URL として参照されるが、supported-browser baseline を WOFF2 対応に固定できるなら WOFF2 だけを残す。参照がない `regular/selection.json` と `browser_wasi_shim/dist/tsconfig.tsbuildinfo` は無条件で外す。 | 約 3 行、最大 6,112,123 bytes |

### P2

| 順位 | 判定 | Finding と方針 | 削減目安 |
| ---: | --- | --- | ---: |
| 19 | `delete` | [tests/utils/svg_compare.py:32](../../tests/utils/svg_compare.py#L32) の未使用 `tolerance`、`ignore_comments`、caller がない `compare_svg_files` と `quick_hash_compare`、空の `tests/utils/__init__.py`、[test_regression.py:569](../../tests/test_regression.py#L569) の skip-only placeholder を削除する。 | 55〜65 行 |
| 20 | `shrink` | clipboard copy が [app-setup.js:2114](../../gbdraw/web/js/app/app-setup.js#L2114)、[orthogroups.js:58](../../gbdraw/web/js/app/orthogroups.js#L58)、[gallery.js:803](../../gbdraw/web/gallery/gallery.js#L803) など 4 箇所にあり、secure context、permission、user activation の扱いも重複する。Clipboard API と非 secure/LAN 用 `execCommand` fallback を持つ helper 1 個にする。 | 50〜60 行 |
| 21 | `delete` | internal function の捨てられる field、state、accessor と同値 branch を削る。主な対象は [radial_layout.py:1705](../../gbdraw/diagrams/circular/radial_layout.py#L1705)、[layout/circular.py:78](../../gbdraw/layout/circular.py#L78)、[layout/linear.py:187](../../gbdraw/layout/linear.py#L187)、[labels/circular.py:4721](../../gbdraw/labels/circular.py#L4721) である。38 個の `ARG` 警告は cache key、互換 callback、P0 projector と重なる `session_io.py:4349,4415` を除外してから call-site ごと縮める。 | 80〜120 行 |
| 22 | `delete` | [gbdraw/web/js/config.js:18](../../gbdraw/web/js/config.js#L18) の `DEBUG` は常に false だが、logger、DI 引数、debug-only DOM query、counter を複数 module に運ぶ。flag と plumbing をまとめて削除する。 | 35〜45 行 |
| 23 | `shrink` | Base64 codec は [file-content-cache.js:34](../../gbdraw/web/js/services/file-content-cache.js#L34)、`gallery-session-migration.js`、`session-request.js`、`config.js`、`session-resources.js` に重複する。bytes/text encode/decode を既存 service 1 箇所に置く。 | 20〜25 行 |
| 24 | `shrink` | pairwise style と collinearity color の値正規化が Python の options、CLI、diagram、legend に重複する。pairwise style は [config/models/objects.py:27](../../gbdraw/config/models/objects.py#L27)、collinearity color は [analysis/collinearity.py:183](../../gbdraw/analysis/collinearity.py#L183) をそれぞれ既存 owner とし、surface ごとの parser を増やさない。 | 約 20 行 |
| 25 | `delete` | [requirements.txt:1](../../requirements.txt#L1) は repository 内に caller がなく、runtime と build dependency を `pyproject.toml` と二重所有する。`pip install -e .[dev]` に統一し、secondary manifest を削除する。 | 8 行 |
### P3

| 順位 | 判定 | Finding と方針 | 削減目安 |
| ---: | --- | --- | ---: |
| 26 | `yagni` | [package.json:11](../../package.json#L11) の `@types/node` は [playwright.config.js:1](../../playwright.config.js#L1) の editor 向け `// @ts-check` を支える可能性がある一方、TypeScript、`tsconfig`、自動 typecheck script はない。editor checking を意図的にやめる場合だけ、直接依存と lockfile の `undici-types` を削除する。 | 条件つき直接依存 1 個 |
| 27 | `shrink` | [index.html:103](../../gbdraw/web/index.html#L103) の未使用 CSS selector と、`text-download.js`、`zip.js`、`session-file.js`、`export.js` の Blob download anchor 手順を削る。download は既存 service の `downloadBlob(blob, filename)` 1 個に寄せる。 | 30〜40 行 |
| 28 | `stdlib` | [render/formats.py:42](../../gbdraw/render/formats.py#L42) の stable dedup loop を `list(dict.fromkeys(formats))` にし、[api/request_render.py:1744](../../gbdraw/api/request_render.py#L1744) の同値 branch を単一 call にする。 | 約 11 行 |
| 29 | `delete` | `missing_inputs.json`、`_reproduced/missing_inputs.json`、`tests/test_inputs/desktop.ini` は consumer がない生成物または OS junk である。削除後、生成先と `desktop.ini` を ignore rule に加える。 | 小、tracked junk 3 個 |

## 実施方針

### 1. まず未到達コードを消す

P0 の interactive v1、旧 orthogroup 推論、collinearity test-only surface を先に削除する。新しい façade、deprecated alias、feature flag は作らない。supported public API snapshot に載らない direct internal import は互換対象にしない。旧 orthogroup 値は persisted reader で `anchor_core_v1` に移すだけにする。

### 2. writer と state owner を 1 個にする

新規 session writer は canonical request からだけ描画設定を作る。raw argv は provenance と released legacy session の replay に限定する。Web の cancellation、history、session import は editor、feature、orthogroup、result の domain capture/apply を共有し、transaction 固有の差分だけを各 caller が持つ。

### 3. テストを実行環境へ戻す

Python packaging test は package、CSP、manifest、bundle topology を検査する。JavaScript の挙動は `tests/web/*.test.mjs`、browser integration は Playwright が所有する。source code の識別子、関数名、control-flow text を Python から固定しない。Node 起動 wrapper と static server は各 1 owner にする。

### 4. Git を履歴として使う

完了した implementation plan と「履歴のため」に残した example は current tree から削除する。再現対象の example、active plan、manifest 参照、release compatibility record は残す。削除する SVG は recipe と manifest allowance も同時に外す。

### 5. 小さな削減に abstraction を追加しない

P2/P3 は既存 owner、標準ライブラリ、browser native API へ寄せる場合だけ実施する。数行の重複を消すための新規 module や generic framework は作らない。

## 検証ゲート

| 変更群 | 最低限の確認 |
| --- | --- |
| Python dead code と公開面 | `ruff check gbdraw/`、`pytest tests/test_public_contract.py tests/test_dead_api_cleanup.py tests/test_collinearity.py tests/test_protein_colinearity.py -q` |
| canonical session writer | `pytest tests/test_session_io.py tests/test_session_request_codec.py tests/test_api_requests.py -q`、session/request の直接 Node tests、session restore Playwright |
| interactive v1 削除 | interactive SVG、standalone interactivity、match popup の focused tests。v1 identifier が production tree に残らないことも `rg` で確認 |
| Web test owner 移行 | 全 `tests/web/*.test.mjs`、対象 Playwright spec、`tests/test_web_packaging.py` の packaging subset |
| docs capture 集約 | capture contract tests、manifest validation、代表 scenario の clean-directory capture |
| example と配布 asset | reproduce-example manifest、Web packaging、offline bundle 検査、配布 tree の file-size inventory |
| 最終 gate | `pytest tests/ -v -m "not slow"`、`pytest tests/test_output_comparison.py::TestOutputComparison -v`。意図しない reference regeneration は行わない |

削減は owner 単位で行い、production、test、documentation、generated artifact の diff を別々に確認する。P0 を一度に混ぜず、各 owner の focused gate が通った後で次へ進む。

## 削除しないもの

- pandas、svgwrite、Biopython、fonttools、bcbio-gff、Python 3.10 用 tomli fallback は production import があるため残す。CairoSVG、Pillow、Playwright、pytest 系、ruff、build、setuptools、wheel にも現行 test/tool/workflow owner がある。
- `setup.py`、`MANIFEST.in`、Conda recipe は別 package format と browser-wheel hook を所有し、packaging tests が使うため残す。
- released session schema の reader/migration、canonical request/resource codec、Web worker boundary は役割が異なる。P0 の raw argv writer projection と混同して削除しない。
- `gbdraw.render.export.save_figure` は 0.16 までの明示済み compatibility path なので、その release までは残す。
- Web の production module は app、Gallery、worker、Cloudflare の entry から到達する。whole orphan module は見つからなかった。standalone runtime も別 execution context として必要であり、削るのは未到達 v1 island と別の dead helper 2 個だけである。
- JSON clone は Vue proxy、`undefined`、JSON-only semantics のため `structuredClone` と同値ではない。native 置換候補にしない。
- active Gallery example、reference output、owner-maintained social preview は今回の削減対象にしない。

## 削減見積り

| 区分 | 保守的な削減 |
| --- | ---: |
| Python production | 約 4,700 行 |
| Web production | 約 1,450 行 |
| tests、docs、capture、CI | 約 5,600 行 |
| 端数・重複の保守的控除 | 約 250 行 |
| 合計 | 約 11,500 行 |
| tracked/deployable assets | 約 13 MiB |
| dependency | editor checking を外す場合に直接 npm 依存最大 1 個 |
| CI | 重複 job 1 個 |

`deps possible` は Finding 26 の条件つき上限を示す。

net: -11500 lines, -1 deps possible.
