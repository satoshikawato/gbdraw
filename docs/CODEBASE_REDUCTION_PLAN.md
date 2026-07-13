# Codebase Reduction Refactoring Plan

- 作成日: 2026-07-13
- 対象: `gbdraw/`, `tests/`, `tools/`, `docs/`, `examples/`
- 目的: 公開機能、出力、互換性、配布形態を変えずに、保守対象のコード量とリポジトリ容量を減らす
- 方針: YAGNI による削除を最優先し、KISS を保てる範囲で DRY と SOLID を適用する

## 1. 成功条件

この計画における「機能に影響を与えない」とは、少なくとも次の条件がすべて維持されることを指す。

1. Python 公開 API の import path、`__all__`、関数・クラスのシグネチャ、既定値が変わらない。
2. Circular/Linear CLI のオプション、別名、既定値、help、生成される `Namespace` が変わらない。
3. 同一入力と同一設定に対する SVG とインタラクティブ SVG の構造・意味・操作結果が変わらない。
4. PNG、PDF、EPS、PS の生成経路と対応形式が変わらない。
5. 既存セッション、旧バージョンのセッション、ギャラリー例を読み込める。
6. Web UI の単一ページ・ビルド不要・オフライン動作という制約を維持する。
7. browser wheel、Cloudflare Pages 用成果物、Python package の内容と動作を維持する。
8. テストケースの削減は、検証する入力・分岐・期待値の削減を意味しない。
9. 各 PR は原則としてテキスト行数または追跡対象バイト数が純減する。

単なるファイル分割やコード移動は削減として数えない。抽象化を追加する場合は、同じ PR 内で重複実装を削除し、差し引きで純減させる。

## 2. 監査スナップショット

数値は 2026-07-13 時点の作業ツリーを静的解析した概算である。生成物や vendor を含む値と、保守対象ソースの値は分けて扱う。

| 区分 | 規模・検出結果 | 解釈 |
|---|---:|---|
| Production Python | 約 55,999 行 | 主な保守対象 |
| Python tests | 約 37,141 行 | 重複削減余地はあるが検証範囲を優先 |
| Web source（vendor 除外） | 約 55,678 行 | HTML、JS、CSS、埋め込み Python を含む |
| Tools | 約 5,998 行 | loader と生成処理に小規模な重複あり |
| Docs | 約 7,622 行 | 完了済み計画書が 1,375 行 |
| 未参照の private Python 定義候補 | 30 定義、762 行 | 動的参照を確認してから削除 |
| Production Python の同一関数本体 | 30 群、重複側 432 行 | 共通化の総量上限であり純減値ではない |
| Web JS の同一関数本体 | 59 群、重複側 702 行 | worker 間重複を含む |
| Python test の同一関数本体 | 27 群、重複側 624 行 | fixture/factory/parameterize の候補 |
| Tools の同一関数本体 | 1 群、重複側 9 行 | 低リスク |
| 完全同一ファイル | 30 群、余剰 35 ファイル、約 50.7 MB | パス互換性と配布内容の確認が必要 |
| 生成 palette SVG | 110 ファイル、約 399 MB | docs の URL と表示を維持できる場合のみ整理 |

### 2.1 構造上の主要所見

- `gbdraw.features.colors` と `gbdraw.features.visibility` に import cycle が 1 件ある。
- `gbdraw/layout/scalar_axis.py` に既存の汎用目盛り関数がある一方、Circular/Linear depth drawer が同等処理を再実装している。
- Circular/Linear CLI には共通名のオプションが 77 個あり、そのうち 55 個は `add_argument` 呼び出しが同一である。
- `cli_utils/common.py` には共通引数 helper が既にあるが、Circular/Linear の parser から十分に再利用されていない。
- Web の全 79 ES module は `app.js` から到達可能であり、ファイル単位で明白な dead module は確認できなかった。
- `run-analysis.js`、`app-setup.js`、`services/config.js`、`python-helpers.js` が大きい。ただし、分割だけでは行数が減らないため、サイズだけを理由に再編しない。
- `tests/test_web_packaging.py` は約 3,888 行、約 832 個の assert を持ち、実装文字列への依存が多い。安全な共通化を妨げるため、外部挙動を検証するテストへ段階的に置換する。

### 2.2 監査時のテスト基準値

実行済みコマンド:

```bash
/home/kawato/micromamba/envs/gbdraw-dev/bin/python -m pytest tests/ -q -m "not slow"
```

結果:

- 1111 passed
- 24 skipped
- 2 failed
- 6 deselected

失敗 2 件はいずれも環境依存だった。

- CairoSVG がないため PNG 変換を完了できない。
- `libcairo.so.2` がないため CLI の smoke render を完了できない。

監査環境には Node と Ruff がなく、Web unit test と lint は未実行である。改修前に依存関係を揃え、これらを含む基準値を固定する。

## 3. 原則の適用ルール

### 3.1 YAGNI

- 到達不能、未参照、完了済み、生成可能なものを先に削除する。
- 将来の拡張だけを理由に残された private helper は保持しない。
- 「いつか使うかもしれない」抽象化は作らない。
- 削除判断は文字列検索だけで行わず、import、callback、monkeypatch、template、Pyodide 埋め込み、配布 manifest からの動的参照を確認する。

### 3.2 KISS

- 既存の明示的 helper を再利用し、registry、DSL、reflection、万能 utility を導入しない。
- 1～2 箇所しか使わず、共通化後に条件分岐が増える処理は統合しない。
- 設定 dataclass の明示的な `from_dict` は、反射的な汎用化より読みやすいため原則維持する。
- 巨大関数は、責務を移動するだけでは分割しない。重複や不要分岐を除去できる単位でのみ縮小する。

### 3.3 DRY

- 完全同一、または入力・出力契約が同一と証明できる処理だけを共通化する。
- 共通化の優先順位は「既存 helper の利用」「同一領域内の小 helper」「新しい共有 module」の順とする。
- Circular/Linear、main thread/worker、Python/Web の違いを隠すための boolean flag は増やさない。
- 共通化後も呼び出し側の意味が名前から読める API にする。

### 3.4 SOLID

- SRP: parsing、normalization、rendering、I/O の境界を維持し、重複処理をその責務の owner に集約する。
- OCP: 公開 API や既存形式を変更せず、内部 helper の差し替えで対応する。
- LSP: drawer、config、track 等の既存差し替え可能性とテスト double を壊さない。
- ISP: 呼び出し側が不要な巨大 context を受け取る共通 helper は作らない。
- DIP: 高位 API から Web/CLI 固有実装への逆依存を増やさず、既存の `web_support` や低位 utility に依存方向を揃える。

## 4. 実施順序

各段階は独立した小さい PR にする。後続段階は前段階の完了を前提とするが、条件を満たさない候補は無理に実施せず保留する。

### Phase 0: 非回帰契約の固定

目的は、削減前後の同値性を機械的に判定できる状態にすることである。

実施内容:

1. `gbdraw.api.__all__`、主要 module の公開 symbol、関数シグネチャを JSON 化して比較可能にする。
2. Circular/Linear CLI について、代表引数と全既定値から生成される `Namespace` と `--help` を保存・比較する。
3. reference SVG、interactive SVG、session、gallery、offline bundle の基準値を固定する。
4. Node/Ruff/Cairo の不足を解消し、全検証コマンドの改修前結果を記録する。
5. 追加の characterization test が必要な場合は、対象削減と同じ PR で導入し、PR 全体では純減させる。

完了条件:

- 後述の共通ゲートが現行実装で通る。
- 環境依存で通せない項目は理由、代替検証、CI 上の実行場所が明記される。

### Phase 1: 削除のみで減らす

#### 1-A. 未参照 private 定義

静的解析で未参照と判定された次の 30 定義を、動的参照確認後に削除する。

| ファイル | 候補 |
|---|---|
| `analysis/collinearity.py` | `_representative_proteins`, `_is_chain_accepted`, `_find_best_chain` |
| `analysis/protein_colinearity.py` | `_feature_hash_inputs`, `_select_distribution_split_orthogroup_edges_from_directional_hits` |
| `api/diagram.py` | `_linear_slots_define_renderer`, `_dinucleotides_from_linear_slots`, `_max_depth_from_dataframes`, `_load_depth_table`, `_load_depth_tables`, `_depth_track_count_from_data`, `_resolve_definition_position` |
| `features/tracks.py` | `_find_best_track_split_overlaps_by_strand` |
| `labels/circular.py` | `_leader_anchor_candidates` |
| `diagrams/circular/assemble.py` | `_label_overlaps_other_labels`, `_feature_track_ratio_factor_from_draw_width`, `_default_numeric_slot_width_px`, `_default_gc_skew_layout_without_depth`, `_default_gc_skew_gap_px`, `_distributed_lane_specs_between_bounds`, `_default_track_center_radius_px`, `_arena_from_center_and_width`, `_resolve_ratio_annulus_center_avoiding_forbidden_bands`, `_resolve_fixed_width_annulus_center_avoiding_forbidden_bands`, `_shrink_fixed_width_annulus_to_avoid_forbidden_bands`, `_resolve_ticks_center_radius_avoiding_feature_band` |
| `diagrams/circular/radial_layout.py` | `_scaled_inside_auto_spacing` |
| `diagrams/linear/assemble.py` | `_clone_skew_config_with_dinucleotide`, `_precalculate_depth_dataframes`, `_apply_shared_depth_axis` |

候補ごとの削除条件:

- Python/JS/template/文字列/manifest から参照されない。
- `getattr`、entry point、callback 登録、test monkeypatch の対象ではない。
- 削除前後で public symbol snapshot と全テストが一致する。

#### 1-B. その他の明白な不要物

- `cli.py` の空の `if ...: pass` を削除する。
- 他文書や navigation から未参照で、実装完了済みであることを確認して次の計画書を削除する。
  - `docs/CLI_TSV_INPUTS_IMPLEMENTATION_PLAN.md`
  - `docs/CLI_TSV_AUDIT_FIX_PLAN.md`
  - `docs/TUTORIALS_AUDIT_IMPLEMENTATION_PLAN.md`
- Web の CSS/HTML/JS と packaging test から未参照であることを確認し、Web vendor 配下の重複 Liberation Sans を削除する。Python 側で必要な font は維持する。

見込み:

- テキスト約 2,100 行の純減。
- font 約 1.65 MB の削減。

### Phase 2: 既存の Python 共通部品へ寄せる

新しい抽象化を作る前に、既に存在する owner へ重複処理を集約する。

#### 2-A. Scalar axis

Circular/Linear depth drawer の次の処理を `layout/scalar_axis.py` の既存実装へ置換する。

- major tick の計算
- minor tick の計算
- fraction scaling
- percent tick formatting
- axis bounds の共通部分

対象となる既存 API:

- `scalar_axis_tick_values`
- `scalar_axis_small_tick_values`
- `scaled_scalar_fraction`
- `format_percent_tick`

出力する値の列、丸め、端点、空入力、log/linear の差を table-driven test で前後比較する。

#### 2-B. 入力・selector helper

次の完全同一または同一契約の処理を、それぞれの責務を持つ既存 module に統合する。

- records table の読み込みと record-major depth file の解決
- coordinate map reader の 4 重複
- feature selector と label filtering の qualifier/selector 正規化
- depth table の読み込み、index 検証、既定 axis の計算
- slot の小規模な正規化処理

#### 2-C. Import cycle

`features.colors` と `features.visibility` の双方が使う低位の正規化処理を、既存の `features.selector_values` 等の適切な owner に移し、循環依存を解消する。循環回避専用の新 module は、純減にならない限り作らない。

見込み:

- Production Python で約 400～650 行の純減。

### Phase 3: CLI parser の重複削減

Circular/Linear CLI の共通引数を `cli_utils/common.py` の既存 helper から構築する。

実施内容:

1. 両 parser の `add_argument` を option 名、順序、action、type、default、choices、help、metavar 単位で機械比較する。
2. 完全一致する 55 呼び出しから、既存 helper で表現できるグループを移す。
3. `add_input_args`、`add_output_args` 等の既存 helper を優先して使う。
4. Circular/Linear 固有の引数は各 parser に明示的に残す。
5. 汎用 option schema、decorator、reflection ベースの parser generator は導入しない。

非回帰条件:

- `--help` の option 順序と文言が一致する。
- 省略時と代表的な全指定時の `Namespace` が型を含めて一致する。
- 既存の別名、deprecated option、エラー終了コードと stderr が一致する。

見込み:

- 約 250～450 行の純減。実測で純減しないグループは共通化しない。

### Phase 4: Web の完全同一 utility を統合する

Web の全 module が到達可能であるため、ファイル削除ではなく関数単位で重複を減らす。

優先候補:

- `concatUint8Arrays` の 4 重複
- `normalizeOptionalText`
- `normalizeStringArray`
- orthogroup membership mode の正規化
- transform parsing
- plot title の正規化
- JSON clone
- DOM class 操作 helper
- text download
- SVG `defs` の確保

ルール:

- `utils.js` のような無制限の置き場は作らず、codec、DOM、SVG、LOSAT など責務別の既存 module に置く。
- 呼び出しが 2 箇所だけで短く、import 行を含めると純減しない処理は残す。
- browser と Node の双方で import できる副作用なしの ES module に限定する。

### Phase 5: LOSAT と worker の重複削減

module worker が ES module を import できる構成を利用し、main thread と worker の同一 runtime 処理を共有する。

候補:

- `runLosatPairDirect`
- `instantiateDirectLosat`
- Uint8Array の連結・変換
- worker RPC の request/result/error/cleanup lifecycle

実施条件:

- threaded/non-threaded、WASI、fallback の各経路を個別にテストできる。
- worker の URL 解決、CSP、offline packaging、Cloudflare bundle が維持される。
- 共通 module が環境判定の条件分岐だらけにならない。

Phase 4～5 の見込み:

- Web JS で合計約 350～700 行の純減。

### Phase 6: Web UI と Python/Web 境界の限定的整理

この段階は、削減量を実測できる項目だけを実施する。

#### 6-A. HTML の反復

- `index.html` に 49 箇所ある `auto-value-field` を、既存の Vue component 方針に沿って限定的に component 化する。
- 反復する modal shell と track field も、slot/props が少なく、生成 markup が同一になる場合だけ component 化する。
- 汎用 form renderer や JSON schema driven UI は導入しない。

#### 6-B. 埋め込み Python

`python-helpers.js` 内の埋め込み Python と、`linear.py`、`web_support/feature_metadata.py`、`web_support/orthogroup_metadata.py` の重複責務を確認する。

- Web 固有の adapter は `web_support` に置く。
- 既存 serializer、ID、coordinate helper を再利用する。
- Pyodide への module 配布に新たな複雑性が生じる場合は統合しない。
- 単なるコード移動は行わず、同一ロジックを 1 実装削除できる場合に限定する。

#### 6-C. 巨大関数

`run-analysis.js`、`app-setup.js`、`api/diagram.py`、Circular/Linear assemble は、次の条件を満たす小単位だけを抽出する。

- 2 箇所以上の重複を消せる。
- 引数が少なく責務が明確になる。
- 抽出後の総行数が減る。

可読性目的の大規模な分割は、本計画の削減 PR とは分離する。

見込み:

- 条件を満たす場合のみ約 200～500 行の純減。

### Phase 7: テストコードの重複削減

検証範囲を維持したまま、fixture と test construction の重複を減らす。

実施内容:

- Circular track/feature test の `fake_add_gc_content_group_on_canvas`、`fake_add_gc_skew_group_on_canvas`、capture helper を `tests` 内の共有 fixture/factory にする。
- 同一 body で入力と期待値だけが異なる test を `pytest.mark.parametrize` にする。
- Web の setup、DOM fixture、worker mock を責務別 helper にする。
- `tests/test_web_packaging.py` の実装文字列 assert は、可能な箇所から Node/Playwright の外部挙動テストへ置換する。

維持するもの:

- 入力値、境界値、例外、期待出力の組み合わせ数。
- CSP、asset inclusion、offline packaging の検証。
- reference output と session compatibility の検証。
- 回帰バグごとの意図が分かる test 名または case id。

見込み:

- Test code で約 350～500 行の純減。

### Phase 8: リポジトリ容量の削減

コード行数の削減とは別 PR とし、配布パスや docs URL を壊さないことを優先する。

#### 8-A. 完全同一 fixture/input

完全同一の 30 群、余剰 35 ファイル、約 50.7 MB について、参照元、package inclusion、docs からの利用を一覧化する。

- 内部 test 専用でパスが契約でない場合は、参照を canonical file に統一して余剰 copy を削除する。
- `examples/` の公開パスや tutorial command が参照するファイルは維持する。
- symlink は Windows、sdist/wheel、GitHub ZIP、PyPI packaging の互換性を確認できない限り使わない。

最大の候補は `examples/NC_010162.gb` と `tests/test_inputs/NC_010162.gb` の重複約 30.8 MB である。公開 example path と test path のどちらも契約として必要なら削除しない。

#### 8-B. 生成 palette SVG

約 399 MB の palette SVG は docs から参照されているため、直ちに削除しない。次のすべてを満たせる場合のみ追跡対象から外す。

- 同一内容を決定的に再生成できる。
- docs の URL と表示結果を維持できる。
- release、gallery、offline bundle の生成工程で自動生成できる。
- CI が欠落と再現性を検証できる。

容量削減の見込み:

- 安全確認済みの完全同一ファイルと font: 最大約 52.4 MB。
- palette SVG: 追加で最大約 399 MB。ただし上記条件を満たす場合に限る。

## 5. 実施しないもの

次の項目は一見削減候補に見えるが、現状では互換性または機能上の役割があるため対象外とする。

- `circular_diagram_components.py`: 後方互換 shim として維持する。
- Circular label の legacy fallback: dense label 配置で現在も使用されるため維持する。
- `svg/`, `layout/`, `tracks/`: 参照中であり、旧構造という理由だけでは削除しない。
- `standalone-interactivity-assets.js`: standalone SVG の runtime asset として維持する。
- offline 用 vendor、Pyodide、WASM: 未使用と証明できた個別 asset 以外は削除しない。
- portable session の旧 field、alias、migration: 互換性維持のため削除しない。
- reference SVG: 出力同値性の oracle として維持する。
- config dataclass の明示的 `from_dict`: reflection 化による見かけ上の削減は行わない。

## 6. PR 分割案

| PR | 内容 | 主なリスク | 必須ゲート |
|---|---|---|---|
| 0 | 基準値と同値性チェックの固定 | 基準不足 | 全ゲートの現行結果を記録 |
| 1 | 未参照 private、no-op、完了文書、未使用 font | 動的参照 | import/public API/full tests/package |
| 2 | scalar axis、coordinate/selector/depth helper | 数値・端点差 | unit/reference SVG/full tests |
| 3 | CLI 共通引数 | default/help 差 | Namespace/help/CLI smoke |
| 4 | Web 小 utility | import/CSP | Node/Playwright/offline bundle |
| 5 | LOSAT/worker runtime | thread/fallback 差 | 全 worker mode/Cloudflare/offline |
| 6 | UI 反復、Python/Web adapter | markup/session 差 | Playwright/session/gallery |
| 7 | test fixture/parameterize | coverage 低下 | case matrix/coverage/full tests |
| 8 | duplicate file、生成物 | path/package/docs 破損 | build/docs/gallery/release artifact |

各 PR は 1 テーマ、レビュー可能な差分量、単独 revert 可能な状態にする。複数 Phase を 1 PR に混在させない。

## 7. 共通の検証ゲート

### 7.1 静的・Python

```bash
ruff check gbdraw/ --select=E,F,W --ignore=E501,W503
pytest tests/ -v -m "not slow"
pytest tests/ -v
python -m build
git diff --check
```

### 7.2 Web

```bash
node --test tests/web/*.test.mjs
node -e "console.log(require.resolve('@playwright/test'))"
playwright test
```

Node の Playwright がない場合は、AGENTS.md の指示に従い Python Playwright で対象ブラウザ確認を実行する。Chromium sandbox error は Playwright 不在として扱わず、同一チェックを必要な sandbox escalation 付きで再実行する。

### 7.3 Packaging・offline

```bash
python tools/prepare_browser_wheel.py
python tools/verify_gui_offline.py
python tools/prepare_cloudflare_pages.py
```

対象 PR に応じて gallery/session の refresh/verification tool も実行する。生成物は検証に使用するが、`dist/`、`gbdraw.egg-info/`、gitignored browser wheel は commit しない。

### 7.4 同値性

- public API snapshot が一致する。
- CLI help、Namespace、終了コード、stderr が一致する。
- reference SVG が一致する。意図しない差分は 1 byte でも不合格とする。
- interactive SVG の表示、選択、drag、download、session round-trip が一致する。
- old/new session と gallery example が読み込める。
- wheel/sdist/offline/Cloudflare artifact 内の必須 file 一覧が一致する。
- 削除対象への参照が `rg`、import graph、package manifest に残っていない。

### 7.5 削減量

PR ごとに次を記録する。

- Production Python の追加・削除行数
- Web source の追加・削除行数
- Test/Tools/Docs の追加・削除行数
- 追跡対象ファイルの増減バイト数
- helper/module 数の増減

原則として `git diff --stat` が純減にならない PR は、この計画の削減成果として merge しない。

## 8. 期待値と打ち切り基準

機械的に検出された削減余地の総量は約 3,904 行である。

- 未参照 private: 762 行
- Production Python の同一本体: 432 行
- Web JS の同一本体: 702 行
- Python test の同一本体: 624 行
- Tools の同一本体: 9 行
- 完了済み計画文書: 1,375 行

これは共通 helper の追加行や互換 adapter を含まない gross ceiling である。現実的な初回目標は次のとおりとする。

| 指標 | 保守的目標 | 条件付き上限 |
|---|---:|---:|
| テキスト純減 | 2,800～3,400 行 | 各共通化が純減する範囲 |
| 追跡対象容量 | 約 1.65 MB | duplicate file を安全に統合できれば約 52.4 MB |
| 生成 palette 容量 | 0 MB | 再生成・URL 維持が成立すれば約 399 MB |

次の場合は、その候補の改修を打ち切る。

- 互換 branch や環境判定が増え、純減しない。
- 出力、help、session、artifact に差が出る。
- 共通 helper の引数や flag が増え、元実装より理解しにくい。
- 動的参照の不存在を証明できない。
- テストを弱めないと統合できない。

削減量より同値性を優先し、条件を満たさない候補は現状維持を正解とする。

## 9. 完了定義

本計画は、以下を満たした時点で完了とする。

1. Phase 1～7 の安全な候補が完了または根拠付きで保留されている。
2. 各 merge 済み PR が純減し、共通検証ゲートを通過している。
3. 公開 API、CLI、図の出力、Web UI、session、offline/package の機能差がない。
4. 実測した削減行数・容量と、保留項目の理由が最終記録に残っている。
5. 新しい循環依存、汎用すぎる helper、巨大な option schema、未使用 compatibility layer を導入していない。

## 10. 実施結果（2026-07-13）

安全性と純減条件を満たす候補は実装し、条件を満たさない候補は根拠を確認して保留した。

| Phase | 結果 | 主な内容 |
|---|---|---|
| 0 | 完了 | 公開 API、Circular/Linear CLI action・既定値・代表 `Namespace`・help の SHA-256 契約を `tests/test_public_contract.py` と `tests/fixtures/public_contract.json` に固定 |
| 1 | 完了 | 未参照 private 定義、no-op、完了文書、未使用 Web font を削除 |
| 2 | 完了 | scalar axis、records/depth/coordinate/selector helper を既存 owner に統合し、`features.colors`/`features.visibility` の循環依存を解消 |
| 3 | 完了 | Circular parser の色・解析引数を既存 `add_color_args` / `add_analysis_args` に統合。契約 snapshot で順序、文言、型、既定値を固定 |
| 4 | 完了 | binary codec、文字列正規化、transform、title、JSON clone、DOM class、text download、SVG serialization/defs の安全な重複を統合 |
| 5 | 完了 | LOSAT direct/WASI runtime を共有。main/worker の dispatch、RPC、cancel 境界と thread budget は維持 |
| 6 | 完了 | `auto-value-field` を component 化し、Pyodide の orthogroup serializer を `web_support` の共通実装へ統合 |
| 7 | 完了 | 数値 track capture、feature factory、同一 label test body を共有。入力・期待値・元の test 名は維持 |
| 8-A | 完了 | canonical example を参照するよう test を変更し、重複 test input と未使用 Web font を削除 |
| 8-B | 保留 | palette SVG は公開 docs/example URL から直接参照される。URL を維持する release 時自動生成工程がないため追跡を継続 |

### 10.1 実測削減量

- 基準 commit `cbfced9` から主要削減 commit `fea52cd` まで、追跡対象は 1,347,675,940 bytes から 1,297,545,546 bytes へ 50,130,394 bytes 純減した。
- 同 commit は 1,355 行追加、698,118 行削除で、重複 fixture を含め 696,763 行純減した。
- その後の契約固定と追加共通化は、本記録を除いて 588 行追加、712 行削除、差し引き 124 行純減した。契約 fixture/test の追加を含む。

### 10.2 意図的に統合しなかった項目

- threaded/non-threaded worker RPC は message shape と lifecycle が異なり、共通化すると環境分岐が増えるため現状を維持した。
- modal shell と巨大関数は、移動だけでは純減せず責務境界も改善しない候補を保留した。
- `analysis/collinearity.py` と `analysis/protein_colinearity.py` の短い row 数値 helper、および closure 固有の test double は、import/factory の追加を含めると純減または単純化にならないため維持した。
- `tests/test_web_packaging.py` の CSP、asset、埋め込み source の assert は artifact 契約そのものを検証するため維持した。Node/Playwright への置換だけを目的とした test infrastructure は追加していない。

### 10.3 性能・並列実行の確認

- LOSAT の `DEFAULT_MAX_WORKERS=4`、total thread budget 16、threads/job 上限 16、auto allocation、threaded 判定、module worker 数は変更していない。
- 共通化したのは 1 pair の標準 WASI setup/start/decode のみで、main/worker は従来どおり別 realm で実行する。main 側の abort check と worker 側の job lifecycle も維持した。
- Pyodide serializer は旧埋め込み実装と出力一致を確認した。合成 10,000 member の比較では旧 24.8 ms、共通実装 30.9 ms（約 6.1 ms 増）であり、LOSAT 計算後の serialization として許容できる範囲だった。

### 10.4 検証記録

- 変更対象の契約・描画・packaging test: 310 passed, 19 skipped。削減前 HEAD と現行の公開 API/CLI 契約も同一 snapshot に完全一致した。
- 全 test: fast 1114 passed / 24 skipped、slow 4 passed / 2 skipped。2 分割の合計 1118 passed / 26 skipped で収集した 1144 件をすべてカバーした。
- `prepare_browser_wheel.py`、offline `check-assets`、最終 wheel の `inspect-wheel`、isolated `python -m build`、`prepare_cloudflare_pages.py` は成功した。
- Node、Node Playwright、Python Playwright、Ruff は実行環境に存在しない。AGENTS.md 記載の全経路を確認し、Web は Python packaging/埋め込み実行 test で代替した。
