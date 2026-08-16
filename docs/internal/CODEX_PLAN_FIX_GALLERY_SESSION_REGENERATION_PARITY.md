# Gallery session Load → Generate parity 修正計画

Status: diagnosis complete; implementation not started
Date: 2026-08-15

## 1. 結論

通常の保存セッションに「レシピ、`renderRequest`、Load 後 UI state、Generate 後 SVG の全てが常に同値」という四層契約は追加しない。

現行セッション設計では、次の分離は意図的である。

- `renderRequest` / `results`: 直近に Generate した committed artifact
- `config` / `ui`: その後の未 Generate 変更も含む active draft

`tests/web/session-draft-authority.test.mjs` もこの不一致を正式な動作として固定している。したがって、両者を常に同一化すると、「設定を変えたが、まだ Generate していない」状態を保存できなくなる。

必要なのは次の一つの追加契約だけである。

> Gallery として公開するセッションは、公開時点に限り clean snapshot であり、active rendering draft は committed `renderRequest` から機械的に再構成できること。

これを Gallery 更新パイプラインの単一 invariant として実装する。四層それぞれに Majanivirus/Hepatoplasmataceae 専用の `1.0` 断言を重複配置しない。

## 2. 予定する authority

| 対象 | authority | 契約 |
| --- | --- | --- |
| 通常の Web session | `renderRequest` と active draft の二つ | 未 Generate の draft 差分を許容する |
| 公開済み Gallery session | canonical `renderRequest` | 公開処理が rendering draft をここから再構成する |
| `editorState` | canonical live-edit overrides | 新しい SVG に既存の一元化済み replay owner が適用する |
| Gallery の表示用 metadata | `GallerySessionExample` と明示 allowlist | rendering 意図と混ぜない |

Majanivirus と Hepatoplasmataceae の Shaft Width Ratio には新しい「レシピ owner」を追加しない。これらは完全な CLI recipe を現時点で持たず、セッション自体の canonical `renderRequest` が再生入力である。対象 3 例の canonical override を意図的な `1.0` に更新し、Gallery 公開 validator と既存の汎用 request projection 契約で後続層への反映を保証する。

## 3. 再現結果

ローカルの現行 Web app で対象セッションを Load し、Generate Diagram を実行して次を確認した。いずれも Generate 自体は成功し、エラー表示は出ない。問題は入力意図が静かにすり替わることである。

### 3.1 `HmmtDNA_ATskew`

- 保存 SVG の `viewBox` は約 `1280.93 × 953.70`、再生後は約 `1315.41 × 958.85` に変わる。
- 保存 request には AT skew slot の `positive_color=#deaf6e` / `negative_color=#7294e3` がある。
- Worker から返った生 SVG にもこの 2 色は残っている。
- mounted SVG に取り込んだ後、palette watcher が AT skew を `#f6a8f1` / `#bd96df` に上書きする。凡例には元の色が残るため、軌道と凡例が不一致になる。
- 保存 request は `objects.definition.circular.interval=30` を持つが、Web からの再生 request はこれを持たず、既定値 `20` に戻る。
- 保存 canonical request を CLI で replay すると、保存 SVG の layout と slot color を再現できる。

### 3.2 `majanivirus_orthogroup`

- 保存 request の LOSATP filter は bitscore `100`、e-value `1e-3`、identity `20`。
- 保存 `config.adv` で bitscore と identity が文字列の `"100"` / `"20"` になっている。
- Load 後 UI state は bitscore `50`、e-value `1e-3`、identity `70` になる。
- 保存 request/resources には custom default-color table と specific-color rules があるが、active `config.colors` / `config.rules` は空で、Load 後は default palette に戻る。
- 保存図は 152 similarity groups / 627 comparison paths を持つ。再生後は変更された filter で再計算され、193 groups / 262 comparison paths に変わる。
- WSSV-like proteins、BIRP、tyrosine recombinase などの custom 色と凡例が失われ、CDS 中心の default palette になる。
- 保存 request の Shaft Width Ratio は `0.5`、Load 後 UI は現行既定値 `1.0`である。

### 3.3 Hepatoplasmataceae 2 例

対象:

- `hepatoplasmataceae_orthogroup.gbdraw-session.json.gz`
- `hepatoplasmataceae_collinear.gbdraw-session.json.gz`

共通する再現結果:

- 保存 request の filter は bitscore `50`、e-value `0.01`、identity `0`。
- active draft はこれらを省略しており、Load 後は bitscore `50`、e-value `1e-5`、identity `70` に戻る。
- 保存 request は 4 本の precomputed comparison と、collinear で 554 groups、orthogroup で 577 groups の結果を持つ。
- active `linearComparisonPlan` は `adjacent` / `losat` / empty edges であり、保存済み comparison intent を表していない。
- 両セッションとも保存 request の Shaft Width Ratio は `0.5`、Load 後 UI は `1.0`である。

## 4. 根因

### Root cause A: Gallery refresh が committed artifact と古い active draft を結合している

`tools/refresh_gallery_sessions.py::_merge_refreshed_gallery_artifacts()` は、再生した `renderRequest`、`results`、`resources` 等を差し替える一方、`config` は `REFRESHED_GALLERY_ARTIFACT_KEYS` に含めず旧セッションから残す。

その後に次の個別同期を行うが、active rendering draft 全体の契約にはなっていない。

- `_sync_legacy_legend_control_with_render_request()`
- `_sync_circular_track_draft_with_render_request()`
- `_restore_rendered_palette_file_binding()`

`_validate_staged_gallery_session()` も schema、resource、geometry、cache、SVG は検査するが、active draft が次の Generate で committed intent を再現するかは検査しない。

### Root cause B: current-writer config の不正な scalar type を静かに既定値へ戻す

`gbdraw/web/js/services/config.js::safeDeepMerge()` は、target と source の primitive type が同じ場合だけ代入する。そのため、number state に保存された string `"100"` / `"20"` を merge すると無視される。

`validateCurrentWriterActiveConfig()` は domain や field name、collection shape は検査するが、これら scalar の型を検査しない。現行 writer の壊れたデータが validation を通過し、merge で静かに捨てられる。

### Root cause C: current-session projection が空の color draft を canonical resource より優先する

Load 時の `projectCanonicalSessionRequest(..., deferResourceContent: true)` は、通常セッションの未 Generate 変更を守るために必要な仕様である。しかし Gallery の active draft が `colors={}` / `rules=[]` のままなので、canonical resources 中の color table を hydrate せず、空の draft を次の Generate に使う。

`deferResourceContent` 自体を廃止するのではなく、Gallery 公開時にだけ canonical resources を hydrate した clean draft を作る必要がある。

### Root cause D: skew slot の色 precedence が post-mount styling で重複実装されている

`gbdraw/web/js/app/run-analysis.js` は Generate 後に `appliedPaletteColors` を更新し、`gbdraw/web/js/app/svg-styles.js::applyPaletteToSvg()` の watcher を起動する。

`applyPaletteToSvg()` は `data-gbdraw-slot-renderer="dinucleotide_skew"` の各 group に global `skew_high` / `skew_low` を上書きし、slot 固有の `positive_color` / `negative_color` を考慮しない。

その一方、`gbdraw/web/js/app/track-slot-colors.js::resolveTrackSlotSkewColorValue()` はすでに「slot 固有色 > palette 継承色」の precedence を所有している。問題は生成側ではなく、mounted SVG の色更新がこの owner を使っていないことである。

### Root cause E: circular definition interval の派生規則が CLI にだけある

`gbdraw/circular.py` は `--definition_font_size 28` に対して `objects.definition.circular.interval=int(28 + 2)=30` を派生させる。Web の canonical request builder は font size を送るが interval は送らないため、Python 側既定値 `20` に戻る。

legacy flat override に対する同じ派生処理は `gbdraw/session_request_codec.py::_migrate_flat_config_overrides()` にもあるが、current canonical path には適用されない。派生値の owner が cross-surface で一元化されていない。

### Root cause F: Shaft Width Ratio は「新しい既定値の regression」ではなく、古い explicit value の化石化である

現行既定値は Python config、Web state、feature-rendering utility のいずれも `1.0` である。既定値は commit `395670cb` で `0.5` から `1.0` に変更された。

しかし対象 3 セッションの canonical request には explicit `0.5` が残っている。Gallery refresh は古い canonical request を `--session` で replay するため、新しい既定値を適用せず `0.5` を再生し続けた。

これは全セッションの `0.5` を自動変換すべき問題ではない。`0.5` を意図的に選ぶ図を壊すため、対象 3 例のみを `1.0` へ明示更新する。

## 5. 実装計画

### Phase 1: 通常 session と Gallery clean snapshot の契約をテストで分離する

1. `tests/web/session-draft-authority.test.mjs` の既存契約を保持する。
   - active draft と committed artifact の不一致を禁止しない。
   - Load は current-writer draft を正しく優先する。
2. Gallery 公開用の純粋な normalization/validation API を `gbdraw/web/js/services/gallery-session-migration.js` に追加する。
   - 例: `normalizeGallerySessionForPublication(session)`
   - canonical `renderRequest` と `resources` から rendering draft を再構成する。
   - `projectCanonicalSessionRequest()` を resource hydration 有効で再利用し、Python に二重の projection を実装しない。
   - `modeProfiles` と `linearComparisonPlan` も、正規化済み active values と canonical comparison descriptors から作る。
   - selected tab、drawer、tutorial metadata などの non-rendering state は明示 allowlist だけを保存する。
   - `editorState` はこの config projection に混ぜず、現行の canonical replay owner に委ねる。
3. `tools/promote_gallery_session.mjs` からこの公開用 owner を呼び出せるようにする。
4. 「セッション全体の JSON 一致」ではなく、render-affecting field の normalized semantic projection だけを比較する。

### Phase 2: Gallery refresh を clean-snapshot publisher に変える

1. `tools/refresh_gallery_sessions.py` で fresh `renderRequest` / `results` / `resources` を merge した後、Phase 1 の公開用 normalization を一度だけ実行する。
2. rendering domain に対する次の ad hoc repair を、汎用 clean projection がカバーした後に削除する。
   - `_sync_legacy_legend_control_with_render_request()`
   - `_sync_circular_track_draft_with_render_request()`
   - `_restore_rendered_palette_file_binding()`
3. `_validate_staged_gallery_session()` から同じ semantic validator を呼び、次を公開前に fail-fast させる。
   - number であるべき filter が string
   - canonical resource はあるが active colors/rules が空
   - comparison descriptors と active plan/mode が不一致
   - arrow geometry やその他の render-affecting override が不一致
4. この validator は `gbdraw/web/gallery/sessions/` の公開物にだけ適用する。ユーザー保存セッションには適用しない。

### Phase 3: current-writer scalar validation を fail-loud にする

1. `gbdraw/web/js/services/config.js::validateCurrentWriterActiveConfig()` に、`form` / `adv` の既知 scalar field の型契約を追加する。
2. current session version が number として書くべき値を string で持つ場合は、`safeDeepMerge()` に曖昧な全般 coercion を追加せず、明示エラーにする。
3. 過去 version の string scalar を受ける必要がある場合は、その version migration 境界でフィールド別に正規化する。current-writer import と混ぜない。
4. テストで、Majanivirus 型の `"100"` / `"20"` が静かに default へ落ちないことを固定する。

### Phase 4: HmmtDNA の 2 つの cross-surface bug を一元化する

#### 4A. skew color precedence

1. `gbdraw/web/js/app/svg-styles.js::applyPaletteToSvg()` が skew group を更新するとき、対応する track slot を解決する。
2. 色の決定は既存の `resolveTrackSlotSkewColorValue()` だけに委ねる。
3. GC skew の palette 継承は維持し、AT skew の explicit 2 色は palette watcher で上書きしない。
4. raw Worker SVG と mounted SVG で slot の色意味が一致するユニットテストを追加する。

#### 4B. circular definition interval

1. `font_size` から `interval=int(font_size + 2)` を派生する規則を、CLI アダプタではなく Python の typed config/request normalization 境界へ移す。
2. canonical config override で circular definition font size があり、interval が省略された場合だけ派生する。
3. explicit interval がある場合は常にそれを優先する。
4. `gbdraw/circular.py` と `gbdraw/session_request_codec.py` に重複している派生処理は、共通の純粋関数または config resolver を通す。
5. CLI、Python canonical replay、Web Worker request が同じ effective interval `30` を得るテストを追加する。

### Phase 5: Shaft Width Ratio を対象 3 例で `1.0` に更新する

対象:

- `majanivirus_orthogroup`
- `hepatoplasmataceae_orthogroup`
- `hepatoplasmataceae_collinear`

手順:

1. 各セッションの canonical request を typed codec で decode する。
2. `objects.features.arrow_geometry.shaft_width_ratio` の explicit override を `1.0` にする。raw JSON path の手編集はしない。
3. 変更後 request を入力に Gallery refresh を実行し、session、SVG、thumbnail を再生する。
4. `tests/test_gallery_session_semantics.py` に 3 ID を一つの parameterized assertion として追加し、canonical value が `1.0` であることを固定する。
5. Load state と Generate request への投影は、Phase 1–2 の汎用 clean-snapshot validator と既存/追加の request projection test に任せる。3 例×複数層の literal `1.0` テストは作らない。

### Phase 6: Gallery assets を generator 経由で更新する

1. `tools/refresh_gallery_sessions.py` で対象 session を再生する。
2. `tools/prepare_interactive_gallery_assets.py` で次を同期する。
   - `gbdraw/web/gallery/sessions/`
   - `gbdraw/web/gallery/sources/`
   - `gbdraw/web/gallery/examples/`
   - `gbdraw/web/gallery/thumbnails/`
   - `gbdraw/web/gallery/examples.json`
3. generated session/SVG/gzip を手編集しない。
4. HmmtDNA は現在の committed artifact が意図した図なので、再生後も視覚的に同一であることを確認する。
5. Shaft Width Ratio 変更で図が変わる Majanivirus/Hepatoplasmataceae は、次の session-based tutorial media を capture script で再比較する。
   - `majanivirus_orthogroup/manual-07-01-orthogroup-preview.webp`
   - `majanivirus_orthogroup/manual-08-01-orthogroup-popup.webp`
   - `hepatoplasmataceae_orthogroup/manual-07-01-orthogroup-overview.webp`
   - `hepatoplasmataceae_orthogroup/manual-08-01-orthogroup-popup.webp`
   - `hepatoplasmataceae_collinear/manual-06-01-collinear-overview.webp`
   - `hepatoplasmataceae_collinear/manual-07-01-collinear-block-popup.webp`
6. HmmtDNA の `manual-09-01-atskew-preview.webp`、`manual-10-01-feature-popup.webp`、`post-02-01-legend-editor.webp` は capture comparison を行い、意図した図が変わらない限り不要に置き換えない。
7. 実際に置き換えた media だけ `docs/internal/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md` の hash/capture 記録を更新する。

## 6. Regression strategy

「四層を全例で全比較」ではなく、次の 3 段構成にする。

### 6.1 単体契約

- canonical request → active config projection
- resource table → colors/rules/whitelist hydration
- current-writer scalar type validation
- explicit skew slot color → mounted SVG color precedence
- circular definition font size → effective interval
- active config → next canonical request の round trip

これらは例固有値ではなく、各 owner の汎用契約としてテストする。

### 6.2 Gallery 公開 gate

全 Gallery session に対し、公開用 normalization の前後で render-affecting semantic signature が安定することを検査する。対象は filter、palette/rules、track slots、geometry、comparison plan/mode、label/legend 設定などであり、timestamp、selected tab、drawer state は除外する。

### 6.3 代表 browser E2E

`tests/web/gallery-session-regeneration.playwright.spec.js` に次の小さな matrix を追加する。

| 例 | この例でしか見えない境界 |
| --- | --- |
| `HmmtDNA_ATskew` | per-slot skew color precedence、derived circular interval |
| `majanivirus_orthogroup` | custom resource hydration、non-default filters、dense cached LOSATP result |
| `hepatoplasmataceae_orthogroup` | orthogroup mode、precomputed comparison reuse |
| `hepatoplasmataceae_collinear` | collinear mode、precomputed comparison reuse |

E2E は SVG の完全文字列 hash ではなく、次の semantic signature を比較する。

- normalized `viewBox` / geometry metadata
- record、track、comparison/group 数
- 必須 label/legend key
- 代表色と explicit slot colors
- LOSATP mode/filter と reusable result identity
- resolved arrow shaft ratio

ID、timestamp、serialization 順のみの差で失敗する brittle な exact-SVG hash は追加しない。

## 7. Acceptance criteria

### Architecture

- 通常 session は dirty draft を保存・Load できる。
- Gallery 公開処理だけが clean snapshot を強制する。
- rendering draft の公開用 projection owner は 1 つである。
- color resource parsing、skew precedence、derived interval に重複 owner を追加しない。
- Majanivirus/Hepatoplasmataceae 専用の新しい recipe/config 層を作らない。

### `HmmtDNA_ATskew`

- Load → Generate 後も AT skew が `#deaf6e` / `#7294e3` である。
- raw Worker SVG と mounted/result SVG の色が一致する。
- definition font size `28` から effective interval `30` が得られる。
- committed 図の track、label、legend、layout と再生図が semantic parity を保つ。

### `majanivirus_orthogroup`

- Load 後が bitscore `100`、e-value `1e-3`、identity `20`、orthogroup mode である。
- custom default colors と specific rules が active draft に hydrate される。
- Generate 後も 152 groups / 627 comparison paths と必須 legend/custom colors を保つ。
- canonical と resolved Shaft Width Ratio が `1.0` である。

### Hepatoplasmataceae

- 両例とも Load 後が bitscore `50`、e-value `0.01`、identity `0` である。
- orthogroup は orthogroup mode / 577 groups、collinear は collinear mode / 554 groups を再生後も保つ。
- precomputed comparison の利用意図と凡例/色が維持される。
- 両例とも canonical と resolved Shaft Width Ratio が `1.0` である。

### Existing live-edit behavior

- feature 色/ラベル/表示性などの即時反映が Generate 無しで動く。
- Save → Load で canonical editor overrides が残る。
- Generate/reflow 後にも editor overrides が新 SVG へ replay される。
- これらの既存契約に Gallery clean-snapshot 処理を混ぜない。

## 8. 検証コマンド

実装時に実行する最小 gate:

```bash
node tests/web/session-draft-authority.test.mjs
node tests/web/session-request.test.mjs
node tests/web/gallery-session-migration.test.mjs
node tests/web/svg-styles-track-groups.test.mjs
pytest tests/test_session_request_codec.py -v
pytest tests/test_circular_feature_width.py -v
pytest tests/test_refresh_gallery_sessions.py -v
pytest tests/test_gallery_session_semantics.py -v
```

browser gate:

```bash
npx playwright test tests/web/gallery-session-regeneration.playwright.spec.js
```

asset/tutorial gate:

```bash
python tools/prepare_interactive_gallery_assets.py
python tools/capture_gallery_tutorial_screenshots.py --example HmmtDNA_ATskew --check
python tools/capture_gallery_tutorial_screenshots.py --example majanivirus_orthogroup --check
python tools/capture_gallery_tutorial_screenshots.py --example hepatoplasmataceae_orthogroup --check
python tools/capture_gallery_tutorial_screenshots.py --example hepatoplasmataceae_collinear --check
```

`prepare_interactive_gallery_assets.py` は生成ツールであり `--check` を持たない。実行後に generated-asset diff を個別にレビューする。最後にリスクに応じて次も実行する。

```bash
pytest tests/ -v -m "not slow"
ruff check gbdraw/
git diff --check
```

## 9. 対象外

- 通常 session を常に clean にする変更
- 全保存セッションの Shaft Width Ratio `0.5` を `1.0` に自動 migration すること
- CLI-only field を無制限に Web config へ passthrough する新層
- Gallery 全例の SVG exact hash 固定
- 本件と無関係な Gallery のデザイン変更

## 10. 実装時のチェックリスト

- [ ] 通常 session の dirty-draft 契約を維持する
- [ ] Gallery publication normalization を 1 owner で実装する
- [ ] refresh の古い `config` 保存を廃止する
- [ ] ad hoc rendering-field sync を削除する
- [ ] current-writer scalar type mismatch を fail-loud にする
- [ ] resource-backed colors/rules を clean draft へ hydrate する
- [ ] comparison plan/mode を canonical intent から再構成する
- [ ] skew slot color precedence を既存 resolver へ一元化する
- [ ] circular definition interval の派生規則を Python boundary で一元化する
- [ ] Majanivirus と Hepatoplasmataceae 2 例の Shaft Width Ratio を `1.0` にする
- [ ] Gallery session/SVG/thumbnail/catalog を generator 経由で更新する
- [ ] 必要な tutorial media だけを capture script で更新する
- [ ] unit、Gallery publication、browser E2E の 3 段 gate を通す
- [ ] generated asset と production/test diff を分けてレビューする
