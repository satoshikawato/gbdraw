# Test Failure Handoff 2026-08-09

## 現在の基準

最終 shared tree で次を実行した。

```bash
env TMPDIR=/tmp \
  PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  /home/kawato/micromamba/bin/python -m pytest \
  tests/ -m "not slow" --timeout=1800 -q
```

結果は 2,925 passed、17 skipped、6 deselected、5 failed、実行時間18分52秒だった。長時間テストは [AGENTS.md](../../AGENTS.md) の規則どおり、30分以上の timeout で実行する。

この文書の5件は Ponytail audit の実装で発生した回帰ではない。うち3件は HEAD archive または HEAD と同一の blob で再現し、残る2件は既存 worktree の stale snapshot/test fixture である。次のセッションでは、下記の owner ごとに直す。

## 1. WSSV の比較配列を session import で復元する

### 失敗

```text
tests/test_web_packaging.py::test_web_javascript_behavior[session-definition-rehydration]
tests/web/session-definition-rehydration.test.mjs:52
expected homology-comparison sources: 20
actual: 0
```

```bash
env TMPDIR=/tmp \
  PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  node tests/web/session-definition-rehydration.test.mjs
```

同じ失敗は HEAD archive でも再現する。`tests/test_wssv_gallery_fastas.py` は 3/3 通過し、file-backed builder 単体では20件を復元できる。

### 原因

[config.js](../../gbdraw/web/js/services/config.js) の `importSession()` は、current schema の `editorState.featureCatalog.schema === 3` を見ただけで `buildRestoredMatchSequenceSources()` を呼ばない。WSSV catalog の `sequenceSources` は circular reference だけで、20個の `homology-comparison` source を持たない。その後 `matchSequenceRegistry` を catalog source だけで reset するため0件になる。

確認箇所:

- `gbdraw/web/js/services/config.js` の `hasAuthoritativeCatalogSequenceSources` と `restoredFileSequenceSources`
- `gbdraw/web/js/services/config.js` の `catalogSequenceSources` と `matchSequenceRegistry.reset()`
- `gbdraw/web/js/app/match-sequences.js::buildRestoredMatchSequenceSources`
- `gbdraw/web/js/services/feature-catalog.js::featureStateFromCatalog`

### 修正条件

- schema-3 catalog を feature identity の authority として維持する。
- catalog にない file-backed `homology-comparison` source は復元する。
- catalog source と file-derived source は既存 registry の同一性規則で重複排除する。
- WSSV は source index 0〜19、各配列非空になる。
- HmmtDNA の通常 session import を壊さない。
- 削除済み match payload v1 や legacy converter は戻さない。

最低限の確認:

```bash
node tests/web/session-definition-rehydration.test.mjs
python -m pytest tests/test_wssv_gallery_fastas.py -q --timeout=1800
node tests/web/feature-catalog.test.mjs
```

## 2. orthogroup sidecar test fixture に record-scoped identity を持たせる

### 失敗

```text
tests/test_session_io.py::test_session_sidecar_saves_complete_orthogroup_state
Orthogroup 'og_1' member must identify one biological feature by a
record-scoped stable feature ID or source feature index and a valid record index/key.
```

### 原因と修正箇所

[test_session_io.py](../../tests/test_session_io.py) の fixture は member に `featureSvgId` と `proteinId` だけを与えている。catalog schema 3 は rendered SVG ID を biological join key として受け付けない。

production validator を緩めず、fixture の2 member を次の形に直す。

```python
{"recordIndex": 0, "stableFeatureSvgId": "feature-1", "proteinId": "protein-1"}
{"recordIndex": 1, "stableFeatureSvgId": "feature-2", "proteinId": "protein-2"}
```

`biological_feature_metadata` 側の `record_idx` と `stable_feature_id` に一致させる。必要なら canonical `recordKey` を使ってもよいが、index/key の競合を作らない。

```bash
python -m pytest \
  tests/test_session_io.py::test_session_sidecar_saves_complete_orthogroup_state \
  --timeout=1800 -q
```

## 3. Linear CLI contract snapshot を意図と照合する

### 失敗

```text
tests/test_public_contract.py::test_public_api_and_cli_contract_snapshot
```

現在の Linear CLI hash:

```json
{
  "actions_sha256": "2e6de596a6e45d3cbde9c359306d0fc49c10f0630185677c7ba6ee80acee85d5",
  "defaults_sha256": "436805006d27849e1a291ab7c948261e25e56b9efd5e381b2fbc4aa9e87a2180",
  "representative_sha256": "51722d3740c6278c268111cf08cafbc75dcb7bc1df1378ab6088a7c1a41665c2"
}
```

期待値は [public_contract.json](../../tests/fixtures/public_contract.json)、収集処理は [test_public_contract.py](../../tests/test_public_contract.py)、parser owner は `gbdraw/linear.py::_get_args()` である。public Python API と Circular CLI の snapshot は一致している。

hash だけを即更新しない。まず action、default、representative namespace のどれが変わったかを出力し、現在の Linear CLI 変更が意図した contract か確認する。

```bash
python -c 'import json; from tests.test_public_contract import build_contract; print(json.dumps(build_contract()["linear_cli"], indent=2))'
python -m pytest tests/test_public_contract.py -vv --timeout=1800
```

変更が意図どおりなら `tests/fixtures/public_contract.json` の Linear CLI 3 hash だけを更新する。意図していない場合は `gbdraw/linear.py` または CLI helper の option/default を直す。API hash と Circular CLI hash は変更しない。

## 4. tutorial manifest の旧 direct-link assertion を削除する

### 失敗

```text
tests/test_cli_comparison_how_to_recipe_contracts.py:85
tests/test_cli_how_to_recipe_contracts.py:278
```

両テストと対象 chapter は HEAD と同じ blob である。commit `a68e24ef` で読者導線を `GETTING_TUTORIAL_DATA.md` に移した後、次の文字列 assertion だけが残った。

```python
assert "../../../gbdraw/web/tutorial-data/manifest.json" in source
```

修正は上の assertion 2行の削除でよい。chapter に内部 manifest の direct link を戻さない。現在の reader-facing owner は [GETTING_TUTORIAL_DATA.md](../GETTING_TUTORIAL_DATA.md) である。

```bash
python -m pytest \
  tests/test_cli_how_to_recipe_contracts.py \
  tests/test_cli_comparison_how_to_recipe_contracts.py \
  --timeout=1800 -q
```

修正前の結果は 17 passed、2 failed である。

## 5. final pytest 外で確認した session preflight failure

これは上の5 failed には含まれないが、対象 Playwright を直接実行すると残る。

```text
tests/web/depth-track-session.playwright.spec.js
Session preflight rejects invalid canonical data without resetting live state
expected legend: bottom
actual legend: right
```

HEAD archive でも同じ結果になる。fixture は canonical request で `output.legend=bottom`、legacy UI で `linearLegendPosition=right` を同時に持つ。現在は [layout-preferences.js](../../gbdraw/web/js/app/layout-preferences.js) の legacy migration が `right` を採用し、test は active canonical output の `bottom` を期待する。

修正前に authority を決める。current-schema canonical request を active render state の authority とするなら、dormant layout preference を保持しつつ active linear legend だけ canonical 値にする。legacy UI を優先する仕様なら、test expectation と session compatibility 文書を合わせる。単に `right` を `bottom` へ固定しない。

確認箇所:

- `gbdraw/web/js/services/config.js` の `restoreLayoutPreferences(ui, { preserveActive: ... })`
- `gbdraw/web/js/app/layout-preferences.js::migrateLegacyLayoutPreferences`
- `tests/web/depth-track-session.playwright.spec.js` の repaired canonical session と v40 layout-preference round trip

```bash
env TMPDIR=/tmp \
  PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  npx playwright test tests/web/depth-track-session.playwright.spec.js \
  --project=chromium \
  --grep "Session preflight rejects invalid canonical data" \
  --workers=1 --timeout=1800000
```

## 今回すでに解消したもの

次の項目は最終 gate で通過したため、上記5件の修正時に戻さない。

- circular preset oracle と linear label indexed helper の stale 引数
- cancellation rollback の raw editor-state 復元
- orthogroup membership scope と collinearity presentation scope の分離
- embedded canonical Web request の comparison metadata
- browser-rendered definition gap 3件。sandbox外、`--timeout=1800` では 3/3 通過
- output comparison 16件、example 再現 14件、CLI/Python の全 recipe

修正後の最終確認:

```bash
ruff check gbdraw/
git diff --check
env TMPDIR=/tmp \
  PATH=/home/kawato/micromamba/bin:/usr/local/bin:/usr/bin:/bin \
  python -m pytest tests/ -m "not slow" --timeout=1800 -q
```
