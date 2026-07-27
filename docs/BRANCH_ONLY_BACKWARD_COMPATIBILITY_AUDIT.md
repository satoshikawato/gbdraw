# Branch-only backward compatibility audit

## 1. 目的と結論

この文書は、`2026-07-20` ブランチを `main` へマージする前に、`main` に一度も存在しなかった中間実装のためだけに残っている後方互換コードを削除するための監査結果である。

監査時点:

- 監査日: 2026-07-27
- `main`: `b05a6bb89e65b6a0b211a1ee10c1a3a402c21165`
- `2026-07-20`: `16762f3d1b997bae0ff8e81ba2af1c03e9ca5f26`
- merge base: `b22b0eb7265f1841bf701e6b161e28d6f8b79ea3`
- 監査対象差分: `b22b0eb..16762f3`

結論は次のとおり。

1. **削除対象は session v34 の受理と、session v35 の LOSAT protein identity/cache 移行層である。**
2. v34 と v35 はこのブランチ内で導入・置換され、`main` の全履歴には存在しない。
3. `main` に実在する session v27–33、canonical request schema 1/2、protein raw schema 2、protein derived schema 1 の移行は維持する。
4. feature underlay、repeat feature shape、gallery request schema 2→3、旧 feature ID の移行も `main` の実在セッションを再現するために必要であり、今回の削除対象ではない。
5. 現行 session version は **36 のまま維持**し、reader の対応集合を **27–33 と 36** にするのが最小かつ安全である。34 に詰め直すと、開発ブランチで既に作られた別仕様の v34 と番号が衝突する。

本セッションでは監査のみを行い、実装コード、テスト、fixture、既存文書は変更していない。

### 作業ツリーに関する注意

監査時の作業ツリーには、多数の tracked file に行末差分だけの変更が見えていた。`git diff --ignore-space-at-eol` では、監査文書を作る前の substantive な tracked diff は確認されなかった。次セッションでは無関係な変更を reset や一括整形せず、対象 hunk だけを編集する。

## 2. 履歴上の境界

### 2.1 実在した形式

| 到達点 | session | canonical request | protein raw | protein derived | identity manifest |
|---|---:|---:|---:|---:|---:|
| `main` | 33 | 2 | 2 | 1 | なし |
| branch-only 中間形 | 34 | 3 | 2 | 1 | なし |
| branch-only 中間形 | 35 | 3 | 3 | 2 | 1 |
| 現行 `HEAD` | 36 | 3 | 4 | 3 | 2 |

注意:

- 上表の schema 数値は artifact ごとに別 namespace である。例えば「schema 2」を一括削除してはいけない。
- 現行 manifest 内の `proteinSets` や `recordAnalyses` の subobject が `schema: 1` を使うことと、削除対象の **top-level identity manifest schema 1** は別物である。
- nucleotide LOSAT raw cache schema 2 も現行契約であり、protein raw schema 2 の移行処理と混同してはいけない。

### 2.2 branch-only 世代の導入履歴

| commit | 内容 | 判定 |
|---|---|---|
| `cdbea26` | session v34、canonical request schema 3、feature underlay を導入 | v34 reader だけ削除。現行機能と main v33 からの移行は維持 |
| `73cf163` | session v35、protein raw 3、derived 2、manifest 1、長形式 transport ID を導入 | `main` 未到達 |
| `16762f3` | session v36、raw 4、derived 3、manifest 2、compact runtime handle を導入 | 現行形式。v35 互換層だけ削除 |

履歴確認の代表コマンド:

```bash
git log main -S'CURRENT_SESSION_VERSION = 34' -- gbdraw/session_io.py
git log main -S'CURRENT_SESSION_VERSION = 35' -- gbdraw/session_io.py
git log main -S'SESSION_VERSION = 34' -- gbdraw/web/js/services/config.js
git log main -S'SESSION_VERSION = 35' -- gbdraw/web/js/services/config.js
git log main -S'PROTEIN_LOSAT_CACHE_SCHEMA = 3' -- gbdraw/analysis/protein_colinearity.py
git log main -S'proteinRawV35Candidates' -- gbdraw
git log --all -S'CURRENT_SESSION_VERSION = 35' -- gbdraw/session_io.py
```

最初の 6 コマンドは結果なしであり、最後のコマンドには `73cf163` と `16762f3` だけが現れる。positive control として、`main` には session v33 と raw schema 2 / `p_r_` identity の履歴が存在する。

## 3. 削除対象の全体像

削除対象は次の 3 群に分かれる。

### A. session v34 の受理

v34 専用の大きな migrator はない。実質的な削除対象は以下である。

- Python/JS の supported-version set から `34` を外す。
- v34 を旧形式の fixture 値として使うテストを、実在する v33 に変更する。
- 「v27–34」を migration source とする文書を「v27–33」に変更する。

feature underlay 自体や、`sourceSessionVersion <= 33` の migration は削除しない。

### B. session v35 の LOSAT 移行層

以下はすべて `73cf163` の中間設計を `16762f3` の現行設計へ移すためのコードであり、`main` 由来データには不要である。

- protein raw schema 3 の識別・検証・quarantine・promotion
- protein derived schema 2 の evidence 保存
- top-level protein identity manifest schema 1 の検証
- 長形式 readable transport ID の生成・解析・書き換え
- `proteinRawV35Candidates` / `proteinDerivedV35Evidence`
- v35 candidate の copy-on-success state machine
- v35 専用の UI reference rewrite、snapshot、rollback

### C. v35 API 名を残す transitional alias

現行値を古いフィールド名でも返すだけの alias も `main` には存在しないため削除する。

- `transport_id`
- `binding_hash`
- `binding_hashes`
- `derived_mapping_hashes`
- `query_binding_hash`
- `subject_binding_hash`
- Web helper response の `binding_hash` / `derived_mapping_hash` fallback

## 4. Python/backend の具体的な削除箇所

行番号は監査時点の `16762f3` の目安である。実装時は symbol 名で検索する。

### 4.1 `gbdraw/analysis/protein_colinearity.py`

| 範囲・symbol | 処置 | 理由・注意 |
|---|---|---|
| `LEGACY_PROTEIN_IDENTITY_MANIFEST_SCHEMA` | 削除 | branch-only top-level manifest schema 1 |
| `V35_PROTEIN_LOSAT_CACHE_SCHEMA` | 削除 | branch-only protein raw schema 3 |
| `_TRANSPORT_FIELD_PATTERN`, `_VALID_LOSAT_TRANSPORT_ID_RE` | 削除 | v35 長形式 transport ID 専用 |
| `CdsProtein.transport_id`, `CdsProtein.binding_hash` | 削除 | v35 名の transitional alias |
| `ProteinLosatPairIdentity.query_binding_hash`, `.subject_binding_hash` | 削除 | docstring でも raw schema 3 migration alias と明記 |
| `V35ProteinRawCachePromotion`, `V35ProteinRawCachePromotionScan` | 削除 | v35 promotion 専用 |
| `ProteinExtractionResult.binding_hashes`, `.derived_mapping_hashes` | 削除 | 現行 canonical は runtime/display binding |
| `build_losat_transport_id()` | 削除 | v35 長形式 ID generator |
| `is_v35_protein_losat_cache_entry()` | 削除 | raw schema 3 discriminator |
| `_v35_protein_losat_cache_key()` | 削除 | raw schema 3 key |
| `_legacy_binding_transport_ids()` | 削除 | manifest 1 の `transportIds` resolver |
| `validate_v35_protein_raw_entry_references()` | 削除 | raw 3 + manifest 1 validator |
| `promote_v35_protein_raw_cache_entries()` | 削除 | raw 3→4 promotion |
| 上記の assignment、import、`__all__` | 削除 | dead API を残さない |

#### manifest validator の要リファクタ

`validate_legacy_protein_identity_manifest()` をそのまま全削除してはいけない。現行 `validate_protein_identity_manifest()` が、`proteinSets` と `recordAnalyses` の共通検証をこの関数経由で再利用している。

実装時は次の順にする。

1. `proteinSets` と `recordAnalyses` の schema-neutral な共通 validator を private helper として抽出する。
2. current manifest schema 2 validator をその helper に直接接続する。
3. top-level schema 1、`recordInstances`、`transportIds`、旧 `bindingHash` の検証だけを削除する。
4. current manifest の hash、membership、runtime/display binding、runtime handle の一意性検査が弱まっていないことをテストする。

#### このファイルで維持するもの

- `PROTEIN_LOSAT_CACHE_SCHEMA = 4`
- `LEGACY_PROTEIN_LOSAT_CACHE_SCHEMA = 2`
- `LOSAT_CACHE_SCHEMA` import alias。これは `main` にも存在する。
- `percent_encode_losat_transport_field()`。現行 TSV export hydration が使用する。
- raw schema 2 / `p_r_` 用の `LegacyProteinRaw*`
- `build_legacy_protein_reference_map()`
- `promote_legacy_protein_raw_cache_entries()`
- legacy schema 2 の token/FASTA/reverse-complement 検証

### 4.2 `gbdraw/session_io.py`

削除対象:

- `SUPPORTED_SESSION_VERSIONS` の `34`, `35`
- v35 raw/derived/manifest/candidate constants
- v35 長形式 transport ID regex
- `validate_session()` の `version == 35` dispatch
- `classify_raw_losat_cache_entry()` の `protein-v35`
- `validate_v35_session_artifacts()`
- current validator の `proteinRawV35Candidates` / `proteinDerivedV35Evidence`
- `normalize_current_session_artifacts()` の v35 引数、manifest 1、raw 3、derived 2 の処理
- `_is_valid_v35_protein_identity_manifest()`
- v35 candidate envelope/evidence の validator と normalizer
- `build_session_json()` の v35 引数
- 関連 export

変更後の supported set:

```python
frozenset({27, 28, 29, 30, 31, 32, 33, CURRENT_SESSION_VERSION})
```

`CURRENT_SESSION_VERSION` は 36 のままにする。

維持対象:

- v27–33 session reader
- canonical request schema 1/2→3
- protein raw schema 2 candidate envelope
- protein derived schema 1 evidence
- nucleotide raw schema 2
- v27–30 repeat rectangle migration

### 4.3 `gbdraw/api/request_render.py`

削除対象:

- v35 constants/functions の import
- `PreparedDiagramRequest`, `RequestRenderResult`, `_PreparedLinearArtifacts` の v35 candidate/evidence fields
- `_V35_PROTEIN_REFERENCE_RE`
- `_v35_candidate_envelope()`
- `_promote_v35_candidates()`
- `_build_v35_protein_reference_map()`
- `_v35_derived_evidence_entries()`
- generic reference rewrite/collection 内の v35 分岐
- `_prepare_linear_artifacts()` 内の raw 3 分類、promotion、reference rewrite、warning、result propagation

維持対象:

- raw schema 2 candidate の promotion
- `p_r_` reference map/rewrite
- derived schema 1 evidence
- legacy candidate の pair-local miss と copy-on-success

### 4.4 `gbdraw/cli_utils/session.py`

削除対象:

- `DiagramRunResult` の v35 candidate/evidence fields
- `build_session_json()` への v35 propagation
- `legacyArtifacts.proteinRawV35Candidates`
- `legacyArtifacts.proteinDerivedV35Evidence`

削除後に `copy` などの import が未使用になった場合のみ整理する。

## 5. Web/JavaScript の具体的な削除箇所

### 5.1 `gbdraw/web/js/app/losat-cache.js`

削除対象:

- `V35_PROTEIN_LOSAT_CACHE_SCHEMA`
- `V35_LOSAT_DERIVED_CACHE_SCHEMA`
- `V35_PROTEIN_IDENTITY_MANIFEST_SCHEMA`
- v35 candidate schema constant
- `isV35ProteinRawLosatCacheEntry()`
- `'protein-v35'` classification
- `validateV35ProteinIdentityManifest()`
- `buildV35ProteinReferenceMap()`
- `validateV35ProteinRawEntryReferences()`
- v35 candidate envelope の create/normalize/transition/serialize 一式
- protein derived schema 2 を migration source として許可する分岐

維持対象:

- current raw 4 / derived 3 / manifest 2
- main 由来 raw 2 / derived 1
- generic legacy candidate envelope と promotion

### 5.2 `gbdraw/web/js/services/config.js`

削除対象:

- `SUPPORTED_SESSION_VERSIONS` の `34`, `35`
- v35 constants/helpers の import
- `validateSessionLosatArtifacts()` の v35 分岐
- `quarantineV35ProteinArtifacts()`
- `normalizeV35DerivedEvidence()`
- `legacyV35RawCandidates`, `legacyV35DerivedEvidence`
- state の import/export/reset/snapshot にある v35 fields
- `proteinRawV35Candidates`, `proteinDerivedV35Evidence`

維持対象:

- `sourceSessionVersion <= 33` の feature rendering migration
- imported linear track slot migration
- canonical request schema 1/2 reader
- `promoteGallerySessionToCanonicalV3()`
- main 由来 raw 2 / derived 1 の `legacyArtifacts`

### 5.3 `gbdraw/web/js/app/run-analysis.js`

削除対象:

- v35 helper imports
- `V35_PROTEIN_REFERENCE_RE` と v35 reference collector
- v35 transport ID rewrite
- v35 candidate/evidence state の destructuring、snapshot、rollback
- Python handle `promoteV35ProteinCache`
- `tryPromoteV35ProteinEntry()`
- v35 source manifest による UI reference migration
- transaction の `kind === 'v35'` 分岐
- v35 candidate/evidence の commit/cleanup/cache-clear
- `migratedFromSchema: 3`

また、Python helper の戻り値について以下の fallback を削除し、現行名だけを使う。

- `res.binding_hash`
- `res.derived_mapping_hash`
- ローカルの `derivedMappingHash`

現行 canonical fields は `runtime_binding_hash` と `display_binding_hash` である。

維持対象:

- `LEGACY_PROTEIN_REFERENCE_RE`
- raw schema 2 / `p_r_` promotion
- generic legacy transaction の copy-on-success、cancel/error rollback

### 5.4 `gbdraw/web/js/app/python-helpers.js`

削除対象:

- `promote_v35_losatp_cache_candidates()`
- serialized `CdsProtein` の `transport_id`, `binding_hash`
- extraction result の `binding_hash`, `derived_mapping_hash`
- reconstruction 時の同名 optional fields

維持対象:

- `promote_legacy_losatp_cache_candidates()`
- `runtime_handle`
- `runtime_binding_hash`
- `display_binding_hash`

既存 docstring に古い変換先 schema が残る場合は、raw schema 2→4 に更新する。

### 5.5 `gbdraw/web/js/state.js`

削除対象:

- `legacyProteinRawV35Candidates`
- `legacyProteinDerivedV35Evidence`
- 関連 export

`colorScopeDialog` の “Backward-compatible alias fields” は今回の対象ではない。この alias は `2553998` で導入され、`main` の祖先に存在する。

### 5.6 v35 ID の display-boundary regex

次のファイルには、v35 長形式 transport ID を UI/derived output に出さないための regex がある。

- `gbdraw/web/js/app/losat-cache.js`
- `gbdraw/web/js/app/feature-utils.js`
- `gbdraw/web/js/services/standalone-interactivity.js`
- `gbdraw/web/js/services/standalone-interactivity-assets.js`

これは旧 session を読み込む migrator ではなく、防御的な漏洩検査でもあるため、機械的に削除しない。

`standalone-interactivity-assets.js` は `gbdraw/render/interactive_svg.py` が standalone SVG に埋め込む runtime の実体である。`standalone-interactivity.js` だけを変更しても、生成 SVG の runtime は変わらない。

推奨方針:

1. v35 session/reference の reader と rewrite は削除する。
2. display-boundary の現行契約を `h_...`、`p_r_...`、feature analysis ID で引き続きテストする。
3. v35 regex を残すなら、「unsupported historical shape の防御」であることをコメントに明記する。
4. dead design の denylist も完全に消す方針なら、v35 case だけを削除し、現行・main 由来 internal ID の非表示テストは残す。

この 4 ファイルの regex は、互換削除に必須の変更ではない。次セッションでは migration 本体の削除と別 commit hunk として判断する。

## 6. 非 LOSAT 領域の監査結果

branch-only の非 LOSAT 後方互換として確認できたものは **session v34 の受理だけ**である。

以下は名前やコメントに `legacy` / `migration` / `alias` があっても削除対象ではない。

### 6.1 feature underlay / repeat shape

維持するもの:

- `FeatureShape`
- `create_feature_dict()`
- `_precalculate_feature_dicts()`
- Python の v27–30 repeat rectangle migration
- request schema 1/2 の repeat rectangle migration
- Web の `sourceSessionVersion <= 33` feature rendering migration
- highlight/feature underlay の現行実装

理由: `main` v33 以前の描画再現性を保つために必要である。v34 の reader を削除することと、v33→current migration を削除することは別である。

### 6.2 gallery session migration

維持するもの:

- `gbdraw/web/js/services/gallery-session-migration.js`
- `preflightSessionImport()` からの `promoteGallerySessionToCanonicalV3()`
- `bf45ce2` の legacy gallery settings 保持
- feature metadata の legacy SVG ID migration
- `makeLinearRenderedFeatureId()`

理由: `main` の tracked gallery session は v33 / renderRequest schema 2、現行 gallery は v36 / schema 3 である。これは branch-only v34/v35 のためではなく、`main` 配布物を直接読むための migration である。

### 6.3 その他の非対象

- feature color rule の `_record_N` 付き/なし hash candidates
- label alias
- diagram assembly の backward-compatible wrapper
- `gbdraw/layout/linear.py` の legacy middle-layout factor
- public orthogroup membership mode aliases

これらは `main` 以前の API/出力契約に由来する。今回の branch-only 監査を根拠に削除してはいけない。

## 7. テストと fixture の整理

### 7.1 完全削除

- `tests/fixtures/sessions/BGC0000708-BGC0000713.schema-v3.gbdraw-session.json.gz`
- `tests/web/session-v35-manifest-preservation.test.mjs`

この schema-v3 fixture は branch-only v35 形式専用である。

### 7.2 v35 部分だけ削除

- `tests/test_protein_colinearity.py`
  - long transport ID generator
  - raw schema 3 promotion
  - transitional alias assertions
- `tests/test_api_request_render.py`
  - v35 manifest helper
  - v35 candidate/reference migration
- `tests/test_session_io.py`
  - v35 constants/helper/validation/round-trip/promotion
- `tests/test_losat_legacy_fixture.py`
  - schema-v3 fixture case
- `tests/web/losat-cache.test.mjs`
  - manifest 1/raw 3/derived 2/v35 candidate state
- `tests/web/session-losat-cache-validation.test.mjs`
  - v35 validation cases
- `tests/web/losat-cache-migration.playwright.spec.js`
  - v35 fixture、save-before-generate、promotion、rollback scenario
- `tests/run_losat_cache_browser_acceptance.py`
  - v35 fixture と v35 assertion/scenario
- `tests/test_web_packaging.py`
  - v35 transport ID 専用 case。display-boundary regex を残す場合はこの case も残せる。
- `tests/web/feature-display-boundaries.test.mjs`
  - v35 transport ID case の扱いは 5.6 の方針に合わせる。
- `tests/web/orthogroups-stable-identity.test.mjs`
  - v35 transport ID case の扱いは 5.6 の方針に合わせる。
- `tests/web/interactive-svg-search-performance.playwright.spec.js`
  - embedded standalone runtime の v35 long-ID display-boundary case。regex を残すなら意図的に維持し、削除するなら current/main internal ID case へ置換する。

### 7.3 fixture 値の変更

- `tests/fixtures/sessions/BGC0000708-BGC0000713.schema-v2.expected.json`
  - nested `v35` block を削除
  - v35 長形式 ID の regex expectation は 5.6 の方針に合わせる
- `tests/test_session_io.py`
  - v34 を使う synthetic legacy case を v33 に変更
  - test 名の `v27_to_v34` を `v27_to_v33` に変更
- `tests/test_refresh_gallery_sessions.py`
  - stale-version の例として使う v34 を v33 など実在値に変更
- `tests/web/session-authority.test.mjs`
- `tests/web/gallery-session-migration.test.mjs`
  - v33 fixture に branch-only top-level manifest schema 1 と `legacyArtifacts` を混在させている assertion を削除するか、current v36 の独立 fixture/test に分離する
  - v33 fixture の manifest だけを current schema 2 に置換してはならない。それも version と artifact contract が不整合になる。
  - authority/gallery migration test 自体は維持する

### 7.4 必ず維持する fixture と acceptance

- `tests/fixtures/sessions/BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz`
- 上記 fixture の raw schema 2 / `p_r_` verified migration
- 25 pair が `25 hits / 0 misses / 0 jobs` になる browser acceptance
- derived schema 1 の復元
- cancel/renderer error 時の generic schema 2 promotion rollback

この schema-v2 fixture は session v33 の実在データであり、`main` 由来互換性の positive control である。

## 8. tool と生成物

### 8.1 `tools/refresh_gallery_sessions.py`

削除対象:

- `proteinRawV35Candidates` の postcondition
- v35 source fixture/candidate を前提にした検査

長形式 transport ID の漏洩検査は 5.6 の方針に合わせる。現行 `h_...`、main 由来 `p_r_...`、feature analysis ID の漏洩検査は維持する。

### 8.2 tracked gallery

監査時点の current gallery session は v36 / request schema 3 / raw 4 / derived 3 / manifest 2 へ更新済みである。v35 専用の tracked **session artifact** はない。

ただし、v35 long-ID の防御 regex は `standalone-interactivity-assets.js` から 11 本すべての `gbdraw/web/gallery/examples/*.svg` に埋め込まれている。5.6 で regex の削除を選ぶ場合は standalone asset と gallery examples を生成ツールで同期する。`gbdraw/web/gallery/sources/*.svg` にはこの embedded runtime はない。

実装削除だけで diagram geometry は変わらないはずである。wheel/cache-bust や session serializer の変更で gallery refresh が必要になった場合も、生成ツールを使い、JSON/SVG を手編集しない。

## 9. 既存ドキュメントの更新対象

以下から v34/v35 を公開済み migration history として扱う記述を削除する。

- `docs/FAQ.md`
- `docs/TUTORIALS/4_Protein_Comparisons.md`
- `docs/TUTORIALS/8_Interactive_SVG_Sessions.md`
- `docs/RELEASE_NOTES_0.14.0b0.md`
- `docs/PYTHON_SESSION_COMPATIBILITY_MATRIX.md`

推奨表現:

```text
Supported session versions are 27–33 and 36.
Versions 34 and 35 were branch-internal development formats and are not supported.
```

特に release note の「session 34 で新 default を記録した」という履歴表現は、現行 v36 の仕様として書き直す。

次の 2 文書は v35 migration を設計の中心にしている。

- `docs/LOSAT_CACHE_IDENTITY_AND_LINEAR_SPACING_REGRESSION_IMPLEMENTATION_PLAN.md`
- `docs/LOSAT_COMPACT_RUNTIME_HANDLE_MIGRATION_IMPLEMENTATION_PLAN.md`

選択肢:

1. 現行 invariant と main schema 2 migration だけを残して全面改稿する。
2. 冒頭に「unmerged, superseded design」と明記して historical record として残す。
3. 他の文書から参照されていなければ削除する。

現行仕様書と誤読される状態のまま v35 acceptance criteria を残してはいけない。

## 10. 次セッションの推奨実装順

### Phase 1: contract を固定

1. `CURRENT_SESSION_VERSION = 36` を維持する。
2. Python/JS の supported set を `27–33, 36` にする。
3. v34/v35 が unsupported として拒否されるテストを先に置く。
4. main schema 2 fixture の positive-control test は残す。

### Phase 2: Python model と validator

1. current manifest の共通 validator を schema-neutral helper へ抽出する。
2. transitional dataclass fields/properties を削除する。
3. v35 manifest/raw detector と promotion を削除する。
4. current と main schema 2 の focused test を通す。

### Phase 3: session/API/CLI propagation

1. `session_io.py` の v35 envelope/evidence を削除する。
2. `request_render.py` の v35 preparation/rewrite/propagation を削除する。
3. `cli_utils/session.py` の v35 result/save fields を削除する。
4. Python export と未使用 import を整理する。

### Phase 4: Web state と orchestration

1. `losat-cache.js` の v35 discriminators/validators/state machine を削除する。
2. `config.js` の quarantine/import/export/reset を削除する。
3. `run-analysis.js` の v35 promotion/reference/transaction を削除する。
4. `python-helpers.js` と `state.js` の aliases/state を削除する。
5. Python worker interface と JS caller を同じ change set で更新する。

### Phase 5: tests、fixture、docs

1. v35 fixture と専用 test を削除する。
2. mixed test から v35 case だけを削除する。
3. v34 synthetic values を v33 に変更する。
4. main schema 2 acceptance を実行する。
5. 文書を `27–33, 36` 契約へ更新する。

### Phase 6: optional defensive regex

5.6 の判断を独立して行う。dead v35 denylist を削除する場合でも、現行 internal ID の display-boundary test を弱めない。削除を選んだ場合は `standalone-interactivity-assets.js`、embedded-runtime Playwright test、11 本の generated gallery example SVG を同時に更新する。

## 11. 検証チェックリスト

### 11.1 source-level gate

production path から次の v35 migration symbol が消えていることを確認する。

```bash
rg -n \
  'V35|v35|proteinRawV35Candidates|proteinDerivedV35Evidence|legacy_protein_.*v35|promote_v35' \
  gbdraw/analysis gbdraw/api gbdraw/cli_utils gbdraw/session_io.py \
  gbdraw/web/js/app gbdraw/web/js/services gbdraw/web/js/state.js tools
```

display-boundary regex を意図的に残す場合、その残存だけを review する。
特に `standalone-interactivity-assets.js` と generated gallery SVG 内の同じ regex を見落とさない。

```bash
rg -l -F '~f_[0-9a-f]{64}' \
  gbdraw/web/js/services/standalone-interactivity-assets.js \
  gbdraw/web/gallery/examples
```

transitional alias の残存も確認する。

```bash
rg -n \
  'build_losat_transport_id|query_binding_hash|subject_binding_hash|binding_hashes|derived_mapping_hashes|derived_mapping_hash' \
  gbdraw tests
```

`binding_hash` という部分文字列を一括禁止してはいけない。現行 `runtime_binding_hash` / `display_binding_hash` は必要である。

### 11.2 focused Python

```bash
python -m pytest \
  tests/test_protein_colinearity.py \
  tests/test_session_io.py \
  tests/test_api_request_render.py \
  tests/test_losat_legacy_fixture.py \
  tests/test_session_request_codec.py \
  tests/test_refresh_gallery_sessions.py -v
```

### 11.3 Web

最低限、次を実行する。

- browser wheel を使う検証の前に `python tools/prepare_browser_wheel.py` を実行する。生成 wheel は commit しない。
- LOSAT cache unit tests
- session LOSAT validation tests
- session request / authority / gallery migration tests
- schema 2 migration Playwright test
- schema 2 browser acceptance
- feature display-boundary tests

期待条件:

- v27–33 と 36 だけを受理
- v34/v35 を拒否
- main schema 2 fixture が raw 4 へ verified promotion
- `25 hits / 0 misses / 0 jobs`
- derived schema 1 を復元
- current raw 4 / derived 3 / manifest 2 round-trip
- failed promotion/cancel/render error で state rollback
- gallery request schema 2→3 promotion を維持

### 11.4 repository-wide

```bash
ruff check gbdraw/
python -m pytest tests/ -v -m "not slow"
python -m build
```

diagram geometry に意図しない変更がないことを、必要に応じて read-only reference comparison でも確認する。

## 12. 完了条件

次のすべてを満たしたら branch-only compatibility cleanup は完了である。

- Python と Web の supported session が 27–33 と 36
- v34/v35 reader、quarantine、promotion、candidate state がない
- v35 top-level manifest 1、protein raw 3、protein derived 2 の reader がない
- v35 長形式 transport ID の generator/rewrite がない
- transitional Python/worker aliases がない
- v35 test fixture と専用 test がない
- main 由来 schema 2 fixture と zero-job acceptance が通る
- request schema 1/2、feature rendering、gallery、旧 feature ID の migration が残る
- current manifest validator の hash/membership/identity 検査が弱まっていない
- user-facing docs が v34/v35 を公開済み互換対象として説明していない
