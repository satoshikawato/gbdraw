# 設計: 個別 override を元の feature の identity で指す（Q4、2026-10-04）

状態: 設計（実行時のコードは含まない）。実装は下の PR 分割に従う。

対象は `dev` の `3bdfbadb`（#757、#759、#760 を含む）。所見の ID（OV-xx）、根本原因（RC-x）、Owner の決定 Q1〜Q4 は
[FINDINGS.md](./FINDINGS.md)（PR #758）による。Feature placement の mode 別の request 投影は PR #762（未 merge）を前提にする。
JS のパスは `gbdraw/web/js/` からの相対パス、Python のパスはリポジトリからの相対パス。

## 0. 決定と範囲

Owner の決定（2026-10-04、Q4 = A、原文）:

> Feature placement と同じく (record_key, 元の feature の ID と instance) を key にし、Python が元の catalog から解決する。
> crop・逆相補・並べ替え・複製でも同じ feature に当たり続ける。Session と生成表の形式が変わるので、設計文書を先に作ってから実装する。
> 色の This feature only は PD-OI-069 のまま。

範囲（Web の popup と Editor で作る feature 単位の編集）:

- **Feature visibility**: On / Off / Exclude from matching
- **Label visibility**: On / Off、および label の文字（"Keep hidden (apply text only)" を含む）
- 選択した feature から作る annotation（`app/annotations.js::addSelectedFeatures`）

範囲外: 色の「This feature only」（§3.5 で推奨を述べる）、stroke、qualifier 単位の editor 規則（product / protein_id の scope）、
手入力の rule 表、bulk の label 置換そのもの（§3.3 で扱いを決める）。

解く所見: OV-01〜OV-03（RC-1。OV-02 は §3.6）、OV-04 と OV-05 の色以外の部分（RC-2。色は #760 で直した）、OV-11、OV-12（RC-5）、
および文字だけの行が label を強制する潜在的な危険。

## 1. identity

### 1.1 型

identity は placement と同じ組 `(record_key, biological_feature_id)` とし、新しい instance 欄は足さない。この組だけで複製も重複も区別できる。

- `record_key`: request の record ごとの key。`gbdraw/api/record_planning.py` が provenance に入れ、ALL の record input は
  `<key>:<source record index + 1>` に展開する。描画した record の `annotations["gbdraw_record_key"]` にも同じ値が入る。
  - Linear: `services/session-request.js::buildRecords` が `seq.uid`（`state.js` が作る `linear-seq-<uuid>`）を使う。
    Session を読むと uid は `recordKey` から戻る（`services/session-request.js` の `uid: String(record.recordKey || …)`）。
    並べ替えても変わらない。同じファイルを 2 回読み込むと、slot ごとに別の uid になる。
  - Circular: 1 record の表示は `circularRecordKey`（`circular-<recordId>-<selector>`、または保存した `recordKey`）。
    multi-record canvas と batch は、ファイル中の全 record をファイル順に並べて `record-<n>` とする
    （`app/record-options.js::resolveCircularRequestRecordSet`）。同じ record ID が 2 回あっても key は別になる。
- `biological_feature_id`: `gbdraw/features/source.py::build_source_feature_catalog` が crop の前の元 record から作る。
  - 値は元座標の hash（`compute_feature_hash_from_location_parts`、record ID・type・座標から計算）。
  - 同じ record に同じ hash の feature が 2 つ以上あると、`gbdraw/features/ids.py::disambiguate_feature_ids` が `<hash>~<source_feature_index>` にする。
  - Web の catalog（`gbdraw/web_support/feature_catalog.py::_normalized_biological_features`）も同じ規則で同じ値を出す。
    rendered feature ごとに `svgId → (recordKey, biologicalFeatureId)` を必ず持たせ、持たせられなければ `GbdrawError` にする（`_normalized_rendered_features`）。
- 確かめた事実（`3bdfbadb` で実行）: 同じ座標の CDS が 2 つある record の source catalog は `fbd0c2a00~0`、`fbd0c2a00~1` を持つ。
  一方、`SourceFeatureIdentity.matches(key="hash", …)` は `stable_feature_id`（`~` の付かない値）としか比べない。
  そのため `hash=fbd0c2a00~1` は 0 件、`hash=fbd0c2a00` は 2 件に当たる。
  - 結果として `--feature_placement_table` は重複の片方を指せない。
  - `app/run-info.js` の Source recipe は placement 行を `feature_selector=<biologicalFeatureId>`（key なし）で書くので、重複の feature では CLI が失敗する。
  - `docs/REFERENCE/input-formats-and-tsv-schemas.md` の「`hash=<biologicalFeatureId>`」は、この場合に成り立たない。PR-Q4-1 で直す。

### 1.2 既知の制約（今回は変えない）

- Circular では、1 record の表示と multi-record canvas / batch で `record_key` が違う。片方で作った編集は、もう片方では当たらない。
  placement も今は同じ。mode の切り替えと同じく、行は消さずに当たらないだけで、戻せば効く（R2、#762 の `requestFeaturePlacements`）。
  §3.4 の通知で利用者に知らせる。key をそろえるには保存済みの placement の key の移行が要るので、別の課題にする。
- `record-display/feature-anchor.js` の Record rotation の anchor と、`gbdraw/api/record_planning.py::_alignment_source_feature` の Similarity alignment の anchor も同じ組を使う。
  ただし照合の規則が違う（source index / stable id への fallback）ので、今回は共有 resolver に入れない。

## 2. 現状（種類ごと）

| | Feature visibility | Label visibility | label の文字 | 選択から作る annotation | Feature placement（参考） |
|---|---|---|---|---|---|
| Web state と key | `state.featureVisibilityOverrides` `{描画 ID: mode}` | `labelVisibilityOverrides` `{描画 ID: on/off}` | `labelTextFeatureOverrides` と `labelTextFeatureOverrideSources` `{描画 ID: …}` | `state.annotationSets[].annotations[].target`（`featureSpan`） | `featurePlacementOverrides` `{JSON[recordKey, bioId]: row}` |
| 行を作る処理 | `state.js::featureVisibilityRules` = `app/feature-visibility.js::deriveFeatureVisibilityRulesForBoundary`。selector cache（`buildFeatureVisibilitySelectorCache`）が一意な qualifier（protein_id、locus_tag など）、元座標の `record_location`、元座標の hash の順に選ぶ | `app/feature-editor/label-override-table.js::buildLabelOverrideRows`: `hash` 行（On = 文字あり、Off = 空） | 同じ関数: `selectStableFeatureKey`（locus_tag / gene / 元座標の `record_location`、最後は `hash`） | `app/annotations/target-actions.js::featureTargetsFromSelection`: `hash=<selector.hash>`（元座標） | `services/feature-placement.js::canonicalFeaturePlacements` |
| request | `diagramOptions.featureVisibilityTableFile`（`services/session-request.js::addGeneratedTableResources`）。列は 5 つで、`featureId` は書かない | `diagramOptions.labelOverrideFile`（`addLabelOverrideResource`。`canonicalLabelOverrideRows` → `generatedLabelOverrideTsv` → 生成の順） | 同左 | `diagramOptions.annotations.sets` | `diagramOptions.featurePlacements`（schema 7 以降） |
| Python の照合 | `gbdraw/features/visibility.py::_rule_matches_feature`: 描画した record の `get_feature_hash` / `location` / `record_location` | `gbdraw/labels/filtering.py::get_label_text`: `hash` 行は filter と scope より先。`has_forced_label_overrides` は文字のある `hash` 行を On とみなす | 同じ。`hash` 以外の行は filter の後 | `gbdraw/annotations/resolve.py::_feature_matches`（描画した record の hash） | `gbdraw/features/placement.py::resolve_placement_inputs`: 元 catalog から `source_feature_index` で解決 |
| Session | `features.featureVisibilityOverrides`（`services/config.js::exportSessionDocument`）。読み込みは `applyFeatureStateData`、それより古い形は `splitLegacyVisibilityRules` | `features.labelVisibilityOverrides` | `features.labelTextFeatureOverrides` と `…Sources`、表は `features.labelOverrideRows` | `config.annotationSets` | `config.featurePlacementOverrides`（`session_io.py::_validate_display_placement_drafts` が key を確かめる） |
| History | feature intent（`services/history-snapshot.js::buildFeatureIntentData`）。selector cache は generated owner set | feature intent | feature intent。`canonicalLabelOverrideRows` は artifact checkpoint だけ | config domain | config domain |
| live の投影 | `app/feature-editor/visibility-actions.js::reconcileFeatureVisibility` は描画 ID で引く。`app/candidate-render.js::compilePlanBundle` は `resultIndexesByRenderedId` で引く | `label-actions.js::applyStoredVisibilityOverridesToSvg` は `data-label-feature-id` で引く。`compilePlanBundle` も同じ | `projectLabelTextIntent` は `data-label-feature-id` で引く | なし（Generate で描く） | なし（Generate で描く） |
| 書き出しと読み込み | `visibility-actions.js::downloadFeatureVisibilityRulesTsv`（取り込みの UI はない） | Export Label TSV（`label-actions.js::downloadLabelOverrideTable`。生きている map から作り直す） | 同左。`loadLabelOverrideTable` は文字だけを map に戻す（visibility は戻さない） | `app/annotations.js`（`encodeAnnotationTable` / `parseAnnotationTableWithNotice`） | Run Info の `--feature_placement_table` |
| 刈り取り（R2） | `app/run-analysis.js` が、元の source を置き換えて成功した Generate の中で `pruneUnmatchedFeatureOverrides` を呼ぶ。描画 ID から得た hash で生存を判定する | 同じ（PD-OI-065） | 同じ（bulk は残す） | なし | `app/record-display-options.js` の source watcher が、置き換えた source の行をファイルの置き換え時に全部消す |

そのほかの事実:

- 色と stroke の override は、すでに `services/feature-catalog.js::stableFeatureOverrideKey`（`recordKey\0biologicalFeatureId`）を key にしている。
  `compilePlanBundle` は `catalogAdmission.renderedTargetsByOverrideKey` で、key から各 Result の描画 ID を引く。identity での live 投影は、この仕組みを使い回す。
- Python の CLI `--session` は `features.*` の map を読まない。個別の編集については、`renderRequest` の表と `featurePlacements` と `annotations` だけを使う
  （`gbdraw/cli_utils/session.py::render_canonical_session_if_present`）。
- 死んでいる export（テストからしか呼ばれない）: `app/feature-visibility.js` の `featureVisibilityRulesFromOverrideCache`、
  `buildExactHashFeatureVisibilityRule`、`buildFeatureVisibilityOverrideCache`、`getEditorFeatureVisibilityMode`、`upsertEditorFeatureVisibilityRule`、
  `buildEditorFeatureVisibilityRule`、`removeEditorFeatureVisibilityRule`（`tests/web/feature-visibility.test.mjs` だけが使う）。
- Exclude from matching は、LOSATP のタンパク質の抽出にも効く。
  - Generate 中の抽出: `gbdraw/analysis/protein_colinearity.py` が `should_include_feature_in_analysis` を呼ぶ。
  - browser の helper: `DIAGRAM_HELPER_OPERATIONS.EXTRACT_CDS_PROTEIN_FASTA` には visibility 表を `role: 'visibility'` で渡し、その TSV 文字列を cache key にする（`app/run-analysis.js`）。
- 潜在的な危険: `label-override-table.js` は、文字だけの編集でも、一意な locus_tag / gene / `record_location` がない feature では `hash` 行を作る。
  Python はそれを On とみなす（`has_forced_label_overrides`）。今は行が照合されないので表に出ないが、照合されるようになると "Keep hidden" でも label が出る。

## 3. 目標の契約

### 3.1 Python: identity の型と resolver を 1 つにする

- `gbdraw/features/source.py`（元 catalog の owner）に次を置く。
  - `FeatureIdentity(record_key, biological_feature_id)`（frozen）: 空でない、NUL を含まない（今の `FeaturePlacementOverride.__post_init__` の検査を移す）。
  - `resolve_feature_identities(*, records, record_keys, source_catalogs, identities) -> Mapping[FeatureIdentity, IdentityBinding]`
    - `IdentityBinding(record_index, source_feature_index | None, status)`
    - `status` は次の 4 つ。
      - `present`: 描画した record にその feature がある。
      - `crop_excluded`: 元 catalog にはあり、描画した record になく、record が crop されている。
      - `absent`: 元 catalog にはあり、描画した record にない（GFF の読み込みで落ちた）。
      - `unresolved`: record key は request にあるが、元 catalog にその ID がない。
    - request にない record key は `ValidationError`（`diagnostic=`、code `FEATURE_IDENTITY`）にする。Web は #762 の投影で、今の request の record の行だけを送る。
  - runtime の feature の source index は placement と同じ方法で求める。`_source_feature_index(feature)` がなければ `_feature_source_index_map(record.features)`。
  - `SourceFeatureIdentity.matches`: key が `hash` のとき、`stable_feature_id` に加えて `biological_feature_id` とも一致を見る。`~n` 付きの値は hash と衝突しないので、元座標の表の `hash=` の意味は広がるだけになる。
- `gbdraw/features/placement.py::resolve_placement_inputs` はこの resolver を使い、中の `known` と `runtime` の組み立てを消す。
  表の行の解決（record selector がちょうど 1 record、feature selector がちょうど 1 feature）は `resolve_identity_table_rows` として切り出し、PR-Q4-3 の表と共有する。
  placement 固有の分類（`hidden` / `underlay` / `foreground`）は、`present` / `absent` の上に今のまま残す。
- 抽象を足す条件（CLAUDE.md）: この resolver は PR-Q4-1 で placement の 2 つの経路（正確な行と表の行）をまとめ、中の重複を同じ PR で消す。
  PR-Q4-2 で visibility、label、annotation が同じ resolver に載る。`resolve_identity_table_rows` は PR-Q4-3 まで消費者が 1 つの私的な分割で、PR-Q4-3 で 2 つ目の表が使う。
- per-record の結果は今の placement と同じ経路で渡す。
  - `PreparedDiagramInputs.placements`（records と整列した `ResolvedPlacementInputs`）を `ResolvedFeatureInputs(record_key, placements, overrides)` に広げる。
  - `overrides` は `source_feature_index → FeatureOverride` の対応。
  - 解決は、今の `gbdraw/api/request_render.py::_materialize_placement_inputs` と同じ場所で、Generate ごとに 1 回だけ行う。元 catalog は provenance の `source_feature_catalog` を使う。

### 3.2 request: 型付きの欄にする（表の新しい selector にはしない）

決定: `diagramOptions.featureOverrides` を足す（request schema 9）。label-override 表と feature-visibility 表に新しい selector を入れる案は採らない。理由:

- 表の照合は正規表現で、`record_id` と `feature_type` の列を持ち、`record_key` を表せない。
- label 表の `hash` 行は「文字がある = On」で、文字と visibility を分けられない。潜在的な危険が残る。
- 表は先に当たった行が勝つので、個別の編集と手入力の規則の優先順位を、行の並びで表すことになる。
- `featurePlacements`（schema 7）と同じ形にでき、Python の型付きの層で値を確かめられる（R7）。

行の形（JSON の key は固定。値がないところは `null`。並びは `recordKey`、次に `biologicalFeatureId` の code point 順）:

```json
{"recordKey": "linear-seq-…", "biologicalFeatureId": "fbd0c2a00~1",
 "featureVisibility": "off", "labelVisibility": null, "labelText": null}
```

- `featureVisibility` は `on` / `off` / `exclude_matching` / `null`、`labelVisibility` は `on` / `off` / `null`、`labelText` は文字列か `null`。
- `labelText` は改行・タブ・NUL を拒否する（今の TSV 正規化と同じく 1 行）。
- 3 つとも `null` の行と、同じ identity の行の重複は拒否する。既定（Default）は行がないことで表す（placement の Auto と同じ）。
- Python の型: `gbdraw/features/overrides.py`（新規、`placement.py` と同じ層）の
  `FeatureOverride(record_key, biological_feature_id, feature_visibility=None, label_visibility=None, label_text=None)`。
  `CircularDiagramOptions` と `LinearDiagramOptions` に `feature_overrides: tuple[FeatureOverride, ...]` を足す。
- codec: `gbdraw/session_request_codec.py::_decode_diagram_options` は schema 9 以上でこの配列を必須にし、8 以下では `()` にする。
  `_encode_diagram_options` は `featurePlacements` の隣に書く。JS は `services/session-request.js` の reader / writer と `promoteCanonicalRenderRequestToCurrent` に同じ規則を入れる。
- annotation の対象に `{"kind": "featureIdentity", "recordKey", "biologicalFeatureId", "envelope", "circularPath"}` を足す。
  Python 側は `gbdraw/annotations/models.py::FeatureIdentitySpan`。1 つの対象は 1 feature（選択から作る対象は今も feature ごとに 1 つ）。
- Web から request への投影は #762 の `requestFeaturePlacements` と同じく、今の request の record（ALL の `<key>:<n>` 展開を含む）に属する行だけを載せる。
  ほかの mode の行は draft に残す（R2）。

### 3.3 Python の消費者と優先順位

identity の行は、その feature についての表の規則より先に見る。どの消費者も、§3.1 の per-record の対応を source index で引く。

| 消費者 | 変更 |
|---|---|
| `features/visibility.py::should_render_feature` | `on` → 描く、`off` → 描かない、`exclude_matching` → 表を見ずに type の選択で決める（今の editor 行が manual 規則より前にあるのと同じ結果） |
| `features/visibility.py::should_include_feature_in_analysis` | `off` と `exclude_matching` → 除く |
| `features/factory.py`（描画）、`web_support/feature_metadata.py::extract_features_from_records_payload`（catalog の描画集合）、`features/placement.py`（状態）、`analysis/protein_colinearity.py`（Generate 中の抽出） | 上の 2 つを per-record の対応つきで呼ぶ |
| `EXTRACT_CDS_PROTEIN_FASTA` helper | その record の行を受け取り、同じ resolver で解決する（R4）。cache key には visibility 表の代わりに正規化した行を入れる（OIPC-C04） |
| GFF の読み込み（`request_render.py` の `resolve_candidate_feature_types`） | GFF の record に `featureVisibility: on` の行が 1 つでもあれば `gff_keep_all_features` にする（今は表の `feature_type` 列で同じ効果を得ている）。読み込み時間の測定を PR に付ける |
| `labels/filtering.py::get_label_text` | `labelVisibility: off` → label なし。`on` → scope・whitelist・blacklist より先に表示する。文字は `labelText`、なければ通常の文字（qualifier の優先順で選び、表の非 `hash` 規則を当てた文字）。それも空なら `<type> <location>`（今 Web の `resolveDefaultLabelText` が作る代替を Python に移す）。`labelVisibility: null` で `labelText` があれば、label が通常どおり出るときだけ文字を変える（表の非 `hash` 規則より先）。**文字だけの行は label を強制しない** |
| `config/models/render_profiles.py` の `forced_labels` | 表の `hash` 行に加えて、`labelVisibility: on` の identity 行も数える |
| `annotations/resolve.py` | `FeatureIdentitySpan` は解決した source index の feature から区間を作る。解決できなければ今の `feature_selector_unmatched` 警告（`missing_count`）で飛ばす。新しい警告の語は足さない |

表の `hash` / `location` / `record_location` の規則は意味を変えない（CLI の契約）。Web はこれらの行を個別の編集のためには作らなくなる（PR-Q4-4）。

bulk の label 置換（`labelTextBulkOverrides`）は source text の規則なので、表に残す。
今は対象の feature ごとの行に展開していて、元座標の selector を使う。この展開の結果を identity 行の `labelText` にする。
同じ identity に個別の文字があれば、そちらが勝つ（今の行の並びと同じ結果）。展開できないときの `* * label ^text$` 行は表に残す。

### 3.4 解決できない identity（OIPC-C03: 黙って無視しない）

- Python は `present` 以外の行を通知として返す。
  - `Diagram.feature_identity_notices`
  - Web の metadata `featureIdentityNotices`: `{recordKey, biologicalFeatureId, status, kinds, resultIndex}` の配列
  - CLI は 1 行ずつ log に出す（`annotation_warnings` と同じ経路: `gbdraw/web_support/request_render.py`、`gbdraw/circular.py`、`gbdraw/linear.py`）
- `crop_excluded` と `absent` は休眠とし、行は残す。crop を戻すと再び効く（placement の今の扱いと同じ）。
  Web は Result の近くに「N 件の feature の編集は今の crop / 表示の外にあります」と出す。
- `unresolved`（Q3 の推奨 A）: Generate は成功させ、通知する。
  - Web は、元の source を置き換えて成功した Generate の中の reconcile（`pruneUnmatchedFeatureOverrides`、R2）で、置き換えた source の record の `unresolved` 行を消す。
  - 同じ Generate で、前の committed request にあって今の request にない同じ mode の record の行も消す。
  - どちらも消した件数を通知に出す。これは PD-OI-065 の結果（もう無い feature への編集だけを外し、残る feature への編集は保つ）を、visibility と placement にも当てはめたものになる。
  - annotation は利用者が作った set なので消さない。解決できない対象は、今の `feature_selector_unmatched` 警告で飛ばす（§3.3）。
  - source を置き換えていない Generate では消さず、"Remove N unmatched feature edits" の操作を通知に付ける（明示の Reset なので R2 に沿う）。
- 判定は Python の状態だけで行う。Web の catalog には crop で外れた feature がないので、Web だけでは `crop_excluded` と `unresolved` を区別できない。

### 3.5 色の「This feature only」: 動かさない（推奨）

- PD-OI-069（`A / STABLE-HASH-ONLY`）と Owner の Q4 の文言に従い、色は今の `hash` の色規則のままにする。
  - 規則（`manualSpecificRules`）は rule 表に見えていて、CLI の色表と同じ形で書き出される。identity 行に移すと、色規則の表し方が 2 つになる。
  - 複製で色を共有するのは、PD-OI-069 が受け入れた残余リスク。
- 残るリスク: crop や逆相補を変えると描画 hash が変わり、色規則が黙って外れる（OV-11 の色の部分）。
  対処の案（新しい Product 判断が要る）: PR-Q4-4 の後に、当たる feature がなくなった単独 feature の色規則（`hash` 行）を、live 照合（`app/rule-matching.js`）の結果で数えて知らせる。
  もし Owner が crop に強い色を望むなら、PD-OI-069 の改訂として identity 行に移す。resolver は共有なので、移すのは小さな変更で済む。
- 複製では、visibility と label は複製ごと（Q1）、色は共有、という差が残る。docs の Feature presentation 節に書く。

### 3.6 OV-02: 手入力の規則の live 照合

OV-02 は identity の編集ではない。手入力の色規則（`location` / `record_location`）の live 照合で、Web が Python に渡す値の座標系が違う。

- 今の動き: `app/rule-matching.js::ruleFeaturePayload` は、`hash` だけを描画した record の値（`getFeatureColorRuleHash`）で渡す。
  `location` と `record_location` は元座標の値（`normalizeFeatureSelectorMetadata`）で渡す。
  Label TSV の取り込み（`label-actions.js::loadLabelOverrideTable`）も同じ payload を使う。
- 決定: Python は catalog の rendered feature に、描画した record の selector の値 `drawnSelector: {hash, location, recordLocation}` を付ける（feature catalog schema 5）。
  - 値は `gbdraw/web_support/feature_metadata.py::_biological_selector_values` が元座標の値で上書きする前の `build_feature_selector_values` の結果で、新しい計算は要らない。
  - `ruleFeaturePayload` はこの 3 つを使う。Web 側で座標を変換しない（R4）。
  - catalog 3 / 4（古い Session）にはこの値がない。その feature の `location` と `record_location` の規則は、次の Generate まで live で照合しない（pending。R4 の「辞退」）。
- catalog の版は Session の版に結び付いている（`services/session-authority.js` の `version === 44 ? 4 : 3`、`session_io.py` の catalog schema の判定）。
  そのため catalog 5 は Session 45 と同じ PR-Q4-4 で入れる。

## 4. 永続化の形式

### 4.1 事実

- `origin/main`（`fe6861f0`）: Session 44、request schema 8（受け付けるのは 1、2、5、6、7、8）、catalog 4（legacy 3）。
  `features.featureVisibilityOverrides` などの描画 ID の map と `diagramOptions.featurePlacements` は main にある
  （`featurePlacements` は first-parent の `4e8c9380` から）。
- 最新の release tag は `0.13.0`（Session 30、canonical request なし）。`0.14.0` の tag はない。
- main に最初に現れた commit: Session 31 は `10d3a3d2`、33 は `6b89c781`、39 は `17e2c9de`、40 は `8228ffab`、44 は `fe6861f0`。
- 既存の fixture の個別 override の map は、すべて空。
  中身があるのは audit の probe の Session だけ（`/home/kawato/gbdraw-baselines/override-precedence-audit-20261004/probes/override-audit-b/`
  の `b-crop.gbdraw-session.json` と `b-circ-dup.gbdraw-session.json`、gzip、v44 / schema 8、`#757` head の dev で保存）。

### 4.2 変える形式

| 名前空間 | 変更 | PR | 互換の reader / migrator | 根拠 |
|---|---|---|---|---|
| canonical render request | 8 → 9: `diagramOptions.featureOverrides`、annotation の `featureIdentity` 対象 | PR-Q4-2 | schema 8 を読み、`featureOverrides: []` を補う（`promoteCanonicalRenderRequestToCurrent` を 5〜8 → 9 に、Python は decode → 再 encode） | schema 8 は main にある |
| Session envelope | 44 → 45: `features.featureOverrides`（`{JSON[recordKey, bioId]: row}`）が 4 つの描画 ID の map を置き換える。`runMetadata.featureIdentityNotices` を足す | PR-Q4-4 | Session 44 → 45 の key の migrator（§4.3） | 40〜44 は main にある |
| Session 31〜33、39〜42 | 今の読み込み経路の後で、同じ migrator を通す | PR-Q4-4 | §4.3 | main にある |
| feature catalog | 4 → 5: rendered feature に `drawnSelector`（§3.6） | PR-Q4-4 | catalog 3 と 4 を読み、`drawnSelector` がなければ live 照合を辞退する（移行で値を作らない） | catalog 3 と 4 は main にある |

- PR-Q4-2 では Session を 44 のままにする。dev では request 7 → 8 を Session 44 のまま変えた前例がある。
  PR-Q4-2 の時点で main に昇格しても、request 9 は完成した形なので、後から読み替えは要らない。
- Session 45 は PR-Q4-4 だけで導入する。PR-Q4-4 を分けるなら、全部が merge するまで dev から main への昇格を止める。
  45 の途中の形には reader を作らない（CLAUDE.md の「Persisted-format compatibility」）。
- PR-Q4-4 の Session 45 の reader は、`config.annotationSets` の `featureIdentity` 対象をすでに受け付ける
  （`services/session-active-config-contract.js` の行の検査）。そうすれば PR-Q4-5 は Session の形を変えずに済む。
- 互換の経路が増えるので（下記）、PR-Q4-2 と PR-Q4-4 は `ARCHITECTURE_EXCEPTION`。正確な head について maintainer の決定を得る。
  - PR-Q4-2: request の名前空間で、schema 8 の reader が 1 増える（`CB +1`）。Web の個別編集の表の経路は、PR-Q4-4 で消すまで残る（`PE` を一時的に +1。削除の条件は PR-Q4-4 の merge）。
  - PR-Q4-4: Session の名前空間で、44 の reader と key の migrator が 1 増える（`CB +1`）。次を消す。
    - 表の経路（`PE -1`）
    - selector cache の owner
    - placement の source watcher

### 4.3 Session の key の migrator（PR-Q4-4）

古い key K（描画 ID）を identity にする規則。保存してある catalog（v39 以降は `editorState.featureCatalog`。catalog 3 を 4 にする移行は anchor profile を足すだけ）を使う。

1. K が catalog の rendered feature の `svgId` にあれば、その ID で描いたすべての identity にする（全 Result）。
   今の live 投影は、同じ描画 ID のすべての Result に当たっていた（`compilePlanBundle` の `resultIndexesByRenderedId`）ので、見えていた結果を保つ。
2. それ以外は K を `<h>[_record_<n>][__instance_<s>_<digest>]` と読む。record `n`（なければ唯一の record）の biological feature のうち、
   `stableFeatureId == h`（`s` があれば `sourceFeatureIndex == s` も）のものがちょうど 1 つなら、その identity にする。
3. どれにも当たらなければ行を捨て、件数を警告に出す。今の `app/session-feature-metadata.js::migrateFeatureOverrideState` と同じ書き方にする。

この 2 つの規則で、保存された Session の実際の形を覆える。

- crop も逆相補もない record では、描画 hash = 元座標の hash で、Python が行を照合していた。隠した feature は catalog にないが、規則 2 で当たる。
- crop した record（OV-01）と `__instance_` の複製（OV-05）では、Python は行を照合していなかった。feature は描かれて catalog にあるので、規則 1 で当たる。
  OV-01 の再現 B で Y に当たっていた行も、K は X の描画 ID なので X に戻る。
- v31〜33（catalog なし）は、今の回復（`buildSessionFeatureRecoveryPlan`）の後で規則 2 を使う。record は `renderRequest.records[record_idx].recordKey` から決める。
  record input がすべて `exactly_one` でないときは捨て、件数を警告に出す。

付け足すこと:

- label の map を 1 つでも移したときは `features.labelOverrideRows` を空にする。これは、その map から最後の Generate で作った表の写しである（map が変わると `app/watchers.js` が空にする不変条件と同じ）。
  map が空で表だけある Session（CLI で作ったものなど）は、今のまま表を送る。
- 必要な fixture（positive。main の writer で作る）:
  - `tests/fixtures/sessions/` に v44 を 2 つ。上の probe の手順を `origin/main` で再生して保存する。
    - Linear で 2 record、crop と逆相補、隠した feature と label の文字と Label Off を含むもの
    - Circular で同じ record を 2 回並べた canvas
  - v33 を 1 つ。`6b89c781` 以降で 33 を書いていた main の commit で作る。作れないときは v31〜33 を捨てる側にし、PR にそう書く。

### 4.4 版の数を上げるときに直すテスト（`tests/` 全体を grep した結果、#762 の前）

- request 8 → 9（PR-Q4-2）:
  - Python のテスト
    - `tests/test_documentation_reference_contracts.py`（集合と docs の文字列）
    - `tests/test_documentation_contracts.py`
    - `tests/test_session_request_codec.py`
    - `tests/test_web_packaging.py`（`BUNDLED_REQUEST_SCHEMAS`）
    - `tests/test_web_runtime_capabilities.py`
    - `tests/test_gui_interactive_capture_contracts.py` と `docs/capture/flows/how_to/interactive_sessions.py` の `CURRENT_RENDER_REQUEST_SCHEMA`
    - `tests/test_gallery_session_semantics.py`
    - `tests/run_losat_cache_browser_acceptance.py`
    - `tests/test_run_info_exact_replay.py`
    - `tests/test_joint_display_placement_surfaces.py`
  - JS のテスト（`tests/web/`）
    - `losat-session-schema-contract.test.mjs`（JS / Python / adapter の一致）
    - `session-request.test.mjs`
    - `gallery-session-publication.test.mjs`
    - `gallery-session-migration.test.mjs`
    - `session-cli-compatibility.test.mjs`
    - `joint-display-placement.test.mjs`
    - `run-info.test.mjs`
    - `imported-comparison-intent.test.mjs`
  - Playwright の spec
    - `interactive-svg-v3`
    - `joint-display-placement`
    - `contracts/session-regenerate-intent`
    - `contracts/active-result-edit-transaction`
  - コード側の直書き
    - `services/session-authority.js` の `schema === 8`
    - `session_request_codec.py` と `session-request.js` の `>= 8` / `< 8`
    - `services/gallery-session-publication.js`（`CURRENT_REQUEST_SCHEMA`、`ACCEPTED_REQUEST_SCHEMAS`）
    - `gbdraw/web_support/capabilities.py`
- Session 44 → 45（PR-Q4-4）:
  - Python のテスト
    - `tests/test_session_io.py`
    - `tests/test_composition_surface_contracts.py`
    - `tests/test_api_session.py`
    - `tests/test_session_compat.py`
    - `tests/test_record_planning.py`
  - JS のテスト（`tests/web/`）
    - `run-analysis-simple-path.test.mjs`
    - `alignment-reset-receipt.test.mjs`
    - `similarity-alignment-actions.test.mjs`
  - Playwright の spec
    - `composite-session-resources`
    - `settings-only-session`
    - `similarity-alignment-ui`
    - `session-active-mode`
    - `linear-typography`
    - `right-drawer`
    - `linear-multi-record`
    - `vibrio-session-save.performance`
    - `contracts/session-custom-skew-colors`
    - `contracts/vibrio-full-generation.serial.spec.js`
  - コード側の直書き
    - `services/session-authority.js` の `[42, 44]`、`=== 44`、`< 44`、`version === 44 ? 4 : 3`
    - `session_io.py` の `(41, 42, CURRENT)`、`(42, CURRENT)`、および catalog schema の判定（今の版を 44 と決め打ちしている）
    - catalog 4 → 5: `services/feature-catalog.js` と `gbdraw/web_support/feature_catalog.py` の `FEATURE_CATALOG_SCHEMA`、
      `session_io.py` の `CURRENT_FEATURE_CATALOG_SCHEMA`
    - `services/gallery-session-migration.js` の `version: 44`
    - `services/config.js` の版の判定
  - docs: `docs/SESSION_COMPATIBILITY.md` の版の表
- `tests/deploy*` はない。`playwright.gallery-publication.config.js` の spec には版の直書きがない。それでも、版を上げる PR は `rg -n "\b44\b|\b8\b" tests` の結果を確かめ、
  main だけ・deploy だけで走る suite も含めて期待値を同じ PR で直す（Owner の指示、2026-10-04）。
- Gallery の Session（10 件）は PR-Q4-4 の後で `tools/refresh_gallery_sessions.py` で作り直す（生成物）。

## 5. CLI と Python API

- 意味を変えないもの:
  - `--label_table` と `--feature_visibility_table` の規則の行。`hash` は描画した record の hash で、`location` と `record_location` も同じ。label 表の `hash` 行は今のまま On / Off を強制する。
  - `--annotation_table` の `feature_selector`。
  - `--feature_placement_table` の列。
- 広げるもの: 元座標の表（`--feature_placement_table` と新しい表）の `hash=` は、`~n` 付きの biological feature ID にも完全一致で当たる（PR-Q4-1）。
  Run Info の placement の行は `feature_selector=hash=<biologicalFeatureId>` と書く。
- Python API に足すもの:
  - `CircularDiagramOptions` / `LinearDiagramOptions.feature_overrides`（`FeatureOverride`）
  - annotation の `FeatureIdentitySpan`
  - `Diagram.feature_identity_notices`
  - `docs/REFERENCE/typed-requests.md` の Feature placement intent の節を、Feature identity overrides に広げる。
- Q2 = A のとき（推奨）は `--feature_override_table` を足す（Circular / Linear）。
  - 列: `record`（任意。record ID か `#index`）、`feature_selector`（必須。placement と同じ完全一致の selector）、`feature_visibility`、`label_visibility`、`label_text`。
  - API の入力は `feature_override_table` と `feature_override_table_file`。`feature_overrides` とは排他にする。
  - 符号化の前に、`resolve_identity_table_rows` で行を解決する（placement の表と同じ）。
  - Run Info は `featureOverrides` をこの表に投影する。投影できなければ、今のとおり `sourceRecipe.unavailableReason` にする（R12）。

## 6. Web

### 6.1 state と key

- `state.featureOverrides`: `{JSON.stringify([recordKey, biologicalFeatureId]): row}`。row は request の行に、Web だけの `labelSourceText`（今の `labelTextFeatureOverrideSources`。bulk の B6 に使う）を足したもの。
  この map が次を置き換える。
  - `featureVisibilityOverrides`
  - `labelVisibilityOverrides`
  - `labelTextFeatureOverrides`
  - `labelTextFeatureOverrideSources`
  - `featureVisibilitySelectorCacheOwner`（と `app/watchers.js::refreshFeatureVisibilitySelectorCache`）
- key の検査と request への投影は `services/feature-placement.js` が持つ（#762 の `canonicalFeaturePlacements` / `requestFeaturePlacements` を、行の検査関数を引数にして共有する）。ファイル名は変えてもよい。
- feature から key を作る処理は 1 つにする: `feature.record_key` と `feature.biological_feature_id`（catalog の取り込みで、描画したすべての feature に付く。`services/feature-catalog.js` の `validateAndProjectCatalogItem`）。
  catalog の取り込みは、`renderedTargetsByOverrideKey` と同じ走査で、`overrideKeyByRenderedId`（Result ごと）も作る（CW-02）。

### 6.2 live の投影（R1、R3）

- `app/candidate-render.js::compilePlanBundle`: visibility と label の操作を、色と同じく `renderedTargetsByOverrideKey.get(biologicalFeatureKey(row))` で作る。
  `resultIndexesByRenderedId` で引いていた 3 か所を消す。これで batch の別の Result でも、同じ identity に当たる（今は同じ描画 ID の文字列のときだけ当たる）。
- `visibility-actions.js::reconcileFeatureVisibility` / `effectiveFeatureVisibility`: `featureOverrides[keyOf(feature)]` を読む。
- `label-actions.js` の `projectLabelTextIntent` と `applyStoredVisibilityOverridesToSvg`: `data-label-feature-id` → `overrideKeyByRenderedId` → 行を読む。
- popup（`app/feature-editor/svg-actions.js::buildClickedFeaturePayload`、`label-actions.js::syncClickedFeatureLabelState`）と
  setter（`setFeatureVisibility`、`applyFeatureVisibilityScope`、`buildSelectedFeaturesVisibilityCommand`、`updateClickedFeatureLabelText`、
  `applyClickedFeatureVisibilityOverride`、`handleHiddenLabelTextChoice`、`handleLabelTextScopeChoice`）: 書く先は identity の行。
  "Keep hidden (apply text only)" は `labelVisibility: null` と `labelText` の行になる。
- `app/run-analysis.js` の `requiredLabelFeatureIds`: 選んだ Result で、`labelVisibility: on` かつ `present` の行の描画 ID。
  Q1 / Q2 のダイアログ（FINDINGS の PR-7）が、描けない On を Apply の前に止める。
- 古い label TSV の取り込み（`loadLabelOverrideTable`）は、行を今の Result の label で評価して文字を決める（今と同じ）。結果は identity の行に書く。

### 6.3 刈り取りと通知（R2）

- `pruneUnmatchedFeatureOverrides` は owner と呼び出し位置を変えない（`app/run-analysis.js` の、元の source を置き換えて成功した Generate の中）。入力を Python の `featureIdentityNotices` に変え、§3.4 の規則で行を消す。
  Q3 = A なら placement の行も同じ関数で扱い、`app/record-display-options.js` の source watcher による削除をなくす（R10）。
- `featureIdentityNotices` は `annotationWarnings` と同じく artifact が持つ値にする（`services/history-snapshot.js` の generated owner set、Session の `runMetadata`、`services/session-authority.js` の検査）。

### 6.4 History、Reset、書き出し

- History: feature intent（`buildFeatureIntentData` / `applyFeatureIntentData`）と artifact checkpoint（`buildFeatureStateData` / `applyFeatureStateData`）の 4 つの map を `featureOverrides` にする。`compactGeneratedArtifactSignature` も同じ。
- Reset: `services/reset.js::resetEditorDraftState` と `services/config.js::resetSessionBaseline`。
- `services/session-authority.js::ARTIFACT_FEATURE_FIELDS` と `gbdraw/session_io.py::CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS`: 古い 4 つの名前を 45 の writer では禁止する。
  Python は 45 の `features.featureOverrides` の key を、`_validate_display_placement_drafts` と同じ方法で確かめる。
- standalone SVG の書き出し（`services/svg-serialization.js::captureSvgExport`）: 埋め込む形は今のまま（label の feature ID を key にした文字）。
  書き出す Result の `renderedTargetsByOverrideKey` で、identity の行を描画 ID に投影する。`services/standalone-interactivity.js` は変えない。
- 書き出しと読み込み（Q2 = A）:
  - "Export feature edits TSV" と "Load feature edits TSV" を足す。読み込みは Worker の helper で `resolve_identity_table_rows` を呼ぶ（R4）。
  - Export Label TSV と Feature visibility の TSV は規則の行だけを書く。
  - annotation の TSV は `featureIdentity` の対象を、今の配置での `FeatureSpan`（`record=#<index>`、`feature_selector=hash=<描画 hash>`）として書き、その旨を通知する。
    annotation には今も CLI の lossless な recipe がない（`app/run-info.js` が unavailable にしている）。

### 6.5 消すもの（PR-Q4-4）

- `app/feature-visibility.js`
  - 死んでいる 7 つの export（§2）
  - `buildFeatureVisibilitySelectorCache`、`preserveFeatureVisibilitySelectorCacheForOverrides`、`featureVisibilityOverridesToRules`、`buildRuleFromSelectorCacheEntry`
  - 描画 ID を key にした set / get の helper（identity 版に置き換える）
  - `pruneUnmatchedFeatureOverrides` の hash による判定
- `app/feature-editor/label-override-table.js`: visibility の行と個別の文字の行を作る分岐、`selectStableFeatureKey`、`getRegexForFeature`。
- `app/feature-utils.js::getFeatureGenerationHash`: ほかに使う所が残らなければ消す。
- `featureSelectorSafetyScope` の一式（state、History の owner set、catalog の投影、Python の `_build_selector_safety_scope` と `selectorSafetyScopeBuildCount`）: selector cache が最後の消費者なら消す（CW-05）。`tests/web/computation-ownership.*` も合わせて直す。
- `gbdraw/web/CLAUDE.md`: R2 の owner の説明、R3 の投影、「Computation ownership」の owner、「Request and session boundary」を直す。

## 7. PR の分割（依存の順。どれも単独で merge でき、CI が通る）

| PR | 内容 | 先に失敗させるテスト | 消すもの |
|---|---|---|---|
| PR-Q4-1（#762 の後） | §3.1 の `FeatureIdentity` と `resolve_feature_identities`。placement をそれに載せる。元座標の表の `hash=` を `~n` まで完全一致にする。Run Info の placement の行を `hash=<id>` にする。永続化の形式は変えない | `tests/test_feature_placement.py`: 同じ座標の CDS 2 つのうち `hash=<h>~1` の行がその 1 つに解決する（今は "matched 0"）。`tests/test_run_info_exact_replay.py`: 重複した feature の placement の recipe を再生できる | `resolve_placement_inputs` の中の identity 解決と、表の行の解決の重複 |
| PR-Q4-2 | request schema 9、`featureOverrides`、`FeatureIdentitySpan`、§3.3 の消費者、`feature_identity_notices`、API の欄、JS の codec と promoter（Web は空の配列を書く）、docs | `tests/test_feature_identity_overrides.py`（新規。`tests/fixtures/web_batch_two_records.gb` を使い、audit の probe の `b_collide.gb` と `b_dup_ids.gb` を `tests/fixtures/` に入れる）: crop した record で X を Off にすると Y は残る、TX の文字は TX だけ、同じ record 2 回で 1 つ目だけ、`~1` の Label On はその instance だけ、文字だけの行は scope 外で label を出さない、identity の annotation は X に付く。`tests/test_session_request_codec.py`: schema 8 → 9 の昇格。`tests/web/session-request.test.mjs`: writer の欄 | なし（Web の表の経路は PR-Q4-4 まで残す。§4.2 の一時的な PE） |
| PR-Q4-3 | `--feature_override_table` と API の表の入力、Run Info の投影（Q2 = A のとき） | `tests/test_run_info_exact_replay.py`: `featureOverrides` を持つ request の recipe を CLI で再生すると、同じ SVG になる（今は欄がないので失敗する） | なし |
| PR-Q4-4（PR-7 の label ダイアログの後） | §6 の全部、Session 45 と migrator、catalog 5 と `ruleFeaturePayload`（§3.6） | Playwright `tests/web/feature-identity-overrides.playwright.spec.js`（新規）: OV-01 再現 B（Generate 後に Y と TY が変わらず、Session を CLI で描き直しても同じ）、OV-02（crop した record の `record_location` 色規則が、live でも Generate でも X を塗る）、OV-04（Generate 後も 1 つ目の複製だけ）、OV-05（重複の片方に Label On → Generate が成功し、その instance だけ）、OV-11（crop を 201 → 101、逆相補を外しても編集が残る）、OV-12（並べ替えの後、popup が Off を示す）。node: v44 の fixture の migrator（§4.3 の規則 1〜3）と、source 置き換えでの刈り取り（`tests/web/feature-visibility-reconciliation.test.mjs`） | §6.5、placement の source watcher、Web の個別編集の表の経路 |
| PR-Q4-5 | 選択から作る annotation を `featureIdentity` の対象にする。Editor では対象を "Selected feature: <caption>" と表示する。annotation の TSV は §6.4 の形 | Playwright: crop した record で X を選んで annotation にし、Generate すると X に付き、crop を変えても X に付いたまま（OV-03） | `featureTargetsFromSelection` の `selector.hash` |

- 分類: PR-Q4-1、PR-Q4-3、PR-Q4-5 は通常の変更（owner と経路の簡潔な証拠）。PR-Q4-2 と PR-Q4-4 は §4.2 の `ARCHITECTURE_EXCEPTION`。
- どの PR の後で main に昇格しても、形式は完結している（§4.2）。ただし PR-Q4-4 を分けたときは、全部が merge するまで昇格しない。
- PR-Q4-4 の検証では、ローカルで次を回す。
  - 変えた spec と、grep で当たる spec（`feature-visibility-*`、`feature-label-visual-unit`、`session-*`、`history-*`、`gui-audit-20260930-editor`、`multi-result-edit-matrix`、`non-edit-state-preservation`）
  - Circular と Linear の両方
  - 全部の functional の spec は PR の CI に任せる。

## 8. Owner に決めてほしいこと

### Q1: 複製した record での Feature visibility の編集は、複製ごとか

- **A（推奨）**: 複製ごと。identity にそのまま従い、live の表示（今も複製ごと）と Generate が一致する。label は PD-OI-069 の「instance 単位の編集」で、もともと複製ごと。
- B: すべての複製に当てる。編集の時に同じ biological feature ID を持つ複製の行も作る。表示と説明を変える必要がある。

### Q2: identity の行の表形式（CLI、Run Info、TSV の書き出し）

- **A（推奨）**: `--feature_override_table` を足す（§5）。
  - Run Info の recipe はそれで個別の編集を表せる（今は label 表と visibility 表で表していて、それより正確になる）。
  - Web の "Export / Load feature edits TSV" は同じ形式を使う。
  - Label TSV と visibility TSV は規則の行だけになる。
  - annotation の TSV は今の配置の形で書く（§6.4）。
- B: CLI の表は足さない。個別の編集がある間、Run Info の recipe は unavailable にする（R12）。書き出しは規則の行だけにする。
- C: 今の表に、今の配置の描画 hash で書く。同じ配置の CLI では使えるが、crop・逆相補・複製で外れる。R12 があるので Run Info には使えない。

### Q3: 解決できなくなった identity（元の source を置き換えた、ファイルを編集した）

- **A（推奨）**: §3.4 のとおり。Python は Generate を止めずに通知する。Web は、元の source を置き換えて成功した Generate で `unresolved` の行だけを消し、件数を知らせる。
  - placement にも同じ規則を当てる（PD-OI-065 の結果を visibility、label、placement でそろえる）。
  - 変わること: 置き換えた後も残る feature の placement は消えなくなる。CLI と API では、正確な placement の行の未知の identity がエラーから警告になる（`docs/REFERENCE/input-formats-and-tsv-schemas.md` を直す）。
- B: placement は今のまま（ファイルを置き換えた時点でその record の行を全部消し、Python は未知の identity をエラーにする）。A は visibility と label にだけ当てる。
