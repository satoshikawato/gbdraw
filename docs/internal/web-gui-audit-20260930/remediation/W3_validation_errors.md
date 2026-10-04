<!-- Raw design report of workstream W3 (sub-agent output, 2026-09-30, base origin/dev 4c89bab1). Evidence for 01_REMEDIATION_PROPOSAL.md, not Product or architecture authority. -->

# W3 設計提案: 入力検証、不正値の黙示変換、エラー正規化、Web↔CLI の recipe と表の互換

対象は DEV `4c89bab1` で、パスはリポジトリ直下から書く。コードの変更、commit、issue 作成はしていない。scratch に置いたのは使い捨ての probe だけ。「確認」は dev で実行または読んで確かめた事実、「推定」は未検証の判断を指す。

## 0. 設計の前提として新たに確認した事実

**JS 側のエラー分類**
- `normalizeUserFacingError` の動き（`gbdraw/web/js/services/error-normalization.js:387-425`）:
  - `code ∈ DEFINITIONS` の object が来たら、その `code/stage/context` をそのまま使う（`:392`）。
  - それ以外の Error は `message` を逆引きする。使う表は `NATIVE_VALIDATIONS`（完全一致 102 件、`:131-263`）と約 20 本の正規表現（`:264-335`）。
- 構造化 error の受け口はすでにある。`app/legend-layout/decoration-continuity.js:8-14` と `services/session-import-client.js:6-8` が使っている。
- それでも Web JS 全体の `throw new Error(` は 743 箇所あり、利用者向けの検証でも大半が plain Error。逆引き表に無い文言は UNKNOWN になる。

**Python 側のエラー分類**
- Python → Worker → JS の境界はすでに構造化されている（`gbdraw/web_support/error_adapter.py:214-264`、`gbdraw/web/js/workers/diagram-generation-worker.js:156,494`）。
- 問題は Python 内部の分類が英文の照合であること（`error_adapter.py:42-127` の `_CONSTRAINTS/_EXACT/_TEMPLATES`）。
- `_web_error_field`（`:165-167`）は隠れ hook で、設定しているのは `gbdraw/web/js/app/python-helpers.js:1722,1728` の 2 箇所だけ。

**Session import**
- Worker はすでに `SESSION_IMPORT_{READ,PARSE,REPLY}_FAILED` と `SESSION_IMPORT_CRASH` を送っている（`workers/session-import-worker.js:36-46`、`session-import-client.js:37-43`）。
- しかしこれらの code は DEFINITIONS に無いので捨てられ、UNKNOWN になる。
- 例外は PARSE で、`services/config.js:4238-4241` が別の code に書き換えている。

**summary に位置情報が出ない**
- context には `inputOrdinal`、`row`、`slotIndex`、`seriesIndex` を保持している。
- それでも summary には出さない（`error-normalization.js:406-416`。例外は DECORATION の Result N だけ）。
- つまり「何番目の Sequence か」「何行目か」が画面から消えている。

**エラーの包み直しで code が失われる**
- `config.js:4028-4030`、`:4073-4075`、`:4190-4193` は、原因の error を code の無い plain message に包み直している。

**UNKNOWN の action**
- UNKNOWN の actions に `save-session` が入っている。
- そのため Save 自体が失敗したときも、エラーパネルに Save Session ボタンがもう一度出る（`gbdraw/web/index.html:4710`）。
- `generate` action を持つ code（EXPORT_INPUT など）には、対応するボタンが無い（`index.html:4704-4718`）。

**CLI の実測（dev、HmmtDNA）**
- `-w 0`、`-w -100`、`-s 0`: exit 0 のまま GC と skew の track が消える（SVG 55,257 → 33,385 byte）。
- `-n G`: IndexError の traceback（`gbdraw/analysis/skew.py:85-87`）。
- `-n XY`: 成功し、平坦な track になる。
- `--block_stroke_width -1`: ValueError の traceback（`gbdraw/configurators/legend.py:73-75`）。
- `--label_font_size -5`: 成功する。
- 比較の閾値（`1e-50x`、`nan`、identity 150/−5、bitscore −1、alignment_length 12.5/−10）は、CLI がすべて拒否する。
- `gbdraw/cli.py:176-178` が整形するのは `GbdrawError` だけ。ValueError と IndexError は traceback のまま出る。

**Generate が draft を書き換える**
- `gbdraw/web/js/app/run-analysis.js` の Generate 経路には、`adv.*=`、`form.*=`、`circularConservation.*=` の代入が 80 箇所ある。
- 例: `:2715-2731`、`:2996-3010`、`:2598`、`:2626-2639`。

**Python → JS の生成データ経路**
- 経路: `gbdraw/mode_profiles.py:175-198`（`mode_profiles_payload`）→ `tools/generate_mode_profiles.py` → `gbdraw/web/js/mode-profiles.generated.js`。
- `tests/test_web_mode_profiles.py` が `--check` で同期を確かめている。
- 比較の閾値の既定値は、すでにこの経路を通っている。

**リリース版（0.13.0）の挙動**
- Web は argv で `--definition_font_size` を CLI に渡していた（`git show 0.13.0:gbdraw/web/js/app/run-analysis.js:1905`）。したがって定義の間隔は `int(font+2)` だった。
- Gallery の全 session はこの派生規則と整合している（HmmtDNA_ATskew と tobacco は font 28 / interval 30、それ以外は font 18）。
- 0.13.0 が書く Session は v30 で、当時の Save に feature catalog の要件は無かった。

---

## X-01（P2）エラー変換の抜け（IN-07、GE-05、TR-05、SE-09、および SE-05/SE-06/PV-04/CO-01/IN-04 の文言）

### 1. 分類
- **IMPLEMENT_EXISTING_AUTHORITY**。根拠は次の 2 つ。
  - PD-OI-046 の Must preserve「すべての移行対象で既知 validation の修正情報を保持する」。
  - issue-601 の MASTER_PLAN の「producer は失敗の意味を所有し、利用者文言を持たない」と「message.includes… を新しい classifier にしない」。
- 新しい Product Decision は不要。
- ただし任意の Phase 4「Show field」ボタンは新しい affordance なので、実装前に Owner の確認が要る。

### 2. 根本原因の再確認
監査の指摘は正しいが、原因は「表が足りない」より一段深い。

- **(a) 意味の owner が 2 箇所ある（OE の余剰 1）。**
  - 利用者向けの失敗の意味を、producer の英文と normalizer の逆引き表の両方で決めている。
  - 文言を 1 字変えると、何も言わずに UNKNOWN になる。
  - 同じ英文が 3 箇所に複製されている例もある（`app/app-setup.js:461`、`services/session-request.js:1097`、`error-normalization.js:153`）。
- **(b) 構造化された code でも語彙に無いと捨てる。** session import が該当する。
- **(c) Python も英文の句で照合している。**
  - `identity must be a finite number in [0, 100].` や `alignment_length must be >= 0 and an integer.` は句が表に無い。
  - 結果として field の無い VALIDATION_UNCLASSIFIED になる。
- **(d) 位置情報が summary に出ない。** 上記 §0 のとおり。
- **(e) wrapper が原因の code を捨てる。** 上記 §0 のとおり。
- **(f) `REASONS.READ` の文が重複する。**
  - `REASONS.READ`（`error-normalization.js:83`）は `INPUT_UNREADABLE` の定義文とほぼ同文。
  - Python の template（`error_adapter.py:122`）が READ を付けるので、案内が 2 回出る。

**監査の修正点**
- fit error で slot id を出さないこと、depth のエラーでファイル名を出さないことは PD-OI-046 のとおりで、正しい。
- 足りないのは、代わりに示すべき slotIndex、使える帯の半径、seriesIndex である。

### 3. 修正案（目標設計）

**責任の分担（SRP）**
- 値が有効かの判断: 各 domain の validator が持つ。
  - JS: `app/current-option-values.js`、`app/track-slot-validation.js`、request の projection。
  - Python: typed options、config、codec。
- 表示: `services/error-normalization.js` だけが持つ。

**JS の producer 契約**
- `error-normalization.js` に約 8 行の export を 1 つ足す。
  ```js
  export const diagnosticError = (code, context = {}, { stage = 'request-validation', operation } = {}) =>
    Object.assign(new Error(`${code}/${context.reason || ''}`), { code, stage, ...(operation && { operation }), context });
  ```
- message は固定の識別子にして、record 名などを補間しない。これで未捕捉時の console にも私的な値が出ない。
- 消費側は既存の `CODES.has(object.code)` 経路（`:392`）と `contextFor`（`:372-384`）をそのまま使うので、変更は要らない。

**Python の producer 契約**
- `gbdraw/exceptions.py` の `GbdrawError` に kw-only の `diagnostic: Mapping | None` を足す。中身は `{"code","reason","field",` 整数 context`}`。
- `serialize_web_error` は、原因の chain 上に `diagnostic` があればそれを最優先で使う。
  - code、reason、field、context key は allowlist で検証する。
  - allowlist は JS と同じ語彙にし、parity test で縛る。
- `diagnostic` が無いときだけ、従来の英文表へ fallback する。
- `_web_error_field` は同じ PR で `diagnostic` に置き換えて削除する。

**Worker → JS の境界**
- diagram Worker と session Worker は変更しない。
- LOSAT worker（`workers/losat-worker.js:169`）は文字列を送っている。利用者向けの 1 件（CO-01）は、呼び出し側の `services/losat.js:793` で typed にする。

**表示の変更**
- summary の locator を data 表で表示する。
  - `inputOrdinal`: code に応じて 'Sequence N'、'Result N'、'Comparison FASTA N'。
  - `row`: 'Line N'。
  - `slotIndex`: 'Track row N'（1 始まりで表示）。
  - `seriesIndex`: 'Depth series N'（1 始まりで表示）。
  - 帯域: 'Available band: 181–209 px.'（CANNOT_FIT のとき）。
- `decorationResult` の特例（`:414-416`）はこの表に吸収して削除する。
- `FIELD_LABELS`（`:126-128`）に data 行を足す（window → 'Window'、evalue → 'E-value'、feature_width_circular → 'Feature Track Width' など）。
- `index.html:4704-4718` にボタンを 1 行足す: `errorDisplay.actions.includes('generate') && errorDisplay.operation !== 'generate'` で Generate を出す。

**UI と field の結合**
- 結合のキーは `errorDisplay.context.field`（有界の識別子）。新しい state は作らない。
- 任意の Phase 4 では、`app-setup.js:3569` の `errorDisplay` から `errorField` computed を 1 つ作り、対象の input に `data-error-field` と `:aria-invalid` を付ける。focus ボタンは Owner の確認後にする。
- 既存の inline error（comparison height、`app/circular-track-slots/measure-input.js`）も、この語彙から文言を取る。上記 (a) の 3 重の英文複製はこれで 1 つになる。

**語彙の追加（すべて data 行）**
- DEFINITIONS:
  - `SESSION_SAVE_REQUIRES_GENERATE`: "Save Session needs current feature metadata for the loaded Result. Generate once, then Save Session."、actions `['generate']`。
  - `THREADED_LOSAT_UNAVAILABLE`: "Threaded LOSAT needs a cross-origin isolated page. Choose Serial LOSAT execution, or open gbdraw.app or `gbdraw gui`."、actions `['edit-comparison','retry']`。
- REASONS: 追加は `SELECT_RECORD_FOR_REGION`、`DISCOVERY_PENDING`、`SESSION_FORMAT`、`LEGEND_NAME_CONFLICT`、`PRESERVED_RECORDS`、`POSITIVE_INTEGER_OR_AUTO`、`DINUCLEOTIDE`、`SPECIFIC_COLOR`、`DEPTH_VALUES`。`READ` は削除する。
- OPERATIONS（JS だけ）: `session-save` と `session-load`。`operationErrorTitle` にも対応する見出しを足す。
- STAGES: `parse` を足す。
- context key: `innerPx`、`outerPx`、`configPath`。
  - `configPath` は canonical な config の path に限る。Python では `canonical_config_override_paths()` の集合、JS では path 形の pattern かつ 80 字以内で検証する。

**producer ごとの対応**

| 出典 | 現状 | 新しい typed |
|---|---|---|
| `services/config.js:2883` | UNKNOWN | REGION_INVALID {inputOrdinal, SELECT_RECORD_FOR_REGION} |
| `config.js:2875` / `:2861, 2865` | UNKNOWN | NO_RECORDS {inputOrdinal} / RECORD_SELECTION {DISCOVERY_PENDING} |
| `app/annotations/record-catalog.js:100, 105-106` の issue 文字列 | UNKNOWN | issues を文字列から `{code, context}` の object に変える（INPUT_REQUIRED {inputOrdinal} など） |
| `record-catalog.js:60-72` と `app/annotations/validation.js:34-60` | 正規表現 7 本（`error-normalization.js:301-309`）に頼っている | ANNOTATION_TARGET {reason} を返す。文字列の prefix 連結をやめる |
| `session-request.js:979-981`（"(automatic)"） | UNKNOWN | RECORD_SELECTION {SELECT_ONE} |
| `session-request.js:1036-1039` | UNKNOWN | RECORD_SELECTION {NO_MATCH, AMBIGUOUS} |
| `session-request.js:641` | UNKNOWN（IN-08） | INPUT_REQUIRED |
| `session-request.js:876` | UNKNOWN | REGION_INVALID {BOTH_ENDPOINTS} |
| `app/circular-track-slots.js:771-777` | UNKNOWN | INPUT_INVALID {field: feature_width_circular など, POSITIVE_OR_AUTO}。ラベル表 `:753-760` を field 表に置き換える |
| `config.js:1189, 1588` | UNKNOWN | INPUT_INVALID {field: schema, SESSION_FORMAT} |
| Session import の READ / PARSE / CRASH / REPLY | UNKNOWN、または書き換え | Worker が最終的な code を送る（READ → INPUT_UNREADABLE {schema}、PARSE → INPUT_INVALID {schema, JSON_FORMAT}、CRASH/REPLY → SESSION_IMPORT_UNAVAILABLE）。`config.js:4238-4241` は削除 |
| `config.js:4019, 4030`（SE-05） | UNKNOWN | SESSION_SAVE_REQUIRES_GENERATE、operation 'session-save' |
| `config.js:4073-4075, 4190-4193` | 原因の code を捨てる | 原因が既知 code の model ならそれを再 throw。そうでなければ INPUT_INVALID {schema} |
| `app/legend/entry-actions.js:735`（PV-04） | UNKNOWN | INPUT_INVALID {field: legend, LEGEND_NAME_CONFLICT} |
| `services/losat.js:793`（CO-01） | UNKNOWN | THREADED_LOSAT_UNAVAILABLE |
| `services/imported-comparison-intent.js:441, 446`（SE-06） | UNKNOWN | COMPARISON_INPUT {PRESERVED_RECORDS} |
| `gbdraw/mode_profiles.py:43, 49, 52` | UNCLASSIFIED | diagnostic {identity, PERCENT} / {alignment_length, NONNEGATIVE_INTEGER} |
| `gbdraw/diagrams/circular/radial_layout.py:934-938, 1193-1198` | reason だけ | {TRACK_LAYOUT, CANNOT_FIT, slotIndex=`intent.slot_index`, innerPx, outerPx}。slot id は出さない |
| `gbdraw/analysis/depth.py:92` | INPUT_UNREADABLE + READ（同じ文が 2 回） | DEPTH_INVALID {DEPTH_VALUES, seriesIndex}。seriesIndex は depth track の読み込みループで付ける |
| `gbdraw/config/modify.py:280-285` | UNCLASSIFIED（例: scale interval 1.5） | INPUT_INVALID {INTEGER, FINITE など, configPath}。表示は 'Setting: objects.scale.interval.' |

**同じ PR で削除する旧経路**
- 移行した producer に対応する `NATIVE_VALIDATIONS` の行と正規表現。
  - `Sequence #` の正規表現（`error-normalization.js:314-321`）は、`app/feature-metadata-extraction.js:238-244` を移行する PR で同時に消す。
- `config.js:4238-4241` の書き換え。
- decoration-continuity の手組みの error と、session-import-client の `importError`（どちらも `diagnosticError` に置き換える）。
- producer が存在しない dead 行 'An exact reference feature is required for alignment.'。
- `REASONS.READ`。
- `_web_error_field`。

**追加しないもの**
- 汎用の validation framework、Error のサブクラス階層。
- raw message の表示。
- Python の範囲規則を JS に手書きで複製すること。
- 新しい session field、新しい Worker。

### 4. 上位の対応
- `gbdraw/web/CLAUDE.md` の Debugging principles に 2 行を足す。
  - 「利用者が直せる失敗は、JS では `diagnosticError`、Python では `diagnostic=` を付けて、code と有界 context で投げる。」
  - 「message を照合する分類器は追加しない。」
- 強制の方法:
  - **(a) 表駆動 test。** 実際の producer を呼び、返る code と context を assert する。既存の方式（`tests/web/error-normalization.test.mjs:106-137`）に合わせる。
  - **(b) ratchet test。** 次の件数を baseline 以下に保つ。
    - `NATIVE_VALIDATIONS.size`。
    - `nativeValidation` 内の正規表現の本数。
    - Python の `_EXACT`、`_TEMPLATES`、`_CONSTRAINTS` の件数。
    - 増加は禁止し、減ったら baseline を下げる。
  - **(c) 語彙の parity test。** 既存の FIELDS の test（`error-normalization.test.mjs:151-160`）を CODES、REASONS、context key に広げる。
  - **(d) privacy の test。** `diagnosticError(` の第 1 引数が文字列リテラルであることを確かめる。
- 文言を文字列として走査する test は採らない。組み立てた文言（16 件）を検出できないため。

### 5. テスト
- **node**（`tests/web/error-normalization.test.mjs`）: 上の表の producer を実際に呼ぶ。
  - 例: `materializeLinearRecordFiles` に 2 record の catalog と region を渡す → REGION_INVALID {inputOrdinal: 1}。
  - 例: `normalizeCircularGeometryShortcuts({featureWidth:0})` → INPUT_INVALID {feature_width_circular}。
  - summary に 'Sequence 1'、'Line 3'、'Track row 2' が出ること。
  - 'Replace or reselect' が 1 回だけ出ること。
- **pytest**（`tests/test_web_error_adapter.py`）:
  - `diagnostic` が優先されること。
  - allowlist に無い key が捨てられること。
  - fit error が slotIndex、innerPx、outerPx を持つこと。
  - identity 150 → PERCENT。
- **Playwright**（`tests/web/error-boundary.playwright.spec.js` を拡張）:
  - Feature Width 0、truncated .gz、v39 Session の Save（Generate ボタンが出ること）、古くなった Circular selector。
  - すべてで code ≠ UNKNOWN、かつ旧 Result が保持されること。
- 今回のバグを捕まえる assertion: 「利用者向け validation の producer の全行が、UNKNOWN も VALIDATION_UNCLASSIFIED も返さない」。

### 6. 規模・PR
- **PR-1（Ordinary、8 files）**
  - 対象: `error-normalization.js`、`config.js`、`session-request.js`、`circular-track-slots.js`、`decoration-continuity.js`、`session-import-client.js`、`workers/session-import-worker.js`、`index.html`。
  - 見積: churn 250〜350、net +60〜90。
  - 直るもの: IN-07 の主要部、GE-05 の JS 部、SE-05 の文言、SE-09、IN-08 の文言。
- **PR-2（Ordinary、net はマイナス）**
  - 対象: `record-catalog.js`、`annotations/validation.js`、`config.js`、`app-setup.js:2827,2843`、`feature-metadata-extraction.js`、`error-normalization.js`（正規表現の削除）。
- **PR-3（Python が主。Web は 2 files）**
  - 対象: exceptions、adapter、mode_profiles、radial_layout、depth、config/modify、および `python-helpers.js` と `error-normalization.js`。
- **PR-4 以降（ratchet の系列、各 Ordinary、net はマイナス）**
  - 対象: `current-option-values.js`（`:24-30`、`:95-112`）、`track-slot-validation.js`。
  - `track-slot-validation.js` では geometry_invalid の英文を再解析している箇所（`error-normalization.js:359-366`、`:284-289`）を issue code にする。
  - 続いて `services/export.js`、`app/similarity-alignment.js`、session binding の validator。
  - 最後の PR で `nativeValidation` と `NATIVE_VALIDATIONS` を削除する。旧経路の削除なので `architecture-change` label を付けてよい。
- 他 workstream の IN-04、SE-06、PV-04、CO-01 は、それぞれの挙動修正 PR の中で上の code を使う。別 PR は作らない。

### 7. 依存・リスク
- PR-1 は X-02、FE-12、他 workstream の前提になる。
- POSITIVE_INTEGER_OR_AUTO を導入すると、既存 test の期待値が変わる（例: candidate limit の reason）。
- catalog の issues を文字列から object に変えるので、表示している箇所を確認する。文言は `normalizeUserFacingError(issue).summary` で得る。
- `GbdrawError.__init__` の拡張について、pickle を使う経路が無いことを確認する（現状は無いと推定）。
- ratchet の文書では OE を 2 → 1（意味の owner）、PE を 2 → 1（message 経路の廃止）と記録する。

---

## X-02（P2）不正な数値の黙示変換（GE-04、CO-09、GE-08、TR-04）

### 1. 分類
- **IMPLEMENT_EXISTING_AUTHORITY**。根拠:
  - OIPC-C01、C03、C04、C06。
  - `gbdraw/web/CLAUDE.md` の "Changing a setting" の手順 3 と 4、および "Python owns request validation"。
- CLI の `-w 0`、`-s 0`、`-n G`/`XY` の拒否も OIPC-C01 による。無効な値で空または平坦な track を黙って出しているため。
- dinucleotide の文字集合は Owner の確認を推奨するが、ブロッキングではない。推奨は大小無視の `^[ACGT]{2}$`。N や U を許す余地は残る。

### 2. 根本原因の再確認
- **(a) JS に coercer が複数ある。**
  - `run-analysis.js:885-895` は既定値に置き換え、その値を draft に書き戻す。
  - `session-request.js:554-564` は不正値を null、つまり Auto にする。
  - `session-request.js:2463-2466` の `Number()` は NaN を作る。
    - JSON 化で NaN は null になり、Python は既定値を使う。
    - 一方、メモリ上の request には NaN が残り、それを Run Info が読んで `--evalue NaN` を出す。
- **(b) Python の typed options が値を検証しない。**
  - `gbdraw/api/options.py:641-645` は window、step、depth_window、depth_step、dinucleotide を検証しない。
  - 既存の `_validate_positive_int`（`:300-311`）があるのに使っていない。そのため CLI も誤る。
- **(c) dinucleotide の扱いが入口で食い違う。**
  - slot の `nt` は短いと黙って既定値になる（`gbdraw/diagrams/circular/assemble.py:1118-1122`）。
  - トップレベルの `-n` は IndexError になる。
- **(d) 構造上、Generate が draft を 80 箇所で書き換える。**
- 監査の「CLI rejects most of these」は、window/step の 0 と負値については誤り。

### 3. 修正案

**誰が検証するか**
- 範囲、整数であること、文字集合は Python の typed layer が唯一の owner になる。
- JS は「JSON で表現できるか」だけを見る。
  - 空欄と 'auto' → null。
  - 有限の数 → そのまま送る（0、−100、500.5、150 も送る）。
  - 数に読めない文字列 → `diagnosticError('INPUT_INVALID', {field, reason:'FINITE'})`。
  - JSON が NaN を運べないので、ここだけは JS に残す。
- 例外として比較の閾値がある。閾値は Python の render より前に、LOSAT の後処理、cache key、protein 変換 helper（`run-analysis.js:3080-3083`、`:4186-4189`、`:4246-4249`）で使われる。
  - そのため Python の定義域を data として生成物に載せ、JS はその data を評価する。手で複製はしない。

**PR-A（Python）**
- `gbdraw/mode_profiles.py`:
  - 定義域の表 `COMPARISON_THRESHOLD_DOMAINS` を置く: evalue {min 0}、bitscore {min 0}、identity {0..100}、alignment_length {整数, min 0}。
  - `:25-52` の 3 つの validator を、この表を読む 1 関数に置き換えて削除する。検証失敗は `diagnostic=` 付きで raise する。
  - `mode_profiles_payload()` に `comparisonDomains` を足し、`mode-profiles.generated.js` を再生成する。
  - 既定値は変わらないので `MODE_PROFILE_VERSION` は据え置く。
  - 同じファイルに `validate_dinucleotide(value, field)` を置く（2 文字で ACGT）。
- `gbdraw/api/options.py` の `_ModeDiagramOptions.__post_init__`:
  - 既存の loop の形（`:1080-1101`）で window、step、depth_window、depth_step を `_validate_positive_int(allow_none=True)` にかける。
  - `plot_title_font_size` を `_validate_positive_real` に、dinucleotide を上の validator にかける。
  - CLI、Python API、Web の typed request がすべてこの 1 箇所を通る。
- slot の `nt`: `gbdraw/tracks/circular.py` と `gbdraw/tracks/linear.py:187-188, 374-377` で同じ validator を使い、`assemble.py:1118-1122` の黙った fallback を削除する。
- CLI adapter の重複した検査（`gbdraw/circular.py:512-513`、`gbdraw/linear.py:793-794`）を削除する。
- `error_adapter.py:46` の「must be a positive integer or None」を POSITIVE_INTEGER_OR_AUTO に変える。今は案内から「整数」が落ちている。

**PR-B（Web）**
- `utils/optional-positive-number.js`:
  - 既存の `DECIMAL_NUMBER_PATTERN` を再利用した `classifyOptionalNumber`（'auto' にも対応）を足す。
  - 既存の `classifyOptionalPositiveNumber` はその上に作る。
  - `projectOptionalNumber(value, field)` は分類して、不正なら `diagnosticError` を投げる。
- `session-request.js`:
  - `optionalNumber` と `optionalPositiveInteger` を削除する。
  - window、step、depth の window/step（`:2455-2458`）と plotTitleFontSize を `projectOptionalNumber` に置き換える。
  - config leaf（`:1131-1248`）は `configPath` context 付きで同じ helper を使う。
  - 比較の閾値（`:2463-2466`）は resolver の結果を使う。
  - `canonicalOptionalPositiveNumber`（`:575-582`）も統合する。
- `mode-profiles.js`: `resolveComparisonThresholds(adv, mode)` を置く。
  - 空欄はその mode の既定値にする。
  - 生成された定義域で検証し、違反は typed error にする。
  - evalue は trim した文字列も返して、cache key の `String(evalue)` の形を変えない。
- `run-analysis.js`:
  - `normalizeBlastThreshold*`（`:885-895`）と draft への書き戻し（`:2715-2731`、`:2996-3010`）を削除する。
  - Generate の冒頭で resolver を 1 回呼び、その snapshot を下流のすべてに渡す。
  - depth の window/step の ad hoc 検査（`:2371-2384` の該当部）を削除し、対応する `NATIVE_VALIDATIONS` の行（`error-normalization.js:142-143`）も削除する。
- `app-setup.js:455-462`: 分類器と語彙を使うように変える。
- `index.html:4526, 4532`: depth の input（`:1532`、`:1538`）に揃えて `min="1" step="1"` を付ける。

**追加しないもの**
- JS 側の汎用 schema や validator の登録簿。
- JS 側の window/step の範囲検査。
- 既定値への黙った置き換えと、draft への書き戻し。

**UI の結果**
- 入力欄の値はそのまま残る。
- summary は例えば 'An input value is invalid. Field: Window. Use Auto or an integer greater than zero.' になる。
- Phase 4 の `errorField` で入力欄を aria-invalid にできる（任意）。

### 4. 上位の対応
- `gbdraw/web/CLAUDE.md` の "Changing a setting" の手順 3 に 1 文を足す: 「projection は値を変換しない（空欄 → null、数 → 数、それ以外は typed INPUT_INVALID）。Generate は draft を書き換えない」。
- ratchet を 2 つ置く。
  - `run-analysis.js` の Generate 経路での draft 代入の件数（現在 80）。増加は禁止する。
  - 下記の parity vector。

### 5. テスト
- 共有の vector を 1 つ作る: `tests/fixtures/option_domain_vectors.json`。項目は field、CLI flag、mode、invalid の値、valid の値、期待する reason。
- **pytest（新規 `tests/test_option_domain_parity.py`）**: 同じ invalid 値について、次の 3 経路で同じ field と reason になることを確かめる。
  - (i) CLI: subprocess で実行し、非 0 で終了し、traceback が無いこと。
  - (ii) codec での typed request の decode: `ValidationError.diagnostic` を見る。
  - (iii) `serialize_web_error`。
  - 値は window 0/−100/500.5、step 0/10.5、depth_window 1.5、dinucleotide G/XY/GCA、evalue nan、identity 150/−5、bitscore −1、alignment_length 12.5/−10。
- **node**（`session-request.test.mjs`、`mode-profiles.test.mjs`、`optional-positive-number.test.mjs`）:
  - 数値が変わらずに request へ渡ること。
  - 数でない文字列は FINITE になること。
  - resolver の reason が Python と一致すること。
  - Generate の後に `adv.*` が変わらないこと。
  - recipe に NaN が出ないこと。
- **Playwright 1 本**: Circular で Window 0 → 失敗し、field 名が表示され、入力欄は 0 のまま、旧 Result は保持される。
- 今回のバグを捕まえる assertion:
  - 同じ不正値が CLI、typed request、Web で同じ field と reason で拒否されること。
  - Generate が draft を書き換えないこと。

### 6. 規模・PR
- **PR-A**: Python が主。Web 側は生成物 1 file だけ。
- **PR-B**: Web Ordinary、約 7 files、churn 300〜450、net はマイナス。
- **PR-C（後続、1〜2 PR）**: 同じ族の残りを同じ方式で直す。
  - multi-record の ratio（`run-analysis.js:873-883`）。
  - plot title（`:2598`、`:2635`）。
  - label spacing（`:2638`、`:3163`）。
  - ring の width/gap（`:2709-2714`）。
  - gencode（`:2740-2746`）。
  - depth_height（`:4535`）。

### 7. 依存・リスク
- X-01 の PR-1 が前提。
- **PR-A を先に入れること。** 逆順だと Web が 0 を送り、Python が受理して空の track を描いてしまう。
- CLI の挙動と文言が変わるので、release notes に書く。
- LOSAT の派生 cache key は値を変えないので、cache は無効にならない。
- 残るリスク: `type=number` の badInput（例 '1e'）は browser が '' を返すため、既定値扱いのまま残る。

---

## SE-05（P2）旧形式の Session を読むと Save が失敗する

### 1. 分類
- **文言: IMPLEMENT_EXISTING_AUTHORITY**（PD-OI-046）。
- **挙動: PRODUCT_DECISION_REQUIRED**。
  - 0.13.0 は v30 を書き、Load → Save が catalog 無しで成立していた。
  - dev は v<40（schema-v2 の fixture は v33）を読んだ後の Save を拒否している（`config.js:4016-4019`）。
  - つまり、リリース済みの継続操作を決定記録なしに退役させている。
- 選択肢:
  - **A / REQUIRE_GENERATE_BEFORE_SAVE（推奨）**
    - v<40、または catalog の無い artifact Session は、1 回 Generate してから Save する。
    - Load 時に notice を出し、Save の失敗は SESSION_SAVE_REQUIRES_GENERATE（Generate ボタン付き）にする。
    - `docs/SESSION_COMPATIBILITY.md` に明記する。format は変えない。
  - **B / SAVE_TIME_CATALOG_RECOVERY**
    - Save 時に、既存の `buildSessionFeatureRecoveryPlan`（`app/session-feature-metadata.js:733`）で旧 Result に整列した catalog を作り、現行 format で保存する。
    - 回復結果が `validateFeatureCatalog`（`config.js:4023-4027`）に受け入れられるかは未確認なので、EVIDENCE_REQUIRED。
  - **C / Load 時に catalog を作る**
    - PD-OI-044 の「preview-only Load で Python を初期化しない」との関係を確かめる必要がある。
    - 監査の v39 Load に 13.7 s かかっており、旧形式の Load ではすでに Worker を使っている可能性がある。
    - B と同じ evidence が要る。
  - **自動 Generate する、または Result を捨てて保存する案は NOT_ALLOWED**。保存が artifact を変えてしまうため（OIPC-C07、PD-OI-045）。
- A を推奨する理由:
  - B と C は、旧 renderer の SVG ID と現行 catalog の対応を heuristics で決めることになる。誤った feature identity を現行 format に固定してしまう危険がある。
  - A は既存の test（v39 は Generate してから Save する）とも一致する。

### 2. 根本原因の再確認
- 監査のとおり。旧形式の Load では `adoptCatalog: currentSchemaSession`（`config.js:4646`）が false なので catalog が null のまま残り、`:4018-4019` で Save が止まる。
- 加えて、エラーパネルが失敗したばかりの Save ボタンをもう一度出す（§0 参照）。

### 3. 修正案（A の場合）
- `config.js:231-232` を削除し、`:4019` と `:4030` を `diagnosticError('SESSION_SAVE_REQUIRES_GENERATE', {}, {operation:'session-save', stage:'result-admission'})` にする。これは X-01 の PR-1 に含める。
- `:4835` で operation を渡す。
- 旧形式の Session に Result があれば、Load の notice に 1 行出す。出す場所は `applySessionFeatureRecoveryPlan`（`:3960-3993`）の状態表示。
- docs に 1 段落足す。

### 4. 上位の対応
- A を採る場合は PD-OI の record を追加する。B を採る場合は evidence の後に別の record にする。

### 5. テスト
- Playwright で v39 と schema-v2 の fixture（`tests/fixtures/sessions/BGC0000708-BGC0000713.v39.gbdraw-session.json.gz`）を使う。
  - Save → SESSION_SAVE_REQUIRES_GENERATE が出て、Generate ボタンが表示される。
  - Generate → Save が成功し、再 Load した内容が一致する。
- v40–v44 の Save が成功することは維持する。

### 6. 規模・PR
- 文言は X-01 の PR-1 に含める。
- A の docs と record は authority の PR として別に出す。

### 7. 依存・リスク
- 文言の修正は判断を待たずに先行できる。

---

## GE-03（P2）Run Info の Source recipe で図を再現できない

### 1. 分類
- **IMPLEMENT_EXISTING_AUTHORITY**。根拠:
  - `docs/REFERENCE/web-app.md:814-818`「unavailable Source recipe includes a reason; it does not silently omit unsupported settings」。
  - `docs/REFERENCE/session-and-request-compatibility.md:248-253`。
- **Circular**: 0.13.0 の挙動の復元になる。今の dev の「間隔 20 固定」を残すなら、0.13.0 の挙動を退役させることになるので、その場合は PRODUCT_DECISION が要る。
- **Linear**: unlinked で Auto の意味は PD-OI-011 で決まっている。したがって recipe を unavailable にするのが正しい。
- CLI に `--ruler_label_font_size auto` を足す案は新しい capability で、別の判断になる（不要と判断）。

### 2. 根本原因の再確認（監査の修正）
- **Circular**: recipe の欠陥ではなく Web の退行である。
  - typed request 化（`session-request.js:1163-1165`）のとき、control の結合規則が落ちた。
  - 結合規則は CLI（`gbdraw/circular.py:951-954`）と legacy の flat key（`gbdraw/session_request_codec.py:1690-1709`）にはまだある。
  - recipe の guard（`run-info.js:965-973`）が働くのは、interval が明示されている場合だけ。
- **Linear**:
  - CLI（`gbdraw/linear.py:1164-1166`）は、ruler を指定しなければ scale の font に追従させる。
  - Web は unlinked のとき、ruler を指定しなければ config の既定値（short 20、long 12。`gbdraw/data/config.toml:67-69`）を使う。
  - short と long で既定値が違うので、CLI の単一の値では表せない。

### 3. 修正案
- **Circular**
  - 派生規則を Python の 1 つの helper（`circular_definition_interval_for_font`）にまとめ、`circular.py:951-954` と codec `:1700` の 2 つのコピーをそれに置き換える。
  - codec の `_migrate_flat_config_overrides`（`:1538`）で、circular mode かつ `objects.definition.circular.font_size` があり、`objects.definition.circular.interval` が無いときだけ interval を補う。
  - JS は変えない。request の bytes、Gallery の request、逆射影（`session-request.js:4502`）を変えずに済むため。
  - JS から派生値を送る案は採らない。逆射影で派生値が明示値に変わり、以後 font を変えても interval が追従しなくなるため。
  - Python API の `config_overrides` は codec を通らないので、leaf が独立である性質は保たれる。
- **Linear**
  - `run-info.js` の `appendConfigOverrides`（`:870-882`、`:944-955`）に条件を足す。
  - 条件: scale の font の override があり、ruler label の font の override が無く、かつ ruler label が描かれる（`objects.scale.style === 'ruler'` または `canvas.linear.ruler_on_axis`）。
  - このとき、理由を付けた `SourceRecipeUnavailable` を投げる。

### 4. 上位の対応
- 不要。

### 5. テスト
- **pytest（codec）**:
  - font 30 → interval 32。
  - font 30 と interval 20 → 20。
  - font 18 → 20。
- CLI、`--session`、Web の typed request の 3 経路で、HmmtDNA の SVG が一致すること。
- **node（`run-info.test.mjs`）**:
  - unlinked、Auto、scale 10 → unavailable。
  - linked → 2 つの flag が出る。
- Gallery の再生成が不要であることを確認する。全 session が font 18 か、interval を明示している。

### 6. 規模・PR
- Python 約 30 行と、`run-info.js` の 1 file。TR-07 と同じ PR にまとめてよい。

### 7. 依存・リスク
- dev で作った session のうち、font≠18 で interval の無いものは、再 Generate で間隔が変わる。未リリースの format なので互換の義務は無いが、fixture を grep して確かめる。
- font 30 では Web も CANNOT_FIT になり得る（0.13.0 と同じ挙動）。PV-08 の fit の案内とあわせて扱う。

---

## TR-07（P3）slot の legend label に "," や " #" があると CLI コマンドが壊れる

### 1. 分類
- **IMPLEMENT_EXISTING_AUTHORITY**（上と同じ Source recipe の契約）。
- CLI の文法に escape を足す案は、公開 grammar の変更なので別の判断になる（推奨しない）。
- Web でこれらの文字を禁止する案は、正しく表示できている label を削ることになるので不可。

### 2. 根本原因の再確認
- 監査のとおり。`circular-track-slots.js:1149` と `linear-track-slots.js:438-439` が値をそのまま連結している。Python は `gbdraw/tracks/parsing.py:30, 141`（strip_inline_comment）と `:118`（split_kv_list）で分割する。
- 追加の点:
  - slot id に `:`、`@`、`,`、` #` があっても同じく壊れる。
  - JS の parser（`circular-track-slots.js:1030-1033`）は `=` の無い token を黙って無視する。そのため JS で読み戻して確かめる方法では検出できない。
  - TSV（`gbdraw/io/cli_tables.py:816`）も同じ分割規則なので、TSV でも表せない。

### 3. 修正案
- `run-info.js` に私的関数 `assertCliSlotTokenLossless` を足す。
  - Python の分割規則だけを再現する: ` #` や `\t#` で切る → `@` で 1 回分割 → head を `:` で分割 → 残りを `,` で分割 → 各 token に `=` が必須。
  - 組み立てた token を分け直し、id、renderer、param が元と一致しなければ、Exact replay を案内する理由付きで unavailable にする。
- 既存の文字列 slot 用の検査（`:1165-1199`）もこの関数に寄せ、重複を消す。

### 4. 上位の対応
- 不要。

### 5. テスト
- node: 'GC skew (1 kb, AT-rich)'、'a #b'、id 'x:y' → unavailable になり、理由が付く。
- 同じ文字列の表を pytest で Python の parser に通し、判定が一致することを確かめる。

### 6. 規模・PR
- `run-info.js` の 1 file、+25/−10。

### 7. 依存・リスク
- 低い。

---

## FE-12（P3、P2 の候補）Specific table の色 `none`

### 1. 分類
- **IMPLEMENT_EXISTING_AUTHORITY**。
  - 値が有効かを決める owner は Python。
  - Web 自身が `none` の規則を作っている（`app/feature-editor/color-actions.js:1530-1531`）のに、Web の parser がそれを拒否するのは自己矛盾。
- 実行環境で判定が変わる挙動は NOT_ALLOWED なので取り除く。
- 4 桁と 8 桁の hex（α 付き）を許す案は新しい capability で、別の判断になる。

### 2. 根本原因の再確認（監査は過小評価）
色の定義域が 3 通りあり、食い違っている。

| 実行環境 | 受理する値 | 根拠 |
|---|---|---|
| Python（CLI と Web の render） | `none`、SVG の色名 147 個、#RGB、#RRGGBB | `gbdraw/features/colors.py:33-37`、`gbdraw/core/color.py:7, 27-35`、`gbdraw/io/colors.py:174-192` |
| browser の JS | CSS の色名（hex に解決）、3/4/6/8 桁の hex。`none` は拒否 | `app/file-imports.js:53-57` |
| DOM の無い JS | 英字の単語なら何でも受理（'notacolor' も通ることを node で確認） | `:54` |

- 実測（CLI）:
  - `#ff000080` と `#f008` は ValueError の traceback で停止する。
  - `notacolor` は 'Unknown color name' で停止する。
- 推定（Web）: `#ff000080` は import 時に受理され、Generate 時に `gbdraw/api/prepared.py:793` で失敗する。
- `read_color_table`（`gbdraw/io/colors.py:270-316`）は色を検証しないので、エラーに行番号が付かない。
- 推定、未検証: Web で `none` の規則を作って保存した Session は、Load 時の projection（`session-request.js:4143`、main thread）で throw する可能性がある。

### 3. 修正案
- **JS**: `app/color-utils.js` に `normalizeSpecificRuleColor` を置く。
  - `none`（大小無視）→ `'none'`。
  - #RGB と #RRGGBB → 小文字にする。
  - 色名 → browser で hex に解決する。
  - それ以外 → null。DOM の無い環境では色名も受理しない。
  - `file-imports.js:53-57` をこの関数に置き換え、失敗は TABLE_INVALID {row, color, SPECIFIC_COLOR} にする。
  - 同じ関数の他の 2 つの throw も typed にし、`error-normalization.js:322-325` と `:330-331` の正規表現を削除する。
- **Python**: `read_color_table` が各行の色を同じ定義域で検証する（`none` は許す）。失敗は `diagnostic` {TABLE_INVALID, color, SPECIFIC_COLOR, row} で raise する。
  - CLI では traceback ではなく、行番号付きの `ERROR:` が出る。
- 追加しないもの: α 付き hex への対応、色名表の JS への生成。

### 4. 上位の対応
- 不要。

### 5. テスト
- 共有の vector を作る: `tests/fixtures/specific_color_domain.json`。
  - valid: none、NONE、red、grey、#abc、#A1B2C3。
  - invalid: notacolor、#f008、#ff000080、rgb(1,2,3)、transparent。
- node（DOM stub あり・なしの両方）と pytest でこの vector を使う。
- 'none' の規則が serialize → parse の往復で壊れないこと。
- Playwright で、none の規則を含む Session の Save → Load → Generate が通ること。

### 6. 規模・PR
- Web 3 files と Python 1〜2 files。Ordinary。

### 7. 依存・リスク
- 色名入りの TSV を Node で読む既存 test は、hex か DOM stub に直す（OIPC-C08）。
- Worker から parseSpecificRules を呼ぶ経路が無いことを grep で確認する。現状はどの Worker も import していない。
- 'rebeccapurple' は Web では hex に変換されて通るが、CLI に直接渡すと拒否される。残るリスクは小さい。

---

## GUI の外で見つかった項目の扱い
- **`-n G` の IndexError と、`-w 0`、`-s 0` で track が空になる件**: X-02 の PR-A に含める。同じ typed options の owner で直る。
- **`--block_stroke_width -1` の traceback と、負の font size や stroke 幅の受理**: 設計は同じ（Python の型付き定義域と `diagnostic`）だが、別の PR にする。
  - config leaf に `gt=0` や `ge=0` を data として持たせ、`gbdraw/config/modify.py:258-290` で評価する。
  - font size >0 と stroke 幅 ≥0 は SVG/CSS の仕様で一意に決まるので、IMPLEMENT_EXISTING_AUTHORITY。
  - offset、spacing、track_axis_gap、label_rotation は、負の値に意味があるかを確かめる必要があるので EVIDENCE_REQUIRED。
  - `configurators/legend.py:67-75` の ValueError も ValidationError に変える。

## 横断的な提案
1. 失敗の意味は producer が code と有界 context で持ち、文言は normalizer が持つ。これを JS と Python の両方で統一する。英文の分類表は、件数が減る一方の ratchet で管理する。
2. context に保持している locator（Sequence、Line、Track row、Depth series、帯の px）は、必ず summary に表示する。
3. Generate は draft を書き換えない。draft 代入の件数を ratchet で 80 から減らしていく。
4. JS の projection は値を変換しない。範囲は Python が判断する。JS が Python より先に値を使う場合だけ、Python から生成した data を評価する。
5. Web と CLI の定義域の parity vector を 1 つ持ち、pytest と node の両方で回す。
6. CLI adapter に意味の規則（interval の派生、ruler の追従、`--depth_window` の検査）を置かない。規則は Python の helper 1 つに集める。
7. Source recipe は一般規則として「組み立てた argv を CLI の分割規則で読み戻し、一致しなければ unavailable」にする。TR-07 の関数を `--definition_line_style` や `--feature_shape` にも使う。
8. browser や第三者のライブラリが投げるエラーは、message で判定しない。operation が分かる呼び出しの境界で typed に包む。
9. error の actions は実際に押せるようにする。Generate ボタンを足し、Save 自身の失敗では Save ボタンを出さない。
10. 検証の test は、実際の producer を呼ぶ表駆動にする。文字列の定数を assert しない。

## 誤分類、またはバグではないと考える指摘
- **TR-05 と X-01**「fit error が slot id を落とす」「depth のエラーにファイル名が無い」: 名前を出さないのは PD-OI-046 のとおり。欠けているのは slotIndex、帯の半径、seriesIndex。
- **depth の値が数値でない場合**: 読めないのではなく内容の検証失敗なので、INPUT_UNREADABLE ではなく DEPTH_INVALID が正しい。案内の重複は REASONS.READ の設計の問題。
- **SE-09 の原因の記述**: 表が足りないのではない。Worker がすでに送っている code を、語彙に無いという理由で捨てている。
- **GE-03 の Circular**: recipe の不具合ではなく、Web の退行。
- **X-02 の「CLI も拒否する」**: window と step の 0 と負の値は、CLI も受理して track が空になる。GE-08（P3）と TR-04（P2）は原因が同じなので、重大度は P2 に揃えるべき。
- **GE-04 の「→ default」と CO-09 の「→ 0」**: 同じ現象。Linear の identity の既定値が 0 なので 0 に見える（Circular の conservation では 70 になる）。
- **FE-12 は過小評価**: 定義域が両方向に食い違っている。Session の往復が壊れる可能性もあるので、要確認の P2 候補。
- **"Canonical resource … is missing"**: 独立した検証の文言ではなく、IN-08（Generate 前の Save）の症状。
- **`NATIVE_VALIDATIONS` の 'An exact reference feature is required for alignment.'**: これを投げる producer が無い dead 行。
