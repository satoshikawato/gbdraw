# gbdraw Web Analytics 実装計画

- 対象: `gbdraw` Web UI
- 想定バージョン: 0.14 系以降
- 文書ステータス: 実装前の設計・受け入れ基準
- 最終更新: 2026-07-16

## 1. 目的

Web UI の訪問数ではなく、次の開発判断に必要な行動を、個人情報やゲノムデータを送らずに測定する。

1. 入力した利用者が、図の生成とエクスポートまで到達しているか。
2. circular / linear、single / multi、各 comparison workflow のどれが使われているか。
3. 各機能が実際に成果物へ反映されているか。
4. 特定機能の利用がエクスポート到達と結び付いているか。
5. 入力、比較、描画、エクスポートのどこで失敗しているか。

イベント数は増やさず、以下の 5 種類に限定する。

- `input_loaded`
- `render_completed`
- `feature_used`
- `export_completed`
- `operation_failed`

## 2. 設計上の決定

### 2.1 採用する方針

- 機能ごとのイベント名は作らず、`feature_used` の `feature_name` で区別する。
- 自由記述値は送らず、イベント名、パラメータ名、値を allowlist で固定する。
- 操作回数ではなく、意味のある状態変化または成果物への反映を記録する。
- `feature_used` は同一ブラウザータブのセッション内で、同じ分析コンテキストごとに 1 回だけ送る。
- 自動再描画や内部 reflow は `render_completed` に含めず、利用者が実行した Generate の成功だけを数える。
- ローカル利用・オフライン利用では analytics を完全な no-op にする。
- GA4 の障害、ブロッカー、CSP、設定ミスが gbdraw の動作へ影響しない構造にする。
- Web UI から GA4 へ生データ、ファイル情報、入力文字列、エラーメッセージを送らない。

### 2.2 初期スコープ

Phase 1 ではメイン SPA の入力、生成、主要機能、エクスポート、失敗を実装する。

Phase 2 で Gallery を追加する。Gallery は別ページであり、初期表示時の自動サンプル読込と利用者による明示的な選択を区別する必要があるため、Phase 1 へ混ぜない。

### 2.3 初期スコープ外

- 個人単位の追跡やログイン ID との結合
- CLI そのもののテレメトリー
- A/B テスト
- マウス移動、hover、pan、zoom、ドラッグ回数などの操作ログ
- 任意文字列を使うエラー分析
- サンプル、accession、feature、label、orthogroup 単位の分析
- BigQuery export を前提とする実装

## 3. KPI と分析目的

### 3.1 主 KPI

#### Output reach rate

入力を完了したセッションのうち、少なくとも 1 回エクスポートまで到達した割合。

```text
output_reach_rate =
  sessions with export_completed after input_loaded
  / sessions with input_loaded
```

この指標を、利用者が gbdraw で成果物を得られたかを表す主 KPI とする。ブラウザーから実際にファイルがディスクへ保存されたことは検知できないため、定義上の到達点は「ファイル生成成功かつダウンロード開始」とする。

### 3.2 ドライバー指標

#### Render completion rate

```text
render_completion_rate =
  sessions with render_completed after input_loaded
  / sessions with input_loaded
```

入力後に図を生成できたかを測る。Output reach が低下した場合に、生成前と生成後のどちらに問題があるかを分解する。

#### Render-to-export rate

```text
render_to_export_rate =
  sessions with export_completed after render_completed
  / sessions with render_completed
```

図の生成には成功したが、出力まで進まないケースを測る。

### 3.3 機能判断用の指標

#### Feature adoption rate

```text
feature_adoption_rate(feature_name) =
  sessions with feature_used(feature_name)
  / sessions with render_completed
```

分母は全訪問ではなく、少なくとも図を生成したセッションとする。入力前に使える一部機能についても、製品成果へ到達可能なセッションを基準に比較できる。

#### Feature-to-export rate

```text
feature_to_export_rate(feature_name) =
  sessions with feature_used(feature_name) followed by export_completed
  / sessions with feature_used(feature_name)
```

因果効果を意味する指標ではない。機能利用者のワークフロー完遂率として、維持・改善・導線見直しの仮説作りに使う。

### 3.4 ガードレール

#### Operation failure rate

```text
operation_failure_rate(stage) =
  operation_failed count
  / corresponding operation attempts
```

GA4 の 5 イベントだけでは全 stage の試行回数を直接記録しないため、初期レポートでは次の近似を使う。

- `input`: `operation_failed(input)` と `input_loaded` の合計を分母にする。
- `render`: `operation_failed(render)` と `render_completed` の合計を分母にする。
- `export`: `operation_failed(export)` と `export_completed` の合計を分母にする。
- `init` と `comparison`: 件数と影響セッション数を監視し、率が必要になった時点で BigQuery または追加設計を検討する。

### 3.5 目標値

初回リリース時には数値目標を置かない。計測の正しさを確認した後、2～4 週間のベースラインを採取し、次の順で目標を設定する。

1. 内部・テストトラフィックを除外する。
2. layout と comparison mode ごとの母数を確認する。
3. 週次変動とリリース前後の差を確認する。
4. Output reach rate の改善目標と Operation failure rate の上限を設定する。

母数が小さい区分は週次で断定せず、28 日ローリング値を併記する。

## 4. イベント契約

### 4.1 イベント一覧

| イベント | 送信条件 | 主な用途 |
|---|---|---|
| `input_loaded` | 描画に必要な主入力がブラウザー内で読み込み可能になった | 入力方式と funnel 起点 |
| `render_completed` | 利用者が実行した Generate が成功し、1 個以上の SVG 結果を得た | layout / workflow 利用と生成成功 |
| `feature_used` | 固定リストの機能が、意味のある状態変化または成功結果へ初めて反映された | 機能 adoption |
| `export_completed` | 出力ファイルの生成に成功し、ダウンロード開始処理まで進んだ | funnel 終点と出力形式 |
| `operation_failed` | 利用者が開始した処理が回復不能なエラーで終了した | stage 別の失敗監視 |

イベント名はこの 5 種類以外を analytics サービスから送信できないようにする。

### 4.2 固定パラメータと値

| パラメータ | 許可値 | 備考 |
|---|---|---|
| `app_version` | リリースバージョン | wheel 名から導出し、コードへ重複して直書きしない |
| `layout` | `circular`, `linear` | 画面上の layout |
| `record_scope` | `single`, `multi` | 生成対象レコード数に基づく |
| `comparison_mode` | `none`, `pairwise`, `orthogroup`, `collinear` | 比較なしを必ず `none` として送る |
| `input_format` | `genbank`, `gff3_fasta`, `session` | 主入力の論理形式 |
| `input_source` | `upload`, `session_restore` | サンプルは Gallery Phase 2 で別途扱う |
| `feature_group` | `input`, `track`, `editing`, `labels`, `workflow`, `discovery` | `feature_name` の固定分類 |
| `feature_name` | 4.6 の固定値 | 任意文字列は禁止 |
| `export_format` | `svg`, `png`, `pdf` | Interactive SVG も形式は `svg` |
| `export_variant` | `static`, `interactive` | PNG/PDF は `static` |
| `stage` | `init`, `input`, `comparison`, `render`, `export` | 失敗箇所 |
| `error_code` | 4.7 の固定値 | 生の例外文字列は禁止 |

未知値、`undefined`、空文字列は送信前に除外する。固定値を追加する場合は、実装、テスト、この文書、GA4 設定を同一変更で更新する。

### 4.3 `input_loaded`

必須パラメータ:

```js
{
  layout: "circular",
  input_format: "genbank",
  input_source: "upload",
  app_version: "0.14.0b0"
}
```

送信条件:

- Circular GenBank: 主 GenBank ファイルの読込が成功した時。
- Circular GFF: GFF3 と FASTA の両方が揃い、組として利用可能になった時。
- Linear: 各入力行について、GenBank または GFF3 + FASTA の論理入力が完成した時。
- Session import: session JSON の検証と状態復元が成功した時。`input_format: "session"`、`input_source: "session_restore"` とする。

送信しないケース:

- ファイル選択ダイアログを開いただけ。
- GFF3 と FASTA の片方しかない。
- 比較用、BLAST、config、style などの補助ファイルを読み込んだだけ。
- 同じ `File` オブジェクトを UI 再描画で再評価しただけ。
- 読込または検証に失敗した。

ファイル名、サイズ、更新日時、MIME type、パスは送らない。ページ内重複判定が必要な場合は `WeakMap<File, number>` のようなメモリー内識別だけを使い、GA4 や storage へ識別子を保存しない。

### 4.4 `render_completed`

必須パラメータ:

```js
{
  layout: "linear",
  record_scope: "multi",
  comparison_mode: "orthogroup",
  app_version: "0.14.0b0"
}
```

送信条件:

- 利用者が Generate を明示的に開始した。
- `runAnalysisInternal()` が `{ status: "ok" }` で終了した。
- 結果に 1 個以上の有効な SVG がある。

送信しないケース:

- window resize、panel resize、設定変更などに伴う内部 reflow。
- stale run として結果を破棄した処理。
- 利用者が cancel した処理。
- 生成途中で結果の一部だけが作られたが、全体は error になった処理。

`layout`、`record_scope`、`comparison_mode` は Generate 開始時に snapshot を取り、非同期処理中の UI 変更による値のずれを防ぐ。

### 4.5 `export_completed`

必須パラメータ:

```js
{
  export_format: "svg",
  export_variant: "interactive",
  layout: "linear",
  record_scope: "multi",
  comparison_mode: "pairwise",
  app_version: "0.14.0b0"
}
```

送信時点:

- SVG: Blob または出力文字列の生成に成功し、download click を開始した後。
- Interactive SVG: 対話機能を含む SVG の生成に成功し、download click を開始した後。
- PNG: canvas 描画と `toBlob` が成功し、download click を開始した後。
- PDF: PDF document の生成と `save()` 呼出しまで成功した後。

ブラウザーは利用者が保存を完了したかを確実には取得できないため、「download initiated」を完了定義とする。この制約を GA4 レポートの注記にも残す。

Interactive SVG は独立した `feature_used` にせず、次の組で分析する。

```js
{
  export_format: "svg",
  export_variant: "interactive"
}
```

これにより「機能を試した」と「成果物として出力した」を混同しない。

### 4.6 `feature_used`

共通例:

```js
{
  feature_group: "track",
  feature_name: "custom_track",
  layout: "linear",
  record_scope: "multi",
  comparison_mode: "pairwise",
  app_version: "0.14.0b0"
}
```

固定 registry:

| `feature_name` | `feature_group` | 初回送信条件 |
|---|---|---|
| `region_selection` | `input` | region crop が成功した生成結果へ反映された |
| `comparison_track` | `track` | 比較 track / ring を含む生成が成功した |
| `custom_track` | `track` | 有効な custom track を含む生成が成功した |
| `gc_content_track` | `track` | `gc_content` group を含む生成が成功した |
| `gc_skew_track` | `track` | `gc_skew` group を含む生成が成功した |
| `feature_selection` | `editing` | 選択 feature 数が 0 から 1 以上になった |
| `feature_editing` | `editing` | 色、stroke、visibility などの編集が実際の history commit になった |
| `label_editing` | `editing` | label text または visibility の編集が実際の history commit になった |
| `orthogroup_top_labels` | `labels` | orthogroup top label を含む生成が成功した |
| `undo` | `editing` | `history.undo()` が `true` を返し、状態が変化した |
| `redo` | `editing` | `history.redo()` が `true` を返し、状態が変化した |
| `cli_command_view` | `workflow` | 成功結果の Run info / CLI command を初めて表示した |
| `cli_command_download` | `workflow` | CLI helper ZIP の生成に成功し、download を開始した |
| `gallery` | `discovery` | Phase 2 で、明示的な sample 選択後の iframe load が成功した |

補足:

- GC content / skew は既定表示の場合があるため、この値は「明示的に ON にした人数」ではなく「その track を含む成果物を生成したセッション」を表す。
- `feature_editing` と `label_editing` は panel を開いた時ではなく、実際の差分が history へ commit された時に送る。
- feature 名は UI ラベルから動的生成しない。

同一セッション内の dedupe key:

```text
app_version + feature_name + layout + record_scope + comparison_mode
```

`sessionStorage` には固定 enum から作った key のみを保存する。`window.gtag` が存在せずイベント送信を試行できないローカル環境では、dedupe 済みとして記録しない。

### 4.7 `operation_failed`

必須パラメータ:

```js
{
  stage: "render",
  error_code: "diagram_generation_failed",
  layout: "circular",
  record_scope: "single",
  comparison_mode: "none",
  app_version: "0.14.0b0"
}
```

初期 error code registry:

| `stage` | `error_code` | 用途 |
|---|---|---|
| `init` | `runtime_init_failed` | Pyodide / runtime / 必須モジュール初期化失敗 |
| `input` | `missing_primary_input` | Generate 時に主入力が不足 |
| `input` | `invalid_input` | 主入力の読込・検証・parse 失敗 |
| `comparison` | `comparison_failed` | 比較データ生成または解析失敗 |
| `render` | `diagram_generation_failed` | diagram 生成本体の失敗 |
| `render` | `postprocess_failed` | SVG 後処理または結果組立ての失敗 |
| `export` | `svg_export_failed` | static SVG の生成失敗 |
| `export` | `interactive_svg_export_failed` | Interactive SVG の生成失敗 |
| `export` | `png_export_failed` | image load、canvas、DPI、`toBlob` を含む PNG 生成失敗 |
| `export` | `pdf_export_failed` | PDF document 生成または save 呼出し失敗 |

実装上の原則:

- `Error.message` の文字列解析で code を決めない。
- `runAnalysisInternal()` 内で `currentStage` を明示的に更新する。
- catch した場所で固定 code を返し、`app-setup.js` から 1 回だけ送る。
- `Error.message`、stack trace、Python traceback、stdout / stderr は送らない。
- 利用者の cancel は失敗として送らない。
- 同一操作について内側と外側の catch から二重送信しない。

### 4.8 記録しない操作

以下は送信しない。

- マウス移動と hover
- pan / zoom の各操作
- label の各ドラッグ
- slider の各 input event
- 色 picker の各変化
- 設定 panel を開いただけ
- 自動再描画と resize reflow
- inactive tab の表示
- undo / redo の no-op
- export button を押したが生成に失敗した場合の `export_completed`

## 5. コンテキスト判定

### 5.1 `record_scope`

- 生成対象レコードが 1 件なら `single`。
- 生成対象レコードが 2 件以上なら `multi`。
- upload slot 数ではなく、Generate 開始時に有効な入力件数で判定する。

### 5.2 `comparison_mode`

| UI 状態 | 値 |
|---|---|
| 比較なし | `none` |
| Circular comparison ring | `pairwise` |
| Linear BLAST upload / blastn / tblastx | `pairwise` |
| LOSATP blastp pairwise | `pairwise` |
| LOSATP orthogroup | `orthogroup` |
| LOSATP collinear | `collinear` |

コンテキスト判定を各 click handler へ重複実装しない。`analytics-context.js` の純粋関数へ集約し、state の組合せを table-driven test で固定する。

## 6. 実装アーキテクチャ

```text
利用者の操作
    │
    ▼
app / service の成功判定
    │
    ▼
analytics-context.js
  - state から固定 enum を生成
  - 結果から反映済み feature を判定
    │
    ▼
analytics.js
  - event / parameter allowlist
  - privacy guard
  - session dedupe
  - gtag がなければ no-op
    │
    ▼
window.gtag("event", ...)
  ※ Cloudflare Pages の hosted bundle のみ
```

analytics 層は UI state を直接 import しない。コンテキストを引数で渡すことで循環依存を避け、unit test を単純にする。

## 7. ファイル別の実装計画

### 7.1 新規: `gbdraw/web/js/services/analytics.js`

責務:

- 5 イベントの送信 API
- イベント・パラメータ・値の allowlist 検証
- `window.gtag` の存在確認
- `feature_used` の session dedupe
- analytics 内部例外の握りつぶしと開発環境向け debug log
- privacy guard

公開 API 案:

```js
configureAnalytics({ getContext, appVersion });
trackInputLoaded(params);
trackRenderCompleted(params);
trackFeatureUsed(featureName, params);
trackExportCompleted(params);
trackOperationFailed(params);
```

実装条件:

- generic な `track(eventName, arbitraryParams)` を app 側へ公開しない。
- `File`、`Blob`、`Error`、配列、任意 object を payload として受け付けない。
- 文字列値も parameter ごとの allowlist または version pattern を通す。
- analytics の例外を呼出元へ throw しない。
- `gtag` がない場合は何も保存せず return する。
- production で `debug_mode` を付けない。

### 7.2 新規: `gbdraw/web/js/app/analytics-context.js`

責務:

```js
buildAnalyticsContext(state);
hasCompletePrimaryInput(state);
detectRenderedFeatures(state, results);
```

要件:

- DOM text や表示 label ではなく state の正規化済み値を参照する。
- feature 検出は state と生成済み SVG の安定した group ID / metadata のみを使う。
- SVG 全体を GA4 へ渡さない。
- 結果検査は小さな文字列検索または DOM query に限定し、再 parse の負荷を避ける。
- Circular / Linear の分岐をこのファイルへ集約する。

### 7.3 `gbdraw/web/js/config.js`

- `GBDRAW_WHEEL_NAME` から `GBDRAW_VERSION` を導出する。
- analytics 用に version literal を別途追加しない。
- wheel 名を解析できない場合はイベント自体を落とすか、固定値 `unknown` を採用するかを実装前に決める。推奨はイベントを落として packaging test で検知すること。
- commit SHA、build path、環境変数は送らない。

### 7.4 `gbdraw/web/js/app/app-setup.js`

このファイルを instrumentation の composition root とする。

実装項目:

1. 起動時に `configureAnalytics()` を呼ぶ。
2. 主入力の完成状態を監視し、`input_loaded` を送る。
3. manual Generate の直前に analytics context を snapshot する。
4. `runAnalysisInternal()` の戻り値に応じて `render_completed` または `operation_failed` を 1 回送る。
5. 成功結果から `detectRenderedFeatures()` を呼び、該当する `feature_used` を送る。
6. undo / redo の戻り値が `true` の時だけ送る。
7. feature / label の history commit callback を接続する。
8. 成功結果後に Run info tab を初めて表示した時、`cli_command_view` を送る。
9. CLI helper ZIP のダウンロード成功時、`cli_command_download` を送る。
10. session import 成功時、`input_loaded` を送る。

注意点:

- generic file uploader へ一括 instrumentation しない。補助ファイルまで `input_loaded` として数えてしまうためである。
- result panel の既定 tab が自動表示されるだけなら `cli_command_view` を送らない。利用者の明示的な tab 選択を条件にする。
- app initialization で復元される UI state と、明示的な session import を区別する。

### 7.5 `gbdraw/web/js/app/run-analysis.js`

戻り値を次の形へ統一する。

```js
{ status: "ok" }
{ status: "cancelled" }
{ status: "stale" }
{
  status: "error",
  failureStage: "render",
  errorCode: "diagram_generation_failed"
}
```

実装項目:

- `currentStage` を `init` → `input` → `comparison` → `render` の順に更新する。
- catch block では固定の `failureStage` と `errorCode` を返す。
- 既存の利用者向け詳細エラー表示は維持するが、analytics payload へ混ぜない。
- stale / cancelled は `operation_failed` に変換しない。
- 成功判定は結果 SVG が 1 個以上確定した後に行う。
- reflow 用の内部呼出しには analytics を追加せず、外側の manual Generate wrapper だけで送信する。

### 7.6 `gbdraw/web/js/services/history.js`

history commit と analytics を疎結合にするため、optional callback を追加する。

API 案:

```js
createHistory({
  onActionCommitted: ({ analyticsFeature }) => {}
});

history.commit(change, {
  analyticsFeature: "feature_editing"
});
```

要件:

- 実際に差分がある commit のみ callback を呼ぶ。
- callback は optional で、analytics 無効時に既存挙動を変えない。
- history layer は GA4 や `window.gtag` を import しない。
- feature visibility など command-style API が boolean を返す場合は、成功した wrapper 側で送る方法でもよい。方式は feature editor 全体で統一する。

### 7.7 `gbdraw/web/js/services/export.js`

各 public export 処理の成功点と失敗点を明示する。

実装項目:

- static SVG、Interactive SVG、PNG、PDF の成功後に `trackExportCompleted()` を呼ぶ。
- export 開始時に生成結果の analytics context を受け取り、現在の UI state から再計算しない。
- 各 export の catch / failure branch で固定 `error_code` を返すか callback する。
- PNG では `img.onerror`、canvas context 不在、DPI 処理失敗、`toBlob(null)` を失敗として扱う。
- anchor の作成前や Blob 生成前に成功イベントを送らない。
- ファイル名を analytics API に渡さない。

`export.js` が analytics service を直接 import すると再利用性が落ちる場合は、`onCompleted` / `onFailed` callback を `app-setup.js` から渡す。依存方向は `app-setup.js` → service のまま維持する。

### 7.8 `tools/prepare_cloudflare_pages.py`

現行の hosted bundle への GA script 注入方式を維持し、次を追加する。

- analytics notice を「訪問」だけでなく「ページ利用と機能操作の集計」を含む表現へ更新する。
- ローカルの `gbdraw/web/index.html` に GA measurement ID を直書きしない。
- hosted bundle にだけ script、CSP、notice を注入する。
- preview deploy は本番 property を汚染しないよう、test data stream を使うか analytics を無効化する。
- marker の置換漏れと二重注入を packaging test で検出する。

### 7.9 新規: `docs/WEB_ANALYTICS.md`

実装 PR で、長期運用用の短い event contract を別文書として作る。本計画書から次だけを抜き出す。

- 5 イベントの定義
- parameter / enum registry
- feature registry
- error code registry
- privacy denylist
- GA4 custom dimension 一覧

本計画書は実装手順と意思決定を残し、`WEB_ANALYTICS.md` を将来の registry の正本とする。

## 8. Gallery Phase 2

対象:

- `gbdraw/web/gallery/index.html`
- `gbdraw/web/gallery/gallery.js`
- `tools/prepare_cloudflare_pages.py`

実装条件:

1. Gallery index に GA script 用 marker と notice 用 marker を追加する。
2. hosted build 時のみタグを注入する。
3. Gallery の CSP に必要な Google Analytics endpoint を hosted bundle で追加する。
4. メイン SPA と同じ generic analytics service を再利用する。
5. 初期表示時に自動読込される sample は `gallery` として数えない。
6. 利用者が sample を明示的に選択し、その後 iframe load が成功した時だけ `feature_used: gallery` を送る。
7. sample ID、accession、title、URL を送らない。
8. Gallery からの SVG / session download はメイン SPA の `export_completed` に含めない。含める場合は funnel 定義の再設計を先に行う。

## 9. プライバシーとコンプライアンス

### 9.1 送信禁止データ

次の値は analytics API の引数へ渡さない。

- ファイル名、ファイルパス、ファイルサイズ、更新日時
- accession、record ID
- organism、species、strain
- feature 名、label text、gene / product 名
- orthogroup 名または ID
- CLI command と引数
- 塩基配列、アミノ酸配列
- SVG 内容、session JSON、config 内容
- 利用者が入力した任意文字列
- 生のエラーメッセージ、traceback、stack、stdout、stderr
- IP、メールアドレス、利用者 ID として使える独自識別子

### 9.2 実装ガード

- analytics service は primitive な allowlist 値だけを受け取る。
- object の再帰 serialize は行わない。
- `Error`、`File`、`Blob` を検出したらイベント全体を破棄する。
- analytics payload を production console に出力しない。
- GA4 で User-ID を設定しない。
- Google Signals、広告 personalization、data retention は運用ポリシーに合わせて最小化する。

### 9.3 リリース前の確認事項

次はコードだけでは決められないため、公開前の release gate とする。

- 対象地域と運用主体に必要な cookie consent / Consent Mode の要否
- privacy notice の文面
- GA4 property の data retention
- internal traffic と開発者トラフィックの除外方法
- preview / staging 用 data stream

## 10. GA4 管理画面の設定

### 10.1 カスタムディメンション

以下を event scope で登録する。

1. `app_version`
2. `layout`
3. `record_scope`
4. `comparison_mode`
5. `input_format`
6. `input_source`
7. `feature_group`
8. `feature_name`
9. `export_format`
10. `export_variant`
11. `stage`
12. `error_code`

パラメータ送信開始前または同時に登録する。登録前の値が標準レポートで遡及利用できると仮定しない。

### 10.2 Key event

`export_completed` を key event として登録する。`render_completed` は診断用イベントとして扱い、初期段階では key event にしない。

### 10.3 Exploration

最低限、次を作成する。

#### Funnel: Input to output

```text
input_loaded
  → render_completed
  → export_completed
```

- closed funnel
- ordered steps
- segment: `layout`, `record_scope`, `comparison_mode`

#### Feature adoption

- row: `feature_name`
- value: Sessions、Total users
- filter: `event_name = feature_used`
- breakdown: `layout` または `comparison_mode`

#### Feature to export

```text
feature_used(feature_name)
  → export_completed
```

- same session
- ordered steps
- breakdown: `feature_name`

#### Failure monitoring

- row: `stage`, `error_code`
- value: Event count、Sessions
- breakdown: `app_version`

### 10.4 データ品質チェック

リリース直後は DebugView と Realtime で次を確認する。

- イベント名が 5 種類だけである。
- parameter 値が registry 内に収まっている。
- 1 回の Generate 成功で `render_completed` が 1 回だけ送られる。
- reflow でイベントが増えない。
- feature dedupe が意図どおり働く。
- local build からイベントが来ない。
- preview traffic が production 集計へ混ざらない。
- payload に禁止データがない。

## 11. テスト計画

### 11.1 Unit test: analytics transport

新規: `tests/web/analytics.test.mjs`

必須ケース:

- `window.gtag` がなくても例外にならない。
- 送信できる event 名が 5 種類だけである。
- 未知の event、parameter、enum 値を破棄する。
- `File`、`Blob`、`Error`、object を payload に含められない。
- version が wheel 名から正しく導出される。
- 同一 feature / context は session 内で 1 回だけ送られる。
- layout または comparison mode が異なる場合は別 context として送られる。
- `gtag` 不在時に dedupe 済みにならない。
- analytics 内部の storage 例外が app へ伝播しない。

### 11.2 Unit test: context mapping

新規: `tests/web/analytics-context.test.mjs`

table-driven test:

- Circular single / multi
- Circular comparison
- Linear comparison none / pairwise / orthogroup / collinear
- GenBank / GFF3 + FASTA / session restore
- incomplete GFF3 + FASTA
- region selection
- custom track
- GC content / GC skew
- orthogroup top labels
- result SVG が空のケース

### 11.3 Packaging test

更新: `tests/test_web_packaging.py`

Phase 1:

- source `index.html` に実 Measurement ID がない。
- source bundle は analytics endpoint へ接続しない。
- hosted main index には GA tag、必要 CSP、notice が 1 回だけある。
- marker が hosted output に残らない。
- preview 設定が production stream を使わない。

Phase 2:

- hosted Gallery にだけ GA tag と必要 CSP が入る。
- Gallery の marker が残らない。
- local Gallery は外部 analytics を読み込まない。

### 11.4 Browser integration test

fake `window.gtag` を注入し、送信 call を配列へ記録する。

正常系:

1. primary input を選択する。
2. `input_loaded` が 1 回送られる。
3. Generate を実行する。
4. `render_completed` と、反映された `feature_used` が送られる。
5. Interactive SVG を export する。
6. `export_completed(svg, interactive)` が送られる。

異常・境界系:

- Generate failure で `operation_failed` が 1 回だけ送られる。
- cancel と stale run では failure を送らない。
- resize reflow で `render_completed` を送らない。
- no-op undo / redo で `feature_used` を送らない。
- 有効な undo / redo は session 中 1 回だけ送る。
- PNG の `toBlob(null)` で `operation_failed(export, png_export_failed)` を送る。
- export failure で `export_completed` を送らない。
- payload に filename、SVG、error message が含まれない。

Node の `@playwright/test` が利用可能なら既存 JS spec の方式を使う。利用できない場合は Python Playwright で同等の targeted browser check を行う。

### 11.5 実行コマンド案

```bash
node --test tests/web/analytics.test.mjs tests/web/analytics-context.test.mjs
python -m pytest tests/test_web_packaging.py -k analytics -q
python -m pytest tests/ -m "not slow"
```

ブラウザーテストは repository の Playwright 利用可否を確認してから実行する。

## 12. 実装順序

### Phase 0: 契約と管理設定の準備

- [ ] `docs/WEB_ANALYTICS.md` を作成する。
- [ ] event、parameter、feature、error code registry を確定する。
- [ ] privacy notice と consent 要件を確認する。
- [ ] production と preview / test の data stream 方針を決める。
- [ ] GA4 custom dimensions を登録する。

完了条件: 実装者とレビュー担当者が「何を送るか／送らないか」を registry だけで判定できる。

### Phase 1A: Analytics 基盤

- [ ] `analytics.js` を追加する。
- [ ] version 導出を追加する。
- [ ] allowlist、privacy guard、no-op、dedupe を実装する。
- [ ] `analytics.test.mjs` を追加する。

完了条件: UI へ接続しなくても、transport の全 guard を unit test で確認できる。

### Phase 1B: コンテキスト

- [ ] `analytics-context.js` を追加する。
- [ ] layout、record scope、comparison mode を実装する。
- [ ] primary input 完成判定を実装する。
- [ ] rendered feature 検出を実装する。
- [ ] table-driven test を追加する。

完了条件: Circular / Linear の主要 state 組合せが固定 enum へ一意に変換される。

### Phase 1C: 入力・生成・失敗

- [ ] `app-setup.js` で primary input を監視する。
- [ ] `run-analysis.js` の structured result を実装する。
- [ ] manual Generate wrapper から成功・失敗イベントを送る。
- [ ] cancel、stale、reflow の除外を実装する。
- [ ] feature result detection を接続する。

完了条件: input → render の funnel と stage 別 failure が重複なく取得できる。

### Phase 1D: 編集・履歴・CLI

- [ ] history commit callback または統一 wrapper を実装する。
- [ ] feature selection / editing / label editing を接続する。
- [ ] undo / redo の成功判定を接続する。
- [ ] CLI command view / download を接続する。

完了条件: panel open や no-op では送信されず、実際の状態変化だけが session 中 1 回記録される。

### Phase 1E: Export

- [ ] static SVG 成功・失敗を接続する。
- [ ] Interactive SVG 成功・失敗を接続する。
- [ ] PNG の全非同期 failure branch を接続する。
- [ ] PDF 成功・失敗を接続する。
- [ ] export context を生成結果へ固定する。

完了条件: 成功イベントはファイル生成成功後だけ送られ、各失敗は固定 code になる。

### Phase 1F: Hosted packaging と検証

- [ ] Cloudflare Pages の注入処理と notice を更新する。
- [ ] packaging test を追加する。
- [ ] fake gtag の browser integration test を追加する。
- [ ] test data stream の DebugView で検証する。
- [ ] production へ main SPA をリリースする。

完了条件: source / local build は外部送信せず、hosted build だけが承認済み設定で送信する。

### Phase 1G: データ品質確認

- [ ] リリース後 24～72 時間の event / parameter を監査する。
- [ ] 5 event 以外がないことを確認する。
- [ ] enum の未知値と `(other)` の兆候を確認する。
- [ ] funnel の step 順序と重複を確認する。
- [ ] 2～4 週間の baseline を採取する。
- [ ] KPI 目標と failure guardrail を設定する。

### Phase 2: Gallery

- [ ] hosted Gallery の tag / CSP / notice 注入を実装する。
- [ ] 明示的な sample 選択と iframe load 成功を接続する。
- [ ] initial auto-load を除外する。
- [ ] sample 識別情報が payload にないことをテストする。
- [ ] Gallery download を main funnel から分離する。

## 13. 受け入れ基準

すべて満たした時に Phase 1 完了とする。

- [ ] analytics event 名は 5 種類だけである。
- [ ] parameter 名と値が固定 registry に限定される。
- [ ] source / local / offline Web UI は `gtag` 不在でも完全に動作する。
- [ ] analytics 障害が Generate、編集、export を失敗させない。
- [ ] ファイル情報、ゲノム情報、任意文字列、生エラーが payload にない。
- [ ] `render_completed` は manual Generate 成功時に 1 回だけ送られる。
- [ ] cancel、stale run、reflow は failure / render 成功として数えない。
- [ ] feature は意味のある状態変化または成功結果への反映時だけ送られる。
- [ ] `feature_used` は定義した context ごとに session 中 1 回へ dedupe される。
- [ ] undo / redo の no-op は記録されない。
- [ ] `export_completed` はファイル生成成功後だけ送られる。
- [ ] Interactive SVG を static SVG と区別できる。
- [ ] feature 利用から同一 session の export 到達を分析できる。
- [ ] production、preview、local の traffic が意図せず混在しない。
- [ ] custom dimensions、key event、4 種類の Exploration が設定されている。
- [ ] unit、packaging、browser integration test が通る。

## 14. レビュー時の重点項目

コードレビューでは次を先に確認する。

1. イベントを送る場所が「操作開始」ではなく、定義した成功点になっているか。
2. instrumentation のために既存 app service が analytics へ強く依存していないか。
3. payload が allowlist から組み立てられ、既存 object を spread していないか。
4. async failure、cancel、stale、reflow の境界がテストされているか。
5. version、feature、error code の registry が複数箇所へ重複していないか。
6. Gallery と main SPA の funnel が混ざっていないか。

## 15. 変更時の運用ルール

新しい機能を計測対象へ追加する場合:

1. 開発判断に使う問いを 1 文で書く。
2. 既存の `feature_used` で表現できることを確認する。
3. 成功または意味のある状態変化を送信点として定義する。
4. 固定 `feature_name` と `feature_group` を registry へ追加する。
5. privacy review と unit test を追加する。
6. GA4 custom dimension の追加が不要であることを確認する。

新しいイベント名は原則追加しない。5 イベントで表現できない場合は、既存 KPI では答えられない具体的な意思決定と、イベント追加による運用コストを設計レビューで示す。

## 16. 参考資料

- [Google Analytics: Best practices to avoid sending Personally Identifiable Information (PII)](https://support.google.com/analytics/answer/6366371?hl=en)
- [GA4: About custom dimensions and metrics](https://support.google.com/analytics/answer/14240153?hl=en)
- [Google tag: Set up event parameters](https://developers.google.com/analytics/devguides/collection/ga4/event-parameters)
- [GA4: Recommended events](https://support.google.com/analytics/answer/9267735?hl=en)
- [Google Analytics privacy controls](https://support.google.com/analytics/answer/9019185?hl=en)
