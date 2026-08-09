# 次セッション改善計画

Status: completed on 2026-08-10  
Scope: current worktreeのfirst-party Python/Webコード、テスト、session/API境界  
Rule: 実装はこの文書を実行する次セッションで行う。この文書の作成ではproduction codeを変更しない。

関連文書: [Ponytail監査](PONYTAIL_AUDIT_2026-08-09.md)、[session互換性](../SESSION_COMPATIBILITY.md)、[Web規約](../../gbdraw/web/CLAUDE.md)

## 実施結果

- current schemaの`renderRequest`、resource、draft、UI、artifactのauthorityを分離し、CLI sidecarは既にresolve済みのrequestを再利用するようにした。公開`build_session_json()`の従来のresolve semanticsは維持した。
- cancellation transactionはgenerated artifactのcapture/applyだけをhistory ownerと共有し、feature identity、cache、worker lifecycleはcaller側に残した。
- Web behavior testはpytest経由のNode wrapperを廃止し、CIでNode 20の直接testとして実行するようにした。stale/重複fixtureは現在のcanonical identityとsession authorityへ移行した。
- endpointを持たないschema 3 homology matchをstandalone interactive SVGで扱えるようにし、record identityを保持する直接Node/Chromium回帰を追加した。
- reference SVG、public export、session/request version、offline runtime URL、生成owner assetに意図した変更はない。

検証結果:

- focused Python contracts: 638 passed、1 skipped
- direct Node: 62 passed
- full Playwright Chromium: 76 passed
- CLI recipe/onboarding contracts: 24 passed
- reference output comparison: 16 passed
- `ruff check gbdraw/`、browser wheel再生成、`verify_gui_offline.py check-assets`: passed
- non-slow pytest全体は修正前のonboarding timeoutを除き2902 passed、17 skippedで、修正後のowner testは5 passed

Phase 6の同一fixture中央値は、HEAD比でrenderがmulti Circular +0.7%、multi Linear -1.83%、single Circular +3.04%、session load -1.9%、write -2.27%、canonical decode +4.22%だった。すべて5%の中止基準内で、request/resource payload sizeはHEADと同一だった。Nodeプロセス全体の簡易値は環境変動と追加assertion時間を含むため、Web性能改善の根拠には採用しない。

`efaef9b2`基準のfirst-party production差分はPython +49行、Web JS +39行（合計+88行）、test差分は-218行、workflowは+14行だった。大幅なLOC削減ではなく、session resolve、generated artifact、behavior testのowner削減を成果とする。

## 目的

行数の大幅削減ではなく、次の変更波及範囲を小さくする。

- CLI、Python、Webが同じtyped request/planner/render経路を使う。
- 描画設定、未生成draft、UI状態、生成artifact、互換データを別ownerとして扱う。
- 同じ値の正規化、変換、保存をsurfaceごとに再実装しない。
- 現在の描画結果、session互換、offline/privacy、公開Python APIを維持する。
- 変更ごとに、production、test、docs、generated artifactの差分を分離してレビューできる。

20〜50%のLOC削減は目標にしない。監査基点ではPython+Web production約161,064行に対して、確実な削減候補は約5,290行（約3.3%）である。Finding #5の約860行はv27〜30 session replay readerであり、current writerの重複として削除してはいけない。したがって、次セッションの成功条件は「小さな純減」と「owner数・変更波及・重複変換の減少」である。

## 守る契約

次を変更しない。変更が必要になった場合は、そのphaseを中止し、別のbreaking-change計画に切り出す。

- package rootの初心者向けPython APIと`gbdraw.api`のtyped integration API
- Circular/LinearのSVG geometry、既存reference output、出力形式とoverwrite/error semantics
- session versions 27–33、39–40のreaderと、canonical request schemas 1、2、5の互換境界
- Webのsame-origin、offline、privacy、sanitization、worker分離、no-build構成
- `renderRequest`、resources、artifacts、draft、UIのauthority分担
- 既存の公開API contract snapshotとsession positive fixtures

## Phase 0: 開始時の再監査とベースライン

### 作業

1. 作業開始時に、既存の利用者変更と前セッションの差分を記録する。

   ```bash
   git status --short --untracked-files=all
   git diff --stat
   git diff --name-status
   ```

2. 監査に記載された対象が既に削除・変更済みかを再確認する。前セッションで消えたdead codeを復元しない。

   ```bash
   rg -n "v1|anchor_core_v1|renderRequest|session_to_cli_args|buildCanonicalRenderRequest" \
     gbdraw gbdraw/web tests docs/internal
   ```

3. 数値をHEADとcurrent worktreeで分けて記録する。generated、vendor、reference outputをproduction LOCに混ぜない。

4. 実装前の最小ゲートを実行する。

   ```bash
   pytest tests/test_public_contract.py tests/test_dead_api_cleanup.py -q
   pytest tests/test_api_requests.py tests/test_api_request_render.py \
     tests/test_session_io.py tests/test_session_request_codec.py -q
   node tests/web/session-authority.test.mjs
   node tests/web/session-request.test.mjs
   ```

### 完了条件

- current worktreeで既に行われた削除と、次セッションの対象が区別できる。
- baselineのテスト失敗、reference outputのdirty状態、ブラウザ実行可否を記録できる。
- 既存の失敗を新しい変更の失敗として扱わない。

## Phase 1: 描画設定のauthorityを固定する（P0）

### 目標

`renderRequest`を「最後に成功した描画意図」の唯一の正本にする。ただし、未生成draft、UI状態、resource bytes、SVG、解析cacheをrequestへ押し込まない。

### authority表

| 状態 | 正本 | してはいけないこと |
| --- | --- | --- |
| committed drawing intent | `renderRequest` / typed `DiagramRequest` | `config.form`やraw argvを再び描画設定の正本にする |
| input bytes / browser files | `resources`、`webFiles`のbinding | requestへbase64やFile objectを二重保存する |
| 未生成・inactive draft | Web editor/config draft | 未生成値を最後のrenderとして扱う |
| zoom、pan、selection、popup | `ui` | rendererのoptionsへ混ぜる |
| SVG、feature catalog、LOSAT cache、manifest | artifact namespace | fresh render requestのoptionsへ混ぜる |
| 旧session | `session_compat`/legacy reader | current writerへlegacy分岐を持ち込む |
| CLI invocation | provenance | v31以降のrender replayの入力にする |

### 作業

1. `gbdraw/api/request_render.py`、`gbdraw/session.py`、`gbdraw/session_request_codec.py`を、現行requestのvalidation・resolve・encode/decodeのownerとして確認する。
2. `gbdraw/web/js/services/session-authority.js`のauthority表と、`gbdraw/web/js/services/session-request.js`のcanonical projectionの不一致だけを列挙する。
3. 現在のsession writerで、committed render値を別のcurrent-state枝へ重ねて書いている箇所が本当にdraftやresource bindingを必要としない場合だけ削る。`config.linearComparisonPlan`、inactive track、resource original name、legacy reader用データは削除しない。
4. load pathは、current schemaではcanonical requestからUIをhydrateする一方向経路にする。旧schemaはreaderでcurrent requestへ移行してから同じ経路に入れる。
5. requestを複数回resolveしたり、rendererがraw configを直接読む経路がないことを確認する。

### 検証

```bash
pytest tests/test_api_requests.py tests/test_api_request_render.py \
  tests/test_session_io.py tests/test_session_request_codec.py \
  tests/test_api_library_usage.py -q
node tests/web/session-authority.test.mjs
node tests/web/session-request.test.mjs
node tests/web/session-resources.test.mjs
```

Playwrightが必要な場合は、Node packageだけで判断せず、Python Playwrightも確認する。

```bash
command -v playwright && playwright --version
python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
```

### 中止条件

- v27〜30の`session_to_cli_args()`と`_gui_session_to_cli_args()`を削除しない。これは互換readerであり、current writerの重複ではない。
- 未生成draft、inactive track、保存済みartifactが失われる場合は、current sessionの重複整理を中止する。
- public contract snapshot、session version、canonical request schemaが意図せず変わった場合は中止する。

### 期待値

純減は数百〜1,000行程度を上限の目安とし、削減量を約860行として計上しない。主成果は、session/requestのauthorityとcurrent/legacy境界の明確化である。

## Phase 2: Web状態とtransaction snapshotの重複を縮める（P1）

### 対象

- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/services/history-snapshot.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/session-authority.js`

### 方針

`run-analysis.js`のcancellation rollback全体をclone型history snapshotに置き換えない。history snapshotとcancellation transactionは、結果のidentity、cache/manifestの参照、asyncのclose、telemetry、worker cancellationという異なる契約を持つ。

次の範囲だけを共通化する。

1. generated artifactのdomain capture/applyを1 ownerにする。
2. cancellation固有のlocal state、result identity、cache/manifest参照、nextTick/nextFrame、transient closeはcaller側に残す。
3. mutation前後の結果identityとresource referenceをテストする。

### 検証

```bash
node tests/web/history-snapshot.test.mjs
node tests/web/run-analysis.test.mjs
node tests/web/session-authority.test.mjs
```

Playwrightでは、Generate → cancel → retry、session load → edit → save、stale result guardを対象にする。単純なsnapshot equalityだけで完了にしない。

## Phase 3: テストのownerを整理する（P1）

### 方針

`tests/test_web_packaging.py`のうち、production source textの識別子やcontrol-flow文字列を固定するテストを、既存の直接Node testまたはPlaywrightへ移す。

残すもの:

- wheel/sdistの内容
- CSP、manifest、entry point、vendor topology
- offline bundleのresource URL
- 直接のbehavior testでまだ覆われていないpackaging contract

削るもの:

- Pythonからinline JavaScriptを生成してNodeへ渡すだけの重複runner
- 既存Node testと同じ条件分岐をsource textで再確認するassertion
- pytestからNodeを呼ぶだけの同型wrapper

### 完了条件

- behaviorはNode/Playwrightがowner、package topologyはPython packaging testがownerになる。
- 失敗時に「どの実行環境の契約が壊れたか」が直接分かる。
- testを減らした箇所ごとに、同じbehaviorの直接testが存在する。

```bash
pytest tests/test_web_packaging.py -q
node --test tests/web/*.test.mjs
npx playwright test tests/web/*.playwright.spec.js
```

Node Playwrightが使えない場合は、対象フローをPython Playwrightで確認し、未検証のまま削除しない。

## Phase 4: Python APIとCLIの境界を保つ（P1）

### 方針

初心者向け`gbdraw` façadeとtyped `gbdraw.api`を同じ名前空間へ強制統合しない。現在はroot façadeがtyped requestへ変換し、CLI/Webもtyped plannerへ収束する構造が正しい。

1. 新機能はtyped APIを正本にする。
2. root APIは既存の簡易利用と`Diagram.save()`の互換adapterとして維持し、新機能を追加しない。
3. `circular.py`、`linear.py`のargparse処理は、parserの責務とrequest constructionの責務を混ぜない。
4. rendererが`argparse.Namespace`、CLI文字列、Web legacy configを受ける新経路を作らない。
5. 公開optionのdefault・validation・downstream consumerを1つのownerで追跡できる表を、変更時のチェックリストとして使う。

### 検証

```bash
pytest tests/test_public_contract.py tests/test_api_library_usage.py \
  tests/test_api_requests.py tests/test_api_request_render.py -q
ruff check gbdraw/
```

公開signature、module owner、signature hashが変わる場合は、先に互換方針を決める。snapshotを機械的に全面更新しない。

## Phase 5: dead code・旧テスト・重複manifestを段階的に削る（P1/P2）

Phase 0でcurrent worktreeを再確認した後、まだ存在するものだけを対象にする。

優先順:

1. live callerがないinteractive payload v1、旧private helper、test-only API
2. 既存Node/Playwright testで置き換え可能なsource-text test
3. capture側がmanifestを再所有しているscenario/fixture metadata
4. callerのないWeb exports、encoder-only path、debug plumbing
5. 参照のないplan文書、example、manifest entry、生成junk

各削除は次の単位で行う。

```text
caller確認 → owner確認 → 対応test確認 → production差分 → test差分 → docs/generated差分
```

公開alias、persisted reader、vendor fallback、active Gallery assetは、単にcallerが少ないという理由で削除しない。

## Phase 6: 性能と変更波及を測る

行数ではなく、次の値を変更前後で記録する。

- 固定fixtureのCircular/Linear render median
- session load/saveとcanonical request decode時間
- WebのGenerate → cancel → retry時間
- worker messageのrequest/resource payloadサイズ
- 同一resourceの読み込み回数とcache hit/miss
- 1つの設定変更で変更が必要になったproduction/test/docsファイル数

render medianが5%以上悪化した場合、またはsession restoreが遅くなった場合は、削減を成果とみなさず原因を調べる。性能改善を主張するには同一fixture・同一形式・同一browser条件で比較する。

## 次セッションの実行順

1. Phase 0のstatus、baseline、既存削除の確認
2. Phase 1のauthority matrixとcurrent writerの狭い整理
3. focused testsとpublic/session contract確認
4. Phase 2のsnapshot共通化を、identity/cancellationテスト付きで検討
5. Phase 3の重複testを1グループだけ移行し、実行時間と失敗診断を比較
6. Phase 4のAPI/CLI境界を確認
7. 残ったdead codeとmanifestをPhase 5で個別削除
8. 全差分をproduction、test、docs、generatedに分けてレビュー
9. 必要な範囲だけ次の小さな計画へ分割

## 完了条件

- current render、Web editor、session replay、CLI、Python APIのfocused testsが通る。
- Circular/Linear、single/grid/batch、offline browser経路に意図しない差分がない。
- canonical request、resource、draft、artifact、legacyのownerが文書とコードで一致する。
- 削除量と、変更波及ファイル数、性能測定値を記録している。
- public contract snapshotとreference outputの差分が意図したものだけである。
- 20%削減を理由に機能、互換性、offline/privacy、テストを無断で捨てていない。

## やらないこと

- `interface.py`を削除して初心者APIをtyped requestへ強制する
- v27〜30 session replayを削除して約860行を数える
- WebをCDN依存または新しいstate frameworkへ移す
- canonical requestへUI、draft、artifact、resource bytesを詰め込む
- 行数だけを目的にminify、コード生成、巨大generic registryを導入する
- generated wheel、reference output、Gallery owner assetを手編集する

この計画の成果は、数万行の見かけ上の削減ではなく、次の設定変更が少数のownerと検証で完了し、既存surface間の不整合を増やさないことである。
