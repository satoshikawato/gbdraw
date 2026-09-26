# Issue #597 — 総合実装計画

## 目的・作業対象

gbdraw は microbial genome の Circular/Linear 図を作る Python tool と、build step のない Web SPA を持つ。
本計画は [Issue #597](https://github.com/satoshikawato/gbdraw/issues/597) の次の3領域を扱う。

| 問題 | 到達すべき状態 |
| --- | --- |
| Circular に一つの GenBank uploader しかない | 複数選択・後から追加し、元ファイル単位で Replace/Remove。各ファイル内の全 record を保持する。 |
| record discovery と crop の適用条件が分かりにくい | 通常 upload は自動探索。適用可能な一件用 controls を自動展開。saved preview は遅延探索を明示し、Inspect/Generate から継続できる。 |
| 大きい Save/Load が主スレッドを停止させる | 既存 streaming export を保ち、import parse を専用 JS Worker に移す。projection/validation/restore の長い処理を分割し、閲覧可能な応答性を保つ。 |

実装対象ブランチ: `fix/issue-597-input-session-20260926`。
計画開始 base: 最新 `origin/dev` の `d457b7189b137185a8dec800819a312c30b969fa`、2026-09-26取得。
対象コードは `gbdraw/web/` と既存 typed Session/request の境界、関連 tests、再生成が必要な artifacts、既存 docs。
本資料の `decisions/`、`sessions/`、`evidence/`、`results/`、`authority-candidates/` は
`docs/internal/issue-597-input-session-implementation-20260926/` 配下。Web module の短いpathsは
`gbdraw/web/js/` からの相対path、`index.html` は `gbdraw/web/index.html` を指す。

各セッションは [SESSION_WORKFLOW.md](./SESSION_WORKFLOW.md) に従って、この remote branch を独立 checkout に取得する。
shared checkout の branch、依存環境、ローカルサーバー、他 session の worktree を操作しない。
各セッション終了時に担当差分・結果を一つの commit にし、同名 remote branch に push する。
実装セッションを同時に走らせない。複数のセッションが同じ shared branch に同時 push する workflow は採用しない。

## 規範資料と正規の authority

本計画は次の規範を参照する。実装担当者は開始時に最新 branch の本文を読み、計画と異なる場合は規範の境界を確認する。

- [Architecture Fitness-function Ratchet](../ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md): owner/path/compatibility の定義、ordinary evidence、例外 packet と exact-head review。
- [Product Impact Ratchet](../PRODUCT_IMPACT_RATCHET.md): Product preflight、authority 検索、承認記録の保存、authority-before-runtime。
- [Product Decision Packet Template](../PRODUCT_DECISION_PACKET_TEMPLATE.md): 未承認の別 concern が生じた場合の独立した判断資料。
- [Web Change Policy](../WEB_CHANGE_POLICY.md): trusted-base check、Gate/Review、checker/evidence/authority/runtime の分離。

正規の実装・authority の配置は repository guidance と上記規範を参照する。
S00 は `tools/web-product-impact-map.json`、`tools/web-product-decisions.json`、
`tools/web-change-policy.json`、`tools/web-architecture-rules.json` と
`docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` の存在・適用範囲を調べる。
候補 receipt を保存した本ディレクトリを、それらに代わる authority として実装しない。

## 用語

- source: 一回の upload に対応する入力 instance。元の file bytes と UID を持つ。
- record: 一 source 内の配列要素。accession は表示名で、同一性は source UID と local index で決める。
- resource: 保存された byte payload。同じ payload を再利用しても source instance は別に保つ。
- active draft: 次の Generate に使う編集中の入力・設定。committed artifact と異なる状態を持てる。
- committed artifact: 最後に成功した request、Result、catalog、関連 provenance。
- semantic mutation: source/config/editor/History/cache/Result など、保存内容または生成意味を変更する操作。
- authority: base branch に統合された支持挙動または認識済み Product/privileged policy。単なる新コードや tests は authority ではない。
- OE/PE/CB: semantic owner excess、canonical path excess、active persisted compatibility path count。Architecture Ratchet の独立した3軸。

## 承認済み Product outcome

以下3件の選択、rationale、must-preserve、may-retire、residual risk、owner/date は文書化済み。
実装担当者は Product 選択の再確認を要求しない。

| Concern / revision | 選択 | 完全な記録 |
| --- | --- | --- |
| `diagram-generation.circular-source-collection` / 1 | `INDEPENDENT_ORDERED_SOURCES` | [source collection](./decisions/01_CIRCULAR_SOURCES.md) |
| `diagram-generation.circular-transform-discoverability` / 1 | `REVEAL_APPLICABLE_SINGLE_RECORD_CONTROLS` | [discovery](./decisions/02_RECORD_DISCOVERY.md) |
| `web.session-operation-consistency` / 1 | `EXCLUSIVE_SEMANTIC_SESSION_OPERATION` | [Session operation](./decisions/03_SESSION_OPERATIONS.md) |

Product approval は failing security/scientific/performance/architecture gate を免除しない。
認識済み authority への統合は S00、dev の trusted base への反映は runtime 開始条件。
未承認の outcome が新たに必要なら、影響する範囲だけ Product preflight/Decision Pack を作る。
既に承認された3件の範囲を広げたり、例外判断をその承認から推定したりしない。

## 確認済み基準動作と残る調査

基準コードでは `files.c_gb` は scalar、Circular uploader は `multiple=false`。
typed request の `records[]` は複数 source を既に表現し、CLI 由来 Session の複数 source は composite file view に復元される。

通常の GenBank upload は watcher → `refreshCircularRecordOrder()` → `record-discovery.js` の JS fast path を自動起動する。
Chromium で重複 accession の2 records、Generate 前の rotation spinbuttons 2、manual Load links 0、Worker constructions 0 を観察した。
BUG-02 を「自動探索コードがない」と診断して parser を追加しない。
saved-preview Load は rows 0/Worker 0 のまま source discovery を遅延するが、status が loading と表示される。
新 native source に交換すると再び自動探索する。

Save は既に256 KiB分割 JSONと `CompressionStream`、8 chunks ごとの task yield。
Load の `JSON.parse(text)`、preflight、restore と Save projection の一部が main thread に残る。
専用 Worker の full-object reply 自体が long task/heap を増やす可能性は未測定。
大規模 Save/Load の性能方式は S01 で計測して一つに確定する。Worker 使用だけを完了条件にしない。
既存5 Nodeファイル・30 checks は pass。これは新 collection/transport の完成証拠ではない。

## 責任・正規状態・実行経路

| 責任 | 最終 owner / boundary | 削除・維持すべき経路 |
| --- | --- | --- |
| Circular sources の add/replace/remove/reset/restore と discovery reconciliation | `app/circular-sources.js` | `run-analysis.js` 内の独立 Circular discovery lifecycle を同じ change で削除 |
| Source の軽量 record metadata と rare-format helper | 既存 `app/record-discovery.js` / diagram helper | parser 全面改修、main-thread Python、並行 parser pipeline を追加しない |
| Source-bound record display intent | 既存 `app/record-display-options.js` | fixed Circular source UID を collection の UID に対応させる |
| Fresh/restore canonical request、record ordering、grouping、resources | 既存 `services/session-request.js` | fresh-only builder、argv generation、合成 GenBank renderer を追加しない |
| Persisted file binding grammar/backing と payload reuse | 既存 `session-resource-backing.js` / `session-resources.js` | current writer は一つ。old reader の normalization だけを境界に残す |
| Save/Load coordination、validation、candidate、atomic adoption | 既存 `services/config.js` | codec と render model の第二 owner にしない |
| JSON/gzip import transport/lifecycle | `services/session-file.js` 境界下の focused import client / `workers/session-import-worker.js` | `config.js` の直接 JSON.parse を削除。transport は一方式のみ |
| JSON/gzip export | 既存 `services/session-file.js` | 一括 stringify や export専用Workerへ置き換えない |
| Current Result admission/sanitization | 既存 `services/svg-result-ingestion.js` | Worker data を直接 mounted SVG/Result に渡さない |
| Composition/availability wiring | `app/app-setup.js`、既存 mutation owners | templates/watchers は adapters。二重 lock/selected/visibility refs を追加しない |

Generate は引き続き `run-analysis.js → session-request.js → diagram-generation.js → diagram Worker → typed renderer`。
Pyodide runtime は diagram Worker の一つだけ。新 Session import Worker は JS のみ。
入力・結果は browser 内、runtime assets は same-origin、SPA に build step を加えない。

### Circular source collection と identity

正規入力は `files.circularGenbankSources = [{ uid, file }]`。
追加は配列末尾、同時選択は FileList の順序、file 内 records は原順。
Replace は同じ card 位置に新 UID、Remove は対象 UID だけを除く。
同名・同一 accession・同じ bytes の再追加を別 source instance として残す。
resource payload の同一化と source/record instance の同一化を混同しない。

per-source discovery catalog に version、status、records、error、operation を持つ。
全体 record list と count は computed view。別の mutable flat list を持たない。
record key は UID/local `#index` から作り、request selector は元 source 内の local index。
UIDを filename/hash/card order だけで決めない。`recordKey` と source UID を同一用途の二つの mutable ID にしない。

position、depth sparse matrix、transform draft、annotation selection、batch output prefixes、comparison endpoints を
同じ ordered record catalog から投影する。record を削除して隣の record に index を流用しない。
対象 source 交換・削除時だけ intent を無効化し、別 source の rows/transforms/depth/selector は保つ。
依存 imported comparison は既存 unresolved-intent owner に渡し、自動的に No comparison にしない。
committed request/Result/download は次の successful Generate まで保持。
add/replace/remove は一つずつ History operation、discovery completion は追加 History entry ではない。
Undo/Redo/restore は source version に対して reconcile し、watcher の実行時機に依存しない。

fresh/reset の grouping は既存 shared canvas、saved explicit opt-out は保持、明示一件選択は single。
single/grid/batch と valid one-record grid を維持。file 数から grouping を推測しない。
crop は一件選択の既存 journey、grid の record rotation は既存 display intent を使う。
既存単一 GFF3/FASTA 入力は維持する。multi-GFF pairing、file reorder、任意 record subset は範囲外。

### Bindings と互換性

new writer は bindings schema 3 の `circularGenbankSources: [{ uid, file: <file binding> }]`。
`c_gb` は current inventory から削除するが released reader は残す。
通常 `file` は leaf binding。真実な元 source/intent mapping が復元不能な legacy combined input は
一 source とし、既存 composite backing/grammar を `file` 境界で使う。Load だけで再エンコードしない。
信頼できる components と既存 request identity で個別元 source を復元できる場合は元順・繰返しを保つ。
推測で合成 bytes を分割しない。

基準 `origin/main` は `4556e04e929a4a85ad28d1833ce7304bd764881c`、Session 42/bindings 1/2。
tag `0.13.0` は Session 30。current dev は Session 44/request 8/bindings 2。
実装開始時に main first-parent/tags を再確認し、Session/request/bindings/catalog/cache の namespace を別に記録する。
未公開の Session 44/bindings 2 artifacts は Session 44/bindings 3 に再生成し、その組だけの reader を作らない。
既に公開された namespace があればその形式を bounded reader と positive fixture で保持する。
Session 44 と request 8 は維持し、必要性のない version 増加や branch-only migration chain は作らない。
reader normalization で source UID と旧 transforms/positions/depth/selector を一度だけ対応付ける。
CLI/Python session validation/materialization が Web binding inventory を使う箇所も更新し、別 decoder を作らない。

### Discovery と一件用 controls

status は `idle / deferred / loading / ready / error`。
new upload/replace は自動 discovery。通常形式は JS fast path、rare format は既存 full-parser helper。
completion で UID/version/mode/input type を検証する。Retry は新 parser pipeline ではない。
source failure は安全な filename/stage/error を示し、別 source と Result を保つ。

saved preview Load は Python Worker を構築せず source bytes の再探索を始めない。
`Records not inspected` と `Inspect source records` を示し、Inspect/Generate の operation を待つ。
active draft と committed artifact が異なる場合、catalog count を current source の検証済み count として使わない。
loading は実際の in-flight operation に限る。

一件 section は `Single-record crop, orientation and titles`。
upload/selector change で一件編集が applicable になったときに自動展開する。
manual close を保ち、無関係な更新では再展開しない。upload/select の focus と scroll anchor を保つ。
grid/batch の適用不可理由と一件選択への導線は section 外に出す。
grouping の自動変更、先頭 record の自動選択、crop条件の緩和はしない。
disclosure/status は transient UI で、追加 Session schema/History entry/Generate trigger にしない。

### Session の排他・応答性・atomic adoption

Save/Load 中は semantic mutations を停止し、閲覧・scroll・pan/zoom・検索を維持する。
既存 `sessionSavePending` / `sessionImportPending` から単一 availability を導出する。
Generate/source/editor/History/Reset/cache/別 Load/Save の owning actions と DOM を同じ predicate に従わせる。
programmatic request も `{ status: 'busy', reason }` 相当の明示結果と retry continuation を持つ。
重複 Save は既存 single-flight join、一度の download。異種操作は join しない。
既に Generate/automatic reflow が動いている場合は Save/Load を開始せず理由を出す。
pending の publish、start、settlement、error/取消、teardown を一つの existing session lifecycle にまとめる。
遅れて完了する file import/auto rerender/cache mutation も同じ consistency checkpoint の前後で処理する。

タイトル確定後の Save 開始から handoff/error/cancel まで同じ document が対象。
adopted immutable payload を再利用し、小さい mutable config/overrides を必要な範囲で固定する。
full graph clone/sign/hash/base64 再変換、独立 snapshot owner、汎用 lock/queue framework は追加しない。
Load は private candidate の validation/normalization/preparation を完了してから既存 admission/transaction で採用。
failed adoption は主データを先に復元し、その後 transient UI を reconcile。旧 Result/History を保つ。
private candidate construction の yield を、部分的な live document の公開に使わない。

S01 は File/Blob → Worker decode/parse → main candidate の transport を実測する。
full-object reply が満たせばそれを一方式として採用。long clone/peak memory が残れば
既知 Session sections、feature batches、巨大 resource/TSV strings の bounded transport を先に比較する。
計測結果から一 transport を選び、通常 runtime に size-based fallback 二経路を残さない。
Vue proxies/DOM/backing internals を postMessage しない。codec error を main JSON.parse fallback で隠さない。
client が operation ID、structured error、stale rejection、settlement termination、teardown を所有する。
`config.js` の validation を Worker に複製せず、返却 data を信頼済み Session と扱わない。

export と重い projection/preflight/restore loops は task/paint/input の機会を持つ時間予算で分割する。
byte count と `nextTick()` / microtask だけを応答性の証拠にしない。
JSON/gzip equivalence、gzip magic、fatal UTF-8、200 MiB file/512 MiB expanded limit、
50 MiB compressed download confirmation と既存 unsupported-browser error を維持する。
DOM SVG sanitization/mount が長い場合も工程別に測り、codec だけの改善を全体完了として報告しない。

## SOLID / KISS / DRY / YAGNI の適用

| 原則 | Code / architecture / workflow の具体的制約 |
| --- | --- |
| SRP | source lifecycle、parser、request projection、codec、adoption、3 Product concerns を別責任にする。authority/evidence/runtime の delivery も分ける。 |
| OCP | 複数 source は既存 typed records の data。新 renderer/argv/CLI surface を増やさない。 |
| LSP | Native File と restored file view は既存 content boundary を通す。同名や同 payload でも source instances を置換可能な一つとして潰さない。 |
| ISP | parser は metadata、renderer は request/resources、codec は plain JSON/bytes。全 Vue state を引数や Worker message にしない。 |
| DIP | source owner は text-reader/helper の既存境界に依存。setup/template は composition/adapters。低水準 transport に Product 判断を埋めない。 |
| KISS | 順次 sessions、一 shared branch、一 import transport、append/target replace/remove。独立 clones で共通 checkout の lock 管理を不要にする。 |
| DRY | fresh/restore projection、identity mapping、availability、writer/reader normalization は各一 owner。session 共通規約は一ファイルを参照する。 |
| YAGNI | multi-GFF pairing/reorder/subset/exportWorker/general RPC を追加しない。性能 evidence のない改修と公開ページ増設をしない。 |

## Authority と review の順序

S00 は approved receipt を既存 authority へ統合するための inert patch を準備する。
この implementation branch に active policy/static authority/checker を混在させない。
Product Contract amendment、privileged preauthorization、必要な mapped evidence 更新は、
maintainer が separate authority-only/evidence-only PR として dev に merge する。
その PR branch は最新 dev から作り、implementation branch の docs/runtime commits を含めない。
この task の同名 implementation branch push の許可を、別 target の push/PR merge 許可に拡張しない。
準備済み patch と exact target で必要な外部手続きだけを行う。

S00 の evidence/patch はこの branch に commit/pushする。S01 の independent measurements は進められる。
runtime session 開始時は `git fetch origin`、既存 authority の dev merge SHA を確認し、
clean implementation checkout で `git merge --no-edit origin/dev` して取り込む。
未反映なら依存 runtime を停止し、完了済み independent evidence を commit/pushして具体的未反映 target を報告する。
承認済み Product を再度問い直したり、candidate authority だけで runtime を開始したりしない。

new Worker constructor/importer、source lifecycle operator が trusted policy に含まれるか S00/S01 で調べる。
必要な preauthorization は実際の検出 subjects と最小の paths に限る。
checker mechanics を変える必要がある場合は checker-only → authority-only → runtime の順序。
全 future paths の一括 permission、policy mode 緩和、detector narrowing は使わない。

bindings reader で new persisted compatibility path が増えるなら Architecture Ratchet の例外条件。
OE/PE 改善で CB 増加を相殺しない。Product receipt から architecture exception approval を推定しない。
S02 で complete before/after sets と positive released fixtures を記録、S08 の final exact head で
maintainer 用の fully populated review packet を用意する。hard gate failure は別 authority が適切に merged するまで残る。
最終 approval は final push 後の exact head に対して行い、agent が承認 comment を投稿しない。

## セッション順序・成果物

| Session | 責任 | 前提 | 終了条件 |
| --- | --- | --- | --- |
| S00 | approved Product serialization、authority/privileged intake | 計画 branch取得 | 3 selected outcomesの候補 patch、外部merge対象と required paths、authority未反映状態を明記してcommit/push |
| S01 | baseline、real large Session、transport feasibility | S00のintake | input/browser/heap/heartbeat/stage evidence、一 transport の理由、budgets、main/tag namespace evidenceをcommit/push |
| S02 | Circular source collection/request/bindings/History を整合した一変更で実装 | 必要な authority merged、S00/S01 | C-01〜C-05、source lifecycle old path除去、released reader/instance identity、draft round trip、architecture sets |
| S03 | disclosure/status/deferred/focus/controls | S02、discovery authority merged | D-01〜D-04、single/grid/batch、native/helper settlement、keyboard/390 px |
| S04 | semantic-operation exclusivity と全 action availability | S03、session policy authority merged | S-02のbusy/cross-operation/late mutation、same-document 保存、rollback、read-only navigation |
| S05 | selected import Worker/transport/lifecycle | S04、S01 transport、privileged authority merged | main parse削除、JSON/gzip/legacy/settings-only、Worker failure/stale/limits、transfer metrics |
| S06 | measured heavy save/preflight/restore work の分割 | S05 | S-01〜S-05の大規模性能・heap/content equivalence、streaming export維持 |
| S07 | branch artifacts再生成、既存 public docs、公開 workflow検証 | S06 | current writerへの更新、旧 branch-only paths除去、generator evidence、screenshots/docsの再現性 |
| S08 | 統合・final exact-head review packet | S07、全 results | focused/full required gates、architecture/product review、残件解消、final commit/push、immutable head packet |

S02 の source/request/writer/UI を異なる intermediate formats の sessions に分割しない。
新writerで読めなくなるbranch-owned Session artifactsがrequired testsの入力なら、S02で同じowner generatorを使って
current formatへ再生成し、そのcommitに含める。S07まで旧bytesを残してtestsをskipしたり一時readerを足したりしない。
S07はその有効evidenceをreuseし、残るpublic docs/capturesとartifact inventoryを完成させる。
S00 の結果が未mergeでも S01 は独立して完了できる。通常は上記順序で、一 session が次の仕事まで実行しない。
各結果は `results/Sxx_RESULT.md` と必要な data/logs。未測定項目は未測定と記し、TODOを完了扱いにしない。

## Acceptance catalog

| ID | 必須の検証 |
| --- | --- |
| C-01 | 2 files/2+3 records →2 cards/5 records。単一 multi-record file →1 card。multi-selectとseparate appendが同じordered universe。 |
| C-02 | 同名、duplicate accession、同 bytes再追加でinstancesを区別。native/restored、Generate、fresh Load、CLI/Python replayで正しいrecord endpoints。 |
| C-03 | add/replace/remove/Undo/Redoは対象sourceだけ。対象外rows/transforms/depth sparse cells/series/selectorを保つ。late discoveryを拒否。 |
| C-04 | fresh shared canvas、explicit saved opt-out、single/grid/batch、one-record grid、batch prefixes、single crop。missing source comparison intentをdiscardしない。 |
| C-05 | before-first-Generate draft、active/inactive modes、committedと異なるdraft、settings-onlyのSave/Load。payloadは二重保存しない。 |
| D-01 | native GenBank/DDBJ、GFF+FASTA選択だけで探索開始。fast-pathはPython Worker 0。helperはoperation settlementを待つ。 |
| D-02 | one/two/duplicate IDs、incomplete pair、invalid、rapid replace/remove、mode/History restoreのtruthful statuses/predicates。 |
| D-03 | preview LoadはPython Worker 0。deferredとloadingを区別。Inspect/Generate/native replacementへ継続。別draftのmetadata混同なし。 |
| D-04 | applicable singleの自動展開、manual close保持、grid/batch理由/導線。focus/scroll anchor、keyboard、390 px。自動topology/selection変更なし。 |
| S-01 | JSON/gzip、supported historical、settings-only、draft/committed mismatch、existing validation/adoption/replay。 |
| S-02 | busy、duplicate Save join、cross-operation、Generate/reflow中、programmatic action、late file/cache/editor completion、取消、crash、unsafe keys、limits、rollback。旧request/resources/Result/Historyを保持。 |
| S-03 | real 10+ records/25,000+ features/full pairwiseと既存VibrioでSave/fresh Loadのstage/wall times、100 ms heartbeat、long tasks、main+Worker/process peak memory、transfer/copy bytes。 |
| S-04 | request/resources/catalog/cache/overrides/manifest/saved SVGの意味を比較。fresh Load → Generate → CLI/Python replay。外部sequence送信なし。 |
| S-05 | 既存budgetsを弱めない。担当工程の目標heartbeat p95≤250 ms/max≤500 msはbaseline実測後にfixture/environmentで固定。DOM mountの未解消停止も報告し、全体未達なら完了扱いしない。 |
| A-01 | single semantic owner、one canonical request/Result path、first-party cyclesなし、privileged deltaはbase authorityで許可、OE/PE/CB scopeと除去pathsをreview。 |
| W-01 | 全sessionが最新同名remote branchを独立checkoutで使い、担当範囲のみcommit/push。failed pushはremote stateを確認してから回復。 |

real large fixture は pinned source IDs/checksums と再生成 recipe を持つ。単なる複製を10 biological genomesの証拠と呼ばない。
synthetic stress fixture は補助に使えるが、S-03 の real dataset の代わりにはならない。
benchmark artifact にsource sequence/全 comparison rowsをログしない。入力fileの扱いは既存privacyに従う。

## Verification と公開資料

各sessionはfocused assertionsとrequired gates。long testsは30分以上を許し、incremental monitoring。
Node Playwrightがない場合もPython pathを調べる。Chromium sandbox failureは同じcheckを適切にescalateする。
browser wheelが必要なら `python tools/prepare_browser_wheel.py`。生成wheelをcommitしない。
主要commandsとsession別targetsは各prompt/common workflowに記載する。

参考出力は通常testsでread-only。意図的geometry変更がなければ更新しない。
必要なgenerated Session/ Gallery artifactsはowner generatorで再生成。social previewを変更しない。
既存 public technical/input/session explanation を更新し、workflowごとに新ページを増やさない。
procedural docs/screenshotsを実際に編集するS07で該当skillsを読む。internal planの作成だけでdocs capture workflowを実行しない。

final acceptance は production/tests/docs/generated diffs を別々にreviewし、未測定/失敗/authority未統合を区別する。
S08のfinal source/evidence commitをpushした後、exact head SHAを入れたarchitecture review packetをrepo外のMarkdownへ生成する。
commit自身のSHAをそのcommit内文書に書こうとせず、レビュー後のdocs commitでheadを変えない。
maintainerがそのheadへ必要なexceptionを手動判断する。PR作成/merge/deployはこの計画書pushとは別の外部操作。
