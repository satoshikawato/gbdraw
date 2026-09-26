# S01 — Session responsiveness and import transfer evidence

日付: 2026-09-26–27。**S01計測完了、transport選定済み。full-pipeline性能はFAIL、runtime受入は未完了。** S01のみ、BUG-02 / BUG-20。Production runtime/defaults、active authority/checker の変更ゼロ。
BUG-01、S02、新 bindings、Circular source collection、authority PR、S03以降の実装は行っていない。

## Checkout / authority / source

専用 checkout `/tmp/gbdraw-issue597-S01.j3R30y`。
Branch `fix/issue-597-input-session-20260926`、upstream `origin/fix/issue-597-input-session-20260926`。
指定の single-branch clone、fetch、pull --ff-only、dev/main/tags 明示 fetch を完了し、開始時 tree は clean。
開始/測定 runtime HEAD `9584ed904c7d1ee93f5a07ea4f3a2a465b061724`。
Fetched dev / PR #610 merge `af5d942af60353dda199aa487da9152a3576b3fe` は dev と HEAD の ancestor。
Fetched main `4556e04e929a4a85ad28d1833ce7304bd764881c`。
[authority evidence](../evidence/S01-authority.json) は OIPC revision 21、PD-OI-044/045 と
human / serialized receipts 全18項目の一致、active OIPC と fetched dev の byte 一致を確認している。
承認済み Product 選択の再承認は不要。

共有 checkout・他 session の branch/environment/server は操作していない。
通常 sandbox は既知の bubblewrap `/mnt/wslg/distro` mount エラーで起動不能だったため、同じ専用 `/tmp` の commands を escalation で実行した。
依存追加や共有 cache の更新をせず、fixture の LOSAT cache と CLI replay cache は専用出力内に限定した。

Environment、source/input/recipe/authority SHA、raw artifacts の SHA と再現 commands は
[S01_REPRODUCTION.md](../evidence/S01_REPRODUCTION.md) と [S01-artifacts.json](../evidence/S01-artifacts.json) を参照。
生成 wheel SHA `e6e338b0c21e03aaa348c60e8314d61b83d262def1494bf85ce29a769aa3816a`。
cache bust は不変。ignored wheel は commit していない。

## Baseline discovery and namespace

[discovery](../evidence/S01-discovery.json) の native GenBank / DDBJ accession / GFF+FASTA は ready、Python Worker 0。
同一ファイル内の duplicate accession 二件は selector #1/#2 と二 controls を保持。
遅延 read の旧入力を新入力に置換した後、旧 settlement は新 record を上書きしない。
Remove は idle/0 records、不完全 GFF pair は idle、invalid native は error と旧 Result を保持。
既存 parser/helper の forced fast-path miss は既存 Python Worker で settlement まで待ち、正しい metadata を返す。
新 parser は追加していない。
保存済み Circular preview は Result 1 / Worker 0 だが discovery loading/0 records のまま。
native replacement は active draft を ready にし、旧 committed request / Result は保存する。
fresh settings-only v42 は idle / Result 0 / Worker 0。page error / external requests はゼロ。
補助 JSON/gzip の Unicode/lone-surrogate/own-key transport 比較は意味一致。
これは real biological fixture の代用や unsafe-key admission の承認ではない。

[namespace evidence](../evidence/S01-namespaces-and-source.json):
main は Session 42/request 7/bindings 2/catalog 3、HEAD は 44/8/2/4。
protein raw cache 4/nucleotide 2/derived 3/identity manifest 2 は main と HEAD で同一。
tag 0.13.0 は Session 30、現在の split namespace modules は存在しない。
v39/v41-bindings1 の witness は main first-parent。
settings-only v42 の元 witness は branch commit だが、同じ fixture bytes の main merge witness `3fd508410d492e20d9963e41d808e9e7078d4d2c` を確認。
v44/schema7 の witness は main first-parent/tag 上になく、main の同一 fixture witness も存在しない。
従って「released」とは扱わず branch-owned current-format evidence として記録した。
S01 は既存 reader/migration/公開説明を変更せず、新 compatibility path を増やしていない。

## Fixtures / measurements

[fixture manifest](../evidence/S01-real-fixture.json): real 12 chromosomes / six assemblies、62,089 biological catalog features（元sourceのtranslated CDS 29,858）。
全132 non-self directed pairs と全12 self pairs を集合一致で確認、raw entries 144。
一つの combined GenBank file を既存 Linear records-table CLI に渡して生成した既存 writer 44/request 8/bindings 2/catalog 4。
compressed 134,471,293 bytes / expanded 466,767,621 bytes。既存200/512 MiB上限内。
Fixture SHA-256 `a60cce6ecb8a0c499f2d801feabc65b122c78793dafcc15644b46c6ef1ecc98f`。
初回coverage assertionはself pairsを誤って除外する件数期待で失敗した。実データを診断し、必要な132組の完全な集合と12自己比較をそれぞれ厳密に検証し直した。CLI成果物の捏造やsearch/cache削除はしていない。

最終 run は生成完了後、pipeline → transport の順。各dataset三回、各回 fresh browser context。
[pipeline metrics](../evidence/S01-final-pipeline.json) / [transport metrics](../evidence/S01-final-transport.json) に全反復の工程時間、bytes、longtasks、main/Worker/process memory、欠測を保存した。

| Dataset / phase | Wall ms | Heartbeat p95 ms | Heartbeat max ms | Long task max ms | 250/500目標 |
| --- | ---: | ---: | ---: | ---: | --- |
| vibrio load | 12,879–16,531 | 12,685.6–16,067.8 | 12,685.6–16,067.8 | 12,451–16,040 | FAIL (0/3) |
| vibrio save | 4,330–8,082 | 186.1–335.9 | 188.2–496.8 | 174–312 | FAIL (2/3) |
| real load | 111,364–129,825 | 100.2–100.2 | 95,246.3–107,418.2 | 95,186–107,358 | FAIL (0/3) |
| real save | 21,728–25,771 | 178.9–204.9 | 1,798.1–2,188.3 | 1,621–1,943 | FAIL (0/3) |
| vibrio whole only | 696–718 | 107.2–119.3 | 107.2–119.3 | 82–92 | PASS (3/3), codec-only |
| vibrio bounded only | 2,749–2,806 | 101.2–101.6 | 101.6–102.4 | 0–0 | PASS (3/3), codec-only |
| real whole only | 2,604–3,088 | 100.2–103.1 | 396.8–410.3 | 335–394 | PASS (3/3), codec-only |
| real bounded only | 12,158–14,105 | 101.4–102.3 | 102.0–145.9 | 0–0 | PASS (3/3), codec-only |

全 pipeline は **FAIL**。real Loadはp95約100 msでも、最大95–107秒の停止が残るためPASSではない。
主工程は real preflight 7,400–11,018 ms（catalog validation 5,120–7,493 ms、raw-cache validation 2,079–3,178 msを内包）、restore gap 7,434–10,896 ms、SVG admission 920–1,315 ms、DOM mount 87,687–98,544 ms。
Vibrio DOM mountも8,876–11,939 ms。real Save projection 378–472 ms / compression 14,952–18,010 ms。
その他のread/decompress/parse/sanitize/encode CPU、nested stagesとobserver spansはJSONに残している。
Vibrio Saveの既存数値budgetは三回すべてPASS。wall最大8,082 ms、observed heap delta最大184,484,549 bytes、heartbeat max最大496.8 ms。追加p95目標は三回目335.9 msでFAIL。元Node Playwright UI specの全assertionsをPASSと呼んでいない。
Load decode bytes: Vibrio 95,989,532 / real 466,767,621。Save re-encode bytes: Vibrio 96,041,953 / real 466,918,752。バイト数は実際にobserverが見たencode/decode出力で、未知の内部copyを0と推定しない。

| Dataset / phase | Main precise observed peak MiB | Worker CDP used / backing peak MiB | Dedicated Chromium RSS peak MiB |
| --- | ---: | ---: | ---: |
| vibrio load | 581.2 | N/A (Worker 0) | 1,461.8 |
| vibrio save | 549.7 | N/A (Worker 0) | 1,360.6 |
| real load | 2,075.0 | 22.5 / 102.1 | 3,792.2 |
| real save | 1,703.2 | 12.6 / 76.5 | 3,068.8 |

real Loadの既存 diagram-generation/Python helper Workerは各fresh contextで一回起動した。Saveは同じcumulative constructor一件を引き継ぎ、新たなconstructor増分0。これは新JS import Workerではなく、preview Worker-free条件の現baseline未達。
`ready.records=1`はsingle multi-record fileのLinear UI entry数、`requestRecords=12`がcanonical record数。両者を混同しない。
busy LoadのCDP timeout、transport終端のclosed-Worker connectionもraw/compact JSONに保存。main/Worker usedSize/backingStorageSizeとprecise/RSSは別定義のsampled lower bounds。完全な瞬間peakや欠測ゼロの証明は未完了。

工程定義、観測区間、sampling / copy / wire-byte の限界は reproduction 文書に固定した。
full pipeline と disposable transport-only は別の判定。stream spans と nested preflight stages は加算しない。
main と Worker の sampled heap、外部 backing storage、Chromium process RSS は定義を分けて報告する。
観測された peak は lower bound。busy-main CDP timeout をゼロや PASS と扱わない。
既存 Vibrio numerical budgets は wall <22,341 ms、observed heap delta <1,898,125,842 bytes、
heartbeat max <1,000 ms または ≤2,875.8 ms のまま。p95≤250/max≤500 ms の追加目標も緩和しない。

初期 BLAST+ generation は途中で自身の process のみを停止し、completed fixture と扱わない。
初回 LOSAT generation は SVG 完了 / Session 未保存で中断し、exit status は取得不能。成功とは扱わない。
初期の narrower codec observer による Vibrio/transport probes は supplemental evidence に限る。
最終判定には complete fixture と final observer による三反復だけを用いる。

## Selected transport and inert permission

**File/Blob → JS Worker read/decompress/JSON.parse → whole-object reply を一方式として選定。**
二dataset・各三回で完全payload semantic hash一致、p95≤250/max≤500を満たした。
real whole wall 2,604–3,088 ms、transfer/assembly 713–812 ms、main longtask 335–394 msを観測。
longtasksはゼロではないが、指定した目標内。Native clone wire/copy bytesは未測定（null）、転送buffer bytesは実測0。
observed real main precise peak約446.8 MiB、Worker CDP used peak約865 MiB、process RSS約2.20 GiB。

bounded比較も目標内、全payload一致だが、real data messages 69,520、UTF-16転送buffer計695,553,402 bytes、clone conservative upper estimate468,463,253 bytes。
データ/転送buffer一messageのupper最大262,144 bytes（metadataを除く）。real wall 12,158–14,105 ms、transfer/assembly 10,202–11,946 ms。
main precise observed peak約789 MiB、Worker CDP used約895 MiB、process RSS約2.24 GiBで、wholeより負担が増えた。
従って通常runtimeにbounded/wholeのsize fallbackや汎用RPCを持ち込まず、whole一方式を実装入口にする。
120 ms sampling pause +120 ms終端flush、ACK等はdisposable observerのための制御で、runtimeの要件ではない。
Workerはreply settlementで終了し、preflight/restore/DOM workに並行して保持しない。
Worker replyはplain未信頼candidate。main全量text/JSON.parse、Vue/DOM proxy、codec fallbackは採用しない。
メモリ値は観測値であり、termination付近のCDP欠測を完全PASSと扱っていない。final実装のmemory/latency再検証はS05/S06に残す。

constructor / client: `services/session-import-client.js`。
Worker: `workers/session-import-worker.js`。codec owner: 既存 `services/session-file.js`。
予定の一方向は config → client → Worker → codec。
codec は client を逆 import しない。config は既存 Save codec import を保持する。
validation、projection、adoption、History、request、SVG owner は既存のままで、Worker が複製・信頼済み化しない。
client が operation ID、error/stale rejection、settlement、termination を所有する。

[permission candidate](../authority-candidates/02-import-worker-permission.patch) は inert のみ更新した。
既存 detector で constructor 一件と Worker → session-file の privileged importer 一件が必要。
候補は policy 一 target、operator permissions 48→49 / importer 54→55、他の許可は不変。
S00 の zero-importer probe を実際の codec 共用 edge に更新し、コピーへの `--inert` 検証と `git apply --check` が PASS。
active policy/checker は source/fetched dev と同一。別 permission-only authority PR の dev merge は未実施。
詳細と将来の別 delivery 手順は [candidate README](../authority-candidates/README.md)。

## Verification / review / acceptance

[semantic evidence](../evidence/S01-semantics.json) と [verification exits](../evidence/S01-verifications.json) に実際の結果を保存した。
全saved六件で request/resources/manifest/derived cache、全source editor fields、catalog、raw cache source fields、saved SVG semantics が一致。
通常semantic verifierは六件ともcore九項目PASS。Vibrio/real各一件でPython reader records 4/12、CLI exit0、SVG 5,128,555/19,621,657 bytes。
fresh browser Loadもcanonical request、catalog/cache、preview readinessは一致し、page error/external requestsは0。
ただしreal fresh-Load verifierは既存Python/helper Worker一件のため **exit1、14項目PASS / Worker-free項目FAIL**。
この失敗はassertionを弱めず残した。runtimeで解消するのはS03/S04/S05以後であり、S01で「全受入PASS」と報告しない。

| Command / check | Result |
| --- | --- |
| `python tools/inspect_issue597_s01_contracts.py` | PASS: merge ancestry、revision21、18/18 receipt fields、active contract bytes。 |
| `python tools/inspect_issue597_s01_namespaces.py` | PASS: namespaceごとのsource/blob SHAとmain/tag/fixture witnesses。unreleased witnessはPASS扱いしない。 |
| `python tools/characterize_issue597_s01.py --output <task-dir>/discovery.json` | PASS: actual baseline assertions、JSON/gzip auxiliary equivalence、page/network errors0。 |
| `python tools/prepare_issue597_s01_fixture.py --output <task-dir>/fixture --threads 8 --verify-only` | PASS after diagnosed coverage-test correction: real sources/unique sequences/counts/full exact pair coverage/schema/200/512 MiB limits。 |
| `python tools/measure_issue597_s01.py --mode pipeline --repetitions 3 --fixture <real.gz> --output <task-dir>/pipeline` | exit0、6 completed observations。既存Vibrio数値budget PASS、追加full-pipeline目標FAIL。 |
| 同recipe `--mode transport` | exit0、各方式×二dataset×三回=12 completed comparisons、payload equivalenceとlatency targets PASS。memory終端欠測は別記。 |
| `python tools/verify_issue597_s01_outputs.py ... [--replay --fresh-load]` | core semantic全六件PASS。replay両dataset PASS。real Worker-free assertionは上記FAIL / exit1。 |
| 共通規約の五Node files | 30 PASS、0 skip。`/tmp/issue597-s01-focused-node-final-j3R30y.log`。 |
| `node --test tests/web/session-resource-backing.test.mjs tests/web/settings-only-session.test.mjs tests/web/session-losat-cache-validation.test.mjs tests/web/session-feature-metadata.test.mjs` | 63 PASS、0 skip。`/tmp/issue597-s01-extra-node-j3R30y.log`。 |
| `pytest tests/test_record_metadata.py` | 14 PASS、0 skip。`/tmp/issue597-s01-python-tests-j3R30y.log`。 |
| `ruff check gbdraw/` | PASS、unchanged production evidenceを再利用。 |
| six recipes + inert verifier Ruff、six recipe format、Python AST / JS syntax | PASS。 |
| scoped Markdown links / Bash blocks / diff checks | PASS。final data/doc更新後にもscoped確認。 |
| inert `git apply --check` / `verify_candidate.py privileged "$PWD" --inert` | PASS、active policyに適用せずコピーでshape/detectorを確認。 |
| trusted-base local `node tools/check-web-change-budget.mjs --base origin/dev` | PASS / Review CLEAR（staged担当成果物をorigin/devと比較）。committed-head gateはcommit後に再確認。trusted CIの結果と混同しない。 |

Node `@playwright/test`未導入のため元のNode browser specsは実行していない。Pythonによる上記targeted browser/measurementを実行した。
全pytest、whole historical browser matrix、fresh Load→Generate、rollback/cancel/crash/busy UI全assertionsは未実施。production変更がない範囲で無目的にbroad testsを繰り返していない。

Production diff: gbdraw runtime と public defaults の変更ゼロ。
Tests/recipes: tools 六 recipes、tests/web/fixtures 二 disposable probes のみ。app routing から import されない。
Docs/evidence: S01 結果、計測値/fingerprints/再現手順、必要な inert permission 候補のみ。
Generated: Gallery、reference outputs、social preview、dist/egg-info、tracked wheel の変更ゼロ。
OE/PE/CB: runtime owner/path/compatibility delta はゼロ。候補の実装前 permissions は別 authority delivery。
例外 sets や新 reader、parallel runtime fallback、general RPC framework は追加していない。

| Acceptance | S01 状態 |
| --- | --- |
| D-01/D-02 | 不足する native/helper、duplicate/replace/remove/incomplete/invalid baseline を固定。新 statuses 実装の受入は未実施。 |
| D-03 | saved preview / draft / settings-only の差を観測。loading の未解消を記録。deferred/Inspect/Generate 実装は未着手。 |
| D-04 | 自動展開、manual-close、390 px/focus/keyboard は未実施。S03 の条件を維持。 |
| S-01 | 既存形式 JSON/gzip/settings-only/namespace を調査。全 historical browser matrix は未実施。 |
| S-02 | exclusive semantic-operation、busy/cross-operation、cancel/crash/rollback の runtime 受入は未実施。 |
| S-03 | real/Vibrio各三回の全工程・codec比較を完了。latencyとmemory欠測を別記。完全peak保証は未完了 |
| S-04 | 全六savedの意味一致と両Python/CLI replay、fresh request/catalog/cache/preview一致PASS。real Worker-free assertionはFAIL。fresh Load→Generate の Generate は未測定、S05/S06 で必要。 |
| S-05 | whole transport-onlyは目標内。両Loadとreal Save、Vibrio Save一回は全体目標FAIL。transport-only 成立と full-pipeline 完了を混同しない。 |
| A-01 | runtime delta 0、既存 detector と最小 inert paths を確認。permission dev merge は待ち。 |
| W-01 | 専用 checkout・明示 stage・一 commit・同名 non-force push と remote/local 一致を完了して終了する。 |

## Remaining conditions / next entry

Product merge は確認済み。Worker permission の別 authority PR review/dev merge SHA、implementation への取り込みと gate PASS が S05 の開始条件。
S01 方式が成立しなければ S05 の import transport 実装は停止する。
S06 は measured bottlenecks を task/paint/input に機会を持つよう分割し、同じ fixtures/environment の full-pipeline target を満たす必要がある。
DOM mount の停止、memory sample gaps、unexecuted acceptance を codec success で閉じない。

次 session は [S03_INSTRUCTION_PROMPT.md](../sessions/S03_INSTRUCTION_PROMPT.md) と本結果を読む。
独立 checkout を取得し、Product ancestry/receipt と source drift を再確認してから、その session の scope だけを開始する。
S01 は S03以降や authority PR に進まず終了する。

English commit title: `Measure issue 597 session responsiveness and transfer costs`

English summary: Add reproducible biological fixtures and isolated browser probes for discovery, Session stage costs and import transport. Preserve measured failures and record the minimal pending Worker permission without changing runtime or active authority.
