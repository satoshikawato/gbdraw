# INSTRUCTION PROMPT — S08: 統合・実待ち時間・回帰検証と最終引き継ぎ

S07.5からS07.8までの追加削減を終了し、S08を実行してください。
調査や計画だけで止めず、採用候補を固定し、実経路での検証、必要最小限の測定、
検出した範囲内の不具合修正、最終handoffまで進めてください。
ユーザーは「今回『もっと速くならないのか』っていったやつはもうやらなくていい」と明示しました。
追加高速化は終了案件です。未完了の宿題・S08の残課題・完了条件として持ち越さないでください。
新しい局所最適化やcache/runtime再設計を始めるセッションではありません。

## 1. 所在と引き継ぎ

Repository:
`/mnt/c/Users/genom/GitHub/gbdraw`

直前worktree:
`/mnt/c/Users/genom/GitHub/gbdraw/.worktrees/collinear-s07-8-20260916`

直前branch:
`perf/collinear-s07-8-local-reductions-20260916`

直前HEAD / S07.7比較基準:
`02ca8f950a7561ae9569c66f3f91d30702b53a17`

**S07.8のproduction、tests、handoff、測定証拠とこのプロンプトは未コミットです。**
HEADだけを引き継いでS07.8を引き継いだことにしないでください。

S07.8最終production SHA-256:
- `gbdraw/analysis/protein_colinearity.py`:
  `816111a76a99755034a3e4e5479230ba6f643cf810c576503f5d0dba267f453a`
- `gbdraw/analysis/collinearity_units.py`:
  `730c1354e98fcc539d7d395cc0c06f1fb3bbd4f71f57f9be211e2d9c8dc7ad8f`

開始時に共有repositoryと直前worktreeのbranch/HEAD/upstream/dirtyを確認してください。
originをfetchし、最新origin/devからupstreamなしの専用branch/worktreeを.worktrees配下に
作成してください。未統合依存をcommit/patch単位で照合し、同等patchを重複適用しないこと。
S07.8未コミット差分は、production・tests・文書・証拠を区別し、hash照合して新worktreeへ
一度だけ引き継いでください。generated wheelやvenvをsourceとしてコピーしないでください。
共有dev、既存worktreeを編集せず、stash/reset/cleanで既存変更を隠さないでください。

最新devとの衝突は現在の契約を保って解消し、S08で検証する統合sourceと依存差分を記録すること。
S07.8は**検証候補として継承**し、性能承認済みとは扱わないでください。
追加の性能採用実験や高速化をS08の開始・完了条件にしてはいけません。
S08の変更判断は統合上の正確性・契約・実在する不具合に基づけてください。

## 2. 必読・確定済み判断

計画ディレクトリ:
`docs/internal/collinear_similarity_performance_plan_2026-09-15/`

- repositoryのAGENTS.md、CLAUDE.md、gbdraw/web/CLAUDE.md。
- 最新architecture/Product ratchet、該当behavior contracts、Session互換契約。
- MASTER_PLAN.mdと本SESSION_08_INTEGRATION.md。
- results/S01.md〜S04.md、S06.md、S07.mdの実在するhandoff、S02のPATH-B判断受領記録。
- results/S07_6.md、S07_7.md、S07_LOCAL_REDUCTIONS.md、S07_8.md。
- results/data/s07-8-start.json、s07-8-stages.json、s07-8-review.json、
  s07-8-baseline-reuse.json、s07-8-benchmark-summary.json。
- 必要な既存runner、browser/CLI証拠と元report。リンク先全証拠を無差別に再実行しないこと。

判断は以下を維持してください。
- S07はユーザー承認済み。過去の合格・留保を変更しない。
- S06の全path列挙除去の主目的は了承済み。将来merge前の互換性Reviewは残る。
  この留保やGalleryの一律高速化をS08開始の妨げにしない。
- S05は却下済み。共通prepared class、cache owner、共有LRU置換、Worker transaction、
  scheduler/runtime改修を復活させない。
- PATH-B、推論ON/OFF、limits、既存のraw再利用・mode別limit記憶は確定済み。
  再質問せず、最新baseの実装とauthorityを照合して回帰対象にする。
- 新しいmaterialな挙動変更が必要な場合だけ該当Product判断に戻す。

S07.8の正確性・browser/CLI検証は完了していますが、実用的な時間・メモリ改善は未確立です。
旧7/21回と新3回のpolicy差で比較guardは全6条件を拒否しています。
Hepの保存値との悪化をコード起因と断定せず、Vibrioの改善値も性能合格へ変換しないでください。
約93%/91%という過去の累積観測を、厳密な全workflow改善率として使わないでください。

ユーザーは追加削減を終了しました。coverage union、数値検証の配列化、top-2化等の
新案を実装せず、その実現やS07.8の性能合格をS08完了の必須条件に追加しないでください。

## 3. 担当範囲

- S03/S04/S06/S07とS07.6〜S07.8候補を含む統合sourceの正確性・回帰検証。
- CLI/Python/Webそれぞれの公開経路、保存互換、実browser lifecycle、offline動作。
- ユーザーが待つWeb時間の区間別確認と、既存native測定との境界の整理。
- 検出したin-scope不具合の最小修正と対応する回帰テスト。
- results/S08.md、MASTER_PLAN.md、必要な既存契約文書の更新と最終監査。

過去に性能改善を立証できなかった事実は測定履歴として記載するだけにしてください。
追加高速化の未達・未完了とは扱わず、その解消のために実装・原因追跡・測定を再開しないこと。
通常のS08計測で支配的処理が見つかった場合も、現状の観測として記録し、
新しい改善課題や次の最適化セッションを自動で作らないでください。

## 4. 統合検証

1. 基準と統合source、authority、依存、互換readerを照合する。
   科学的結果のoracleは契約・独立oracle・適切な凍結実装を用いる。
   latest devに正当な契約変更があれば旧sourceの無条件一致を優先しない。
2. HSP、sparse/dense support、lossless paths、strict conflict、exact count、merge順、
   singleton、unit全mode、全ID/順序を確認する。固定coreと拡張後membershipを混同しない。
   通常経路での全path materializationがないことと明示的tuple APIを区別する。
3. Similarity有限/無制限、Collinear adjacent/all・推論ON/OFF・unit auto/cds/locus、
   公開されたmember/anchor limits、multi-record sourcesを実経路で確認する。
   source job数、directional table数、LOSAT起動数を区別し、必要なself/reverse evidenceを保持する。
4. 実browserでGenerate、warm Generate、色/block/member/filter変更、向き・順序・source変更を確認。
   cache key、実行段階、provenance、結果が一致すること。49/64/81表の既存cache境界は
   現行契約の検証に使い、S05で却下したcache redesignの要求に変えない。
5. raw完了後cancel→member変更→retry、raw設定変更、Clear Cache、Session/History置換、
   stale完了、Worker再作成を確認する。既存契約どおりlast successful Resultが維持されること。
   mockだけで最終合格にせず、実Worker/解析を通し、必要なcancel境界を確実に観測する。
6. supported old Sessionのload/current save/fresh load/regenerate、typed resource・metadataを検証。
   両modeのCLI replayと、公開Python経路を確認する。native/Pyodideの一致は各runtime内で比較する。
7. 最終sourceからwheelを生成して中身のPython sourceを照合する。installed CLIをcheckout外から実行。
   browser-offline-qa skillを適用し、local assetsのみの実行、lazy initialization、繰り返し実行時の
   保持とrelease、共有Worker/cacheを使うCircular smokeを確認する。
   source SPAと配布用bundleのどちらを検証したか明記し、未確認の配布物へ結果を拡張しない。
8. 不具合を修正した場合はfocused gateから影響範囲を再確認し、測定sourceと最終sourceを一致させる。
   科学的意味、limits、public validation、float計算順、Session authorityを性能都合で変えない。

## 5. 測定 — 必要な境界だけ、各3回まで

まず既存reportの入力・明示設定・到達source・依存を照合し、再利用可能な結果を使ってください。
S01から全benchmarkをやり直さず、HEADやrunner全体の差だけで同じstageを再測定しないこと。
基準sourceは各測定の目的に合わせて実在handoffから特定する。
S07.7→S07.8比較とS03以前→統合版比較を混ぜないでください。

native代表caseはHep Collinear ON/OFF、Vibrio Collinear ON/OFF、Hep Similarity member=5/None。
既存の`tools/benchmark_protein_comparison.py`を使用し、無関係なstageを測らないこと。
同条件比較が必要なのに旧reportを再利用できない場合は、その理由と必要caseを先に固定し、
該当baselineと統合sourceを各3回まで測定する。比較guardを迂回しない。

S08での実browser時間の確認境界は、**Generate操作→最終SVGの画面反映**です。
これは統合版の現状確認であり、「もっと速く」の再検討や新たな秒単位短縮目標ではありません。
既存browser runner・計測hookを再利用し、必要ならそこへ最小限追加する。
新しい並行benchmark frameworkや常駐profiling機構は作らないでください。

- 代表caseと設定を測定前に固定する。少なくともHep Collinear ON/OFF、Hep Similarity、
  大きいVibrio Collinear ONを含め、正確性の全組合せを時間測定の全組合せへ膨らませない。
- raw未生成で実LOSAT検索を行う経路、保存raw再利用で再解析する経路、derived結果再利用のwarm経路を分ける。
  「新browser contextだが保存rawあり」をraw coldと呼ばない。
- runtime初期化・入力準備、raw検索、post-search、typed結果/転送、描画・DOM反映について、
  観測できる境界の所要秒数を記録する。並行処理や親子の累積時間を足し合わせない。
  分離できない区間は合算のまま明示する。
- 旧版との同条件browser比較がない場合、現在の内訳を報告し、改善率を捏造しない。
  現在の支配的区間の特定を、過去全revisionのbrowser再測定へ拡張しない。
- 各case・source・測定境界につき3 samplesまで、warmup最大1回。全sample、中央値、ばらつきを記録。
  cold/warmの状態を毎回規定どおり復元し、warmupによってcold測定をwarmへ変えない。
- noiseだけを理由に自動追加・7回・21回への増加をしない。3回で結論が出なければ未確立とする。
- profile/counter、tracemallocは時間測定と分離し、必要caseのみ各1回。
  Python allocation、RSS、browser/Wasm memory、retained/peakを混同しない。
- build、重いtests、自分の別benchmarkと時間測定を重ねない。外部processを停止・変更しない。
  競合時は記録し、長時間待たず独立作業へ進む。
- 時間は秒単位と短縮した絶対秒数で示す。処理回数削減、局所時間、native post-search、
  Web全体時間を分ける。新しい速度目標を捏造して未達の作業を増やさない。

## 6. 必須gate

focused testsから関連native、typed/Session、全体gateへ進む。
同じ最終sourceのfull suiteに含まれる検証を別commandで機械的に重複実行しないこと。

```bash
pytest tests/test_protein_colinearity.py tests/test_collinearity.py tests/test_collinearity_units.py -v
pytest tests/test_web_feature_catalog.py tests/test_session_request_codec.py tests/test_session_compat.py -v
pytest tests/ -v -m "not slow"
pytest tests/test_output_comparison.py::TestOutputComparison -v
ruff check gbdraw/
```

最新inventoryのcache、derived/raw identity、Session、Worker、architecture/ProductのNode gateも実行。
Python/Node両方のPlaywright導入経路を確認し、NodeがなければPython Playwrightを使用する。
sandbox制約の実browser起動失敗を未導入と取り違えず、同じcheckを必要な権限で再実行する。
各test commandには最低30分を許容し、incrementalに監視すること。
参照SVG、test-owned timeoutを変更して失敗を回避せず、pytestにwall-clock assertionを追加しない。

## 7. 最終成果物と終了条件

`results/S08.md`と`MASTER_PLAN.md`に以下を記録する。

- base/依存/未コミット差分の継承、最終source hash、wheel一致、環境、fixture、再現command。
- 正確性、実browser lifecycle、offline、Python/CLI、保存互換の結果と未実施・失敗。
- 同条件で立証した改善、過去の観測値、未確立の性能、code/host原因未分離を区別した結果。
- native/Webの時間境界、絶対秒数、転送量、memory、operation countと現状の時間内訳。
- S07.8候補を含む最終sourceと、統合上必要だった修正・取り下げとその理由。
  追加高速化の証明は不要。候補の継承を性能合格と表現しない。
- S07の既承認、S06の将来merge前Review、S05却下を維持した達成範囲。
- production / tests / benchmark tooling / documentation / generated evidenceを別々に監査。
- rollback、統合検証上の未実施・実在する不具合、English proposed commit titleと短いsummary。
  終了した追加高速化を残課題へ戻さない。

技術的な統合検証の完了と性能の立証範囲を分けて報告する。
必須gateの未実施・未解決を成功扱いしない一方、追加高速化を終了条件へ持ち込まないこと。
PATH-Bの通常経路の全列挙除去と、明示的全量取得に残る不可避なコストを区別する。
検証・必要な修正・記録が済んだらS08を終了し、際限なく追加改善へ進まないでください。

push、PR作成、統合merge、tag、deploy、外部への送信は認可しません。
共有ツリーや既存worktreeの変更、public showcaseの無関係な更新も行わないでください。
