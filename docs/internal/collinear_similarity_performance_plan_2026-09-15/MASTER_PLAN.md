# Collinear / Similarity groups 計算量改修 — 総合計画書

- 作成日: 2026-09-15
- 状態: S01の基準・再現資産とS02の調査・テスト専用prototype・同値検証・設計が完了。S02のProduct選択はPATH-Bとして受領済み（decision revision 1、authority文書化済み、review・base統合は未完了）。S03のHSP集計実装・同値検証・広いgate・観測測定と引き継ぎが完了。独立再測定はユーザーの最終指示で省略し、他benchmarkとの重複を測定条件の未達として記録。S04の疎な所属探索・metadata索引の本番改修、同値検証、広いgate、測定と引き継ぎが完了。最終時間判定は22 pass / 2 inconclusive / 0 regressionで、時間gateの全項目合格ではない。S05以降の本番改修は未実施。
- 対象: Similarity groups（内部トークン `orthogroup`）と Collinear の解析、Web 中間キャッシュ、結果・保存メタデータ。

## 1. 目的と完了の意味

共通の解析処理を再利用できる構造に整理し、総当たりを evidence の索引参照へ置き換え、全経路の列挙を通常の解析・描画・保存から切り離す。

目標は次の四つである。

1. 同じ入力を繰り返し転送・検証・parse・集計しない。
2. 疎な入力の所属候補探索を、全タンパク質と全グループの積から実在 evidence に基づく処理へ変える。
3. 有限のキャッシュ予算で、設定変更時に再利用できる解析段階を保持する。
4. 経路情報を失わないグラフ表現によって、通常処理が全経路数に比例して増える問題を解消する。

2026-09-15にProduct Decision OwnerがPATH-B / lossless-graphを選択した。
根拠は[完全なhuman receiptとauthority引き継ぎ](results/S02_PRODUCT_DECISION.md)であり、
この計画書の推奨ではない。authorityのbase統合前に依存runtimeを実装しない。

## 2. 調査の基準と更新手順

### 2.1 計画時に確認した状態

- 発端: 共有作業ツリーの未追跡 `docs/internal/report.md`（このworktreeには含めない）。監査対象は旧作業ツリー `51d2ae362ac154a3e360727ad086ad1d27bbc989` と未コミット変更。
- 文書保存時に fetch して確認した `origin/dev`: `65f231af175c0dbbbbf9b4674566e84ce5dddac5`。
- 現在の共有作業ツリーには対象・対象外の未コミット変更がある。旧作業ツリーをそのまま実装ベースにしない。
- report の `/tmp/gbdraw-mode-audit/` は計画時に存在しなかった。再現コードを利用可能と仮定しない。
- report の秒数、239 passed / 1 skipped、Node import 失敗は過去の監査結果であり、新しいベースの測定・テスト結果ではない。

| report の論点 | ベース確認結果 | 扱い |
|---|---|---|
| 全 manifest のジョブごとの転送・検証 | `47f5cebd` でキー生成をバッチ化済み | 再実装せず、検証回数・キー・転送量を再確認 |
| HSP ペアごとの DataFrame / 行イテレータ処理 | 構造が残存 | S03 |
| 全経路の列挙 | 構造が残存 | S02 で契約調査、S06 で承認済み結果を実装 |
| 未所属タンパク質と全 core group の総当たり | 構造が残存 | S04 |
| filtered / converted の64件共有LRU | 構造が残存 | S05 |
| 不要な表示 projection、メタデータ線形検索 | 改善余地あり。selector に `comparison_pairs=()` が追加済み | S04・S05で既存境界を利用 |
| Collinear cluster 結合 | 構造が残存 | S07、再測定で必要性を判定 |
| downstream cancel 後の raw 検索再実行 | `7db0539a` で対応済み | 回帰を防ぐ |
| optional Collinear inference、mode別limit記憶 | `65f231af` でruntimeも統合済み | ON/OFF・旧Session・mode切替を回帰検証 |

参照: [バッチキー生成](https://github.com/satoshikawato/gbdraw/commit/47f5cebd)、[キャンセル後の再利用](https://github.com/satoshikawato/gbdraw/commit/7db0539a)、[optional inference と limit 記憶](https://github.com/satoshikawato/gbdraw/commit/65f231af)。統合済みはコードの確認結果であり、本書作成時に当該runtimeのテストを実行したという意味ではない。

### 2.2 Product authority と別作業への依存

S01 の確認時には `origin/dev` が `9a4f7e29ab1b99676bb783321f8a9f40f969d04c`
へ進んでいた。追加差分は legacy LOSAT candidate の消費を成功時まで確定しない
retry 修正（PR #534）。計画済みの性能改修へ重複実装しない。
更新された入力・authority・consumer の一覧は
[S01 inventory](results/S01_INVENTORY.md)、実測と再現コマンドは
[S01 handoff](results/S01.md)を参照する。以後のセッションはこの基準との
差分を確認する。S02 の経路表現の選択と S03 以降の本番改修は別セッションである。

S02開始時のfetchでもbaseは`9a4f7e29`のままで、S01成果は未統合だったため
新規worktreeへcherry-pickした。調査・検証結果は[S02 handoff](results/S02.md)、
全consumerと公開履歴は[S02 inventory](results/S02_INVENTORY.md)、具体的な型・
順位・保存・互換性は[S02 contract](results/S02_PATH_CONTRACT.md)を参照する。
[Decision Pack](results/PATH_DECISION_PACK.md)は判断時の比較記録でありauthorityではない。
受領したPATH-Bの文書化・review・base統合状態は[判断引き継ぎ](results/S02_PRODUCT_DECISION.md)を参照する。
S06にはS05完了と、選択済みoutcomeを認可するbase authorityが引き続き必要。

S03開始時のfetchでもbaseは`9a4f7e29`。未統合のS01/S02成果を依存順に
新規worktreeへ引き継ぎ、authority-only commitは適用しなかった。
[S03 handoff](results/S03.md)に処理owner、同値結果、測定・環境条件、S04の
実在APIを記録する。S03の実装・検証・観測測定と引き継ぎは完了。ユーザーの
最終指示で独立再測定を省略したため、性能の数値判定24 passと、他benchmarkとの
重複による測定条件未達を分けて報告する。独立測定の性能合格とは扱わない。

S04開始時もfetch後のbaseは`9a4f7e29`。未統合のS01〜S03を一度ずつ引き継いだ
`d7695a62`を固定baselineとして、所属evidenceとmetadataの索引を実装した。
[S04 handoff](results/S04.md)にsnapshotの寿命、削除した走査、28 cases / 51 stagesの一致、
4,295 Python testsと214 Node tests、実browser、時間・操作回数・メモリの証拠を記録する。
最終の21-sample測定は22 pass / 2 inconclusive / 0 regression。sparse-200とunrelated-200は
MAD/中央値が5%を超え、数値上の性能合格にはしない。S03の省略指示や観測値は流用していない。
authority-only commitは未統合のままS04へ混ぜていない。S05は未実施。

[Option Integrity Product Contract revision 5](https://github.com/satoshikawato/gbdraw/blob/9e28581a3d0e7fb9117da0fbc11d46a945b0e59e/docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md) は、以下を既に選択している。

- `PD-OI-001/002/004`: Web の raw/member limit のモード別初期値・記憶、明示した有限値と無制限の維持。
- `PD-OI-018`: 必要な全レコードの比較、入力ソース単位の検索、実際の検索DB範囲を含む raw identity。
- `PD-OI-021`: Web Collinear の **Infer orthogroups with self-comparisons** は fresh/reset で OFF。ON は既存推論を維持し、古いセッションで値が欠ける場合は歴史的 ON を維持。
- `PD-OI-022`: 完了した raw 検索は downstream cancel 後も再利用できる。未完了 batch を完了扱いしない。

先行する口頭計画で未反映とした optional inference と limit retention は、文書保存時の `65f231af` でruntimeへ統合された。S01 はこの既存実装を含むbaseで検証する。性能改修内に別の UI / defaults / migration 実装を作らず、ONの推論とOFFの直接block構築をそれぞれ維持する。既に承認された outcome に、同じ Product 選択を再要求しない。

単純な一行一レコード・推論 ON の隣接比較は方向付き evidence が `3R-2`、Similarity groups / Collinear All records 推論 ON は `R²`。推論 OFF の Collinear は self を除く。実際の LOSAT invocation 数は source batching、選択ペア、互換な検索設定によって異なる。レコード比較数・ソースジョブ数・プロセス起動数を別々に記録する。

S01 および各実装開始時に最新の `origin/dev` を確認し、既に解決した項目、改名した owner、更新された契約をこの表と handoff に反映する。計画時の schema 番号や古い行番号を新規実装の根拠に固定しない。

## 3. 科学的結果と互換性の不変条件

- 選択された scope に必要な self / reverse / 非隣接 evidence を保つ。OFF で不要になる evidence は承認済み契約に従う。
- `max_hsps=1`、隠れた candidate/member cap、経路の打ち切りを性能対策として導入しない。
- HSP union coverage、代表 HSP の同点順位、ペア順・行順、無効入力の扱いを維持する。
- support / diagnostic evidence、domain-only、best/second-best、role、confidence、assignment reason、関連辺を維持する。
- 固定 core snapshot と、所属拡張後に record-local 競合判定が使う集合を区別する。
- block の anchor、singleton、merge 順、境界条件、orientation、ID を維持する。
- group/path ID、出力順、型、空値・省略・エラーの契約を同値最適化の対象に含める。
- raw identity に実際の検索DB範囲・`searchContext`・入力 bindings を保つ。表示設定変更による raw 再利用を壊さない。
- failed/canceled/stale Generate は最後に成功した Result と committed request を置換しない。
- 公開 API と各保存 namespace の変更は S02 の調査・Product 判断・互換性方針に従う。

現在のコード・テストは観測事実であり、自動的に Product authority にはならない。既存の明確な authority と食い違う期待値は、性能改修で凍結せず S01/S02 で分類する。

## 4. 目標アーキテクチャ

```mermaid
flowchart TD
    A[CLI / Python / Web 入力] --> B[既存の型付き要求・検索計画]
    B --> C[検証済み raw evidence]
    C --> D[フィルタ・member選択・必要な正規化]
    D --> E[共通 group inference]
    E --> F[所属・関係グラフ]
    F --> G[表示リンク / Collinear block projection]
    D -->|Collinear inference OFF| G
    G --> H[描画・metadata・保存]
```

この図は責務を表す。各箱のために新しい公開 class や framework を作る指示ではない。

| 責務 | 現在の主な owner / 方針 |
|---|---|
| protein identity、HSP 集計、所属判定 | `gbdraw/analysis/protein_colinearity.py`。Python に意味を集約 |
| block 形成 | `gbdraw/analysis/collinearity.py`。推論の有無で block algorithm を複製しない |
| 生物学的単位・座標 | 既存の unit / record planning owner を使用 |
| Web 検索の orchestration | `gbdraw/web/js/app/run-analysis.js` と既存 source batching / LOSAT services |
| helper transport | `python-helpers.js`、diagram worker protocol / worker。新しい推論ロジックは Python package へ置く |
| 中間キャッシュの寿命 | 一つの diagram Worker 内 owner。既存 parsed biological input cache と役割を混同しない |
| metadata / catalog | `gbdraw/web_support/orthogroup_metadata.py`、`feature_catalog.py` |
| 型付き保存 | `gbdraw/session_request_codec.py` と既存 resource/session owner |

両モードと各 surface は同じ解析実装を呼ぶ。キャッシュの有無で別の推論経路を持たない。私有関数への分解は必要な範囲にとどめ、移した旧実装は同じ変更で削除する。

### 4.1 原則の適用

| 原則 | 具体的な判断基準 |
|---|---|
| SRP | 推論、projection、キャッシュ寿命、保存形式の変更理由を分離 |
| OCP | 既存の二つ以上の consumer が共通結果を利用。未知の mode のための拡張点は作らない |
| LSP | 同値改修では返却型・順序・例外を維持。tuple を暗黙に lazy object にしない |
| ISP | block consumer に不要な表示リンク・全経路の生成を要求しない |
| DIP | 解析関数は検証済みデータと設定を受け取り、Vue / Worker / cache に依存しない |
| KISS | 一回走査、dict の逆引き、直近一解析の容量管理から始める |
| DRY | identity、スコア、同点順位、所属判定、形式変換を一つの owner が定義 |
| YAGNI | GPU、別言語移植、汎用 graph/cache framework、未要求の経路閲覧 UI を先行導入しない |

## 5. 計算量の改修設計

記号: `M` は manifest の転送・検証対象量、`J` はキー数、`H` は HSP 行数、`h_p` はペア p の HSP 数、`U` は未所属数、`A` は core member 総数、`E` は保持 evidence 数、`V` はグラフノード数、`P` は経路数、`L` は全経路を展開した総要素数、`N` は block anchor 数。

### 5.1 バッチキー生成（既存改修を検証）

全 manifest の扱いを `J` 回からバッチにつき一回へ減らし、キー生成部分を概ね `O(M+J)` にする。これは key helper に限定した評価であり、他の validation が存在しないという主張ではない。順序、direction、searchContext、失敗、キャンセルを検証する。

### 5.2 HSP 集計（S03）

テーブルの行を一回走査してペア別 accumulator へ振り分ける。代表行、coverage 区間、HSP 数、alignment length 合計だけを保持し、不要な DataFrame / tuple の全量コピーを避ける。

目標は概ね `O(H + Σ h_p log h_p)`。coverage union の整列は必要であり、全体を線形と呼ばない。member limit は集計より前の現在の意味を維持し、選ばれたペアの HSP を落とさない。

### 5.3 疎な所属候補と metadata（S04）

incoming/outgoing evidence と protein→group の逆引きから候補を得る。候補グループだけを絞ってから全 member を再走査する構成で終えず、実際に接続する evidence から support / diagnostic を還元する。

固定 core snapshot 用と拡張後集合用の索引は段階境界で構築する。既存 scoring と tie-break の owner を使う。candidate sort を含む仕事を実在 evidence と接続先グループに基づく量へ近づけるが、密な入力や selector 全体の線形性は主張しない。

protein→RBH group、protein→member、endpoint→edge metadata の逆引きも必要な consumer で一度作る。多対多、重複、最初に選ばれる edge と出力順を保持する。最新 selector の `comparison_pairs=()` が要件を満たす場合、不要な表示生成の抑止に利用する。

### 5.4 解析段階とキャッシュ（S05）

| 段階 | identity に必要な主な情報 | 再利用可能な例 |
|---|---|---|
| raw | 実入力・protein/record bindings・方向・検索引数・検索DB範囲 | member limit、色、block条件の変更 |
| filtered | raw identity と filter thresholds | member limit、block条件の変更 |
| normalized | filtered identity、member limit、protein長など消費する情報 | 推論に無関係な block条件の変更 |
| inferred | normalized evidence、所属・順位・命名が読む情報、推論設定 | 色や独立した block条件の変更 |
| block/display | anchor、unit、gap、drift、conflict、display pairs、座標変換 | 下流 appearance の変更 |

この表は実装前に consumer の実際の read-set で具体化する。各行のために独立した cache を増設する指示ではない。raw と persisted derived identity の既存 owner・形式を維持する。

`_protein_sort_key()` は座標を読み、group/path ID の並びに影響する。向き・record順・注釈を無条件に inferred key から除かない。向き変更では raw/normalized evidence を再利用し、必要な下流を再計算する。独立性を証明した段階だけ reuse を広げる。

第一案は、直近一解析の必要な中間データ一式を単一 owner がバイト予算内で保持する方式。予算は native/browser の実測と既存 lifecycle に基づいて決め、上限値と算定方法を S05 handoff に記載する。件数上限だけには依存しない。

- 一式が予算内なら保持し、64表前後でも sequential eviction を起こさない。
- 超過時は保持しない。通常の計算経路で正確な結果を作り、理由を診断情報に残す。
- 保持予算と実行中の peak memory を区別する。新旧 snapshot、転送コピー、DataFrame、Wasm heap も測定する。
- 段階結果は不変に扱い、変更する境界でのみコピーする。
- 新しい cache hit は必要な入力・identity validation を迂回しない。
- キャッシュ公開の確定点と failure/cancel/stale 時の破棄を既存契約に合わせる。完了 raw retry の保存は既存 owner に任せる。
- Worker 終了、Clear Cache、入力 / Session / History 置換による失効を明示する。
- Python converted JSON は既存 JS derived cache と consumer を照合して整理し、同じ完成 payload の保持 owner を追加しない。

### 5.5 経路表現（S02・S06）

現在の `_build_ortholog_paths()` は全経路を列挙し、protein path で重複除去し、edge列の優先順位を決め、整列後の連番を path ID とする。各 edge には最初に含まれる path ID が入り、sharedProteinIds は複数経路への登場に依存する。

推奨は、同じノード・辺・選択規則を表す lossless graph を標準の中間・保存表現にすることである。DAG の構造保持は `O(V+E)`。非循環性と重複除去規則を証明した後、動的計画法で正確な経路数等を計算する。sort、大整数、旧 path ID の順位再現まで一律に線形とはしない。

全経路の明示的な列挙を consumer が要求する場合は、少なくとも `Ω(L)` の時間・出力量が必要。generator 化だけで、その後の JSON 配列化の指数コストは解消しない。

S02 の決定対象:

| 安定 choice code / outcome ID | 内容 | 性能上の限界 |
|---|---|---|
| `PATH-A` / `exhaustive-current` | 現在の全経路配列の公開・保存契約を維持 | 全量要求の指数コストは残る |
| `PATH-B` / `lossless-graph` | 関係グラフを通常の公開・保存経路に採用。必要な既存 consumer の経路取得契約を明示 | 通常経路は全列挙を回避。明示的な全量取得は出力量に比例 |

`PATH-B` はscenario revision 1の完全なhuman receiptで選択済み。S02 は API、exact count、path ID、shared情報、互換読込、保存形式、失敗時継続を具体化した。authorityのreview・base統合は未完了。DAG でない到達可能入力を勝手に切り落としたり、edge 方向を変えたりしない。循環への対応が未決なら影響部分を判断へ戻す。

`OrthogroupResult`、`OrthologPath`、`OrthologEdge`、typed resource、SVG属性、catalog、popup、persisted raw/derived/session を consumer ごとに追跡する。catalog が最後に配列を count へ縮約しても、前段の encode が全量展開していれば未解決である。

JavaScript の安全な整数範囲を超える count も扱う。数値表現を暗黙に丸めず、必要な契約を決定する。未要求の paging UI や汎用 graph export は追加しない。

### 5.6 Cluster merge（S07、条件付き）

query order を一度整列し、二分探索で conflict 判定区間を狭める。merge 可否の内部判定だけなら max_conflicts 超過で終了する。exact count の consumer に打ち切り値を返さない。

順序済み endpoint を再利用し、成長 cluster の全コピーを確定時まで遅らせる案を比較する。sort と copy の最適化が merge 順・block ID・singleton 保存を変えないことを検証する。`O(log N + k)` の区間候補取得は、cluster 全体の線形性を意味しない。高度な多次元索引はこの段階の実測で必要と確認された場合だけ検討する。

## 6. セッションと依存関係

各ファイルをそのセッションの instruction prompt として使用する。共通規則は本書の §7〜§9 を参照する。全セッションを一度に実行する指示ではない。

| Session | Prompt | 必須前提 | 成果 / 次への条件 |
|---|---|---|---|
| S01 | [基準・再現・契約確認](SESSION_01_BASELINE.md) | なし | baseline、入力hash、観測とauthorityの対応、測定スクリプト |
| S02 | [経路契約と設計判断](SESSION_02_PATH_CONTRACT.md) | S01 | consumer inventory、比較設計、必要なら Product Decision Pack |
| S03 | [HSP集計](SESSION_03_HSP_AGGREGATION.md) | S01 | [完了報告](results/S03.md)：同値な一回走査集計、時間・メモリの観測比較。独立再測定は最終指示で省略、測定条件未達を明記 |
| S04 | [疎な候補探索とmetadata](SESSION_04_SPARSE_SUPPORT_AND_METADATA.md) | S03 | [完了報告](results/S04.md)：evidence・metadata索引、全結果一致。時間22 pass / 2 inconclusive、全項目合格は未達 |
| S05 | [共通解析境界と容量付きcache](SESSION_05_SHARED_ANALYSIS_AND_CACHE.md) | S04 | 段階依存表、bounded reuse、Web lifecycle 検証 |
| S06 | [経路表現の実装](SESSION_06_PATH_REPRESENTATION.md) | S05、S02の選択を認可するbase authority | 選択済み契約の runtime / reader / writer、性能検証 |
| S07 | [Cluster merge評価・必要な改修](SESSION_07_CLUSTER_MERGE.md) | S05。S06完了時はその結果も含む | 改修または実測に基づく見送り判断 |
| S08 | [統合・性能・回帰検証](SESSION_08_INTEGRATION.md) | S03〜S05、S07。全体完了にはS06の扱い確定 | 正確な達成範囲、最終 gate、残課題 |

推奨順序は S01 → S02 → S03 → S04 → S05 → S06 → S07 → S08。S02 が判断待ちでも S03〜S05 を進められる。S06 が判断待ちなら S07 と限定的な統合検証は進められるが、指数増加対策を完了扱いにしない。これらは依存関係の説明であり、sub-agent の自動起動を要求しない。

optional inference と limit retention は統合済みの回帰対象である。開始時のbaseにその変更が欠ける場合はcheckout/依存revisionを特定する。OFF の期待動作を仮実装した test stub で最終合格にせず、統合された実経路で確認する。

## 7. 全セッション共通の実行規則

1. リポジトリの `AGENTS.md`、`CLAUDE.md`、Webに触る場合は `gbdraw/web/CLAUDE.md` を読む。instruction prompt の実行だけから `$execute-plan-with-evidence` の利用を推定しない。
2. 編集前に working tree を確認し、対象・対象外の既存変更を区別する。必要な新規 work branch は fetch 後の最新 `origin/dev` から `git switch --no-track -c <branch-name> origin/dev` で作る。汚れた共有ツリーで安全に切り替えられない場合は isolated worktree を利用し、ユーザーの作業を stash/reset で隠さない。
3. 同じセッションの再開は既存 work branch と handoff を確認して継続する。前セッションの必要な runtime 変更が base に無い場合、混ぜて実装を複製せず、依存差分を特定する。承認済みの継続方法があればそれに従う。
4. 本書は push、merge、公開、deploy、tag の許可ではない。コミットを行う場合は branch/upstream を確認し、`main` / `dev` に直接作成しない。
5. 指定セッションの担当範囲を実装・検証・handoff まで完了する。単にコード案を示して実行済みとしない。対象外の大規模整理や並列 owner を追加しない。
6. [Architecture Fitness Function Ratchet](../ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md) に従い、owner/path の前後、削除した旧経路、振る舞い検証、rollback を記録する。OE/PE/CB の完全な集合と例外判断は同 policy の例外条件が成立する場合だけ作る。
7. [Product Impact Ratchet](../PRODUCT_IMPACT_RATCHET.md) の該当範囲・developer preflight を確認する。`IMPLEMENT_EXISTING_AUTHORITY` / `EVIDENCE_REQUIRED` / `PRODUCT_DECISION_REQUIRED` / `NOT_ALLOWED` を理由付きで分類する。選択未決の依存 runtime だけを止め、独立した作業は継続する。
8. 新しい Product 選択は [Decision Pack template](../PRODUCT_DECISION_PACKET_TEMPLATE.md) を用いる。機械表現は人の明示的選択だけから作り、同じ候補の authority で同じ候補の runtime を自己認可しない。存在しない `BD-###` を作らない。
9. 公開の reader/migrator 追加は `main` first-parent または release tag の契約と代表fixtureで裏付ける。session/request/resource/cache/catalog の namespace を別々に確認する。branch-only 中間形式の migration chain を残さない。
10. production / tests / documentation / generated artifacts の差分を分けて監査する。通常検証で reference_outputs を更新しない。生成wheelは必要時に作成しコミットしない。owner-maintained social preview は変更しない。
11. benchmark の旧実装は必要ならテスト・比較用に隔離し、本番 fallback として残さない。cache miss / over-budget は明示した同じ計算経路で処理する。

## 8. 検証設計

### 8.1 再現資産

S01 は既存の `tools/benchmark_diagram_layout.py` の測定・source-root 比較方式を参考に、今回の解析用に必要な最小の runner を用意する。既存の適切な runner があれば拡張する。候補名は `tools/benchmark_protein_comparison.py`。実際に採用した一つのコマンドと引数を S01 handoff に固定し、後続は同じ runner を使う。

- 入力: Gallery の GenBank / 保存raw、合成fixtureのseed・generator・設定。
- 同一性: ファイルhash、protein/runtime identity、record/source順、scope、inference、各limit、filter、thread/search設定。
- 出力: stageごとの一致、group/member/edge/path/block metadata、typed resource、SVG geometryと必要なsemantic属性。
- 比較基準: 同じ設定を与えた変更前後。fresh default の変化と algorithm の効果を混ぜない。
- 保存: 実行コマンド、exit code、source commit/diff、依存バージョン、各sampleを機械可読結果へ記録。

再現コードとfixture生成方法はリポジトリで保持する。大きな一時生成物は作業ディレクトリに置けるが、handoffを一時パスだけに依存させない。report が未追跡なら、引用元としての扱いを明記し、最新コードで再現を構築する。

### 8.2 ケースと合格条件

| Case | 必要な検証 |
|---|---|
| manifest | バッチ1回検証、各キー・順序・searchContext一致、改ざん拒否、warm時転送量 |
| HSP | overlap、disjoint、reverse、clamp、missing/unknown ID、NaN/非有限、同点、重複、空入力 |
| sparse support | groups/unassigned 200/200、400/400、800/800。無関係group追加で候補数が総当たり増加しない |
| support correctness | incoming-only、same-record、cross-record、domain-only、best/second、tie、snapshot、record-local競合、dense入力 |
| cache | 49/64/81表、予算内外、cold/warm、色・block・member・filter・向き・順序・入力変更 |
| paths | 小規模で全path/ID/shared一致。R=8/12/16で既存式を照合。R=24以上はcompact経路で全列挙を避ける |
| path count | JS安全整数範囲を超えるcaseを、全経路展開せず正確に処理 |
| merge | 300/600/1200 anchors、strict boundary、reverse、singleton、max_conflicts、chain merge |
| cross-surface | CLI/Python/Webがそれぞれ公開する範囲で、Similarity、Collinear adjacent/all、inference ON/OFF、有限/無制限limit、multi-record sources |
| lifecycle | raw完了→downstream cancel→member変更→retry、raw設定変更、Clear Cache、Session/History置換、stale完了、Worker再作成 |

各caseは実装上必要な規模で実行する。旧実装の巨大な全列挙で故意にメモリを枯渇させず、大規模旧経路の理論値と実測を区別する。

時間測定は warmup と複数sampleの中央値・ばらつきを記録する。S01 で測定回数と許容変動を確定し、変更後の結果を見て合格基準を緩めない。profiling、tracemalloc、RSS、browser/Wasm memory は別runで測る。native の改善率をそのまま Web 全体の改善率としない。

安定した回帰テストは検証回数、走査行数、候補数、parse数、保持量、全列挙の有無を中心にする。wall-clock gate は同条件の反復測定とnoiseに基づく。小さい入力の負担増、cache oversize時の性能、dense入力も報告する。

### 8.3 必要な既存 gate

変更に対応するfocused testsを先に実行し、最後に共通解析・Web変更の広さに応じて以下を実行する。

```bash
pytest tests/test_protein_colinearity.py tests/test_collinearity.py tests/test_collinearity_units.py -v
pytest tests/test_web_feature_catalog.py tests/test_session_request_codec.py tests/test_session_compat.py -v
pytest tests/ -v -m "not slow"
pytest tests/test_output_comparison.py::TestOutputComparison -v
ruff check gbdraw/
```

Node testsは最新treeの実在inventoryに合わせる。主要対象は `losat-cache.test.mjs`、`run-analysis-derived-cache.test.mjs`、`run-analysis-simple-path.test.mjs`、raw/session identity、Worker lifecycle、architecture/Product contracts。report の Node import 失敗を現在も起きると断定せず、失敗時は所有する境界を診断する。

Browser検証では以下を確認する。

```bash
command -v playwright && playwright --version
python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
node -e "console.log(require.resolve('@playwright/test'))"
```

Node版が無ければ Python Playwright で該当チェックを行う。Chromium sandbox制約なら同じlocal checkを必要な権限で再実行する。テストcommandは少なくとも30分を許容してincrementalに監視し、短いtest-owned timeoutを性能問題の回避として変更しない。

実browserでlocal assetsのみの動作、実helper/typed render経路、warm reuse、cancel、保存と再生成を確認する。共有render/cacheに触れた場合はCircularもsmoke確認する。public figureを更新する必要が生じた場合に限り対応skillと生成ownerを使用する。

## 9. Handoff と完了判定

各sessionは `results/S01.md` など一つの短いhandoffを作る。未実行の結果ファイルを今から成功状態で用意しない。handoffには次を含める。

1. base/head、branch/upstream、依存sessionのrevision。
2. 完了事項、未完了事項、該当authorityとpreflight分類。
3. owner/pathの変更、旧経路の削除、scope外の既存変更。
4. 再現コマンド、fixture/hash、tests/benchmarksの実測結果とlimits。
5. cache/保存形式/互換性への影響とrollback方法。
6. 次sessionが使う実在API・資産・前提。既承認outcomeを再質問しないための根拠。
7. 英語の proposed commit title と短い英語summary。実commit/pushの有無を明記。

S08は「同値な性能改善」「経路契約」「Web lifecycle」「architecture」「未達」を分けて総括する。S02の選択待ち、依存revisionの不足、未実施browser checkを成功扱いにしない。PATH-Aを選択した場合は、指数出力量の残存を最終報告に明記する。
