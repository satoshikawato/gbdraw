# S07.5 — Gallery主要処理の追加改修計画

作成日: 2026-09-16。**計画完了。本番実装・S08は未実施。**

最初に **S07.6: Collinear unit索引の重複生成を削減**し、次に
**S07.7: fit用の一時DataFrameを必要な4列へ絞る**。
unit索引はCollinear推論OFFの最大の残存ownerで、ONにも適用できる。
fit行選択はCollinear ONとSimilarity groupsの共通経路にあり、局所prototypeで
選択順・fit係数の完全一致と一時領域の減少を確認した。
数値変換境界の統合、member/RBHのアルゴリズム置換、group全体のrank索引は
今回の実装順に入れない。大きな親関数の割合だけを根拠に広い改修を始めない。

S07の性能評価はユーザーの「合格とします」「合格でいいです」により**合格済み**。
14 pass / 6 inconclusive / 0 confirmed regression、+10.92%のnoisy microcase、
外部作業・歴史的baselineの制約、Gallery isolated mergeのpeak増加をそのまま保持する。
再測定はS07.5や後続改修の開始条件ではない。
S06の通常経路から全path列挙を除いた主目的はユーザー了承済み。
S06の互換性Reviewは将来のmerge前に必要だが、ローカル計画・実装を止める条件ではない。
S05は2026-09-15却下済み。共通解析class、cache owner、Worker transaction、
scheduler、runtime変更を依存に含めない。

## 1. 基準revisionとauthority

| 項目 | 確認した状態 |
|---|---|
| fetch後の`origin/dev` | `5e0cb0fa1d9920592da431b764fa2c9133b5bfdf` |
| S07 branch / HEAD | `perf/collinear-s07-cluster-merge-20260916` / `f52ed68e15e527b08e4d9247b990684040fe4270` |
| S07開始時状態 | clean、upstreamなし |
| S07直前baseline | `ec16110cf08b894b6dce16aa973e1d6112554ac9` |
| 元のS06固定commit | `9e6c6e740e5886b7cc627386c303973325a785a9`。上記baselineとtree同値 |
| 今回のworktree | `/mnt/c/Users/genom/GitHub/gbdraw/.worktrees/collinear-s07-5-20260916` |
| 今回のbranch | `perf/collinear-s07-5-gallery-plan-20260916`、upstreamなし |
| 調査する本番コード | S07 HEADのまま。共有dirty `dev`のコードを混ぜない |

fetch後、専用worktreeで`git switch --no-track -c <branch> origin/dev`を実行した。
`git cherry`で未統合・patch-uniqueと確認した直系依存
`3a3b29c7 → e3f93172 → 603e0787 → f35f2e63 → dc79915a → ec16110c → f52ed68e`
を、fast-forwardで一度だけ引き継いだ。新しいmerge commitは作成していない。
Git管理領域のread-only制約は同じローカル操作をsandbox escalationで再実行して解消した。
共有ツリーの既存変更に書き込まず、stash/reset/cleanを使用していない。

参照したauthorityは、最新baseの
[Product ratchet](../../PRODUCT_IMPACT_RATCHET.md)、
[architecture ratchet](../../ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)、
[OIPC revision 6](../../OPTION_INTEGRITY_PRODUCT_CONTRACT.md)
（PD-OI-001/002/004/005/006/007/014/016/018/021/022/023）、
[comparison semantics](../../../REFERENCE/comparison-programs-thresholds-and-results.md)、
[Session互換性](../../../SESSION_COMPATIBILITY.md)、Web CLAUDE、
Product Impact map/decisions、[S02契約](S02_PATH_CONTRACT.md)と
[receipt](S02_PRODUCT_DECISION.md)、[S06実装の具体化](S06.md)。
S02文書の「authority未統合」「S05完了を前提」は歴史的記録であり、
PR #536でのPATH-B統合・S05却下を反映したMASTER/S06/S07と本依頼が現在の前提である。
存在しない`BD-###`は引用しない。

今回の計画と二つの局所改修のpreflightは **IMPLEMENT_EXISTING_AUTHORITY**。
既存の科学的結果、型、順序、失敗、保存、次の操作を保つ実装だけを対象とする。
選択すべき新outcomeがないためPRODUCT_DECISION_REQUIREDではなく、禁止された案を
採用しないためNOT_ALLOWEDでもない。後述の不足は性能・同値性の検証課題で、
現時点のProduct選択を未決にするEVIDENCE_REQUIREDではない。
同値にできずmaterialな効果変更が必要になった部分だけ、証拠を整理して既存の
Product判断手順に戻す。PATH-B・推論ON/OFF・limitsは再質問しない。

## 2. 再利用した証拠と追加調査の範囲

[s07-5-evidence-reuse.json](data/s07-5-evidence-reuse.json)に、report hash、
source照合、関数単位のrunner比較、既存sample・MAD・profile・counterをまとめた。
これは測定の再実行ではない。S07のreview台帳の全entryも変更前にhash照合した。

| 再利用対象 | 用途と適用範囲 |
|---|---|
| `s07-start.json`, `s07-acceptance.json`, `s07-review.json` | revision、明示的合格、source/evidence対応 |
| `s07-baseline-reuse.json`, `s07-comparison.json`, `s07-summarize.py` | 既存判定・測定境界・集計方式を維持。S07集計を上書きしない |
| `s07-current-probe.json.gz`, `s07-vibrio-current-probe.json.gz` | Collinear ON/OFFの親子関係・残存処理の診断 |
| `s07-current-gallery-timing.json.gz`, `s07-vibrio-current-timing.json.gz` | 変更前post-search baseline。Hepは21、Vibrioは7 samples |
| `s07-current-memory.json.gz`, `s07-vibrio-current-memory.json.gz` | post-searchのPython peak/retained、出力hash |
| `s06-current-probe.json.gz`, `s06-current-timing-final.json.gz`, `s06-current-memory-final.json.gz` | **Similarity自身**のprofile、21-sample時間、memory。S07はこの解析経路を変更していない |
| [S03](S03.md)、[S04](S04.md)、[S06](S06.md)、[S07](S07.md)と小規模oracle | 既に削除した仕事、同値性、既存の留保・回帰保護 |

Similarityのprotein inference、lossless graph、metadata、typed serializationの
source hashはS06 reportと一致する。S07のCollinearity merge変更と、S06 probe後の
Session promotion修正はSimilarityのnative `post_search`から呼ばれない。
runnerはファイル全体が変わったが、`gallery_input`、`measure_stage`、canonical/hash、
既存caseのparse/filter/analyze関数はAST一致。追加counterはmerge用である。
従ってSimilarityの既存結果を無効にして再測定する理由はない。

既存証拠に不足していたのは次の二点に限定した。

1. unit builderの重複オブジェクト・member sort・alias一時setの実数量。
   各Galleryの抽出済みproteinからcall-local counterで確認した。
2. fit行選択で24列すべてを運ぶ必要があるか、4列投影の追加copyを含めても効果が
   ありそうか。既存入力から作った同一のfit直前frameに対し、計画専用prototypeを比較した。

追加コードは[data/s07-5-diagnose.py](data/s07-5-diagnose.py)だけで、
本番からimportしない。inventory、timing、tracemallocは別process。
S03/S04/S06/S07の全体benchmark、pytest全体、wheel、browserは再実行していない。
新しいraw検索、描画、保存も実行していない。

診断コード内の独立した小規模list oracleは48 cases（size 0/1/3/4/7/8/15/16/31/32/33/65、
fraction 0.1/0.25/0.5/1.0）で、full tie・bin境界・0 score・入力不変を確認した。
Gallery三入力では全fit点の順序とモデル係数が完全一致。これはprototypeの局所証拠で、
本番変更後の科学的出力・Session・geometry検証を済ませたという意味ではない。

## 3. 測定境界と現在のコスト

native `post_search`は準備済みthreshold-qualified tablesとprotein extractionを入力とする。
Collinearはunit構築、scope選別、member選択、ONなら正規化・推論、anchor/block構築まで。
Similarityはmember選択、正規化・推論、selector内部のdisplay-edge選択まで。
raw検索、parse/filter、protein抽出、下流のgenomic display変換、metadata/typed encode、
最終SVG描画、ファイル保存、browser/Worker起動を含まない。独立stageを加算しない。

### 3.1 既存native全体時間とmemory

| ケース | median ms | MAD/median | samples | peak / retained bytes |
|---|---:|---:|---:|---:|
| Hep Collinear ON | 507.203 | 1.50% | 21 | 19,017,096 / 10,252,071 |
| Hep Collinear OFF | 114.958 | 1.76% | 21 | 6,582,680 / 4,853,226 |
| Vibrio Collinear ON | 3,044.608 | 0.91% | 7 | 115,083,765 / 31,999,765 |
| Vibrio Collinear OFF | 458.673 | 0.44% | 7 | 26,246,539 / 9,677,006 |
| Hep Similarity | 879.741 | 5.052% | 21 | 29,024,530 / 9,244,243 |

Similarityの時間はS06の**inconclusiveな観測**のまま再利用する。精度のよいbaselineへ
黙って読み替えない。peakは別runのtracemalloc一回値で、RSS・Wasm容量ではない。
これらはWeb全体待ち時間でもS07.5の改修結果でもない。

### 3.2 計測器入りprofileの内訳

下表は各caseの`analyze`累積時間に対する割合。太字の親の子を別行に示しているため、
**縦に合計しない**。instrumentation自体の負担も含む。

| owner / 関数 | Hep ON | Vibrio ON | Hep Similarity | Hep OFF | Vibrio OFF |
|---|---:|---:|---:|---:|---:|
| **`_normalize_directional_hit_tables`** | **48.55%** | **42.27%** | **49.92%** | 不通過 | 不通過 |
| その子: S03 HSP aggregation | 20.93% | 26.99% | 20.74% | 不通過 | 不通過 |
| その子: `_select_normalized_fit_rows` | 17.93% | 7.45% | 19.30% | 不通過 | 不通過 |
| **`_build_anchor_core_orthogroups`** | **18.88%** | **27.12%** | **19.30%** | 不通過 | 不通過 |
| `_anchor_core_hit_rank`（複数phaseにまたがる） | 4.89% | 6.79% | 6.29% | 不通過 | 不通過 |
| `_select_member_candidate_hits_per_query` | 9.79% | 6.57% | 10.88% | 28.98% | 17.82% |
| `_select_orthogroup_edges`（OFFのRBH等） | 不通過※ | 不通過※ | 不通過 | 24.90% | 15.15% |
| `build_collinearity_unit_index` | 5.13% | 8.80% | 不通過 | 24.63% | 54.18% |
| `_coerce_outfmt6_numeric_columns`（複数親の子） | 8.96% | 6.84% | 9.84% | 23.07% | 16.50% |

※このON caseは`rbh`。ONの`all/one_to_one`ではこの関数も通る。
OFFのmember選択とedge選択は兄弟phaseだが、numeric coercionはその中に含まれる。
`_fit_expected_bitscore_model`の比率18.01% / 7.51% / 19.37%のほぼ全部が
fit行選択であり、係数計算そのものを高速化する根拠は薄い。
rankは25,828 / 200,406 / 53,696 callsだが、group builder全体の割合を
rank削減の利益と見なすことはできない。

### 3.3 入力と実行設定を混同しない

| 入力 | proteins | 保存raw tables / rows | filter後rows（scope適用前） | 推論で正規化するtables / member後HSP rows |
|---|---:|---:|---:|---:|
| Hep Collinear | 2,829 | 13 / 98,947 | 13,720 | 13 / 11,201 |
| Vibrio Collinear | 24,027 | 59 / 657,743 | 170,223 | 31 / 80,888 |
| Hep Similarity | 2,829 | 25 / 183,661 | 21,935 | 25 / 18,699 |

Collinear OFFの選別対象はselfを除く8 / 20 tables。
Vibrioの59保存tables全部を31正規化tablesと混同しない。
sourceCountはHep 5、Vibrio 11。raw job数とdirectional table数も同じではなく、
保存cacheから歴史的raw検索時間・実invocation数は復元していない。

既存runnerが**明示的に実行する**設定はmember=5、bitscore=50、evalue=0.01、
identity=0、alignment_length=0、Collinear adjacent/auto/rbh、related=2、
`LosslessCollinearityParameters()`（min_anchors=1、gap=0、drift=0、conflicts=1、either）。
`-off`以外は推論ON。Similarityは全方向evidenceを使い、displayは既存の隣接projection。
これは保存Session設定やfresh defaultsを一律に実行した測定ではない。

Hep Collinear Sessionは保存ResultがCollinearでもactive draftのmodeは`orthogroup`。
Vibrio保存draftのblock設定はmin_anchors=3、gap=2、drift=2で、上記runnerと異なる。
Similarity fresh member limitは無制限だが、この既存caseは明示的5。
今後は「既存benchmark比較」と「保存recipeを実modeで再生成」を別case名・hash・境界で記録する。
Vibrio Similarityの独立profileは存在しないため、その速度寄与は推定しない。

## 4. 最初に実装するS07.6 — unit索引

**owner:** `gbdraw/analysis/collinearity_units.py::build_collinearity_unit_index`。
callerは`collinearity.py::build_orthogroup_collinearity_blocks_from_hits`でON/OFF共通。
native検索builderと`build_collinearity_blocks_from_hits`、Web
`python-helpers.js::convert_losatp_blastp_pairs_to_genomic_payload`がこの既存経路へ収束する。
Similarity selectorはこのownerを呼ばないため、Similarity高速化は主張しない。

### 4.1 実装する三つの変更

1. aliasesを既存`_unit_aliases`で先に作り、`CollinearityUnit`を一回だけ生成する。
   aliases空の中間instanceと、その`__dict__`を展開した再構築を削除する。
2. membersは、同一keyでsort済みの`sorted_proteins`から順にappendされた部分列。
   CDS singletonもこの順を満たす。`cds_members`作成時の同じkeyによる再sortを削除し、
   元のmembers順をtupleにする。最初のprotein sortとgroup sortは維持する。
3. recordごとの`alias_targets: alias -> set[unit_id]`を、返却するunique-alias dictと
   ambiguous-alias setの一回走査に置き換える。初出unitを保持し、別unitへの衝突で
   dictから除外してambiguousへ移す。一度ambiguousになったaliasは三回目以降に復活させない。
   同じunit内の重複はambiguousにしない。生き残るdict keyの初出順を維持する。

別の索引class、lazy unit、永続cache、モード別fast pathは作らない。
3は返却される既存二containerだけで完結する局所実装であり、別のalias解決ownerではない。
`strong_locus_id`、代表protein順位、座標/strand、display_name、警告と例外は変更しない。

### 4.2 処理量とmemory

`P`: proteins、`U`: units、`m_u`: unit内members、`A`: unit別の重複除去済みalias総数。

| 構造 | Hep | Vibrio | 改修後 |
|---|---:|---:|---|
| unit object生成 | 5,658 | 48,054 | 2,829 / 24,027（2U→U） |
| membersの再sort calls | 2,829 | 24,027 | 0（元の順序は維持） |
| 再sortへのmember参照 | 2,829 | 24,027 | 0 |
| `alias_targets.setdefault(alias, set())` | 15,769 | 151,181 | 0（A個のset生成と一時dictを削除） |
| record内distinct aliases | 15,761 | 150,911 | 不変 |
| ambiguous aliases | 8 | 229 | 不変 |

両Galleryはこの抽出でcollapsed unitが0。従って再sort削除は実Galleryでは主に
singleton sortの固定費削減であり、大きいlocusの実測効果とは言わない。
instance生成とalias追加試行はwrapperで実countし、member sort・alias setdefault数は
その出力のU/Aと現行sourceの一回ずつのloopを照合した値である。
Similarity inventoryにも同じ抽出入力のunit診断を載せたが、Similarity本番がunitを
構築するという測定ではない。
一般入力では`Σ m_u log m_u`の再sortを除くが、record protein sort、unit group sort、
代表選別、座標集計、alias正規化は残る。全体をO(P)とは呼ばない。
`_unit_sort_key`の四つの独立min（特にmin end）は変えない。

返却indexの内容・保持量は不変。一時instance、dict展開、member list、alias setの
割当量は減る。最終outputが支配すればpost-search peakの差は小さい。
中間instanceはunitごと、alias_targetsはrecordごとに作られるので、全割当数を同時保持量と
解釈しない。alias部分の余分なpeakは主に最大record内alias数に依存し、全recordのAとは異なる。
今回unit改修の時間・peak比較は未測定なので、何%速くなるという下限保証はしない。
OFFのprofile 24.63% / 54.18%はownerの機会を示すだけで、削減可能割合ではない。
ONの寄与は同ownerの5.13% / 8.80%の範囲内の一部、Similarityは0。
この大きなOFF owner、確実に不要な生成、狭い変更境界から最優先とした。

## 5. 次に実装するS07.7 — fit行選択の列投影

**owner:** `protein_colinearity.py::_select_normalized_fit_rows`。
共通callerは`_normalize_directional_hit_table → _fit_expected_bitscore_model`。
その上流は`select_rbh_orthogroup_edges_from_directional_hits`で、
Collinear ONとSimilarityが利用する。OFFは通らない。

現行は各tableのaggregation/coverage判定後の24列を
`length_product, query, subject`順にsortし、最大8 binsを作り、各binを
`bitscore desc, query, subject`順にsortしてtop fractionを選ぶ。
Hep ON 13 tables / 104 bins / 532 fit点、Similarity 25 / 200 / 937、
Vibrio ON 31 / 232 / 4,010。bin数が8未満の小表もある。

採用する変更は、関数冒頭の非empty境界で
`length_product, query, subject, bitscore`の**4列だけ**を一時投影し、
その後は同じpandas stable sort、bin境界、ceil/head、scalar `_row_float`と
`math.log10`を使うもの。一時frame以外の正規化出力を4列に削ってはいけない。
score計算式、fit係数計算、sum/zip順、fallback、coverage判定位置は変えない。
numpy回帰・vectorized power・heap/nlargest・全table共通fitは採用しない。

計算量のsort項は`O(N log N + Σ b log b)`のまま。
削減するのはsort/slice/itertuplesが運ぶ列とnamedtuple生成の固定費。
投影自体に4N列要素分の新しい仕事があるため、極小表の回帰は別検証する。
元の入力frameは呼出中保持され、深いbyte推定は共有stringの物理copy量ではない。

### 5.1 計画専用prototypeの観測

[inventory](data/s07-5-inventory.json.gz)、[timing](data/s07-5-fit-timing.json.gz)、
[memory](data/s07-5-fit-memory.json.gz)、[timing command/host](data/s07-5-fit-timing-command.json)。
対象は**準備済み全fit直前frames→順序付きfit点list**だけで、候補の4列投影を含む。
full normalizationやpost-searchの比較ではない。

| 入力 | original→4列 median ms | MAD/median original→4列 | 変化 | Python peak bytes original→4列 |
|---|---:|---:|---:|---:|
| Hep ON | 96.779→66.209 | 0.88%→1.87% | −31.59% | 859,823→466,763 |
| Hep Similarity | 185.371→127.135 | 0.43%→0.85% | −31.42% | 1,021,742→604,291 |
| Vibrio ON | 266.284→194.016 | 0.97%→0.58% | −27.14% | 2,237,692→1,126,115 |

warmup 1、双方7 samples、CPU 3、PYTHONHASHSEED=0。各sampleと完全な出力hashを保存。
同じcase内でoriginal→候補を連続測定し、準備とhash計算はtiming外。
全体の21-sample baselineとの比較ではないため、ここは7 samples。
数値判定は3 pass、中央値10%超悪化・MAD/median 5%超の基準は維持した。
ただしホスト前後snapshotで別のLOSAT処理とGit作業を観測した。
**独立した性能合格としては採用しない**。外部processに介入せず、都合のよい再測定もしない。
このセッション自身のbuild/checkout/test/別benchmarkをtimingへ重ねていない。

tracemallocは別process、一回ずつ。retainedは427,915→292,233、652,738→386,418、
1,027,748→802,551 bytesで、同一の返却fit点を保持した時点の残留割当を含む。
これは永続cache削減量ではない。各caseの元frame群と4列frame群のdeep bytesは
3,657,316→1,637,116、6,306,740→2,823,260、29,614,452→13,250,832。
これらを同時に全て保持する本番変更を指示していない。実装は現在どおりtableごとに処理する。

期待はON/Similarityのfit行選択を数十ms単位で軽くできる可能性。
上の局所差を別runの全体medianで割った「全体改善率」は作らない。
Pyodideのpandas、JS transfer、typed encoding、renderが支配すればWeb待ち時間への寄与は
小さくなる。OFFの改善は0。実装受入には最終sourceでの検証が別途必要である。

## 6. 見送る対象と再着手条件

| 対象 / owner | 残る具体的仕事とmemory | 今回見送る理由 / 再着手条件 |
|---|---|---|
| numeric coercion / `_coerce_outfmt6_numeric_columns` | callごとのfull-frame copy、10 numeric columnsの`to_numeric`、numeric matrix、Python finite scan。ONはHep 26/Vibrio 62、OFFは16/40、Similarityはdisplay分込み54 calls | validationはpublic parse/filter/best-hit等でも必要。有限member選択はcoerceしたrankでpairを選び**元HSP rows**を返すため、下流coercionを単に削れない。呼出間共有・検証済みclassはS05復活になる。4列fit/unit改修後も残存支配が確認され、整数/nullable/object/overflow/NaNと同一errorを保つ局所finite検証prototypeに時間・memory根拠が得られた場合のみ再提案 |
| member候補 / `_select_member_candidate_hits_per_query`, `select_top_hits_per_query` | 有限はcoerce、6-key stable sort、distinct pair→query head、二つのMultiIndexと元行mask。概ねO(H log H)、一時frame/index O(H)。無制限はvalidate+copy/resetでsortしない | member limitの適用位置と全HSP retentionを維持する必要がある。native OFFの割合は高いが、文字列/カテゴリ/数値ID、同点、元順序とdtype、public selectorとの契約を崩す危険に比べ、残存sort/copyの削減量が未分離。有限/無制限別のsort/copy行数と実時間が不足。coercionと一緒に大改修しない |
| best-hit/RBH / `_select_orthogroup_edges`, `select_best_hits_per_query`, reciprocal selectors | OFFの同じmember tableを再coerce・再sort。forward/reverse reciprocityは必須。`one_to_one`はsubject順sortも必要。projectionにはoriginal rowsの保持が必要 | top-K順位を横流しする新しいprepared-result境界を作らない。public関数は単独呼出を維持する。S07.6後にOFFがなおこのownerで支配されることと、局所の冗長drop_duplicates等だけの効果を証明してから。今回に新しい高速path/flagを追加しない |
| group rank / `_anchor_core_hit_rank`, `_best_rows_by_query_target_record`, `_select_anchor_core_edges`, `_select_record_local_paralog_edges`, `_build_anchor_core_orthogroups` | cross-record bucket sort、direction全体の複数sort、rank再計算、edge/member/annotation objects。rank自身4.89–6.79%で親の19–27%とは違う。仮の全rank dictはO(E)の9要素tuple等を追加保持 | stable sortをmin/top2へ替えるには各consumerが読む範囲とtie証明が必要。親group builderにはname metadata・graph等も含まれる。rankのcall-site別重複数、rank一回保持のpeak、巨大groupのsupport tuple再copy量は未測定。現在のscopeへは採用しない |
| S03 HSP aggregation | 既にper-pair DataFrameを除去。残る一回のHSP走査、pair accumulator、区間union sort、代表rank、output materializationは意味上必要。Vibrio ONでは約27%のprofileを占める | 比率が大きいことはS03未実施を意味しない。interval/rankの追加変更は別の正確性領域。重複仕事の新しい実証がある場合のみ |
| S04 support/metadata | coreと拡張membershipの二snapshotに対して別索引を構築する仕事は必要。Hep ONは2回/13,584 input visits/423 evidence reduction visits、Vibrio ONは2回/107,814/20,752 | 二snapshotを統合すると所属連鎖が変わる。S04を巻き戻して総当たりへ戻さない。current graph/name/edge生成の残存負担とは分ける |
| unit内alias list membership、locus取得 | Hep 36,777 / Vibrio 312,351 alias追加試行。list長合計141,951 / 1,224,403は比較回数の上限的診断で、実比較数ではない | 通常unitのalias数は小さく、set追加が得か未確認。代表/優先順位と文字列正規化を動かす追加変更は不要。大locus多alias入力でここが支配された場合に検討 |
| S06 graph / S07 merge追加索引 | graphのexact DPと明示tuple API、strict rectangle・exact count・merge順・singletonは維持。既存Galleryではaccepted merge/conflict callが0 | 全列挙除去とchain改善は完了・了承済み。dense mergeのexclusion setが問題になる代表設定が新たに観測されるまで後回し |
| raw検索、render/save/transfer/cache | 現在のprofile境界の外。S06 browser helper観測には起動・最終描画を含まない | このnative profileから律速や改善率を断定できない。S08で実境界を分けて観測。LOSAT runtime変更、S05 cache/Worker変更は本計画の依存にしない |

見送りは永久に最適化不能という結論ではない。新しい根拠を得て別の小さな計画として扱う。
今回の二sessionを終えても効果が小さい場合、残りの全候補を自動で実装しない。

groupの仕事量を具体化すると、Hep ON / Vibrio ON / Similarityのdedupe後direction evidenceは
6,792 / 53,907 / 13,760件。現行sourceはanchor選択、record-local選択、related-edge選択で
それぞれこの全体をrank順にsortする（3E件のkey評価と3回のO(E log E) sortに、
bucketやtie判定のrank評価が加わる）。同じimmutable rowの順位再計算は削減候補だが、
phaseごとの走査・membership判定は必要であり、shared rank保持によるO(E)追加memoryと
既存pandas row objectの寿命延長をまだ比較していない。field数を減らした一時frameと異なり、
新しい全rank索引を今すぐ採用する根拠にはしない。

## 7. 同値性を証明する方法

### 7.1 S07.6の独立oracle

小規模の読みやすいreferenceをtests/prototypesに隔離する。入力順を意図的にshuffleし、
locus grouping、representative、alias→unit集合を素直なlist/dict/setで計算する。
本番で削除するsort・中間instanceの再現だけをoracleの判断根拠にしない。
f52ed68eの旧builderをtest-onlyで凍結した完全differentialも併用する。
単なるunit数ではなく`CollinearityUnitIndex`全field、dict順、tuple順、ambiguity、
例外type/message、警告内容、入力不変性を比較する。

- auto/cds/locus、empty record、records省略/短いrecords列、複数recordの同名alias。
- gene_parent_id > locus_tag > gene_id > GeneIDの優先、空白/欠損、同じgene labelで異なるlocus。
- 一locusに多数CDS、同座標・同長・source_protein_id有無・feature_index・protein_idの各tie。
- sorted member順とgroupの四独立minの違いが出る交差座標、reverse/mixed/unknown strand。
- aliasのunit内重複、二unit衝突、三unit衝突、曖昧化後の再出現、空値。
- global unit counter、全属性、最後のprotein mapping更新順、代表SVG IDを含む。
- structural guard: U回だけのunit生成、member再sortなし、per-alias set構築なし。
  wall-clock assertionをpytestへ入れない。

既存`tests/test_collinearity_units.py`の4 testsだけではtie/alias衝突の上記全体を
証明できないため、これらは実装sessionで追加する。Galleryはcollapsed unit 0なので、
multi-CDS locus syntheticを省略してはいけない。

### 7.2 S07.7の独立oracle

小規模list oracleでglobal length sort→bin→score sort→採択行ordinalを計算し、
同じ順でscalar log/sumを行う。前後のfit点tuple、slope/interceptのfloat値、
fallback有無、normalized_score、coverage/domain metadataを**完全一致**で比較する。
許容誤差でfloat計算順の変化を隠さない。fit行に使わない列を様々に変えても
fit結果が不変で、最終normalized DataFrameにはそれらの列が元どおり残ることを確認する。

- fit開始の最小行数前後、4/8 bin境界と端数、fraction ceilの境界、full tieと元ordinal。
- 等しい長さ・9桁round後の同一x・3点未満・denominator epsilon以下・非有限modelのfallback。
- string/categorical/numeric IDを持つ許容入力、extra columns、非連続/重複index、empty。
- 不明/欠損protein ID、zero length、reverse/overlap/disjoint HSP、coverage閾値直前直後。
- public numeric rejection（nonnumeric/NaN/±inf）、有限/無制限member選択、全HSP保持。
- 同じmodel下の`10 ** (slope * math.log10(length_product) + intercept)`、sqrt fallbackと
  `_row_float`呼出順は元のまま。vector演算への置換はscope外。

今回の48 checksは正規化済み小入力の局所oracle。上記全boundaryを既に実施したとは扱わない。

### 7.3 結果全体・既存契約

両sessionで、独立した小規模caseから全group/member/role/confidence/support/edge属性、
group/path/block ID・順序、DF columns/dtypes/index/values、typed JSON、semantic attrsと
geometryを旧sourceと比較する。core snapshotと拡張後membershipを別々に保つ。
S06の小規模exhaustive oracleと通常経路のmaterialization禁止を維持し、
S07のstrict境界・exact count・merge順・singletonの既存testsもそのまま通す。
member/raw limits、search DB scope・identityは変更しない。

## 8. 実装時の検証と測定

### 8.1 focused gateと広いgate

S07.6の最初のgate:

```bash
PYTHONPATH=. pytest tests/test_collinearity_units.py tests/test_collinearity.py -v
```

S07.7の最初のgate（新しいfit oracleを同じfileまたは専用の実在test fileへ追加）:

```bash
PYTHONPATH=. pytest tests/test_protein_colinearity.py tests/test_sparse_support.py -v
```

最終native candidateに対する共通gate:

```bash
PYTHONPATH=. pytest tests/test_protein_colinearity.py tests/test_collinearity.py tests/test_collinearity_units.py tests/test_lossless_ortholog_paths.py tests/test_ortholog_path_contract.py tests/test_protein_comparison_benchmark.py -v
PYTHONPATH=. pytest tests/test_web_feature_catalog.py tests/test_session_request_codec.py tests/test_session_compat.py -v
pytest tests/test_output_comparison.py::TestOutputComparison -v
ruff check gbdraw/
```

新規prototype/test/runnerもruff対象に含める。図のgeometryは変わらない前提なので、
referenceを更新して差を解消しない。test commandは最低30分を許容しincrementalに監視、
test-owned timeoutは変更しない。
全体not-slow gateは原則S08の統合revisionで一度実施する。
局所改修がshared caller/public validationへ拡大した場合は、そのsessionでも広いgateが必要。
通った無関係なsuiteをsession名だけで再実行しない。

### 8.2 native測定の受入条件

- input/settings/source/dependency hashを照合してから既存baselineを選ぶ。
  full-file差だけで無効化せず、測定callableと到達依存の差を確認する。
- ユーザーの後続指示により、今後の時間測定は全caseでwarmup 1 / **3 samples**を標準とする。
  旧手順のHep 21、Vibrio/unit 7という要求は適用しない。3は測定予算であり、分散・
  検定力・必要精度から導いた数ではない。各値・中央値・MADを報告し、統計的確証を約束しない。
  過去の7/21回の実測値と当時の判定は書き換えない。既存reportの比較は保存されたpolicyを使う。
  noiseだけを理由に自動で回数を増やさない。profile/counter、memoryは各case 1回。
  unit/fit単体は原因帰属が必要な場合だけ測り、prototype値を本番合格baselineに流用しない。
- timing、profile/counter、tracemallocは別run。timingにこの作業のbuild、checkout、重いtests、
  別benchmarkを重ねない。開始前に実hostを確認し、外部benchmarkが稼働中ならtimingの受入測定を
  保留して独立作業を進める。外部processを停止・変更しない。実行後に競合が判明した観測は残す。
- 全sample、median、MAD、warmup、CPU affinity、PYTHONHASHSEED、依存version、開始終了、
  コマンド、exit、host観測、input/output hashを保存する。median悪化**>10%**をflag、
  どちらかのMAD/median **>5%**ならinconclusive。都合のよいsamplesを選ばず、未確立のまま報告する。
- Similarity既存baselineの5.052% noiseは残る。独立した受入判断に追加測定が必要なら、
  別途ユーザーが追加を指示した場合だけ、そのcase/同じ境界の両sourceを一組測る。S03〜S07全体は再測定しない。
- output差は速度に関係なく不合格。memoryはinputの外でstartし、return outputを保持した
  peak/retainedを記録。減らない場合も実値を残し、RSSやbrowser memoryへ読み替えない。

基本の既存runnerコマンド（post-search境界だけ）:

```bash
taskset -c 3 env PYTHONHASHSEED=0 python tools/benchmark_protein_comparison.py run --source-root . --cases gallery-collinear gallery-collinear-off gallery-orthogroup --stages post_search --measure timing --warmups 1 --samples 3 --output <new-current-report.json.gz>
taskset -c 3 env PYTHONHASHSEED=0 python tools/benchmark_protein_comparison.py run --source-root . --cases gallery-collinear-vibrio gallery-collinear-vibrio-off --stages post_search --measure timing --warmups 1 --samples 3 --output <new-vibrio-current-report.json.gz>
```

S07.6だけならSimilarityの時間を新規測定する必要はない。
S07.7だけならOFFの不通過はsource/counterとfocused testで確認でき、unchanged OFFの
21-sample再測定を機械的に追加しない。必要なstageを同じrunnerへ追加し、独立した
別benchmark frameworkを作らない。baseline artifactは新規結果と別名にする。

representative correctness matrixは両Collinear scope、ON/OFF、member=1/5/None、
unit=auto/cds/locus、anchor=rbh/one_to_one/all、forward/reverse、small/sparse/dense。
全直積を時間測定する必要はない。Gallery default比較のほか、正確性に必要な小規模caseと
下記の不足する実recipe境界を選ぶ。

### 8.3 Web/CLIをどこまで検証するか

今回の計画では実行不要。Python実装後は生成wheelのsource一致を確認し、既存
`tests/run_cluster_merge_browser_acceptance.py`（明示Collinear）、
`tests/run_hsp_browser_acceptance.py`（raw reuse・失敗/retry）を、対象modeに合う範囲で使う。
Gallery名だけでmodeを断定せず、各Generateのrequest/helper/provenanceでmode・推論・
scope・raw/member limit・block paramsをassertする。SimilarityのOFFを作ってはならない。

二sessionを同じ未公開candidateへ続ける場合、S07.6はCollinear ON/OFFのreal helperと
typed result確認、S07.7完了時に両modeの最終wheelでGenerate/save/fresh-load/regenerateを確認。
最終S08でdesktop/narrow、保存recipeのVibrio 3/2/2とHepの明示Collinear、
Similarity明示5/無制限を含め、offline・geometry・cancel後raw reuseまで統合確認する。
各sessionを独立公開単位にするならその都度必要なreal-browser gateを省略しない。

CLIはsourceとは別の専用installから、checkout外の作業directoryでPYTHONPATHを外し、
現在のgraph/typed payloadを含むSessionと代表旧Sessionをreplayする。
保存済みanalysisのreplayだけでは新しいunit/fit計算を通らないため、prepared rawを使う
新規Python/typed生成も併用する。Session replay時間は解析時間として報告しない。

PlaywrightはCLIとPythonの両方を確認し、Node specsには`@playwright/test`を確認する。
NodeがなければPythonでtargeted check、Chromium/socketのsandbox拒否は同じcheckを
必要な権限で再実行する。wheelはignoredの生成物、public Gallery/mediaは更新しない。
共有render/Workerを変更しない二局所改修ではCircular再測定を毎回増やさず、S08で共通gateを行う。

Web性能を主張するなら、起動、raw検索、raw reuse後helper解析、最終render/preview、
saveを別区間で測る。cold/既存converted-cache hit/derived missを区別する。
既存S06 browser timingは異なるsource・境界・競合条件なので新しいWeb改善率の分母にしない。

## 9. セッション分割・終了条件・rollback

| session | 範囲 / 依存 | 終了条件 |
|---|---|---|
| S07.5（今回） | S07証拠とsource調査、二局所案を選択 | 本計画・MASTER・局所診断・handoffが完成。本番変更0 |
| S07.6 | unit builderの三変更。S07 sourceを基準 | 独立unit oracle、全index一致、focused/関連gate、counter、isolatedとCollinear全体time/memory、適切なreal helper確認、結果report |
| S07.7 | 4列fit投影。unit変更とは機能上独立、通常はS07.6後に実施 | fit点/係数/最終結果完全一致、boundary/gate、最終source測定、両modeの実browser/CLI範囲、結果report |
| S08 | S03/S04/S06/S07と、採用したS07.6/S07.7の統合 | 既存S08の回帰・性能・lifecycle検証、実際の待ち時間への寄与と留保を報告 |

S07.6/S07.7は一session一commit相当の変更単位とし、先行のsource/測定hashをhandoffへ固定。
同じownerを二重に作らず、superseded production経路はそのdiffで削除する。
owner/path/public API/schema/compatibility branchは増えないので通常の非増加ratchet証拠で足りる。
S06の互換性例外を今回の局所改修で解消済みと宣言しない。

時間だけでなく構造削減と同値性を満たすこと。安定した代表caseで回帰した変更は原因を
限定して修正し、解決できなければ当該sessionの本番差分を戻して見送りreportを残す。
noise/競合で性能合格を判断できない場合は実装完了と性能判定未確立を分ける。
S07の合格はこの結果によって変わらない。
採用/見送り/未確立を明示すればS08へ渡せるが、未確立を性能合格として数えない。
次sessionを開始する許可や外部公開は本計画そのものから推定しない。

rollbackはunit三変更、fit一変更を各session単位で戻す。保存形式・raw identityに変更がなく、
data migrationやcache key更新は不要。戻したsourceで必要なwheel/専用CLIを再生成する。
共有ツリーのstash/reset/cleanや他sessionのrevertを行わない。

## 10. S08への受け渡しと実施しない作業

S08は追加最適化のsessionではなく統合検証である。S07.6/7のdispositionと測定に使った
sourceを揃え、変わった依存・欠けた境界だけを新たに測る。
S08旧promptの「S03〜S05」や「共有LRUが残っていないか」はS05却下以前の記述であり、
S05実装・LRU置換を前提にしない。本書と更新MASTERの依存表を現在の指示として扱う。
49/64/81 cacheやWorker lifecycleは既存挙動の回帰確認であって、新cache設計の受入試験ではない。

S08が確認する性能仮説は、(a) OFFのunit仕事削減がnative post-searchに現れるか、
(b) ON/Similarityのfit削減が全体にも残るか、(c) 同じ設定・cache状態のWeb Generate待ち時間に
寄与するか、の三段階。局所の27–32%を(c)の数字として報告しない。
raw検索が支配する場合も既存raw設定を変えて成功扱いにしない。

本セッションでは本番コード、tests/reference_outputs、Gallery assets、authority、
S08の実行を変更しない。push、PR、mergeによる統合、tag、deployを行わない。

## 11. 次の実装session用プロンプト

以下をそのまま次sessionへ渡す。対象は**S07.6だけ**で、S07.7/S08の一括実行ではない。

```text
gbdrawのCollinear / Similarity groups性能改修について、S07.6「unit索引の重複生成削減」を実装してください。

Repository: /mnt/c/Users/genom/GitHub/gbdraw
計画worktree: /mnt/c/Users/genom/GitHub/gbdraw/.worktrees/collinear-s07-5-20260916
計画: docs/internal/collinear_similarity_performance_plan_2026-09-15/results/S07_5_PLAN.md
基準本番commit: f52ed68e15e527b08e4d9247b990684040fe4270
S07 worktree: /mnt/c/Users/genom/GitHub/gbdraw/.worktrees/collinear-s07-20260916

最初にAGENTS.md/CLAUDE.md、Web CLAUDE、最新ratchets、MASTER_PLAN、S07/S07.5を読み、
branch/HEAD/upstream/dirtyを確認してください。共有変更をstash/reset/cleanで隠さず保持してください。
originをfetchし最新origin/devからupstreamなし専用branch/worktreeを.worktrees配下に作り、
必要な未統合依存だけをpatch同値性確認後に順番に引き継いでください。
計画成果が未commitなら上記計画worktreeの文書・証拠を読んで使い、旧/tmpを必須依存にしないでください。

S07はユーザーにより性能合格済みです。6 inconclusive等の数値・条件は残し、合格を撤回したり
独立再測定を開始条件にしたりしないでください。S06の全path列挙除去の主目的は了承済みで、
互換性Reviewは将来のmerge前の課題です。S05のclass/cache/Worker/runtime変更は復活させないでください。

実装範囲はgbdraw/analysis/collinearity_units.pyの既存ownerに限定します。
1. aliases空の中間CollinearityUnitと__dict__再構築を除き、aliases込みで一回生成する。
2. sorted_proteinsの部分列であるmembers順を使い、cds_members作成時の同一key再sortを除く。
3. alias_targetsのper-alias setを除き、既存unique dict/ambiguous setへ一回で還元する。
   二回目の別unit衝突で曖昧化し、三回目以降も復活させない。unique dictの順序を保つ。
最初のprotein/group sort、四つの独立min、locus優先、代表順位、全属性、例外・警告は維持してください。
Similarityはこのunit builderを呼ばないため、Similarityの性能改善と報告しないでください。

計画§7.1の独立小規模oracleと凍結differentialで全indexを比較し、auto/cds/locus、
multi-CDS locus、tie、reverse、欠損、unit内重複、二/三unit alias衝突、dict順を検証してください。
計画§8に従ってfocused/関連gateと必要なreal helper検証を実施してください。
test commandには最低30分を許容し、test-owned timeoutやreference SVGを変更しないでください。

測定は既存S07 post_searchを参照し、新規の時間測定は3 samplesを標準にしてください。
unit isolatedは必要な場合だけbaseline/currentをwarmup1・3 samplesで測り、counterとmemoryは別run各1回。
過去の保存policyによる判定を変更せず、noiseによる自動再測定や回数増加はしないでください。
入力設定・source依存hash・sample・median/MAD・retained/peak・完全出力hashを保存してください。
中央値悪化>10%、MAD/median>5%の基準は維持。timingにbuild/checkout/tests/別benchmarkを重ねず、
外部作業を停止しないでください。既存native設定と保存Gallery実recipeを区別してください。

S07.7のfit改修、numeric/member/RBH/rankの追加最適化、S08は実行しないでください。
新しいmaterialな効果変更が必要な部分だけProduct判断手順へ戻し、既承認事項は再質問しないでください。
結果・変更owner・同値性・測定の限界・rollbackをresults/S07_6.mdへ記録し、MASTERを更新してください。
production/tests/docs/generated evidenceを別々に監査し、English proposed commit titleとsummaryを付けてください。
push/PR/merge/tag/deployは行わないでください。
```

S07.7を別途依頼する場合の実行プロンプト:

```text
gbdrawのS07.7「fit行選択の一時DataFrame削減」を実装してください。
Repository: /mnt/c/Users/genom/GitHub/gbdraw
S07.5計画worktree: /mnt/c/Users/genom/GitHub/gbdraw/.worktrees/collinear-s07-5-20260916
計画: docs/internal/collinear_similarity_performance_plan_2026-09-15/results/S07_5_PLAN.md

上記計画、MASTER、S07/S07.6の実在するhandoff、repository規則とratchetsを読んでください。
S07.6未採用ならそのdispositionを記録し、unit変更を仮定せずS07 sourceから進めて構いません。
fetch後の最新origin/dev由来の専用worktreeを.worktrees配下へ作り、必要な未統合依存だけ引き継ぎ、
既存dirtyを保持してください。S07合格、S06了承とmerge前Review、S05却下は維持してください。

本番範囲はprotein_colinearity.py::_select_normalized_fit_rowsだけです。
非empty入力からlength_product/query/subject/bitscoreの4列を一時投影し、
同じpandas sort/bin/ceil/head/itertuples/scalar logを使ってfit点を返してください。
元の正規化table全列を維持し、fit係数、sum順、score計算、fallback、HSP、member選択を変えないでください。
計画専用prototypeは参照証拠であり本番からimportしないでください。

計画§7.2/7.3の独立oracle・完全float一致・境界検証、§8のfocused/関連gateと最終sourceの
両mode real-browser/CLI検証を行ってください。test-owned timeout・referenceを変更しないでください。
Hep ON/Similarity/Vibrio ONとも、新規測定は3 samplesを標準として互換baselineを参照してください。
noiseによる自動再測定や回数増加はせず、過去の判定と今回の参考値を分けてください。
Similarityの既存noise5.052%とprototype timingの外部LOSAT競合は残し、性能合格へ読み替えないでください。
必要な追加測定はその境界の両sourceだけに限定し、timing/profile/memoryを分離してください。
OFFはこのownerを通らず、無変更のOFF benchmarkを理由なく再測定しないでください。

numeric/member/RBH/rank改修やS05/S08へ進まないでください。
results/S07_7.mdとMASTERへ結果・同値性・counter/time/memory・限界・S08引き継ぎ・rollbackを記録し、
English proposed commit title/summaryを付けてください。push/PR/merge/tag/deployは行わないでください。
```

## 12. 今回の検証・引き継ぎ

今回の差分は計画文書二つ、計画専用診断とその証拠だけ。
本番・本番tests・public generated artifactsの差分は0。既存S07 evidenceは変更していない。
診断の48-case self-check、Gallery fit点/model一致、sample hash一致、別run memory、
source/evidence照合、診断ruffと文書リンク/diff検証を実施した。
実行していない後続gateを合格と記載しない。
最終診断scriptは測定時から未使用の`ExitStack` importを一行削除しただけで、
測定関数・fixture・処理境界は同じ。report内の測定時script hashは書き換えず、
[最終照合記録](data/s07-5-review.json)に両hashと再構成検証を残す。

再現（worktree rootから。timingは実hostの競合を確認してから実行）:

```bash
python docs/internal/collinear_similarity_performance_plan_2026-09-15/results/data/s07-5-diagnose.py --measure inventory --output <inventory.json.gz>
taskset -c 3 env PYTHONHASHSEED=0 python docs/internal/collinear_similarity_performance_plan_2026-09-15/results/data/s07-5-diagnose.py --measure timing --output <fit-timing.json.gz>
python docs/internal/collinear_similarity_performance_plan_2026-09-15/results/data/s07-5-diagnose.py --measure memory --output <fit-memory.json.gz>
```

今回の成果は専用worktreeの未コミット差分として引き継ぐ。依存S07まではcommit済み。
新規commit、push、PR、統合merge、tag、deployは行っていない。

**Proposed commit title:** `Plan targeted Gallery analysis optimizations after S07`

**English summary:** Prioritize equivalent unit-index and fit-row improvements using
archived profiles and focused diagnostics. Define independent oracles, measurement
boundaries, and implementation handoffs while preserving S07 acceptance.
