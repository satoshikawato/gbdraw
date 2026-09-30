# Web GUI 監査 65 件の修正提案（dev 4c89bab1）

Status: 提案（判断は確定済み）。runtime、Product Contract、CI、checker は変更していない。この文書は Product authority ではない。

> 判断は [02_DECISION_PACK.md](02_DECISION_PACK.md) で確定した（2026-09-30）。実装の手順と PR の分け方は [03_IMPLEMENTATION_REFERENCE.md](03_IMPLEMENTATION_REFERENCE.md) に従う。
> この文書の推奨と 02 の決定が異なる項目では、02 に従う。異なる項目は次の 3 つ。
> - CO-04: B。ファイル単位の database を維持し、要求していない自己検索を除く。
> - IN-06: B。record ごとに自分の推定値で定義を付ける。
> - Dinucleotide: U を許し、T として扱う。
>
> 第 5 節と第 6 節の表は、判断を求めたときの記録として残す。
作成日: 2026-09-30。基準: `origin/dev` `4c89bab1`。行番号はすべてこの commit のもの。
読者: gbdraw の Owner（Product Decision Owner）と、修正を実装するセッション。

- 監査結果: [README.md](README.md)（65 件、ID は同じものを使う）
- 各 workstream の詳細と根拠: [remediation/](remediation/)。この文書の各項目の末尾に対応するファイルを示す。
- 付録は設計担当サブエージェントの報告そのもので、行番号や実験の記録を含む。authority ではない。

JS のパスは `gbdraw/web/js/` からの相対パスで書く。それ以外はリポジトリ直下からのパス。
分類の略号は次のとおり。これは [Product Impact Ratchet](../PRODUCT_IMPACT_RATCHET.md) の procedural classification である。

| 略号 | 意味 |
|---|---|
| IEA | IMPLEMENT_EXISTING_AUTHORITY。既存の authority が結果を 1 つに決めている。判断を待たずに実装できる |
| PDR | PRODUCT_DECISION_REQUIRED。製品として妥当な結果が 2 つ以上ある。Owner の決定まで実装を止める |
| ER | EVIDENCE_REQUIRED。判断の前に決定的な証拠が必要 |
| 決定済 | この会話で Owner が選んだもの（第 5 節） |

## 1. 原則の読み替え

SOLID、KISS、DRY、YAGNI は、このリポジトリの既存の規約（`CLAUDE.md` の Key Architecture、`gbdraw/web/CLAUDE.md`、[architecture ratchet](../ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)）に沿って次のように読む。

| 原則 | 読み替え | 主な規則 |
|---|---|---|
| SRP | 値と遷移ごとに semantic owner を 1 つにする。修正は呼び出し側ではなく owner で行う | R2、R3、R6、R10 |
| OCP | 能力は既存の境界に data として足す（表の行、plan の op、共有ベクタ）。分岐は足さない | R6、R7、R9 |
| LSP | 同じ要求は、Circular/Linear、single/grid/batch、Web/CLI で同じ結果になる | R4、R5、R8 |
| ISP | helper は必要最小限の入力（committed request や Result）だけを受け取る | R1 |
| DIP | draft ではなく、canonical request と committed Result に依存する | R1、R5 |
| KISS | 最小の変更にする。guard を足すより、分岐した経路を消す | IN-01、GE-02 |
| DRY | 同じ事実の複製（parser、評価器、reader、翻訳）を 1 つにし、旧経路は同じ PR で消す | R4、R9 |
| YAGNI | 新しい framework、schema field、option は、2 つ以上の実経路を統合し、旧経路を消すときだけ足す | 全体 |

## 2. 要旨

- **65 件の分類**（第 7 節の各項目の分類を集計したもの。複数の部分に分かれる項目は主な部分で数える）
  - IEA: 48 件。判断を待たずに実装できる。
  - PDR: 13 件。うち 4 件は決定済（IN-01、GE-02、FE-06、PV-04）。
    - IEA の項目の一部として、TR-10(b) と D-TR02 も決定済。
  - ER: 1 件（PV-10）。
  - 欠陥ではない、または受容済みの残余リスク: 3 件（PV-01、PV-03、PV-07。PD-OI-052 による）。
    - PV-03 と PV-07 は、機能要望として第 6 節の判断事項にも挙げる。
- **原因の共通化。** 65 件の原因は 10 の設計規則（第 4 節）にまとまる。各規則の強制手段（テスト、checker、CLAUDE.md の不変条件）で、同じクラスの再発を防ぐ。
- **CI が緑だった理由**（[W8](remediation/W8_verification_workflow.md)）。
  - 既存テストの多くは「Result・表示・出力が互いに一致するか」を見るだけで、fresh Generate や CLI と比べていない。
  - fixture が単一 Result、単一 record、正常な入力に偏っている。
  - 契約（OIPC-C01/C03/C04）を汎用に実行するテストがない。
- **設計段階で見つかった重大な訂正**（第 3 節）。
  - PV-08 の真因は、Multi-Record Canvas が depth 入力なしでも空の depth スロットを確保すること。Web 既定の Circular 出力はすべて、成功した図も含めて単一 record 経路と幾何が違う。私が dev の CLI で確認した。
  - CO-05 は conservation の reader にも及ぶ。
- **進め方**（第 8 節）。
  - Phase 0: 判断を一括で求め、guard と検査基盤を用意する。
  - Phase 1: P1。
  - Phase 2: P2。
  - Phase 3: P3。
  - Phase 4: 仕上げ。
  - 合計約 50 PR。ホットなファイル（`app/run-analysis.js`、`app/app-setup.js`、`services/config.js`、`services/session-request.js`）を触る PR は直列にする。

## 3. 監査結果の訂正と、新たに見つかった不具合

### 3.1 訂正

| ID | 監査の記述 | 訂正 | 出典 |
|---|---|---|---|
| PV-08 | grid が定義文字列のために中央を確保する（推定） | 定義用の予約半径は単一 record と同じ。真因は `gbdraw/api/diagram.py:2897, 2994` が depth 入力なしでも `_precomputed_depth_tracks=[]` を渡し、`:2264` の `is not None` で `show_depth` が True になること。37 px の空の depth スロットが入る（`021b4c49` からの退行。main にも含まれる）。Salmonella だけは、tuckin と separate strands の組み合わせで定義の円が入りきらないという第 2 の原因もある。「CLI の既定では成功」は不正確で、単一キャンバスでも `--separate_strands` で失敗する | W6、私も確認 |
| CO-05 | `gbdraw/analysis/conservation.py:147-149` は対処済み | 対処は DataFrame 入力だけ。ファイル入力（`:175-180`）も同じようにずれる | W6、W8 |
| IN-02 | 原因は `app/record-discovery.js:137-138` | 同じ正規表現が `app/run-analysis.js:772-774`（LOSAT の FASTA 抽出）と `app/match-sequences.js:820-822` にもある。record-discovery だけを直すと、Prokka 形式の入力で LOSAT が壊れる | W2 |
| IN-03 | 原因は JS の高速経路 | Worker の helper `list_gff_fasta_records`（`app/python-helpers.js:1645-1660`）も、`09c99fb8` から FASTA を列挙するだけになっている。落ちるのは「GFF の行を 1 つも持たない FASTA 配列」で、feature が 0 でも埋め込み `##FASTA` や region 行があれば残る | W2 |
| IN-06 | `app/linear-sources.js:247-255, 271-276` | 行番号が誤り（ファイルは 235 行）。該当箇所は `:107-145` | W2 |
| FE-05 | 回避策は target record を消すこと | batch では Python が record 未指定を拒否するので、回避策にならない。有効な設定が 1 つもない | W2 |
| FE-06 | IUPAC 文字だけの遺伝子名が原因 | All が配列と `/translation` を含むこと自体が原因。IUPAC の展開は増幅にすぎない。Interactive SVG は同じ実装ではなく、手で移植した重複 | W7 |
| PV-01 | P2（P1 に近い） | PD-OI-052 が受容した残余リスク（clipping は自動 clamp しない）。欠陥ではない | W1b |
| PV-03、PV-07 | 欠陥 | PD-OI-052 が「padding と legend 順の継承は保証しない」と明記している。機能要望 | W1b |
| TR-09 | 初回の stack 有効化で ticks が入る | 初回の件は UI 契約どおり（「saved custom stack を使う」）。欠陥は preset Reset の件だけ | W7 |
| TR-10(b) | help tip が hover 専用 | `js/components.js:12-16` に明記された意図的な設計。変更は判断事項（決定済） | W7 |
| TR-12 | 生成器を特定する | 生成器はなく、手書き。相対パスの形式自体がホスト上では無効 | W7 |
| CO-01 | serial へ fallback しない | PD-OI-018 が「threaded は既定で、support check は残る」と決めている。欠陥は UNKNOWN 表示と stage の誤り | W5 |
| CO-06 | ID を照合しない | 位置での割り当ては設計どおり。欠陥は metadata の矛盾。ID 照合をするかは判断事項 | W5 |
| GE-03 | recipe の不具合 | Circular は Web の退行。0.13.0 は `--definition_font_size` を渡して間隔を `int(font+2)` にしていた | W3 |
| X-02 | CLI はほとんど拒否する | window/step の 0 と負値は CLI も受理し、空の track を描く。GE-08（P3）と TR-04（P2）は同じ原因なので P2 に揃える | W3 |
| FE-12 | P3 | 色の定義域が 3 通りある（Python、browser JS、DOM なし JS）。Session の往復が壊れる可能性もあり、P2 の候補 | W3 |
| IN-08 | Save の不具合 | Save の拒否は仕様（`docs/REFERENCE/web-app.md:210`）。欠陥は案内だけ | W2 |
| GE-02 | 失敗した Generate が「保持した」と表示する | その表示は正しい。欠陥は不正値を受け付けること | W1b |
| SE-05 | Save の失敗 | 0.13.0 では Load → Save が catalog なしで成立していた。リリース済みの継続操作を、決定記録なしに退役させている | W3 |
| SE-06 | `services/imported-comparison-intent.js:191` が原因 | ここを緩めると EDITABLE になり、Generate が保存された BLAST ではなく Web 既定の LOSAT を実行した（match 8715 → 112）。読み取り専用の扱いは PD-OI-008 どおりで正しい。本当の原因は 2 つ。(A) CLI の binding uid `cli-seq-N` と request の `record-N` の食い違い。(B) Inherit が schema 7 の比較を schema 8 に昇格しないこと（N-17） | W4 |
| SE-08 | tobacco の Session に限る | CLI encoder が派生値（`whitelist_map`、`priority_map`、`label_override_rules`）を `diagramOptions.config` に書くことが原因。Linear の CLI Session はすべて該当する。読み込みの遅さは lazy Worker の規則違反ではなく、PD-OI-009 の型付き検証が Python を要するため | W4 |
| SE-02 | ラベルのクリック | キーボード（checkbox の Space、radio の矢印）も同じ原因で Undo できない | W4 |
| SE-04 | select と range | range の input は存在しない | W4 |
| GE-06 | `state.js:793-802` が原因 | Generate 中の編集は仕様で許されている。破られているのは History の不変条件（transaction の間は before が current）。`services/history.js:1000-1049` の `undo`/`redo` が、開いている artifact transaction を見ない | W4 |

### 3.2 新たに見つかった不具合

| ID | 内容 | 重大度（案） | 担当 |
|---|---|---|---|
| N-01 | Web 既定（Multi-Record Canvas）の Circular 出力すべてに空の depth 帯が入る。GC content の幅が 74.1 → 57.4 px（HmmtDNA、tuckin）。PV-08 と同じ原因 | P1 | W6 |
| N-02 | `app/record-discovery.js:84` の ORGANISM 正規表現も改行をまたぎ、空の ORGANISM が `Unclassified.` になる | P3 | W2 |
| N-03 | `load_comparisons` が存在しないファイルや解析できないファイルを warning だけ出して飛ばすので、後ろの比較が 1 つ前の record 対に描かれる（`gbdraw/io/comparisons.py:66-71, 80-87`） | P1 | W6 |
| N-04 | Multi-Record Canvas で、Region annotation と depth slot の `legend_label` も凡例から消える | P2 | W6 |
| N-05 | GFF3 の phase を読まないので、LOSATP の protein が読み枠のずれた配列になる | P2 | W6 |
| N-06 | 色ルールの caption が既定の type 名と同じだと、CLI の凡例がルールの色を落とす（`gbdraw/legend/table.py:175-230`） | P2 | W1b |
| N-07 | crop したとき、範囲外の比較行が record の外に描かれる（`gbdraw/linear_comparison.py:25-33`、CLI で確認） | P2 | W5 |
| N-08 | アップロードした表と逆相補の組み合わせで ribbon を誤る（CLI で確認） | P1 | W5 |
| N-09 | `reverseComplementOverride` がカードの checkbox を上書きする（コードからの判断） | P3 | W5 |
| N-10 | LOSAT の job 数の見積もりが、実行の batch 計画と別に実装されている（`app/losat-settings.js:48-90`） | P3 | W5 |
| N-11 | Interactive SVG でも FE-11 と同じ座標の誤りがある（`services/standalone-interactivity-assets.js:2621-2624, 4475-4482, 5765-5770`、コードからの判断） | P3 | W2 |
| N-12 | UNKNOWN のエラーパネルは、Save が失敗したときにも Save Session ボタンを出す。`generate` action にはボタンがない（`index.html:4704-4718`） | P3 | W3 |
| N-13 | CLI の `-n XY` が受理され、平坦な track を描く | P3 | W3 |
| N-14 | `tests/web/helpers/session-regeneration-contract.cjs:10` の `number()` が、ない metric を 0 と読む。`exact(key, 0)` が metric なしでも通る（要確認） | 検査基盤 | W8 |
| N-15 | 文書と実装のずれ: `docs/internal/SELECTIVE_CI.md:76, 124` の smoke 件数、ratchet の「at most three rules」と実装の上限 4 | 文書 | W8 |
| N-16 | `runLabelReflow` が draft から request を組み、adopt せずに Result を置き換える（`app/run-analysis.js:1898, 4579, 4770`）。「Generate で反映」の draft がラベルの reflow で Result に漏れる可能性（未検証の仮説） | ER | W1b |
| N-17 | Inherit が、確定済みの比較を schema 8 に昇格せずに候補へコピーする（`services/imported-comparison-intent.js:449`）。main の CLI が書く v42 の sidecar（schema 7、`settings.alignOrthogroupFeature` が残る）を Python が拒否する。`origin/main` の CLI 出力で再現した | P2 | W4 |
| N-18 | legend と diagram の drag（`app/legend/drag-actions.js:87`、`app/legend-layout/diagram-drag.js:88`）も、text 欄にフォーカスがあると SE-03 と同じく transaction を共有しうる（コードからの判断） | P3 | W4 |
| N-19 | Session 読み込みの rollback（`services/config.js:3579, 3636`）でも catalog が複製され、admit されない。SE-01 と同じ原因の潜在経路（コードからの判断） | P2 | W4 |
| N-20 | checkpoint のたびに feature catalog 全体を JSON で複製して署名している。HmmtDNA の色変更 1 件で 1,397,698 bytes。参照で持つと 24% 減る（prototype で計測） | 性能 | W4 |

## 4. 設計規則

各規則は、1 つの semantic owner と 1 つの強制手段を持つ。新しい framework は作らず、既存の境界と既存のテスト基盤を使う。

### R1: Result の書き手を限る

即時の編集が現在の Result を変えてよいのは、次の 3 つの場合だけとする。

- (a) Generate の compiler（`app/candidate-render.js` の `compilePlanBundle` → `services/svg-result-ingestion.js`）が同じ executor で再適用する editor intent。
- (b) Python oracle との parity が保たれた composition owner の編集。
- (c) committed Session と宣言済みの投影 field から自動で再描画する場合。

どれにも当たらない設定は「Applies on Generate」と表示する。draft や helper による部分的な再描画は、Result に書かない。

- **対象:** IN-01、GE-02、PV-02、PV-10、N-16。
- **強制手段:**
  - privileged capability の "Mounted SVG/Result replacement"（`tools/web-change-policy.json`）を "Result content commit" と "SVG serialization" に分ける。前者の allowlist は縮小方向にしか動かさない。
  - live 編集と fresh Generate の parity テスト（第 8.2 節の G-A）。
  - SVG 断片を返す Python helper（`app/python-helpers.js:1855-1872` の登録表）を足すときは、architecture review を求める。
- **根拠:** web CLAUDE.md の live-edit invariants、OIC-024（表示の真実性）、PD-OI-052 の「同じ候補境界」。

### R2: editor intent の寿命

editor override は、表示の変化（Result 選択、mount、record 選択、mode、非表示、reflow）で作成、剪定、削除しない。削除してよいのは次の場合だけとする。

- 明示の操作（Reset、Import）
- Undo/Redo
- Session の置き換え
- 成功した source 置換 Generate の中での owner reconcile

これは OIPC-C05/C06 の Web 版である。

- **対象:** FE-01、IN-05、TR-06、SE-10。
- **強制手段:**
  - Node contract: feature 集合が重ならない SVG を順に bind しても、override の map が変わらない。
  - G-C（編集ではない操作の前後で、利用者が持つ状態が変わらない）。

### R3: domain ごとに投影は 1 つ

live action、History apply、Result mount は、同じ投影関数を呼ぶ。表示の母集合は、mounted Result の committed metadata から導く。同じ概念の選択 ref を 2 つ持たない。dialog の reactive object は表示用の値だけを持つ。

- **対象:** FE-02/PV-09、FE-03、FE-04、FE-10。
- **強制手段:** 各 action の直後に reconcile を呼んでも DOM が変わらないことを確かめる parity テスト。

### R4: fast path は「完全一致か decline」

Python が持つ事実を JS が先に出す場合は、canonical loader と完全に一致するか、判断を Worker に回す（decline）。

- 1 つの共有ベクタ（`tests/fixtures/record_metadata_inference_cases.json` を拡張）を、loader、Worker helper、JS の 3 者で実行する。
- 「Generate が読むもの」に答える Worker helper は canonical loader を呼ぶ。
- JS の GenBank 見出し parser は 1 つにする。
- JS が Python の照合を先回りする箇所（大文字小文字、色の定義域）は、Python と同じ同値関係にする。

- **対象:** IN-02、IN-03、N-02、FE-08、FE-12。
- **強制手段:**
  - 共有ベクタのテスト。
  - `tests/web/architecture-contracts.test.mjs` に、`/^VERSION` と `/^ACCESSION` の正規表現が owner 以外に現れないことを確かめる assertion。

### R5: 比較 evidence の座標系と再利用

座標系を 3 つに分ける。

| 記号 | 座標系 | 内容 |
|---|---|---|
| F | 探索座標系 | 選択と crop の後の配列。元の鎖の向き。LOSAT の出力そのもの |
| V | 表示座標系 | F に実効の逆相補を適用したもの |
| Src | source 座標系 | 入力ファイル上の座標 |

規則は次のとおり。

- raw と Save Raw は F のまま。
- 保存する evidence とアップロードの表を F に統一するか V に統一するかは、Q-FRAME の決定による（推奨は F）。
- 人が読む座標（popup、FASTA ヘッダ）は Src にする。
- Linear の向きの owner は `seq.region_reverse` 1 つにする。
- 再利用の判定は、抽出しなくても分かる入力をすべて data として比べる。invalidate イベントは最適化にとどめる。
- LOSAT の job 計画は、実行と見積もりで共通の純関数にする。

- **対象:** CO-02、CO-03、CO-04、CO-07、CO-10、N-07、N-08、N-09、N-10。
- **強制手段:** metamorphic テスト（G-F）。
  - ファイルのまとめ方を変えても結果が同じ。
  - rotate+orient と逆相補が同じ。
  - visibility を変えた再利用と新規実行が同じ。
  - Save Raw → 再アップロードが同じ。

### R6: 失敗の意味はエラーを出す側が持つ

利用者が直せる失敗は、code と範囲の決まった context で投げる。

- JS: `services/error-normalization.js` に `diagnosticError(code, context)` を 1 つ置く。
- Python: `GbdrawError(..., diagnostic=)` を使う。

文言は normalizer だけが持つ。

- message を照合する分類器は新しく足さない。
- context にある locator（Sequence N、Line N、Track row N、Depth series N、帯の px）は summary に出す。
- エラーの action は実際に押せるものだけにする。

- **対象:** X-01（IN-07、GE-05、TR-05、SE-09）と、SE-05、SE-06、PV-04、CO-01、IN-04、IN-08 の文言、N-12。
- **強制手段:**
  - 実際の producer を呼ぶ表駆動テスト。
  - `NATIVE_VALIDATIONS`、`nativeValidation` の正規表現、Python の `_EXACT`/`_TEMPLATES`/`_CONSTRAINTS` の件数を減る方向にしか動かさない ratchet。
  - 語彙の JS/Python parity テスト。
- **根拠:** PD-OI-046。

### R7: 値の検証は Python の typed layer が持つ

- JS の projection は値を変換しない。空欄は null、数は数、それ以外は typed の INPUT_INVALID にする。
- JS が Python より先に値を使う比較の閾値は、Python から生成した定義域（`gbdraw/mode_profiles.py` → `mode-profiles.generated.js`）で評価する。
- Generate は draft を書き換えない。

- **対象:** X-02（GE-04、GE-08、TR-04、CO-09）、N-13、GUI の外の CLI の不具合。
- **強制手段:**
  - 共有ベクタ `tests/fixtures/option_domain_vectors.json` を、CLI、typed request、Web で同じ reason になるかで検査する。
  - Generate 経路の draft 代入の件数（現在 80）を ratchet にする。
- **根拠:** OIPC-C01、C03。

### R8: Multi-Record Canvas は配置だけを持つ

- 凡例、depth の有無、slot の構成は、単一 record 経路と同じ関数から得る。
- 「入力がない（None）」と「入力はあるが空（[]）」を区別する。
- **対象:** TR-01、PV-08、N-01、N-04。
- **強制手段:** 1 record の Multi-Record Canvas と単一 record の経路で、幾何、凡例、定義が一致することを確かめる parametrize した pytest。

### R9: 入力形式ごとに reader は 1 つ

- BLAST 表は `gbdraw/io/comparisons.py` の 1 つの reader と 1 つの正規化関数だけで読む。
- 読めない入力は、飛ばさずに ValidationError にする。
- **対象:** CO-05、N-03。
- **強制手段:** `names=COMPARISON_COLUMNS` が `gbdraw/` 内で 1 か所だけであることを確かめる静的テスト（前例: `tests/test_dead_api_cleanup.py`）。

### R10: watcher で状態を修復しない

- slot の reconcile は、入力を変える遷移から owner の関数を明示的に呼ぶ。
- 1 つの概念には 1 つの builder だけを置く。
- **対象:** TR-02、TR-03、TR-06、TR-09。
- **強制手段:** Circular と Linear に同じ表を流す table-driven の parity テスト。
- **根拠:** web CLAUDE.md の「watcher execution is not an invariant mechanism」。

### R11: History のトランザクション境界

- transaction は 1 つの owner（control または gesture）に属する。別の owner の transaction を始めるときは、開いているものを先に確定する。
- 「開いている intent を確定する」処理（`settlePendingIntent`）は 1 つにまとめる。`begin`、`beginCheckpoint`、`beginArtifactReplacement`、`runUndoableCommand`、`undo`、`redo` がこれを共有する。
- discrete control の transaction は、値を確定するイベントの capture phase で始める（checkbox と radio は `change`、button は `click`）。text 系は focus で始める。操作の手段（pointer、ラベル、キーボード）に左右されない。
- History の利用可否は reactive な判定 1 つにまとめ、`undo`、`redo`、`canUndo`、`canRedo` で共有する。
- Generate が所有する artifact（feature catalog）は、履歴項目でも復元経路でも参照で持つ。JSON で複製しない。`state.featureCatalog` は null か admit 済みの catalog だけにする。

- **対象:** SE-01、SE-02、SE-03、SE-04、GE-06、FE-04、N-18、N-19、N-20。
- **強制手段:** 入力イベントの組み合わせテスト（control の種類 × 操作の手段 × フォーカスの状態 × 処理中）。確認するのは、変更した control ごとに step がちょうど 1 つ増え、Undo が後入れ先出しで戻ること。
- **根拠:** web CLAUDE.md の History 規則、`docs/REFERENCE/web-app.md:706-707`。

### R12: 補助の規則（Web↔CLI、アクセシビリティ）

- **Source recipe:** 組み立てた argv を CLI の分割規則で読み戻し、元と一致しなければ理由付きで unavailable にする（GE-03、TR-07）。
- **アクセシビリティ:** 入力には、可視ラベルと一致する accessible name を付ける。placeholder や状態で変わる title は名前にしない（TR-10(a)）。help tip はキーボードとタップで開ける（TR-10(b)、決定済）。

## 5. Owner の決定（2026-09-30、Owner: satoshikawato）

この会話で Owner が選んだ結果を、選ばれた案の記述のまま記録する。
理由、保持すべき効果、退役範囲、受容リスクは与えられていないので空欄にする。テンプレートがこれらの推測での記入を禁じているため（[PRODUCT_DECISION_PACKET_TEMPLATE.md](../PRODUCT_DECISION_PACKET_TEMPLATE.md)）。
Product Contract に PD-OI-056 以降として登録するには、Owner がこの 4 項目を記入する必要がある。登録は authority だけの PR として出すか、実装と reviewed co-change で同時に出す。

| 項目 | Owner の発言 | 選ばれた結果 |
|---|---|---|
| FE-06 | 「FE-06: 案A」 | 検索の All から配列の内容（Nucleotide sequence、Amino acid sequence、`/translation` の値）を外す。配列の検索は専用の Nucleotide / Amino acid field で行い、IUPAC の展開はそこに残す。module 版と Interactive SVG 版を同時に変える |
| TR-10(b) | 「TR-10(b): 推奨案でいきましょう」 | A。すべての help tip を、キーボードとタップで開ける disclosure button にする。id を自動で振り、hover 専用の分岐を削除し、tip を label の外に移し、対象の control から `aria-describedby` で参照する |
| D-TR02 | 「D-TR02: 足す」 | A。Circular と Linear の両方で、論理 depth series が最初の source を得たとき、その index を参照する行（有効・無効を問わない）がなければ managed 行を 1 つ足す。series が source を失ったら managed 行を除く |
| IN-01 | 「IN-01: Aでいきましょう」「IN-01（P1）：推奨案でいきます」 | A。Circular の Species、Strain、Plot title、Title position、Title font、Default font size、Keep Full Definition を「Applies on Generate」にする。定義行を即時に作り直す経路は削除する |
| GE-02 | 「GE-02: 推奨案でいきましょう」 | C。全体の stroke 設定の即時反映をやめ、「Applies on Generate」にする（`applyStylesToSvg` の live 経路を削除） |
| PV-04 | 「PV-04: 推奨案でいきましょう」 | A。feature のある凡例項目を既存の名前に変えるときも、既存の衝突ダイアログ（Merge / Suffix / Cancel）を出す。衝突先が色ルールの caption なら、PD-OI-042 の既存の区別処理に任せる |

receipt の下書き（1 件分の例。残りも同じ形）:

```text
PRODUCT_DECISION
Concern: <新しい concern key。登録時に決める>
Scenario revision: 1
Choice: A / <上の「選ばれた結果」>
Rationale: （Owner 未記入）
Must preserve: （Owner 未記入）
May retire: （Owner 未記入）
Accepted residual risk: （Owner 未記入）
Owner: satoshikawato
Decision date: 2026-09-30
```

## 6. 判断を求めた事項（記録。結果は 02 の「判断の記録」）

判断は 5 つの Pack にまとめて一度に提示する。推奨は engineering の推奨で、Product authority ではない。

### DP-1 即時編集と装飾の継承

| 項目 | 選択肢 | 推奨 |
|---|---|---|
| FE-02/PV-09 batch の編集範囲 | A: 編集時に全 Result へ伝播する / B: 表示時に正本から投影し直す / C: 現状を明示する | B。label は構造上 B しか取れないので、方式が 1 つで済む。既存の History after-apply の投影と統合できる |
| PV-03 legend 順 | A: 現状を表示する / B: 既存の candidate plan に `legendOrder` op を足す / C: request field にする | B。順序の意図は `editorState.legend.entries` に保存済みなので、新しい schema は要らない |
| PV-07 padding | A: 現状維持 / B: Session の `ui.canvasPadding` を候補の公開前に全出力へ適用する / C: request field にする | B。PD-OI-052 が示す緩和策（padding で調整）を持続させられる |
| PV-01 side 変更と delta | A: 現状維持（受容済み） / B: side 変更で delta を 0 にする / C: bbox が viewBox と交わらないときは公開前に止める | A。改善するなら PV-07 の B を先に入れる |
| FE-01 source 置換の細部（任意） | A: 対象単位で剪定する / B: 現状の all-or-nothing / C: Reset まで残す | A |

### DP-2 入力と record

| 項目 | 選択肢 | 推奨 |
|---|---|---|
| IN-06 1 ファイル内の複数 organism | A: 全 record で同じときだけ file 既定にする / B: record ごとに定義を入れる / C: 推定値の層を足す | A |
| IN-03 FASTA だけにある配列を描くか | 描く / 現状どおり落とす | 落とす（CLI の科学的出力を変えない） |
| IN-04 単一 record → 単一 record の置換 | A: 新規 upload と同じにする / B: 残す | B |
| FE-09 重複 record ID の「This feature only」 | A: 常に stable hash を出す / B: Python で instance id を照合する / C: 無効にして理由を示す | A |
| SE-10 Reset の範囲 | A: Linear の record ごとの表示状態と alignment plan を初期化する / B: Circular も残す / C: 明記する | A |
| IN-02 BOM を Python で受け付けるか | 受け付ける / 今は decline で揃える | 今は decline で揃え、受け付けるかは別途判断 |

### DP-3 比較

| 項目 | 選択肢 | 推奨 |
|---|---|---|
| Q-FRAME（CO-07、CO-03(b)、CO-10、N-08） | A: 表は V、Save Raw は V に変換 / B: V と F を明示 / C: すべての比較表を F にし、planner が投影する / D: source 絶対座標 | C。JS 側の F→V 変換と helper を削除できる。main の Session には reader で一度だけ変換する |
| CO-04 ファイルのまとめ方（PD-OI-018 の改訂） | A: DB を subject の 1 record に限る / B: batch を維持し自己検索を除く / C: limit で切り替える | A。CLI と一致する。native の計測で時間は 1.1〜2.5 倍 |
| CO-06 ID の照合 | A: 厳格 / B: 矛盾はエラー、不明は警告 / C: metadata だけ直す | B（Python の planner に置く） |
| CO-08 消えた group の名前 | A: member 集合の完全一致で付け替え、残りは dormant として保存する / B: 通知付きで削除する / C: 重なりで引き継ぐ | A |
| CO-10 popup の座標 | A: Src だけ / B: Src を主にし、違うときだけ表の座標も出す / C: ラベルだけ付ける | B |

### DP-4 Python 本体、CLI、Session

| 項目 | 選択肢 | 推奨 |
|---|---|---|
| CO-05 12 列を超える表 | A: ちょうど 12 列だけを受け付ける / B: 先頭 12 列が型どおりなら受け付ける | B。決定までは A で先に出せる |
| PV-08 定義文字列が入らないとき | A: エラーだけ改善する / B: 種名の行を単語で折り返して配置し直す / C: フォントを縮小する | B。決定までは A |
| SE-05 旧形式 Session の Save | A: Generate してから Save（Load 時に notice） / B: Save 時に catalog を回復する（ER） / C: Load 時に作る | A |
| Dinucleotide の文字集合 | `^[ACGT]{2}$`（大小無視） / N・U を許す | `^[ACGT]{2}$` |
| 負の offset、spacing、gap、rotation | ER（負値の意味を確かめる） | — |
| SE-08(b) CLI Session の読み込み時の Python 検証（tobacco で 9.9 秒、そのうち Worker の起動が 9.2 秒） | A: 読み込み時の検証を続ける / B: 最初に Python を使う操作まで遅らせる / C: CLI が full config ではなく差分の configOverrides を書く | 既存のファイルには A。C は、replay で結果が変わらないことを証明したうえで別の architecture PR にする。B は推奨しない。判断の前に証拠を用意する（ER） |
| SE-06 の比較を EDITABLE にするか | する / しない | しない。読み取り専用のまま、Inherit で使い続ける（PD-OI-008） |

### DP-5 Generate 中の操作

| 項目 | 選択肢 | 推奨 |
|---|---|---|
| GE-06 Generate 中の Undo/Redo | A: History の artifact replacement か checkpoint が開いている間は busy として拒否する / B: Undo/Redo の前に Cancel し、直前の Result を保ってから Undo する | A。ボタンにも同じ判定を使う。overlay にはすでに理由が出ている |
| overlay の後ろのキーボード編集 | 一貫して modal にする（inert） / 現状のまま編集を許す | modal にする。GE-06 と一緒に決める |

## 7. バグ別の修正案

各項目の書式は次のとおり。

- **分類:** 分類と、それを決める根拠。
- **修正:** owner と変更の内容。
- **削除:** 同じ PR で消す旧経路。
- **テスト:** 検出できる assertion。
- **PR:** 規模と、同梱する項目。

詳細な行番号と実験の記録は、括弧内の付録にある。

### 7.1 P1

**IN-01（=GE-01）Circular の定義編集が Result を誤った内容で書き換える**（[W1b](remediation/W1b_live_svg_legend.md)）
- 分類: 決定済（A）。現状の維持は OIC-009/010 と科学的出力の整合性に反する。
- 修正: `app/results.js` の即時の定義再生成を削除する。Species、Strain、Titles & Record Labels に、既存の「Applies on Generate」表記（`index.html:1472` と同じ span）を付ける。Linear は 2026-04 からこの形なので、両モードが揃う。
- 削除:
  - `app/results.js:14-55, 132-138, 177-232, 234-492`
  - `app/watchers.js:186-196, 574-585`
  - `app/python-helpers.js:1693-1827, 1866`
  - `workers/diagram-generation-worker.js:745-773` の helper op
  - 呼び出し元がなくなる `app/legend-layout/reposition-actions.js:155-173`、`app/legend-layout/composition-actions.js:850-894`
- テスト: Region と Record label を付けて Generate し、定義系の field を編集する。`results[0].content` と保存した Result がバイト一致すること。Generate 後は crop の長さ、GC%、label を保持すること。
- PR: 約 10 files、net 約 −600。Architecture profile、`architecture-change`。
  - 先に guard だけの PR が必要: `tests/web/architecture-contracts.test.mjs:307-314, 593` の厳密な件数を上限に変える。

**IN-02 空の ACCESSION/VERSION で record ID が `KEYWORDS` になる**（[W2](remediation/W2_records_topology.md)）
- 分類: IEA（PD-OI-044 の exact source-bound identity。loading は Python が持つ）。
- 修正:
  - `app/record-discovery.js` の見出し正規表現を `[ \t]+` にする。
  - `genbankHeaderIds(chunk)` を export し、`app/run-analysis.js:772-775` と `app/match-sequences.js:820-825` もこれを使う。
  - BOM で始まるファイルは decline して Worker に回す。
  - N-02（ORGANISM）も同じ PR で直す。
- テスト: 共有ベクタ（Prokka 形式、ACCESSION のみ、GI 付き VERSION、空の ORGANISM、BOM）を loader、helper、JS で実行する。Prokka 形式の単一 record で `contig_1.svg` が出ること。LOSAT の抽出も通ること。
- PR: IN-03 と同梱。4 files、churn 約 110。Ordinary。main への backport は判断事項。

**IN-03 GFF3+FASTA で、feature のない FASTA 配列が NO_MATCH になる**（W2）
- 分類: 検出を loader に合わせる部分は IEA。FASTA だけの配列を描くかは PDR（DP-2）。
- 修正:
  - `list_gff_fasta_records` を `load_gff_fasta` の結果に戻す（`09c99fb8` より前の意味）。
  - JS の高速経路は、GFF の seqid と FASTA の id の積を FASTA 順で返す。判断できない入力（8 列未満、id の重複など）では decline する。
- テスト: gffpair、sequence-region だけ、埋め込み `##FASTA`、region 型、GFF の seqid が FASTA にない、の各ベクタ。
- PR: IN-02 と同梱。

**PV-08（と N-01）Web 既定で長い /organism を描けない。成功した図の幾何も違う**（[W6](remediation/W6_python_core_export.md)）
- 分類: 主な原因とエラー文言は IEA（`docs/RELEASE_NOTES_0.14.0b0.md` の sparse depth の定義）。自動で回避するかは PDR（DP-4）。
- 修正:
  - `gbdraw/api/diagram.py:2994-2995` で、depth 入力がないときは `_precomputed_depth_tracks=None` を渡す。`[]` は「入力はあるが、この record のセルがない」だけを表す。
  - `gbdraw/diagrams/circular/radial_layout.py` の 2 つの raise を 1 つの生成関数にまとめ、実際に置けなかった slot と、定義帯が原因であることを示す。
  - `error_adapter.py` に TRACK_LAYOUT/DEFINITION_RESERVED を足す。
- テスト: depth 入力なしの Multi-Record Canvas に depth slot がないこと。1 record の Multi-Record Canvas と単一 record で、gc_content と gc_skew の widthPx が一致すること（R8）。Mycobacterium の名前が Web 既定で成功すること。
- PR: Python の数行とテスト。Review REQUIRED（Web 既定の Circular 出力がすべて変わる）。Gallery の Vnig_TUMSAT-TG-2018 を generator から作り直す。

**FE-01 ラベル文字と表示の編集が黙って消える**（[W1a](remediation/W1a_multi_result_editor.md)）
- 分類: IEA（OIPC-C05/C06、web CLAUDE.md の live-edit）。
- 修正:
  - `app/feature-editor/label-actions.js` の `syncLabelEditor` は、`sourceReplaced` が真のときだけ消す。述語は既存のもの（`app/run-analysis.js:1952-1957`）を `bindingOptions` で渡す。
  - `hasSourceBoundEditorIntent` に label の map を加える。
- 削除: `labelOverrideContextKey` の全経路（state、history-snapshot、config、session-authority、reset、app-setup）。writer から field を外し、reader は既存の `copyFields` で無視する。
- テスト:
  - Node: feature 集合が重ならない SVG を順に同期しても、3 つの map が変わらないこと。
  - Playwright: batch の Result 往復 → Generate → Save で編集が残ること。非表示 → Generate → 表示 → Generate で label が戻ること。
- PR: 単独。8 files、net はマイナス。Ordinary。Review REQUIRED（compatibility-path）。

**TR-01（と N-04）Multi-Record Canvas の凡例が custom track slot を無視する**（W6）
- 分類: IEA（`docs/CLI_Reference.md:758, 781, 1348-1350`）。
- 修正: `gbdraw/diagrams/circular/assemble.py:2729-2788` の凡例の組み立てを `build_circular_legend_table(records, ...)` に抜き出す。単一 record と Multi-Record Canvas の両方がこれを呼ぶ。
- 削除: `gbdraw/api/diagram.py:3101-3158` の独自の組み立て（conservation の重複を含む）。
- テスト: 1 record で Multi-Record Canvas と単一 record の凡例が一致すること（slot の色、legend_label、AT skew、annotation、depth、conservation）。
- PR: PV-08 の後。Python +40/−70 行。Review REQUIRED（凡例の owner の移動）。

**CO-02 visibility を変えても古い LOSATP 結果を再利用する**（[W5](remediation/W5_comparisons.md)）
- 分類: IEA（OIPC-C04、PD-OI-022）。
- 修正: `app/run-analysis.js:468-541` の `canReuseResolvedProteinArtifacts` に `featureVisibility` を加え、committed Session の `featureVisibilityTableFile` と正規化してから比べる。比較の helper は `services/session-request.js` が export する。
- テスト: Similarity、Pairwise、Collinear の各モードで、visibility を変えて Generate した結果が、cache を消した新規実行と一致すること。
- PR: 2 files、約 25 行。Review REQUIRED（科学的出力）。

**CO-03 回転と Orient forward で LOSATN の ribbon が誤る**（W5）
- 分類: IEA（PD-OI-032、OIC-021 AC-14/15）。
- 修正:
  - (a) Linear の向きの owner を `seq.region_reverse` だけにする。rotate も alignment と同じく `region_reverse` に書き、override を null にする。`resolveLinearRecordReverse` を 1 つ export し、`app/run-analysis.js:492, 3177` と `app/match-sequences.js:868` はこれを使う。
  - (b) Q-FRAME が決まるまでの暫定策: 向きが変わるとき、生成済みの nucleotide 行を投影し直す。アップロードした比較を持つ record では、Orient forward を理由付きで不可にする。
  - Q-FRAME が C に決まれば、(b) は削除する。
- テスト: 4 経路（カードの逆相補、popup の Orient forward、その後の新規 Generate、Save→Load→Generate）で ribbon の端点が一致し、LOSAT の job が 0 であること。
- PR: 4〜5 files、約 150 行。Review REQUIRED。main には override がないので、互換 reader は不要。

**CO-04 ファイルのまとめ方でリンクが変わる**（W5）
- 分類: 決定済（02 の D-19 = B。PD-OI-018 revision 4）。
- 修正（B）:
  - `app/linear-sources.js:151-212` の source ファイル単位の batch は維持する。
  - 要求していない自己検索（record 自身への検索）は、どのモードでも実行しない。現在は、Collinear の inference が OFF のときだけ除いている。
  - 同じ source の中の record どうしの比較は、query の record を除いた database で検索する。
  - 純関数 `planLosatSourceJobs` を作り、実行と `app/losat-settings.js:48-90` の見積もりの両方で使う（N-10）。
  - E-value の database の範囲を、docs と Run Info に書く。
  - CLI は変えない（D-40）。Web と E-value が違うことを docs に書く。
- テスト:
  - 1 ファイルに入れた 2 record で、Max target seqs = 1 のとき、別ファイルと同じリンク数（19）になること。
  - 自己検索の job がないこと。
  - job 数の見積もりと実行が一致すること。
  - 既存の `tests/web/linear-sources.test.mjs:190-238` を新しい規定に合わせて更新する。
- PR: 3〜4 files。

**CO-05（と N-03）13 列以上の outfmt 6 を誤読する**（W6）
- 分類: 誤読と読み飛ばしをなくすのは IEA（`docs/REFERENCE/input-formats-and-tsv-schemas.md:28`）。12 列を超える表を受け付けるかは PDR（DP-4）。
- 修正:
  - `gbdraw/io/comparisons.py` に `read_comparison_table` と `normalize_comparison_dataframe` を置く。後者は `conservation.py` の `_coerce_comparison_dataframe` を移したもの。
  - 4 か所の reader（`io/comparisons.py:73-78`、`session_request_codec.py:3330-3335`、`api/record_planning.py:1426-1431`、`analysis/conservation.py:175-180`）をこれに置き換える。
  - `load_comparisons` の読み飛ばしを ValidationError にする。
- テスト: 12/13/14 列、11 列、途中の行だけ列が多い・少ない、数値でない値の各入力を、4 つの入口で確かめる。存在しないファイルで CLI が非ゼロで終了すること。
- PR: Python +60/−45 行。Review REQUIRED（科学的出力、CLI の挙動の変更を CHANGELOG に書く）。

### 7.2 P2

#### 入力と Session

**IN-04 cropped Linear ファイルの置換で古い crop が残る**（W2）
- 分類: 複数 record への置換は IEA（PD-OI-018）。単一 → 単一は DP-2。
- 修正: `app/app-setup.js:1044-1062` の `expandDiscoveredLinearRecords` は、`region_record_id` があるときだけ早期に return する。展開する行は source 単位の値だけで作る。
- テスト: crop と定義を付けてから 2 record のファイルに置き換えると、2 行に展開され、region と定義が空になること。
- PR: 1 file、約 12 行。IN-06 と同梱。

**IN-05 モード往復で Multi-Record Canvas の並び順が消える**（W2）
- 分類: IEA（OIPC-C06）。
- 修正: `runCircularRecordRefresh` は mode が circular でなければ return する。削除するのは入力がなくなったときだけにする。
- テスト: Node の `refreshCircularRecordOrder` の harness で、並び順が往復後も残ること。
- PR: 1 file、net −2。

**IN-06 1 ファイル内の全 record が 1 番目の organism になる**（W2）
- 分類: 決定済（02 の D-12 = B）。
- 修正（B）:
  - `app/app-setup.js:1064-1081` で、`records[0].inferredDefinition` を file の既定値に入れる処理をやめる。
  - 定義の優先順位を「record に入力した値 → file に入力した値 → その record 自身の推定値」にする。対象は `app/linear-sources.js:107-145` の実効値の計算。
  - record ごとの推定値は、discovery の結果（`records[i].inferredDefinition`）から取る。
  - Session に推定値を保存する必要がある場合は、current writer の項目として足す。旧 Session は読み込み時に推定し直す。
- テスト: organism が異なる 3 record のファイルで、各 record の定義がそれぞれの organism になること。file に入力した定義は全 record に適用されること。Reset で推定値に戻ること。Session の往復で変わらないこと。
- PR: IN-04 と同梱。

**IN-08 最初の Generate の前の Save が UNKNOWN になる**（W2）
- 分類: IEA（拒否は仕様。案内は PD-OI-046）。
- 修正: active mode の入力の有無を確かめる関数を `services/config.js` に 1 つ置き、Save と Generate の両方から呼ぶ（`app/run-analysis.js:2692-2696` を置き換える）。
- PR: 2 files。X-01 の PR-1 の後。

**SE-05 旧形式の Session を読むと Save が失敗する**（[W3](remediation/W3_validation_errors.md)）
- 分類: 文言は IEA。挙動は PDR（DP-4、推奨 A）。
- 修正（A の場合）: `services/config.js:4019, 4030` を `SESSION_SAVE_REQUIRES_GENERATE`（Generate ボタン付き）にする。旧形式の Session を読んだときに notice を出す。`docs/SESSION_COMPATIBILITY.md` に明記する。
- テスト: v39 と schema-v2 の fixture で Save → 案内とボタンが出ること。Generate → Save が成功すること。

**SE-06 CLI の Linear+BLAST Session で Generate できない**（[W4](remediation/W4_history_session.md)）
- 分類: Inherit の不具合は IEA（PD-OI-008 の「実行可能だが投影できない比較は、明示的な inheritance で使い続けられる」）。読み取り専用の扱いは正しい。
- 修正:
  - `services/session-request.js:3785-3799` の CLI sidecar 用の分岐で、Linear の binding uid を、確定済みの request のファイル単位の recordKey に置き換える。
  - `app/run-analysis.js:4590-4593` の Inherit は、`promoteCanonicalRenderRequestToCurrent` を通してから渡す（N-17）。
  - 読み込み側を直す。CLI の書き出し側（`cli-seq-N`）は変えない。0.12.0 と main に存在する形式で、main の CLI が書く v42 は Web で読めるため。
- テスト: `tests/web/session-cli-compatibility.playwright.spec.js` に `linear blast` を追加する。main で作った v42 の sidecar を来歴付きで fixture に固定する。
  - PRESERVED_READ_ONLY のままであること。
  - Inherit → Generate が ok で、comparisons と svgSemantics が確定済みのものと一致すること。
- PR: 2 files、約 30 行。Review REQUIRED（互換経路）。prototype で、v42 と v44 の両方で Generate できることを確かめた。SE-07 と同じ分岐を触るので順に進める。

**SE-07 CLI Session で凡例の位置が失われる**（W4）
- 分類: IEA（`docs/REFERENCE/session-and-request-compatibility.md:50-52`、PD-OI-019、OIPC-C05）。
- 修正:
  - projection が、確定済みの (mode, grouping) の `layoutPreferences` を返す。作成には既存の `createDefaultLayoutPreferences` と `updateActiveLayoutPreference` を使う。
  - アクセサ経由の `legend` と `plot_title_position` を projection から削除する。
  - `restoreLayoutPreferences` は、保存された値がなければ projection の値を使う。`preserveActive` の分岐は削除する。
  - merge の順序で解決する細工は入れない。
- テスト: CLI 往復に `--legend upper_left`、`--multi_record_canvas --legend upper_right`、Linear の `--legend left` を追加する。読み込み後と最初の Generate の後で、凡例の位置が変わらないこと。
- PR: 2 files、約 60 行。Review REQUIRED。Gallery の publication への影響を確かめる。

**SE-08 Qualifier Priority の編集が無視される**（W4）
- 分類: (a) は IEA（PD-OI-009、OIPC-C03）。(b) の読み込みの遅さは DP-4。
- 修正（Python のみ）:
  - `gbdraw/web_support/config_overrides.py` で、派生した key を比較と raw の対象から除く。派生 key の集合は `gbdraw/labels/filtering.py` に 1 か所で定義する。
  - 防御として、`gbdraw/api/diagram.py:283-291` で table を付けるときに古い map を捨てる。
- テスト: tobacco と通常の CLI Linear の config で、保持する設定が空になること。古い `priority_map` がある config に table を付けると、table のとおりのラベルになること。
- PR: Python のみ（Web の規模には入らない）。main で保存された Session のための互換 fixture を置く。

#### History

**SE-01 Undo/Redo 後の catalog が admit されない**（W4）
- 分類: IEA（`docs/REFERENCE/web-app.md:706-707`）。
- 修正（R11）:
  - `services/history-snapshot.js` の `buildArtifactCheckpoint` は、catalog を JSON から除き、参照を checkpoint を key とする WeakMap に持つ。
  - `applyArtifactCheckpoint` で参照を戻す。
  - `services/config.js` の `applyEditorStateData` と `normalizeEditorStateData` は catalog を複製せず、admit 済みのものだけを state に入れる。
- 削除: `adoptCatalog` option と複製の分岐（`:1169-1171`）、使われていない `applyGeneratedArtifactSnapshot`。
- テスト: checkpoint → Undo → Redo の後も同じ catalog object であること。Linear → Circular で pageerror がなく 37 features であること。
- PR: 2 files、約 60 行、net は負。import rollback の経路（N-19）と checkpoint のサイズ（N-20）も同時に直る。

**SE-02、SE-03、SE-04 入力の History の境界**（W4）
- 分類: IEA。
- 修正（R11）:
  - `app/history-inputs.js`: checkbox と radio は capture-phase の `change` で transaction を始め、pointerdown の分岐から外す。button は click の capture だけで始める。
  - `services/history.js`: `begin(label, { source, owner })` に owner を持たせる。重複する確定処理 3 か所を `settlePendingIntent()` にまとめる。
  - `app/history-shortcuts.js:7`: select を text 系から外す。
  - `onKeyDown`: 修飾キーを押したままのキーは無視する。
  - `undo()` と `redo()`: 先頭で `settlePendingIntent()` を呼ぶ。
  - drag の 2 ファイル（N-18）でも owner を渡す。
- テスト: 入力イベントの組み合わせテスト（`tests/web/history-inputs.playwright.spec.js`）。text 欄では Ctrl+Z がブラウザ標準の undo のままであること。
- PR: 5 files、churn 約 100。Ordinary。Undo の step が増えるのは意図した変化として明記する。prototype で確かめた。

**GE-06 Generate 中の Ctrl+Z**（W4）
- 分類: PDR（DP-5、推奨 A）。現状は History の不変条件に反するので維持できない。
- 修正（A の場合）: `services/history.js` に `activeReplacement` を足す。`historyAvailability()` を `undo`、`redo`、`canUndo`、`canRedo` で共有する。開始時と解除時に `touchTransaction()` を呼び、判定を reactive にする。
- テスト: `beforeDiagramGenerationResponse` の hook で Generate を止めて Ctrl+Z/Ctrl+Y を押す。止めている間は確定済みの request と履歴の数が変わらないこと。既存の `tests/web/session-operation-consistency.test.mjs:28-33` を書き換える。
- PR: 約 20 行。SE-02〜04 の PR にまとめてよい。

**FE-04 Exact product での非表示が Redo のあと再び表示される**（W1a）
- 分類: IEA（History は owner の action と同じ操作を記録する）。
- 修正: `app/feature-visibility.js` の `resolveEffectiveFeatureVisibility` が、editor の exact-qualifier rule を Python と同じく大小無視で評価するようにする。action と reconcile は同じ matcher を使う（`getMatchingQualifierFeatures` を置き換える）。
- テスト: product rule のあと `reconcileFeatureVisibility()` が `off` を出すこと。Undo → Redo のあとも hidden のままであること。
- PR: 2 files、約 60 行。FE-10 と同梱。

#### Generate とオプション

**GE-02 stroke の即時反映が不正な値を Result に書く**（W1b）
- 分類: 決定済（C）。
- 修正: `app/svg-styles.js:409-548`（`applyStylesToSvg`）、`:650-667`（watcher）、`:681`（export）を削除する。`index.html:4329, 4458` の表示を「Applies on Generate」にする。
- 残すもの: 個々の feature の stroke を Auto に戻すための `originalSvgStroke`（`app/feature-editor/canvas-actions.js:87-97`）は残す。
- テスト: Generate 後に stroke を 5 → 空 → −1 と変えても、Result がバイト一致すること。`tests/web/history-generated-authority.playwright.spec.js:180-200` を書き換える。
- PR: 2 files、net 約 −140。

**GE-03 Run Info の Source recipe で図を再現できない**（W3）
- 分類: IEA（`docs/REFERENCE/web-app.md:814-818`。Circular は 0.13.0 の挙動の復元）。
- 修正:
  - Circular: `circular_definition_interval_for_font` を Python に 1 つ置き、`gbdraw/circular.py:951-954` と codec の複製を置き換える。codec の移行で、font があって interval がないときだけ interval を補う。
  - Linear: scale の font の override があって ruler の override がなく、ruler が描かれるときは、recipe を理由付きで unavailable にする。
- PR: Python 約 30 行と `app/run-info.js` 1 file。TR-07 と同梱。

**X-01 エラー変換の抜け**（W3）
- 分類: IEA（PD-OI-046）。
- 修正（R6）:
  - producer の契約（`diagnosticError`、`diagnostic=`）を定める。
  - locator を summary に出す（`decorationResult` の特例は削除する）。
  - `FIELD_LABELS` と REASONS に data の行を足す。
  - session import の Worker が送っている code を捨てずに使う。
  - wrapper が原因の code を捨てないようにする。
  - Generate ボタンを足す。
- 削除: 移行した producer に対応する `NATIVE_VALIDATIONS` の行と正規表現、`services/config.js:4238-4241` の書き換え、`REASONS.READ`、`_web_error_field`、producer がない dead 行。
- PR:
  - PR-1（8 files、net +60〜90）
  - PR-2（annotation、net はマイナス）
  - PR-3（Python）
  - PR-4 以降は ratchet の系列で、最後の PR で `nativeValidation` を削除する。

**X-02 不正な数値が黙って既定値になる**（W3）
- 分類: IEA（OIPC-C01、C03）。
- 修正（R7）:
  - PR-A（Python が先）: `gbdraw/api/options.py` の typed options で window、step、depth window/step、dinucleotide を検証する（dinucleotide は ACGTU の 2 文字で大小を区別せず、U は T として扱う。02 の D-26）。閾値の定義域を `gbdraw/mode_profiles.py` の表にし、生成物に載せる。slot の `nt` も同じ validator を使う。
  - PR-B（Web）: `optionalNumber` と `optionalPositiveInteger` を削除し、`projectOptionalNumber` に置き換える。`normalizeBlastThreshold*` と draft への書き戻し（`app/run-analysis.js:885-895, 2715-2731, 2996-3010`）を削除し、`resolveComparisonThresholds` を 1 回だけ呼ぶ。
  - PR-C: 同じ族の残り（multi-record ratio、plot title、label spacing など）。
- 順序: PR-A を先に入れる。逆順だと、Web が 0 を送り、Python がそれを受理して空の track を描く。
- テスト: 共有ベクタで、CLI、typed request、Web が同じ field と reason で拒否すること。Generate の後に `adv.*` が変わらないこと。

#### 比較

**CO-01 Threaded の既定値が、分離されていない環境で UNKNOWN になる**（W5）
- 分類: IEA（PD-OI-018 と `docs/REFERENCE/web-app.md:19-22`）。threaded は厳格なままにする。
- 修正:
  - `services/losat.js:101-109` の判定を `losatThreadingPrecondition()` に切り出し、`THREADED_LOSAT_UNAVAILABLE` を投げる。判定は、キャッシュされていない job を dispatch する時点で行う。
  - LOSAT の実行中の失敗を request-validation ではなく LOSAT 用の stage で報告する（`app/run-analysis.js:1903`）。
  - 設定欄に「この環境では不可」と表示する。
- PR: 2〜3 files、約 40 行。

**CO-06 アップロードした表の ID を照合しない**（W5）
- 分類: 照合するかは PDR（DP-3、推奨 B）。metadata の矛盾は IEA。
- 修正（B の場合）: Python の planner（`project_source_bound_comparisons` の隣）で照合する。`data-*-record-id` は常に端点の record から取る。
- PR: Python 2 files、約 80 行。CO-05 の後。

**CO-07 逆相補の Save Raw が CLI で再現しない**（W5）
- 分類: PDR（DP-3 の Q-FRAME）。
- 修正（C の場合）: planner（`gbdraw/api/record_planning.py:216-300`）が、feature binding を持たない表を F から投影し、1..L を検証する（N-07）。JS の F→V 変換と `convert_losat_nucleotide_to_display_tsv` を削除する。main の Session の nucleotide 行は、reader で一度だけ V→F に変換する。
- PR: Architecture profile、約 350 行、net ≈0。Review REQUIRED（科学的出力と互換）。

**CO-08 Similarity group の名前が別の group へ移る**（W5）
- 分類: 移る点は IEA。消えた group の扱いは PDR（DP-3、推奨 A）。
- 修正（A の場合）: `app/run-analysis.js:1452-1461` の `pruneOrthogroupOverrides` を、member 集合の完全一致で付け替える `rekeyOrthogroupOverrides` に置き換える。
- テスト: 現在の誤りを正しいとしている既存テスト `tests/web/orthogroup-computation-cache.test.mjs:43-45` を書き換える（OIPC-C08）。
- PR: 3〜4 files、約 100 行。Session に項目を 1 つ足す（absent は空）。

#### Feature 編集

**FE-02（+PV-09）batch で範囲指定の編集が表示中の Result にしか反映されない**（W1a）
- 分類: PDR（DP-1、推奨 B）。
- 修正（B の場合）: History after-apply の処理を `reconcileMountedEditorIntent` にまとめ、Result を mount するときにも呼ぶ。凡例は、既存の candidate plan の legend 操作を mounted root に idempotent に適用する。
- PR: 5〜7 files。Architecture profile。FE-01、FE-03、FE-04 の後。

**FE-03 batch で drawer が別の record を並べ、Edit が効かない**（W1a）
- 分類: IEA（one canonical selected value、同じ availability predicate）。
- 修正: `state.js` の `filteredFeatures` を、mounted Result の `renderedFeatureIdentities` に限る。record picker は、mounted Result に record が 2 つ以上あるときだけ出す。
- PR: 3 files。FE-02 より前。

**FE-05 Circular 複数 record で「Selected features」が Generate を止める**（W2）
- 分類: IEA（Python の typed request、PD-OI-041）。
- 修正: `services/session-request.js:1017-1036` の Circular の record 集合の解決を純関数として export し、annotation の catalog はそれに依存する。
- 削除: `allowExplicitSelectors` の概念（`app/annotations/record-catalog.js:128, 145, 152`、`app/annotations/validation.js:53-55`、`app/annotations/record-selector.js` の該当部分）。
- PR: 6 files、net 約 −10。

**FE-06 既定の All 検索で遺伝子名が配列に一致する**（[W7](remediation/W7_tracks_search_shell.md)）
- 分類: 決定済（A）。
- 修正: `app/feature-search/search-core.js:392-397` で、All から配列の内容と `/translation` を外す。`services/standalone-interactivity-assets.js:2661-2664` も同じ PR で変える。
- テスト:
  - All で `CYTB`、`dnaA`、`gyrA` が名前の一致だけを返すこと。専用 field の配列検索は従来どおりであること。
  - 埋め込み版と module 版で、一致する id の集合が等しいこと（既存の `embeddedFunctionSource` を使う）。
- PR: 2 files。Gallery の `gallery/examples/*.svg` を generator で作り直す。

**FE-07（と N-05）/translation のない CDS で開始コドンが M にならない**（W6）
- 分類: IEA（決定的な科学の規則）。LOSATP と orthogroup の出力が変わるので、preflight を記録し、Review REQUIRED とする。
- 修正:
  - 翻訳の関数を `gbdraw/core/sequence.py` に 1 つ置く。5′ 端が完全で、読み枠が 1 で、table の start codon なら M にする。GFF3 の phase を読み枠に使う。
  - `gbdraw/web_support/feature_metadata.py:262-295` と `gbdraw/analysis/protein_colinearity.py:2761-2775` をこれに置き換える。
- テスト: GTG、TTG で M になること。部分的な CDS はそのまま訳すこと。GFF3 の phase=1。MG1655 の `/translation` と一致すること。
- PR: Python +45/−25 行。protein cache の昇格経路（`app/run-analysis.js:3516`）が古い結果を昇格しないか確かめる（ER）。

#### 凡例・Preview・出力

**PV-01** — 欠陥ではない（PD-OI-052 の受容済みリスク）。DP-1 で A（現状維持）を推奨する。

**PV-02 feature のない凡例項目の名前変更が Generate で戻る**（W1b）
- 分類: IEA（`docs/REFERENCE/web-app.md:101-102`、`legendRenames` の契約）。
- 修正: `app/feature-editor/color-actions.js` の `renameLegendEntryInSvg` は、renderer が生成した行では caption だけを更新し、`originalCaption` と original 系の値は変えない。
- テスト: rename の後、`compileDirectEditorMutationPlan` が `legendRenames` を出すこと。
- PR: PV-04、PV-12 と同梱（2 files、約 60 行）。横書きの凡例で sort、rename、Generate を組み合わせた Playwright を必須にする。

**PV-03** — 欠陥ではない（PD-OI-052）。機能要望として DP-1（推奨 B）。

**PV-04 既存名への rename が UNKNOWN になる**（W1b）
- 分類: 決定済（A）。
- 修正: 衝突判定（`app/feature-editor/color-actions.js:787-798`）を、feature/rule の分岐（`:780-785`）より前に移す。文言は X-01 の `LEGEND_NAME_CONFLICT`。
- 別件: CLI の凡例衝突（N-06）は Python の owner で別に扱う。

**PV-05 Web の PDF が px を pt として扱う**（W6）
- 分類: IEA（CSS の単位、同じファイルの PNG の `dpi/96`）。
- 修正: `services/export.js:382-396` で 0.75 を掛ける。定数は PNG の処理と共有する。
- PR: PV-06 と同梱。CHANGELOG に書く。

#### トラック

**TR-02 削除・無効化・移動した depth 行が戻る**（W7）
- 分類: 中核は IEA。Linear の追加の挙動は決定済（D-TR02 = A）。
- 修正:
  - 純関数 `managedDepthSlotAdditions` を `app/depth-track-state.js` に置く。除去は既存の `dropInvalidManagedDepthSlots` を使う。
  - 各 editor に commit wrapper（`reconcileCircularDepthSlots`、`reconcileLinearDepthSlots`）を置き、depth 入力を変える遷移（`setCircularDepthFile`、`removeCircularDepthTrack`、`addLinearDepthTrack`、`removeLinearDepthTrack`、`setLinearDepthFiles`）から明示的に呼ぶ。
- 削除: watcher `app/app-setup.js:2185-2206`、ensure の claim・再有効化・付け替え（`app/circular-track-slots.js:1733-1766`、`app/linear-track-slots.js:1349-1377`）。
- テスト: 無効な行は reconcile しても変わらないこと。削除した行は戻らないこと。同じ表を Linear にも流すこと。既存の `tests/web/depth-track-session.playwright.spec.js:469-553` を書き換える。
- PR: 4 files、churn 200〜300、net ≤0。`architecture-change`。TR-03 と同梱。

**TR-03 アップローダーの Remove で depth 行が残る**（W7）
- 分類: IEA。
- 修正: TR-02 の reconcile を `setCircularDepthFile` から無条件に呼ぶ。
- 残る課題: Linear 側は PD-OI-025 により結論が異なり、ER。

### 7.3 P3

| ID | 分類 | 修正 | 付録 |
|---|---|---|---|
| FE-08 | IEA（Python は IGNORECASE。`docs/REFERENCE/web-app.md:735`） | `app/feature-selector.js` の一意性キーの value を大小無視にする。出力する値と正規表現は変えない | W6 |
| FE-09 | PDR（DP-2、推奨 A） | `app/feature-editor/rule-actions.js:447-470` の instance id の分岐を削除し、常に stable hash を出す | W2 |
| FE-10 | IEA | `resetColorDialog.defaultColor` を削除し、確定の時点で palette から得る。Cancel は owner の遷移を呼ぶ。キャンセルしたときに限らず、dialog を一度使った後にも起きる | W1a |
| FE-11（と N-11） | IEA（INSDC の 1 始まり、Python の `location_parts[].display`） | `buildFeatureLocation` を export して共有する。検索の Start/End の項目は削除する。standalone は別 PR で直し、Gallery を作り直す | W2 |
| FE-12 | IEA（P2 候補） | JS の `normalizeSpecificRuleColor` と Python の行ごとの検証を、共有ベクタで拘束する | W3 |
| GE-07 | IEA | `app/legend-layout/reposition-actions.js:116-119` の throw を `return false` にする | W1b |
| GE-09 | IEA（lazy Worker の規則「後の操作は Worker を再利用する」） | `services/diagram-generation.js:605-616` で、`hadActiveRequest` が偽なら Worker を終了させない（2 行）。Cancel の後に要求が送られないことは `app/run-analysis.js:2012-2016` が保証している | W4 |
| PV-06 | IEA | UTF-16 の index で繰り返し、空白も出力して `xml:space="preserve"` を付ける | W6 |
| PV-07 | 欠陥ではない（PD-OI-052）。DP-1（推奨 B） | B の場合、Session の padding を候補に適用し、`app/app-setup.js:2298-2301` の初期化を削除する | W1b |
| PV-10 | ER | live と Generate の凡例 bbox を計測し、差が reflow にあれば oracle にケースを足して JS を合わせる。補正係数は足さない | W1b |
| PV-11 | IEA（PD-OI-038） | `closeRightDrawer` で、focus が drawer の中にあれば Editor toggle に戻す | W1b |
| PV-12 | IEA（OIC-024） | `index.html:6512` の文言を実際の操作に合わせる | W1b |
| CO-10 | FASTA ヘッダは IEA、popup は PDR（DP-3） | renderer が source への写像を属性として出し、ヘッダと popup で 1 つの関数を使う | W5 |
| SE-04 | IEA | 第 7.2 節の SE-02、SE-03、SE-04 と同じ PR で直す | W4 |
| SE-10 | PDR（DP-2、推奨 A） | A の場合、Reset の checkpoint 内で Linear の record ごとの表示状態を初期化し、alignment は `clearCommittedPlan` で消す | W2 |
| TR-06 | IEA | suppress の印を form から計算する双方向の純関数にする。`disableCircularTrackSlotsForSuppress` と `restoreCircularTrackSlotsForSuppress` を削除する | W7 |
| TR-07 | IEA | `app/run-info.js` に `assertCliSlotTokenLossless` を置き、CLI の分割規則で読み戻して一致しなければ unavailable にする。GE-03 と同梱 | W3 |
| TR-08 | IEA | `findTrackSlotGeometry` を slotId での一致だけにし、根拠のない index fallback を削除する | W7 |
| TR-09 | preset Reset は IEA。初回は欠陥ではない | `resetCircularTrackSlotsFromSimpleControls` を preset 版への委譲にし、`showTicks` を渡す | W7 |
| TR-10 | (a) IEA、(b) 決定済（A） | (a) `index.html` に可視ラベルと同じ `aria-label` を付ける（証拠では 26 件）。(b) help tip を disclosure button にする。churn 500〜800 なのでパネル単位で分ける | W7 |
| TR-11 | IEA | 行を見出し（3 列）と操作の行に分ける CSS を足し、両モードで使う | W7 |
| TR-12 | IEA | リンクを `https://github.com/satoshikawato/gbdraw/blob/main/docs/REFERENCE/web-app.md#record-selection-and-layout` にする。`tests/test_web_packaging.py` に tutorial JSON の link checker を足す | W7 |

## 8. 進め方

### 8.1 フェーズ

**Phase 0（判断、guard、検査基盤。runtime の挙動は変えない）**

- 0-1 判断: DP-1〜DP-5 を一度に Owner に示す。決定済の 6 件は、第 5 節の未記入の項目を Owner が埋めてから PD-OI-056 以降として登録する。
- 0-2 guard だけの PR:
  - `tests/web/architecture-contracts.test.mjs` の厳密な件数を上限に変える（IN-01 の削除に必要）。
  - privileged capability の分割は、authority の PR を出してから checker の PR を出す（checker と authority の分離規則）。
- 0-3 検査基盤の PR（tests-only）:
  - 既知の不具合に印を付ける（Playwright の `test.fail(true, 'IN-01')`、pytest の `xfail(strict=True)`、ベクタの `knownDefect`）。修正 PR が印を外す。
  - 共有ベクタ（record 検出、色の定義域、option の定義域、比較表の列数）。
  - helper での自動検査: `generateAndWaitForResult` で UNKNOWN、pageerror、Run Info の `NaN` がないことを確かめる。15 の spec がそのまま検査網になる。
  - 2 record の batch fixture。
- 0-4 文書のずれを直す（N-15）。

**Phase 1（P1。科学的出力と CLI も正しくなる順）**

1. PV-08 と N-01（主な原因と文言）
2. TR-01 と N-04
3. CO-05 と N-03（A で先に出す）
4. FE-07 と N-05
5. IN-02、IN-03、N-02
6. FE-01
7. IN-01（決定済。0-2 の後）
8. CO-02
9. CO-03(a)
10. CO-04（DP-3 の後）

参照 SVG が変わる PR は、`--update-reference-outputs` で作り直して Review REQUIRED にする。Gallery に影響するものは generator で作り直す。

**Phase 2（P2）**

| まとまり | 項目 |
|---|---|
| 基盤 | X-01 PR-1 → X-02 PR-A（Python）→ X-02 PR-B（Web） |
| History | GE-09 → SE-01 → SE-02〜04（GE-06 は DP-5 の後に同じ PR へ入れてよい）→ FE-04 と FE-10 |
| Feature 編集 | FE-03 → FE-02（DP-1 の後）、FE-05、FE-06（決定済） |
| 即時編集 | GE-02（決定済）、PV-02 と PV-04 と PV-12 |
| 入力 | IN-04 と IN-06、IN-05、IN-08、SE-05 |
| Session と CLI | SE-06 → SE-07（同じ分岐なので直列）、SE-08(a)（Python）、GE-03 と TR-07 |
| 比較 | CO-01、CO-06、CO-07（DP-3 の後）、CO-08 |
| 出力 | PV-05 と PV-06 |
| トラック | TR-02 と TR-03（決定済） |

**Phase 3（P3）** owner ファイル単位でまとめる。

- TR-06、TR-08、TR-09
- TR-10(a) と TR-11
- TR-10(b)（パネル単位）
- TR-12
- FE-08、FE-09、FE-11
- GE-07（PV-10 の証拠の後）
- PV-11、CO-10（SE-04 と GE-09 は Phase 2 の History の PR に含める）

**Phase 4（仕上げ）**

- "Result content commit" の allowlist を最終形まで狭める。
- 行列テスト（G-A、G-B）を、条件を満たせば Product Impact map の contract に昇格させる。
- 定期監査の手順を定着させる（第 8.3 節）。

### 8.2 再発防止のガード（[W8](remediation/W8_verification_workflow.md)）

| ID | ガード | 段階 | 検出できた項目 |
|---|---|---|---|
| G-A | 即時編集 ≡ Generate の metamorphic テスト。fixture は crop、逆相補、record ラベル、grid を含む | dev staging | IN-01、GE-02、PV-01〜03、PV-07、PV-10、CO-08、SE-07 |
| G-B | Result topology の行列（single/grid/batch × 編集の種類） | dev staging | FE-01〜03、FE-05、PV-09 |
| G-C | 編集ではない操作の前後で、利用者が持つ状態が変わらないこと | dev staging（Node の部分は PR） | FE-01、IN-05、TR-02、TR-06、TR-09、PV-07、SE-07、FE-04 |
| G-D | JS/Python の共有ベクタ | PR | IN-02、IN-03、IN-06、FE-08、FE-12、CO-05 |
| G-E | Source recipe と CLI の一致（代表 10 probe。57 probe の完全版は tool として残す） | dev staging | GE-03、TR-07、SE-06、SE-07、CO-07 |
| G-F | 比較の metamorphic テスト | PR と dev staging | CO-02〜04、CO-06、CO-07、CO-10 |
| G-G | エラーと数値の網羅（producer の網羅、helper での自動検査、数値の表） | PR | X-01、X-02、SE-05、PV-04、IN-04、IN-08、GE-07 |
| G-H | History の境界の行列（入力方法 × 処理中） | dev staging | SE-01〜04、GE-06、FE-04 |
| G-I | Python の性質テスト（1 record の grid ≡ single、Web 既定 × 長い organism、GTG の翻訳、数値の制約） | PR | TR-01、PV-08、FE-07 |
| G-J | 既存の仕組みを使う静的ガード（allowlist の縮小、reader が 1 つ、Gallery のリンク、a11y） | PR と dev staging | クラス 1/2、CO-05、TR-12、TR-10 |

- **予算:** PR smoke は 19/19 で埋まっているので追加しない。PR 段階に足すのは Node と pytest だけで、合計 15 秒未満を目安にする。
- **現在の誤りを正しいとしている既存テスト:** 修正と同時に書き換える（OIPC-C08）。
  - `tests/web/composition-layout.test.mjs:796`（PV-01。PD-OI-052 の範囲なら維持）
  - `tests/web/session-operation-consistency.test.mjs:28-33`（GE-06）
  - `tests/web/orthogroup-computation-cache.test.mjs:43-45`（CO-08）

### 8.3 ワークフロー

- **PR の規模:** Web の Ordinary profile（8 files / churn 800 / 純増 100）に収める。owner の移動や旧経路の削除は `architecture-change` にする。
- **並べ方:** ホットなファイルを触る PR は直列にする。dev の staging が緑になってから、同じクラスの次の修正をマージする。1 回の staging の周期でマージする runtime の PR は 2〜3 本までにする。
- **自動マージ:** Owner の許可（2026-09-29）に従い、dev への自動マージはしてよい。ただし、次の PR は Owner のレビューを待つことを推奨する。Review REQUIRED は Gate を失敗させないため。
  - Product Contract の co-change
  - 参照出力の変更
  - size の超過
  - `architecture-change`
  - policy の文書
- **推奨案で進める運用:** 残りの判断に「推奨案で進める」の standing instruction を使う場合は、receipt に Owner-delegated と明記する。AGENTS.md の「製品の選択肢を自律的に選ばない」と緊張するので、Pack を一度提示して確認する運用を推奨する。
- **定期監査:** 新しい gate や workflow は作らず、dev → main の昇格 PR のチェック項目にする。
  - 手順: 前回の昇格以降に変わった領域を、時間を区切って探索的に監査する。
  - 出力: この README と同じ形式の文書。確認した不具合のクラスごとに「ガードを足す」か「足さない判断」を記録する。
  - 終了条件: P1 がすべて直っているか、Owner が明示的に waiver を出していること。

## 9. 未検証と残るリスク

- 設計担当のエージェントの確認の範囲はそれぞれ異なる。
  - W4 は History の 7 項目を、複製した JS への試作 patch で確かめた。
  - それ以外の多くは、コードを読んだ結果と監査の証拠による。
  - 各付録に「確認」と「推定」の区別がある。実装時には、各 PR の最初に監査の spec を assert 付きに直し、dev で失敗することを確かめてから修正する。
- SE-07 は仕組みの確認だけで、修正は試していない。SE-08(a) の Python の修正も試していない。
- ER の項目:
  - SE-08(b): CLI Session の読み込み時の Python 検証をどう扱うか。
  - PV-10: live と Generate の差の出所。
  - N-16: ラベルの reflow で draft が Result に漏れるかどうか。
  - FE-07: protein cache の昇格経路が古い結果を昇格しないか。
  - TR-03: Linear 側の扱い。
  - SE-05 の B: Save 時の catalog 回復が成り立つか。
  - 負の offset や spacing に意味があるか。
- 修正すると Web 既定の Circular 出力がすべて変わる項目（PV-08 と N-01）がある。Gallery、チュートリアルの画像、docs の capture の作り直しが必要かを、Phase 1 の前に確かめる。
- main（リリース版）にも含まれる不具合がある（IN-02 と IN-03 の起源、PV-08 の `021b4c49`、CO-05、TR-01）。backport するかは Owner の判断。
