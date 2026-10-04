# Decision Pack — Web GUI 監査の修正（2026-09-30）

Status: 確定（2026-09-30）。結果は末尾の「判断の記録」にある。
Owner: satoshikawato。基準: `origin/dev` `4c89bab1`。
関連文書: [01_REMEDIATION_PROPOSAL.md](01_REMEDIATION_PROPOSAL.md)（各項目の背景と修正の中身）。

この Pack の推奨は engineering の推奨で、Product authority ではない。
承認された receipt は、次のセッションが [OPTION_INTEGRITY_PRODUCT_CONTRACT.md](../OPTION_INTEGRITY_PRODUCT_CONTRACT.md) に PD-OI-056 以降として登録する。登録するとき、文面は翻訳も追加もせずに再現する。

## 回答の仕方

- 「推奨で承認」と答えると、その項目の receipt の全文（Choice、Rationale、Must preserve、May retire、Accepted residual risk）を、ここに書かれたとおりに承認したことになる。
- 別の選択肢を選んだ項目は、その案の receipt をこちらで書き起こし、もう一度確認をお願いする。
- 文面の一部だけを直したい場合は、直したい項目と文言を示す。
- A〜C の区分:
  - **A:** 決定済みの 6 件。receipt の残りの項目の承認だけをお願いする。
  - **B:** Product の判断。
  - **C:** 現状を維持する項目。新しい record は作らない。
  - **W:** 作業の進め方。Product authority ではない。

## A. 決定済みの 6 件（receipt の残りの項目の承認）

### D-01 FE-06 検索の All の範囲
- Concern: `web.feature-search.all-field-scope`

```text
PRODUCT_DECISION
Concern: web.feature-search.all-field-scope
Scenario revision: 1
Choice: A / ALL-EXCLUDES-SEQUENCE-CONTENT
Rationale: 既定の All で遺伝子名を検索したとき、名前の一致だけが返るようにする。配列への偶然の一致で結果が埋まらないようにする。
Must preserve: 専用の Nucleotide sequence と Amino acid sequence の field による配列検索（IUPAC の展開を含む）、Label・qualifier・Location など他の field、編集後のラベルの検索、Preview と Interactive SVG の検索結果の一致、件数の表示。
May retire: All が配列の内容と /translation の値に一致する動作。
Accepted residual risk: 配列の断片を All で検索していた利用者は、専用の field を選ぶ必要がある。作り直す前の Interactive SVG は旧挙動のまま。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-02 TR-10(b) help tip の到達性
- Concern: `web.help-tip.reachability`

```text
PRODUCT_DECISION
Concern: web.help-tip.reachability
Scenario revision: 1
Choice: A / FOCUSABLE-DISCLOSURE-TIPS
Rationale: キーボードとタッチの利用者が、hover と同じ説明に届くようにする。常時表示の説明を help tip に移した（PD-OI-054）ので、tip に届くことが必要になった。
Must preserve: 各 tip の文言、hover での表示、周囲の label と control の accessible name、Escape で閉じられること、390 px での操作、既存の id 付き tip の挙動。
May retire: id のない tip を hover 専用の aria-hidden の icon にする設計（js/components.js の意図のコメント）。
Accepted residual risk: tab stop が最大 175 個増え、キーボードでの移動が長くなる。Gallery と docs の capture を撮り直すことがある。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-03 D-TR02 managed depth 行の扱い
- Concern: `web.depth.managed-slot-lifecycle`

```text
PRODUCT_DECISION
Concern: web.depth.managed-slot-lifecycle
Scenario revision: 1
Choice: A / ADD-ON-FIRST-SOURCE-BOTH-MODES
Rationale: Circular と Linear で、Depth ファイルを設定したときの track 行の振る舞いを揃える。利用者が消した行や無効にした行を勝手に戻さない。
Must preserve: 明示の track slot が有効なときの authority、利用者が削除・無効化・移動した行、行の params と legend_label、Reset による再生成、Undo/Redo、Session の往復、PD-OI-025 の論理 series の範囲。論理 series が最初の source を得たとき、その index を参照する行（有効・無効を問わない）がなければ managed 行を 1 つ足し、series が source を失ったら managed 行を除く。
May retire: 無関係な切り替えのたびに Circular の watcher が depth 行を作り直す・再び有効にする・付け替える動作。Linear の "Add Depth TSV series" が無効な行を再び有効にする動作。
Accepted residual risk: 無効な stack に depth 行を持たない旧 Session は、読み込んで stack を有効にしても行が自動では足されない（Reset で作れる）。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-04 IN-01 Circular の定義系の設定の反映
- Concern: `web.circular.definition-application`

```text
PRODUCT_DECISION
Concern: web.circular.definition-application
Scenario revision: 1
Choice: A / CIRCULAR-DEFINITION-APPLIES-ON-GENERATE
Rationale: Result の定義行に、crop の長さ・GC%・record label と食い違う値が書き込まれないようにする。Linear と同じ「Applies on Generate」に揃える。
Must preserve: Species、Strain、Plot title、Title position、Title font、Default font size、Keep Full Definition の編集・保存・History、Generate 後の正しい定義（region の長さと GC%、record label と subtitle、逆相補、grid の順序）、Linear の現在の挙動、他の即時編集（色、ラベル、凡例など）。各設定には Applies on Generate の表示を付ける。
May retire: 上の Circular の設定の即時反映と、そのための helper（regenerate_definition_svgs）。
Accepted residual risk: これらを変えたとき、Generate するまで preview が変わらない。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-05 GE-02 全体の stroke 設定の反映
- Concern: `web.stroke.application`

```text
PRODUCT_DECISION
Concern: web.stroke.application
Scenario revision: 1
Choice: C / STROKE-APPLIES-ON-GENERATE
Rationale: 空欄・不正な値・Auto への戻しが、Generate と違う stroke として Result に残らないようにする。個々の feature の stroke 指定を全体の設定で上書きしないようにする。
Must preserve: 全体の stroke 設定（block、line、axis、scale の幅と色）の編集・保存・Generate での適用、個々の feature の stroke 編集の即時反映と Auto への復元、不正な値の検証。各設定には Applies on Generate の表示を付ける。
May retire: 全体の stroke 設定の即時反映。
Accepted residual risk: 全体の stroke を変えたとき、Generate するまで preview が変わらない。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-06 PV-04 凡例の名前の衝突
- Concern: `web.legend.rename-collision`

```text
PRODUCT_DECISION
Concern: web.legend.rename-collision
Scenario revision: 1
Choice: A / DIALOG-FOR-ALL-ENTRIES
Rationale: 凡例の名前を既存の名前に変えたとき、feature の有無にかかわらず同じ Merge / Suffix / Cancel の選択を示す。原因の分からないエラーで止めない。
Must preserve: 衝突しない rename の即時反映、feature のない項目の既存ダイアログ、衝突先が色ルールの caption のときの PD-OI-042 の区別、Undo/Redo、Generate と Session での保持。
May retire: feature のある項目の衝突で、UNKNOWN のエラーを出して何もしない動作。
Accepted residual risk: Merge を選ぶと 2 つの凡例項目が 1 つの色と名前にまとまる。
Owner: satoshikawato
Decision date: 2026-09-30
```

## B. Product の判断

### D-07 FE-02/PV-09 batch の即時編集の範囲
- 選択肢:
  - A: 編集時に全 Result へ伝播する。
  - **B（推奨）:** Result を表示するときに正本から投影する。
  - C: 表示中の Result だけに反映されると明示する。

```text
PRODUCT_DECISION
Concern: web.batch.live-edit-projection
Scenario revision: 1
Choice: B / PROJECT-ON-MOUNT
Rationale: batch の各 Result を選んだとき、その Result の preview と出力に、すでに行った色・非表示・凡例・ラベルの編集が反映されているようにする。編集のたびに全 Result を処理し直すことは避ける。
Must preserve: 表示中の Result への即時反映、Undo/Redo で全 Result の見え方が戻ること、Save → Load → Result 選択、各 Result の export、編集がないときの zero fast path、stale・cancel のときの旧 Result の保持、label の DOM identity。
May retire: 表示していない Result に、Generate まで古い色・非表示・凡例が残る動作。
Accepted residual risk: 一度も表示していない Result は、Session に保存される SVG の bytes が、次に表示するか Generate するまで古い（Load 後に選べば投影される）。大きな batch では表示のたびに投影のコストがかかる。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-08 PV-03 凡例の順序の継承
- 選択肢:
  - A: 現状（Generate で既定の順序に戻る）を表示する。
  - **B（推奨）:** 既存の candidate plan に順序の操作を足す。
  - C: request の field にする。

```text
PRODUCT_DECISION
Concern: web.legend.order-continuity
Scenario revision: 1
Choice: B / CARRY-LEGEND-ORDER
Rationale: Sort や Move で整えた凡例の順序を、色や font を直して Generate するたびにやり直さなくてよいようにする（PD-OI-052 と同じ負担をなくす）。
Must preserve: 編集がないときの既定の順序、Sort と Move の即時反映、Undo/Redo、Session の往復、batch の全出力への適用、新しく現れた項目の表示、PD-OI-052 の装飾 delta。
May retire: Generate が凡例の順序を既定に戻す動作。
Accepted residual risk: 並べ替えの後に新しく現れた項目は末尾に置かれる。消えた項目の順序の情報は捨てる。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-09 PV-07 canvas padding の継承
- 選択肢:
  - A: 現状（Generate で 0 に戻る）。
  - **B（推奨）:** Session の padding を、候補を公開する前に全出力へ適用する。
  - C: request の field にする。

```text
PRODUCT_DECISION
Concern: web.canvas.padding-continuity
Scenario revision: 1
Choice: B / CARRY-CANVAS-PADDING
Rationale: PD-OI-052 が clipping の緩和策として示す padding を、Generate のたびに入れ直さなくてよいようにする。
Must preserve: padding の編集と即時反映、Reset、Undo/Redo、Session の往復、batch の全出力、export。padding を二重に適用しないこと。
May retire: Generate が canvas padding を 0 に戻す動作。
Accepted residual risk: 図の大きさが大きく変わる設定変更の後も同じ padding が残るので、余白が合わないことがある。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-10 FE-01 source を置き換えたときの label 編集の扱い
- 選択肢:
  - **A（推奨）:** もう存在しない target の編集だけを外す。
  - B: 現状どおり、どれか 1 つでも target が消えたら全部消す。
  - C: Reset するまで全部残す。

```text
PRODUCT_DECISION
Concern: web.labels.source-replacement-reconciliation
Scenario revision: 1
Choice: A / PRUNE-UNMATCHED-TARGETS
Rationale: 別のゲノムに置き換えて Generate したとき、もう存在しない feature への label の編集だけを外し、残る feature への編集は保つ。
Must preserve: 表示の変化（Result 選択、mount、record 選択、mode、非表示、reflow）では label の override を作成・削除しないこと、bulk の label override（matcher として残す）、Undo による復元、Session の往復。
May retire: source の置き換えのとき、target が 1 つでも消えると label の override をすべて消す動作。
Accepted residual risk: 置き換えた後のゲノムに同じ identity の feature があれば、その override はそのまま適用される。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-11 即時編集と Generate の一致を契約にする
- 選択肢:
  - **A（推奨）:** 契約にする。新しい acceptance contract OIC-027 を足す。
  - B: 契約にしない。

```text
PRODUCT_DECISION
Concern: web.live-edit.regeneration-parity
Scenario revision: 1
Choice: A / LIVE-EDIT-EQUALS-REGENERATION
Rationale: 即時の編集で見えている図が、次の Generate、Session の読み込み、export でも同じになることを保証する。
Must preserve: 各即時編集の応答の速さ、Live edit と Applies on Generate の表示の正確さ（OIC-024）、PD-OI-052 の装飾 delta。即時に編集した Result は、同じ draft から新しく Generate した Result と、対象要素の意味（位置、色、文字、表示）で一致する。Applies on Generate の設定は、Generate の前に Result を変えない。
May retire: 即時の編集と Generate で結果が違ってよいという暗黙の扱い。Generate の compiler で再現できない即時編集は、Applies on Generate に切り替える。
Accepted residual risk: 即時に反映される設定が減ることがある（D-04 と D-05 と同じ方向）。parity を取れない即時編集は退役しうる（D-30）。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-12 IN-06 1 ファイルに複数の生物があるときの定義（Linear）
- 選択肢:
  - A: 全 record で同じときだけ file の既定値にする（当初の推奨）。
  - **B（Owner の選択）:** record ごとに、その record 自身の推定値で定義を付ける。
  - C: 推定値の層を足す。
- 補足: Circular はすでに、record ごとに自分の生物名で定義を作っている（Species 欄に入力した場合だけ、全 record に同じ値が入る）。

```text
PRODUCT_DECISION
Concern: web.linear.file-default-definition
Scenario revision: 1
Choice: B / PER-RECORD-INFERRED-DEFINITION
Rationale: 1 つのファイルに別の生物の record が入っていても、各 record に自分の生物名の定義を付ける。
Must preserve: record ごとの Definition の編集、利用者が file に入力した Definition を全 record に適用すること、全 record が同じ生物のときの現在の表示、Reset で推定値に戻ること、Session の往復、Circular の現在の挙動。定義の優先順位は、record に入力した値 → file に入力した値 → その record 自身の推定値。
May retire: 1 番目の record の推定値を file の既定値として全 record に使う動作。
Accepted residual risk: file 欄の「Using file default」は、利用者が file に入力した値だけを指すようになる。推定値を保存するために Session の形式が変わる場合がある（旧 Session は読み込み時に推定し直す）。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-13 IN-03 FASTA だけにある配列
- 選択肢:
  - **A（推奨）:** GFF の行を持つ record だけを扱う（CLI と同じ）。
  - B: FASTA だけの配列も feature のない record として描く（CLI の出力が変わる）。

```text
PRODUCT_DECISION
Concern: diagram-generation.gff-fasta-record-universe
Scenario revision: 1
Choice: A / GFF-ANNOTATED-RECORDS-ONLY
Rationale: GFF3+FASTA の record の集合を CLI と同じにし、科学的な出力を変えない。
Must preserve: GFF の行を持つ record の表示（feature が 0 でも region 行や埋め込み ##FASTA を持つものを含む）、FASTA の順序、CLI の出力、PD-OI-018 の他の項目。
May retire: GFF の行を 1 つも持たない FASTA の配列を、record の候補として一覧に出す動作（Generate できない候補）。
Accepted residual risk: FASTA だけにある配列は描けない。描くには GFF に region 行を足す必要がある。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-14 FE-09 record ID が重複するときの「This feature only」
- 選択肢:
  - **A（推奨）:** 常に stable hash を出す。
  - B: Python が instance id も照合できるようにする。
  - C: 無効にして理由を示す。

```text
PRODUCT_DECISION
Concern: web.feature-color.duplicate-record-instance
Scenario revision: 1
Choice: A / STABLE-HASH-ONLY
Rationale: Web が作る色ルールを、Python が必ず照合できる値にする（OIPC-C03）。
Must preserve: 重複しない feature の「This feature only」、label の instance 単位の編集、Undo/Redo、Session。
May retire: record ID が重複するとき、rendered instance id を色ルールに書く動作。
Accepted residual risk: 同じ record ID を持つ同一の複製がある場合、「This feature only」の色は両方の複製に適用される。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-15 SE-10 Reset Settings の範囲
- 選択肢:
  - **A（推奨）:** Linear の record ごとの表示状態と alignment plan も初期化する。
  - B: Circular 側も残すようにして揃える。
  - C: 現状を明記する。

```text
PRODUCT_DECISION
Concern: web.reset.linear-record-display
Scenario revision: 1
Choice: A / RESET-LINEAR-RECORD-DISPLAY
Rationale: Reset Settings の範囲を Circular と Linear で揃える。
Must preserve: ファイル、展開された行の record 選択（region_record_id）、file の既定値、depth の割り当て、Undo による復元。
May retire: Reset Settings の後も、Linear の record ごとの Definition・Subtitle・region・逆相補と alignment plan が残る動作。
Accepted residual risk: Reset で Linear の record ごとの表示設定が消える（Undo で戻せる）。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-16 FE-11 検索の Start/End の項目
- 選択肢:
  - **A（推奨）:** 削除し、Location で検索する。
  - B: 残して 1 始まりに直す。

```text
PRODUCT_DECISION
Concern: web.feature-search.location-fields
Scenario revision: 1
Choice: A / LOCATION-FIELD-ONLY
Rationale: 位置の検索と表示を 1 始まりの INSDC 形式に揃え、誤った一致をなくす。
Must preserve: Location での検索（原点をまたぐ feature と分割された feature を含む）、drawer と popup の位置と長さの表示。
May retire: 検索の Start と End の項目（0 始まりの生の値）。
Accepted residual risk: 開始位置の数値だけで検索していた場合は、Location の値で検索し直す必要がある。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-17 PV-05 Web の PDF の物理的な大きさ
- 選択肢:
  - **A（推奨）:** CSS の px を pt に換算する（CLI と PNG に合わせる）。
  - B: 現状（1 px = 1 pt）を維持する。

```text
PRODUCT_DECISION
Concern: web.export.pdf-physical-size
Scenario revision: 1
Choice: A / CSS-PX-TO-PT
Rationale: Web の PDF の物理的な大きさを、CLI（CairoSVG）の PDF と、PNG の DPI に揃える。
Must preserve: PDF の見た目、文字の抽出、ページが 1 枚であること、Web の PNG と SVG の大きさ。
May retire: Web の PDF を 1 px = 1 pt で作る動作。
Accepted residual risk: Web で作る PDF の物理サイズは、これまでの 75% になる。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-18 Q-FRAME 比較表の座標系（CO-07、CO-03(b)、CO-10、N-08、N-07）
- 座標系の定義:
  - F（探索座標系）: record 選択と crop の後の配列。1 始まりで、元の鎖の向き。LOSAT の出力そのもの。
  - V（表示座標系）: F に実効の逆相補を適用したもの。
- 選択肢:
  - A: 表は V のまま。Save Raw は V に変換して出す。
  - B: 表は V、raw は F のまま。座標系を明示する。
  - **C（推奨）:** すべての比較表を F にし、planner が向きを投影する。
  - D: source の絶対座標にする。

```text
PRODUCT_DECISION
Concern: comparison.table-coordinate-frame
Scenario revision: 1
Choice: C / SEARCH-FRAME-EVERYWHERE
Rationale: 比較表の座標を表示の向きに関係なく同じ意味にし、逆相補・回転・再アップロード・CLI で、同じ表が同じ相同領域を指すようにする。
Must preserve: LOSAT の raw cache と Save Raw の内容、crop の意味（表は crop 後の record 内の座標）、feature binding を持つ protein 比較の投影、向きを変えない場合の既存の BLAST 表の結果、main で保存された Session の読み込み（読み込み時に一度だけ変換する）。表の座標が record の範囲外なら、検証で止めるか警告する。
May retire: アップロードした表と CLI の -b を表示座標（逆相補の後）として読む動作、JS 側の探索座標から表示座標への変換、範囲外の行を record の外に描く動作。
Accepted residual risk: 逆相補と -b を組み合わせていた CLI の利用者にとって、表の座標の意味が変わる（release note に書く）。crop より前の全長の座標で作った表は、crop すると範囲外として止まる。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-19 CO-04 LOSAT の database の範囲（PD-OI-018 の revision 4）
- 選択肢:
  - A: database を subject の 1 record に限る（当初の推奨）。
  - **B（Owner の選択）:** source ファイル単位の batch を維持し、要求していない自己検索を除く。
  - C: 上限の有無で方式を切り替える。
- Owner の指摘: 今の PD-OI-018 は B に当たる。複数 replicon のゲノムでは、染色体ごとではなく、ゲノム（＝ファイル）単位の E-value が望ましい。
- CO-04 でリンクが 0 本になった原因は、要求していない自己一致（record 自身への hit）が Max target seqs を埋めたこと。

```text
PRODUCT_DECISION
Concern: diagram-generation.linear-record-universe-and-search-scope
Scenario revision: 4
Supersedes: PD-OI-018, scenario revision 3
Choice: B / FILE-DATABASE-WITHOUT-UNREQUESTED-SELF
Rationale: 1 つのファイルを 1 つのゲノムとして扱い、複数 replicon のゲノムでもゲノム単位の E-value で比較する。そのうえで、要求していない自己一致が Max target seqs を埋めてリンクが消えることを防ぐ。
Must preserve: revision 3 の item 1〜7（source ファイル単位の batch、database の範囲を raw cache の identity に含めることを含む）。表示に使わない自己検索（record 自身への検索）は、どのモードでも実行しない。同じ source の中の record どうしの比較は、query の record を除いた database で検索する。E-value の database は subject 側の source ファイル（query と同じ source なら query の record を除いたもの）であることを、docs と Run Info に明記する。job 数の見積もりは、実行の計画と同じ関数から出す。
May retire: 自己検索を除くのが Collinear の inference が OFF のときだけという限定と、要求していない自己一致で Max target seqs が埋まりリンクが消える動作。
Accepted residual risk: 同じ record でも、ファイルの分け方を変えると E-value が変わる（1 ファイル = 1 ゲノムという前提）。record の対ごとに検索する CLI とは、複数 record のファイルで E-value が一致しない。同じ source の中の比較は別の job になり、job 数が増える。
Owner: satoshikawato
Decision date: 2026-09-30
```

D-19 に関連して、CLI は今回変えない（D-40）。

### D-20 CO-06 アップロードした表の record ID の照合
- 選択肢:
  - A: 厳格に照合する。
  - **B（推奨）:** 矛盾はエラー、不明な ID は警告にする。
  - C: metadata だけ直す。

```text
PRODUCT_DECISION
Concern: comparison.uploaded-table-record-binding
Scenario revision: 1
Choice: B / CONTRADICTION-ERROR-UNKNOWN-WARN
Rationale: query と subject を取り違えた表で、誤った相同領域を黙って描かないようにする。
Must preserve: record ID と一致する表の結果、version 接尾辞の違い（.1 など）を許すこと、ID が record と無関係な表の位置による割り当て（警告付き）、CLI と Web で同じ結果。
May retire: 端点と逆の record を指す行や、ほかの record を指す行を、位置のまま描く動作と、metadata の ID と index の矛盾。
Accepted residual risk: ID が record と無関係な表は、今と同じく位置で割り当てられる（警告は出る）。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-21 CO-08 Similarity group の名前の継承
- 選択肢:
  - **A（推奨）:** member 集合の完全一致で付け替え、一致しないものは dormant として保存する。
  - B: 一致しないものは通知付きで削除する。
  - C: 重なりで引き継ぐ。

```text
PRODUCT_DECISION
Concern: web.similarity-group.override-identity
Scenario revision: 1
Choice: A / REKEY-BY-MEMBERSET-KEEP-DORMANT
Rationale: 利用者が付けた group の名前と説明を同じ member の group に付け続け、別の group へ移したり黙って消したりしない。
Must preserve: group の名前と説明の編集、Session の往復、Undo/Redo、Interactive SVG への出力、og_* の ID の表示。
May retire: ID の文字列だけを頼りに名前を残す動作と、残らない名前を黙って消す動作。
Accepted residual risk: member が 1 つでも変わった group には名前が付かない（dormant として保存し、一覧と Clear から扱える）。Session に項目が 1 つ増える。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-22 CO-10 match popup の座標
- 選択肢:
  - A: source の座標だけを出す。
  - **B（推奨）:** source の座標を主にし、違うときだけ表の座標も出す。
  - C: ラベルだけ付ける。

```text
PRODUCT_DECISION
Concern: web.match-popup.coordinates
Scenario revision: 1
Choice: B / SOURCE-PRIMARY-WITH-TABLE
Rationale: match popup と FASTA ヘッダの座標を、feature popup と同じ入力ファイルの座標にする。
Must preserve: 取り出す配列そのもの、逆鎖の扱い、表の座標の参照（違うときだけ併記）、Interactive SVG の popup。
May retire: match popup と FASTA ヘッダが crop 後の表示座標だけを出す動作。
Accepted residual risk: popup の行が 1 行増えることがある。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-23 CO-05 12 列を超える outfmt 6
- 選択肢:
  - A: ちょうど 12 列だけを受け付ける。
  - **B（推奨）:** 先頭 12 列が型どおりなら受け付ける。

```text
PRODUCT_DECISION
Concern: comparison.outfmt6-extra-columns
Scenario revision: 1
Choice: B / FIRST-12-COLUMNS
Rationale: -outfmt "6 std qlen slen" のように列を足した表を、CLI と Web でそのまま使えるようにする。
Must preserve: 12 列の表の結果、outfmt 7 のコメント行、空のファイル、先頭 12 列の型の検証、CLI・Web・conservation で同じ規則。
May retire: 列を足した表を黙って誤読する動作。存在しないファイルや読めないファイルを飛ばして、後ろの比較をずらす動作。
Accepted residual risk: 13 列目以降は使わずに捨てる（INFO ログを出す）。存在しないファイルを渡していた CLI の実行は失敗に変わる。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-24 PV-08 定義文字列が入らないとき
- 選択肢:
  - A: エラーの案内だけを改善する。
  - **B（推奨）:** 配置に失敗したときだけ、種名の行を単語の区切りで折り返して配置し直す。
  - C: フォントを縮小する。

```text
PRODUCT_DECISION
Concern: diagram-generation.circular-definition-fit
Scenario revision: 1
Choice: B / WRAP-DEFINITION-ON-FIT-FAILURE
Rationale: よくある細菌の長い学名（subsp.、serovar、str. などを含むもの）でも、Web の既定の設定で図を作れるようにする。
Must preserve: これまで成功していた出力（折り返さない）、定義の文字そのもの、center_reserved_radius と definition_font_size を明示した場合の扱い（折り返しを適用しない）、失敗したときの案内。
May retire: 定義の円が入りきらないとき、配置し直さずに失敗する動作。
Accepted residual risk: 折り返しても入らない場合は今までどおり失敗し、改善した案内を出す。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-25 SE-05 旧形式の Session の Save
- 選択肢:
  - **A（推奨）:** Generate してから Save する（読み込み時と Save 時に案内する）。
  - B: Save のときに catalog を回復する（証拠が必要）。
  - C: 読み込みのときに catalog を作る。

```text
PRODUCT_DECISION
Concern: web.session.legacy-save
Scenario revision: 1
Choice: A / REQUIRE-GENERATE-BEFORE-SAVE
Rationale: 旧形式の Session から、feature の identity が確かでない状態のまま現行の形式を書き出さない。
Must preserve: 旧形式の Session の読み込みと preview、Generate 後の Save、v40 以降の Session の Save、PD-OI-045 の Session 操作。Save が必要とする Generate を案内し、エラーパネルから Generate を実行できるようにする。
May retire: 0.13.0 で可能だった「旧形式の Session を読み込んで、そのまま Save する」操作。
Accepted residual risk: 旧形式の Session を保存し直すには 1 回 Generate が必要。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-26 Dinucleotide の文字集合
- 選択肢:
  - A: ACGT の 2 文字（大小を区別しない）に限る（当初の推奨）。
  - **B（Owner の選択）:** U も許し、T と同じ塩基として扱う。

```text
PRODUCT_DECISION
Concern: options.dinucleotide-alphabet
Scenario revision: 1
Choice: B / ACGTU-PAIRS
Rationale: 無効な指定で空や平坦な track を黙って描かないようにし、CLI の traceback もなくす。RNA の表記（U）でも指定できるようにする。
Must preserve: ACGT の 2 文字の指定（大小を区別しない）、slot の nt、CLI と Web で同じ検証。U は T と同じ塩基として扱う（AU は AT と同じ結果になり、配列中の U も T として数える）。
May retire: 2 文字でない指定を黙って既定値に戻す動作、XY のような塩基でない文字の受理、G での IndexError。
Accepted residual risk: N などの曖昧な塩基の記号は指定できない。凡例などの表示名は入力した文字（AU）のまま出す。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-27 負の値の検証の範囲
- 選択肢:
  - **A（推奨）:** 仕様で意味が決まる値だけを検証する（フォントサイズ > 0、stroke 幅 ≥ 0）。
  - B: offset、spacing なども調べて検証する。

```text
PRODUCT_DECISION
Concern: options.nonnegative-style-values
Scenario revision: 1
Choice: A / SPEC-DETERMINED-ONLY
Rationale: SVG と CSS の仕様で意味が決まる値だけを検証し、意味を確かめていない値は変えない。
Must preserve: offset、spacing、track_axis_gap、label_rotation の現在の受理範囲。CLI、Python API、Web、Session で同じ検証とエラー。
May retire: 0 以下のフォントサイズと負の stroke 幅の受理、描画の途中の ValueError による traceback。
Accepted residual risk: 負の offset などの意味は、今回は確かめない。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-28 GE-06 Generate 中の Undo/Redo
- 選択肢:
  - **A（推奨）:** History の artifact の置き換えか checkpoint が開いている間は、busy として拒否する。
  - B: 先に Cancel してから Undo する。
- 推奨 A の範囲: PD-OI-051 が保証する Generate 中の draft の編集は、そのまま許す。

```text
PRODUCT_DECISION
Concern: web.generation.in-flight-history
Scenario revision: 1
Choice: A / BUSY-UNDO-DURING-GENERATE
Rationale: Generate 中に Undo や Redo を押しても、確定済みの request と Result が古いものに戻ったり、実行中の Generate が黙って捨てられたりしないようにする。
Must preserve: Generate の Cancel、処理中の表示、Generate が終わった後の Undo/Redo、Save と Load の拒否（既存）、PD-OI-051 による Generate 中の draft の編集。Undo と Redo のボタンとショートカットは同じ判定を使い、拒否した理由を示す。
May retire: Generate 中の Undo と Redo。
Accepted residual risk: 長い LOSAT の実行中は、Cancel するか終わるまで Undo できない。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-29 TR-03 Linear で source のない series を手動の行が参照するとき
- 選択肢:
  - **A（推奨）:** Circular と同じ row issue（source なし）を示し、Generate の前に止める。
  - B: 空の track を描く。

```text
PRODUCT_DECISION
Concern: web.depth.linear-sourceless-series
Scenario revision: 1
Choice: A / ROW-ISSUE-BEFORE-GENERATE
Rationale: Linear で File の Depth を消した後、source を持たない series を有効な手動の行が参照していても、原因の分からない失敗にしない。
Must preserve: PD-OI-025（論理 series の保持、File 単位の apply と clear が 1 つの undoable 操作）、OIC-020、Circular の row issue と同じ文言。
May retire: この状態の Generate が汎用のエラーで失敗する動作。
Accepted residual risk: 利用者が行を無効にするか削除する必要がある。
Owner: satoshikawato
Decision date: 2026-09-30
```

### D-30 PV-10 Linear の凡例の side の即時移動（条件付き）
- 選択肢:
  - **A（推奨）:** 次のセッションで、即時移動と Generate の差を計測する。JS の reflow を Python に合わせて一致させ、一致させられない場合だけ Linear の即時移動を Applies on Generate にする。
  - B: 常に Generate で再描画する。

```text
PRODUCT_DECISION
Concern: web.legend.live-side-move
Scenario revision: 1
Choice: A / PARITY-OR-APPLY-ON-GENERATE
Rationale: 即時に見えた凡例の配置と、Generate 後の配置が食い違わないようにする（D-11 の契約）。
Must preserve: Circular の凡例の side の即時移動、凡例の drag、PD-OI-052 の装飾 delta、Undo/Redo。Linear で parity を取れる場合は、即時移動も残す。
May retire: parity を取れない場合に限り、Linear の凡例の side の即時移動（凡例を持たない図での side 変更の例外を含む）。
Accepted residual risk: 退役した場合、Linear では side を変えても Generate まで preview が変わらない。
Owner: satoshikawato
Decision date: 2026-09-30
```

## C. 現状を維持する項目（承認すれば新しい record は作らない）

| ID | 項目 | 推奨 | 理由 |
|---|---|---|---|
| D-31 | PV-01 side を変えたときの drag の delta | 現状維持 | PD-OI-052 が受容済み。はみ出しは padding（D-09）と Reset で直す |
| D-32 | IN-04 単一 record から単一 record への置き換えで crop と定義を残す | 現状維持 | 同じ genome の版を差し替える操作を妨げない |
| D-33 | IN-02 BOM で始まる GenBank | JS は Python に合わせて判断を Worker に回す。Python での受理は足さない | CLI と同じ失敗になり、案内も出る |
| D-34 | TR-09 最初に stack を有効にしたとき、保存済みの既定 stack を使う | 現状維持 | UI の契約どおり |
| D-35 | SE-08(b) CLI Session の読み込み時の Python 検証 | 現状維持。CLI が差分の configOverrides を書く案は、replay で結果が変わらないことを証明してから別の課題にする | PD-OI-009 の検証を保つ |
| D-36 | SE-06 CLI Session の比較を EDITABLE にするか | しない（読み取り専用のまま、Inherit で使う） | PD-OI-008 |
| D-37 | Generate 中の overlay の後ろでのキーボード編集 | 現状維持（許す） | PD-OI-051 が Generate 中の draft の編集を保証している |
| D-38 | エラーの field へ移動するボタン | 足さない（`aria-invalid` による表示だけ） | 新しい操作は不要 |
| D-39 | 色の 4 桁・8 桁の hex（α 付き） | 足さない（`none` の受理だけを直す） | 新しい能力は不要 |
| D-40 | CLI の LOSAT の database（D-19 の関連） | 今回は変えない（record の対ごとの検索のまま）。複数 record のファイルで Web と E-value が違うことを docs に書く | CLI の科学的出力を変えない |

## W. 作業の進め方（Product authority ではない）

| ID | 項目 | 選択肢 | 推奨 |
|---|---|---|---|
| W-1 | main にも含まれる不具合 | — | **決定済み:** dev で直し、次の dev → main の昇格で main に入れる。backport はしない |
| W-2 | 範囲 | A: 65 件と N-01〜N-20 のすべて / B: 65 件だけ | A |
| W-3 | PR の粒度 | A: 約 50 本の小さな PR（size profile の範囲内） / B: 約 15 本の workstream 単位の PR（size の Review REQUIRED を受け入れる） | B。1 セッションで終えるには A は遅い。size の超過は Gate を失敗させない |
| W-4 | authority の登録 | A: 最初に 1 本の authority だけの PR で PD-OI-056 以降と OIC-027 をまとめて登録する / B: 実装の PR ごとに co-change する | A |
| W-5 | dev への自動マージ | A: すべての PR を、必須の check が通ったら自動でマージする（main は対象外） / B: Product Contract、参照出力、architecture-change の PR は Owner のレビューを待つ | A。receipt はこのセッションで承認済みになるため。main への昇格は別に判断する |
| W-6 | 参照出力と Gallery の作り直し | 許可する / しない | 許可する。PV-08 の修正で Web 既定の Circular 出力がすべて変わり、FE-06 と FE-11 で Interactive SVG の runtime が変わるため。目視の確認を記録する |
| W-7 | 次のセッションでの並列化 | A: Workflow（複数のサブエージェント）で、ファイルが重ならない workstream を並列に進める / B: 1 つのエージェントで順に進める | A。ホットなファイルを触る PR は直列にする。トークンの消費は大きい |
| W-8 | 再発防止のガード（G-A〜G-J） | A: すべて入れる / B: 修正に付随するテストだけ | A |

## 判断の記録

**Status: 確定（2026-09-30、Owner: satoshikawato）。**

Owner の発言は次の 2 つで、翻訳も要約もせずに記す。

1. D-01〜D-39 と W-1〜W-8 を提示したときの回答:
   > だいたいそのままで承認。けど、D-12: 1ファイルに複数生物あるときは、Definitionはそれぞれにつけてくれるとうれしいかも。Circularでしょ？ D-19: 今まではBだったんじゃないの？たとえばmultiple repliconのゲノムだったら、染色体ごとじゃなくてゲノム=ファイル単位のE-valueが欲しいんじゃないかな？ D-26: UはOKにして。
2. D-12、D-19、D-26 の書き直した receipt と、D-19 の CLI の扱い（推奨 A）を提示したときの回答:
   > OK,これで承認します。

| 区分 | 項目 | 結果 |
|---|---|---|
| A | D-01〜D-06 | 上の receipt の文面のとおり承認 |
| B | D-07〜D-11、D-13〜D-18、D-20〜D-25、D-27〜D-30 | 推奨の receipt の文面のとおり承認 |
| B | D-12 | B / PER-RECORD-INFERRED-DEFINITION（上の文面） |
| B | D-19 | B / FILE-DATABASE-WITHOUT-UNREQUESTED-SELF（上の文面。PD-OI-018 revision 4） |
| B | D-26 | B / ACGTU-PAIRS（上の文面） |
| C | D-31〜D-39 | 推奨のとおり現状を維持。新しい record は作らない |
| C | D-40 | CLI は変えず、docs に違いを書く。2 つ目の回答で承認されたものとして扱う |
| W | W-1 | dev で直し、次の dev → main の昇格で main に入れる（backport なし） |
| W | W-2〜W-8 | 推奨のとおり（範囲は N-01〜N-20 を含む全件、PR は workstream 単位、authority を先に 1 本、dev への自動マージ、参照出力と Gallery の作り直し、Workflow による並列化、ガードをすべて入れる） |

### 登録の方法（次のセッション）

- **登録する record:** 区分 A と B の receipt 30 件。
  - D-19 は、PD-OI-018 の本文を revision 4 に置き換える。supersession の書式と、前の revision を Git 履歴に残す規則に従う。
  - それ以外の 29 件は、PD-OI-056 から番号を振る。
  - D-11 に対応して、acceptance contract catalog に OIC-027 を足す。
- **receipt の書き方:** 既存の record（PD-OI-046、PD-OI-052 など）と同じ形にする。
  - 各 receipt の text block を翻訳も追加もせずに再現する。
  - JSON の表現と、UTF-8 の SHA-256（最後の改行を除く）を付ける。
  - Decision source には、この文書のパスと、上の 2 つの発言を書く。
- **PR の出し方:** authority だけの PR（W-4）として出し、`Authority metadata` の revision の一覧も更新する。
- **C 区分:** D-31〜D-40 は record を作らない。関係する docs（D-40 の CLI と Web の E-value の違いなど）は、実装の PR で更新する。
