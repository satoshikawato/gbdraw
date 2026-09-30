<!-- Raw design report of workstream W6 (sub-agent output, 2026-09-30, base origin/dev 4c89bab1). Evidence for 01_REMEDIATION_PROPOSAL.md, not Product or architecture authority. -->

# W6 修正提案：Python core、描画、Web PDF

前提として、行番号はすべて `4c89bab1` のもので、パスはリポジトリ直下からの相対パスです。dev worktree は読み取り専用で扱い、変更、commit、issue 作成はしていません。検証には dev の CLI と、プロセス内の monkeypatch による probe を使いました。probe は `…/scratchpad/design/w6/` の probe*.py と fixprobe.py にあります。

**主な発見は 3 つで、いずれも監査の記述と異なります。**
- **PV-08:** 主な原因は定義文字列ではなく、Multi-Record Canvas（MRC）が depth 入力なしでも空の depth スロットを確保していることです。この回帰は 021b4c49 から入っています。影響は Web 既定の circular 出力すべてに及び、Generate に成功している図でも幾何が単一 record の経路と違います。
- **CO-05:** conservation の BLAST ファイル読み込みにも同じ不具合があります。監査は対処済みとしていましたが誤りです。
- **FE-07:** GFF3 の phase を読んでいないため、LOSATP の protein 入力が読み枠のずれた配列になります。これは監査にない新しい不具合です。

---

## CO-05（P1）13 列以上の outfmt 6 が誤読される

**1. 分類**
- **誤読をなくすこと自体:** IMPLEMENT_EXISTING_AUTHORITY です。根拠は `docs/REFERENCE/input-formats-and-tsv-schemas.md:28` の「BLAST-compatible input uses the 12 outfmt 6 columns … qseqid … bitscore」です。列がずれた ribbon を黙って描く現状は、科学的出力の整合性の点で残せません。
- **12 列を超える表の扱い:** PRODUCT_DECISION_REQUIRED です。受け付ける入力の範囲を広げるかどうかの判断になるためです。
  - A: ちょうど 12 列だけを受け付ける。12 列を超える表も足りない表も、ファイル名、行、列数を示す ValidationError にする。
  - B: 先頭 12 列が std の型として読める場合に限り、13 列以上も受け付ける。末尾の列は INFO ログを出して無視する。docs の :28 は「先頭 12 列」に改める。
  - 推奨は B です。`-outfmt "6 std qlen slen"` はよく使われます。conservation の DataFrame 経路（`gbdraw/analysis/conservation.py:147-149`）はすでに先頭 12 列を採っているので、B ならファイルと DataFrame が 1 つの規則になります。決定が出るまでは、既存の権限だけで実装できる A で先に出せます。
- **付随する不具合:** `load_comparisons` は、存在しないファイルや解析できないファイルを warning だけ出して飛ばします。その結果、後ろの比較が 1 つ前の record 対に描かれます。誤った対に ribbon を描くことは許されないので、これも IMPLEMENT_EXISTING_AUTHORITY です。なお `--comparison_table` の経路（`gbdraw/api/record_planning.py:1432-1436`）は、すでに ValidationError を出します。

**2. 根本原因の再確認**
- pandas で再現しました。`read_csv(names=12列)` に 13 列の表を渡すと、先頭列が index になります。その結果 query は `SB`、identity は 200、qstart は 300 になり、bitscore には 13 列目の 3000 が入ります。14 列では先頭 2 列が MultiIndex になり、evalue に 3000 が入って閾値ですべて落ちます（ribbon 0 本）。ここまでは監査のとおりです。
- **監査の誤り:** 「`conservation.py:147-149` は対処済み」は誤りです。ファイル入力は `conservation.py:175-180` で同じ `names=COMPARISON_COLUMNS` で読まれ、列がずれた状態で `_coerce_comparison_dataframe` に渡ります。13 列と 14 列のファイルを `_load_conservation_file` に渡すと、閾値を evalue≤10 に緩めても 0 hits でした。147-149 が効くのは、Python API から位置列の DataFrame を渡した場合だけです。
- 同じ読み方をしている箇所は 4 つです。
  - `gbdraw/io/comparisons.py:73-78`（CLI の `-b`、schema<2 の request）
  - `gbdraw/session_request_codec.py:3330-3335`（Web の nucleotideBlast、schema≥2）
  - `gbdraw/api/record_planning.py:1426-1431`（CLI の `--comparison_table`）
  - `gbdraw/analysis/conservation.py:175-180`（CLI と Web の similarity ring）
- 影響しない箇所もあります。
  - `gbdraw/web/js/app/python-helpers.js:453-457` も同じ書き方ですが、入力は内部の LOSAT 出力で常に 12 列です。
  - LOSATP の parser（`protein_colinearity.py:3070-3103`）は 12 列を厳格に要求します。
  - JS はアップロードされた表を解析せず、Python に resource として渡すだけです。
- 付随する不具合も再現しました。`linear --gbk SA SB SC -b missing.tsv SA_SB.tsv` は warning だけ出して成功し、SVG には `comparison1` しかありません（`io/comparisons.py:66-71, 80-87`）。空ファイルは `names=` を指定していれば空の DataFrame になり、位置はずれません。ただし `header=None` に変えると EmptyDataError になるので、処理が必要です。

**3. 修正案（コード）**
- 担当を `gbdraw/io/comparisons.py` の 1 か所にまとめます。
  - `read_comparison_table(source, *, label)` を追加します。パスとテキストの両方を受け付けます。
  - `normalize_comparison_dataframe(df)` を追加します。これは `conservation.py:139-158` の `_coerce_comparison_dataframe` を移したものです。
- 読み込みは `pd.read_csv(..., sep="\t", comment="#", header=None, dtype=str)` とし、次のように扱います。

| 入力 | 扱い |
|---|---|
| 空ファイル（EmptyDataError） | 12 列の空の DataFrame |
| 途中の行だけ列が多い（ParserError） | 行番号つきの ValidationError |
| 列数が 12 未満 | ValidationError |
| 列数が 12 超 | A ならエラー、B なら先頭 12 列を採って INFO |
| 先頭 12 列に欠損（短い行） | 行番号つきの ValidationError |
| 数値列に数値でない値 | `pd.to_numeric(errors="raise")` の失敗を、列名と行を示す ValidationError にする |

  返す DataFrame は RangeIndex で、列は COMPARISON_COLUMNS の 12 列に固定します。
- 上の 4 か所の読み込みを、この関数の呼び出しに置き換えます。conservation からは `_coerce_comparison_dataframe` を削除します。
- `load_comparisons` の 2 つの skip（存在しない :66-71、解析失敗 :80-87）と、例外を握りつぶす `except Exception` を ValidationError に変えます。
- codec は reader のメッセージを CanonicalRequestDecodingError に残すようにします。現状は "Could not decode BLAST resource" だけになり、詳細が消えます。
- `gbdraw/web_support/error_adapter.py` に列数エラーの template を 1 行追加します（COMPARISON_INPUT/COLUMN_COUNT）。JS 側の文言は W3（X-01）の catalog と合わせます。
- 追加しないもの: outfmt 7 の "# Fields:" 行を使った列名の対応付け、新しいオプション、JS 側での事前検証。
- 任意の追加: `python-helpers.js:453` をテキスト入力版の関数に置き換えると、5 つ目の複製もなくなります。

**4. 上位の対応**
- 「BLAST 表を読むのは `io/comparisons.py` の関数だけ」を担当の約束にします。
- `names=COMPARISON_COLUMNS` が `gbdraw/` 内でこの 1 か所にしかないことを、静的テストで守ります。`tests/test_dead_api_cleanup.py` と同じ種類のテストです。

**5. テスト**
- `tests/test_comparisons.py` に、次の入力を並べた parametrize テストを追加します。
  - 12、13、14 列の表
  - outfmt 7 に列を足した表
  - 11 列の表
  - 途中の行だけ列が多い表、少ない表
  - 数値列に数値でない値がある表
- 今回の不具合を捕まえる assert は次の 2 つです。B なら 13 列と 14 列の入力で成り立ちます。A なら、同じ入力で列数を示すエラーになることを確かめます。
  - `df.loc[0, ["query","identity","qstart","evalue"]].tolist() == ["SA", 95.0, 101, 1e-50]`
  - `isinstance(df.index, pd.RangeIndex)`
- 同じ fixture を 4 つの入口で確かめます。
  - `load_comparisons`
  - `_decode_comparisons`（既存の codec テストの仕組みを使う）
  - `--comparison_table`
  - `_load_conservation_file`
- CLI では、存在しないファイルを渡すと非ゼロで終了することを確かめます。
- Web は Python の経路なので、新しい Playwright spec は不要です。wheel を作り直したうえで、既存の比較 spec で回帰を確認します。
- `tests/reference_outputs/` の再生成は不要です（12 列を超える参照出力はありません）。

**6. 規模と PR**
- Python の本体は +60/−45 行、テストは +100 行程度です。
- 1 PR にまとめます。科学的出力が変わるので Review REQUIRED です。

**7. 依存とリスク**
- B を選ぶ場合は docs の改訂が必要です。
- 13 列の表を含む既存の Session は、次の Generate で ribbon が変わります（正しい表示になります）。保存形式は変わらないので移行は不要です。
- CLI で存在しないファイルを渡していた利用者は、成功から失敗に変わります。意図した変更として CHANGELOG に書きます。

---

## TR-01（P1）MRC の凡例が custom track slot を無視する

**1. 分類**
- IMPLEMENT_EXISTING_AUTHORITY です。
- 根拠は `docs/CLI_Reference.md` の次の箇所です。
  - :758 と :781: `legend_label`、`positive_color`、`negative_color` は slot の renderer パラメータ。
  - :1348-1350: custom stack が有効なときは slot が凡例の権威。
- PD-OI-019 は Web の既定を MRC にしましたが、MRC で凡例の意味を変える決定はありません。

**2. 根本原因の再確認**
- 監査のとおりでした。HmmtDNA に次の 3 つの slot を設定して、CLI で再現しました。
  - gc_content: `legend_label=MY GC`
  - gc_skew: 正を `#ff0000`、負を `#0000ff`
  - at_skew（`nt=AT`）を追加

| 経路 | 凡例に出た内容 |
|---|---|
| 単一 record | `MY GC`、GC skew ±（#ff0000/#0000ff）、AT skew ± |
| `--multi_record_canvas` | `GC content`、GC skew ±（#6dded3/#ad72e3）。AT skew はない |

- 原因の流れは次のとおりです。
  - MRC は各 record を `legend="none"` で組み立てます（`gbdraw/api/diagram.py:2962-2972`）。
  - そのうえで `diagram.py:3101-3158` で凡例を別に組み立てます。
  - この経路は `_sync_legend_table_for_circular_slots`（`gbdraw/diagrams/circular/assemble.py:1277-1414`）も、`sync_annotation_legend_entries`（`assemble.py:2783`）も呼びません。
- **監査にない、同じ原因の欠落:**
  - Region annotation の `legend_label` も MRC で消えます。bracket の行に MY REGION を付けると、単一 record では出て、MRC では出ませんでした（実験で確認）。
  - depth slot の `legend_label` も同じ経路で無視されます（コードからの判断）。
- `diagram.py:3134-3158` の conservation gradient の生成は、`assemble.py:1392-1413` と重複しています。

**3. 修正案（コード）**
- `assemble.py:2729-2788` にある凡例の組み立てを、同じファイルの 1 つの関数に抜き出します。
  - 組み立ての順序は check_feature_presence → precompute_used_color_rules → prepare_legend_table → depth の同期 → slot の同期 → annotation の同期です。
  - 関数は `build_circular_legend_table(records: Sequence[SeqRecord], *, feature_config, gc_config, skew_config, depth_config, depth_df, depth_tracks, slots, annotations, conservation_tracks, conservation_min_identity, cfg, profile)` の形にします。
  - check_feature_presence と precompute_used_color_rules は、すでに record のリストを受け付けます。
- 単一 record の経路は、`[gb_record]` と `effective_circular_track_slots` でこの関数を呼びます。振る舞いは変わりません。
- MRC は `records` で呼び、残りの引数は次のものを渡します。
  - slot: `prepare_annotation_track_slots(resolved_annotations, records, parsed_circular_track_slots, mode="circular", …)`（`gbdraw/annotations/planning.py:22`、複数 record に対応済み）で得た実際の slot
  - depth: `representative_depth_tracks(record_depth_track_data)`
  - conservation: `first_record_conservation_tracks`
- `diagram.py:3101-3158` の独自の組み立ては削除します。depth の同期と conservation の重複もこれでなくなります。
- 追加しないもの: record ごとに作った凡例を後で結合する仕組み（並び順が変わるため）、新しい凡例モデル。

**4. 上位の対応**
- 「circular の凡例表は `build_circular_legend_table` だけが作る」を担当の約束にします。
- MRC は配置だけを自分で行い、凡例の意味は単一 record と同じ関数から得ます。

**5. テスト**
- `tests/test_circular_multi_canvas.py` に、1 record での一致テストを parametrize で追加します。
  - `assemble_circular_diagram_from_records([rec], …)` と `assemble_circular_diagram_from_record(rec, …)` の凡例（テキスト列と fill 列）が一致すること。
  - ケースは、上の 3 slot、annotation の legend_label、depth slot の legend_label、色付きの conservation です。
  - 現状では最初のケースで失敗します。
- 2 record でも slot 由来の項目が出ることを確かめます。既存の `test_assemble_circular_diagram_from_records_shared_legend_and_unique_ids` の隣に置きます。
- CLI では `--multi_record_canvas --circular_track_slot …` の凡例を確かめます。
- Web では Playwright に 1 ケース追加し、既定（MRC ON）で `MY GC` と AT skew が凡例に出ることを確かめます。監査の spec `tests/web/audit-tracks/mrc-legend-disabled` を assert 付きに直して使います。
- `tests/reference_outputs/` の再生成は不要です（16 件に MRC はありません）。
- MRC と conservation を組み合わせた凡例は、並び順が単一 record と同じになるので変わります。`test_circular_conservation.py` の MRC 凡例テストを確認します。

**6. 規模と PR**
- Python は +40/−70 行、テストは +120 行程度です。
- 1 PR にまとめ、PV-08 の後に出します。
- 凡例の担当を移すので architecture-bearing で、Review REQUIRED です。ratchet は ordinary non-increasing（重複した経路を削除する）として、担当と経路の根拠を記録します。

**7. 依存とリスク**
- PV-08 と同じ関数（`assemble_circular_diagram_from_records`）を変更します。
- Gallery の Vnig_TUMSAT-TG-2018 は MRC ですが、slot も annotation もないので、凡例は変わらない見込みです（要確認）。

---

## PV-08（P1）既定設定では、長い /organism を持つ record を Generate できない

**1. 分類**
- **主な原因（MRC の空の depth リング）:** IMPLEMENT_EXISTING_AUTHORITY です。
  - `docs/RELEASE_NOTES_0.14.0b0.md` の sparse depth の節は、帯を確保するのは「depth 入力がある中で、ある record のセルが欠けている場合」だとしています。
  - 同じ節で、どの record にもファイルがないグループは無効、dense な入力の幾何は変わらない、とも述べています。
  - 021b4c49（2026-07-19）による回帰です。
- **エラーメッセージの改善:** IMPLEMENT_EXISTING_AUTHORITY です。既存の対処法を事実どおり示すだけです。
- **定義文字列が本当に収まらない場合の自動回避:** PRODUCT_DECISION_REQUIRED です。
  - A: 幾何は変えず、エラーだけ改善する。
  - B: 配置に失敗したときだけ、種名の行を単語の区切りで折り返して配置し直す。文字は変えず、これまで成功していた出力も変わらない。`center_reserved_radius` が明示されていれば適用しない。
  - C: 配置に失敗したときだけ、定義のフォントを縮小する。`definition_font_size` が明示されていれば適用しない。
  - 推奨は B です。表示内容を保ち、参照出力も動かさずに、Web 既定でよく見る細菌名（「… subsp. … serovar … str. …」）を表示できます。決定が出るまでは A にします。

**2. 根本原因の再確認（監査の仮説は誤り）**
- 定義文字列のために中央に確保する半径は、単一 record でも MRC でも同じでした（Salmonella 241.9 px、E. coli 168.7 px）。したがって「grid が定義文字列のために中央を確保する」という監査の仮説は誤りです。
- **主な原因:**
  - MRC は depth 入力がなくても、`_precomputed_depth_tracks=[]` を各 record に渡します（`diagram.py:2897` の `[[] for _ in records]` と :2994）。
  - 受け取る側の `diagram.py:2264` の `precomputed_depth_tracks_provided = _precomputed_depth_tracks is not None` が True になります。
  - その結果 :2273-2279 で show_depth が True になり、preset が 37 px の depth slot を内側に置きます。
  - probe で確かめると、MRC の slot は features, ticks, **depth**, gc_content, gc_skew の順で、単一 record には depth がありません。
  - **成功している図にも影響があります。** HmmtDNA（Homo sapiens）では次のとおりでした。つまり Web 既定の circular 出力は、すべて単一 record と幾何が違います。

| track_type | 空の depth 帯 | GC content の幅（単一 → MRC） |
|---|---|---|
| tuckin | 239.7–268.4 px | 74.1 → 57.4 px |
| middle | 274.5–308.7 px | 74.0 → 68.6 px |

  - depth 入力がないときに `[]` を None に置き換える monkeypatch を当てると、Mycobacterium … BCG と Candidatus Hepatoplasma crinochetorum は Web 既定で成功しました。
- **2 つ目の原因（Salmonella だけに該当）:**
  - tuckin では数値トラックがすべて内側に並びます。
  - 中心の定義文字列は、半径が「最も長い行の幅の半分 + max(8, 0.02×radius)」の円として確保されます（`assemble.py:1143-1196`）。
  - GC content と GC skew の最小幅の合計（67.1 px）を除くと、features と ticks に残るのは [309.0, 386.1] の範囲です。separate strands の features（78 px）と ticks、その内側のラベル帯はここに入りません。
  - 単一キャンバスでも `--separate_strands` だけで失敗します（309.0 px）。
  - `--definition_font_size 16` では失敗し（範囲 103 px）、15 なら成功します。`--plot_title_position top` や `--center_reserved_radius 200` でも成功します。
- 監査の「grid、separate strands、GC、skew のどれかを off にすれば成功する」も、Salmonella では成り立ちません。grid を off にしても separate strands が残るためです。
- **エラーが 'ticks' を名指しする理由:**
  - stack group の失敗では `first_unplaced = intents[-1].slot_id`（`gbdraw/diagrams/circular/radial_layout.py:1192`）が使われます。これはグループの末尾の slot で、実際に置けなかった slot でも、原因の定義帯でもありません。
  - Web は `error_adapter.py:118` でこれを TRACK_LAYOUT/CANNOT_FIT にまとめ、`error-normalization.js:82, 109` の汎用の文に置き換えるので、Python の詳細も消えます。

**3. 修正案（コード）**
- **主な原因（1〜2 行の修正）:**
  - `diagram.py:2994-2995` で、depth 入力がないとき（`cfg.canvas.show_depth` が False）は `_precomputed_depth_tracks=None` と `_precomputed_depth_track_count=None` を渡します。
  - `[]` は「depth 入力はあるが、この record のセルがない」だけを意味する、と :2194 の private 引数の docstring に明記します。
  - 部分的な depth で帯を確保する動作は残します（`tests/test_depth_track.py:1619, 1657, 1697`）。
- **エラーメッセージ:**
  - `radial_layout.py` の 2 つの raise（:935 と :1192-1198）を、1 つのエラー生成関数にまとめます。
  - 実際に置けなかった slot を名指しします。
  - 置ける範囲の内側の限界を定義帯が決めているときは、次の文を付けます。これは既存の `_inside_stack_failure_hint`（:1124）が occupied を受け取るように広げて実装します。「The center definition text reserves {R:.1f}px; shorten the species/record label, reduce the definition font size, set center_reserved_radius, or place tracks outside.」
  - `error_adapter.py` の :118 より前に、定義帯用の template（TRACK_LAYOUT/DEFINITION_RESERVED）を置きます。
  - JS の catalog に 1 行追加します。文言は Web のラベル名（index.html:3294 の "Center Reserved Radius" など）に合わせ、W3 の X-01 と調整します。
- **B が承認された場合:**
  - circular の DefinitionDrawer に最大行幅を渡します。
  - `_definition_reserved_radius_px` と実際の描画が、同じ折り返し結果を使うようにします。
  - 配置できず（CANNOT_FIT）、原因が定義帯のときだけ、1 回配置し直します。
- 追加しないもの: track を外側へ自動で移すこと、新しいオプション。

**4. 上位の対応**
- MRC が record ごとの組み立てに渡す値は、「入力がない（None）」と「入力はあるが空（[]）」を区別する約束にします。
- 1 record の MRC と単一 record で、幾何と凡例が一致することを不変条件にします（TR-01 と共通）。

**5. テスト**
- `tests/test_circular_multi_canvas.py`:
  - depth 入力なしの MRC で、`canvas._gbdraw_track_slot_geometry["records"][*]["slots"]` に slotId "depth" がないこと。
  - 1 record の MRC と単一 record で、gc_content と gc_skew の widthPx と radiusFactor が一致すること。今回の不具合はこの assert で捕まります。
- `tests/test_depth_track.py` の sparse 系のテストが引き続き通ること。
- CLI:
  - Mycobacterium の名前で `--multi_record_canvas --track_type tuckin --gc --skew --separate_strands` が成功すること。
  - Salmonella の名前では、メッセージに "center definition text reserves" が含まれ、adapter が DEFINITION_RESERVED を返すこと。
- Web では Playwright を 1 本追加し、既定設定で Mycobacterium の record を Generate できることを確かめます。
- `tests/reference_outputs/` の再生成は不要です。
- Gallery の **Vnig_TUMSAT-TG-2018** は MRC で depth がないため幾何が変わります。generator から作り直し、Gallery の publication と capture の期待値も更新します。

**6. 規模と PR**
- Python の数行と、エラー生成関数 +25 行、adapter +1 行、JS +1 行、テスト +80 行です。
- 1 PR にまとめ、TR-01 より先に出します。MRC の出力がすべて変わるので Review REQUIRED です。
- B は、決定後に別の PR にします。

**7. 依存とリスク**
- 021b4c49 は origin/main にも入っています。
- 修正すると Web 既定の circular 出力が全体に変わります。Gallery とチュートリアルの画像を作り直す必要があるか確認が要ります。
- browser wheel の再生成が必要です。

---

## FE-07（P2）/translation のない CDS で GTG/TTG が M にならない

**1. 分類**
- IMPLEMENT_EXISTING_AUTHORITY です。根拠は決定的な科学の規則です。5' 端が完全な CDS で、遺伝暗号表の開始コドンから始まる場合、最初のアミノ酸は M になります（INSDC の /translation と Biopython の `cds=True` の規則）。
- LOSATP と orthogroup の出力が変わるので、preflight を記録し、Review REQUIRED とします。

**2. 根本原因の再確認**
- 再現しました。監査の fixture `single.gb` の TESTA_0006（GTG、table 11）で、次の 2 つがどちらも `VRSKRQLSRT` を返しました。
  - `_extract_amino_acid_sequence`（`gbdraw/web_support/feature_metadata.py:290`）
  - `_translate_cds_feature`（`gbdraw/analysis/protein_colinearity.py:2767`）
- したがって、CLI と Web の LOSATP の protein 入力、および orthogroup も影響を受けます（監査では「コードからの判断」でしたが、関数を実行して確認しました）。
- 規模は次のとおりです。MG1655 の codon_start=1 の CDS 4,300 本のうち、開始コドンは GTG 338、TTG 80、ATT 4、CTG 2 本でした。そのまま訳した結果が /translation と違うのは 427 本（約 10%）です。GFF3 の入力では、これらの N 末端がすべて V、L、I になります。
- 比較への影響は N 末端の 1 残基です。/translation のある GenBank と GFF3 を混ぜた比較や、閾値に近い判定で差が出ます。
- **監査にない隣の不具合（実験で確認）:** GFF3 の phase が読まれません。
  - `partial=true;start_range=.,1`、phase=1 の CDS で試しました。
  - metadata は「3 で割り切れない」として翻訳を飛ばします。
  - LOSATP は読み枠がずれた `HETRV` を入力に使います。
  - `gbdraw/io/genome.py:147-157` は phase を qualifiers に残していますが、どちらの翻訳も codon_start しか見ていません。

**3. 修正案（コード）**
- `gbdraw/core/sequence.py` に翻訳の関数を 1 つ追加します。
  - 遺伝暗号表は transl_table で、なければ 1 にします。
  - 読み枠は codon_start を使います。なければ、GFF3 の 5' 端のパーツの phase に 1 を足した値を使います。
  - 5' 端が部分的かどうかは、次のどちらかで判定します。
    - fuzzy な 5' 位置（+ 鎖は start の BeforePosition、− 鎖は end の AfterPosition）
    - GFF3 の `start_range`（+ 鎖）か `end_range`（− 鎖）
  - 5' 端が完全で、読み枠が 1 で、最初のコドンが `CodonTable.unambiguous_dna_by_id[table].start_codons` に含まれるなら、最初の残基を M にします。それ以外はそのまま訳します（to_stop=False）。
- 次の 2 か所をこの関数に置き換えます。
  - `feature_metadata.py:262-295` の翻訳部分
  - `protein_colinearity.py:2761-2775` の本体
- `protein_colinearity.py:2651-2669` の `_translation_table` と `_codon_start` は、ほかで使われていないので削除します。
- どの CDS を含めるか除くかの規則は、呼び出し側ごとに今のまま残します（metadata の pseudo/fuzzy の除外、LOSATP の内部 stop の除外）。
- 追加しないもの: `/transl_except`（Sec/Pyl）への対応、`cds=True` の厳格なモード（終止コドンのない CDS を落としてしまう）。

**4. 上位の対応**
- 「CDS のアミノ酸配列は core の 1 関数だけが作る」を約束にします。/translation がある場合は、これまでどおりそれを優先します。

**5. テスト**
- `tests/test_web_feature_metadata.py`:
  - table 11 の GTG と TTG、table 1 の TTG が M で始まること。
  - 5' 端が部分的な CDS と codon_start=2 の CDS はそのまま訳すこと。
  - GFF3 の phase=1 で正しい読み枠になること。
- `tests/test_protein_colinearity.py`:
  - GFF3 の GTG 開始の CDS で、LOSATP の入力が M で始まること。
  - 正解の基準として、MG1655 の codon_start=1 で transl_except のない CDS について、新しい関数の結果が /translation と一致することを確かめます。数件を抜き出す fast 版と、全件の slow 版を用意します。
- `tests/reference_outputs/` の再生成は不要です（protein 比較はありません）。

**6. 規模と PR**
- Python は +45/−25 行、テストは +80 行程度です。
- 1 PR にまとめます。Review REQUIRED です。

**7. 依存とリスク**
- Web の protein cache のキーは `queryProteinSetHash` を含みます（`run-analysis.js:225-262`）。配列が変われば別のキーになる見込みです。
- ただし `tryPromoteLegacyProteinEntry`（`run-analysis.js:3516`）の昇格の経路が、配列を見ずに古い結果を昇格しないかは、実装時に確かめる必要があります（EVIDENCE）。
- rich-v1 の Session は、/translation と違う amino_acid_sequence を保存します（`session_io.py:232-237, 303-308`）。古い Session を読み込むと、Generate するまでは古い値が表示されます。保存形式は変わらないので、移行は不要です。
- Gallery の orthogroup の Session が /translation のない入力を含むかは、実装時に確認します。

---

## FE-08（P3）大文字小文字だけが違う feature にラベル編集が広がる

**1. 分類**
- IMPLEMENT_EXISTING_AUTHORITY で、**直すのは JS 側**です。
- 根拠は次のとおりです。
  - `docs/REFERENCE/web-app.md:735`「Color rules and Label TSV selectors use case-insensitive Python regular expressions」
  - `tests/test_label_overrides.py:232, :288`
  - Python の照合はすべて IGNORECASE です（`gbdraw/features/visibility.py:217`、`gbdraw/features/colors.py:79`、`gbdraw/labels/filtering.py:292, :328`）。
- Python 側を大文字小文字を区別するように変えると、公開されている契約と、CLI 利用者の表が壊れます。
- JS が `(?-i:^orfA$)` を出す案も検討しましたが、採りません。TSV に新しい書き方が増え、手書きの表の意味とずれるためです。

**2. 根本原因の再確認**
- `makeSelectorUniquenessKey`（`gbdraw/web/js/app/feature-selector.js:219-221`、:236、:253）は、qualifier だけを小文字にし、value はそのまま使います。
- record_id と feature_type は Python 側も完全一致で比べる（`gbdraw/features/selector_values.py:257-262`）ので、ずれているのは value だけです。
- この一意性の判定は共有されているので、ラベルだけでなく、色と表示の "This feature only" も同じ条件で広がります（コードからの判断）。

**3. 修正案（コード）**
- value の比較キーを 1 つ定義し、`makeSelectorUniquenessKey` の value の部分だけに使います。キーは `normalizeSelectorText(value).toUpperCase().toLowerCase()` です。
  - これは Python の IGNORECASE より広く同じとみなしますが、広すぎても locus_tag か hash に切り替わるだけなので安全です。
  - 出力する値と regex は変えません。
- bulk の編集（`app/feature-editor/label-override-table.js:88-99`、:124 の `*\t*\tlabel\t^text$`）も、Python では大文字小文字を無視して当たります。件数表示のグループ分けに同じキーを使います（コードからの判断、要確認）。

**4. 上位の対応**
- 不要です。「JS の一意性の判定は、Python の照合と同じ同値関係を使う」ことを、テストで固定します。

**5. テスト**
- `tests/web/feature-selector.test.mjs`: gene が orfA と ORFA の 2 つの feature で、orfA の selector が locus_tag になること。
- label-override の行が `locus_tag ^CASE1_0001$` になること。
- pytest で、その行を `labels/filtering` に通すと 1 つの feature にだけ当たること。

**6. 規模と PR**
- JS は +4/−2 行、テストは +40 行です。Web のサイズ判定は CLEAR の範囲です。

**7. 依存とリスク**
- Session を保存し直すと、行が gene から locus_tag に変わることがあります。意味は同じか、より正確になります。

---

## PV-05（P2）Web の PDF が px を pt として扱う

**1. 分類**
- IMPLEMENT_EXISTING_AUTHORITY です。根拠は次の 4 つです。
  - CSS の単位の定義（1 in = 96 px = 72 pt）
  - 同じファイルの PNG が `dpi / 96` を使っていること（`gbdraw/web/js/services/export.js:325`）
  - `docs/REFERENCE/output-formats-and-export.md:59-61`
  - CLI（CairoSVG）の PDF と揃うこと
- Web の PDF の物理的な大きさを既存の約束とみなすかどうかは、Owner に確認します。根拠となる記述は見つかっていません。

**2. 根本原因の再確認**
- `export.js:382-396` は `unit:'pt', format:[w,h]` とし、`doc.svg` にも同じ w、h を渡しています。
- 監査の証拠では、Web は 1491.15 pt、CairoSVG は 1118.36 pt（×0.75）、PNG は 72 DPI で 1118 px でした。

**3. 修正案（コード）**
- `CSS_PX_PER_INCH = 96` と `PT_PER_CSS_PX = 72/96` を定義し、PNG の処理と共有します。
- PDF の format と `doc.svg` の width と height に 0.75 を掛けます。
- viewBox がない場合は、`prepareSvgForPdf` で `0 0 w h` を補います。
- jsPDF の `px_scaling` hotfix には頼りません。

**4. 上位の対応**
- 不要です。

**5. テスト**
- Playwright で、MediaBox が SVG の幅と高さ × 0.75（±0.01）になることを確かめます。
- `tests/reference_outputs/` の再生成は不要です。

**6. 規模と PR**
- JS は +6/−3 行です。PV-06 と同じ PR にします。

**7. 依存とリスク**
- Web で作る PDF は 75% の大きさになり、CLI と同じになります。CHANGELOG に書きます。

---

## PV-06（P3）PDF のテキスト層でスペースが失われる

**1. 分類**
- IMPLEMENT_EXISTING_AUTHORITY です。PDF のテキストが正しいことは、既存のテスト `tests/web/gui-audit-regressions.playwright.spec.js:201` でも期待されています。

**2. 根本原因の再確認（コードから）**
- `flattenTextPathsForPdf`（`export.js:163-209`）は、1 文字ごとに `<text>` を作ります。空白だけの `<text>` は、svg2pdf の既定の空白処理で捨てられます。
- 同梱の svg2pdf は textPath に対応していませんが、`xml:space` 属性と CSS の `white-space` は読みます。
- 文字の列挙は `Array.from`（コードポイント単位）で、位置の取得は `getStartPositionOfChar(i)`（UTF-16 単位）です。このため BMP 外の文字があると、それ以降の位置がずれ、末尾の文字も落ちます。ブラウザでの単位は実装時に確かめます。

**3. 修正案（コード）**
- UTF-16 の index で `getNumberOfChars()` まで繰り返します。サロゲートペアは 1 文字にまとめ、次の位置を飛ばします。
- 空白の文字も出力し、`xml:space="preserve"`（`setAttributeNS(XML_NS)`）を付けます。
- 新しい仕組みは作りません。

**4. 上位の対応**
- 不要です。

**5. テスト**
- HmmtDNA_basic_circular の PDF を `readPdfText`（`tests/web/helpers/pdf-text.cjs`）で読み、"cytochrome c oxidase subunit I" と "1 kbp" が含まれることを確かめます。
- BMP 外の文字を含む曲線ラベルで、位置が崩れないことも確かめます。

**6. 規模と PR**
- JS は +10/−5 行で、PV-05 と 1 PR にします。変更は `export.js` だけで、Web のサイズ判定は CLEAR の範囲です。

**7. 依存とリスク**
- PDF のテキストオブジェクトが、空白の分だけ増えます。

---

## GUI の外で見つかったもの（利用者向けの検証は W3 の担当。ここでは core のどこで防ぐべきかだけを示す）

- **`--block_stroke_width -1` の traceback:**
  - 描画の途中で `gbdraw/configurators/legend.py:73, 75` が ValueError を投げ、CLI がこれを捕まえないために traceback になります。ValidationError ではないためです。
  - 防ぐ場所は config model です（`gbdraw/config/models/objects.py`）。既存の :203 と同じ形で、stroke の幅を「有限で 0 以上」に検証します。そうすれば `apply_config_overrides` で、CLI、API、Web、Session がすべて同じエラーになります。
- **負のフォントサイズや stroke 幅がそのまま受け付けられる:**
  - 上と同じく config model で防ぎます（`gbdraw/config/models/objects.py` と `gbdraw/config/models/labels.py`）。
- **`-w 0` と `-s 0` で GC トラックが空になる:**
  - `gbdraw/analysis/skew.py:65, 94` と `gbdraw/analysis/gc.py:91` は、空の表を返すだけの低いレベルの処理です。
  - 防ぐ場所は、単一 record と MRC の両方が通る `_resolve_circular_window_step`（`gbdraw/api/diagram.py:653`）です。window と step に、既存の `_validate_positive_optional`（:737）を適用します。Linear の同じ箇所にも適用します。
- **Dinucleotide に `G` を指定すると IndexError:**
  - `gbdraw/analysis/gc.py:86` と `gbdraw/analysis/skew.py:87` が、`nt_list[1]` を無条件に参照しています。
  - 防ぐ場所は、GcContentConfigurator と GcSkewConfigurator を作る option の境界と、`_slot_dinucleotide`（`assemble.py:1122-1126`）です。`_slot_dinucleotide` は今、2 文字未満の値を黙って既定値に戻しています。
  - 検証の関数は 1 つにして共有します。`XY` のような塩基でない文字を拒否するかどうかは W3 が決めます。

---

## 横断的な提案

1. 「1 record の MRC と単一 record の経路で、幾何、凡例、定義が一致する」ことを、parametrize した 1 つの pytest で不変条件にします。PV-08 と TR-01 は、どちらもこの一致が破れた問題（LSP 違反）です。
2. MRC が自分で行うのは配置だけにします。凡例、depth の有無、slot の実際の構成は、単一 record の経路の関数から得ます。
3. BLAST 表は、`io/comparisons.py` の 1 つの読み込み関数と 1 つの正規化関数だけで扱います。`names=COMPARISON_COLUMNS` が再び増えないよう、静的テストで禁止します。
4. 「入力がない（None）」と「入力はあるが空（[]）」を区別します。021b4c49 の回帰は、これを混同したことが原因でした。
5. 入力の検証は、config model か option の境界で行います。描画の途中で ValueError を出すことをやめ、CLI の traceback を防ぎます。
6. JS が Python の照合を先回りして判定する箇所は、Python と同じ同値関係に合わせます。JS と Python の間の契約テストを 1 本置きます。
7. 新しいエラー（列数、定義帯）は、Python のメッセージ、error_adapter の template、JS の catalog の 3 つを同じ PR で揃えます。W3 の X-01 と調整します。
8. 科学的出力が変わる PR（CO-05、FE-07）では preflight を記録します。Gallery の出力は generator が作るものなので、影響がある Session だけを作り直します。
9. PR の順序は次のとおりです。PV-08-B と CO-05-B は Owner の決定後に進めます。
   1. PV-08（主な原因とメッセージ）
   2. TR-01
   3. CO-05（A で先に出せる）
   4. FE-07
   5. Web（FE-08、PV-05 と PV-06）

## 監査の分類誤り、事実の誤り、範囲の不足

- **PV-08:** 原因の仮説が誤りです。主な原因は MRC の空の depth slot で、2 つ目の原因は tuckin と separate strands の組み合わせで定義の円が入りきらないことです。「CLI の既定では成功する」も不正確で、単一キャンバスでも `--separate_strands` だけで Salmonella は失敗します。影響範囲は監査の記述より広く、Web 既定の circular 出力すべてに及びます。
- **CO-05:** conservation が対処済みという記述は誤りです。ファイルを飛ばすことで比較の位置がずれる不具合も記載されていません。
- **TR-01:** Region annotation の `legend_label` と depth slot の `legend_label` の欠落が記載されていません。
- **FE-07:** LOSATP への影響は、実行して確認しました。GFF3 の phase が読まれない不具合が記載されていません。
- **FE-08:** 影響はラベルに限らず、色と表示の "This feature only" にも及びます。どちらの側を直すかは docs で決まっていて、直すのは JS 側です。
- **`--block_stroke_width -1`:** 直接の原因は `legend.py` の ValueError で、根本は config の検証がないことです。
