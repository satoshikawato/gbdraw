# Issue #288 Record Selector Implementation Plan

- 作成日: 2026-07-15
- 対象 Issue: [#288 Record id input should be a select box](https://github.com/satoshikawato/gbdraw/issues/288)
- 対象: Web UI の Linear input にある record selector
- 状態: 実装済み (2026-07-16)
- 方針: 既存の selector/session contract を維持し、SOLID、KISS、DRY を満たす最小の Web UI 改修とする

## 1. 結論

実現可能性は高い。CLI、Python loader、canonical session request はすでに record ID と
`#index` selector を扱えるため、diagram の assembly/rendering contract は変更しない。

実装の中心は次の三点とする。

1. upload 済み GenBank、または GFF3 と組み合わせた FASTA から record ID と配列長を列挙する。
2. 現在の自由入力を native `<select>` に置き換え、`record id (X,XXX bp)` と表示する。
3. record 一覧を session に保存せず、file から再生成できる一時的な UI state として管理する。

ユーザーが選んだ値は既存の `seq.region_record_id` にそのまま保存する。session version、
canonical request schema、CLI option、Python public API は変更しない。

## 2. 現状と変更境界

### 2.1 現状

- `gbdraw/web/index.html` の Linear input は `seq.region_record_id` を text input で編集する。
- `gbdraw/web/js/app/run-analysis.js` は同じ値を record selector として既存 CLI args に渡す。
- `gbdraw/web/js/services/session-request.js` は `record ID` と `#index` を canonical selector に変換する。
- `gbdraw/web/js/app/python-helpers.js` の `list_genbank_records()` は、Circular record order 用に
  `selector`、`record_id`、`record_length` をすでに返している。
- Python の `parse_record_selector()` / `select_record()` は、record ID と `#index` の双方を
  validation 付きで処理する。

### 2.2 変更する範囲

- Linear input の record selector UI。
- GenBank/FASTA record 一覧の読み込みと正規化。
- file upload、file 差し替え、mode 切り替え、Pyodide 初期化、session 復元に伴う一覧更新。
- 上記の Web unit/browser/packaging test。

### 2.3 変更しない範囲

- `seq.region_record_id` の field 名と文字列型。
- `--record_id`、`--region`、`#index` の CLI semantics。
- region start/end と reverse-complement の適用順。
- session version と canonical request schema。
- Circular record order のユーザー向け挙動。
- SVG geometry、rendering、reference output。

## 3. UI 仕様

### 3.1 表示

ラベルを `Record ID (optional)` から `Record (optional)` に変更する。record ID が重複した
場合には内部値として `#index` を使うため、`Record` の方が実際の意味に合う。

通常の option は次の形式で表示する。

```text
NC_000913.3 (4,641,652 bp)
```

- 長さは省略せず、`en-US` の桁区切りを使った整数 bp とする。
- `kbp` / `Mbp` には丸めない。record の識別と region 座標の判断に正確な長さが必要なためである。
- 長い record ID は値を改変せず表示する。

### 3.2 option の値

option の表示と保存値を分ける。

| 条件 | 表示例 | `seq.region_record_id` に保存する値 |
|---|---|---|
| ID が file 内で一意 | `RecA (1,234 bp)` | `RecA` |
| ID が重複 | `RecA (1,234 bp) [#1]` | `#1` |
| ID が重複 | `RecA (1,100 bp) [#2]` | `#2` |

一意な ID は record の並び替えに強いため ID を保存する。重複 ID は Python loader が曖昧として
拒否するため、その場合だけ既存の `#index` contract を使用する。

### 3.3 selector 未指定

先頭に空値の option を置く。

```text
Automatic (no explicit selector)
```

空値は現在と同じく selector を指定しないことを意味する。upload 時に先頭 record を自動選択
しない。自動選択すると、単一の multi-record file を selector なしで読み込む既存動作が変わる
可能性があるためである。

### 3.4 loading、error、stale selection

- file 未選択時は select を disabled にし、`Upload a sequence file first` を表示する。
- Pyodide 初期化中または record 読み込み中は disabled にし、`Loading records...` を表示する。
- 読み込み失敗時は短い inline error を表示し、詳細は console に残す。
- session 復元値または file 差し替え前の値が新しい一覧に存在しない場合、値を黙って削除しない。
  synthetic option として `RecA (not found in current file)` を表示し、ユーザーに再選択を促す。
- record 変更時に region start/end を自動消去または自動補正しない。既存 validation を維持する。

## 4. Target architecture

```text
file upload / session restore / mode change / Pyodide ready
                            |
                            v
                 linear selector watcher
                            |
                            v
              focused record-selector controller
                            |
                            v
                shared record discovery helper
                            |
                            v
          Python list_sequence_records(path, format)
                            |
                            v
       transient options keyed by linear sequence uid
                            |
                            v
     native <select> --v-model--> seq.region_record_id
                            |
                            v
       existing run-analysis and session-request paths
```

### 4.1 Python helper contract

`list_genbank_records()` を format-aware な内部 helper に一般化する。

```python
list_sequence_records(path, format)
```

対応 format は今回必要な二つだけに限定する。

- `genbank`: Linear GenBank と Circular record order に使用する。
- `fasta`: GFF3+FASTA input の selector 候補に使用する。

戻り値は現在の最小 contract を維持する。

```json
{
  "records": [
    {
      "selector": "#1",
      "record_id": "NC_000913.3",
      "record_length": 4641652
    }
  ]
}
```

GFF3 mode では paired FASTA の record ID と長さを候補として表示する。GFF3 と FASTA の整合性は
既存の authoritative loader が generation 時に引き続き検証する。この UI のためだけに別の
GFF parser や record matching engine は作らない。

### 4.2 Web record discovery

小さい共通 module `gbdraw/web/js/app/record-discovery.js` を追加し、次だけを所有させる。

- file を一意な Pyodide FS path に stage する。
- `list_sequence_records()` を呼ぶ。
- payload を `{selector, recordId, recordLength}` の配列へ正規化する。
- `finally` で一時 file を削除する。

Circular と Linear はこの helper を共有する。format、file、temporary path を渡すだけにし、
mode-specific state や UI 文言は持たせない。

### 4.3 Linear selector controller

`gbdraw/web/js/app/linear-record-selector.js` を追加し、Linear UI 固有の責務を持たせる。

- `linearSeqs` の `uid` ごとの loading/ready/error state。
- 一意 ID と重複 ID の selector value 決定。
- exact bp label の生成。
- current value が一覧にない場合の synthetic option。
- inactive uid の state 削除。
- 非同期 refresh の generation token 管理。

state の概念形は次のとおりとする。

```javascript
{
  [uid]: {
    status: 'idle' | 'loading' | 'ready' | 'error',
    records: [],
    error: ''
  }
}
```

この state は derived/transient data であり、`linearSeqs`、history snapshot、session JSON、
canonical request に追加しない。

### 4.4 Watcher

watch 対象は値全体の deep watch ではなく、record 一覧に影響する最小の dependency とする。

- `mode`
- `lInputType`
- `pyodideReady`
- 各 `linearSeqs` の `uid`
- GenBank mode の `seq.gb` file reference
- GFF3 mode の `seq.fasta` file reference

refresh は KISS を優先して sequence row ごとに直列実行する。各 row は `uid` を含む一意な
temporary path を使う。refresh generation が変わった場合、完了済みの古い結果を state に
commit しない。

## 5. SOLID、KISS、DRY の適用

### 5.1 SOLID

- SRP: record discovery、Linear selector state、template rendering、diagram generation を分離する。
  `run-analysis.js` に新しい file parsing/UI state を直接追加しない。
- OCP: 既存 selector/session contract を変更せず、record 候補の供給を追加する。将来 format を
  増やす場合も `list_sequence_records()` の明示的 format map を拡張できる。
- LSP: `seq.region_record_id` を受け取る既存 run/session path から見た入力 contract を変えない。
- ISP: record discovery helper に app 全体の state object を渡さない。必要な file、format、
  Pyodide access、temporary path だけを渡す。
- DIP: Linear selector controller は Pyodide global を直接探索せず、composition root から渡された
  record reader に依存する。専用 class/interface hierarchy は作らない。

### 5.2 KISS

- native `<select>` を使い、searchable combobox、virtual list、外部 dependency を追加しない。
- record option を session に保存せず、file から再生成する。
- selector の新しい型や schema を作らず、既存の string contract を使う。
- file record 数が多い場合の検索 UI は、実際の利用上の問題が確認されるまで追加しない。
- region 座標の auto-clamp、record 自動選択、GFF の独自 parser を実装しない。

### 5.3 DRY

- Circular 専用の `list_genbank_records()` と新しい Linear parser を並存させず、
  `list_sequence_records()` に統合する。
- record payload normalization と temporary-file lifecycle は `record-discovery.js` の一箇所に置く。
- option label/value の決定は pure helper にし、template、watcher、test に同じロジックを書かない。
- CLI/session の既存 selector normalization を Web UI 側に再実装しない。UI は既存 contract に合う
  `record ID` または `#index` を生成するだけにする。

## 6. 実装フェーズ

### Phase 0: Characterization test

実装前に次の非回帰 contract を test で固定する。

1. 空の `region_record_id` は selector 未指定として canonical request に変換される。
2. `RecA` は `recordId` selector、`#2` は `recordIndex: 1` に変換される。
3. session round-trip で `region_record_id` が保持される。
4. record selector と region start/end/reverse の既存 CLI args が変わらない。
5. Circular record order の既存 record list と position merge が変わらない。

完了条件:

- selector 値の生成を変えずに test が通る。
- SVG reference update を必要としない。

### Phase 1: Shared record discovery

1. Python helper を `list_sequence_records(path, format)` に一般化する。
2. GenBank と FASTA のみを明示的に許可し、unknown format は structured error にする。
3. 空 file、parse error、0 records、0-length record を扱う。
4. `record-discovery.js` に staging、Python call、payload normalization、cleanup を実装する。
5. Circular record order を shared helper の GenBank path へ移す。

完了条件:

- Circular の record list、order、row assignment が改修前と一致する。
- GenBank と FASTA で ID、index、length が取得できる。
- temporary file が成功時と失敗時の双方で残らない。

### Phase 2: Linear selector state

1. `linear-record-selector.js` に transient state と pure option builder を追加する。
2. record ID の出現回数を一度だけ集計し、一意 ID と重複 ID の value を決める。
3. `uid` を identity として、sequence row の reorder 後も正しい options を維持する。
4. 削除された row の state を消す。
5. generation token と file reference check で stale async result を拒否する。
6. unmatched current value を synthetic option と warning に変換する。

完了条件:

- file 差し替えを連続実行しても、最後に選んだ file の options だけが表示される。
- row reorder/remove/add 後に options が別 row へ移動しない。
- session restore value を勝手に空にしない。

### Phase 3: Template and lifecycle wiring

1. text input を native `<select>` に置き換える。
2. loading、disabled、error、unmatched warning を表示する。
3. watcher に最小 dependency を追加する。
4. upload、clear、mode change、Pyodide ready、session restore から同じ refresh path を使う。
5. select 操作を既存 history/undo contract に合わせる。必要なら select change 全体を一つの
   undoable action とするが、record options 自体は history に含めない。

完了条件:

- GenBank と GFF3+FASTA の両 mode で record list が表示される。
- 選択後の Generate、LOSAT、feature metadata extraction が同じ selector を使う。
- save/load 後も選択値が復元され、options は embedded file から再生成される。

### Phase 4: Verification and documentation

1. unit、packaging、session、CLI selector test を実行する。
2. browser で multi-record GenBank upload と session restore を確認する。
3. offline mode で外部 request が増えていないことを確認する。
4. Issue #288 の acceptance criteria と実装結果を照合する。
5. user-facing release note が必要な release 単位なら一行追加する。

## 7. Test plan

### 7.1 Web unit tests

新規 `tests/web/record-selector.test.mjs` で pure behavior を検証する。

- exact bp formatting: `4641652 -> 4,641,652 bp`。
- 一意 ID は ID value を使う。
- 重複 ID は `#index` value と suffix を使う。
- record 順を変えても一意 ID value は同じ。
- blank selector option。
- unmatched current value の synthetic option。
- invalid/missing length の fallback 表示。
- state purge と stale generation reject。

### 7.2 Python/Web packaging tests

- embedded helper に `list_sequence_records()` が含まれる。
- helper が GenBank/FASTA format を受理する。
- template が text input ではなく select を使用する。
- Web module が package/wheel に含まれる。
- CSP または CDN dependency の変更がない。

### 7.3 Session and CLI tests

- unique record ID の save/load round-trip。
- duplicate ID 用 `#index` の canonical request round-trip。
- selector 付き whole-record render。
- selector と region を組み合わせた render。
- reverse complement と selector の組み合わせ。
- legacy session の unmatched selector を load しても値を失わない。

### 7.4 Browser tests

最低限、次の user flow を実ブラウザで確認する。

1. 2 records の GenBank を Linear row に upload する。
2. `RecA (1,234 bp)` と `RecB (567 bp)` が表示される。
3. `RecB` を選択し、region を指定せず Generate する。
4. `RecB` のみが対象になる。
5. region start/end を追加し、同じ record の region が描画される。
6. session を保存して読み直し、`RecB` が選択状態に戻る。
7. 同じ ID を持つ 2 records では `[#1]` / `[#2]` が表示され、選択できる。
8. GFF3+FASTA mode で FASTA record IDs が表示される。

Playwright 検証前に Node と Python の両 installation path を確認する。Node の
`@playwright/test` がなければ Python Playwright で同じ flow を検証する。

### 7.5 Suggested commands

```bash
node tests/web/record-selector.test.mjs
python -m pytest tests/test_web_packaging.py -k "record or python_helpers or index" -q
python -m pytest tests/test_session_io.py tests/test_linear_selectors.py -q
python -m pytest tests/ -q -m "not slow"
```

Web source の変更で SVG geometry は変わらないため、reference output は更新しない。比較 test が
差分を報告した場合は、selector の対象 record が意図せず変わっていないかを先に調査する。

## 8. Acceptance criteria

すべて満たしたとき Issue #288 を完了とする。

1. Linear の record selector が text input ではなく select である。
2. 各 option が `record id (X,XXX bp)` 形式で表示される。
3. GenBank と GFF3+FASTA の両 input mode で候補が読み込まれる。
4. selector 未指定の既存動作が変わらない。
5. 一意 ID は record ID、重複 ID は `#index` で正しく選択できる。
6. file 差し替え、row reorder/remove、mode change で stale options が残らない。
7. session save/load と canonical request round-trip で選択値が失われない。
8. unmatched legacy value を黙って削除しない。
9. Circular record order に回帰がない。
10. offline SPA のままで、build step、CDN dependency、server API を追加しない。
11. fast test、targeted Web test、targeted browser test が通る。
12. reference SVG に意図しない差分がない。

## 9. Risks and mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| 古い非同期 parse 結果が新しい file state を上書きする | 誤った record を選択する | generation token、file reference check、一意 temporary path |
| duplicate record ID | ID selector が曖昧で generation が失敗する | duplicate のみ `#index` を保存し、label に index を表示 |
| session の selector が file と一致しない | select が空に見える、値を失う | synthetic option と inline warning、値は保持 |
| Pyodide 初期化前に session/file が設定される | options が読み込まれない | `pyodideReady` を watcher dependency に含める |
| GFF3 と FASTA の record 集合が異なる | 無効な候補を選べる | 候補は FASTA から表示し、既存 GFF/FASTA loader を最終 authority として error を明示 |
| 多数 record の native select が使いにくい | 選択に時間がかかる | 今回は native select。実測上必要になった場合だけ検索 UI を別 Issue にする |
| shared helper 化で Circular が回帰する | multi-record canvas の order が壊れる | Phase 0 characterization と targeted browser test を先に追加 |

## 10. Non-goals and follow-ups

今回実装しないもの:

- searchable combobox、typeahead、virtual scrolling。
- record description、organism、topology 等の追加 metadata 表示。
- selected record length に基づく region end の auto-fill、auto-clamp。
- record 選択時の organism/strain/subtitle 自動入力。
- multi-record file を複数 Linear rows に自動展開する機能。
- selector/schema の typed object 化や session version bump。
- CLI の default record-loading semantics の変更。

native select の実用性、record 数、GFF/FASTA mismatch の発生頻度を実装後に確認し、必要なら
それぞれ独立した Issue として扱う。

## 11. Definition of done

- acceptance criteria をすべて満たす。
- record discovery と selector option logic がそれぞれ単一 owner にある。
- `run-analysis.js` に file parsing と selector UI state を追加していない。
- session/canonical request に derived record options を保存していない。
- 新規 abstraction が GenBank/FASTA record discovery 以外の責務を持たない。
- public API、CLI、session schema、SVG output の互換性が維持される。
- 実行した test と未実行 test、その理由を PR に記録する。

## 12. Implementation result

- `list_sequence_records(path, format)` と shared Web record discovery を追加し、GenBank と FASTA を
  同じ contract で列挙するようにした。
- Linear selector は UID ごとの transient state と generation token を使い、file 差し替え、row
  reorder/remove、mode 切り替え、session 復元時の stale result を拒否する。
- native `<select>`、exact bp label、duplicate ID の `#index`、unmatched value の synthetic option、
  loading/error state を実装した。
- Circular record order も shared discovery を使用し、既存の position merge contract を維持した。
- Web unit 17 files、fast Python suite (`1411 passed, 16 skipped`)、targeted packaging/build、Python
  Playwright の GenBank/GFF3+FASTA/session flow を実行し、すべて成功した。reference SVG の更新はない。
