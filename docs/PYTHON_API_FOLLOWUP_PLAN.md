# Python API Follow-up Plan

- 作成日: 2026-07-14
- 対象: `docs/PYTHON_API_IMPROVEMENT_PLAN.md` の「13. 次の作業」
- 対象バージョン: release readiness は `0.14.0b0`、session bridge は設計ゲート通過後のリリース
- 状態: Workstream A～D、E0、E1 完了。Workstream E2 以降は canonical request 設計を継続
- 目的: Python API 改修後に残ったテスト安全性、リリース説明、lint debt、option 監査、session 公開境界を、既存契約を壊さず段階的に完了する

## 0. 結論

後続作業は次の順序で進める。

| 順序 | 作業 | 完了時の成果 |
|---:|---|---|
| 1 | reference SVG 生成の opt-in 化 | 通常の full suite が tracked fixture を更新しない |
| 2 | `0.14.0b0` release note | behavior correction、追加 API、session 非公開判断を利用者へ説明する |
| 3 | Ruff debt の baseline 化と解消 | production Python の lint が 0 件になり、CI が blocking gate になる |
| 4 | `DiagramOptions` 監査 | 70 field の owner、mode、test、重複候補が一覧化される |
| 5 | session bridge の設計・実装 | CLI 引数文字列を公開契約にせず、検証・再実行・保存を公開 API だけで行える |

依存関係は次のとおりである。

```text
reference safety
    ├── reliable full test evidence
    │       ├── release readiness
    │       ├── Ruff cleanup
    │       └── DiagramOptions audit
    └── fixture update workflow

DiagramOptions audit
    └── typed request model ADR
            └── session validation/materialization
                    └── session conversion/render/save
```

reference safety、release note、Ruff、`DiagramOptions` 監査を
`0.14.0b0` の release-readiness milestone とする。session bridge は現行 ADR の
再開条件を満たしてから実施し、`0.14.0b0` のリリースをブロックしない。

### 0.1 実施結果

2026-07-14 に A～D、2026-07-15 に validation follow-up と E0 を完了した。E は公開実装へ
進めず、canonical request model の設計・実装へ移行した。

| Workstream | 状態 | 結果 |
|---|---|---|
| A | 完了 | reference 生成 13 case を明示 opt-in にし、比較失敗時の actual SVG を `tmp_path` へ隔離した |
| B | 完了 | [`RELEASE_NOTES_0.14.0b0.md`](RELEASE_NOTES_0.14.0b0.md) を追加し、`DOCS.md` からリンクした |
| C | 完了 | Ruff 0.15.12 の 85 diagnostics を 0 にし、CI lint を blocking にした |
| D | 完了 | [70-field audit](DIAGRAM_OPTIONS_AUDIT.md) と non-default forwarding test を追加した |
| E | E0 完了 / public gated | [compatibility matrix](PYTHON_SESSION_COMPATIBILITY_MATRIX.md) で legacy lossless conversion 不可を確認し、version 31 canonical request を次の境界とした |

監査では 70 field すべてに owner があり、dead field はなかった。分類は shared 39、
circular-only 10、linear-only 21 で、shared のうち 2 field は互換 alias である。
この結果から、`DiagramOptions` は `0.14.0b0` では分割しない。

検証結果:

- Python API targeted suite: `547 passed`
- full non-slow suite: `1343 passed, 16 skipped, 6 deselected`
- reference/output comparison: `13 passed, 13 skipped`
- `ruff check gbdraw/`: 成功
- reference SVG: 13 件の SHA-256 が着手時から不変
- public Markdown local target check: 成功

## 1. 着手時の基準値

2026-07-14 の監査値は次のとおりである。

| 項目 | 現状 |
|---|---|
| Python API targeted suite | 308 passed |
| `gbdraw.api.__all__` | 87 symbol |
| `DiagramOptions` | 70 field |
| `CircularMultiRecordOptions` | 5 field |
| reference 生成ケース | Circular 8、Linear 5、合計 13 |
| Ruff | `ruff 0.15.12`、85 diagnostics |
| Ruff 内訳 | `W291` 22、`W293` 55、`F401` 6、`F841` 2 |
| session schema | current 30、supported 27～30 |
| session replay owner | `gbdraw.session_io.session_to_cli_args` |
| session 公開状態 | `gbdraw.api` からは非公開 |

`DiagramOptions` の 70 field はすべて `gbdraw/api/diagram.py` から少なくとも
1 回参照される。したがって、監査では単純な文字列検索で dead field と判定せず、最終 consumer
まで値が届くかを確認する。

Ruff の既存 85 件は production Python に限った値である。現在の共通コマンドに含まれる
`W503` は Ruff 0.15.12 では有効な rule selector ではない。また CI は
`|| true` と `continue-on-error` の両方を持つため、現状の lint job は
失敗を release gate にできていない。

## 2. 共通原則

### 2.1 既存作業と改行コードを混ぜない

- 実装開始時に `git status --short` と reference SVG の差分を記録する。
- 既存の未コミット変更を変更対象へ混ぜない。
- CRLF/LF の全ファイル置換を機能差分や lint cleanup と同じ PR に含めない。
- 広範囲の機械修正前後で `git diff --ignore-cr-at-eol` も確認する。

### 2.2 公開契約を先に固定する

- `gbdraw.api.__all__`、signature、dataclass field/default の snapshot を維持する。
- behavior correction 以外では既定 SVG と CLI behavior を変更しない。
- option の存在ではなく、正規化後の値が最終 owner へ届くことを test する。
- public symbol を追加する PR では、同じ PR で契約 snapshot と文書を更新する。

### 2.3 薄い公開層を維持する

- parsing、normalization、I/O、rendering の owner を共有し、`gbdraw.api` は
  re-export または小さな adapter に限定する。
- CLI parser、`argparse.Namespace`、CLI option 名、argument index を Python API の
  型へ持ち込まない。
- option class や request model は、caller と実装の合計複雑度が下がる場合だけ追加する。

### 2.4 検証結果を汚染しない

- reference の生成と比較を別操作にする。
- full suite 実行後に `tests/reference_outputs/` の差分を必ず確認する。
- binary export、Playwright、LOSAT の環境依存失敗は、機能回帰と分けて報告する。
- generated browser wheel と `dist/` は追跡対象へ追加しない。

## 3. Workstream A: reference SVG 生成を明示 opt-in にする

状態: 完了。

### 3.1 問題

`tests/test_output_comparison.py::TestGenerateReferences` は通常の pytest collection に
含まれ、13 個の tracked SVG を直接上書きできる。比較 test より先に生成 test が動くと、
意図しない geometry change を新しい正解として受け入れてしまう。

さらに `TestOutputComparison` は比較失敗時の `*.actual.svg` を
`tests/reference_outputs/` に書くため、生成を止めても通常 test が tracked directory を
変更する余地がある。

### 3.2 採用する契約

1. reference 更新には pytest option `--update-reference-outputs` を必須にする。
2. 生成 class に `reference_generation` marker を付ける。
3. option なしでは生成 test を skip し、write path へ到達させない。
4. 比較 test の actual output は `tmp_path` 配下へ保存し、tracked directory を
   読み取り専用の fixture 集合として扱う。
5. reference 更新は atomic replacement にし、途中失敗で既存 fixture を壊さない。
6. CI の通常 job から reference 更新 option を渡さない。

環境変数だけによる opt-in は採用しない。shell profile や CI secret に残った値で意図せず
生成できるためである。

### 3.3 変更対象

- `tests/conftest.py`
  - `pytest_addoption` で `--update-reference-outputs` を登録する。
  - marker と option を確認する collection hook または fixture を追加する。
- `tests/test_output_comparison.py`
  - `TestGenerateReferences` を marker 対象にする。
  - option なしの write を禁止する。
  - temporary file を同じ directory に作ってから `Path.replace()` する。
  - comparison failure の actual path を `tmp_path` へ移す。
- `pyproject.toml`
  - `reference_generation` marker を `--strict-markers` 用に登録する。
- `AGENTS.md`、`CLAUDE.md`
  - comparison-only と明示的 update のコマンドを分けて記載する。
- `.github/workflows/test.yml`
  - fast test 後に reference directory が clean であることを確認する。

### 3.4 実装手順

1. 現在の tracked reference 13 件について path と SHA-256 を記録する。
2. option なしで generator class を実行し、全 case が skip される characterization test を
   追加する。
3. comparison class を実行し、reference directory の hash が不変であることを確認する。
4. actual output の保存先を temporary directory へ移す。
5. explicit option 付き生成を 1 case で実行し、その case だけ更新されることを確認する。
6. 13 case 全体の生成は geometry 更新が承認された場合だけ行う。
7. CI に `git diff --exit-code -- tests/reference_outputs/` を追加する。失敗時にも
   差分を観察できるよう、この確認 step は `if: always()` で実行する。

### 3.5 検証コマンド

```bash
# 通常の比較。tracked reference は変更しない。
python -m pytest \
  tests/test_output_comparison.py::TestOutputComparison \
  tests/test_output_comparison.py::TestQuickValidation \
  -q

# opt-in がないため全 case が skip されることを確認する。
python -m pytest \
  tests/test_output_comparison.py::TestGenerateReferences \
  -q

# geometry 更新を意図した場合だけ実行する。
python -m pytest \
  tests/test_output_comparison.py::TestGenerateReferences \
  --update-reference-outputs \
  -q

git diff --exit-code -- tests/reference_outputs/
```

### 3.6 完了条件

- 通常の full non-slow suite が reference SVG を上書きしない。
- comparison failure でも `tests/reference_outputs/` に新規 file を作らない。
- option なしでは生成 path へ到達できない。
- explicit update の使用方法と review 手順が repository guidance に記載される。
- 更新後の reference を同じ comparison test が読み、期待どおり pass する。

## 4. Workstream B: `0.14.0b0` release note

状態: 完了。

### 4.1 配置

repository に changelog の継続運用がまだないため、今回だけで新しい release 管理制度を
作らず、versioned document `docs/RELEASE_NOTES_0.14.0b0.md` を canonical source とする。
`docs/DOCS.md` の Reference からリンクし、GitHub Release を作る場合は同文書を
要約して使用する。

### 4.2 必須内容

#### Behavior corrections

- `DiagramOptions.collinearity_anchor_mode` の `all`、
  `one_to_one`、`rbh` と alias が正規化後の collinearity owner まで届く。
- 既定値 `rbh` は変わらない。
- `save_figure_to` は明示要求した形式を生成できない場合に
  `ExportError` を送出する。
- 戻り値には関数終了時に実在する path だけを含める。
- CLI の optional converter 不在時の warning-and-skip behavior は維持する。

#### Added public capabilities

- `CircularMultiRecordOptions` と `build_circular_multi_diagram`。
- `InteractiveSvgContext`、`build_interactive_svg_context`、
  `enrich_svg`。
- records、conservation、circular-track table の model/reader。
- label whitelist、qualifier priority、label override reader。
- `DiagramOptions` の label DataFrame/file fields。
- `render_to_bytes(..., interactive_context=...)`。

#### Compatibility and migration

- 既存 `assemble_*` import path と既定値は維持する。
- binary export caller は `gbdraw.exceptions.ExportError` または
  `GbdrawError` を捕捉する。
- DataFrame と file path を同時指定した場合は `ValidationError` になる。
- session replay はこの release では公開 Python API ではなく、CLI/GUI 内部機能のままである。
- session 非公開判断から公開 bridge の再開条件へ ADR のリンクを付ける。

### 4.3 実装手順

1. `8f952d4..67cb66d` の public contract diff を release note の added/changed
   一覧と照合する。
2. behavior correction ごとに「以前」「現在」「移行方法」を 3～6 行で記述する。
3. strict export の短い before/after example を追加する。
4. `docs/PYTHON_API.md` と矛盾する capability 名がないことを確認する。
5. version、project metadata、CLI help の `0.14.0b0` 表記を照合する。
6. 文書内 code block を docs test の既存 execution model で実行するか、実行対象外なら
   擬似コードであることを明記する。

### 4.4 完了条件

- behavior correction、added API、互換性、session 境界が 1 文書から把握できる。
- `save_figure_to` の既存 caller が必要な migration を判断できる。
- release note の symbol 名が `gbdraw.api` から実際に import できる。
- `docs/DOCS.md` から到達できる。
- release note に未実装の session bridge を利用可能と書かない。

## 5. Workstream C: Ruff debt を 0 にして CI gate にする

状態: 完了。

### 5.1 固定 baseline

`ruff 0.15.12` と次の設定相当で 85 件である。

```text
W291  22  trailing whitespace
W293  55  blank line contains whitespace
F401   6  unused import
F841   2  unused local variable
```

file 別内訳:

| Rule | Files |
|---|---|
| `W291` | `gbdraw/cli.py` 5、`gbdraw/crop_genbank.py` 17 |
| `W293` | `gbdraw/cli.py` 3、`gbdraw/crop_genbank.py` 18、`features/objects.py` 2、`features/tracks.py` 26、`svg/circular_features.py` 6 |
| `F401` | `crop_genbank.py` 2、`features/factory.py` 1、`features/objects.py` 1、`render/groups/linear/length_bar.py` 1、`svg/text_path.py` 1 |
| `F841` | `gbdraw/layout/linear.py` 2 |

### 5.2 Tool contract の修正

1. `ruff==0.15.12` を development dependency に追加する。
2. `pyproject.toml` に次を置き、CLI option の重複をなくす。

```toml
[tool.ruff]
required-version = "==0.15.12"

[tool.ruff.lint]
select = ["E", "F", "W"]
ignore = ["E501"]
```

3. `AGENTS.md`、`CLAUDE.md`、CI のコマンドを
   `ruff check gbdraw/` に統一する。
4. `W503` を command line から削除する。
5. cleanup が 0 件になるまでは既存 CI behavior を一時維持できるが、command 自体は
   必ず実行して diagnostics を表示する。

version 更新は通常 dependency update と同様に独立 PR で行い、Ruff を毎回 unpinned latest
へ上げない。

### 5.3 Cleanup の分割

#### PR C1: whitespace 77 件

- `W291` と `W293` だけを safe fix する。
- source 全体の format、quote、import sort は同時に変更しない。
- 改行コード変更が混ざっていないことを確認する。
- SVG/reference output と public contract が不変であることを確認する。

#### PR C2: unused import 6 件

- import が side effect、registration、type-only reference、public re-export を担っていないか
  1 件ずつ確認する。
- `__all__` による public export なら単純削除せず、意図を明示する。
- optional dependency の import timing を変えない。

#### PR C3: unused variable 2 件

- 代入右辺に副作用がある場合は呼び出し自体を残す。
- tuple unpacking や callback contract のため必要なら、意味のある underscore 名へ変える。
- layout output の snapshot/reference test を実行する。

#### PR C4: blocking CI

- `ruff check gbdraw/` が 0 件であることを確認する。
- workflow の `|| true` を削除する。
- lint step/job の `continue-on-error` を削除する。
- test source や Web source の新しい lint scope はこの PR に含めない。

### 5.4 検証

各 cleanup PR で次を行う。

```bash
ruff check gbdraw/
python -m pytest tests/test_public_contract.py -q
python -m pytest tests/test_output_comparison.py::TestOutputComparison -q
python -m pytest tests/ -v -m "not slow"
git diff --exit-code -- tests/reference_outputs/
```

### 5.5 完了条件

- Ruff version と rule selection が local/CI/documentation で一致する。
- production Python の `E`、`F`、`W`
  （`E501` を除く）が 0 件になる。
- CI は新しい violation で失敗する。
- broad formatter、import sorter、line-ending rewrite が cleanup に混ざらない。
- reference SVG と public contract に差分がない。

## 6. Workstream D: `DiagramOptions` 70 field の監査

状態: 完了。結果は [`DIAGRAM_OPTIONS_AUDIT.md`](DIAGRAM_OPTIONS_AUDIT.md) に記録した。

### 6.1 目的

目的は dataclass を分割することではなく、各 field の実利用を証明し、重複、mode 固有、
互換用、未到達の候補を識別することである。監査 PR では field の削除や rename を行わない。

### 6.2 成果物

`docs/DIAGRAM_OPTIONS_AUDIT.md` を作り、70 field すべてについて次を記録する。

| 列 | 内容 |
|---|---|
| Field/default | public contract 上の型と既定値 |
| Category | config、feature、label、scalar track、comparison、layout、output 等 |
| Mode | Circular single、Circular multi、Linear |
| High-level reader | どの `build_*` が読むか |
| Forwarding path | adapter → `assemble_*` の引数 |
| Final owner | parser/normalizer/configurator/renderer の実 consumer |
| Non-default test | 指定値が owner まで届く test |
| Documentation | recipe または reference へのリンク |
| Classification | live、mode-specific、compatibility、alias、duplicate candidate、dead candidate |
| Action | keep、document、test、deprecate proposal、group proposal |

### 6.3 監査手順

1. `dataclasses.fields(DiagramOptions)` を source of truth に field list を生成する。
2. AST/`rg` で `options.<field>` の直接参照を機械収集する。
3. 各参照を `build_circular_diagram`、`build_circular_multi_diagram`、
   `build_linear_diagram` から最終 owner まで手動で追跡する。
4. assembler が引数を受けるだけで使っていない場合は live と判定しない。
5. alias/normalization を持つ field は、alias、canonical value、default の三つを test する。
6. DataFrame/file pair は入力形態の違いとして記録し、同時指定時の
   `ValidationError` test を確認する。
7. singular/plural depth field は Circular/Linear の意味と precedence を明記する。
8. docs にない live field と、docs にあるが owner へ届かない field を別々に列挙する。
9. `tests/fixtures/public_contract.json` との差分がないことを確認する。

### 6.4 判定ルール

- **keep**: 複数 mode で共有され、現在の bundle が caller code を減らす。
- **mode-specific**: 1 mode だけで live。直ちに移動せず、誤用時の validation/documentation を
  先に検討する。
- **duplicate candidate**: 同一意味・同一 precedence・同一 owner であることを test で
  証明できる場合だけ候補にする。
- **dead candidate**: public builder、low-level assembler、docs、session migration のいずれでも
  意味を持たない場合に限る。
- **deprecation candidate**: replacement、warning period、major release までの維持方法を
  提示できる場合だけ提案する。

### 6.5 分割提案の gate

`DiagramOptions` の分割は次をすべて満たす場合だけ別計画にする。

1. 既存 `DiagramOptions` を互換層として維持できる。
2. caller と implementation の合計 code が増えない。
3. Circular/Linear の boolean mode switch を増やさない。
4. option forwarding test の重複が減る。
5. public contract、既定 SVG、CLI behavior を維持できる。

条件を満たさない場合、監査結果は「70 field を維持し、mode/capability 表を文書化する」を
正式な結論とする。分割しない判断も完了結果である。

### 6.6 完了条件

- 70 field すべてに final owner または未到達の証拠がある。
- non-default forwarding が未検証の field は test backlog として明示される。
- field を削除・分割せずに監査文書を review できる。
- 変更提案には net code量、互換性、deprecation path がある。

## 7. Workstream E: CLI 非依存の session bridge

状態: E0、E1、request render foundation 完了、public bridge は gated。version 27～30 から
lossless typed request を復元できないため公開 symbol は追加せず、version 31 canonical
request の前提実装を進める。

### 7.1 現状と設計ゲート

現行 session は version 30、supported version は 27～30 で、次の二つの形を持つ。

1. GUI-authored session: `config`、`ui`、`files`、
   `results`、optional `losatCache`。
2. CLI sidecar: 上記に加え `cliInvocation.args`、
   `fileBindings`、`renderFormats` を replay の authoritative data として持つ。

現行 `session_to_cli_args()` は embedded file を temporary directory に展開し、
CLI argument string を返す。この関数をそのまま re-export しない。

session bridge は、`Circular/Linear` 共通の typed input request と pure conversion が
設計できた場合だけ再開する。次のいずれかが残る場合は ADR を更新して停止する。

- public model が CLI option 名または argument index を保持する。
- `gbdraw.api` が `gbdraw.circular`、`gbdraw.linear`、
  parser helper を import する。
- supported legacy session を lossless に変換できる証拠がない。
- temporary file の cleanup owner が不明である。

### 7.2 Phase E0: compatibility matrix と ADR 更新

状態: 完了。結果は
[`PYTHON_SESSION_COMPATIBILITY_MATRIX.md`](PYTHON_SESSION_COMPATIBILITY_MATRIX.md) と
[`ADR_PYTHON_SESSION_API.md`](ADR_PYTHON_SESSION_API.md) に記録した。

session version 27～30 について、GUI/CLI sidecar × Circular/Linear の matrix を作る。

各 cell で次を調べる。

- GenBank、GFF3+FASTA、records table。
- selector、region、reverse、multi-record position。
- color、visibility、shape、label table。
- depth、conservation、custom track slot。
- nucleotide comparison、protein comparison、collinearity。
- output formats、interactive metadata、plot title、legend。
- LOSAT cache、record metadata、saved SVG result。

`config` と `files` だけで typed request を再現できるかを確認し、
`cliInvocation.args` だけに存在する情報を gap として列挙する。

ADR には次を決定事項として追記する。

1. request model の owner module。
2. legacy session の対応範囲。
3. canonical typed request payload を session schema に追加するか。
4. schema version を上げる条件と Web/Python の migration policy。
5. public symbol を追加する release boundary。

lossless conversion が不可能なら、session public replay を実装せず、必要な canonical payload
を新しい session schema に追加する設計を先に行う。

### 7.3 Phase E1: typed request model

状態: internal incubation model として完了。実装と次の normalization/rendering phase は
[`PYTHON_SESSION_CANONICAL_REQUEST_PLAN.md`](PYTHON_SESSION_CANONICAL_REQUEST_PLAN.md) に記録した。

名前は ADR で確定するが、責務は次のように分ける。

- **Record input spec**
  - GenBank または GFF3+FASTA source。
  - selector、region、reverse、label、subtitle、grid position。
- **Circular render request**
  - record inputs、`DiagramOptions`、optional
    `CircularMultiRecordOptions`。
- **Linear render request**
  - record inputs、comparison inputs、`DiagramOptions`。
- **Output request**
  - format、output directory/prefix、overwrite、interactive context policy。
- **Render result**
  - 実際に生成した path/bytes、mode、warnings、feature/match metadata、
    optional LOSAT cache。

Circular と Linear は別 dataclass にし、巨大な mode switch を作らない。共通 field は
小さな immutable value object で共有する。既存 table model と reader を再利用し、
同じ validation を再実装しない。

この phase の完了条件:

- in-memory caller と file caller の両方を表現できる。
- CLI parser を使わずに既存 high-level builder を呼べる。
- request → render の pure/side-effect boundary が明確である。
- request model 単体で session に依存しない。
- `DiagramOptions` 監査で mode-specific と判定した field の扱いが決まっている。

### 7.4 Phase E2: validated document と materialization

候補 public contract:

```python
document = load_session_document(path)

with materialize_session(document) as materialized:
    request = session_to_request(materialized)
```

必要な型:

- validated `SessionDocument`。
- context manager である `MaterializedSession`。
- materialized path、source document、warnings、optional cache metadata。

必要な例外:

- format/envelope error。
- unsupported/future version error。
- embedded file/base64/size/path error。
- replay conversion error。

これらは既存 caller が `GbdrawError` または `ValidationError` で
まとめて捕捉できる hierarchy にする。raw `JSONDecodeError`、
`binascii.Error`、`OSError` を公開境界から漏らさない。

materialization の要件:

- basename sanitization と directory containment を維持する。
- declared size と decoded size を検証する。
- depth codec schema と integer overflow を検証する。
- context exit 後に全 temporary file が削除される。
- caller が path を lifetime 外へ保持できないことを文書化する。
- duplicate filename、partial failure、cleanup failure を test する。

### 7.5 Phase E3: session → typed request conversion

conversion は mode ごとの pure function とし、CLI arg list を返さない。

```python
request = circular_request_from_session(materialized)
# または
request = linear_request_from_session(materialized)
```

実施順:

1. current GUI Circular。
2. current GUI Linear。
3. canonical typed payload を持つ CLI sidecar。
4. legacy version 27～30 migration。
5. unsupported legacy combination の明示 error。

同じ session に GUI state、canonical request、`cliInvocation` が併存する場合の
precedence を ADR で固定する。推奨 precedence は canonical typed request、GUI config、
legacy CLI invocation の順だが、legacy CLI invocation を public conversion に使用するには
CLI 名を contract にしない migration owner が必要である。

conversion test は、現行 CLI replay の最終 SVG group/metadata と typed request replay を
比較する。ただし CLI args 自体の一致を期待値にしない。

### 7.6 Phase E4: render bridge

候補 contract:

```python
with open_session(path) as session:
    result = render_session(
        session,
        output_directory=output_dir,
        output_prefix="restored",
        formats=("svg", "interactive_svg"),
    )
```

result は次を返す。

- 実際に生成した path。
- optional in-memory bytes。
- mode と normalized request。
- warnings。
- feature/match metadata。
- optional LOSAT cache metadata。

`save_figure_to` と同様、存在しない path を成功結果に含めない。部分失敗は
warning で成功扱いせず、型付き exception と cause を返す。render 中の temporary input は
context の lifetime 内に保持する。

### 7.7 Phase E5: session build/save と LOSAT cache

session 保存は CLI args ではなく、typed request と structured render result を入力にする。

必要条件:

- GUI が読める config/files/results を生成する。
- canonical typed request payload を保存する。
- source path をそのまま JSON に埋めず、embedded file と binding を安全に構築する。
- LOSAT cache は optional で、欠落しても再描画できる。
- cache が存在する場合は schema、pair identity、content integrity を検証する。
- atomic JSON write を維持する。
- save → load → request → render の round-trip test がある。

session schema を更新する場合は、Python の
`CURRENT_SESSION_VERSION` と Web の `SESSION_VERSION` を同じ PR で更新し、
旧 version の read compatibility を維持する。

### 7.8 Phase E6: public export と文書

public symbol は E0～E5 の gate が通った後に一括して追加する。

- `gbdraw.api.__all__` と contract fixture を更新する。
- `docs/PYTHON_API.md` に load/inspect/render/save recipe を追加する。
- optional cache、temporary lifetime、version compatibility、error handling を説明する。
- `docs/ADR_PYTHON_SESSION_API.md` の status を superseded または amended にする。
- CLI/GUI tutorial は既存 behavior を維持し、Python API recipe と混同しない。

### 7.9 session test matrix

最低限、次を test する。

| 領域 | Cases |
|---|---|
| Envelope | invalid JSON、wrong format、missing version、future version、unsupported version |
| Embedded files | invalid base64、size mismatch、path traversal、duplicate name、depth codec error |
| Lifetime | success cleanup、conversion failure cleanup、render failure cleanup |
| Modes | Circular single/multi、Linear single/multi |
| Inputs | GenBank、GFF3+FASTA、records table、selector、region、reverse |
| Tracks | depth、conservation、custom slots、label/color/visibility table |
| Comparisons | BLAST、precomputed protein、orthogroup/collinearity metadata |
| Output | SVG、interactive SVG、binary failure、overwrite、multiple formats |
| Schema | GUI/CLI shape、versions 27～30、新 canonical request version |
| Cache | absent、valid、invalid、round-trip |
| Dependency | `gbdraw.api` import が CLI modules を import しない |

### 7.10 session 完了条件

- GUI session と CLI sidecar を public API だけで検証・再実行できる。
- public request/result に CLI arg string、option name、arg index がない。
- temporary file lifetime が context manager で保証される。
- version、embedded file、conversion、render failure を型で識別できる。
- cache なしでも描画でき、cache ありの round-trip が可能である。
- CLI/GUI の現行 replay behavior が変わらない。
- public contract、Python docs、session ADR、Web schema test が一致する。

この条件を満たせない場合は、部分的な replay API を公開せず、ADR に blocker と再開条件を
記録して Workstream E を gated のまま維持する。

## 8. PR 分割

| PR | 内容 | 規模 | 依存 | 主なリスク |
|---|---|---:|---|---|
| A1 | reference generation opt-in と actual output 隔離 | S | なし | fixture 更新手順の混乱 |
| A2 | CI reference cleanliness gate と guidance | S | A1 | CI の false positive |
| B1 | `0.14.0b0` release note と docs link | S | A1 | 実装との差異 |
| C0 | Ruff version/config/command の固定 | S | A1 | local/CI version 不一致 |
| C1 | `W291/W293` 77 件 | S | C0 | 改行差分の混入 |
| C2 | `F401/F841` 8 件 | S | C1 | side effect の誤削除 |
| C3 | Ruff blocking CI | S | C2 | 未確認 scope の gate 化 |
| D1 | 70-field inventory と owner trace | M | A1 | syntactic use の誤判定 |
| D2 | forwarding test/documentation gap | M | D1 | test の過剰 mock 化 |
| E0 | compatibility matrix と ADR 更新 | M | D2 | legacy 情報不足 |
| E1 | CLI 非依存 typed request model | M–L | E0 | model の肥大化 |
| E2 | document validation/materialization | M | E1 | security、cleanup |
| E3 | GUI/CLI session conversion | L | E2 | schema/precedence |
| E4 | render bridge と structured result | M–L | E3 | export 部分失敗 |
| E5 | save/cache round-trip、Web schema | L | E4 | Python/Web 乖離 |
| E6 | public export、docs、contract | M | E5 | beta API 固定化 |

A～D は release-readiness series、E は session-bridge series として分ける。E0 の gate で
停止した場合も、A～D の完了と `0.14.0b0` release を妨げない。

## 8.1 次期計画

`DiagramOptions` 監査後の validation と E0 compatibility matrix は
[`PYTHON_API_VALIDATION_PLAN.md`](PYTHON_API_VALIDATION_PLAN.md) で完了した。E0 の結論に従い、
次は [`PYTHON_SESSION_CANONICAL_REQUEST_PLAN.md`](PYTHON_SESSION_CANONICAL_REQUEST_PLAN.md) の
Phase R3 として canonical payload schema 1 の codec を実装する。legacy replay の公開や
version 31 bump は request/schema contract が固まるまで行わない。

## 9. 共通検証ゲート

### 9.1 Python API targeted

```bash
python -m pytest \
  tests/test_api_requests.py \
  tests/test_api_request_render.py \
  tests/test_api_library_usage.py \
  tests/test_public_contract.py \
  tests/test_comparisons.py \
  tests/test_feature_shapes.py \
  tests/test_gc_content_percent.py \
  tests/test_circular_conservation.py \
  tests/test_circular_track_slots.py \
  tests/test_linear_track_slots.py \
  tests/test_depth_track.py \
  tests/test_interactive_svg_cli_format.py \
  tests/test_linear_selectors.py \
  tests/test_cli_tables.py \
  -q -m "not slow"
```

### 9.2 Session targeted

```bash
python -m pytest \
  tests/test_session_io.py \
  tests/test_api_library_usage.py \
  tests/test_public_contract.py \
  tests/test_interactive_svg_cli_format.py \
  -q -m "not slow"
```

### 9.3 Full and lint

```bash
python -m pytest tests/ -v -m "not slow"
ruff check gbdraw/
git diff --exit-code -- tests/reference_outputs/
```

### 9.4 Web/session schema

Web session schema を変更した PR では次も行う。

```bash
node --test \
  tests/web/file-imports.test.mjs \
  tests/web/history-inputs.test.mjs \
  tests/web/history.test.mjs \
  tests/web/session-feature-metadata.test.mjs

npx playwright test tests/web/depth-track-session.playwright.spec.js
python tools/prepare_browser_wheel.py
python -m build
```

Node の `@playwright/test` がない場合は Python Playwright で同等の targeted browser
check を行う。Chromium sandbox error は同じ test を必要な sandbox escalation 付きで再実行する。

## 10. リスクと対策

| リスク | 対策 |
|---|---|
| generator が regression fixture を先に更新する | opt-in option、marker、CI clean check |
| comparison failure が tracked directory を汚す | actual output を `tmp_path` へ隔離 |
| Ruff cleanup が全ファイル format になる | rule 単位 PR、version pin、diff/line-ending review |
| unused import 削除で registration が消える | side effect、`__all__`、optional import を個別監査 |
| option audit が分割ありきになる | owner trace を成果物にし、分割は net code gate 後だけ |
| session model が CLI parser の写像になる | mode-specific typed request、CLI import prohibition test |
| legacy session が情報を失う | version/mode compatibility matrix、lossless でなければ停止 |
| embedded file が path traversal や resource abuse を起こす | containment、size/schema validation、context cleanup |
| Python/Web session version がずれる | version bump と両側 test を同一 PR にする |
| partial export を成功扱いする | structured result は実在 artifact のみ、failure は型付き例外 |
| session bridge が `0.14.0b0` を遅らせる | release-readiness と session milestone を分離 |

## 11. 完了の定義

### 11.1 Release-readiness milestone

状態: 完了。targeted、full non-slow、reference comparison、lint の各 gate は
0.1 の実施結果どおり成功した。

次をすべて満たしたとき A～D を完了とする。

1. 通常 test は reference SVG を生成・更新しない。
2. reference 更新には `--update-reference-outputs` が必要である。
3. `0.14.0b0` release note が behavior correction と migration を説明する。
4. Ruff diagnostics が 0 件で、CI が blocking になる。
5. `DiagramOptions` 70 field の owner/mode/test/action が一覧化される。
6. public contract と既定 SVG に意図しない差分がない。
7. targeted、full non-slow、reference cleanliness、lint gate が通る。

### 11.2 Session-bridge milestone

次をすべて満たしたとき E を完了とする。

1. ADR の再開条件を満たす typed request model がある。
2. GUI/CLI session の supported path を CLI parser なしで request へ変換できる。
3. validation、materialization、render、save の lifetime/error contract が型付きである。
4. LOSAT cache が optional に round-trip する。
5. Python/Web schema、public contract、docs、tests が一致する。
6. unsupported legacy path は silent fallback せず、明示 error または ADR の非対応判断になる。

## 12. 実施 checklist

- [x] 作業ツリーと reference SVG の既存差分を記録した。
- [x] A1 の opt-in protection を他作業より先に実装した。
- [x] `0.14.0b0` release note の実装 diff を再確認した。
- [x] Ruff 85 件の baseline を同じ version で再現した。
- [x] `DiagramOptions` contract snapshot を保存した。
- [x] session ADR の再開条件を review した。
- [x] session version 27～30 の compatibility matrix を作り、E0 gate を判定した。
- [ ] browser wheel、Playwright、CairoSVG の個別確認は Web/session schema を変更する E phase で行う。
