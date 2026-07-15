# Python Session Compatibility Matrix

- 調査日: 2026-07-15
- 対象: session version 27～30、GUI-authored / CLI sidecar、Circular / Linear
- 結論: internal replay は維持できるが、public typed request への lossless conversion は未達
- 関連 ADR: [`ADR_PYTHON_SESSION_API.md`](ADR_PYTHON_SESSION_API.md)

## 1. 判定基準

この文書では、次の三つを別の互換性として扱う。

| 判定 | 意味 |
|---|---|
| Envelope | `format`、`version`、`files` を `validate_session()` が受理する |
| Internal replay | `session_to_cli_args()` と現行 CLI owner が再実行できる |
| Typed conversion | CLI option 名、argument index、parser module に依存せず、公開 request へ情報を欠落なく変換できる |

`SUPPORTED_SESSION_VERSIONS` に含まれることは Envelope の保証にすぎない。version 27～30 の
すべてについて Typed conversion が可能である、という意味ではない。

## 2. version history

| Version | schema change | Python/Web の現状 |
|---:|---|---|
| 27 | Python の session validation/materialization と `cliInvocation` schema 1 を導入した基準 version | current reader が受理する。実 history shape を使った Python replay fixture はない |
| 28 | optional `editorState` とその Web normalization を追加 | Web は欠落時に default を補う。Python は editor state を typed input に変換しない |
| 29 | optional `losatDerivedCache` を追加 | Web は保存・復元する。Python の public request model には cache owner がない |
| 30 | feature visibility を manual rules と per-feature overrides に分離 | Web は legacy rules を分割する。Python は visibility edit state を typed input に変換しない |

version 30 の導入後にも、`config.cliOptions` と `cliInvocation` からの
`multiRecordPositions` 補完が version bump なしで追加されている。このため `version == 30`
だけでは document shape を一意に決められず、将来の typed payload 追加には version 31 への
bump が必要である。

## 3. version × origin × mode matrix

記号:

- **CLI**: `cliInvocation.schema == 1`、string `args`、完全な `fileBindings` / `files.cliTables`
  がある場合、stored CLI invocation を internal replay できる。
- **GUI-partial**: `config` / `files` から現行 `_gui_session_to_cli_args()` が扱う subset を
  best-effort で再構成できる。lossless ではない。
- **No typed**: CLI 非依存の lossless typed conversion を証明できない。

| Version | GUI Circular | GUI Linear | CLI Circular | CLI Linear |
|---:|---|---|---|---|
| 27 | GUI-partial / No typed | GUI-partial / No typed | CLI / No typed | CLI / No typed |
| 28 | GUI-partial / No typed | GUI-partial / No typed | CLI / No typed | CLI / No typed |
| 29 | GUI-partial / No typed | GUI-partial / No typed | CLI / No typed | CLI / No typed |
| 30 | GUI-partial / No typed | GUI-partial / No typed | CLI / No typed | CLI / No typed |

この表の CLI 判定は version から導いたものではなく、`cliInvocation` の payload が完全である
場合の条件付き判定である。`cliInvocation` が欠落・不完全なら同じ version でも GUI-partial
または明示 error になる。

## 4. capability matrix

| Capability | GUI `config` / `files` fallback | CLI sidecar | Typed conversion blocker |
|---|---|---|---|
| GenBank | Circular/Linear とも materialize 可能 | binding により replay | typed record source と lifetime がない |
| GFF3 + FASTA | Circular/Linear とも materialize 可能 | binding により replay | pair identity と lifetime がない |
| records table | 専用復元なし | `files.cliTables` と dependency rewrite | row model、embedded dependency の typed owner がない |
| selector / region / reverse | Linear row field の安全な subset を復元 | stored args で replay | per-record input spec がない |
| multi-record position | GUI config の current fieldは保持されるが Python fallback は CLI arg にしない | stored args。Web は不足時に `cliInvocation` から補完 | version 30 内にも shape 差がある |
| default / feature color table | GUI fallback に専用復元なし | binding/CLI table で replay | `files` と `config.colors` の precedence が未定義 |
| visibility / shape / label edit state | Web state には存在するが Python fallback は再現しない | CLI に表現された範囲だけ replay | editor state と render request の境界が未定義 |
| blacklist / whitelist / qualifier priority / label override | Python fallback に専用復元なし | stored args と binding で replay | GUI text/rule state と file input の正規形がない |
| depth | basic Circular/Linear file と主要 option を復元 | binding で replay | multi-track metadata と custom slot を完全には表現しない |
| conservation | Circular BLAST file の基本 subset を復元 | stored args で replay | FASTA-derived data、labels/colors/table の完全な復元がない |
| custom circular/linear track slots | Python fallback は slot array を再現しない | stored args で replay | renderer-specific typed slot model がない |
| nucleotide comparison | Linear GUI の基本 BLAST path/threshold subset を復元 | stored args で replay | table/dependency と comparison source の統一 model がない |
| protein / orthogroup / collinearity | Linear GUI の主要 option subset を復元 | stored args と source session で replay | generated/precomputed input と cache metadata の owner がない |
| plot title / legend / common layout | 多くを復元 | stored args で replay | GUI-only の一部 advanced state は欠落する |
| output formats / interactive metadata | GUI fallback は override がなければ SVG に固定 | args と `renderFormats` に保持 | saved result、format request、interactive policy の typed 分離がない |
| `losatCache` / `losatDerivedCache` | CLI argument listだけでは表現されない | Linear internal owner が source session を直接参照 | cache schema、identity、optional lifetime が public contract にない |
| saved SVG `results` | artifact として保持するが replay input にはしない | sidecar に保存 | render request と previous result を区別する型がない |

## 5. 実証できた範囲と test gap

- current version の synthetic tests は、GUI Circular、GUI Linear、CLI binding、CLI table
  dependency、Circular save/replay、Linear protein/orthogroup option を cover する。
- version 29 の repository fixture では、GUI-only Linear session と CLI Circular sidecar の
  `session_to_cli_args()` smoke conversion が成功した。
- version 27 は Envelope acceptance test があるが、当時の GUI/CLI document shape を固定した
  replay fixture はない。
- version 28 の historical replay fixture はない。
- version 29 の smoke conversion は typed losslessness を証明しない。
- version 30 の後発 field が version bump なしで追加されているため、current synthetic shape を
  version 27～29 に relabel する test は migration evidence として採用しない。

## 6. E0 decision

1. request model の owner は `gbdraw.api.requests` とする。ただし E1～E5 完了までは
   `gbdraw.api.__all__` から export しない。
2. version 27～30 の `cliInvocation` replay は CLI/GUI の internal compatibility path として
   維持する。CLI argument string から public request を作ることは保証しない。
3. GUI-only legacy session は current best-effort replay を維持するが、lossless public
   conversion の対象にしない。unsupported field を推測で default 化しない。
4. 新しい session schema は CLI 非依存の canonical typed render request payload を持つ。
   payload の mode、record inputs、options、comparisons、output policy を schema で検証する。
5. canonical payload 追加時は session version を 31 に上げ、Python と Web の version、writer、
   migration test を同じ change で更新する。
6. temporary materialization は context manager が所有し、typed request 内の path はその
   lifetime 外へ持ち出せない契約にする。
7. public session symbol は version 31 canonical payload の load → materialize → request →
   render → save round-trip が両 mode で通った release でのみ追加する。

したがって Workstream E の public bridge gate は閉じたままである。次に進める作業は、
session replay の公開ではなく、session から独立した request model と version 31 canonical
payload の schema 設計である。

[Validation plan](PYTHON_API_VALIDATION_PLAN.md) | [Follow-up plan](PYTHON_API_FOLLOWUP_PLAN.md) | [Python API](PYTHON_API.md)
