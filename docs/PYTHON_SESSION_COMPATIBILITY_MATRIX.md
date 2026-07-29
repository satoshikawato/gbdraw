# Python Session Compatibility Matrix

- 調査日: 2026-07-15（最終更新: 2026-07-29）
- 対象: session version 27～33、36～39、GUI-authored / CLI sidecar、Circular / Linear
- 結論: version 27～30 は internal replay のみ。version 31～33、36～39 は canonical public typed conversion を提供
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
| 31 | CLI 非依存の typed `renderRequest` と canonical resource map を追加 | Python/Web とも canonical request を復元し、public typed bridge が受理する |
| 32 | materialized annotation set、target、style、track binding を canonical request に追加 | version 31 を migrate し、annotation を typed request で保持する |
| 33 | Linear custom-track geometry を schema 2 へ更新 | version 32 の feature-slot no-op 値を compatibility rule に従って migrate する |
| 36 | compact session-global runtime handle と runtime/display binding の分離を導入 | reader は protein raw schema 4、derived schema 3、identity manifest schema 2 を current artifact として検証する。nucleotide raw schema 2 は変更しない |
| 37 | canonical request schema 4 で重複した assembly-level output prefix を削除 | version 37 writer は `renderRequest.output.prefix` だけを書き、`diagramOptions.output.outputPrefix` は書かない |
| 38 | canonical request schema 5 で明示的な grouping と Circular batch output array を追加 | current writer は Circular `single` / `grid` / `batch`、Linear `single` を保存する。batch は record ごとの output array、それ以外は単一 output object を使う |
| 39 | Web の重複した active/per-mode legend・title-position fields を `ui.layoutPreferences` へ統合 | current writer は一つの per-mode preference tree を保存し、active 値は mode と Circular grouping から計算する。canonical request schema は 5 のまま |

Version 34 と 35 は未マージの開発ブランチ内だけで使われた形式であり、reader の対応対象ではない。
Version 39 の reader は version 27～33、36～38 から直接移行する。Version
36～39 は同じ protein cache/manifest artifact 契約を使う。

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
| 31 | Canonical / Typed | Canonical / Typed | Canonical / Typed | Canonical / Typed |
| 32 | Canonical / Typed | Canonical / Typed | Canonical / Typed | Canonical / Typed |
| 33 | Canonical / Typed | Canonical / Typed | Canonical / Typed | Canonical / Typed |
| 36 | Canonical / Typed | Canonical / Typed | Canonical / Typed | Canonical / Typed |
| 37 | Canonical / Typed | Canonical / Typed | Canonical / Typed | Canonical / Typed |
| 38 | Canonical / Typed | Canonical / Typed | Canonical / Typed | Canonical / Typed |
| 39 | Canonical / Typed | Canonical / Typed | Canonical / Typed | Canonical / Typed |

この表のversion 27～30における CLI 判定は version から導いたものではなく、`cliInvocation` の payload が完全である
場合の条件付き判定である。`cliInvocation` が欠落・不完全なら同じ version でも GUI-partial
または明示 error になる。

## 4. capability matrix

| Capability | GUI `config` / `files` fallback | CLI sidecar | Legacy blocker / canonical status |
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
| protein / orthogroup / collinearity | Linear GUI の主要 option subset を復元 | stored args と source session で replay | version 31～33、36～39 は canonical comparison options/resources、version 36～39 は compact runtime handle と current protein identity/cache artifact を保持 |
| plot title / legend / common layout | 多くを復元 | stored args で replay | GUI-only の一部 advanced state は欠落する |
| output formats / interactive metadata | GUI fallback は override がなければ SVG に固定 | args と `renderFormats` に保持 | saved result、format request、interactive policy の typed 分離がない |
| `losatCache` / `losatDerivedCache` | CLI argument listだけでは表現されない | Linear owner が source session の artifact を参照 | version 36～39 は raw schema 4、derived schema 3、manifest schema 2 を current artifact として検証し、version 27～33 の protein raw schema 2 と derived schema 1 を migration artifact として隔離する |
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
- version 31～33、36～39 は canonical request の load/materialize/typed conversion と save round-trip を
  focused session/API test で固定する。current version 39 writer は `renderRequest.schema == 5` のみを出力する。
- schema 1～3 の `diagramOptions.output.outputPrefix` は compatibility
  reader が受理して無視する。schema 4 と 5 では同 field を拒否し、
  `renderRequest.output.prefix` だけを output policy の authority とする。
- supported session の `cliInvocation.args` に旧 `--show_gc`、
  `--show_skew`、`--suppress_gc`、`--suppress_skew`、`--depth`、
  `--show_depth` が残る場合は、resource materialization 後に canonical
  option へ変換して binding index を再計算する。fresh CLI parser は旧名を
  受理しない。
- 同じ reader-only 境界で、保存済み `--depth_tick_interval`、
  `--feature_table`、`--collinear_max_gene_gap`、Circular の
  `multi_record_size_mode=sqrt`、Linear の `label_placement=on_feature` と
  `track_layout=spreadout|tuckin` を、それぞれ
  `--depth_large_tick_interval`、`--feature_visibility_table`、
  `--collinear_max_unit_gap`、`auto`、`above_feature`、
  `above|below` へ変換する。fresh CLI/API は旧名・旧値を受理しない。
- label の canonical override は Circular の
  `labels.circular.scope` / `labels.circular.placement`、Linear の
  `labels.linear.scope` / `labels.linear.placement` /
  `labels.linear.rotation`、両 mode 共通の rendering policy である
  `labels.rendering` である。保存済みの flat `show_labels` /
  `allow_inner_labels`、`canvas.show_labels`、mode-nested canvas
  `show_labels`、`canvas.circular.allow_inner_labels` は reader が canonical
  leaf へ移行する。fresh config override validator と current writer は旧
  field を受理・出力しない。
- 保存済み Circular slot の `spacing` は schema 1～3 migration 内だけで
  legacy gap として運び、`strict`、`compress`、`reserve` は旧 no-op/
  derived state として除去する。fresh slot parser は四つの key を拒否し、
  gap は `inner_gap_px` と `outer_gap_px` で指定する。
- この reader-only migration が内部で使う
  `__gbdraw_legacy_spacing` transport は canonical request schema 1～3 と
  old-session reader の境界外には出さない。schema 4 と 5 の decoder は同 key を
  拒否し、current writer も出力しない。pixel 単位の legacy spacing は
  `innerGapPx` / `outerGapPx` に変換できるが、radius に対する factor 単位の
  spacing は replay だけを保証する。factor 値は current schema 5 へ無損失に
  再保存できないため、保存前に明示的な inner/outer pixel gap へ移行する。
- typed output の明示的 prefix は `sample.v1` のような dot を保持する。
  Separate-diagram Circular batch の implicit record ID が重複した場合は、
  input order で最初に空いている `_2`, `_3`, ... suffix を付ける。
  `CircularBatchRequest` は explicit `batch` grouping と record ごとの
  `RenderOutputRequest` を保持し、schema 5 session で round-trip する。
- `--annotation-table` は現行 `--annotation_table` の active alias、
  `--gc_content_tick_interval` は現行
  `--gc_content_large_tick_interval` の active alias である。これらは
  reader-only migration や削除対象には含めない。
- version 34 と 35 は branch-internal development format として明示的に拒否する。
- version 27～33 の protein raw schema 2 は引き続き
  `legacyArtifacts.proteinRawCandidates` に隔離する。完全な FASTA、program、arguments、direction、
  feature mapping を検証できた candidate だけを raw schema 4 へ昇格する。derived schema 1 は
  migration evidence として復元する。nucleotide raw entry は schema 2 のまま保持する。
- version 36～39 の generated protein FASTA、raw QUERY/SUBJECT、protein map、derived reference は
  compact runtime handle を使う。完全な feature analysis identity と display metadata は manifest
  schema 2 が一度だけ保持する。**Save Raw LOSAT TSV** は manifest で QUERY/SUBJECT を readable
  alias へ hydrate し、内部 handle を current session の外へ出さない。

## 6. E0 decision と現在の到達点

1. request model の owner は `gbdraw.api.requests` とする。version 31 gate の完了後は
   public session/request symbols を `gbdraw.api` から export する。
2. version 27～30 の `cliInvocation` replay は CLI/GUI の internal compatibility path として
   維持する。CLI argument string から public request を作ることは保証しない。
3. GUI-only legacy session は current best-effort replay を維持するが、lossless public
   conversion の対象にしない。unsupported field を推測で default 化しない。
4. 新しい session schema は CLI 非依存の canonical typed render request payload を持つ。
   payload の mode、record inputs、options、comparisons、output policy を schema で検証する。
5. canonical payload 追加時は session version を 31 に上げ、Python と Web の version、writer、
   migration test を同じ change で更新する。その後も schema に意味のある変更は version 32、33、36～39 で同じ原則を守る。
6. temporary materialization は context manager が所有し、typed request 内の path はその
   lifetime 外へ持ち出せない契約にする。
7. public session symbol は version 31 canonical payload の load → materialize → request →
   render → save round-trip が両 mode で固定されたため公開されている。version 32、33、36～39 もこの typed boundary を維持する。

したがって Workstream E の public bridge gate は version 31 で開放済みである。
current writer は version 39 と canonical request schema 5 を出力し、reader は version 27～33、36～39 を受理する。
version 27～30 は引き続き internal replay のみで、public typed conversion の境界は version 31 のままである。

[Validation plan](PYTHON_API_VALIDATION_PLAN.md) | [Follow-up plan](PYTHON_API_FOLLOWUP_PLAN.md) | [Python API](PYTHON_API.md)
