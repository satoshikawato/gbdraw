# Changelog

All notable changes to gbdraw are documented here, in the style of
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). This project uses
[Semantic Versioning](https://semver.org/); pre-1.0 minor versions
(`0.MINOR.0`) may include breaking changes.

Detailed, per-release notes (migration steps, session/schema compatibility,
and full feature descriptions) live under `docs/RELEASE_NOTES_*.md`. This
file is the short, chronological index; follow the links below for the full
write-up of a release.

## [Unreleased]

Fixes from the 2026-09-30 Web GUI audit of `dev`. The plan and the approved
decisions are in
[`docs/internal/web-gui-audit-20260930/`](./docs/internal/web-gui-audit-20260930/03_IMPLEMENTATION_REFERENCE.md).

- Comparison tables: every BLAST outfmt 6/7 reader (CLI `-b` and
  `--comparisons_table`, Web uploads, Circular similarity rings, and the LOSATP
  parser) now reads the first 12 columns by position and validates their types.
  Tables with extra columns, such as `-outfmt "6 std qlen slen"`, are no longer
  misread. Only lines that start with `#` are comments, so a `#` inside an ID
  no longer truncates the row. The web app reports a malformed table as a
  comparison-input error with its line number instead of an unclassified
  validation error (CO-05).
- **CLI behavior change:** a missing, unreadable, or malformed `-b` file now
  stops the run with a non-zero exit status instead of being skipped, which
  shifted later tables onto the wrong record pair (N-03).
- Linear comparisons now reject a table whose query or subject IDs name the other
  endpoint or another displayed record: the CLI names the conflicting ID, and
  the web app reports a comparison-endpoint error and keeps the previous
  Result. Unknown IDs keep the positional pair (the CLI logs a warning), and
  SVG record-ID metadata always names the endpoint records (CO-06).

- Circular Multi-Record Canvas output (the Web default) no longer reserves an
  empty depth slot when there is no depth input; a one-record canvas now has the
  same slot geometry, legend and center definition as a single-record diagram,
  so every Web-default Circular figure changes (PV-08, N-01).
- A Multi-Record Canvas legend now uses the same builder as a single-record
  diagram: custom slot labels and colors, added skew slots, region annotation
  and depth slot `legend_label` values now appear (TR-01, N-04).
- When the center definition blocks the inside tracks, the species line is
  wrapped at word boundaries and the tracks are placed again. An explicit
  `center_reserved_radius` keeps one line. `definition_font_size` 18 is the
  default, so 18 counts as not explicit and may wrap; any other value is
  explicit and never wraps. The remaining failure names the slot and the reserved definition
  radius, and the Web shows it as `TRACK_LAYOUT` / `DEFINITION_RESERVED`; a
  failure caused by an explicit `center_reserved_radius` names that radius
  instead (`CENTER_RESERVED`) (PV-08, PD-OI-078).
- A specific-color rule whose caption equals a generated legend row, such as a
  rule captioned `rRNA` on tRNA features, keeps its own color row with a hex
  suffix (`rRNA [#ff0000]`) instead of being dropped (N-06).

- Non-pseudo CDS translated without `/translation` now start with `M` when the
  5' end is complete, the reading frame starts at the first base, and the first
  codon is a start codon of `transl_table` (for example `GTG` and `TTG` in
  table 11), as in INSDC `/translation`. This changes **Copy aa FASTA**,
  Interactive SVG feature metadata, and LOSATP and similarity-group protein
  inputs for GFF3 input and for GenBank CDS without `/translation` (FE-07).
- GFF3 CDS phase now sets the reading frame. A CDS with phase 1 or 2 was
  previously translated out of frame, or skipped by the feature popup when its
  length was not a multiple of 3 (N-05).
- LOSATP now skips a CDS whose `codon_start` is invalid instead of translating
  it from the first base. Saved LOSATP rows for a record whose proteins changed
  are not reused; the next search recomputes them.

- The CLI, the Python API, and typed Web or Session requests reject the same
  invalid values and name the option or setting. `--window`, `--step`,
  `--depth_window`, and `--depth_step` take positive integers; before, `0` drew
  empty GC content and skew tracks (X-02).
- A dinucleotide (`-n/--nt`, a track slot's `nt`) is two letters from `A`, `C`,
  `G`, `T`, and `U` in any case. `U` is counted as `T`, so `AU` matches `AT`.
  `-n XY` no longer draws flat tracks and `-n G` no longer raises `IndexError`
  (D-26, N-13).
- Font sizes must be greater than zero and stroke widths zero or greater;
  `--block_stroke_width -1` no longer ends in a traceback. Offsets, spacing,
  `track_axis_gap`, and label rotation keep their current ranges (D-27).
- A Circular definition font set without an interval in a Web or Session
  request, or in Python API `config_overrides`, uses the font size plus 2 as the
  line interval, as the CLI and 0.13.0 do (GE-03).
- Qualifier Priority and whitelist edits apply to Sessions written by the CLI
  or saved from `main`. Label maps compiled by an earlier render are no longer
  kept as preserved settings or reused when a table is attached (SE-08).
- Python validation errors reach the Web with a code, field, reason, and
  Track row, Depth series, table line, or setting path instead of an
  unclassified error (X-01).

- Web validation failures report a recognized code with the field and the
  Sequence, Line, Track row, Depth series, setting, or available band that
  applies, instead of an unknown error. Session import keeps the Worker's
  diagnosis: a file that is not JSON or not a Session names the format problem
  (X-01, SE-09).
- Generate sends numeric settings as typed and no longer replaces a rejected
  value with Auto or a default: a GC window of 0, a decimal step, a non-numeric
  e-value, or a comparison filter outside its range is rejected with the field
  named, and Generate no longer rewrites the settings it reads (X-02: TR-04,
  GE-08, GE-04, CO-09).
- Save before the first Generate explains the missing input with the same check
  Generate uses (IN-08). After loading a Session older than Session 40, Save asks
  for one Generate and offers it; Generate then Save writes the current format
  (D-25, SE-05). The error panel offers **Generate** when that is the correction
  and no **Save Session** after a failed Save (N-12).

- Each Linear record gets a **Definition** inferred from its own `/organism`
  and `/strain`; the first record's value is no longer the file default for all
  records. A definition is the record's own value, then the file default you
  typed, then the inferred value, and **Using file default** marks only a value
  you typed. Choosing another record in a row infers that record's definition.
  Sessions save it as `inferred_definition`; older and CLI Sessions are not
  re-inferred when loaded (IN-06, D-12).
- **Reset Settings** also clears each Linear record's **Definition**,
  **Subtitle**, crop, and reverse complement, and the comparison alignment plan;
  Files, selected records, file defaults, and Depth assignments stay, and
  **Undo** restores all of it, alignment plan included (SE-10, D-15).
- A GenBank record whose ACCESSION or VERSION line is empty (Prokka format) is
  named by its LOCUS instead of `KEYWORDS` in the record choices, the LOSAT FASTA
  extraction, and the match sequences, and an empty ORGANISM line no longer gives
  `Unclassified.` as the organism. The browser reads the header in one place and
  hands any file it cannot read exactly to the Python reader; a file that starts
  with a UTF-8 byte-order mark is rejected as the command line rejects it instead
  of being counted as one record (IN-02, N-02, D-33).
- GFF3 + FASTA record discovery lists only the sequences that have GFF3 rows, in
  FASTA order, like the command line; a FASTA-only sequence is no longer offered
  as a record that Generate cannot draw (IN-03, D-13).
- Replacing a cropped Linear File with a multi-record File expands it into one
  row per record carrying only File-level values, instead of keeping one row
  with the old crop and definition (IN-04).
- A Circular source turned inactive keeps its Multi-Record Canvas record order
  across a switch to Linear and back (IN-05).
- Circular **Region Annotations** offer the records the next Generate draws and
  require a target record whenever it draws several, whether on one canvas or
  as separate diagrams; before, separate diagrams accepted a target the Python
  request could not bind (FE-05).
- **This feature only** always writes the stable feature hash, also when two
  records share a record ID; before, it wrote a rendered ID Python could not
  match. For identical duplicate features the color applies to both (FE-09, D-14).

- History records one Undo step for each checkbox, radio button, or button
  change, also when it is made with its label text or the keyboard, or while a
  text field has focus. Such changes were previously dropped or merged into
  another step, so edits can now produce more Undo steps (SE-02, SE-03, N-18).
- Ctrl+Z, Ctrl+Shift+Z, and Ctrl+Y (Cmd on macOS) now work while a select has
  focus; text fields keep the browser's text undo (SE-04).
- **Undo** and **Redo** are unavailable while a diagram is generating, and the
  header names the reason. They previously reverted the draft behind a running
  Generate. Settings edits remain available (GE-06).
- Undo or Redo of a step such as a feature color that adds a legend entry, or
  **Reset Settings**, no longer empties the feature list after a mode switch;
  these History steps no longer copy the feature metadata (SE-01, N-19, N-20).
- Cancel during Generate preparation no longer stops the loaded diagram
  engine, so the next Generate reuses it (GE-09).
- Undo, Redo, and the rollback of a failed Session Load restore settings
  exactly. Before, a round trip set unset track slot sides, the Features lane
  direction, and both track axis indexes, pinned an unset Circular
  multi-record legend position, and rewrote the Result's named stroke color as
  hex (F-1).

- A CLI Linear BLAST Session keeps its comparison read-only, and **Inherit
  saved comparison** then Generate succeeds: each Linear file takes the record
  identity of the committed request instead of the CLI binding uid, and the
  saved comparison of an older CLI Session (version 42) is promoted to the
  current request schema before reuse (SE-06, N-17, D-36).
- A CLI Session keeps its `--legend` position through loading and the first
  Generate. The legend position was written into the wrong layout slot and
  replaced by the Web default (SE-07).
- **Source recipe** is unavailable, with a reason, when a track slot legend
  label contains `,` or ` #` (the CLI would cut or reject it), and for a Linear
  scale font without a ruler-label font while ruler labels are drawn (the CLI
  ruler labels would follow the scale font) (TR-07, GE-03).

- Circular **Species**, **Strain**, plot title text, position, and font size,
  **Keep Full Definition with Plot Title**, and **Default font size** now apply
  on Generate, as in Linear. Editing them after Generate no longer rebuilds the
  record definitions from the whole input file, which replaced the region
  length, GC%, and record label in the Result, Save Session, and exports
  (IN-01).
- The global block, line, axis, and scale stroke colors and widths now apply on
  Generate. A cleared or invalid width no longer reaches the current Result,
  and the global setting no longer overwrites a feature's own stroke before
  Generate. Stroke edits on selected features and legend entries stay live
  (GE-02).

- Label text and label visibility edits are no longer cleared when another
  Result, record, or hidden feature changes the displayed diagram. After
  Generate replaces a source file, only the edits of features that no longer
  exist are removed (FE-01).
- With one Result per record, the Features list and its **Edit** actions follow
  the displayed Result, and its color, visibility, legend, and label edits reach
  the other Results when they are displayed, including their export and saved
  Session (FE-02, FE-03, PV-09).
- Redo, or an unrelated Undo, keeps a feature hidden by **Exact product** or
  **Exact protein ID**; the editor matches these rules like Generate, ignoring
  case (FE-04).
- **Reset fill color** uses the default color of the feature it resets, also
  after the Reset dialog of another feature was canceled (FE-10).
- A label rerender no longer applies settings that wait for Generate, such as
  **Species** or a block stroke width (N-16).

- A legend order made with **Sort** or the move buttons, and a renamed legend
  entry without features such as **GC content**, now survive Generate and Save
  Session; entries that appear later follow the ordered ones (PV-02, PV-03,
  D-08).
- Canvas padding now survives Generate and reaches every Result of a batch,
  without being applied twice (PV-07, D-09).
- Renaming a legend entry with features to the caption of another entry of a
  different color offers **Merge**, **Suffix**, or **Cancel** instead of
  failing with an unrecognized error (PV-04, D-06).
- A Linear **Legend position** change now applies on Generate, like Circular;
  the in-place move did not match the generated layout. Changing it after a
  Generate without a legend no longer raises an error (PV-10, D-30, GE-07).
- Closing the **Editor** with Close or Escape returns keyboard focus to the
  Editor toggle (PV-11), and the **Legend** tab help names the edits it offers
  (PV-12).

- Web LOSAT keeps one job per source file pair, but a record no longer
  searches a database that contains itself unless that search was requested.
  Two records in one file now give the same links as two separate files with
  **Max target seqs** `1`, where self hits previously removed every link.
  Comparisons within one file need more jobs. E-values for multi-record files
  still differ from the CLI, which searches each record pair (CO-04, D-40).
- The Settings LOSAT job count now comes from the plan that Generate runs
  (N-10). Run Info states the E-value database.
- A LOSATP Generate after a Feature visibility change searches again instead of
  reusing comparisons that still included hidden proteins (CO-02).
- Threaded LOSAT without cross-origin isolation fails with a LOSAT diagnostic
  instead of UNKNOWN, and the option reads **Threaded (unavailable here)**
  (CO-01).
- Rotating a Linear record with **Orient feature forward** now sets the File
  card's reverse-complement checkbox, so LOSATN ribbons stay on the homologous
  region and the checkbox controls the orientation again (CO-03, N-09).
- A renamed Similarity group keeps its name only on the same members after
  regrouping; other names are saved until those members return (CO-08).
- The web app now shows **Comparison table rows were placed by position.**
  when an uploaded comparison table has sequence IDs that match no displayed
  record. The CLI already logged this warning; the browser discarded it
  (CO-06 carry-over, PD-OI-074). Sessions store the notice in
  `runMetadata.comparisonWarnings`.

- **Changed (command line and Python API):** comparison tables (`-b/--blast`,
  `--comparisons_table`, `linear_comparisons`) are now read in the search
  frame, the selected and cropped record on its source strand. Before, a table
  was read after `--reverse_complement` or a reversed region was applied, so a
  forward BLAST table drew the wrong region on a reversed record. A table
  written for the reversed display needs its coordinates of that record
  converted (`L + 1 - x`). Output changes only for comparisons that touch a
  reverse-complemented record (N-08, CO-07, PD-OI-073).
- A comparison row outside its cropped record stops the run with a
  comparison-input error instead of being drawn beyond the record (N-07).
- **Save Raw LOSAT TSV** of a reverse-complemented record uploads again to the
  same ribbons; the web app no longer converts LOSAT rows to the displayed
  orientation (CO-07).
- Linear match popups and match FASTA headers report input-file coordinates,
  as feature popups do, with the table interval added for a cropped record. The
  SVG record group of a cropped or reverse-complemented Linear record carries
  `data-gbdraw-record-source-start`, `-end`, and `-step` (CO-10, PD-OI-076).
- Sessions saved from `main` with a reverse-complemented Linear record load
  with the same ribbons: the CLI and the web app convert their stored
  comparison rows to the search frame once at Load, and the web app rewrites
  the stored table bytes. A CLI sidecar of `-b` with `--reverse_complement`
  replays the ribbons of the original run. An empty Similarity alignment
  target (`""`) in a `main` Web Session is read as no target. A comparison row
  outside its record reports the new `SEARCH_FRAME` reason of
  `COMPARISON_INPUT` (CO-07, PD-OI-073).

- Web Custom Track Slots: a Depth row is added only when a logical Depth series
  gets its first file and removed when the series loses its last file, in both
  Circular and Linear. Unrelated toggles and **Add Depth TSV series** no longer
  re-create, re-enable, or move Depth rows, and removing a Depth file with the
  uploader **Remove** leaves a stack that generates (TR-02, TR-03).
- Linear: an enabled Depth row whose series has no file shows the same row issue
  as Circular and stops Generate on that row (TR-03).
- Turning off **Hide GC Content** or **Hide GC Skew** while the custom stack is
  off restores the rows that Hide disabled; rows you disabled stay disabled
  (TR-06).
- A disabled track row no longer shows another row's resolved **(auto)**
  geometry (TR-08).
- **Reset to Tuckin**, **Reset to Middle**, and **Reset to Spreadout** follow
  **Show Coordinate Scale**, like **Reset** (TR-09).

- Feature Search and Interactive SVG search: **All** no longer matches
  nucleotide or amino-acid sequences or `/translation` values; use the
  **Nucleotide** and **Amino acid** fields for sequence search. The raw 0-based
  **Start**/**End** search values are removed; **Location** matches the 1-based
  INSDC location (FE-06, FE-11).
- The feature list, feature popup, hover summary, match popup feature
  sections, and Interactive SVG popups show every part of split and
  origin-spanning locations, 1-based, with the summed length (FE-11, N-11).
- A one-feature label, color, or visibility rule no longer spreads to a
  feature whose qualifier value differs only in case (FE-08).
- Specific-color tables accept `none` in any case in the Web app and the CLI.
  Both reject hex colors with alpha and unknown names; the CLI reports the line
  and no longer reads `None`, `NA`, or `null` cells as blank (FE-12).
- Web PDF pages convert CSS px to pt, so they are 75% of their former size and
  match the CLI PDF. Curved and tick labels keep their spaces in the PDF text
  layer (PV-05, PV-06).

- Every visible Web form control has an accessible name. Controls with a
  visible label use that label, including **Window**, **Step**,
  **Dinucleotide**, **GC Content Mode**, and the depth and color-rule fields;
  Custom Track Slots row controls are named by slot id in both modes.
  Placeholders and state-dependent titles are no longer names (TR-10).
- Every **?** help tip is a **Help** button that hover, keyboard focus, a click,
  or a tap opens and **Escape** closes. Tips sit outside labels, and the
  control a tip explains references its text as the accessible description
  (TR-10, PD-OI-057).
- Custom Track Slots rows put the move, duplicate, and remove buttons on their
  own line in both modes, so the slot id and renderer stay readable (TR-11).
- The Vibrio harveyi group Gallery tutorial links to the current Linear record
  layout reference, and a packaging test checks every tutorial link (TR-12).

<!-- web-gui-audit-20260930 P19 -->

<!-- web-gui-audit-20260930 P20 -->

## [0.14.0](./docs/RELEASE_NOTES_0.14.0.md)

- Added circular-record display-start rotation and manual feature lane placement.
- Improved Run Info / Exact replay, Save Session resource preservation, preview
  navigation, and Generate processing-stage feedback.
- Includes the beta's package-root Python API and current session/request
  compatibility, plus isolated wheel/sdist installation and local GUI packaging fixes.
- See the [full release notes](./docs/RELEASE_NOTES_0.14.0.md) for migration,
  compatibility, and installation availability. Publication dates are recorded
  in [GitHub Releases](https://github.com/satoshikawato/gbdraw/releases).

## [0.14.0b0](./docs/RELEASE_NOTES_0.14.0b0.md) — unreleased (beta)

- Added a small top-level Python interface (`read_genbank`, `read_gff`,
  `draw_circular`, `draw_linear`, mode-specific `CircularOptions` /
  `LinearOptions`, and a first-party `Diagram` result) alongside the
  existing typed `gbdraw.api` request/session/table contracts.
- One Circular function now handles both single- and multi-record input.
- Removed obsolete low-level convenience re-exports from `gbdraw.api`.
- See [the full release notes](./docs/RELEASE_NOTES_0.14.0b0.md) for the
  complete list of changes, including architecture/API and session-format
  updates.

## Earlier releases

Releases before 0.14.0b0 predate this changelog and were not recorded with
per-version release notes. Their tags and dates are listed below for
reference; see `git log <tag>` or the
[GitHub tag list](https://github.com/satoshikawato/gbdraw/tags) for the
commits each one contains.

| Version | Date |
| --- | --- |
| 0.13.0 | 2026-07-05 |
| 0.12.1 | 2026-06-27 |
| 0.12.0 | 2026-06-26 |
| 0.11.0 | 2026-05-07 |
| 0.10.0 | 2026-04-29 |
| 0.9.2 | 2026-04-08 |
| 0.9.1 | 2026-04-06 |
| 0.9.0 | 2026-04-06 |
| 0.8.0 | 2025-12-18 |
| 0.7.0 | 2025-10-26 |
| 0.6.0 | 2025-10-09 |
| 0.5.3 | 2025-09-29 |
| 0.5.2 | 2025-09-09 |
| 0.5.1 | 2025-09-08 |
| 0.5.0 | 2025-09-01 |
| 0.4.0 | 2025-08-07 |
| 0.3.0 | 2025-07-24 |
| 0.2.0 | 2025-05-25 |
| 0.1.1 | 2025-05-18 |
| 0.1.0 | 2025-05-14 |
