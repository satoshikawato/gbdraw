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

<!-- web-gui-audit-20260930 P03 -->

<!-- web-gui-audit-20260930 P04 -->

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

- **Changed:** The CLI, the Python API, and typed Web or Session requests
  reject the same invalid values with a named option or setting. `--window`, `--step`,
  `--depth_window`, and `--depth_step` must be positive integers; `-n XY`,
  `-n G`, and other values that are not two of `A`, `C`, `G`, `T`, `U` are
  rejected instead of drawing flat tracks or raising `IndexError`. `U` is
  counted as `T`, so `AU` matches `AT`. Font sizes must be greater than zero
  and stroke widths zero or greater; `--block_stroke_width -1` no longer ends
  in a traceback. Offsets, spacing, `track_axis_gap`, and label rotation are
  unchanged.
- **Fixed:** A Circular definition font set in a Web or Session request
  without an interval again uses the font size plus 2 as the line interval,
  as in 0.13.0 and the CLI. Python API `config_overrides` follow the same
  rule.
- **Fixed:** Qualifier Priority and label whitelist edits apply to Sessions
  written by the CLI or saved from `main`. Label maps compiled by an earlier
  render are no longer kept as preserved settings or reused when a table is
  attached.
- **Changed:** Python validation failures report field, reason, Track row,
  Depth series, and canonical setting identifiers to the Web instead of an
  unclassified error. A test keeps unclassified validation sites from growing.

<!-- web-gui-audit-20260930 P07 -->

<!-- web-gui-audit-20260930 P08 -->

<!-- web-gui-audit-20260930 P09 -->

<!-- web-gui-audit-20260930 P10 -->

<!-- web-gui-audit-20260930 P11 -->

<!-- web-gui-audit-20260930 P12 -->

<!-- web-gui-audit-20260930 P13 -->

<!-- web-gui-audit-20260930 P14 -->

<!-- web-gui-audit-20260930 P15 -->

<!-- web-gui-audit-20260930 P16 -->

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

<!-- web-gui-audit-20260930 P18 -->

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
