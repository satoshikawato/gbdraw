[Home](./DOCS.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [CLI Reference](./CLI_Reference.md)

# Tutorials Audit Implementation Plan

Date: 2026-07-10

This plan turns the `docs/TUTORIALS` audit into concrete documentation work.
The goal is to keep the tutorials aligned with the current CLI behavior while
adding coverage for important features that are currently only documented in
the CLI reference, recipes, gallery tutorials, or implementation tests.

## Goals

- Fix tutorial text that is inconsistent with current CLI validation or can be
  read as supporting invalid option combinations.
- Make each command block reproducible by either creating every referenced
  local input file or explicitly pointing to the source of that input.
- Add tutorial coverage for high-value current features that users are likely
  to miss if they only read `docs/TUTORIALS`.
- Keep `docs/TUTORIALS` task-oriented. Use the CLI reference for exhaustive
  option listings.
- Prefer canonical public spellings in new tutorial text, especially
  `interactive-svg` instead of the accepted compatibility alias
  `interactive_svg`.
- Avoid broad rewrites and line-ending churn in existing tutorial files.

## Non-Goals

- Do not change CLI behavior as part of this documentation task.
- Do not duplicate the full CLI reference inside tutorials.
- Do not regenerate example SVGs unless a tutorial command is intentionally
  changed and the existing image no longer represents the command.
- Do not alter web gallery tutorial JSON in this pass.

## Confirmed Documentation Fixes

### Precomputed BLAST vs Protein Blastp

Problem:

- `Tutorial 2` currently says comparative plots need one BLAST result file for
  each adjacent comparison.
- Current CLI also supports generated protein comparisons with
  `--protein_blastp_mode pairwise`, `orthogroup`, and `collinear`.
- Current CLI rejects `--protein_blastp_mode` combined with `-b/--blast`.

Implementation:

- Change the Tutorial 2 input section to distinguish:
  - precomputed nucleotide/protein BLAST tables supplied with `-b/--blast`
  - generated protein blastp workflows supplied with `--protein_blastp_mode`
- Add an explicit warning that `-b/--blast` and `--protein_blastp_mode` are
  mutually exclusive.
- Replace the sentence "Files supplied with `-b/--blast` are not used to infer
  orthogroups or collinear blocks; they remain a separate comparison input."
  with wording that says users must choose one comparison source per run.
- Add at least one minimal command for:
  - `--protein_blastp_mode pairwise`
  - `--protein_blastp_mode orthogroup`
  - `--protein_blastp_mode collinear`

Acceptance checks:

- `python -m gbdraw.cli linear --help` still matches the documented option
  names.
- The tutorial no longer implies that `-b/--blast` can be used together with
  `--protein_blastp_mode`.

### Linear Label Rendering Wording

Problem:

- `Tutorial 3` says to use `--label_placement above_feature` "without
  `--label_rendering`."
- Current behavior allows the default `--label_rendering auto`; only
  `embedded_only` and `external_only` conflict with `above_feature`.

Implementation:

- Reword the section to say:
  - `--label_placement above_feature` works with the default
    `--label_rendering auto`
  - do not combine it with `--label_rendering embedded_only` or
    `--label_rendering external_only`

Acceptance checks:

- The text matches CLI validation in `gbdraw/linear.py`.

### Reproducible Input Files

Problem:

- Some Tutorial 3 command blocks reference local files without showing how to
  obtain or create them in the same tutorial flow.
- Examples include `O157_H7.gbk`, `HmmtDNA.gbk`, and `label_override.tsv`.

Implementation:

- Add download commands or cross-links for each referenced GenBank input.
- For `label_override.tsv`, either:
  - write the example table to `label_override.tsv` before the command, or
  - change the command to use an already existing example file and explain its
    source.
- Keep table examples small and focused.

Acceptance checks:

- A reader can follow Tutorial 3 from a clean directory without guessing where
  local inputs come from.

## New Tutorial Coverage

### Tutorial 4: Protein Comparisons Without Precomputed BLAST

Purpose:

- Teach the generated protein blastp workflows that are now central to linear
  comparison plots.

Content:

- Required inputs: two or more annotated GenBank or GFF3+FASTA records with CDS
  translations or translatable CDS features.
- Runtime selection:
  - bundled LOSAT on Linux x86_64 when available
  - `losat` on `PATH`
  - NCBI BLAST+ `blastp` fallback
  - explicit `--losatp_bin`, `--ncbi_blastp_bin`, and `--losatp_threads`
- Commands for:
  - `--protein_blastp_mode pairwise`
  - `--protein_blastp_mode orthogroup`
  - `--protein_blastp_mode collinear`
  - `--collinear_min_anchors 2`
  - `--collinear_color_mode orientation_identity`
  - `--show_labels orthogroup_top`
- A short explanation of when to prefer precomputed `-b/--blast` instead.

Acceptance checks:

- Commands do not include `-b/--blast`.
- The tutorial states that NCBI BLAST+ output may not exactly match LOSAT.

### Tutorial 5: Table-Driven CLI Inputs

Purpose:

- Promote row-coupled TSV workflows out of the CLI reference into a practical
  tutorial.

Content:

- `--records_table` for linear inputs:
  - `gbk`
  - `record_label`
  - `record_subtitle`
  - `record_id`
  - `region`
  - `reverse_complement`
  - `order`
- `--records_table` for GFF3+FASTA rows.
- `--records_table` for circular multi-record placement:
  - `row`
  - `column`
  - `--multi_record_canvas`
- `--conservation_table` for circular conservation rings.
- `--circular_track_table` for circular slot order and AT skew.
- A short note that table-relative paths resolve against the table file.

Acceptance checks:

- The tutorial states the relevant mutual exclusions:
  - `--records_table` cannot be combined with `--gbk`, `--gff`, or `--fasta`
  - circular table placement should not be combined with
    `--multi_record_position`
  - `--conservation_table` cannot be combined with `--conservation_blast`,
    `--conservation_labels`, or `--conservation_colors`
  - `--circular_track_table` cannot be combined with inline circular track slot
    options

### Tutorial 6: Coverage Depth and Quantitative Tracks

Purpose:

- Document coverage/depth rendering, GC percent mode, additional skew tracks,
  and axis controls in one workflow.

Content:

- Depth TSV format: `reference`, `position`, `depth`.
- `--depth` for a single logical depth track.
- Repeatable `--depth_track` for multiple logical tracks.
- Track labels and colors:
  - `--depth_track_label`
  - `--depth_track_color`
- Axis controls:
  - `--show_depth_axis`
  - `--show_depth_ticks`
  - `--depth_large_tick_interval`
  - `--depth_small_tick_interval`
  - `--share_depth_axis`
- Scaling:
  - `--depth_log_scale`
  - `--no_depth_log_scale`
  - `--depth_min`
  - `--depth_max`
- GC percent mode:
  - `--gc_content_mode percent`
  - `--gc_content_min_percent`
  - `--gc_content_max_percent`
  - GC content tick options

Acceptance checks:

- The tutorial explains that `--depth` and `--depth_track` are alternatives.
- Circular and linear examples use the mode-specific dimensions:
  `--depth_width` for circular and `--depth_height` for linear.

### Tutorial 7: Linear Layout, Definitions, and Rulers

Purpose:

- Give linear plotting its own customization tutorial instead of scattering the
  details across circular-focused pages.

Content:

- `--track_layout above|middle|below`
- `--track_axis_gap`
- `--scale_style ruler`
- `--ruler_on_axis`
- `--record_label`
- `--record_subtitle`
- `--plot_title`
- `--plot_title_position center|top|bottom`
- `--definition_line_style`
- `--show_replicon`
- `--hide_accession`
- `--hide_length`
- `--keep_definition_left_aligned`
- `--linear_track_order` and a minimal `--linear_track_slot` example

Acceptance checks:

- Examples use options that are accepted by `gbdraw linear --help`.
- The tutorial explains that `--record_label` and `--record_subtitle` are
  repeatable and order-sensitive unless a `--records_table` is used.

### Tutorial 8: Interactive SVG and Session Round Trips

Purpose:

- Show users how to move between CLI output, standalone interactive SVG, saved
  GUI sessions, and regenerated diagrams.

Content:

- `-f interactive-svg`
- Expected outputs:
  - `output.svg`
  - `output.interactive.svg`
- What interactive SVG includes at a high level:
  - feature metadata
  - pairwise match metadata where available
  - orthogroup metadata where available
- `--save_session`
- `--session_output`
- `--session`
- `gbdraw gui` as local web UI entry point

Acceptance checks:

- New tutorial uses canonical `interactive-svg` spelling.
- Commands do not require CairoSVG for interactive SVG output.

### Tutorial 9: Feature Visibility and Shapes

Purpose:

- Cover per-feature visibility and shape overrides, which are powerful but
  currently easy to miss.

Content:

- `--feature_shape TYPE=arrow|rectangle`
- `--feature_visibility_table`
- Visibility table actions:
  - `show`
  - `off`
  - `exclude_matching`
- Selectors such as:
  - record ID
  - feature type
  - qualifier
  - regex value
- Relationship to `-k/--features`, color tables, and label tables.

Acceptance checks:

- The visibility table uses the current action names and avoids deprecated
  `suppress` wording.

## Index and Cross-Link Updates

Update `docs/TUTORIALS/TUTORIALS.md` after adding or renaming tutorials:

- Keep Tutorials 1-3 as the introductory path.
- Add an "Advanced workflows" group for Tutorials 4-9.
- Link to gallery tutorials as examples of full web UI reproductions.
- Link back to `CLI_Reference.md` for exhaustive syntax.

Update related pages only where needed:

- `docs/QUICKSTART.md`: mention that protein comparisons and table manifests
  have dedicated tutorials.
- `docs/RECIPES.md`: keep short copy-paste commands and link to the relevant
  tutorial for explanation.
- `docs/CLI_Reference.md`: avoid duplicating tutorial prose, but ensure any new
  tutorial names are linked from the related guides section if useful.

## Suggested Implementation Order

1. Correct existing Tutorials 2 and 3 wording.
2. Add missing input setup commands to Tutorial 3.
3. Add Tutorial 4 for generated protein comparisons.
4. Add Tutorial 5 for table-driven inputs.
5. Add Tutorial 6 for depth and quantitative tracks.
6. Add Tutorial 7 for linear layout and definition text.
7. Add Tutorial 8 for interactive SVG and sessions.
8. Add Tutorial 9 for feature visibility and shapes.
9. Update `TUTORIALS.md`, `QUICKSTART.md`, and `RECIPES.md` links.
10. Run documentation sanity checks.

## Validation Plan

Minimum checks:

```bash
python -m gbdraw.cli circular --help
python -m gbdraw.cli linear --help
rg -n "interactive_svg|--show_labels \\[\\{all,first,none\\}\\]|Files supplied with `-b/--blast`" docs
```

Recommended command smoke checks:

```bash
python -m gbdraw.cli linear \
  --gbk tests/test_inputs/MjeNMV.gb tests/test_inputs/MelaMJNV.gb \
  --protein_blastp_mode pairwise \
  -o /tmp/gbdraw_tutorial_pairwise \
  -f svg

python -m gbdraw.cli linear \
  --gbk tests/test_inputs/MjeNMV.gb tests/test_inputs/MelaMJNV.gb \
  -b tests/test_inputs/MjeNMV.MelaMJNV.tblastx.out \
  --protein_blastp_mode orthogroup \
  -o /tmp/gbdraw_invalid_mix \
  -f svg
```

Expected result for the second command:

- It must fail with `--protein_blastp_mode cannot be used with -b/--blast`.

Optional broader checks:

```bash
pytest tests/test_cli_tables.py -v
pytest tests/test_linear_label_placement.py -v
pytest tests/test_interactive_svg_cli_format.py -v
```

## Risks and Mitigations

- Risk: tutorials become too long.
  Mitigation: split advanced workflows into separate pages and keep each page
  focused on one practical outcome.
- Risk: tutorial commands depend on external NCBI availability.
  Mitigation: prefer `tests/test_inputs` or committed `examples` where possible,
  and use download commands only when the source is stable and already used by
  existing docs.
- Risk: LOSAT and NCBI BLAST+ availability differs by platform.
  Mitigation: document runtime resolution explicitly and avoid claiming byte-for
  byte identical output between runtimes.
- Risk: example images drift from commands.
  Mitigation: only change image-linked commands when necessary, and regenerate
  images in a separate, reviewable step if needed.

## Completion Criteria

- Existing tutorial mismatches are corrected.
- New tutorial pages cover the major CLI features currently absent from
  `docs/TUTORIALS`.
- `docs/TUTORIALS/TUTORIALS.md` links all tutorial pages.
- All new examples use accepted current option names.
- The validation checks above pass or have documented reasons for not running.
