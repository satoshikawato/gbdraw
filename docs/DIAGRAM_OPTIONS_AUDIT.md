[Python API](./PYTHON_API.md) | [API improvement plan](./PYTHON_API_IMPROVEMENT_PLAN.md) | [Follow-up plan](./PYTHON_API_FOLLOWUP_PLAN.md)

# `DiagramOptions` field audit

- Date: 2026-07-14
- Version: `0.14.0b0`
- Contract: 70 fields in `gbdraw.api.DiagramOptions`

## All 70 fields reach an owner

Every field is read by at least one public `build_*` adapter and has a non-default
forwarding test. No dead field was found. The audit found 31 mode-specific fields,
two compatibility aliases, and two depth input forms whose behavior changes by
mode. Those findings do not justify splitting the dataclass in this release.

The main problem is discoverability. A caller can set a field for the wrong mode,
and the current adapter may ignore it. The action for `0.14.0b0` is to keep the
contract, publish the mode table below, and record validation work separately.

## Reading the tables

Mode abbreviations:

- `CS`: `build_circular_diagram`
- `CM`: `build_circular_multi_diagram`
- `L`: `build_linear_diagram`

Owner codes identify the function or subsystem that consumes the normalized value:

| Code | Final owner |
|---|---|
| `C` | `apply_config_overrides`, `GbdrawConfig`, canvas/configurator constructors |
| `F` | feature table readers, `compile_feature_visibility_rules`, `FeatureDrawingConfigurator` |
| `LB` | `LabelsFilteringConfig` and circular/linear label policy |
| `G` | GC content/skew configurators and scalar track renderers |
| `D` | `normalize_depth_tracks`, `build_depth_track_dataframes`, depth configurator/renderers |
| `CV` | conservation source loader, record normalization, circular conservation slots |
| `T` | plot-title resolvers, `DefinitionGroup`, canvas placement |
| `P` | BLAST/protein comparison loaders, `BlastMatchConfigurator`, linear comparison groups |
| `K` | protein grouping and collinearity builders/normalizers |

Forwarding test codes:

| Code | Test |
|---|---|
| `H-S` | `test_diagram_options_forward_shared_non_default_values` |
| `H-B` | `test_diagram_option_bundles_forward_non_default_values` |
| `H-LB` | `test_diagram_options_attach_explicit_label_tables` and `test_diagram_options_attach_explicit_label_files` |
| `H-D` | `test_plural_depth_options_forward_with_mode_specific_shape` |
| `H-C` | `test_diagram_options_forward_circular_only_non_default_values` |
| `H-X` | `test_diagram_options_forward_linear_only_non_default_values` |

Owner tests are grouped by subsystem: `O-C` is
`tests/test_api_library_usage.py`; `O-F` is `tests/test_feature_visibility.py` and
`tests/test_feature_shapes.py`; `O-LB` is `tests/test_label_rendering_policy.py` and
`tests/test_label_overrides.py`; `O-G` is `tests/test_gc_content_percent.py` and the
track-slot suites; `O-D` is `tests/test_depth_track.py`; `O-CV` is
`tests/test_circular_conservation.py`; `O-T` is
`tests/test_circular_multi_canvas.py` and `tests/test_linear_track_layout.py`;
`O-P` is `tests/test_comparisons.py`; `O-K` is `tests/test_collinearity.py` and
`tests/test_protein_colinearity.py`.

`API` in the Docs column means `docs/PYTHON_API.md` names the field or demonstrates
the capability. `Audit` means this file is the only field-level Python API
reference. CLI documentation may still describe the equivalent CLI option.

## Configuration, feature, and label fields

| # | Field, type, default | Mode and reader | Forwarding and final owner | Test | Docs | Classification and action |
|---:|---|---|---|---|---|---|
| 1 | `config`; config or dict; `None` | CS, CM, L | `_resolve_diagram_options_config` to `config_dict` or `cfg` to `C` | H-S, O-C | API | live; keep |
| 2 | `config_overrides`; mapping; `None` | CS, CM, L | config resolver to `config_overrides` or merged `cfg` to `C` | H-S, O-C | Audit | live; keep |
| 3 | `colors`; `ColorOptions`; `None` | CS, CM, L | bundle fields to color readers and `F` | H-B, O-F | Audit | live bundle; keep |
| 4 | `tracks`; `TrackOptions`; `None` | CS, CM, L | mode-specific slot parser/axis normalizer to `F`, `G`, or `D` | H-B, O-F, O-G, O-D | API | live bundle; keep and document subfield modes |
| 5 | `output`; `OutputOptions`; `None` | CS, CM, L | prefix/legend/title position to canvas and `T` | H-B, O-T | API | live bundle; keep |
| 6 | `selected_features_set`; sequence; `None` | CS, CM, L | same-name assembler argument to `F` | H-S, O-F | API | live; keep |
| 7 | `feature_table`; DataFrame; `None` | CS, CM, L | `feature_table` to compatibility resolver to `F` | H-S, O-F | Audit | compatibility alias; keep |
| 8 | `feature_table_file`; path; `None` | CS, CM, L | `feature_table_file` to visibility reader to `F` | H-S, O-F | Audit | compatibility alias; keep |
| 9 | `feature_visibility_table`; DataFrame; `None` | CS, CM, L | same-name argument to visibility resolver to `F` | H-S, O-F | Audit | canonical input; keep |
| 10 | `feature_visibility_table_file`; path; `None` | CS, CM, L | same-name argument to visibility reader to `F` | H-S, O-F | Audit | canonical input; keep |
| 11 | `label_whitelist_table`; DataFrame; `None` | CS, CM, L | config resolver to `whitelist_df` to `LB` | H-LB, O-LB | API | DataFrame/file pair; keep |
| 12 | `label_whitelist_file`; path; `None` | CS, CM, L | whitelist reader to `whitelist_df` to `LB` | H-LB, O-LB | API | DataFrame/file pair; keep |
| 13 | `qualifier_priority_table`; DataFrame; `None` | CS, CM, L | config resolver to `qualifier_priority_df` to `LB` | H-LB, O-LB | API | DataFrame/file pair; keep |
| 14 | `qualifier_priority_file`; path; `None` | CS, CM, L | priority reader to `qualifier_priority_df` to `LB` | H-LB, O-LB | API | DataFrame/file pair; keep |
| 15 | `label_override_table`; DataFrame; `None` | CS, CM, L | config resolver to `label_override_df` to `LB` | H-LB, O-LB | API | DataFrame/file pair; keep |
| 16 | `label_override_file`; path; `None` | CS, CM, L | override reader to `label_override_df` to `LB` | H-LB, O-LB | API | DataFrame/file pair; keep |
| 17 | `feature_shapes`; mapping; `None` | CS, CM, L | same-name assembler argument to `F` and feature drawers | H-S, O-F | Audit | live; keep |

`feature_table*` and `feature_visibility_table*` have the same final owner. The
former names are compatibility aliases. Passing both table forms, or both file
forms, raises `ValidationError`; no deprecation is proposed for the beta release.

## Scalar and depth fields

| # | Field, type, default | Mode and reader | Forwarding and final owner | Test | Docs | Classification and action |
|---:|---|---|---|---|---|---|
| 18 | `dinucleotide`; string; `"GC"` | CS, CM, L | `dinucleotide` to GC/skew configurators and `G` | H-S, O-G | Audit | live; keep |
| 19 | `window`; integer; `None` | CS, CM, L | `window` to mode resolver, then `G` and `D` | H-S, O-G, O-D | Audit | live; keep |
| 20 | `step`; integer; `None` | CS, CM, L | `step` to mode resolver, then `G` and `D` | H-S, O-G, O-D | Audit | live; keep |
| 21 | `depth_window`; integer; `None` | CS, CM, L | depth window resolver to `D` | H-S, O-D | Audit | live; keep |
| 22 | `depth_step`; integer; `None` | CS, CM, L | depth step resolver to `D` | H-S, O-D | Audit | live; keep |
| 23 | `depth_table`; DataFrame; `None` | CS, CM, L | singular depth source to `D` | H-S, O-D | Audit | input form; keep |
| 24 | `depth_file`; path; `None` | CS, CM, L | singular depth reader to `D` | H-S, O-D | API | input form; keep |
| 25 | `depth_tables`; DataFrame sequence; `None` | CS, CM, L | CS uses element 0 as singular fallback; CM/L pass the sequence to `D` | H-D, O-D | Audit | mode-dependent input form; keep and validate CS length later |
| 26 | `depth_files`; path sequence; `None` | CS, CM, L | CS uses element 0 as singular fallback; CM/L pass the sequence to `D` | H-D, O-D | Audit | mode-dependent input form; keep and validate CS length later |
| 27 | `depth_track_tables`; nested DataFrame sequence; `None` | CS, CM, L | per-record/per-track source matrix to `D` | H-S, O-D | Audit | live; keep |
| 28 | `depth_track_files`; nested path sequence; `None` | CS, CM, L | per-record/per-track file matrix to `D` | H-S, O-D | Audit | live; keep |
| 29 | `depth_track_labels`; string sequence; `None` | CS, CM, L | track metadata normalization to `D` and legend | H-S, O-D | Audit | live; keep |
| 30 | `depth_track_colors`; string sequence; `None` | CS, CM, L | track metadata normalization to `D` and legend | H-S, O-D | Audit | live; keep |
| 31 | `depth_track_heights`; scalar sequence; `None` | L | linear adapter to normalized depth specs and linear canvas | H-X, O-D | Audit | linear-only; keep and document |
| 32 | `depth_track_large_tick_intervals`; scalar sequence; `None` | CS, CM, L | track metadata normalization to depth axes | H-S, O-D | Audit | live; keep |
| 33 | `depth_track_small_tick_intervals`; scalar sequence; `None` | CS, CM, L | track metadata normalization to depth axes | H-S, O-D | Audit | live; keep |
| 34 | `depth_track_tick_font_sizes`; scalar sequence; `None` | CS, CM, L | track metadata normalization to depth axes | H-S, O-D | Audit | live; keep |

For CS, `depth_table` and `depth_file` take precedence over their plural forms.
When the singular value is absent, the adapter uses only the first plural entry.
CM and L preserve the plural sequence. A later behavior change should reject more
than one plural entry in CS instead of ignoring entries after index 0.

## Conservation and title fields

| # | Field, type, default | Mode and reader | Forwarding and final owner | Test | Docs | Classification and action |
|---:|---|---|---|---|---|---|
| 35 | `conservation_blast_files`; path sequence; `None` | CS, CM | conservation source loader to `CV` | H-C, O-CV | API | circular-only; keep |
| 36 | `conservation_dataframes`; DataFrame sequence; `None` | CS, CM | conservation source loader to `CV` | H-C, O-CV | Audit | circular-only; keep |
| 37 | `conservation_reference`; literal; `"auto"` | CS, CM | reference normalizer to `CV` | H-C, O-CV | API | circular-only; keep |
| 38 | `conservation_labels`; string sequence; `None` | CS, CM | track normalization and legend to `CV` | H-C, O-CV | API | circular-only; keep |
| 39 | `conservation_colors`; string sequence; `None` | CS, CM | track normalization and gradient legend to `CV` | H-C, O-CV | Audit | circular-only; keep |
| 40 | `conservation_ring_width`; float; `None` | CS, CM | slot width override to `CV` layout | H-C, O-CV | Audit | circular-only; keep |
| 41 | `conservation_ring_gap`; float; `None` | CS, CM | slot spacing override to `CV` layout | H-C, O-CV | Audit | circular-only; keep |
| 42 | `plot_title`; string; `None` | CS, CM, L | normalized title to `T` | H-S, O-T | Audit | live; keep |
| 43 | `plot_title_font_size`; float; `None` | CS, CM, L | title size resolver/config override to `T` | H-S, O-T | Audit | live; keep |
| 44 | `keep_full_definition_with_plot_title`; boolean; `False` | CS, CM | circular definition profile selection in `T` | H-C, O-T | Audit | circular-only; keep |
| 45 | `species`; string; `None` | CS, CM | circular `DefinitionGroup` override in `T` | H-C, O-T | API | circular-only; keep |
| 46 | `strain`; string; `None` | CS, CM | circular `DefinitionGroup` override in `T` | H-C, O-T | Audit | circular-only; keep |

Linear diagrams derive definition text from each record and do not read `species`,
`strain`, or `keep_full_definition_with_plot_title` from `DiagramOptions`.

## Linear comparison fields

| # | Field, type, default | Mode and reader | Forwarding and final owner | Test | Docs | Classification and action |
|---:|---|---|---|---|---|---|
| 47 | `blast_files`; path sequence; `None` | L | nucleotide BLAST input to `P` | H-X, O-P | API | linear-only; keep |
| 48 | `protein_comparisons`; DataFrame sequence; `None` | L | precomputed comparison input to `P` | H-X, O-P | API | linear-only; keep |
| 49 | `orthogroups`; result object; `None` | L | precomputed group metadata to linear assembly/alignment | H-X, O-K | Audit | linear-only; keep |
| 50 | `protein_blastp_mode`; literal; `"none"` | L | mode normalizer to protein search or collinearity owner | H-X, O-K | API | linear-only; keep |
| 51 | `pairwise_match_style`; literal; `"ribbon"` | L | style resolver/config override to comparison groups | H-X, O-P | Audit | linear-only; keep |
| 52 | `collinearity_blocks`; result or block sequence; `None` | L | block conversion to comparisons and `K` metadata | H-X, O-K | Audit | linear-only; keep |
| 53 | `collinearity_params`; parameter object; `None` | L | collinearity builder `params` in `K` | H-X, O-K | Audit | linear-only; keep |
| 54 | `collinearity_unit_mode`; string; `"auto"` | L | collinearity builder `unit_mode` in `K` | H-X, O-K | Audit | linear-only; keep |
| 55 | `collinearity_anchor_mode`; string; `"rbh"` | L | anchor normalizer to collinearity builder `edge_mode` | H-X, O-K | API | linear-only; keep |
| 56 | `collinearity_search_scope`; string; `"adjacent"` | L | scope normalizer to collinearity builder `search_scope` | H-X, O-K | Audit | linear-only; keep |
| 57 | `collinearity_color_mode`; string; `"orientation"` | L | color normalizer to block conversion and comparison legend | H-X, O-K | Audit | linear-only; keep |
| 58 | `losatp_bin`; string; `"losat"` | L | protein search executable in `K` | H-X, O-K | Audit | linear-only; keep |
| 59 | `ncbi_blastp_bin`; string; `None` | L | fallback protein search executable in `K` | H-X, O-K | Audit | linear-only; keep |
| 60 | `losatp_threads`; integer; `None` | L | protein search thread count in `K` | H-X, O-K | Audit | linear-only; keep |
| 61 | `protein_blastp_max_hits`; integer; `5` | L | pairwise protein search `max_hits` in `K` | H-X, O-K | Audit | linear-only; keep |
| 62 | `protein_blastp_candidate_limit`; integer; `None` | L | protein search `candidate_limit` in `K` | H-X, O-K | Audit | linear-only; keep |
| 63 | `orthogroup_membership_mode`; string; `"anchor_core_v1"` | L | membership normalizer to grouping/collinearity in `K` | H-X, O-K | Audit | linear-only; keep |
| 64 | `orthogroup_member_max_hits`; integer; `5` | L | grouping/collinearity member limit in `K` | H-X, O-K | Audit | linear-only; keep |
| 65 | `collinear_max_paralog_links_per_orthogroup`; integer; `2` | L | grouping and collinearity paralog-link limits in `K` | H-X, O-K | Audit | linear-only; keep |
| 66 | `align_orthogroup_feature`; string; `None` | L | linear assembler orthogroup alignment target | H-X, O-K | Audit | linear-only; keep |
| 67 | `evalue`; float; `1e-5` | CS, CM, L | conservation or comparison `BlastMatchConfigurator`; protein search filters | H-S, O-CV, O-P, O-K | API | shared threshold; keep |
| 68 | `bitscore`; float; `50.0` | CS, CM, L | conservation or comparison `BlastMatchConfigurator`; protein search filters | H-S, O-CV, O-P, O-K | API | shared threshold; keep |
| 69 | `identity`; float; `70.0` | CS, CM, L | conservation or comparison filtering and legends; protein search filters | H-S, O-CV, O-P, O-K | API | shared threshold; keep |
| 70 | `alignment_length`; integer; `0` | CS, CM, L | conservation or comparison filtering; protein search filters | H-S, O-CV, O-P, O-K | API | shared threshold; keep |

The comparison fields are mutually constrained. For example, precomputed
`collinearity_blocks` cannot be combined with search mode, protein DataFrames, or
nucleotide BLAST files. These constraints belong to the linear assembler and are
covered by the comparison and collinearity owner suites.

## Classification totals

| Classification | Count | Fields |
|---|---:|---|
| Shared live field | 39 | configuration, feature, label, scalar/depth fields, title, thresholds |
| Circular-only | 10 | conservation fields, circular definition controls |
| Linear-only | 21 | `depth_track_heights` and comparison/protein/collinearity fields |
| Compatibility alias | 2 | `feature_table`, `feature_table_file` (also counted as shared) |
| Dead candidate | 0 | none |

## Decision and backlog

Keep `DiagramOptions` intact for `0.14.0b0`. A split would add conversion code while
39 fields remain shared, and there is no dead field to remove. The public
`assemble_*` functions remain the lower-level escape hatch.

The audit leaves three follow-up items:

1. Reject a non-default mode-specific field when it is passed to the wrong builder,
   or add support where that makes sense. This applies to circular conservation and
   definition fields, linear comparison fields, and `depth_track_heights`.
2. Reject CS `depth_tables` or `depth_files` sequences longer than one. The current
   adapter uses index 0 and ignores later entries.
3. Keep the aliases in fields 7 and 8 for the beta series. Any later deprecation
   needs a warning period and a major-release removal target.

These are validation and compatibility tasks. They do not require a new option
class, and no field should be removed as part of them.

[Python API](./PYTHON_API.md) | [API improvement plan](./PYTHON_API_IMPROVEMENT_PLAN.md) | [Follow-up plan](./PYTHON_API_FOLLOWUP_PLAN.md)
