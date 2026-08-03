[Documentation home](../DOCS.md) | [Python Tutorial](../TUTORIALS/PYTHON/README.md) | [Python how-to guides](../HOW_TO/PYTHON/README.md) | [Typed requests](typed-requests.md)

# Python API reference

The package-root API is the stable, beginner-facing interface for reading annotated sequences, drawing Circular or Linear diagrams, and returning or saving the result. Integrations that need explicit input cardinality, planning, output preflight, current analysis artifacts, or session conversion use the typed [`gbdraw.api` request interface](typed-requests.md).

## Public functions

```text
read_genbank(paths: str | PathLike[str] | Sequence[str | PathLike[str]]) -> list[SeqRecord]
read_gff(gff_paths: str | PathLike[str] | Sequence[str | PathLike[str]], fasta_paths: str | PathLike[str] | Sequence[str | PathLike[str]], *, features: Sequence[str] | None = None) -> list[SeqRecord]
draw_circular(records: RecordCollection, *, options: CircularOptions | None = None, layout: CircularLayout | None = None) -> Diagram
draw_linear(records: RecordCollection, *, options: LinearOptions | None = None, layout: LinearLayout | None = None) -> Diagram
```

`read_genbank()` returns every record from every supplied file. `read_gff()` requires equally sized GFF3 and FASTA path lists and matching sequence IDs. Both drawing functions reject an empty record collection, non-`SeqRecord` members, and option or layout objects for the wrong mode.

## Package-root exports

The `gbdraw` package exports the four functions above, `__version__`, and these public types:

| Area | Types |
|---|---|
| Result | `Diagram` |
| Shared presentation | `FeatureOptions`, `LabelOptions`, `TitleOptions`, `Thresholds`, `DepthTrackOptions` |
| Circular | `CircularOptions`, `CircularLayout`, `CircularTrackOptions`, `ComparisonRingOptions`, `ComparisonRingTrackOptions` |
| Linear | `LinearOptions`, `LinearLayout`, `LinearTrackOptions`, `LinearComparisonOptions` |
| Compatibility aliases | `ConservationOptions`, `ConservationTrackOptions` |

`ConservationOptions` and `ConservationTrackOptions` are identity aliases for `ComparisonRingOptions` and `ComparisonRingTrackOptions`. New code should use the comparison-ring names.

## Shared option defaults

| Type and field | Type | Default |
|---|---|---|
| `FeatureOptions.types` | sequence of feature names or `None` | `CDS`, `rRNA`, `tRNA`, `tmRNA`, `ncRNA`, `misc_RNA`, `repeat_region` |
| `FeatureOptions.color_table` | path or `DataFrame` | `None` |
| `FeatureOptions.default_colors` | path or `DataFrame` | `None` |
| `FeatureOptions.palette` | string | `default` |
| `FeatureOptions.visibility` | path or `DataFrame` | `None` |
| `FeatureOptions.shapes` | mapping | `None` |
| `LabelOptions.whitelist` | path or `DataFrame` | `None` |
| `LabelOptions.qualifier_priority` | path or `DataFrame` | `None` |
| `LabelOptions.overrides` | path or `DataFrame` | `None` |
| `TitleOptions.text` | string or `None` | `None` |
| `TitleOptions.position` | `none`, `center`, `top`, `bottom`, or `None` | `None` |
| `TitleOptions.font_size` | positive number or `None` | `None` |
| `CircularOptions.legend`, `LinearOptions.legend` | legend position token | `right` |
| `CircularOptions.dinucleotide`, `LinearOptions.dinucleotide` | two-base string | `GC` |
| `window`, `step`, `depth_window`, `depth_step` | integer or `None` | `None` |
| `annotations`, `config`, `config_overrides` | mode-appropriate value or `None` | `None` |
| `depth_tracks` | sequence of `DepthTrackOptions` | empty |

Unset `Thresholds` values resolve through the mode profile:

| Mode | E-value | Bitscore | Identity | Alignment length | GC / skew on a fresh request |
|---|---:|---:|---:|---:|---|
| Circular | `1e-5` | `50` | `70` | `0` | on / on |
| Linear | `1e-2` | `50` | `0` | `0` | off / off |

E-value, bitscore, and identity must be finite and non-negative; identity is limited to `0`–`100`. Alignment length must be a non-negative integer.

## Layout and track options

`CircularLayout` defaults to `size="auto"`, `min_radius_ratio=0.55`, `column_gap_ratio=0.10`, `row_gap_ratio=0.05`, and automatic positions. `size` accepts `linear`, `auto`, or `equal`.

`LinearLayout` defaults to `record_gap=24.0` and automatic positions. Position strings identify a displayed record and its grid row, for example `#1@1`.

`CircularTrackOptions` and `LinearTrackOptions` accept an ordered `slots` sequence and a zero-based `axis_index`. Circular tracks also accept `center_reserved_radius`. An explicit slot sequence is authoritative; the axis index must agree with the selected mode and slot types.

Each `DepthTrackOptions` represents one logical series. `source` is one path or `DataFrame` for one displayed record, or one path, `DataFrame`, or `None` per record. `label`, `color`, tick intervals, and tick font size default to `None`. `height` is supported only by Linear diagrams; Circular options reject it.

## Circular options

`CircularOptions.tracks` defaults to an empty `CircularTrackOptions`. `comparison_rings` defaults to an empty `ComparisonRingOptions`; `species` and `strain` default to `None`; `keep_full_definition_with_title` defaults to `False`. Circular title position accepts `None`, `none`, `top`, or `bottom`.

`ComparisonRingOptions` defaults to no tracks, `reference="auto"`, and automatic ring width and gap. `reference` accepts `query`, `subject`, or `auto`. Every `ComparisonRingTrackOptions` requires a BLAST/LOSAT table path or `DataFrame`; its label, color, and comparison sequence source are optional. All tracks in one diagram must use the same source kind. If any ring has a label or color, every ring must provide that field.

`CircularOptions.comparison_rings` is the canonical field. The `conservation` constructor and attribute alias remains available for compatible code, but passing both names is an error.

## Linear options

`LinearOptions.tracks` defaults to an empty `LinearTrackOptions`; `comparisons` defaults to `LinearComparisonOptions()`. Linear title position accepts `None`, `center`, `top`, or `bottom`.

Important `LinearComparisonOptions` defaults are:

| Field | Default |
|---|---|
| `protein_mode` | `none` |
| `match_style` | `ribbon` |
| `collinearity_unit` | `auto` |
| `collinearity_anchor` | `rbh` |
| `collinearity_scope` | `adjacent` |
| `collinearity_color` | `orientation` |
| `losat_executable` | `losat` |
| `blastp_executable` | `None` |
| `threads` | `None` |
| `max_hits` | `5` |
| `candidate_limit` | `None` |
| `orthogroup_membership` | `anchor_core_v1` |
| `orthogroup_member_max_hits` | `5` |
| `max_paralog_links` | `2` |

`protein_mode` accepts `none`, `pairwise`, `orthogroup`, or `collinear`. The `orthogroup` token means gbdraw Similarity groups; it does not claim phylogenetic orthology.

## Canonical label overrides

`config_overrides` uses canonical dotted leaf paths:

| Path | Accepted values |
|---|---|
| `labels.circular.scope` | `none`, `outer`, `both` |
| `labels.circular.placement` | `horizontal`, `radial` |
| `labels.linear.scope` | `none`, `all`, `first`, `orthogroup_top` |
| `labels.linear.placement` | `auto`, `above_feature` |
| `labels.linear.rotation` | finite degrees |
| `labels.rendering` | `auto`, `embedded_only`, `external_only` |

Fresh Python requests reject the retired flat label aliases and old `canvas.*show_labels` paths. Supported persisted-session readers perform compatibility migration; new code should not emit retired names.

## `Diagram` output

`Diagram.mode` is `circular` or `linear`, and `Diagram.records` is the rendered record tuple. `Diagram.save()`, `Diagram.to_svg()`, and `Diagram.to_bytes()` provide four output forms:

| Method | Result |
|---|---|
| `to_svg()` | static SVG text |
| `to_svg(interactive=True)` | interactive SVG text with supported metadata |
| `to_bytes(format="svg")` | SVG, interactive SVG, PNG, PDF, EPS, or PS bytes |
| `save(path, *, format=None, overwrite=False)` | writes exactly one file and returns its `Path` |

PNG, PDF, EPS, and PS require CairoSVG. `save()` infers a known format from the path unless `format` is explicit and refuses to replace an existing file unless `overwrite=True`.

## Errors

Catch `gbdraw.exceptions.GbdrawError` for expected gbdraw failures and `ValidationError` for invalid records, paths, options, or mode combinations. Pin a gbdraw version and compare representative outputs after upgrading; stable Python calls do not imply byte-identical SVG geometry.

## Related

- [Draw and save your first genome diagram from Python](../TUTORIALS/PYTHON/first-genome-diagram.md)
- [Python how-to guides](../HOW_TO/PYTHON/README.md)
- [Typed request reference](typed-requests.md)
- [Session and request compatibility](session-and-request-compatibility.md)
- [Output format and export reference](output-formats-and-export.md)
