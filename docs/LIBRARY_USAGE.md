[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)

# Library Usage (Python API)

This chapter explains how to use gbdraw as a Python library (not just the CLI).
It is useful for pipelines and notebooks.

## 1. Core workflow

A minimal flow looks like this:

1. Load input data into `SeqRecord`
2. (Optional) Load config (`config.toml`) to override defaults
3. (Optional) Load colors (palette / TSV) to override defaults
4. Assemble the SVG via the API
5. Save the SVG

Minimum required arguments:
- Circular: `gb_record`
- Linear: `records` (BLAST is optional; use `blast_files=None` or omit)

If `selected_features_set` is omitted, it defaults to:
`["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]`.
You can also import `DEFAULT_SELECTED_FEATURES` from `gbdraw.api`.

## 2. Minimal circular example

```python
from Bio import SeqIO

from gbdraw.api import assemble_circular_diagram_from_record
from gbdraw.api.render import save_figure

# 1) Load GenBank
record = next(SeqIO.parse("NC_000913.gbk", "genbank"))

# 2) Assemble (common options; defaults exist for config/colors/features)
canvas = assemble_circular_diagram_from_record(
    record,
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    output_prefix="ecoli_circular",
    legend="right",
)

# 3) Save (SVG/PNG/PDF)
save_figure(canvas, ["svg"])
```

### 2.1 Strict minimum (required args only)

```python
from Bio import SeqIO

from gbdraw.api import assemble_circular_diagram_from_record
from gbdraw.api.render import save_figure

record = next(SeqIO.parse("NC_000913.gbk", "genbank"))

canvas = assemble_circular_diagram_from_record(record)

save_figure(canvas, ["svg"])
```

### 2.2 Control tracks with `track_specs`

You can fine-tune visibility and placement of GC/skew/legend/etc. using
`track_specs` (list of strings or `TrackSpec`). Track specs are experimental
and currently only applied to circular diagrams.

```python
from gbdraw.api import parse_track_specs

track_specs = parse_track_specs(
    [
        "gc_skew@show=false",
        "legend@show=false",
        "features@w=0.12",
    ],
    mode="circular",
)

canvas = assemble_circular_diagram_from_record(
    record,
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    output_prefix="ecoli_circular_tracks",
    legend="right",
    track_specs=track_specs,
)
```

`features@w=...` accepts both factor and px units:
- `features@w=0.12` uses a factor of the circular radius.
- `features@w=48px` uses an absolute width in pixels.

### 2.3 Control labels/ticks with `track_specs`

```python
track_specs = parse_track_specs(
    [
        "labels@show=true",
        "ticks@show=false",
        "axis@show=true",
    ],
    mode="circular",
)

canvas = assemble_circular_diagram_from_record(
    record,
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    output_prefix="ecoli_circular_labels",
    legend="right",
    track_specs=track_specs,
)
```

## 3. Minimal linear example

```python
from gbdraw.api import assemble_linear_diagram_from_records
from gbdraw.api.io import load_gbks
from gbdraw.api.render import save_figure

records = load_gbks(["Genome1.gbk", "Genome2.gbk"], mode="linear", load_comparison=True)

canvas = assemble_linear_diagram_from_records(
    records,
    blast_files=["Genome1_Genome2.blast.out"],  # use None if you do not have BLAST
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    output_prefix="pairwise_linear",
    legend="right",
)

save_figure(canvas, ["svg"])
```

### 3.1 Strict minimum (required args only)

```python
from gbdraw.api import assemble_linear_diagram_from_records
from gbdraw.api.io import load_gbks
from gbdraw.api.render import save_figure

records = load_gbks(["Genome1.gbk"], mode="linear", load_comparison=False)

canvas = assemble_linear_diagram_from_records(
    records,
    blast_files=None,
)

save_figure(canvas, ["svg"])
```

### 3.2 Linear without BLAST

If you do not want comparison ribbons, pass `blast_files=None`.

```python
records = load_gbks(["Genome1.gbk"], mode="linear", load_comparison=False)

canvas = assemble_linear_diagram_from_records(
    records,
    blast_files=None,
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    output_prefix="linear_no_blast",
    legend="right",
)

save_figure(canvas, ["svg"])
```

### 3.3 GFF3 + FASTA input

Use `load_gff_fasta()` to load paired GFF3 + FASTA files.

```python
from gbdraw.api.io import load_gff_fasta

records = load_gff_fasta(
    gff_list=["Genome1.gff"],
    fasta_list=["Genome1.fna"],
    mode="linear",
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    load_comparison=False,
)

canvas = assemble_linear_diagram_from_records(
    records,
    blast_files=None,
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    output_prefix="linear_gff_fasta",
    legend="right",
)

save_figure(canvas, ["svg"])
```

### 3.4 Multiple records + region cropping

If your input files contain multiple records, you can select specific records
and crop regions.

```python
from gbdraw.api.io import parse_region_specs, apply_region_specs

records = load_gbks(
    ["Genome1.gbk", "Genome2.gbk"],
    mode="linear",
    load_comparison=False,
    record_selectors=["RecordA", "#0"],  # first by ID, second by index
)

region_specs = parse_region_specs([
    "RecordA:100000-250000",
    "#0:50000-120000:rc",
])

records = apply_region_specs(records, region_specs)

canvas = assemble_linear_diagram_from_records(
    records,
    blast_files=None,
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    output_prefix="linear_multi_region",
    legend="right",
)

save_figure(canvas, ["svg"])
```

## 4. Modify config programmatically

If you want CLI-like behavior, use `modify_config_dict()`.

```python
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml

config_dict = load_config_toml("gbdraw.data", "config.toml")
config_dict = modify_config_dict(
    config_dict,
    show_labels=True,
    track_type="middle",
    legend_font_size=18,
)
```

You can also pass `config_overrides` directly to the API (kwargs to `modify_config_dict`):

```python
canvas = assemble_circular_diagram_from_record(
    record,
    config_overrides={"strandedness": True, "show_labels": True},
)
```

Direct edits to `config_dict` are also fine.

## 5. Color tables

- `load_default_colors(...)`: load a palette
- `read_color_table(...)`: rule-based color overrides (`-t` equivalent)

```python
default_colors = load_default_colors("", palette="orchid")
color_table = read_color_table("custom_colors.tsv")
```

You can also pass palette/TSV paths directly to the API:

```python
canvas = assemble_circular_diagram_from_record(
    record,
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    default_colors_palette="orchid",
    default_colors_file="default_colors.tsv",
    color_table_file="custom_colors.tsv",
)
```

The same options are available for the linear API; `load_comparison` is inferred from `blast_files`.

## 6. API Reference

### 6.1 `assemble_circular_diagram_from_record(...)`

Required parameters:
- `gb_record` (`SeqRecord`): input record

Optional parameters:

| Name | Type | Default | Notes |
| --- | --- | --- | --- |
| `selected_features_set` | `Sequence[str] \| None` | CLI default list | Feature types to draw. Defaults to `["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]`. |
| `config_dict` | `dict \| None` | `None` | If `None`, loads built-in `config.toml`. |
| `config_overrides` | `Mapping[str, object] \| None` | `None` | Passed to `modify_config_dict`. Not compatible with `cfg`. |
| `color_table` | `DataFrame \| None` | `None` | Preloaded color table. Takes precedence over `color_table_file`. |
| `color_table_file` | `str \| None` | `None` | TSV color table path, loaded when `color_table` is `None`. |
| `default_colors` | `DataFrame \| None` | `None` | Preloaded palette. Takes precedence over `default_colors_*`. |
| `default_colors_palette` | `str` | `"default"` | Palette name in `color_palettes.toml`. |
| `default_colors_file` | `str \| None` | `None` | TSV overrides for the default palette. |
| `output_prefix` | `str` | `"out"` | Output name prefix (used in exported files). |
| `legend` | `str` | `"right"` | Legend placement (e.g. `"right"`, `"none"`). |
| `dinucleotide` | `str` | `"GC"` | GC/AT content target. |
| `window` | `int \| None` | `None` | Sliding window size; auto if `None`. |
| `step` | `int \| None` | `None` | Sliding step size; auto if `None`. |
| `species` | `str \| None` | `None` | Display label. |
| `strain` | `str \| None` | `None` | Display label. |
| `track_specs` | `Sequence[str \| TrackSpec] \| None` | `None` | Track controls. |
| `cfg` | `GbdrawConfig \| None` | `None` | Prebuilt config object. If provided, keep it consistent with `config_dict`. |

### 6.2 `assemble_linear_diagram_from_records(...)`

Required parameters:
- `records` (`Sequence[SeqRecord]`): input records (non-empty)
- `blast_files` (`Sequence[str] \| None`): BLAST outfmt 6/7 files, or `None`

Optional parameters:

| Name | Type | Default | Notes |
| --- | --- | --- | --- |
| `selected_features_set` | `Sequence[str] \| None` | CLI default list | Feature types to draw. Defaults to `["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]`. |
| `config_dict` | `dict \| None` | `None` | If `None`, loads built-in `config.toml`. |
| `config_overrides` | `Mapping[str, object] \| None` | `None` | Passed to `modify_config_dict`. Not compatible with `cfg`. |
| `color_table` | `DataFrame \| None` | `None` | Preloaded color table. Takes precedence over `color_table_file`. |
| `color_table_file` | `str \| None` | `None` | TSV color table path, loaded when `color_table` is `None`. |
| `default_colors` | `DataFrame \| None` | `None` | Preloaded palette. Takes precedence over `default_colors_*`. |
| `default_colors_palette` | `str` | `"default"` | Palette name in `color_palettes.toml`. |
| `default_colors_file` | `str \| None` | `None` | TSV overrides for the default palette. |
| `output_prefix` | `str` | `"out"` | Output name prefix (used in exported files). |
| `legend` | `str` | `"right"` | Legend placement (e.g. `"right"`, `"none"`). |
| `dinucleotide` | `str` | `"GC"` | GC/AT content target. |
| `window` | `int \| None` | `None` | Sliding window size; auto if `None`. |
| `step` | `int \| None` | `None` | Sliding step size; auto if `None`. |
| `evalue` | `float` | `1e-5` | BLAST filter. |
| `bitscore` | `float` | `50.0` | BLAST filter. |
| `identity` | `float` | `70.0` | BLAST filter (percent). |
| `cfg` | `GbdrawConfig \| None` | `None` | Prebuilt config object. If provided, keep it consistent with `config_dict`. |

Notes:
- If you pass both `color_table` and `color_table_file`, `color_table` is used.
- If you pass both `default_colors` and `default_colors_*`, `default_colors` is used.
- For linear diagrams, palette loading uses `load_comparison=bool(blast_files)`.

## 7. Exceptions (recommended pattern)

The library raises exceptions instead of calling `sys.exit()`.

```python
from gbdraw.exceptions import GbdrawError

try:
    # load / assemble / save
    ...
except GbdrawError as exc:
    print(f"gbdraw error: {exc}")
```

Common exceptions:
- `ConfigError`: config.toml missing or invalid
- `InputFileError`: input file missing or unreadable
- `ParseError`: input file parse errors
- `ValidationError`: validation errors

## 8. References

- CLI options: [CLI Reference](./CLI_Reference.md)
- Color palette examples: [color_palette_examples](../examples/color_palette_examples.md)

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)
