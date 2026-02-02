[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)

# Library Usage (Python API)

This chapter explains how to use gbdraw as a Python library (not just the CLI).
It is useful for pipelines and notebooks.

## 1. Core workflow

A minimal flow looks like this:

1. Load input data into `SeqRecord`
2. Load config (`config.toml`)
3. Load colors (palette / TSV)
4. Assemble the SVG via the API
5. Save the SVG

## 2. Minimal circular example

```python
from Bio import SeqIO

from gbdraw.api import assemble_circular_diagram_from_record
from gbdraw.config.toml import load_config_toml
from gbdraw.io.colors import load_default_colors, read_color_table
from gbdraw.render.export import save_figure

# 1) Load GenBank
record = next(SeqIO.parse("NC_000913.gbk", "genbank"))

# 2) Config and colors
config_dict = load_config_toml("gbdraw.data", "config.toml")
default_colors = load_default_colors(user_defined_default_colors="", palette="default")
color_table = read_color_table("")  # empty string means no table

# 3) Assemble
canvas = assemble_circular_diagram_from_record(
    record,
    config_dict=config_dict,
    color_table=color_table,
    default_colors=default_colors,
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    output_prefix="ecoli_circular",
    legend="right",
)

# 4) Save (SVG/PNG/PDF)
save_figure(canvas, ["svg"])
```

### 2.1 Control tracks with `track_specs`

You can fine-tune visibility and placement of GC/skew/legend/etc. using
`track_specs` (list of strings or `TrackSpec`).

```python
from gbdraw.api import parse_track_specs

track_specs = parse_track_specs(
    [
        "gc_skew@show=false",
        "legend@show=false",
        "features@r=0.82,w=0.12",
    ],
    mode="circular",
)

canvas = assemble_circular_diagram_from_record(
    record,
    config_dict=config_dict,
    color_table=color_table,
    default_colors=default_colors,
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    output_prefix="ecoli_circular_tracks",
    legend="right",
    track_specs=track_specs,
)
```

### 2.2 Control labels/ticks with `track_specs`

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
    config_dict=config_dict,
    color_table=color_table,
    default_colors=default_colors,
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    output_prefix="ecoli_circular_labels",
    legend="right",
    track_specs=track_specs,
)
```

## 3. Minimal linear example

```python
from gbdraw.api import assemble_linear_diagram_from_records
from gbdraw.config.toml import load_config_toml
from gbdraw.io.colors import load_default_colors, read_color_table
from gbdraw.io.genome import load_gbks
from gbdraw.render.export import save_figure

records = load_gbks(["Genome1.gbk", "Genome2.gbk"], mode="linear", load_comparison=True)
config_dict = load_config_toml("gbdraw.data", "config.toml")
default_colors = load_default_colors(user_defined_default_colors="", palette="default", load_comparison=True)
color_table = read_color_table("")

canvas = assemble_linear_diagram_from_records(
    records,
    blast_files=["Genome1_Genome2.blast.out"],  # use None if you do not have BLAST
    config_dict=config_dict,
    color_table=color_table,
    default_colors=default_colors,
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    output_prefix="pairwise_linear",
    legend="right",
)

save_figure(canvas, ["svg"])
```

### 3.1 Linear without BLAST

If you do not want comparison ribbons, pass `blast_files=None`.

```python
records = load_gbks(["Genome1.gbk"], mode="linear", load_comparison=False)
config_dict = load_config_toml("gbdraw.data", "config.toml")
default_colors = load_default_colors(user_defined_default_colors="", palette="default", load_comparison=False)
color_table = read_color_table("")

canvas = assemble_linear_diagram_from_records(
    records,
    blast_files=None,
    config_dict=config_dict,
    color_table=color_table,
    default_colors=default_colors,
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    output_prefix="linear_no_blast",
    legend="right",
)

save_figure(canvas, ["svg"])
```

### 3.2 GFF3 + FASTA input

Use `load_gff_fasta()` to load paired GFF3 + FASTA files.

```python
from gbdraw.io.genome import load_gff_fasta

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
    config_dict=config_dict,
    color_table=color_table,
    default_colors=default_colors,
    selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"],
    output_prefix="linear_gff_fasta",
    legend="right",
)

save_figure(canvas, ["svg"])
```

### 3.3 Multiple records + region cropping

If your input files contain multiple records, you can select specific records
and crop regions.

```python
from gbdraw.io.regions import parse_region_specs, apply_region_specs

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
    config_dict=config_dict,
    color_table=color_table,
    default_colors=default_colors,
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

config_dict = load_config_toml("gbdraw.data", "config.toml")
config_dict = modify_config_dict(
    config_dict,
    show_labels=True,
    track_type="middle",
    legend_font_size=18,
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

## 6. Exceptions (recommended pattern)

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

## 7. References

- CLI options: [CLI Reference](./CLI_Reference.md)
- Color palette examples: [color_palette_examples](../examples/color_palette_examples.md)

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)
