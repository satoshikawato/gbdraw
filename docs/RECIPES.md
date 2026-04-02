[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | **Recipes** | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)

# Recipes

Quick copy-paste solutions for common tasks. For explanations and screenshots, see the [Tutorials](./TUTORIALS/TUTORIALS.md).

## Circular Diagrams

### Basic circular plot

```bash
gbdraw circular --gbk genome.gb -o output -f svg
```

### Separate forward and reverse strands

```bash
gbdraw circular --gbk genome.gb -o output -f svg --separate_strands
```

### Use a built-in palette

```bash
gbdraw circular --gbk genome.gb -o output -f svg --separate_strands -p orchid
```

See [color_palette_examples.md](../examples/color_palette_examples.md) for the available palettes.

### Add centered organism text

```bash
gbdraw circular \
  --gbk genome.gb \
  -o output \
  -f svg \
  --species "<i>Escherichia coli</i>" \
  --strain "K-12"
```

### Add a top or bottom plot title

```bash
gbdraw circular \
  --gbk genome.gb \
  -o output \
  -f svg \
  --plot_title "Comparison overview" \
  --plot_title_position top
```

### Show labels

```bash
gbdraw circular --gbk genome.gb -o output -f svg --track_type middle --labels
```

### Show outer and inner labels

```bash
gbdraw circular --gbk genome.gb -o output -f svg --track_type middle --labels both
```

### Hide GC, skew, and legend

```bash
gbdraw circular --gbk genome.gb -o output -f svg --suppress_gc --suppress_skew --legend none
```

### Resolve overlapping features

```bash
gbdraw circular --gbk plasmid.gb -o output -f svg --resolve_overlaps
```

### Place multiple records on one shared canvas

```bash
gbdraw circular \
  --gbk file_with_multiple_records.gb \
  -o output \
  -f svg \
  --multi_record_canvas \
  --multi_record_size_mode auto
```

## Linear Diagrams

### Basic linear plot

```bash
gbdraw linear --gbk genome.gb -o output -f svg
```

### Show labels for all records

```bash
gbdraw linear --gbk genome.gb -o output -f svg --show_labels all
```

### Show labels only on the first record

```bash
gbdraw linear --gbk genome1.gb genome2.gb -o output -f svg --show_labels first
```

### Put tracks below the axis

```bash
gbdraw linear --gbk genome.gb -o output -f svg --track_layout below
```

### Use a ruler on the axis

```bash
gbdraw linear \
  --gbk genome.gb \
  -o output \
  -f svg \
  --track_layout below \
  --scale_style ruler \
  --ruler_on_axis
```

### Crop a region

```bash
gbdraw linear \
  --gbk genome.gb \
  --region Record1:101-300 \
  -o output \
  -f svg
```

### Reverse-complement selected inputs

```bash
gbdraw linear \
  --gbk genome1.gb genome2.gb \
  --reverse_complement false true \
  -o output \
  -f svg
```

## Comparative Genomics

### Two-genome comparison

```bash
blastn -query genome1.fasta -subject genome2.fasta -outfmt 7 -out blast.out

gbdraw linear \
  --gbk genome1.gb genome2.gb \
  -b blast.out \
  --align_center \
  --separate_strands \
  -o comparison \
  -f svg
```

### Multi-genome comparison

```bash
gbdraw linear \
  --gbk A.gb B.gb C.gb \
  -b A_vs_B.blast.out B_vs_C.blast.out \
  --align_center \
  --separate_strands \
  -o multi_comparison \
  -f svg
```

### Filter BLAST ribbons

```bash
gbdraw linear \
  --gbk genome1.gb genome2.gb \
  -b blast.out \
  --evalue 1e-99 \
  --bitscore 5000 \
  --identity 90 \
  --alignment_length 1000 \
  -o filtered \
  -f svg
```

## Color and Label Tables

### Override default feature colors

Create `modified_colors.tsv`:

```tsv
CDS	#d3d3d3
```

Then:

```bash
gbdraw circular --gbk genome.gb -d modified_colors.tsv -o output -f svg
```

### Highlight specific genes

Create `highlight.tsv`:

```tsv
CDS	product	polymerase	red	Polymerase
CDS	product	helicase	blue	Helicase
```

Then:

```bash
gbdraw circular --gbk genome.gb -t highlight.tsv -o output -f svg
```

### Prefer gene names over product labels

Create `priority.tsv`:

```tsv
CDS	gene
```

Then:

```bash
gbdraw circular --gbk genome.gb --labels --qualifier_priority priority.tsv -o output -f svg
```

### Override label text after filtering

Create `label_override.tsv`:

```tsv
*	*	label	^hypothetical protein$	HP
```

Then:

```bash
gbdraw circular --gbk genome.gb --labels --label_table label_override.tsv -o output -f svg
```

## GFF3 Input

### Circular plot from GFF3 + FASTA

```bash
gbdraw circular --gff annotations.gff --fasta sequence.fasta -o output -f svg
```

### Linear plot from GFF3 + FASTA

```bash
gbdraw linear --gff annotations.gff --fasta sequence.fasta -o output -f svg
```

## Output Formats

### Export to multiple formats

```bash
gbdraw circular --gbk genome.gb -o output -f svg,png,pdf
```

PNG, PDF, EPS, and PS require CairoSVG to be installed.

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | **Recipes** | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)
