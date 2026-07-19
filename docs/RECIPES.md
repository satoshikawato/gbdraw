[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | **Recipes** | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)

# Recipes

Copy-paste command examples for common tasks. For explanations and screenshots, see the [command-line guides](./TUTORIALS/TUTORIALS.md).

## Circular diagrams

### Basic circular genome diagram

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

### Add a center label

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

`--track_type` selects a circular preset. Custom track slots can reorder tracks while inheriting omitted geometry from that preset; the circular axis stays fixed and is not configurable as a slot.

### Reorder circular tracks without calculating radii

```bash
gbdraw circular \
  --gbk genome.gb \
  -o output \
  -f svg \
  --track_type middle \
  --circular_track_order features,ticks,gc_skew,gc_content
```

The order list above swaps GC skew and GC content while keeping the record-specific `middle` preset radii, widths, spacing, and tick defaults unless a slot explicitly overrides them.

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

## Linear diagrams

### Basic linear genome diagram

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
  --reverse_complement false \
  --reverse_complement true \
  -o output \
  -f svg
```

### Use a records table for row-specific inputs

Create `records.tsv`:

```tsv
gbk	record_label	record_id	region	reverse_complement	order
genome1.gb	Genome 1	#1		0	1
genome2.gb	Genome 2	#1	101-20000	1	2
```

Then run:

```bash
gbdraw linear \
  --records_table records.tsv \
  -o output \
  -f svg
```

See [Use TSV manifests for CLI inputs](./TUTORIALS/5_Table_Driven_Inputs.md) for examples covering GenBank, GFF3+FASTA, circular placement, BLAST similarity rings, and track slots.

### Highlight named regions

Create `annotations.tsv`:

```tsv
set_id	id	mark	start	end	label	fill	fill_opacity	legend_label
review	window	band	1000	5000	Review window	#f59e0b	0.25	Review region
```

```bash
gbdraw circular \
  --gbk genome.gb \
  --annotation_table annotations.tsv \
  --circular_track_slot review:annotations@set_id=review,side=outside,w=28px \
  --circular_track_slot features:features@side=overlay \
  --circular_track_slot ticks:ticks@side=inside \
  -o annotated \
  -f svg
```

For a Linear row, replace the circular slots with `--linear_track_slot review:annotations@set_id=review,side=above,h=28px` and a feature slot.

## Comparative genomics

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

### Selected comparisons across multi-record rows

```bash
gbdraw linear \
  --records_table examples/linear_multi_records.tsv \
  --comparisons_table examples/linear_multi_comparisons.tsv \
  --linear_record_gap 28 \
  --identity 97 \
  --alignment_length 500 \
  --scale_style ruler \
  --ruler_on_axis \
  -o linear_multi_record \
  -f svg
```

The records table assigns `row` and `column`; the comparisons table declares `blast`, `query`, and `subject`. Use this form whenever a row contains more than one record. `-b/--blast` remains the shorter adjacency-based form for one-record-per-row layouts.

### Run a protein comparison from CDS annotations

```bash
gbdraw linear \
  --gbk genome1.gb genome2.gb genome3.gb \
  --protein_blastp_mode orthogroup \
  --show_labels orthogroup_top \
  --pairwise_match_style curve \
  -o protein_orthogroup \
  -f svg
```

Use `--protein_blastp_mode pairwise`, `orthogroup`, or `collinear`. The `orthogroup` mode creates gbdraw similarity groups for visualization; it does not infer phylogeny-based orthogroups. Do not combine these modes with `-b/--blast`. See [Draw protein matches from annotated CDS features](./TUTORIALS/4_Protein_Comparisons.md) for runtime selection and collinear examples.

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

### Circular BLAST similarity rings

```bash
blastn -query comparison.fasta -subject reference.fasta -outfmt 7 -out comparison_vs_reference.blast.out

gbdraw circular \
  --gbk reference.gb \
  --conservation_blast comparison_vs_reference.blast.out \
  --conservation_reference subject \
  --conservation_labels "Comparison" \
  --identity 75 \
  --alignment_length 500 \
  -o circular_conservation \
  -f svg
```

Pass more `--conservation_blast` files to add rings. Each ring shows raw BLAST HSP spans, not an inferred measure of evolutionary conservation. Reverse-coordinate BLAST rows are drawn as reverse hits rather than circular wraparound hits.

For a standalone interactive SVG with both matched spans available, add one comparison FASTA per BLAST source:

```bash
gbdraw circular \
  --gbk reference.gbk \
  --conservation_blast comparison_vs_reference.blast.out \
  --conservation_fasta comparison.fna \
  --conservation_reference subject \
  -f interactive_svg \
  -o circular_comparison
```

## Color and label tables

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

## GFF3 input

### Circular plot from GFF3 + FASTA

```bash
gbdraw circular --gff annotations.gff --fasta sequence.fasta -o output -f svg
```

### Linear plot from GFF3 + FASTA

```bash
gbdraw linear --gff annotations.gff --fasta sequence.fasta -o output -f svg
```

## Output formats

### Export to multiple formats

```bash
gbdraw circular --gbk genome.gb -o output -f svg,png,pdf
```

PNG, PDF, EPS, and PS require CairoSVG to be installed.

### Export standalone interactive SVG

```bash
gbdraw circular --gbk genome.gb -o output -f interactive_svg
gbdraw linear --gbk genome1.gb genome2.gb -o output -f svg,interactive_svg
```

`interactive_svg` writes the normal `output.svg` plus `output.interactive.svg`. The older spelling `interactive-svg` remains accepted as a compatibility alias.
It does not require CairoSVG, Node.js, Playwright, Chromium, or a web build step.
Open the interactive file in a browser; some desktop SVG viewers block embedded scripts.

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | **Recipes** | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)
