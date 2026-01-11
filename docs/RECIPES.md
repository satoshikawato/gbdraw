[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | **Recipes** | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)

# Recipes

Quick copy-paste solutions for common tasks. For detailed explanations, see the [Tutorials](./TUTORIALS/TUTORIALS.md).

---

## Circular Diagrams

### Basic Circular Plot
```bash
gbdraw circular --gbk genome.gb -o output -f svg
```

### Circular Plot with Separated Strands
```bash
gbdraw circular --gbk genome.gb -o output -f svg --separate_strands
```

### Circular Plot with Custom Color Palette
```bash
gbdraw circular --gbk genome.gb -o output -f svg --separate_strands -p orchid
```

See [color_palette_examples.md](../examples/color_palette_examples.md) for all 55+ available palettes.

### Circular Plot with Title
```bash
gbdraw circular --gbk genome.gb -o output -f svg --separate_strands \
  --species "<i>Escherichia coli</i>" --strain "K-12"
```

### Circular Plot with Labels (Small Genomes)
```bash
gbdraw circular --gbk genome.gb -o output -f svg --track_type middle --show_labels
```

### Circular Plot with Inner Labels (Mitochondria/Plasmids)
```bash
gbdraw circular --gbk genome.gb -o output -f svg --track_type middle \
  --show_labels --allow_inner_labels
```

### Minimal Circular Plot (No GC/Skew/Legend)
```bash
gbdraw circular --gbk genome.gb -o output -f svg \
  --suppress_gc --suppress_skew -l none
```

### Resolve Overlapping Features (Plasmids)
```bash
gbdraw circular --gbk plasmid.gb -o output -f svg --resolve_overlaps
```

---

## Linear Diagrams

### Basic Linear Plot
```bash
gbdraw linear --gbk genome.gb -o output -f svg
```

### Linear Plot with Separated Strands
```bash
gbdraw linear --gbk genome.gb -o output -f svg --separate_strands
```

### Linear Plot with GC Content
```bash
gbdraw linear --gbk genome.gb -o output -f svg --show_gc
```

### Linear Plot with GC Content and Skew
```bash
gbdraw linear --gbk genome.gb -o output -f svg --show_gc --show_skew
```

### Linear Plot with Ruler Scale
```bash
gbdraw linear --gbk genome.gb -o output -f svg --scale_style ruler
```

---

## Comparative Genomics

### Two-Genome Comparison
```bash
# Step 1: Run BLAST
blastn -query genome1.fasta -subject genome2.fasta -outfmt 6 -out blast.out

# Step 2: Generate comparison plot
gbdraw linear --gbk genome1.gb genome2.gb -b blast.out -o comparison -f svg \
  --align_center --separate_strands
```

### Multi-Genome Comparison
```bash
# For genomes A, B, C - need BLAST for consecutive pairs: A-B and B-C
gbdraw linear --gbk A.gb B.gb C.gb \
  -b A_vs_B.blast.out B_vs_C.blast.out \
  -o multi_comparison -f svg --align_center --separate_strands
```

### Filter BLAST Matches
```bash
gbdraw linear --gbk genome1.gb genome2.gb -b blast.out -o filtered -f svg \
  --evalue 1e-99 --bitscore 5000 --identity 90
```

---

## Color Customization

### Override Default CDS Color
Create `modified_colors.tsv`:
```
CDS	#d3d3d3
```
Then:
```bash
gbdraw circular --gbk genome.gb -d modified_colors.tsv -o output -f svg
```

### Highlight Specific Genes
Create `highlight.tsv`:
```
CDS	product	polymerase	red	Polymerase
CDS	product	helicase	blue	Helicase
```
Then:
```bash
gbdraw circular --gbk genome.gb -t highlight.tsv -o output -f svg
```

---

## Label Control

### Hide "Hypothetical Protein" Labels
```bash
gbdraw circular --gbk genome.gb --show_labels \
  --label_blacklist "hypothetical protein" -o output -f svg
```

### Show Only Specific Genes
Create `whitelist.tsv`:
```
CDS	gene	dnaA
CDS	gene	dnaB
CDS	product	DNA polymerase
```
Then:
```bash
gbdraw circular --gbk genome.gb --show_labels \
  --label_whitelist whitelist.tsv -o output -f svg
```

### Use Gene Names Instead of Product
Create `priority.tsv`:
```
CDS	gene
```
Then:
```bash
gbdraw circular --gbk genome.gb --show_labels \
  --qualifier_priority priority.tsv -o output -f svg
```

---

## Output Formats

### Export to Multiple Formats
```bash
gbdraw circular --gbk genome.gb -o output -f svg,png,pdf
```

> **Note:** PNG, PDF, EPS, and PS formats require CairoSVG to be installed.

---

## GFF3 Input

### Using GFF3 Instead of GenBank
```bash
gbdraw circular --gff annotations.gff --fasta sequence.fasta -o output -f svg
```

---

## Track Layout Options

### Tuckin (Default)
Features inside the axis circle, GC/skew tracks in the center.
```bash
gbdraw circular --gbk genome.gb --track_type tuckin -o output -f svg
```

### Middle
Features on the axis circle (ideal for labels).
```bash
gbdraw circular --gbk genome.gb --track_type middle -o output -f svg
```

### Spreadout
Features outside the axis circle.
```bash
gbdraw circular --gbk genome.gb --track_type spreadout -o output -f svg
```

---

## Large Genomes

### Bacterial Genome (1-10 Mb)
Window/step sizes auto-adjust. Override if needed:
```bash
gbdraw circular --gbk bacteria.gb -o output -f svg \
  -w 10000 -s 1000
```

### Large Genome (>10 Mb)
```bash
gbdraw circular --gbk large.gb -o output -f svg \
  -w 100000 -s 10000
```

---

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | **Recipes** | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)
