[Home](./DOCS.md) | [Installation](./INSTALL.md) | **Quickstart** | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)

[< Back to Installation](./INSTALL.md) | [Go to Tutorials >](./TUTORIALS/TUTORIALS.md)

# Quickstart

This quickstart shows how to draw a circular genome diagram with the CLI. If you prefer a GUI, use [https://gbdraw.app/](https://gbdraw.app/) or run `gbdraw gui` locally after installation.

## 1. Confirm the installation

Make sure `gbdraw` is available in your environment:

```bash
gbdraw -h
```

## 2. Download a sample GenBank file

This example uses the *Escherichia coli* K-12 reference genome.

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=gbwithparts&retmode=text" -O NC_000913.gbk
```

## 3. Generate a circular genome diagram

```bash
gbdraw circular \
  --gbk NC_000913.gbk \
  -o ecoli_k12_plot \
  -f svg \
  --separate_strands
```

This command:

- reads `NC_000913.gbk`
- writes `ecoli_k12_plot.svg`
- draws a circular genome diagram
- separates forward and reverse strands into different feature tracks

## 4. Inspect the output

Open the resulting SVG in a browser or vector editor.

![Circular E. coli K-12 genome diagram with separate forward and reverse strand tracks](../examples/ecoli_k12_plot.svg)

## 5. Optional: add feature labels or a center label

Show labels for a focused set of RNA features. Restricting this large *E. coli* genome to rRNA and tRNA features keeps the labels readable:

```bash
gbdraw circular \
  --gbk NC_000913.gbk \
  -o ecoli_k12_labeled \
  -f svg \
  --track_type middle \
  --features rRNA,tRNA \
  --labels
```

![Circular E. coli K-12 genome with labeled rRNA and tRNA features](../examples/quickstart-labeled-rna-features.svg)

Add a species and strain to the center of the diagram:

```bash
gbdraw circular \
  --gbk NC_000913.gbk \
  -o ecoli_with_title \
  -f svg \
  --separate_strands \
  --species "<i>Escherichia coli</i>" \
  --strain "K-12"
```

![Circular E. coli K-12 genome diagram with an italic species name and strain in the center](../examples/ecoli_with_title.svg)

## 6. Optional: linear mode selectors

Linear mode can target specific records or regions:

- `--record_id`: select a record by ID or `#index`
- `--reverse_complement`: reverse-complement per input file
- `--region`: crop a region with `record_id:start-end[:rc]`

Example:

```bash
gbdraw linear \
  --gbk NC_000913.gbk \
  --record_id NC_000913.3 \
  --region NC_000913.3:100000-250000 \
  -o ecoli_linear_region \
  -f svg
```

![Linear E. coli K-12 diagram cropped to positions 100000 through 250000](../examples/quickstart-linear-region.svg)

If you need an index selector in the shell, quote it:

```bash
gbdraw linear \
  --gbk Genome1.gbk Genome2.gbk \
  --record_id Genome1_Chr1 \
  --record_id '#0' \
  --reverse_complement false \
  --reverse_complement true \
  -o genome_pair_selected \
  -f svg
```

## Next steps

- Continue with [Style a circular genome diagram](./TUTORIALS/1_Customizing_Plots.md)
- Use [Draw protein matches from annotated CDS features](./TUTORIALS/4_Protein_Comparisons.md) to run protein searches during diagram generation
- Use [TSV manifests for CLI inputs](./TUTORIALS/5_Table_Driven_Inputs.md) when records need separate labels, selectors, crops, or orientations
- Browse [Recipes](./RECIPES.md) for common command patterns
- Use the [CLI Reference](./CLI_Reference.md) for the full option list
- Explore more figures in the [Gallery](./GALLERY.md)

[< Back to Installation](./INSTALL.md) | [Go to Tutorials >](./TUTORIALS/TUTORIALS.md)

[Home](./DOCS.md) | [Installation](./INSTALL.md) | **Quickstart** | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)
