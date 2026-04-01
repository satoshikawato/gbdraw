[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 1](./1_Customizing_Plots.md) | [Go to Tutorial 3 >](./3_Advanced_Customization.md)

# Tutorial 2: Comparative Genomics with BLAST

**Goal:** use `gbdraw linear` to visualize sequence similarity between genomes with BLAST result files.

## 1. Required Inputs

Comparative plots need:

1. Two or more annotated genomes in GenBank or GFF3 + FASTA format
2. One BLAST result file for each adjacent comparison

Accepted BLAST table formats are outfmt 6 and outfmt 7.

## 2. Prepare a Pairwise Example

This example compares *Escherichia coli* and *Shigella dysenteriae*.

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=gbwithparts&retmode=text" -O Escherichia_coli.gbk
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text" -O Escherichia_coli.fasta

wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NZ_CP055055.1&rettype=gbwithparts&retmode=text" -O Shigella_dysenteriae.gbk
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NZ_CP055055.1&rettype=fasta&retmode=text" -O Shigella_dysenteriae.fasta
```

Run BLAST:

```bash
blastn \
  -query Escherichia_coli.fasta \
  -subject Shigella_dysenteriae.fasta \
  -outfmt 7 \
  -out Escherichia_coli-Shigella_dysenteriae.blastn.out
```

## 3. Generate the Pairwise Plot

```bash
gbdraw linear \
  --gbk Escherichia_coli.gbk Shigella_dysenteriae.gbk \
  -b Escherichia_coli-Shigella_dysenteriae.blastn.out \
  --align_center \
  --separate_strands \
  -o Escherichia_Shigella_pair \
  -f svg
```

![Escherichia_Shigella_pair.svg](../../examples/Escherichia_Shigella_pair.svg)

## 4. Select Records or Regions

Linear mode supports three selectors:

- `--record_id`: choose a record by ID or `#index` for each input file
- `--reverse_complement`: reverse-complement selected inputs
- `--region`: crop a region with `record_id:start-end[:rc]`

Example: select one record and crop a region.

```bash
gbdraw linear \
  --gbk Escherichia_coli.gbk \
  --record_id NC_000913.3 \
  --region NC_000913.3:100000-250000 \
  -o Escherichia_coli_region \
  -f svg
```

Example: reverse-complement the second input.

```bash
gbdraw linear \
  --gbk Genome1.gbk Genome2.gbk \
  --reverse_complement false true \
  -o Genome1_Genome2_rc \
  -f svg
```

If you use `--region` or `--reverse_complement` together with BLAST, make sure the BLAST coordinates still match the selected inputs.

## 5. Compare More Than Two Genomes

For `A -> B -> C -> D`, provide BLAST files for `A vs B`, `B vs C`, and `C vs D`.

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_004337.2&rettype=gbwithparts&retmode=text" -O Shigella_flexneri.gbk
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_004337.2&rettype=fasta&retmode=text" -O Shigella_flexneri.fasta

wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NZ_CP026802.1&rettype=gbwithparts&retmode=text" -O Shigella_sonnei.gbk
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NZ_CP026802.1&rettype=fasta&retmode=text" -O Shigella_sonnei.fasta
```

```bash
blastn -query Shigella_dysenteriae.fasta -subject Shigella_flexneri.fasta -outfmt 7 -out Shigella_dysenteriae-Shigella_flexneri.blastn.out
blastn -query Shigella_flexneri.fasta -subject Shigella_sonnei.fasta -outfmt 7 -out Shigella_flexneri-Shigella_sonnei.blastn.out
```

```bash
gbdraw linear \
  --gbk Escherichia_coli.gbk Shigella_dysenteriae.gbk Shigella_flexneri.gbk Shigella_sonnei.gbk \
  -b Escherichia_coli-Shigella_dysenteriae.blastn.out Shigella_dysenteriae-Shigella_flexneri.blastn.out Shigella_flexneri-Shigella_sonnei.blastn.out \
  --align_center \
  --separate_strands \
  --evalue 1e-99 \
  --bitscore 5000 \
  -o Escherichia_Shigella_multi \
  -f svg
```

![Escherichia_Shigella_multi.svg](../../examples/Escherichia_Shigella_multi.svg)

> [!IMPORTANT]
> BLAST files must follow the same order as the genome list.

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 1](./1_Customizing_Plots.md) | [Go to Tutorial 3 >](./3_Advanced_Customization.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
