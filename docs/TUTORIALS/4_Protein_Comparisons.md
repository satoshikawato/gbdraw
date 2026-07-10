[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 3](./3_Advanced_Customization.md) | [Go to Tutorial 5 >](./5_Table_Driven_Inputs.md)

# Tutorial 4: Protein Comparisons Without Precomputed BLAST

**Goal:** build linear comparison plots from CDS-derived proteins without preparing BLAST tables yourself.

## 1. Prepare Annotated GenBank Inputs

The generated protein workflows need two or more annotated GenBank or GFF3 + FASTA records with CDS translations, or CDS features that can be translated.

This tutorial uses three majanivirus GenBank records:

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738868.1&rettype=gbwithparts&retmode=text" -O MjeNMV.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738874.1&rettype=gbwithparts&retmode=text" -O MelaMJNV.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738870.1&rettype=gbwithparts&retmode=text" -O PemoMJNVA.gb
```

If you are working from a source checkout, the same files are also available under `examples/`.

## 2. Runtime Selection

For `--protein_blastp_mode pairwise`, `orthogroup`, and `collinear`, gbdraw resolves the protein search runtime in this order:

1. bundled native LOSAT on Linux x86_64 when available
2. `losat` on `PATH`
3. NCBI BLAST+ `blastp` on `PATH`

Use explicit paths when you need to control the runtime:

```bash
gbdraw linear \
  --gbk MjeNMV.gb MelaMJNV.gb \
  --protein_blastp_mode pairwise \
  --losatp_bin /path/to/losat \
  --losatp_threads 4 \
  -o majani_pairwise_losat \
  -f svg
```

Or force NCBI BLAST+:

```bash
gbdraw linear \
  --gbk MjeNMV.gb MelaMJNV.gb \
  --protein_blastp_mode pairwise \
  --ncbi_blastp_bin /path/to/blastp \
  -o majani_pairwise_blastp \
  -f svg
```

NCBI BLAST+ output is compatible with the workflow, but it may not produce exactly the same hit set as LOSAT.

## 3. Pairwise Protein Ribbons

`pairwise` runs adjacent protein searches and draws pairwise ribbons from the resulting matches.

```bash
gbdraw linear \
  --gbk MjeNMV.gb MelaMJNV.gb \
  --protein_blastp_mode pairwise \
  --align_center \
  --pairwise_match_style curve \
  -o MjeNMV_MelaMJNV_pairwise_protein \
  -f svg
```

## 4. Orthogroup Ribbons

`orthogroup` groups related CDS-derived proteins across all input records before drawing adjacent display ribbons.

```bash
gbdraw linear \
  --gbk MjeNMV.gb MelaMJNV.gb PemoMJNVA.gb \
  --protein_blastp_mode orthogroup \
  --show_labels orthogroup_top \
  --pairwise_match_style curve \
  --align_center \
  -o majani_orthogroup \
  -f svg
```

`--show_labels orthogroup_top` labels the topmost displayed member of each orthogroup, which is useful when the same group appears in multiple records.

## 5. Collinear Blocks

`collinear` keeps protein-supported matches that occur in compatible local order.

```bash
gbdraw linear \
  --gbk MjeNMV.gb MelaMJNV.gb PemoMJNVA.gb \
  --protein_blastp_mode collinear \
  --collinear_min_anchors 2 \
  --collinear_color_mode orientation_identity \
  --pairwise_match_style curve \
  --align_center \
  -o majani_collinear \
  -f svg
```

`--collinear_min_anchors 2` removes singleton blocks. `--collinear_color_mode orientation_identity` separates forward and inverted blocks while still encoding identity.

## 6. When to Prefer Precomputed `-b/--blast`

Use precomputed BLAST tables when you need exact reproducibility from a specific BLAST version, custom database settings, nucleotide comparisons, translated nucleotide searches, or a workflow that has already filtered hits upstream.

Do not combine `-b/--blast` with `--protein_blastp_mode`. The CLI rejects that combination because the two options define different comparison sources.

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 3](./3_Advanced_Customization.md) | [Go to Tutorial 5 >](./5_Table_Driven_Inputs.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
