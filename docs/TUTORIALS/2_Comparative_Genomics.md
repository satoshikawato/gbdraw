[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 1](./1_Customizing_Plots.md) | [Go to Tutorial 3 >](./3_Advanced_Customization.md)

# Tutorial 2: Comparative Genomics with BLAST

**Goal:** use `gbdraw linear` to visualize sequence similarity between genomes with precomputed BLAST tables or generated protein blastp comparisons.

## 1. Required Inputs

Comparative plots need two or more annotated genomes in GenBank or GFF3 + FASTA format.

For precomputed comparisons, supply one nucleotide or protein BLAST result file for each adjacent comparison with `-b/--blast`. Accepted BLAST table formats are outfmt 6 and outfmt 7.

For generated protein comparisons, do not supply BLAST tables. Use `--protein_blastp_mode pairwise`, `orthogroup`, or `collinear`; gbdraw extracts proteins from CDS features and runs the selected blastp workflow.

> [!IMPORTANT]
> `-b/--blast` and `--protein_blastp_mode` are mutually exclusive. Choose one comparison source per run: either precomputed BLAST tables or generated protein blastp comparisons.

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

The default pairwise link style is `ribbon`, which draws straight filled ribbons and is best when the exact BLAST alignment span is the main signal. Use `--pairwise_match_style curve` for curved filled ribbons in dense synteny-style views; the curve style still preserves each match span from `qstart/qend` and `sstart/send`.

For generated collinear protein comparisons, `--collinear_color_mode orientation_identity` uses separate forward and inverted identity gradients. Collinear blocks use RBH anchors. In the web UI, Evidence scope controls which record pairs provide collinearity evidence; Adjacent pairs produces local collinear gene groups, while All records can provide collinearity-backed global metadata. Ribbons are still drawn for adjacent display pairs.

## 4. Generate Protein Comparisons Without BLAST Tables

The same two GenBank files can be compared by generated protein blastp workflows. These commands intentionally omit `-b/--blast`.

Pairwise protein ribbons:

```bash
gbdraw linear \
  --gbk Escherichia_coli.gbk Shigella_dysenteriae.gbk \
  --protein_blastp_mode pairwise \
  --align_center \
  -o Escherichia_Shigella_protein_pairwise \
  -f svg
```

Orthogroup-supported ribbons:

```bash
gbdraw linear \
  --gbk Escherichia_coli.gbk Shigella_dysenteriae.gbk \
  --protein_blastp_mode orthogroup \
  --show_labels orthogroup_top \
  --pairwise_match_style curve \
  -o Escherichia_Shigella_orthogroup \
  -f svg
```

Collinear blocks:

```bash
gbdraw linear \
  --gbk Escherichia_coli.gbk Shigella_dysenteriae.gbk \
  --protein_blastp_mode collinear \
  --collinear_min_anchors 2 \
  --collinear_color_mode orientation_identity \
  --pairwise_match_style curve \
  -o Escherichia_Shigella_collinear \
  -f svg
```

See [Tutorial 4](./4_Protein_Comparisons.md) for runtime selection, labels, and collinear tuning.

## 5. Interpreting Orthogroups and Collinear Blocks

`orthogroup` groups related CDS-derived proteins across the input records, while `collinear` groups those protein-supported matches when their genes occur in a compatible order.

Both modes run protein `blastp` searches on proteins obtained from CDS features. gbdraw uses bundled native LOSAT when available, then a `losat` executable on `PATH`, then NCBI BLAST+ `blastp` on `PATH`. Linux x86_64 can use the bundled LOSAT binary. macOS and Windows do not currently ship bundled LOSAT binaries, so install NCBI BLAST+ and make `blastp` available on `PATH`, or pass it explicitly with `--ncbi_blastp_bin /path/to/blastp`. You can still force a native LOSAT executable on any platform with `--losatp_bin /path/to/losat`. NCBI BLAST+ fallback provides compatible outfmt 6 protein comparisons, but it is not guaranteed to produce exactly the same hit set as LOSAT.

Choose one comparison source per run. Precomputed `-b/--blast` tables are used directly as pairwise comparison input; generated `--protein_blastp_mode` runs infer pairwise protein matches, orthogroups, or collinear blocks from CDS proteins.

### Orthogroups

In evolutionary genomics, an orthogroup is a set of genes descended from a single gene in the last common ancestor of the taxa being compared. Different genomes may contribute one gene, multiple genes, or no detectable member. An orthogroup is therefore not necessarily a set of one-to-one orthologs.

> [!NOTE]
> gbdraw infers similarity-based groups for visualization. Reciprocal and near-reciprocal protein matches define the group cores, and additional proteins are assigned when the evidence supports membership in the same family. Because gbdraw does not infer gene or species trees, these groups may not necessarily reflect true phylogenetic orthology. Use a dedicated orthology workflow such as OrthoFinder when orthology inference itself is the analysis goal.

### Collinear Blocks

An anchor is a protein-supported gene pair associated with the same gbdraw orthogroup. A collinear block is a run of anchors with compatible order in two records. For multi-anchor blocks, `plus` means that the anchors occur in the same order in both records; `minus` means that their order is reversed and suggest an inversion.

Multi-anchor blocks combine protein similarity with conserved local gene order. They can highlight conserved gene neighborhoods, including candidate operons or gene clusters, but they do not by themselves establish shared function or cotranscription. Genes lying between anchors are not automatically homologous or members of the same orthogroup.

The default `--collinear_min_anchors 1` retains singleton links. Set it to `2` to require at least two anchors per rendered block; higher values require more anchors.

## 6. Select Records or Regions

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

For larger comparisons, put row-specific selectors and crops in a records table instead of repeating several order-sensitive CLI options.

```tsv
gbk	record_label	record_id	region	reverse_complement	order
Genome1.gbk	Genome 1	#1		0	1
Genome2.gbk	Genome 2	#1	50000-180000	1	2
```

```bash
gbdraw linear \
  --records_table records.tsv \
  -b Genome1_vs_Genome2.blast.outfmt7.txt \
  -o Genome1_Genome2_table \
  -f svg
```

`--records_table` is an alternative input source, so do not combine it with `--gbk`, `--gff`, or `--fasta`.

## 7. Compare More Than Two Genomes

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
  --alignment_length 1000 \
  -o Escherichia_Shigella_multi \
  -f svg
```

Use `--alignment_length` when you want to hide very short BLAST hits and keep only longer ribbons.

![Escherichia_Shigella_multi.svg](../../examples/Escherichia_Shigella_multi.svg)

> [!IMPORTANT]
> BLAST files must follow the same order as the genome list.

## 8. Circular Conservation Rings

Circular mode can draw BLAST HSPs as concentric conservation rings around each displayed record. This is useful when one annotated reference genome should be compared against one or more unannotated FASTA sequences.

Generate BLAST with the displayed circular genome as the subject:

```bash
blastn \
  -query Shigella_dysenteriae.fasta \
  -subject Escherichia_coli.fasta \
  -outfmt 7 \
  -out Shigella_vs_Escherichia.blastn.out
```

Then render the circular reference with a conservation ring:

```bash
gbdraw circular \
  --gbk Escherichia_coli.gbk \
  --conservation_blast Shigella_vs_Escherichia.blastn.out \
  --conservation_reference subject \
  --conservation_labels "Shigella dysenteriae" \
  --identity 75 \
  --alignment_length 500 \
  -o Escherichia_conservation \
  -f svg
```

Use one `--conservation_blast` file per ring. If gbdraw can identify the displayed reference IDs on only one BLAST side, `--conservation_reference auto` works; otherwise set `query` or `subject` explicitly. BLAST rows with `start > end` on the selected reference side are drawn as reverse-orientation hits, not as circular wraparound hits.

In the web app, circular mode also supports Conservation Rings from uploaded BLAST outfmt 6/7 files or browser LOSAT `blastn`. LOSAT mode uses each comparison FASTA as the query and the displayed circular genome as the subject, so the generated BLAST results are passed to gbdraw with `conservation_reference=subject`.

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 1](./1_Customizing_Plots.md) | [Go to Tutorial 3 >](./3_Advanced_Customization.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
