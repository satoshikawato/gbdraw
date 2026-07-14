[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 1](./1_Customizing_Plots.md) | [Go to Tutorial 3 >](./3_Advanced_Customization.md)

# Tutorial 2: Compare genomes with BLAST and protein searches

Use precomputed BLAST tables or protein similarity searches from translated CDS features to draw linear links. The final section compares LOSATN and TLOSATX hits in circular mode.

## 1. Required inputs

Linear comparative plots require two or more annotated genomes in GenBank or GFF3 + FASTA format.

For precomputed searches, supply one tabular nucleotide or protein BLAST output file for each adjacent record pair with `-b/--blast`. gbdraw accepts outfmt 6 and outfmt 7.

To run a protein search during diagram generation, omit `-b/--blast` and select `--protein_blastp_mode pairwise`, `orthogroup`, or `collinear`. gbdraw obtains amino acid sequences from CDS features and runs the selected workflow.

> [!IMPORTANT]
> `-b/--blast` and `--protein_blastp_mode` are mutually exclusive. Choose one source per run: precomputed tabular BLAST output or a protein search run by gbdraw.

## 2. Prepare a pairwise example

This example uses the first two sequences from the [Hepatoplasmataceae five-genome comparison](../GALLERY.md#hepatoplasmataceae-five-genome-comparison): *Candidatus Tyloplasma litorale* (AP027078.1) and *Candidatus Hepatoplasma vulgare* (AP027131.1).

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AP027078.1&rettype=gbwithparts&retmode=text" -O AP027078.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AP027078.1&rettype=fasta&retmode=text" -O AP027078.fasta

wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AP027131.1&rettype=gbwithparts&retmode=text" -O AP027131.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AP027131.1&rettype=fasta&retmode=text" -O AP027131.fasta
```

Install NCBI BLAST+ so that `tblastx` is available on `PATH`, then run a translated nucleotide search and save the outfmt 7 results:

```bash
tblastx \
  -query AP027078.fasta \
  -subject AP027131.fasta \
  -outfmt 7 \
  -out AP027078_AP027131.tblastx.out
```

## 3. Generate the pairwise plot

```bash
gbdraw linear \
  --gbk AP027078.gb AP027131.gb \
  -b AP027078_AP027131.tblastx.out \
  --align_center \
  --separate_strands \
  -o hepatoplasmataceae_pair \
  -f svg
```

![Two aligned Hepatoplasmataceae genomes connected by precomputed TBLASTX ribbons](../../examples/tutorial-2-pairwise-blast.svg)

The default `ribbon` style draws each retained BLAST high-scoring pair (HSP) as a straight filled ribbon between its query and subject coordinate spans. The `--pairwise_match_style curve` option bends the same spans, which can make dense tracks easier to follow; it does not merge or resize HSPs.

## 4. Run protein searches during diagram generation

gbdraw can instead generate protein matches directly from CDS-derived amino acid sequences. These commands omit `-b/--blast`.

All three modes compare CDS-derived proteins. Unless you pass an explicit executable path, gbdraw chooses the protein-search runtime in this order: bundled native LOSAT, `losat` on `PATH`, then NCBI BLAST+ `blastp` on `PATH`.

The package currently bundles LOSAT only for Linux x86_64. On macOS and Windows, install NCBI BLAST+ and make `blastp` available on `PATH`, or pass it with `--ncbi_blastp_bin /path/to/blastp`. Use `--losatp_bin /path/to/losat` to force a native LOSAT executable on any platform.

The NCBI BLAST+ fallback produces compatible 12-column outfmt 6 output, but its hit set can differ from LOSAT.

Pairwise protein ribbons:

```bash
gbdraw linear \
  --gbk AP027078.gb AP027131.gb \
  --protein_blastp_mode pairwise \
  --align_center \
  -o hepatoplasmataceae_protein_pairwise \
  -f svg
```

![Pairwise protein comparison between two Hepatoplasmataceae genomes](../../examples/tutorial-2-protein-pairwise.svg)

Ribbons between proteins in the same gbdraw orthogroup:

```bash
gbdraw linear \
  --gbk AP027078.gb AP027131.gb \
  --protein_blastp_mode orthogroup \
  --show_labels orthogroup_top \
  --pairwise_match_style curve \
  -o hepatoplasmataceae_orthogroup \
  -f svg
```

![Protein ribbons based on gbdraw similarity groups between two Hepatoplasmataceae genomes](../../examples/tutorial-2-protein-orthogroup.svg)

Collinear blocks:

```bash
gbdraw linear \
  --gbk AP027078.gb AP027131.gb \
  --protein_blastp_mode collinear \
  --collinear_min_anchors 2 \
  --collinear_color_mode orientation_identity \
  --pairwise_match_style curve \
  -o hepatoplasmataceae_collinear \
  -f svg
```

![Collinear protein blocks between two Hepatoplasmataceae genomes](../../examples/tutorial-2-protein-collinear.svg)

See [Tutorial 4](./4_Protein_Comparisons.md) for runtime selection, labels, and collinear tuning.

## 5. Interpreting orthogroups and collinear blocks

`pairwise` draws filtered protein matches between adjacent records. `orthogroup` assigns CDS-derived proteins to similarity-based groups and draws links between members of the same group. `collinear` combines compatible runs of those protein-match anchors into blocks.

### Orthogroups

In evolutionary genomics, an orthogroup is a set of genes descended from a single gene in the last common ancestor of the taxa being compared. Different genomes may contribute one gene, multiple genes, or no detectable member. An orthogroup is therefore not necessarily a set of one-to-one orthologs.

> [!NOTE]
> gbdraw builds similarity-based groups for visualization; it does not infer gene or species trees. These groups are not a substitute for phylogenetically informed orthology inference. Use a dedicated workflow such as OrthoFinder when orthology assignments are the analysis goal.

### Collinear blocks

An anchor links a pair of CDS-derived proteins assigned to the same gbdraw orthogroup. A collinear block is a run of anchors with compatible order in two records. For multi-anchor blocks, `plus` means that the anchors occur in the same order in both records; `minus` means that their order is reversed, consistent with an inversion between the displayed regions.

Multi-anchor blocks combine protein similarity with conserved local gene order. They can highlight conserved gene neighborhoods, including candidate operons or gene clusters, but they do not by themselves establish shared function or cotranscription. Genes lying between anchors are not automatically homologous or members of the same orthogroup.

The default `--collinear_min_anchors 1` retains singleton links. Set it to `2` to require at least two anchors per rendered block; higher values require more anchors.

Blocks are built from reciprocal best-hit (RBH) anchors. The `--collinear_color_mode orientation_identity` option assigns separate identity gradients to forward and inverted blocks. In the web app, `Evidence scope` controls which record pairs are searched for collinearity evidence. `Adjacent pairs` limits searches to neighboring records. `All records` searches every record pair for grouping and block support. Both settings render blocks only between adjacent displayed records.

## 6. Select records or regions

Linear mode supports three selectors:

- `--record_id`: choose one record by ID or `#index`; repeat the option for multiple input files
- `--reverse_complement`: set one Boolean orientation value per occurrence; repeat the option in input-file order
- `--region`: crop a region with `record_id:start-end[:rc]`; repeat the option for multiple regions

Example: select one record and crop a region.

```bash
gbdraw linear \
  --gbk AP027078.gb \
  --record_id AP027078.1 \
  --region AP027078.1:100000-250000 \
  -o AP027078_region \
  -f svg
```

Example: reverse-complement the second input.

```bash
gbdraw linear \
  --gbk AP027078.gb AP027131.gb \
  --reverse_complement false \
  --reverse_complement true \
  -o hepatoplasmataceae_second_rc \
  -f svg
```

Each occurrence of `--record_id` and `--reverse_complement` accepts one value. Repeat the option rather than placing several values after a single occurrence.

If you use `--region` or `--reverse_complement` with BLAST results, make sure the coordinates match the selected inputs.

For larger comparisons, put row-specific selectors and crops in a records table instead of repeating several order-sensitive CLI options. Create `records.tsv`:

```tsv
gbk	record_label	record_id	region	reverse_complement	order
AP027078.gb	T. litorale	AP027078.1	100000-250000	0	1
AP027131.gb	H. vulgare	AP027131.1	50000-180000	1	2
```

```bash
gbdraw linear \
  --records_table records.tsv \
  -o hepatoplasmataceae_regions \
  -f svg
```

![Two cropped Hepatoplasmataceae regions with the H. vulgare record reverse-complemented](../../examples/tutorial-2-record-selectors.svg)

This selector-only example omits `-b/--blast`; gbdraw does not remap full-record BLAST coordinates to cropped or reverse-complemented inputs.

`--records_table` is an alternative input source, so do not combine it with `--gbk`, `--gff`, or `--fasta`.

## 7. Compare more than two genomes

The Gallery comparison extends the pair above to five genomes. Download the remaining three GenBank and FASTA records:

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AP027133.1&rettype=gbwithparts&retmode=text" -O AP027133.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AP027133.1&rettype=fasta&retmode=text" -O AP027133.fasta

wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AP027132.1&rettype=gbwithparts&retmode=text" -O AP027132.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AP027132.1&rettype=fasta&retmode=text" -O AP027132.fasta

wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NZ_CP006932.1&rettype=gbwithparts&retmode=text" -O NZ_CP006932.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NZ_CP006932.1&rettype=fasta&retmode=text" -O NZ_CP006932.fasta
```

The first adjacent BLAST output file was created in Section 2. Generate the remaining three in the same order as the genome list:

```bash
tblastx -query AP027131.fasta -subject AP027133.fasta -outfmt 7 -out AP027131_AP027133.tblastx.out
tblastx -query AP027133.fasta -subject AP027132.fasta -outfmt 7 -out AP027133_AP027132.tblastx.out
tblastx -query AP027132.fasta -subject NZ_CP006932.fasta -outfmt 7 -out AP027132_NZ_CP006932.tblastx.out
```

```bash
gbdraw linear \
  --gbk AP027078.gb AP027131.gb AP027133.gb AP027132.gb NZ_CP006932.gb \
  -b AP027078_AP027131.tblastx.out AP027131_AP027133.tblastx.out AP027133_AP027132.tblastx.out AP027132_NZ_CP006932.tblastx.out \
  --align_center \
  --separate_strands \
  --show_gc \
  --show_skew \
  --palette default \
  -o hepatoplasmataceae_default \
  -f svg
```

Add `--evalue`, `--bitscore`, or `--alignment_length` to filter retained BLAST HSPs before drawing ribbons.

The output should have four ribbon bands, each connecting one adjacent pair in the five-genome input order.

![Five aligned Hepatoplasmataceae genomes connected by four adjacent BLAST ribbon tracks](../../examples/hepatoplasmataceae_default.svg)

> [!IMPORTANT]
> BLAST output files must follow the same order as the genome list.

## 8. Circular search-hit rings: LOSATN and TLOSATX

Circular mode can draw one ring of pairwise search hits per query sequence around an annotated reference. This example uses MjeNMV as the reference and the other nine majaniviruses as queries for LOSATN nucleotide searches and TLOSATX translated nucleotide searches. Comparing the diagrams shows where nucleotide similarity is sparse but translated-sequence similarity remains detectable.

### 8.1 Download the reference and queries

Download MjeNMV as the annotated GenBank reference:

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738868.1&rettype=gbwithparts&retmode=text" -O MjeNMV.gb
```

Download the nine comparison genomes as FASTA queries:

```bash
while read -r name accession; do
  wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=fasta&retmode=text" -O "${name}.fasta"
done <<'EOF'
MelaMJNV LC738874.1
PemoMJNVA LC738870.1
PeseMJNV LC738873.1
PemoMJNVB LC738871.1
LvMJNV LC738872.1
TrcuMJNV LC738879.1
MellatMJNV AP027153.1
MeenMJNV LC738876.1
MejoMJNV LC738878.1
EOF
```

The same GenBank and FASTA files are available under `examples/` in a source checkout.

### 8.2 Configure the shared circular plot

Open the [gbdraw web app](https://gbdraw.app/), select `Circular`, and upload `MjeNMV.gb` as the input genome. Use these settings for both comparisons:

| Setting | Value |
|---|---|
| Track Preset | Spreadout |
| Hide GC Content | Selected |
| Hide GC Skew | Selected |
| Legend Position | Left |
| Pairwise Comparisons | Run LOSAT |
| Ring Width | 8 |
| Ring Gap | 2 |
| Bitscore | 50 |
| E-value | `1e-5` |

In `Pairwise Comparisons`, click `Add Seq` and add the nine FASTA files in the order shown in the download block. Each FASTA is a query; the displayed MjeNMV genome is the subject and the coordinate reference for every ring.

### 8.3 Generate the nucleotide-level LOSATN rings

Select `LOSATN`, set `Task` to `blastn`, `Minimum Identity` to `50`, and `Minimum Length` to `100`. Click `Generate Diagram`.

![MjeNMV reference with nine majanivirus LOSATN rings; most nucleotide-level rings are sparse](../../examples/tutorial-2-majanivirus-losatn.svg)

Most rings contain only isolated nucleotide HSPs. MelaMJNV is the exception: its outer ring remains dense because it is more similar to MjeNMV at the nucleotide level than the other queries.

### 8.4 Generate the translated TLOSATX rings

Keep the same reference, query order, colors, and layout. Change `LOSAT Mode` to `TLOSATX`, then set `Minimum Identity` to `30` and `Minimum Length` to `30`. Keep `Reference gencode` and each query `Subject gencode` at `1`, then click `Generate Diagram` again.

![MjeNMV reference with nine majanivirus TLOSATX rings showing broader translated-sequence similarity](../../examples/tutorial-2-majanivirus-tlosatx.svg)

The eight queries with few retained nucleotide hits produce denser rings in the TLOSATX plot, while MelaMJNV remains dense in both. At these thresholds, translated-sequence similarity remains detectable in many regions that have no retained LOSATN HSP.

At the displayed thresholds, each of the eight divergent queries yields 108–177 LOSATN HSPs and 769–2,000 TLOSATX HSPs.

> [!NOTE]
> MjeNMV is annotated as linear, so gbdraw reports a topology warning. This tutorial intentionally uses a circular view to compare similarity rings; use linear mode when the displayed topology must match the record annotation.

The web app can save raw LOSAT results with `Save Raw LOSAT TSV`. To draw the same rings from the CLI, pass one saved outfmt 6/7 file per ring with `--conservation_blast` and select `--conservation_reference subject`. Rows with `start > end` on the selected reference side are drawn as reverse-orientation hits, not as circular wraparound hits.

[< Back to the Tutorials Index](./TUTORIALS.md)
[< Back to Tutorial 1](./1_Customizing_Plots.md) | [Go to Tutorial 3 >](./3_Advanced_Customization.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
