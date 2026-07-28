[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)

[< Back to the guide index](./TUTORIALS.md)
[< Previous: Set feature colors and labels](./3_Advanced_Customization.md) | [Next: Use TSV manifests >](./5_Table_Driven_Inputs.md)

# Draw protein matches from annotated CDS features

Compare annotated CDS translations as individual pairwise hits, gbdraw similarity groups across records, or locally ordered collinear blocks. Each example builds a linear diagram by translating annotated CDS features before running the protein search. The runtime is configurable, and the similarity groups are search-derived visualization groups rather than phylogeny-based orthogroups.

## 1. Prepare annotated GenBank inputs

The protein-search modes need two or more annotated GenBank or GFF3 + FASTA records with CDS translations, or CDS features that can be translated.

This tutorial uses three majanivirus GenBank records:

```bash
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738868.1&rettype=gbwithparts&retmode=text" -O MjeNMV.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738874.1&rettype=gbwithparts&retmode=text" -O MelaMJNV.gb
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=LC738870.1&rettype=gbwithparts&retmode=text" -O PemoMJNVA.gb
```

If you are working from a source checkout, the same files are also available under `examples/`.

## 2. Runtime selection

Unless you pass an explicit executable path, `--protein_blastp_mode pairwise`, `orthogroup`, and `collinear` choose the protein search program in this order:

1. bundled native LOSAT on Linux x86_64 when available
2. `losat` on `PATH`
3. NCBI BLAST+ `blastp` on `PATH`

Use an explicit path when you need to control the runtime:

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

Pass only one of `--losatp_bin` and `--ncbi_blastp_bin` in a command.

## 3. Pairwise protein matches

`pairwise` searches each adjacent record pair and draws a ribbon for every retained protein match.

```bash
gbdraw linear \
  --gbk MjeNMV.gb MelaMJNV.gb \
  --protein_blastp_mode pairwise \
  --align_center \
  --pairwise_match_style curve \
  -o tutorial-protein-pairwise \
  -f svg
```

This writes `tutorial-protein-pairwise.svg`. The curved ribbons connect CDS-derived protein hits between the two adjacent records.

![Pairwise majanivirus protein comparison with curved ribbons between two linear records](../../examples/tutorial-protein-pairwise.svg)

## 4. gbdraw similarity-group ribbons (`orthogroup` mode)

`orthogroup` assigns CDS-derived proteins to gbdraw similarity groups across all input records before drawing ribbons between adjacent records. The mode does not perform phylogeny-based orthology inference.

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

This writes `majani_orthogroup.svg`.

`--show_labels orthogroup_top` labels the topmost displayed member of each gbdraw similarity group, which is useful when the same group appears in multiple records.

![Protein ribbons based on gbdraw similarity groups across three majanivirus records](../../examples/majani_orthogroup.svg)

## 5. Collinear blocks

`collinear` combines protein-match anchors that occur in compatible local order.

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

This writes `majani_collinear.svg`.

`--collinear_min_anchors 2` removes singleton blocks. `--collinear_color_mode orientation_identity` separates forward and inverted blocks while still encoding identity.

In the web preview or standalone interactive SVG, click a collinear block to export its complete query and subject block envelopes as nucleotide FASTA. The envelopes may include intergenic sequence and genes that are not anchors. This remains distinct from the orthogroup-member nucleotide and amino-acid actions in the same popup.

For a multi-record layout, add `--collinear_search_scope all` to search every record pair. gbdraw renders the accepted blocks only between adjacent layout rows and omits same-row ribbons. This supports all-vs-all comparisons between two groups of replicons, including chromosome-to-noncorresponding-chromosome blocks. See the [two-strain *Vibrio* example](./7_Linear_Layout.md#two-strains-with-multiple-replicons).

![Collinear protein blocks across three majanivirus records](../../examples/majani_collinear.svg)

## 6. When to prefer precomputed `-b/--blast`

Use precomputed BLAST tables when you need to preserve an existing result, use custom database settings, compare nucleotide or translated nucleotide sequences, or draw hits that were filtered by an upstream workflow.

Do not combine `-b/--blast` with `--protein_blastp_mode`. The CLI rejects that combination because the two options define different comparison sources.

For Python workflows with multi-record rows, use `LinearComparisonOptions(protein_mode="pairwise", pairs=((0, 2), (1, 3)))` to run only the declared record pairs. Pair indices are zero-based and must connect adjacent layout rows. Omitting `pairs` preserves adjacent-record behavior. See the [Python API linear example](../PYTHON_API.md#linear-diagrams-and-comparisons).

## 7. When a saved protein search is reused

Current sessions identify each CDS from its complete location, strand, and stable record-instance binding. The protein identity manifest keeps that full machine identity and its display metadata once. Generated protein FASTA, the raw LOSAT QUERY and SUBJECT fields, protein maps, and derived comparison references use deterministic session-global runtime handles with the form `h_[a-z2-7]{26}`. These compact IDs avoid repeating a long feature hash and alias in every hit row.

Changing an upload filename, file modification time, resource name, display alias, or saving and loading the same biological inputs does not invalidate a compatible raw protein-search result. A display-only change can rebuild derived display metadata without rerunning LOSAT. Changing an amino-acid sequence, the selected protein set, record-instance binding, feature location or strand, or a search argument does invalidate the affected result. Query/subject direction is significant.

Session versions 36 and 37 keep current protein raw results in cache schema 4,
derived protein comparisons in schema 3, and the protein identity manifest in
schema 2. Nucleotide LOSAT results remain schema 2, so protein and nucleotide
entries with different schemas may appear in one valid session.

Sessions 27–33 retain the verified migration path for schema-2 protein candidates and derived schema-1 evidence. On import, gbdraw keeps those artifacts separate from current hits. Generation verifies the complete FASTA content, program and arguments, direction, and one-to-one feature mapping before writing a schema-4 copy and rebuilding derived schema 3 without rerunning LOSAT. A candidate that cannot be verified is ignored only for that pair, and LOSAT runs normally. Session versions 34 and 35 were branch-internal development formats and are not supported.

The internal handles do not appear in a normal **Save Raw LOSAT TSV** download. Immediately before download, gbdraw uses the manifest to replace only QUERY and SUBJECT with the readable, percent-encoded form of a `protein_id`, `locus_tag`, GFF `ID`, or location fallback. Duplicate aliases within one record instance receive a deterministic short ordinal. The hydrator preserves comments, row order, columns 3–12, numeric spelling, and line endings, and it aborts the whole download if any handle cannot be resolved. User-uploaded comparison TSV is left unchanged.

[< Back to the guide index](./TUTORIALS.md)
[< Previous: Set feature colors and labels](./3_Advanced_Customization.md) | [Next: Use TSV manifests >](./5_Table_Driven_Inputs.md)

[Home](../DOCS.md) | [Installation](../INSTALL.md) | [Quickstart](../QUICKSTART.md) | [Tutorials](./TUTORIALS.md) | [Recipes](../RECIPES.md) | [CLI Reference](../CLI_Reference.md) | [Gallery](../GALLERY.md) | [FAQ](../FAQ.md) | [About](../ABOUT.md)
