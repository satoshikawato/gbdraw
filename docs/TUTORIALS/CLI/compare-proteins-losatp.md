[Documentation home](../../DOCS.md) | [All Tutorials](../README.md) | [Comparison reference](../../REFERENCE/comparison-programs-thresholds-and-results.md)

# Create protein Similarity groups with LOSATP from the command line

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web app](../GUI/compare-proteins-losatp.md) | **Command line** | [Python API](../PYTHON/compare-proteins-losatp.md) |

This recipe compares five complete aminoglycoside biosynthetic gene clusters,
builds the same 23 Similarity groups as the browser Tutorial, and aligns the
records to the group containing `CAG38695.1` (`og_1`).

## Step 1: Prepare the inputs

Copy the five BGC GenBank files, the two color tables, and
`cds_gene_qualifier_priority.tsv` into an empty directory.

## Step 2: Run LOSATP and draw the figure

<!-- executable:T-CLI-08:start -->
```bash
gbdraw linear \
  --gbk BGC0000708.gbk BGC0000709.gbk BGC0000711.gbk BGC0000712.gbk BGC0000713.gbk \
  --record_label '<i>Streptomyces lividus</i> CBS 844.73' \
  --record_label '<i>Streptomyces fradiae</i> ATCC 10745' \
  --record_label '<i>Streptomyces fradiae</i> MCIMB 8233' \
  --record_label '<i>Streptomyces rimosus</i> subsp. <i>paromomycinus</i> NRRL 2455' \
  --record_label '<i>Streptomyces ribosidificus</i> ATCC 21294' \
  --record_subtitle 'Lividomycin biosynthetic gene cluster' \
  --record_subtitle 'Neomycin biosynthetic gene cluster' \
  --record_subtitle 'Neomycin biosynthetic gene cluster' \
  --record_subtitle 'Paromomycin biosynthetic gene cluster' \
  --record_subtitle 'Ribostamycin biosynthetic gene' \
  --reverse_complement false \
  --reverse_complement false \
  --reverse_complement false \
  --reverse_complement false \
  --reverse_complement true \
  --protein_blastp_mode orthogroup \
  --losatp_threads 1 \
  --bitscore 50 \
  --evalue 0.01 \
  --identity 30 \
  --alignment_length 0 \
  --align_orthogroup_feature CAG38695.1 \
  --palette orange \
  --default_colors BGC0000708-BGC0000713_default_colors.tsv \
  --table BGC0000708-BGC0000713_specific_colors.tsv \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --show_labels first \
  --label_font_size 18 \
  --label_placement above_feature \
  --label_rotation 45 \
  --feature_height 75 \
  --block_stroke_color '#262626' \
  --block_stroke_width 2 \
  --line_stroke_width 2 \
  --axis_stroke_width 5 \
  --scale_style ruler \
  --track_layout middle \
  --keep_definition_left_aligned \
  --definition_line_style 'name:size=20,weight=bold' \
  --definition_line_style 'subtitle:size=20' \
  --definition_line_style 'accession:size=20,color=#7b7c7d' \
  --definition_line_style 'length:size=20,color=#7b7c7d' \
  --pairwise_match_style curve \
  --plot_title 'LOSATP Similarity groups across five whole BGC records' \
  --plot_title_position bottom \
  --legend bottom \
  -o bgc_losatp_groups \
  -f svg
```
<!-- executable:T-CLI-08:end -->

The bundled LOSAT runtime performs four adjacent searches and writes
`bgc_losatp_groups.svg`.

![Five BGC records aligned to Similarity group og_1](../../images/t-cli-08/bgc_losatp_groups.svg)

## Step 3: Verify the groups

The recipe pins the bundled LOSAT version, checks 232 qualified raw rows, 23
groups, 77 adjacent links, the reversed fifth record, and the final
presentation:

```bash
python docs/recipes/run_cli_scenarios.py --scenario T-CLI-08 --check
```

## What you built

You built the same `og_1`-aligned five-record comparison as the browser
Tutorial. The Similarity groups describe graph components in this dataset; they
are not a phylogenetic orthology claim.
