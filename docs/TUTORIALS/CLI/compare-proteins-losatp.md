[Documentation home](../../DOCS.md) | [All Tutorials](../README.md) | [Get tutorial files](../../GETTING_TUTORIAL_DATA.md) | [Comparison reference](../../REFERENCE/comparison-programs-thresholds-and-results.md)

# Create protein Similarity groups with LOSATP from the command line

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web app](../GUI/compare-proteins-losatp.md) | **Command line** | [Python API](../PYTHON/compare-proteins-losatp.md) |

This recipe compares five complete aminoglycoside biosynthetic gene clusters,
builds the same 23 Similarity groups as the browser Tutorial, and aligns the
records to the group containing `CAG38695.1` (`og_1`).

## Files used in this Tutorial

| Kind | Filename | What to do |
| --- | --- | --- |
| Download | [`BGC0000708.gbk`](https://mibig.secondarymetabolites.org/repository/BGC0000708.5/BGC0000708.gbk) | Download MIBiG entry `BGC0000708.5` and save it as `BGC0000708.gbk`. |
| Download | [`BGC0000709.gbk`](https://mibig.secondarymetabolites.org/repository/BGC0000709.5/BGC0000709.gbk) | Download MIBiG entry `BGC0000709.5` and save it as `BGC0000709.gbk`. |
| Download | [`BGC0000711.gbk`](https://mibig.secondarymetabolites.org/repository/BGC0000711.5/BGC0000711.gbk) | Download MIBiG entry `BGC0000711.5` and save it as `BGC0000711.gbk`. |
| Download | [`BGC0000712.gbk`](https://mibig.secondarymetabolites.org/repository/BGC0000712.5/BGC0000712.gbk) | Download MIBiG entry `BGC0000712.5` and save it as `BGC0000712.gbk`. |
| Download | [`BGC0000713.gbk`](https://mibig.secondarymetabolites.org/repository/BGC0000713.5/BGC0000713.gbk) | Download MIBiG entry `BGC0000713.5` and save it as `BGC0000713.gbk`. |
| Download | [`BGC0000708-BGC0000713_default_colors.tsv`](../../../gbdraw/web/tutorial-data/aminoglycoside-bgc-five/BGC0000708-BGC0000713_default_colors.tsv) | Save this repository-hosted base color table with the exact filename. |
| Download | [`BGC0000708-BGC0000713_specific_colors.tsv`](../../../gbdraw/web/tutorial-data/aminoglycoside-bgc-five/BGC0000708-BGC0000713_specific_colors.tsv) | Save this repository-hosted feature color table with the exact filename. |
| Download | [`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv) | Save this repository-hosted label rule as `cds_gene_qualifier_priority.tsv`. |
| Generated | `bgc_losatp_groups.svg` | The command writes this SVG after running LOSATP. |
| Reference result | [`bgc_losatp_groups.svg`](../../images/t-cli-08/bgc_losatp_groups.svg) | Compare your generated SVG with this versioned result. |

This Tutorial has no Create files.

## Step 1: Prepare the working directory

Create and enter an empty directory:

```bash
mkdir gbdraw-cli-losatp-groups
cd gbdraw-cli-losatp-groups
```

The five sequence links download the versioned BGC entries directly from
MIBiG. The other links are repository-hosted support tables; select **Download
raw file** for those. Save every file with the exact name in the table. See
[Get the tutorial files](../../GETTING_TUTORIAL_DATA.md) for browser,
PowerShell, and identity-check instructions.

On macOS, Linux, or WSL, run:

```bash
gbdraw_data_base="https://raw.githubusercontent.com/satoshikawato/gbdraw/main/gbdraw/web/tutorial-data"
curl -L "https://mibig.secondarymetabolites.org/repository/BGC0000708.5/BGC0000708.gbk" -o BGC0000708.gbk
curl -L "https://mibig.secondarymetabolites.org/repository/BGC0000709.5/BGC0000709.gbk" -o BGC0000709.gbk
curl -L "https://mibig.secondarymetabolites.org/repository/BGC0000711.5/BGC0000711.gbk" -o BGC0000711.gbk
curl -L "https://mibig.secondarymetabolites.org/repository/BGC0000712.5/BGC0000712.gbk" -o BGC0000712.gbk
curl -L "https://mibig.secondarymetabolites.org/repository/BGC0000713.5/BGC0000713.gbk" -o BGC0000713.gbk
curl -L "$gbdraw_data_base/aminoglycoside-bgc-five/BGC0000708-BGC0000713_default_colors.tsv" -o BGC0000708-BGC0000713_default_colors.tsv
curl -L "$gbdraw_data_base/aminoglycoside-bgc-five/BGC0000708-BGC0000713_specific_colors.tsv" -o BGC0000708-BGC0000713_specific_colors.tsv
curl -L "$gbdraw_data_base/shared/cds_gene_qualifier_priority.tsv" -o cds_gene_qualifier_priority.tsv
```

Confirm that the five downloaded records report the expected BGC accessions:

```bash
grep -H '^VERSION' BGC0000708.gbk BGC0000709.gbk BGC0000711.gbk BGC0000712.gbk BGC0000713.gbk
```

The working directory should now contain:

```text
gbdraw-cli-losatp-groups/
├── BGC0000708.gbk
├── BGC0000709.gbk
├── BGC0000711.gbk
├── BGC0000712.gbk
├── BGC0000713.gbk
├── BGC0000708-BGC0000713_default_colors.tsv
├── BGC0000708-BGC0000713_specific_colors.tsv
└── cds_gene_qualifier_priority.tsv
```

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

Expected output: the bundled LOSAT runtime performs four adjacent searches and
writes the Generated `bgc_losatp_groups.svg` in the working directory.

Open `bgc_losatp_groups.svg` and compare its record order and link layout with
the image below.

![Five BGC records aligned to Similarity group og_1](../../images/t-cli-08/bgc_losatp_groups.svg)

The image above is the Reference result. Verify 232 raw rows, 23 Similarity
groups, and 77 adjacent links. The fifth record should remain reversed, matching
the browser Tutorial's alignment.

## Next steps

- [Create protein Similarity groups for another CLI project](../../HOW_TO/CLI/create-protein-similarity-groups.md)
- [Choose a genome-comparison method](../../EXPLANATION/choose-a-genome-comparison-method.md)
