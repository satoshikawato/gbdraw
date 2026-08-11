[Home](../../DOCS.md) | [Tutorials](../README.md) | [Get tutorial files](../../GETTING_TUTORIAL_DATA.md) | [CLI reference](../../REFERENCE/command-line.md)

# Recreate the Gallery chloroplast map from the command line

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web-app workflow](../GUI/build-an-annotated-chloroplast-map.md) | **This page** | [Python workflow](../PYTHON/build-an-annotated-chloroplast-map.md) |

All three workflows produce the same complete tobacco plastome map. Only the
interface changes.

This command-line workflow reproduces the Interactive SVG Gallery chloroplast
map with functional gene colors, radial labels on both strands, one inner
LSC/IRb/SSC/IRa bracket lane, GC content, and an upper-left legend.

## Files used in this Tutorial

| Kind | Filename | What to do |
| --- | --- | --- |
| Download | [`NC_001879.gbk`](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001879.2&rettype=gbwithparts&retmode=text) | Download NCBI accession [`NC_001879.2`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001879.2) in full GenBank format and save it as `NC_001879.gbk`. |
| Download | [`nicotiana-tabacum-regions.tsv`](../../../gbdraw/web/tutorial-data/tobacco-plastome-regions/nicotiana-tabacum-regions.tsv) | Save this repository-hosted region table with the exact filename. |
| Download | [`chloroplast_specific_table.tsv`](../../../gbdraw/web/tutorial-data/tobacco-plastome-regions/chloroplast_specific_table.tsv) | Save this repository-hosted feature color table with the exact filename. |
| Download | [`qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/tobacco-plastome-regions/qualifier_priority.tsv) | Save this repository-hosted label rule with the exact filename. |
| Generated | `cli_annotated_chloroplast.svg` | The command writes this SVG. |
| Reference result | [`cli_annotated_chloroplast.svg`](../../images/t-cli-06/cli_annotated_chloroplast.svg) | Compare your generated SVG with this versioned result. |

This Tutorial has no Create files.

## Step 1: Prepare the working directory

Create and enter an empty directory:

```bash
mkdir gbdraw-cli-chloroplast
cd gbdraw-cli-chloroplast
```

The sequence link downloads accession `NC_001879.2` directly from NCBI in full
GenBank format. The other links are repository-hosted support tables; select
**Download raw file** for those. Save every file with the exact name in the
table. See [Get the tutorial files](../../GETTING_TUTORIAL_DATA.md) for
browser, PowerShell, and identity-check instructions.

On macOS, Linux, or WSL, run:

```bash
ncbi_efetch="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
gbdraw_data_base="https://raw.githubusercontent.com/satoshikawato/gbdraw/main/gbdraw/web/tutorial-data"
curl -L "${ncbi_efetch}?db=nuccore&id=NC_001879.2&rettype=gbwithparts&retmode=text" -o NC_001879.gbk
curl -L "$gbdraw_data_base/tobacco-plastome-regions/nicotiana-tabacum-regions.tsv" -o nicotiana-tabacum-regions.tsv
curl -L "$gbdraw_data_base/tobacco-plastome-regions/chloroplast_specific_table.tsv" -o chloroplast_specific_table.tsv
curl -L "$gbdraw_data_base/tobacco-plastome-regions/qualifier_priority.tsv" -o qualifier_priority.tsv
```

Confirm that the source record reports `VERSION     NC_001879.2`:

```bash
grep '^VERSION' NC_001879.gbk
```

The working directory should now contain:

```text
gbdraw-cli-chloroplast/
├── NC_001879.gbk
├── nicotiana-tabacum-regions.tsv
├── chloroplast_specific_table.tsv
└── qualifier_priority.tsv
```

## Step 2: Run the documented command

<!-- executable:T-CLI-06:start -->
```bash
gbdraw circular \
  --gbk NC_001879.gbk \
  -t chloroplast_specific_table.tsv \
  -k CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,rep_origin \
  --species '<i>Nicotiana tabacum</i>' \
  --track_type tuckin \
  --separate_strands \
  --gc \
  --no-skew \
  --labels both \
  --label_placement radial \
  --outer_label_x_radius_offset 0.9 \
  --outer_label_y_radius_offset 0.9 \
  --inner_label_x_radius_offset 0.975 \
  --inner_label_y_radius_offset 0.975 \
  --qualifier_priority qualifier_priority.tsv \
  --annotation_table nicotiana-tabacum-regions.tsv \
  --circular_track_slot 'features:features@side=overlay,lane_direction=split' \
  --circular_track_slot 'plastome_regions:annotations@set_id=plastome_regions,side=inside,r=0.65,w=20px,inner_gap_px=1,outer_gap_px=1,show_labels=true,padding_px=1,overflow=compress' \
  --circular_track_slot 'gc_content:dinucleotide_content@side=inside,r=0.56,w=0.08,nt=GC,legend_label=GC content' \
  --block_stroke_color black \
  --block_stroke_width 1 \
  --line_stroke_width 2 \
  --axis_stroke_width 3 \
  --definition_font_size 28 \
  --legend upper_left \
  -o cli_annotated_chloroplast \
  -f svg
```
<!-- executable:T-CLI-06:end -->

Expected output: gbdraw writes the Generated
`cli_annotated_chloroplast.svg` in the working directory.

The explicit slots are the same three-slot stack used by the GUI and Python
workflows. Because no tick or skew slot is declared, neither appears in the
finished figure.

## Step 3: Inspect the result

Open `cli_annotated_chloroplast.svg`. Confirm `NC_001879.2`, 147 logical
features, radial labels, all four structural-region brackets, GC content, and
the functional-color legend.

![Gallery-style tobacco plastome drawn from the command line](../../images/t-cli-06/cli_annotated_chloroplast.svg)

The image above is the Reference result. Your SVG should have the same complete
record, structural-region brackets, track order, labels, and functional colors.

The recipe runner executes the literal command above in a clean directory and
checks the complete record, feature count, annotation identities, three-slot
order, labels, colors, and parity with the Python API rendering.

## Next steps

- [Annotation table fields](../../REFERENCE/input-formats-and-tsv-schemas.md#annotation-table-fields)
- [Review track and annotation schemas](../../REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md)
