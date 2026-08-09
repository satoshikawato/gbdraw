[Documentation home](../../DOCS.md) | [Tutorials](../README.md) | [CLI tutorials](README.md) | [Get tutorial files](../../GETTING_TUTORIAL_DATA.md) | [CLI reference](../../REFERENCE/command-line.md)

# Highlight mitochondrial features without editing the GenBank file

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web-app workflow](../GUI/highlight-mitochondrial-features.md) | **This page** | [Python workflow](../PYTHON/highlight-mitochondrial-features.md) |

This project starts with the complete 16,569 bp human mitochondrial reference
and turns it into a focused presentation without changing its GenBank
annotations.

## What you will build

You will make a baseline map, then a second SVG in which selected CDS and rRNA
features have deliberate colors and labels. All 13 mitochondrial CDS remain
visible. The D-loop is added as a named, origin-spanning region annotation,
using the same bracket semantics as the chloroplast Gallery example. Feature
shapes and strokes are presentation settings, not edits to the biological
record.

## Files used in this Tutorial

| Kind | Filename | What to do |
| --- | --- | --- |
| Download | [`HmmtDNA.gbk`](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=gbwithparts&retmode=text) | Download NCBI accession [`NC_012920.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1) in full GenBank format and save it as `HmmtDNA.gbk`. |
| Download | [`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv) | Save this repository-hosted label rule as `cds_gene_qualifier_priority.tsv`. |
| Create | `tables/presentation_colors.tsv` | Create this color-rule table under Step 1. |
| Create | `tables/presentation_labels.tsv` | Create this label whitelist under Step 1. |
| Create | `tables/presentation_label_overrides.tsv` | Create this label-replacement table under Step 1. |
| Create | `tables/mitochondrial_regions.tsv` | Create this D-loop annotation table under Step 1. |
| Generated | `mitochondrial_features_baseline.svg` | The first command writes this baseline SVG. |
| Generated | `mitochondrial_features_highlighted.svg` | The second command writes the finished SVG. |
| Reference result | [`mitochondrial_features_baseline.svg`](../../images/t-cli-03/mitochondrial_features_baseline.svg) | Compare the Generated baseline with this versioned result. |
| Reference result | [`mitochondrial_features_highlighted.svg`](../../images/t-cli-03/mitochondrial_features_highlighted.svg) | Compare the Generated final SVG with this versioned result. |

Install gbdraw so that `gbdraw circular -h` succeeds.

## Step 1: Prepare the source record and presentation tables

### Create the working directory and download the source inputs

Create the project directory and its `tables` subdirectory:

```bash
mkdir gbdraw-cli-mitochondrial-features
cd gbdraw-cli-mitochondrial-features
mkdir tables
```

The sequence link downloads accession `NC_012920.1` directly from NCBI in full
GenBank format. For the repository-hosted label rule, select **Download raw
file**. Save both files with the exact names in the table. See [Get the
tutorial files](../../GETTING_TUTORIAL_DATA.md) for browser, PowerShell, and
identity-check instructions.

On macOS, Linux, or WSL, run:

```bash
ncbi_efetch="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
gbdraw_data_base="https://raw.githubusercontent.com/satoshikawato/gbdraw/main/gbdraw/web/tutorial-data"
curl -L "${ncbi_efetch}?db=nuccore&id=NC_012920.1&rettype=gbwithparts&retmode=text" -o HmmtDNA.gbk
curl -L "$gbdraw_data_base/shared/cds_gene_qualifier_priority.tsv" -o cds_gene_qualifier_priority.tsv
```

Confirm that the source record reports `VERSION     NC_012920.1`:

```bash
grep '^VERSION' HmmtDNA.gbk
```

### Define color rules

Save this as `tables/presentation_colors.tsv`. Rows are evaluated in order;
the first matching feature-specific rule wins over the base palette.

```tsv
CDS	gene	^ND(4L|[1-6])$	#3B82F6	NADH dehydrogenase
CDS	gene	^COX[1-3]$	#EF4444	Cytochrome c oxidase
CDS	gene	^ATP[68]$	#F59E0B	ATP synthase
CDS	gene	^CYTB$	#8B5CF6	Cytochrome b
rRNA	gene	^RNR[12]$	#10B981	Ribosomal RNA
```

### Select and rename labels

Save this whitelist as `tables/presentation_labels.tsv`:

```tsv
CDS	gene	^(ND[1-6]|ND4L|COX[1-3]|ATP[68]|CYTB)$
rRNA	gene	^RNR[12]$
```

Then save these resolved-label replacements as
`tables/presentation_label_overrides.tsv`:

```tsv
record_id	feature_type	qualifier	value	label_text
NC_012920.1	rRNA	label	^s-rRNA$	12S rRNA
NC_012920.1	rRNA	label	^l-rRNA$	16S rRNA
```

The qualifier-priority file selects `gene` for all 13 CDS labels. The whitelist
keeps those CDS and the two rRNAs in scope, and the override table renames the
two rRNA labels.

### Add the D-loop as a region annotation

Save this as `tables/mitochondrial_regions.tsv`:

```tsv
set_id	id	mark	record	start	end	coordinate_space	wraps_origin	label	lane	stroke	stroke_width	line_cap	label_color	label_font_size	label_orientation	label_offset
mitochondrial_regions	d_loop	bracket	NC_012920.1	16024	576	source	true	D-loop	0	#202020	3	tick	#202020	14	tangent	7
```

The source D-loop joins bases 16,024–16,569 and 1–576. One row with
`wraps_origin=true` preserves that biology and draws a single named bracket.
The finished command draws every CDS and tRNA as an arrow and rRNA as a
rectangle.

### Check the working directory

After downloading the source record and support file and creating the four
tables, the working directory should contain:

```text
gbdraw-cli-mitochondrial-features/
├── HmmtDNA.gbk
├── cds_gene_qualifier_priority.tsv
└── tables/
    ├── presentation_colors.tsv
    ├── presentation_labels.tsv
    ├── presentation_label_overrides.tsv
    └── mitochondrial_regions.tsv
```

## Step 2: Run both reproducible commands

Run the block from the directory containing the two inputs and the `tables`
directory.

<!-- executable:T-CLI-03:start -->
```bash
gbdraw circular \
  --gbk HmmtDNA.gbk \
  --labels none \
  --legend right \
  -o mitochondrial_features_baseline \
  -f svg

gbdraw circular \
  --gbk HmmtDNA.gbk \
  -k CDS,rRNA,tRNA \
  --table tables/presentation_colors.tsv \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --label_whitelist tables/presentation_labels.tsv \
  --label_table tables/presentation_label_overrides.tsv \
  --annotation_table tables/mitochondrial_regions.tsv \
  --feature_shape CDS=arrow \
  --feature_shape rRNA=rectangle \
  --feature_shape tRNA=arrow \
  --arrow_head_length_ratio auto \
  --arrow_shaft_width_ratio 0.72 \
  --track_type middle \
  --circular_track_slot 'ticks:ticks@side=outside,tick_label_layout=label_out_tick_in' \
  --circular_track_slot 'features:features@side=overlay,lane_direction=split' \
  --circular_track_slot 'mitochondrial_regions:annotations@set_id=mitochondrial_regions,side=inside,w=24px,show_labels=true,padding_px=1,overflow=compress' \
  --labels both \
  --label_rendering auto \
  --block_stroke_color '#1F2937' \
  --block_stroke_width 1.5 \
  --axis_stroke_color '#374151' \
  --axis_stroke_width 4 \
  --line_stroke_color '#9CA3AF' \
  --line_stroke_width 1.5 \
  --species '<i>Homo sapiens</i>' \
  --plot_title 'Human mitochondrial feature presentation' \
  --plot_title_position top \
  --legend right \
  -o mitochondrial_features_highlighted \
  -f svg
```
<!-- executable:T-CLI-03:end -->

## Step 3: Verify the two outputs

Expected output: the first command writes the Generated
`mitochondrial_features_baseline.svg`, and the second writes the Generated
`mitochondrial_features_highlighted.svg`.

Open `mitochondrial_features_baseline.svg` first. Verify that its definition
names `NC_012920.1`, reports `16,569 bp`, and shows the complete circular
record.

![Baseline human mitochondrial map before presentation overrides](../../images/t-cli-03/mitochondrial_features_baseline.svg)

Open `mitochondrial_features_highlighted.svg`. Check that all 13 CDS labels use
their gene names, including `COX1`, and that `12S rRNA` and `16S rRNA` are also
present. The blue, red, amber, violet, and green rules should appear in the
legend. rRNA blocks are rectangular, directional features remain arrows, and
the D-loop appears as one origin-spanning inner bracket.

![Human mitochondrial map with all CDS, explicit feature colors, and a D-loop region bracket](../../images/t-cli-03/mitochondrial_features_highlighted.svg)

Both images above are Reference results. Compare their record definition,
labels, feature shapes, annotation bracket, colors, and legend with your two
Generated SVGs.

## Check the source record

Only presentation inputs changed. `HmmtDNA.gbk` stayed byte-for-byte identical;
the SVG still contains the same complete `NC_012920.1` sequence context.

## Next steps

- [Set colors, labels, visibility, shapes, and strokes](../../HOW_TO/CLI/set-colors-labels-visibility-shapes-and-strokes.md)
- [Review feature-rule and label schemas](../../REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md)
