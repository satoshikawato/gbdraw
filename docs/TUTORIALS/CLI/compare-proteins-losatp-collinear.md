[Documentation home](../../DOCS.md) | [All Tutorials](../README.md) | [Get tutorial files](../../GETTING_TUTORIAL_DATA.md) | [Comparison reference](../../REFERENCE/comparison-programs-thresholds-and-results.md)

# Reproduce the Hepatoplasmataceae Collinear map from the command line

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web app](../GUI/compare-proteins-losatp-collinear.md) | **Command line** | [Python API](../PYTHON/compare-proteins-losatp-collinear.md) |

Run an adjacent LOSATP search across five complete genomes and draw
the same Collinear blocks used by the GUI Tutorial.

## Files used in this Tutorial

| Kind | Filename | What to do |
| --- | --- | --- |
| Download | [`AP027078.gb`](https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=AP027078.1&sat=3&satkey=69902295) | Download the official NCBI Revision History snapshot of `AP027078.1` and save it as `AP027078.gb`. |
| Download | [`AP027131.gb`](https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=AP027131.1&sat=3&satkey=69902296) | Download the official NCBI Revision History snapshot of `AP027131.1` and save it as `AP027131.gb`. |
| Download | [`AP027133.gb`](https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=AP027133.1&sat=3&satkey=69902298) | Download the official NCBI Revision History snapshot of `AP027133.1` and save it as `AP027133.gb`. |
| Download | [`AP027132.gb`](https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=AP027132.1&sat=3&satkey=69902297) | Download the official NCBI Revision History snapshot of `AP027132.1` and save it as `AP027132.gb`. |
| Download | [`NZ_CP006932.gb`](https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=NZ_CP006932.1&sat=60&satkey=39275474) | Download the official NCBI Revision History snapshot of `NZ_CP006932.1` in full GenBank format and save it as `NZ_CP006932.gb`. |
| Generated | `losatp_collinear.svg` | The command writes this SVG after running LOSATP. |
| Reference result | [`losatp_collinear.svg`](../../images/t-cli-10/losatp_collinear.svg) | Compare your generated SVG with this versioned result. |

This Tutorial has no Create files.

Use the five Revision History links, not the live accession downloads. NCBI can
update a record's annotation without changing its sequence accession version.
These pinned revisions preserve the exact feature tables used to reproduce this
Tutorial's 2,994 displayed features and 500 Collinear matches.

## Step 1: Prepare the working directory

Create and enter an empty directory:

```bash
mkdir gbdraw-cli-losatp-collinear
cd gbdraw-cli-losatp-collinear
```

Each sequence link downloads an annotation revision directly from NCBI in full
GenBank format. Save the records with the exact names and in the order shown in
the table. See [Get the tutorial
files](../../GETTING_TUTORIAL_DATA.md) for browser, PowerShell, and
identity-check instructions.

On macOS, Linux, or WSL, run:

```bash
curl -L "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=AP027078.1&sat=3&satkey=69902295" -o AP027078.gb
curl -L "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=AP027131.1&sat=3&satkey=69902296" -o AP027131.gb
curl -L "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=AP027133.1&sat=3&satkey=69902298" -o AP027133.gb
curl -L "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=AP027132.1&sat=3&satkey=69902297" -o AP027132.gb
curl -L "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&db=nuccore&report=gbwithparts&id=NZ_CP006932.1&sat=60&satkey=39275474" -o NZ_CP006932.gb
```

Confirm that the records report the five pinned versions:

```bash
grep -H '^VERSION' AP027078.gb AP027131.gb AP027133.gb AP027132.gb NZ_CP006932.gb
```

The working directory should now contain:

```text
gbdraw-cli-losatp-collinear/
├── AP027078.gb
├── AP027131.gb
├── AP027133.gb
├── AP027132.gb
└── NZ_CP006932.gb
```

## Step 2: Render a standard SVG

<!-- executable:T-CLI-10:start -->
```bash
gbdraw linear \
  --gbk AP027078.gb AP027131.gb AP027133.gb AP027132.gb NZ_CP006932.gb \
  --protein_blastp_mode collinear \
  --losatp_threads 32 \
  --collinear_search_scope adjacent \
  --collinear_max_unit_gap 0 \
  --collinear_min_anchors 1 \
  --collinear_max_diagonal_drift 0 \
  --collinear_max_conflicts_in_merge_gap 1 \
  --collinear_color_mode orientation_identity \
  --bitscore 50 \
  --evalue 0.01 \
  --identity 0 \
  --alignment_length 0 \
  --pairwise_match_style curve \
  --track_layout middle \
  --align_center \
  --separate_strands \
  --gc \
  --skew \
  --scale_style ruler \
  --palette ajisai \
  --plot_title 'LOSATP Collinear blocks across Hepatoplasmataceae' \
  --plot_title_position top \
  --legend right \
  -o losatp_collinear \
  -f svg
```
<!-- executable:T-CLI-10:end -->

Expected output: LOSATP performs four adjacent protein searches,
and gbdraw writes the Generated `losatp_collinear.svg` in the working directory.

Open `losatp_collinear.svg` and compare its record layout and Collinear blocks
with the image below.

![Five Hepatoplasmataceae genomes with adjacent Collinear blocks](../../images/t-cli-10/losatp_collinear.svg)

The image above is the Reference result. Verify five complete records in the
documented order, centered alignment, rulers, GC content and skew, and 500
rendered Collinear match elements colored by orientation and identity.

## Next steps

- [Draw Collinear blocks for another CLI project](../../HOW_TO/CLI/draw-collinear-protein-blocks.md)
- [Choose a genome-comparison method](../../EXPLANATION/choose-a-genome-comparison-method.md)
