[Documentation home](../../DOCS.md) | [All Tutorials](../README.md) | [Comparison reference](../../REFERENCE/comparison-programs-thresholds-and-results.md)

# Reproduce the Hepatoplasmataceae Collinear map from the command line

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web app](../GUI/compare-proteins-losatp-collinear.md) | **Command line** | [Python API](../PYTHON/compare-proteins-losatp-collinear.md) |

Run the qualified adjacent LOSATP search across five complete genomes and draw
the same Collinear blocks used by the GUI Tutorial.

## Step 1: Prepare the records

Copy `AP027078.gb`, `AP027131.gb`, `AP027133.gb`, `AP027132.gb`, and
`NZ_CP006932.gb` into an empty working directory.

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

![Five Hepatoplasmataceae genomes with adjacent Collinear blocks](../../images/t-cli-10/losatp_collinear.svg)

## Step 3: Verify the replay

```bash
python docs/recipes/run_cli_scenarios.py --scenario T-CLI-10 --check
```

The validator checks all five record IDs and lengths, the adjacent four-pair
scope, 500 displayed Collinear matches, curved ribbons,
orientation-and-identity colors, GC content, GC skew, ruler, title, and legend.

## What you built

You reproduced the Gallery's adjacent Collinear project from the five source
records with the Gallery project's fixed 32-thread search setting.
