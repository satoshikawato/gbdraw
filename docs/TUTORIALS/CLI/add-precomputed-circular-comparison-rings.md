[Documentation home](../../DOCS.md) | [All Tutorials](../README.md) | [Comparison reference](../../REFERENCE/comparison-programs-thresholds-and-results.md)

# Add precomputed circular comparison rings from the command line

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web app](../GUI/add-precomputed-circular-comparison-rings.md) | **Command line** | [Python API](../PYTHON/add-precomputed-circular-comparison-rings.md) |

Use three frozen TLOSATX tables to place Danio, Drosophila, and
Caenorhabditis rings inside the complete human mitochondrial reference.

## Step 1: Prepare the inputs

Copy `HmmtDNA.gbk`, the three `.tlosatx.tsv` files, their three companion
`.fna` files, and `cds_gene_qualifier_priority.tsv` into an empty directory.

## Step 2: Draw the three rings

<!-- executable:T-CLI-09:start -->
```bash
gbdraw circular \
  --gbk HmmtDNA.gbk \
  --conservation_blast danio-human.tlosatx.tsv drosophila-human.tlosatx.tsv caenorhabditis-human.tlosatx.tsv \
  --conservation_fasta NC_002333.2.fna NC_024511.2.fna NC_001328.1.fna \
  --conservation_reference subject \
  --conservation_labels 'Danio rerio (NC_002333.2)' 'Drosophila melanogaster (NC_024511.2)' 'Caenorhabditis elegans (NC_001328.1)' \
  --conservation_colors '#4E79A7' '#F28E2B' '#59A14F' \
  --bitscore 50 \
  --evalue 1e-5 \
  --identity 40 \
  --alignment_length 50 \
  --conservation_ring_width 18 \
  --conservation_ring_gap 4 \
  --species '<i>Homo sapiens</i>' \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --track_type middle \
  --labels out \
  --definition_font_size 18 \
  --plot_title 'Precomputed TLOSATX rings around Homo sapiens mtDNA' \
  --plot_title_position top \
  --legend right \
  -o precomputed_circular_rings \
  -f svg
```
<!-- executable:T-CLI-09:end -->

![Human mitochondrial reference with three TLOSATX rings](../../images/t-cli-09/precomputed_circular_rings.svg)

## Step 3: Verify the result

```bash
python docs/recipes/run_cli_scenarios.py --scenario T-CLI-09 --check
```

The validator checks the subject-reference direction, ring order, labels,
comparison FASTA identities, 106 retained HSPs, and standard-SVG safety.

## What you built

You reproduced the browser ring map from fixed evidence without rerunning a
sequence search. Every ring position is a coordinate on `NC_012920.1`.
