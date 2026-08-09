[Documentation home](../../DOCS.md) | [All Tutorials](../README.md) | [Get tutorial files](../../GETTING_TUTORIAL_DATA.md) | [Command-line reference](../../REFERENCE/command-line.md)

# Compare Lambda and DE3 from the command line

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web app](../GUI/compare-genomes-losatn.md) | **Command line** | [Python API](../PYTHON/compare-genomes-losatn.md) |

This variant uses the six-row LOSATN result produced by the browser Tutorial.
It draws both complete phage records in the same order and keeps the 120 px
comparison band used by the GUI.

## Files used in this Tutorial

| Kind | Filename | What to do |
| --- | --- | --- |
| Download | [`NC_001416.gb`](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001416.1&rettype=gbwithparts&retmode=text) | Download NCBI accession [`NC_001416.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1) in full GenBank format and save it as `NC_001416.gb`. |
| Download | [`NC_042057.1.gb`](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_042057.1&rettype=gbwithparts&retmode=text) | Download NCBI accession [`NC_042057.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_042057.1) in full GenBank format and save it as `NC_042057.1.gb`. |
| Download | [`lambda-de3.losatn.tsv`](../../../gbdraw/web/tutorial-data/lambda-de3-comparison/lambda-de3.losatn.tsv) | Save this repository-hosted LOSATN result as `lambda-de3.losatn.tsv`. |
| Generated | `lambda-de3-losatn.svg` | The command writes this comparison SVG. |
| Reference result | [`lambda-de3-losatn.svg`](../../images/t-cli-07/lambda-de3-losatn.svg) | Compare your generated SVG with this versioned result. |

This Tutorial has no Create files.

## Step 1: Prepare the working directory

Create and enter an empty directory:

```bash
mkdir gbdraw-cli-losatn
cd gbdraw-cli-losatn
```

The two sequence links download accession-pinned full GenBank records directly
from NCBI. The LOSATN link is a repository-hosted support table; select
**Download raw file** for it. Save every file with the exact name in the table.
See [Get the tutorial files](../../GETTING_TUTORIAL_DATA.md) for browser,
PowerShell, and identity-check instructions.

On macOS, Linux, or WSL, run:

```bash
ncbi_efetch="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
gbdraw_data_base="https://raw.githubusercontent.com/satoshikawato/gbdraw/main/gbdraw/web/tutorial-data"
curl -L "${ncbi_efetch}?db=nuccore&id=NC_001416.1&rettype=gbwithparts&retmode=text" -o NC_001416.gb
curl -L "${ncbi_efetch}?db=nuccore&id=NC_042057.1&rettype=gbwithparts&retmode=text" -o NC_042057.1.gb
curl -L "$gbdraw_data_base/lambda-de3-comparison/lambda-de3.losatn.tsv" -o lambda-de3.losatn.tsv
```

Confirm that the two records report the pinned versions:

```bash
grep -H '^VERSION' NC_001416.gb NC_042057.1.gb
```

Expect `NC_001416.1` and `NC_042057.1`. The working directory should now
contain:

```text
gbdraw-cli-losatn/
├── NC_001416.gb
├── NC_042057.1.gb
└── lambda-de3.losatn.tsv
```

## Step 2: Draw the finished comparison

<!-- executable:T-CLI-07:start -->
```bash
gbdraw linear \
  --gbk NC_001416.gb NC_042057.1.gb \
  --record_id NC_001416.1 \
  --record_id NC_042057.1 \
  --blast lambda-de3.losatn.tsv \
  --bitscore 50 \
  --evalue 0.01 \
  --identity 0 \
  --alignment_length 0 \
  --comparison_height 120 \
  -o lambda-de3-losatn \
  -f svg
```
<!-- executable:T-CLI-07:end -->

Expected output: gbdraw writes the Generated `lambda-de3-losatn.svg` in the
working directory.

Open that SVG. It should contain six ribbons between the complete Lambda and
DE3 records.

![Lambda and DE3 LOSATN comparison from the CLI](../../images/t-cli-07/lambda-de3-losatn.svg)

The image above is the Reference result. Verify that both accessions and all
six endpoint pairs match the TSV and that the ribbons use the 120 px comparison
band set above.

## Next steps

- [Draw comparisons from prepared evidence](../../HOW_TO/CLI/draw-precomputed-comparisons.md)
- [Review comparison thresholds and result semantics](../../REFERENCE/comparison-programs-thresholds-and-results.md)
