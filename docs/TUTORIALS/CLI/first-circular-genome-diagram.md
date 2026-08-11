[Home](../../DOCS.md) | [Tutorials](../README.md) | [Get tutorial files](../../GETTING_TUTORIAL_DATA.md) | [CLI reference](../../REFERENCE/command-line.md) | [Export](../../REFERENCE/output-formats-and-export.md)

# Create a reproducible circular diagram from the command line

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web-app workflow](../GUI/first-circular-genome-diagram.md) | **This page** | [Python workflow](../PYTHON/first-genome-diagram.md) |

All three workflows build the same labeled human mitochondrial figure. Only
the interface changes.

You will draw the 16,569 bp human mitochondrial reference genome as a standard
SVG. The finished figure includes 37 displayed features, coordinate ticks, GC
content, GC skew, and concise feature labels. CDS labels come from the `gene`
qualifier rather than the longer `product` text.

## What you'll need

| Kind | Filename | What to do |
| --- | --- | --- |
| Download | [`HmmtDNA.gbk`](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=gbwithparts&retmode=text) | Download NCBI accession [`NC_012920.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1) in full GenBank format and save it as `HmmtDNA.gbk`. |
| Download | [`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv) | Save this repository-hosted label rule as `cds_gene_qualifier_priority.tsv`. |
| Generated | `human_mitochondrion.svg` | The command writes this SVG. |
| Reference result | [`human_mitochondrion.svg`](../../images/t-cli-01/human_mitochondrion.svg) | Compare your generated SVG with this versioned result. |

This Tutorial has no Create files. You need gbdraw installed so that
`gbdraw -h` succeeds and an empty working directory.

Verify the sequence input against the versioned NCBI accession shown above.
The executable recipe separately pins its offline regression copy; that copy
is not a reader input.

## Step 1: Prepare the working directory

Create and enter a new directory:

```bash
mkdir gbdraw-cli-circular
cd gbdraw-cli-circular
```

The sequence link downloads accession `NC_012920.1` directly from NCBI in full
GenBank format. For the repository-hosted label rule, select **Download raw
file**. Save both files with the exact names in the table. See [Get the
tutorial files](../../GETTING_TUTORIAL_DATA.md) for browser, PowerShell, and
identity-check instructions.

On macOS, Linux, or WSL, download both files with `curl`:

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

Before running gbdraw, the working directory should contain:

```text
gbdraw-cli-circular/
├── HmmtDNA.gbk
└── cds_gene_qualifier_priority.tsv
```

## Step 2: Generate the first diagram

Run this command from the directory containing `HmmtDNA.gbk`:

<!-- executable:T-CLI-01:start -->
```bash
gbdraw circular \
  --gbk HmmtDNA.gbk \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --separate_strands \
  --track_type middle \
  --labels out \
  --species "<i>Homo sapiens</i>" \
  --legend right \
  -o human_mitochondrion \
  -f svg
```
<!-- executable:T-CLI-01:end -->

Expected output: the command prints `Generated SVG: human_mitochondrion.svg`
and writes the Generated file in the current directory.

## Step 3: Inspect the SVG

Open the Generated `human_mitochondrion.svg` in a browser or vector editor.
Check the center definition for `NC_012920.1` and `16,569 bp`, then follow the
two inner plots for GC content and GC skew. CDS labels should use short gene
symbols such as `ND1`, `COX1`, and `CYTB`, not product descriptions such as
"NADH dehydrogenase subunit 1."

![Labeled circular human mitochondrial genome with gene-symbol CDS labels, coordinate ticks, GC content, GC skew, and a right legend](../../images/t-cli-01/human_mitochondrion.svg)

The image above is the Reference result. Your SVG should have the same record,
track order, labels, and legend even if metadata or XML formatting differs.

The committed recipe checks the XML structure, exact record metadata, 37 stable
feature IDs, both GC tracks, ticks, and the absence of scripts, event handlers,
and external links. It also checks that all 13 CDS gene symbols are present and
that the longer CDS product descriptions are absent from label text.

## If the command fails

- `gbdraw: command not found`: activate the environment where gbdraw is installed.
- `Output file already exists`: return to a new empty directory or choose a new
  output prefix. gbdraw does not overwrite an existing file by default.
