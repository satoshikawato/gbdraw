[Home](../../DOCS.md) | [Tutorials](../README.md) | [CLI reference](../../REFERENCE/command-line.md) | [Export](../../REFERENCE/output-formats-and-export.md)

# Create a reproducible circular diagram from the command line

You will draw the 16,569 bp human mitochondrial reference genome as a standard
SVG. The finished figure includes 37 displayed features, coordinate ticks, GC
content, GC skew, and concise feature labels. CDS labels come from the `gene`
qualifier rather than the longer `product` text.

## What you'll need

- gbdraw installed so that `gbdraw -h` succeeds;
- the packaged [`HmmtDNA.gbk`](../../../gbdraw/web/tutorial-data/human-mitochondrion/HmmtDNA.gbk)
  fixture, saved with that filename;
- the packaged [`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv)
  label rule, saved with that filename;
- an empty working directory.

The [fixture manifest](../../../gbdraw/web/tutorial-data/manifest.json) records
the accession, sequence length, feature count, and SHA-256 checksum used by the
executable recipe.

## Step 1: Prepare the working directory

Create a new directory, place `HmmtDNA.gbk` and
`cds_gene_qualifier_priority.tsv` inside it, and enter the directory:

```bash
mkdir gbdraw-cli-circular
cd gbdraw-cli-circular
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

The command prints `Generated SVG: human_mitochondrion.svg` and writes the
diagram in the current directory.

## Step 3: Inspect the SVG

Open `human_mitochondrion.svg` in a browser or vector editor. Check the center
definition for `NC_012920.1` and `16,569 bp`, then follow the two inner plots for
GC content and GC skew. CDS labels should use short gene symbols such as `ND1`,
`COX1`, and `CYTB`, not product descriptions such as “NADH dehydrogenase
subunit 1.”

![Labeled circular human mitochondrial genome with gene-symbol CDS labels, coordinate ticks, GC content, GC skew, and a right legend](../../images/t-cli-01/human_mitochondrion.svg)

The committed recipe checks the XML structure, exact record metadata, 37 stable
feature IDs, both GC tracks, ticks, and the absence of scripts, event handlers,
and external links. It also checks that all 13 CDS gene symbols are present and
that the longer CDS product descriptions are absent from label text.

## If the command fails

- `gbdraw: command not found`: activate the environment where gbdraw is installed.
- `Output file already exists`: return to a new empty directory or choose a new
  output prefix. gbdraw does not overwrite an existing file by default.

## What you built

You created a labeled, editable standard SVG from a fixed public GenBank record.
Its qualifier-priority rule makes the label source explicit and reproducible. Run
`python docs/recipes/run_cli_scenarios.py --scenario T-CLI-01` from a repository
checkout to regenerate the published figure and repeat its semantic checks.
Use the [CLI reference](../../REFERENCE/command-line.md) for other Circular
options and the [export guide](../../REFERENCE/output-formats-and-export.md)
when a journal requires another format. For figure-review and regeneration
guidance, see [Prepare publication and reproducible figures](../../EXPLANATION/prepare-publication-and-reproducible-figures.md).
