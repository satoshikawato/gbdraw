[Home](../../DOCS.md) | [Tutorials](../README.md) | [CLI reference](../../REFERENCE/command-line.md) | [Export](../../REFERENCE/output-formats-and-export.md)

# Create a reproducible linear diagram from the command line

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web-app workflow](../GUI/first-linear-genome-diagram.md) | **This page** | Not yet migrated |

The implemented variants build the same complete Lambda map. A Python variant
must reproduce that figure rather than introduce a different example.

You will draw the complete 48,502 bp Lambda reference genome with concise gene
labels and a ruler on the record axis. The result is one standard SVG containing
all 73 displayed CDS features.

## What you'll need

- gbdraw installed so that `gbdraw -h` succeeds;
- the packaged [`NC_001416.gb`](../../../gbdraw/web/tutorial-data/lambda/NC_001416.gb)
  fixture;
- the packaged [`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv)
  label rule;
- an empty working directory.

The label rule prefers short `gene` values over long product descriptions. The
[fixture manifest](../../../gbdraw/web/tutorial-data/manifest.json) records both
files and their checksums.

## Step 1: Prepare the working directory

Create a new directory, place both files inside it, and enter the directory:

```bash
mkdir gbdraw-cli-linear
cd gbdraw-cli-linear
```

## Step 2: Generate the first diagram

Run this command from the directory containing the two fixture files:

<!-- executable:T-CLI-02:start -->
```bash
gbdraw linear \
  --gbk NC_001416.gb \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --show_labels all \
  --scale_style ruler \
  --track_layout above \
  --ruler_on_axis \
  -o lambda_linear \
  -f svg
```
<!-- executable:T-CLI-02:end -->

The command prints `Generated SVG: lambda_linear.svg` and writes the diagram in
the current directory.

## Step 3: Inspect the SVG

Open `lambda_linear.svg`. The definition at the left identifies `NC_001416.1`
and `48,502 bp`. The axis is marked every 5 kbp, and short labels such as `A`,
`B`, `J`, and `int` remain legible near their CDS features.

![Linear Lambda genome with concise gene labels and a ruler](../../images/t-cli-02/lambda_linear.svg)

The committed recipe checks the XML structure, record metadata, all 73 stable
feature IDs, representative gene labels, and ruler labels. It also rejects
active or external content in this standard SVG.

## If the command fails

- `--ruler_on_axis is ignored`: keep `--scale_style ruler` and
  `--track_layout above` in the same command.
- `Output file already exists`: use a new empty directory or a new output
  prefix. gbdraw refuses to replace the existing SVG by default.

## What you built

You created a whole-record Lambda map with compact biological labels and an
axis ruler. Run `python docs/recipes/run_cli_scenarios.py --scenario T-CLI-02`
from a repository checkout to regenerate the published figure and repeat its
semantic checks. Use the [CLI reference](../../REFERENCE/command-line.md) for
other Linear layout options. See [Choose circular or linear](../../EXPLANATION/choose-circular-or-linear.md)
when deciding which layout suits another record.
