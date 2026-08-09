[Home](../../DOCS.md) | [Tutorials](../README.md) | [Get tutorial files](../../GETTING_TUTORIAL_DATA.md) | [Session compatibility](../../REFERENCE/session-and-request-compatibility.md)

# Create an interactive figure and reproduce it from a CLI session

## Choose how to build this figure

| GUI | CLI | Python API |
| --- | --- | --- |
| [Web-app workflow](../GUI/create-and-resume-an-interactive-figure.md) | **This page** | [Python workflow](../PYTHON/create-and-resume-an-interactive-figure.md) |

Create the same human mitochondrial figure as the web tutorial, export its
self-contained Interactive SVG, save the request as a compressed session, and
replay the session under a new output name.

## Files used in this Tutorial

| Kind | Filename | What to do |
| --- | --- | --- |
| Download | [`HmmtDNA.gbk`](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=gbwithparts&retmode=text) | Download NCBI accession [`NC_012920.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1) in full GenBank format and save it as `HmmtDNA.gbk`. |
| Download | [`cds_gene_qualifier_priority.tsv`](../../../gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv) | Save this repository-hosted label rule as `cds_gene_qualifier_priority.tsv`. |
| Generated | `interactive_human_mitochondrion.svg` | The first command writes the static master. |
| Generated | `interactive_human_mitochondrion.interactive.svg` | The first command writes the offline Interactive SVG. |
| Generated | `interactive_handoff.gbdraw-session.json.gz` | The first command writes the compressed session used by the second command. |
| Generated | `restored_interactive_figure.svg` | The second command replays the session into this static SVG. |
| Reference result | [`restored_interactive_figure.svg`](../../images/t-cli-11/restored_interactive_figure.svg) | Compare the restored Generated SVG with this versioned result. |

This Tutorial has no Create files. You need gbdraw installed and an empty
working directory.

## Step 1: Prepare the inputs, then export the figure and session

### Prepare the working directory

Create and enter an empty directory:

```bash
mkdir gbdraw-cli-interactive-session
cd gbdraw-cli-interactive-session
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

The working directory should now contain:

```text
gbdraw-cli-interactive-session/
├── HmmtDNA.gbk
└── cds_gene_qualifier_priority.tsv
```

### Export the figure and session

<!-- executable:T-CLI-11:start -->
```bash
gbdraw circular \
  --gbk HmmtDNA.gbk \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --separate_strands \
  --track_type middle \
  --labels out \
  --species '<i>Homo sapiens</i>' \
  --legend right \
  --session_output interactive_handoff.gbdraw-session.json.gz \
  -o interactive_human_mitochondrion \
  -f svg,interactive_svg

gbdraw circular \
  --session interactive_handoff.gbdraw-session.json.gz \
  -o restored_interactive_figure \
  -f svg
```
<!-- executable:T-CLI-11:end -->

Expected output: the first command writes
`interactive_human_mitochondrion.svg`,
`interactive_human_mitochondrion.interactive.svg`, and
`interactive_handoff.gbdraw-session.json.gz`. The second command reads only
the Generated session and writes `restored_interactive_figure.svg`.

## Step 2: Verify the restored figure

Open the Generated `restored_interactive_figure.svg`. It should match the
original static master exactly: the same 37 feature IDs, the same `COX1` search
metadata, and the same visible content, reproduced entirely from the saved
session.

![Human mitochondrial map restored from the CLI session](../../images/t-cli-11/restored_interactive_figure.svg)

The image above is the Reference result. Use it to verify the record identity,
visible labels and tracks, and legend before testing the Interactive SVG in a
browser.

## Next steps

- [Review session and request compatibility](../../REFERENCE/session-and-request-compatibility.md)
- [Review output formats and export](../../REFERENCE/output-formats-and-export.md)
