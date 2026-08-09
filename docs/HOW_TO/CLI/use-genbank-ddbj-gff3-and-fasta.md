[Documentation home](../../DOCS.md) | [CLI how-to guides](README.md) | [Input formats](../../REFERENCE/input-formats-and-tsv-schemas.md) | [CLI reference](../../REFERENCE/command-line.md)

# How to use GenBank, DDBJ, GFF3, and FASTA inputs

Use either an annotated flat file or a matched GFF3 and FASTA pair. The two
commands below draw the same complete 48,502 bp Lambda reference record, so you
can check that changing the input representation does not change the biological
record being shown.

## Prerequisites

- Install gbdraw so that `gbdraw -h` succeeds.
- Start in an empty working directory.
- Download accession `NC_001416.1` from NCBI EFetch twice:
  - save the
    [complete GenBank flat file](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001416.1&rettype=gbwithparts&retmode=text)
    as `NC_001416.gb`;
  - save the
    [complete FASTA sequence](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001416.1&rettype=fasta&retmode=text)
    as `NC_001416.fna`.
- Copy the repository support annotation
  [`NC_001416.gff3`](../../../gbdraw/web/tutorial-data/lambda-gff3/NC_001416.gff3)
  into the same directory.

See [Get the tutorial files](../../GETTING_TUTORIAL_DATA.md) for the
authoritative-download and accession-check procedure.

All three files represent the whole `NC_001416.1` sequence. The GFF3 and FASTA
files were derived from the same versioned record without cropping or
splitting it.

## Draw the GenBank or DDBJ flat file

Files downloaded from GenBank or DDBJ in GenBank flat-file format use `--gbk`.
Run the first command below to create `lambda_genbank.svg`.

## Draw the matched GFF3 and FASTA pair

GFF3 contains annotations and coordinates; FASTA supplies the nucleotide
sequence. Pass the files together with `--gff` and `--fasta`. Each GFF3 sequence
ID must exactly match a FASTA header ID.

<!-- executable:H-CLI-01:start -->
```bash
gbdraw linear \
  --gbk NC_001416.gb \
  --show_labels none \
  --scale_style ruler \
  --track_layout above \
  --ruler_on_axis \
  -o lambda_genbank \
  -f svg

gbdraw linear \
  --gff NC_001416.gff3 \
  --fasta NC_001416.fna \
  --show_labels none \
  --scale_style ruler \
  --track_layout above \
  --ruler_on_axis \
  -o lambda_gff3 \
  -f svg
```
<!-- executable:H-CLI-01:end -->

Both diagrams identify `NC_001416.1`, show its complete 48,502 bp coordinate
range, and contain the same 73 CDS features on their original strands.

![Complete Lambda record loaded from a GenBank flat file](../../images/h-cli-01/lambda_genbank.svg)

![The same complete Lambda record loaded from matched GFF3 and FASTA files](../../images/h-cli-01/lambda_gff3.svg)

## Diagnose a GFF3 and FASTA ID mismatch

Inspect the first column of the GFF3 and the identifier immediately after `>` in
each FASTA header. For this fixture, both are `NC_001416.1`. If the FASTA header
uses another ID, gbdraw stops without writing a diagram and reports:

```text
No matching FASTA record found for GFF record NC_001416.1.
```

The executable recipe creates a temporary mismatched header, checks this failure
message, and confirms that no SVG is written. That negative input is used only
inside the clean test directory; it is not a public tutorial fixture.

## Use more than one biological record

Supply multiple complete files after `--gbk`, or one matched FASTA file for each
GFF3 file. Do not divide one genome into artificial contigs just to demonstrate
a multi-record layout. The [multi-record CLI layout guide](arrange-multiple-circular-records-and-tracks.md)
uses separate complete biological records and covers placement choices.

Run `python docs/recipes/run_cli_scenarios.py --scenario H-CLI-01` from a
repository checkout to regenerate both figures and repeat the success and
failure checks. See the [input-format reference](../../REFERENCE/input-formats-and-tsv-schemas.md)
for supported flat-file variants and GFF3 identity rules.
