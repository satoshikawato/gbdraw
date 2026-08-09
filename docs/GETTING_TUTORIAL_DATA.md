[Documentation home](./DOCS.md) | [Start here](./TUTORIALS/README.md) | [Task guides](./HOW_TO/README.md) | [Gallery](./GALLERY.md)

# Get the tutorial inputs

Every gbdraw Tutorial names its input records, accessions, local filenames, and
supporting files. Download sequence records from the authoritative record
database named in the Tutorial. Do not substitute a sequence copy bundled in
this repository, and do not start a Tutorial from a prebuilt Gallery session.

Repository-hosted files are used only for gbdraw-specific support data such as
label rules, depth tables, and precomputed comparison results. A session may be
reloaded only after you create it in an earlier step of the same Tutorial.

## Download an NCBI or DDBJ sequence in a browser

Each sequence row links its versioned accession, for example
[`NC_012920.1`](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1). Keep the
accession version: `.1`, `.2`, and later suffixes identify a specific nucleotide
sequence version.

1. Open the accession link in the Tutorial.
2. On the NCBI record page, select **Send to**, then **File**.
3. Choose **GenBank (full)** for an annotated GenBank input, or **FASTA** when
   the Tutorial explicitly asks for FASTA.
4. Select **Create File**.
5. Rename the download to the exact local filename shown in the Tutorial and
   move it into that Tutorial's empty working directory.

For records whose primary entry is at DDBJ, the Tutorial still links a
versioned public record page and states which format and filename to use.

Database curators can update annotations without changing the nucleotide
sequence version. When a Tutorial supplies an official NCBI Revision History
download, use that exact link: it pins the feature table needed to reproduce
the documented figure. Do not replace it with the current record page or a
repository copy.

## Download an NCBI sequence from the command line

NCBI EFetch can return an accession-pinned GenBank flatfile directly. This
example saves the human mitochondrial record as the filename expected by the
Tutorials:

```bash
mkdir gbdraw-tutorial
cd gbdraw-tutorial

ncbi_efetch="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
curl -L --fail \
  "${ncbi_efetch}?db=nuccore&id=NC_012920.1&rettype=gbwithparts&retmode=text" \
  -o HmmtDNA.gbk
```

For a FASTA input, change `rettype=gbwithparts` to `rettype=fasta`. The
equivalent PowerShell download is:

```powershell
New-Item -ItemType Directory -Path gbdraw-tutorial
Set-Location gbdraw-tutorial

$ncbiEfetch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
$recordUri = "${ncbiEfetch}?db=nuccore&id=NC_012920.1&rettype=gbwithparts&retmode=text"
Invoke-WebRequest -Uri $recordUri -OutFile "HmmtDNA.gbk"
```

The CLI Tutorials provide complete EFetch commands for every NCBI record they
use, except where a Tutorial provides a pinned NCBI Revision History URL for an
annotation-sensitive record. Use that URL verbatim. MIBiG Tutorials instead
use accession-specific GenBank download URLs from the
[MIBiG repository](https://mibig.secondarymetabolites.org/repository/).

## Download repository support data

Some workflows need gbdraw-specific rule tables, quantitative tracks, or
precomputed comparison tables. These are not substitutes for the original
sequence records. Open the support-file link in the Tutorial, select
**Download raw file**, and save it with the exact filename shown.

On Windows, choose **All files** if the save dialog tries to append `.txt`.
Do not copy the rendered GitHub page into a text file.

If you plan to follow several Tutorials, you may clone the repository once to
obtain only those support files:

```bash
git clone --depth 1 https://github.com/satoshikawato/gbdraw.git gbdraw-tutorial-support
```

Support data live under
`gbdraw-tutorial-support/gbdraw/web/tutorial-data/`. Copy only the support
files listed by the current Tutorial into a new working directory. Continue to
download every GenBank or FASTA sequence from its accession-linked database.

## Check what you downloaded

After downloading `NC_012920.1` and one label rule for a human mitochondrial
Tutorial, the working directory should contain:

```text
gbdraw-tutorial/
├── HmmtDNA.gbk
└── cds_gene_qualifier_priority.tsv
```

On macOS, Linux, or WSL:

```bash
ls -l
head -n 1 HmmtDNA.gbk
grep -m 1 '^VERSION' HmmtDNA.gbk
head -n 1 cds_gene_qualifier_priority.tsv
```

On PowerShell:

```powershell
Get-ChildItem
Get-Content HmmtDNA.gbk -TotalCount 4
Get-Content cds_gene_qualifier_priority.tsv -TotalCount 1
```

The GenBank file must start with `LOCUS` and contain the requested versioned
accession, such as `VERSION     NC_012920.1`. It must not start with HTML such
as `<!DOCTYPE`. A successful download alone is not proof: each Tutorial also
names the expected output and shows a regenerated reference image.

## What the file labels mean

| Label used in the guides | Meaning |
| --- | --- |
| Authoritative download | A versioned sequence record downloaded from NCBI, DDBJ, MIBiG, or another database named by the Tutorial |
| Support download | A repository-maintained rule, track, annotation, or precomputed result used with the authoritative sequence |
| Create | A text, TSV, JSON, or Python file whose complete contents appear in the guide |
| Generated | An output or session created by following the current Tutorial |
| Reference result | A result used to compare with a new run; do not use it as an input unless the guide explicitly says to |

The machine-readable [tutorial-data manifest](../gbdraw/web/tutorial-data/manifest.json)
records the accession and checksum of the frozen copies used by automated
offline tests. Those internal copies make regression tests deterministic; they
are not the reader-facing sequence acquisition path.

[Documentation home](./DOCS.md) | [Start here](./TUTORIALS/README.md) | [Task guides](./HOW_TO/README.md) | [Gallery](./GALLERY.md)
