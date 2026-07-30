[Home](./DOCS.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | **Workflow guide** | [Python API](./PYTHON_API.md) | [Export](./EXPORT.md)

# Choose a gbdraw workflow

Use this page before collecting inputs. It separates the biological question from the interface and output format.

## Diagram layout

| Need | Choose | Why |
|---|---|---|
| One complete microbial, organelle, or viral genome | Circular | Shows features and quantitative tracks around a closed coordinate axis. |
| Regions, multiple records, or adjacent genome comparisons | Linear | Keeps record order, orientation, coordinate spans, and comparison links explicit. |
| Several replicons from one assembly on one canvas | Circular multi-record | Preserves circular layout while arranging chromosomes and plasmids together. |

Start with the [Beginner circular](https://gbdraw.app/gallery/#HmmtDNA_basic_circular) or [Beginner linear](https://gbdraw.app/gallery/#lambda_basic_linear) web tutorial.

## Annotation input

| Available input | Use | Check first |
|---|---|---|
| GenBank/GBFF | `--gbk` | Feature qualifiers, record IDs, and CDS translations are already bundled. |
| GFF3 plus FASTA | `--gff` with `--fasta` | Record IDs must match; CDS strand and phase must be valid. See [GFF3 + FASTA](./GFF3_FASTA.md). |
| Many records with row-specific settings | `--records_table` | Store paths, labels, selectors, crop, orientation, and order in one TSV. |

## Comparison evidence

| Question | Input or mode | What the drawing means |
|---|---|---|
| Where are retained local nucleotide or translated-nucleotide hits? | Precomputed BLAST outfmt 6/7 | Each ribbon/ring is an HSP coordinate span. |
| Which adjacent CDS translations have retained protein hits? | `--protein_blastp_mode pairwise` | Individual filtered protein-search matches. |
| Which CDS translations belong to the same search-derived group? | `--protein_blastp_mode orthogroup` | gbdraw similarity groups; not phylogeny-based orthogroups. |
| Which protein-match anchors occur in compatible local order? | `--protein_blastp_mode collinear` | Collinear blocks supported by ordered anchors. |

Use [precomputed BLAST comparisons](./TUTORIALS/2_Comparative_Genomics.md) when search results already exist. Use [protein comparisons](./TUTORIALS/4_Protein_Comparisons.md) when gbdraw should translate annotated CDS features and run LOSATP or BLASTP.

## Interface

| Interface | Best for |
|---|---|
| Hosted web app | Interactive exploration with no installation; genomic data stays in the browser. |
| Local `gbdraw gui` | The same browser workflow offline or with local assets. |
| CLI | Reproducible commands, batch work, and all static export formats. |
| `gbdraw` Python API | Python pipelines that need records and output in memory. See [Python API](./PYTHON_API.md). |

### Web UX profile v1

The browser interface has an explicit UX profile. These defaults are chosen for
interactive editing and intentionally differ from the CLI and Python defaults:

| Setting | Browser default |
|---|---|
| Feature strands | Separate strands on |
| Circular legend | Left |
| Linear legend | Bottom |
| One Circular record | Single diagram |
| Several records in one Circular input | Separate results; select **Multi-Record Canvas** for a grid |

A loaded session restores its saved values instead of applying this profile again.
The CLI keeps separate Circular outputs by default and uses
`--multi_record_canvas` for a grid. The package-root Python API automatically uses a
grid when `draw_circular()` receives several records; the typed API selects the
single- or multi-record builder explicitly.

Circular Web input accepts one GenBank/GBFF file or one GFF3 + FASTA pair per run. One file
may contain several records. Use the CLI or Python API when records must come from
several source files. Circular legends support Left, Right, Top, Bottom, all four
corner positions, and None.

The annotation editor covers common styles. Use an annotation TSV for the complete
style contract, including fields that are not shown in the form; see
[region annotation tables](./TUTORIALS/5_Table_Driven_Inputs.md#7-region-annotation-tables).

## Output and reproducibility

- Use SVG or PDF for editable vector output; use PNG only when a raster format is required.
- Use `interactive_svg` when feature or match popups are part of the deliverable.
- Save a `.gbdraw-session.json.gz` file from the web app when browser state, inputs, and post-generation edits must be restored. The CLI reads both forms, writes `.gbdraw-session.json` by default, and writes gzip when an explicit `--session_output` path ends in `.gz`.
- Store the command, input accessions/checksums, search settings, and gbdraw version with publication figures.

See [Export for publication](./EXPORT.md) and [interactive SVG/session workflows](./TUTORIALS/8_Interactive_SVG_Sessions.md).

[Home](./DOCS.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | **Workflow guide**
