[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [How-to guides](../HOW_TO/README.md) | [Reference](../REFERENCE/README.md)

# Choose a genome-comparison method

Pick the evidence type before choosing colors or link shapes. Each method answers a different question.

| Method | Evidence | Use it when | Result shown |
|---|---|---|---|
| Uploaded BLAST table | Existing outfmt 6 or 7 rows | The search was run elsewhere or must be frozen | One retained HSP per link or ring span |
| LOSATN | Nucleotide local alignments | Closely related nucleotide sequences should be compared directly | Pairwise nucleotide spans |
| TLOSATX | Translated nucleotide alignments | Coding similarity may remain after nucleotide divergence | Translated query and subject spans |
| LOSATP Pairwise | Annotated CDS translations | Individual protein matches between records matter | Filtered protein-match links |
| LOSATP Similarity groups | Search-derived protein groups across records | Shared membership matters more than local order | Group membership and group links |
| LOSATP Collinear | Compatible ordered protein anchors | Conserved local order or rearrangement is the focus | Blocks supported by ordered anchors |
| Circular similarity rings | Prepared or generated similarity evidence around one reference | Coverage around a Circular reference is the intended view | One ordered ring per evidence set |

Similarity groups are search-derived connected groups. They are not phylogenetic orthogroups. Collinear blocks depend on retained anchors, gap limits, and order; they are not whole-genome synteny calls. A link means that a row passed the configured filters, not that two complete genes or genomes are equivalent.

LOSATN and TLOSATX run in the web app. The CLI reads prepared nucleotide or translated-nucleotide evidence through `--blast`, `--comparisons_table`, or `--conservation_blast`. The CLI can run Pairwise, Similarity-group, and Collinear protein searches through `--protein_blastp_mode`.

Exact programs, thresholds, coordinate conventions, and result fields are in [Comparison programs, thresholds, and result semantics](../REFERENCE/comparison-programs-thresholds-and-results.md).
