[Documentation home](../DOCS.md) | [LOSATN Tutorial](../TUTORIALS/GUI/compare-genomes-losatn.md) | [Protein-comparison Tutorial](../TUTORIALS/GUI/compare-proteins-losatp.md) | [Choose a method](../EXPLANATION/choose-a-genome-comparison-method.md)

# Comparison programs, thresholds, and result semantics

## Capability matrix

| Evidence or mode | Web app | CLI search | Prepared input in CLI/API |
|---|---:|---:|---:|
| BLAST outfmt 6/7 | Upload | No | Yes |
| LOSATN nucleotide search | Yes | No | Yes, after exporting the result table |
| TLOSATX translated-nucleotide search | Yes | No | Yes, after exporting the result table |
| LOSATP Pairwise proteins | Yes | Yes | Yes |
| LOSATP Similarity groups | Yes | Yes | Session/current artifacts |
| LOSATP Collinear blocks | Yes | Yes | Session/current artifacts |
| Circular similarity rings | Yes | Uses prepared conservation evidence | Yes |

## Filters

Search results can be filtered by e-value, bit score, identity, and alignment length. Protein modes may add per-query candidate limits, reciprocal or grouping rules, and Collinear anchor/gap constraints. A documented result count is valid only for the recorded program version, inputs and order, task, thresholds, scheduling, and thread count.

Prepared BLAST-style rows retain query, subject, identity, alignment length, mismatch and gap counts, query span, subject span, e-value, and bit score. Commented outfmt 7 headers are accepted where the reader supports them. Query and subject direction are meaningful and must match the displayed record mapping.

## Frozen Lambda–DE3 example

The nucleotide-comparison guides use two complete RefSeq records in a fixed order:

1. query: Lambda `NC_001416.1`, 48,502 bp;
2. subject: Enterobacteria phage DE3 `NC_042057.1`, 42,925 bp.

They do not use cropped records or artificial contigs. Pinned LOSAT 0.1.0 browser
WASM runs use serial scheduling and one thread. LOSATN `megablast` produces six
raw rows, all six of which pass the documented default display filters. TLOSATX
produces 397 raw rows; 266 pass the default filters, while the How-to guide uses
a 1,000 bp minimum alignment length to show seven links. The exact arguments,
coordinate ranges, orientations, byte sizes, and checksums are recorded in the
[fixture manifest](../../gbdraw/web/tutorial-data/manifest.json).

## Result meanings

- Pairwise links represent retained individual search rows.
- Similarity-group membership is derived from retained protein-search relationships. It is not a phylogenetic orthology call.
- A Collinear block is a compatible ordered run of retained protein anchors under the selected unit-gap and diagonal-drift rules.
- A Circular ring paints retained spans against the selected reference side; ring order follows the evidence order.

An empty drawing can mean that the search returned no rows, that rows failed thresholds, that endpoints did not map to displayed records or features, or that CDS translations were unusable. Inspect exported raw evidence and warnings before changing thresholds.

## Cache identity

A reusable comparison result binds sequence content, selected protein set, record and feature identity, query/subject direction, program, and meaningful search arguments. Display labels, filenames, and modification times do not substitute for content identity. Export raw evidence for a durable record; a cache is not the provenance artifact.
